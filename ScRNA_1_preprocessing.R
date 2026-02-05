rm(list=ls())


setwd(this.path::this.dir())
getwd()


library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(harmony)
library(SingleR)
library(clustree)
library(RColorBrewer)
library(ggsci)
library(cowplot)
library(DoubletFinder) 
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(GSVA)
library(ggprism)
library(limma)
library(ggthemes)






############################################################################################################################################################
################################################################     1. Data Propressing 
############################################################################################################################################################
SCdata<-list()
dir_file<-dir("../seekonetools_readcount")    ###  51 samples

########################    step 1.  Create Seurat Object  &  merge
for(i in 1:length(dir_file)){
  counts <- Read10X(paste0("../seekonetools_readcount/",dir_file[i]))
  SCdata[[i]]<-CreateSeuratObject(counts = counts,project = dir_file[i],min.cells = 0,min.features = 0)
  SCdata[[i]]<-RenameCells(SCdata[[i]],add.cell.id=dir_file[i])
  SCdata[[i]][["percent.mt"]] <- PercentageFeatureSet(SCdata[[i]], pattern = "^MT-")
}
datasets <-as.list(dir_file)
SCdata <- merge(SCdata[[1]], y = SCdata[2:length(SCdata)], add.cell.ids = datasets, project = "allsc")
dim(SCdata);table(SCdata$orig.ident)      ##  33538 414198
a = data.frame(table(SCdata$orig.ident))
colnames(a) <- c('sampleID','origincellNum')


########################    step 2.  QC
selected_c <- WhichCells(SCdata, expression = nFeature_RNA >= 500 & percent.mt < 10 & nFeature_RNA <= 10000)
selected_f <- rownames(SCdata)[Matrix::rowSums(SCdata@assays$RNA@counts > 0 ) >= 0.0001*dim(SCdata)[2]]  ###  41.4198 cell/gene 
SCdata <- subset(SCdata, features = selected_f, cells = selected_c)
table(SCdata$orig.ident);dim(SCdata)   #  21879 genes * 149516   cells 
a1 = data.frame(table(SCdata$orig.ident))
a$Filtersamplenum <- sapply(as.character(a$sampleID),function(x) a1$Freq[which(a1$Var1 == x)])
rm(a1)


########################    step 3.  Check Double
SCdata@meta.data$doubletFinder <- 0;SCdata@meta.data$doubletFinder1 <- 0
sampleID <- unique(SCdata@meta.data$orig.ident)
dim_usage=30;RunPCAdim <- NULL
for(i in 1:length(sampleID)){
  choose_data <- subset(SCdata,orig.ident == sampleID[i])
  choose_data <- NormalizeData(choose_data, normalization.method = "LogNormalize", scale.factor = 10000)
  all.genes <- rownames(choose_data)
  choose_data <- ScaleData(object = choose_data,features = all.genes)
  choose_data <- FindVariableFeatures(object = choose_data,nfeatures = 2000)
  choose_data <- RunPCA(choose_data, features = VariableFeatures(object = choose_data))
  data_use<-Stdev(object = choose_data,reduction = 'pca')
  RunPCAdim[[i]] = c(sum(data_use[1:40]),sum(data_use[1:30]),sum(data_use[1:20]),sum(data_use[1:10]))
  names(RunPCAdim)[i] = sampleID[i]
  saveRDS(RunPCAdim,"./RunPCAdim.rds")
  choose_data <- FindNeighbors(choose_data, dims = 1:dim_usage)
  choose_data <- FindClusters(choose_data, resolution = 0.4) 
  choose_data <- RunUMAP(choose_data, dims = 1:dim_usage, label = T)
  #########     Check Double1
  ## pK Identification 
  sweep.res.list_kidney <- paramSweep_v3(choose_data, PCs = 1:dim_usage, sct = FALSE)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  mpK<-as.numeric(as.vector(bcmvn_kidney$pK[which.max(bcmvn_kidney$BCmetric)]))
  ## Homotypic Doublet Proportion Estimate 
  homotypic.prop <- modelHomotypic(choose_data@meta.data$seurat_clusters)
  DoubletRate =  ncol(choose_data)*8*1e-6 
  nExp_poi <- round(DoubletRate*nrow(choose_data@meta.data))  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  ## Run DoubletFinder with varying classification stringencies 
  choose_data <- doubletFinder_v3(choose_data, PCs = 1:dim_usage, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  choose_data@meta.data$doubletFinder <- choose_data@meta.data[which(grepl('DF.classifications',colnames(choose_data@meta.data)))][,1]
  a = choose_data@meta.data[which(grepl('DF.classifications',colnames(choose_data@meta.data)))]
  for(j in 1:dim(choose_data)[2]){
    SCdata@meta.data$doubletFinder[which(colnames(SCdata@assays$RNA@counts) == rownames(a)[j])] <- a[j,1]
  }
  table(SCdata@meta.data$doubletFinder)
  choose_data$doubletFinder1 <- choose_data@meta.data[,grepl('DF.classifications',colnames(choose_data@meta.data))]
  SCdata$doubletFinder1[names(choose_data$doubletFinder1)] <- choose_data$doubletFinder1
  table(SCdata@meta.data$doubletFinder1)
}


########################    step 4.  Seurat & harmony
SCdata <- NormalizeData(SCdata, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(SCdata)
SCdata <- ScaleData(object = SCdata,features = all.genes)
SCdata <- FindVariableFeatures(object = SCdata,nfeatures = 2000)
SCdata <- RunPCA(SCdata, features = VariableFeatures(object = SCdata))
data_use<-Stdev(object = SCdata,reduction = 'pca')
sum(data_use[1:30]);sum(data_use[1:25]);sum(data_use[1:20]);sum(data_use[1:15]);sum(data_use[1:10])  ##  dims=1:30  94.85%    DIMS=30 95.97726%

table(SCdata@meta.data$orig.ident);dim(SCdata)
SCdata_harmony <- SCdata %>%  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SCdata_harmony, 'harmony')
res.used <- c(seq(0.2,6,by=0.2),0.1,0.3,0.5,0.7,0.9);res.used
SCdata_harmony <- SCdata_harmony %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = res.used) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  RunTSNE(reduction = "harmony", dims = 1:30) %>%
  identity()
SCdata_clus.tree.out <- clustree(SCdata_harmony) +theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low = "grey80", high = "red")
print(SCdata_clus.tree.out)


sel.clust = "RNA_snn_res.1"
SCdata_harmony_choose<- SetIdent(SCdata_harmony, value = sel.clust)
SCdata_harmony_choose@meta.data$seurat_clusters <- SCdata_harmony_choose@meta.data[[sel.clust]]
p1 <- DimPlot(SCdata_harmony_choose, reduction = "umap", group.by = "orig.ident")  
p2 <- DimPlot(SCdata_harmony_choose, reduction = "umap",label.size = 5,label = T);p2
table(SCdata_harmony_choose@meta.data$seurat_clusters)
# saveRDS(SCdata_harmony_choose,'SCdata_Round1_finish.rds')


########################    step 5. Figure & saved
umap_coords <- Embeddings(SCdata_harmony_choose,reduction = "umap")
cluster_info <- SCdata_harmony_choose@meta.data
allcell_info <- data.frame(cluster_info,UMAP_1=umap_coords[,1],UMAP_2=umap_coords[,2])
allcell_info <- allcell_info[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","RNA_snn_res.1","seurat_clusters","doubletFinder","UMAP_1","UMAP_2")]
allcell_info$cellbarcode = rownames(allcell_info)
allcell_info$cellbarcode = sapply(allcell_info$cellbarcode,function(x) paste0(unlist(strsplit(x,split = '[-_]'))[3:5],collapse = '_'))
write.csv(allcell_info,"1_Processing_umap_cluster_info.csv", row.names = FALSE)
allcell_info$seurat_clusters <- factor(allcell_info$seurat_clusters)
library(ggplot2)
ggplot(allcell_info,aes(x=UMAP_1,y=UMAP_2,color=seurat_clusters))+geom_point(size=1.5,alpha=0.7)+theme_minimal()+labs(x="UMAP 1",y="UMAP 2",color="Cluster")
