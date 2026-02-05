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
library(export) 
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(GSVA)
library(ggprism)
library(limma)
library(ggthemes)
library(scibet)   
library(tidyverse)
library(viridis)
library(ggsci)
library(scater)
library(devtools)
library(SingleR)    



################################################################        Cell Type Annoation
########################    step 1. SingleR
hpca.se <- celldex::HumanPrimaryCellAtlasData()
clusters=SCdata_harmony_choose@meta.data$seurat_clusters
pred.hesc <- SingleR(SCdata_harmony_choose@assays$RNA@data, ref = hpca.se, labels = hpca.se$label.main,clusters = clusters)
table(pred.hesc$labels)
celltype = data.frame(ClusterID=rownames(pred.hesc), celltype=pred.hesc$labels, stringsAsFactors = F)
SCdata_harmony_choose@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
p1<-DimPlot(SCdata_harmony_choose,reduction = "umap",label = TRUE,label.size = 4)
p2<-DimPlot(SCdata_harmony_choose, reduction = "umap", group.by = "singleR",label = T,label.size = 5)
plot_grid(p1,p2)


########################    step 2. SciBet
exp_tpm <- calculateTPM(SCdata_harmony_choose@assays$RNA@counts)
exp_tpm <- t(exp_tpm)
model <- readr::read_csv("major_human_cell_types.csv")
model <- pro.core(model)
prd <- LoadModel(model)
scibettype <- prd(exp_tpm)
SCdata_harmony_choose$scibet <- scibettype
p3<-DimPlot(SCdata_harmony_choose, reduction = "umap", group.by = "scibet",label = T,label.size = 5)
plot_grid(p1,p3)


########################    step 3. Diff-Genes
DEGS<-FindAllMarkers(SCdata_harmony_choose,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


########################    step 4. Marker-gene Expression
FeaturePlot(object = SCdata_harmony_choose, features = c("PTPRC",'MKI67'), cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0,label.size = 4)
FeaturePlot(object = SCdata_harmony_choose, features = c("CD3D","CD3E","CD4","CD8A","GNLY","NKG7","KLRD1","KLRF1"),cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0,label.size = 4)
FeaturePlot(object = SCdata_harmony_choose, features = c("MS4A1","CD19","CD79A","MZB1"),cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0,label.size = 4)
FeaturePlot(object = SCdata_harmony_choose, features = c("LYZ","CD14","CD68","CD33",'CD163','C1QB','TPSAB1','TPSB2','CLEC9A','CD1C','LAMP3'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0,label.size = 4)
FeaturePlot(object = SCdata_harmony_choose, features = c("PECAM1","CD34","CDH5","ENG","VWF","PLVAP"),cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0,label.size = 4)
FeaturePlot(object = SCdata_harmony_choose, features = c("ACTA2","PDGFRB","COL1A1","PDGFRB"),cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0,label.size = 4)
FeaturePlot(object = SCdata_harmony_choose, features = c("ACTA2",'COL1A2','BGN','DCN'),cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0,label.size = 4)
FeaturePlot(object =SCdata_harmony_choose, features = c("ALB",'AFP',"SOX2",'KRT8'), cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0,label.size = 4)


########################    step 5. Finish cell type anno
SCdata_harmony_choose@meta.data$cell<-'others'
SCdata_harmony_choose@meta.data$cell[which(SCdata_harmony_choose@meta.data[["seurat_clusters"]] %in% c(2,3,4,7,14,17,22,29))]<-"Tcell"
SCdata_harmony_choose@meta.data$cell[which(SCdata_harmony_choose@meta.data[["seurat_clusters"]] %in% c(10))]<-"NK"
SCdata_harmony_choose@meta.data$cell[which(SCdata_harmony_choose@meta.data[["seurat_clusters"]] %in% c(12,19,28))]<-"Bcell"
SCdata_harmony_choose@meta.data$cell[which(SCdata_harmony_choose@meta.data[["seurat_clusters"]] %in% c(1,6,8,18,20,21,23))]<-"myeloid"
SCdata_harmony_choose@meta.data$cell[which(SCdata_harmony_choose@meta.data[["seurat_clusters"]] %in% c(0,9,11,13,25,26))]<-"Hepatocytes"
SCdata_harmony_choose@meta.data$cell[which(SCdata_harmony_choose@meta.data[["seurat_clusters"]] %in% c(5,24,27,15))]<-"Endothelial"
SCdata_harmony_choose@meta.data$cell[which(SCdata_harmony_choose@meta.data[["seurat_clusters"]] %in% c(16))]<-"CAFHSC"
table(SCdata_harmony_choose@meta.data$cell)
DimPlot(SCdata_harmony_choose,reduction = "umap",group.by = 'cell',label = TRUE,label.size = 5)
# saveRDS(SCdata_harmony_choose,'SCdata_Round1_finish.rds')


########################    step 6. Figure & saved
umap_coords <- Embeddings(SCdata_harmony_choose,reduction = "umap")
cluster_info <- SCdata_harmony_choose@meta.data
allcell_info <- data.frame(cluster_info,UMAP_1=umap_coords[,1],UMAP_2=umap_coords[,2])
allcell_info <- allcell_info[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","RNA_snn_res.1","seurat_clusters","doubletFinder", "doubletFinder1","singleR","scibet","cell","UMAP_1","UMAP_2")]
allcell_info$cellbarcode = rownames(allcell_info)
allcell_info$cellbarcode = sapply(allcell_info$cellbarcode,function(x) paste0(unlist(strsplit(x,split = '[-_]'))[3:5],collapse = '_'))
write.csv(allcell_info,"2_CellAnno_Round1_info.csv", row.names = FALSE)
library(ggplot2)
allcell_info$cell <- factor(allcell_info$cell)
ggplot(allcell_info,aes(x=UMAP_1,y=UMAP_2,color=cell))+geom_point(size=1.5,alpha=0.7)+theme_minimal()+labs(x="UMAP 1",y="UMAP 2",color="Cluster")
ggplot(allcell_info,aes(x=UMAP_1,y=UMAP_2,color=singleR))+geom_point(size=1.5,alpha=0.7)+theme_minimal()+labs(x="UMAP 1",y="UMAP 2",color="Cluster")








################################################################     InferCNV_Round1
###  each sample: all immune as reference and stromal and Hepa predicted
library(infercnv)
sampleID <- unique(SCdata_harmony_choose$orig.ident)
SCdata_harmony_choose$cellR1 = 'immune'
SCdata_harmony_choose$cellR1[which(SCdata_harmony_choose$cell %in% c('CAFHSC','Endothelial'))] <-  'stromal'
SCdata_harmony_choose$cellR1[which(SCdata_harmony_choose$cell == 'Hepatocytes')] <-  'Hepatocytes'
table(SCdata_harmony_choose$cellR1)
for(i in 1:length(sampleID)){
  count_sample <- subset(SCdata_harmony_choose,orig.ident == sampleID[i])
  count <- GetAssayData(count_sample,slot = 'counts')
  count <- as.data.frame(count)
  
  geneAnno<-read.table("geneAnno.txt",header = F,sep="\t");geneAnno<-data.frame(geneAnno)
  colnames(geneAnno)<-c("SYMBOL","chr","start","end")
  geneInfor<-geneAnno[which(geneAnno$SYMBOL %in% rownames(count)),]
  geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
  write.table(geneInfor,file = 'geneFileNormal1.txt',sep = '\t',quote = F,col.names = F,row.names = F)
  
  groupinfo<-data.frame(v1=colnames(count),v2<-count_sample$cellR1)
  write.table(groupinfo,file = 'groupinfoNormal1.txt',sep = '\t',quote = F,col.names = F,row.names = F)
  write.table(count,file = 'expFileNormal1.txt',sep = '\t',quote = F)
  
  gene <- c("chrMT","chrGL000009.2","chrGL000194.1","chrGL000195.1","chrGL000219.1","chrKI270734.1","chrGL000218.1","chrKI270721.1","chrKI270726.1","chrKI270711.1","chrY")
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix= 'expFileNormal1.txt',annotations_file='groupinfoNormal1.txt',
                                      delim="\t",gene_order_file= 'geneFileNormal1.txt',ref_group_names=c("immune"),chr_exclude = gene)
  infercnv_obj2 = infercnv::run(infercnv_obj,cutoff=0.1,out_dir=paste0("InferCNV_",sampleID[i],"_ImmuneRef"),cluster_by_groups=T,denoise=FALSE,HMM=FALSE,num_threads = 20)
  rm(infercnv_obj,infercnv_obj2);gc()
  gc()
}