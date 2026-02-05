rm(list=ls())

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
library(tidyverse)
library(viridis)
library(ggsci)
library(scater)
library(devtools)
library(SingleR)   


setwd(this.path::this.dir())
getwd()


SCdata_All <- readRDS("SCdata_Round1_finish.rds")
dim(SCdata_All);table(SCdata_All$cell)
SCdata <- subset(SCdata_All,cell %in% c('CAFHSC','Endothelial'))
dim(SCdata);table(SCdata@meta.data$cell)
table(SCdata@meta.data$orig.ident);table(SCdata@meta.data$seurat_clusters)




################################################################     1. Round1   
SCdata <- NormalizeData(SCdata, normalization.method = "LogNormalize", scale.factor = 10000)
SCdata <- FindVariableFeatures(object = SCdata,nfeatures = 2000)
all.genes <- rownames(SCdata)
SCdata <- ScaleData(object = SCdata,features = all.genes)
SCdata <- RunPCA(SCdata, features = VariableFeatures(object = SCdata))
data_use<-Stdev(object = SCdata,reduction = 'pca')
sum(data_use[1:50]);sum(data_use[1:40]);sum(data_use[1:30]);sum(data_use[1:20]);sum(data_use[1:10]) 
ElbowPlot(SCdata,ndims = 50)


dim_usage=40
SCdata_harmony <- SCdata %>%  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SCdata_harmony, 'harmony')
res.used <- c(seq(0.2,6,by=0.2),0.1,0.3,0.5,0.7,0.9);res.used
SCdata_harmony <- SCdata_harmony %>%  FindNeighbors(reduction = "harmony", dims = 1:dim_usage) %>% FindClusters(resolution = res.used) %>%
  RunUMAP(reduction = "harmony", dims = 1:dim_usage) %>% RunTSNE(reduction = "harmony", dims = 1:dim_usage) %>%  identity()
SCdata_clus.tree.out <- clustree(SCdata_harmony) +theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low = "grey80", high = "red")
print(SCdata_clus.tree.out)


sel.clust = "RNA_snn_res.0.6"
SCdata_choose<- SetIdent(SCdata_harmony, value = sel.clust)
SCdata_choose@meta.data$seurat_clusters <- SCdata_choose@meta.data[[sel.clust]]
p1 <- DimPlot(SCdata_choose, reduction = "umap", group.by = "orig.ident");p1
p2 <- DimPlot(SCdata_choose, reduction = "umap",label.size = 4,label = T);p2
table(SCdata_choose@meta.data$seurat_clusters)

double <- data.frame(table(SCdata_choose@meta.data$seurat_clusters,SCdata_choose@meta.data$doubletFinder))
colnames(double) <- c('cluster','doubleres','number')
double$sum <- sapply(double$cluster,function(x) length(which(SCdata_choose@meta.data$seurat_clusters == x))) 
double$ratio <- double$number/double$sum
double[which(double$doubleres == 'Doublet' & double$ratio > 0.1),]


FeaturePlot(object =SCdata_choose, features = c("PTPRC",'MKI67'), cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0,label.size = 4)
FeaturePlot(object = SCdata_choose, features = c("GNLY","NKG7","KLRD1","KLRF1","CD3D","CD3E","CD4","CD8A","MS4A1","CD19","CD79A",'IGHG3',"MZB1","TRBC2",'TRAC'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_choose, features = c("CD14","CD68",'CD163','C1QB','LYZ','TPSAB1','CLEC9A','S100A8','S100A9'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_choose, features = c("PECAM1","CD34","CDH5","ENG",'PLVAP',"ACTA2","PDGFRB","COL1A1","PDGFRB",'COL2A1','BGN','DCN'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_choose, features = c("ALB",'AFP',"SOX2",'KRT8',"AMBP","TTR"),cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0,label.size = 4)
DotPlot(SCdata_choose,features = c("CD3D",'CD3E',"CD4","CD8A",'TRBC2','TRAC',"GNLY","NKG7",'CD79A','IGHG3','MS4A1',"CD14","CD68",'CD163','C1QB','LYZ','MKI67','STMN1','PECAM1','CD34','COL1A1','ALB',"KRT8",'TTR','AMBP','KRT18','KRT19','AFP'),assay='RNA')


hpca.se <- celldex::HumanPrimaryCellAtlasData()
clusters=SCdata_choose@meta.data$seurat_clusters
pred.hesc <- SingleR(SCdata_choose@assays$RNA@data, ref = hpca.se, labels = hpca.se$label.main,clusters = clusters)
table(pred.hesc$labels)
celltype = data.frame(ClusterID=rownames(pred.hesc), celltype=pred.hesc$labels, stringsAsFactors = F)
SCdata_choose@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
p1<-DimPlot(SCdata_choose,reduction = "umap",label = TRUE,label.size = 4)
p2<-DimPlot(SCdata_choose, reduction = "umap", group.by = "singleR",label = T,label.size = 5)
plot_grid(p1,p2)


DEGS<-FindAllMarkers(SCdata_choose,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)


SCdata_choose@meta.data$cellR1 = 'stromal'
SCdata_choose@meta.data$cellR1[which(SCdata_choose$seurat_clusters %in% c(10))] <- 'doublet_TNKandEC'
SCdata_choose@meta.data$cellR1[which(SCdata_choose$seurat_clusters %in% c(12))] <- 'doublet_myeloidandEC'
table(SCdata_choose@meta.data$cellR1)
#saveRDS(SCdata_choose,"SCdata_stromal_Round1_seurat.rds")






################################################################     2. Round2   Seurat   Process
round1Hepa_stromal <- readRDS('SCdata_Hepa_Round1_stromal.rds');dim(round1Hepa_stromal)
SCdata_stromal <- subset(SCdata_choose,cellR1 == 'stromal')
SCdata  <- merge(SCdata_stromal,round1Hepa_stromal)
dim(SCdata)  


SCdata <- NormalizeData(SCdata, normalization.method = "LogNormalize", scale.factor = 10000)
SCdata <- FindVariableFeatures(object = SCdata,nfeatures = 2000)
all.genes <- rownames(SCdata)
SCdata <- ScaleData(object = SCdata,features = all.genes)
SCdata <- RunPCA(SCdata, features = VariableFeatures(object = SCdata))
data_use<-Stdev(object = SCdata,reduction = 'pca')
sum(data_use[1:50]);sum(data_use[1:40]);sum(data_use[1:30]);sum(data_use[1:20]);sum(data_use[1:10]) 
ElbowPlot(SCdata,ndims = 50)


dim_usage=40
SCdata_harmony <- SCdata %>%  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SCdata_harmony, 'harmony')
res.used <- c(seq(0.2,6,by=0.2),0.1,0.3,0.5,0.7,0.9);res.used
SCdata_harmony <- SCdata_harmony %>%  FindNeighbors(reduction = "harmony", dims = 1:dim_usage) %>% FindClusters(resolution = res.used) %>%
  RunUMAP(reduction = "harmony", dims = 1:dim_usage) %>% RunTSNE(reduction = "harmony", dims = 1:dim_usage) %>%  identity()
SCdata_clus.tree.out <- clustree(SCdata_harmony) +theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low = "grey80", high = "red")
print(SCdata_clus.tree.out)


sel.clust = "RNA_snn_res.1.2"
SCdata_choose<- SetIdent(SCdata_harmony, value = sel.clust)
SCdata_choose@meta.data$seurat_clusters <- SCdata_choose@meta.data[[sel.clust]]
p1 <- DimPlot(SCdata_choose, reduction = "umap", group.by = "orig.ident");p1
p2 <- DimPlot(SCdata_choose, reduction = "umap",label.size = 4,label = T);p2
table(SCdata_choose@meta.data$seurat_clusters)


double <- data.frame(table(SCdata_choose@meta.data$seurat_clusters,SCdata_choose@meta.data$doubletFinder))
colnames(double) <- c('cluster','doubleres','number')
double$sum <- sapply(double$cluster,function(x) length(which(SCdata_choose@meta.data$seurat_clusters == x))) 
double$ratio <- double$number/double$sum
double[which(double$doubleres == 'Doublet' & double$ratio > 0.1),]


FeaturePlot(object = SCdata_choose, features = c("PTPRC",'MKI67'), cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0,label.size = 4)
FeaturePlot(object = SCdata_choose, features = c("GNLY","NKG7","KLRD1","KLRF1","CD3D","CD3E",'CD3G',"CD4","CD8A","MS4A1","CD19","CD79A","MZB1","TRBC2",'TRAC'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_choose, features = c("CD14","CD68",'CD163','C1QB','LYZ','TPSAB1','CLEC9A','S100A8','S100A9'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_choose, features = c("PECAM1","CD34","CDH5","ENG",'PLVAP',"ACTA2","PDGFRB","COL1A1","PDGFRB",'COL2A1','BGN','DCN',"ALB",'AFP',"SOX2",'KRT8'), cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
DotPlot(SCdata_choose,features = c('PTPRC',"CD3D",'CD3E','CD3G',"CD4","CD8A",'TRBC2','TRAC',"GNLY","NKG7",'CD79A','IGHG3','MS4A1',"CD14","CD68",'CD163','C1QB','LYZ','S100A8','MKI67','STMN1','PECAM1','CD34','COL1A1','ACTA2','ALB',"KRT8",'TTR','AMBP','KRT18','KRT19','AFP'),assay='RNA')


DEGS<-FindAllMarkers(SCdata_choose,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)


SCdata_choose@meta.data$cellR1 = 'stromal'
SCdata_choose@meta.data$cellR1[which(SCdata_choose$seurat_clusters %in% c(16,20))] <- 'doublet_delete'
SCdata_choose@meta.data$cellR1[which(SCdata_choose$seurat_clusters %in% c(21))] <- 'delete_like'
table(SCdata_choose@meta.data$cellR1)
#saveRDS(SCdata_choose,"SCdata_stromal_Round2_seurat.rds")



################################################################     3. Round3
SCdata <- subset(SCdata_choose,cellR1 == 'stromal')
dim(SCdata)

SCdata <- NormalizeData(SCdata, normalization.method = "LogNormalize", scale.factor = 10000)
SCdata <- FindVariableFeatures(object = SCdata,nfeatures = 2000)
all.genes <- rownames(SCdata)
SCdata <- ScaleData(object = SCdata,features = all.genes)
SCdata <- RunPCA(SCdata, features = VariableFeatures(object = SCdata))
data_use<-Stdev(object = SCdata,reduction = 'pca')
sum(data_use[1:50]);sum(data_use[1:40]);sum(data_use[1:30]);sum(data_use[1:20]);sum(data_use[1:10]) 
ElbowPlot(SCdata,ndims = 50)


dim_usage=40
SCdata_harmony <- SCdata %>%  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SCdata_harmony, 'harmony')
res.used <- c(seq(0.2,8,by=0.2),0.1,0.3,0.5,0.7,0.9);res.used
SCdata_harmony <- SCdata_harmony %>%  FindNeighbors(reduction = "harmony", dims = 1:dim_usage) %>% FindClusters(resolution = res.used) %>%
  RunUMAP(reduction = "harmony", dims = 1:dim_usage) %>% RunTSNE(reduction = "harmony", dims = 1:dim_usage) %>%  identity()
SCdata_clus.tree.out <- clustree(SCdata_harmony) +theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low = "grey80", high = "red")
print(SCdata_clus.tree.out)


sel.clust = "RNA_snn_res.0.9"
SCdata_choose <- SetIdent(SCdata_harmony, value = sel.clust)
SCdata_choose@meta.data$seurat_clusters <- SCdata_choose@meta.data[[sel.clust]]
p1 <- DimPlot(SCdata_choose, reduction = "umap", group.by = "orig.ident");p1
p2 <- DimPlot(SCdata_choose, reduction = "umap",label.size = 4,label = T);p2
table(SCdata_choose@meta.data$seurat_clusters)

double <- data.frame(table(SCdata_choose@meta.data$seurat_clusters,SCdata_choose@meta.data$doubletFinder))
colnames(double) <- c('cluster','doubleres','number')
double$sum <- sapply(double$cluster,function(x) length(which(SCdata_choose@meta.data$seurat_clusters == x))) 
double$ratio <- double$number/double$sum
double[which(double$doubleres == 'Doublet' & double$ratio > 0.1),]


FeaturePlot(object = SCdata_choose, features = c("PTPRC",'MKI67'), cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0,label.size = 4)
FeaturePlot(object = SCdata_choose, features = c("GNLY","NKG7","KLRD1","KLRF1","CD3D","CD3E",'CD3G',"CD4","CD8A","MS4A1","CD19","CD79A",'IGHG3',"MZB1","TRBC2",'TRAC'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_choose, features = c("CD14","CD68",'CD163','C1QB','LYZ','TPSAB1','CLEC9A','S100A8','S100A9'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_choose, features = c("PECAM1","CD34","CDH5","ENG",'PLVAP',"ACTA2","PDGFRB","COL1A1","PDGFRB",'COL2A1','BGN','DCN',"ALB",'AFP',"SOX2",'KRT8'), cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
DotPlot(SCdata_choose,features = c('PTPRC',"CD3D",'CD3E',"CD4","CD8A",'TRBC2','TRAC',"GNLY","NKG7",'CD79A','IGHG3','MS4A1',"CD14","CD68",'CD163','C1QB','LYZ','S100A8','MKI67','STMN1','PECAM1','CD34','COL1A1','ACTA2','ALB',"KRT8",'TTR','AMBP','KRT18','KRT19','AFP'),assay='RNA')


DEGS<-FindAllMarkers(SCdata_choose,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)


SCdata_choose@meta.data$cellR1 = 'stromal'
table(SCdata_choose@meta.data$cellR1)
SCdata_choose$label <- 0
SCdata_choose$label[which(SCdata_choose$seurat_clusters == 0)] <- 'EC_KDR_KLF4'
SCdata_choose$label[which(SCdata_choose$seurat_clusters == 1)] <- 'EC_KDR_ESM1'
SCdata_choose$label[which(SCdata_choose$seurat_clusters == 3)] <- 'EC_IGFBP3'
SCdata_choose$label[which(SCdata_choose$seurat_clusters == 4)] <- 'EC_IGFBP3_KLF4'
SCdata_choose$label[which(SCdata_choose$seurat_clusters == 5)] <- 'EC_KDR_CXCR4'
SCdata_choose$label[which(SCdata_choose$seurat_clusters == 6)] <- 'EC_MALL'
SCdata_choose$label[which(SCdata_choose$seurat_clusters == 7)] <- 'EC_IGFBP3_CXCL2'
SCdata_choose$label[which(SCdata_choose$seurat_clusters == 9)] <- 'EC_ANXA1'
SCdata_choose$label[which(SCdata_choose$seurat_clusters == 11)] <- 'EC_APOA2'
SCdata_choose$label[which(SCdata_choose$seurat_clusters == 12)] <- 'EC_STMN1'
SCdata_choose$label[which(SCdata_choose$seurat_clusters == 13)] <- 'EC_CLEC4G'
SCdata_choose$label[which(SCdata_choose$seurat_clusters == 14)] <- 'EC_CXCL10'
SCdata_choose$label[which(SCdata_choose$seurat_clusters == 2)] <- 'Me_CD36_STEAP4'
SCdata_choose$label[which(SCdata_choose$seurat_clusters %in% c(8,16))] <- 'Me_COL1A1_C7'
SCdata_choose$label[which(SCdata_choose$seurat_clusters == 10)] <- 'Me_MYH11'
SCdata_choose$label[which(SCdata_choose$seurat_clusters == 15)] <- 'Me_FABP3_DAPL1'
table(SCdata_choose$label)
#saveRDS(SCdata_choose,"SCdata_stromal_seurat.rds")





# label <- unique(SCdata_choose$label);label
# DiffDegs <- NULL
# for(i in 1:length(label)){
#   SCdata_choose$labeluse <- 'con'
#   SCdata_choose$labeluse[which(SCdata_choose$label == label[i])] <- 'Treat'
#   table(SCdata_choose$labeluse)
#   DEGS <-FindMarkers(SCdata_choose,ident.1 = 'Treat',ident.2 ='con',group.by = 'labeluse',min.pct = 0.25,logfc.threshold = 0.25,only.pos = TRUE)
#   DEGS$genes <- rownames(DEGS)
#   DEGS$clusterlabel = label[i]
#   DiffDegs <- rbind(DiffDegs,DEGS)
# }
# #saveRDS(DiffDegs,"Stromal_DEGs_OnlyUP.rds")


