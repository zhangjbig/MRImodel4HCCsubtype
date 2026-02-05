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
SCdata <- subset(SCdata_All,cell %in% c('Hepatocytes'))
dim(SCdata);table(SCdata@meta.data$cell);table(SCdata@meta.data$orig.ident)



################################################################     1. Round1   Seurat   Process
SCdata <- NormalizeData(SCdata, normalization.method = "LogNormalize", scale.factor = 10000)
SCdata <- FindVariableFeatures(object = SCdata,nfeatures = 2000)
all.genes <- rownames(SCdata)
SCdata <- ScaleData(object = SCdata,features = all.genes)
SCdata <- RunPCA(SCdata, features = VariableFeatures(object = SCdata))
data_use<-Stdev(object = SCdata,reduction = 'pca')
sum(data_use[1:50]);sum(data_use[1:40]);sum(data_use[1:30]);sum(data_use[1:20]);sum(data_use[1:10]) 
ElbowPlot(SCdata,ndims = 50)


dim_usage=40;sum(data_use[1:40])
SCdata_harmony <- SCdata %>%  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SCdata_harmony, 'harmony')
res.used <- c(seq(0.2,6,by=0.2),0.1,0.3,0.5,0.7,0.9);res.used
SCdata_harmony <- SCdata_harmony %>%  FindNeighbors(reduction = "harmony", dims = 1:dim_usage) %>% FindClusters(resolution = res.used) %>%
  RunUMAP(reduction = "harmony", dims = 1:dim_usage) %>% RunTSNE(reduction = "harmony", dims = 1:dim_usage) %>%  identity()
SCdata_clus.tree.out <- clustree(SCdata_harmony) +theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low = "grey80", high = "red")
print(SCdata_clus.tree.out)


sel.clust = "RNA_snn_res.2"
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
FeaturePlot(object = SCdata_choose, features = c("GNLY","NKG7","KLRD1","KLRF1","CD3D","CD3E","CD4","CD8A","MS4A1","CD19","CD79A","MZB1","TRBC2",'TRAC'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_choose, features = c("CD14","CD68",'CD163','C1QB','LYZ','TPSAB1','CLEC9A','S100A8','S100A9'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_choose, features = c("PECAM1","CD34","CDH5","ENG",'PLVAP',"ACTA2","PDGFRB","COL1A1","PDGFRB",'COL2A1','BGN','DCN'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_choose, features = c("ALB",'AFP',"SOX2",'KRT8',"AMBP","TTR"),cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0,label.size = 4)
DotPlot(SCdata_choose,features = c('PTPRC',"CD3D",'CD3E','CD3G',"CD4","CD8A",'TRBC2','TRAC',"GNLY","NKG7",'CD79A','IGHG3','MS4A1',"CD14","CD68",'CD163','C1QB','LYZ','S100A8','MKI67','STMN1','PECAM1','CD34','COL1A1','ACTA2','ALB',"KRT8",'TTR','AMBP','KRT18','KRT19','AFP'),assay='RNA')
p=DotPlot(SCdata_choose,features = c('PTPRC',"CD3D",'CD3E','CD3G',"CD4","CD8A",'TRBC2','GNLY','NKG7','MS4A1','CD79A','CD68','CD163','LYZ','C1QB','CD14','S100A8','PECAM1','CD34','COL1A1','ACTA2','ALB','KRT8','AMBP','TTR','AFP'),assay='RNA');p
p=p$data;p=p[which(p$avg.exp > 1 & p$pct.exp > 25),]


hpca.se <- celldex::HumanPrimaryCellAtlasData()
clusters=SCdata_choose@meta.data$seurat_clusters
pred.hesc <- SingleR(SCdata_choose@assays$RNA@data, ref = hpca.se, labels = hpca.se$label.main,clusters = clusters)
table(pred.hesc$labels)
celltype = data.frame(ClusterID=rownames(pred.hesc), celltype=pred.hesc$labels, stringsAsFactors = F)
SCdata_choose@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
p1<-DimPlot(SCdata_choose,reduction = "umap",label = TRUE,label.size = 4)
p2<-DimPlot(SCdata_choose, reduction = "umap", group.by = "singleR",label = T,label.size = 5)
p1+p2


DEGS<-FindAllMarkers(SCdata_choose,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)


SCdata_choose$cellR1 <- 'hepatocytes'
SCdata_choose$cellR1[which(SCdata_choose$seurat_clusters == 15)] <- 'doublet-delete_TNKandKRT8'
SCdata_choose$cellR1[which(SCdata_choose$seurat_clusters %in% c(20,25))] <- 'doublet-delete_MyeloidandKRT8'
SCdata_choose$cellR1[which(SCdata_choose$seurat_clusters == 28)] <- 'stromal_like'
table(SCdata_choose$cellR1)
# saveRDS(SCdata_choose,"SCdata_hepatocytes_Round1_seurat.rds")
# SCdata_stromal <- subset(SCdata_choose,cellR1 == 'stromal_like')
# saveRDS(SCdata_stromal,"SCdata_Hepa_Round1_stromal.rds")




################################################################     2. Round2   Seurat   Process
SCdata_Tcell_Round1_hepatocytes <- readRDS("SCdata_Tcell_Round1_hepatocytes.rds")
SCdata_Bcell_Round1_Hepa <- readRDS("SCdata_Bcell_Round1_Hepa.rds")
SCdata_hepa <- subset(SCdata_choose,cellR1 %in% c('hepatocytes'))
SCdata <- merge(SCdata_hepa,y=c(SCdata_Tcell_Round1_hepatocytes,SCdata_Bcell_Round1_Hepa))
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


sel.clust = "RNA_snn_res.1.6"
SCdata_choose<- SetIdent(SCdata_harmony, value = sel.clust)
SCdata_choose@meta.data$seurat_clusters <- SCdata_choose@meta.data[[sel.clust]]
p1 <- DimPlot(SCdata_choose, reduction = "umap", group.by = "orig.ident") ;p1
p2 <- DimPlot(SCdata_choose, reduction = "umap",label.size = 4,label = T) ;p2
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
DotPlot(SCdata_choose,features = c('PTPRC',"CD3D",'CD3E',"CD4","CD8A",'TRBC2','TRAC',"GNLY","NKG7",'CD79A','IGHG3','MS4A1',"CD14","CD68",'CD163','C1QB','LYZ','S100A8','MKI67','STMN1','PECAM1','CD34','COL1A1','ACTA2','ALB',"KRT8",'TTR','AMBP','KRT18','KRT19','AFP'),assay='RNA')
p=DotPlot(SCdata_choose,features = c('PTPRC',"CD3D",'CD3E','CD3G',"CD4","CD8A",'TRBC2','TRAC','GNLY','NKG7','CD160','MZB1','MS4A1','CD79A','IGHG3','CD68','CD163','LYZ','C1QB','PECAM1','COL1A1','ACTA2','ALB','KRT8','AMBP','TTR','AFP'),assay='RNA');p
p=p$data;p=p[which(p$avg.exp > 1 & p$pct.exp > 25),]


DEGS<-FindAllMarkers(SCdata_choose,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)


SCdata_choose$cellR1 <- 'hepatocytes'
SCdata_choose$cellR1[which(SCdata_choose$seurat_clusters == 18)] <- 'Tcell_like'
SCdata_choose$cellR1[which(SCdata_choose$seurat_clusters == 28)] <- 'Bcell_like'
table(SCdata_choose$cellR1)
SCdata_T <- subset(SCdata_choose,cellR1 == 'Tcell_like')
SCdata_Bcell <- subset(SCdata_choose,cellR1 == 'Bcell_like')
# saveRDS(SCdata_T,"SCdata_Hepa_Round2_Tcell.rds")
# saveRDS(SCdata_Bcell,"SCdata_Hepa_Round2_Bcell.rds")
# saveRDS(SCdata_choose,"SCdata_hepatocytes_Round2_seurat.rds")




################################################################    3.    Round3   Seurat   Process
SCdata1 <- subset(SCdata_choose,cellR1 == 'hepatocytes')
SCdata_Tcell_Round2_Hepa <- readRDS("SCdata_Tcell_Round2_Hepa.rds");dim(SCdata_Tcell_Round2_Hepa)
SCdata_CD4T_Round1_Hepa <- readRDS("SCdata_CD4T_Round1_Hepa.rds");dim(SCdata_CD4T_Round1_Hepa)
SCdata <- merge(SCdata1,y=c(SCdata_Tcell_Round2_Hepa,SCdata_CD4T_Round1_Hepa))
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


sel.clust = "RNA_snn_res.0.6"
SCdata_choose<- SetIdent(SCdata_harmony, value = sel.clust)
SCdata_choose@meta.data$seurat_clusters <- SCdata_choose@meta.data[[sel.clust]]
p1 <- DimPlot(SCdata_choose, reduction = "umap", group.by = "orig.ident") ;p1
p2 <- DimPlot(SCdata_choose, reduction = "umap",label.size = 4,label = T) ;p2
table(SCdata_choose@meta.data$seurat_clusters)


double <- data.frame(table(SCdata_choose@meta.data$seurat_clusters,SCdata_choose@meta.data$doubletFinder))
colnames(double) <- c('cluster','doubleres','number')
double$sum <- sapply(double$cluster,function(x) length(which(SCdata_choose@meta.data$seurat_clusters == x))) 
double$ratio <- double$number/double$sum
double[which(double$doubleres == 'Doublet' & double$ratio > 0.1),]


FeaturePlot(object = SCdata_choose, features = c("PTPRC",'MKI67',"GNLY","NKG7","KLRD1","KLRF1","CD3D","CD3E","CD4","CD8A","MS4A1","CD19","CD79A","IGHG3","MZB1","TRBC2",'TRAC'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_choose, features = c("CD14","CD68",'CD163','C1QB','LYZ','TPSAB1','CLEC9A','S100A8','CLEC9A','CLEC10A','TPSAB1','VCAN'),cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_choose, features = c("PECAM1","CD34","CDH5","ENG",'PLVAP',"ACTA2","PDGFRB","COL1A1",'COL2A1','BGN','DCN',"ALB",'AFP',"SOX2",'KRT8'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
DotPlot(SCdata_choose,features = c('PTPRC',"CD3D",'CD3E','CD3G',"CD4","CD8A",'TRBC2','TRAC',"GNLY","NKG7",'CD79A','IGHG3','MS4A1',"CD14","CD68",'CD163','C1QB','LYZ','CD1C','CLEC9A','S100A8','MKI67','STMN1','PECAM1','CD34','COL1A1','ACTA2','ALB',"KRT8",'TTR','AMBP','KRT18','KRT19','AFP'),assay='RNA')
p=DotPlot(SCdata_choose,features = c('PTPRC',"CD3D",'CD3E','CD3G',"CD4","CD8A",'TRBC2','TRAC','GNLY','NKG7','CD160','MZB1','MS4A1','CD79A','IGHG3','CD68','CD163','LYZ','C1QB','S100A8','PECAM1','COL1A1','ACTA2','ALB','KRT8','AMBP','TTR','AFP'),assay='RNA');p
p=p$data;p=p[which(p$avg.exp > 1 & p$pct.exp > 25),]


DEGS<-FindAllMarkers(SCdata_choose,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)


SCdata_choose$cellR1 <- 'hepatocytes'
SCdata_choose$cellR1[which(SCdata_choose$seurat_clusters == c(13))] <- 'delete_PTPRCandHepa'
SCdata_choose$cellR1[which(SCdata_choose$seurat_clusters == c(9))] <-  'delete_CD4CD68andHepa'
SCdata_choose$cellR1[which(SCdata_choose$seurat_clusters == c(14))] <- 'doublet_delete'
table(SCdata_choose$cellR1)
#saveRDS(SCdata_choose,"SCdata_hepatocytes_Round3_seurat.rds")
# SCdata_T <- subset(SCdata_choose,cellR1 == 'Tcell_delete');dim(SCdata_T)
# SCdata_CD68 <- subset(SCdata_choose,cellR1 == 'CD68_delete');dim(SCdata_CD68)
# saveRDS(SCdata_T,"SCdata_Hepa_Round3_CD4Tcell_delete.rds")
# saveRDS(SCdata_CD68,"SCdata_Hepa_Round3_CD68_delete.rds")



################################################################    4.   Round4   Seurat   Process
SCdata <- subset(SCdata_choose,cellR1 == 'hepatocytes')
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
res.used <- c(seq(0.2,8,by=0.2),0.1,0.3,0.5,0.7);res.used
SCdata_harmony <- SCdata_harmony %>%  FindNeighbors(reduction = "harmony", dims = 1:dim_usage) %>% FindClusters(resolution = res.used) %>%
  RunUMAP(reduction = "harmony", dims = 1:dim_usage) %>% RunTSNE(reduction = "harmony", dims = 1:dim_usage) %>%  identity()
SCdata_clus.tree.out <- clustree(SCdata_harmony) +theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low = "grey80", high = "red")
print(SCdata_clus.tree.out)


sel.clust = "RNA_snn_res.0.3"
SCdata_choose<- SetIdent(SCdata_harmony, value = sel.clust)
SCdata_choose@meta.data$seurat_clusters <- SCdata_choose@meta.data[[sel.clust]]
p1 <- DimPlot(SCdata_choose, reduction = "umap", group.by = "orig.ident")  
p2 <- DimPlot(SCdata_choose, reduction = "umap",label.size = 4,label = T);p2
table(SCdata_choose@meta.data$seurat_clusters)

double <- data.frame(table(SCdata_choose@meta.data$seurat_clusters,SCdata_choose@meta.data$doubletFinder))
colnames(double) <- c('cluster','doubleres','number')
double$sum <- sapply(double$cluster,function(x) length(which(SCdata_choose@meta.data$seurat_clusters == x))) 
double$ratio <- double$number/double$sum
double[which(double$doubleres == 'Doublet' & double$ratio > 0.1),]


FeaturePlot(object = SCdata_choose, features = c("PTPRC",'MKI67',"GNLY","NKG7","KLRD1","KLRF1","CD3D","CD3E","CD4","CD8A","MS4A1","CD19","CD79A","IGHG3","MZB1","TRBC2",'TRAC'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_choose, features = c("CD14","CD68",'CD163','C1QB','LYZ','TPSAB1','CLEC9A','S100A8','CLEC9A','CLEC10A','TPSAB1','VCAN'),cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_choose, features = c("PECAM1","CD34","CDH5","ENG",'PLVAP',"ACTA2","PDGFRB","COL1A1",'COL2A1','BGN','DCN',"ALB",'AFP',"SOX2",'KRT8'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
DotPlot(SCdata_choose,features = c('PTPRC',"CD3D",'CD3E','CD3G',"CD4","CD8A",'TRBC2','TRAC',"GNLY","NKG7",'CD79A','IGHG3','MS4A1',"CD14","CD68",'CD163','C1QB','LYZ','CD1C','CLEC9A','S100A8','MKI67','STMN1','PECAM1','CD34','COL1A1','ACTA2','ALB',"KRT8",'TTR','AMBP','KRT18','KRT19','AFP'),assay='RNA')
VlnPlot(SCdata_choose,features = c('PTPRC',"CD3D","CD4","CD8A",'TRBC2','GNLY','NKG7','CD79A','CD68','CD163','LYZ','C1QB','S100A8','ALB','KRT8','AFP'),assay='RNA')
p = DotPlot(SCdata_choose,features = c('PTPRC',"CD3D",'CD3E','CD3G',"CD4","CD8A",'TRBC2','TRAC',"GNLY","NKG7",'CD79A','IGHG3','MS4A1',"CD14","CD68",'CD163','C1QB','LYZ','CD1C','CLEC9A','S100A8','MKI67','STMN1','PECAM1','CD34','COL1A1','ACTA2','ALB',"KRT8",'TTR','AMBP','KRT18','KRT19','AFP'),assay='RNA')
p1 = p$data;p1=p1[which(p1$avg.exp > 1 & p1$pct.exp > 25),]


DEGS<-FindAllMarkers(SCdata_choose,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)


SCdata_choose$cellR1 <- 'hepatocytes'
SCdata_choose$cellR1[which(SCdata_choose$seurat_clusters == 9)] <- 'delete_CD68andHepa'
table(SCdata_choose$cellR1)
# saveRDS(SCdata_choose,"SCdata_hepatocytes_Round4_seurat.rds")




################################################################    5.   Round5   Seurat   Process
SCdata <- subset(SCdata_choose,cellR1 == 'hepatocytes')
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
res.used <- c(seq(0.2,8,by=0.2),0.1,0.3,0.5,0.7);res.used
SCdata_harmony <- SCdata_harmony %>%  FindNeighbors(reduction = "harmony", dims = 1:dim_usage) %>% FindClusters(resolution = res.used) %>%
  RunUMAP(reduction = "harmony", dims = 1:dim_usage) %>% RunTSNE(reduction = "harmony", dims = 1:dim_usage) %>%  identity()
SCdata_clus.tree.out <- clustree(SCdata_harmony) +theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low = "grey80", high = "red")
print(SCdata_clus.tree.out)


sel.clust = "RNA_snn_res.0.4"
SCdata_choose<- SetIdent(SCdata_harmony, value = sel.clust)
SCdata_choose@meta.data$seurat_clusters <- SCdata_choose@meta.data[[sel.clust]]
p1 <- DimPlot(SCdata_choose, reduction = "umap", group.by = "orig.ident")  
p2 <- DimPlot(SCdata_choose, reduction = "umap",label.size = 4,label = T);p2
table(SCdata_choose@meta.data$seurat_clusters)


double <- data.frame(table(SCdata_choose@meta.data$seurat_clusters,SCdata_choose@meta.data$doubletFinder))
colnames(double) <- c('cluster','doubleres','number')
double$sum <- sapply(double$cluster,function(x) length(which(SCdata_choose@meta.data$seurat_clusters == x))) 
double$ratio <- double$number/double$sum
double[which(double$doubleres == 'Doublet' & double$ratio > 0.1),]


FeaturePlot(object = SCdata_choose, features = c("PTPRC",'MKI67',"GNLY","NKG7","KLRD1","KLRF1","CD3D","CD3E","CD4","CD8A","MS4A1","CD19","CD79A","IGHG3","MZB1","TRBC2",'TRAC'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_choose, features = c("CD14","CD68",'CD163','C1QB','LYZ','TPSAB1','CLEC9A','S100A8','CLEC9A','CLEC10A','TPSAB1','VCAN'),cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_choose, features = c("PECAM1","CD34","CDH5","ENG",'PLVAP',"ACTA2","PDGFRB","COL1A1",'COL2A1','BGN','DCN',"ALB",'AFP',"SOX2",'KRT8'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
DotPlot(SCdata_choose,features = c('PTPRC',"CD3D",'CD3E','CD3G',"CD4","CD8A",'TRBC2','TRAC',"GNLY","NKG7",'CD79A','IGHG3','MS4A1',"CD14","CD68",'CD163','C1QB','LYZ','CD1C','CLEC9A','S100A8','MKI67','STMN1','PECAM1','CD34','COL1A1','ACTA2','ALB',"KRT8",'TTR','AMBP','KRT18','KRT19','AFP'),assay='RNA')
VlnPlot(SCdata_choose,features = c('PTPRC',"CD3D","CD4","CD8A",'TRBC2','GNLY','NKG7','CD79A','CD68','CD163','LYZ','C1QB','S100A8','ALB','KRT8','AFP'),assay='RNA')


DEGS<-FindAllMarkers(SCdata_choose,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)


SCdata_choose$cellR1 <- 'hepatocytes'
table(SCdata_choose$cellR1)
# saveRDS(SCdata_choose,"SCdata_hepatocytes_Round5_seurat.rds")


################################################################    6. Hepa - cell clustering
SCdata_choose$label <- 0
SCdata_choose$label[which(SCdata_choose$seurat_clusters %in% c(0))] <- 'Hepa_SERPINC1_HP'
SCdata_choose$label[which(SCdata_choose$seurat_clusters %in% c(1))] <- 'Hepa_IGFBP1_CP'
SCdata_choose$label[which(SCdata_choose$seurat_clusters %in% c(2))] <- 'Hepa_cl2'
SCdata_choose$label[which(SCdata_choose$seurat_clusters %in% c(4,6))] <- 'Hepa_STMN1'
SCdata_choose$label[which(SCdata_choose$seurat_clusters %in% c(3))] <- 'Hepa_EGR1_ATF3'
SCdata_choose$label[which(SCdata_choose$seurat_clusters %in% c(5))] <- 'Hepa_PKM_TM4SF'
SCdata_choose$label[which(SCdata_choose$seurat_clusters %in% c(7))] <- 'Hepa_NEAT1_CAPN12'
SCdata_choose$label[which(SCdata_choose$seurat_clusters %in% c(8))] <- 'Hepa_IGFBP7'
SCdata_choose$label[which(SCdata_choose$seurat_clusters %in% c(9))] <- 'Hepa_DKK3'
#saveRDS(SCdata_choose,"SCdata_hepatocytes_seurat.rds")

