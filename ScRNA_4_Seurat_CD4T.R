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
library(tidyverse)
library(viridis)
library(ggsci)
library(scater)
library(devtools)
library(SingleR)  



SCdata <- readRDS("SCdata_Tcell_Round2_CD4.rds")
dim(SCdata)



################################################################     1. Round1   Seurat   Process
SCdata <- NormalizeData(SCdata, normalization.method = "LogNormalize", scale.factor = 10000)
SCdata <- FindVariableFeatures(object = SCdata,nfeatures = 2000)
all.genes <- rownames(SCdata)
SCdata <- ScaleData(object = SCdata,features = all.genes)
SCdata <- RunPCA(SCdata, features = VariableFeatures(object = SCdata))
data_use<-Stdev(object = SCdata,reduction = 'pca')
sum(data_use[1:50]);sum(data_use[1:40]);sum(data_use[1:30]);sum(data_use[1:20]);sum(data_use[1:10]) 
ElbowPlot(SCdata,ndims = 50)


dim_usage=50
SCdata_harmony <- SCdata %>%  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SCdata_harmony, 'harmony')
res.used <- c(seq(0.2,8,by=0.2),0.1,0.3,0.5,0.7,0.9);res.used
SCdata_harmony <- SCdata_harmony %>% FindNeighbors(reduction = "harmony", dims = 1:dim_usage) %>% FindClusters(resolution = res.used) %>%
  RunUMAP(reduction = "harmony", dims = 1:dim_usage) %>% RunTSNE(reduction = "harmony", dims = 1:dim_usage) %>% identity()
SCdata_clus.tree.out <- clustree(SCdata_harmony) +theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low = "grey80", high = "red")
print(SCdata_clus.tree.out)


sel.clust = 'RNA_snn_res.1.4'
SCdata_harmony_choose<- SetIdent(SCdata_harmony, value = sel.clust)
SCdata_harmony_choose@meta.data$seurat_clusters <- SCdata_harmony_choose@meta.data[[sel.clust]]
p1 <- DimPlot(SCdata_harmony_choose, reduction = "umap", group.by = "orig.ident")  
p2 <- DimPlot(SCdata_harmony_choose, reduction = "umap",label.size = 4,label = T);p2 
table(SCdata_harmony_choose@meta.data$seurat_clusters)

double <- data.frame(table(SCdata_harmony_choose@meta.data$seurat_clusters,SCdata_harmony_choose@meta.data$doubletFinder))
colnames(double) <- c('cluster','doubleres','number')
double$sum <- sapply(double$cluster,function(x) length(which(SCdata_harmony_choose@meta.data$seurat_clusters == x))) 
double$ratio <- double$number/double$sum
double[which(double$doubleres == 'Doublet' & double$ratio > 0.1),]


FeaturePlot(object = SCdata_harmony_choose, features = c("PTPRC",'MKI67',"GNLY","NKG7","KLRD1","KLRF1","CD3D","CD3E",'CD3G',"CD4","CD8A",'CD8B',"MS4A1","CD19","CD79A","MZB1","TRBC2",'TRAC'),cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_harmony_choose, features = c("CD14","CD68",'CD163','C1QB','LYZ','TPSAB1','CLEC9A','S100A8','S100A9','CD1C','PTPRC','VCAN'), cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_harmony_choose, features = c("PECAM1","CD34","CDH5","ENG",'PLVAP',"ACTA2","PDGFRB","COL1A1","PDGFRB",'COL2A1','BGN','DCN',"ALB",'AFP',"SOX2",'KRT8'), cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
DotPlot(SCdata_harmony_choose,features = c('PTPRC',"CD3D",'CD3E','CD3G',"CD4","CD8A","CD8B",'TRBC2','TRAC',"GNLY","NKG7",'CD79A','IGHG3','MS4A1',"CD14","CD68",'CD163','C1QB','LYZ','S100A8','MKI67','STMN1','PECAM1','CD34','COL1A1','ACTA2','ALB',"KRT8",'TTR','AMBP','KRT18','KRT19','AFP'),assay='RNA')
p=DotPlot(SCdata_harmony_choose,features = c('PTPRC',"CD3D",'CD3E','CD3G',"CD4","CD8A",'TRBC2','GNLY','NKG7','MS4A1','CD79A','IGHG3','CD68','CD163','LYZ','C1QB','S100A8','PECAM1','COL1A1','ACTA1','ALB','KRT8','AMBP','TTR','AFP'),assay='RNA');p
p=p$data;p=p[which(p$avg.exp > 1 &  p$pct.exp > 25),]


DEGS<-FindAllMarkers(SCdata_harmony_choose,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)


SCdata_harmony_choose$cellR1 = 'CD4T'
SCdata_harmony_choose$cellR1[which(SCdata_harmony_choose$seurat_clusters %in% c(14,17))]<-"CD8T"
SCdata_harmony_choose$cellR1[which(SCdata_harmony_choose$seurat_clusters %in% c(7))]<-"doublet_HepaandPTPRC"
table(SCdata_harmony_choose$cellR1)



SCdata_subset <- subset(SCdata_harmony_choose,seurat_clusters %in% c(7));dim(SCdata_subset)
SCdata <- NormalizeData(SCdata_subset, normalization.method = "LogNormalize", scale.factor = 10000)
SCdata <- FindVariableFeatures(object = SCdata,nfeatures = 2000)
all.genes <- rownames(SCdata)
SCdata <- ScaleData(object = SCdata,features = all.genes)
SCdata <- RunPCA(SCdata, features = VariableFeatures(object = SCdata))
data_use<-Stdev(object = SCdata,reduction = 'pca')
sum(data_use[1:50]);sum(data_use[1:40]);sum(data_use[1:30]);sum(data_use[1:20]);sum(data_use[1:10]) 
dim_usage=50
SCdata_harmony_subset <- SCdata %>%  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SCdata_harmony_subset, 'harmony')
res.used <- c(seq(0.2,8,by=0.2),0.1,0.3,0.5,0.7,0.9);res.used
SCdata_harmony_subset <- SCdata_harmony_subset %>% FindNeighbors(reduction = "harmony", dims = 1:dim_usage) %>% FindClusters(resolution = res.used) %>%
  RunUMAP(reduction = "harmony", dims = 1:dim_usage) %>% RunTSNE(reduction = "harmony", dims = 1:dim_usage) %>% identity()
SCdata_clus.tree.out <- clustree(SCdata_harmony_subset) +theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low = "grey80", high = "red")
print(SCdata_clus.tree.out)
sel.clust = "RNA_snn_res.1"
SCdata_harmony_choose_subset <- SetIdent(SCdata_harmony_subset, value = sel.clust)
SCdata_harmony_choose_subset@meta.data$seurat_clusters <- SCdata_harmony_choose_subset@meta.data[[sel.clust]]
p2 <- DimPlot(SCdata_harmony_choose_subset, reduction = "umap",label.size = 4,label = T);p2 
table(SCdata_harmony_choose_subset@meta.data$seurat_clusters)
double <- data.frame(table(SCdata_harmony_choose_subset@meta.data$seurat_clusters,SCdata_harmony_choose_subset@meta.data$doubletFinder))
colnames(double) <- c('cluster','doubleres','number')
double$sum <- sapply(double$cluster,function(x) length(which(SCdata_harmony_choose_subset@meta.data$seurat_clusters == x))) 
double$ratio <- double$number/double$sum
double[which(double$doubleres == 'Doublet' & double$ratio > 0.1),]
FeaturePlot(object = SCdata_harmony_choose_subset, features = c("PTPRC",'MKI67',"GNLY","NKG7","KLRD1","KLRF1","CD3D","CD3E",'CD3G',"CD4","CD8A","MS4A1","CD19","CD79A","IGHG3","MZB1","TRBC2",'TRAC'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_harmony_choose_subset, features = c("CD14","CD68",'CD163','C1QB','LYZ','TPSAB1','CLEC9A','S100A8','CLEC9A','CLEC10A','TPSAB1','VCAN'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_harmony_choose_subset, features = c("PECAM1","CD34","CDH5","ENG",'PLVAP',"ACTA2","PDGFRB","COL1A1","PDGFRB",'COL2A1','BGN','DCN',"ALB",'AFP',"SOX2",'KRT8'),cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
DotPlot(SCdata_harmony_choose_subset,features = c('PTPRC',"CD3D",'CD3E','CD3G',"CD4","CD8A",'TRBC2','TRAC',"GNLY","NKG7",'CD79A','IGHG3','MS4A1',"CD14","CD68",'CD163','C1QB','LYZ','S100A8','MKI67','STMN1','PECAM1','CD34','COL1A1','ACTA2','ALB',"KRT8",'TTR','AMBP','KRT18','KRT19','AFP'),assay='RNA')
SCdata_harmony_choose_subset$cellR1subet <- 'CD4T'
SCdata_harmony_choose_subset$cellR1subet[which(SCdata_harmony_choose_subset$seurat_clusters %in% c(2))] <- 'Hepa_like'
SCdata_harmony_choose_subset$cellR1subet[which(SCdata_harmony_choose_subset$seurat_clusters %in% c(7,8))] <- 'CD8T'
table(SCdata_harmony_choose_subset$cellR1subet)
#saveRDS(SCdata_harmony_choose_subset,"SCdata_harmony_choose_subset.rds")


SCdata_harmony_choose$cellR1 <- 'CD4T'
SCdata_harmony_choose$cellR1[which(SCdata_harmony_choose$seurat_clusters %in% c(14,17))] <- 'CD8T'
SCdata_harmony_choose$cellR1[names(SCdata_harmony_choose_subset$cellR1subet)] <- SCdata_harmony_choose_subset$cellR1subet
table(SCdata_harmony_choose$cellR1)
DimPlot(SCdata_harmony_choose,group.by = 'cellR1',reduction = "umap",label.size = 4,label = T)
# saveRDS(SCdata_harmony_choose,"SCdata_CD4T_Round1_seurat.rds")
# SCdata_CD8 <- subset(SCdata_harmony_choose,cellR1 %in% c('CD8T'))
# saveRDS(SCdata_CD8,".SCdata_CD4T_Round1_CD8T.rds")
# SCdata_Hepa <- subset(SCdata_harmony_choose,cellR1 %in% c('Hepa_like'))
# saveRDS(SCdata_Hepa,"SCdata_CD4T_Round1_Hepa.rds")




################################################################     2. Round2   Seurat   Process
SCdata_CD4T <- subset(SCdata_harmony_choose,cellR1 == 'CD4T');dim(SCdata_CD4T )
SCdata_CD8T_Round1_CD4T <- readRDS("SCdata_CD8T_Round1_CD4T.rds");dim(SCdata_CD8T_Round1_CD4T)
SCdata <- merge(SCdata_CD4T,SCdata_CD8T_Round1_CD4T)
dim(SCdata)


SCdata <- NormalizeData(SCdata, normalization.method = "LogNormalize", scale.factor = 10000)
SCdata <- FindVariableFeatures(object = SCdata,nfeatures = 2000)
all.genes <- rownames(SCdata)
SCdata <- ScaleData(object = SCdata,features = all.genes)
SCdata <- RunPCA(SCdata, features = VariableFeatures(object = SCdata))
data_use<-Stdev(object = SCdata,reduction = 'pca')
sum(data_use[1:50]);sum(data_use[1:40]);sum(data_use[1:30]);sum(data_use[1:20]);sum(data_use[1:10]) 
ElbowPlot(SCdata,ndims = 50)

dim_usage=50
SCdata_harmony <- SCdata %>%  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SCdata_harmony, 'harmony')
res.used <- c(seq(0.2,8,by=0.2),0.1,0.3,0.5,0.7,0.9);res.used
SCdata_harmony <- SCdata_harmony %>% FindNeighbors(reduction = "harmony", dims = 1:dim_usage) %>% FindClusters(resolution = res.used) %>%
  RunUMAP(reduction = "harmony", dims = 1:dim_usage) %>% RunTSNE(reduction = "harmony", dims = 1:dim_usage) %>% identity()
SCdata_clus.tree.out <- clustree(SCdata_harmony) +theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low = "grey80", high = "red")
print(SCdata_clus.tree.out)

sel.clust = "RNA_snn_res.1"
SCdata_harmony_choose<- SetIdent(SCdata_harmony, value = sel.clust)
SCdata_harmony_choose@meta.data$seurat_clusters <- SCdata_harmony_choose@meta.data[[sel.clust]]
p1 <- DimPlot(SCdata_harmony_choose, reduction = "umap", group.by = "orig.ident")  
p2 <- DimPlot(SCdata_harmony_choose, reduction = "umap",label.size = 4,label = T);p2 
table(SCdata_harmony_choose@meta.data$seurat_clusters)
table(SCdata_harmony_choose$seurat_clusters[which(colnames(SCdata_harmony_choose) %in% colnames(SCdata_CD8T_Round1_CD4T))])

double <- data.frame(table(SCdata_harmony_choose@meta.data$seurat_clusters,SCdata_harmony_choose@meta.data$doubletFinder))
colnames(double) <- c('cluster','doubleres','number')
double$sum <- sapply(double$cluster,function(x) length(which(SCdata_harmony_choose@meta.data$seurat_clusters == x))) 
double$ratio <- double$number/double$sum
double[which(double$doubleres == 'Doublet' & double$ratio > 0.1),]


FeaturePlot(object = SCdata_harmony_choose, features = c("PTPRC",'MKI67',"GNLY","NKG7","KLRD1","KLRF1","CD3D","CD3E",'CD3G',"CD4","CD8A","MS4A1","CD19","CD79A","IGHG3","MZB1","TRBC2",'TRAC'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_harmony_choose, features = c("CD14","CD68",'CD163','C1QB','LYZ','TPSAB1','CLEC9A','S100A8','S100A9','VCAN','CD1C'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_harmony_choose, features = c("PECAM1","CD34","CDH5","ENG",'PLVAP',"ACTA2","PDGFRB","COL1A1","PDGFRB",'COL2A1','BGN','DCN',"ALB",'AFP',"SOX2",'KRT8'), cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
DotPlot(SCdata_harmony_choose,features = c('PTPRC',"CD3D",'CD3E','CD3G',"CD4","CD8A",'CD8B','TRBC2','TRAC',"GNLY","NKG7",'CD79A','IGHG3','MS4A1',"CD14","CD68",'CD163','C1QB','LYZ','S100A8','MKI67','STMN1','PECAM1','CD34','COL1A1','ACTA2','ALB',"KRT8",'TTR','AMBP','KRT18','KRT19','AFP'),assay='RNA')
p=DotPlot(SCdata_harmony_choose,features = c('PTPRC',"CD3D","CD4","CD8A",'TRBC2','GNLY','NKG7','MS4A1','CD79A','IGHG3','CD68','CD163','LYZ','C1QB','PECAM1','COL1A1','ACTA1','ALB','KRT8','AMBP','TTR','AFP'),assay='RNA');p
p=p$data;p=p[which(p$avg.exp > 1 & p$pct.exp > 25),]


DEGS<-FindAllMarkers(SCdata_harmony_choose,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)


SCdata_harmony_choose$cellR1 = 'CD4T'
table(SCdata_harmony_choose$cellR1)
#saveRDS(SCdata_harmony_choose,"SCdata_CD4T_Round2_CD4T.rds")


SCdata_subset <- subset(SCdata_harmony_choose,seurat_clusters %in% c(8,11,12,13,9));dim(SCdata_subset)
SCdata <- NormalizeData(SCdata_subset, normalization.method = "LogNormalize", scale.factor = 10000)
SCdata <- FindVariableFeatures(object = SCdata,nfeatures = 2000)
all.genes <- rownames(SCdata)
SCdata <- ScaleData(object = SCdata,features = all.genes)
SCdata <- RunPCA(SCdata, features = VariableFeatures(object = SCdata))
data_use<-Stdev(object = SCdata,reduction = 'pca')
sum(data_use[1:50]);sum(data_use[1:40]);sum(data_use[1:30]);sum(data_use[1:20]);sum(data_use[1:10]) 
dim_usage=50
SCdata_harmony_subset <- SCdata %>%  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SCdata_harmony_subset, 'harmony')
res.used <- c(seq(0.2,8,by=0.2),0.1,0.3,0.5,0.7,0.9);res.used
SCdata_harmony_subset <- SCdata_harmony_subset %>% FindNeighbors(reduction = "harmony", dims = 1:dim_usage) %>% FindClusters(resolution = res.used) %>%
  RunUMAP(reduction = "harmony", dims = 1:dim_usage) %>% RunTSNE(reduction = "harmony", dims = 1:dim_usage) %>% identity()
SCdata_clus.tree.out <- clustree(SCdata_harmony_subset) +theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low = "grey80", high = "red")
print(SCdata_clus.tree.out)
sel.clust = "RNA_snn_res.5.6"
SCdata_harmony_choose_subset <- SetIdent(SCdata_harmony_subset, value = sel.clust)
SCdata_harmony_choose_subset@meta.data$seurat_clusters <- SCdata_harmony_choose_subset@meta.data[[sel.clust]]
p2 <- DimPlot(SCdata_harmony_choose_subset, reduction = "umap",label.size = 4,label = T);p2 
table(SCdata_harmony_choose_subset@meta.data$seurat_clusters)
double <- data.frame(table(SCdata_harmony_choose_subset@meta.data$seurat_clusters,SCdata_harmony_choose_subset@meta.data$doubletFinder))
colnames(double) <- c('cluster','doubleres','number')
double$sum <- sapply(double$cluster,function(x) length(which(SCdata_harmony_choose_subset@meta.data$seurat_clusters == x))) 
double$ratio <- double$number/double$sum
double[which(double$doubleres == 'Doublet' & double$ratio > 0.1),]
p=DotPlot(SCdata_harmony_choose_subset,features = c('PTPRC',"CD3D",'CD3E','CD3G',"CD4","CD8A",'TRBC2','GNLY','NKG7','MS4A1','CD79A','IGHG3','CD68','CD163','LYZ','C1QB','CLEC9A','CLEC10A','TPSAB1','VCAN','PECAM1','COL1A1','ACTA1','ALB','KRT8','AMBP','TTR','AFP'),assay='RNA');p
p=p$data;p=p[which(p$avg.exp > 1 & p$pct.exp > 25),]
FeaturePlot(object = SCdata_harmony_choose_subset,features = c("PTPRC",'MKI67',"CD3D",'CD3E','CD3G','CD4','CD8A','TRBC2',"GNLY","NKG7","MS4A1","CD79A",'IGHG3','CD68','C1QB','TPSAB1',
                                                               'CLEC9A','S100A8','CD1C','CD163',"PECAM1","CD34","ACTA2","COL1A1","ALB",'AFP',"SOX2",'KRT8'),cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0,label.size = 4)
SCdata_harmony_choose_subset$cellR1subet <- 'CD4T'
SCdata_harmony_choose_subset$cellR1subet[which(SCdata_harmony_choose_subset$seurat_clusters %in% c(28,33))] <- 'CD8T'
SCdata_harmony_choose_subset$cellR1subet[which(SCdata_harmony_choose_subset$seurat_clusters %in% c(29,36))] <- 'doublet-delete'
table(SCdata_harmony_choose_subset$cellR1subet)

SCdata_harmony_choose$cellR1 <- 'CD4T'
SCdata_harmony_choose$cellR1[names(SCdata_harmony_choose_subset$cellR1subet)] <- SCdata_harmony_choose_subset$cellR1subet
table(SCdata_harmony_choose$cellR1)
DimPlot(SCdata_harmony_choose,group.by = 'cellR1',reduction = "umap",label.size = 4,label = T)
#saveRDS(SCdata_harmony_choose,"SCdata_CD4T_Round2_seurat.rds")
# SCdata_CD8 <- subset(SCdata_harmony_choose,cellR1 %in% c('CD8T'))
# saveRDS(SCdata_CD8,"SCdata_CD4T_Round2_CD8T.rds")





################################################################     3. Round4  Seurat   Process
SCdata_CD4T <- subset(SCdata_harmony_choose,cellR1 %in% c('CD4T'));dim(SCdata_CD4T )
SCdata_CD8T_Round2_CD4T <- readRDS("SCdata_CD8T_Round2_CD4T.rds");dim(SCdata_CD8T_Round2_CD4T)
SCdata <- merge(SCdata_CD4T,SCdata_CD8T_Round2_CD4T)
dim(SCdata)

SCdata <- NormalizeData(SCdata, normalization.method = "LogNormalize", scale.factor = 10000)
SCdata <- FindVariableFeatures(object = SCdata,nfeatures = 2000)
all.genes <- rownames(SCdata)
SCdata <- ScaleData(object = SCdata,features = all.genes)
SCdata <- RunPCA(SCdata, features = VariableFeatures(object = SCdata))
data_use<-Stdev(object = SCdata,reduction = 'pca')
sum(data_use[1:50]);sum(data_use[1:40]);sum(data_use[1:30]);sum(data_use[1:20]);sum(data_use[1:10]) 
ElbowPlot(SCdata,ndims = 50)

dim_usage=50
SCdata_harmony <- SCdata %>%  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SCdata_harmony, 'harmony')
res.used <- c(seq(0.2,8,by=0.2),0.1,0.3,0.5,0.7,0.9);res.used
SCdata_harmony <- SCdata_harmony %>% FindNeighbors(reduction = "harmony", dims = 1:dim_usage) %>% FindClusters(resolution = res.used) %>%
  RunUMAP(reduction = "harmony", dims = 1:dim_usage) %>% RunTSNE(reduction = "harmony", dims = 1:dim_usage) %>% identity()
SCdata_clus.tree.out <- clustree(SCdata_harmony) +theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low = "grey80", high = "red")
print(SCdata_clus.tree.out)

sel.clust = "RNA_snn_res.0.6"
SCdata_harmony_choose<- SetIdent(SCdata_harmony, value = sel.clust)
SCdata_harmony_choose@meta.data$seurat_clusters <- SCdata_harmony_choose@meta.data[[sel.clust]]
p1 <- DimPlot(SCdata_harmony_choose, reduction = "umap", group.by = "orig.ident")  
p2 <- DimPlot(SCdata_harmony_choose, reduction = "umap",label.size = 4,label = T);p2 
table(SCdata_harmony_choose@meta.data$seurat_clusters)

double <- data.frame(table(SCdata_harmony_choose@meta.data$seurat_clusters,SCdata_harmony_choose@meta.data$doubletFinder))
colnames(double) <- c('cluster','doubleres','number')
double$sum <- sapply(double$cluster,function(x) length(which(SCdata_harmony_choose@meta.data$seurat_clusters == x))) 
double$ratio <- double$number/double$sum
double[which(double$doubleres == 'Doublet' & double$ratio > 0.1),]

DotPlot(SCdata_harmony_choose,features = c('PTPRC',"CD3D",'CD3E','CD3G',"CD4","CD8A",'CD8B','TRBC2','TRAC',"GNLY","NKG7",'CD79A','IGHG3','MS4A1',"CD14","CD68",'CD163','C1QB','LYZ','S100A8','MKI67','STMN1','PECAM1','CD34','COL1A1','ACTA2','ALB',"KRT8",'TTR','AMBP','KRT18','KRT19','AFP'),assay='RNA')
p=DotPlot(SCdata_harmony_choose,features = c('PTPRC',"CD3D","CD4","CD8A",'TRBC2','GNLY','NKG7','MS4A1','CD79A','IGHG3','CD68','CD163','LYZ','C1QB','PECAM1','COL1A1','ACTA1','ALB','KRT8','AMBP','TTR','AFP'),assay='RNA');p
p=p$data;p=p[which(p$avg.exp > 1 & p$pct.exp > 25),]
FeaturePlot(object = SCdata_harmony_choose,features = c("PTPRC",'MKI67',"CD3D",'CD3E','CD3G','CD4','CD8A','TRBC2',"GNLY","NKG7","MS4A1","CD79A",'IGHG3','CD68','C1QB','TPSAB1','CLEC9A','S100A8','CD1C','CD163',"PECAM1","CD34","ACTA2","COL1A1","ALB",'AFP',"SOX2",'KRT8'),cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0,label.size = 4)

DEGS<-FindAllMarkers(SCdata_harmony_choose,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

SCdata_harmony_choose$cellR1 = 'CD4T'
table(SCdata_harmony_choose$cellR1)


SCdata_harmony_choose$label <- 0
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(0))] <- 'CD4T_NR4A1_BAG3'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(1))] <- 'CD4T_FOXP3'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(2))] <- 'CD4T_GPR183'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(3))] <- 'CD4T_CD40LG'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(4))] <- 'CD4T_GZMA'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(5))] <- 'CD4T_CXCL13'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(6))] <- 'CD4T_SERPINA1'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(7,8,13))] <- 'CD4T_STMN1'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(9))] <- 'CD4T_GZMB'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(10))] <- 'CD4T_IFIT3'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(11))] <- 'CD4T_FOXP3_SERPINA1'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(12))] <- 'CD4T_MT1E'
#saveRDS(SCdata_harmony_choose,"SCdata_CD4T_seurat.rds")



