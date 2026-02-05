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
SCdata <- subset(SCdata_All,cell %in% c('Tcell'))
dim(SCdata);table(SCdata@meta.data$cell)
table(SCdata@meta.data$orig.ident)



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
res.used <- c(seq(0.2,6,by=0.2),0.1,0.3,0.5,0.7,0.9);res.used
SCdata_harmony <- SCdata_harmony %>% FindNeighbors(reduction = "harmony", dims = 1:dim_usage) %>% FindClusters(resolution = res.used) %>%
  RunUMAP(reduction = "harmony", dims = 1:dim_usage) %>% RunTSNE(reduction = "harmony", dims = 1:dim_usage) %>% identity()
SCdata_clus.tree.out <- clustree(SCdata_harmony) +theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low = "grey80", high = "red")
print(SCdata_clus.tree.out)


sel.clust = 'RNA_snn_res.3'
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


FeaturePlot(object =SCdata_harmony_choose, features = c("PTPRC",'MKI67'), cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0,label.size = 4)
FeaturePlot(object = SCdata_harmony_choose, features = c("GNLY","NKG7","KLRD1","KLRF1","CD3D","CD3E","CD4","CD8A","MS4A1","CD19","CD79A","MZB1","TRBC2",'TRAC'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_harmony_choose, features = c("CD14","CD68",'CD163','C1QB','LYZ','TPSAB1','CLEC9A','S100A8','S100A9'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_harmony_choose, features = c("PECAM1","CD34","CDH5","ENG",'PLVAP',"ACTA2","PDGFRB","COL1A1","PDGFRB",'COL2A1','BGN','DCN'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_harmony_choose, features = c("ALB",'AFP',"SOX2",'KRT8',"AMBP","TTR"),cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0,label.size = 4)
DotPlot(SCdata_harmony_choose,features = c('PTPRC',"CD3D",'CD3E',"CD4","CD8A",'TRBC2','TRAC',"GNLY","NKG7",'CD79A','IGHG3','MS4A1',"CD14","CD68",'CD163','C1QB','LYZ','S100A8','MKI67','STMN1','PECAM1','CD34','COL1A1','ACTA2','ALB',"KRT8",'TTR','AMBP','KRT18','KRT19','AFP'),assay='RNA')


hpca.se <- celldex::HumanPrimaryCellAtlasData()
clusters=SCdata_harmony_choose@meta.data$seurat_clusters
pred.hesc <- SingleR(SCdata_harmony_choose@assays$RNA@data, ref = hpca.se, labels = hpca.se$label.main,clusters = clusters)
table(pred.hesc$labels)
celltype = data.frame(ClusterID=rownames(pred.hesc), celltype=pred.hesc$labels, stringsAsFactors = F)
SCdata_harmony_choose@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
p1<-DimPlot(SCdata_harmony_choose,reduction = "umap",label = TRUE,label.size = 8)
p2<-DimPlot(SCdata_harmony_choose, reduction = "umap", group.by = "singleR",label = T,label.size = 5)
plot_grid(p1,p2)


DEGS<-FindAllMarkers(SCdata_harmony_choose,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)


SCdata_harmony_choose$cellR1 = 'Tcell'
SCdata_harmony_choose$cellR1[which(SCdata_harmony_choose$seurat_clusters %in% c(7,40))]<-"NK_like"
SCdata_harmony_choose$cellR1[which(SCdata_harmony_choose$seurat_clusters %in% c(33))]<-"myeloid_like"
SCdata_harmony_choose$cellR1[which(SCdata_harmony_choose$seurat_clusters %in% c(37))]<-"Bcell_like"
SCdata_harmony_choose$cellR1[which(SCdata_harmony_choose$seurat_clusters %in% c(35))]<-"Hepa_like"
SCdata_harmony_choose$cellR1[which(SCdata_harmony_choose$seurat_clusters %in% c(34))]<-"unknow"
table(SCdata_harmony_choose$cellR1)
# SCdata_NKcell <- subset(SCdata_harmony_choose,cellR1 == 'NK')
# SCdata_myeloid <- subset(SCdata_harmony_choose,cellR1 == 'myeloid')
# SCdata_Bcell <- subset(SCdata_harmony_choose,cellR1 == 'Bcell')
# SCdata_hepatocytes <- subset(SCdata_harmony_choose,cellR1 == 'Hepa')
# saveRDS(SCdata_NKcell,"SCdata_Tcell_Round1_NK.rds")
# saveRDS(SCdata_myeloid,"SCdata_Tcell_Round1_myeloid.rds")
# saveRDS(SCdata_Bcell,"SCdata_Tcell_Round1_Bcell.rds")
# saveRDS(SCdata_hepatocytes,"SCdata_Tcell_Round1_hepatocytes.rds")
# saveRDS(SCdata_uncertain,"SCdata_Tcell_Round1_uncertain.rds")




################################################################     2. Round2   Seurat   Process
SCdata_Tcell <- subset(SCdata_harmony_choose,cellR1 == 'Tcell');dim(SCdata_Tcell)
SCdata_NKcell_Round1_Tcell <- readRDS("SCdata_NKcell_Round1_Tcell.rds");dim(SCdata_NKcell_Round1_Tcell )
SCdata_NKcell_Round2_Tcell <- readRDS("SCdata_NKcell_Round2_Tcell.rds");dim(SCdata_NKcell_Round2_Tcell)
SCdata_myeloid_Round2_Tcells <- readRDS("SCdata_myeloid_Round2_Tcells.rds");dim(SCdata_myeloid_Round2_Tcells)
SCdata_Hepa_Round2_Tcell <- readRDS("SCdata_Hepa_Round2_Tcell.rds");dim(SCdata_Hepa_Round2_Tcell)
SCdata <- merge(SCdata_Tcell,y=c(SCdata_NKcell_Round1_Tcell,SCdata_NKcell_Round2_Tcell,SCdata_myeloid_Round2_Tcells,SCdata_Hepa_Round2_Tcell))
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
res.used <- c(seq(0.2,6,by=0.2),0.1,0.3,0.5,0.7,0.9);res.used
SCdata_harmony <- SCdata_harmony %>% FindNeighbors(reduction = "harmony", dims = 1:dim_usage) %>% FindClusters(resolution = res.used) %>%
  RunUMAP(reduction = "harmony", dims = 1:dim_usage) %>% RunTSNE(reduction = "harmony", dims = 1:dim_usage) %>% identity()
SCdata_clus.tree.out <- clustree(SCdata_harmony) +theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low = "grey80", high = "red")
print(SCdata_clus.tree.out)


sel.clust = "RNA_snn_res.0.8"
SCdata_harmony_choose<- SetIdent(SCdata_harmony, value = sel.clust)
SCdata_harmony_choose@meta.data$seurat_clusters <- SCdata_harmony_choose@meta.data[[sel.clust]]
p1 <- DimPlot(SCdata_harmony_choose, reduction = "umap", group.by = "orig.ident")  
p2 <- DimPlot(SCdata_harmony_choose, reduction = "umap",label.size = 4,label = T);p2 
table(SCdata_harmony_choose@meta.data$seurat_clusters)
table(SCdata_harmony_choose@meta.data$seurat_clusters[which(colnames(SCdata_harmony_choose) %in% colnames(SCdata_myeloid_Round2_Tcells))])

double <- data.frame(table(SCdata_harmony_choose@meta.data$seurat_clusters,SCdata_harmony_choose@meta.data$doubletFinder))
colnames(double) <- c('cluster','doubleres','number')
double$sum <- sapply(double$cluster,function(x) length(which(SCdata_harmony_choose@meta.data$seurat_clusters == x))) 
double$ratio <- double$number/double$sum
double[which(double$doubleres == 'Doublet' & double$ratio > 0.1),]


FeaturePlot(object =SCdata_harmony_choose, features = c("PTPRC",'MKI67'), cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0,label.size = 4)
FeaturePlot(object = SCdata_harmony_choose, features = c("GNLY","NKG7","KLRD1","KLRF1","CD3D","CD3E",'CD3G',"CD4","CD8A","MS4A1","CD19","CD79A","IGHG3","MZB1","TRBC2",'TRAC'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_harmony_choose, features = c("CD14","CD68",'CD163','C1QB','LYZ','TPSAB1','CLEC9A','S100A8','S100A9'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_harmony_choose, features = c("PECAM1","CD34","CDH5","ENG",'PLVAP',"ACTA2","PDGFRB","COL1A1","PDGFRB",'COL2A1','BGN','DCN',"ALB",'AFP',"SOX2",'KRT8'), cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
DotPlot(SCdata_harmony_choose,features = c('PTPRC',"CD3D",'CD3E','CD3G',"CD4","CD8A",'TRBC2','TRAC',"GNLY","NKG7",'CD79A','IGHG3','MS4A1',"CD14","CD68",'CD163','C1QB','LYZ','S100A8','MKI67','STMN1','PECAM1','CD34','COL1A1','ACTA2','ALB',"KRT8",'TTR','AMBP','KRT18','KRT19','AFP'),assay='RNA')

SCdata_harmony_choose$cellR1 = 'Tcell'
table(SCdata_harmony_choose$cellR1)
# saveRDS(SCdata_harmony_choose,"SCdata_Tcell_Round2_seurat_beforesubset.rds")


####  CD4T/CD8T-part1
SCdata_subset <- subset(SCdata_harmony_choose,seurat_clusters %in% c(10,13,16,9,14,15))
dim(SCdata_subset)
SCdata <- NormalizeData(SCdata_subset, normalization.method = "LogNormalize", scale.factor = 10000)
SCdata <- FindVariableFeatures(object = SCdata,nfeatures = 2000)
all.genes <- rownames(SCdata)
SCdata <- ScaleData(object = SCdata,features = all.genes)
SCdata <- RunPCA(SCdata, features = VariableFeatures(object = SCdata))
data_use<-Stdev(object = SCdata,reduction = 'pca')
sum(data_use[1:50]);sum(data_use[1:40]);sum(data_use[1:30]);sum(data_use[1:20]);sum(data_use[1:10]) 

dim_usage=40
SCdata_harmony_subset <- SCdata %>%  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SCdata_harmony_subset, 'harmony')
res.used <- c(seq(0.2,6,by=0.2),0.1,0.3,0.5,0.7,0.9);res.used
SCdata_harmony_subset <- SCdata_harmony_subset %>% FindNeighbors(reduction = "harmony", dims = 1:dim_usage) %>% FindClusters(resolution = res.used) %>%
  RunUMAP(reduction = "harmony", dims = 1:dim_usage) %>% RunTSNE(reduction = "harmony", dims = 1:dim_usage) %>% identity()
SCdata_clus.tree.out <- clustree(SCdata_harmony_subset) +theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low = "grey80", high = "red")
print(SCdata_clus.tree.out)
sel.clust = "RNA_snn_res.5"
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
SCdata_harmony_choose_subset$cellR1subet[which(SCdata_harmony_choose_subset$seurat_clusters %in% c(46))] <- 'delete_like_Tcellneg'
SCdata_harmony_choose_subset$cellR1subet[which(SCdata_harmony_choose_subset$seurat_clusters %in% c(24))] <- 'Hepa_like'
SCdata_harmony_choose_subset$cellR1subet[which(SCdata_harmony_choose_subset$seurat_clusters %in% c(41))] <- 'B_like'
SCdata_harmony_choose_subset$cellR1subet[which(SCdata_harmony_choose_subset$seurat_clusters %in% c(12,36,43))] <- 'NK_like'
SCdata_harmony_choose_subset$cellR1subet[which(SCdata_harmony_choose_subset$seurat_clusters %in% c(9,22,34,42))] <- 'myeloid_like'
SCdata_harmony_choose_subset$cellR1subet[which(SCdata_harmony_choose_subset$seurat_clusters %in% c(1,7,8,20,26,30,38,40,44,45))] <- 'CD8T'
table(SCdata_harmony_choose_subset$cellR1subet)
table(SCdata_harmony_choose_subset$cellR1subet[which(colnames(SCdata_harmony_choose_subset) %in% colnames(SCdata_myeloid_Round2_Tcells))])
#saveRDS(SCdata_harmony_choose_subset,"SCdata_harmony_choose_subset1.rds")


####  CD4T/CD8T-part2
SCdata_subset <- subset(SCdata_harmony_choose,seurat_clusters %in% c(8,12))
dim(SCdata_subset)
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
res.used <- c(seq(0.2,6,by=0.2),0.1,0.3,0.5,0.7,0.9);res.used
SCdata_harmony_subset <- SCdata_harmony_subset %>% FindNeighbors(reduction = "harmony", dims = 1:dim_usage) %>% FindClusters(resolution = res.used) %>%
  RunUMAP(reduction = "harmony", dims = 1:dim_usage) %>% RunTSNE(reduction = "harmony", dims = 1:dim_usage) %>% identity()
SCdata_clus.tree.out <- clustree(SCdata_harmony_subset) +theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low = "grey80", high = "red")
print(SCdata_clus.tree.out)

sel.clust = "RNA_snn_res.1.4"
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
p=DotPlot(SCdata_harmony_choose_subset,features = c('PTPRC',"CD3D",'CD3E','CD3G',"CD4","CD8A",'TRBC2','GNLY','NKG7','MS4A1','CD79A','IGHG3','CD68','CD163','LYZ','C1QB','CLEC9A','CLEC10A','TPSAB1','VCAN','PECAM1','COL1A1','ACTA1','ALB','KRT8','AMBP','TTR','AFP'),assay='RNA');p
p=p$data;p=p[which(p$avg.exp > 1 & p$pct.exp > 25),]

SCdata_harmony_choose_subset$cellR1subet <- 'CD8T'
SCdata_harmony_choose_subset$cellR1subet[which(SCdata_harmony_choose_subset$seurat_clusters %in% c(8))] <- 'NK_like'
SCdata_harmony_choose_subset$cellR1subet[which(SCdata_harmony_choose_subset$seurat_clusters %in% c(10,11))] <- 'CD4T'
table(SCdata_harmony_choose_subset$cellR1subet)
#saveRDS(SCdata_harmony_choose_subset,"SCdata_harmony_choose_subset2.rds")


####  CD4T/CD8T-part3
SCdata_harmony_choose_subset1 <- readRDS("SCdata_harmony_choose_subset1.rds")
SCdata_harmony_choose_subset2 <- readRDS("SCdata_harmony_choose_subset2.rds")
table(SCdata_harmony_choose_subset1$cellR1subet);table(SCdata_harmony_choose_subset2$cellR1subet)

SCdata_harmony_choose$cellR1 <- 'Tcells'
SCdata_harmony_choose$cellR1[which(SCdata_harmony_choose$seurat_clusters %in% c(0,1,2,3,11))] <- 'CD4T'
SCdata_harmony_choose$cellR1[which(SCdata_harmony_choose$seurat_clusters %in% c(4,5,6,7))] <- 'CD8T'
SCdata_harmony_choose$cellR1[names(SCdata_harmony_choose_subset1$cellR1subet)] <- SCdata_harmony_choose_subset1$cellR1subet
SCdata_harmony_choose$cellR1[names(SCdata_harmony_choose_subset2$cellR1subet)] <- SCdata_harmony_choose_subset2$cellR1subet
table(SCdata_harmony_choose$cellR1)
DimPlot(SCdata_harmony_choose,group.by = 'cellR1',reduction = "umap",label.size = 4,label = T)
#saveRDS(SCdata_harmony_choose,"SCdata_Tcell_Round2_seurat_Aftersubset.rds")

# SCdata_CD4 <- subset(SCdata_harmony_choose,cellR1 == 'CD4T')
# SCdata_CD8 <- subset(SCdata_harmony_choose,cellR1 %in% c('CD8T'))
# SCdata_B <- subset(SCdata_harmony_choose,cellR1 == 'B_like')
# SCdata_NK <- subset(SCdata_harmony_choose,cellR1 == 'NK_like')
# SCdata_myeloid <- subset(SCdata_harmony_choose,cellR1 == 'myeloid_like')
# SCdata_Hepa <- subset(SCdata_harmony_choose,cellR1 == 'Hepa_like')
# saveRDS(SCdata_B,"SCdata_Tcell_Round2_Bcell.rds")
# saveRDS(SCdata_NK,"SCdata_Tcell_Round2_NKcell.rds")
# saveRDS(SCdata_myeloid,"SCdata_Tcell_Round2_myeloid.rds")
# saveRDS(SCdata_Hepa,"SCdata_Tcell_Round2_Hepa.rds")
# saveRDS(SCdata_CD4,"SCdata_Tcell_Round2_CD4.rds")
# saveRDS(SCdata_CD8,"SCdata_Tcell_Round2_CD8.rds")

