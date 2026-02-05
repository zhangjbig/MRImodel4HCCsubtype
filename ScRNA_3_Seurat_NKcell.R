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


# SCdata_All <- readRDS("SCdata_Round1_finish.rds")
dim(SCdata_All);table(SCdata_All$cell)
SCdata <- subset(SCdata_All,cell %in% c('NK'))
dim(SCdata);table(SCdata@meta.data$cell)
table(SCdata@meta.data$orig.ident)



################################################################     1. Round1 
SCdata <- NormalizeData(SCdata, normalization.method = "LogNormalize", scale.factor = 10000)
SCdata <- FindVariableFeatures(object = SCdata,nfeatures = 2000)
all.genes <- rownames(SCdata)
SCdata <- ScaleData(object = SCdata,features = all.genes)
SCdata <- RunPCA(SCdata, features = VariableFeatures(object = SCdata))
data_use<-Stdev(object = SCdata,reduction = 'pca')
sum(data_use[1:50]);sum(data_use[1:40]);sum(data_use[1:30]);sum(data_use[1:20]);sum(data_use[1:10]) 


dim_usage=50
SCdata_harmony <- SCdata %>%  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SCdata_harmony, 'harmony')
res.used <- c(seq(0.2,6,by=0.2),0.1,0.3,0.5,0.7,0.9);res.used
SCdata_harmony <- SCdata_harmony %>% FindNeighbors(reduction = "harmony", dims = 1:dim_usage) %>% FindClusters(resolution = res.used) %>%
  RunUMAP(reduction = "harmony", dims = 1:dim_usage) %>%  RunTSNE(reduction = "harmony", dims = 1:dim_usage) %>%  identity()
SCdata_clus.tree.out <- clustree(SCdata_harmony) +theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low = "grey80", high = "red")
print(SCdata_clus.tree.out)


sel.clust = "RNA_snn_res.1.6"
SCdata_harmony_choose<- SetIdent(SCdata_harmony, value = sel.clust)
SCdata_harmony_choose@meta.data$seurat_clusters <- SCdata_harmony_choose@meta.data[[sel.clust]]
p1 <- DimPlot(SCdata_harmony_choose, reduction = "umap", group.by = "orig.ident")  
p2 <- DimPlot(SCdata_harmony_choose, reduction = "umap",label.size = 4,label = T);p2
table(SCdata_harmony_choose@meta.data$seurat_clusters)

double <- data.frame(table(SCdata_harmony_choose@meta.data$seurat_clusters,SCdata_harmony_choose@meta.data$doubletFinder))
colnames(double) <- c('cluster','doubleres','number')
double$sum <- sapply(double$cluster,function(x) length(which(SCdata_harmony_choose@meta.data$seurat_clusters == x))) 
double$ratio <- double$number/double$sum
double[which(double$doubleres == 'Doublet' & double$ratio >0.1),]


FeaturePlot(object =SCdata_harmony_choose, features = c("PTPRC",'MKI67'), cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0,label.size = 4)
FeaturePlot(object = SCdata_harmony_choose, features = c("GNLY","NKG7","KLRD1","KLRF1","CD3D","CD3E","CD4","CD8A","MS4A1","CD19","CD79A","IGHG3","MZB1","TRBC2",'TRAC'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_harmony_choose, features = c("CD14","CD68",'CD163','C1QB','LYZ','TPSAB1','CLEC9A','S100A8','S100A9'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_harmony_choose, features = c("PECAM1","CD34","CDH5","ENG",'PLVAP',"ACTA2","PDGFRB","COL1A1","PDGFRB",'COL2A1','BGN','DCN'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_harmony_choose, features = c("ALB",'AFP',"SOX2",'KRT8',"AMBP","TTR"),cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0,label.size = 4)
DotPlot(SCdata_harmony_choose,features = c('PTPRC',"CD3D",'CD3E','CD3G',"CD4","CD8A",'TRBC2','TRAC',"GNLY","NKG7",'CD79A','IGHG3','MS4A1',"CD14","CD68",'CD163','C1QB','LYZ','S100A8','MKI67','STMN1','PECAM1','CD34','COL1A1','ACTA2','ALB',"KRT8",'TTR','AMBP','KRT18','KRT19','AFP'),assay='RNA')


hpca.se <- celldex::HumanPrimaryCellAtlasData()
clusters=SCdata_harmony_choose@meta.data$seurat_clusters
pred.hesc <- SingleR(SCdata_harmony_choose@assays$RNA@data, ref = hpca.se, labels = hpca.se$label.main,clusters = clusters)
table(pred.hesc$labels)
celltype = data.frame(ClusterID=rownames(pred.hesc), celltype=pred.hesc$labels, stringsAsFactors = F)
SCdata_harmony_choose@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
p1<-DimPlot(SCdata_harmony_choose,reduction = "umap",label = TRUE,label.size = 4)
p2<-DimPlot(SCdata_harmony_choose, reduction = "umap", group.by = "singleR",label = T,label.size = 5)
plot_grid(p1,p2)


DEGS<-FindAllMarkers(SCdata_harmony_choose,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)


SCdata_harmony_choose$cellR1 <- 'NK'
SCdata_harmony_choose$cellR1[which(SCdata_harmony_choose$seurat_clusters %in% c(0,6,13))] <- 'Tcell_like'
table(SCdata_harmony_choose$cellR1)
SCdata_Tcell <- subset(SCdata_harmony_choose,cellR1 == 'Tcell_like')
# saveRDS(SCdata_Tcell ,"SCdata_NKcell_Round1_Tcell.rds")



################################################################     2. Round2   Seurat   Process
SCdata_Tcell_Round1_NK <- readRDS("SCdata_Tcell_Round1_NK.rds");dim(SCdata_Tcell_Round1_NK)
SCdata_NK <- subset(SCdata_harmony_choose,cellR1 %in% c('NK'))
SCdata <- merge(SCdata_NK,SCdata_Tcell_Round1_NK )
dim(SCdata)

SCdata <- NormalizeData(SCdata, normalization.method = "LogNormalize", scale.factor = 10000)
SCdata <- FindVariableFeatures(object = SCdata,nfeatures = 2000)
all.genes <- rownames(SCdata)
SCdata <- ScaleData(object = SCdata,features = all.genes)
SCdata <- RunPCA(SCdata, features = VariableFeatures(object = SCdata))
data_use<-Stdev(object = SCdata,reduction = 'pca')
sum(data_use[1:50]);sum(data_use[1:40]);sum(data_use[1:30]);sum(data_use[1:20]);sum(data_use[1:10]) 


dim_usage=50
SCdata_harmony <- SCdata %>%  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SCdata_harmony, 'harmony')
res.used <- c(seq(0.2,6,by=0.2),0.1,0.3,0.5,0.7,0.9);res.used
SCdata_harmony <- SCdata_harmony %>% FindNeighbors(reduction = "harmony", dims = 1:dim_usage) %>% FindClusters(resolution = res.used) %>%
  RunUMAP(reduction = "harmony", dims = 1:dim_usage) %>%  RunTSNE(reduction = "harmony", dims = 1:dim_usage) %>%  identity()
SCdata_clus.tree.out <- clustree(SCdata_harmony) +theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low = "grey80", high = "red")
print(SCdata_clus.tree.out)


sel.clust = "RNA_snn_res.1.8"
SCdata_harmony_choose<- SetIdent(SCdata_harmony, value = sel.clust)
SCdata_harmony_choose@meta.data$seurat_clusters <- SCdata_harmony_choose@meta.data[[sel.clust]]
p1 <- DimPlot(SCdata_harmony_choose, reduction = "umap", group.by = "orig.ident")  
p2 <- DimPlot(SCdata_harmony_choose, reduction = "umap",label.size = 4,label = T);p2
table(SCdata_harmony_choose@meta.data$seurat_clusters)

double <- data.frame(table(SCdata_harmony_choose@meta.data$seurat_clusters,SCdata_harmony_choose@meta.data$doubletFinder))
colnames(double) <- c('cluster','doubleres','number')
double$sum <- sapply(double$cluster,function(x) length(which(SCdata_harmony_choose@meta.data$seurat_clusters == x))) 
double$ratio <- double$number/double$sum
double[which(double$doubleres == 'Doublet' & double$ratio >0.1),]


FeaturePlot(object =SCdata_harmony_choose, features = c("PTPRC",'MKI67'), cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0,label.size = 4)
FeaturePlot(object = SCdata_harmony_choose, features = c("GNLY","NKG7","KLRD1","KLRF1","CD3D","CD3E","CD3G","CD4","CD8A","MS4A1","CD19","CD79A","MZB1","TRBC2",'TRAC'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_harmony_choose, features = c("CD14","CD68",'CD163','C1QB','LYZ','TPSAB1','CLEC9A','S100A8','S100A9'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_harmony_choose, features = c("PECAM1","CD34","CDH5","ENG",'PLVAP',"ACTA2","PDGFRB","COL1A1","PDGFRB",'COL2A1','BGN','DCN',"ALB",'AFP',"SOX2",'KRT8'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
DotPlot(SCdata_harmony_choose,features = c('PTPRC',"CD3D",'CD3E',"CD4","CD8A",'TRBC2','TRAC',"GNLY","NKG7","KLRD1","KLRF1",'CD79A','IGHG3','MS4A1',"CD14","CD68",'CD163','C1QB','LYZ','S100A8','MKI67','STMN1','PECAM1','CD34','COL1A1','ACTA2','ALB',"KRT8",'TTR','AMBP','KRT18','KRT19','AFP'),assay='RNA')


DEGS<-FindAllMarkers(SCdata_harmony_choose,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)

SCdata_harmony_choose$cellR1 <- 'NK'
SCdata_harmony_choose$cellR1[which(SCdata_harmony_choose$seurat_clusters %in% c(7,16,17))] <- 'Tcell_like'
table(SCdata_harmony_choose$cellR1)
SCdata_Tcell <- subset(SCdata_harmony_choose,cellR1 == 'Tcell_like')
# saveRDS(SCdata_Tcell ,"SCdata_NKcell_Round2_Tcell.rds")




################################################################     3. Round3   Seurat   Process
SCdata_NKcell <- subset(SCdata_harmony_choose,cellR1 %in% c('NK'));dim(SCdata_NKcell)
SCdata_Tcell_Round2_NKcell <- readRDS("SCdata_Tcell_Round2_NKcell.rds");dim(SCdata_Tcell_Round2_NKcell)
SCdata <- merge(SCdata_NKcell,SCdata_Tcell_Round2_NKcell)
dim(SCdata)

SCdata <- NormalizeData(SCdata, normalization.method = "LogNormalize", scale.factor = 10000)
SCdata <- FindVariableFeatures(object = SCdata,nfeatures = 2000)
all.genes <- rownames(SCdata)
SCdata <- ScaleData(object = SCdata,features = all.genes)
SCdata <- RunPCA(SCdata, features = VariableFeatures(object = SCdata))
data_use<-Stdev(object = SCdata,reduction = 'pca')
sum(data_use[1:50]);sum(data_use[1:40]);sum(data_use[1:30]);sum(data_use[1:20]);sum(data_use[1:10]) 


dim_usage=50
SCdata_harmony <- SCdata %>%  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SCdata_harmony, 'harmony')
res.used <- c(seq(0.2,8,by=0.2),0.1,0.3,0.5,0.7,0.9);res.used
SCdata_harmony <- SCdata_harmony %>% FindNeighbors(reduction = "harmony", dims = 1:dim_usage) %>% FindClusters(resolution = res.used) %>%
  RunUMAP(reduction = "harmony", dims = 1:dim_usage) %>%  RunTSNE(reduction = "harmony", dims = 1:dim_usage) %>%  identity()
SCdata_clus.tree.out <- clustree(SCdata_harmony) +theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low = "grey80", high = "red")
print(SCdata_clus.tree.out)


sel.clust = "RNA_snn_res.2.4"
SCdata_harmony_choose<- SetIdent(SCdata_harmony, value = sel.clust)
SCdata_harmony_choose@meta.data$seurat_clusters <- SCdata_harmony_choose@meta.data[[sel.clust]]
p1 <- DimPlot(SCdata_harmony_choose, reduction = "umap", group.by = "orig.ident")  
p2 <- DimPlot(SCdata_harmony_choose, reduction = "umap",label.size = 4,label = T);p2
table(SCdata_harmony_choose@meta.data$seurat_clusters)
table(SCdata_harmony_choose@meta.data$seurat_clusters[which(colnames(SCdata_harmony_choose) %in% colnames(SCdata_Tcell_Round2_NKcell))])

double <- data.frame(table(SCdata_harmony_choose@meta.data$seurat_clusters,SCdata_harmony_choose@meta.data$doubletFinder))
colnames(double) <- c('cluster','doubleres','number')
double$sum <- sapply(double$cluster,function(x) length(which(SCdata_harmony_choose@meta.data$seurat_clusters == x))) 
double$ratio <- double$number/double$sum
double[which(double$doubleres == 'Doublet' & double$ratio >0.1),]


FeaturePlot(object =SCdata_harmony_choose, features = c("PTPRC",'MKI67'), cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0,label.size = 4)
FeaturePlot(object = SCdata_harmony_choose, features = c("GNLY","NKG7","KLRD1","KLRF1","CD3D","CD3E","CD3G","CD4","CD8A","MS4A1","CD19","CD79A","MZB1","TRBC2",'TRAC'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_harmony_choose, features = c("CD14","CD68",'CD163','C1QB','LYZ','TPSAB1','CLEC9A','S100A8','S100A9'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_harmony_choose, features = c("PECAM1","CD34","CDH5","ENG",'PLVAP',"ACTA2","PDGFRB","COL1A1","PDGFRB",'COL2A1','BGN','DCN',"ALB",'AFP',"SOX2",'KRT8'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
DotPlot(SCdata_harmony_choose,features = c('PTPRC',"CD3D",'CD3E',"CD4","CD8A",'TRBC2','TRAC',"GNLY","NKG7","KLRD1","KLRF1",'CD79A','IGHG3','MS4A1',"CD14","CD68",'CD163','C1QB','LYZ','S100A8','MKI67','STMN1','PECAM1','CD34','COL1A1','ACTA2','ALB',"KRT8",'TTR','AMBP','KRT18','KRT19','AFP'),assay='RNA')


DEGS<-FindAllMarkers(SCdata_harmony_choose,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)


SCdata_harmony_choose$cellR1 <- 'NK'
SCdata_harmony_choose$cellR1[which(SCdata_harmony_choose$seurat_clusters == 20)] <- 'Tcells-like'
table(SCdata_harmony_choose$cellR1)
SCdata_Tcell <- subset(SCdata_harmony_choose,cellR1 == 'Tcells-like')
# saveRDS(SCdata_Tcell ,"Cdata_NKcell_Round3_Tcell.rds")



################################################################     4. Round4  
SCdata_NK <- subset(SCdata_harmony_choose,cellR1 =='NK');dim(SCdata_NK)
SCdata <- SCdata_NK
dim(SCdata)

SCdata <- NormalizeData(SCdata, normalization.method = "LogNormalize", scale.factor = 10000)
SCdata <- FindVariableFeatures(object = SCdata,nfeatures = 2000)
all.genes <- rownames(SCdata)
SCdata <- ScaleData(object = SCdata,features = all.genes)
SCdata <- RunPCA(SCdata, features = VariableFeatures(object = SCdata))
data_use<-Stdev(object = SCdata,reduction = 'pca')
sum(data_use[1:50]);sum(data_use[1:40]);sum(data_use[1:30]);sum(data_use[1:20]);sum(data_use[1:10]) 


dim_usage=50
SCdata_harmony <- SCdata %>%  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SCdata_harmony, 'harmony')
res.used <- c(seq(0.2,8,by=0.2),0.1,0.3,0.5,0.7,0.9);res.used
SCdata_harmony <- SCdata_harmony %>% FindNeighbors(reduction = "harmony", dims = 1:dim_usage) %>% FindClusters(resolution = res.used) %>%
  RunUMAP(reduction = "harmony", dims = 1:dim_usage) %>%  RunTSNE(reduction = "harmony", dims = 1:dim_usage) %>%  identity()
SCdata_clus.tree.out <- clustree(SCdata_harmony) +theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low = "grey80", high = "red")
print(SCdata_clus.tree.out)


sel.clust = "RNA_snn_res.0.4"
SCdata_harmony_choose<- SetIdent(SCdata_harmony, value = sel.clust)
SCdata_harmony_choose@meta.data$seurat_clusters <- SCdata_harmony_choose@meta.data[[sel.clust]]
p1 <- DimPlot(SCdata_harmony_choose, reduction = "umap", group.by = "orig.ident")  
p2 <- DimPlot(SCdata_harmony_choose, reduction = "umap",label.size = 4,label = T);p2
table(SCdata_harmony_choose@meta.data$seurat_clusters)

double <- data.frame(table(SCdata_harmony_choose@meta.data$seurat_clusters,SCdata_harmony_choose@meta.data$doubletFinder))
colnames(double) <- c('cluster','doubleres','number')
double$sum <- sapply(double$cluster,function(x) length(which(SCdata_harmony_choose@meta.data$seurat_clusters == x))) 
double$ratio <- double$number/double$sum
double[which(double$doubleres == 'Doublet' & double$ratio >0.1),]


FeaturePlot(object =SCdata_harmony_choose, features = c("PTPRC",'MKI67'), cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0,label.size = 4)
FeaturePlot(object = SCdata_harmony_choose, features = c("GNLY","NKG7","KLRD1","KLRF1","CD3D","CD3E","CD4","CD8A","MS4A1","CD19","CD79A","MZB1",'IGHG3',"TRBC2",'TRAC'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_harmony_choose, features = c("CD14","CD68",'CD163','C1QB','LYZ','TPSAB1','CLEC9A','S100A8','S100A9'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_harmony_choose, features = c("PECAM1","CD34","CDH5","ENG",'PLVAP',"ACTA2","PDGFRB","COL1A1","PDGFRB",'COL2A1','BGN','DCN',"ALB",'AFP',"SOX2",'KRT8'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
DotPlot(SCdata_harmony_choose,features = c('PTPRC',"CD3D",'CD3E',"CD4","CD8A",'TRBC2','TRAC',"GNLY","NKG7",'CD79A','IGHG3','MS4A1',"CD14","CD68",'CD163','C1QB','LYZ','S100A8','MKI67','STMN1','PECAM1','CD34','COL1A1','ACTA2','ALB',"KRT8",'TTR','AMBP','KRT18','KRT19','AFP'),assay='RNA')
VlnPlot(SCdata_harmony_choose,features = c('PTPRC',"CD3D","CD4","CD8A",'TRBC2','GNLY','NKG7','IGHG3','CD68','CD163','LYZ','C1QB','S100A8','ALB','KRT8','AFP'),assay='RNA')


DEGS<-FindAllMarkers(SCdata_harmony_choose,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)


SCdata_harmony_choose$cellR1 <- 'NK'
table(SCdata_harmony_choose$cellR1)
SCdata_harmony_choose$label <- 0
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(0))] <- 'NK_FCGR3A_CXCR4'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(1))] <- 'NK_CD160_B3GNT7'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(2))] <- 'NK_IFNG'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(3))] <- 'NK_CD160_HSPA6'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(4))] <- 'NK_CD160_IKZF3'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(5))] <- 'NK_FCGR3A_PLAC8'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(6))] <- 'NK_STMN1'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(7))] <- 'NK_SERPINA1'
# saveRDS(SCdata_harmony_choose,"SCdata_NK_seurat.rds")


