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






################################################################     1. Round1   Seurat   Process
SCdata <- readRDS("SCdata_Tcell_Round2_CD8.rds")
dim(SCdata)
SCdata_NKcell_Round3_Tcell <- readRDS("SCdata_NKcell_Round3_Tcell.rds")
SCdata <- merge(SCdata,SCdata_NKcell_Round3_Tcell)
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

sel.clust = 'RNA_snn_res.1.2'
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


FeaturePlot(object = SCdata_harmony_choose, features = c("PTPRC",'MKI67',"GNLY","NKG7","KLRD1","KLRF1","CD3D","CD3E",'CD3G',"CD4","CD8A","MS4A1","CD19","CD79A",'IGHG3',"MZB1","TRBC2",'TRAC'), cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_harmony_choose, features = c("CD14","CD68",'CD163','C1QB','LYZ','TPSAB1','CLEC9A','S100A8','S100A9','CD1C','VCAN'), cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_harmony_choose, features = c("PECAM1","CD34","CDH5","ENG",'PLVAP',"ACTA2","PDGFRB","COL1A1","PDGFRB",'COL2A1','BGN','DCN',"ALB",'AFP',"SOX2",'KRT8'), cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
DotPlot(SCdata_harmony_choose,features = c('PTPRC',"CD3D","CD3E",'CD3G',"CD4","CD8A",'TRBC2','TRAC',"GNLY","NKG7",'CD79A','IGHG3','MS4A1',"CD14","CD68",'CD163','C1QB','LYZ','S100A8','MKI67','STMN1','PECAM1','CD34','COL1A1','ACTA2','ALB',"KRT8",'TTR','AMBP','KRT18','KRT19','AFP'),assay='RNA')
p=DotPlot(SCdata_harmony_choose,features = c('PTPRC',"CD3D",'CD3E','CD3G',"CD4","CD8A",'TRBC2','GNLY','NKG7','MS4A1','CD79A','IGHG3','CD68','CD163','LYZ','C1QB','S100A8','PECAM1','COL1A1','ACTA1','ALB','KRT8','AMBP','TTR','AFP'),assay='RNA');p
p=p$data;p=p[which(p$avg.exp > 1 &  p$pct.exp > 25),]


hpca.se <- celldex::HumanPrimaryCellAtlasData()
clusters=SCdata_harmony_choose@meta.data$seurat_clusters
pred.hesc <- SingleR(SCdata_harmony_choose@assays$RNA@data, ref = hpca.se, labels = hpca.se$label.main,clusters = clusters)
table(pred.hesc$labels)
celltype = data.frame(ClusterID=rownames(pred.hesc), celltype=pred.hesc$labels, stringsAsFactors = F)
SCdata_harmony_choose@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
p1<-DimPlot(SCdata_harmony_choose,reduction = "umap",label = TRUE,label.size = 8)
p2<-DimPlot(SCdata_harmony_choose, reduction = "umap", group.by = "singleR",label = T,label.size = 5)
p1+p2


SCdata_subset <- subset(SCdata_harmony_choose,seurat_clusters %in% c(3,4,9,10,11));dim(SCdata_subset)
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
sel.clust = "RNA_snn_res.2"
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
SCdata_harmony_choose_subset$cellR1subet[which(SCdata_harmony_choose_subset$seurat_clusters %in% c(6,10))] <- 'CD4T'
table(SCdata_harmony_choose_subset$cellR1subet)


SCdata_harmony_choose$cellR1 <- 'CD8T'
SCdata_harmony_choose$cellR1[names(SCdata_harmony_choose_subset$cellR1subet)] <- SCdata_harmony_choose_subset$cellR1subet
table(SCdata_harmony_choose$cellR1)
DimPlot(SCdata_harmony_choose,group.by = 'cellR1',reduction = "umap",label.size = 4,label = T)
# saveRDS(SCdata_harmony_choose,"SCdata_CD8T_Round1_seurat.rds")
# SCdata_CD4 <- subset(SCdata_harmony_choose,cellR1 %in% c('CD4T'))
# saveRDS(SCdata_CD4,"SCdata_CD8T_Round1_CD4T.rds")






################################################################     2. Round2   Seurat   Process
SCdata_CD8 <- subset(SCdata_harmony_choose,cellR1 == 'CD8T');dim(SCdata_CD8)
SCdata_CD4T_Round1_CD8T <- readRDS("SCdata_CD4T_Round1_CD8T.rds")
SCdata <- merge(SCdata_CD8,SCdata_CD4T_Round1_CD8T)
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

sel.clust = "RNA_snn_res.0.5"
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


FeaturePlot(object = SCdata_harmony_choose, features = c("PTPRC",'MKI67',"GNLY","NKG7","KLRD1","KLRF1","CD3D","CD3E","CD4","CD8A","MS4A1","CD19","CD79A","IGHG3","MZB1","TRBC2",'TRAC'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_harmony_choose, features = c("CD14","CD68",'CD163','C1QB','LYZ','TPSAB1','CLEC9A','S100A8','S100A9'),  cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
FeaturePlot(object = SCdata_harmony_choose, features = c("PECAM1","CD34","CDH5","ENG",'PLVAP',"ACTA2","PDGFRB","COL1A1","PDGFRB",'COL2A1','BGN','DCN',"ALB",'AFP',"SOX2",'KRT8'), cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0)
DotPlot(SCdata_harmony_choose,features = c('PTPRC',"CD3D",'CD3E','CD3G',"CD4","CD8A",'TRBC2','TRAC',"GNLY","NKG7",'CD79A','IGHG3','MS4A1',"CD14","CD68",'CD163','C1QB','LYZ','S100A8','MKI67','STMN1','PECAM1','CD34','COL1A1','ACTA2','ALB',"KRT8",'TTR','AMBP','KRT18','KRT19','AFP'),assay='RNA')
p=DotPlot(SCdata_harmony_choose,features = c('PTPRC',"CD3D","CD4","CD8A",'TRBC2','GNLY','NKG7','MS4A1','CD79A','IGHG3','CD68','CD163','LYZ','C1QB','PECAM1','COL1A1','ACTA1','ALB','KRT8','AMBP','TTR','AFP'),assay='RNA');p
p=p$data;p=p[which(p$avg.exp > 1 & p$pct.exp > 25),]


DEGS<-FindAllMarkers(SCdata_harmony_choose,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)


SCdata_subset <- subset(SCdata_harmony_choose,seurat_clusters %in% c(3,9));dim(SCdata_subset)
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
sel.clust = "RNA_snn_res.7"
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
SCdata_harmony_choose_subset$cellR1subet <- 'CD8T'
SCdata_harmony_choose_subset$cellR1subet[which(SCdata_harmony_choose_subset$seurat_clusters %in% c(8,12,14,15))] <- 'CD4T'
table(SCdata_harmony_choose_subset$cellR1subet)


SCdata_harmony_choose$cellR1 <- 'CD8T'
SCdata_harmony_choose$cellR1[names(SCdata_harmony_choose_subset$cellR1subet)] <- SCdata_harmony_choose_subset$cellR1subet
table(SCdata_harmony_choose$cellR1)
DimPlot(SCdata_harmony_choose,group.by = 'cellR1',reduction = "umap",label.size = 4,label = T)
#saveRDS(SCdata_harmony_choose,"SCdata_CD8T_Round2_seurat_aftersubset.rds")
# SCdata_CD4 <- subset(SCdata_harmony_choose,cellR1 %in% c('CD4T'))
# saveRDS(SCdata_CD4,"SCdata_CD8T_Round2_CD4T.rds")




################################################################     3. Round3   Seurat   Process
SCdata_CD8 <- subset(SCdata_harmony_choose,cellR1 == 'CD8T');dim(SCdata_CD8)
SCdata_CD4T_Round2_CD8T <- readRDS("SCdata_CD4T_Round2_CD8T.rds");dim(SCdata_CD4T_Round2_CD8T)
SCdata <- merge(SCdata_CD8,SCdata_CD4T_Round2_CD8T)
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

sel.clust = "RNA_snn_res.0.5"
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


DotPlot(SCdata_harmony_choose,features = c('PTPRC',"CD3D",'CD3E','CD3G',"CD4","CD8A",'TRBC2','TRAC',"GNLY","NKG7",'CD79A','IGHG3','MS4A1',"CD14","CD68",'CD163','C1QB','LYZ','S100A8','MKI67','STMN1','PECAM1','CD34','COL1A1','ACTA2','ALB',"KRT8",'TTR','AMBP','KRT18','KRT19','AFP'),assay='RNA')
p=DotPlot(SCdata_harmony_choose,features = c('PTPRC',"CD3D","CD4","CD8A",'TRBC2','GNLY','NKG7','MS4A1','CD79A','IGHG3','CD68','CD163','LYZ','C1QB','PECAM1','COL1A1','ACTA1','ALB','KRT8','AMBP','TTR','AFP'),assay='RNA');p
p=p$data;p=p[which(p$avg.exp > 1 & p$pct.exp > 25),]
FeaturePlot(object =SCdata_harmony_choose,features = c("PTPRC",'MKI67',"CD3D",'CD3E','CD3G','CD4','CD8A','TRBC2',"GNLY","NKG7","MS4A1","CD79A",'IGHG3','CD68','C1QB','TPSAB1',
                                                       'CLEC9A','S100A8','CD1C','CD163',"PECAM1","CD34","ACTA2","COL1A1","ALB",'AFP',"SOX2",'KRT8'),cols = c("gray", "blue"),max.cutoff = 2,min.cutoff = 0,label.size = 4)


DEGS<-FindAllMarkers(SCdata_harmony_choose,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

SCdata_harmony_choose$cellR1 = 'CD8T'
table(SCdata_harmony_choose$cellR1)


SCdata_harmony_choose$label <- 0
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(0))] <- 'CD8T_NR4A1'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(1))] <- 'CD8T_CXCR4_RGCC'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(2))] <- 'CD8T_IKZF3'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(3))] <- 'CD8T_FGFBP2_FCGR3A'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(4))] <- 'CD8T_SLC4A10'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(5))] <- 'CD8T_CTLA4_LAG3'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(6))] <- 'CD8T_TRDC'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(7))] <- 'CD8T_STMN1'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(8))] <- 'CD8T_SERPINA1'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(9))] <- 'CD8T_CD40LG'
SCdata_harmony_choose$label[which(SCdata_harmony_choose$seurat_clusters %in% c(10))] <- 'CD8T_IFIT3'
#saveRDS(SCdata_harmony_choose,"SCdata_CD8T_seurat.rds")


