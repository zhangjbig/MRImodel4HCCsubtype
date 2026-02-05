rm(list=ls())


library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(harmony)
library(cowplot)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(limma)
library(DOSE)
library(dendextend)


setwd(this.path::this.dir())
getwd()




###########################################################################      Round2-Merge Cells with Anno-label
##################################      step 1. Merge Data-Round2merge
CD4Tcell <- readRDS("SCdata_CD4T_seurat.rds")   
CD8Tcell <- readRDS("SCdata_CD8T_seurat.rds")   
Bcell <- readRDS("SCdata_Bcell_seurat.rds")    
NKcell <- readRDS("SCdata_NK_seurat.rds")   
myeloid <- readRDS("SCdata_myeloid_seurat.rds")  
stromal<- readRDS("SCdata_stromal_seurat.rds")
tumorcell <- readRDS("SCdata_hepatocytes_seurat.rds")  
allfinish1 <- merge(CD4Tcell,y=c(CD8Tcell,Bcell,NKcell,myeloid ,stromal,tumorcell))  
dim(allfinish1);table(allfinish1$cell);table(allfinish1$cellR1)


allfinish1 <- NormalizeData(allfinish1, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(allfinish1)
allfinish1 <- ScaleData(object = allfinish1,features = all.genes)
allfinish1 <- FindVariableFeatures(object = allfinish1,nfeatures = 2000)
allfinish1 <- RunPCA(allfinish1, features = VariableFeatures(object = allfinish1))
#saveRDS(allfinish1,"allfinish1_round2merge_allcells.rds")


allcell_info <- allfinish1@meta.data
allcell_info <- allcell_info[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","doubletFinder","cell","cellR1","sample","tissue")]
allcell_info$cellbarcode = rownames(allcell_info)
allcell_info$cellbarcode = sapply(allcell_info$cellbarcode,function(x) paste0(unlist(strsplit(x,split = '[-_]'))[3:5],collapse = '_'))
# write.csv(allcell_info,"5_Round2_MergeCellWithAnno_information.csv", row.names = FALSE)





###########################################################################      Round2-InferCV : Stromal Ref
table(allfinish1$cellR1)
NonImmune_cell <- subset(allfinish1,cellR1 %in% c('hepatocytes','stromal'))
dim(NonImmune_cell);table(NonImmune_cell$cellR1)


count <- GetAssayData(NonImmune_cell,slot = 'counts')
count <- as.data.frame(count)
dim(count)
write.table(count,file = 'expFileNormal.txt',sep = '\t',quote = F)

geneAnno<-read.table("geneAnno.txt",header = F,sep="\t");geneAnno<-data.frame(geneAnno)
colnames(geneAnno)<-c("SYMBOL","chr","start","end")
geneInfor<-geneAnno[which(geneAnno$SYMBOL %in% rownames(count)),]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
write.table(geneInfor,file = 'geneFileNormal.txt',sep = '\t',quote = F,col.names = F,row.names = F)

groupinfo<-data.frame(v1=colnames(count),v2<-NonImmune_cell@meta.data$cellR1)   
write.table(groupinfo,file = 'groupinfoNormal.txt',sep = '\t',quote = F,col.names = F,row.names = F)

gene <- c("chrMT","chrGL000009.2","chrGL000194.1","chrGL000195.1","chrGL000219.1","chrKI270734.1","chrGL000218.1","chrKI270721.1","chrKI270726.1","chrKI270711.1","chrY")
infercnv_obj = CreateInfercnvObject(raw_counts_matrix= 'expFileNormal.txt',annotations_file='groupinfoNormal.txt',
                                    delim="\t",gene_order_file= 'geneFileNormal.txt',ref_group_names=c("stromal"),chr_exclude = gene )
infercnv_obj2 = infercnv::run(infercnv_obj,cutoff=0.1,out_dir="Result",cluster_by_groups=T,denoise=FALSE,HMM=FALSE,num_threads = 15)



