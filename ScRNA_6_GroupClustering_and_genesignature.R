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




###########################################################################      Diff-Genes From ScRNA-seq
allfinish1 <- readRDS(allfinish1,"allfinish1_round2merge_allcells.rds")
dim(allfinish1);table(allfinish1$cell);table(allfinish1$cellR1)


##################################      step 1. First DEGs 
allfinish1$group <- '0'
allfinish1$group[which(allfinish1$tissue == 'pang' & allfinish1$cellR1 == 'stromal')]  <- 'pang'
allfinish1$group[which(allfinish1$tissue == 'tumor' & allfinish1$cellR1 == 'hepatocytes')]  <- 'tumor'
table(allfinish1$group)
DEGS_group1a <- FindMarkers(allfinish1,ident.1 ='tumor',ident.2 = 'pang',group.by = 'group',logfc.threshold = 0.25,min.pct = 0.25)
DEGS_group1chooseballup <- rownames(DEGS_group1a)[which(DEGS_group1a$p_val_adj <0.05 &  DEGS_group1a$avg_log2FC > log2(1.5))]
DEGS_group1chooseballdown <- rownames(DEGS_group1a)[which(DEGS_group1a$p_val_adj <0.05 &  DEGS_group1a$avg_log2FC < log2(1/1.5))]



##################################      step 2. Itera-subsampling DEGs 
sampleresult <- NULL;iteration=300;DEG <- NULL
table(allfinish1$group)
pang1 <- colnames(allfinish1)[which(allfinish1$group == 'pang')]
tumor1 <- colnames(allfinish1)[which(allfinish1$group == 'tumor')]
iteration = 300
for(i in 1:iteration){
  pang_choose <- sample(c(1:length(pang1)),size =  round(0.8*length(pang1)))
  tumor_choose <- sample(c(1:length(tumor1)),size = round(0.8*length(tumor1)))
  allfinish1$groupchoose <- '0'
  allfinish1$groupchoose[which(colnames(allfinish1) %in% pang1[pang_choose])] <- 'pang_choose'
  allfinish1$groupchoose[which(colnames(allfinish1) %in% tumor1[tumor_choose])] <- 'tumor_choose'
  table(allfinish1$groupchoose)
  sampleresult[[2*i-1]] <- c(pang1[pang_choose],tumor1[tumor_choose])
  sampleresult[[2*i]] <- allfinish1$groupchoose
}
for(i in 1:iteration){
  cat(identical(names(sampleresult[[2*i]]),colnames(allfinish1)),sep="\n")   
  allfinish1$groupchoose <- sampleresult[[2*i]]
  table(allfinish1$groupchoose)
  
  DEGS_group1a <- FindMarkers(allfinish1,ident.1 ='tumor_choose',ident.2 = 'pang_choose',group.by = 'groupchoose',logfc.threshold = 0.25,min.pct = 0.25)
  DEGS_group1chooseb <- rownames(DEGS_group1a)[which(DEGS_group1a$p_val_adj <0.05 &  DEGS_group1a$avg_log2FC > log2(1.5))]
  DEG[[i]] <- DEGS_group1a
}
# saveRDS(DEG,'FindMarkers_DEG_300iter_new.rds')
allgeneup <- NULL
for(i in 1:iteration){
  gene <-  DEG[[i]]
  allgeneup <- c(allgeneup,rownames(gene)[which(gene$p_val_adj< 0.05 & gene$avg_log2FC > log2(1.5))])
}
a = data.frame(table(allgeneup))
achoose <- a$allgene[which(a$Freq >= 0.95*iteration)];length(achoose) 
finish_geneup <- intersect(DEGS_group1chooseballup,achoose)
#saveRDS(finish_geneup,'Diffgene_finish_from_ScRNA.rds')




###########################################################################     Round1-clustering
tumorcell <- readRDS("SCdata_hepatocytes_seurat.rds") 
dim(tumorcell)
table(allfinish1$cellR1)
## Protein Data 
data_IBAQ_format <- readRDS("data_IBAQ_formatused.rds")
## sample Info
tumor = unique(as.character(allfinish1$orig.ident[which(allfinish1$tissue == 'tumor')]))
group_result <- data.frame(table(allfinish1$orig.ident,allfinish1$tissue))
group_result <- group_result[which(group_result$Freq > 0),]
colnames(group_result) <- c('sample','tissue','cellnumber')
group_result$sample <- as.character(group_result$sample)
group_result$tissue <- as.character(group_result$tissue)
group_result$sampleID <- sapply(group_result$sample,function(x) unlist(strsplit(x,split = '[-_]'))[2])
group_result$sampleID1 <- paste0("X",group_result$sampleID)
head(group_result)


choosegene <- finish_geneup
expr<-AverageExpression(tumorcell,assays = 'RNA',slot='data',group.by = 'orig.ident')[[1]]    
expr<-expr[rowSums(expr)>0,]
expr<-as.matrix(expr);dim(expr)
choose <- expr[choosegene,tumor]
choose_data_scale<-scale(t(choose));dim(choose_data_scale)

d<-dist(choose_data_scale)
fit1<-hclust(d,method = 'ward.D2')
plot(fit1,hang = -1,cex=.8)
groups <- cutree(fit1, k=2)
rect.hclust(fit1, k=2, border="red")
groups_use <- data.frame(groups)
groups_use$sampleID <-   rownames(groups_use)
groups_use$sampleuse<- sapply(groups_use$sampleID,function(x) unlist(strsplit(x,split = '[-_]'))[2])
table(groups_use$groups)

color = c(seq(-2,0,length=200),seq(0,3,length=200),seq(3,6,length=200))
bk = colorRampPalette(c('green','black','red'))(n=length(unique(color)))
pheatmap::pheatmap(as.matrix(t(choose_data_scale)),clustering_method = "ward.D2",color = bk,breaks = unique(color))


group_result$groupD <- 'pang'
group_result$groupD[which(group_result$sampleID %in% groups_use$sampleuse[which(groups_use$groups == 1)])] <- 'G1'
group_result$groupD[which(group_result$sampleID %in% groups_use$sampleuse[which(groups_use$groups == 2)])] <- 'G2'
table(group_result$groupD)


tumorcell$groupchoose <- '0'
tumorcell$groupchoose[which(tumorcell$sample %in% group_result$sampleID[which(group_result$groupD == 'G1')])] <- 'G1tumor'
tumorcell$groupchoose[which(tumorcell$sample %in% group_result$sampleID[which(group_result$groupD == 'G2')])] <- 'G2tumor'
table(tumorcell$groupchoose)    
DEGS_groupPub1 <- FindMarkers(tumorcell,ident.1 ='G2tumor',ident.2 = 'G1tumor',group.by = 'groupchoose',logfc.threshold = 0.25,min.pct = 0.25)
choosegene_scRNA  <- c(rownames(DEGS_groupPub1)[which(DEGS_groupPub1$p_val_adj < 0.05 & abs(DEGS_groupPub1$avg_log2FC) > log2(1.5))]);log2(1.5);log2(1/1.5)  


tumorsample = group_result$sampleID[which(group_result$groupD %in% c('G1','G2'))]
tumorsample <- intersect(tumorsample,colnames(data_IBAQ_format))
data_IBAQ_choose <- data_IBAQ_format[which(data_IBAQ_format$PG.Genes %in% choosegene_scRNA),c('PG.Genes',tumorsample)]   
data_IBAQ_choose <- data.frame(data_IBAQ_choose);dim(data_IBAQ_choose)   
rownames(data_IBAQ_choose) <- data_IBAQ_choose$PG.Genes
data_IBAQ_choose <- subset(data_IBAQ_choose,select=-c(PG.Genes))
colnames(data_IBAQ_choose);dim(data_IBAQ_choose)  
label <- sapply(colnames(data_IBAQ_choose),function(x) group_result$groupD[which(group_result$sampleID1 == x)]);table(label)
Sig_choose <- data.frame();sum=NULL
for(i in 1:dim(data_IBAQ_choose)[1]){
  data = data_IBAQ_choose[i,]
  if(length(which(data[1,] == 'NaN')) > length(tumorsample)/2){sum = c(sum,rownames(data_IBAQ_choose)[i]);next}
  if(length(which(data[1,] == 'NaN')) < length(tumorsample)/2 & length(which(data[1,] == 'NaN')) > 0){
    data <- data[,-which(data[1,] == 'NaN')]
  }
  label_choose <- sapply(colnames(data),function(x) group_result$groupD[which(group_result$sampleID1 == x)]);table(label_choose)
  p <-t.test(as.numeric(data[1,which(label_choose == 'G1')]),as.numeric(data[1,which(label_choose == 'G2')]))
  a <- c(rownames(data),p$p.value,mean(as.numeric(data[1,which(label_choose == 'G1')])),mean(as.numeric(data[1,which(label_choose == 'G2')])))
  Sig_choose <- rbind(Sig_choose,a)
  colnames(Sig_choose) <- c('Gene','Pvalue','MeanG1','MeanG2')
}
Sig_choose$FC <- as.numeric(Sig_choose$MeanG2)/as.numeric(Sig_choose$MeanG1)
sig <- Sig_choose[which(Sig_choose$Pvalue <0.05),]     
Pub1_finish <- intersect(intersect(sig$Gene,choosegene_scRNA),finish_geneup)
# saveRDS(Pub1_finish,'Protein_SigGene_R1.rds')




###########################################################################     Round2-clustering
choosegene <- Pub1_finish 
expr<-AverageExpression(tumorcell,assays = 'RNA',slot='data',group.by = 'orig.ident')[[1]]    
expr<-expr[rowSums(expr)>0,]
expr<-as.matrix(expr);dim(expr)
choose <- expr[choosegene,tumor]
choose_data_scale<-scale(t(choose));dim(choose_data_scale)

d<-dist(choose_data_scale)
fit1<-hclust(d,method = 'ward.D2')
plot(fit1,hang = -1,cex=.8)
groups <- cutree(fit1, k=2)
rect.hclust(fit1, k=2, border="red")
groups_use <- data.frame(groups)
groups_use$sampleID <-   rownames(groups_use)
groups_use$sampleuse<- sapply(groups_use$sampleID,function(x) unlist(strsplit(x,split = '[-_]'))[2])
table(groups_use$groups)

color = c(seq(-2,0,length=200),seq(0,3,length=200),seq(3,6,length=200))
bk = colorRampPalette(c('green','black','red'))(n=length(unique(color)))
pheatmap::pheatmap(as.matrix(t(choose_data_scale)),clustering_method = "ward.D2",color = bk,breaks = unique(color))


group_result$groupD <- 'pang'
group_result$groupD[which(group_result$sampleID %in% groups_use$sampleuse[which(groups_use$groups == 1)])] <- 'G1'
group_result$groupD[which(group_result$sampleID %in% groups_use$sampleuse[which(groups_use$groups == 2)])] <- 'G2'
table(group_result$groupD)


tumorcell$groupchoose <- '0'
tumorcell$groupchoose[which(tumorcell$sample %in% group_result$sampleID[which(group_result$groupD == 'G1')])] <- 'G1tumor'
tumorcell$groupchoose[which(tumorcell$sample %in% group_result$sampleID[which(group_result$groupD == 'G2')])] <- 'G2tumor'
table(tumorcell$groupchoose)    
DEGS_groupPub1 <- FindMarkers(tumorcell,ident.1 ='G2tumor',ident.2 = 'G1tumor',group.by = 'groupchoose',logfc.threshold = 0.25,min.pct = 0.25)
choosegene_scRNA  <- c(rownames(DEGS_groupPub1)[which(DEGS_groupPub1$p_val_adj < 0.05 & abs(DEGS_groupPub1$avg_log2FC) > log2(1.5))]);log2(1.5);log2(1/1.5)  


tumorsample = group_result$sampleID[which(group_result$groupD %in% c('G1','G2'))]
tumorsample <- intersect(tumorsample,colnames(data_IBAQ_format))
data_IBAQ_choose <- data_IBAQ_format[which(data_IBAQ_format$PG.Genes %in% choosegene_scRNA),c('PG.Genes',tumorsample)]   
data_IBAQ_choose <- data.frame(data_IBAQ_choose);dim(data_IBAQ_choose)   
rownames(data_IBAQ_choose) <- data_IBAQ_choose$PG.Genes
data_IBAQ_choose <- subset(data_IBAQ_choose,select=-c(PG.Genes))
colnames(data_IBAQ_choose);dim(data_IBAQ_choose)  
label <- sapply(colnames(data_IBAQ_choose),function(x) group_result$groupD[which(group_result$sampleID1 == x)]);table(label)
Sig_choose <- data.frame();sum=NULL
for(i in 1:dim(data_IBAQ_choose)[1]){
  data = data_IBAQ_choose[i,]
  if(length(which(data[1,] == 'NaN')) > length(tumorsample)/2){sum = c(sum,rownames(data_IBAQ_choose)[i]);next}
  if(length(which(data[1,] == 'NaN')) < length(tumorsample)/2 & length(which(data[1,] == 'NaN')) > 0){
    data <- data[,-which(data[1,] == 'NaN')]
  }
  label_choose <- sapply(colnames(data),function(x) group_result$groupD[which(group_result$sampleID1 == x)]);table(label_choose)
  p <-t.test(as.numeric(data[1,which(label_choose == 'G1')]),as.numeric(data[1,which(label_choose == 'G2')]))
  a <- c(rownames(data),p$p.value,mean(as.numeric(data[1,which(label_choose == 'G1')])),mean(as.numeric(data[1,which(label_choose == 'G2')])))
  Sig_choose <- rbind(Sig_choose,a)
  colnames(Sig_choose) <- c('Gene','Pvalue','MeanG1','MeanG2')
}
Sig_choose$FC <- as.numeric(Sig_choose$MeanG2)/as.numeric(Sig_choose$MeanG1)
sig <- Sig_choose[which(Sig_choose$Pvalue <0.05),]   
Pub2_finish <- intersect(intersect(sig$Gene,choosegene_scRNA),finish_geneup)
# saveRDS(Pub2_finish,'Protein_SigGene_R2.rds')




###########################################################################     Round3-clustering
choosegene <- Pub2_finish 
expr<-AverageExpression(tumorcell,assays = 'RNA',slot='data',group.by = 'orig.ident')[[1]]    
expr<-expr[rowSums(expr)>0,]
expr<-as.matrix(expr);dim(expr)
choose <- expr[choosegene,tumor]
choose_data_scale<-scale(t(choose));dim(choose_data_scale)

d<-dist(choose_data_scale)
fit1<-hclust(d,method = 'ward.D2')
plot(fit1,hang = -1,cex=.8)
groups <- cutree(fit1, k=2)
rect.hclust(fit1, k=2, border="red")
groups_use <- data.frame(groups)
groups_use$sampleID <-   rownames(groups_use)
groups_use$sampleuse<- sapply(groups_use$sampleID,function(x) unlist(strsplit(x,split = '[-_]'))[2])
table(groups_use$groups)

color = c(seq(-2,0,length=200),seq(0,3,length=200),seq(3,6,length=200))
bk = colorRampPalette(c('green','black','red'))(n=length(unique(color)))
pheatmap::pheatmap(as.matrix(t(choose_data_scale)),clustering_method = "ward.D2",color = bk,breaks = unique(color))


group_result$groupD <- 'pang'
group_result$groupD[which(group_result$sampleID %in% groups_use$sampleuse[which(groups_use$groups == 1)])] <- 'G1'
group_result$groupD[which(group_result$sampleID %in% groups_use$sampleuse[which(groups_use$groups == 2)])] <- 'G2'
table(group_result$groupD)


tumorcell$groupchoose <- '0'
tumorcell$groupchoose[which(tumorcell$sample %in% group_result$sampleID[which(group_result$groupD == 'G1')])] <- 'G1tumor'
tumorcell$groupchoose[which(tumorcell$sample %in% group_result$sampleID[which(group_result$groupD == 'G2')])] <- 'G2tumor'
table(tumorcell$groupchoose)    
DEGS_groupPub1 <- FindMarkers(tumorcell,ident.1 ='G2tumor',ident.2 = 'G1tumor',group.by = 'groupchoose',logfc.threshold = 0.25,min.pct = 0.25)
choosegene_scRNA  <- c(rownames(DEGS_groupPub1)[which(DEGS_groupPub1$p_val_adj < 0.05 & abs(DEGS_groupPub1$avg_log2FC) > log2(1.5))]);log2(1.5);log2(1/1.5)  


tumorsample = group_result$sampleID[which(group_result$groupD %in% c('G1','G2'))]
tumorsample <- intersect(tumorsample,colnames(data_IBAQ_format))
data_IBAQ_choose <- data_IBAQ_format[which(data_IBAQ_format$PG.Genes %in% choosegene_scRNA),c('PG.Genes',tumorsample)]   
data_IBAQ_choose <- data.frame(data_IBAQ_choose);dim(data_IBAQ_choose)   
rownames(data_IBAQ_choose) <- data_IBAQ_choose$PG.Genes
data_IBAQ_choose <- subset(data_IBAQ_choose,select=-c(PG.Genes))
colnames(data_IBAQ_choose);dim(data_IBAQ_choose)  
label <- sapply(colnames(data_IBAQ_choose),function(x) group_result$groupD[which(group_result$sampleID1 == x)]);table(label)
Sig_choose <- data.frame();sum=NULL
for(i in 1:dim(data_IBAQ_choose)[1]){
  data = data_IBAQ_choose[i,]
  if(length(which(data[1,] == 'NaN')) > length(tumorsample)/2){sum = c(sum,rownames(data_IBAQ_choose)[i]);next}
  if(length(which(data[1,] == 'NaN')) < length(tumorsample)/2 & length(which(data[1,] == 'NaN')) > 0){
    data <- data[,-which(data[1,] == 'NaN')]
  }
  label_choose <- sapply(colnames(data),function(x) group_result$groupD[which(group_result$sampleID1 == x)]);table(label_choose)
  p <-t.test(as.numeric(data[1,which(label_choose == 'G1')]),as.numeric(data[1,which(label_choose == 'G2')]))
  a <- c(rownames(data),p$p.value,mean(as.numeric(data[1,which(label_choose == 'G1')])),mean(as.numeric(data[1,which(label_choose == 'G2')])))
  Sig_choose <- rbind(Sig_choose,a)
  colnames(Sig_choose) <- c('Gene','Pvalue','MeanG1','MeanG2')
}
Sig_choose$FC <- as.numeric(Sig_choose$MeanG2)/as.numeric(Sig_choose$MeanG1)
sig <- Sig_choose[which(Sig_choose$Pvalue <0.05),]    
Pub3_finish <- intersect(intersect(sig$Gene,choosegene_scRNA),finish_geneup)
# saveRDS(Pub3_finish,'Protein_SigGene_R3.rds')


identical(Pub2_finish,Pub3_finish)
###  TRUE
# write.csv(group_result,'6_group_result_finish.csv')





