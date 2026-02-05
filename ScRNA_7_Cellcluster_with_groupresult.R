rm(list=ls())


library(Seurat)
library(ggplot2)
library(cowplot)
library(ggpubr)


setwd(this.path::this.dir())
getwd()


allfinish1 <- readRDS(allfinish1,"allfinish1_round2merge_allcells.rds")
dim(allfinish1);table(allfinish1$cell);table(allfinish1$cellR1)

group_result <- read.csv('6_group_result_finish.csv')
group_result <- data.frame(group_result)



##############################################################################      calculate cell ratio
sample=unique(allfinish1$orig.ident)
label = as.character(unique(allfinish1$label))
Result_2 <- as.data.frame(matrix(0,nrow = length(label),ncol = length(sample)))
rownames(Result_2) <- label;colnames(Result_2) <- sample
table(allfinish1$cellR1)
allfinish1$cellgroup <- allfinish1$cellR1
allfinish1$cellgroup[which(allfinish1$cellR1 %in% c('Bcell','CD4T','CD8T','NK'))] <- 'TBNK'
table(allfinish1$cellgroup)
for(i in 1:length(sample)){
  for(c in 1:length(label)){
    cellgroup_1 <- unique(allfinish1$cellgroup[which(allfinish1$label == label[c])])
    sum <- length(which(allfinish1$orig.ident == sample[i] & allfinish1$cellgroup == cellgroup_1))
    Result_2[c,i] <- length(which(allfinish1$orig.ident == sample[i] & allfinish1$label == label[c]))/sum
  }
}
Result_2 <- t(Result_2)
Result_2 <- data.frame(Result_2)
Result_2$sampleID <- sapply(rownames(Result_2),function(x) unlist(strsplit(x,split = '[-_]'))[2])
Result_2$group <- sapply(Result_2$sampleID,function(x) group_result$groupD[which(group_result$sampleID == x)])
Result_2 <- Result_2[which(Result_2$group %in% c('G1','G2')),]
# write.csv(Result_2,"7_samplelabelRatio_samplecellType_Result.csv")


data_merge <- Result_2
label = as.character(unique(allfinish1$label))
count1 <- data.frame(label= as.character(),G1G2p = as.numeric(),G1ratio=as.numeric(),G2ratio=as.numeric())
for(i in 1:length(label)){
  pt = t.test(data_merge[which(data_merge$group == 'G1'),colnames(data_merge) == label[i]],data_merge[which(data_merge$group == 'G2'),colnames(data_merge) == label[i]])
  count1[i,] <- c(label[i],pt$p.value,mean(data_merge[which(data_merge$group == 'G1'),colnames(data_merge) == label[i]]),mean(data_merge[which(data_merge$group == 'G2'),colnames(data_merge) == label[i]]))
}
count1$G1G2p <- as.numeric(count1$G1G2p)


source('geom_uperrorbar.R')
for(i in 1:length(label)){
  data <- data_merge[,c(i,81,82)]
  data$group <- factor(data$group,levels = c('G1','G2'))
  colnames(data)[1] <- 'celllabel'
  comparisons <- list(c('G1','G2'))
  p1 = ggplot(data,aes(x=group,y=celllabel,fill=group))+geom_boxplot(outlier.shape = NA)+geom_jitter(width=0.1,size=0.8)+stat_compare_means(comparisons=comparisons,method='t.test')+
    scale_fill_manual(values=c(G2='#C74D26',G1='#308192'))+ylab(label[i])+theme_bw();print(p1)
  p3=ggplot(data,aes(x=group,y=celllabel,fill=group))+geom_bar(stat="summary",fun="mean",position=position_dodge())+stat_summary(aes(col=group),fun.data='mean_sd',geom = "uperrorbar",colour="black",width=0.15)+
    geom_jitter(width = 0.1,size=.8)+ylab(label[i])+stat_compare_means(comparisons=comparisons,method="t.test")+scale_fill_manual(values=c(pang="#CCCCCC",G2='#C74D26',G1='#308192'))+theme_bw();print(p3)
}



data_mergetumor <- data_merge[which(data_merge$group %in% c('G1','G2')),c("Hepa_SERPINC1_HP","Hepa_STMN1","group")]
data_mergetumor$HPSTMN1ratio <- log2(1+data_mergetumor$Hepa_SERPINC1_HP)-log2(data_mergetumor$Hepa_STMN1+1)
mean(data_mergetumor$HPSTMN1ratio[which(data_mergetumor$group == 'G1')]);mean(data_mergetumor$HPSTMN1ratio[which(data_mergetumor$group == 'G2')])
wilcox.test(data_mergetumor$HPSTMN1ratio[which(data_mergetumor$group == 'G1')],data_mergetumor$HPSTMN1ratio[which(data_mergetumor$group == 'G2')])
ggplot(data_mergetumor,aes(x=group,y=HPSTMN1ratio,fill=group))+geom_bar(stat="summary",fun="mean",position=position_dodge(0.5))+stat_summary(aes(col=group),fun.data='mean_sd',geom="uperrorbar",colour="black",width=0.3)+
  geom_jitter(width=0.15,size=2)+ylab('HPSTMN1Ratio')+stat_compare_means(method="wilcox.test")+scale_fill_manual(values=c(G2='#C74D26',G1='#308192'))+theme_bw()


HPSTMN1ratioAUC <- pROC::roc(data_mergetumor$group,data_mergetumor$HPSTMN1ratio);HPSTMN1ratioAUC
plot(HPSTMN1ratioAUC, print.thres=TRUE,print.auc=TRUE)






##############################################################################     ROE Result
divMatrix <- function(m1,m2){
  dim_m1 <- dim(m1)
  dim_m2 <- dim(m2)
  if(sum(dim_m1 == dim_m2) ==2){
    div.result <- matrix(rep(0,dim_m1[1]*dim_m1[2]),nrow  =dim_m1[1])
    rownames(div.result) <- rownames(m1)
    colnames(div.result) <- colnames(m1)
    for(i in 1:dim_m1[1]){
      for(j in 1:dim_m1[2]){
        div.result[i,j] <- m1[i,j]/m2[i,j]
      }
    }
    return(div.result)
  }
}
ROIE <- function(crosstab){
  ## Calculate the Ro/e value from the given crosstab
  ## Args:
  #' @crosstab: the contingency table of given distribution
  ## Return:
  ## The Ro/e matrix
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  return(roie)
}


allfinish1$groupD <- sapply(allfinish1$sample,function(x) group_result$groupD[which(group_result$sampleID == x)])
table(allfinish1$groupD)
tumor = allfinish1@meta.data
tumor <- tumor[,c('cellR1','groupD','orig.ident')]
tumor <- table(tumor$cellR1,tumor$groupD);tumor
roe <- ROIE(tumor)
roe
pheatmap::pheatmap(roe,cluster_rows = T,cluster_cols = T,scale = 'none',display_numbers=T)


