rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()




library(broom)
library(readxl)
library(pROC)
library(irr)
library(dplyr)





##########################################################     Part1:  T1-weight
clinical <- read.csv('ScRNA_6_group_result_finish.csv')
clinical <- data.frame(clinical)
head(clinical)


RFErf_features <- c(T1_RFErf_featurechoose,T2_RFErf_featurechoose,T1CE_RFErf_featurechoose)
RF_features <- c(T1_RF_featurechoose,T2_RF_featurechoose,T1CE_RF_featurechoose)
Pvalue_featurechoose <- c(T1_Pvalue_featurechoose,T2_Pvalue_featurechoose,T1CE_Pvalue_featurechoose)
library(VennDiagram)
venn.plot <- venn.diagram(x = list(RFE=RFErf_features,RF=RF_features,Pchoose = Pvalue_featurechoose),filename = NULL,
                          col = "black",wd = 2,fill = c("#BF1D2D","#397FC7","#F1B656"),cat.col = c("#BF1D2D","#397FC7","#F1B656"),cat.cex = 1.5,cat.fontface = "bold")
grid.draw(venn.plot)
Featureschoose =  intersect(RFErf_features,intersect(RF_features,Pvalue_featurechoose));Featureschoose





#######################------------------------------------------------------
datainput_T1 = readRDS('T1datause_choose.rds')
datainput_T2 = readRDS('T2datause_choose.rds')
datainput_T2 = datainput_T2[,-c(1,2)]
datainput_T1CE = readRDS('T1CEdatause_choose.rds')
datainput_T1CE = datainput_T1CE[,-c(1,2)]
datainput = bind_cols(datainput_T1,datainput_T2,datainput_T1CE)  
dim(datainput)


T1_ChooseResult_origin <- readRDS('T1_ChooseResult_origin.rds')
T2_ChooseResult_origin <- readRDS('T2_ChooseResult_origin.rds')
T1CE_ChooseResult_origin <- readRDS('T1CE_ChooseResult_origin.rds')
ChooseResult_origin <- bind_rows(T1_ChooseResult_origin,T2_ChooseResult_origin,T1CE_ChooseResult_origin)


###   Choose Sig features and calculate AUC value
data_sample <- data.frame(datainput[,c('sample','group',Featureschoose)])
data_sample$ID <- paste0(data_sample$sample,":",data_sample$group)
rownames(data_sample) <- data_sample$ID
data_sample$groupuse <- factor(ifelse(data_sample$group =='G1',1,0))
data_sample <- subset(data_sample,select=-c(sample,group,ID))

combine_predicted_result <- NULL;error <- 0
for(i in 2:length(Featureschoose)){
  combine_i <- combn(Featureschoose,i)
  for(j in 1:dim(combine_i)[2]){
    datause <- data_sample[,c(combine_i[,j],'groupuse')]
    coef_choose <- ChooseResult_origin[which(ChooseResult_origin$features %in% combine_i[,j]),c("features","LogCoef")]
    if(identical(sort(coef_choose$features),sort(combine_i[,j]))){
      datause <- datause[,c('groupuse',coef_choose$features)]
      datause$score <- apply(datause[,-1],1,function(x) coef_choose$LogCoef  %*% x)
      rocresult <- roc(datause$groupuse,datause$score,quiet = TRUE);rocresult
    }else{
      cat(paste0(i,"_",j),sep='\n');error <- error+1
    }
    result <- c(i,paste0(combine_i[,j],collapse = ' ; '),i,rocresult$auc)
    combine_predicted_result  <- rbind(combine_predicted_result,result)
    colnames(combine_predicted_result) <- c('com_Num','com_features','com_Num1','Auc')
  }
}
#write.csv(combine_predicted_result,'combine_predicted_result.csv')
rm(i,combine_i,j,datause,rocresult,result)


data <- combine_predicted_result;data <- data.frame(data)
data <- data[order(data$Auc,decreasing = T),]
data$Auc <- as.numeric(data$Auc)
data$label <- rownames(data);data$label <- factor(data$label,levels = data$label)
ggplot2::ggplot(data=data,aes(x=label,y=Auc))+geom_bar(stat='identity',col='grey',fill='grey',width=0.6)+coord_cartesian(ylim = c(0.7,1))+theme_bw()




#######################------------------------------------------------------
###   choose MAX AUC value and Calculated Score with Features
data$com_features[which(data$Auc == max(data$Auc))]
features <- c('T1_wavelet.HHL_glcm_Imc2','T2_original_shape_Sphericity','T2_wavelet.HLH_firstorder_RootMeanSquared','T1ce_wavelet.HLL_firstorder_Skewness')
data_sample1 <- data_sample[,c(features,'groupuse')]
data_sample1$score <- 0
for(j in 1:length(features)){
  for(i in 1:34){
    datause <- data_sample1[,c(features[j],'groupuse')]
    datause$groupuse <- as.numeric(datause$groupuse)
    model1 <- glm(groupuse ~.  , data=datause, family = binomial())
    summary(model1)
    value= tidy(model1);coef <- value$estimate[-1];coef <- round(coef,2)
    data_sample1$score[i] <- data_sample1$score[i] + coef*datause[i,1]
  }
  cat(identical(round(as.numeric(ChooseResult_origin$LogCoef[which(ChooseResult_origin$features == features[j])]),2),round(coef,2)),sep="\n")
  cat(coef,sep="\n")
}
coef_choose <- ChooseResult_origin[which(ChooseResult_origin$features %in% features),c("features","LogCoef")]
data_sample1$score1 <- apply(data_sample1[,1:4],1,function(x) round(coef_choose$LogCoef,2)  %*% x)
identical(data_sample1$score,data_sample1$score1)
###   Roc
roc(data_sample1$groupuse,data_sample1$score1)
rocobj <- roc(data_sample1$groupuse,data_sample1$score1);rocobj 
plot(rocobj,legacy.axes = TRUE,thresholds="best",print.thres="best",print.auc=TRUE) 
threshold <- 2213.609
rm(i,j,data,datause,model1,value,coef,rocobj)


