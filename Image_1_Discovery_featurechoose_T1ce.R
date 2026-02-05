rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()




library(broom)
library(readxl)
library(pROC)
library(irr)
library(dplyr)





##########################################################     Part1:  T1CE-weight
clinical <- read.csv('ScRNA_6_group_result_finish.csv')
clinical <- data.frame(clinical)
head(clinical)


#######################------------------------------------------------------
T1CEdata <- read_excel('Discovery Cohort_upload.xlsx',sheet = 'Radiomics-FeaturesT1ce')
T1CEdata <- data.frame(T1CEdata)
colnames(T1CEdata)[-1] <- paste0("T1ce_",colnames(T1CEdata)[-1])
sample <- T1CEdata$PatientID
group <- sapply(sample,function(x) clinical$groupD[which(clinical$PatientID == x)]);group

T1CEdatause <- data.frame(sample= sample,group=group,T1CEdata)
T1CEdatause <- subset(T1CEdatause,select=-c(PatientID))
rownames(T1CEdatause) <- T1CEdatause$sample
colnames(T1CEdatause)[1:50]
T1CEdatause <- T1CEdatause[,c(1:2,40:dim(T1CEdatause)[2])]
features <- colnames(T1CEdatause)[-c(1,2)]
means <- colMeans(T1CEdatause[,-c(1,2)])
sds <- apply(T1CEdatause[,-c(1,2)], 2, sd)
cv_values <- abs(sds/means);min(cv_values);max(cv_values)
T1CEcv_report <- data.frame(Feature = features, CV = cv_values,means=means,sd=sds)   
hist(T1CEcv_report$CV[which(T1CEcv_report$CV <= 50)], breaks = 100, main = "CV值分布", xlab = "CV值")
hist(T1CEcv_report$CV[which(T1CEcv_report$CV <= 10)], breaks = 20, main = "CV值分布", xlab = "CV值")
T1CEdatause_choose <- T1CEdatause[,c("sample","group",T1CEcv_report$Feature[which(T1CEcv_report$CV < 2.5)])]  
dim(T1CEdatause_choose)
saveRDS(T1CEdatause_choose,"T1CEdatause_choose.rds")


#######################-------- Step1:   RandomForest choose features  
library(randomForest)
set.seed(124)
T1CEdatauseRF <- subset(T1CEdatause_choose,select=-(sample))
T1CEdatauseRF$group[which(T1CEdatauseRF$group == 'G1')] <- 1
T1CEdatauseRF$group[which(T1CEdatauseRF$group == 'G2')] <- 0
T1CEdatauseRF$group <- as.factor(T1CEdatauseRF$group)
table(T1CEdatauseRF$group)
dim(T1CEdatauseRF)

T1CEfit <- randomForest(group~., data = T1CEdatauseRF)
plot(T1CEfit)
T1CEimportance <-data.frame(importance(T1CEfit))
T1CEimportance$features <- rownames(T1CEimportance)
T1CEimportance1 <- T1CEimportance %>% top_n(ceiling(n() * 0.05), wt = MeanDecreaseGini)


itera=5000
RF_result <- NULL
for(i in 1:itera){
  T1CEfit <- randomForest(group~., data = T1CEdatauseRF)
  # plot(T1CEfit)
  T1CEimportance <-data.frame(importance(T1CEfit))
  T1CEimportance$features <- rownames(T1CEimportance)
  topn <- T1CEimportance %>% top_n(ceiling(n() * 0.05), wt = MeanDecreaseGini)
  T1CEimportance$Top5per <- ifelse(T1CEimportance$features %in% topn$features,'Top','Untop')
  T1CEimportance$run <- paste0('itera_',i)
  RF_result[[i]] <- list(all=T1CEimportance,topn=topn)
  names(RF_result)[i] <- paste0('itera_',i)
  # saveRDS(T1CEfit,paste0("T1CE_RF_model_itera",i,".rds"))
}

pos = seq(1,length(RF_result),1)
RF_summary = lapply(pos,function(x) RF_result[[x]][['topn']])
RF_summary = do.call(rbind,RF_summary)
T1CE_RF_featureProb = data.frame(table(RF_summary$features))
colnames(T1CE_RF_featureProb) <- c('featuresT1CE','choosenum')
T1CE_RF_featureProb$selPorb = T1CE_RF_featureProb$choosenum/itera
T1CE_RF_featurechoose <- as.character(T1CE_RF_featureProb$featuresT1CE[which(T1CE_RF_featureProb$selPorb >= 0.75)])
rm(T1CEdatauseRF,T1CEfit,itera,i,topn,T1CEimportance,pos)



#######################-------- Step2: RFE
library(caret)
dataMatrix <- T1CEdatause_choose
colnames(dataMatrix)[1:20]
dataMatrix <- subset(dataMatrix,select = -c(sample,group))
dim(dataMatrix) 
group <- sapply(rownames(dataMatrix) ,function(x) clinical$groupD[which(clinical$PatientID == x)])
group <- as.factor(as.character(group));table(group)  


control <- rfeControl(functions = rfFuncs,method = "cv",number = 5)
T1CE_rfresults <- rfe(dataMatrix,group, sizes = seq(1,dim(dataMatrix)[2],1),rfeControl = control)
T1CE_rfresults
length(T1CE_rfresults$optVariables)  
plot(T1CE_rfresults, type = c("g", "o"))
T1CERFE_rffeatures = predictors(T1CE_rfresults)
T1CERFE_rffeatures = T1CE_rfresults$optVariables

itera=500
RFErf_result <- NULL
for(i in 1:itera){
  rfeControl = rfeControl(functions=rfFuncs,method="cv",saveDetails=T, number=5,allowParallel=T)
  rfProfile <- rfe(dataMatrix,group,sizes = seq(1,dim(dataMatrix)[2],1),rfeControl = rfeControl)
  p1 = ggplot(data = rfProfile, metric = "Accuracy") + theme_bw()
  p2 = ggplot(data = rfProfile, metric = "AccuracySD") + theme_bw()
  cowplot::plot_grid(p1,p2)
  # saveRDS(rfProfile,paste0("T1CE_RFErf_model_itera",i,".rds"))
  
  choosefeaturesRFErf = varImp(rfProfile)
  choosefeaturesRFErf = data.frame(features=rownames(choosefeaturesRFErf),Overall=choosefeaturesRFErf$Overall)
  rferfchoose = predictors(rfProfile)
  choosefeaturesRFErf$choose <- ifelse(choosefeaturesRFErf$features %in% rferfchoose,1,0)
  choosefeaturesRFErf$run <- paste0('itera_',i)
  RFErf_result[[i]] <- list(all=choosefeaturesRFErf,topn=rferfchoose)
  names(RFErf_result)[i] <- paste0('itera_',i)
  cat(i,sep="\n")
}

pos = seq(1,length(RFErf_result),1)
RFErf_summary = lapply(pos,function(x) RFErf_result[[x]][['topn']])
RFErf_summary = unlist(RFErf_summary)
T1CE_RFErf_featureProb = data.frame(table(RFErf_summary))
colnames(T1CE_RFErf_featureProb) <- c('featuresT1CE','choosenum')
T1CE_RFErf_featureProb$selPorb = T1CE_RFErf_featureProb$choosenum/itera
T1CE_RFErf_featurechoose <- as.character(T1CE_RFErf_featureProb$featuresT1CE[which(T1CE_RFErf_featureProb$selPorb >= 0.6)])
rm(dataMatrix,control,T1CE_rfresults,T1CERFE_rffeatures,itera,i,rfeControl,rfProfile,p1,p2,choosefeaturesRFErf,rferfchoose,pos)



#######################-------- Step3:   Logistic regression and Wilcoxtest  choose features  
features <- colnames(T1CEdatause_choose)[-c(1,2)]
ChooseResult_origin <- data.frame()
for(i in 1:length(features)){
  a = T1CEdatause_choose[,c(1,2,which(colnames(T1CEdatause_choose) == features[i]))]
  a$groupuse <- ifelse(a$group =='G1',1,0)
  a$groupuse <- factor(a$groupuse)
  colnames(a) <- c('sample','group','featureValue','groupuse')
  p1 = wilcox.test(a$featureValue[which(group == 'G1')],a$featureValue[which(group == 'G2')])
  FC = mean(a$featureValue[which(a$group == 'G1')])/mean(a$featureValue[which(a$group == 'G2')])
  model1 <- glm(groupuse ~ featureValue, data= a, family = binomial())
  value= tidy(model1);p2=value$p.value[2];logcoef = value$estimate[2]
  result = c(mean(a$featureValue[which(a$group == 'G1')]),mean(a$featureValue[which(a$group == 'G2')]),p1$p.value,FC,p2,logcoef)
  ChooseResult_origin <- rbind(ChooseResult_origin,result)
  rownames(ChooseResult_origin)[i] <- features[i] 
  colnames(ChooseResult_origin) <- c('MeanG1','MeanG2','Wilpvalue','FC','LogPvalue','LogCoef')
}
ChooseResult_origin$features <- rownames(ChooseResult_origin)
ChooseResult_origin$Wilqvalue <- p.adjust(ChooseResult_origin$Wilpvalue,method = 'fdr')
ChooseResult_origin$Logqvalue <- p.adjust(ChooseResult_origin$LogPvalue,method = 'fdr')
length(which(ChooseResult_origin$Wilpvalue < 0.05 & ChooseResult_origin$LogPvalue < 0.05))
T1CE_Pvalue_featurechoose <- ChooseResult_origin$features[which(ChooseResult_origin$LogPvalue < 0.05 & ChooseResult_origin$Wilpvalue < 0.05)]
T1CE_Pvalue_featurechoose
saveRDS(ChooseResult_origin,"T1CE_ChooseResult_origin.rds")
rm(features,i,a,p1,FC,model1,value,result)
