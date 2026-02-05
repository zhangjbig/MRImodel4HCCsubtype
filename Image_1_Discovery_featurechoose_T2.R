rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()




library(broom)
library(readxl)
library(pROC)
library(irr)
library(dplyr)





##########################################################     Part1:  T2-weight
clinical <- read.csv('ScRNA_6_group_result_finish.csv')
clinical <- data.frame(clinical)
head(clinical)


#######################------------------------------------------------------
T2data <- read_excel('Discovery Cohort_upload.xlsx',sheet = 'Radiomics-FeaturesT2')
T2data <- data.frame(T2data)
colnames(T2data)[-1] <- paste0("T2_",colnames(T2data)[-1])
sample <- T2data$PatientID
group <- sapply(sample,function(x) clinical$groupD[which(clinical$PatientID == x)]);group

T2datause <- data.frame(sample= sample,group=group,T2data)
T2datause <- subset(T2datause,select=-c(PatientID))
rownames(T2datause) <- T2datause$sample
colnames(T2datause)[1:50]
T2datause <- T2datause[,c(1:2,40:dim(T2datause)[2])]
features <- colnames(T2datause)[-c(1,2)]
means <- colMeans(T2datause[,-c(1,2)])
sds <- apply(T2datause[,-c(1,2)], 2, sd)
cv_values <- abs(sds/means);min(cv_values);max(cv_values)
T2cv_report <- data.frame(Feature = features, CV = cv_values,means=means,sd=sds)   
hist(T2cv_report$CV[which(T2cv_report$CV <= 50)], breaks = 100, main = "CV值分布", xlab = "CV值")
hist(T2cv_report$CV[which(T2cv_report$CV <= 10)], breaks = 20, main = "CV值分布", xlab = "CV值")
T2datause_choose <- T2datause[,c("sample","group",T2cv_report$Feature[which(T2cv_report$CV < 1)])]  
dim(T2datause_choose)
saveRDS(T2datause_choose,"T2datause_choose.rds")


#######################-------- Step1:   RandomForest choose features  
library(randomForest)
set.seed(124)
T2datauseRF <- subset(T2datause_choose,select=-(sample))
T2datauseRF$group[which(T2datauseRF$group == 'G1')] <- 1
T2datauseRF$group[which(T2datauseRF$group == 'G2')] <- 0
T2datauseRF$group <- as.factor(T2datauseRF$group)
table(T2datauseRF$group)
dim(T2datauseRF)

T2fit <- randomForest(group~., data = T2datauseRF)
plot(T2fit)
T2importance <-data.frame(importance(T2fit))
T2importance$features <- rownames(T2importance)
T2importance1 <- T2importance %>% top_n(ceiling(n() * 0.05), wt = MeanDecreaseGini)


itera=5000
RF_result <- NULL
for(i in 1:itera){
  T2fit <- randomForest(group~., data = T2datauseRF)
  # plot(T2fit)
  T2importance <-data.frame(importance(T2fit))
  T2importance$features <- rownames(T2importance)
  topn <- T2importance %>% top_n(ceiling(n() * 0.05), wt = MeanDecreaseGini)
  T2importance$Top5per <- ifelse(T2importance$features %in% topn$features,'Top','Untop')
  T2importance$run <- paste0('itera_',i)
  RF_result[[i]] <- list(all=T2importance,topn=topn)
  names(RF_result)[i] <- paste0('itera_',i)
  # saveRDS(T2fit,paste0("T2_RF_model_itera",i,".rds"))
}

pos = seq(1,length(RF_result),1)
RF_summary = lapply(pos,function(x) RF_result[[x]][['topn']])
RF_summary = do.call(rbind,RF_summary)
T2_RF_featureProb = data.frame(table(RF_summary$features))
colnames(T2_RF_featureProb) <- c('featuresT2','choosenum')
T2_RF_featureProb$selPorb = T2_RF_featureProb$choosenum/itera
T2_RF_featurechoose <- as.character(T2_RF_featureProb$featuresT2[which(T2_RF_featureProb$selPorb >= 0.75)])
rm(T2datauseRF,T2fit,itera,i,topn,T2importance,pos)



#######################-------- Step2: RFE
library(caret)
dataMatrix <- T2datause_choose
colnames(dataMatrix)[1:20]
dataMatrix <- subset(dataMatrix,select = -c(sample,group))
dim(dataMatrix) 
group <- sapply(rownames(dataMatrix) ,function(x) clinical$groupD[which(clinical$PatientID == x)])
group <- as.factor(as.character(group));table(group)  


control <- rfeControl(functions = rfFuncs,method = "cv",number = 5)
T2_rfresults <- rfe(dataMatrix,group, sizes = seq(1,dim(dataMatrix)[2],1),rfeControl = control)
T2_rfresults
length(T2_rfresults$optVariables)  
plot(T2_rfresults, type = c("g", "o"))
T2RFE_rffeatures = predictors(T2_rfresults)
T2RFE_rffeatures = T2_rfresults$optVariables

itera=500
RFErf_result <- NULL
for(i in 1:itera){
  rfeControl = rfeControl(functions=rfFuncs,method="cv",saveDetails=T, number=5,allowParallel=T)
  rfProfile <- rfe(dataMatrix,group,sizes = seq(1,dim(dataMatrix)[2],1),rfeControl = rfeControl)
  p1 = ggplot(data = rfProfile, metric = "Accuracy") + theme_bw()
  p2 = ggplot(data = rfProfile, metric = "AccuracySD") + theme_bw()
  cowplot::plot_grid(p1,p2)
  # saveRDS(rfProfile,paste0("T2_RFErf_model_itera",i,".rds"))
  
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
T2_RFErf_featureProb = data.frame(table(RFErf_summary))
colnames(T2_RFErf_featureProb) <- c('featuresT2','choosenum')
T2_RFErf_featureProb$selPorb = T2_RFErf_featureProb$choosenum/itera
T2_RFErf_featurechoose <- as.character(T2_RFErf_featureProb$featuresT2[which(T2_RFErf_featureProb$selPorb >= 0.6)])
rm(dataMatrix,control,T2_rfresults,T2RFE_rffeatures,itera,i,rfeControl,rfProfile,p1,p2,choosefeaturesRFErf,rferfchoose,pos)



#######################-------- Step3:   Logistic regression and Wilcoxtest  choose features  
features <- colnames(T2datause_choose)[-c(1,2)]
ChooseResult_origin <- data.frame()
for(i in 1:length(features)){
  a = T2datause_choose[,c(1,2,which(colnames(T2datause_choose) == features[i]))]
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
T2_Pvalue_featurechoose <- ChooseResult_origin$features[which(ChooseResult_origin$LogPvalue < 0.05 & ChooseResult_origin$Wilpvalue < 0.05)]
T2_Pvalue_featurechoose
saveRDS(ChooseResult_origin,"T2_ChooseResult_origin.rds")
rm(features,i,a,p1,FC,model1,value,result)
