rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()




library(broom)
library(readxl)
library(pROC)
library(irr)
library(dplyr)





##########################################################     Part1:  T1-weight
clinical <- read.csv('./datafile/Discovery_group_result.csv')
clinical <- data.frame(clinical)
head(clinical)


#######################------------------------------------------------------
T1data <- read_excel('./datafile/Discovery Cohort_upload.xlsx',sheet = 'Radiomics-FeaturesT1')
T1data <- data.frame(T1data)
colnames(T1data)[-1] <- paste0("T1_",colnames(T1data)[-1])
sample <- T1data$PatientID
group <- sapply(sample,function(x) clinical$Group_result[which(clinical$PatientID == x)]);group

T1datause <- data.frame(sample= sample,group=group,T1data)
T1datause <- subset(T1datause,select=-c(PatientID))
rownames(T1datause) <- T1datause$sample
colnames(T1datause)[1:50]
T1datause <- T1datause[,c(1:2,40:dim(T1datause)[2])]
features <- colnames(T1datause)[-c(1,2)]
means <- colMeans(T1datause[,-c(1,2)])
sds <- apply(T1datause[,-c(1,2)], 2, sd)
cv_values <- abs(sds/means);min(cv_values);max(cv_values)
T1cv_report <- data.frame(Feature = features, CV = cv_values,means=means,sd=sds)   
hist(T1cv_report$CV[which(T1cv_report$CV <= 50)], breaks = 100, main = "CV值分布", xlab = "CV值")
hist(T1cv_report$CV[which(T1cv_report$CV <= 10)], breaks = 20, main = "CV值分布", xlab = "CV值")
T1datause_choose <- T1datause[,c("sample","group",T1cv_report$Feature[which(T1cv_report$CV < 1.5)])]  
dim(T1datause_choose)
# saveRDS(T1datause_choose,"T1datause_choose.rds")


#######################-------- Step1:   RandomForest choose features  
library(randomForest)
set.seed(124)
T1datauseRF <- subset(T1datause_choose,select=-(sample))
T1datauseRF$group[which(T1datauseRF$group == 'G1')] <- 1
T1datauseRF$group[which(T1datauseRF$group == 'G2')] <- 0
T1datauseRF$group <- as.factor(T1datauseRF$group)
table(T1datauseRF$group)
dim(T1datauseRF)

T1fit <- randomForest(group~., data = T1datauseRF)
plot(T1fit)
T1importance <-data.frame(importance(T1fit))
T1importance$features <- rownames(T1importance)
T1importance1 <- T1importance %>% top_n(ceiling(n() * 0.05), wt = MeanDecreaseGini)


itera=5000
RF_result <- NULL
for(i in 1:itera){
  T1fit <- randomForest(group~., data = T1datauseRF)
  # plot(T1fit)
  T1importance <-data.frame(importance(T1fit))
  T1importance$features <- rownames(T1importance)
  topn <- T1importance %>% top_n(ceiling(n() * 0.05), wt = MeanDecreaseGini)
  T1importance$Top5per <- ifelse(T1importance$features %in% topn$features,'Top','Untop')
  T1importance$run <- paste0('itera_',i)
  RF_result[[i]] <- list(all=T1importance,topn=topn)
  names(RF_result)[i] <- paste0('itera_',i)
  # saveRDS(T1fit,paste0("T1_RF_model_itera",i,".rds"))
}

pos = seq(1,length(RF_result),1)
RF_summary = lapply(pos,function(x) RF_result[[x]][['topn']])
RF_summary = do.call(rbind,RF_summary)
T1_RF_featureProb = data.frame(table(RF_summary$features))
colnames(T1_RF_featureProb) <- c('featuresT1','choosenum')
T1_RF_featureProb$selPorb = T1_RF_featureProb$choosenum/itera
T1_RF_featurechoose <- as.character(T1_RF_featureProb$featuresT1[which(T1_RF_featureProb$selPorb >= 0.75)])
rm(T1datauseRF,T1fit,itera,i,topn,T1importance,pos)



#######################-------- Step2: RFE
library(caret)
dataMatrix <- T1datause_choose
colnames(dataMatrix)[1:20]
dataMatrix <- subset(dataMatrix,select = -c(sample,group))
dim(dataMatrix) 
group <- sapply(rownames(dataMatrix) ,function(x) clinical$Group_result[which(clinical$PatientID == x)])
group <- as.factor(as.character(group));table(group)  


control <- rfeControl(functions = rfFuncs,method = "cv",number = 5)
T1_rfresults <- rfe(dataMatrix,group, sizes = seq(1,dim(dataMatrix)[2],1),rfeControl = control)
T1_rfresults
length(T1_rfresults$optVariables)  
plot(T1_rfresults, type = c("g", "o"))
T1RFE_rffeatures = predictors(T1_rfresults)
T1RFE_rffeatures = T1_rfresults$optVariables

itera=500
RFErf_result <- NULL
for(i in 1:itera){
  rfeControl = rfeControl(functions=rfFuncs,method="cv",saveDetails=T, number=5,allowParallel=T)
  rfProfile <- rfe(dataMatrix,group,sizes = seq(1,dim(dataMatrix)[2],1),rfeControl = rfeControl)
  p1 = ggplot(data = rfProfile, metric = "Accuracy") + theme_bw()
  p2 = ggplot(data = rfProfile, metric = "AccuracySD") + theme_bw()
  cowplot::plot_grid(p1,p2)
  # saveRDS(rfProfile,paste0("T1_RFErf_model_itera",i,".rds"))
  
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
T1_RFErf_featureProb = data.frame(table(RFErf_summary))
colnames(T1_RFErf_featureProb) <- c('featuresT1','choosenum')
T1_RFErf_featureProb$selPorb = T1_RFErf_featureProb$choosenum/itera
T1_RFErf_featurechoose <- as.character(T1_RFErf_featureProb$featuresT1[which(T1_RFErf_featureProb$selPorb >= 0.6)])
rm(dataMatrix,control,T1_rfresults,T1RFE_rffeatures,itera,i,rfeControl,rfProfile,p1,p2,choosefeaturesRFErf,rferfchoose,pos)



#######################-------- Step3:   Logistic regression and Wilcoxtest  choose features  
features <- colnames(T1datause_choose)[-c(1,2)]
ChooseResult_origin <- data.frame()
for(i in 1:length(features)){
  a = T1datause_choose[,c(1,2,which(colnames(T1datause_choose) == features[i]))]
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
T1_Pvalue_featurechoose <- ChooseResult_origin$features[which(ChooseResult_origin$LogPvalue < 0.05 & ChooseResult_origin$Wilpvalue < 0.05)]
T1_Pvalue_featurechoose
saveRDS(ChooseResult_origin,"T1_ChooseResult_origin.rds")
rm(features,i,a,p1,FC,model1,value,result)
