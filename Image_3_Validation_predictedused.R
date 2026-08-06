rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()




library(broom)
library(readxl)
library(pROC)
library(dplyr)





##########################################################    
features <- c('T1_wavelet.HHL_glcm_Imc2','T2_original_shape_Sphericity','T2_wavelet.HLH_firstorder_RootMeanSquared','T1ce_wavelet.HLL_firstorder_Skewness')
choosecoef = c(6.32,30.75,7.29,1.51)
threshold = 2213.609




DataVali <- read_excel('./datafile/example.xlsx')
DataVali <- data.frame(DataVali)
PatientID <- DataVali$PatientID
rownames(DataVali) <- DataVali$PatientID
DataVali <- subset(DataVali,select = -c(PatientID))
head(DataVali)

identical(colnames(DataVali),features)

DataVali$score <- apply(DataVali,1,function(x) choosecoef  %*% x)
DataVali$threshold <- 'low'
DataVali$threshold[which(DataVali$score >= threshold)] <- 'high'
table(DataVali$threshold)
rownames(DataVali) <- PatientID

write.csv(DataVali,"./datafile/example_result.csv")
