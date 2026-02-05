rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()




library(broom)
library(readxl)
library(pROC)
library(irr)
library(dplyr)





##########################################################    
features <- c('T1_wavelet.HHL_glcm_Imc2','T2_original_shape_Sphericity','T2_wavelet.HLH_firstorder_RootMeanSquared','T1ce_wavelet.HLL_firstorder_Skewness')
choosecoef = c(6.32,30.75,7.29,1.51)
threshold = 2213.609




DataVali_T1 <- read_excel('Validation Cohort 7_upload.xlsx',sheet = 'Radiomics-FeaturesT1_1')
DataVali_T1 <- data.frame(DataVali_T1)
colnames(DataVali_T1)[-1] <- paste0('T1_',colnames(DataVali_T1)[-1])
DataVali_T1choose <- DataVali_T1[,c(1,which(colnames(DataVali_T1) == features[1]))]
rownames(DataVali_T1choose) <- DataVali_T1choose$PatientID

DataVali_T2 <- read_excel('Validation Cohort 7_upload.xlsx',sheet = 'Radiomics-FeaturesT2_1')
DataVali_T2 <- data.frame(DataVali_T2)
colnames(DataVali_T2)[-1] <- paste0('T2_',colnames(DataVali_T2)[-1])
DataVali_T2choose <- DataVali_T2[,c(1,which(colnames(DataVali_T2) %in% features[2:3]))]
rownames(DataVali_T2choose) <- DataVali_T2choose$PatientID

DataVali_T1ce <- read_excel('Validation Cohort 7_upload.xlsx',sheet = 'Radiomics-FeaturesT1ce_1')
DataVali_T1ce <- data.frame(DataVali_T1ce)
colnames(DataVali_T1ce)[-1] <- paste0('T1ce_',colnames(DataVali_T1ce)[-1])
DataVali_T1cechoose <- DataVali_T1ce[,c(1,which(colnames(DataVali_T1ce) %in% features[4]))]
rownames(DataVali_T1cechoose) <- DataVali_T1cechoose$PatientID

DataVali <- bind_cols(DataVali_T1choose,DataVali_T2choose,DataVali_T1cechoose)  
DataVali <- DataVali[,which(grepl('^T1_|T2_|T1ce_',colnames(DataVali)))]
DataVali$score <- apply(DataVali,1,function(x) choosecoef  %*% x)
DataVali$threshold <- 'low'
DataVali$threshold[which(DataVali$score >= threshold)] <- 'high'
table(DataVali$threshold)


