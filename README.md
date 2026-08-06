###  MRImodel4HCCsubtype

## This is a resource for the article: "Novel HBV-HCC Molecular Subtypes Coupled with a Non-Invasive MRI Model", 

## It includes all code used for the analysis of single cell RNA sequencing data, MRI image data and final prediction model.




###  Environment:
R version: 4.3.0 (2023-04-21)

RStudio : 2023.03.1

Operating system tested: Windows 11 ; Ubuntu 20.04

Required R packages: Seurat v4.4.0 ; DoubletFinder v2.0.3 ; harmony v1.1.0 ; clustree v0.5.0 ; infercnv v1.16.0 ; Cellchat v2.1.2 ; randomForest v4.7-1.1 ; caret v7.0-1; 

Required Python software/package: cellphonedb v5

Note: No non-standard hardware is required.

Install: Install the required R packages from CRAN, Bioconductor, or their
official repositories according to the package documentation.
 Example: install.packages("Seurat")

CellPhoneDB v5 should be installed following its official installation
instructions.


###  File description -Code.R:

## Single Cell RNA Sequencing Data Analysis Code:

ScRNA_1. Preprocessing of scRNASeq Data. Including sample merge,quality control,doublet detection, and cell clustering.

ScRNA_2. Cell Cluster Annotation for Main Cell Types. Identify and label the primary mainly cell types in the dataset.

ScRNA_3. Clustering of Each Main Cell Type. Further segment the major cell types into subpopulations.

ScRNA_4. CD4+ T and CD8+ T Cell Clustering.

ScRNA_5. Merging All Cells and Tumor Cell Validation. 

ScRNA_6.Identification of Molecular Subtypes. Identify gene signatures associated with molecularsubtypes and assign each sample to its corresponding molecular subtype based on these signature genes.

ScRNA_7. Comparison of Cell Subtype Ratios Across Molecular Subtypes. Examine how the distribution of cell subtypes changes based on molecular classification.


## MRI Image Data Analysis Code:

Image_1.Discovery Cohort MRI Data and Feature Selection.  Perform feature selection using machine learning techniques to identify the most relevant features in the MRI data.

Image_2. Building the Optimal Prediction Model.  Use the selected MRI features to train and build thepredictive model.

Image_3. Prediction on New Validation Dataset. Test the model on a validation dataset.

geom_uperrorbar.R: Function for creating error bars in ggplot visualizations.


###  File description - ./datafile/:

major_human_cell_types.csv: Input file for SciBet cell annotation, containing a list of known major human cell types.

geneAnno.txt: Input file used for the inferCNV analysis to annotate gene information.

Discovery_group_result.csv: G1/G2 molecular subtypes for each sample in discovery cohort.

Discovery Cohort_upload.xlsx: The T1, T2 and T1ce-weighted MRI features identified for feature selection and model development in discovery cohort.

example.xlsx: Sample data. used for result prediction.

example_result.csv: The running results of the sample data.


### Prediction

Usage: Open Image_3_Validation_predictedused.R in RStudio. Then run the script or source("Image_3_Validation_predictedused.R").

-- Input: Example input file(. xlsx)  illustrating the required image features data format. Example:"./datafile/example.xlsx"

-- Output: The model-predicted probability of HCC subtype for each sample. Example:"./datafile/example_result.csv"