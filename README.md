# MRImodel4HCCsubtype
This is a resource for the article: 
"Novel HBV-HCC Molecular Subtypes Coupled with a Non-Invasive MRI Model", 
It includes all code used for the analysis of single cell RNA sequencing data, MRI image data and final prediction model.


Single Cell RNA Sequencing Data Analysis Code:
ScRNA_1. Preprocessing of scRNASeq Data. Including sample merge,quality control,doublet detection, and cell clustering.
ScRNA_2. Cell Cluster Annotation for Main Cell Types. Identify and label the primary mainly cell types in the dataset.
ScRNA_3. Clustering of Each Main Cell Type. Further segment the major cell types into subpopulations.
ScRNA_4. CD4+ T and CD8+ T Cell Clustering.
ScRNA_5. Merging All Cells and Tumor Cell Validation. 
ScRNA_6.Identification of Molecular Subtypes. Identify gene signatures associated with molecular subtypes and assign each sample to its corresponding molecular subtype based on these signature genes.
ScRNA_7. Comparison of Cell Subtype Ratios Across Molecular Subtypes. Examine how the distribution of cell subtypes changes based on molecular classification.


MRI Image Data Analysis Code:
Image_1.Discovery Cohort MRI Data and Feature Selection.  Perform feature selection using machine learning techniques to identify the most relevant features in the MRI data.
Image_2. Building the Optimal Prediction Model.  Use the selected MRI features to train and build thepredictive model.
Image_3. Prediction on New Validation Dataset. Test the model on a validation dataset.


Files Used: including several essential files.
geom_uperrorbar: Function for creating error bars in ggplot visualizations.
major_human_cell_types.csv: Input file for SciBet cell annotation, containing a list of known major human cell types.
geneAnno.txt: Input file used for the inferCNV analysis to annotate gene information.
