### Required libraries and utility functions

# Set parent directory as working directory:
getwd()
setwd("..")
getwd() # make sure working directory is the parent directory

# Directory structure:
data_dir <- "./data"
dir.create(data_dir, recursive = T, showWarnings = F)
save_dir <- "./results"
dir.create(save_dir, recursive = T, showWarnings = F)

#### Load libraries ####

library(TCGAbiolinks) # to download data
library(readxl) # te read clinical data
library(miRBaseConverter) # to convert miRNA IDs
library(edgeR) # to normalize expression
library(survival) # to coxPH model
library(survminer) # to coxPH model
library(biospear) # to build lasso model
library(gridExtra) # to export kaplan meier plots
library(forestploter) # to coxPH forestPlots
library(grid) # to coxPH forestPlots
library(modelr) # to CrossValidation patient distribution
library(matrixStats) # to C Index final calculation
library(data.table) # to long format for C index plot 
library(survivalROC) # to AUC ROC curves
library(GEOquery) # to download GEO data for QC
library(sigQC) # to Quality Control of signature
library(ggplot2) # to functional enrichment plots
library(multienrichjam) # to functional enrichment network, install from remotes::install_github("jmw86069/multienrichjam")
library(enrichplot) # to functional enrichment network

source("./src/utils.R")

#----




