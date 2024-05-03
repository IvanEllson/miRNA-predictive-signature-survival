# miRNA Expression Meta-Analysis

setwd('~/Desktop/Review code/miRNA_meta_analysis')
data_folder <- "/data"
dir.create(paste0(".", data_folder))



### Import data ###

### TARGET-AML


library(TCGAbiolinks)


#### Clinical data OLD: ####
#query_clinical_TARGET_AML <- GDCquery(project = "TARGET-AML",
#                                  data.category = "Clinical",
#                                  data.type = "Clinical Supplement",
#                                  data.format = "xlsx")
#GDCdownload(query_clinical_TARGET_AML)
#
#install.packages("rematch")
#library(readxl)
#
#clinfiles_path <- "./GDCdata/TARGET-AML/harmonized/Clinical/Clinical_Supplement"
#clinfiles <- list.files(path = clinfiles_path, recursive = T)
#
#TARGET_AML_ClinicalData_Validation <- read_excel(paste0(clinfiles_path, "/", clinfiles[1]))
#TARGET_AML_ClinicalData_AML1031 <- read_excel(paste0(clinfiles_path, "/", clinfiles[2]))
#TARGET_AML_ClinicalData_CDE <- read_excel(paste0(clinfiles_path, "/", clinfiles[3])) # data set de diccionario de variables
#TARGET_AML_ClinicalData_Discovery <- read_excel(paste0(clinfiles_path, "/", clinfiles[4]))
#TARGET_AML_ClinicalData_Discovery_and_Validation_Tumor_Content_and_RIN_Supplement <- read_excel(paste0(clinfiles_path, "/", clinfiles[5])) # info de muestras: anatomic site, contenido de tumor, etc
#
#
#library(dplyr)
#dim(TARGET_AML_ClinicalData_Discovery)
#dim(TARGET_AML_ClinicalData_Validation)
#dim(TARGET_AML_ClinicalData_AML1031)
#
#sample_info_clinical <- rbind(TARGET_AML_ClinicalData_Discovery, # juntamos data frames clinicos de pacientes
#                              TARGET_AML_ClinicalData_Validation,
#                              TARGET_AML_ClinicalData_AML1031)
#dim(sample_info_clinical)
#
## Vamos a ver si los duplicados son muestras diferentes del mismo paciente:
#table(duplicated(sample_info_clinical$`TARGET USI`)|duplicated(sample_info_clinical$`TARGET USI`, fromLast = T))
#alldups <- sample_info_clinical[duplicated(sample_info_clinical$`TARGET USI`)|duplicated(sample_info_clinical$`TARGET USI`, fromLast = T),]
#alldups <- alldups[order(alldups$`TARGET USI`),]
#dim(alldups)
#alldups <- alldups[!duplicated.data.frame(alldups),]
#dim(alldups)
#alldups[alldups=="N/A"] <- NA
#alldups <- alldups[!duplicated.data.frame(alldups),]
#dim(alldups)
#alldups <- alldups[duplicated(alldups$`TARGET USI`)|duplicated(alldups$`TARGET USI`, fromLast = T),]
## utils::View(alldups)
#
## Vemos que dentro de los duplicados no hay pacientes duplicados con diferentes datos (muestras diferentes), asi que los quitamos todos:
#sample_info_clinical <- sample_info_clinical[!duplicated(sample_info_clinical$`TARGET USI`),]
#
#TARGET_AML_ClinicalData_Discovery_and_Validation_Tumor_Content_and_RIN_Supplement$`TARGET USI`%in%sample_info_clinical$`TARGET USI`
#table(duplicated(TARGET_AML_ClinicalData_Discovery_and_Validation_Tumor_Content_and_RIN_Supplement$`TARGET USI`)) # dentro del data set de info de muestras no hay pacientes duplicados.
#
## Añadimos info de muestras a nuestro data frame clinico (solo anatomic site):
##sample_info_clinical$'Anatomic Site' <- NA
#dim(sample_info_clinical)
#all_samples_clin_data <- merge(sample_info_clinical, 
#                              TARGET_AML_ClinicalData_Discovery_and_Validation_Tumor_Content_and_RIN_Supplement,
#                              by = "TARGET USI",
#                              all = T)
#dim(all_samples_clin_data)
##utils::View(all_samples_clin_data)
#
#names(all_samples_clin_data)[names(all_samples_clin_data)=='TARGET USI'] <- 'TARGET_USI'
#all_samples_clin_data$TARGET_USI <- gsub("-", "_", all_samples_clin_data$TARGET_USI)


#----


#### Clinical data: ####

### Download

query_clinical_TARGET_AML <- GDCquery(project = "TARGET-AML",
                                      data.category = "Clinical",
                                      data.type = "Clinical Supplement",
                                      data.format = "xlsx")
GDCdownload(query_clinical_TARGET_AML)

###


### Load ###

library(readxl)
clinfiles_path <- "./GDCdata/TARGET-AML/harmonized/Clinical/Clinical_Supplement"
clinfiles <- list.files(path = clinfiles_path, recursive = T)
TARGET_AML_ClinicalData_LowDepthRNAseq <- read_excel(paste0(clinfiles_path, "/", clinfiles[1]))
TARGET_AML_ClinicalData_Validation <- read_excel(paste0(clinfiles_path, "/", clinfiles[2]))
TARGET_AML_ClinicalData_AML1031 <- read_excel(paste0(clinfiles_path, "/", clinfiles[3]))
TARGET_AML_ClinicalData_CDE <- read_excel(paste0(clinfiles_path, "/", clinfiles[4])) # data set de diccionario de variables
TARGET_AML_ClinicalData_Discovery <- read_excel(paste0(clinfiles_path, "/", clinfiles[5]))
TARGET_AML_ClinicalData_AAML1031_AAML0631_additional <- read_excel(paste0(clinfiles_path, "/", clinfiles[6]))
TARGET_AML_ClinicalData_Discovery_and_Validation_Tumor_Content_and_RIN_Supplement <- read_excel(paste0(clinfiles_path, "/", clinfiles[7])) # info de muestras: anatomic site, contenido de tumor, etc

###


### Explore ###

dim(TARGET_AML_ClinicalData_LowDepthRNAseq)
dim(TARGET_AML_ClinicalData_Validation)
dim(TARGET_AML_ClinicalData_AML1031)
dim(TARGET_AML_ClinicalData_Discovery)
dim(TARGET_AML_ClinicalData_AAML1031_AAML0631_additional)
dim(TARGET_AML_ClinicalData_Discovery_and_Validation_Tumor_Content_and_RIN_Supplement)

###


### Modify ###

colnames(TARGET_AML_ClinicalData_LowDepthRNAseq)[colnames(TARGET_AML_ClinicalData_LowDepthRNAseq)=="Gene Fusion...67"] <- "Gene Fusion"

###


### Join data ###

common_cols <- Reduce(intersect, list(colnames(TARGET_AML_ClinicalData_LowDepthRNAseq), # intersect several vectors
                                      colnames(TARGET_AML_ClinicalData_Validation),
                                      colnames(TARGET_AML_ClinicalData_AML1031),
                                      colnames(TARGET_AML_ClinicalData_Discovery),
                                      colnames(TARGET_AML_ClinicalData_AAML1031_AAML0631_additional)))

TARGET_AML_ClinicalData_LowDepthRNAseq <- as.data.frame(subset(TARGET_AML_ClinicalData_LowDepthRNAseq, select = common_cols))
TARGET_AML_ClinicalData_Validation <- as.data.frame(subset(TARGET_AML_ClinicalData_Validation, select = common_cols))
TARGET_AML_ClinicalData_AML1031 <- as.data.frame(subset(TARGET_AML_ClinicalData_AML1031, select = common_cols))
TARGET_AML_ClinicalData_Discovery <- as.data.frame(subset(TARGET_AML_ClinicalData_Discovery, select = common_cols))
TARGET_AML_ClinicalData_AAML1031_AAML0631_additional <- as.data.frame(subset(TARGET_AML_ClinicalData_AAML1031_AAML0631_additional, select = common_cols))

# Create list of data frames:
df_list <- list(TARGET_AML_ClinicalData_LowDepthRNAseq, # juntamos data frames clinicos de pacientes
                TARGET_AML_ClinicalData_Validation,
                TARGET_AML_ClinicalData_AML1031,
                TARGET_AML_ClinicalData_Discovery,
                TARGET_AML_ClinicalData_AAML1031_AAML0631_additional)

# Merge all data frames in list:
sample_info_clinical <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
#utils::View(sample_info_clinical)

###


### Transform

#table(sample_info_clinical$`Gene Fusion`, exclude = NULL)
sample_info_clinical[sample_info_clinical=="NA"] <- NA # Convert "NA" character to NA

###


### Reduce and collapse


# Reduce data frame eliminating duplicated elements based on one column but collapsing different values into one row:
source("~/Desktop/Utility/Useful R codes.R")
sample_info_clinical <- df.reduce_by_id(sample_info_clinical, id_col = "TARGET USI", collapse_by = " ||| ")
#utils::View(sample_info_clinical)
dim(sample_info_clinical)

###



## Vamos a ver si los duplicados muestras diferentes del mismo paciente:
#table(duplicated(sample_info_clinical$`TARGET USI`)|duplicated(sample_info_clinical$`TARGET USI`, fromLast = T))
#alldups <- sample_info_clinical[duplicated(sample_info_clinical$`TARGET USI`)|duplicated(sample_info_clinical$`TARGET USI`, fromLast = T),]
#alldups <- alldups[order(alldups$`TARGET USI`),]
#dim(alldups)
#alldups <- alldups[!duplicated.data.frame(alldups),]
#dim(alldups)
#alldups[alldups=="N/A"] <- NA
#alldups <- alldups[!duplicated.data.frame(alldups),]
#dim(alldups)
#alldups <- alldups[duplicated(alldups$`TARGET USI`)|duplicated(alldups$`TARGET USI`, fromLast = T),]
#utils::View(alldups)
#
## Vemos que dentro de los duplicados no hay pacientes duplicados con diferentes datos (muestras diferentes), asi que los quitamos todos:
#sample_info_clinical <- sample_info_clinical[!duplicated(sample_info_clinical$`TARGET USI`),]
#
#TARGET_AML_ClinicalData_Discovery_and_Validation_Tumor_Content_and_RIN_Supplement$`TARGET USI`%in%sample_info_clinical$`TARGET USI`
#table(duplicated(TARGET_AML_ClinicalData_Discovery_and_Validation_Tumor_Content_and_RIN_Supplement$`TARGET USI`)) # dentro del data set de info de muestras no hay pacientes duplicados.


### Add extra data ###

# Añadimos info de muestras a nuestro data frame clinico (solo analomic site):
#sample_info_clinical$'Anatomic Site' <- NA
dim(sample_info_clinical)
dim(TARGET_AML_ClinicalData_Discovery_and_Validation_Tumor_Content_and_RIN_Supplement)
intersect(sample_info_clinical$`TARGET USI`, TARGET_AML_ClinicalData_Discovery_and_Validation_Tumor_Content_and_RIN_Supplement$`TARGET USI`)
all_samples_clin_data <- merge(sample_info_clinical, 
                               TARGET_AML_ClinicalData_Discovery_and_Validation_Tumor_Content_and_RIN_Supplement,
                               by = "TARGET USI",
                               all = T)
dim(all_samples_clin_data)
#utils::View(all_samples_clin_data)

###


### Modify ###

names(all_samples_clin_data)[names(all_samples_clin_data)=='TARGET USI'] <- 'TARGET_USI'
all_samples_clin_data$TARGET_USI <- gsub("-", "_", all_samples_clin_data$TARGET_USI)

###


### Save ###

write.table(all_samples_clin_data,
            file = paste0(".", data_folder, "/", "clinical_dataset.txt"),
            sep = ";")

###

#----



#### Expression data (immature miRNAs): ####

#query_TARGET_AML <- GDCquery(project = "TARGET-AML",
#                             data.category = "Transcriptome Profiling",
#                             data.type = "miRNA Expression Quantification",
#                             workflow.type = "BCGSC miRNA Profiling")
#GDCdownload(query_TARGET_AML)
#TARGET_miRNA_expr_data <- GDCprepare(query = query_TARGET_AML,
#                                     summarizedExperiment = T)
## *some miRNAs are cross-mapped: miRNA sequence originating from one genomic region is mapped to another location.
#
## Organize data frame:
#counts <- TARGET_miRNA_expr_data[,startsWith(colnames(TARGET_miRNA_expr_data),"read_count_")]
#colnames(counts) <- gsub("read_count_", "", colnames(counts))
#counts_TARGET_AML <- data.frame(miRNA_ID = TARGET_miRNA_expr_data$miRNA_ID,
#                                counts)
#
#rownames(counts_TARGET_AML) <- counts_TARGET_AML$miRNA_ID
#counts_TARGET_AML <- counts_TARGET_AML[,-1]
#
## save:
#write.table(counts_TARGET_AML, paste0(".", data_folder, "/miRNA_expr_counts_immature.txt"), sep = "\t", row.names = T, quote = F)
#
#
## Check TARGET codes at https://ocg.cancer.gov/sites/default/files/OCG-Project-Codes_Tissues-and-Samples_03-02-2021.pdf
## Tenemos muestras de todo tipo: AML, induction failure AML, de Bone marrow, de Peripheral Blood, etc.
## Vamos a hacer un data frame con estos datos:
#
#
#strsplit(colnames(counts_TARGET_AML), '\\.')
#
#seq_samples_clin <- data.frame(sample_code = colnames(counts_TARGET_AML),
#                               project = vapply(strsplit(colnames(counts_TARGET_AML), '\\.'), `[`, 1, FUN.VALUE=character(1)),
#                               tumor_code = vapply(strsplit(colnames(counts_TARGET_AML), '\\.'), `[`, 2, FUN.VALUE=character(1)),
#                               case_identifier = vapply(strsplit(colnames(counts_TARGET_AML), '\\.'), `[`, 3, FUN.VALUE=character(1)),
#                               tissue_code = vapply(strsplit(colnames(counts_TARGET_AML), '\\.'), `[`, 4, FUN.VALUE=character(1)),
#                               nucleic_acid = vapply(strsplit(colnames(counts_TARGET_AML), '\\.'), `[`, 5, FUN.VALUE=character(1)),
#                               TARGET_USI = vapply(strsplit(colnames(counts_TARGET_AML), '\\.'), function(x)
#                                 paste(x[seq.int(3)], collapse='_'), character(1L)))
#
#table(seq_samples_clin$project)
#table(seq_samples_clin$tumor_code)
#table(seq_samples_clin$tissue_code)
#table(seq_samples_clin$nucleic_acid)
#
#
## Vamos a filtrar muestras:
## Nos quedamos con las muestras de de novo AML (tumor_code = 20) y tejido no canceroso (tumor_code = 00):
#seq_samples_clin <- seq_samples_clin[seq_samples_clin$tumor_code == "20" | seq_samples_clin$tumor_code == "00",]
#table(seq_samples_clin$tissue_code)
#
## Sabemos que los tissue_code que tenemos son:
## 03 = Primary blood derived cancer – peripheral blood
## 04 = Recurrent blood derived cancer – bone marrow
## 09 = Primary blood derived cancer – bone marrow
## 14 = Bone marrow normal
## 40 = Recurrent blood derived cancer – peripheral blood
## 50 = Cell line from patient tumor
#
## Nos quedamos con las muestras de Bone Marrow (04A, 09A) y Bone Marrow Normal (14A):
#seq_samples_clin <- seq_samples_clin[seq_samples_clin$tissue_code == "04A" | seq_samples_clin$tissue_code == "09A" | seq_samples_clin$tissue_code == "14A",]
#
#table(seq_samples_clin$tumor_code)
#table(seq_samples_clin$tissue_code)
#
## Tenemos 1509 muestras (42 Recurrent Blood Derived Cancer y 1467 Primary Blood Derived Cancer) y 62 controles.
#
#
## Vamos a ver si hay pacientes duplicados aqui (case_identifier):
#table(duplicated(seq_samples_clin$case_identifier)) # tenemos 34 duplicados
#
## Vamos a ver los duplicados y sus originales:
#
#alldups <- seq_samples_clin[duplicated(seq_samples_clin$case_identifier) | duplicated(seq_samples_clin$case_identifier, fromLast = T),]
#alldups <- alldups[order(alldups$case_identifier),]
##utils::View(alldups)
## Vemos que hay pacientes repetidos porque:
## Tienen muestras de Recurrent blood derived cancer (04) y Primary blood derived cancer (09).
## O bien tienen repetidas las amplificaciones de la misma alicuota (nucleic_acid: 01R, 02R, 03R...)
#
## Vamos a quedarnos con las muestras de Primary blood derived cancer (09) y vamos a quitarnos las Recurrent (04), que solo hay 42 muestras
#
#seq_samples_clin <- seq_samples_clin[!seq_samples_clin$tissue_code == "04A",]
#
#table(seq_samples_clin$tumor_code)
#table(seq_samples_clin$tissue_code)
#table(seq_samples_clin$nucleic_acid)
#
## Vamos a volver a ver duplicados:
#table(duplicated(seq_samples_clin$case_identifier)) # Ya hay solo 4 duplicados
#alldups <- seq_samples_clin[duplicated(seq_samples_clin$case_identifier) | duplicated(seq_samples_clin$case_identifier, fromLast = T),]
#alldups <- alldups[order(alldups$case_identifier),]
##utils::View(alldups)
#
## Vemos que los duplicados son alicuotas repetidas.
## Nos quedamos con la ultima (nucleic_acid mas alto):
#alldups
#
#c <- unique(alldups$case_identifier)
#
#dup_to_remove <- list()
#for (i in 1:length(c)) {
#  dup_to_remove[i] <- alldups[alldups$case_identifier==c[i],][order(alldups[alldups$case_identifier==c[i],]$nucleic_acid, decreasing = T),][2,]$sample_code
#}
#
#dup_to_remove
#dim(seq_samples_clin)
#
#seq_samples_clin <- seq_samples_clin[!seq_samples_clin$sample_code%in%dup_to_remove,]
#dim(seq_samples_clin)
#
#any(duplicated(seq_samples_clin$case_identifier)) # ya no hay pacientes duplicados
#
#
#write.table(seq_samples_clin,
#            file = paste0(".", data_folder, "/", "clinical_dataset_2.txt"),
#            sep = ";")


#----




#### Expression data (mature miRNAs): ####

library(TCGAbiolinks)
query_TARGET_AML <- GDCquery(project = "TARGET-AML",
                             data.category = "Transcriptome Profiling",
                             data.type = "Isoform Expression Quantification",
                             workflow.type = "BCGSC miRNA Profiling")
GDCdownload(query_TARGET_AML)
TARGET_miRNA_expr_data_iso <- GDCprepare(query = query_TARGET_AML,
                                     summarizedExperiment = T)
# *some miRNAs are cross-mapped: miRNA sequence originating from one genomic region is mapped to another location.

#write.table(TARGET_miRNA_expr_data_iso, paste0(".", data_folder, "/TARGET_miRNA_expr_data_iso.txt"), sep = "\t", row.names = F, quote = F)


miRNA_expr_long <- TARGET_miRNA_expr_data_iso
#table(gsub(",MIMAT.*", "", miRNA_expr_long$miRNA_region), exclude = NULL) # check types of miRNAs
miRNA_expr_long <- miRNA_expr_long[startsWith(miRNA_expr_long$miRNA_region, "mature"),] # Get rid of immature miRNAs
miRNA_expr_long$accession <- gsub("mature,", "", miRNA_expr_long$miRNA_region)
#table(miRNA_expr_long$barcode) # different miRNA isoforms per sample

miRNA_expr_counts <- data.frame("accession" = na.omit(unique(miRNA_expr_long$accession)))
library("miRBaseConverter")
miRversion <- checkMiRNAVersion(unique(miRNA_expr_long$miRNA_ID)) # Check miRbase version = v21
miRNA_expr_counts$miRNA_ID_mature <- miRNA_AccessionToName(miRNA_expr_counts$accession, targetVersion = miRversion)[,2]
library(tibble)
miRNA_expr_counts <- column_to_rownames(miRNA_expr_counts, "accession")


for (i in 1:length(unique(miRNA_expr_long$barcode))) {
  patient <- miRNA_expr_long[miRNA_expr_long$barcode==unique(miRNA_expr_long$barcode)[i],c("accession", "read_count")]
  patient_uniqueAccession <- aggregate(read_count~accession, patient, sum) # collapse to unique accession by sum
  miRNA_expr_counts[,unique(miRNA_expr_long$barcode)[i]] <- 0 # Not measured as 0 counts
  miRNA_expr_counts[patient_uniqueAccession$accession,unique(miRNA_expr_long$barcode)[i]] <- patient_uniqueAccession$read_count
}

rownames(miRNA_expr_counts) <- miRNA_expr_counts$miRNA_ID_mature
miRNA_expr_counts <- miRNA_expr_counts[-1]
miRNA_expr_counts[1:10,1:10]

any(apply(miRNA_expr_counts, 1, function(x) any(is.na(x)))) # check nas values

write.table(miRNA_expr_counts, paste0(".", data_folder, "/miRNA_expr_counts_mature.txt"), sep = "\t", row.names = T, quote = F)


# How to collapse miRNA isoforms to one mature miRNA:
# Option 1: keep accession of max counts (also collapses _1, _2 (functionally identical miRNAs expressed from different loci), etc)
# Option 2: official version (https://github.com/bcgsc/mirna): Sum of all counts by MIMAT accession


seq_samples_clin <- data.frame(sample_code = colnames(miRNA_expr_counts),
                               project = vapply(strsplit(colnames(miRNA_expr_counts), '\\-'), `[`, 1, FUN.VALUE=character(1)),
                               tumor_code = vapply(strsplit(colnames(miRNA_expr_counts), '\\-'), `[`, 2, FUN.VALUE=character(1)),
                               case_identifier = vapply(strsplit(colnames(miRNA_expr_counts), '\\-'), `[`, 3, FUN.VALUE=character(1)),
                               tissue_code = vapply(strsplit(colnames(miRNA_expr_counts), '\\-'), `[`, 4, FUN.VALUE=character(1)),
                               nucleic_acid = vapply(strsplit(colnames(miRNA_expr_counts), '\\-'), `[`, 5, FUN.VALUE=character(1)),
                               TARGET_USI = vapply(strsplit(colnames(miRNA_expr_counts), '\\-'), function(x)
                                 paste(x[seq.int(3)], collapse='_'), character(1L)))





# Check TARGET codes at https://ocg.cancer.gov/sites/default/files/OCG-Project-Codes_Tissues-and-Samples_03-02-2021.pdf
# Tenemos muestras de todo tipo: AML, induction failure AML, de Bone marrow, de Peripheral Blood, etc.
# Vamos a hacer un data frame con estos datos:


table(seq_samples_clin$project)
table(seq_samples_clin$tumor_code)
table(seq_samples_clin$tissue_code)
table(seq_samples_clin$nucleic_acid)


# Vamos a filtrar muestras:
# Nos quedamos con las muestras de de novo AML (tumor_code = 20) y tejido no canceroso (tumor_code = 00):
# Nos quitamos las Induction failure AML (tumor_code = 21)
seq_samples_clin <- seq_samples_clin[seq_samples_clin$tumor_code == "20" | seq_samples_clin$tumor_code == "00",]
table(seq_samples_clin$tissue_code)

# Sabemos que los tissue_code que tenemos son:
# 03 = Primary blood derived cancer – peripheral blood
# 04 = Recurrent blood derived cancer – bone marrow
# 09 = Primary blood derived cancer – bone marrow
# 14 = Bone marrow normal
# 40 = Recurrent blood derived cancer – peripheral blood
# 50 = Cell line from patient tumor

# Nos quedamos con las muestras de Bone Marrow (04A, 09A) y Bone Marrow Normal (14A):
seq_samples_clin <- seq_samples_clin[seq_samples_clin$tissue_code == "04A" | seq_samples_clin$tissue_code == "09A" | seq_samples_clin$tissue_code == "14A",]

table(seq_samples_clin$tumor_code)
table(seq_samples_clin$tissue_code)

# Tenemos 1509 muestras (42 Recurrent Blood Derived Cancer y 1467 Primary Blood Derived Cancer) y 62 controles.


# Vamos a ver si hay pacientes duplicados aqui (case_identifier):
table(duplicated(seq_samples_clin$case_identifier)) # tenemos 34 duplicados

# Vamos a ver los duplicados y sus originales:

alldups <- seq_samples_clin[duplicated(seq_samples_clin$case_identifier) | duplicated(seq_samples_clin$case_identifier, fromLast = T),]
alldups <- alldups[order(alldups$case_identifier),]
#utils::View(alldups)
# Vemos que hay pacientes repetidos porque:
# Tienen muestras de Recurrent blood derived cancer (04) y Primary blood derived cancer (09).
# O bien tienen repetidas las amplificaciones de la misma alicuota (nucleic_acid: 01R, 02R, 03R...)

# Vamos a quedarnos con las muestras de Primary blood derived cancer (09) y vamos a quitarnos las Recurrent (04), que solo hay 42 muestras

seq_samples_clin <- seq_samples_clin[!seq_samples_clin$tissue_code == "04A",]

table(seq_samples_clin$tumor_code)
table(seq_samples_clin$tissue_code)
table(seq_samples_clin$nucleic_acid)

# Vamos a volver a ver duplicados:
table(duplicated(seq_samples_clin$case_identifier)) # Ya hay solo 4 duplicados
alldups <- seq_samples_clin[duplicated(seq_samples_clin$case_identifier) | duplicated(seq_samples_clin$case_identifier, fromLast = T),]
alldups <- alldups[order(alldups$case_identifier),]
#utils::View(alldups)

# Vemos que los duplicados son alicuotas repetidas.
# Nos quedamos con la ultima (nucleic_acid mas alto):
alldups

c <- unique(alldups$case_identifier)
dup_to_remove <- list()
for (i in 1:length(c)) {
  dup_to_remove[i] <- alldups[alldups$case_identifier==c[i],][order(alldups[alldups$case_identifier==c[i],]$nucleic_acid, decreasing = T),][2,]$sample_code
}

dup_to_remove
dim(seq_samples_clin)

seq_samples_clin <- seq_samples_clin[!seq_samples_clin$sample_code%in%dup_to_remove,]
dim(seq_samples_clin)

any(duplicated(seq_samples_clin$case_identifier)) # ya no hay pacientes duplicados


write.table(seq_samples_clin,
            file = paste0(".", data_folder, "/clinical_dataset_2.txt"),
            sep = ";")


#----


#### Compare counts mature and immature ####

mature <- read.table(paste0(".", data_folder, "/miRNA_expr_counts_mature.txt"), sep = "\t")
immature <- read.table(paste0(".", data_folder, "/miRNA_expr_counts_immature.txt"), sep = "\t")
iso <- read.delim(paste0(".", data_folder, "/TARGET_miRNA_expr_data_iso.txt"))

intersect(colnames(mature), colnames(immature))[1]
"TARGET.20.PATELT.03A.05R"

rownames(mature)
mir <- "-7a"
mir <- "-7b"
mir <- "-4470"
mir <- "-629"
mature_1 <- mature[grep(mir, rownames(mature)),"TARGET.20.PATELT.03A.05R", drop = F]
immature_1 <- immature[grep(mir, rownames(immature)),"TARGET.20.PATELT.03A.05R", drop = F]
iso_1 <- iso[iso$barcode=="TARGET-20-PATELT-03A-05R"&grepl(mir, iso$miRNA_ID),c("miRNA_ID", "read_count", "miRNA_region")]

sum(immature_1$TARGET.20.PATELT.03A.05R)
sum(mature_1$TARGET.20.PATELT.03A.05R)
sum(iso_1$read_count)

immature_1$TARGET.20.PATELT.03A.05R==sum(mature_1$TARGET.20.PATELT.03A.05R) 
immature_1$TARGET.20.PATELT.03A.05R==sum(iso_1$read_count) # immature is the sum of all isoforms (matures and precursor sequence)

# Ususally precursor, unnanotated, etc isoform are of very low counts. 
# --> It is ok to get the immature miRNA counts from summing the two mature miRNA counts.


#----



#### Juntar clinica y seleccionar pacientes ####

# Ahora vamos a juntar los dos data sets clinicos (los de la clinica de pacientes y los de las muestras secuenciadas):
# El factor comun de ambos datasets es el TARGET_USI:



all_samples_clin_data <- read.delim(paste0(".", data_folder, "/", "clinical_dataset.txt"),
                                    sep = ";")
seq_samples_clin <- read.delim(paste0(".", data_folder, "/", "clinical_dataset_2.txt"),
                               sep = ";")
counts <- read.table(paste0(".", data_folder, "/miRNA_expr_counts_mature.txt"), sep = "\t", check.names = F)


dim(seq_samples_clin)
dim(all_samples_clin_data)
table(seq_samples_clin$TARGET_USI%in%all_samples_clin_data$TARGET_USI) # 64 pacientes no tienen datos clínicos (posiblemente aqui esten los controles sanos)

not_in_clin <- seq_samples_clin[!seq_samples_clin$TARGET_USI%in%all_samples_clin_data$TARGET_USI,]
table(not_in_clin$tumor_code) 
to_remove <- not_in_clin[not_in_clin$tumor_code!="00",]$TARGET_USI
# efectivamente aqui tenemos los 62 sanos y tambien 2 pacientes (estos los quitamos):
seq_samples_clin <- seq_samples_clin[!seq_samples_clin$TARGET_USI%in%to_remove,]


# Fusionamos todos los datos clinicos:
full_clin <- merge(seq_samples_clin, all_samples_clin_data, by = "TARGET_USI", all.x = T)
#utils::View(full_clin)
rownames(full_clin) <- full_clin$sample_code

dim(full_clin) #1523 pacientes despues de seleccionar muestras de interes y quitar duplicados

# Reducimos el dataset de counts para coincidir con la clinica filtrada:
counts <- counts[,rownames(full_clin)]
dim(full_clin)
dim(counts) # 1523 pacientes y 1881 miRNAs

all(colnames(counts)==rownames(full_clin))

colnames(full_clin) <- gsub("\\.", "_", colnames(full_clin))
colnames(full_clin) <- gsub(" ", "_", colnames(full_clin)) # Importante para posterior analisis

#----


#### Filter miRNAs ####


library(edgeR)
keep <- filterByExpr(counts, min.count = 10)
counts_filt <- counts[keep,]
dim(counts)
dim(counts_filt) # We go from 2280 to 270 miRNAs


#----


#### Save ####

setwd('~/Desktop/Review code/miRNA_meta_analysis')
data_dir <- "./data"
dir.create(data_dir)

write.table(counts, file = paste0(data_dir, "/", "miRNA_expression_data.txt"), sep = "\t")
write.table(counts_filt, file = paste0(data_dir, "/", "miRNA_expression_data_filt.txt"), sep = "\t")
write.table(full_clin, file = paste0(data_dir, "/", "clin.txt"), sep = "\t")

#----




#### ------------------Differential Expression Analysis:---------------------- ####

#### Load Data and parameters for analysis: ####

setwd('~/Desktop/Review code/miRNA_meta_analysis')
data_dir <- "./data"
save_dir <- "./results"
dir.create(save_dir)

### LOAD ###
# Clinical data:
clin <- read.table(file = paste0(data_dir, "/", "clin.txt"), sep = "\t") # 1523 pacientes: 62 Sanos y 1461 AML
# Expression data:
expr <- read.table(file = paste0(data_dir, "/", "miRNA_expression_data_filt.txt"), sep = "\t") # 215 miRNAs (filtrados)
colnames(expr) <- gsub("\\.", "-", colnames(expr))

all(rownames(clin)==colnames(expr))


### PARAMETERS ###
#utils::View(clin) # We can see all the clinical data set and choose our variables of interest
variables <- c("tumor_code"
               #"Event_Free_Survival_Time_in_Days",
               #"Overall_Survival_Time_in_Days"
) # Is important that no variable name has characters like "-", ".", ";", " ", "(" or ")". Only "_" is valid
tmm_normalization <- FALSE # if we want to normalize counts by TMM method. Recommended when using limma or edgeR method, not for voomlimma.
DE_method = "voomLimma" # options: "voomLimma" / "limma / "edgeR"
DE_adjustMethod <- "BH" # options: "bonferroni" / "BH"...
pval <- 0.05 # Cut p-value
logFC <- 2 # Cut logFC
outliers <- c("TARGET.20.PASPLH.09A.01R") # choose outlier samples to quit (based on OutlierPCAHplot results)
dev.off() # Clean graphics just in case

#----


#### Outliers: ####

### PLOT OUTLIERS ###
library(rrcov)
pcaH <- PcaHubert(t(expr))
which(pcaH@flag=='FALSE')
sort(pcaH$sd)
pdf(paste0(save_dir, "/", "OutlierPCAHplot.pdf"),width = 10,height = 10) #export as pdf
plot(pcaH, xlim = c(0,max(pcaH$sd)+5))
dev.off()

### REMOVE OUTLIERS ###
clin <- clin[!rownames(clin)%in%outliers,]
expr <- expr[,rownames(clin)]
all(rownames(clin)==colnames(expr))


#----


#### Create Design Matrix: ####

# We first create a dgList object with edgeR package:
library(edgeR)
dgList <- DGEList(counts=expr, genes=rownames(expr))

clin$tumor_code <- relevel(factor(clin$tumor_code), ref = "0")

designMat <- model.matrix(~0+. , data = clin[,variables, drop=F]) # No intercept. We will need to specify the contrast matrix later.
#designMat <- drop.coef(designMat) # We remove redundant variables to keep the matrix full rank
designMat

#----


#### MDS plot: ####

set.seed(12345)
vartoplot <- variables[1] # we can select here the variable to color samples
factornum <- as.numeric(as.factor(clin[,vartoplot]))
library(RColorBrewer)
colors <- brewer.pal(n = max(factornum), name = 'Accent') # automatic coloring
colorVector <- colors[factornum]
lcpm <- cpm(dgList, log=TRUE)

library(limma)
pdf(paste0(save_dir, "/", "MDSplot.pdf"), width = 10, height = 10) #export as pdf
par(mar=c(17, 4.1, 4.1, 17), xpd=TRUE)
plotMDS(lcpm, 
        col = colorVector, 
        #xaxt = "n", 
        #yaxt = "n",
        #top = 10000,
        pch = 16,
        main = vartoplot, # distances between gene expression profiles
        #dim = c(3,4)
)
legend("bottomright", pch=20, legend = unique(clin[,vartoplot]), col = unique(colorVector),inset = c(0,-0.6))
dev.off()

#----


#### Normalization: ####

if (tmm_normalization==TRUE) {
  dgList <- calcNormFactors(dgList, method="TMM")
  tmm_normalized_counts <- cpm(dgList)
}

dgList <- estimateDisp(dgList, design = designMat, robust = TRUE)
plotBCV(dgList)

#----


#### Differential Expression: ####

### CONTRAST MATRIX ###
colnames(designMat)
contr.matrix <- makeContrasts(diff= tumor_code20 - tumor_code0, levels=designMat) # Important step: Select the contrasts to make.
contr.matrix

### DIFFERENTIAL EXPRESSION ###
# Several differential expression methods:

if (DE_method == "voomLimma") {
  v <- voom(dgList, designMat, plot = T) # voom normalization method (logCPM)
  vfit <- lmFit(v, designMat)
  vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
  efit <- eBayes(vfit)
  
  summary(decideTests(efit, adjust.method = DE_adjustMethod, p.value = pval, lfc = logFC))
  deg <- topTable(efit, coef=1, adjust.method = DE_adjustMethod, number = nrow(efit))
  write.table(deg, paste0(save_dir, "/", "deg_table.txt"), sep = "\t", quote = F)
}

if (DE_method == "limma") {
  fit <- lmFit(tmm_normalized_counts, designMat)
  fit2 <- contrasts.fit(fit, contr.matrix)
  efit <- eBayes(fit2)
  
  summary(decideTests(efit, adjust.method = DE_adjustMethod, p.value = pval, lfc = logFC))
  deg <- topTable(efit, coef=1, adjust.method = DE_adjustMethod, number = nrow(efit))
  write.table(deg, paste0(save_dir, "/", "deg_table.txt"), sep = "\t", quote = F)
}

if (DE_method == "edgeR") {
  fit <- glmQLFit(dgList, designMat, robust = TRUE)
  plotQLDisp(fit)
  qlf <- glmQLFTest(fit, contrast = contr.matrix)
  
  summary(decideTests(qlf, adjust.method = DE_adjustMethod))
  deg_tT <- topTags(qlf, n = nrow(qlf), adjust.method = DE_adjustMethod)
  deg <- deg_tT$table
  colnames(deg)[6] <- "adj.P.Val"
  write.table(deg, paste0(save_dir, "/", "deg_table.txt"), sep = "\t", quote = F)
}


#----

#### Get DE Genes: ####

### LOAD ###
deg <- read.table(paste0(save_dir, "/", "deg_table.txt"), sep = "\t")

### CUT ###
# We apply the cutoff to get the DE gene lists:
deg_cut <- deg[deg$adj.P.Val<=pval,]
deg_cut <- deg_cut[abs(deg_cut$logFC)>=logFC,]
up <- deg_cut[deg_cut$logFC>=0,]
down <- deg_cut[deg_cut$logFC<0,]

nrow(deg_cut)
nrow(up)
nrow(down)

### SAVE ###
write.table(deg_cut$genes, paste0(save_dir, "/", "deg_list.txt"), col.names = F, row.names = F, quote = F) # export deg list
write.table(up$genes, paste0(save_dir, "/", "Up_deg_list.txt"), col.names = F, row.names = F, quote = F) # export up deg list
write.table(down$genes, paste0(save_dir, "/", "Down_deg_list.txt"), col.names = F, row.names = F, quote = F) # export down deg list

#----


#### Heatmap: ####

vartoplot <- variables[1]
df_simple <- as.data.frame(clin[,vartoplot, drop=F])
df_simple[,1] <- as.factor(df_simple[,1])
df_simple.ordered <- df_simple[order(df_simple[,1], decreasing = T), , drop = FALSE]

expression_matrix_DE <- as.data.frame(v$E) # Use v$E or tmm_normalized_counts (depending of the DE method used)
expression_matrix_DE <- expression_matrix_DE[rownames(deg_cut),rownames(df_simple.ordered)]
rownames(df_simple.ordered) == colnames(expression_matrix_DE)

library(pheatmap)
dev.off()
pdf(paste0(save_dir, "/", "Heatmap.pdf"),width = 10,height = 10) #export as pdf
pheatmap(expression_matrix_DE,
         cluster_rows=T, 
         show_rownames=TRUE, 
         show_colnames =F, 
         cluster_cols=F, 
         annotation_col=df_simple.ordered,
         scale = 'row',
         breaks = seq(-3,3,length.out = 100),
         main = "Differential Expression ",
         fontsize_row = 4,
         fontsize_col = 3,
         color = greenred(100)
)
dev.off()


#----











#### ----------------Survival Analysis:-------------------- ####



#### Choose variables ####

data_subset <- "All_Risks"
folds <- 10

#data_subset <- "High_Risk"
#folds <- 5

#data_subset <- "Standard_Risk"
#folds <- 10

#data_subset <- "Low_Risk"
#folds <- 10


CIcovars <- c(#"Risk_group"#,
              #"WBC_at_Diagnosis",
              #"SCT_in_1st_CR",
              #"FLT3.ITD_positive.",
              #"CNS_disease"
              )
#stopifnot(all(CIcovars%in%colnames(surv_table_extended)))
#colnames(surv_table_extended)
#colnames(clin)



#+ Age_at_Diagnosis_in_Days 
#+ WBC_at_Diagnosis 
#+ SCT_in_1st_CR        # post transplant variable. Severan NAs
#+ FAB_Category        # FAB outdated
#+ MLL                 # mutations not included: correlated to risk group
#+ FLT3.ITD_positive.
#+ NPM_mutation
#+ CEBPA_mutation
#+ Primary_Cytogenetic_Code
#+ Bone_marrow_leukemic_blast_percentage_...
#+ Peripheral_blasts_...
#+ CNS_disease
#+ Chloroma

#----

#### Load data and variables ####

setwd('~/Desktop/Review code/miRNA_meta_analysis')
data_dir <- "./data"
dir.create(data_dir, recursive = T, showWarnings = F)
#save_dir <- "./results"
save_dir <- "./results_peerReview"
dir.create(save_dir, recursive = T, showWarnings = F)

### LOAD 
# Clinical data:
clin <- read.table(file = paste0(data_dir, "/", "clin.txt"), sep = "\t") # 1523 pacientes: 62 Sanos y 1461 AML
#utils::View(clin)
# Expression data:
expr <- read.table(file = paste0(data_dir, "/", "miRNA_expression_data.txt"), sep = "\t") # 1881 miRNAs (no filtrados)


### MODIFY 
rownames(expr) <- gsub("-", "_", rownames(expr))
colnames(expr) <- gsub("\\.", "_", colnames(expr))
rownames(clin) <- gsub("\\.", "_", rownames(clin))
rownames(clin) <- gsub("\\-", "_", rownames(clin))
all(colnames(expr)==rownames(clin))

### FIX VARIABLES

table(clin$Risk_group, exclude = NULL)
clin$Risk_group[!(clin$Risk_group=="High Risk" | clin$Risk_group=="Standard Risk" | clin$Risk_group=="Low Risk")] <- NA
table(clin$Risk_group, exclude = NULL)


table(clin$FLT3_ITD_positive_, exclude = NULL)
clin$FLT3_ITD_positive_[grepl("No", clin$FLT3_ITD_positive_, ignore.case = T)] <- "No"
clin$FLT3_ITD_positive_[grepl("Yes", clin$FLT3_ITD_positive_, ignore.case = T)] <- "Yes"
clin$FLT3_ITD_positive_[grepl("Unknown", clin$FLT3_ITD_positive_, ignore.case = T)] <- NA
table(clin$FLT3_ITD_positive_, exclude = NULL)

table(clin$SCT_in_1st_CR, exclude = NULL)
clin$SCT_in_1st_CR[clin$SCT_in_1st_CR=="Unknown"] <- NA
table(clin$SCT_in_1st_CR, exclude = NULL)

clin$WBC_at_Diagnosis_status <- cut(clin$WBC_at_Diagnosis,
                                    breaks = c(-Inf, 100, Inf),
                                    labels = c("low", "high")) # Los valores de expresion igual al cut se asignan a low.
table(clin$WBC_at_Diagnosis_status, exclude = NULL)


table(clin$MRD_at_end_of_course_1, exclude = NULL)
table(clin$MRD_at_end_of_course_2, exclude = NULL)
summary(as.numeric(clin$MRD_at_end_of_course_1_perc))
summary(as.numeric(clin$MRD_at_end_of_course_2_perc))

colnames(clin)[colnames(clin)=="MRD___at_end_of_course_1"] <- "MRD_at_end_of_course_1_perc"
colnames(clin)[colnames(clin)=="MRD___at_end_of_course_2"] <- "MRD_at_end_of_course_2_perc"
clin$MRD_at_end_of_course_1[clin$MRD_at_end_of_course_1 %in% c("Unknown")] <- NA
clin$MRD_at_end_of_course_2[clin$MRD_at_end_of_course_2 %in% c("Unknown")] <- NA
clin$MRD_at_end_of_course_1_perc <- gsub(" \\|\\|\\|.+", "", clin$MRD_at_end_of_course_1_perc)
clin$MRD_at_end_of_course_2_perc <- gsub(" \\|\\|\\|.+", "", clin$MRD_at_end_of_course_2_perc)
clin$MRD_at_end_of_course_1_perc <- as.numeric(clin$MRD_at_end_of_course_1_perc)
clin$MRD_at_end_of_course_2_perc <- as.numeric(clin$MRD_at_end_of_course_2_perc)

table(clin$MRD_at_end_of_course_1, exclude = NULL)
table(clin$MRD_at_end_of_course_2, exclude = NULL)
summary(as.numeric(clin$MRD_at_end_of_course_1_perc))
summary(as.numeric(clin$MRD_at_end_of_course_2_perc))


table(clin$Primary_Cytogenetic_Code, exclude = NULL)
clin$Primary_Cytogenetic_Code[grepl("inv\\(16\\)", clin$Primary_Cytogenetic_Code)] <- "inv(16)"
clin$Primary_Cytogenetic_Code[grepl("MLL", clin$Primary_Cytogenetic_Code)] <- "MLL"
clin$Primary_Cytogenetic_Code[grepl("Normal", clin$Primary_Cytogenetic_Code)] <- "Normal"
clin$Primary_Cytogenetic_Code[grepl("Unknown", clin$Primary_Cytogenetic_Code)] <- NA

table(clin$NPM_mutation, exclude = NULL)
clin$NPM_mutation[grepl("no", clin$NPM_mutation, ignore.case = T) & !grepl("yes", clin$NPM_mutation, ignore.case = T)] <- "No"
clin$NPM_mutation[!clin$NPM_mutation %in% c("No", "Yes")] <- NA

table(clin$CEBPA_mutation, exclude = NULL)
clin$CEBPA_mutation[grepl("no", clin$CEBPA_mutation, ignore.case = T) & clin$CEBPA_mutation!="Unknown"] <- "No"
clin$CEBPA_mutation[!clin$CEBPA_mutation %in% c("No", "Yes")] <- NA

table(clin$WT1_mutation, exclude = NULL)
clin$WT1_mutation[grepl("no", clin$WT1_mutation, ignore.case = T) & !grepl("yes", clin$WT1_mutation, ignore.case = T) & clin$WT1_mutation!="Unknown"] <- "No"
clin$WT1_mutation[grepl("yes", clin$WT1_mutation, ignore.case = T) & !grepl("no", clin$WT1_mutation, ignore.case = T) & clin$WT1_mutation!="Unknown"] <- "Yes"
clin$WT1_mutation[!clin$WT1_mutation %in% c("No", "Yes")] <- NA

table(clin$c_Kit_Mutation_Exon_8, exclude = NULL)
clin$c_Kit_Mutation_Exon_8[grepl("no", clin$c_Kit_Mutation_Exon_8, ignore.case = T) & !grepl("yes", clin$c_Kit_Mutation_Exon_8, ignore.case = T) & !grepl("Not done", clin$c_Kit_Mutation_Exon_8, ignore.case = T)] <- "No"
clin$c_Kit_Mutation_Exon_8[clin$c_Kit_Mutation_Exon_8 == "Not done ||| No"] <- "No"
clin$c_Kit_Mutation_Exon_8[grepl("yes", clin$c_Kit_Mutation_Exon_8, ignore.case = T)] <- "Yes"
clin$c_Kit_Mutation_Exon_8[grepl("Not done", clin$c_Kit_Mutation_Exon_8, ignore.case = T)] <- NA

table(clin$c_Kit_Mutation_Exon_17, exclude = NULL)
clin$c_Kit_Mutation_Exon_17[grepl("no", clin$c_Kit_Mutation_Exon_17, ignore.case = T) & !grepl("yes", clin$c_Kit_Mutation_Exon_17, ignore.case = T) & !grepl("Not done", clin$c_Kit_Mutation_Exon_17, ignore.case = T)] <- "No"
clin$c_Kit_Mutation_Exon_17[clin$c_Kit_Mutation_Exon_17 == "Not done ||| No"] <- "No"
clin$c_Kit_Mutation_Exon_17[grepl("Not done", clin$c_Kit_Mutation_Exon_17, ignore.case = T)] <- NA

table(clin$FLT3_PM, exclude = NULL)
clin$FLT3_PM[grepl("no", clin$FLT3_PM, ignore.case = T) & !grepl("yes", clin$FLT3_PM, ignore.case = T) & clin$FLT3_PM!="Unknown"] <- "No"
clin$FLT3_PM[grepl("yes", clin$FLT3_PM, ignore.case = T) & !grepl("no", clin$FLT3_PM, ignore.case = T) & clin$FLT3_PM!="Unknown"] <- "Yes"
clin$FLT3_PM[!clin$FLT3_PM %in% c("No", "Yes")] <- NA



# cytogenetic abnormalities:
for (i in 25:42) {
  print(table(clin[,i], exclude = NULL))
}
for (i in 25:42) {
  clin[,i][clin[,i]=="Unknown"] <- NA
}


# Fix survival variables: choose the max survival of several values (last measures)
clin$Overall_Survival_Time_in_Days <- sapply(strsplit(as.character(clin$Overall_Survival_Time_in_Days)," \\|\\|\\| "), function(x) max(as.numeric(x)) )
clin$Event_Free_Survival_Time_in_Days <- sapply(strsplit(as.character(clin$Event_Free_Survival_Time_in_Days)," \\|\\|\\| "), function(x) max(as.numeric(x)) )



# Explore Low Risk group characteristics:
table(clin$Risk_group, exclude = NULL)
# Explore CEBPA and risk groups:
table(clin$CEBPA_mutation, clin$Risk_group, exclude = NULL)
# Explore NPM1 and risk groups:
table(clin$NPM_mutation, clin$Risk_group, exclude = NULL)
# Explore t(8;21) and risk groups:
table(clin$t_8_21, clin$Risk_group, exclude = NULL)
# Explore inv(16) and risk groups:
table(clin$inv_16, clin$Risk_group, exclude = NULL)

# Explore High Risk group characteristics:
# Explore FLT3 and risk groups:
table(clin$FLT3_ITD_positive_, clin$Risk_group, exclude = NULL)
# Explore FLT3-ITD allelic ratio:
tapply(as.numeric(clin$FLT3_ITD_allelic_ratio), clin$Risk_group, summary)
# Explore monosomy 5 and risk groups:
table(clin$monosomy_5, clin$Risk_group, exclude = NULL)
# Explore monosomy 7 and risk groups:
table(clin$monosomy_7, clin$Risk_group, exclude = NULL)
# Explore del(5q) and risk groups:
table(clin$del5q, clin$Risk_group, exclude = NULL)



all_miRNAs_mature <- rownames(expr) # save only mature miRNAs
all_miRNAs <- rownames(expr) # this will expand later


#----

#### Patient filter ####
# Seleccionamos solo los enfermos y que tengan datos de supervivencia:
clin <- clin[clin$tumor_code==20 &
               !is.na(clin$Overall_Survival_Time_in_Days) &
               !is.na(clin$Event_Free_Survival_Time_in_Days) &
               !is.na(clin$Vital_Status),] # from 1523 to 1441 patients (-82)

any(is.na(clin$Overall_Survival_Time_in_Days)) # No hay NAs
any(is.na(clin$Vital_Status)) # No hay NAs]

nrow(clin)

# FILTRO EDAD
max(clin$Age_at_Diagnosis_in_Days/365)
#clin <- clin[clin$Age_at_Diagnosis_in_Days/365<19,]
clin <- clin[clin$Age_at_Diagnosis_in_Days/365<21,] # <21 years

nrow(clin)

# Filter covariates NAs
if (length(CIcovars)>=1) {
  for (i in 1:length(CIcovars)) {
    nas <- nrow(clin[is.na(clin[,CIcovars[i]]),])
    print(paste0(CIcovars[i], " has ", nas, " NAs")) # print how many NAs in each covariate
  }
  n_original <- nrow(clin)
  for (i in 1:length(CIcovars)) {
    clin <- clin[!is.na(clin[,CIcovars[i]]),] # remove all NAs in any covariate
  }
  print(paste0(n_original-nrow(clin), " patients removed due to NAs in covariates"))
}


# Apply filter to expr dataset:
expr <- expr[,rownames(clin)]

stopifnot(all(colnames(expr)==rownames(clin)))


#----

#### Patient selection by risk ####

### Get variables info:

table(clin$Protocol, exclude = NULL) # Diferentes protocolos
table(clin$Gemtuzumab_ozogamicin_treatment, exclude = NULL)
table(clin$Protocol[!is.na(clin$Gemtuzumab_ozogamicin_treatment)])
#utils::View(clin)
table(clin$Protocol[clin$FLT3.ITD_positive.=="No"])
table(clin$Risk_group, exclude = NULL)
table(clin$Risk_group, clin$Protocol,  exclude = NULL)

clinfiles_path <- "./GDCdata/TARGET-AML/harmonized/Clinical/Clinical_Supplement"
clinfiles <- list.files(path = clinfiles_path, recursive = T)
library(readxl)
TARGET_AML_ClinicalData_CDE <- read_excel(paste0(clinfiles_path, "/", clinfiles[4])) # data set de diccionario de variables
#utils::View(TARGET_AML_ClinicalData_CDE)



### Seleccionar subset de datos elegido:
if (data_subset=="High_Risk") {
  clin <- clin[clin$Risk_group=="High Risk"&!is.na(clin$Risk_group),]
  expr <- expr[, rownames(clin)]
} 
if (data_subset=="Standard_Risk") {
  clin <- clin[clin$Risk_group=="Standard Risk"&!is.na(clin$Risk_group),]
  expr <- expr[, rownames(clin)]
} 
if (data_subset=="Low_Risk") {
  clin <- clin[clin$Risk_group=="Low Risk"&!is.na(clin$Risk_group),]
  expr <- expr[, rownames(clin)]
} 

#----

#### Discovery - Validation set selection ####

nrow(clin)
table(clin$Protocol, exclude = NULL)
Discovery_set <- rownames(clin)[clin$Protocol=="AAML0531"|clin$Protocol=="AAML03P1"|clin$Protocol=="CCG2961"]
Validation_set <- rownames(clin)[clin$Protocol=="AAML1031"]
#Discovery_set <- rownames(clin)[clin$Protocol=="AAML1031"|clin$Protocol=="AAML03P1"]
#Validation_set <- rownames(clin)[clin$Protocol=="AAML0531"|clin$Protocol=="CCG2961"]
length(Discovery_set)
length(Validation_set)
clin$PatientSet <- "Discovery"
clin[Validation_set,]$PatientSet <- "Validation"

#table(clin[Validation_set,]$Risk_group)
#table(clin[Discovery_set,]$Risk_group)


#----

#### Individual MiRNA selection ####

search <- "hsa_mir_146a"
all_miRNAs[grep(search, all_miRNAs, ignore.case = T)]
sel_miRNAs <- c("hsa_mir_199a", 
                "hsa_mir_381",
                "hsa_miR_206", # already mature
                "hsa_miR_193b_3p", # already mature
                "hsa_miR_139_5p", # already mature
                "hsa_mir_100",
                "hsa_mir_195",
                "hsa_miR_375", # already mature
                "hsa_mir_335",
                "hsa_mir_370",
                "hsa_mir_29a",
                "hsa_mir_122",
                "hsa_mir_192",
                "hsa_mir_155",
                "hsa_mir_196b",
                "hsa_mir_34b")
miR3_miRNAs <- c("hsa_mir_146b",
                 "hsa_mir_181c",
                 "hsa_mir_4786")
miR4_miRNAs <- c("hsa_mir_509",
                 "hsa_mir_542",
                 "hsa_mir_3667",
                 "hsa_mir_146a")


# Calculate immature miRNAs counts (by sum of matures counts (unnormalized))
# Create new miRNAs:


newmir_df <- data.frame(row.names = colnames(expr),
                        "hsa_mir_199a" = colSums(expr[c("hsa_miR_199a_3p", "hsa_miR_199a_5p"),]),
                        "hsa_mir_381" = colSums(expr[c("hsa_miR_381_3p", "hsa_miR_381_5p"),]),
                        "hsa_mir_100" = colSums(expr[c("hsa_miR_100_3p", "hsa_miR_100_5p"),]),
                        "hsa_mir_195" = colSums(expr[c("hsa_miR_195_3p", "hsa_miR_195_5p"),]),
                        "hsa_mir_335" = colSums(expr[c("hsa_miR_335_3p", "hsa_miR_335_5p"),]),
                        "hsa_mir_370" = colSums(expr[c("hsa_miR_370_3p", "hsa_miR_370_5p"),]),
                        "hsa_mir_29a" = colSums(expr[c("hsa_miR_29a_3p", "hsa_miR_29a_5p"),]),
                        "hsa_mir_122" = colSums(expr[c("hsa_miR_122_3p", "hsa_miR_122_5p"),]),
                        "hsa_mir_192" = colSums(expr[c("hsa_miR_192_3p", "hsa_miR_192_5p"),]),
                        "hsa_mir_155" = colSums(expr[c("hsa_miR_155_3p", "hsa_miR_155_5p"),]),
                        "hsa_mir_196b" = colSums(expr[c("hsa_miR_196b_3p", "hsa_miR_196b_5p"),]),
                        "hsa_mir_34b" = colSums(expr[c("hsa_miR_34b_3p", "hsa_miR_34b_5p"),]),
                        "hsa_mir_146b" = colSums(expr[c("hsa_miR_146b_3p", "hsa_miR_146b_5p"),]), # from miR3
                        "hsa_mir_181c" = colSums(expr[c("hsa_miR_181c_3p", "hsa_miR_181c_5p"),]), # from miR3
                        "hsa_mir_4786" = colSums(expr[c("hsa_miR_4786_3p", "hsa_miR_4786_5p"),]), # from miR3
                        "hsa_mir_509" = colSums(expr[c("hsa_miR_509_3p", "hsa_miR_509_5p"),]), # from miR4
                        "hsa_mir_542" = colSums(expr[c("hsa_miR_542_3p", "hsa_miR_542_5p"),]), # from miR4
                        "hsa_mir_3667" = colSums(expr[c("hsa_miR_3667_3p", "hsa_miR_3667_5p"),]), # from miR4
                        "hsa_mir_146a" = colSums(expr[c("hsa_miR_146a_3p", "hsa_miR_146a_5p"),])) # from miR4

stopifnot(all(colnames(expr)==colnames(t(newmir_df)))) # check

expr_notnorm_mature <- expr # Save expr of matures before adding immatures

expr <- rbind(expr, t(newmir_df)) # add new immature mirnas to expr

all_miRNAs <- rownames(expr) # overwrite


#----

#### NORMALIZE EXPRESSION ####

# Normalization to log2(CPM), also called log2(RPM)

head(expr[1:5,1:5])

expr_notnorm <- expr # save unnormalized counts

library(edgeR)
d <- DGEList(counts=expr, genes=rownames(expr))
d <- calcNormFactors(d
                     , method="TMM" # If we want TMM normalization
                     )
expr_tmm <- cpm(d, normalized.lib.sizes = T, log = T)
#expr_tmm <- cpm(d, normalized.lib.sizes = T, log = F)

head(expr_tmm[1:5,1:5])

expr <- expr_tmm

#----

#### Crear data frame con variables de interes ####

# Association between single miRNAs and survival (Univariate Cox Proportional Hazards Model)
# http://www.sthda.com/english/wiki/cox-proportional-hazards-model

library("survival")
library("survminer")

texpr <- as.data.frame(t(expr)) # t of normalized expression
sort(colnames(texpr))

surv_vars <- data.frame(tumor_code = clin$tumor_code,
                        OS = clin$Overall_Survival_Time_in_Days,
                        EFS = clin$Event_Free_Survival_Time_in_Days,
                        Vital_Status = clin$Vital_Status,
                        row.names = rownames(clin))
surv_vars$status <- 0
surv_vars$status[surv_vars$Vital_Status=="Dead"] <- 1



surv_table <- merge(surv_vars,
      #texpr[,miRNAs], #si quisieramos hacerlo solo con los miRNAs seleccionados
      texpr,
      by = "row.names")

#utils::View(surv_table)

rownames(surv_table) <- surv_table$Row.names

#----

#### Create surv_dataset with all interesting variables: ####


### Select co-variables:
clin$Gender
clin$Race
clin$Ethnicity
clin$Age_at_Diagnosis_in_Days
clin$Risk_group
clin$FAB_Category
clin$WBC_at_Diagnosis
clin$Bone_marrow_leukemic_blast_percentage____
clin$Peripheral_blasts____
clin$CNS_disease
clin$Chloroma
clin$FLT3_ITD_positive_
clin$NPM_mutation
clin$CEBPA_mutation
clin$WT1_mutation
clin$c_Kit_Mutation_Exon_8
clin$c_Kit_Mutation_Exon_17
clin$t_8_21_
clin$inv_16_
clin$SCT_in_1st_CR
clin$PatientSet
clin$FLT3_PM
colnames(clin[,25:42])


# Creamos data frame extendido:

surv_vars_extended <- data.frame(row.names = rownames(clin),
                                 tumor_code = clin$tumor_code,
                                 OS = clin$Overall_Survival_Time_in_Days,
                                 EFS = clin$Event_Free_Survival_Time_in_Days,
                                 Vital_Status = clin$Vital_Status,
                                 Gender = clin$Gender,
                                 Race = clin$Race,
                                 Ethnicity = clin$Ethnicity,
                                 Age_at_Diagnosis_in_Days = clin$Age_at_Diagnosis_in_Days,
                                 Risk_group = clin$Risk_group,
                                 FAB_Category = clin$FAB_Category,
                                 WBC_at_Diagnosis = clin$WBC_at_Diagnosis,
                                 WBC_at_Diagnosis_status = clin$WBC_at_Diagnosis_status,
                                 Bone_marrow_leukemic_blast_percentage____ = clin$Bone_marrow_leukemic_blast_percentage____,
                                 Peripheral_blasts____ = clin$Peripheral_blasts____,
                                 CNS_disease = clin$CNS_disease,
                                 Chloroma = clin$Chloroma,
                                 FLT3.ITD_positive. = clin$FLT3_ITD_positive_,
                                 NPM_mutation = clin$NPM_mutation,
                                 CEBPA_mutation = clin$CEBPA_mutation,
                                 WT1_mutation = clin$WT1_mutation,
                                 FLT3_PM = clin$FLT3_PM,
                                 c_Kit_Mutation_Exon_8 = clin$c_Kit_Mutation_Exon_8,
                                 c_Kit_Mutation_Exon_17 = clin$c_Kit_Mutation_Exon_17,
                                 SCT_in_1st_CR = clin$SCT_in_1st_CR,
                                 Primary_Cytogenetic_Code = clin$Primary_Cytogenetic_Code,
                                 Protocol = clin$Protocol,
                                 PatientSet = clin$PatientSet,
                                 MRD_at_end_of_course_1 = clin$MRD_at_end_of_course_1,
                                 MRD_at_end_of_course_2 = clin$MRD_at_end_of_course_2,
                                 MRD_at_end_of_course_1_perc = clin$MRD_at_end_of_course_1_perc,
                                 MRD_at_end_of_course_2_perc = clin$MRD_at_end_of_course_2_perc)
surv_vars_extended <- merge(surv_vars_extended,
                            clin[,25:42],
                            by = "row.names")
rownames(surv_vars_extended) <- surv_vars_extended$Row.names
dim(surv_vars_extended)

# Modify
table(surv_vars_extended$Primary_Cytogenetic_Code, exclude = NULL)
surv_vars_extended$Primary_Cytogenetic_Code[grepl("\\|\\|\\|", surv_vars_extended$Primary_Cytogenetic_Code)] <- "Other"
table(surv_vars_extended$FAB_Category, exclude = NULL)
surv_vars_extended$FAB_Category[grepl("\\|\\|\\|", surv_vars_extended$FAB_Category)] <- "Unknown"
table(surv_vars_extended$SCT_in_1st_CR, exclude = NULL)
surv_vars_extended$SCT_in_1st_CR[is.na(surv_vars_extended$SCT_in_1st_CR)] <- "Unknown"


# Create status variable
surv_vars_extended$status <- 0
surv_vars_extended$status[surv_vars_extended$Vital_Status=="Dead"] <- 1


miRNAs <- all_miRNAs

surv_table_extended <- merge(surv_vars_extended,
                    #texpr[,miRNAs], #si quisieramos hacerlo solo con los miRNAs seleccionados
                    surv_table[,miRNAs],
                    by = "row.names")
rownames(surv_table_extended) <- surv_table_extended$Row.names
#utils::View(surv_table_extended)

#----

#### AMLmiR36 Scores ####

search <- "30c"
all_miRNAs[grep(search, all_miRNAs)]

AMLmiR36_miRNAs <- c("hsa_miR_106a_3p",
                     "hsa_miR_584_5p",
                     "hsa_miR_1247_3p",
                     "hsa_miR_155_5p",
                     "hsa_miR_130b_3p",
                     "hsa_miR_320a",
                     "hsa_miR_34c_5p",
                     "hsa_miR_30c_2_3p",
                     "hsa_miR_30e_3p",
                     "hsa_miR_450a_5p",
                     "hsa_miR_296_5p",
                     "hsa_miR_935",
                     "hsa_miR_502_3p",
                     "hsa_miR_181b_3p",
                     "hsa_miR_363_3p",
                     "hsa_miR_362_3p",
                     "hsa_miR_100_5p",
                     "hsa_miR_132_3p",
                     "hsa_miR_340_3p",
                     "hsa_miR_181c_5p",
                     "hsa_miR_148b_3p",
                     "hsa_miR_4662a_5p",
                     "hsa_miR_1287_3p",
                     "hsa_miR_664b_5p",
                     "hsa_miR_539_5p",
                     "hsa_miR_217",
                     "hsa_miR_181c_3p",
                     "hsa_miR_335_3p",
                     "hsa_miR_1180_3p",
                     "hsa_miR_202_5p",
                     "hsa_miR_146a_5p",
                     "hsa_miR_2110",
                     "hsa_miR_375",
                     "hsa_let_7g_5p",
                     "hsa_miR_139_5p",
                     "hsa_miR_409_5p")

AMLmiR36_coef <- c(0.244530267,
                   0.155286578,
                   0.070234817,
                   0.06433791,
                   0.042520575,
                   0.03910599,
                   0.027657901,
                   0.018074572,
                   0.017215965,
                   0.015253859,
                   0.015230088,
                   0.013849354,
                   0.013641305,
                   0.010667128,
                   0.006792342,
                   0.005215312,
                   -0.001863252,
                   -0.002882947,
                   -0.003799344,
                   -0.011324553,
                   -0.018062943,
                   -0.019307849,
                   -0.02419678,
                   -0.032959784,
                   -0.037325674,
                   -0.040306357,
                   -0.040989372,
                   -0.04528897,
                   -0.047275703,
                   -0.048514515,
                   -0.050372946,
                   -0.059026016,
                   -0.062192126,
                   -0.062530666,
                   -0.063621683,
                   -0.071435852)

stopifnot(all(AMLmiR36_miRNAs %in% all_miRNAs)) # check


### Getting AMLmiR36 Scores:

AMLmiR36 <- data.frame(miRNA = AMLmiR36_miRNAs,
                       Coefficient = AMLmiR36_coef)

#write.table(AMLmiR36, paste0(save_dir, "/signature_AMLmiR36.txt"), sep = "\t", row.names = F)

#AMLmiR36 <- read.delim(paste0(save_dir, "/signature_AMLmiR36.txt"))
#AMLmiR36_miRNAs <- AMLmiR36$miRNA
#AMLmiR36_coef <- AMLmiR36$Coefficient

# Getting Score for each patient: sum of all coef * miRNA:
AMLmiR36_scores <- as.data.frame(colSums(expr[AMLmiR36$miRNA,]*AMLmiR36$Coefficient)) # using norm expresion (log(cpm))
colnames(AMLmiR36_scores) <- "AMLmiR36_score"


### Adding AMLmiR36 score to surv_table:
stopifnot(all(rownames(surv_table)==rownames(AMLmiR36_scores))) # check
stopifnot(all(rownames(surv_table_extended)==rownames(AMLmiR36_scores))) # check
surv_table$AMLmiR36_score <- AMLmiR36_scores$AMLmiR36_score
surv_table_extended$AMLmiR36_score <- AMLmiR36_scores$AMLmiR36_score

#----

#### miR3 scores ####


search <- "146b"
all_miRNAs[grep(search, all_miRNAs)]
miR3_miRNAs <- c("hsa_mir_146b",
                 "hsa_mir_181c",
                 "hsa_mir_4786")
miR3_coef <- c(1.652,
               -1.838,
               -1.455)

stopifnot(all(miR3_miRNAs%in%all_miRNAs)) # check

### Getting miR3 Scores:

miR3 <- data.frame(miRNA = miR3_miRNAs,
                   Coefficient = miR3_coef)

#write.table(miR3, paste0(save_dir, "/signature_miR3.txt"), sep = "\t", row.names = F)

#miR3 <- read.delim(paste0(save_dir, "/signature_miR3.txt"))
#miR3_miRNAs <- miR3$miRNA
#miR3_coef <- miR3$Coefficient

# Getting Score for each patient: sum of all coef * miRNA:
miR3_scores <- as.data.frame(colSums(expr[miR3$miRNA,]*miR3$Coefficient)) # using norm expresion (log(cpm))
colnames(miR3_scores) <- "miR3_score"


### Adding miR3 score to surv_table:
stopifnot(all(rownames(surv_table)==rownames(miR3_scores))) # check
stopifnot(all(rownames(surv_table_extended)==rownames(miR3_scores))) # check
surv_table$miR3_score <- miR3_scores$miR3_score
surv_table_extended$miR3_score <- miR3_scores$miR3_score

#----

#### miR4 scores ####

search <- "146a"
all_miRNAs[grep(search, all_miRNAs)]
miR4_miRNAs <- c("hsa_mir_509",
                 "hsa_mir_542",
                 "hsa_mir_3667",
                 "hsa_mir_146a")
miR4_coef <- c(0.914,
               0.759,
               -0.837,
               -0.856)


stopifnot(all(miR3_miRNAs%in%all_miRNAs)) # check

### Getting miR4 Scores:

miR4 <- data.frame(miRNA = miR4_miRNAs,
                   Coefficient = miR4_coef)

#write.table(miR4, paste0(save_dir, "/signature_miR4.txt"), sep = "\t", row.names = F)

#miR4 <- read.delim(paste0(save_dir, "/signature_miR4.txt"))
#miR4_miRNAs <- miR4$miRNA
#miR4_coef <- miR4$Coefficient

# Getting Score for each patient: sum of all coef * miRNA:
miR4_scores <- as.data.frame(colSums(expr[miR4$miRNA,]*miR4$Coefficient)) # using norm expresion (log(cpm))
colnames(miR4_scores) <- "miR4_score"


### Adding miR4 score to surv_table:
stopifnot(all(rownames(surv_table)==rownames(miR4_scores))) # check
stopifnot(all(rownames(surv_table_extended)==rownames(miR4_scores))) # check
surv_table$miR4_score <- miR4_scores$miR4_score
surv_table_extended$miR4_score <- miR4_scores$miR4_score

#----

#### miR24 scores ####

#search <- "hsa_miR_20b_5p"
#all_miRNAs[grep(search, all_miRNAs)]
#miR24_miRNAs <- c("hsa_miR_20b_5p",
#                  "hsa_miR_223_3p",
#                  "hsa_miR_193a_3p",
#                  "hsa_miR_24_3p",
#                  "hsa_miR_128_3p",
#                  "hsa_miR_17_5p",
#                  "hsa_miR_199b_5p",
#                  "hsa_miR_181c_5p",
#                  "hsa_miR_181a_5p",
#                  "hsa_miR_181b_5p",
#                  "hsa_miR_21_5p",
#                  "hsa_miR_222_5p",
#                  "hsa_miR_331_5p",
#                  "hsa_miR_373_3p",
#                  "hsa_miR_708_5p",
#                  "hsa_miR_34b_5p",
#                  "hsa_miR_195_5p",
#                  "hsa_miR_151a_5p",
#                  "hsa_miR_30b_5p",
#                  "hsa_miR_22_3p",
#                  "hsa_let_7g_5p",
#                  "hsa_let_7i_5p",
#                  "hsa_miR_1290",
#                  "hsa_miR_9_5p")
#
#stopifnot(all(miR24_miRNAs%in%all_miRNAs)) # check
#
#### Getting mir24 coefs (using OS) IN DISCOVERY DATA: 
#surv <- Surv(time = surv_table_extended[Discovery_set,]$OS, event = surv_table_extended[Discovery_set,]$status) # OS
#cox_reg_model_uni_mir24 <- coxph(formula = surv ~ ., data = surv_table_extended[Discovery_set,][,miR24_miRNAs]) # univariate Cox to obtain coefs
#miR24_coef <- cox_reg_model_uni_mir24$coefficients
#stopifnot(all(names(miR24_coef) == miR24_miRNAs)) # check
#miR24_coef <- unname(miR24_coef)
#
#summary(cox_reg_model_uni_mir24) #we can check model of 24 mirnas
#
#miR24 <- data.frame(miRNA = miR24_miRNAs,
#                    Coefficient = miR24_coef)
#
#write.table(miR24, paste0(save_dir, "/signature_miR24.txt"), sep = "\t", row.names = F)

miR24 <- read.delim(paste0(save_dir, "/signature_miR24.txt"))
miR24_miRNAs <- miR24$miRNA
miR24_coef <- miR24$Coefficient



# Getting Score for each patient: sum of all coef * miRNA:
miR24_scores <- as.data.frame(colSums(expr[miR24$miRNA,]*miR24$Coefficient))
colnames(miR24_scores) <- "miR24_score"


### Adding miR24 score to surv_table:
stopifnot(all(rownames(surv_table)==rownames(miR24_scores))) # check
stopifnot(all(rownames(surv_table_extended)==rownames(miR24_scores))) # check
surv_table$miR24_score <- miR24_scores$miR24_score
surv_table_extended$miR24_score <- miR24_scores$miR24_score


#----

#### Develop a new signature ####

### Lasso method
## covariates (lasso only allows numeric cvrts)
#cvrts <- c(#"Risk_group"
#  "Age_at_Diagnosis_in_Days"
#  #, "WBC_at_Diagnosis" # Has 1 NA
#  #, "SCT_in_1st_CR"
#  , "FLT3.ITD_positive."
#  #, "CNS_disease"
#)
#
#
#stopifnot(all(all_miRNAs_mature%in%colnames(surv_table_extended[Discovery_set,])))
#
#library(edgeR)
#keep <- filterByExpr(expr_notnorm_mature, min.count = 10)
#expr_notnorm_mature_filt <- expr_notnorm_mature[keep,]
#dim(expr_notnorm_mature)
#dim(expr_notnorm_mature_filt) # We go from 2280 to 270 miRNAs
#
#data <- surv_table_extended[Discovery_set, c("OS", "status", rownames(expr_notnorm_mature_filt))]
#biomarkers <- rownames(expr_notnorm_mature_filt)
##data <- surv_table_extended[Discovery_set,]
##biomarkers <- rownames(expr_notnorm_mature)
#
#BiocManager::install("biospear", force = T)
#
#library(biospear)
#fit_lasso <- BMsel(data = data, # use Discovery set and filtered miRNAs
#                   x = biomarkers,
#                   y = c("OS", "status"),
#                   #z = cvrts, # no covariates
#                   inter = FALSE,
#                   method = "lasso",
#                   #dfmax = 50,
#                   folds = 5
#)
#
#length(fit_lasso$lasso)
#if (length(fit_lasso$lasso)<=0) {stop("No predictors from lasso model")}

#----

#### new score ####


#miR37_miRNAs <- names(fit_lasso$lasso)
#miR37_coef <- unname(fit_lasso$lasso)
#
#
#### Getting miR37 Scores:
#
#miR37 <- data.frame(miRNA = miR37_miRNAs,
#                   Coefficient = miR37_coef)
#
#
#write.table(miR37, paste0(save_dir, "/signature_new_miRNA.txt"), sep = "\t", row.names = F)


miR37 <- read.delim(paste0(save_dir, "/signature_new_miRNA.txt"))
miR37_miRNAs <- miR37$miRNA
miR37_coef <- miR37$Coefficient

# Getting Score for each patient: sum of all coef * miRNA:
miR37_scores <- as.data.frame(colSums(expr[miR37$miRNA,]*miR37$Coefficient))
colnames(miR37_scores) <- "miR37_score"


### Adding miR37 score to surv_table:
stopifnot(all(rownames(surv_table_extended)==rownames(miR37_scores)))
stopifnot(all(rownames(surv_table)==rownames(miR37_scores)))
surv_table$miR37_score <- miR37_scores$miR37_score
surv_table_extended$miR37_score <- miR37_scores$miR37_score


#----

#### Random 37 miRNA signature ####

#set.seed(1029)
#random_signature_miRNAs <- sample(all_miRNAs, 37)
#
#### Getting coefs (using OS) IN DISCOVERY DATA: 
#surv <- Surv(time = surv_table_extended[Discovery_set,]$OS, event = surv_table_extended[Discovery_set,]$status) # OS
#cox_reg_model_uni_random_signature <- coxph(formula = surv ~ ., data = surv_table_extended[Discovery_set,][,random_signature_miRNAs]) # univariate Cox to obtain coefs
#random_signature_coef <- cox_reg_model_uni_random_signature$coefficients
#stopifnot(all(names(random_signature_coef) == random_signature_miRNAs)) # check
#random_signature_coef <- unname(random_signature_coef)
#
#summary(cox_reg_model_uni_random_signature) #we can check model of 37 random mirnas
#
#random_signature <- data.frame(miRNA = random_signature_miRNAs,
#                         Coefficient = random_signature_coef)
#
#write.table(random_signature, paste0(save_dir, "/random_signature.txt"), sep = "\t", row.names = F)

random_signature <- read.delim(paste0(save_dir, "/random_signature.txt"))
random_signature_miRNAs <- random_signature$miRNA
random_signature_coef <- random_signature$Coefficient



# Getting Score for each patient: sum of all coef * miRNA:
random_signature_scores <- as.data.frame(colSums(expr[random_signature$miRNA,]*random_signature$Coefficient))
colnames(random_signature_scores) <- "random_signature_score"


### Adding miR24 score to surv_table:
stopifnot(all(rownames(surv_table)==rownames(random_signature_scores))) # check
stopifnot(all(rownames(surv_table_extended)==rownames(random_signature_scores))) # check
surv_table$random_signature_score <- random_signature_scores$random_signature_score
surv_table_extended$random_signature_score <- random_signature_scores$random_signature_score


#----


#### Select all miRNAs and scores to analyze ####

miRNAs <- c(sel_miRNAs, "AMLmiR36_score", "miR3_score", "miR4_score", "miR24_score", "miR37_score", "random_signature_score")

#----

#### Categorizamos expresion (maxstat) ####

## Podemos pasar la expresión a categórica ("low" and "high"):
## Separamos los grupos por método estadístico de maxstat:
## https://cran.r-project.org/web/packages/maxstat/vignettes/maxstat.pdf
#library(maxstat)
#library(survival)
#
##miRNAs <- c(sel_miRNAs, "AMLmiR36_score", "miR3_score", "miR4_score", "miR24_score")
#
#miRNA_results <- data.frame(row.names = miRNAs)
#miRNA_results$cut_pval <- NA
#miRNA_results$LogRankTest_pval_OS <- NA
#miRNA_results$LogRankTest_pval_EFS <- NA
#
#surv_table_categorical <- surv_table
#for (i in 1:length(miRNAs)) {
#  mt <- maxstat.test(Surv(OS, status) ~ surv_table[,miRNAs[i]], data=surv_table, smethod="LogRank", pmethod="condMC") # Cortamos con el método de p-value aproximation "conditional Monte-Carlo. Basandonos en el OS.
#  surv_table_categorical[,miRNAs[i]] <- cut(surv_table_categorical[,miRNAs[i]],
#                                            breaks = c(-Inf, mt$estimate, Inf),
#                                            labels = c("low", "high")) # Los valores de expresion igual al cut se asignan a low.
#  miRNA_results$cut_pval[i] <- mt$p.value
#}
#
#
#surv_table_categorical[,miRNAs] 
#miRNA_results
#
#
#
#
#
#### Create categorical extended data set:
#surv_table_categorical_extended <- merge(surv_vars_extended,
#                                         #texpr[,miRNAs], #si quisieramos hacerlo solo con los miRNAs seleccionados
#                                         surv_table_categorical[,sel_miRNAs],
#                                         by = "row.names")
#rownames(surv_table_categorical_extended) <- surv_table_categorical_extended$Row.names
##utils::View(surv_table_categorical_extended)
#
#
#
##Check
#m <- 6
#l <- surv_table_categorical[surv_table_categorical[,sel_miRNAs[m]]=="low",c("OS", "Vital_Status", sel_miRNAs[m])]
#h <- surv_table_categorical[surv_table_categorical[,sel_miRNAs[m]]=="high",c("OS", "Vital_Status", sel_miRNAs[m])]
#table(l$Vital_Status)
#table(h$Vital_Status)
#summary(l)
#summary(h)

#----


#### Categorizamos expresion (median) ####

# Podemos pasar la expresión a categórica ("low" and "high"):
# Separamos los grupos por mediana:


miRNA_results <- data.frame(row.names = miRNAs)
miRNA_results$LogRankTest_pval_OS <- NA
miRNA_results$LogRankTest_pval_EFS <- NA

surv_table_categorical <- surv_table
for (i in 1:length(miRNAs)) {
  median <- median(surv_table[,miRNAs[i]], na.rm = T)
  surv_table_categorical[,miRNAs[i]] <- cut(surv_table_categorical[,miRNAs[i]],
                                            breaks = c(-Inf, median, Inf),
                                            labels = c("low", "high")) # Los valores de expresion igual al cut se asignan a low.
}

surv_table_categorical[,miRNAs] 
miRNA_results





### Create categorical extended data set:
surv_table_categorical_extended <- merge(surv_vars_extended,
                                         #texpr[,miRNAs], #si quisieramos hacerlo solo con los miRNAs seleccionados
                                         surv_table_categorical[,sel_miRNAs],
                                         by = "row.names")
rownames(surv_table_categorical_extended) <- surv_table_categorical_extended$Row.names
#utils::View(surv_table_categorical_extended)



#Check
m <- 6
l <- surv_table_categorical[surv_table_categorical[,sel_miRNAs[m]]=="low",c("OS", "Vital_Status", sel_miRNAs[m])]
h <- surv_table_categorical[surv_table_categorical[,sel_miRNAs[m]]=="high",c("OS", "Vital_Status", sel_miRNAs[m])]
table(l$Vital_Status)
table(h$Vital_Status)
summary(l)
summary(h)

#----


#### Compare outcomes across patient sets ####

surv_table_extended$OS
surv_table_extended$EFS
surv_table_extended$status
table(surv_table_extended$Protocol)
table(surv_table_extended$Risk_group)

x <- data.frame(timeOS = surv_table_extended$OS,
                timeEFS = surv_table_extended$EFS,
                status = surv_table_extended$status,
                Protocol = surv_table_extended$Protocol,
                Risk_Group = surv_table_extended$Risk_group,
                PatientSet = surv_table_extended$PatientSet)
x$timeOS <- (x$timeOS+1)/365
x$timeEFS <- (x$timeEFS+1)/365
fitOSprotocol <- survfit(Surv(timeOS, status) ~ Protocol, data = x)
fitEFSprotocol <- survfit(Surv(timeEFS, status) ~ Protocol, data = x)
fitOSpatientSet <- survfit(Surv(timeOS, status) ~ PatientSet, data = x)
fitEFSpatientSet <- survfit(Surv(timeEFS, status) ~ PatientSet, data = x)
fitOSriskGroup <- survfit(Surv(timeOS, status) ~ Risk_Group, data = x)
fitEFSriskGroup <- survfit(Surv(timeEFS, status) ~ Risk_Group, data = x)

custom_theme <- function() { # Custom theme para poner el titulo en el centro
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5, vjust=4, family = "ArialMT"),
      plot.margin = unit(c(0.4,0.4,0.4,0.4), "cm")
    )
}
plot_protocolOS <- ggsurvplot(fitOSprotocol,
                   risk.table = TRUE, 
                   pval = TRUE,
                   pval.coord = c(1, 0.25),
                   conf.int = TRUE,
                   break.x.by = 2,
                   surv.median.line = "none", # o "hv"
                   #palette = c(),
                   axes.offset = F,
                   xlim = c(0,max(x$timeOS)),
                   xlab = "Time (years)",
                   ylab = "OS",
                   ggtheme = custom_theme(),
                   tables.theme = theme_cleantable(),
                   surv.scale = "percent",
                   risk.table.y.text = FALSE,
                   tables.y.text = FALSE, 
                   legend = "right",
                   legend.title="",
                   risk.table.title = "",
                   title = "Survival by protocol"
)
plot_protocolEFS <- ggsurvplot(fitEFSprotocol,
                              risk.table = TRUE, 
                              pval = TRUE,
                              pval.coord = c(1, 0.25),
                              conf.int = TRUE,
                              break.x.by = 2,
                              surv.median.line = "none", # o "hv"
                              #palette = c(),
                              axes.offset = F,
                              xlim = c(0,max(x$timeEFS)),
                              xlab = "Time (years)",
                              ylab = "EFS",
                              ggtheme = custom_theme(),
                              tables.theme = theme_cleantable(),
                              surv.scale = "percent",
                              risk.table.y.text = FALSE,
                              tables.y.text = FALSE, 
                              legend = "right",
                              legend.title="",
                              risk.table.title = "",
                              title = "Survival by protocol"
)
plot_patientSetOS <- ggsurvplot(fitOSpatientSet,
                               risk.table = TRUE, 
                               pval = TRUE,
                               pval.coord = c(1, 0.25),
                               conf.int = TRUE,
                               break.x.by = 2,
                               surv.median.line = "none", # o "hv"
                               #palette = c(),
                               axes.offset = F,
                               xlim = c(0,max(x$timeOS)),
                               xlab = "Time (years)",
                               ylab = "OS",
                               ggtheme = custom_theme(),
                               tables.theme = theme_cleantable(),
                               surv.scale = "percent",
                               risk.table.y.text = FALSE,
                               tables.y.text = FALSE, 
                               legend = "right",
                               legend.title="",
                               risk.table.title = "",
                               title = "Survival by patient set"
)
plot_patientSetEFS <- ggsurvplot(fitEFSpatientSet,
                                risk.table = TRUE, 
                                pval = TRUE,
                                pval.coord = c(1, 0.25),
                                conf.int = TRUE,
                                break.x.by = 2,
                                surv.median.line = "none", # o "hv"
                                #palette = c(),
                                axes.offset = F,
                                xlim = c(0,max(x$timeEFS)),
                                xlab = "Time (years)",
                                ylab = "EFS",
                                ggtheme = custom_theme(),
                                tables.theme = theme_cleantable(),
                                surv.scale = "percent",
                                risk.table.y.text = FALSE,
                                tables.y.text = FALSE, 
                                legend = "right",
                                legend.title="",
                                risk.table.title = "",
                                title = "Survival by patient set"
)
plot_RiskGroupOS <- ggsurvplot(fitOSriskGroup,
                              risk.table = TRUE, 
                              pval = TRUE,
                              pval.coord = c(1, 0.25),
                              conf.int = TRUE,
                              break.x.by = 2,
                              surv.median.line = "none", # o "hv"
                              #palette = c(),
                              axes.offset = F,
                              xlim = c(0,max(x$timeOS)),
                              xlab = "Time (years)",
                              ylab = "OS",
                              ggtheme = custom_theme(),
                              tables.theme = theme_cleantable(),
                              surv.scale = "percent",
                              risk.table.y.text = FALSE,
                              tables.y.text = FALSE, 
                              legend = "right",
                              legend.title="",
                              risk.table.title = "",
                              title = "Survival by Risk Group"
)
plot_RiskGroupEFS <- ggsurvplot(fitEFSriskGroup,
                               risk.table = TRUE, 
                               pval = TRUE,
                               pval.coord = c(1, 0.25),
                               conf.int = TRUE,
                               break.x.by = 2,
                               surv.median.line = "none", # o "hv"
                               #palette = c(),
                               axes.offset = F,
                               xlim = c(0,max(x$timeEFS)),
                               xlab = "Time (years)",
                               ylab = "EFS",
                               ggtheme = custom_theme(),
                               tables.theme = theme_cleantable(),
                               surv.scale = "percent",
                               risk.table.y.text = FALSE,
                               tables.y.text = FALSE, 
                               legend = "right",
                               legend.title="",
                               risk.table.title = "",
                               title = "Survival by Risk Group"
)
while (!is.null(dev.list()))  dev.off()
pdf(paste0(save_dir, "/Groups_SurvivalPlots.pdf"), width = 8, height = 5, family="ArialMT")
plot_protocolOS
plot_protocolEFS
plot_patientSetOS
plot_patientSetEFS
plot_RiskGroupOS
plot_RiskGroupEFS
while (!is.null(dev.list()))  dev.off()


#----

#### Kaplan Meier plot ####


#miRNAs <- c(sel_miRNAs, "AMLmiR36_score", "miR3_score", "miR4_score", "miR24_score")


custom_theme <- function() { # Custom theme para poner el titulo en el centro
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5, vjust=4, family = "ArialMT"),
      plot.margin = unit(c(0.4,0.4,0.4,0.4), "cm")
    )
}


# OS

p <- list()
for (i in 1:length(miRNAs)) {
  x <- data.frame(time = surv_table_categorical$OS,
                  status = surv_table_categorical$status,
                  expression = surv_table_categorical[,miRNAs[i]])
  x$time <- x$time+1
  x$time <- x$time/365
  fit <- survfit(Surv(time, status) ~ expression, data = x)
  #print(fit)
  #summary(fit)
  #summary(fit)$table
  
  # Log Rank Test:
  log_rank_test <- survdiff(Surv(time, status) ~ expression, data = x)
  log_rank_test
  p.val <- 1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1)
  miRNA_results$LogRankTest_pval_OS[i] <- p.val # Mismo pval que aparece en el plot
  
  ### PLOT ###
  plot <- ggsurvplot(fit,
                     risk.table = TRUE, 
                     pval = TRUE,
                     pval.coord = c(1, 0.25),
                     conf.int = TRUE,
                     #legend.labs=c("", ""),  
                     legend = "none", # No legend
                     #size=0.7,
                     break.x.by = 2,
                     surv.median.line = "none", # o "hv"
                     palette = c("#E7B800", "#2E9FDF"),
                     axes.offset = F,
                     xlim = c(0,max(x$time)),
                     xlab = "",
                     ylab = "",
                     ggtheme = custom_theme(),
                     tables.theme = theme_cleantable(),
                     surv.scale = "percent",
                     #tables.col = "strata",
                     #risk.table.col = "strata",
                     risk.table.y.text = FALSE,
                     tables.y.text = FALSE, 
                     #legend.title="",
                     risk.table.title = "",
                     title = gsub("_", "-", miRNAs[i])
  )
  p[[i]] <- plot$plot
}


while (!is.null(dev.list()))  dev.off()
pdf(paste0(save_dir, "/", data_subset, "_", "SurvivalPlots_OS.pdf"), width = 12, height = 17, family="ArialMT"
) #export as pdf
library(gridExtra)
args <- c(p, list(ncol = 3, left = "Overall Survival", bottom = "Time (Years)"))
do.call(grid.arrange, args)
while (!is.null(dev.list()))  dev.off()


#while (!is.null(dev.list()))  dev.off()
#library(gridExtra)
#args <- c(p, list(ncol = 3, left = "Overall Survival", bottom = "Time (Years)"))
#plot <- do.call(grid.arrange, args)
#ggsave(paste0(save_dir, "/", data_subset, "_", "SurvivalPlots_OS.eps"),
#       plot = plot,
#       width = 12,
#       height = 20) #export as eps
#while (!is.null(dev.list()))  dev.off()



# EFS

p <- list()
for (i in 1:length(miRNAs)) {
  x <- data.frame(time = surv_table_categorical$EFS,
                  status = surv_table_categorical$status,
                  expression = surv_table_categorical[,miRNAs[i]])
  x$time <- x$time+1
  x$time <- x$time/365
  fit <- survfit(Surv(time, status) ~ expression, data = x)
  #print(fit)
  #summary(fit)
  #summary(fit)$table
  
  # Log Rank Test:
  log_rank_test <- survdiff(Surv(time, status) ~ expression, data = x)
  log_rank_test
  p.val <- 1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1)
  miRNA_results$LogRankTest_pval_EFS[i] <- p.val # Mismo pval que aparece en el plot
  
  ### PLOT ###
  plot <- ggsurvplot(fit,
                     risk.table = TRUE, 
                     pval = TRUE,
                     pval.coord = c(1, 0.25),
                     conf.int = TRUE,
                     #legend.labs=c("", ""),  
                     legend = "none", # No legend
                     break.x.by = 2,
                     surv.median.line = "none", # o "hv
                     palette = c("#E7B800", "#2E9FDF"),
                     axes.offset = F,
                     xlim = c(0,max(x$time)),
                     xlab = "",
                     ylab = "",
                     ggtheme = custom_theme(),
                     tables.theme = theme_cleantable(),
                     surv.scale = "percent",
                     #tables.col = "strata",
                     #risk.table.col = "strata",
                     risk.table.y.text = FALSE,
                     tables.y.text = FALSE, 
                     legend.title="",
                     risk.table.title = "",
                     title = gsub("_", "-", miRNAs[i])
  )
  p[[i]] <- plot$plot
}


while (!is.null(dev.list()))  dev.off()
pdf(paste0(save_dir, "/", data_subset, "_", "SurvivalPlots_EFS.pdf"), width = 12, height = 17, family="ArialMT"
) #export as pdf
library(gridExtra)
args <- c(p, list(ncol = 3, left = "Event Free Survival", bottom = "Time (Years)"))
do.call(grid.arrange, args)
while (!is.null(dev.list()))  dev.off()


#while (!is.null(dev.list()))  dev.off()
#library(gridExtra)
#args <- c(p, list(ncol = 3, left = "Event Free Survival", bottom = "Time (Years)"))
#plot <- do.call(grid.arrange, args)
#ggsave(paste0(save_dir, "/", data_subset, "_", "SurvivalPlots_OS.eps"),
#       plot = plot,
#       width = 12,
#       height = 20) #export as eps
#while (!is.null(dev.list()))  dev.off()



miRNA_results$LogRankTest_pval_OS <- round(miRNA_results$LogRankTest_pval_OS, 3)
miRNA_results$LogRankTest_pval_EFS <- round(miRNA_results$LogRankTest_pval_EFS, 3)
miRNA_results

#----

#### Cox Proportional Hazards model: ####


# Function to calculate and plot univariate and multivariate coxph
perform_uni_multi_coxph <- function(data, time_var, status_var, predictors, covars, save_table = F, save_table_file = NULL){
  
  ### Fit coxph
  
  surv <- Surv(time = data[,time_var], event = data[,status_var])
  cox <- data.frame("predictors"=predictors,
                    "beta_uni"=NA,
                    "HR_uni"=NA,
                    "HR.CI_lower_uni"=NA,
                    "HR.CI_upper_uni"=NA,
                    "wald.test_uni"=NA,
                    "p.value_uni"=NA,
                    "beta_multi"=NA,
                    "HR_multi"=NA,
                    "HR.CI_lower_multi"=NA,
                    "HR.CI_upper_multi"=NA,
                    "wald.test_multi"=NA,
                    "p.value_multi"=NA)
  for (i in 1:length(predictors)) {
    # Univariate cox analysis:
    cox_reg_model_uni <- coxph(formula = surv ~ data[,predictors[i]], data = data)
    # Multivariate cox analysis:
    cox_reg_model_multi <- coxph(formula = as.formula(paste("surv ~ data[,predictors[i]] +", paste(covars, collapse = " + "))),
                                 data = data) 
    
    unimod <- summary(cox_reg_model_uni)
    cox$beta_uni[i] <- signif(unimod$coef[1], digits=2);#coeficient beta
    cox$HR_uni[i] <- signif(unimod$coef[2], digits=2);#exp(beta)
    cox$HR.CI_lower_uni[i] <- signif(unimod$conf.int[,"lower .95"], 2)
    cox$HR.CI_upper_uni[i] <- signif(unimod$conf.int[,"upper .95"],2)
    cox$wald.test_uni[i] <- signif(unimod$wald["test"], digits=2)
    cox$p.value_uni[i] <- signif(unimod$coef[,"Pr(>|z|)"], digits=2)
    
    multimod <- summary(cox_reg_model_multi)
    cox$beta_multi[i] <- signif(multimod$coef[1,1], digits=2);#coeficient beta
    cox$HR_multi[i] <- signif(multimod$coef[1,2], digits=2);#exp(beta)
    cox$HR.CI_lower_multi[i] <- signif(multimod$conf.int[1,"lower .95"], 2)
    cox$HR.CI_upper_multi[i] <- signif(multimod$conf.int[1,"upper .95"],2)
    cox$wald.test_multi[i] <- signif(multimod$wald["test"], digits=2)
    cox$p.value_multi[i] <- signif(multimod$coef[1,"Pr(>|z|)"], digits=2)
  }
  
  #cox[,unlist(lapply(cox, is.numeric))] <- round(cox[,unlist(lapply(cox, is.numeric))], 3) #round numeric cols
  cox
  
  # Coef exlanation: Positive coef means that an increase in the predictor variable is associated with an increased hazard of the event
  # Positive coef: more expression more risk. 
  # Negative coef: less expression more risk.
  
  # Save
  if (save_table) {
    write.table(cox, file = save_table_file, sep = "\t", row.names = F)
  }
  
  
  
  ### Plot
  
  library(forestploter)
  library(grid)
  
  cox.fig <- data.frame(gsub("_", "-", cox$predictors), # Sustituimos _ solo para aplicarlo a miRNAs !!!
                        paste(rep(" ", 20), collapse = " "),
                        paste0(cox$HR_uni, " [", 
                               cox$HR.CI_lower_uni, "-", 
                               cox$HR.CI_upper_uni, "]"),
                        cox$p.value_uni,
                        "  ",
                        paste(rep(" ", 20), collapse = " "),
                        paste0(cox$HR_multi, " [", 
                               cox$HR.CI_lower_multi, "-", 
                               cox$HR.CI_upper_multi, "]"),
                        cox$p.value_multi)
  
  cox.fig[,4][cox.fig[,4]>=0.001] <- round(cox.fig[,4][cox.fig[,4]>=0.001], 3)
  cox.fig[,4][cox.fig[,4]<0.001] <- "<0.001"
  cox.fig[,8][cox.fig[,8]>=0.001] <- round(cox.fig[,8][cox.fig[,8]>=0.001], 3)
  cox.fig[,8][cox.fig[,8]<0.001] <- "<0.001"
  
  
  colnames(cox.fig) <- c("", "", "HR [95% CI]", "p-value", "", "", "HR [95% CI]", "p-value")
  
  p <- forest(cox.fig,
              est = list(cox$HR_uni, cox$HR_multi),
              lower = list(cox$HR.CI_lower_uni, cox$HR.CI_lower_multi), 
              upper = list(cox$HR.CI_upper_uni, cox$HR.CI_upper_multi),
              ci_column = c(2, 6),
              ref_line = 1,
              xlim = c(0, max(cox$HR.CI_upper_uni)+0.5),
              ticks_at = c(0, 0.5, 1, 2),
              xlab = "Hazard Ratio",
              theme = forest_theme(core = list(bg_params=list(fill = c("white"))),
                                   arrow_label_just = "end",
                                   arrow_type = "closed"),
              nudge_y = 0.2)
  p <- insert_text(p, text = "Univariate",
                   part = "header",
                   row = 1,
                   col = 2:4,
                   gp = gpar(fontface = "bold"))
  p <- add_text(p, text = "Multivariate",
                part = "header",
                row = 1,
                col = 6:8,
                gp = gpar(fontface = "bold"))
  p <- add_border(p, 
                  part = "header", 
                  row = 1,
                  col = 2:4,
                  gp = gpar(lwd = 1))
  p <- add_border(p, 
                  part = "header", 
                  row = 1,
                  col = 6:8,
                  gp = gpar(lwd = 1))
  p <- add_border(p, 
                  part = "header", 
                  row = 2,
                  gp = gpar(lwd = 1))
  p <- edit_plot(p, col = c(3,4,7,8), 
                 which = "text",
                 hjust = unit(0.5, "npc"),
                 x = unit(0.5, "npc"))
  p <- insert_text(p, text = time_var,
                   part = "header",
                   row = 1,
                   col = 2:8,
                   gp = gpar(fontface = "bold"))
  p <- add_border(p, 
                  part = "header", 
                  row = 1,
                  col = 2:8,
                  gp = gpar(lwd = 1))
  return(p)
}


# Multivariate variables to choose:
#colnames(surv_vars_extended)
#colnames(surv_table_extended[1:28])
#table(surv_table_extended$status, surv_table_extended$CNS_disease)
#table(surv_table_extended$status, surv_table_extended$Chloroma)
#table(surv_table_extended$status, surv_table_extended$FAB_Category)
#table(surv_table_extended$status, surv_table_extended$SCT_in_1st_CR)
#table(surv_table_extended$status, surv_table_extended$FLT3.ITD_positive.)

coxcovars <- c("WBC_at_Diagnosis",
               "SCT_in_1st_CR",
               "FLT3.ITD_positive.",
               "CNS_disease")

predictors <- c("miR37_score", "AMLmiR36_score", "miR24_score", "miR4_score", "miR3_score",# "random_signature_score", 
                "hsa_mir_199a", "hsa_mir_335", "hsa_mir_196b", "hsa_mir_100", "hsa_mir_192", 
                "hsa_miR_139_5p", "hsa_mir_34b", "hsa_mir_195", "hsa_mir_381", "hsa_mir_370", "hsa_miR_193b_3p", 
                "hsa_miR_206", "hsa_mir_122", "hsa_mir_29a", "hsa_miR_375", "hsa_mir_155")

# Apply function:


# OS
p <- perform_uni_multi_coxph(
  data = surv_table_extended,
  time_var = "OS",
  status_var = "status",
  predictors = predictors,
  covars = coxcovars,
  save_table = T,
  save_table_file = paste0(save_dir, "/", data_subset, "_coxph_OS.txt")
)

while (!is.null(dev.list()))  dev.off()
pdf(paste0(save_dir, "/", data_subset, "_coxph_OS.pdf"), width = 12, height = 7, family="ArialMT")
p
while (!is.null(dev.list()))  dev.off()



# EFS
p <- perform_uni_multi_coxph(
  data = surv_table_extended,
  time_var = "EFS",
  status_var = "status",
  predictors = predictors,
  covars = coxcovars,
  save_table = T,
  save_table_file = paste0(save_dir, "/", data_subset, "_coxph_EFS.txt")
)

while (!is.null(dev.list()))  dev.off()
pdf(paste0(save_dir, "/", data_subset, "_coxph_EFS.pdf"), width = 12, height = 7, family="ArialMT")
p
while (!is.null(dev.list()))  dev.off()
















### EFS
survival_type <- "EFS"
surv_vars_extended[,survival_type]
#miRNAs <- c(sel_miRNAs, "AMLmiR36_score", "miR3_score", "miR4_score", "miR24_score")
surv <- Surv(time = surv_vars_extended$EFS, event = surv_vars_extended$status) # EFS
miRNA_cox <- data.frame("miRNA"=miRNAs,
                        "beta_uni"=NA,
                        "HR_uni"=NA,
                        "HR.CI_lower_uni"=NA,
                        "HR.CI_upper_uni"=NA,
                        "wald.test_uni"=NA,
                        "p.value_uni"=NA,
                        "beta_multi"=NA,
                        "HR_multi"=NA,
                        "HR.CI_lower_multi"=NA,
                        "HR.CI_upper_multi"=NA,
                        "wald.test_multi"=NA,
                        "p.value_multi"=NA)
for (i in 1:length(miRNAs)) {
  # Univariate cox analysis:
  cox_reg_model_uni <- coxph(formula = surv ~ surv_table_extended[,miRNAs[i]], data = surv_table_extended)
  # Multivariate cox analysis:
  cox_reg_model_multi <- coxph(formula = surv ~ surv_table_extended[,miRNAs[i]] 
                               #+ Risk_group # remove when selecting a risk group
                               #+ Age_at_Diagnosis_in_Days 
                               #+ Gender              # not related to outcome
                               + WBC_at_Diagnosis 
                               + SCT_in_1st_CR
                               #+ FAB_Category        # FAB outdated
                               #+ MLL                 # mutations not included: correlated to risk group
                               + FLT3.ITD_positive.
                               #+ NPM_mutation
                               #+ CEBPA_mutation
                               #+ Primary_Cytogenetic_Code
                               #+ Bone_marrow_leukemic_blast_percentage_...
                               #+ Peripheral_blasts_...
                               + CNS_disease
                               #+ Chloroma
                               ,data = surv_table_extended) 
  
  unimod <- summary(cox_reg_model_uni)
  miRNA_cox$beta_uni[i] <- signif(unimod$coef[1], digits=2);#coeficient beta
  miRNA_cox$HR_uni[i] <- signif(unimod$coef[2], digits=2);#exp(beta)
  miRNA_cox$HR.CI_lower_uni[i] <- signif(unimod$conf.int[,"lower .95"], 2)
  miRNA_cox$HR.CI_upper_uni[i] <- signif(unimod$conf.int[,"upper .95"],2)
  miRNA_cox$wald.test_uni[i] <- signif(unimod$wald["test"], digits=2)
  miRNA_cox$p.value_uni[i] <- signif(unimod$wald["pvalue"], digits=2)
  
  multimod <- summary(cox_reg_model_multi)
  miRNA_cox$beta_multi[i] <- signif(multimod$coef[1,1], digits=2);#coeficient beta
  miRNA_cox$HR_multi[i] <- signif(multimod$coef[1,2], digits=2);#exp(beta)
  miRNA_cox$HR.CI_lower_multi[i] <- signif(multimod$conf.int[1,"lower .95"], 2)
  miRNA_cox$HR.CI_upper_multi[i] <- signif(multimod$conf.int[1,"upper .95"],2)
  miRNA_cox$wald.test_multi[i] <- signif(multimod$wald["test"], digits=2)
  miRNA_cox$p.value_multi[i] <- signif(multimod$wald["pvalue"], digits=2)
}
#miRNA_cox[,unlist(lapply(miRNA_cox, is.numeric))] <- round(miRNA_cox[,unlist(lapply(miRNA_cox, is.numeric))], 3) #round numeric cols
miRNA_cox

# Coef exlanation: Positive coef means that an increase in the predictor variable is associated with an increased hazard of the event
# Positive coef: more expression more risk. 
# Negative coef: less expression more risk.

# Save
write.table(miRNA_cox, file = paste0(save_dir, "/", data_subset, "_coxph_EFS.txt"), sep = "\t", row.names = F)


#----


#### Cox Proportional Hazards model (explore covariables): ####


# CoxPH to single predictor:

predictors <- c("miR37_score", "AMLmiR36_score", "miR24_score", "miR4_score", "miR3_score", "random_signature_score", 
                "hsa_mir_199a", "hsa_mir_335", "hsa_mir_196b", "hsa_mir_100", "hsa_mir_192", 
                "hsa_miR_139_5p", "hsa_mir_34b", "hsa_mir_195", "hsa_mir_381", "hsa_mir_370", "hsa_miR_193b_3p", 
                "hsa_miR_206", "hsa_mir_122", "hsa_mir_29a", "hsa_miR_375", "hsa_mir_155")
#pred <- predictors[1]
data <- surv_table_extended
dim(data)
time_var <- "EFS"#"OS"#  # CHANGE FOR RE-ANALYSIS
status_var <- "status"

covars <- c("WBC_at_Diagnosis", # Numeric - 0% NAs
            "SCT_in_1st_CR",
            "FLT3.ITD_positive.",
            "CNS_disease",
            #"MRD_at_end_of_course_1_perc", # Numeric - 10% NAs
            "MRD_at_end_of_course_1", # 10% NAs
            #"MRD_at_end_of_course_2_perc", # Numeric -  32% NAs (too many) --> quit
            #"MRD_at_end_of_course_2", # 25% NAs (too many) --> quit
            #"Primary_Cytogenetic_Code", # 2% NAs but redundant with other variables --> quit
            "NPM_mutation",
            "CEBPA_mutation",
            #"WT1_mutation", # 65% NAs --> quit
            #"FLT3_PM", # 65% NAs --> quit
            #"c_Kit_Mutation_Exon_8", # 90% NAs --> quit
            #"c_Kit_Mutation_Exon_17", # 90% NAs --> quit
            colnames(clin[,25:42]) #cytogenetic abnormalities
            )


# Explore NAs in data:
for (i in 1:length(covars)) {
  cat(paste0(covars[i], ":"))
  if (!is.numeric(data[,covars[i]])) {
    print(table(data[,covars[i]], exclude = NULL))
  }
  print(paste0("Prop NAs: ", mean(is.na(data[,covars[i]]))))
  cat("\n")
}

hist(data$MRD_at_end_of_course_1_perc)
table(data$MRD_at_end_of_course_1_perc==0, exclude = NULL)
table(data$MRD_at_end_of_course_1, exclude = NULL)

## Change NA to Unknown:
#for (i in 1:length(covars)) {
#  if (!is.numeric(data[,covars[i]])) {
#    data[,covars[i]][is.na(data[,covars[i]])] <- "Unknown"
#  }
#}

# Change Unknown to NA:
for (i in 1:length(covars)) {
  if (!is.numeric(data[,covars[i]])) {
    data[,covars[i]][data[,covars[i]] == "Unknown"] <- NA
  }
}

# Remove covars with less than 20 patients with Yes:
newcovars <- covars
for (i in 1:length(covars)) {
  if (!is.numeric(data[,covars[i]])) {
    if (sum(data[,covars[i]] == "Yes", na.rm = T) < 20) {
      print(paste0(covars[i], " removed"))
      newcovars <- newcovars[newcovars != covars[i]]
    }
  }
}
covars <- newcovars


surv <- Surv(time = data[,time_var], event = data[,status_var])

# Multivariate cox analysis:
#cox_reg_model_multi <- coxph(formula = as.formula(paste("surv ~ ", pred, "+", paste(covars, collapse = " + "))),
#                              data = data) 
#cox_reg_model_multi_allPreds <- coxph(formula = as.formula(paste("surv ~ ", paste(c(predictors, covars), collapse = " + "))),
#                             data = data) 
#cox_reg_model_uni_allPreds_noCovars <- coxph(formula = as.formula(paste("surv ~ ", paste(predictors, collapse = " + "))),
#                                      data = data[Validation_set,]) 
cox_reg_model_onlyCovars <- coxph(formula = as.formula(paste("surv ~ ", paste(covars, collapse = " + "))),
                                             data = data) 

#summary(cox_reg_model_multi)
#summary(cox_reg_model_multi_allPreds)
#summary(cox_reg_model_uni_allPreds_noCovars)
multimod <- summary(cox_reg_model_onlyCovars)

cox <- data.frame("predictors"=covars,
                  "beta_multi"=NA,
                  "HR_multi"=NA,
                  "HR.CI_lower_multi"=NA,
                  "HR.CI_upper_multi"=NA,
                  "wald.test_multi"=NA,
                  "p.value_multi"=NA)
cox$beta_multi <- signif(multimod$coef[,1], digits=2);#coeficient beta
cox$HR_multi <- signif(multimod$coef[,2], digits=2);#exp(beta)
cox$HR.CI_lower_multi <- signif(multimod$conf.int[,"lower .95"], 2)
cox$HR.CI_upper_multi <- signif(multimod$conf.int[,"upper .95"],2)
cox$p.value_multi <- signif(multimod$coef[,5], digits=2)

# Fix covars names:
cox$predictors <- c("WBC at Diagnosis",
                    "SCT in 1st CR",
                    "FLT3 ITD positive",
                    "CNS disease",
                    "MRD at end of course 1",
                    "NPM mutation",
                    "CEBPA mutation",
                    "t(6;9)",
                    "t(8;21)",
                    "t(6;11)(q27;q23)",
                    "t(9;11)(p22;q23)",
                    "t(10;11)(p11.2;q23)",
                    "t(11:19)(q23:p13.1)",
                    "inv(16)",
                    "del7q",
                    "del9q",
                    "monosomy 7",
                    "trisomy 8",
                    "trisomy 21",
                    "MLL",
                    "Minus Y",
                    "Minus X")

# Plot:

library(forestploter)
library(grid)

cox.fig <- data.frame(cox$predictors,
                      paste(rep(" ", 20), collapse = " "),
                      paste0(cox$HR_multi, " [", 
                             cox$HR.CI_lower_multi, "-", 
                             cox$HR.CI_upper_multi, "]"),
                      cox$p.value_multi)

cox.fig[,4][cox.fig[,4]>=0.001] <- round(cox.fig[,4][cox.fig[,4]>=0.001], 3)
cox.fig[,4][cox.fig[,4]<0.001] <- "<0.001"

colnames(cox.fig) <- c("", "", "HR [95% CI]", "p-value")

p <- forest(cox.fig,
            est = cox$HR_multi,
            lower = cox$HR.CI_lower_multi, 
            upper = cox$HR.CI_upper_multi,
            ci_column = 2,
            ref_line = 1,
            xlim = c(0, max(cox$HR.CI_upper_multi)+0.5),
            ticks_at = c(0, 0.5, 1, 2),
            xlab = "Hazard Ratio",
            theme = forest_theme(core = list(bg_params=list(fill = c("white"))),
                                 arrow_label_just = "end",
                                 arrow_type = "closed"),
            nudge_y = 0.2
            )

p <- add_border(p, 
                part = "header", 
                row = 1,
                col = 1:4,
                gp = gpar(lwd = 1))
p <- edit_plot(p, col = c(3,4), 
               which = "text",
               hjust = unit(0.5, "npc"),
               x = unit(0.5, "npc"))
p <- insert_text(p, text = time_var,
                 part = "header",
                 row = 1,
                 col = 2:4,
                 gp = gpar(fontface = "bold"))
p <- add_border(p, 
                part = "header", 
                row = 1,
                col = 2:4,
                gp = gpar(lwd = 1))
p


while (!is.null(dev.list()))  dev.off()
pdf(paste0(save_dir, "/", data_subset, "_coxph_", time_var, "_onlyCovars.pdf"), width = 12, height = 7, family="ArialMT")
p
while (!is.null(dev.list()))  dev.off()



#----


#### Cox Proportional Hazards model (adding all covariables): ####


# Function to calculate and plot univariate and multivariate coxph
perform_uni_multi_coxph <- function(data, time_var, status_var, predictors, covars, save_table = F, save_table_file = NULL){
  
  ### Fit coxph
  
  surv <- Surv(time = data[,time_var], event = data[,status_var])
  cox <- data.frame("predictors"=predictors,
                    "beta_uni"=NA,
                    "HR_uni"=NA,
                    "HR.CI_lower_uni"=NA,
                    "HR.CI_upper_uni"=NA,
                    "wald.test_uni"=NA,
                    "p.value_uni"=NA,
                    "beta_multi"=NA,
                    "HR_multi"=NA,
                    "HR.CI_lower_multi"=NA,
                    "HR.CI_upper_multi"=NA,
                    "wald.test_multi"=NA,
                    "p.value_multi"=NA)
  for (i in 1:length(predictors)) {
    # Univariate cox analysis:
    cox_reg_model_uni <- coxph(formula = surv ~ data[,predictors[i]], data = data)
    # Multivariate cox analysis:
    cox_reg_model_multi <- coxph(formula = as.formula(paste("surv ~ data[,predictors[i]] +", paste(covars, collapse = " + "))),
                                 data = data) 
    
    unimod <- summary(cox_reg_model_uni)
    cox$beta_uni[i] <- signif(unimod$coef[1], digits=2);#coeficient beta
    cox$HR_uni[i] <- signif(unimod$coef[2], digits=2);#exp(beta)
    cox$HR.CI_lower_uni[i] <- signif(unimod$conf.int[,"lower .95"], 2)
    cox$HR.CI_upper_uni[i] <- signif(unimod$conf.int[,"upper .95"],2)
    cox$wald.test_uni[i] <- signif(unimod$wald["test"], digits=2)
    cox$p.value_uni[i] <- signif(unimod$coef[,"Pr(>|z|)"], digits=2)
    
    multimod <- summary(cox_reg_model_multi)
    cox$beta_multi[i] <- signif(multimod$coef[1,1], digits=2);#coeficient beta
    cox$HR_multi[i] <- signif(multimod$coef[1,2], digits=2);#exp(beta)
    cox$HR.CI_lower_multi[i] <- signif(multimod$conf.int[1,"lower .95"], 2)
    cox$HR.CI_upper_multi[i] <- signif(multimod$conf.int[1,"upper .95"],2)
    cox$wald.test_multi[i] <- signif(multimod$wald["test"], digits=2)
    cox$p.value_multi[i] <- signif(multimod$coef[1,"Pr(>|z|)"], digits=2)
  }
  
  #cox[,unlist(lapply(cox, is.numeric))] <- round(cox[,unlist(lapply(cox, is.numeric))], 3) #round numeric cols
  cox
  
  # Coef exlanation: Positive coef means that an increase in the predictor variable is associated with an increased hazard of the event
  # Positive coef: more expression more risk. 
  # Negative coef: less expression more risk.
  
  # Save
  if (save_table) {
    write.table(cox, file = save_table_file, sep = "\t", row.names = F)
  }
  
  
  
  ### Plot
  
  library(forestploter)
  library(grid)
  
  cox.fig <- data.frame(gsub("_", "-", cox$predictors), # Sustituimos _ solo para aplicarlo a miRNAs !!!
                        paste(rep(" ", 20), collapse = " "),
                        paste0(cox$HR_uni, " [", 
                               cox$HR.CI_lower_uni, "-", 
                               cox$HR.CI_upper_uni, "]"),
                        cox$p.value_uni,
                        "  ",
                        paste(rep(" ", 20), collapse = " "),
                        paste0(cox$HR_multi, " [", 
                               cox$HR.CI_lower_multi, "-", 
                               cox$HR.CI_upper_multi, "]"),
                        cox$p.value_multi)
  
  cox.fig[,4][cox.fig[,4]>=0.001] <- round(cox.fig[,4][cox.fig[,4]>=0.001], 3)
  cox.fig[,4][cox.fig[,4]<0.001] <- "<0.001"
  cox.fig[,8][cox.fig[,8]>=0.001] <- round(cox.fig[,8][cox.fig[,8]>=0.001], 3)
  cox.fig[,8][cox.fig[,8]<0.001] <- "<0.001"
  
  
  colnames(cox.fig) <- c("", "", "HR [95% CI]", "p-value", "", "", "HR [95% CI]", "p-value")
  
  p <- forest(cox.fig,
              est = list(cox$HR_uni, cox$HR_multi),
              lower = list(cox$HR.CI_lower_uni, cox$HR.CI_lower_multi), 
              upper = list(cox$HR.CI_upper_uni, cox$HR.CI_upper_multi),
              ci_column = c(2, 6),
              ref_line = 1,
              xlim = c(0, max(cox$HR.CI_upper_uni)+0.5),
              ticks_at = c(0, 0.5, 1, 2),
              xlab = "Hazard Ratio",
              theme = forest_theme(core = list(bg_params=list(fill = c("white"))),
                                   arrow_label_just = "end",
                                   arrow_type = "closed"),
              nudge_y = 0.2)
  p <- insert_text(p, text = "Univariate",
                   part = "header",
                   row = 1,
                   col = 2:4,
                   gp = gpar(fontface = "bold"))
  p <- add_text(p, text = "Multivariate",
                part = "header",
                row = 1,
                col = 6:8,
                gp = gpar(fontface = "bold"))
  p <- add_border(p, 
                  part = "header", 
                  row = 1,
                  col = 2:4,
                  gp = gpar(lwd = 1))
  p <- add_border(p, 
                  part = "header", 
                  row = 1,
                  col = 6:8,
                  gp = gpar(lwd = 1))
  p <- add_border(p, 
                  part = "header", 
                  row = 2,
                  gp = gpar(lwd = 1))
  p <- edit_plot(p, col = c(3,4,7,8), 
                 which = "text",
                 hjust = unit(0.5, "npc"),
                 x = unit(0.5, "npc"))
  p <- insert_text(p, text = time_var,
                   part = "header",
                   row = 1,
                   col = 2:8,
                   gp = gpar(fontface = "bold"))
  p <- add_border(p, 
                  part = "header", 
                  row = 1,
                  col = 2:8,
                  gp = gpar(lwd = 1))
  return(p)
}



# Multivariate variables to choose:
#colnames(surv_vars_extended)
#colnames(surv_table_extended[1:28])
#table(surv_table_extended$status, surv_table_extended$CNS_disease)
#table(surv_table_extended$status, surv_table_extended$Chloroma)
#table(surv_table_extended$status, surv_table_extended$FAB_Category)
#table(surv_table_extended$status, surv_table_extended$SCT_in_1st_CR)
#table(surv_table_extended$status, surv_table_extended$FLT3.ITD_positive.)

coxcovars <- c("WBC_at_Diagnosis", # Numeric - 0% NAs
            "SCT_in_1st_CR",
            "FLT3.ITD_positive.",
            "CNS_disease",
            "MRD_at_end_of_course_1", # Numeric - 10% NAs
            #"MRD_at_end_of_course_2_perc", # Numeric -  32% NAs (too many) --> quit
            #"Primary_Cytogenetic_Code", # 2% NAs but redundant with other variables --> quit
            "NPM_mutation",
            "CEBPA_mutation",
            #"WT1_mutation", # 65% NAs --> quit
            #"FLT3_PM", # 65% NAs --> quit
            #"c_Kit_Mutation_Exon_8", # 90% NAs --> quit
            #"c_Kit_Mutation_Exon_17", # 90% NAs --> quit
            colnames(clin[,25:42]) #cytogenetic abnormalities
)

predictors <- c("miR37_score", "AMLmiR36_score", "miR24_score", "miR4_score", "miR3_score",# "random_signature_score", 
                "hsa_mir_199a", "hsa_mir_335", "hsa_mir_196b", "hsa_mir_100", "hsa_mir_192", 
                "hsa_miR_139_5p", "hsa_mir_34b", "hsa_mir_195", "hsa_mir_381", "hsa_mir_370", "hsa_miR_193b_3p", 
                "hsa_miR_206", "hsa_mir_122", "hsa_mir_29a", "hsa_miR_375", "hsa_mir_155")


# Apply function:

# OS
p <- perform_uni_multi_coxph(
  data = surv_table_extended,
  time_var = "OS",
  status_var = "status",
  predictors = predictors,
  covars = coxcovars,
  save_table = F
)

while (!is.null(dev.list()))  dev.off()
pdf(paste0(save_dir, "/", data_subset, "_coxph_OS_withAllCovars.pdf"), width = 12, height = 7, family="ArialMT")
p
while (!is.null(dev.list()))  dev.off()



# EFS
p <- perform_uni_multi_coxph(
  data = surv_table_extended,
  time_var = "EFS",
  status_var = "status",
  predictors = predictors,
  covars = coxcovars,
  save_table = F
)

while (!is.null(dev.list()))  dev.off()
pdf(paste0(save_dir, "/", data_subset, "_coxph_EFS_withAllCovars.pdf"), width = 12, height = 7, family="ArialMT")
p
while (!is.null(dev.list()))  dev.off()



#----


#### Cox Proportional Hazards model (adding significant covariables p<0.005): ####


# Function to calculate and plot univariate and multivariate coxph
perform_uni_multi_coxph <- function(data, time_var, status_var, predictors, covars, save_table = F, save_table_file = NULL){
  
  ### Fit coxph
  
  surv <- Surv(time = data[,time_var], event = data[,status_var])
  cox <- data.frame("predictors"=predictors,
                    "beta_uni"=NA,
                    "HR_uni"=NA,
                    "HR.CI_lower_uni"=NA,
                    "HR.CI_upper_uni"=NA,
                    "wald.test_uni"=NA,
                    "p.value_uni"=NA,
                    "beta_multi"=NA,
                    "HR_multi"=NA,
                    "HR.CI_lower_multi"=NA,
                    "HR.CI_upper_multi"=NA,
                    "wald.test_multi"=NA,
                    "p.value_multi"=NA)
  for (i in 1:length(predictors)) {
    # Univariate cox analysis:
    cox_reg_model_uni <- coxph(formula = surv ~ data[,predictors[i]], data = data)
    # Multivariate cox analysis:
    cox_reg_model_multi <- coxph(formula = as.formula(paste("surv ~ data[,predictors[i]] +", paste(covars, collapse = " + "))),
                                 data = data) 
    
    unimod <- summary(cox_reg_model_uni)
    cox$beta_uni[i] <- signif(unimod$coef[1], digits=2);#coeficient beta
    cox$HR_uni[i] <- signif(unimod$coef[2], digits=2);#exp(beta)
    cox$HR.CI_lower_uni[i] <- signif(unimod$conf.int[,"lower .95"], 2)
    cox$HR.CI_upper_uni[i] <- signif(unimod$conf.int[,"upper .95"],2)
    cox$wald.test_uni[i] <- signif(unimod$wald["test"], digits=2)
    cox$p.value_uni[i] <- signif(unimod$coef[,"Pr(>|z|)"], digits=2)
    
    multimod <- summary(cox_reg_model_multi)
    cox$beta_multi[i] <- signif(multimod$coef[1,1], digits=2);#coeficient beta
    cox$HR_multi[i] <- signif(multimod$coef[1,2], digits=2);#exp(beta)
    cox$HR.CI_lower_multi[i] <- signif(multimod$conf.int[1,"lower .95"], 2)
    cox$HR.CI_upper_multi[i] <- signif(multimod$conf.int[1,"upper .95"],2)
    cox$wald.test_multi[i] <- signif(multimod$wald["test"], digits=2)
    cox$p.value_multi[i] <- signif(multimod$coef[1,"Pr(>|z|)"], digits=2)
  }
  
  #cox[,unlist(lapply(cox, is.numeric))] <- round(cox[,unlist(lapply(cox, is.numeric))], 3) #round numeric cols
  cox
  
  # Coef exlanation: Positive coef means that an increase in the predictor variable is associated with an increased hazard of the event
  # Positive coef: more expression more risk. 
  # Negative coef: less expression more risk.
  
  # Save
  if (save_table) {
    write.table(cox, file = save_table_file, sep = "\t", row.names = F)
  }
  
  
  
  ### Plot
  
  library(forestploter)
  library(grid)
  
  cox.fig <- data.frame(gsub("_", "-", cox$predictors), # Sustituimos _ solo para aplicarlo a miRNAs !!!
                        paste(rep(" ", 20), collapse = " "),
                        paste0(cox$HR_uni, " [", 
                               cox$HR.CI_lower_uni, "-", 
                               cox$HR.CI_upper_uni, "]"),
                        cox$p.value_uni,
                        "  ",
                        paste(rep(" ", 20), collapse = " "),
                        paste0(cox$HR_multi, " [", 
                               cox$HR.CI_lower_multi, "-", 
                               cox$HR.CI_upper_multi, "]"),
                        cox$p.value_multi)
  
  cox.fig[,4][cox.fig[,4]>=0.001] <- round(cox.fig[,4][cox.fig[,4]>=0.001], 3)
  cox.fig[,4][cox.fig[,4]<0.001] <- "<0.001"
  cox.fig[,8][cox.fig[,8]>=0.001] <- round(cox.fig[,8][cox.fig[,8]>=0.001], 3)
  cox.fig[,8][cox.fig[,8]<0.001] <- "<0.001"
  
  
  colnames(cox.fig) <- c("", "", "HR [95% CI]", "p-value", "", "", "HR [95% CI]", "p-value")
  
  p <- forest(cox.fig,
              est = list(cox$HR_uni, cox$HR_multi),
              lower = list(cox$HR.CI_lower_uni, cox$HR.CI_lower_multi), 
              upper = list(cox$HR.CI_upper_uni, cox$HR.CI_upper_multi),
              ci_column = c(2, 6),
              ref_line = 1,
              xlim = c(0, max(cox$HR.CI_upper_uni)+0.5),
              ticks_at = c(0, 0.5, 1, 2),
              xlab = "Hazard Ratio",
              theme = forest_theme(core = list(bg_params=list(fill = c("white"))),
                                   arrow_label_just = "end",
                                   arrow_type = "closed"),
              nudge_y = 0.2)
  p <- insert_text(p, text = "Univariate",
                   part = "header",
                   row = 1,
                   col = 2:4,
                   gp = gpar(fontface = "bold"))
  p <- add_text(p, text = "Multivariate",
                part = "header",
                row = 1,
                col = 6:8,
                gp = gpar(fontface = "bold"))
  p <- add_border(p, 
                  part = "header", 
                  row = 1,
                  col = 2:4,
                  gp = gpar(lwd = 1))
  p <- add_border(p, 
                  part = "header", 
                  row = 1,
                  col = 6:8,
                  gp = gpar(lwd = 1))
  p <- add_border(p, 
                  part = "header", 
                  row = 2,
                  gp = gpar(lwd = 1))
  p <- edit_plot(p, col = c(3,4,7,8), 
                 which = "text",
                 hjust = unit(0.5, "npc"),
                 x = unit(0.5, "npc"))
  p <- insert_text(p, text = time_var,
                   part = "header",
                   row = 1,
                   col = 2:8,
                   gp = gpar(fontface = "bold"))
  p <- add_border(p, 
                  part = "header", 
                  row = 1,
                  col = 2:8,
                  gp = gpar(lwd = 1))
  return(p)
}



# Multivariate variables to choose:
#colnames(surv_vars_extended)
#colnames(surv_table_extended[1:28])
#table(surv_table_extended$status, surv_table_extended$CNS_disease)
#table(surv_table_extended$status, surv_table_extended$Chloroma)
#table(surv_table_extended$status, surv_table_extended$FAB_Category)
#table(surv_table_extended$status, surv_table_extended$SCT_in_1st_CR)
#table(surv_table_extended$status, surv_table_extended$FLT3.ITD_positive.)

coxcovars <- c(#"WBC_at_Diagnosis", # Numeric - 0% NAs
               "SCT_in_1st_CR",
               #"FLT3.ITD_positive.",
               #"CNS_disease",
               "MRD_at_end_of_course_1_perc", # Numeric - 10% NAs
               #"MRD_at_end_of_course_2_perc", # Numeric -  32% NAs (too many) --> quit
               #"Primary_Cytogenetic_Code", # 2% NAs but redundant with other variables --> quit
               "NPM_mutation",
               "CEBPA_mutation",
               #"WT1_mutation", # 65% NAs --> quit
               #"FLT3_PM", # 65% NAs --> quit
               #"c_Kit_Mutation_Exon_8", # 90% NAs --> quit
               #"c_Kit_Mutation_Exon_17", # 90% NAs --> quit
               "t_8_21_",
               "inv_16_"
               #colnames(clin[,25:42]) #cytogenetic abnormalities
)

predictors <- c("miR37_score", "AMLmiR36_score", "miR24_score", "miR4_score", "miR3_score",# "random_signature_score", 
                "hsa_mir_199a", "hsa_mir_335", "hsa_mir_196b", "hsa_mir_100", "hsa_mir_192", 
                "hsa_miR_139_5p", "hsa_mir_34b", "hsa_mir_195", "hsa_mir_381", "hsa_mir_370", "hsa_miR_193b_3p", 
                "hsa_miR_206", "hsa_mir_122", "hsa_mir_29a", "hsa_miR_375", "hsa_mir_155")


# Apply function:

# OS
p <- perform_uni_multi_coxph(
  data = surv_table_extended,
  time_var = "OS",
  status_var = "status",
  predictors = predictors,
  covars = coxcovars,
  save_table = F
)



while (!is.null(dev.list()))  dev.off()
pdf(paste0(save_dir, "/", data_subset, "_coxph_OS_withSignifCovars.pdf"), width = 12, height = 7, family="ArialMT")
p
while (!is.null(dev.list()))  dev.off()



# EFS
p <- perform_uni_multi_coxph(
  data = surv_table_extended,
  time_var = "EFS",
  status_var = "status",
  predictors = predictors,
  covars = coxcovars,
  save_table = F
)

while (!is.null(dev.list()))  dev.off()
pdf(paste0(save_dir, "/", data_subset, "_coxph_EFS_withSignifCovars.pdf"), width = 12, height = 7, family="ArialMT")
p
while (!is.null(dev.list()))  dev.off()



#----



#### Train Test patient selection (Validation set) ####

### Division of datasets in train and test by protocol:

#nrow(clin)
#table(clin$Protocol, exclude = NULL)
##trainy_set <- rownames(clin)[clin$Protocol=="AAML0531"|clin$Protocol=="AAML03P1"|clin$Protocol=="CCG2961"]
##test_set <- rownames(clin)[clin$Protocol=="AAML1031"]
#train_set <- rownames(clin)[clin$Protocol=="AAML1031"|clin$Protocol=="AAML03P1"]
#test_set <- rownames(clin)[clin$Protocol=="AAML0531"|clin$Protocol=="CCG2961"]
#length(train_set)
#length(test_set)
#train <- rownames(clin)%in%train_set
#test <- rownames(clin)%in%test_set


### Division of datasets in train and test by random X fold Cross Validation:

library(modelr)
set.seed(1994)
CV <- crossv_kfold(clin[Validation_set,], k=folds) # In Validation Set
CV
#rownames(clin[Validation_set,])[CV[[1]][[1]]$idx] # train CV1
#rownames(clin[Validation_set,])[CV[[1]][[2]]$idx] # train CV2
#rownames(clin[Validation_set,])[CV[[2]][[1]]$idx] # test CV1

#rownames(clin[Validation_set,])[CV[[1]][[1]]$idx] %in% rownames(clin[Validation_set,])[CV[[1]][[2]]$idx] # train CV sets are different but share patients
#rownames(clin[Validation_set,])[CV[[1]][[1]]$idx] %in% rownames(clin[Validation_set,])[CV[[2]][[2]]$idx] # train Fold 1 and test Fold 2 share some patients
#rownames(clin[Validation_set,])[CV[[2]][[1]]$idx] %in% rownames(clin[Validation_set,])[CV[[2]][[2]]$idx] # test CV sets are different and dont share patients.



#----

#### C index ####


#source("~/Desktop/Utility/Useful R codes.R")
#### Apply function
##allvars_to_test <- c("AMLmiR36_score", "miR3_score", "miR4_score", "miR24_score", sel_miRNAs)
#allvars_to_test <- miRNAs
#
#Cindex <- data.frame(var = allvars_to_test,
#                     OS_C_index_train = NA,
#                     OS_C_index_test = NA,
#                     EFS_C_index_train = NA,
#                     EFS_C_index_test = NA)
#for (i in 1:length(allvars_to_test)) {
#  os <- get_coxph_C(dataset = surv_table,
#                         var = allvars_to_test[i], # una o varias,
#                         status_var = "status",
#                         time_var = "OS",
#                         train = train,
#                         test = test)
#  efs <- get_coxph_C(dataset = surv_table,
#                                var = allvars_to_test[i], # una o varias,
#                                status_var = "status",
#                                time_var = "EFS",
#                                train = train,
#                                test = test)
#  Cindex[Cindex$var==allvars_to_test[i],]$OS_C_index_train <- os$C_index_train
#  Cindex[Cindex$var==allvars_to_test[i],]$OS_C_index_test <- os$C_index_test
#  Cindex[Cindex$var==allvars_to_test[i],]$EFS_C_index_train <- efs$C_index_train
#  Cindex[Cindex$var==allvars_to_test[i],]$EFS_C_index_test <- efs$C_index_test
#}
#Cindex
#
#write.table(Cindex, paste0(save_dir, "/", data_subset, "_", "Cindex.txt"))
#
#
#
#
#### Plot separated (bar plot)
#
#Cindex <- read.table(paste0(save_dir, "/", data_subset, "_", "Cindex.txt"))
#Cindex$var <- gsub("_", "-", Cindex$var)
#
### OS
#ci <- data.frame(var = Cindex$var,
#                 OS = Cindex$OS_C_index_test)
#ci <- ci[order(ci$OS, decreasing = T),] # order
#ci$var <- factor(ci$var, levels = ci$var) # para que queden ordenados en ggplot
#while (!is.null(dev.list()))  dev.off()
#pdf(paste0(save_dir, "/", data_subset, "_", "Cindex_OS.pdf"), width = 8, height = 3.5) #export as pdf
#ggplot(data = ci, aes(x = var, y = OS)) +
#  geom_bar(stat = "identity", position=position_dodge()) + 
#  theme_minimal() + 
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#  coord_cartesian(ylim = c(0.35, 0.65)) + 
#  ggtitle(data_subset) + 
#  theme(plot.title = element_text(hjust = 0.5)) +
#  xlab("") +
#  ylab("C Index")
#while (!is.null(dev.list()))  dev.off()
#
### EFS
#ci <- data.frame(var = Cindex$var,
#                 EFS = Cindex$EFS_C_index_test)
#ci <- ci[order(ci$EFS, decreasing = T),] # order
#ci$var <- factor(ci$var, levels = ci$var) # para que queden ordenados en ggplot
#while (!is.null(dev.list()))  dev.off()
#pdf(paste0(save_dir, "/", data_subset, "_", "Cindex_EFS.pdf"), width = 8, height = 3.5) #export as pdf
#ggplot(data = ci, aes(x = var, y = EFS)) +
#  geom_bar(stat = "identity", position=position_dodge()) + 
#  theme_minimal() + 
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#  coord_cartesian(ylim = c(0.35, 0.65)) + 
#  ggtitle(data_subset) + 
#  theme(plot.title = element_text(hjust = 0.5)) +
#  xlab("") +
#  ylab("C Index")
#while (!is.null(dev.list()))  dev.off()
#
#
#
### OS and EFS:
#ci <- data.frame(var = Cindex$var, # create new dataset of non-splitted variables
#                 OS = Cindex$OS_C_index_test,
#                 EFS = Cindex$EFS_C_index_test)
#ci <- ci[order(ci$OS, decreasing = T),] # order by OS
#library(magrittr)
#ci <- ci %>%
#  tidyr::pivot_longer(cols = c(OS, EFS),
#                      names_to = "survival_type",
#                      values_to = "value")
#ci$var <- factor(ci$var, levels = unique(ci$var)) # para que queden ordenados en ggplot
#while (!is.null(dev.list()))  dev.off()
#pdf(paste0(save_dir, "/", data_subset, "_", "Cindex.pdf")) #export as pdf
#ggplot(data = ci, aes(x = var, y = value, fill = survival_type)) +
#  geom_bar(stat = "identity", position=position_dodge()) + 
#  theme_minimal() + 
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#  scale_y_continuous(limits = c(0, 1)) +
#  ggtitle(data_subset) + 
#  theme(plot.title = element_text(hjust = 0.5)) +
#  xlab("") +
#  ylab("Concordance Index")
#while (!is.null(dev.list()))  dev.off()


#----

#### C index with CrossValidation ####


stopifnot(all(rownames(clin[Validation_set,])==rownames(surv_table_extended[Validation_set,]))) #check

source("~/Desktop/Utility/Useful R codes.R")
### Apply function
#allvars_to_test <- c("AMLmiR36_score", "miR3_score", "miR4_score", "miR24_score", sel_miRNAs)
allvars_to_test <- miRNAs

CindexByFold <- list()

for (fold in 1:nrow(CV)) {
  Cindex <- data.frame(var = allvars_to_test,
                     OS_C_index_train = NA,
                     OS_C_index_test = NA,
                     EFS_C_index_train = NA,
                     EFS_C_index_test = NA)
  for (i in 1:length(allvars_to_test)) {
    os <- get_coxph_C(dataset = surv_table_extended[Validation_set,], # Using Validation set
                      var = c(allvars_to_test[i], CIcovars), # una o varias,
                      status_var = "status",
                      time_var = "OS",
                      train = rownames(clin[Validation_set,])[CV[[1]][[fold]]$idx],
                      test = rownames(clin[Validation_set,])[CV[[2]][[fold]]$idx])
    efs <- get_coxph_C(dataset = surv_table_extended[Validation_set,],
                       var = c(allvars_to_test[i], CIcovars), # una o varias,
                       status_var = "status",
                       time_var = "EFS",
                       train = rownames(clin[Validation_set,])[CV[[1]][[fold]]$idx],
                       test = rownames(clin[Validation_set,])[CV[[2]][[fold]]$idx])
    Cindex[Cindex$var==allvars_to_test[i],]$OS_C_index_train <- os$C_index_train
    Cindex[Cindex$var==allvars_to_test[i],]$OS_C_index_test <- os$C_index_test
    Cindex[Cindex$var==allvars_to_test[i],]$EFS_C_index_train <- efs$C_index_train
    Cindex[Cindex$var==allvars_to_test[i],]$EFS_C_index_test <- efs$C_index_test
  }
  colnames(Cindex)[-1] <- paste0(colnames(Cindex)[-1], "_fold", fold)
  CindexByFold[[fold]] <- Cindex
}



write.table(CindexByFold, paste0(save_dir, "/", data_subset, "_", "Cindex_", fold, "foldCV.txt"), row.names = F)


# Calculate median CI
Cindex_tab <- read.csv(paste0(save_dir, "/", data_subset, "_", "Cindex_", fold, "foldCV.txt"), sep="") # 10 fold CV
Cindex_OS <- data.frame(row.names = Cindex_tab$var, "miRNA" = Cindex_tab$var,
                            Cindex_tab[,startsWith(colnames(Cindex_tab), "OS_C_index_test")])
Cindex_EFS <- data.frame(row.names = Cindex_tab$var, "miRNA" = Cindex_tab$var,
                         Cindex_tab[,startsWith(colnames(Cindex_tab), "EFS_C_index_test")])


Cindex_table <- data.frame(row.names = rownames(Cindex_OS),
                           "miRNA" = rownames(Cindex_OS))
library(matrixStats)
Cindex_table$CI_median_OS <- rowMedians(as.matrix(Cindex_OS[,-1]))
Cindex_table$CI_median_EFS <- rowMedians(as.matrix(Cindex_EFS[,-1]))
Cindex_table

write.table(Cindex_table, paste0(save_dir, "/", data_subset, "_", "Cindex_table_", fold, "foldCV.txt"), row.names = F)


#----




#### C index with CrossValidation (no scores, all signature miRNAs) ####


#stopifnot(all(rownames(clin[Validation_set,])==rownames(surv_table_extended[Validation_set,]))) #check
#
#source("~/Desktop/Utility/Useful R codes.R")
#### Apply function
##allvars_to_test <- c("AMLmiR36_score", "miR3_score", "miR4_score", "miR24_score", sel_miRNAs)
#allvars_to_test <- miRNAs
#
#
#
#CindexByFold <- list()
#for (fold in 1:nrow(CV)) {
#  Cindex <- data.frame(var = allvars_to_test,
#                       OS_C_index_train = NA,
#                       OS_C_index_test = NA,
#                       EFS_C_index_train = NA,
#                       EFS_C_index_test = NA)
#  for (i in 1:length(allvars_to_test)) {
#    if (grepl("score", allvars_to_test[i])) {
#      stopifnot(all(c(eval(as.name(gsub("_score", "_miRNAs", allvars_to_test[i])))) %in% colnames(surv_table_extended[Validation_set,])))
#      os <- get_coxph_C(dataset = surv_table_extended[Validation_set,], # Using Validation set
#                        var = c(eval(as.name(gsub("_score", "_miRNAs", allvars_to_test[i]))), CIcovars), # all miRNAs of the signature
#                        status_var = "status",
#                        time_var = "OS",
#                        train = rownames(clin[Validation_set,])[CV[[1]][[fold]]$idx],
#                        test = rownames(clin[Validation_set,])[CV[[2]][[fold]]$idx])
#      efs <- get_coxph_C(dataset = surv_table_extended[Validation_set,],
#                         var = c(eval(as.name(gsub("_score", "_miRNAs", allvars_to_test[i]))), CIcovars), # all miRNAs of the signature
#                         status_var = "status",
#                         time_var = "EFS",
#                         train = rownames(clin[Validation_set,])[CV[[1]][[fold]]$idx],
#                         test = rownames(clin[Validation_set,])[CV[[2]][[fold]]$idx])
#    } else {
#      os <- get_coxph_C(dataset = surv_table_extended[Validation_set,], # Using Validation set
#                        var = c(allvars_to_test[i], CIcovars), # una o varias,
#                        status_var = "status",
#                        time_var = "OS",
#                        train = rownames(clin[Validation_set,])[CV[[1]][[fold]]$idx],
#                        test = rownames(clin[Validation_set,])[CV[[2]][[fold]]$idx])
#      efs <- get_coxph_C(dataset = surv_table_extended[Validation_set,],
#                         var = c(allvars_to_test[i], CIcovars), # una o varias,
#                         status_var = "status",
#                         time_var = "EFS",
#                         train = rownames(clin[Validation_set,])[CV[[1]][[fold]]$idx],
#                         test = rownames(clin[Validation_set,])[CV[[2]][[fold]]$idx])
#    }
#    
#    Cindex[Cindex$var==allvars_to_test[i],]$OS_C_index_train <- os$C_index_train
#    Cindex[Cindex$var==allvars_to_test[i],]$OS_C_index_test <- os$C_index_test
#    Cindex[Cindex$var==allvars_to_test[i],]$EFS_C_index_train <- efs$C_index_train
#    Cindex[Cindex$var==allvars_to_test[i],]$EFS_C_index_test <- efs$C_index_test
#  }
#  colnames(Cindex)[-1] <- paste0(colnames(Cindex)[-1], "_fold", fold)
#  CindexByFold[[fold]] <- Cindex
#}
#
#
#write.table(CindexByFold, paste0(save_dir, "/", data_subset, "_", "Cindex_", fold, "foldCV_allSigMiRNAs.txt"), row.names = F)

#----

#### Time dependent ROC curve (all risks only) ####

sdata <- surv_table_extended[Validation_set,]

library(survivalROC)
#roc <- survivalROC(Stime = sdata$OS, status = sdata$status, marker = sdata$miR37_score, 
#            entry = NULL, 
#            predict.time = 100, 
#            cut.values = NULL, method = "KM", lambda = NULL, span = NULL, window = "symmetric")
#
#plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),
#     xlab=paste( "FP", "\n", "AUC = ",round(roc$AUC,3)),
#     ylab="TP",main="ROC, Method = KM \n Year = 1")
#abline(0,1)


### OS

auc_miR37 <- NULL
auc_AMLmiR36 <- NULL
auc_miR24 <- NULL
auc_miR4 <- NULL
auc_miR3 <- NULL
auc_random <- NULL
for (t in 30:max(sdata$OS)) {
  auc_miR37 <- c(auc_miR37, survivalROC.C(Stime = sdata$OS, status = sdata$status, marker = sdata$miR37_score, predict.time = t, span = 0.05)$AUC)
  auc_AMLmiR36 <- c(auc_AMLmiR36, survivalROC.C(Stime = sdata$OS, status = sdata$status, marker = sdata$AMLmiR36_score, predict.time = t, span = 0.05)$AUC)
  auc_miR24 <- c(auc_miR24, survivalROC.C(Stime = sdata$OS, status = sdata$status, marker = sdata$miR24_score, predict.time = t, span = 0.05)$AUC)
  auc_miR4 <- c(auc_miR4, survivalROC.C(Stime = sdata$OS, status = sdata$status, marker = sdata$miR4_score, predict.time = t, span = 0.05)$AUC)
  auc_miR3 <- c(auc_miR3, survivalROC.C(Stime = sdata$OS, status = sdata$status, marker = sdata$miR3_score, predict.time = t, span = 0.05)$AUC)
  auc_random <- c(auc_random, survivalROC.C(Stime = sdata$OS, status = sdata$status, marker = sdata$random_signature_score, predict.time = t, span = 0.05)$AUC)
}


while (!is.null(dev.list()))  dev.off()
#pdf(paste0(save_dir, "/All_risks_AUC_OS.pdf"), width = 5, height = 5.5)
pdf(paste0(save_dir, "/", data_subset, "_AUC_OS.pdf"), width = 5, height = 5.5)

plot(seq(min(sdata$OS)/365, max(sdata$OS)/365, length.out = length(auc_miR37)), auc_miR37, type = "l", xlab = "Time (years)", ylab = "AUC(t)", col = "green",
     xlim = c(0, 6), ylim = c(0.3, 0.8), main = "Overall Survival")
lines(seq(min(sdata$OS)/365, max(sdata$OS)/365, length.out = length(auc_AMLmiR36)), auc_AMLmiR36, col = "red")
lines(seq(min(sdata$OS)/365, max(sdata$OS)/365, length.out = length(auc_miR24)), auc_miR24, col = "blue")
lines(seq(min(sdata$OS)/365, max(sdata$OS)/365, length.out = length(auc_miR4)), auc_miR4, col = "orange")
lines(seq(min(sdata$OS)/365, max(sdata$OS)/365, length.out = length(auc_miR3)), auc_miR3, col = "yellow")
lines(seq(min(sdata$OS)/365, max(sdata$OS)/365, length.out = length(auc_random)), auc_random, col = "grey")
legend(3.65, 0.49, legend = c("miR37", "AMLmiR36", "miR24", "miR4", "miR3", "random"),
       col = c("green", "red", "blue", "orange", "yellow", "grey"), lty = c(1, 1, 1, 1, 1, 1), bty = "n")

while (!is.null(dev.list()))  dev.off()


### EFS

auc_miR37 <- NULL
auc_AMLmiR36 <- NULL
auc_miR24 <- NULL
auc_miR4 <- NULL
auc_miR3 <- NULL
auc_random <- NULL
for (t in 30:max(sdata$EFS)) {
  auc_miR37 <- c(auc_miR37, survivalROC.C(Stime = sdata$EFS, status = sdata$status, marker = sdata$miR37_score, predict.time = t, span = 0.05)$AUC)
  auc_AMLmiR36 <- c(auc_AMLmiR36, survivalROC.C(Stime = sdata$EFS, status = sdata$status, marker = sdata$AMLmiR36_score, predict.time = t, span = 0.05)$AUC)
  auc_miR24 <- c(auc_miR24, survivalROC.C(Stime = sdata$EFS, status = sdata$status, marker = sdata$miR24_score, predict.time = t, span = 0.05)$AUC)
  auc_miR4 <- c(auc_miR4, survivalROC.C(Stime = sdata$EFS, status = sdata$status, marker = sdata$miR4_score, predict.time = t, span = 0.05)$AUC)
  auc_miR3 <- c(auc_miR3, survivalROC.C(Stime = sdata$EFS, status = sdata$status, marker = sdata$miR3_score, predict.time = t, span = 0.05)$AUC)
  auc_random <- c(auc_random, survivalROC.C(Stime = sdata$EFS, status = sdata$status, marker = sdata$random_signature_score, predict.time = t, span = 0.05)$AUC)
}


while (!is.null(dev.list()))  dev.off()
#pdf(paste0(save_dir, "/All_risks_AUC_EFS.pdf"), width = 5, height = 5.5)
pdf(paste0(save_dir, "/", data_subset, "_AUC_EFS.pdf"), width = 5, height = 5.5)

plot(seq(min(sdata$EFS)/365, max(sdata$EFS)/365, length.out = length(auc_miR37)), auc_miR37, type = "l", xlab = "Time (years)", ylab = "AUC(t)", col = "green",
     xlim = c(0, 6), ylim = c(0.3, 0.8), main = "Event Free Survival")
lines(seq(min(sdata$EFS)/365, max(sdata$EFS)/365, length.out = length(auc_AMLmiR36)), auc_AMLmiR36, col = "red")
lines(seq(min(sdata$EFS)/365, max(sdata$EFS)/365, length.out = length(auc_miR24)), auc_miR24, col = "blue")
lines(seq(min(sdata$EFS)/365, max(sdata$EFS)/365, length.out = length(auc_miR4)), auc_miR4, col = "orange")
lines(seq(min(sdata$EFS)/365, max(sdata$EFS)/365, length.out = length(auc_miR3)), auc_miR3, col = "yellow")
lines(seq(min(sdata$EFS)/365, max(sdata$EFS)/365, length.out = length(auc_random)), auc_random, col = "grey")
legend(3.65, 0.49, legend = c("miR37", "AMLmiR36", "miR24", "miR4", "miR3", "random"),
       col = c("green", "red", "blue", "orange", "yellow", "grey"), lty = c(1, 1, 1, 1, 1, 1), bty = "n")

while (!is.null(dev.list()))  dev.off()




#----







#### Plot C index all risks together ####


#### Plot together (line plot)
#
#Cindex_All <- read.table("./results/All_Risks_Cindex.txt")
#Cindex_Low <- read.table("./results/Low_Risk_Cindex.txt")
#Cindex_Standard <- read.table("./results/Standard_Risk_Cindex.txt")
#Cindex_High <- read.table("./results/High_Risk_Cindex.txt")
#
#
#C_OS <- data.frame(miRNA = Cindex_All$var,
#                   All_Risks = Cindex_All$OS_C_index_test,
#                   Low_Risk = Cindex_Low$OS_C_index_test,
#                   Standard_Risk = Cindex_Standard$OS_C_index_test,
#                   High_Risk = Cindex_High$OS_C_index_test)
#C_EFS <- data.frame(miRNA = Cindex_All$var,
#                    All_Risks = Cindex_All$EFS_C_index_test,
#                    Low_Risk = Cindex_Low$EFS_C_index_test,
#                    Standard_Risk = Cindex_Standard$EFS_C_index_test,
#                    High_Risk = Cindex_High$EFS_C_index_test)
#
#C_OS <- C_OS[order(C_OS$All_Risks, decreasing = T),]
#C_OS[,2:5] <- round(C_OS[,2:5], digits = 3)
#C_EFS <- C_EFS[order(C_EFS$All_Risks, decreasing = T),]
#C_EFS[,2:5] <- round(C_EFS[,2:5], digits = 3)
#C_OS
#C_EFS
#
#
#
## Long format:
#library(data.table)
#C_OS.long <- as.data.frame(melt(setDT(C_OS), id.vars = "miRNA", variable.name = "Risk"))
#C_EFS.long <- as.data.frame(melt(setDT(C_EFS), id.vars = "miRNA", variable.name = "Risk"))
#
#
#C_OS.long$miRNA <- factor(C_OS.long$miRNA, levels = unique(C_OS.long$miRNA)) # para que queden ordenados en ggplot
#C_EFS.long$miRNA <- factor(C_EFS.long$miRNA, levels = unique(C_EFS.long$miRNA)) # para que queden ordenados en ggplot
#
#
#plot_Cindex <- function(data, title, xlab, ylab){
#  plot_C <- ggplot(data, aes(x=miRNA, y=value, group = Risk)) +
#    geom_line(aes(color=Risk,
#                  #linetype=Risk,
#                  alpha=Risk)) +
#    geom_point(aes(color=Risk,
#                   alpha=Risk)) +
#    scale_color_manual(values = c("#04050D", "#03A62C", "#F2CB05", "#F20707")) +
#    #scale_linetype_manual(values = c("solid", rep("dashed", 3))) +
#    scale_alpha_manual(values = c(1, rep(0.7, 3))) +
#    theme_minimal() + 
#    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#    scale_y_continuous(limits = c(0.25, 0.75)) +
#    ggtitle(title) + 
#    theme(plot.title = element_text(hjust = 0.5)) +
#    xlab(xlab) +
#    ylab(ylab)
#  return(plot_C)
#}
#
#while (!is.null(dev.list()))  dev.off()
#pdf(paste0(save_dir, "/All_data_Cindex.pdf"), width = 6, height = 5)
#plot_Cindex(C_OS.long, title = "OS", xlab = "", ylab = "C Index")
#plot_Cindex(C_EFS.long, title = "EFS", xlab = "", ylab = "C Index")
#while (!is.null(dev.list()))  dev.off()



#----

#### Plot C index (crossValidation) all risks together ####


#### Plot together (boxplot)
#
#Cindex_All <- read.csv(paste0("./results/All_Risks_Cindex_", "5", "foldCV.txt"), sep="") # 5 fold CV
#Cindex_Low <- read.csv(paste0("./results/Low_Risk_Cindex_", "3", "foldCV.txt"), sep="") # 3 fold CV
#Cindex_Standard <- read.csv(paste0("./results/Standard_Risk_Cindex_", "5", "foldCV.txt"), sep="") # 5 fold CV
#Cindex_High <- read.csv(paste0("./results/High_Risk_Cindex_", "3", "foldCV.txt"), sep="") # 3 fold CV
#
#colnames(Cindex_All)[-1] <- paste0("All_Risks_", colnames(Cindex_All)[-1])
#colnames(Cindex_Low)[-1] <- paste0("Low_Risk_", colnames(Cindex_Low)[-1])
#colnames(Cindex_Standard)[-1] <- paste0("Standard_Risk_", colnames(Cindex_Standard)[-1])
#colnames(Cindex_High)[-1] <- paste0("High_Risk_", colnames(Cindex_High)[-1])
#
#
#
#C_OS <- data.frame(miRNA = Cindex_All$var) # Dataset con todos los C index
#rownames(C_OS) <- C_OS$miRNA
#C_OS <- cbind(C_OS, Cindex_All[,startsWith(colnames(Cindex_All), "All_Risks_OS_C_index_test")])
#C_OS <- cbind(C_OS, Cindex_Low[,startsWith(colnames(Cindex_Low), "Low_Risk_OS_C_index_test")])
#C_OS <- cbind(C_OS, Cindex_Standard[,startsWith(colnames(Cindex_Standard), "Standard_Risk_OS_C_index_test")])
#C_OS <- cbind(C_OS, Cindex_High[,startsWith(colnames(Cindex_High), "High_Risk_OS_C_index_test")])
#
#C_EFS <- data.frame(miRNA = Cindex_All$var)
#rownames(C_EFS) <- C_EFS$miRNA
#C_EFS <- cbind(C_EFS, Cindex_All[,startsWith(colnames(Cindex_All), "All_Risks_EFS_C_index_test")])
#C_EFS <- cbind(C_EFS, Cindex_Low[,startsWith(colnames(Cindex_Low), "Low_Risk_EFS_C_index_test")])
#C_EFS <- cbind(C_EFS, Cindex_Standard[,startsWith(colnames(Cindex_Standard), "Standard_Risk_EFS_C_index_test")])
#C_EFS <- cbind(C_EFS, Cindex_High[,startsWith(colnames(Cindex_High), "High_Risk_EFS_C_index_test")])
#
#
## Ordenar los miRNAs por C index (primero signatures y luego miRNAs individuales):
#signatures <- rownames(C_OS)[endsWith(rownames(C_OS), "score")]
#indiv <- rownames(C_OS)[!endsWith(rownames(C_OS), "score")]
#
#signatures_ord <- signatures[order(rowMeans(C_OS[signatures, startsWith(colnames(C_OS), "All_Risks")]), decreasing = T)] # ordenar signatures por media de C index de los All risks
#indiv_ord <- indiv[order(rowMeans(C_OS[indiv, startsWith(colnames(C_OS), "All_Risks")]), decreasing = T)]
##C_OS <- C_OS[c(signatures_ord, indiv_ord),]
#
#
## Orden personalizado:
#C_OS <- C_OS[c("miR37_score", "AMLmiR36_score", "miR24_score", "miR4_score", "miR3_score", 
#               "hsa_mir_199a", "hsa_mir_335", "hsa_mir_196b", "hsa_mir_100", "hsa_mir_192", 
#               "hsa_miR_139_5p", "hsa_mir_34b", "hsa_mir_195", "hsa_mir_381", "hsa_mir_370", "hsa_miR_193b_3p", 
#               "hsa_miR_206", "hsa_mir_122", "hsa_mir_29a", "hsa_miR_375", "hsa_mir_155"),]
#
#
## Mismo orden para EFS:
#C_EFS <- C_EFS[rownames(C_OS),]
#
#
## Long format:
#library(data.table)
#C_OS.long <- as.data.frame(melt(setDT(C_OS), id.vars = "miRNA", variable.name = "Risk"))
#C_OS.long$fold <- gsub(".*fold", "", C_OS.long$Risk)
#C_OS.long$Risk <- gsub("_OS.*", "", C_OS.long$Risk)
#
#C_EFS.long <- as.data.frame(melt(setDT(C_EFS), id.vars = "miRNA", variable.name = "Risk"))
#C_EFS.long$fold <- gsub(".*fold", "", C_EFS.long$Risk)
#C_EFS.long$Risk <- gsub("_EFS.*", "", C_EFS.long$Risk)
#
## Fix order for ggplot
#C_OS.long$miRNA <- factor(C_OS.long$miRNA, levels = unique(C_OS.long$miRNA))
#C_OS.long$Risk <- factor(C_OS.long$Risk, levels = unique(C_OS.long$Risk))
#
#C_EFS.long$miRNA <- factor(C_EFS.long$miRNA, levels = unique(C_EFS.long$miRNA))
#C_EFS.long$Risk <- factor(C_EFS.long$Risk, levels = unique(C_EFS.long$Risk))
#
##utils::View(C_OS.long)
#
#
## Plot
#
#plot_Cindex <- function(data, title, xlab, ylab){
#  plot_C <- ggplot(data, aes(x = miRNA, y = value, fill = Risk)) +
#    geom_boxplot(outlier.size = 0.5) +
#    scale_fill_manual(values = c("#878787", "#03A62C", "#F2CB05", "#F20707")) +
#    theme_minimal() + 
#    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#    scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
#    ggtitle(title) + 
#    theme(plot.title = element_text(hjust = 0.5)) +
#    xlab(xlab) +
#    ylab(ylab)
#  return(plot_C)
#}
#
#
#
#while (!is.null(dev.list()))  dev.off()
##pdf(paste0(save_dir, "/All_data_Cindex_", folds, "foldCV.pdf"), width = 10, height = 7)
#pdf(paste0(save_dir, "/All_data_Cindex_XfoldCV.pdf"), width = 10, height = 7)
#plot_Cindex(C_OS.long, title = "OS", xlab = "", ylab = "C Index")
#plot_Cindex(C_EFS.long, title = "EFS", xlab = "", ylab = "C Index")
#while (!is.null(dev.list()))  dev.off()





#----

#### Plot C index (crossValidation) all risks separated ####


#### Plot together (boxplot)
#
#Cindex_All <- read.csv(paste0("./results/All_Risks_Cindex_", "5", "foldCV.txt"), sep="") # 5 fold CV
#Cindex_Low <- read.csv(paste0("./results/Low_Risk_Cindex_", "5", "foldCV.txt"), sep="") # 5 fold CV
#Cindex_Standard <- read.csv(paste0("./results/Standard_Risk_Cindex_", "5", "foldCV.txt"), sep="") # 5 fold CV
#Cindex_High <- read.csv(paste0("./results/High_Risk_Cindex_", "5", "foldCV.txt"), sep="") # 5 fold CV
#
## Separate OS and EFS
#Cindex_All_OS <- data.frame(row.names = Cindex_All$var, "miRNA" = Cindex_All$var,
#                            Cindex_All[,startsWith(colnames(Cindex_All), "OS_C_index_test")])
#Cindex_Low_OS <- data.frame(row.names = Cindex_Low$var, "miRNA" = Cindex_Low$var,
#                            Cindex_Low[,startsWith(colnames(Cindex_Low), "OS_C_index_test")])
#Cindex_Standard_OS <- data.frame(row.names = Cindex_Standard$var, "miRNA" = Cindex_Standard$var,
#                                 Cindex_Standard[,startsWith(colnames(Cindex_Standard), "OS_C_index_test")])
#Cindex_High_OS <- data.frame(row.names = Cindex_High$var, "miRNA" = Cindex_High$var,
#                             Cindex_High[,startsWith(colnames(Cindex_High), "OS_C_index_test")])
#Cindex_All_EFS <- data.frame(row.names = Cindex_All$var, "miRNA" = Cindex_All$var,
#                             Cindex_All[,startsWith(colnames(Cindex_All), "EFS_C_index_test")])
#Cindex_Low_EFS <- data.frame(row.names = Cindex_Low$var, "miRNA" = Cindex_Low$var,
#                             Cindex_Low[,startsWith(colnames(Cindex_Low), "EFS_C_index_test")])
#Cindex_Standard_EFS <- data.frame(row.names = Cindex_Standard$var, "miRNA" = Cindex_Standard$var,
#                                  Cindex_Standard[,startsWith(colnames(Cindex_Standard), "EFS_C_index_test")])
#Cindex_High_EFS <- data.frame(row.names = Cindex_High$var, "miRNA" = Cindex_High$var,
#                              Cindex_High[,startsWith(colnames(Cindex_High), "EFS_C_index_test")])
#
#ordered_miRNAs <- c("miR37_score", "AMLmiR36_score", "miR24_score", "miR4_score", "miR3_score", 
#                    "hsa_mir_199a", "hsa_mir_335", "hsa_mir_196b", "hsa_mir_100", "hsa_mir_192", 
#                    "hsa_miR_139_5p", "hsa_mir_34b", "hsa_mir_195", "hsa_mir_381", "hsa_mir_370", "hsa_miR_193b_3p", 
#                    "hsa_miR_206", "hsa_mir_122", "hsa_mir_29a", "hsa_miR_375", "hsa_mir_155")
#
#mRlist <- list(Cindex_All_OS = Cindex_All_OS,
#               Cindex_Low_OS = Cindex_Low_OS,
#               Cindex_Standard_OS = Cindex_Standard_OS,
#               Cindex_High_OS = Cindex_High_OS,
#               Cindex_All_EFS = Cindex_All_EFS,
#               Cindex_Low_EFS = Cindex_Low_EFS,
#               Cindex_Standard_EFS = Cindex_Standard_EFS,
#               Cindex_High_EFS = Cindex_High_EFS)
#
#p <- list()
#
#colors <- rep(c("#878787", "#03A62C", "#F2CB05", "#F20707"), 2)
#
#plot_Cindex_ind <- function(data, title, xlab, ylab, color){
#  plot_C <- ggplot(data, aes(x = miRNA, y = value)) +
#    geom_boxplot(outlier.size = 0.5, fill = color, alpha = 0.5) +
#    theme_minimal() + 
#    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#    scale_y_continuous(breaks = seq(0, 1, by = 0.1)#, limits = c(0, 1)
#    ) +
#    ggtitle(title) + 
#    theme(plot.title = element_text(hjust = 0.5)) +
#    xlab(xlab) +
#    ylab(ylab)
#  return(plot_C)
#}
#
#for (i in 1:length(mRlist)) {
#  mRlist[[i]] <- mRlist[[i]][ordered_miRNAs,] # order miRNAs
#  mRlist[[i]]$miRNA <- gsub("_", "-", mRlist[[i]]$miRNA)
#  library(data.table)
#  mRlist[[i]] <- as.data.frame(melt(setDT(mRlist[[i]]), id.vars = "miRNA", variable.name = "fold")) # long format
#  mRlist[[i]]$fold <- gsub(".*fold", "", mRlist[[i]]$fold)
#  mRlist[[i]]$miRNA <- factor(mRlist[[i]]$miRNA, levels = unique(mRlist[[i]]$miRNA)) # fix order for ggplot
#  # Plot
#  p[[i]] <- plot_Cindex_ind(mRlist[[i]], title = names(mRlist)[i], xlab = "", ylab = "C Index", color = colors[i])
#}
#
#while (!is.null(dev.list()))  dev.off()
#pdf(paste0(save_dir, "/All_data_Cindex_5foldCV_sep.pdf"), width = 16, height = 8, family="ArialMT"
#) #export as pdf
#library(gridExtra)
#args <- c(p, list(ncol = 4))
#do.call(grid.arrange, args)
#while (!is.null(dev.list()))  dev.off()



#----

#### Plot C index (crossValidation) all risks separated (10Fold) ####


### Plot together (boxplot)

Cindex_All <- read.csv(paste0(save_dir, "/All_Risks_Cindex_", "10", "foldCV.txt"), sep="") # 10 fold CV
Cindex_Low <- read.csv(paste0(save_dir, "/Low_Risk_Cindex_", "10", "foldCV.txt"), sep="") # 10 fold CV
Cindex_Standard <- read.csv(paste0(save_dir, "/Standard_Risk_Cindex_", "10", "foldCV.txt"), sep="") # 10 fold CV
Cindex_High <- read.csv(paste0(save_dir, "/High_Risk_Cindex_", "5", "foldCV.txt"), sep="") # 5 fold CV

# Separate OS and EFS
Cindex_All_OS <- data.frame(row.names = Cindex_All$var, "miRNA" = Cindex_All$var,
                            Cindex_All[,startsWith(colnames(Cindex_All), "OS_C_index_test")])
Cindex_Low_OS <- data.frame(row.names = Cindex_Low$var, "miRNA" = Cindex_Low$var,
                            Cindex_Low[,startsWith(colnames(Cindex_Low), "OS_C_index_test")])
Cindex_Standard_OS <- data.frame(row.names = Cindex_Standard$var, "miRNA" = Cindex_Standard$var,
                                 Cindex_Standard[,startsWith(colnames(Cindex_Standard), "OS_C_index_test")])
Cindex_High_OS <- data.frame(row.names = Cindex_High$var, "miRNA" = Cindex_High$var,
                             Cindex_High[,startsWith(colnames(Cindex_High), "OS_C_index_test")])
Cindex_All_EFS <- data.frame(row.names = Cindex_All$var, "miRNA" = Cindex_All$var,
                             Cindex_All[,startsWith(colnames(Cindex_All), "EFS_C_index_test")])
Cindex_Low_EFS <- data.frame(row.names = Cindex_Low$var, "miRNA" = Cindex_Low$var,
                             Cindex_Low[,startsWith(colnames(Cindex_Low), "EFS_C_index_test")])
Cindex_Standard_EFS <- data.frame(row.names = Cindex_Standard$var, "miRNA" = Cindex_Standard$var,
                                  Cindex_Standard[,startsWith(colnames(Cindex_Standard), "EFS_C_index_test")])
Cindex_High_EFS <- data.frame(row.names = Cindex_High$var, "miRNA" = Cindex_High$var,
                              Cindex_High[,startsWith(colnames(Cindex_High), "EFS_C_index_test")])

ordered_miRNAs <- c("miR37_score", "AMLmiR36_score", "miR24_score", "miR4_score", "miR3_score", "random_signature_score", 
                    "hsa_mir_199a", "hsa_mir_335", "hsa_mir_196b", "hsa_mir_100", "hsa_mir_192", 
                    "hsa_miR_139_5p", "hsa_mir_34b", "hsa_mir_195", "hsa_mir_381", "hsa_mir_370", "hsa_miR_193b_3p", 
                    "hsa_miR_206", "hsa_mir_122", "hsa_mir_29a", "hsa_miR_375", "hsa_mir_155")

mRlist <- list(Cindex_All_OS = Cindex_All_OS,
               Cindex_Low_OS = Cindex_Low_OS,
               Cindex_Standard_OS = Cindex_Standard_OS,
               Cindex_High_OS = Cindex_High_OS,
               Cindex_All_EFS = Cindex_All_EFS,
               Cindex_Low_EFS = Cindex_Low_EFS,
               Cindex_Standard_EFS = Cindex_Standard_EFS,
               Cindex_High_EFS = Cindex_High_EFS)

p <- list()

colors <- rep(c("#878787", "#03A62C", "#F2CB05", "#F20707"), 2)

plot_Cindex_ind <- function(data, title, xlab, ylab, color){
  plot_C <- ggplot(data, aes(x = miRNA, y = value)) +
    geom_boxplot(outlier.size = 0.5, fill = color, alpha = 0.5) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1)#, limits = c(0, 1)
    ) +
    ggtitle(title) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(xlab) +
    ylab(ylab)
  return(plot_C)
}

for (i in 1:length(mRlist)) {
  mRlist[[i]] <- mRlist[[i]][ordered_miRNAs,] # order miRNAs
  mRlist[[i]]$miRNA <- gsub("_", "-", mRlist[[i]]$miRNA)
  library(data.table)
  mRlist[[i]] <- as.data.frame(melt(setDT(mRlist[[i]]), id.vars = "miRNA", variable.name = "fold")) # long format
  mRlist[[i]]$fold <- gsub(".*fold", "", mRlist[[i]]$fold)
  mRlist[[i]]$miRNA <- factor(mRlist[[i]]$miRNA, levels = unique(mRlist[[i]]$miRNA)) # fix order for ggplot
  # Plot
  p[[i]] <- plot_Cindex_ind(mRlist[[i]], title = names(mRlist)[i], xlab = "", ylab = "C Index", color = colors[i])
}

while (!is.null(dev.list()))  dev.off()
pdf(paste0(save_dir, "/All_data_Cindex_10foldCV_sep.pdf"), width = 16, height = 8, family="ArialMT"
) #export as pdf
library(gridExtra)
args <- c(p, list(ncol = 4))
do.call(grid.arrange, args)
while (!is.null(dev.list()))  dev.off()



#----




#### Compare Cindex from signature score vs miRNAs of signature ####

#vars <- "miR37_score"
#vars <- miR37_miRNAs
#vars <- "miR24_score"
#vars <- miR24_miRNAs
#vars <- "AMLmiR36_score"
#vars <- AMLmiR36_miRNAs # surprisingly good (0.70)
#vars <- "miR3_score"
#vars <- miR3_miRNAs
#vars <- "miR4_score"
#vars <- miR4_miRNAs
#
#
#CindexByFold <- data.frame("CItrain" = NA,
#                           "CItest" = NA)
#for (fold in 1:nrow(CV)) {
#  ci <- get_coxph_C(dataset = surv_table_extended[Validation_set,], # Using Validation set
#                    var = c(vars, CIcovars), # una o varias,
#                    status_var = "status",
#                    time_var = "OS",
#                    train = rownames(clin[Validation_set,])[CV[[1]][[fold]]$idx],
#                    test = rownames(clin[Validation_set,])[CV[[2]][[fold]]$idx])
#  CindexByFold[fold,"CItrain"] <- ci$C_index_train
#  CindexByFold[fold,"CItest"] <- ci$C_index_test
#}
#colMeans(CindexByFold)


#----

#### Plot C index (crossValidation) all risks together (all signature miRNAs) ####


#### Plot together (boxplot)
#
#Cindex_All <- read.csv(paste0("./results/All_Risks_Cindex_", "5", "foldCV_allSigMiRNAs.txt"), sep="") # 5 fold CV
#Cindex_Low <- read.csv(paste0("./results/Low_Risk_Cindex_", "5", "foldCV_allSigMiRNAs.txt"), sep="") # 5 fold CV
#Cindex_Standard <- read.csv(paste0("./results/Standard_Risk_Cindex_", "5", "foldCV_allSigMiRNAs.txt"), sep="") # 5 fold CV
#Cindex_High <- read.csv(paste0("./results/High_Risk_Cindex_", "5", "foldCV_allSigMiRNAs.txt"), sep="") # 5 fold CV
#
#colnames(Cindex_All)[-1] <- paste0("All_Risks_", colnames(Cindex_All)[-1])
#colnames(Cindex_Low)[-1] <- paste0("Low_Risk_", colnames(Cindex_Low)[-1])
#colnames(Cindex_Standard)[-1] <- paste0("Standard_Risk_", colnames(Cindex_Standard)[-1])
#colnames(Cindex_High)[-1] <- paste0("High_Risk_", colnames(Cindex_High)[-1])
#
#
#
#C_OS <- data.frame(miRNA = Cindex_All$var) # Dataset con todos los C index
#rownames(C_OS) <- C_OS$miRNA
#C_OS <- cbind(C_OS, Cindex_All[,startsWith(colnames(Cindex_All), "All_Risks_OS_C_index_test")])
#C_OS <- cbind(C_OS, Cindex_Low[,startsWith(colnames(Cindex_Low), "Low_Risk_OS_C_index_test")])
#C_OS <- cbind(C_OS, Cindex_Standard[,startsWith(colnames(Cindex_Standard), "Standard_Risk_OS_C_index_test")])
#C_OS <- cbind(C_OS, Cindex_High[,startsWith(colnames(Cindex_High), "High_Risk_OS_C_index_test")])
#
#C_EFS <- data.frame(miRNA = Cindex_All$var)
#rownames(C_EFS) <- C_EFS$miRNA
#C_EFS <- cbind(C_EFS, Cindex_All[,startsWith(colnames(Cindex_All), "All_Risks_EFS_C_index_test")])
#C_EFS <- cbind(C_EFS, Cindex_Low[,startsWith(colnames(Cindex_Low), "Low_Risk_EFS_C_index_test")])
#C_EFS <- cbind(C_EFS, Cindex_Standard[,startsWith(colnames(Cindex_Standard), "Standard_Risk_EFS_C_index_test")])
#C_EFS <- cbind(C_EFS, Cindex_High[,startsWith(colnames(Cindex_High), "High_Risk_EFS_C_index_test")])
#
#
## Ordenar los miRNAs por C index (primero signatures y luego miRNAs individuales):
#signatures <- rownames(C_OS)[endsWith(rownames(C_OS), "score")]
#indiv <- rownames(C_OS)[!endsWith(rownames(C_OS), "score")]
#
#signatures_ord <- signatures[order(rowMeans(C_OS[signatures, startsWith(colnames(C_OS), "All_Risks")]), decreasing = T)] # ordenar signatures por media de C index de los All risks
#indiv_ord <- indiv[order(rowMeans(C_OS[indiv, startsWith(colnames(C_OS), "All_Risks")]), decreasing = T)]
##C_OS <- C_OS[c(signatures_ord, indiv_ord),]
#
#
## Orden personalizado:
#C_OS <- C_OS[c("miR37_score", "AMLmiR36_score", "miR24_score", "miR4_score", "miR3_score", 
#               "hsa_mir_199a", "hsa_mir_335", "hsa_mir_196b", "hsa_mir_100", "hsa_mir_192", 
#               "hsa_miR_139_5p", "hsa_mir_34b", "hsa_mir_195", "hsa_mir_381", "hsa_mir_370", "hsa_miR_193b_3p", 
#               "hsa_miR_206", "hsa_mir_122", "hsa_mir_29a", "hsa_miR_375", "hsa_mir_155"),]
#
#
## Mismo orden para EFS:
#C_EFS <- C_EFS[rownames(C_OS),]
#
#
## Long format:
#library(data.table)
#C_OS.long <- as.data.frame(melt(setDT(C_OS), id.vars = "miRNA", variable.name = "Risk"))
#C_OS.long$fold <- gsub(".*fold", "", C_OS.long$Risk)
#C_OS.long$Risk <- gsub("_OS.*", "", C_OS.long$Risk)
#
#C_EFS.long <- as.data.frame(melt(setDT(C_EFS), id.vars = "miRNA", variable.name = "Risk"))
#C_EFS.long$fold <- gsub(".*fold", "", C_EFS.long$Risk)
#C_EFS.long$Risk <- gsub("_EFS.*", "", C_EFS.long$Risk)
#
## Fix order for ggplot
#C_OS.long$miRNA <- factor(C_OS.long$miRNA, levels = unique(C_OS.long$miRNA))
#C_OS.long$Risk <- factor(C_OS.long$Risk, levels = unique(C_OS.long$Risk))
#
#C_EFS.long$miRNA <- factor(C_EFS.long$miRNA, levels = unique(C_EFS.long$miRNA))
#C_EFS.long$Risk <- factor(C_EFS.long$Risk, levels = unique(C_EFS.long$Risk))
#
##utils::View(C_OS.long)
#
#
## Plot
#
#plot_Cindex <- function(data, title, xlab, ylab){
#  plot_C <- ggplot(data, aes(x = miRNA, y = value, fill = Risk)) +
#    geom_boxplot(outlier.size = 0.5) +
#    scale_fill_manual(values = c("#878787", "#03A62C", "#F2CB05", "#F20707")) +
#    theme_minimal() + 
#    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#    scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
#    ggtitle(title) + 
#    theme(plot.title = element_text(hjust = 0.5)) +
#    xlab(xlab) +
#    ylab(ylab)
#  return(plot_C)
#}
#
#
#
#while (!is.null(dev.list()))  dev.off()
##pdf(paste0(save_dir, "/All_data_Cindex_", folds, "foldCV.pdf"), width = 10, height = 7)
#pdf(paste0(save_dir, "/All_data_Cindex_XfoldCV_allSigMiRNAs.pdf"), width = 10, height = 7)
#plot_Cindex(C_OS.long, title = "OS", xlab = "", ylab = "C Index")
#plot_Cindex(C_EFS.long, title = "EFS", xlab = "", ylab = "C Index")
#while (!is.null(dev.list()))  dev.off()





#----







#### Functional enrichment analysis of new signature ####


setwd('~/Desktop/Review code/miRNA_meta_analysis')

Enrichment.TAM2_result <- read.delim("./results/Functional analysis TAM2/Results table.txt")


### Functions

functions.TAM2 <- Enrichment.TAM2_result[Enrichment.TAM2_result$Category=="Function",]
functions.TAM2 <- functions.TAM2[order(functions.TAM2$Term, decreasing = F),]
functions.TAM2 <- functions.TAM2[order(functions.TAM2$P.value, decreasing = F),]
functions.TAM2_top <- functions.TAM2[1:10,]


## Barplot
functions.TAM2_top$TermCount <- paste0(functions.TAM2_top$Term, " (", functions.TAM2_top$Count, ")")
functions.TAM2_top$logpval <- -log10(functions.TAM2_top$P.value)
functions.TAM2_top$TermCount <- factor(functions.TAM2_top$TermCount, levels = unique(functions.TAM2_top$TermCount))
library(ggplot2)
plotbarF <- ggplot(functions.TAM2_top, aes(x=logpval, y=TermCount)) +
  geom_bar(stat = "identity") +
  scale_y_discrete(limits = rev(levels(functions.TAM2_top$TermCount))) +
  xlab("-log10(P.value)") +
  ylab("Function") +
  theme_minimal()

# Save
while (!is.null(dev.list()))  dev.off()
pdf("./results/Functional analysis TAM2/Bar_functions.pdf", width = 4, height = 4, family="ArialMT")
plotbarF
while (!is.null(dev.list()))  dev.off()


## Network plot

#remotes::install_github("jmw86069/multienrichjam")
library(multienrichjam)
functions.TAM2_top.enr <- enrichDF2enrichResult(enrichDF = functions.TAM2_top, keyColname = "Term", 
                                                geneColname = "miRNA", pvalueColname = "P.value", 
                                                pvalueCutoff = 0.05, descriptionColname = "Term")
#BiocManager::install("enrichplot")
library(enrichplot)
plotnetF <- cnetplot(functions.TAM2_top.enr, 
                  showCategory = 10,
                  node_label = "category",
                  cex_category = 1,
                  cex_gene = 1,
                  cex_label_category = 1,
                  cex_label_gene = 0.7,
                  color_category = "#E5C494",
                  color_gene = "#B3B3B3"#,
                  #hilight.params = list(category = c("Hematopoiesis"), alpha_hilight = 1, alpha_no_hilight = 0.3)
                  )

# Save
while (!is.null(dev.list()))  dev.off()
pdf("./results/Functional analysis TAM2/Net_functions.pdf", width = 7, height = 7, family="ArialMT")
plotnetF
while (!is.null(dev.list()))  dev.off()



### Diseases

diseases.TAM2 <- Enrichment.TAM2_result[Enrichment.TAM2_result$Category=="Disease",]
diseases.TAM2 <- diseases.TAM2[order(diseases.TAM2$Term, decreasing = F),]
diseases.TAM2 <- diseases.TAM2[order(diseases.TAM2$P.value, decreasing = F),]
diseases.TAM2_top <- diseases.TAM2[1:10,]

## Barplot
diseases.TAM2_top$TermCount <- paste0(diseases.TAM2_top$Term, " (", diseases.TAM2_top$Count, ")")
diseases.TAM2_top$logpval <- -log10(diseases.TAM2_top$P.value)
diseases.TAM2_top$TermCount <- factor(diseases.TAM2_top$TermCount, levels = unique(diseases.TAM2_top$TermCount))
library(ggplot2)
plotbarF <- ggplot(diseases.TAM2_top, aes(x=logpval, y=TermCount)) +
  geom_bar(stat = "identity") +
  scale_y_discrete(limits = rev(levels(diseases.TAM2_top$TermCount))) +
  xlab("-log10(P.value)") +
  ylab("Disease") +
  theme_minimal()

# Save
while (!is.null(dev.list()))  dev.off()
pdf("./results/Functional analysis TAM2/Bar_diseases.pdf", width = 4.8, height = 4, family="ArialMT")
plotbarF
while (!is.null(dev.list()))  dev.off()


## Network plot
library(multienrichjam)
diseases.TAM2_top.enr <- enrichDF2enrichResult(enrichDF = diseases.TAM2_top, keyColname = "Term", 
                                               geneColname = "miRNA", pvalueColname = "P.value", 
                                               pvalueCutoff = 0.05, descriptionColname = "Term")
library(enrichplot)
plotnetD <- cnetplot(diseases.TAM2_top.enr, 
                  showCategory = 10,
                  node_label = "category",
                  cex_category = 1,
                  cex_gene = 1,
                  cex_label_category = 1,
                  cex_label_gene = 0.7,
                  color_category = "#E5C494",
                  color_gene = "#B3B3B3"#,
                  #hilight.params = list(category = c("Leukemia, Myeloid, Acute"), alpha_hilight = 1, alpha_no_hilight = 0.3)
                  )

# Save
while (!is.null(dev.list()))  dev.off()
pdf("./results/Functional analysis TAM2/Net_diseases.pdf", width = 7, height = 7, family="ArialMT")
plotnetD
while (!is.null(dev.list()))  dev.off()




## Barplot together


#library(dplyr)
#FD.TAM2_top <- bind_rows(functions.TAM2_top, diseases.TAM2_top)
#
#library(ggplot2)
#plotbarF <- ggplot(diseases.TAM2_top, aes(x=logpval, y=TermCount)) +
#  geom_bar(stat = "identity") +
#  scale_y_discrete(limits = rev(levels(diseases.TAM2_top$TermCount))) +
#  xlab("-log10(P.value)") +
#  ylab("Disease") +
#  theme_minimal()
#
## Save
#while (!is.null(dev.list()))  dev.off()
#pdf("./results/Functional analysis TAM2/Bar_diseases.pdf", width = 4, height = 4, family="ArialMT")
#plotbarF
#while (!is.null(dev.list()))  dev.off()






#----







#### Common miRNAs between signatures: ####


miR37 <- read.delim(paste0(save_dir, "/signature_new_miRNA.txt"))
AMLmiR36 <- read.delim(paste0(save_dir, "/signature_AMLmiR36.txt"))
miR24 <- read.delim(paste0(save_dir, "/signature_miR24.txt"))
miR4 <- read.delim(paste0(save_dir, "/signature_miR4.txt"))
miR3 <- read.delim(paste0(save_dir, "/signature_miR3.txt"))
miR37_miRNAs <- miR37$miRNA
AMLmiR36_miRNAs <- AMLmiR36$miRNA
miR24_miRNAs <- miR24$miRNA
miR4_miRNAs <- miR4$miRNA
miR3_miRNAs <- miR3$miRNA
sel_miRNAs


# Get common immature mirnas:
sets <- list(miR37_miRNAs = miR37_miRNAs,
             AMLmiR36_miRNAs = AMLmiR36_miRNAs,
             miR24_miRNAs = miR24_miRNAs,
             miR4_miRNAs = miR4_miRNAs,
             miR3_miRNAs = miR3_miRNAs,
             sel_miRNAs = sel_miRNAs)

for (i in 1:length(sets)) {
  sets[[i]] <- gsub("_3p", "", sets[[i]])
  sets[[i]] <- gsub("_5p", "", sets[[i]])
  sets[[i]] <- gsub("miR", "mir", sets[[i]])
}


# Venn diagram plot

library(gplots)
venn(sets[c(1:3)])
library(VennDiagram)
v <- venn.diagram(sets[1:5],
             fill = c("red", "green","blue","yellow","orange"),
             alpha = c(0.5, 0.5, 0.5, 0.5, 0.5), cat.cex = 1.5, cex=1.5,
             filename=NULL)
grid.newpage()
grid.draw(v)



miR37_miRNAs
AMLmiR36_miRNAs
miR24_miRNAs
miR4_miRNAs
miR3_miRNAs
sel_miRNAs



intersect(sets$miR24_miRNAs, sets$sel_miRNAs)
intersect(miR37_miRNAs, AMLmiR36_miRNAs)
intersect(miR24_miRNAs, AMLmiR36_miRNAs)
intersect(miR24_miRNAs, miR37_miRNAs)
intersect(intersect(miR24_miRNAs, AMLmiR36_miRNAs), miR37_miRNAs)
intersect(sets$miR24_miRNAs, sets$sel_miRNAs)
intersect(AMLmiR36_miRNAs, miR4_miRNAs)
intersect(AMLmiR36_miRNAs, sel_miRNAs)
intersect(sets$sel_miRNAs, sets$miR37_miRNAs)
intersect(sets$sel_miRNAs, sets$AMLmiR36_miRNAs)
intersect(sets$sel_miRNAs, sets$miR24_miRNAs)
intersect(sets$sel_miRNAs, sets$miR3_miRNAs)
intersect(sets$sel_miRNAs, sets$miR4_miRNAs)




#----















#### --------------------------------Review fix:------------------------------------- ####

#### SigQC ####



### Load expr dataset
#expr_sigQC <- expr_notnorm # save unnormalized counts
#
#library(edgeR)
#d <- DGEList(counts=expr_sigQC, genes=rownames(expr_sigQC))
#d <- calcNormFactors(d
#                     , method="TMM" # If we want TMM normalization
#)
#expr_tmm <- cpm(d, normalized.lib.sizes = T, log = F)
##expr_tmm <- cpm(d, normalized.lib.sizes = T, log = F)
#
#head(expr_tmm[1:5,1:5])
#expr_sigQC <- expr_tmm


expr_sigQC <- expr




### Load new GEO dataset

omicsdataLocation <- "~/Desktop/Data"
library(GEOquery)
GSE97135 <- getGEO("GSE97135", destdir = paste0(omicsdataLocation,"/OmicsData/GEO"))
GSE97135_metadata <- GSE97135[["GSE97135_series_matrix.txt.gz"]]@phenoData@data
GSE97135_expr <- data.frame(GSE97135[["GSE97135_series_matrix.txt.gz"]]@assayData[["exprs"]])
GSE97135_geneAnnot <- GSE97135[["GSE97135_series_matrix.txt.gz"]]@featureData@data
all(rownames(GSE97135_expr)==rownames(GSE97135_geneAnnot))
rownames(GSE97135_expr) <- GSE97135_geneAnnot$miRNA_ID
nrow(GSE97135_expr)
ncol(GSE97135_expr)
rownames(GSE97135_expr) <- gsub("-", "_", rownames(GSE97135_expr))

## Normalize: (NO)
#library(limma)
#plotDensities(GSE97135_expr[,1:5])
#plotDensities(expr_sigQC[,1:5])
#head(GSE97135_expr[1:5,1:5])
#
#library(edgeR)
#d <- DGEList(counts=GSE97135_expr, genes=rownames(GSE97135_expr))
#d <- calcNormFactors(d, method="TMM") # If we want TMM normalization
##GSE97135_expr_tmm <- cpm(d, normalized.lib.sizes = T, log = T)
#GSE97135_expr_tmm <- cpm(d, normalized.lib.sizes = T, log = F)
#head(GSE97135_expr_tmm[1:5,1:5])
#GSE97135_expr <- GSE97135_expr_tmm

nrow(GSE97135_expr)
nrow(expr_sigQC)
length(intersect(rownames(GSE97135_expr), rownames(expr_sigQC)))
#check if all signature miRNAs are in both expression datasets:
table(c(miR37$miRNA, miR24$miRNA, AMLmiR36$miRNA, miR3$miRNA, miR4$miRNA, random_signature$miRNA) %in% intersect(rownames(GSE97135_expr), rownames(expr_sigQC)))
table(c(miR3$miRNA, miR4$miRNA) %in% intersect(rownames(GSE97135_expr), rownames(expr_sigQC)))



# Input of sigQC
sig_list <- list()
sig_list[["miR37_signature"]] <- as.matrix(miR37$miRNA)
#sig_list[["AMLmiR36_signature"]] <- as.matrix(AMLmiR36$miRNA)
#sig_list[["miR24_signature"]] <- as.matrix(miR24$miRNA)
#sig_list[["miR3_signature"]] <- c("hsa_miR_146b_5p", "hsa_miR_181c_5p", "hsa_miR_4786_5p") # no immature in our datasets so pick the dominant mature
#sig_list[["miR4_signature"]] <- c("hsa_miR_509_3p", "hsa_miR_542_3p", "hsa_miR_3667_5p", "hsa_miR_146a_5p") # no immature in our datasets so pick the dominant mature
#sig_list[["random_signature"]] <- as.matrix(random_signature$miRNA)
expr_list <- list()
expr_list[["TARGET_AML_Validation_dataset"]] <- expr_sigQC[,Validation_set]
expr_list[["GSE97135_dataset"]] <- as.matrix(GSE97135_expr)



# Run sigQC
library(sigQC)
set.seed(1994)
make_all_plots(sig_list, expr_list, out_dir = paste0(save_dir, "/", "sigQC_results_onlymiR37"), showResults = F)
dev.off()

rowMeans(expr_sigQC[miR37$miRNA,Validation_set])
rowMedians(expr_sigQC[miR37$miRNA,Validation_set])
rowMeans(expr_sigQC[random_signature$miRNA,Validation_set])
rowMedians(expr_sigQC[random_signature$miRNA,Validation_set])
colMeans(expr_sigQC[random_signature$miRNA,Validation_set])
matrixStats::colMedians(expr_sigQC[random_signature$miRNA,Validation_set])
summary(colMeans(expr_sigQC[random_signature$miRNA,Validation_set]))
summary(matrixStats::colMedians(expr_sigQC[random_signature$miRNA,Validation_set]))
table(matrixStats::colMedians(expr_sigQC[random_signature$miRNA,Validation_set]))
table(matrixStats::colMedians(expr_notnorm[random_signature$miRNA,Validation_set]))


#library(GEOquery)
#GSE13159 <- getGEO("GSE13159", destdir = paste0(omicsdataLocation,"/OmicsData/GEO"))
#GSE13159_metadata <- GSE13159[["GSE13159_series_matrix.txt.gz"]]@phenoData@data
#GSE13159_expr <- data.frame(GSE13159[["GSE13159_series_matrix.txt.gz"]]@assayData[["exprs"]])
#GSE13159_geneAnnot <- GSE13159[["GSE13159_series_matrix.txt.gz"]]@featureData@data


#----



#### Alternatives to LASSO ####


#### Develop signatures ####

stopifnot(all(all_miRNAs_mature%in%colnames(surv_table_extended[Discovery_set,])))

library(edgeR)
keep <- filterByExpr(expr_notnorm_mature, min.count = 10)
expr_notnorm_mature_filt <- expr_notnorm_mature[keep,]
dim(expr_notnorm_mature)
dim(expr_notnorm_mature_filt) # We go from 2280 to 270 miRNAs

data <- surv_table_extended[Discovery_set, c("OS", "status", rownames(expr_notnorm_mature_filt))]
biomarkers <- rownames(expr_notnorm_mature_filt)

library(biospear)
set.seed(918273645)
fit_enet <- BMsel(data = data, # use Discovery set and filtered miRNAs
                     x = biomarkers,
                     y = c("OS", "status"),
                     #z = cvrts, # no covariates
                     inter = FALSE,
                     method = "enet",
                     dfmax = 50,
                     folds = 5
)
fit_gboost <- BMsel(data = data, # use Discovery set and filtered miRNAs
                     x = biomarkers,
                     y = c("OS", "status"),
                     #z = cvrts, # no covariates
                     inter = FALSE,
                     method = "gboost",
                     dfmax = 50,
                     folds = 5
)


length(fit_enet$enet)
length(fit_gboost$gboost)


##enet:
enet <- data.frame(miRNA = names(fit_enet$enet),
                   Coefficient = unname(fit_enet$enet))
enet_miRNAs <- enet$miRNA
enet_coef <- enet$Coefficient
# Getting Score for each patient: sum of all coef * miRNA:
enet_scores <- as.data.frame(colSums(expr[enet$miRNA,]*enet$Coefficient))
colnames(enet_scores) <- "enet_score"

##gboost:
gboost <- data.frame(miRNA = names(fit_gboost$gboost),
                     Coefficient = unname(fit_gboost$gboost))
gboost_miRNAs <- gboost$miRNA
gboost_coef <- gboost$Coefficient
# Getting Score for each patient: sum of all coef * miRNA:
gboost_scores <- as.data.frame(colSums(expr[gboost$miRNA,]*gboost$Coefficient))
colnames(gboost_scores) <- "gboost_score"



##lasso miR37
miR37 <- read.delim(paste0(save_dir, "/signature_new_miRNA.txt"))
miR37_miRNAs <- miR37$miRNA
miR37_coef <- miR37$Coefficient

# Getting Score for each patient: sum of all coef * miRNA:
miR37_scores <- as.data.frame(colSums(expr[miR37$miRNA,]*miR37$Coefficient))
colnames(miR37_scores) <- "miR37_score"



## Random
random_signature <- read.delim(paste0(save_dir, "/random_signature.txt"))
random_signature_miRNAs <- random_signature$miRNA
random_signature_coef <- random_signature$Coefficient



# Getting Score for each patient: sum of all coef * miRNA:
random_signature_scores <- as.data.frame(colSums(expr[random_signature$miRNA,]*random_signature$Coefficient))
colnames(random_signature_scores) <- "random_signature_score"


### Adding scores to surv_table:
stopifnot(all(rownames(surv_table_extended)==rownames(miR37_scores)))
stopifnot(all(rownames(surv_table_extended)==rownames(enet_scores)))
stopifnot(all(rownames(surv_table_extended)==rownames(gboost_scores)))
stopifnot(all(rownames(surv_table_extended)==rownames(random_signature_scores)))
stopifnot(all(rownames(surv_table)==rownames(miR37_scores)))
stopifnot(all(rownames(surv_table)==rownames(enet_scores)))
stopifnot(all(rownames(surv_table)==rownames(gboost_scores)))
stopifnot(all(rownames(surv_table)==rownames(random_signature_scores)))
surv_table$miR37_score <- miR37_scores$miR37_score
surv_table$enet_score <- enet_scores$enet_score
surv_table$gboost_score <- gboost_scores$gboost_score
surv_table$random_signature_score <- random_signature_scores$random_signature_score
surv_table_extended$miR37_score <- miR37_scores$miR37_score
surv_table_extended$enet_score <- enet_scores$enet_score
surv_table_extended$gboost_score <- gboost_scores$gboost_score
surv_table_extended$random_signature_score <- random_signature_scores$random_signature_score

#----


#### Train Test patient selection (Validation set) ####

## Division of datasets in train and test by random X fold Cross Validation:

library(modelr)
set.seed(1994)
CV <- crossv_kfold(clin[Validation_set,], k=folds) # In Validation Set
CV

#----


#### C index with CrossValidation ####

stopifnot(all(rownames(clin[Validation_set,])==rownames(surv_table_extended[Validation_set,]))) #check

source("~/Desktop/Utility/Useful R codes.R")
### Apply function
allvars_to_test <- c("miR37_score", "enet_score", "gboost_score", "random_signature_score")

CindexByFold <- list()

for (fold in 1:nrow(CV)) {
  Cindex <- data.frame(var = allvars_to_test,
                       OS_C_index_train = NA,
                       OS_C_index_test = NA,
                       EFS_C_index_train = NA,
                       EFS_C_index_test = NA)
  for (i in 1:length(allvars_to_test)) {
    os <- get_coxph_C(dataset = surv_table_extended[Validation_set,], # Using Validation set
                      var = c(allvars_to_test[i], CIcovars), # una o varias,
                      status_var = "status",
                      time_var = "OS",
                      train = rownames(clin[Validation_set,])[CV[[1]][[fold]]$idx],
                      test = rownames(clin[Validation_set,])[CV[[2]][[fold]]$idx])
    efs <- get_coxph_C(dataset = surv_table_extended[Validation_set,],
                       var = c(allvars_to_test[i], CIcovars), # una o varias,
                       status_var = "status",
                       time_var = "EFS",
                       train = rownames(clin[Validation_set,])[CV[[1]][[fold]]$idx],
                       test = rownames(clin[Validation_set,])[CV[[2]][[fold]]$idx])
    Cindex[Cindex$var==allvars_to_test[i],]$OS_C_index_train <- os$C_index_train
    Cindex[Cindex$var==allvars_to_test[i],]$OS_C_index_test <- os$C_index_test
    Cindex[Cindex$var==allvars_to_test[i],]$EFS_C_index_train <- efs$C_index_train
    Cindex[Cindex$var==allvars_to_test[i],]$EFS_C_index_test <- efs$C_index_test
  }
  colnames(Cindex)[-1] <- paste0(colnames(Cindex)[-1], "_fold", fold)
  CindexByFold[[fold]] <- Cindex
}


write.table(CindexByFold, paste0(save_dir, "/", "Model_Comparison_", data_subset, "_", "Cindex_", fold, "foldCV.txt"), row.names = F)


# Calculate median CI
Cindex_tab <- read.csv(paste0(save_dir, "/", "Model_Comparison_", data_subset, "_", "Cindex_", fold, "foldCV.txt"), sep="") # 10 fold CV
Cindex_OS <- data.frame(row.names = Cindex_tab$var, "miRNA" = Cindex_tab$var,
                        Cindex_tab[,startsWith(colnames(Cindex_tab), "OS_C_index_test")])
Cindex_EFS <- data.frame(row.names = Cindex_tab$var, "miRNA" = Cindex_tab$var,
                         Cindex_tab[,startsWith(colnames(Cindex_tab), "EFS_C_index_test")])


Cindex_table <- data.frame(row.names = rownames(Cindex_OS),
                           "miRNA" = rownames(Cindex_OS))
library(matrixStats)
Cindex_table$CI_median_OS <- rowMedians(as.matrix(Cindex_OS[,-1]))
Cindex_table$CI_median_EFS <- rowMedians(as.matrix(Cindex_EFS[,-1]))
Cindex_table

write.table(Cindex_table, paste0(save_dir, "/", "Model_Comparison_",  data_subset, "_", "Cindex_table_", fold, "foldCV.txt"), row.names = F)


#----


#### Time dependent ROC curve (all risks only) ####

sdata <- surv_table_extended[Validation_set,]
library(survivalROC)

## OS

auc_miR37 <- NULL
auc_enet <- NULL
auc_gboost <- NULL
auc_random <- NULL
for (t in 30:max(sdata$OS)) {
  auc_miR37 <- c(auc_miR37, survivalROC.C(Stime = sdata$OS, status = sdata$status, marker = sdata$miR37_score, predict.time = t, span = 0.05)$AUC)
  auc_enet <- c(auc_enet, survivalROC.C(Stime = sdata$OS, status = sdata$status, marker = sdata$enet_score, predict.time = t, span = 0.05)$AUC)
  auc_gboost <- c(auc_gboost, survivalROC.C(Stime = sdata$OS, status = sdata$status, marker = sdata$gboost_score, predict.time = t, span = 0.05)$AUC)
  auc_random <- c(auc_random, survivalROC.C(Stime = sdata$OS, status = sdata$status, marker = sdata$random_signature_score, predict.time = t, span = 0.05)$AUC)
}


while (!is.null(dev.list()))  dev.off()
pdf(paste0(save_dir, "/Model_Comparison_All_risks_AUC_OS.pdf"), width = 5, height = 5.5)

plot(seq(min(sdata$OS)/365, max(sdata$OS)/365, length.out = length(auc_miR37)), auc_miR37, type = "l", xlab = "Time (years)", ylab = "AUC(t)", col = "green",
     xlim = c(0, 6), ylim = c(0.3, 0.8), main = "Overall Survival")
lines(seq(min(sdata$OS)/365, max(sdata$OS)/365, length.out = length(auc_enet)), auc_enet, col = "blue")
lines(seq(min(sdata$OS)/365, max(sdata$OS)/365, length.out = length(auc_gboost)), auc_gboost, col = "orange")
lines(seq(min(sdata$OS)/365, max(sdata$OS)/365, length.out = length(auc_random)), auc_random, col = "grey")
legend(3.65, 0.49, legend = c("miR37_score", "enet_score", "gboost_score", "random_signature_score"),
       col = c("green", "blue", "orange", "grey"), lty = c(1, 1, 1, 1), bty = "n")

while (!is.null(dev.list()))  dev.off()


## EFS

auc_miR37 <- NULL
auc_enet <- NULL
auc_gboost <- NULL
auc_random <- NULL
for (t in 30:max(sdata$EFS)) {
  auc_miR37 <- c(auc_miR37, survivalROC.C(Stime = sdata$EFS, status = sdata$status, marker = sdata$miR37_score, predict.time = t, span = 0.05)$AUC)
  auc_enet <- c(auc_enet, survivalROC.C(Stime = sdata$EFS, status = sdata$status, marker = sdata$enet_score, predict.time = t, span = 0.05)$AUC)
  auc_gboost <- c(auc_gboost, survivalROC.C(Stime = sdata$EFS, status = sdata$status, marker = sdata$gboost_score, predict.time = t, span = 0.05)$AUC)
  auc_random <- c(auc_random, survivalROC.C(Stime = sdata$EFS, status = sdata$status, marker = sdata$random_signature_score, predict.time = t, span = 0.05)$AUC)
}


while (!is.null(dev.list()))  dev.off()
pdf(paste0(save_dir, "/Model_Comparison_All_risks_AUC_EFS.pdf"), width = 5, height = 5.5)

plot(seq(min(sdata$EFS)/365, max(sdata$EFS)/365, length.out = length(auc_miR37)), auc_miR37, type = "l", xlab = "Time (years)", ylab = "AUC(t)", col = "green",
     xlim = c(0, 6), ylim = c(0.3, 0.8), main = "Event Free Survival")
lines(seq(min(sdata$EFS)/365, max(sdata$EFS)/365, length.out = length(auc_enet)), auc_enet, col = "blue")
lines(seq(min(sdata$EFS)/365, max(sdata$EFS)/365, length.out = length(auc_gboost)), auc_gboost, col = "orange")
lines(seq(min(sdata$EFS)/365, max(sdata$EFS)/365, length.out = length(auc_random)), auc_random, col = "grey")
legend(3.65, 0.49, legend = c("miR37_score", "enet_score", "gboost_score", "random_signature_score"),
       col = c("green","blue", "orange", "grey"), lty = c(1, 1, 1, 1), bty = "n")

while (!is.null(dev.list()))  dev.off()




#----


#### Plot C index (crossValidation) all risks separated (10Fold) ####


### Plot together (boxplot)

Cindex_All <- read.csv(paste0(save_dir, "/Model_Comparison_All_Risks_Cindex_", "10", "foldCV.txt"), sep="") # 10 fold CV


# Separate OS and EFS
Cindex_All_OS <- data.frame(row.names = Cindex_All$var, "miRNA" = Cindex_All$var,
                            Cindex_All[,startsWith(colnames(Cindex_All), "OS_C_index_test")])

Cindex_All_EFS <- data.frame(row.names = Cindex_All$var, "miRNA" = Cindex_All$var,
                             Cindex_All[,startsWith(colnames(Cindex_All), "EFS_C_index_test")])


ordered_miRNAs <- c("miR37_score", "enet_score", "gboost_score", "random_signature_score")

mRlist <- list(Cindex_All_OS = Cindex_All_OS,
               Cindex_All_EFS = Cindex_All_EFS)

p <- list()

colors <- rep("#878787", 2)

plot_Cindex_ind <- function(data, title, xlab, ylab, color){
  plot_C <- ggplot(data, aes(x = miRNA, y = value)) +
    geom_boxplot(outlier.size = 0.5, fill = color, alpha = 0.5) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1)#, limits = c(0, 1)
    ) +
    ggtitle(title) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(xlab) +
    ylab(ylab)
  return(plot_C)
}

for (i in 1:length(mRlist)) {
  mRlist[[i]] <- mRlist[[i]][ordered_miRNAs,] # order miRNAs
  mRlist[[i]]$miRNA <- gsub("_", "-", mRlist[[i]]$miRNA)
  library(data.table)
  mRlist[[i]] <- as.data.frame(melt(setDT(mRlist[[i]]), id.vars = "miRNA", variable.name = "fold")) # long format
  mRlist[[i]]$fold <- gsub(".*fold", "", mRlist[[i]]$fold)
  mRlist[[i]]$miRNA <- factor(mRlist[[i]]$miRNA, levels = unique(mRlist[[i]]$miRNA)) # fix order for ggplot
  # Plot
  p[[i]] <- plot_Cindex_ind(mRlist[[i]], title = names(mRlist)[i], xlab = "", ylab = "C Index", color = colors[i])
}

while (!is.null(dev.list()))  dev.off()
pdf(paste0(save_dir, "/Model_Comparison_All_data_Cindex_10foldCV_sep.pdf"), width = 16, height = 8, family="ArialMT"
) #export as pdf
library(gridExtra)
args <- c(p, list(ncol = 4))
do.call(grid.arrange, args)
while (!is.null(dev.list()))  dev.off()



#----


#----







#### extra code ####

univ_formulas <- sapply(miRNAs,
                        function(x) as.formula(paste('Surv(Overall_Survival_Time_in_Days, status)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = surv_table, id = rownames(surv_table))})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res_continuous <- as.data.frame(res)



univ_models <- lapply( univ_formulas, function(x){coxph(x, data = surv_table_categorical, id = rownames(surv_table_categorical))})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res_categorical <- as.data.frame(res)



res_continuous
res_categorical







#----






















#### Otros datos ####


setwd('~/Desktop/Review code/miRNA_meta_analysis')
data_dir <- "./data"
dir.create(data_dir)
GEOdata_dir <- paste0(data_dir, "/", "GEOdata")
dir.create(GEOdata_dir)



library(GEOquery)

GSE97135 <- getGEO(GEO = "GSE97135",
                   destdir = GEOdata_dir,
                   AnnotGPL = T)

GSE35320 <- getGEO(GEO = "GSE35320", # 2 different GPL
                   destdir = GEOdata_dir,
                   AnnotGPL = T)

GSE98697 <- getGEO(GEO = "GSE98697", # lncRNAs
                   destdir = GEOdata_dir,
                   AnnotGPL = T)

GSE94066 <- getGEO(GEO = "GSE94066", # RT-PCR of pAML
                   destdir = GEOdata_dir,
                   AnnotGPL = T) 


utils::View(pData(GSE97135[[1]]))
exprs(GSE97135[[1]])

utils::View(pData(GSE35320[[1]]))
exprs(GSE35320[[1]])

utils::View(pData(GSE35320[[2]]))
exprs(GSE35320[[2]])

utils::View(pData(GSE98697[[1]]))
exprs(GSE98697[[1]])

utils::View(pData(GSE94066[[1]]))
exprs(GSE94066[[1]])


#----


#### Create our own model (protocol lasso) ####

# Get all variables together

surv_table_extended[1:10,1:10]

source("~/Desktop/Dream_challenge_FINRISK/DREAM-FINRISK-master/src/jlb_m1/requirements.r")
source("~/Desktop/Dream_challenge_FINRISK/DREAM-FINRISK-master/src/utils/importPseq.r")
source("~/Desktop/Dream_challenge_FINRISK/DREAM-FINRISK-master/src/utils/prepro_functions.r")
source("~/Desktop/Dream_challenge_FINRISK/DREAM-FINRISK-master/src/utils/co-abundances.r")
source("~/Desktop/Dream_challenge_FINRISK/DREAM-FINRISK-master/src/utils/get_scores.r")
source("~/Desktop/Dream_challenge_FINRISK/DREAM-FINRISK-master/src/jlb_m1/fit_model.r")
source("~/Desktop/Dream_challenge_FINRISK/DREAM-FINRISK-master/src/jlb_m1/predRes_helper.r")
source("~/Desktop/Dream_challenge_FINRISK/DREAM-FINRISK-master/src/preprocessing/preprocessing.r")



surv_table_extended$WBC_at_Diagnosis
summary(surv_table_extended$WBC_at_Diagnosis)
summary(surv_table_extended$Age_at_Diagnosis_in_Days)
summary(surv_table_extended$SCT_in_1st_CR)
table(surv_table_extended$SCT_in_1st_CR, exclude = NULL)
summary(surv_table_extended$FLT3.ITD_positive.)
table(surv_table_extended$FLT3.ITD_positive., exclude = NULL)
summary(surv_table_extended$CNS_disease)
table(surv_table_extended$CNS_disease, exclude = NULL)

cvrts <- c(#"Risk_group"
  "Age_at_Diagnosis_in_Days"
  #, "WBC_at_Diagnosis" # Has 1 NA
  #, "SCT_in_1st_CR"
  , "FLT3.ITD_positive."
  #, "CNS_disease"
)

library(biospear)
fit_lasso <- BMsel(data = surv_table_extended[train,],
                   x = all_miRNAs_mature,
                   y = c("OS", "status"),
                   #z = cvrts,
                   inter = FALSE,
                   method = "lasso",
                   dfmax = 50
)
length(fit_lasso$lasso)
names(fit_lasso$lasso)
unname(fit_lasso$lasso)

#prediction <- predRes(res = fit,
#                       method = "lasso",
#                       traindata = surv_table_extended[train,],
#                       newdata = surv_table_extended[test,],
#                       int.cv = FALSE,
#                       int.cv.nfold = 5,
#                       time = c(2, 365, 365*2, 365*3, 365*4, 365*5, 365*6, 365*7, 365*8),
#                       trace = TRUE,
#                       ncores = 20)
#predRes()
#
#prediction[["time = 1"]][["External validation"]]
#prediction[["time = 2"]][["External validation"]]
#prediction[["time = 100"]][["External validation"]]
#prediction[["time = 1000"]][["External validation"]]
#prediction[["time = 2000"]][["External validation"]]
#prediction[["time = 2730"]][["External validation"]]
#prediction[["time = 3000"]][["External validation"]]



# Build final model cox with all train data

f <- as.formula(paste0("Surv(OS, status) ~", paste(names(fit$lasso),collapse='+')))
model <- coxph(f, surv_table_extended[train,])




# Predict risk

risk = function(model, newdata, time) {
  as.numeric(1-summary(
    survfit(model, newdata = newdata, se.fit = F, conf.int = F), 
    times = time, extend = TRUE)$surv)
}
pred <- risk(model, surv_table_extended[test,], time = 365*5)


# Final scores

scores <- data.frame(
  SampleID = rownames(surv_table_extended[test,]),
  Score = pred
)
print(head(scores))


library(Hmisc)
cindex <- Hmisc::rcorr.cens(probs, Surv(OS, status), outx = FALSE)


source("~/Desktop/Utility/Useful R codes.R")
get_coxph_C(dataset = surv_table_extended,
            var = names(fit$lasso),
            time_var = "OS",
            status_var = "status",
            train = train,
            test = test)



# Method 2
get_coxph_C_2 <- function(dataset, var, status_var, time_var, train, test){ # funcion para obtener el concordance index de los modelos coxph a partir de los datos de training y test
  data <- data.frame(time = dataset[,time_var],
                     status = dataset[,status_var],
                     dataset[,var])
  data_train <- data[train,]
  data_test <- data[test,]
  
  
  collect.harrelc <- data.frame()
  
  # Train:
  model <- coxph(Surv(time, status)~., data = data_train#, x = TRUE
  )
  scores <- as.data.frame(predict(model, 
                                  newdata=data_test,     
                                  se.type="expected"))
  library(tibble)
  library(magrittr)
  scores = scores %>%
    set_colnames(c("Score")) %>%
    rownames_to_column(var = "SampleID")
  
  # Align the user provided scores with the true event times
  
  # Calculate Harrell's C statistics
  HarrellC <- Hmisc::rcorr.cens(exp(-scores$Score), Surv(data_test$time, data_test$status), outx=FALSE)
  collect.harrelc <- rbind(collect.harrelc, HarrellC["C Index"])
  colnames(collect.harrelc) <- "HarrellC"
  res <- collect.harrelc
  return(res)
}

get_coxph_C_2(dataset = surv_table_extended,
              var = names(fit$lasso),
              time_var = "OS",
              status_var = "status",
              train = train,
              test = test)






# Feature Selection with biospear


# Fit model

method = "lasso"
cvrts <- c("Age", "BodyMassIndex",
           "SystolicBP", "NonHDLcholesterol",
           "Sex", "disbiosis")
models <- fit_biospear(data = train,
                       biomarkers = setdiff(colnames(train),
                                            c(cvrts, colnames(y_train))),
                       surv = c("Event_time", "Event"),
                       cvrts = cvrts,
                       inter = FALSE,
                       methods = method)
print(models)


# Build final model

# Build cox with all train data

selected_betas <- models$lasso
train <- train[, match(
  c("Event", "Event_time", names(selected_betas)),
  colnames(train))]
surv <- "Surv(Event_time, Event) ~"
f <- as.formula(paste0(surv, paste(names(selected_betas),collapse='+')))
model <- coxph(f, train)


# Predict risk

# Risk Prediction

risk = function(model, newdata, time) {
  as.numeric(1-summary(
    survfit(model, newdata = newdata, se.fit = F, conf.int = F), 
    times = time, extend = TRUE)$surv)
}
pred <- risk(model, test, time = 15)


# Final scores

# Save scores
scores <- data.frame(
  SampleID = rownames(test),
  Score = pred
)
print(head(scores))





#----



#### mir24 Cox Proportional Hazards model: ####

# Multivriate variables to choose:
colnames(surv_vars_extended)
colnames(surv_table_extended[1:28])
table(surv_table_extended$status, surv_table_extended$CNS_disease)
table(surv_table_extended$status, surv_table_extended$Chloroma)
table(surv_table_extended$status, surv_table_extended$FAB_Category)
table(surv_table_extended$status, surv_table_extended$SCT_in_1st_CR)
table(surv_table_extended$status, surv_table_extended$FLT3.ITD_positive.)


### OS

surv <- Surv(time = surv_vars_extended$OS, event = surv_vars_extended$status) # OS

any(colSums(surv_table_extended[,miR24_miRNAs])==0)

# Multivariate cox analysis:
cox_reg_model_multi <- coxph(formula = surv ~ hsa_mir_20b+
                               hsa_mir_223+
                               hsa_mir_193a+
                               hsa_mir_24_1+
                               hsa_mir_128_1+
                               hsa_mir_17+
                               hsa_mir_199b+
                               hsa_mir_181c+
                               hsa_mir_181a_1+
                               hsa_mir_181b_1+
                               hsa_mir_21+
                               hsa_mir_222+
                               hsa_mir_331+
                               hsa_mir_373+
                               hsa_mir_708+
                               hsa_mir_34b+
                               hsa_mir_195+
                               hsa_mir_151a+
                               hsa_mir_30b+
                               hsa_mir_22+
                               hsa_let_7g+
                               hsa_let_7i+
                               hsa_mir_1290+
                               hsa_mir_9_1
                             ,data = surv_table_extended) 

summary(cox_reg_model_multi)


# Coef exlanation: Positive coef means that an increase in the predictor variable is associated with an increased hazard of the event
# Positive coef: more expression more risk. 
# Negative coef: less expression more risk.



### C index 

#surv_table_extended <- surv_table_extended[!is.na(surv_table_extended$Risk_group),]

source("~/Desktop/Utility/Useful R codes.R")
os <- get_coxph_C(dataset = surv_table_extended,
                  var = c(miR24_miRNAs
                          #, "Risk_group" # remove when selecting a risk group
                          #, Age_at_Diagnosis_in_Days 
                          #, Gender              # not related to outcome
                          #, "WBC_at_Diagnosis" 
                          #, "SCT_in_1st_CR"
                          #, FAB_Category        # FAB outdated
                          #, "MLL"                 # mutations not included: correlated to risk group
                          #, "FLT3.ITD_positive."
                          #, NPM_mutation
                          #, CEBPA_mutation
                          #, Primary_Cytogenetic_Code
                          #, Bone_marrow_leukemic_blast_percentage_...
                          #, Peripheral_blasts_...
                          #   , "CNS_disease"
                          #, Chloroma, # una o varias,
                  ),
                  status_var = "status",
                  time_var = "OS",
                  train = train,
                  test = test)
os

os2 <- get_coxph_C(dataset = surv_table_extended,
                   var = "miR24_score",
                   status_var = "status",
                   time_var = "OS",
                   train = train,
                   test = test)
os2






get_coxph_C <- function(dataset,
                        var, 
                        status_var, 
                        time_var, 
                        train, 
                        test){ 
  
  #' Fits a Cox model and calculate the Harrell's C Index
  #' 
  #' @description 
  #' Fits a Cox model for survival data and calculate the Harrell's Concordance Index for the train and test.
  #' C index of test data indicates the predictive power of model.
  #' https://cran.r-project.org/web/packages/SurvMetrics/vignettes/SurvMetrics-vignette.html
  #' 
  #' @param dataset Data frame with patients in rows all variables in columns: time, status and variables (clinical, gene expression, etc)
  #' @param var One or more variables to use in the Cox model.
  #' @param status_var Status variable: Dead or alive for survival data
  #' @param time_var Time variable: Survival time. Can be Overall Survival of Event Free Survival, for example.
  #' @param train Training vector. Can be a TRUE/FALSE vector or a vector of patient names.
  #' @param test Test vector. Can be a TRUE/FALSE vector or a vector of patient names.
  #' 
  
  data <- data.frame(time = dataset[,time_var],
                     status = dataset[,status_var],
                     dataset[,var])
  data_train <- data[train,]
  data_test <- data[test,]
  
  # Train:
  cox_reg_model <- coxph(Surv(time, status)~., data = data_train, x = TRUE)
  C_index_train <- concordance(cox_reg_model)$concordance # or summary(cox_reg_model)$concordance
  dis_time = sort(unique(data_train$time[data_train$status==1]))
  # Test:
  library(pec)
  mat_cox = predictSurvProb(cox_reg_model, data_test, dis_time)
  #View(surv_table_extended[test,])
  med_index = median(1:length(dis_time))
  surv_obj = Surv(time = data_test$time, event = data_test$status) 
  #C index for Cox:
  library(SurvMetrics)
  C_index_test = Cindex(surv_obj, predicted = mat_cox[, med_index])
  res <- data.frame(C_index_train = C_index_train,
                    C_index_test = C_index_test)
  return(res)
}











#https://www.reneshbedre.com/blog/survival-analysis.html#:~:text=and%20drug_2%20treatments.-,Cox's%20proportional%20hazards%20(CPH)%20model,not%20consider%20additional%20independent%20variables.


## Aplicar modelo:
##miRNAs <- sel_miRNAs
#surv <- Surv(time = surv_vars_extended$OS, event = surv_vars_extended$status) # OS
###cox_reg_model <- coxph(formula = surv ~ surv_table_extended[,miRNAs[1]] + surv_table_extended[,miRNAs[2]] + Risk_group + Gender + Race + Ethnicity + Age_at_Diagnosis_in_Days, data = surv_table_extended)  
#cox_reg_model <- coxph(formula = surv ~ AMLmiR36_score + Risk_group + WBC_at_Diagnosis_status + SCT_in_1st_CR + FLT3.ITD_positive. + Gender + Age_at_Diagnosis_in_Days,
#                       data = surv_table_extended)  
#summary(cox_reg_model)
### Podemos acceder a los datos asi:
#library(broom)
#results_cox <- as.data.frame(tidy(cox_reg_model))
### Hazard ratios (= exp(coef)):
#results_cox$HR <- exp(results_cox$estimate)
#results_cox
#
#ggforest(cox_reg_model, data=surv_vars_extended)
#
## ATENCION: ANALISIS ES MULTIVARIANTE (UNAS VARIABLES AFECTAN A LAS OTRAS)





### Predictive evaluation of model:
# We use Concordance Index (CI or C index):
# https://cran.r-project.org/web/packages/SurvMetrics/vignettes/SurvMetrics-vignette.html
source("~/Desktop/Utility/Useful R codes.R")



### Examples:

## Example AMLmiR36_score OS
#get_coxph_C(dataset = surv_table,
#                       var = c("AMLmiR36_score"), # una o varias,
#                       status_var = "status",
#                       time_var = "OS",
#                       train = train,
#                       test = test)
#
## Example AMLmiR36_score OS
#get_coxph_C(dataset = surv_table_extended,
#                       var = c("AMLmiR36_score", "Gender", "Age_at_Diagnosis_in_Days"), # una o varias,
#                       status_var = "status",
#                       time_var = "OS",
#                       train = train,
#                       test = test)
#
## Example AMLmiR36_score OS
#ds <- surv_table_extended
#ds <- ds[!is.na(ds$Risk_group),]
#dim(ds)
#s3 <- sample(rownames(ds), round(length(rownames(ds))/3)) # random
#s3 <- 
#train_3 <- !rownames(ds)%in%s3
#test_3 <- rownames(ds)%in%s3
#table(ds$Risk_group, exclude = NULL)
#
#get_coxph_C(dataset = ds
#                       ,
#                       var = c("AMLmiR36_score"
#                               #, "miR4_score"
#                               #, "Gender"
#                               #, "Age_at_Diagnosis_in_Days"
#                               , "Risk_group"
#                               ), # una o varias,
#                       status_var = "status",
#                       time_var = "OS",
#                       train = train_3,
#                       test = test_3)
#
#allvars_to_test <- c("AMLmiR36_score", "miR3_score", "miR4_score", sel_miRNAs)
#
#ds <- surv_table_extended
#ds <- ds[!is.na(ds$Risk_group),]
#dim(ds)
##sample_test_names <- sample(rownames(ds), round(length(rownames(ds))/3))
#sample_test_names <- rownames(clin)[clin$Protocol=="AAML1031"]
#train_3 <- !rownames(ds)%in%sample_test_names
#test_3 <- rownames(ds)%in%sample_test_names
#table(ds$Risk_group, exclude = NULL)
#
#Ci <- data.frame(miRNA= allvars_to_test,
#                 Cindex=NA)
#
#for (i in 1:length(allvars_to_test)) {
#  Ci_res <- get_coxph_C(dataset = ds,
#                                   var = c(allvars_to_test[i]
#                                           #, "miR4_score"
#                                           #, "Gender"
#                                           #, "Age_at_Diagnosis_in_Days"
#                                           , "Risk_group"
#                                   ),
#                                   status_var = "status",
#                                   time_var = "OS",
#                                   train = train_3,
#                                   test = test_3)
#  Ci$Cindex[i] <- Ci_res$C_index_test
#}
#Ci <- Ci[order(Ci$Cindex, decreasing = T),]
#Ci





### Apply function
allvars_to_test <- c("AMLmiR36_score", "miR3_score", "miR4_score", sel_miRNAs)
Cindex <- data.frame(var = allvars_to_test,
                     OS_C_index_train = NA,
                     OS_C_index_test = NA,
                     EFS_C_index_train = NA,
                     EFS_C_index_test = NA)
for (i in 1:length(allvars_to_test)) {
  os <- get_coxph_C(dataset = surv_table,
                    var = allvars_to_test[i], # una o varias,
                    status_var = "status",
                    time_var = "OS",
                    train = train,
                    test = test)
  efs <- get_coxph_C(dataset = surv_table,
                     var = allvars_to_test[i], # una o varias,
                     status_var = "status",
                     time_var = "EFS",
                     train = train,
                     test = test)
  Cindex[Cindex$var==allvars_to_test[i],]$OS_C_index_train <- os$C_index_train
  Cindex[Cindex$var==allvars_to_test[i],]$OS_C_index_test <- os$C_index_test
  Cindex[Cindex$var==allvars_to_test[i],]$EFS_C_index_train <- efs$C_index_train
  Cindex[Cindex$var==allvars_to_test[i],]$EFS_C_index_test <- efs$C_index_test
}
Cindex

write.table(Cindex, paste0(save_dir, "/", data_subset, "_", "Cindex.txt"))




### Plot separated (bar plot)

Cindex <- read.table(paste0(save_dir, "/", data_subset, "_", "Cindex.txt"))
Cindex$var <- gsub("_", "-", Cindex$var)

## OS
ci <- data.frame(var = Cindex$var,
                 OS = Cindex$OS_C_index_test)
ci <- ci[order(ci$OS, decreasing = T),] # order
ci$var <- factor(ci$var, levels = ci$var) # para que queden ordenados en ggplot
while (!is.null(dev.list()))  dev.off()
pdf(paste0(save_dir, "/", data_subset, "_", "Cindex_OS.pdf"), width = 8, height = 3.5) #export as pdf
ggplot(data = ci, aes(x = var, y = OS)) +
  geom_bar(stat = "identity", position=position_dodge()) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  coord_cartesian(ylim = c(0.35, 0.65)) + 
  ggtitle(data_subset) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("C Index")
while (!is.null(dev.list()))  dev.off()



## OS and EFS:
ci <- data.frame(var = Cindex$var, # create new dataset of non-splitted variables
                 OS = Cindex$OS_C_index_test,
                 EFS = Cindex$EFS_C_index_test)
ci <- ci[order(ci$OS, decreasing = T),] # order by OS
library(magrittr)
ci <- ci %>%
  tidyr::pivot_longer(cols = c(OS, EFS),
                      names_to = "survival_type",
                      values_to = "value")
ci$var <- factor(ci$var, levels = unique(ci$var)) # para que queden ordenados en ggplot
while (!is.null(dev.list()))  dev.off()
pdf(paste0(save_dir, "/", data_subset, "_", "Cindex.pdf")) #export as pdf
ggplot(data = ci, aes(x = var, y = value, fill = survival_type)) +
  geom_bar(stat = "identity", position=position_dodge()) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  ggtitle(data_subset) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("Concordance Index")
while (!is.null(dev.list()))  dev.off()









### Plot together (line plot)

Cindex_All <- read.table("./results/All_Risks_Cindex.txt")


C_OS <- data.frame(miRNA = Cindex_All$var,
                   All_Risks = Cindex_All$OS_C_index_test,
                   Low_Risk = Cindex_Low$OS_C_index_test,
                   Standard_Risk = Cindex_Standard$OS_C_index_test,
                   High_Risk = Cindex_High$OS_C_index_test)

C_OS <- C_OS[order(C_OS$All_Risks, decreasing = T),]
C_OS[,2:5] <- round(C_OS[,2:5], digits = 3)
C_OS
C_EFS



# Long format:
library(data.table)
C_OS.long <- as.data.frame(melt(setDT(C_OS), id.vars = "miRNA", variable.name = "Risk"))
C_EFS.long <- as.data.frame(melt(setDT(C_EFS), id.vars = "miRNA", variable.name = "Risk"))


C_OS.long$miRNA <- factor(C_OS.long$miRNA, levels = unique(C_OS.long$miRNA)) # para que queden ordenados en ggplot
C_EFS.long$miRNA <- factor(C_EFS.long$miRNA, levels = unique(C_EFS.long$miRNA)) # para que queden ordenados en ggplot



#----






#### Functional anaysis miRNAs ####

miRNAs <- sel_miRNAs
miRNAs <- gsub("_", "-", miRNAs)
miRNAs[2]

#BiocManager::install("multiMiR")
library(multiMiR)

example1 <- get_multimir(mirna = "hsa-miR-381-3p", summary = TRUE)
table(example1@data$support_type)
table(example1@data$experiment)
example1@data$target_symbol

unique(example1@data$target_symbol[example1@data$support_type=="positive"|example1@data$support_type=="Functional MTI"])




library()
goana("hsa-miR-381-3p")

topGO()





#----







#### Extra ####


#### AMLmiR36 Scores (immature miRNAs) ####


#search <- "30c"
#all_miRNAs[grep(search, all_miRNAs)]
#AMLmiR36_miRNAs <- c("hsa_mir_106a",
#                     "hsa_mir_584",
#                     "hsa_mir_1247",
#                     "hsa_mir_155",
#                     "hsa_mir_130b",
#                     "hsa_mir_320a",
#                     "hsa_mir_34c",
#                     "hsa_mir_30c_2",
#                     "hsa_mir_30e",
#                     "hsa_mir_450a_1", # Nos quedamos con el _1
#                     "hsa_mir_296",
#                     "hsa_mir_935",
#                     "hsa_mir_502",
#                     "hsa_mir_181b_1", # Nos quedamos con el _1
#                     "hsa_mir_363",
#                     "hsa_mir_362",
#                     "hsa_mir_100",
#                     "hsa_mir_132",
#                     "hsa_mir_340",
#                     "hsa_mir_181c",
#                     "hsa_mir_148b",
#                     "hsa_mir_4662a",
#                     "hsa_mir_1287",
#                     "hsa_mir_664b",
#                     "hsa_mir_539",
#                     "hsa_mir_217",
#                     "hsa_mir_181c", # El 181c se repite 2 veces porque en la firma del paper hay 181c-3p y 181c-5p, pero nosotros solo tenemos un 181c
#                     "hsa_mir_335",
#                     "hsa_mir_1180",
#                     "hsa_mir_202",
#                     "hsa_mir_146a",
#                     "hsa_mir_2110",
#                     "hsa_mir_375",
#                     "hsa_let_7g",
#                     "hsa_mir_139",
#                     "hsa_mir_409")
#
## Problem repeated miRNA: 5p dominant but both present in the model
## Deep dequencing extracted from miR base:
#seq_reads_mir_181c_5p <- 86561
#seq_reads_mir_181c_3p <- 26397
#prop_mir_181c_5p <- seq_reads_mir_181c_5p/(seq_reads_mir_181c_5p+seq_reads_mir_181c_3p)
#prop_mir_181c_3p <- seq_reads_mir_181c_3p/(seq_reads_mir_181c_5p+seq_reads_mir_181c_3p)
#
#AMLmiR36_coef <- c(0.244530267,
#                   0.155286578,
#                   0.070234817,
#                   0.06433791,
#                   0.042520575,
#                   0.03910599,
#                   0.027657901,
#                   0.018074572,
#                   0.017215965,
#                   0.015253859,
#                   0.015230088,
#                   0.013849354,
#                   0.013641305,
#                   0.010667128,
#                   0.006792342,
#                   0.005215312,
#                   -0.001863252,
#                   -0.002882947,
#                   -0.003799344,
#                   -0.011324553*prop_mir_181c_5p,
#                   -0.018062943,
#                   -0.019307849,
#                   -0.02419678,
#                   -0.032959784,
#                   -0.037325674,
#                   -0.040306357,
#                   -0.040989372*prop_mir_181c_3p,
#                   -0.04528897,
#                   -0.047275703,
#                   -0.048514515,
#                   -0.050372946,
#                   -0.059026016,
#                   -0.062192126,
#                   -0.062530666,
#                   -0.063621683,
#                   -0.071435852)
#
#
#
#### Getting AMLmiR36 Scores:
#
#AMLmiR36 <- data.frame(miRNA = AMLmiR36_miRNAs,
#                       Coefficient = AMLmiR36_coef)
#
## Getting Score for each patient: sum of all coef * miRNA:
#AMLmiR36_scores <- as.data.frame(colSums(expr[AMLmiR36$miRNA,]*AMLmiR36$Coefficient))
#colnames(AMLmiR36_scores) <- "AMLmiR36_score"
#
#
#### Adding AMLmiR36 score to surv_table:
#all(rownames(surv_table)==rownames(AMLmiR36_scores))
#surv_table$AMLmiR36_score <- AMLmiR36_scores$AMLmiR36_score
#surv_table_extended$AMLmiR36_score <- AMLmiR36_scores$AMLmiR36_score

#----

#### miR3 scores (immature) ####


#search <- "146b"
#all_miRNAs[grep(search, all_miRNAs)]
#miR3_miRNAs <- c("hsa_mir_146b",
#                 "hsa_mir_181c",
#                 "hsa_mir_4786")
#miR3_coef <- c(1.652,
#               -1.838,
#               -1.455)
#
#
#### Getting miR3 Scores:
#
#miR3 <- data.frame(miRNA = miR3_miRNAs,
#                   Coefficient = miR3_coef)
#
## Getting Score for each patient: sum of all coef * miRNA:
#miR3_scores <- as.data.frame(colSums(expr[miR3$miRNA,]*miR3$Coefficient))
#colnames(miR3_scores) <- "miR3_score"
#
#
#### Adding miR3 score to surv_table:
#all(rownames(surv_table)==rownames(miR3_scores))
#surv_table$miR3_score <- miR3_scores$miR3_score
#surv_table_extended$miR3_score <- miR3_scores$miR3_score

#----

#### miR4 scores (immature) ####


#search <- "146a"
#all_miRNAs[grep(search, all_miRNAs)]
#miR4_miRNAs <- c("hsa_mir_509_1",
#                 "hsa_mir_542",
#                 "hsa_mir_3667",
#                 "hsa_mir_146a")
#miR4_coef <- c(0.914,
#               0.759,
#               -0.837,
#               -0.856)
#
#
#### Getting miR4 Scores:
#
#miR4 <- data.frame(miRNA = miR4_miRNAs,
#                   Coefficient = miR4_coef)
#
## Getting Score for each patient: sum of all coef * miRNA:
#miR4_scores <- as.data.frame(colSums(expr[miR4$miRNA,]*miR4$Coefficient))
#colnames(miR4_scores) <- "miR4_score"
#
#
#### Adding miR4 score to surv_table:
#all(rownames(surv_table)==rownames(miR4_scores))
#surv_table$miR4_score <- miR4_scores$miR4_score
#surv_table_extended$miR4_score <- miR4_scores$miR4_score

#----

#### miR24 scores (immature miRNAs) ####

#search <- "_9"
#all_miRNAs[grep(search, all_miRNAs)]
#miR24_miRNAs <- c("hsa_mir_20b",
#                  "hsa_mir_223",
#                  "hsa_mir_193a",
#                  "hsa_mir_24_1",
#                  "hsa_mir_128_1",
#                  "hsa_mir_17",
#                  "hsa_mir_199b",
#                  "hsa_mir_181c",
#                  "hsa_mir_181a_1",
#                  "hsa_mir_181b_1",
#                  "hsa_mir_21",
#                  "hsa_mir_222",
#                  "hsa_mir_331",
#                  "hsa_mir_373",
#                  "hsa_mir_708",
#                  "hsa_mir_34b",
#                  "hsa_mir_195",
#                  "hsa_mir_151a",
#                  "hsa_mir_30b",
#                  "hsa_mir_22",
#                  "hsa_let_7g",
#                  "hsa_let_7i",
#                  "hsa_mir_1290",
#                  "hsa_mir_9_1")
#
#### Getting mir24 coefs (using OS):
#surv <- Surv(time = surv_vars_extended$OS, event = surv_vars_extended$status) # OS
#cox_reg_model_uni_mir24 <- coxph(formula = surv ~ ., data = surv_table_extended[,miR24_miRNAs])
#miR24_coef <- cox_reg_model_uni_mir24$coefficients
#all(names(miR24_coef) == miR24_miRNAs)#check
#miR24_coef <- unname(miR24_coef)
#
#summary(cox_reg_model_uni_mir24) #we can check model of 24 mirnas
#
#miR24 <- data.frame(miRNA = miR24_miRNAs,
#                   Coefficient = miR24_coef)
#
## Getting Score for each patient: sum of all coef * miRNA:
#miR24_scores <- as.data.frame(colSums(expr[miR24$miRNA,]*miR24$Coefficient))
#colnames(miR24_scores) <- "miR24_score"
#
#
#### Adding miR24 score to surv_table:
#all(rownames(surv_table)==rownames(miR24_scores))
#surv_table$miR24_score <- miR24_scores$miR24_score
#surv_table_extended$miR24_score <- miR24_scores$miR24_score


#----

#### C index ####

#https://www.reneshbedre.com/blog/survival-analysis.html#:~:text=and%20drug_2%20treatments.-,Cox's%20proportional%20hazards%20(CPH)%20model,not%20consider%20additional%20independent%20variables.


## Aplicar modelo:
##miRNAs <- sel_miRNAs
#surv <- Surv(time = surv_vars_extended$OS, event = surv_vars_extended$status) # OS
###cox_reg_model <- coxph(formula = surv ~ surv_table_extended[,miRNAs[1]] + surv_table_extended[,miRNAs[2]] + Risk_group + Gender + Race + Ethnicity + Age_at_Diagnosis_in_Days, data = surv_table_extended)  
#cox_reg_model <- coxph(formula = surv ~ AMLmiR36_score + Risk_group + WBC_at_Diagnosis_status + SCT_in_1st_CR + FLT3.ITD_positive. + Gender + Age_at_Diagnosis_in_Days,
#                       data = surv_table_extended)  
#summary(cox_reg_model)
### Podemos acceder a los datos asi:
#library(broom)
#results_cox <- as.data.frame(tidy(cox_reg_model))
### Hazard ratios (= exp(coef)):
#results_cox$HR <- exp(results_cox$estimate)
#results_cox
#
#ggforest(cox_reg_model, data=surv_vars_extended)
#
## ATENCION: ANALISIS ES MULTIVARIANTE (UNAS VARIABLES AFECTAN A LAS OTRAS)





### Predictive evaluation of model:
# We use Concordance Index (CI or C index):
# https://cran.r-project.org/web/packages/SurvMetrics/vignettes/SurvMetrics-vignette.html
# source("~/Desktop/Utility/Useful R codes.R")



### Examples:

## Example AMLmiR36_score OS
#get_coxph_C(dataset = surv_table,
#                       var = c("AMLmiR36_score"), # una o varias,
#                       status_var = "status",
#                       time_var = "OS",
#                       train = train,
#                       test = test)
#
## Example AMLmiR36_score OS
#get_coxph_C(dataset = surv_table_extended,
#                       var = c("AMLmiR36_score", "Gender", "Age_at_Diagnosis_in_Days"), # una o varias,
#                       status_var = "status",
#                       time_var = "OS",
#                       train = train,
#                       test = test)
#
## Example AMLmiR36_score OS
#ds <- surv_table_extended
#ds <- ds[!is.na(ds$Risk_group),]
#dim(ds)
#s3 <- sample(rownames(ds), round(length(rownames(ds))/3)) # random
#s3 <- 
#train_3 <- !rownames(ds)%in%s3
#test_3 <- rownames(ds)%in%s3
#table(ds$Risk_group, exclude = NULL)
#
#get_coxph_C(dataset = ds
#                       ,
#                       var = c("AMLmiR36_score"
#                               #, "miR4_score"
#                               #, "Gender"
#                               #, "Age_at_Diagnosis_in_Days"
#                               , "Risk_group"
#                               ), # una o varias,
#                       status_var = "status",
#                       time_var = "OS",
#                       train = train_3,
#                       test = test_3)
#
#allvars_to_test <- c("AMLmiR36_score", "miR3_score", "miR4_score", sel_miRNAs)
#
#ds <- surv_table_extended
#ds <- ds[!is.na(ds$Risk_group),]
#dim(ds)
##sample_test_names <- sample(rownames(ds), round(length(rownames(ds))/3))
#sample_test_names <- rownames(clin)[clin$Protocol=="AAML1031"]
#train_3 <- !rownames(ds)%in%sample_test_names
#test_3 <- rownames(ds)%in%sample_test_names
#table(ds$Risk_group, exclude = NULL)
#
#Ci <- data.frame(miRNA= allvars_to_test,
#                 Cindex=NA)
#
#for (i in 1:length(allvars_to_test)) {
#  Ci_res <- get_coxph_C(dataset = ds,
#                                   var = c(allvars_to_test[i]
#                                           #, "miR4_score"
#                                           #, "Gender"
#                                           #, "Age_at_Diagnosis_in_Days"
#                                           , "Risk_group"
#                                   ),
#                                   status_var = "status",
#                                   time_var = "OS",
#                                   train = train_3,
#                                   test = test_3)
#  Ci$Cindex[i] <- Ci_res$C_index_test
#}
#Ci <- Ci[order(Ci$Cindex, decreasing = T),]
#Ci


#-----



#----




