### Download the data

#### Expression data (mature miRNAs) from TARGET-AML: ####

### Download expression data

query_TARGET_AML <- GDCquery(project = "TARGET-AML",
                             data.category = "Transcriptome Profiling",
                             data.type = "Isoform Expression Quantification",
                             workflow.type = "BCGSC miRNA Profiling")
GDCdownload(query_TARGET_AML, directory = "./data/GDCdata")
miRNA_expr_long <- GDCprepare(query = query_TARGET_AML, 
                                         directory = "./data/GDCdata",
                                         summarizedExperiment = T)

# Explore:
table(gsub(",MIMAT.*", "", miRNA_expr_long$miRNA_region), exclude = NULL) # check types of miRNAs
miRNA_expr_long <- miRNA_expr_long[startsWith(miRNA_expr_long$miRNA_region, "mature"),] # Get rid of immature miRNAs
miRNA_expr_long$accession <- gsub("mature,", "", miRNA_expr_long$miRNA_region)
table(miRNA_expr_long$barcode)[1:5] # different miRNA isoforms per sample


### Generate mature miRNA counts:

expr <- data.frame("accession" = na.omit(unique(miRNA_expr_long$accession)))
miRversion <- checkMiRNAVersion(unique(miRNA_expr_long$miRNA_ID)) # Check miRbase version
expr$miRNA_ID_mature <- miRNA_AccessionToName(expr$accession, targetVersion = miRversion)[,2] # get miRNA mature ID
rownames(expr) <- expr$accession
expr <- expr[,-1,drop=F]

# To collapse miRNA isoforms to one mature miRNA: Sum of all counts by MIMAT accession (https://github.com/bcgsc/mirna):
for (i in 1:length(unique(miRNA_expr_long$barcode))) {
  patient <- miRNA_expr_long[miRNA_expr_long$barcode==unique(miRNA_expr_long$barcode)[i],c("accession", "read_count")]
  patient_uniqueAccession <- aggregate(read_count~accession, patient, sum) # collapse to unique accession by sum
  expr[,unique(miRNA_expr_long$barcode)[i]] <- 0 # Not measured as 0 counts
  expr[patient_uniqueAccession$accession,unique(miRNA_expr_long$barcode)[i]] <- patient_uniqueAccession$read_count
}

rownames(expr) <- expr$miRNA_ID_mature
expr <- expr[,-1]
any(apply(expr, 1, function(x) any(is.na(x)))) # check nas values

colnames(expr) <- gsub("-", "_", colnames(expr)) # Sample IDs

#----

#### Clinical data from TARGET-AML: ####

### Download
query_clinical_TARGET_AML <- GDCquery(project = "TARGET-AML",
                                      data.category = "Clinical",
                                      data.type = "Clinical Supplement",
                                      data.format = "xlsx")
GDCdownload(query_clinical_TARGET_AML, directory = "./data/GDCdata")


### Load
clinfiles_path <- "./data/GDCdata/TARGET-AML/Clinical/Clinical_Supplement"
clinfiles <- list.files(path = clinfiles_path, recursive = T)
TARGET_AML_ClinicalData_LowDepthRNAseq <- read_excel(paste0(clinfiles_path, "/", clinfiles[1]))
TARGET_AML_ClinicalData_Discovery <- read_excel(paste0(clinfiles_path, "/", clinfiles[2]))
TARGET_AML_ClinicalData_CDE <- read_excel(paste0(clinfiles_path, "/", clinfiles[3])) # variable dictionary
TARGET_AML_ClinicalData_Validation <- read_excel(paste0(clinfiles_path, "/", clinfiles[4]))
TARGET_AML_ClinicalData_AML1031 <- read_excel(paste0(clinfiles_path, "/", clinfiles[5]))
TARGET_AML_ClinicalData_AAML1031_AAML0631_additional <- read_excel(paste0(clinfiles_path, "/", clinfiles[6]))
TARGET_AML_ClinicalData_Discovery_and_Validation_Tumor_Content_and_RIN_Supplement <- read_excel(paste0(clinfiles_path, "/", clinfiles[7])) # Sample info: anatomic site, tumor content, etc

colnames(TARGET_AML_ClinicalData_LowDepthRNAseq)[colnames(TARGET_AML_ClinicalData_LowDepthRNAseq)=="Gene Fusion...67"] <- "Gene Fusion" # Modify


### Join data
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
df_list <- list(TARGET_AML_ClinicalData_LowDepthRNAseq, # join clinical data frames
                TARGET_AML_ClinicalData_Validation,
                TARGET_AML_ClinicalData_AML1031,
                TARGET_AML_ClinicalData_Discovery,
                TARGET_AML_ClinicalData_AAML1031_AAML0631_additional)

# Merge all data frames in list:
sample_info_clinical <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)

sample_info_clinical[sample_info_clinical=="NA"] <- NA # Convert "NA" character to NA



# Reduce and collapse: Reduce data frame eliminating duplicated elements based on one column but collapsing different values into one row:
sample_info_clinical <- df.reduce_by_id(sample_info_clinical, id_col = "TARGET USI", collapse_by = " ||| ")

# Add extra data
all_samples_clin_data <- merge(sample_info_clinical, 
                               TARGET_AML_ClinicalData_Discovery_and_Validation_Tumor_Content_and_RIN_Supplement,
                               by = "TARGET USI",
                               all = T)

names(all_samples_clin_data)[names(all_samples_clin_data)=='TARGET USI'] <- 'TARGET_USI'
all_samples_clin_data$TARGET_USI <- gsub("-", "_", all_samples_clin_data$TARGET_USI) # Modify sample names



### Clinical data coded in TARGET IDs

# Check TARGET codes at https://ocg.cancer.gov/sites/default/files/OCG-Project-Codes_Tissues-and-Samples_03-02-2021.pdf:
seq_samples_clin <- data.frame(sample_code = colnames(expr),
                               project = vapply(strsplit(colnames(expr), '_'), `[`, 1, FUN.VALUE=character(1)),
                               tumor_code = vapply(strsplit(colnames(expr), '_'), `[`, 2, FUN.VALUE=character(1)),
                               case_identifier = vapply(strsplit(colnames(expr), '_'), `[`, 3, FUN.VALUE=character(1)),
                               tissue_code = vapply(strsplit(colnames(expr), '_'), `[`, 4, FUN.VALUE=character(1)),
                               nucleic_acid = vapply(strsplit(colnames(expr), '_'), `[`, 5, FUN.VALUE=character(1)),
                               TARGET_USI = vapply(strsplit(colnames(expr), '_'), function(x)
                                 paste(x[seq.int(3)], collapse='_'), character(1L)))

# Filter samples in order to keep interesting groups and remove patients with duplicated samples:
# Keep de novo AML (tumor_code = 20) and non-cancerous tissue (tumor_code = 00). Remove Induction failure AML (tumor_code = 21):
seq_samples_clin <- seq_samples_clin[seq_samples_clin$tumor_code == "20" | seq_samples_clin$tumor_code == "00",]
# Keep Bone Marrow (04A, 09A) y Bone Marrow Normal (14A):
seq_samples_clin <- seq_samples_clin[seq_samples_clin$tissue_code == "04A" | seq_samples_clin$tissue_code == "09A" | seq_samples_clin$tissue_code == "14A",]
# Keep Primary blood derived cancer (09). Remove Recurrent (04):
seq_samples_clin <- seq_samples_clin[!seq_samples_clin$tissue_code == "04A",]
# Keep last aliquot (highest nucleic_acid):
alldups <- seq_samples_clin[duplicated(seq_samples_clin$case_identifier) | duplicated(seq_samples_clin$case_identifier, fromLast = T),] # see duplicates
alldups <- alldups[order(alldups$case_identifier),]
c <- unique(alldups$case_identifier)
dup_to_remove <- list()
for (i in 1:length(c)) {
  dup_to_remove[i] <- alldups[alldups$case_identifier==c[i],][order(alldups[alldups$case_identifier==c[i],]$nucleic_acid, decreasing = T),][2,]$sample_code
}
seq_samples_clin <- seq_samples_clin[!seq_samples_clin$sample_code%in%dup_to_remove,]

# Explore:
table(seq_samples_clin$tumor_code, seq_samples_clin$tissue_code)
table(duplicated(seq_samples_clin$case_identifier)) # no duplicated patients

### Join clinical datasets and select patients

table(seq_samples_clin$TARGET_USI%in%all_samples_clin_data$TARGET_USI, seq_samples_clin$tumor_code) # healthy samples and 2 patients without clinical data
not_in_clin <- seq_samples_clin[!seq_samples_clin$TARGET_USI%in%all_samples_clin_data$TARGET_USI,]
to_remove <- not_in_clin[not_in_clin$tumor_code!="00",]$TARGET_USI
seq_samples_clin <- seq_samples_clin[!seq_samples_clin$TARGET_USI%in%to_remove,] # Remove 2 patients without clinical data

clin <- merge(seq_samples_clin, all_samples_clin_data, by = "TARGET_USI", all.x = T) # Merge all clinical data
rownames(clin) <- clin$sample_code

#----

#### Save unprocessed data ####

write.table(expr, file = paste0(data_dir, "/miRNA_expression_data_TARGET.AML_raw.txt"), sep = "\t") # save unprocessed expr
write.table(clin, file = paste0(data_dir, "/clinical_data_TARGET.AML_raw.txt"), sep = "\t") # save unprocessed expr

#----



