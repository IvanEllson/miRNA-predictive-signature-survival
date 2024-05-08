### Quality Control of new signature

#### Download new GEO dataset GSE97135 ####

GSE97135 <- getGEO("GSE97135")
GSE97135_metadata <- GSE97135[["GSE97135_series_matrix.txt.gz"]]@phenoData@data
GSE97135_expr <- data.frame(GSE97135[["GSE97135_series_matrix.txt.gz"]]@assayData[["exprs"]])
GSE97135_geneAnnot <- GSE97135[["GSE97135_series_matrix.txt.gz"]]@featureData@data
stopifnot(all(rownames(GSE97135_expr)==rownames(GSE97135_geneAnnot))) # check
rownames(GSE97135_expr) <- GSE97135_geneAnnot$miRNA_ID # annotate miRNAs
rownames(GSE97135_expr) <- gsub("-", "_", rownames(GSE97135_expr)) # fix miRNA IDs
write.table(GSE97135_expr, paste0(data_dir, "/miRNA_expression_data_GSE97135.txt"), sep = "\t") # save

#----

#### Load ####

newSig <- read.delim(paste0(save_dir, "/signature_new.txt"))
randomSig <- read.delim(paste0(save_dir, "/signature_random.txt"))
GSE97135_expr <- read.table(file = paste0(data_dir, "/miRNA_expression_data_GSE97135.txt"), sep = "\t")# Expression data
TARGET_AML_expr <- read.table(file = paste0(data_dir, "/fulldata_TARGET.AML.txt"), sep = "\t") # Load data
TARGET_AML_expr <- TARGET_AML_expr[TARGET_AML_expr$PatientSet=="Validation",] # Only Validation set
TARGET_AML_expr <- TARGET_AML_expr[,grepl("hsa_miR", colnames(TARGET_AML_expr), ignore.case = F) | grepl("hsa_let", colnames(TARGET_AML_expr), ignore.case = F)] # retrieve mature miRNA expression data
TARGET_AML_expr <- t(TARGET_AML_expr)

#----

#### SigQC pipeline ####

# Input of sigQC
sig_list <- list()
sig_list[["new_signature"]] <- as.matrix(newSig$miRNA)
sig_list[["random_signature"]] <- as.matrix(randomSig$miRNA)
expr_list <- list()
expr_list[["TARGET_AML_Validation_dataset"]] <- as.matrix(TARGET_AML_expr)
expr_list[["GSE97135_dataset"]] <- as.matrix(GSE97135_expr)

# Run sigQC
set.seed(1994)
make_all_plots(sig_list, expr_list, out_dir = paste0(save_dir, "/", "sigQC_results"), showResults = F)
while (!is.null(dev.list()))  dev.off()

#----

