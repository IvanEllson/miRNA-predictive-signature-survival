### Process and normalize data 

#### Load data  ####

clin <- read.table(file = paste0(data_dir, "/", "clinical_data_TARGET.AML.txt"), sep = "\t")# Clinical data
expr <- read.table(file = paste0(data_dir, "/", "miRNA_expression_data_TARGET.AML.txt"), sep = "\t")# Expression data
stopifnot(all(colnames(expr)==rownames(clin))) # check

#----

#### Patient selection by risk (optional) ####

data_subset <- "All_Risks"
#data_subset <- "High_Risk"
#data_subset <- "Standard_Risk"
#data_subset <- "Low_Risk"

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

#### Discovery - Validation patient distribution ####

nrow(clin)
table(clin$Protocol, exclude = NULL)
Discovery_set <- rownames(clin)[clin$Protocol=="AAML0531"|clin$Protocol=="AAML03P1"|clin$Protocol=="CCG2961"]
Validation_set <- rownames(clin)[clin$Protocol=="AAML1031"]
clin$PatientSet <- "Discovery"
clin[Validation_set,]$PatientSet <- "Validation"

#----

#### Immature MiRNA expression calculation ####

# Calculate immature miRNAs counts (by sum of matures counts (unnormalized)):
newmir_df <- data.frame(row.names = colnames(expr),
                        "hsa_mir_199a" = colSums(expr[c("hsa_miR_199a_3p", "hsa_miR_199a_5p"),]), # individual immature miRNAs
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
expr <- rbind(expr, t(newmir_df)) # add new immature mirnas to expr

#----

#### Normalize expression to log2(CPM) ####

# TMM normalization:
d <- DGEList(counts=expr, genes=rownames(expr))
d <- calcNormFactors(d, method="TMM")

# log(CPM) transformation:
expr <- cpm(d, normalized.lib.sizes = T, log = T)

#----

#### Create data frame with interesting variables and miRNA expression ####

# Select covariates to latter analysis:
covars <- c("WBC_at_Diagnosis", # Numeric - 0% NAs
            "SCT_in_1st_CR",
            "FLT3_ITD_positive_",
            "CNS_disease",
            "MRD_at_end_of_course_1", # 10% NAs
            "NPM_mutation",
            "CEBPA_mutation",
            colnames(clin[,25:42])) #cytogenetic abnormalities

clin_selected <- clin[,c("OS",
                         "EFS",
                         "status",
                         "Vital_Status",
                         "Protocol",
                         "PatientSet",
                         "Risk_group",
                         covars)]

# Add miRNA expression:
fulldata <- merge(clin_selected,
                  as.data.frame(t(expr)),
                  by = "row.names")

rownames(fulldata) <- fulldata$Row.names



# Explore NAs in fulldata:
for (i in 1:length(covars)) {
  cat(paste0(covars[i], ":"))
  if (!is.numeric(fulldata[,covars[i]])) {
    print(table(fulldata[,covars[i]], exclude = NULL))
  }
  print(paste0("Prop NAs: ", mean(is.na(fulldata[,covars[i]]))))
  cat("\n")
}


# Change Unknown to NA:
for (i in 1:length(covars)) {
  if (!is.numeric(fulldata[,covars[i]])) {
    fulldata[,covars[i]][fulldata[,covars[i]] == "Unknown"] <- NA
  }
}

# Remove covars with less than 20 patients with Yes:
newcovars <- covars
for (i in 1:length(covars)) {
  if (!is.numeric(fulldata[,covars[i]])) {
    if (sum(fulldata[,covars[i]] == "Yes", na.rm = T) < 20) {
      print(paste0(covars[i], " removed"))
      newcovars <- newcovars[newcovars != covars[i]]
    }
  }
}
covars <- newcovars



#----

#### Save ####

write.table(fulldata, file = paste0(data_dir, "/processed_data_TARGET.AML.txt"), sep = "\t") # save
write.csv(covars, file = paste0(data_dir, "/covars_TARGET.AML.csv"), row.names = F) # save

#----



