### Preprocess data

#### Load data ####

expr <- read.delim(paste0(data_dir, "/miRNA_expression_data_TARGET.AML_raw.txt"))
clin <- read.delim(paste0(data_dir, "/clinical_data_TARGET.AML_raw.txt"), check.names = F)

#----

#### Preprocess data ####

# Reduce expression to samples with clinical data:
expr <- expr[,rownames(clin)]

# Change clinical variable names
colnames(clin) <- gsub("\\.", "_", colnames(clin))
colnames(clin) <- gsub("\\/", "_", colnames(clin))
colnames(clin) <- gsub("\\?", "_", colnames(clin))
colnames(clin) <- gsub("%", "perc", colnames(clin)) 
colnames(clin) <- gsub(" ", "_", colnames(clin)) 


### Fix clin variables

# Risk group
table(clin$Risk_group, exclude = NULL)
clin$Risk_group[!(clin$Risk_group=="High Risk" | clin$Risk_group=="Standard Risk" | clin$Risk_group=="Low Risk")] <- NA
table(clin$Risk_group, exclude = NULL)

# FLT3_ITD_positive
table(clin$FLT3_ITD_positive_, exclude = NULL)
clin$FLT3_ITD_positive_[grepl("No", clin$FLT3_ITD_positive_, ignore.case = T)] <- "No"
clin$FLT3_ITD_positive_[grepl("Yes", clin$FLT3_ITD_positive_, ignore.case = T)] <- "Yes"
clin$FLT3_ITD_positive_[grepl("Unknown", clin$FLT3_ITD_positive_, ignore.case = T)] <- NA
table(clin$FLT3_ITD_positive_, exclude = NULL)

# SCT_in_1st_CR
table(clin$SCT_in_1st_CR, exclude = NULL)
clin$SCT_in_1st_CR[clin$SCT_in_1st_CR=="Unknown"] <- NA
table(clin$SCT_in_1st_CR, exclude = NULL)

# CNS_disease
table(clin$CNS_disease, exclude = NULL)
clin$CNS_disease[clin$CNS_disease=="Unknown"] <- NA
table(clin$CNS_disease, exclude = NULL)

# WBC_at_Diagnosis
clin$WBC_at_Diagnosis
summary(as.numeric(clin$WBC_at_Diagnosis))

# MRD_at_end_of_course_1 and 2
table(clin$MRD_at_end_of_course_1, exclude = NULL)
table(clin$MRD_at_end_of_course_2, exclude = NULL)
clin$MRD_at_end_of_course_1[clin$MRD_at_end_of_course_1 %in% c("Unknown")] <- NA
clin$MRD_at_end_of_course_2[clin$MRD_at_end_of_course_2 %in% c("Unknown")] <- NA
clin$MRD___at_end_of_course_1 <- gsub(" \\|\\|\\|.+", "", clin$MRD_perc_at_end_of_course_1)
clin$MRD_perc_at_end_of_course_2 <- gsub(" \\|\\|\\|.+", "", clin$MRD_perc_at_end_of_course_2)
clin$MRD_perc_at_end_of_course_1 <- as.numeric(clin$MRD_perc_at_end_of_course_1)
clin$MRD_perc_at_end_of_course_2 <- as.numeric(clin$MRD_perc_at_end_of_course_2)
table(clin$MRD_at_end_of_course_1, exclude = NULL)
table(clin$MRD_at_end_of_course_2, exclude = NULL)
summary(as.numeric(clin$MRD_perc_at_end_of_course_1))
summary(as.numeric(clin$MRD_perc_at_end_of_course_2))

# NPM_mutation
table(clin$NPM_mutation, exclude = NULL)
clin$NPM_mutation[grepl("no", clin$NPM_mutation, ignore.case = T) & !grepl("yes", clin$NPM_mutation, ignore.case = T)] <- "No"
clin$NPM_mutation[!clin$NPM_mutation %in% c("No", "Yes")] <- NA
table(clin$NPM_mutation, exclude = NULL)

# CEBPA_mutation
table(clin$CEBPA_mutation, exclude = NULL)
clin$CEBPA_mutation[grepl("no", clin$CEBPA_mutation, ignore.case = T) & clin$CEBPA_mutation!="Unknown"] <- "No"
clin$CEBPA_mutation[!clin$CEBPA_mutation %in% c("No", "Yes")] <- NA
table(clin$CEBPA_mutation, exclude = NULL)

# cytogenetic abnormalities:
for (i in 25:42) {
  print(table(clin[,i], exclude = NULL))
}
for (i in 25:42) {
  clin[,i][clin[,i]=="Unknown"] <- NA
}

# Fix survival variables: choose the max survival of several values (last measures)
clin$Overall_Survival_Time_in_Days <- sapply(strsplit(as.character(clin$Overall_Survival_Time_in_Days)," \\|\\|\\| "), function(x) max(as.numeric(x)))
clin$Event_Free_Survival_Time_in_Days <- sapply(strsplit(as.character(clin$Event_Free_Survival_Time_in_Days)," \\|\\|\\| "), function(x) max(as.numeric(x)))

# Survival variables in years
clin$OS <- (clin$Overall_Survival_Time_in_Days+1)/365
clin$EFS <- (clin$Event_Free_Survival_Time_in_Days+1)/365

# Create vital status variable 0/1
clin$status <- 0 
clin$status[clin$Vital_Status=="Dead"] <- 1

# Age in years
clin$Age <- clin$Age_at_Diagnosis_in_Days/365

#----

#### Patient filter ####

# Select only AML patients with survival data and and age < 21 years:
clin <- clin[clin$tumor_code==20 &
               !is.na(clin$OS) &
               !is.na(clin$EFS) &
               !is.na(clin$Vital_Status),]

any(is.na(clin$OS)) # No NAs
any(is.na(clin$Vital_Status)) # No NAs
any(is.na(clin$Age)) # No NAs

clin <- clin[clin$Age<21,] # <21 years

# Apply filter to expr dataset:
expr <- expr[,rownames(clin)]

rownames(expr) <- gsub("-", "_", rownames(expr)) # fix miRNA IDs

stopifnot(all(colnames(expr)==rownames(clin))) # check

#----

#### Save ####

stopifnot(all(colnames(expr)==rownames(clin))) # check
write.table(expr, file = paste0(data_dir, "/miRNA_expression_data_TARGET.AML.txt"), sep = "\t")
write.table(clin, file = paste0(data_dir, "/clinical_data_TARGET.AML.txt"), sep = "\t")

#----




