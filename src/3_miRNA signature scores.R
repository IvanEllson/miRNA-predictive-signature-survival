### miRNA signature scores

#### Load ####

fulldata <- read.table(file = paste0(data_dir, "/processed_data_TARGET.AML.txt"), sep = "\t") # Load data
AMLmiR36 <- read.delim(paste0(save_dir, "/signature_AMLmiR36.txt"))
miR3 <- read.delim(paste0(save_dir, "/signature_miR3.txt"))
miR4 <- read.delim(paste0(save_dir, "/signature_miR4.txt"))
miR24 <- read.delim(paste0(save_dir, "/signature_miR24.txt"))
newSig <- read.delim(paste0(save_dir, "/signature_new.txt"))
randomSig <- read.delim(paste0(save_dir, "/signature_random.txt"))

#----

#### Calculate all miRNA signature scores: ####

# Getting Score for each patient: sum of all coef * miRNA log(cpm) normalized expresion 

# AMLmiR36 Scores
AMLmiR36_scores <- as.data.frame(colSums(expr[AMLmiR36$miRNA,]*AMLmiR36$Coefficient)) # 
colnames(AMLmiR36_scores) <- "AMLmiR36_score"
fulldata$AMLmiR36_score <- AMLmiR36_scores$AMLmiR36_score # Adding AMLmiR36 score to fulldata

# miR3 scores
miR3_scores <- as.data.frame(colSums(expr[miR3$miRNA,]*miR3$Coefficient))
colnames(miR3_scores) <- "miR3_score"
fulldata$miR3_score <- miR3_scores$miR3_score # Adding AMLmiR36 score to fulldata

# miR4 scores
miR4_scores <- as.data.frame(colSums(expr[miR4$miRNA,]*miR4$Coefficient))
colnames(miR4_scores) <- "miR4_score"
fulldata$miR4_score <- miR4_scores$miR4_score # Adding AMLmiR36 score to fulldata

# miR24 scores
miR24_scores <- as.data.frame(colSums(expr[miR24$miRNA,]*miR24$Coefficient))
colnames(miR24_scores) <- "miR24_score"
fulldata$miR24_score <- miR24_scores$miR24_score # Adding miR24 score to fulldata

# new signature scores
newSig_scores <- as.data.frame(colSums(expr[newSig$miRNA,]*newSig$Coefficient))
colnames(newSig_scores) <- "newSig_score"
fulldata$newSig_score <- newSig_scores$newSig_score # Adding newSig score to fulldata

# random signature scores
random_signature_scores <- as.data.frame(colSums(expr[randomSig$miRNA,]*randomSig$Coefficient))
colnames(random_signature_scores) <- "random_signature_score"
fulldata$random_signature_score <- random_signature_scores$random_signature_score # Adding random signature score to fulldata

#----

#### Select all miRNAs and scores to analyze ####

# Individually selected miRNAs:
sel_miRNAs <- c("hsa_mir_199a", 
                "hsa_mir_381",
                "hsa_miR_206",
                "hsa_miR_193b_3p",
                "hsa_miR_139_5p",
                "hsa_mir_100",
                "hsa_mir_195",
                "hsa_miR_375",
                "hsa_mir_335",
                "hsa_mir_370",
                "hsa_mir_29a",
                "hsa_mir_122",
                "hsa_mir_192",
                "hsa_mir_155",
                "hsa_mir_196b",
                "hsa_mir_34b")

# Signatures:
sel_signatures <- c("newSig_score", 
                    "AMLmiR36_score", 
                    "miR24_score", 
                    "miR4_score", 
                    "miR3_score", 
                    "random_signature_score")

miRNAs <- c(sel_signatures, sel_miRNAs)

#----

#### Save ####

write.table(fulldata, file = paste0(data_dir, "/fulldata_TARGET.AML.txt"), sep = "\t") # save full data
write.csv(miRNAs, file = paste0(data_dir, "/selected_miRNAs_TARGET.AML.csv"), row.names = F) # save selected miRNAs and signatures

#----



