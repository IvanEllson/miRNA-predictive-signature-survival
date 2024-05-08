### Build new model and collect previous signatures

#### Load ####

fulldata <- read.table(file = paste0(data_dir, "/processed_data_TARGET.AML.txt"), sep = "\t") # Load data
covars <- read.csv(file = paste0(data_dir, "/covars_TARGET.AML.csv"))[,1] # Load covars

#----

#### Develop a new signature using Lasso CoxPH ####

expr_notnorm_mature <- read.table(file = paste0(data_dir, "/", "miRNA_expression_data_TARGET.AML.txt"), sep = "\t") # Retrieve unnormalized all-mature miRNA expression
expr_notnorm_mature <- expr_notnorm_mature[,rownames(fulldata)]
stopifnot(all(rownames(expr_notnorm_mature)%in%colnames(fulldata[fulldata$PatientSet=="Discovery",]))) # check

# Filter miRNAs with very low counts:
keep <- filterByExpr(expr_notnorm_mature, min.count = 10)
table(keep)
expr_notnorm_mature_filt <- expr_notnorm_mature[keep,]

# Select miRNAs and compute penalized Lasso CoxPH model using biospear package:
set.seed(1994)
fit_lasso <- BMsel(data = fulldata[fulldata$PatientSet=="Discovery", c("OS", "status", rownames(expr_notnorm_mature_filt))], # use Discovery set and filtered miRNAs
                   x = rownames(expr_notnorm_mature_filt), # All mature miRNAs as biomarkers
                   y = c("OS", "status"), # OS based
                   inter = FALSE,
                   method = "lasso",
                   folds = 5)

length(fit_lasso$lasso)
if (length(fit_lasso$lasso)<=0) {stop("No predictors from lasso model")}

# Getting new signature scores:
newSig <- data.frame(miRNA = names(fit_lasso$lasso),
                     Coefficient = unname(fit_lasso$lasso))

#----

#### AMLmiR36 signature ####

# miRNAs and scores obtained from article https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5721230/
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
stopifnot(all(AMLmiR36_miRNAs %in% rownames(expr))) # check
AMLmiR36 <- data.frame(miRNA = AMLmiR36_miRNAs,
                       Coefficient = AMLmiR36_coef)

#----

#### miR3 signature ####

# miRNAs and scores obtained from article https://pubmed.ncbi.nlm.nih.gov/28473658/
miR3_miRNAs <- c("hsa_mir_146b",
                 "hsa_mir_181c",
                 "hsa_mir_4786")
miR3_coef <- c(1.652,
               -1.838,
               -1.455)
stopifnot(all(miR3_miRNAs%in%rownames(expr))) # check
miR3 <- data.frame(miRNA = miR3_miRNAs,
                   Coefficient = miR3_coef)

#----

#### miR4 signature ####

# miRNAs and scores obtained from article https://pubmed.ncbi.nlm.nih.gov/30242879/
miR4_miRNAs <- c("hsa_mir_509",
                 "hsa_mir_542",
                 "hsa_mir_3667",
                 "hsa_mir_146a")
miR4_coef <- c(0.914,
               0.759,
               -0.837,
               -0.856)
stopifnot(all(miR3_miRNAs%in%rownames(expr))) # check
miR4 <- data.frame(miRNA = miR4_miRNAs,
                   Coefficient = miR4_coef)

#----

#### miR24 signature ####

# miRNAs obtained from article https://pubmed.ncbi.nlm.nih.gov/36951259/
miR24_miRNAs <- c("hsa_miR_20b_5p",
                  "hsa_miR_223_3p",
                  "hsa_miR_193a_3p",
                  "hsa_miR_24_3p",
                  "hsa_miR_128_3p",
                  "hsa_miR_17_5p",
                  "hsa_miR_199b_5p",
                  "hsa_miR_181c_5p",
                  "hsa_miR_181a_5p",
                  "hsa_miR_181b_5p",
                  "hsa_miR_21_5p",
                  "hsa_miR_222_5p",
                  "hsa_miR_331_5p",
                  "hsa_miR_373_3p",
                  "hsa_miR_708_5p",
                  "hsa_miR_34b_5p",
                  "hsa_miR_195_5p",
                  "hsa_miR_151a_5p",
                  "hsa_miR_30b_5p",
                  "hsa_miR_22_3p",
                  "hsa_let_7g_5p",
                  "hsa_let_7i_5p",
                  "hsa_miR_1290",
                  "hsa_miR_9_5p")
stopifnot(all(miR24_miRNAs%in%rownames(expr))) # check

# Getting mir24 coefs (using OS) in Discovery data: 
set.seed(1994)
surv <- Surv(time = fulldata[fulldata$PatientSet=="Discovery",]$OS, event = fulldata[fulldata$PatientSet=="Discovery",]$status) # OS based
cox_model_uni_mir24 <- coxph(formula = surv ~ ., data = fulldata[fulldata$PatientSet=="Discovery",][,miR24_miRNAs]) # Univariate CoxPH model to obtain coefs
summary(cox_model_uni_mir24) #we can check model of 24 mirnas
miR24_coef <- cox_model_uni_mir24$coefficients
stopifnot(all(names(miR24_coef) == miR24_miRNAs)) # check
miR24_coef <- unname(miR24_coef)
miR24 <- data.frame(miRNA = miR24_miRNAs,
                    Coefficient = miR24_coef)

#----

#### Random miRNA signature ####

set.seed(1994)
sig_length <- nrow(newSig) # random signature of same length as our new signature
random_signature_miRNAs <- sample(rownames(expr_notnorm_mature_filt), sig_length) # selected from all mature miRNAs filtered

### Getting coefs (using OS) in Discovery Data: 
surv <- Surv(time = fulldata[fulldata$PatientSet=="Discovery",]$OS, event = fulldata[fulldata$PatientSet=="Discovery",]$status) # OS based
cox_model_uni_random_signature <- coxph(formula = surv ~ ., data = fulldata[fulldata$PatientSet=="Discovery",random_signature_miRNAs]) # univariate Cox to obtain coefs
summary(cox_model_uni_random_signature) # we can check model of random mirnas
stopifnot(all(names(cox_model_uni_random_signature$coefficients) == random_signature_miRNAs)) # check
randomSig <- data.frame(miRNA = random_signature_miRNAs,
                        Coefficient = unname(cox_model_uni_random_signature$coefficients))

#----

#### Save ####

write.table(newSig, paste0(save_dir, "/signature_new.txt"), sep = "\t", row.names = F) # save
write.table(AMLmiR36, paste0(save_dir, "/signature_AMLmiR36.txt"), sep = "\t", row.names = F) # save
write.table(miR3, paste0(save_dir, "/signature_miR3.txt"), sep = "\t", row.names = F) # save
write.table(miR4, paste0(save_dir, "/signature_miR4.txt"), sep = "\t", row.names = F) # save
write.table(miR24, paste0(save_dir, "/signature_miR24.txt"), sep = "\t", row.names = F) # save
write.table(randomSig, paste0(save_dir, "/signature_random.txt"), sep = "\t", row.names = F) # save

#----


