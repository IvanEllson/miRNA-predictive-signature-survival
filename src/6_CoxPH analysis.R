### CoxPH analysis

#### Load ####

fulldata <- read.table(file = paste0(data_dir, "/fulldata_TARGET.AML.txt"), sep = "\t") # Load full data
miRNAs <- read.csv(file = paste0(data_dir, "/selected_miRNAs_TARGET.AML.csv"))[,1] # Load selected miRNAs and signatures
covars <- read.csv(file = paste0(data_dir, "/covars_TARGET.AML.csv"))[,1] # Load covars
covars_trueNames <- c("WBC at Diagnosis",
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

time_var <- "OS"#"EFS"  # Change for re-analysis

#----

#### Multivariate Cox PH (explore covariates): ####

surv <- Surv(time = fulldata[,time_var], event = fulldata$status)

# Multivariate cox analysis:
cox_model_onlyCovars <- coxph(formula = as.formula(paste("surv ~ ", paste(covars, collapse = " + "))),
                              data = fulldata) 

multimod <- summary(cox_model_onlyCovars)

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
cox$predictors <- covars_trueNames

# Plot:
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
pdf(paste0(save_dir, "/coxph_onlyCovars_", time_var, ".pdf"), width = 12, height = 7, family="ArialMT") # Save pdf
p
while (!is.null(dev.list()))  dev.off()

#----

#### Cox Proportional Hazards model (adding all covariates): ####

vars <- miRNAs[!miRNAs%in%c("random_signature_score")]

# Apply custom function:
p <- perform_uni_multi_coxph(
  data = fulldata,
  time_var = time_var,
  status_var = "status",
  predictors = vars,
  covars = covars,
  save_table = F
)

while (!is.null(dev.list()))  dev.off()
pdf(paste0(save_dir, "/coxph_withAllCovars_", time_var, ".pdf"), width = 12, height = 7, family="ArialMT") # Save pdf
p
while (!is.null(dev.list()))  dev.off()

#----

#### Cox Proportional Hazards model (adding significant covariates p<0.005): ####

vars <- miRNAs[!miRNAs%in%c("random_signature_score")]

sig.covars <- c("SCT_in_1st_CR",
                "MRD_at_end_of_course_1",
                "NPM_mutation",
                "CEBPA_mutation",
                "t.8.21.",
                "inv.16.")

# Apply custom function:
p <- perform_uni_multi_coxph(
  data = fulldata,
  time_var = time_var,
  status_var = "status",
  predictors = vars,
  covars = sig.covars,
  save_table = F
)

while (!is.null(dev.list()))  dev.off()
pdf(paste0(save_dir, "/coxph_withSignifCovars_", time_var, ".pdf"), width = 12, height = 7, family="ArialMT") # Save pdf
p
while (!is.null(dev.list()))  dev.off()


#----

