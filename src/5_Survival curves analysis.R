### Survival curves analysis

#### Load ####

fulldata <- read.table(file = paste0(data_dir, "/fulldata_TARGET.AML.txt"), sep = "\t") # Load full data
miRNAs <- read.csv(file = paste0(data_dir, "/selected_miRNAs_TARGET.AML.csv"))[,1] # Load selected miRNAs and signatures

time_var <- "OS"#"EFS"  # Change for re-analysis

#----

#### Categorize expression using median ####

fulldata_categorical <- fulldata[,c(1:which(startsWith(colnames(fulldata), "hsa"))[1]-1, which(colnames(fulldata)%in%miRNAs))] # Clin vars and selected miRNAs
for (i in 1:length(miRNAs)) {
  median <- median(fulldata[,miRNAs[i]], na.rm = T)
  fulldata_categorical[,miRNAs[i]] <- cut(fulldata_categorical[,miRNAs[i]],
                                          breaks = c(-Inf, median, Inf),
                                          labels = c("low", "high"))} # Expression values equal to the median are assigned as low

#----

#### Kaplan Meier plots ####

vars <- miRNAs[!miRNAs%in%c("random_signature_score")]

custom_theme <- function() { # Custom theme to center title
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5, vjust=4, family = "ArialMT"),
      plot.margin = unit(c(0.4,0.4,0.4,0.4), "cm")
    )
}
p <- list()
for (i in 1:length(vars)) {
  x <- data.frame(time = fulldata_categorical[,time_var],
                  status = fulldata_categorical$status,
                  expression = fulldata_categorical[,vars[i]])
  x$expression <- as.numeric(x$expression)
  fit <- survfit(Surv(time, status) ~ expression, data = x)
  
  ## Log Rank Test:
  #log_rank_test <- survdiff(Surv(time, status) ~ expression, data = x)
  #log_rank_test
  #log_rank_test_p.val <- 1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1)
  
  # Plot 
  plot <- ggsurvplot(fit,
                     risk.table = TRUE, 
                     pval = TRUE,
                     pval.coord = c(1, 0.25),
                     conf.int = TRUE,
                     legend = "none",
                     break.x.by = 2,
                     surv.median.line = "none",
                     palette = c("#E7B800", "#2E9FDF"),
                     axes.offset = F,
                     xlim = c(0,max(x$time)),
                     xlab = "",
                     ylab = "",
                     ggtheme = custom_theme(),
                     tables.theme = theme_cleantable(),
                     surv.scale = "percent",
                     risk.table.y.text = FALSE,
                     tables.y.text = FALSE,
                     risk.table.title = "",
                     title = gsub("_", "-", vars[i]))
  p[[i]] <- plot$plot
}

# Save plots:
while (!is.null(dev.list()))  dev.off()
pdf(paste0(save_dir, "/SurvivalPlots_", time_var, ".pdf"), width = 12, height = 17, family="ArialMT") #export as pdf
args <- c(p, list(ncol = 3, left = time_var, bottom = "Time (Years)"))
do.call(grid.arrange, args)
while (!is.null(dev.list()))  dev.off()

#----

