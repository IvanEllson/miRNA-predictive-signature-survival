### Explore survival by patient groups 

#### Load ####

fulldata <- read.table(file = paste0(data_dir, "/fulldata_TARGET.AML.txt"), sep = "\t") # Load data

time_var <- "OS"#"EFS"  # Change for re-analysis

#----

#### Compare outcomes across patient sets ####

custom_theme <- function() { # Custom theme to center title
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5, vjust=4, family = "ArialMT"),
      plot.margin = unit(c(0.4,0.4,0.4,0.4), "cm"))}
plot_protocol <- ggsurvplot(survfit(Surv(fulldata[,time_var], fulldata$status) ~ Protocol, data = fulldata),
                            risk.table = TRUE, 
                            pval = TRUE,
                            pval.coord = c(1, 0.25),
                            conf.int = TRUE,
                            break.x.by = 2,
                            surv.median.line = "none",
                            axes.offset = F,
                            xlim = c(0,max(fulldata[,time_var])),
                            xlab = "Time (years)",
                            ylab = time_var,
                            ggtheme = custom_theme(),
                            tables.theme = theme_cleantable(),
                            surv.scale = "percent",
                            risk.table.y.text = FALSE,
                            tables.y.text = FALSE, 
                            legend = "right",
                            legend.title="",
                            risk.table.title = "",
                            title = "Survival by protocol")
plot_patientSet <- ggsurvplot(survfit(Surv(fulldata[,time_var], fulldata$status) ~ PatientSet, data = fulldata),
                              risk.table = TRUE, 
                              pval = TRUE,
                              pval.coord = c(1, 0.25),
                              conf.int = TRUE,
                              break.x.by = 2,
                              surv.median.line = "none",
                              axes.offset = F,
                              xlim = c(0,max(fulldata[,time_var])),
                              xlab = "Time (years)",
                              ylab = time_var,
                              ggtheme = custom_theme(),
                              tables.theme = theme_cleantable(),
                              surv.scale = "percent",
                              risk.table.y.text = FALSE,
                              tables.y.text = FALSE, 
                              legend = "right",
                              legend.title="",
                              risk.table.title = "",
                              title = "Survival by patient set")

while (!is.null(dev.list()))  dev.off()
pdf(paste0(save_dir, "/Groups_SurvivalPlots_", time_var, ".pdf"), width = 8, height = 5, family="ArialMT")
plot_protocol
plot_patientSet
while (!is.null(dev.list()))  dev.off()

#----

