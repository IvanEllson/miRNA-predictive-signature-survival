### Predictive accuracy evaluation

#### Load ####

fulldata <- read.table(file = paste0(data_dir, "/fulldata_TARGET.AML.txt"), sep = "\t") # Load full data
miRNAs <- read.csv(file = paste0(data_dir, "/selected_miRNAs_TARGET.AML.csv"))[,1] # Load selected miRNAs and signatures

time_var <- "OS"#"EFS"  # Change for re-analysis
folds <- 10 #5 #select number of folds

#----

#### Train Test patient selection in Validation set (Cross Validation) ####

# Division of datasets in train and test by random X fold Cross Validation:
set.seed(1994)
CV <- crossv_kfold(fulldata[fulldata$PatientSet=="Validation",], k=folds) # In Validation Set
CV

#----

#### C index with CrossValidation ####

vars <- miRNAs

# Apply custom function:
CindexByFold <- list()
for (fold in 1:nrow(CV)) {
  Cindex <- data.frame(var = vars,
                       C_index_train = NA,
                       C_index_test = NA)
  for (i in 1:length(vars)) {
    C <- get_coxph_C(dataset = fulldata[fulldata$PatientSet=="Validation",], # Using Validation set
                     var = c(vars[i]),
                     status_var = "status",
                     time_var = time_var,
                     train = rownames(fulldata[fulldata$PatientSet=="Validation",])[CV[[1]][[fold]]$idx],
                     test = rownames(fulldata[fulldata$PatientSet=="Validation",])[CV[[2]][[fold]]$idx])
    Cindex[Cindex$var==vars[i],]$C_index_train <- C$C_index_train
    Cindex[Cindex$var==vars[i],]$C_index_test <- C$C_index_test
  }
  colnames(Cindex)[-1] <- paste0(colnames(Cindex)[-1], "_fold", fold)
  CindexByFold[[fold]] <- Cindex
}
write.table(CindexByFold, paste0(save_dir, "/Cindex_", folds, "foldCV_", time_var, ".txt"), row.names = F) # save

# Calculate median CI
Cindex_tab <- read.csv(paste0(save_dir, "/Cindex_", folds, "foldCV_", time_var, ".txt"), sep="") # load
Cindex_result <- data.frame(row.names = Cindex_tab$var, "miRNA" = Cindex_tab$var,
                            Cindex_tab[,startsWith(colnames(Cindex_tab), "C_index_test")])
Cindex_result$CI_median <- rowMedians(as.matrix(Cindex_result[,-1]))

# Order individual miRNAs by median keeping signatures first
CIindMiRNAs <- Cindex_result[!grepl("score", rownames(Cindex_result)),] # individual miRNAs
ordered_miRNAs <- c(grep("score", rownames(Cindex_result), value = T),
                    rownames(CIindMiRNAs)[order(CIindMiRNAs$CI_median, decreasing = T)])
Cindex_result <- Cindex_result[ordered_miRNAs,] 
Cindex_result$miRNA <- gsub("_", "-", Cindex_result$miRNA)

write.table(Cindex_result, paste0(save_dir, "/Cindex_result_", folds, "foldCV_", time_var, ".txt"), row.names = F) # save

# Plot
Cindex_result.long <- as.data.frame(melt(setDT(Cindex_result[,!colnames(Cindex_result)%in%"CI_median"]), id.vars = "miRNA", variable.name = "fold")) # long format
Cindex_result.long$fold <- gsub(".*fold", "", Cindex_result.long$fold)
Cindex_result.long$miRNA <- factor(Cindex_result.long$miRNA, levels = unique(Cindex_result.long$miRNA)) # fix order for ggplot
plot_C <- ggplot(Cindex_result.long, aes(x = miRNA, y = value)) +
  geom_boxplot(outlier.size = 0.5, fill = "#878787", alpha = 0.5) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  ggtitle(paste0("C-index ", time_var)) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("C-Index")

while (!is.null(dev.list()))  dev.off()
pdf(paste0(save_dir, "/Cindex_", folds, "foldCV_", time_var, ".pdf"), width = 5, height = 5, family="ArialMT") #export as pdf
plot_C
while (!is.null(dev.list()))  dev.off()

#----

#### Time dependent ROC curve (only signatures) ####

vars <- miRNAs

sdata <- fulldata[fulldata$PatientSet=="Validation",]
time_days <- sdata[,time_var]*365-1 # return to days to better definition
sigs <- grep("score", vars, value = T)

auc_list <- vector(mode='list', length=length(sigs)) # empty list to store AUCs
names(auc_list) <- paste0("AUC_", sigs)
# Calculate an AUC for each time point (each day):
for (i in 1:length(sigs)) {
  for (t in 30:max(time_days)) {
    auc_list[[i]] <- c(auc_list[[i]], survivalROC.C(Stime = time_days, status = sdata$status, marker = sdata[,sigs[i]], predict.time = t, span = 0.05)$AUC)
  }
}


# Plot

while (!is.null(dev.list()))  dev.off()
pdf(paste0(save_dir, "/AUC_", time_var, ".pdf"), width = 5, height = 5) # save pdf

plot(seq(min(sdata[,time_var]), max(sdata[,time_var]), length.out = length(auc_list[[1]])), auc_list[[1]], type = "l", xlab = "Time (years)", ylab = "AUC(t)", col = "green",
     xlim = c(0, 6), ylim = c(0.35, 0.75), 
     main = time_var)
lines(seq(min(sdata[,time_var]), max(sdata[,time_var]), length.out = length(auc_list[[2]])), auc_list[[2]], col = "red")
lines(seq(min(sdata[,time_var]), max(sdata[,time_var]), length.out = length(auc_list[[3]])), auc_list[[3]], col = "blue")
lines(seq(min(sdata[,time_var]), max(sdata[,time_var]), length.out = length(auc_list[[4]])), auc_list[[4]], col = "orange")
lines(seq(min(sdata[,time_var]), max(sdata[,time_var]), length.out = length(auc_list[[5]])), auc_list[[5]], col = "yellow")
lines(seq(min(sdata[,time_var]), max(sdata[,time_var]), length.out = length(auc_list[[6]])), auc_list[[6]], col = "grey")
legend(x = "bottomright", legend = gsub("_score", "", sigs),
       col = c("green", "red", "blue", "orange", "yellow", "grey"), 
       lty = c(1, 1, 1, 1, 1, 1), 
       bty = "n",
       cex = 0.7)

while (!is.null(dev.list()))  dev.off()

#----

