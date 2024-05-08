# Utility functions

df.reduce_by_id <- function(df, 
                            id_col = 1, 
                            collapse_by = ":")
{
  #'Collapse duplicates instead of removing
  #'
  #'@description 
  #'Reduce data frame eliminating duplicated elements based on one column but 
  #'collapsing all different values from other columns into one row.
  #'
  #'@param df Data frame
  #'@param id_col Column of data frame to check for duplicates. Default os the first column.
  #'@param collapse_by Character to separate collapsed items. Default ":"
  #'
  
  dup_ids <- sum(duplicated(df[,id_col]))
  for (i in 1:length(unique(df[,id_col]))) {
    for (j in 1:ncol(df)) {
      if (length(unique(df[df[,id_col]==unique(df[,id_col])[i],j]))!=1) {
        df[df[,id_col]==unique(df[,id_col])[i],j] <- 
          paste(na.omit(df[df[,id_col]==
                             unique(df[,id_col])[i],j])
                , collapse = collapse_by)
      }
    }
  }
  df <- df[!duplicated.data.frame(df),]
  message(paste0(dup_ids, " rows removed. Different values collapsed."))
  return(df)
}


perform_uni_multi_coxph <- function(data, 
                                    time_var, 
                                    status_var, 
                                    predictors, 
                                    covars, 
                                    save_table = F, 
                                    save_table_file = NULL)
{
  #'Perform univariate and multivariate coxph analysis for multiple predictors
  #'
  #'@description 
  #'Fit univariate and multivariate cox proportional hazards models for each 
  #'predictor. A forest plot of resulting metrics is returned and a result table
  #' is optionally saved
  #'
  #'
  #'@param data Data frame
  #'@param time_var Time column of data frame.
  #'@param status_var Status column of data frame.
  #'@param predictors Vector of predictors to analyze.
  #'@param covars Vector of covariates to include in multivariate analysis. If 
  #'empty, only univariate analysis will be performed.
  #'@param save_table Logical. Save result table.
  #'@param save_table_file File name to save result table.
  #'
  
  
  
  
  ### Fit coxph
  
  surv <- Surv(time = data[,time_var], event = data[,status_var])
  cox <- data.frame("predictors"=predictors,
                    "beta_uni"=NA,
                    "HR_uni"=NA,
                    "HR.CI_lower_uni"=NA,
                    "HR.CI_upper_uni"=NA,
                    "wald.test_uni"=NA,
                    "p.value_uni"=NA,
                    "beta_multi"=NA,
                    "HR_multi"=NA,
                    "HR.CI_lower_multi"=NA,
                    "HR.CI_upper_multi"=NA,
                    "wald.test_multi"=NA,
                    "p.value_multi"=NA)
  for (i in 1:length(predictors)) {
    # Univariate cox analysis:
    cox_reg_model_uni <- coxph(formula = surv ~ data[,predictors[i]], data = data)
    # Multivariate cox analysis:
    cox_reg_model_multi <- coxph(formula = as.formula(paste("surv ~ data[,predictors[i]] +", paste(covars, collapse = " + "))),
                                 data = data) 
    
    unimod <- summary(cox_reg_model_uni)
    cox$beta_uni[i] <- signif(unimod$coef[1], digits=2);#coeficient beta
    cox$HR_uni[i] <- signif(unimod$coef[2], digits=2);#exp(beta)
    cox$HR.CI_lower_uni[i] <- signif(unimod$conf.int[,"lower .95"], 2)
    cox$HR.CI_upper_uni[i] <- signif(unimod$conf.int[,"upper .95"],2)
    cox$wald.test_uni[i] <- signif(unimod$wald["test"], digits=2)
    cox$p.value_uni[i] <- signif(unimod$coef[,"Pr(>|z|)"], digits=2)
    
    multimod <- summary(cox_reg_model_multi)
    cox$beta_multi[i] <- signif(multimod$coef[1,1], digits=2);#coeficient beta
    cox$HR_multi[i] <- signif(multimod$coef[1,2], digits=2);#exp(beta)
    cox$HR.CI_lower_multi[i] <- signif(multimod$conf.int[1,"lower .95"], 2)
    cox$HR.CI_upper_multi[i] <- signif(multimod$conf.int[1,"upper .95"],2)
    cox$wald.test_multi[i] <- signif(multimod$wald["test"], digits=2)
    cox$p.value_multi[i] <- signif(multimod$coef[1,"Pr(>|z|)"], digits=2)
  }
  
  #cox[,unlist(lapply(cox, is.numeric))] <- round(cox[,unlist(lapply(cox, is.numeric))], 3) #round numeric cols
  cox
  
  # Coef exlanation: Positive coef means that an increase in the predictor variable is associated with an increased hazard of the event
  # Positive coef: more expression more risk. 
  # Negative coef: less expression more risk.
  
  # Save
  if (save_table) {
    write.table(cox, file = save_table_file, sep = "\t", row.names = F)
  }
  
  
  
  ### Plot
  
  library(forestploter)
  library(grid)
  
  cox.fig <- data.frame(gsub("_", "-", cox$predictors), # Sustituimos _ solo para aplicarlo a miRNAs !!!
                        paste(rep(" ", 20), collapse = " "),
                        paste0(cox$HR_uni, " [", 
                               cox$HR.CI_lower_uni, "-", 
                               cox$HR.CI_upper_uni, "]"),
                        cox$p.value_uni,
                        "  ",
                        paste(rep(" ", 20), collapse = " "),
                        paste0(cox$HR_multi, " [", 
                               cox$HR.CI_lower_multi, "-", 
                               cox$HR.CI_upper_multi, "]"),
                        cox$p.value_multi)
  
  cox.fig[,4][cox.fig[,4]>=0.001] <- round(cox.fig[,4][cox.fig[,4]>=0.001], 3)
  cox.fig[,4][cox.fig[,4]<0.001] <- "<0.001"
  cox.fig[,8][cox.fig[,8]>=0.001] <- round(cox.fig[,8][cox.fig[,8]>=0.001], 3)
  cox.fig[,8][cox.fig[,8]<0.001] <- "<0.001"
  
  
  colnames(cox.fig) <- c("", "", "HR [95% CI]", "p-value", "", "", "HR [95% CI]", "p-value")
  
  p <- forest(cox.fig,
              est = list(cox$HR_uni, cox$HR_multi),
              lower = list(cox$HR.CI_lower_uni, cox$HR.CI_lower_multi), 
              upper = list(cox$HR.CI_upper_uni, cox$HR.CI_upper_multi),
              ci_column = c(2, 6),
              ref_line = 1,
              xlim = c(0, max(cox$HR.CI_upper_uni)+0.5),
              ticks_at = c(0, 0.5, 1, 2),
              xlab = "Hazard Ratio",
              theme = forest_theme(core = list(bg_params=list(fill = c("white"))),
                                   arrow_label_just = "end",
                                   arrow_type = "closed"),
              nudge_y = 0.2)
  p <- insert_text(p, text = "Univariate",
                   part = "header",
                   row = 1,
                   col = 2:4,
                   gp = gpar(fontface = "bold"))
  p <- add_text(p, text = "Multivariate",
                part = "header",
                row = 1,
                col = 6:8,
                gp = gpar(fontface = "bold"))
  p <- add_border(p, 
                  part = "header", 
                  row = 1,
                  col = 2:4,
                  gp = gpar(lwd = 1))
  p <- add_border(p, 
                  part = "header", 
                  row = 1,
                  col = 6:8,
                  gp = gpar(lwd = 1))
  p <- add_border(p, 
                  part = "header", 
                  row = 2,
                  gp = gpar(lwd = 1))
  p <- edit_plot(p, col = c(3,4,7,8), 
                 which = "text",
                 hjust = unit(0.5, "npc"),
                 x = unit(0.5, "npc"))
  p <- insert_text(p, text = time_var,
                   part = "header",
                   row = 1,
                   col = 2:8,
                   gp = gpar(fontface = "bold"))
  p <- add_border(p, 
                  part = "header", 
                  row = 1,
                  col = 2:8,
                  gp = gpar(lwd = 1))
  return(p)
}


get_coxph_C <- function(dataset,
                        var, 
                        status_var, 
                        time_var, 
                        train, 
                        test)
{ 
  
  #' Fits a Cox model and calculate the Harrell's C Index
  #' 
  #' @description 
  #' Fits a Cox model for survival data and calculate the Harrell's Concordance Index for the train and test.
  #' C index of test data indicates the predictive power of model.
  #' https://cran.r-project.org/web/packages/SurvMetrics/vignettes/SurvMetrics-vignette.html
  #' 
  #' @param dataset Data frame with patients in rows and variables in columns: time, status and variables (clinical, gene expression, etc)
  #' @param var One or more variables to use in the Cox model.
  #' @param status_var Status variable: Dead or alive for survival data
  #' @param time_var Time variable: Survival time. Can be Overall Survival of Event Free Survival, for example.
  #' @param train Training vector. Can be a TRUE/FALSE vector or a vector of patient names.
  #' @param test Test vector. Can be a TRUE/FALSE vector or a vector of patient names.
  #' 
  
  data <- data.frame(time = dataset[,time_var],
                     status = dataset[,status_var],
                     dataset[,var],
                     row.names = rownames(dataset))
  data_train <- data[train,]
  data_test <- data[test,]
  
  # Train:
  cox_reg_model <- coxph(Surv(time, status)~., data = data_train, x = TRUE) # fit cox model in train
  C_index_train <- concordance(cox_reg_model)$concordance # or summary(cox_reg_model)$concordance 
  dis_time = sort(unique(data_train$time[data_train$status==1])) # times of death
  # Test:
  library(pec)
  mat_cox = predictSurvProb(cox_reg_model, data_test, dis_time#c(25:50) #
  ) # survival probability predicted for each test patient in each time
  med_index = median(1:length(dis_time)) # Get the median time point
  surv_obj = Surv(time = data_test$time, event = data_test$status) # surv object of test
  #C index for Cox:
  library(SurvMetrics)
  C_index_test = Cindex(surv_obj, predicted = mat_cox[, med_index]) # Calculate the Concordance index between the survival time and the predicted result at the median time point.
  res <- data.frame(C_index_train = C_index_train,
                    C_index_test = C_index_test)
  return(res)
}










