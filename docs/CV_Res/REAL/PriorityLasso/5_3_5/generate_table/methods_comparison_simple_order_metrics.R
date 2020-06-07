# load the results
results <- readRDS("../data/results_different_blocks_2020_05_28.Rds")
results_ddspls <- readRDS("../data/results_different_blocks_ddsPLS_weight_2020_05_20.Rds")
results_different_block_combination_1 <- readRDS("../data/results_permute_blocks_small_1_2020_05_28.Rds")
results_different_block_combination_2 <- readRDS("../data/results_permute_blocks_small_2_2020_05_28.Rds")

true_values_list <- lapply(results, function(x) x$test_y)
true_values <- do.call("c", true_values_list)

true_values_list <- lapply(results_different_block_combination_1, function(x) x$test_y)
true_values_block_comb <- do.call("c", true_values_list)

generate_metrics <- function(object = results,
                             object_2 = results_ddspls,
                             object_3 = results_different_block_combination_1,
                             object_4 = results_different_block_combination_2) {
  
  
  # get the indices depending on the method
  # blocks 1-6 and blocks 1-3
  indices <- c(21:24, 9:12)
  
  # ddsPLS
  predicted_values_ddspls_list <- lapply(object_2, function(x) {
    as.vector(x$pred_value_ddspls[, "1"])
  })
  
  # make a list with the predictions of all methods
  predicted_list <- list()
  predicted_list[[1]] <- do.call("c", predicted_values_ddspls_list)
  
  # blocks 1 - 6
  for (i in 1:length(indices)) {
    predicted_values_list <- lapply(object, function(x) {
      x[["pred_value_list"]][[indices[i]]]
      
    })
    predicted_list[[i + 1]] <- do.call("c", predicted_values_list)
  }
  
  # new indices for predictions with 4, 2, 1, 3
  indices <- 13:16
  
  # blocks 4, 2, 1, 3, 5 (but only predictions with 4, 2, 1, 3)
  start_list <- length(predicted_list)
  for (i in 1:length(indices)) {
    predicted_values_list <- lapply(object_3, function(x) {
      x[["pred_value_list"]][[indices[i]]]
      
    })
    predicted_list[[i + start_list]] <- do.call("c", predicted_values_list)
  }
  
  # blocks 4, 2, 1, 3, 6 (but only predictions with 4, 2, 1, 3)
  start_list <- length(predicted_list)
  for (i in 1:length(indices)) {
    predicted_values_list <- lapply(object_4, function(x) {
      x[["pred_value_list"]][[indices[i]]]
      
    })
    predicted_list[[i + start_list]] <- do.call("c", predicted_values_list)
  }
  
  # calculate the metrics
  metrics_list <- lapply(seq_along(predicted_list), function(i) {
    x <- predicted_list[[i]]
    # do a prediction with the cut-off >0.5 => 1
    x_pred <- ifelse(x > 0.5, 1, 0)
    
    # use the correct true values (either for blocks 1-6 or for blocks 4, 2, 1, 3)
    if (i < 10) {
    confmat <- caret::confusionMatrix(data      = as.factor(x_pred), 
                                      reference = as.factor(true_values))
    } else {
      confmat <- caret::confusionMatrix(data      = as.factor(x_pred), 
                                        reference = as.factor(true_values_block_comb))
    }
    
    mcc <- mcc_metric(conf_matrix = confmat)
    bal_acc <- confmat$byClass["Balanced Accuracy"]
    F_1 <- confmat$byClass["F1"]
    
    list(mcc = mcc,
         bal_acc = bal_acc,
         F_1 = F_1)
  })
  
  # return the metrics list
  metrics_list
  
}

mcc_metric <- function(conf_matrix) {
  "Calculate the MCC Metric [Matthews correlation coefficient]!
   Works only for binary cases! 
   If the Conf_Matrix has more than 2 classes it will return NA 
   instead of the MCC!
       
   Definition of the Metric:
      MCC takes into account true and false positives and negatives and is 
      generally regarded as a balanced measure which can be used even if the 
      classes are of very different sizes.
       
    Args: 
      - conf_matrix (confusionMatrix) : Confusion Matrix created with the 
                                        'caret'-Package!
    Return: 
      - Matthews correlation coefficient (numeric): '1' is best & '-1' is worst
  "
  # [0] Check Inputs -----------------------------------------------------------
  # 0-1 check amount of classes:
  if (nrow(conf_matrix$table) != 2) {
    warning("Can not calc the MCC-Metric! Return NULL")
    return(NA)
  }
  
  # 0-2 Check class of conf_matrix
  if (class(conf_matrix) != "confusionMatrix") {
    stop("conf_matrix not of class 'confusionMatrix'")
  }
  
  # [1] Calc the Score ---------------------------------------------------------
  # 1-1 Get the TruePositives, FalsePositives, FalseNegatives & TrueNegatives
  TP <- conf_matrix$table[1,1]
  TN <- conf_matrix$table[2,2]
  FP <- conf_matrix$table[1,2]
  FN <- conf_matrix$table[2,1]
  
  # 1-2 Calc single parts of MMC [bases on TP, TN, FP, FN]
  mcc_num <- (TP * TN - FP * FN)
  mcc_den <- as.double((TP + FP)) * as.double((TP + FN)) * as.double((TN + FP)) * as.double((TN + FN))
  
  # 1-3 Put together the single parts!
  mcc_final <- mcc_num/sqrt(mcc_den)
  
  # [2] Return  ----------------------------------------------------------------
  return(mcc_final)
}

metrics_pl_ddspls <- generate_metrics()

# load the results from Frederik
load("../data/SingleBlockApproach.Rda")
load("../data/BlockWise_Approach.Rda")
load("../data/FoldWise_Approach.Rda")

# for the comparison, use the F1 score

# Frederik's results are per fold of the CV, so use the mean of all folds
# from block 1
f1_single_block <- mean(unlist(lapply(single_block_res$df1, function(x) {
  x$F1
})))

# block-wise use the F1 weighting
f1_block_wise <- mean(unlist(lapply(blockwise_res, function(x) {
  x$f1_weight$F1
})))

# fold-wise
f1_fold_wise <- mean(unlist(lapply(foldwise_res, function(x) {
  x$f1_weight$F1
})))

# create the overview table
overview_metric <- data.frame(ign_zero_16 = metrics_pl_ddspls[[2]]$F_1,
                              ign_incp_16 = metrics_pl_ddspls[[3]]$F_1,
                              imp_maxb_16 = metrics_pl_ddspls[[4]]$F_1,
                              imp_maxn_16 = metrics_pl_ddspls[[5]]$F_1,
                              ign_zero_13 = metrics_pl_ddspls[[6]]$F_1,
                              ign_incp_13 = metrics_pl_ddspls[[7]]$F_1,
                              imp_maxb_13 = metrics_pl_ddspls[[8]]$F_1,
                              imp_maxn_13 = metrics_pl_ddspls[[9]]$F_1,
                              ign_zero_43_5 = metrics_pl_ddspls[[10]]$F_1,
                              ign_incp_43_5 = metrics_pl_ddspls[[11]]$F_1,
                              imp_maxb_43_5 = metrics_pl_ddspls[[12]]$F_1,
                              imp_maxn_43_5 = metrics_pl_ddspls[[13]]$F_1,
                              ign_zero_43_6 = metrics_pl_ddspls[[14]]$F_1,
                              ign_incp_43_6 = metrics_pl_ddspls[[15]]$F_1,
                              imp_maxb_43_6 = metrics_pl_ddspls[[16]]$F_1,
                              imp_maxn_43_6 = metrics_pl_ddspls[[17]]$F_1,
                              ddspls = metrics_pl_ddspls[[1]]$F_1,
                              singl_bl1 = f1_single_block,
                              blockwise = f1_block_wise,
                              foldwise = f1_fold_wise)

overview_metric_final_table <- data.frame(blocks_16_13_other = 
                                            unlist(overview_metric[1, c(1:8, 17:20)]),
                                          blocks_43 = 
                                            c(unlist(overview_metric[1, c(9:16)]), rep(NA, 4)))

# print the table as latex code for copy & paste into the latex file
library(xtable)
print(xtable(overview_metric_final_table))
