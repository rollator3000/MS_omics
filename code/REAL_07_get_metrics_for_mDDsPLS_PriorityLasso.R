"Script to calculate the Metrics of the Priority-Lasso / mdd-sPLS
  - only recived the predictions of the models!
  Save the metrics as files and use them in an other script to create plots!
"
# [0] Librarys and Functions
library(checkmate)

mcc_metric      <- function(conf_matrix) {
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
get_metrics     <- function(preds, truth) {
  "Function to calculate the metrics based on 2 lists 'preds' 'truth'.
   Both lists need the same length & each of the list entrances needs the same
   length aswell!
   
   Args:
    - preds (list) : list of length k, with predicted probabilities!
    - truth (list) : list of length k, with true class labels!  
  
   Return:
    list of length k, where each entrance holds the metrics for 
    preds[[i]] & truth[[i]]
  "
  # [0] Check Input
  # 0-1 Arguments are lists
  assert_list(preds)
  assert_list(truth)
  
  # 0-2 Same length?
  if (length(preds) != length(truth)) {
    stop("'preds' & 'truth' have different lengths")
  }
  
  # 0-3 Each list entrance of the same length
  same_len <- sapply(1:length(preds), FUN = function(i_) {
    length(preds[[i_]]) != length(truth[[i_]])
  })
  if (any(same_len)) {
    stop("entrances of 'preds' & 'truth' have a different length!")
  }
  
  # [1] Calculate the Confusion Matrix
  # 1-1 Loop over each list entrance and compare the predicitions with the truth!
  metrics <- sapply(1:length(preds), FUN = function(i_){
    
    # Convert probabilities to class predictions
    predicted_classes <- ifelse(preds[[i_]] > 0.5,  1, 0)
    
    # Get the true responses
    true_responses <- truth[[i_]]
    
    # Get the Confusion Matrix
    confmat <- caret::confusionMatrix(as.factor(predicted_classes), 
                                      as.factor(true_responses))
    
    # Area under the ROC Curve
    roc1 <- pROC::auc(pROC::roc(as.numeric(predicted_classes), 
                                as.numeric(as.character(true_responses)),
                                levels = levels(as.factor(as.character(true_responses))),
                                direction = "<"))
    
    roc2 <- pROC::auc(pROC::roc(as.numeric(predicted_classes), 
                                as.numeric(as.character(true_responses)),
                                levels = levels(as.factor(as.character(true_responses))),
                                direction = ">"))
    
    # MCC Matthews correlation coefficient [only for binary cases!]
    mcc <- mcc_metric(conf_matrix = confmat)
    
    # Collect all metrics in a list & replace the not defined values
    res <- list("Accuracy"    = confmat$overall["Accuracy"],
                "Kappa"       = confmat$overall["Kappa"],
                "Sensitifity" = confmat$byClass["Sensitivity"],
                "Specificity" = confmat$byClass["Specificity"],
                "Precision"   = confmat$byClass["Precision"],
                "Recall"      = confmat$byClass["Recall"],
                "F1"          = confmat$byClass["F1"],
                "Balance_Acc" = confmat$byClass["Balanced Accuracy"],
                "Pos_Pred_Value" =  confmat$byClass["Pos Pred Value"],
                "Neg_Pred_Value" =  confmat$byClass["Neg Pred Value"],
                "Prevalence"  = confmat$byClass["Prevalence"],      
                "AUC1"        = as.numeric(roc1),
                "AUC2"        = as.numeric(roc2),
                "MCC"         = mcc)
    
    #   If the F1-Score/ Precision/ Recall is NA, then we set it to 0 [-1 for MCC]! 
    if (is.na(res$F1))             res$F1             <- 0
    if (is.na(res$Precision))      res$Precision      <- 0
    if (is.na(res$Recall))         res$Recall         <- 0
    if (is.na(res$MCC))            res$MCC            <- -1
    if (is.na(res$Pos_Pred_Value)) res$Pos_Pred_Value <- 0
    if (is.na(res$Neg_Pred_Value)) res$Neg_Pred_Value <- 0
    
    res
  })
  
  return(metrics)
}

# [1] Section 5.3.1    ---------------------------------------------------------
# 1-1 Load the CV results ['Priority Lasso' & 'mdd-sPLS']
res_PL <- readRDS("./docs/CV_Res/REAL/PriorityLasso/5_3_1/data/results_different_blocks_2020_05_28.Rds")
res_MD <- readRDS("./docs/CV_Res/REAL/PriorityLasso/5_3_1/data/results_different_blocks_ddsPLS_weight_2020_05_20.Rds")

# 1-2 Get the true responses
true_values_list <- lapply(res_PL, function(x) x$test_y)
true_values      <- do.call("c", true_values_list)

# 1-3 ---------- Get Metrics for mdd-sPLS 
# 1-3-1 MD Predictions
predicted_values_ddspls_list <- lapply(res_MD, function(x) {
  as.vector(x$pred_value_ddspls[, "1"])
})

# 1-3-2 Metrics
MD_SPLSS <- get_metrics(predicted_values_ddspls_list, true_values_list)

# 1-4  ---------- Get the Metrics for the different PL methods!
# 1-4-1 Indices to read out results [given by Hagenberg]
indices <- 21:24

# 1-4-2 List to save the results
list_names      <- c("ignore, zero", "ignore, intercept",
                     "impute maximise blocks", "impute, maximise n")
all_res         <- vector("list", length(list_names))
names(all_res)  <- list_names

# 1-4-2 Loop over the indices and get the metrics
j <- 1
for (i_ in indices) {
  
  # Extract predicitons
  predicted_values_list <- lapply(res_PL, function(x) {
    x[["pred_value_list"]][[i_]]
  })
  
  # Calculate the metric
  curr_metric <- get_metrics(predicted_values_list, true_values_list)
  
  all_res[[list_names[j]]] <- curr_metric
  
  # Count up index
  j <- j + 1
}

# 1-5 Add the mdd-sPLS results to the 'all_res' list 
all_res[["mdd_sPLS"]] <- MD_SPLSS

# 1-6 Save the list with the results of all metrics
save(all_res, file = "./docs/CV_Res/REAL/Hagenberg_5_3_1.RData")

# [2] Section 5.3.2    ---------------------------------------------------------
# 2-1 Load the CV results ['Priority Lasso' & 'mdd-sPLS']
res_PL <- readRDS("./docs/CV_Res/REAL/PriorityLasso/5_3_2/data/results_different_blocks_2020_05_28.Rds")
res_MD <- readRDS("./docs/CV_Res/REAL/PriorityLasso/5_3_2/data/results_different_blocks_ddsPLS_weight_2020_05_20.Rds")

# 2-2 Get the true responses
true_values_list <- lapply(res_PL, function(x) x$test_y)
true_values      <- do.call("c", true_values_list)

# 2-3 ---------- Get Metrics for mdd-sPLS 
# 2-3-1 MD Predictions
predicted_values_ddspls_list <- lapply(res_MD, function(x) {
  as.vector(x$pred_value_ddspls[, "1"])
})

# 2-3-2 Metrics
MD_SPLSS <- get_metrics(predicted_values_ddspls_list, true_values_list)

# 2-4  ---------- Get the Metrics for the different PL methods with ALL BLOCKS!
# 2-4-1 List to save the results
list_names      <- c("ignore, zero", "ignore, intercept",
                     "impute maximise blocks", "impute, maximise n")
all_res         <- vector("list", length(list_names))
names(all_res)  <- list_names

# 2-4-2 Loop over the different PL Methods & get the metrics
for (curr_method in 1:4) {
  
  # Get Indices based on current method
  indices <- seq(from = curr_method, to = 20 + curr_method, by = 4)
  
  # Get the name of the current method
  curr_method_name <- names(all_res)[curr_method]
  
  # [1] ----- Predicitions w/ one Block! 
  #     --> Calculate the Metrics & add to 'all_res'
  curr_block1_pred <- lapply(res_PL, function(x) {
    x[["pred_value_list"]][[indices[1]]]
  })
  
  all_res[[curr_method_name]]$block_1 <- get_metrics(preds = curr_block1_pred,
                                                     truth = true_values_list)
  
  # [2] ----- Predicitions w/ two Blocks! 
  #     --> Calculate the Metrics & add to 'all_res'
  curr_block2_pred <- lapply(res_PL, function(x) {
    x[["pred_value_list"]][[indices[2]]]
  })
  
  all_res[[curr_method_name]]$block_1_2 <- get_metrics(preds = curr_block2_pred,
                                                       truth = true_values_list)
  
  # [3] ----- Predicitions w/ three Blocks! 
  #     --> Calculate the Metrics & add to 'all_res'
  curr_block3_pred <- lapply(res_PL, function(x) {
    x[["pred_value_list"]][[indices[3]]]
  })
  
  all_res[[curr_method_name]]$block_1_2_3 <- get_metrics(preds = curr_block3_pred,
                                                         truth = true_values_list)
  
  # [4] ----- Predicitions w/ four Blocks! 
  #     --> Calculate the Metrics & add to 'all_res'
  curr_block4_pred <- lapply(res_PL, function(x) {
    x[["pred_value_list"]][[indices[4]]]
  })
  
  all_res[[curr_method_name]]$block_1_2_3_4 <- get_metrics(preds = curr_block4_pred,
                                                           truth = true_values_list)
  
  # [5] ----- Predicitions w/ five Blocks! 
  #     --> Calculate the Metrics & add to 'all_res'
  curr_block5_pred <- lapply(res_PL, function(x) {
    x[["pred_value_list"]][[indices[5]]]
  })
  
  all_res[[curr_method_name]]$block_1_2_3_4_5 <- get_metrics(preds = curr_block5_pred,
                                                           truth = true_values_list)
  
  # [6] ----- Predicitions w/ all (six) Blocks! 
  #     --> Calculate the Metrics & add to 'all_res'
  curr_block6_pred <- lapply(res_PL, function(x) {
    x[["pred_value_list"]][[indices[6]]]
  })
  
  all_res[[curr_method_name]]$block_1_2_3_4_5_6 <- get_metrics(preds = curr_block6_pred,
                                                               truth = true_values_list)
}

# 2-5  ---------- Get the Metrics for the different PL methods with PRED BLOCKS!
# 2-5-1 List to save the results
list_names       <- c("ignore, zero", "ignore, intercept",
                      "impute maximise blocks", "impute, maximise n")
all_res2         <- vector("list", length(list_names))
names(all_res2)  <- list_names

# 2-5-2 Loop over the different PL Methods & get the metrics
for (curr_method in 1:4) {
  
  # Get Indices based on current method
  indices <- seq(from = curr_method, to = 20 + curr_method, by = 4)
  
  # Get the name of the current method
  curr_method_name <- names(all_res2)[curr_method]
  
  # [1] ----- Predicitions w/ one Block! 
  #     --> Calculate the Metrics & add to 'all_res'
  curr_block1_pred <- lapply(res_PL, function(x) {
    x[["pred_value_single"]][[indices[1]]]
  })
  
  all_res2[[curr_method_name]]$block_1 <- get_metrics(preds = curr_block1_pred,
                                                      truth = true_values_list)
  
  # [2] ----- Predicitions w/ two Blocks! 
  #     --> Calculate the Metrics & add to 'all_res'
  curr_block2_pred <- lapply(res_PL, function(x) {
    x[["pred_value_single"]][[indices[2]]]
  })
  
  all_res2[[curr_method_name]]$block_1_2 <- get_metrics(preds = curr_block2_pred,
                                                        truth = true_values_list)
  
  # [3] ----- Predicitions w/ three Blocks! 
  #     --> Calculate the Metrics & add to 'all_res'
  curr_block3_pred <- lapply(res_PL, function(x) {
    x[["pred_value_single"]][[indices[3]]]
  })
  
  all_res2[[curr_method_name]]$block_1_2_3 <- get_metrics(preds = curr_block3_pred,
                                                          truth = true_values_list)
  
  # [4] ----- Predicitions w/ four Blocks! 
  #     --> Calculate the Metrics & add to 'all_res'
  curr_block4_pred <- lapply(res_PL, function(x) {
    x[["pred_value_single"]][[indices[4]]]
  })
  
  all_res2[[curr_method_name]]$block_1_2_3_4 <- get_metrics(preds = curr_block4_pred,
                                                            truth = true_values_list)
  
  # [5] ----- Predicitions w/ five Blocks! 
  #     --> Calculate the Metrics & add to 'all_res'
  curr_block5_pred <- lapply(res_PL, function(x) {
    x[["pred_value_single"]][[indices[5]]]
  })
  
  all_res2[[curr_method_name]]$block_1_2_3_4_5 <- get_metrics(preds = curr_block5_pred,
                                                              truth = true_values_list)
  
  # [6] ----- Predicitions w/ all (six) Blocks! 
  #           For all blocks, "pred_value_single" doesn't have its own model
  #     --> Calculate the Metrics & add to 'all_res'
  curr_block6_pred <- lapply(res_PL, function(x) {
    x[["pred_value_list"]][[indices[6]]]
  })
  
  all_res2[[curr_method_name]]$block_1_2_3_4_5_6 <- get_metrics(preds = curr_block6_pred,
                                                                truth = true_values_list)
}


# 2-6 Paste the 2 lists 'all_res' & 'all_res2'  to a single list!
all_res_532 <- list("all_block"   = all_res,
                    "pred_blocks" = all_res2,
                    "mdd_splss"   = MD_SPLSS)
save(all_res_532, file = "./docs/CV_Res/REAL/Hagenberg_5_3_2.RData")

# [3] Section 5.3.4    ---------------------------------------------------------
# Setting 1 [4, 2, 1, 3, 5]  ---------------------------------------------------
# 3-1 Load the CV results ['Priority Lasso' & 'mdd-sPLS']
res_PL <- readRDS("./docs/CV_Res/REAL/PriorityLasso/5_3_4/data/results_permute_blocks_small_1_2020_05_28.Rds")

# 3-2 Get the true responses
true_values_list <- lapply(res_PL, function(x) x$test_y)
true_values      <- do.call("c", true_values_list)

# 3-3 Get the indices (all combinations from block 1-5)
indices <- 1:20

# 3-4 loop over the remaining indices and get PL results
# 3-4-1 Initalize list to save results
list_names      <- c("ignore, zero", "ignore, intercept",
                     "impute maximise blocks", "impute, maximise n")
all_res         <- vector("list", length(list_names))
names(all_res)  <- list_names

# 3-4-2 Define the different situations [which blocks used]
all_sits <- c("block_4", "block_4_2", "block_4_2_1", 
              "block_4_2_1_3", "block_4_2_1_3_5")

# 3-4-3 Calculate the metrics for the different methods and blocks for predictions
# 3-4-3-1 IGNORE ZERO
j <- 1
for (i in seq(from = 1, to = 20, by = 4)) {
  
  # Get current predictions
  curr_predicted_values_list <- lapply(res_PL, function(x) {
    x[["pred_value_list"]][[i]]
  })
  
  # Calculate metric
  curr_metric <- get_metrics(curr_predicted_values_list,
                             true_values_list)
  
  # current situatiuon
  curr_sit <- all_sits[j]
  
  # Add the result to the list
  all_res$`ignore, zero`[[curr_sit]] <- curr_metric
  
  # Count up Index
  j <- j + 1
}

# 3-4-3-2 IGNORE INTERCEPT
j <- 1
for (i in seq(from = 2, to = 20, by = 4)) {
  
  # Get current predictions
  curr_predicted_values_list <- lapply(res_PL, function(x) {
    x[["pred_value_list"]][[i]]
  })
  
  # Calculate metric
  curr_metric <- get_metrics(curr_predicted_values_list,
                             true_values_list)
  
  # current situatiuon
  curr_sit <- all_sits[j]
  
  # Add the result to the list
  all_res$`ignore, intercept`[[curr_sit]] <- curr_metric
  
  # Count up Index
  j <- j + 1
}

# 3-4-3-3 IMPUTE MAXIMISE BLOCKS
j <- 1
for (i in seq(from = 3, to = 20, by = 4)) {
  
  # Get current predictions
  curr_predicted_values_list <- lapply(res_PL, function(x) {
    x[["pred_value_list"]][[i]]
  })
  
  # Calculate metric
  curr_metric <- get_metrics(curr_predicted_values_list,
                             true_values_list)
  
  # current situatiuon
  curr_sit <- all_sits[j]
  
  # Add the result to the list
  all_res$`impute maximise blocks`[[curr_sit]] <- curr_metric
  
  # Count up Index
  j <- j + 1
}

# 3-4-3-3 IMPUTE MAXIMISE N
j <- 1
for (i in seq(from = 3, to = 20, by = 4)) {
  
  # Get current predictions
  curr_predicted_values_list <- lapply(res_PL, function(x) {
    x[["pred_value_list"]][[i]]
  })
  
  # Calculate metric
  curr_metric <- get_metrics(curr_predicted_values_list,
                             true_values_list)
  
  # current situatiuon
  curr_sit <- all_sits[j]
  
  # Add the result to the list
  all_res$`impute, maximise n`[[curr_sit]] <- curr_metric
  
  # Count up Index
  j <- j + 1
}

# 3-4-3-4 MDD_SPLS
predicted_values_ddspls_list <- lapply(res_PL, function(x) {
  as.vector(x$pred_ddspls[, "1"])
})

all_res[["mdd_spls"]] <- get_metrics(predicted_values_ddspls_list, truth = true_values_list)

# 3-5 Save the metrics
save(all_res, file = "./docs/CV_Res/REAL/Hagenberg_5_3_4__Setting1.R")

# Setting 2 [4, 2, 1, 3, 6]  ---------------------------------------------------
# 3-1 Load the CV results ['Priority Lasso' & 'mdd-sPLS']
res_PL <- readRDS("./docs/CV_Res/REAL/PriorityLasso/5_3_4/data/results_permute_blocks_small_2_2020_05_28.Rds")

# 3-2 Get the true responses
true_values_list <- lapply(res_PL, function(x) x$test_y)
true_values      <- do.call("c", true_values_list)

# 3-3 Get the indices (all combinations from block 1-5)
indices <- 1:20

# 3-4 loop over the remaining indices and get PL results
# 3-4-1 Initalize list to save results
list_names      <- c("ignore, zero", "ignore, intercept",
                     "impute maximise blocks", "impute, maximise n")
all_res         <- vector("list", length(list_names))
names(all_res)  <- list_names

# 3-4-2 Define the different situations [which blocks used]
all_sits <- c("block_4", "block_4_2", "block_4_2_1", 
              "block_4_2_1_3", "block_4_2_1_3_6")

# 3-4-3 Calculate the metrics for the different methods and blocks for predictions
# 3-4-3-1 IGNORE ZERO
j <- 1
for (i in seq(from = 1, to = 20, by = 4)) {
  
  # Get current predictions
  curr_predicted_values_list <- lapply(res_PL, function(x) {
    x[["pred_value_list"]][[i]]
  })
  
  # Calculate metric
  curr_metric <- get_metrics(curr_predicted_values_list,
                             true_values_list)
  
  # current situatiuon
  curr_sit <- all_sits[j]
  
  # Add the result to the list
  all_res$`ignore, zero`[[curr_sit]] <- curr_metric
  
  # Count up Index
  j <- j + 1
}

# 3-4-3-2 IGNORE INTERCEPT
j <- 1
for (i in seq(from = 2, to = 20, by = 4)) {
  
  # Get current predictions
  curr_predicted_values_list <- lapply(res_PL, function(x) {
    x[["pred_value_list"]][[i]]
  })
  
  # Calculate metric
  curr_metric <- get_metrics(curr_predicted_values_list,
                             true_values_list)
  
  # current situatiuon
  curr_sit <- all_sits[j]
  
  # Add the result to the list
  all_res$`ignore, intercept`[[curr_sit]] <- curr_metric
  
  # Count up Index
  j <- j + 1
}

# 3-4-3-3 IMPUTE MAXIMISE BLOCKS
j <- 1
for (i in seq(from = 3, to = 20, by = 4)) {
  
  # Get current predictions
  curr_predicted_values_list <- lapply(res_PL, function(x) {
    x[["pred_value_list"]][[i]]
  })
  
  # Calculate metric
  curr_metric <- get_metrics(curr_predicted_values_list,
                             true_values_list)
  
  # current situatiuon
  curr_sit <- all_sits[j]
  
  # Add the result to the list
  all_res$`impute maximise blocks`[[curr_sit]] <- curr_metric
  
  # Count up Index
  j <- j + 1
}

# 3-4-3-3 IMPUTE MAXIMISE N
j <- 1
for (i in seq(from = 3, to = 20, by = 4)) {
  
  # Get current predictions
  curr_predicted_values_list <- lapply(res_PL, function(x) {
    x[["pred_value_list"]][[i]]
  })
  
  # Calculate metric
  curr_metric <- get_metrics(curr_predicted_values_list,
                             true_values_list)
  
  # current situatiuon
  curr_sit <- all_sits[j]
  
  # Add the result to the list
  all_res$`impute, maximise n`[[curr_sit]] <- curr_metric
  
  # Count up Index
  j <- j + 1
}

# 3-4-3-4 MDD_SPLS
predicted_values_ddspls_list <- lapply(res_PL, function(x) {
  as.vector(x$pred_ddspls[, "1"])
})

all_res[["mdd_spls"]] <- get_metrics(predicted_values_ddspls_list, truth = true_values_list)

# 3-5 Save the metrics
save(all_res, file = "./docs/CV_Res/REAL/Hagenberg_5_3_4__Setting2.R")

# Setting 3 [4, 2, 1, 5, 3]  ---------------------------------------------------
# 3-1 Load the CV results ['Priority Lasso' & 'mdd-sPLS']
res_PL <- readRDS("./docs/CV_Res/REAL/PriorityLasso/5_3_4/data/results_permute_blocks_small_3_2020_05_28.Rds")

# 3-2 Get the true responses
true_values_list <- lapply(res_PL, function(x) x$test_y)
true_values      <- do.call("c", true_values_list)

# 3-3 Get the indices (all combinations from block 1-5)
indices <- 1:20

# 3-4 loop over the remaining indices and get PL results
# 3-4-1 Initalize list to save results
list_names      <- c("ignore, zero", "ignore, intercept",
                     "impute maximise blocks", "impute, maximise n")
all_res         <- vector("list", length(list_names))
names(all_res)  <- list_names

# 3-4-2 Define the different situations [which blocks used]
all_sits <- c("block_4", "block_4_2", "block_4_2_1", 
              "block_4_2_1_5", "block_4_2_1_5_3")

# 3-4-3 Calculate the metrics for the different methods and blocks for predictions
# 3-4-3-1 IGNORE ZERO
j <- 1
for (i in seq(from = 1, to = 20, by = 4)) {
  
  # Get current predictions
  curr_predicted_values_list <- lapply(res_PL, function(x) {
    x[["pred_value_list"]][[i]]
  })
  
  # Calculate metric
  curr_metric <- get_metrics(curr_predicted_values_list,
                             true_values_list)
  
  # current situatiuon
  curr_sit <- all_sits[j]
  
  # Add the result to the list
  all_res$`ignore, zero`[[curr_sit]] <- curr_metric
  
  # Count up Index
  j <- j + 1
}

# 3-4-3-2 IGNORE INTERCEPT
j <- 1
for (i in seq(from = 2, to = 20, by = 4)) {
  
  # Get current predictions
  curr_predicted_values_list <- lapply(res_PL, function(x) {
    x[["pred_value_list"]][[i]]
  })
  
  # Calculate metric
  curr_metric <- get_metrics(curr_predicted_values_list,
                             true_values_list)
  
  # current situatiuon
  curr_sit <- all_sits[j]
  
  # Add the result to the list
  all_res$`ignore, intercept`[[curr_sit]] <- curr_metric
  
  # Count up Index
  j <- j + 1
}

# 3-4-3-3 IMPUTE MAXIMISE BLOCKS
j <- 1
for (i in seq(from = 3, to = 20, by = 4)) {
  
  # Get current predictions
  curr_predicted_values_list <- lapply(res_PL, function(x) {
    x[["pred_value_list"]][[i]]
  })
  
  # Calculate metric
  curr_metric <- get_metrics(curr_predicted_values_list,
                             true_values_list)
  
  # current situatiuon
  curr_sit <- all_sits[j]
  
  # Add the result to the list
  all_res$`impute maximise blocks`[[curr_sit]] <- curr_metric
  
  # Count up Index
  j <- j + 1
}

# 3-4-3-3 IMPUTE MAXIMISE N
j <- 1
for (i in seq(from = 3, to = 20, by = 4)) {
  
  # Get current predictions
  curr_predicted_values_list <- lapply(res_PL, function(x) {
    x[["pred_value_list"]][[i]]
  })
  
  # Calculate metric
  curr_metric <- get_metrics(curr_predicted_values_list,
                             true_values_list)
  
  # current situatiuon
  curr_sit <- all_sits[j]
  
  # Add the result to the list
  all_res$`impute, maximise n`[[curr_sit]] <- curr_metric
  
  # Count up Index
  j <- j + 1
}

# 3-4-3-4 MDD_SPLS
predicted_values_ddspls_list <- lapply(res_PL, function(x) {
  as.vector(x$pred_ddspls[, "1"])
})

all_res[["mdd_spls"]] <- get_metrics(predicted_values_ddspls_list, truth = true_values_list)

# 3-5 Save the metrics
save(all_res, file = "./docs/CV_Res/REAL/Hagenberg_5_3_4__Setting3.R")

# Setting 4 [4, 2, 1, 6, 3]  ---------------------------------------------------
# 3-1 Load the CV results ['Priority Lasso' & 'mdd-sPLS']
res_PL <- readRDS("./docs/CV_Res/REAL/PriorityLasso/5_3_4/data/results_permute_blocks_small_4_2020_05_28.Rds")

# 3-2 Get the true responses
true_values_list <- lapply(res_PL, function(x) x$test_y)
true_values      <- do.call("c", true_values_list)

# 3-3 Get the indices (all combinations from block 1-5)
indices <- 1:20

# 3-4 loop over the remaining indices and get PL results
# 3-4-1 Initalize list to save results
list_names      <- c("ignore, zero", "ignore, intercept",
                     "impute maximise blocks", "impute, maximise n")
all_res         <- vector("list", length(list_names))
names(all_res)  <- list_names

# 3-4-2 Define the different situations [which blocks used]
all_sits <- c("block_4", "block_4_2", "block_4_2_1", 
              "block_4_2_1_6", "block_4_2_1_6_3")

# 3-4-3 Calculate the metrics for the different methods and blocks for predictions
# 3-4-3-1 IGNORE ZERO
j <- 1
for (i in seq(from = 1, to = 20, by = 4)) {
  
  # Get current predictions
  curr_predicted_values_list <- lapply(res_PL, function(x) {
    x[["pred_value_list"]][[i]]
  })
  
  # Calculate metric
  curr_metric <- get_metrics(curr_predicted_values_list,
                             true_values_list)
  
  # current situatiuon
  curr_sit <- all_sits[j]
  
  # Add the result to the list
  all_res$`ignore, zero`[[curr_sit]] <- curr_metric
  
  # Count up Index
  j <- j + 1
}

# 3-4-3-2 IGNORE INTERCEPT
j <- 1
for (i in seq(from = 2, to = 20, by = 4)) {
  
  # Get current predictions
  curr_predicted_values_list <- lapply(res_PL, function(x) {
    x[["pred_value_list"]][[i]]
  })
  
  # Calculate metric
  curr_metric <- get_metrics(curr_predicted_values_list,
                             true_values_list)
  
  # current situatiuon
  curr_sit <- all_sits[j]
  
  # Add the result to the list
  all_res$`ignore, intercept`[[curr_sit]] <- curr_metric
  
  # Count up Index
  j <- j + 1
}

# 3-4-3-3 IMPUTE MAXIMISE BLOCKS
j <- 1
for (i in seq(from = 3, to = 20, by = 4)) {
  
  # Get current predictions
  curr_predicted_values_list <- lapply(res_PL, function(x) {
    x[["pred_value_list"]][[i]]
  })
  
  # Calculate metric
  curr_metric <- get_metrics(curr_predicted_values_list,
                             true_values_list)
  
  # current situatiuon
  curr_sit <- all_sits[j]
  
  # Add the result to the list
  all_res$`impute maximise blocks`[[curr_sit]] <- curr_metric
  
  # Count up Index
  j <- j + 1
}

# 3-4-3-3 IMPUTE MAXIMISE N
j <- 1
for (i in seq(from = 3, to = 20, by = 4)) {
  
  # Get current predictions
  curr_predicted_values_list <- lapply(res_PL, function(x) {
    x[["pred_value_list"]][[i]]
  })
  
  # Calculate metric
  curr_metric <- get_metrics(curr_predicted_values_list,
                             true_values_list)
  
  # current situatiuon
  curr_sit <- all_sits[j]
  
  # Add the result to the list
  all_res$`impute, maximise n`[[curr_sit]] <- curr_metric
  
  # Count up Index
  j <- j + 1
}

# 3-4-3-4 MDD_SPLS
predicted_values_ddspls_list <- lapply(res_PL, function(x) {
  as.vector(x$pred_ddspls[, "1"])
})

all_res[["mdd_spls"]] <- get_metrics(predicted_values_ddspls_list, truth = true_values_list)

# 3-5 Save the metrics
save(all_res, file = "./docs/CV_Res/REAL/Hagenberg_5_3_4__Setting4.R")