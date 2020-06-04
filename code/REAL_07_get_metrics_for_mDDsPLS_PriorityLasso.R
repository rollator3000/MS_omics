"Script to calculate the Metrics of the Priority-Lasso / mdd-sPLS
  - only recived the predictions of the models!
"
# [0] Librarys and Functions
library(assertthat)

mcc_metric  <- function(conf_matrix) {
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
get_metrics <- function(preds, truth) {
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

# [1] Section 5.3.1
# 1-1 Load the CV results ['Priority Lasso' & 'mdd-sPLS']
res_PL <- readRDS("./docs/CV_Res/REAL/PriorityLasso/5_3_1/data/results_different_blocks_2020_05_06.Rds")
res_MD <- readRDS("./docs/CV_Res/REAL/PriorityLasso/5_3_1/data/results_different_blocks_ddsPLS_weight_2020_05_20.Rds")

# 1-2 Get the true responses
true_values_list <- lapply(res_PL, function(x) x$test_y)
true_values      <- do.call("c", true_values_list)

# 1-3 Get Metrics for mdd-sPLS
# 1-3-1 MD Predictions
predicted_values_ddspls_list <- lapply(res_MD, function(x) {
  as.vector(x$pred_value_ddspls[, "1"])
})

MD_SPLSS <- get_metrics(predicted_values_ddspls_list, true_values_list)


# 1-3 Get the predicted Values
# 1-3-1 Indices depending on the method
indices <- 21:24

# 1-3-2 PL predictions
predicted_values_list_PL <- lapply(res_PL, function(x) {
  x[["pred_value_list"]][[indices[1]]]
})



# Calculate the Metrics
a1 <- get_metrics(predicted_values_list_PL, true_values_list)



mean(unlist(a1["AUC1",]))


i_ = 1
predicted_classes_list_PL <- ifelse(predicted_values_list_PL[[i_]] > 0.5,  1, 0)
confmat                   <- caret::confusionMatrix(as.factor(predicted_classes_list_PL), 
                                                    as.factor(true_values_list[[i_]]))

# --5-2 Are under the ROC Curve
roc1 <- pROC::auc(pROC::roc(as.numeric(predicted_classes_list_PL), 
                            as.numeric(as.character(true_values_list[[i_]])),
                            levels = levels(as.factor(as.character(true_values_list[[i_]]))),
                            direction = "<"))

roc2 <- pROC::auc(pROC::roc(as.numeric(predicted_classes_list_PL), 
                            as.numeric(as.character(true_values_list[[i_]])),
                            levels = levels(as.factor(as.character(true_values_list[[i_]]))),
                            direction = ">"))

# --5-3 MCC Matthews correlation coefficient [only for binary cases!]
mcc <- mcc_metric(conf_matrix = confmat)

# --5-4 Collect all metrics in a list & replace the not defined values
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
