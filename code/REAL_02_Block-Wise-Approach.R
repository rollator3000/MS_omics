"
Script to analyse the real multi-omics data set with the block-wise approach
"
# [0] Set WD, load packages & define functions
setwd("C:/Users/kuche_000/Desktop/MS_omics/")
library(randomForestSRC)
library(checkmate)
library(caret)

mcc_metric            <- function(conf_matrix) {
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
get_blockwise_pred    <- function(Forest, testdata, weights) {
  " For a collection of blockwise fitted RFs [each RF fitted on diff. block] 
    get the aggregated predicition from all trees for a single/ multiple obs.!
    Important these need the same observed features!
    Then each RF creates a prediciton - if possible. These predicitions are
    then averaged in a weighted way and used to find the final predicted class!
  
     Args:
      - Forest (list)         : list filled with the objects of class 'rfsrc'
      - testdata (data.frame) : testdata we want to get predicitons for!
                                Must conatain the response we learned the Forest
                                on [1st. column]
      - weights (vector)      : The weights for the different BlockWise RFs!
                                As some of the blockwise RFs (in Forest) have a 
                                higher Accuracy/ F1/ ... than others, they have 
                                a higher predictive power!
                                  We can weight the different predicitons based
                                  on the oob performance, so that RFs w/ good 
                                  pred. performance have higher influence!
                                  --> if we want equally weighted blocks, pass
                                      a vector filled w/ '1'
                                  --> needs same length as Forest!
     Return:
      - Class witht he highest probabilitiy! 
  "
  # [0] Check Inputs -----------------------------------------------------------
  # 0-1 Reponse Class in the testdata
  assert_data_frame(testdata, min.rows = 1)
  if (!(Forest[[1]]$yvar.names %in% colnames(testdata))) stop("testdata is missing response column!")
  if (colnames(testdata)[1] != Forest[[1]]$yvar.names) stop("'testdata' needs the response column of the trees as first feature")
  
  # 0-2 All Elements of Forest should be of class "rfsrc"
  if (any(sapply(Forest, FUN = function(x) !("rfsrc" %in% class(x))))) {
    stop("not all elements in 'Forest' are of class 'rfsrc'")
  }
  
  # 0-3 weights should be vector of exactly the length of 'Forest'
  assert_vector(weights, len = length(Forest))
  if (any(sapply(weights, function(x) !is.numeric(x)))) {
    print("'weights:\n")
    print(weights)
    stop("'weights' is not only filled with numerics!")
  }
  
  # [1] Remove RFs, that use split variables not avaible in the testdata!  ----
  # 1-1 Get the ColNames we have in our testdata
  test_cols <- colnames(testdata[which(!is.na(testdata))])
  
  # 1-2 Check for each RF whether it uses a variable not existent in testdata!
  # 1-2-1 Before that recored of how many blockwise RFs it originally consist of
  forest_orignal_length <- length(Forest)
  
  # 1-2-2 Get Index of the Forests that need to be removed
  forrest_to_rm <- sapply(Forest, FUN = function(x) any(!(x$xvar.names %in% test_cols)))
  
  # 1-3 If any has to be removed, also remove the corresponding OOB Metric of the tree!
  if (any(forrest_to_rm)) {
    Forest  <- Forest[-c(which(forrest_to_rm))]
    weights <- weights[-c(which(forrest_to_rm))]
  } 
  
  # 1-3-1 Check whether the remaning 'weights' are not only 0's
  #       --> this leads to an error in 'weighted.mean()' calculation
  if (all(weights == 0)) {
    weights <- rep(1, times = length(weights))
  }
  
  # 1-4 Check whether there are any forrests left to do predicitons with 
  #     & if so, print the amount of usable trees!
  if (length(Forest) < 1) {
    print("Forest can not predicit on TestData, as all trees use split vars not avaible in 'testdata'")
    return(NA)
  }
  
  # 1-5 Print Info, how many of the blockwise fitted RF were removed!
  print(paste0(sum(forrest_to_rm), "/", forest_orignal_length , " blockwise RFs of totally had to be removed from 'Forest', as these use splitvars, not in 'testdata'"))
  
  # [2] Use the remaining RFs to create predicitons for the testdata -----------
  predicitions <- lapply(Forest, FUN = function(x) predict(x, testdata))
  
  # [3] Aggregate the Predicitons! ---------------------------------------------
  prob_class0 <- c()
  for (j in 1:nrow(testdata)) {
    prob_class0 <- c(prob_class0,
                     weighted.mean(sapply(1:length(Forest), 
                                          FUN = function(x) predicitions[[x]]$predicted[j, 1]),
                                   w = weights, na.rm = TRUE))
  }
  
  prob_class0 <- ifelse(prob_class0 >= 0.5, 0, 1)
  return(prob_class0)
}
get_oob_weight_metric <- function(blockwise_RF) {
  "Funciton to evaluate how good the predicitive performance of 
   a single blockwise fitted RF is! For this we let the already fitted model 
   do predicitons on its OOB observations and get metrics! 
   [use all trees w/ same OOB observation to do a prediciton!]
   
    Args:
      - blockwise_RF (rfsrc) : RF that was alredy fit on data! 
      
    Return:
      - return the OOB_F1_Score and the OOB_Accuracy
  "
  # [0] Check Inputs  ----------------------------------------------------------
  if (!("rfsrc" %in% class(blockwise_RF))) {
    stop("'blockwise_RF' is not of class 'rfsrc'!")
  }
  
  # [1] Calculate the metrics based on OOB observations!  ----------------------
  # 1-1 Get probabilites for class '0' and remove NAs [ID was in no tree oob!]
  preds_class_0                            <- blockwise_RF$predicted.oob[,1]
  to_rm_NA                                 <- which(is.na(preds_class_0))
  if (length(to_rm_NA) >= 1) preds_class_0 <- preds_class_0[-to_rm_NA] 
  preds_class                              <- ifelse(preds_class_0 > 0.5, 0, 1)
  
  # 1-2 Get true classes and remove the ID, which hasn't been OOB at all!
  true_classes <- blockwise_RF$yvar
  if (length(to_rm_NA) >= 1) true_classes <- true_classes[-to_rm_NA] 
  
  # 1-3 Convert 'preds_class' to factor w/ same levels!
  preds_class <- factor(preds_class, levels = levels(true_classes))
  
  # 1-4 Get Predicitons
  confmat_ <- caret::confusionMatrix(data      = true_classes, 
                                     reference = preds_class)
  
  # [2] Return the F1-Score and the Accuracy!  ---------------------------------
  # 2-1 Extract F1-Score and Accuracy from 'confmat_'
  F1_curr  <- confmat_$byClass["F1"]
  Acc_curr <- confmat_$overall["Accuracy"]
  
  # 2-2 If F1 is not defined replace it by 0!
  if (is.na(F1_curr)) F1_curr <- 0 
  
  # 2-3 Create Vector with names of the metrics!
  oob_results <- c(F1_curr, Acc_curr)
  names(oob_results) <- c("F1", "Acc")
  
  # 2-4 Return the named vector!
  return(oob_results)
}
eval_predicitons      <- function(predicitons, reference) {
  " Calculate the Metrics for a given list of prediciotns and a given list of 
    true responses - response must be categorical!
  
  Args:
    - predicitons (vector) : Vector filled with integers - predicitons from 
                             a model! Must be of class factor! 
    - reference   (vector) : Vector filled with integers - true responses!
                             Must be of class factor!
  
  Return:
    - list of metrics with F1, Accuracy, MMC, AUC, ...
  "
  # [0] Check Inputs  ----------------------------------------------------------
  # 0-1 'predicitons' & 'reference' must be both of class factor!
  assert_factor(predicitons)
  assert_factor(reference)
  
  # 0-2 Must have the same lenght
  if (length(predicitons) != length(reference)) {
    stop("predicitons' & 'reference' are of different length!")
  }
  
  # [1] Evaluate the predicitons  ----------------------------------------------
  # 1-1 Confusion Matrix - from which we can calc/ extract most metrics
  confmat <- caret::confusionMatrix(data      = predicitons, 
                                    reference = reference)
  
  # 1-2 Area under the ROC Curve
  roc1 <- tryCatch(pROC::auc(pROC::roc(testdata[,1], all_forrest_preds_probs_class_0, 
                                       levels = levels(Forest[[1]][[1]]$data$data[,1]), 
                                       direction = "<")),
                   error = function(e) "not defined!")
  if (is.numeric(roc1)) roc1 <- as.numeric(roc1)
  
  roc2 <- tryCatch(pROC::auc(pROC::roc(testdata[,1], all_forrest_preds_probs_class_0, 
                                       levels = levels(Forest[[1]][[1]]$data$data[,1]), 
                                       direction = ">")),
                   error = function(e) "not defined!")
  if (is.numeric(roc2)) roc2 <- as.numeric(roc2)
  
  # 1-3 MCC Matthews correlation coefficient [only for binary cases!]
  mcc <- mcc_metric(conf_matrix = confmat)
  
  # [2] Create a list to collect the results!  ---------------------------------
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
              "AUC1"        = roc1,
              "AUC2"        = roc2,
              "MCC"         = mcc)
  
  # 6-1 If the F1-Score/ Precision/ Recall is NA, then we set it to 0 [-1 for MCC]! 
  if (is.na(res$F1))             res$F1             <- 0
  if (is.na(res$Precision))      res$Precision      <- 0
  if (is.na(res$Recall))         res$Recall         <- 0
  if (is.na(res$MCC))            res$MCC            <- -1
  if (is.na(res$Pos_Pred_Value)) res$Pos_Pred_Value <- 0
  if (is.na(res$Neg_Pred_Value)) res$Neg_Pred_Value <- 0
  
  return(res)
}

# [1] Read in the data and some preprocessing                               ----
# 1-1 Read in the data with Block-Wise missingness & only keep the merged DF!
load("./data/processed/real_data/data 05052020.RData")
rm(df1, df1_mod, df2, df2_mod, df3, df3_mod, 
   df4, df4_mod, df51, df51_mod, df53, df53_mod)

# 1-2 Create DF with information to the outcomes for the CV
#     Completlty copied from Hagenberg!
#         > needed for CV 
index_df <- data.frame(index = 1:521,
                       outcome = y$outcome,
                       group = 0)

#         > weight the obs. such that there is a 50/50 mix of 0 & 1
set.seed(8274)
index_outcome_0 <- sample(rep(1:5, each = 53))
index_outcome_1 <- sample(c(rep(1:5, each = 51), 1))
index_df[index_df$outcome == 0, "group"] <- index_outcome_0
index_df[index_df$outcome == 1, "group"] <- index_outcome_1


# 1-3 Read in the Block-Structure of the data. Tells us which observations
#     have which observed feature blocks
missing_str <- read.csv("./data/processed/real_data/Block Structure.csv", 
                        sep = ";")
missing_str$outcome <- NULL

# [2] Start the 5-fold CV - BLOCK-WISE Approach                             ----
# 2-1 Set the arguments for the fold-wise RandomForest Models!
num_trees         = 25
mtry              = NULL
min_node_size     = 5

# 2-2 Create a list to save the metrics
blockwise_res <- list()

# 2-3 Start the 5 fold CV
for (i in 1:5) {
  
  # 2-3-1 Print Current Fold Status
  print(paste0("FOLD ", i, "/ 5 -------------------------------"))
  
  # 2-3-2 Get the current Test- & Train-Set
  index_test  <- index_df[index_df$group == i, "index"]
  index_train <- index_df[index_df$group != i, "index"]
  train_x <- data_merged[index_train, ]
  train_x <- cbind("outcome_Y" = as.factor(index_df[index_train, "outcome"]), 
                   train_x)
  test_x <- data_merged[index_test, ]
  test_x <- cbind("outcome_Y" = as.factor(index_df[index_test, "outcome"]), 
                  test_x)
  
  # 2-3-3 Grow a Block-Wise RF - each colname in'missing_str' stands for a block
  Forest <- list()
  i_     <- 1
  for (block_ in colnames(missing_str)) {
    
    # --1 Get the Observations that have the features of 'block_'
    rows_ <- row.names(missing_str)[which((row.names(missing_str) %in% row.names(train_x)) & 
                                            missing_str[block_])]
    rows_ <- as.integer(rows_)
    
    # --2 Get the columns of the current Block
    cols_ <- grep(paste0(block_, "_"), colnames(data_merged))
    
    # --2-1 add the response as column aswell!
    cols_ <- c(1, cols_)
    
    # --3 Fit a RF on this fully observed (fold-)subdata!
    # --3-1 Define formula
    response    <- 'outcome_Y'
    formula_all <- as.formula(paste(response, " ~ ."))
    
    # --3-2 Only keep the usable rows & columns and add the response
    curr_fold_train_data <- train_x[which(row.names(train_x) %in% rows_), cols_]
    
    # --4 Fit a Tree on the block data
    print(paste0("Fit FoldWise RF on current block: '", block_, "'"))
    blockwise_rf <- rfsrc(formula = formula_all, data = curr_fold_train_data, 
                          ntree = num_trees, mtry = mtry, nodesize = min_node_size, 
                          samptype = "swr", seed = 12345678, var.used = 'all.trees')
    
    Forest[[i_]] <- blockwise_rf
    rm(blockwise_rf)
    
    # 2-3-5 Count up i_ - used to fill Forest list!
    i_  = i_ + 1
  }
  
  # 2-3-4 Get the OOB_Metrics for the block-wise fitted models
  f1_weights  <- c()
  acc_weights <- c()
  no_weights  <- rep(1, times = length(Forest))
  for (block_RF in Forest) {
    
    curr_weights <- get_oob_weight_metric(block_RF)
    
    f1_weights  <- c(f1_weights, curr_weights["F1"])
    acc_weights <- c(acc_weights, curr_weights["Acc"])
  }
  
  # 2-3-5 Get the predicitions for each observation for all different weightings!
  # 2-3-5-1 Define the vectors to save results
  predicted_no_weight  <- c()
  predicted_f1_weight  <- c()
  predicted_acc_weight <- c()
  true_reponse         <- test_x$outcome_Y
  
  # 2-3-5-2 Start Looping over the different test-folds
  for (index in 1:nrow(test_x)) { 
    
    # --0 Print progress
    print(paste("Evaluation:", index, "/", nrow(test_x), "----"))
    
    predicted_f1_weight <- c(predicted_f1_weight,
                             get_blockwise_pred(Forest = Forest, testdata = test_x[index,], 
                                                weights = f1_weights))
    
    predicted_acc_weight <- c(predicted_acc_weight, 
                              get_blockwise_pred(Forest = Forest, testdata = test_x[index,], 
                                                 weights = acc_weights))
    
    predicted_no_weight <- c(predicted_no_weight,
                             get_blockwise_pred(Forest = Forest, testdata = test_x[index,], 
                                                weights = no_weights))
  }
  
  # 2-3-6 Compare the predicted classes with the true response
  no_weight  <- eval_predicitons(as.factor(predicted_no_weight), 
                                 true_reponse)
  acc_weight <- eval_predicitons(as.factor(predicted_acc_weight), 
                                 true_reponse)
  f1_weight  <- eval_predicitons(as.factor(predicted_f1_weight), 
                                 true_reponse)
  
  # 2-3-7-4 Add the results of the evaluation to the 'foldwise_res' list!
  blockwise_res[[i]] <- list("no_weight"  = no_weight,
                            "acc_weight" = acc_weight,
                            "f1_weight"  = f1_weight)
}

# 2-4 Save the results of the CV
save(blockwise_res, file = "./docs/CV_Res/REAL/BlockWise_Approach.R")
