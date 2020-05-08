"
Script to analyse the real multi-omics data set.
"
# [0] Set WD, load packages & define functions
setwd("C:/Users/kuche_000/Desktop/MS_omics/")
source("./code/GENERAL_simpleRF_adaption.R")
library(pROC)
library(assertthat)
library(checkmate)
library(doParallel)
library(e1071)

detectCores()
registerDoParallel(cores = 2)

all_trees_grown_correctly <- function(trees) {
  " Check, whether 'trees', were grown correctly & if not grow these 
    trees again, as long, as they are grown correctly! 
        --> Growning not correctly: No 'childNodeIDs', no split variables etc...
  
  Args:
   - trees (list) : list filled with object of the class 'Tree'! 
                    For each object in there check, whether it was grown 
                    correctly (has childNodeIDs) - if not grown it again!
  
  Return: 
    - list of trees, where all of these trees were grown correctly
        --> each tree has at least one 'childNodeID'
  "
  # [0] Check Inputs  ----------------------------------------------------------
  # 0-1 All Objects in 'trees' of class 'Tree'
  trees_classes <- sapply(trees, function(x) class(x))
  if (any(!grepl("Tree", trees_classes))) {
    stop("Not all Objects in 'trees' are of class 'Tree'")
  }
  
  # [1] Get the entrance of the objects, that miss child node IDs  -------------
  wrong_trees  <- unlist(lapply(1:length(trees), 
                                FUN = function(x) {
                                  if (length(trees[[x]]$child_nodeIDs) == 0) x
                                }))
  
  # [2] Regrow the trees, that were not grown correctly  -----------------------
  #     If there are any trees not grown correctly, grow them again until all
  #     of the trees were grown correctly!
  while (length(wrong_trees) > 0) {
    
    # grow the errours trees again
    trees[wrong_trees] <- lapply(trees[wrong_trees], 
                                 function(x) {
                                   x$grow(replace = TRUE)
                                   x
                                 })
    
    # check whether any of the trees is not grown correctly!
    wrong_trees  <- unlist(lapply(1:length(trees), 
                                  FUN = function(x) {
                                    if (length(trees[[x]]$child_nodeIDs) == 0) x
                                  }))
  }
  
  # [3] Return the correclty grown trees  --------------------------------------
  return(trees)
}
get_oob_weight_metric     <- function(trees) {
  " Calculate OOB Metrirc ['F1' & 'Acc'] of a list of trees! 
    For this we go  through all OOB Predictions and obtain aggregated 
    predicitons from all trees, that have the same observation as out-of-bag!
    In case the metrics are 'NA' [not defined] they are repalced by their worst
    possible value
  
  Args: 
    - trees (list) : list filled w/ objects of class 'Tree'
  
  Return:
    - Average oob-Acc & oob-F1 metric for 'trees'!
  "
  # [0] Check Input ------------------------------------------------------------
  # 0-1 Make sure 'trees' is a list filled with 'Trees'
  assert_list(trees, min.len = 1, any.missing = FALSE)
  if (any(sapply(trees, function(x) 'Tree' %in% class(x)))) {
    stop("not all elements in 'trees' are of class 'Tree'")
  }
  
  # [1] Get the OOB Predicitons ------------------------------------------------
  # 1-1 Get the trees, that are usable - not pruned in the first split_variable!
  usable_trees <- sapply(1:length(trees), function(x) {
    
    # Check whether the first split_var was pruned!
    if (trees[[x]]$child_nodeIDs[[1]][1] != "pruned") {
      return(x)
    } 
  })
  
  # 1-1-1 If all the trees were pruned in the first node, we can not do
  #       any OOB predictions --> OOB-Accuracy = 0
  if (length(usable_trees) < 1) {
    return(0)
  } else {
    usable_trees <- unlist(usable_trees)
  }
  
  # 1-2 For all usable trees [not pruned in first split variable] get the IDs
  #     of the OOB observations. Then collect all and only keep unique IDs!
  oob_ids_all_trees <- sapply(usable_trees, function(x) trees[[x]]$oob_sampleIDs)
  unique_oob_ids    <- unique(unlist(oob_ids_all_trees))
  
  # 1-3 Loop over all 'unique_oob_ids' and get the oob predictions from all 
  #     trees, that have the same OOB observation!
  all_oob_preds_class0 <- c()
  for (curr_oob in unique_oob_ids) {
    
    # 1-3-1 Get all trees that have 'curr_oob' as OOB observation!
    trees_same_oob <- unlist(sapply(usable_trees, function(x) {
      if (curr_oob %in% trees[[x]]$oob_sampleIDs) x
    }))
    
    # 1-3-2 Get the feas of the observation that is OOB for 'trees_same_oob' 
    curr_oob_feas <- Data$new(data = trees[[1]]$data$data[curr_oob,])
    
    # 1-3-3 Get a Prediciton for the 'curr_oob' from all trees!
    predicted_probs <- sapply(trees_same_oob, 
                              function(x) trees[[x]]$predict(curr_oob_feas))
    
    # 1-3-4 Aggregate the predictions from the different trees and
    #       Get the probability for class 0! [First Row is class 0]
    all_oob_preds_class0 <- c(all_oob_preds_class0, 
                              sum(predicted_probs[1,]) / length(trees_same_oob))
  }
  
  # 1-4 Convert the probs to classes
  predicted_oob_classes <- ifelse(all_oob_preds_class0 > 0.5, 0, 1)
  predicted_oob_classes <- factor(predicted_oob_classes, 
                                  levels = levels(trees[[1]]$data$data[unique_oob_ids, 1]))
  
  # [2] Compare predicted Classes with the true classes & get the 'weight_metric'!
  conf_mat_oob <- caret::confusionMatrix(predicted_oob_classes, 
                                         trees[[1]]$data$data[unique_oob_ids, 1])
  
  # 2-1 Check, that neither Acc nor F1 is not defined, if replace it with worst value
  Acc_ <- conf_mat_oob$overall["Accuracy"]
  f1_  <- conf_mat_oob$byClass["F1"]
  
  if (is.na(Acc_)) Acc_ <- 0 
  if (is.na(f1_))  f1_  <- 0 
  
  # [3] Return it as named vector!
  return(setNames(c(Acc_, f1_),c("Acc", "F1")))
}
get_predicition           <- function(Forest, testdata) {
  " Get the aggregated predicition from all trees for each of the observations
    in testdata! Evaluate the aggregated predicitons & return metrics!
    CAUTION: The observation in 'testdata' need the same observed features,
             as the trees of the fold-wise forests need to pruned!
  
  Args:
    - Forest (list)         : list filled with the objects of class 'Tree'
    - testdata (data.frame) : Testdata we want to get predicitons for!
                              Must conatain the response we learned the
                              Forest on!
  Return:
    - list w/ metrics [accuracy, f1, mcc, roc, ...] 
      > Metrics w/ NA are replaced by their worst possible value!
  "
  # [0] Check Inputs -----------------------------------------------------------
  # 0-1 All Elements in Forest of class 'Tree'
  for (i in 1:length(Forest)) {
    tree_classes <- sapply(1:length(Forest[[i]]), function(x) {
      grepl("Tree", class(Forest[[i]][[x]]))
    })
    
    if (!all(tree_classes)) {
      stop("Forest contains a Object not of class 'Tree'")
    }
  }
  
  # 0-2 Check testdata for DF 
  assert_data_frame(testdata, min.rows = 1)
  
  # [1] Prepare TestData  ------------------------------------------------------
  #     Convert TestData to same Format as the data the trees were originally
  #     trained with, to ensure factor levels/ features are the same ....
  tree_testsets <- list()
  for (i in 1:length(Forest)) {
    tree_testsets[[i]] <- process_test_data(tree = Forest[[i]][[1]], 
                                            test_data = testdata)
  }
  
  # [2] Get Predicitons  ------------------------------------------------------- 
  #     Get a prediction for every observation in TestData from all foldwise
  #     fitted RandomForests ['unique_folds' predicitons per testobs.] 
  tree_preds_all <- list()
  tree_preds_all <- foreach(i = 1:length(Forest)) %do% {
    
    # save the predictions as 'treeX_pred'
    return(get_pruned_prediction(trees = Forest[[i]], 
                                 test_set = tree_testsets[[i]]))
  }
  
  # Check whether any of the RFs is not usable [only NA predicitons]
  not_usable <- sapply(seq_len(length(tree_preds_all)), function(i) {
    all(is.na(tree_preds_all[[i]]$Class))
  })
  
  # 2-1 Check that there are still trees existing, else not preds possible!
  if (all(not_usable)) {
    print("None of the trees are usable for predictions!")
    return("No Predictions possible as all trees are pruned at the 1. splitvar!")
  }
  
  # 2-2 Remove the SubRFs, that are not usable from Forest, tree_preds_all & tree_testsets
  if (any(not_usable)) {
    Forest         <- Forest[-c(which(not_usable))]
    tree_preds_all <- tree_preds_all[-c(which(not_usable))]
    tree_testsets  <- tree_testsets[-c(which(not_usable))]
  }
  
  # [3] If we want to create weighted ensemble of the predicitons, we need to ----
  #     calc the OOB-Accuracy/ -F1-Score per foldwise fitted RF. These are used
  #     as weights when assembeling the predicitons [low Metric -> low weight]!
  #     ATTENTION: No need to calc weights when only a single foldwise RF left!
  
  # 3-1 Firstly prune the trees according to the current testset & then get the oob 
  #     performance, so the pruned RF is rated and not the fully grown one 
  #     [where test set misses the used split variables of the trees]!
  #     --> child nodeIDs to 'pruned', if split_var not in test
  
  #    Loop over all trees and prune them according to the testdata!
  #    [--> had to be added, as the pruning is not saved when running in parallel]
  for (i_ in 1:length(Forest)) {
    curr_test_set <- tree_testsets[[i_]]
    tmp <- sapply(Forest[[i_]], FUN = function(x) x$prune(curr_test_set))
  }
  
  # 3-2 Get the oob performance of the pruned trees!
  print("Calculation OOB Metrics...")
  if (length(Forest) > 1) {
    
    # Get the ACC & F1 weights for each of the foldwise fitted Forests!
    tree_weights <- foreach(l = seq_len(length(Forest))) %do% { # par
      get_oob_weight_metric(trees = Forest[[l]])
    }
    
    # Save them in seperate vectors
    F1_weights <- sapply(seq_len(length(tree_weights)), 
                         function(l) tree_weights[[l]]["F1"])
    Acc_weights <- sapply(seq_len(length(tree_weights)), 
                          function(l) tree_weights[[l]]["Acc"])
    No_weights  <- rep(1, times = length(tree_weights))
  } else {
    
    # Fill the Vectors with 1 - if the whole 'Forest' only consitis of 
    # 1 single foldwise fitted RF!
    No_weights  <- rep(1, times = length(Forest))
    Acc_weights <- rep(1, times = length(Forest))
    F1_weights  <- rep(1, times = length(Forest))
  }
  
  # [4] Aggregate Predictions from the different foldwise fitted RFs!  ---------
  # --- 4-1 Get the probabilities of all test obs to be of class '0' - NO WEIGHTS
  probs_class_0_no_weight <- sapply(1:nrow(testdata), FUN = function(x) {
    
    # Get a probability prediciton from each [still usable] tree!
    preds_all <- sapply(seq_len(length(tree_preds_all)), function(i) {
      tree_preds_all[[i]]$Probs[[x]][1]
    })
    
    # Combine the preditions of the different trees!
    prob_class0 <- weighted.mean(preds_all, w = No_weights, na.rm = TRUE)
    prob_class0
  })
  # 4-1-1 Convert probabilities to class predicitons!
  preds_class_0_no_weight <- ifelse(probs_class_0_no_weight >= 0.5, 0, 1)
  preds_class_0_no_weight <- factor(preds_class_0_no_weight, 
                                    levels = levels(Forest[[1]][[1]]$data$data[,1]))
  
  
  # --- 4-2 Get the probabilities of all test obs to be of class '0' - ACC WEIGHTS
  probs_class_0_acc_weight <- sapply(1:nrow(testdata), FUN = function(x) {
    
    # Get a probability prediciton from each [still usable] tree!
    preds_all <- sapply(seq_len(length(tree_preds_all)), function(i) {
      tree_preds_all[[i]]$Probs[[x]][1]
    })
    
    # Combine the preditions of the different trees!
    prob_class0 <- weighted.mean(preds_all, w = Acc_weights, na.rm = TRUE)
    prob_class0
  })
  # 4-2-1 Convert probabilities to class predicitons!
  preds_class_0_acc_weight <- ifelse(probs_class_0_acc_weight >= 0.5, 0, 1)
  preds_class_0_acc_weight <- factor(preds_class_0_acc_weight, 
                                     levels = levels(Forest[[1]][[1]]$data$data[,1]))
  
  
  # --- 4-3 Get the probabilities of all test obs to be of class '0' - F-1 WEIGHTS
  probs_class_0_f1_weight <- sapply(1:nrow(testdata), FUN = function(x) {
    
    # Get a probability prediciton from each [still usable] tree!
    preds_all <- sapply(seq_len(length(tree_preds_all)), function(i) {
      tree_preds_all[[i]]$Probs[[x]][1]
    })
    
    # Combine the preditions of the different trees!
    prob_class0 <- weighted.mean(preds_all, w = F1_weights, na.rm = TRUE)
    prob_class0
  })
  # 4-3-1 Convert probabilities to class predicitons!
  preds_class_0_f1_weight <- ifelse(probs_class_0_f1_weight >= 0.5, 0, 1)
  preds_class_0_f1_weight <- factor(preds_class_0_f1_weight, 
                                    levels = levels(Forest[[1]][[1]]$data$data[,1]))
  
  
  # [5] Create an vector with the predicted classes for the current observation
  res_all <- c("No_Weighting"  = as.character(preds_class_0_no_weight),
               "Acc_Weighting" = as.character(preds_class_0_acc_weight),
               "F1_Weighting"  = as.character(preds_class_0_f1_weight)) 
  return(res_all)
}
mcc_metric                <- function(conf_matrix) {
  " Calculate the MCC Metric [Matthews correlation coefficient]!
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
eval_predicitons          <- function(predicitons, reference) {
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

# [2] Start the 5-fold CV - FOLD-WISE Approach                              ----
# 2-1 Set the arguments for the fold-wise RandomForest Models!
num_trees         = 300
mtry              = NULL
min_node_size     = 5

# 2-2 Create a list to save the metrics
foldwise_res <- list()

# 2-3 Start the 5 fold CV
for (i in 1:5) {
  
  # 2-3-1 Print Current Fold Status
  print(paste0("FOLD ", i, "/ 5 -------------------------------"))
  
  # 2-3-2 Get the current Test- & Train-Set
  index_test  <- index_df[index_df$group == i, "index"]
  index_train <- index_df[index_df$group != i, "index"]
  train_x <- data_merged[index_train, ]
  train_y <- index_df[index_train, "outcome"]
  test_x <- data_merged[index_test, ]
  test_y <- index_df[index_test, "outcome"]
  
  #      - Paste the response to into the data frame
  train_x$outcome_Y <- as.factor(train_y)
  test_x$outcome_Y  <- as.factor(test_y)
  
  # 2-3-3 Get the Observations that belong to the same fold [same feature space]
  #       - Get for each obs. the index of the observed feas
  observed_feas <- foreach(x = seq_len(nrow(train_x)), .combine = 'c') %dopar% {
    paste0(which(!(is.na(train_x[x,]))), collapse = "_")
  }
  
  #       - Keep the unique observed feas [equals the different folds]
  #         That we use to assign obs. to the differnt folds!
  observed_folds <- unique(observed_feas)
  print(paste0("Found ", length(observed_folds), " unique folds!"))
  
  # 2-3-4 Train foldwise RFs for each fold seperatly! For this loop over all 
  #       folds [observed_folds contain all observed features for diff folds]!
  #       --> Results in a Forest of length 'length(observed_folds)' 
  #           & each entrance consits of 'num_trees' foldwise fitted trees!
  Forest <- list()
  Forest <- foreach(j_ = 1:length(observed_folds)) %do% {
    
    # --1 Get the observed feas for the current fold
    fold_ = observed_folds[j_]
    
    # --2 Get all Obs. with the feture space as in 'fold_'
    fold_obs_ <- which(observed_feas == fold_)
    
    # --3 Get all the indeces of the columns observed with this fold!
    obs_columns_ <- as.numeric(unlist(strsplit(fold_, split = "_")))
    
    # --4 Get all Trainpoints from the obs. w/ same features + 
    #     only keep observed features of these! --> fully observed subdata!
    curr_fold_train_data <- train_x[fold_obs_, obs_columns_]
    
    if (length(unique(curr_fold_train_data$outcome_Y)) == 1) {
      print(paste("Growing on Fold", j_, "failed. TrainFold has steady response"))
      return(NULL)
    } 
    
    # --5 Fit a RF on this fully observed (fold-)subdata!
    formula_all <- as.formula("outcome_Y ~ .")
    
    # --6 Define the foldwise RF and fit it on the foldwise data
    fold_RF <- simpleRF(formula           = formula_all, 
                        data              = curr_fold_train_data, 
                        num_trees         = num_trees, 
                        mtry              = mtry, 
                        min_node_size     = as.integer(min_node_size),
                        replace           = TRUE,  # always TRUE, as we need OOB!
                        splitrule         = NULL,  # always NULL!
                        unordered_factors = "ignore")
    
    fold_RF <- lapply(fold_RF, function(x) {
      x$grow(replace = TRUE)
      x
    })
    
    # --7 Check that all trees were grown correctly
    #     --> none w/o 'child_node_ID' after that!
    fold_RF <- all_trees_grown_correctly(fold_RF)
    
    print(paste("Growing on Fold", j_, "successful"))
    return(fold_RF)
  }
  
  # 2-3-5 Remove the Forests without any entrance & check that length is > 1
  Forest <- Forest[lengths(Forest) != 0]
  if (length(Forest) == 0) stop("Forest ist not existent")
  
  # 2-3-6 Evaluate the fold-wise grown forest!
  # 2-3-6-1 Get the test-observations that belong to the same fold, and get 
  #         predicitions for all observations from the same fold!
  #        [important, as the trees are prunde dpending on the observed feas in test]
  observed_feas_test <- foreach(x = seq_len(nrow(test_x)), .combine = 'c') %dopar% {
    paste0(which(!(is.na(test_x[x,]))), collapse = "_")
  }
  
  # 2-3-6-2 Keep the unique observed test feas [equals the different folds]
  observed_test_folds <- unique(observed_feas_test)
  print(paste0("Found ", length(observed_test_folds), " unique test folds!"))
  
  # 2-3-7 Get the predicted classes from the RF for each test-fold!
  # 2-3-7-1 Define Vectors to save predicitions & true classes
  predicted_no_weight  <- c()
  predicted_f1_weight  <- c()
  predicted_acc_weight <- c()
  true_reponse         <- c()
  
  # 2-3-7-2 Start Looping over the different test-folds
  for (index in 1:length(observed_test_folds)) { 
    
    # --0 Print progress
    print(paste("Evaluation:", index, "/", length(observed_test_folds), "----"))
    
    # --1 Get the current TestFold
    fold_ = observed_test_folds[index]
    
    # --2 Get all Obs. with the feture space as in 'fold_'
    fold_obs_ <- which(observed_feas_test == fold_)
    
    # --3 Get all the indeces of the columns that were observed with this fold!
    obs_columns_ <- as.numeric(unlist(strsplit(fold_, split = "_")))
  
    # --4 Get the current fold - all observations and corresponding features
    curr_fold_test_data <- test_x[fold_obs_, obs_columns_]

    # --5 Get the predicitions [all weightings] on the current fold data []
    curr_preds <- get_predicition(Forest = Forest, 
                                  testdata = curr_fold_test_data)
    
    # --6 Exrtact the predicted classes for the different weighting techniques
    predicted_no_weight <- c(predicted_no_weight, 
                             curr_preds[grep("No", names(curr_preds))])
    
    predicted_f1_weight <- c(predicted_f1_weight,
                             curr_preds[grep("F1", names(curr_preds))])
    
    predicted_acc_weight <- c(predicted_acc_weight,
                              curr_preds[grep("Acc", names(curr_preds))])
    # --7 Extract the true response
    true_reponse <- c(true_reponse, 
                      as.character(curr_fold_test_data$outcome_Y))
    
  }
  
  # 2-3-7-3 Compare the predicted classes with the true response
  no_weight  <- eval_predicitons(as.factor(predicted_no_weight), 
                                 as.factor(true_reponse))
  acc_weight <- eval_predicitons(as.factor(predicted_acc_weight), 
                                 as.factor(true_reponse))
  f1_weight  <- eval_predicitons(as.factor(predicted_f1_weight), 
                                 as.factor(true_reponse))
  
  # 2-3-7-4 Add the results of the evaluation to the 'foldwise_res' list!
  foldwise_res[[i]] <- list("no_weight"  = no_weight,
                            "acc_weight" = acc_weight,
                            "f1_weight"  = f1_weight)
}

# 2-4 Save the results of the CV
save(foldwise_res, file = "./docs/CV_Res/REAL/FoldWise_Approach.R")
