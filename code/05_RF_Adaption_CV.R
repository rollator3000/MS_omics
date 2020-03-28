"Script to CrossValidate the Random Forest adaption of Roman Hornung!

 For each fold [set of observations w/ the same observed features], a seperate RF
 is trained [--> one seperate RF for each fold then]. For predicitons with these 
 foldwise-fitted-RFs it might be, that the foldwise grown trees need to be adjusted
 - depending on the observed features in testdata! Foldwise fitted trees need to be
 pruned, when they use a split variable, that is not avaible in the testset!
 ----- Pruning ----- 
  1. Select a foldwise fitted RandomForest.
  2. For each tree the RF consists of it is checked whether any of these trees uses
     a variable for splitting that is not avaible in the testdata.
  3. Each tree that uses a splitvariable not existent in the testset needs to be 
     pruned [cut off tree, before it splits w/ this variable]
  4. If a tree was pruned:
        - at its first split it can not be used for predicitons anymore!
        - anywhere else than the first split variable,  it can still be used for 
          predictions! Pass a test-obs. down the tree until it reaches a terminal-/
          pruned-node. 
          --> Prediciton equals disribution of the response in terminal-/ pruned-node!

For a final prediction by a RF the predicitons from the different foldwise-fitted-RFs 
are averaged. To get the predicition from a single foldwise-fitted-RF, get the predicition
from every single tree [if not pruned at 1st splitvariable] and average these! 

!!! Attention with running this code on Windows !!!
  - Paralelization leads to errors on Windows! 
    --> Avoid this, by:   - replacing '%dopar%' with '%do%'
                                         OR
                          - run 'registerDoParallel(cores = 1)'
"
# Load Functions, Classes & Librarys!
source("./code/04_simpleRF_adaption.R")
library(pROC)
library(assertthat)
library(checkmate)
library(doParallel)
library(e1071)

detectCores()
registerDoParallel(cores = 8)

load_CV_data        <- function(path) {
  "Load the subsetted, test-train splitted data, with blockwise missingness 
   induced already into the train split!
   This function does a detailed check on the data in 'path', so that all 
   approaches can deal with the data!
  
  Args:
    path (str) : Path to the data we want for CV!
                 Path must point to a list, that consits of two more lists 
                 'data' & 'block_names'!
                    - $data (list) must contain k entrances of 'train' & 'test'
                                   where each is a dataframe!
                    - $block_names (list) must contain at least 3 names!
  Return:
    A list, filled with 'data' & 'block_names' that can be used for CV!
  "
  # [0] Check Inputs  ----------------------------------------------------------
  # 0-1 'path' of type string with '.RData' inside!
  assert_string(path, fixed = ".RData")
  
  # [1] Load the Data & check it  ----------------------------------------------
  # 1-1 Load DF & put data to 'data_'-variable
  data_ <- load(path)
  data_ <- eval(as.symbol(data_))
  
  # 1-2 Check that 'data_' has what we need!
  # 1-2-1 List with 2 entrances
  assert_list(data_, len = 2)
  
  # 1-2-2 'data' and 'block_names' as entrances!
  if (!('data' %in% names(data_))) stop("'path' should lead to list w/ 'data' entry")
  if (!('block_names' %in% names(data_))) stop("'path' should lead to list w/ 'block_names' entry")
  
  # 1-2-3 'data' & 'block_names' should be lists aswell! 
  assert_list(data_$data, min.len = 2)
  assert_list(data_$block_names, min.len = 3)
  
  # 1-2-4 'data' must contain 'train' & 'test' as list names!
  res <- sapply(seq_len(length(data_$data)), function(x) {
    ("train" %in% names(data_$data[[x]])) & ('test' %in% names(data_$data[[x]]))
  })
  
  if (any(!res)) stop("'data' section doesn't contain 'test' & 'train' in each split")
  
  # 1-2-5 Contained 'train' & 'test' in 'data' need to be dataframes!
  res <- sapply(seq_len(length(data_$data)), function(x) {
    test_data_frame(data_$data[[x]]$train) & test_data_frame(data_$data[[x]]$test)
  })
  
  if (any(!res)) stop("test' & 'train' in 'data' are not only dataframes!")
  
  # 1-2-6 'block_names' must contain at least 3 blocknames!
  if (length(names(data_$block_names)) < 3) stop("block_names consist less than 3 names!")
  
  
  # [2] Return the data  -------------------------------------------------------
  return(data_)
}
mcc_metric          <- function(conf_matrix) {
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
copy_forrest        <- function(Forest) {
  "Funciton to copy all trees in a forest!
   Needed, as our trees can be pruned when doing predicitons!
   --> Need to copy the fitted forest, before using it to do predicitons 
       on different testsets!
       
    Args:
      - Forest (list) : list filled lists filled with objects of the class 'tree'
                        e.g. [tree1, tree2, tree3, tree4], where 
                             tree1, ..., tree4 are also lists filled with objects
                             of class 'tree' 
    Return:
      - Forest (list) : Same object as passed, but deeply copied!
  "
  # [0] Check Arguments  -------------------------------------------------------
  # 0-1 Check that all elemts are of class 'Tree'
  for (i in 1:length(Forest)) {
    classes_in_list <- sapply(Forest[[i]], function(x) "Tree" %in% class(x))
    if (any(classes_in_list)) {
      msg <- paste("List", i, "in Forest contains at least one object not of class 'Tree'")
      stop(msg)
    }
  }
  
  # [1] Copy the trees!  -------------------------------------------------------
  Forest_copy <- list() 
  for (i in 1:length(Forest)) {
    assign(paste0("treescopy", as.character(i)), lapply(Forest[[i]], FUN = function(x) x$copy()))
    Forest_copy[[i]] <- eval(parse(text = paste0("treescopy", as.character(i))))
  }
  
  # [2] Return deeply copied Forest  -------------------------------------------
  return(Forest_copy)
}
get_split_vars      <- function(Forest_) {
  "Get the used split variables from a foldwise fitted Forest!
  
  Args:
    - Forest_ (list) : List filled with objects of class 'Trees'
    
  Return:
    - vector of the used split variables & how often they have been 
      selected in the 'Forest_'
  "
  # [0] Check Inputs  ----------------------------------------------------------
  # 0-1 Forest is a list filled with elements of class 'Tree'
  assert_list(Forest_)
  tree_classes <- sapply(1:length(Forest_), function(x) {
    grepl("Tree", class(Forest_[[x]]))
  })
  
  if (!all(tree_classes)) {
    stop("Forest contains a Object not of class 'Tree'")
  }
  
  # 0-2 All trees must have at least 1 element in 'split_varIDs'
  wrong_trees <- unlist(lapply(Forest_, FUN = function(x) {
    length(x$child_nodeIDs) == 0
  })) 
  
  if (any(wrong_trees)) stop("Not all Trees were grown correctly!")
  
  # [1] Extract the used splitVariables  ---------------------------------------
  # 1-1 Get the split variable IDs of each Tree!
  used_split_vars_all <- split_var_ids_all <- sapply(1:length(Forest_), function(x) {
    split_var_ids <- Forest_[[x]]$split_varIDs
    split_var_ids <- split_var_ids[-which(is.na(split_var_ids))]
  })
  
  # 1-2 Unlist the result!
  used_split_vars_all <- unlist(used_split_vars_all)
  
  # 1-3 Convert the Variable IDs to featurenames!
  # 1-3-1 Extract all possible feature names the trees were fit on!
  poss_fea_names <- colnames(Forest_[[1]]$data$data)
  
  # 1-3-2 Convert VariableIDs to names!
  used_split_var_names <- sapply(used_split_vars_all, function(x) {
    poss_fea_names[x]
  })
  
  # [2] Return the names of the used split variables  --------------------------
  return(used_split_var_names)
}
get_oob_weight_metric     <- function(trees) {
  "Calculate OOB Metrirc ['F1' & 'Acc'] of a list of trees! 
   For this we go  through all OOB Predictions and obtain aggregated predicitons
   from all trees, that have the same observation as out-of-bag!
   In case the metrics are 'NA' [not defined] they are repalced by 0 [worst value]
    
    Args: 
      - trees (list)        : list filled w/ objects of class 'Tree'
      
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
all_trees_grown_correctly <- function(trees) {
  "Check, whether 'trees', were grown correctly & if not grow these trees again,
   as long, as they are grown correctly! 
      --> Growning not correctly: No 'childNodeIDs', no split variables etc...
  
   Args:
      trees (list) : list filled with object of the class 'Tree'! 
                     For each object in there check, whether it was grown 
                     correctly (has childNodeIDs) - if not grown it again!
                     
   Return: 
      list of trees, where all of these trees were grown correctly
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
eval_predicitons          <- function(predicitons, reference, used_split_vars) {
  "Calculate the Metrics for a given list of prediciotns and a given list of 
   true responses - response must be categorical!
   
   
  Args:
    - predicitons (vector) : Vector filled with integers - predicitons from 
                             a model! Must be of class factor!
    - reference   (vector) : Vector filled with integers - true responses!
                             Must be of class factor!
    - used_split_vars (vector) : vevctor full of strings - containing the used
                                 split variables
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
  
  # 0-3 'used_split_vars' must be vector filled with strings!
  assert_vector(used_split_vars)
  
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
              "MCC"         = mcc,
              "Selected_Vars" = used_split_vars)
  
  # 6-1 If the F1-Score/ Precision/ Recall is NA, then we set it to 0 [-1 for MCC]! 
  if (is.na(res$F1))             res$F1             <- 0
  if (is.na(res$Precision))      res$Precision      <- 0
  if (is.na(res$Recall))         res$Recall         <- 0
  if (is.na(res$MCC))            res$MCC            <- -1
  if (is.na(res$Pos_Pred_Value)) res$Pos_Pred_Value <- 0
  if (is.na(res$Neg_Pred_Value)) res$Neg_Pred_Value <- 0
  
  return(res)
}
do_evaluation             <- function(Forest, testdata) {
  " Get the aggregated predicition from all trees for each of the observations
    in testdata! Evaluate the aggregated predicitons & return metrics!
  
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
  assert_data_frame(testdata, any.missing = F, min.rows = 1)
  
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
  tree_preds_all <- foreach(i = 1:length(Forest)) %dopar% { # par
    
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
  if (length(Forest) > 1) {
    
    # Get the ACC & F1 weights for each of the foldwise fitted Forests!
    tree_weights <- foreach(l = seq_len(length(Forest))) %dopar% { # par
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
  
  
  # [5] Get Metrics for the current setting - each weighting seperate!  --------
  # 5-1 Firstly extract the used split variables from the fitted trees!
  used_split_vars <- unlist(sapply(1:length(Forest), function(x) {
    get_split_vars(Forest_ = Forest[[x]])
  }))
    
  # 5-2 Get the performance with NO_weighting!
  no_weighting <- eval_predicitons(predicitons = preds_class_0_no_weight,
                                   reference = testdata[,1],
                                   used_split_vars = used_split_vars)
  
  # 5-3 Get the performance with ACC_weighting!
  acc_weighting <- eval_predicitons(predicitons = preds_class_0_acc_weight,
                                    reference = testdata[,1],
                                    used_split_vars = used_split_vars)
  
  # 5-4 Get the performance with F1_weighting!
  f1_weighting <- eval_predicitons(predicitons = preds_class_0_f1_weight,
                                   reference = testdata[,1],
                                   used_split_vars = used_split_vars)
  
  # [6] Return the Results  ----------------------------------------------------
  res_all <- list("no_weighting"  = no_weighting,
                  "acc_weighting" = acc_weighting,
                  "f1_weighting"  = f1_weighting)
  return(as.vector(res_all))
}

do_CV_5_blocks <- function(path = "data/processed/RH_subsetted_12345/missingness_1234/BLCA_1.RData",
                           num_trees = 300, mtry = NULL, min_node_size = 5,
                           unorderd_factors = "ignore") {
  "CrossValidate the Foldwise-Approach when the Traindata has blockwise 
   missingness according to the scenario 1, 2 or 3!
   
   'path' must lead to a list with 2 entrances: 'data' & 'block_names'
     - 'data' is a list filled with 'k' test-train-splits
        --> k-fold-Validation on this test-train-splits!
     - 'block_names' is a list filled with the names of the single blocks 
        & must be ['A', 'B', 'C', 'D', 'clin_block']!
        (Attention: In Scenario2 the order is different, but this is wanted!)
      
   Based on the 'k' test-train-splits in 'data', we will fit a foldwise RFs to the
   train data (that has blockwise missingness in it). Then we ensemble the 
   predicitons from foldwise fitted RFs to single predicitons & rate these with 
   Accuracy, Precision, Specifity, F1-Socre,...
   
   The TestingSituations are different, as we can test the models on fully 
   observed testdata, on testdata w/ 1 missing block, etc...
    --> Results is list with all results from the k test-train splits for all 
        possible testsituations - 20 in total!
   
  Args:
      - path (char)         : path to data w/ blockwise missingness for the CV.
                              Must end in '1.RData', '2.RData' or '3.RData'
                              --> List with 2 entrances: 'data' & 'block_names'
                                  - 'data' consitis of 'k' test-train-splits, 
                                     where train has missingness induced and the 
                                     test-set is fully observed!
                                  - 'block_names' contains all colnames of the 
                                     different blocks!
      - num_trees (int)     : Amount of trees, we shall grow on each foldwise 
                              fitted random Forest must be int > 10!
      - mtry (int)          : Amount of split-variables we try, when looking for 
                              a split variable! 
                                - If 'NULL': mtry = sqrt(p)
      - min_node_size (int) : Amount of Observations a node must at least 
                              contain, so the model keeps on trying to split 
                              them! Must be int >= 1.
      - unorderd_factors (chr) : How to handle non numeric features!
                                ['ignore', 'order_once', 'order_split', 'partition']
  Return:
      - list filled w/:
        * 'res_all' [the CV Results on different testsets]:
           >> A = CNV-Block, B = RNA-Block, C = Mutation-Block, D = Mirna-Block
            - full    : CV Results for each fold on the fully observed testdata!
            - miss1_A : CV Results for each fold on the testdata, w/ missing cnv
                        block [or in scenario2, what was sampled to be block 'A']
            - miss1_B : CV Results for each fold on the testdata, w/ missing rna 
                        block [or in scenario2, what was sampled to be block 'B'] 
                .
                .
            - miss2_AC: CV Results for each fold on the testdata, w/ missing cnv 
                        & mutation block! [or in scenario2, what was sampled to 
                                           be block 'A' & 'C'] 
                .
                .
            - single_D: CV-Results for each fold on the testdata w/ minra block 
                        only! [or in scenario2, what was sampled to be block 'D']
                         
        * 'settings' [settings used to do the CV - all arguments!]
            - datapath, response, mtry, time for CV,.... 
  "
  # [0] Check Inputs  ----------------------------------------------------------
  # 0-1 path must be numeric and have '1.RData' | '2.RData' | '3.RData' in it!
  assert_string(path)
  if (!grepl("1.RData", path) & !grepl("2.RData", path) & !grepl("3.RData", path)) {
    stop("'path' must end in '1.RData' | '2.RData' | '3.RData'")
  }
  
  # 0-2 'num_trees', 'min_node_size' must be intgers > 0  & 'mtry' only if NOT NULL
  assert_int(num_trees, lower = 10)
  assert_int(min_node_size, lower = 1)
  if (!is.null(mtry)) assert_int(mtry)
  
  # 0-3 unorderd factors must be a legit value
  if (!(unorderd_factors %in% c("ignore", "order_once", "order_split", "partition"))) {
    stop("Unknown value for argument 'unordered_factors'")
  }
  
  # [1] Prepare CV  ------------------------------------------------------------
  # 1-1 Load CV-Data [already splitted - data checked in 'load_CV_data' itself]
  curr_data <- load_CV_data(path = path)
  
  # 1-1-1 Must contain 'A', 'B', 'C', 'D' & 'clin_block' as block_names
  corr_block_names <- ("A" %in% names(curr_data$block_names) & 
                       "B" %in% names(curr_data$block_names) &
                       "C" %in% names(curr_data$block_names) & 
                       "D" %in% names(curr_data$block_names) &
                       "clin_block" %in% names(curr_data$block_names))
  
  if (!corr_block_names) stop("'path' lead to a file without 'A', 'B', 'C', 'D' & 'clin_block' as blocknames!")
  
  # 1-2 Create empty lists to store results in!
  # 1-2-1 Full TestSet
  full <- list()
  # 1-2-2 TestSet with 1 missing omics-block
  miss1_A <- list(); miss1_B <- list(); miss1_C <- list(); miss1_D <- list()
  # 1-2-3 TestSet with 2 missing omics-blocks
  miss2_CD <- list(); miss2_BD <- list(); miss2_BC <- list(); miss2_AD <- list()
  miss2_AC <- list(); miss2_AB <- list()
  # 1-2-4 TestSet with 3 missing omics-blocks
  miss3_ABC <- list(); miss3_ABD <- list(); miss3_ACD <- list(); miss3_BCD <- list()
  # 1-2-5 Single BlockTestSet [4 missing blocks!]
  single_A <- list(); single_B <- list(); single_CL <- list(); single_C <- list(); single_D <- list()
  
  # 1-3 Get the amount of test-train splits in data 
  k_splits <- length(curr_data$data)
  
  # 1-4 Start the Timer, so we know how long CV took
  start_time <- Sys.time()
  
  # [2] Start the CV and loop over all test-train splits in data  --------------
  for (i in seq_len(k_splits)) {
    
    # 2-1 Current Fold Status
    print(paste0("FOLD ", i, "/", k_splits, " -------------------------------"))
    
    # 2-2 Extract the test and train set from 'curr_data'
    train <- curr_data$data[[i]]$train
    test  <- curr_data$data[[i]]$test
    
    # 2-3 Get the Observations that belong to the same fold [same feature space]
    # 2-3-1 Get for each obs. the index of the observed feas
    observed_feas <- foreach(x = seq_len(nrow(train)), .combine = 'c') %dopar% {
      paste0(which(!(is.na(train[x,]))), collapse = "_")
    }
    
    # 2-3-2 Keep the unique observed feas [equals the different folds]
    #       That we use to assign obs. to the differnt folds!
    observed_folds <- unique(observed_feas)
    print(paste0("Found ", length(observed_folds), " unique folds!"))
    
    # 2-4 Train foldwise RFs for each fold seperatly! For this loop over all folds
    #     [observed_folds contains all unique features observed for diff folds]!
    #     --> Results in a Forest of length 'lenght(observed_folds)' & each 
    #         entrance consits of 'num_trees' foldwise fitted trees!
    Forest <- list()
    Forest <- foreach(j_ = 1:length(observed_folds)) %do% {
      
      fold_ = observed_folds[j_]
      # 2-4-1 Get all Obs. with the feture space as in 'fold_'
      fold_obs_ <- which(observed_feas == fold_)
      
      # 2-4-2 Get all the indeces of the columns that were observed with this fold!
      obs_columns_ <- as.numeric(unlist(strsplit(fold_, split = "_")))
      
      # 2-4-2 Get all Trainpoints from the obs. w/ same features + 
      #       only keep observed features of these!
      #       --> fully observed subdata!
      curr_fold_train_data <- train[fold_obs_, obs_columns_]
      
      # 2-4-3 Fit a RF on this fully observed (fold-)subdata!
      # 2-4-3-1 Define formula
      response    <- colnames(train)[1]
      formula_all <- as.formula(paste(response, " ~ ."))
      
      # 2-4-3-2 Define the foldwise RF and fit it on the foldwise data
      fold_RF <- simpleRF(formula           = formula_all, 
                          data              = curr_fold_train_data, 
                          num_trees         = num_trees, 
                          mtry              = mtry, 
                          min_node_size     = as.integer(min_node_size),
                          replace           = TRUE,  # always TRUE, as we need OOB!
                          splitrule         = NULL,  # always NULL!
                          unordered_factors = unorderd_factors)
      
      fold_RF <- lapply(fold_RF, function(x) {
        x$grow(replace = TRUE)
        x
      })
      
      # 2-4-3-3 Check that all trees were grown correctly
      #         --> none w/o 'child_node_ID' after that!
      fold_RF <- all_trees_grown_correctly(fold_RF)
      
      return(fold_RF)
    }
    
    # 2-5 Evaluate the RF on the different Testsets! Fromfull TestSet w/o any 
    #     missing blocks to TestSets w/ only one observed Block!
    # 2-5-1 Full TestSet!
    print("Evaluation full TestSet -------------------------------------------")
    curr_Forest <- copy_forrest(Forest)
    full[[i]]   <- do_evaluation(Forest = curr_Forest, testdata = test)
    
    # 2-5-2 TestSet with 1 missing block!
    print("Evaluation TestSet w/ 1 missing omics block------------------------")
    curr_Forest  <- copy_forrest(Forest)
    miss1_A[[i]] <- do_evaluation(Forest = curr_Forest,
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$A)])
    curr_Forest  <- copy_forrest(Forest)
    miss1_B[[i]] <- do_evaluation(Forest = curr_Forest,
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$B)])
    curr_Forest  <- copy_forrest(Forest)
    miss1_C[[i]] <- do_evaluation(Forest = curr_Forest,
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$C)])
    curr_Forest  <- copy_forrest(Forest)
    miss1_D[[i]] <- do_evaluation(Forest = curr_Forest,
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$D)])
    
    # 2-5-3 TestSet with 2 missing blocks!
    print("Evaluation TestSet w/ 2 missing omics blocks-----------------------")
    curr_Forest   <- copy_forrest(Forest)
    miss2_CD[[i]] <- do_evaluation(Forest = curr_Forest,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                 curr_data$block_names$D))])
    curr_Forest   <- copy_forrest(Forest)
    miss2_BD[[i]] <- do_evaluation(Forest = curr_Forest,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                 curr_data$block_names$D))])
    curr_Forest   <- copy_forrest(Forest)
    miss2_BC[[i]] <- do_evaluation(Forest = curr_Forest,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                 curr_data$block_names$C))])
    curr_Forest   <- copy_forrest(Forest)
    miss2_AD[[i]] <- do_evaluation(Forest = curr_Forest,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$D))])
    curr_Forest   <- copy_forrest(Forest)
    miss2_AC[[i]] <- do_evaluation(Forest = curr_Forest,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                 curr_data$block_names$A))])
    
    curr_Forest   <- copy_forrest(Forest)
    miss2_AB[[i]] <- do_evaluation(Forest = curr_Forest,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$B))])
    # 2-5-4 Testset with 3 missing blocks!
    print("Evaluation TestSet w/ 3 missing omics blocks-----------------------")
    curr_Forest    <- copy_forrest(Forest)
    miss3_ABC[[i]] <- do_evaluation(Forest = curr_Forest,
                                    testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                  curr_data$block_names$B,
                                                                                  curr_data$block_names$C))])
    curr_Forest    <- copy_forrest(Forest)
    miss3_ACD[[i]] <- do_evaluation(Forest = curr_Forest,
                                    testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                  curr_data$block_names$C,
                                                                                  curr_data$block_names$D))])
    curr_Forest    <- copy_forrest(Forest)
    miss3_ABD[[i]] <- do_evaluation(Forest = curr_Forest,
                                    testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                  curr_data$block_names$A,
                                                                                  curr_data$block_names$D))])
    
    curr_Forest    <- copy_forrest(Forest)
    miss3_BCD[[i]] <- do_evaluation(Forest = curr_Forest,
                                    testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                  curr_data$block_names$C,
                                                                                  curr_data$block_names$D))])
    # 2-5-5 Evaluation on single Block Testdata
    print("Evaluation TestSet w/ only 1 observed Block -----------------------")
    curr_Forest   <- copy_forrest(Forest)
    single_A[[i]] <- do_evaluation(Forest = curr_Forest,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                 curr_data$block_names$C,
                                                                                 curr_data$block_names$D,
                                                                                 curr_data$block_names$clin_block))])
    curr_Forest   <- copy_forrest(Forest)
    single_B[[i]] <- do_evaluation(Forest = curr_Forest,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$C,
                                                                                 curr_data$block_names$D,
                                                                                 curr_data$block_names$clin_block))])
    
    curr_Forest   <- copy_forrest(Forest)
    single_C[[i]] <- do_evaluation(Forest = curr_Forest,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$B,
                                                                                 curr_data$block_names$D,
                                                                                 curr_data$block_names$clin_block))])
    curr_Forest   <- copy_forrest(Forest)
    single_D[[i]] <- do_evaluation(Forest = curr_Forest,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$B,
                                                                                 curr_data$block_names$C,
                                                                                 curr_data$block_names$clin_block))])
    curr_Forest     <- copy_forrest(Forest)
    single_CL[[i]]  <- do_evaluation(Forest = curr_Forest,
                                     testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                   curr_data$block_names$B,
                                                                                   curr_data$block_names$D,
                                                                                   curr_data$block_names$C))])
  }
  # 2-6 Take the time difference for the k fold CV!
  time_for_CV <- difftime(Sys.time(), start_time, units = "mins")
  
  # [3] Return the metric & settings of the fitting!  --------------------------
  # 3-1 Collect all CV Results in a list!
  res_all <- list("full" = full,
                  "miss1_A"   = miss1_A,   "miss1_B"   = miss1_B,
                  "miss1_C"   = miss1_C,   "miss1_D"   = miss1_D,
                  "miss2_CD"  = miss2_CD,  "miss2_BD"  = miss2_BD,
                  "miss2_BC"  = miss2_BC,  "miss2_AD"  = miss2_AD,
                  "miss2_AC"  = miss2_AC,  "miss2_AB"  = miss2_AB,
                  "miss3_ABC" = miss3_ABC, "miss3_ABD" = miss3_ABD,
                  "miss3_ACD" = miss3_ACD, "miss3_BCD" = miss3_BCD,
                  "single_A"  = single_A,  "single_B"  = single_B,
                  "single_C"  = single_C,  "single_D"  = single_D,
                  "single_CL" = single_CL)
  
  # 3-2 Collect the Settings, used to do the CV!
  settings <- list("data_path"     = path,
                   "num_folds"     = k_splits,
                   "response"      = response,
                   "num_trees"     = num_trees,
                   "mtry"          = mtry, 
                   "min_node_size" = min_node_size,
                   "unorderd_factors" = unorderd_factors,
                   "time_for_CV"      = time_for_CV)
  
  # 3-3 Return both lists!
  return(list("res_all"  = res_all, 
              "settings" = settings))
}

do_CV_2_blocks <- function(path = "data/processed/RH_subsetted_12345/missingness_1234/BLCA_4.RData", 
                           num_trees = 300, mtry = NULL, min_node_size = 5,
                           unorderd_factors = "ignore") {
  "CrossValidate the Approach when the Traindata has blockwise missingness
   according to scenario 4 --> 'path' must end in '4.RData'!
   
   'path' must lead to a list with 2 entrances: 'data'  &  'block_names'
       - 'data' is a list filled with 'k' test-train-splits
          --> k-fold-Validation on this test-train-splits!
       - 'block_names' is a list filled with the names of the single blocks 
          & must be 'A', 'B' & 'clin_block' for this scenario!
      
   Based on the 'k' test-train-splits in 'data', we will fit foldwise RFs to the
   train data (that has blockwise missingness in it). Then we ensemble the 
   predicitons from foldwise fitted RFs to single predicitons & rate these wth 
   Accuracy, Precision, Specifity, F1-Socre,...
   
   The TestingSituations are different, as we can test the models on fully 
   observed testdata, on testdata w/ 1 missing block, etc...
    --> Results is list with all results from the k test-train splits for all 
        possible testsituations - 6 in total!
   
  Args:
      - path (char)         : path to data w/ blockwise missingness for the CV.
                              Must end in '4.RData'
                              --> List with 2 entrances: 'data' & 'block_names'
                                  - 'data' consitis of 'k' test-train-splits, 
                                     where train has missingness induced and the 
                                     test-set is fully observed!
                                  - 'block_names' contains all colnames of the 
                                     different blocks!
      - num_trees (int)     : Amount of trees, we shall grow on each foldwise 
                              fitted random Forest must be int > 10!
      - mtry (int)          : Amount of split-variables we try, when looking for 
                              a split variable! 
                                - If 'NULL': mtry = sqrt(p)
      - min_node_size (int) : Amount of Observations a node must at least 
                              contain, so the model keeps on trying to split 
                              them! Must be int >= 1.
      - unorderd_factors (chr) : How to handle non numeric features!
                                ['ignore', 'order_once', 'order_split', 'partition']
  Return:
      - list filled w/:
        * 'res_all' [the CV Results on different testsets]:
           >> A = CNV-Block, B = RNA-Block, C = Mutation-Block, D = Mirna-Block
            - full    : CV Results for each fold on the fully observed testdata!
            - miss1_A : CV Results for each fold on the testdata, w/ missing A
                        block [actually 2 blocks randomly thrown together]
            - miss1_B : CV Results for each fold on the testdata, w/ missing B 
                        block [actually 2 blocks randomly thrown together]
            - single_A: CV Results for each fold on the testdata, w/ only block
                        A as feature [actually 2 blocks randomly thrown together]
            - single_B: CV Results for each fold on the testdata, w/ only block
                        A as feature [actually 2 blocks randomly thrown together]
            - single_clin: CV Results for each fold on the testdata, w/ only block
                           clinical as feature
        * 'settings' [settings used to do the CV - all arguments!]
            - datapath, response, mtry, time for CV, ...
  "
  # [0] Check Inputs  ----------------------------------------------------------
  # 0-1 path must be numeric and have '1.RData' in it!
  assert_string(path, fixed = "4.RData")
  
  # 0-2 'num_trees', 'min_node_size' must be intgers > 0  & 'mtry' only if NOT NULL
  assert_int(num_trees, lower = 10)
  assert_int(min_node_size, lower = 1)
  if (!is.null(mtry)) assert_int(mtry)
  
  # 0-3 unorderd factors must be a legit value
  if (!(unorderd_factors %in% c("ignore", "order_once", "order_split", "partition"))) {
    stop("Unknown value for argument 'unordered_factors'")
  }
  
  # [1] Prepare CV  ------------------------------------------------------------
  # 1-1 Load CV-Data [already splitted - data checked in 'load_CV_data' itself]
  curr_data <- load_CV_data(path = path)
  
  # 1-1-1 Must contain 'A', 'B', $ 'clin_block' as block_names
  corr_block_names <- ("A" %in% names(curr_data$block_names) & 
                         "B" %in% names(curr_data$block_names) &
                         "clin_block" %in% names(curr_data$block_names))
  
  if (!corr_block_names) stop("'path' lead to a file without 'A', 'B' & 'clin_block' as blocknames!")
  
  # 1-2 Create empty lists to store results in!
  # 1-2-1 Full TestSet
  full <- list()
  # 1-2-2 TestSet with 1 missing omics-block
  miss1_A <- list(); miss1_B <- list(); miss1_C <- list(); miss1_D <- list()
  # 1-2-3 TestSet with 2 missing omics-blocks
  miss2_CD <- list(); miss2_BD <- list(); miss2_BC <- list(); miss2_AD <- list()
  miss2_AC <- list(); miss2_AB <- list()
  # 1-2-4 TestSet with 3 missing omics-blocks
  miss3_ABC <- list(); miss3_ABD <- list(); miss3_ACD <- list(); miss3_BCD <- list()
  # 1-2-5 Single BlockTestSet [4 missing blocks!]
  single_A <- list(); single_B <- list(); single_CL <- list(); single_C <- list(); single_D <- list()
  
  # 1-3 Get the amount of test-train splits in data 
  k_splits <- length(curr_data$data)
  
  # 1-4 Start the Time so we know how long CV took
  start_time <- Sys.time()
  
  # [2] Start the CV and loop over all test-train splits in data  --------------
  for (i in seq_len(k_splits)) {
    
    # 2-1 Current Fold Status
    print(paste0("FOLD ", i, "/", k_splits, " -------------------------------"))
    
    # 2-2 Extract the test and train set from 'curr_data'
    train <- curr_data$data[[i]]$train
    test  <- curr_data$data[[i]]$test
    
    # 2-3 Get the Observations that belong to the same fold [same feature space]
    # 2-3-1 Get for each obs. the index of the observed feas
    observed_feas <- foreach(x = seq_len(nrow(train)), .combine = 'c') %dopar% {
      paste0(which(!(is.na(train[x,]))), collapse = "_")
    }
    
    # 2-3-2 Keep the unique observed feas [equals the different folds]
    #       That we use to assign obs. to the differnt folds!
    observed_folds <- unique(observed_feas)
    print(paste0("Found ", length(observed_folds), " unique folds!"))
    
    # 2-4 Train foldwise RFs for each fold seperatly! For this loop over all folds
    #     [observed_folds contains all unique features observed for diff folds]!
    #     --> Results in a Forest of length 'lenght(observed_folds)' & each 
    #         entrance consits of 'num_trees' foldwise fitted trees!
    Forest <- list()
    Forest <- foreach(j_ = 1:length(observed_folds)) %do% {
      
      fold_ = observed_folds[j_]
      # 2-4-1 Get all Obs. with the feture space as in 'fold_'
      fold_obs_ <- which(observed_feas == fold_)
      
      # 2-4-2 Get all the indeces of the columns that were observed with this fold!
      obs_columns_ <- as.numeric(unlist(strsplit(fold_, split = "_")))
      
      # 2-4-2 Get all Trainpoints from the obs. w/ same features + 
      #       only keep observed features of these!
      #       --> fully observed subdata!
      curr_fold_train_data <- train[fold_obs_, obs_columns_]
      
      # 2-4-3 Fit a RF on this fully observed (fold-)subdata!
      # 2-4-3-1 Define formula
      response    <- colnames(train)[1]
      formula_all <- as.formula(paste(response, " ~ ."))
      
      # 2-4-3-2 Define the foldwise RF and fit it on the foldwise data
      fold_RF <- simpleRF(formula           = formula_all, 
                          data              = curr_fold_train_data, 
                          num_trees         = num_trees, 
                          mtry              = mtry, 
                          min_node_size     = as.integer(min_node_size),
                          replace           = TRUE,  # always TRUE, as we need OOB!
                          splitrule         = NULL,  # always NULL!
                          unordered_factors = unorderd_factors)
      
      fold_RF <- lapply(fold_RF, function(x) {
        x$grow(replace = TRUE)
        x
      })
      
      # 2-4-3-3 Check that all trees were grown correctly
      #         --> none w/o 'child_node_ID' after that!
      fold_RF <- all_trees_grown_correctly(fold_RF)
      
      return(fold_RF)
    }
    
    # 2-5 Evaluate the RF on the different Testsets! Fromfull TestSet w/o any 
    #     missing blocks to TestSets w/ only one observed Block!
    # 2-5-1 Full TestSet!
    print("Evaluation full TestSet -------------------------------------------")
    curr_Forest <- copy_forrest(Forest)
    full[[i]]   <- do_evaluation(Forest = curr_Forest, testdata = test)
    
    # 2-5-2 TestSet with 1 missing block!
    print("Evaluation TestSet w/ 1 missing omics block------------------------")
    curr_Forest  <- copy_forrest(Forest)
    miss1_A[[i]] <- do_evaluation(Forest = curr_Forest,
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$A)])
    curr_Forest  <- copy_forrest(Forest)
    miss1_B[[i]] <- do_evaluation(Forest = curr_Forest,  
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$B)])
    
    # 2-5-3 Evaluation on single Block Testdata
    print("Evaluation TestSet w/ only 1 observed Block -----------------------")
    curr_Forest   <- copy_forrest(Forest)
    single_A[[i]] <- do_evaluation(Forest = curr_Forest, 
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                 curr_data$block_names$clin_block))])
    curr_Forest   <- copy_forrest(Forest)
    single_B[[i]] <- do_evaluation(Forest = curr_Forest, 
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$clin_block))])
    
    curr_Forest     <- copy_forrest(Forest)
    single_CL[[i]]  <- do_evaluation(Forest = curr_Forest, 
                                     testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                   curr_data$block_names$B))])
  }
  # 2-6 Take the time difference for the k fold CV!
  time_for_CV <- difftime(Sys.time(), start_time, units = "mins")
  
  # [3] Return the metric & settings of the fitting!  --------------------------
  # 3-1 Collect all CV Results in a list!
  res_all <- list("full" = full, "miss1_A" = miss1_A, "miss1_B" = miss1_B,
                  "single_A" = single_A, "single_B" = single_B, "single_CL" = single_CL)
  
  # 3-2 Collect the Settings, used to do the CV!
  settings <- list("data_path"     = path,
                   "num_folds"     = k_splits,
                   "response"      = response,
                   "num_trees"     = num_trees,
                   "mtry"          = mtry, 
                   "min_node_size" = min_node_size,
                   "unorderd_factors" = unorderd_factors,
                   "time_for_CV"      = time_for_CV)
  
  # 3-3 Return both lists!
  return(list("res_all"  = res_all, 
              "settings" = settings))
}

# Run Main                                                                  ----
"Run the CV for all DFs from the missForest paper ------------------------------"
DFs_w_gender <- c("COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC","LGG", "BLCA",
                  "LUAD", "LUSC", "PAAD", "SARC", "SKCM", "STAD")

# ----- Situation 1
for (DF in DFs_w_gender) {
  
  print(paste0("----- Situation 1 for DF: '", DF, "' -----"))
  
  # Create the path for the current DF
  curr_path <- paste0("data/processed/RH_subsetted_12345/missingness_1234/", DF, "_1.RData")
  
  sit1 <- do_CV_5_blocks(path = curr_path, num_trees = 300, mtry = NULL, 
                         min_node_size = 5, unorderd_factors = "ignore")
  save(sit1, file = paste0("./docs/CV_Res/gender/Roman_final_subsets/setting1/", DF, ".RData"))
}

# ----- Situation 2
for (DF in DFs_w_gender) {
  
  print(paste0("----- Situation 2 for DF: '", DF, "' -----"))
  
  # Create the path for the current DF
  curr_path <- paste0("data/processed/RH_subsetted_12345/missingness_1234/", DF, "_2.RData")
  
  sit2 <- do_CV_5_blocks(path = curr_path, num_trees = 300, mtry = NULL, 
                         min_node_size = 5, unorderd_factors = "ignore")
  save(sit2, file = paste0("./docs/CV_Res/gender/Roman_final_subsets/setting2/", DF, ".RData"))
}

# ----- Situation 3
for (DF in DFs_w_gender) {
  
  print(paste0("----- Situation 3 for DF: '", DF, "' -----"))
  
  # Create the path for the current DF
  curr_path <- paste0("data/processed/RH_subsetted_12345/missingness_1234/", DF, "_3.RData")
  
  sit3 <- do_CV_5_blocks(path = curr_path, num_trees = 300, mtry = NULL, 
                         min_node_size = 5, unorderd_factors = "ignore")
  save(sit3, file = paste0("./docs/CV_Res/gender/Roman_final_subsets/setting3/", DF, ".RData"))
}

# ----- Situation 4
for (DF in DFs_w_gender) {
  
  print(paste0("----- Situation 4 for DF: '", DF, "' -----"))
  
  # Create the path for the current DF
  curr_path <- paste0("data/processed/RH_subsetted_12345/missingness_1234/", DF, "_4.RData")
  
  sit4 <- do_CV_5_blocks(path = curr_path, num_trees = 300, mtry = NULL, 
                         min_node_size = 5, unorderd_factors = "ignore")
  save(sit4, file = paste0("./docs/CV_Res/gender/Roman_final_subsets/setting4/", DF, ".RData"))
}