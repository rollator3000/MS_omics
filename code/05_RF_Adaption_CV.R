"Script to CrossValidate the Random Forest adaption of Roman Hornung!
 Here for each fold [set of observations w/ the same observed feature space], 
 we train a seperate RF [--> seperate RFs for each fold then]. When we want to 
 do predicitons with these foldwise fitted RFs, we need to adjust the foldwise
 trees, depending on the observed features in the testdata! Foldwise fitted trees
 need to pruned, when they use a split variable, that is not avaible in testset!
 Pruning: 1. Select a single foldwise RandomForest.
          2. For each tree of the RF, check whether any tree uses a variable for
             splitting that is not avaible in the testdata
          3. Each tree, that uses splitvariables, that are not existent in test, 
             need to be pruned [cut off tree, before it splits w/ this variable]
          4. If a tree was pruned at its first (inital) split it can not be used
             for prediciton anymore!
          5. If a tree was pruned anywhere else than the first split variable, 
             it can still be used for predictions, by passing the test obs. down
             the tree until it reaches a terminal/ pruned node.
             --> Class Prediciton is the disribution in the terminal/ pruned node!
  
To obtain a final prediction by a foldwise RF we combine the predicitons from the
different Trees, that were not pruned in the first splitvar., and combine them!
For this we obtain the predicition from each single tree and combine these!
"
# Load Funcitons, Classes, librarys & set the WD!
setwd("C:/Users/kuche_000/Desktop/MS-Thesis/")
source("./code/04_simpleRF_adaption.R")
library(pROC)
library(assertthat)

load_CV_data              <- function(path) {
  "Load the subsetted, test-train splitted data, with blockwise missingness 
   induced already into the train split!
  
  Args:
    path (str) : Path to the data we want for CV!
                 Path must point to a list, that consits of two more lists 
                 'data' & 'block_names'! 
                 These two lists are further checked in this function!
                 - $data (list) must contain n entrances of 'train' & 'test'
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
  
  # 1-2-2 'data' and 'block_names'
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
all_trees_grown_correctly <- function(trees) {
  "Check, whether 'trees', were grown correctly & if not grow these trees again,
   as long, as they are grown correctly! 
      --> Growning not correctly: No 'childNodeIDs' 
  
   Args:
      trees (list) : list filled with object of the class 'Tree'! 
                     For each object in there check, whether it was grown 
                     correctly (has childNodeIDs) - if not grown it again!
                     
   Return: 
      list of trees, where all of these trees were grown correctly
      --> each tree has at least one 'childNodeID'
  "
  # [0] Check Inputs -----------------------------------------------------------
  # 0-1 All Objects in 'trees' of class 'Tree'
  trees_classes <- sapply(trees, function(x) class(x))
  if (any(!grepl("Tree", trees_classes))) {
    stop("Not all Objects in 'trees' are of class 'Tree'")
  }
  
  # [1] Get the entrance of the objects, that miss child node IDs --------------
  wrong_trees  <- unlist(lapply(1:length(trees), 
                                FUN = function(x) {
                                  if (length(trees[[x]]$child_nodeIDs) == 0) x
                                }))
  
  # [2] Regrow the trees, that were not grown correctly ------------------------
  #     If there are any trees not grown correctly, grow them again until all
  #     of the trees were grown correctly!
  while (length(wrong_trees) > 0) {
    
    # grow the errours trees again
    trees[wrong_trees] <- mclapply(trees[wrong_trees], 
                                   function(x) {
                                     x$grow(replace = TRUE)
                                     x
                                   }, mc.cores = 1)
    
    # check whether any of the trees is not grown correctly!
    wrong_trees  <- unlist(lapply(1:length(trees), 
                                  FUN = function(x) {
                                    if (length(trees[[x]]$child_nodeIDs) == 0) x
                                  }))
  }
  
  # [3] Return the correclty grown trees ---------------------------------------
  return(trees)
}
copy_forrest              <- function(Forest) {
  "Funciton to copy all trees in a forest!
   This is needed, as we change our trees when doing predicitons on data,
   where we need to prune our trees! 
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
  # [0] Check Arguments --------------------------------------------------------
  # 0-1 Check that all elemts are of class 'Tree'
  for (i in 1:length(Forest)) {
    classes_in_list <- sapply(Forest[[i]], function(x) "Tree" %in% class(x))
    if (any(classes_in_list)) {
      msg <- paste("List", i, "in Forest contains at least one object not of class 'Tree'")
      stop(msg)
    }
  }
  
  # [1] Copy the trees! --------------------------------------------------------
  Forest_copy <- list() 
  for (i in 1:length(Forest)) {
    assign(paste0("treescopy", as.character(i)), lapply(Forest[[i]], FUN = function(x) x$copy()))
    Forest_copy[[i]] <- eval(parse(text = paste0("treescopy", as.character(i))))
  }
  
  # [2] Return deeply copied Forest --------------------------------------------
  return(Forest_copy)
}
get_oob_weight_metric     <- function(trees, weight_metric) {
  "Calculate OOB Metrirc ['F1' or 'Acc'] of a list of trees! For this we go 
   through all OOB Predictions and obtain aggregated predicitons from all 
   trees, that have the same observation as out-of-bag!
    
    Args: 
      - trees (list)        : list filled w/ objects of class 'Tree'
      - weight_metric (chr) : When assigning weights to the different predictions
                              which metric to use to calc the weight?!
                              [- must be 'Acc' or 'F1']
    Return:
      - average oob-Acc // oob-F1 over all trees in 'trees'-list!
  "
  # [0] Check Input ------------------------------------------------------------
  # 0-1 Make sure 'trees' is a list filled with 'Trees'
  assert_list(trees, min.len = 1, any.missing = FALSE)
  if (any(sapply(trees, function(x) 'Tree' %in% class(x)))) {
    stop("not all elements in 'trees' are of class 'Tree'")
  }
  
  # 0-2 Check 'weight_metric' to be one character w/ F1 or Acc
  #     check that it is character
  assert_character(weight_metric, len = 1)
  
  #     check it has a valid value!
  if (!(weight_metric %in% c("Acc", "F1"))) stop("'weight_metric' must be 'Acc' or 'F1'!")
  
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
  
  if (weight_metric == "Acc") return(conf_mat_oob$overall["Accuracy"])
  if (weight_metric == "F1")  return(conf_mat_oob$byClass["F1"])
}
mcc_metric                <- function(conf_matrix) {
  "Function to calculate the MCC [Matthews correlation coefficient] Metric
   --> only for binary cases! If the Conf_Matrix has more than 2 classes 
       it will return NULL instead of the MCC!
       
    Definition of the Metric:
      MCC takes into account true and false positives and negatives and is 
      generally regarded as a balanced measure which can be used even if the 
      classes are of very different sizes.
       
    Args: 
      - conf_matrix (confusionMatrix) : Confusion Matrix created with the 
                                        'caret'-Package!
    Return: 
      Matthews correlation coefficient [1 is best, -1 is worst!]
  "
  # [0] Check Inputs -----------------------------------------------------------
  # 0-1 check amount of classes:
  if (nrow(conf_matrix$table) != 2) {
    warning("Can not calc the MCC-Metric! Return NULL")
    return(NULL)
  }
  
  # 0-2 Check class of conf_matrix
  if (class(conf_matrix) != "confusionMatrix") {
    stop("conf_matrix not of class 'confusionMatrix'")
  }
  
  # [1] Calc the Score ---------------------------------------------------------
  TP <- conf_matrix$table[1,1]
  TN <- conf_matrix$table[2,2]
  FP <- conf_matrix$table[1,2]
  FN <- conf_matrix$table[2,1]
  
  mcc_num <- (TP*TN - FP*FN)
  mcc_den <- 
    as.double((TP+FP))*as.double((TP+FN))*as.double((TN+FP))*as.double((TN+FN))
  
  mcc_final <- mcc_num/sqrt(mcc_den)
  
  # [2] Return it --------------------------------------------------------------
  return(mcc_final)
}
do_evaluation             <- function(Forest, testdata, weighted, weight_metric) {
  " Get the aggregated predicition from all trees for each of the observations
    in testdata! Evaluate the aggregated predicitons & return metrics!
  
     Args:
      - Forest (list)         : list filled with the objects of class 'Tree'
      - testdata (data.frame) : Testdata we want to get predicitons for!
                                Must conatain the response we learned the
                                Forest on!
      - weighted (boolean)    : Shall the predicitons be weighted when aggregted
                                --> Blocks with higher 'weight_metric' 
                                    [Acc, F1, ...], recieve higher weight!
      - weight_metric (chr)   : Which Metric to use, to weight the predicitons
                                of the different trees in 'Forest'
                                Has to be 'F1' or 'Acc'
                                [will be ignored, when 'weighted' = FALSE]
     Return:
      - list w/ metrics [accuracy, f1, mcc, roc, ...] 
        if F-1 Score is NA [as precision or recall = 0], we set it to 0 [not NA]
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
  
  # 0-2 Check testdata for DF & weighted to be boolean
  assert_data_frame(testdata, any.missing = F, min.rows = 1)
  assert_logical(weighted)
  
  # 0-3 Check weight_metric to be meaningful [if weighted is true at all]!
  if (weighted) {
    
    # check that it is character
    assert_character(weight_metric, len = 1)
    
    # check it has a valid value!
    if (!(weight_metric %in% c("Acc", "F1"))) {
      stop("'weight_metric' must be 'Acc' or 'F1'!")
    }
  }
  
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
  #     fitted RandomForests [4 predicitons per testobs.]      
  tree_preds_all <- list()
  not_usable     <- c()
  for (i in 1:length(Forest)) {
    
    # save the predictions as 'treeX_pred'
    tree_preds_all[[i]] <- get_pruned_prediction(trees = Forest[[i]], 
                                                 test_set = tree_testsets[[i]])
    
    # check whether all of the predicted values are NA 
    # [if so the tree can not be used when ensembling prediciton results!]
    not_usable <- c(not_usable, all(is.na(tree_preds_all[[i]]$Class)))
  }
  
  # 2-1 Check that there are still trees existing, else not preds possible!
  if (all(not_usable)) {
    print("None of the trees are usable for predictions!")
    return("No Predictions possible as all trees are pruned at the 1. splitvar!")
  }
  
  # 2-2 Remove the SubRFs, that are not usable from Forest % tree_preds_all
  if (any(not_usable)) {
    Forest         <- Forest[-c(which(not_usable))]
    tree_preds_all <- tree_preds_all[-c(which(not_usable))]
  }
  
  # [3] If we want to create weighted ensemble of the predicitons, we need to
  #     calc the OOB-Accuracy/ OOB-F1--Score per 'trees' [foldwise RF] and use 
  #     these as weights when assembeling the predicitons!
  #     [lower ACC --> lower weight | lower F1-Score --> lower weight]
  #     ATTENTION: No need to calc weights when there is only one foldwise RF left!
  if (weighted & length(Forest) > 1) {
    tree_weights <- c()
    
    # Think about putting this in parallel!
    tree_weights <- sapply(seq_len(length(Forest)), function(l) {
      get_oob_weight_metric(trees = Forest[[l]], 
                            weight_metric = weight_metric)
    })
    
    # Norm the weights
    tree_weights <- tree_weights /  sum(tree_weights)
    
  } else {
    tree_weights <- rep(1, times = length(Forest))
  }
  
  # [4] Aggregate Predictions from the different trees!
  # 4-1 Get the probabilities of all test obs to be of class '0'
  all_forrest_preds_probs_class_0 <- sapply(1:nrow(testdata), FUN = function(x) {
    
    # Get a probability prediciton from each [still usable] tree!
    preds_all <- sapply(seq_len(length(tree_preds_all)), function(i) {
      tree_preds_all[[i]]$Probs[[x]][1]
    })
    
    # Combine the preditions of the different trees!
    prob_class0 <- weighted.mean(preds_all, w = tree_weights, na.rm = TRUE)
    prob_class0
  })
  
  # 4-2 Convert the probabilities to class predicitons and convert it to 'class'
  all_forrest_preds_class <- ifelse(all_forrest_preds_probs_class_0 >= 0.5, 0, 1)
  all_forrest_preds_class <- factor(all_forrest_preds_class, 
                                    levels = levels(Forest[[1]][[1]]$data$data[,1]))
  
  # [5] Get Metrics for the current setting!
  # 5-1 Confusion Matrix - from which we can calc/ extract most metrics
  confmat <- caret::confusionMatrix(data      = all_forrest_preds_class, 
                                    reference = testdata[,1])
  
  # 5-2 Are under the ROC Curve
  roc1 <- pROC::auc(pROC::roc(testdata[,1], all_forrest_preds_probs_class_0, 
                              levels = levels(Forest[[1]][[1]]$data$data[,1]), 
                              direction = "<"))
  
  roc2 <- pROC::auc(pROC::roc(testdata[,1], all_forrest_preds_probs_class_0, 
                              levels = levels(Forest[[1]][[1]]$data$data[,1]), 
                              direction = ">"))
  
  # 5-3 MCC Matthews correlation coefficient [only for binary cases!]
  mcc <- mcc_metric(conf_matrix = confmat)
  
  # [6] Create a list to collect the results!
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
  
  # 6-1 If the F1-Score/ Precision/ Recall is NA, then we set it to 0 [-1 for MCC]! 
  if (is.na(res$F1))             res$F1             <- 0
  if (is.na(res$Precision))      res$Precision      <- 0
  if (is.na(res$Recall))         res$Recall         <- 0
  if (is.na(res$MCC))            res$MCC            <- -1
  if (is.na(res$Pos_Pred_Value)) res$Pos_Pred_Value <- 0
  if (is.na(res$Neg_Pred_Value)) res$Neg_Pred_Value <- 0
  
  return(as.vector(res))
}

do_CV_1 <- function(path = "data/external/Dr_Hornung/subsetted_12345/missingness_1234/BLCA_1.RData",
                    weighted = TRUE, weight_metric = "Acc", 
                    num_trees = 10, mtry = NULL, min_node_size = NULL,
                    unorderd_factors = "ignore") {
  "CrossValidate the Approach when the Traindata has blockwise missingness
   according to scenario 1 --> 'path' must end in '1.RData', else error!
   
   path must lead to a list with 2 entrances: 'data'  &  'block_names'
   - 'data' is a list filled with 'k' test-train-splits
      --> k-fold-Validation on this test-train-splits!
   - 'block_names' is a list filled with the names of the single blocks 
      & must be ['A', 'B', 'C', 'D', 'clin_block'] for this scenario!
      [A = cnv; B = rna; C = mirna; D = mutation]
      
   Based on the 'k' test-train-splits in 'data', we will fit foldwise RFs to the
   train data (that has blockwise missingness in it). 
   Then we ensemble the predicitons from foldwise fitted RFs to single  
   predicitons & rate these w/ Accuracy, Precision, Specifity, F1-Socre,...
   
   The TestingSituations are different, as we can test the models on fully 
   observed testdata, on testdata w/ 1 missing block, etc...
   
   Args:
      - path (char)         : Path to the data that is already split to test train
                              Must end in '1.RData'
      - weighted (bool)     : When ensembling the prediciton from the fodlwise
                              RFs shall we weight the predictions by the OOB 
                              performance of the foldwise fitted RFs? The higher
                              the OOB-Metric for a foldwise fitted RF, the higher
                              its contribution to the prediction!
      - weight_metric (chr) : Which metric to use to assigning weights to the 
                              different predictions? 
                                - must be 'Acc' or 'F1' 
                                - if weighted = FALSE, it will be ignored!
      - num_trees (int)     : Amount of trees, we shall grow on each foldwise
                              fitted random Forest!
      - mtry (int)          : Amount of split-variables we try, when looking for 
                              a split variable! 
                                - If 'NULL': mtry = sqrt(p)
      - min_node_size (int) : Amount of Observations a node must at least 
                              contain, so the model keeps on trying to split 
                              them!  
                                - If 'NULL: It is automatically set to 10 in the
                                            'simpleRF()' function!
      - unorderd_factors (chr) : How to handle non numeric features!
                                ['ignore', 'order_once', 'order_split', 'partition']
    Return:
      - list filled w/:
        * 'res_all' [the CV Results on different testsets]:
           >> A = CNV-Block, B = RNA-Block, C = Mutation-Block, D = Mirna-Block
            - full    : CV Results for each fold on the fully observed testdata!
            - miss1_A : CV Results for each fold on the testdata, w/ missing cnv block!
            - miss1_B : CV Results for each fold on the testdata, w/ missing rna block!
                .
                . --> done for all possible permutations
                .
            - miss2_AC: CV Results for each fold on the testdata, w/ missing cnv & mutation block!
                .
                .
            - single_D:  CV-Results for each fold on the testdata w/ minra block only!
                         
        * 'settings' [settings used to do the CV - all arguments!]
            - datapath, seed, response, mtry, time for CV,.... 
  "
  # [0] Check Inputs  ----------------------------------------------------------
  # 0-0 mtry, min_node_size & num_trees are all checked within 'simpleRF()'
  # 0-1 path must be numeric and have '1.RData' in it!
  assert_string(path, fixed = "1.RData")
  
  # 0-2 weighted must be a single boolean
  assert_logical(weighted, len = 1)
  
  # 0-3 unorderd factors must be a legit value
  if (!(unorderd_factors %in% c("ignore", "order_once", "order_split", "partition"))) {
    stop("Unknown value for argument 'unordered_factors'")
  }
  
  # 0-4 Check weight_metric to be meaningful [if weighted is true at all]!
  if (weighted) {
    
    # 0-4-1 Check that it is character
    assert_character(weight_metric, len = 1)
    
    # 0-4-2 Check it has a valid value!
    if (!(weight_metric %in% c("Acc", "F1"))) {
      stop("'weight_metric' must be 'Acc' or 'F1'!")
    }
  }
  
  # [1] Prepare CV  ------------------------------------------------------------
  # 1-1 Load CV-Data [already splitted - data checked in 'load_CV_data' itself]
  curr_data <- load_CV_data(path = path)
  
  # 1-1-1 Must contain 'A', 'B', 'C' & 'D' as block_names
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
  
  # [2] Start the CV and loop over all test-train splits in data  --------------
  for (i in seq_len(k_splits)) {
    
    # 2-1 Current Fold Status
    print(paste0("FOLD ", i, "/", k_splits, " -------------------------------"))
    
    # 2-2 Extract the test and train set from 'curr_data'
    train <- curr_data$data[[i]]$train
    test  <- curr_data$data[[i]]$test
    
    # 2-3 Get the Observations that belong to the same fold [same feature space]
    # 2-3-1 Get for each obs. the index of the observed feas
    observed_feas <- sapply(seq_len(nrow(train)), function(x) {
      paste0(which(!(is.na(train[x,]))), collapse = "_")
    })
    
    # 2-3-2 Keep the unique observed feas [equals the different folds]
    #       That we use to assign obs. to the differnt folds!
    observed_folds <- unique(observed_feas)
    
    print(paste0("Found ", length(observed_folds), " unique folds!"))
    
    # 2-4 Train foldwise RFs for each fold seperatly! For this loop over all folds
    #     [observed_folds contains all unique features observed for diff folds]!
    #     --> Results in a Forest of length 'lenght(observed_folds)' & each 
    #         entrance consits of 'num_trees' foldwise fitted trees!
    Forest <- list(); i_  <- 1; start_time <- Sys.time()
    
    for (fold_ in observed_folds) {
      
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
                          min_node_size     = min_node_size,
                          replace           = TRUE,  # always TRUE, as we need OOB!
                          splitrule         = NULL,  # always NULL!
                          unordered_factors = unorderd_factors)
      
      print(paste0("Fit FoldWise RF on current fold: ", i_))
      
      fold_RF <- mclapply(fold_RF, function(x) {
        x$grow(replace = TRUE)
        x
      }, mc.cores = 1)
      
      # 2-4-3-3 Check that all trees were grown correctly
      #         --> none w/o 'child_node_ID' after that!
      fold_RF <- all_trees_grown_correctly(fold_RF)
      
      # 2-4-4 Add the fitted tree to the Forest!
      # 2-4-4-1 Copy the tree first, so we do no overwriting!
      fold_RF_copy <- sapply(fold_RF, function(x) x$copy())
      
      # 2-4-4-2 Remove the fold_RF and add the copied fold_RF_copy to 'Forest'
      rm(fold_RF)
      Forest[[i_]] <- fold_RF_copy
      
      # 2-4-5 Count up i_ - used to fill Forest list!
      i_  = i_ + 1
    }
    
    # 2-5 Evaluate the RF on the different Testsets! Fromfull TestSet w/o any 
    #     missing blocks to TestSets w/ only one observed Block!
    # 2-5-1 Full TestSet!
    print("Evaluation full TestSet -------------------------------------------")
    curr_Forest <- copy_forrest(Forest)
    full[[i]]   <- do_evaluation(Forest = curr_Forest, testdata = test, 
                                 weighted = weighted, weight_metric = weight_metric)
    
    # 2-5-2 TestSet with 1 missing block!
    print("Evaluation TestSet w/ 1 missing omics block------------------------")
    curr_Forest  <- copy_forrest(Forest)
    miss1_A[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                  weight_metric = weight_metric,
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$A)])
    curr_Forest  <- copy_forrest(Forest)
    miss1_B[[i]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                  weight_metric = weight_metric,
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$B)])
    curr_Forest  <- copy_forrest(Forest)
    miss1_C[[i]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                  weight_metric = weight_metric,
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$C)])
    curr_Forest  <- copy_forrest(Forest)
    miss1_D[[i]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                  weight_metric = weight_metric,
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$D)])
    
    # 2-5-3 TestSet with 2 missing blocks!
    print("Evaluation TestSet w/ 2 missing omics blocks-----------------------")
    curr_Forest   <- copy_forrest(Forest)
    miss2_CD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                 curr_data$block_names$D))])
    curr_Forest   <- copy_forrest(Forest)
    miss2_BD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                 curr_data$block_names$D))])
    curr_Forest   <- copy_forrest(Forest)
    miss2_BC[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted, 
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                 curr_data$block_names$c))])
    curr_Forest   <- copy_forrest(Forest)
    miss2_AD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$D))])
    curr_Forest   <- copy_forrest(Forest)
    miss2_AC[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                 curr_data$block_names$A))])
    
    curr_Forest   <- copy_forrest(Forest)
    miss2_AB[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$B))])
    # 2-5-4 Testset with 3 missing blocks!
    print("Evaluation TestSet w/ 3 missing omics blocks-----------------------")
    curr_Forest    <- copy_forrest(Forest)
    miss3_ABC[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                    weight_metric = weight_metric,
                                    testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                  curr_data$block_names$B,
                                                                                  curr_data$block_names$C))])
    curr_Forest    <- copy_forrest(Forest)
    miss3_ACD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                    weight_metric = weight_metric,
                                    testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                  curr_data$block_names$C,
                                                                                  curr_data$block_names$D))])
    curr_Forest    <- copy_forrest(Forest)
    miss3_ABD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted, 
                                    weight_metric = weight_metric,
                                    testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                  curr_data$block_names$A,
                                                                                  curr_data$block_names$D))])
    
    curr_Forest    <- copy_forrest(Forest)
    miss3_BCD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                    weight_metric = weight_metric,
                                    testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                  curr_data$block_names$C,
                                                                                  curr_data$block_names$D))])
    # 2-5-5 Evaluation on single Block Testdata
    print("Evaluation TestSet w/ only 1 observed Block -----------------------")
    curr_Forest   <- copy_forrest(Forest)
    single_A[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                 curr_data$block_names$C,
                                                                                 curr_data$block_names$D,
                                                                                 curr_data$block_names$clin_block))])
    curr_Forest   <- copy_forrest(Forest)
    single_B[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$C,
                                                                                 curr_data$block_names$D,
                                                                                 curr_data$block_names$clin_block))])
    
    curr_Forest   <- copy_forrest(Forest)
    single_C[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$B,
                                                                                 curr_data$block_names$D,
                                                                                 curr_data$block_names$clin_block))])
    curr_Forest   <- copy_forrest(Forest)
    single_D[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$B,
                                                                                 curr_data$block_names$C,
                                                                                 curr_data$block_names$clin_block))])
    curr_Forest     <- copy_forrest(Forest)
    single_CL[[i]]  <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                     weight_metric = weight_metric,
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
  
  # 4-2 Collect the Settings, used to do the CV!
  settings <- list("data_path"     = path,
                   "num_folds"     = k_splits,
                   "response"      = response,
                   "weighted"      = weighted,
                   "weight_metric" = weight_metric,
                   "num_trees"     = num_trees,
                   "mtry"          = mtry, 
                   "min_node_size" = min_node_size,
                   "unorderd_factors" = unorderd_factors,
                   "time_for_CV"      = time_for_CV)
  
  # 4-3 Return both lists!
  return(list("res_all"  = res_all, 
              "settings" = settings))
}

do_CV_2 <- function(path = "data/external/Dr_Hornung/subsetted_12345/missingness_1234/BLCA_2.RData",
                    weighted = TRUE, weight_metric = "Acc", 
                    num_trees = 10, mtry = NULL, min_node_size = NULL,
                    unorderd_factors = "ignore") {
  "CrossValidate the Approach when the Traindata has blockwise missingness
   according to scenario 2 --> 'path' must end in '2.RData', else error!
   
   path must lead to a list with 2 entrances: 'data'  &  'block_names'
   - 'data' is a list filled with 'k' test-train-splits
      --> k-fold-Validation on this test-train-splits!
   - 'block_names' is a list filled with the names of the single blocks 
      & must be ['A', 'B', 'C', 'D', 'clin_block'] for this scenario!
      
   Based on the 'k' test-train-splits in 'data', we will fit foldwise RFs to the
   train data (that has blockwise missingness in it). 
   Then we ensemble the predicitons from foldwise fitted RFs to single  
   predicitons & rate these w/ Accuracy, Precision, Specifity, F1-Socre,...
   
   The TestingSituations are different, as we can test the models on fully 
   observed testdata, on testdata w/ 1 missing block, etc...
   
   Args:
      - path (char)         : Path to the data that is already split to test train
                              Must end in '2.RData'
      - weighted (bool)     : When ensembling the prediciton from the fodlwise
                              RFs shall we weight the predictions by the OOB 
                              performance of the foldwise fitted RFs? The higher
                              the OOB-Metric for a foldwise fitted RF, the higher
                              its contribution to the prediction!
      - weight_metric (chr) : Which metric to use to assigning weights to the 
                              different predictions? 
                                - must be 'Acc' or 'F1' 
                                - if weighted = FALSE, it will be ignored!
      - num_trees (int)     : Amount of trees, we shall grow on each foldwise
                              fitted random Forest!
      - mtry (int)          : Amount of split-variables we try, when looking for 
                              a split variable! 
                                - If 'NULL': mtry = sqrt(p)
      - min_node_size (int) : Amount of Observations a node must at least 
                              contain, so the model keeps on trying to split 
                              them!  
                                - If 'NULL: It is automatically set to 10 in the
                                            'simpleRF()' function!
      - unorderd_factors (chr) : How to handle non numeric features!
                                ['ignore', 'order_once', 'order_split', 'partition']
    Return:
      - list filled w/:
        * 'res_all' [the CV Results on different testsets]:
           >> A = CNV-Block, B = RNA-Block, C = Mutation-Block, D = Mirna-Block
            - full    : CV Results for each fold on the fully observed testdata!
            - miss1_A : CV Results for each fold on the testdata, w/ missing A block!
            - miss1_B : CV Results for each fold on the testdata, w/ missing B block!
                .
                . --> done for all possible permutations
                .
            - miss2_AC: CV Results for each fold on the testdata, w/ missing A & C block!
                .
                .
            - single_D:  CV-Results for each fold on the testdata w/ D block only!
                         
        * 'settings' [settings used to do the CV - all arguments!]
            - datapath, seed, response, mtry, time for CV,.... 
  "
  # [0] Check Inputs  ----------------------------------------------------------
  # 0-0 mtry, min_node_size & num_trees are all checked within 'simpleRF()'
  # 0-1 path must be numeric and have '1.RData' in it!
  assert_string(path, fixed = "2.RData")
  
  # 0-2 weighted must be a single boolean
  assert_logical(weighted, len = 1)
  
  # 0-3 unorderd factors must be a legit value
  if (!(unorderd_factors %in% c("ignore", "order_once", "order_split", "partition"))) {
    stop("Unknown value for argument 'unordered_factors'")
  }
  
  # 0-4 Check weight_metric to be meaningful [if weighted is true at all]!
  if (weighted) {
    
    # 0-4-1 Check that it is character
    assert_character(weight_metric, len = 1)
    
    # 0-4-2 Check it has a valid value!
    if (!(weight_metric %in% c("Acc", "F1"))) {
      stop("'weight_metric' must be 'Acc' or 'F1'!")
    }
  }
  
  # [1] Prepare CV  ------------------------------------------------------------
  # 1-1 Load CV-Data [already splitted - data checked in 'load_CV_data' itself]
  curr_data <- load_CV_data(path = path)
  
  # 1-1-1 Must contain 'A', 'B', 'C' & 'D' as block_names
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
  
  # [2] Start the CV and loop over all test-train splits in data  --------------
  for (i in seq_len(k_splits)) {
    
    # 2-1 Current Fold Status
    print(paste0("FOLD ", i, "/", k_splits, " -------------------------------"))
    
    # 2-2 Extract the test and train set from 'curr_data'
    train <- curr_data$data[[i]]$train
    test  <- curr_data$data[[i]]$test
    
    # 2-3 Get the Observations that belong to the same fold [same feature space]
    # 2-3-1 Get for each obs. the index of the observed feas
    observed_feas <- sapply(seq_len(nrow(train)), function(x) {
      paste0(which(!(is.na(train[x,]))), collapse = "_")
    })
    
    # 2-3-2 Keep the unique observed feas [equals the different folds]
    #       That we use to assign obs. to the differnt folds!
    observed_folds <- unique(observed_feas)
    
    print(paste0("Found ", length(observed_folds), " unique folds!"))
    
    # 2-4 Train foldwise RFs for each fold seperatly! For this loop over all folds
    #     [observed_folds contains all unique features observed for diff folds]!
    #     --> Results in a Forest of length 'lenght(observed_folds)' & each 
    #         entrance consits of 'num_trees' foldwise fitted trees!
    Forest <- list(); i_  <- 1; start_time <- Sys.time()
    
    for (fold_ in observed_folds) {
      
      # 2-4-1 Get all Obs. with the feture space as in 'fold_'
      fold_obs_ <- which(observed_feas == fold_)
      
      # 2-4-2 Get all the indeces of the columns that were observed with this fold!
      obs_columns_ <- as.numeric(unlist(strsplit(fold_, split = "_")))
      
      # 2-4-2 Get all Trainpoints from the obs. w/ same features + 
      #       only keep observed features of these!
      #       --> fully observed subdata!
      curr_fold_train_data <- train[fold_obs_, obs_columns_]
      
      cat(paste("Ncol fold", i_, ":", ncol(curr_fold_train_data), "\n"))
      
      # 2-4-3 Fit a RF on this fully observed (fold-)subdata!
      # 2-4-3-1 Define formula
      response    <- colnames(train)[1]
      formula_all <- as.formula(paste(response, " ~ ."))
      
      # 2-4-3-2 Define the foldwise RF and fit it on the foldwise data
      fold_RF <- simpleRF(formula           = formula_all, 
                          data              = curr_fold_train_data, 
                          num_trees         = num_trees, 
                          mtry              = mtry, 
                          min_node_size     = min_node_size,
                          replace           = TRUE,  # always TRUE, as we need OOB!
                          splitrule         = NULL,  # always NULL!
                          unordered_factors = unorderd_factors)
      
      print(paste0("Fit FoldWise RF on current fold: ", i_))
      
      fold_RF <- mclapply(fold_RF, function(x) {
        x$grow(replace = TRUE)
        x
      }, mc.cores = 1)
      
      # 2-4-3-3 Check that all trees were grown correctly
      #         --> none w/o 'child_node_ID' after that!
      fold_RF <- all_trees_grown_correctly(fold_RF)
      
      # 2-4-4 Add the fitted tree to the Forest!
      # 2-4-4-1 Copy the tree first, so we do no overwriting!
      fold_RF_copy <- sapply(fold_RF, function(x) x$copy())
      
      # 2-4-4-2 Remove the fold_RF and add the copied fold_RF_copy to 'Forest'
      rm(fold_RF)
      Forest[[i_]] <- fold_RF_copy
      
      # 2-4-5 Count up i_ - used to fill Forest list!
      i_  = i_ + 1
    }
    
    # 2-5 Evaluate the RF on the different Testsets! Fromfull TestSet w/o any 
    #     missing blocks to TestSets w/ only one observed Block!
    # 2-5-1 Full TestSet!
    print("Evaluation full TestSet -------------------------------------------")
    curr_Forest <- copy_forrest(Forest)
    full[[i]]   <- do_evaluation(Forest = curr_Forest, testdata = test, 
                                 weighted = weighted, weight_metric = weight_metric)
    
    # 2-5-2 TestSet with 1 missing block!
    print("Evaluation TestSet w/ 1 missing omics block------------------------")
    curr_Forest  <- copy_forrest(Forest)
    miss1_A[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                  weight_metric = weight_metric,
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$A)])
    curr_Forest  <- copy_forrest(Forest)
    miss1_B[[i]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                  weight_metric = weight_metric,
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$B)])
    curr_Forest  <- copy_forrest(Forest)
    miss1_C[[i]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                  weight_metric = weight_metric,
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$C)])
    curr_Forest  <- copy_forrest(Forest)
    miss1_D[[i]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                  weight_metric = weight_metric,
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$D)])
    
    # 2-5-3 TestSet with 2 missing blocks!
    print("Evaluation TestSet w/ 2 missing omics blocks-----------------------")
    curr_Forest   <- copy_forrest(Forest)
    miss2_CD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                 curr_data$block_names$D))])
    curr_Forest   <- copy_forrest(Forest)
    miss2_BD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                 curr_data$block_names$D))])
    curr_Forest   <- copy_forrest(Forest)
    miss2_BC[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted, 
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                 curr_data$block_names$c))])
    curr_Forest   <- copy_forrest(Forest)
    miss2_AD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$D))])
    curr_Forest   <- copy_forrest(Forest)
    miss2_AC[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                 curr_data$block_names$A))])
    
    curr_Forest   <- copy_forrest(Forest)
    miss2_AB[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$B))])
    # 2-5-4 Testset with 3 missing blocks!
    print("Evaluation TestSet w/ 3 missing omics blocks-----------------------")
    curr_Forest    <- copy_forrest(Forest)
    miss3_ABC[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                    weight_metric = weight_metric,
                                    testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                  curr_data$block_names$B,
                                                                                  curr_data$block_names$C))])
    curr_Forest    <- copy_forrest(Forest)
    miss3_ACD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                    weight_metric = weight_metric,
                                    testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                  curr_data$block_names$C,
                                                                                  curr_data$block_names$D))])
    curr_Forest    <- copy_forrest(Forest)
    miss3_ABD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted, 
                                    weight_metric = weight_metric,
                                    testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                  curr_data$block_names$A,
                                                                                  curr_data$block_names$D))])
    
    curr_Forest    <- copy_forrest(Forest)
    miss3_BCD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                    weight_metric = weight_metric,
                                    testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                  curr_data$block_names$C,
                                                                                  curr_data$block_names$D))])
    # 2-5-5 Evaluation on single Block Testdata
    print("Evaluation TestSet w/ only 1 observed Block -----------------------")
    curr_Forest   <- copy_forrest(Forest)
    single_A[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                 curr_data$block_names$C,
                                                                                 curr_data$block_names$D,
                                                                                 curr_data$block_names$clin_block))])
    curr_Forest   <- copy_forrest(Forest)
    single_B[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$C,
                                                                                 curr_data$block_names$D,
                                                                                 curr_data$block_names$clin_block))])
    
    curr_Forest   <- copy_forrest(Forest)
    single_C[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$B,
                                                                                 curr_data$block_names$D,
                                                                                 curr_data$block_names$clin_block))])
    curr_Forest   <- copy_forrest(Forest)
    single_D[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$B,
                                                                                 curr_data$block_names$C,
                                                                                 curr_data$block_names$clin_block))])
    curr_Forest     <- copy_forrest(Forest)
    single_CL[[i]]  <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                     weight_metric = weight_metric,
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
  
  # 4-2 Collect the Settings, used to do the CV!
  settings <- list("data_path"     = path,
                   "num_folds"     = k_splits,
                   "response"      = response,
                   "weighted"      = weighted,
                   "weight_metric" = weight_metric,
                   "num_trees"     = num_trees,
                   "mtry"          = mtry, 
                   "min_node_size" = min_node_size,
                   "unorderd_factors" = unorderd_factors,
                   "time_for_CV"      = time_for_CV)
  
  # 4-3 Return both lists!
  return(list("res_all"  = res_all, 
              "settings" = settings))
}

path = "data/external/Dr_Hornung/subsetted_12345/missingness_1312/COAD_3.RData"
weighted = TRUE
weight_metric = "Acc"
num_trees = 10
mtry = NULL
min_node_size = NULL
unorderd_factors = "ignore"

do_CV_3 <- function(path = "data/external/Dr_Hornung/subsetted_12345/missingness_1234/BLCA_3.RData",
                    weighted = TRUE, weight_metric = "Acc", 
                    num_trees = 10, mtry = NULL, min_node_size = NULL,
                    unorderd_factors = "ignore") {
  "CrossValidate the Approach when the Traindata has blockwise missingness
   according to scenario 3 --> 'path' must end in '3.RData', else error!
   
   path must lead to a list with 2 entrances: 'data'  &  'block_names'
   - 'data' is a list filled with 'k' test-train-splits
      --> k-fold-Validation on this test-train-splits!
   - 'block_names' is a list filled with the names of the single blocks 
      & must be ['A', 'B', 'C', 'D', 'clin_block'] for this scenario!
      
   Based on the 'k' test-train-splits in 'data', we will fit foldwise RFs to the
   train data (that has blockwise missingness in it). 
   Then we ensemble the predicitons from foldwise fitted RFs to single  
   predicitons & rate these w/ Accuracy, Precision, Specifity, F1-Socre,...
   
   The TestingSituations are different, as we can test the models on fully 
   observed testdata, on testdata w/ 1 missing block, etc...
   
   Args:
      - path (char)         : Path to the data that is already split to test train
                              Must end in '3.RData'
      - weighted (bool)     : When ensembling the prediciton from the fodlwise
                              RFs shall we weight the predictions by the OOB 
                              performance of the foldwise fitted RFs? The higher
                              the OOB-Metric for a foldwise fitted RF, the higher
                              its contribution to the prediction!
      - weight_metric (chr) : Which metric to use to assigning weights to the 
                              different predictions? 
                                - must be 'Acc' or 'F1' 
                                - if weighted = FALSE, it will be ignored!
      - num_trees (int)     : Amount of trees, we shall grow on each foldwise
                              fitted random Forest!
      - mtry (int)          : Amount of split-variables we try, when looking for 
                              a split variable! 
                                - If 'NULL': mtry = sqrt(p)
      - min_node_size (int) : Amount of Observations a node must at least 
                              contain, so the model keeps on trying to split 
                              them!  
                                - If 'NULL: It is automatically set to 10 in the
                                            'simpleRF()' function!
      - unorderd_factors (chr) : How to handle non numeric features!
                                ['ignore', 'order_once', 'order_split', 'partition']
    Return:
      - list filled w/:
        * 'res_all' [the CV Results on different testsets]:
           >> A = CNV-Block, B = RNA-Block, C = Mutation-Block, D = Mirna-Block
            - full    : CV Results for each fold on the fully observed testdata!
            - miss1_A : CV Results for each fold on the testdata, w/ missing A block!
            - miss1_B : CV Results for each fold on the testdata, w/ missing B block!
                .
                . --> done for all possible permutations
                .
            - miss2_AC: CV Results for each fold on the testdata, w/ missing A & C block!
                .
                .
            - single_D:  CV-Results for each fold on the testdata w/ D block only!
                         
        * 'settings' [settings used to do the CV - all arguments!]
            - datapath, seed, response, mtry, time for CV,.... 
  "
  # [0] Check Inputs  ----------------------------------------------------------
  # 0-0 mtry, min_node_size & num_trees are all checked within 'simpleRF()'
  # 0-1 path must be numeric and have '1.RData' in it!
  assert_string(path, fixed = "3.RData")
  
  # 0-2 weighted must be a single boolean
  assert_logical(weighted, len = 1)
  
  # 0-3 unorderd factors must be a legit value
  if (!(unorderd_factors %in% c("ignore", "order_once", "order_split", "partition"))) {
    stop("Unknown value for argument 'unordered_factors'")
  }
  
  # 0-4 Check weight_metric to be meaningful [if weighted is true at all]!
  if (weighted) {
    
    # 0-4-1 Check that it is character
    assert_character(weight_metric, len = 1)
    
    # 0-4-2 Check it has a valid value!
    if (!(weight_metric %in% c("Acc", "F1"))) {
      stop("'weight_metric' must be 'Acc' or 'F1'!")
    }
  }
  
  # [1] Prepare CV  ------------------------------------------------------------
  # 1-1 Load CV-Data [already splitted - data checked in 'load_CV_data' itself]
  curr_data <- load_CV_data(path = path)
  
  # 1-1-1 Must contain 'A', 'B', 'C' & 'D' as block_names
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
  
  # [2] Start the CV and loop over all test-train splits in data  --------------
  for (i in seq_len(k_splits)) {
    
    # 2-1 Current Fold Status
    print(paste0("FOLD ", i, "/", k_splits, " -------------------------------"))
    
    # 2-2 Extract the test and train set from 'curr_data'
    train <- curr_data$data[[i]]$train
    test  <- curr_data$data[[i]]$test
    
    # 2-3 Get the Observations that belong to the same fold [same feature space]
    # 2-3-1 Get for each obs. the index of the observed feas
    observed_feas <- sapply(seq_len(nrow(train)), function(x) {
      paste0(which(!(is.na(train[x,]))), collapse = "_")
    })
    
    # 2-3-2 Keep the unique observed feas [equals the different folds]
    #       That we use to assign obs. to the differnt folds!
    observed_folds <- unique(observed_feas)
    
    print(paste0("Found ", length(observed_folds), " unique folds!"))
    
    # 2-4 Train foldwise RFs for each fold seperatly! For this loop over all folds
    #     [observed_folds contains all unique features observed for diff folds]!
    #     --> Results in a Forest of length 'lenght(observed_folds)' & each 
    #         entrance consits of 'num_trees' foldwise fitted trees!
    Forest <- list(); i_  <- 1; start_time <- Sys.time()
    
    for (fold_ in observed_folds) {
      
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
                          min_node_size     = min_node_size,
                          replace           = TRUE,  # always TRUE, as we need OOB!
                          splitrule         = NULL,  # always NULL!
                          unordered_factors = unorderd_factors)
      
      print(paste0("Fit FoldWise RF on current fold: ", i_))
      
      fold_RF <- mclapply(fold_RF, function(x) {
        x$grow(replace = TRUE)
        x
      }, mc.cores = 1)
      
      # 2-4-3-3 Check that all trees were grown correctly
      #         --> none w/o 'child_node_ID' after that!
      fold_RF <- all_trees_grown_correctly(fold_RF)
      
      # 2-4-4 Add the fitted tree to the Forest!
      # 2-4-4-1 Copy the tree first, so we do no overwriting!
      fold_RF_copy <- sapply(fold_RF, function(x) x$copy())
      
      # 2-4-4-2 Remove the fold_RF and add the copied fold_RF_copy to 'Forest'
      rm(fold_RF)
      Forest[[i_]] <- fold_RF_copy
      
      # 2-4-5 Count up i_ - used to fill Forest list!
      i_  = i_ + 1
    }
    
    # 2-5 Evaluate the RF on the different Testsets! Fromfull TestSet w/o any 
    #     missing blocks to TestSets w/ only one observed Block!
    # 2-5-1 Full TestSet!
    print("Evaluation full TestSet -------------------------------------------")
    curr_Forest <- copy_forrest(Forest)
    full[[i]]   <- do_evaluation(Forest = curr_Forest, testdata = test, 
                                 weighted = weighted, weight_metric = weight_metric)
    
    # 2-5-2 TestSet with 1 missing block!
    print("Evaluation TestSet w/ 1 missing omics block------------------------")
    curr_Forest  <- copy_forrest(Forest)
    miss1_A[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                  weight_metric = weight_metric,
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$A)])
    curr_Forest  <- copy_forrest(Forest)
    miss1_B[[i]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                  weight_metric = weight_metric,
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$B)])
    curr_Forest  <- copy_forrest(Forest)
    miss1_C[[i]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                  weight_metric = weight_metric,
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$C)])
    curr_Forest  <- copy_forrest(Forest)
    miss1_D[[i]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                  weight_metric = weight_metric,
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$D)])
    
    # 2-5-3 TestSet with 2 missing blocks!
    print("Evaluation TestSet w/ 2 missing omics blocks-----------------------")
    curr_Forest   <- copy_forrest(Forest)
    miss2_CD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                 curr_data$block_names$D))])
    curr_Forest   <- copy_forrest(Forest)
    miss2_BD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                 curr_data$block_names$D))])
    curr_Forest   <- copy_forrest(Forest)
    miss2_BC[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted, 
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                 curr_data$block_names$c))])
    curr_Forest   <- copy_forrest(Forest)
    miss2_AD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$D))])
    curr_Forest   <- copy_forrest(Forest)
    miss2_AC[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                 curr_data$block_names$A))])
    
    curr_Forest   <- copy_forrest(Forest)
    miss2_AB[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$B))])
    # 2-5-4 Testset with 3 missing blocks!
    print("Evaluation TestSet w/ 3 missing omics blocks-----------------------")
    curr_Forest    <- copy_forrest(Forest)
    miss3_ABC[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                    weight_metric = weight_metric,
                                    testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                  curr_data$block_names$B,
                                                                                  curr_data$block_names$C))])
    curr_Forest    <- copy_forrest(Forest)
    miss3_ACD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                    weight_metric = weight_metric,
                                    testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                  curr_data$block_names$C,
                                                                                  curr_data$block_names$D))])
    curr_Forest    <- copy_forrest(Forest)
    miss3_ABD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted, 
                                    weight_metric = weight_metric,
                                    testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                  curr_data$block_names$A,
                                                                                  curr_data$block_names$D))])
    
    curr_Forest    <- copy_forrest(Forest)
    miss3_BCD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                    weight_metric = weight_metric,
                                    testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                  curr_data$block_names$C,
                                                                                  curr_data$block_names$D))])
    # 2-5-5 Evaluation on single Block Testdata
    print("Evaluation TestSet w/ only 1 observed Block -----------------------")
    curr_Forest   <- copy_forrest(Forest)
    single_A[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                 curr_data$block_names$C,
                                                                                 curr_data$block_names$D,
                                                                                 curr_data$block_names$clin_block))])
    curr_Forest   <- copy_forrest(Forest)
    single_B[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$C,
                                                                                 curr_data$block_names$D,
                                                                                 curr_data$block_names$clin_block))])
    
    curr_Forest   <- copy_forrest(Forest)
    single_C[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$B,
                                                                                 curr_data$block_names$D,
                                                                                 curr_data$block_names$clin_block))])
    curr_Forest   <- copy_forrest(Forest)
    single_D[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$B,
                                                                                 curr_data$block_names$C,
                                                                                 curr_data$block_names$clin_block))])
    curr_Forest     <- copy_forrest(Forest)
    single_CL[[i]]  <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                     weight_metric = weight_metric,
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
  
  # 4-2 Collect the Settings, used to do the CV!
  settings <- list("data_path"     = path,
                   "num_folds"     = k_splits,
                   "response"      = response,
                   "weighted"      = weighted,
                   "weight_metric" = weight_metric,
                   "num_trees"     = num_trees,
                   "mtry"          = mtry, 
                   "min_node_size" = min_node_size,
                   "unorderd_factors" = unorderd_factors,
                   "time_for_CV"      = time_for_CV)
  
  # 4-3 Return both lists!
  return(list("res_all"  = res_all, 
              "settings" = settings))
}

do_CV_4 <- function(path = "data/external/Dr_Hornung/subsetted_12345/missingness_1234/BLCA_3.RData",
                    weighted = TRUE, weight_metric = "Acc", 
                    num_trees = 10, mtry = NULL, min_node_size = NULL,
                    unorderd_factors = "ignore") {
  "CrossValidate the Approach when the Traindata has blockwise missingness
   according to scenario 4 --> 'path' must end in '4.RData', else error!
   
   path must lead to a list with 2 entrances: 'data'  &  'block_names'
   - 'data' is a list filled with 'k' test-train-splits
      --> k-fold-Validation on this test-train-splits!
   - 'block_names' is a list filled with the names of the single blocks 
      & must be ['A', 'B', 'clin_block'] for this scenario!
      
   Based on the 'k' test-train-splits in 'data', we will fit foldwise RFs to the
   train data (that has blockwise missingness in it). 
   Then we ensemble the predicitons from foldwise fitted RFs to single  
   predicitons & rate these w/ Accuracy, Precision, Specifity, F1-Socre,...
   
   The TestingSituations are different, as we can test the models on fully 
   observed testdata, on testdata w/ 1 missing block, etc...
   
   Args:
      - path (char)         : Path to the data that is already split to test train
                              Must end in '3.RData'
      - weighted (bool)     : When ensembling the prediciton from the fodlwise
                              RFs shall we weight the predictions by the OOB 
                              performance of the foldwise fitted RFs? The higher
                              the OOB-Metric for a foldwise fitted RF, the higher
                              its contribution to the prediction!
      - weight_metric (chr) : Which metric to use to assigning weights to the 
                              different predictions? 
                                - must be 'Acc' or 'F1' 
                                - if weighted = FALSE, it will be ignored!
      - num_trees (int)     : Amount of trees, we shall grow on each foldwise
                              fitted random Forest!
      - mtry (int)          : Amount of split-variables we try, when looking for 
                              a split variable! 
                                - If 'NULL': mtry = sqrt(p)
      - min_node_size (int) : Amount of Observations a node must at least 
                              contain, so the model keeps on trying to split 
                              them!  
                                - If 'NULL: It is automatically set to 10 in the
                                            'simpleRF()' function!
      - unorderd_factors (chr) : How to handle non numeric features!
                                ['ignore', 'order_once', 'order_split', 'partition']
    Return:
      - list filled w/:
        * 'res_all' [the CV Results on different testsets]:
           >> A = CNV-Block, B = RNA-Block, C = Mutation-Block, D = Mirna-Block
            - full    : CV Results for each fold on the fully observed testdata!
            - miss1_A : CV Results for each fold on the testdata, w/ missing A block!
            - miss1_B : CV Results for each fold on the testdata, w/ missing B block!
                .
                . --> done for all possible permutations
                .
            - miss2_AC: CV Results for each fold on the testdata, w/ missing A & C block!
                .
                .
            - single_D:  CV-Results for each fold on the testdata w/ D block only!
                         
        * 'settings' [settings used to do the CV - all arguments!]
            - datapath, seed, response, mtry, time for CV,.... 
  "
  # [0] Check Inputs  ----------------------------------------------------------
  # 0-0 mtry, min_node_size & num_trees are all checked within 'simpleRF()'
  # 0-1 path must be numeric and have '1.RData' in it!
  assert_string(path, fixed = "4.RData")
  
  # 0-2 weighted must be a single boolean
  assert_logical(weighted, len = 1)
  
  # 0-3 unorderd factors must be a legit value
  if (!(unorderd_factors %in% c("ignore", "order_once", "order_split", "partition"))) {
    stop("Unknown value for argument 'unordered_factors'")
  }
  
  # 0-4 Check weight_metric to be meaningful [if weighted is true at all]!
  if (weighted) {
    
    # 0-4-1 Check that it is character
    assert_character(weight_metric, len = 1)
    
    # 0-4-2 Check it has a valid value!
    if (!(weight_metric %in% c("Acc", "F1"))) {
      stop("'weight_metric' must be 'Acc' or 'F1'!")
    }
  }
  
  # [1] Prepare CV  ------------------------------------------------------------
  # 1-1 Load CV-Data [already splitted - data checked in 'load_CV_data' itself]
  curr_data <- load_CV_data(path = path)
  
  # 1-1-1 Must contain 'A', 'B', 'C' & 'D' as block_names
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
  
  # [2] Start the CV and loop over all test-train splits in data  --------------
  for (i in seq_len(k_splits)) {
    
    # 2-1 Current Fold Status
    print(paste0("FOLD ", i, "/", k_splits, " -------------------------------"))
    
    # 2-2 Extract the test and train set from 'curr_data'
    train <- curr_data$data[[i]]$train
    test  <- curr_data$data[[i]]$test
    
    # 2-3 Get the Observations that belong to the same fold [same feature space]
    # 2-3-1 Get for each obs. the index of the observed feas
    observed_feas <- sapply(seq_len(nrow(train)), function(x) {
      paste0(which(!(is.na(train[x,]))), collapse = "_")
    })
    
    # 2-3-2 Keep the unique observed feas [equals the different folds]
    #       That we use to assign obs. to the differnt folds!
    observed_folds <- unique(observed_feas)
    
    print(paste0("Found ", length(observed_folds), " unique folds!"))
    
    # 2-4 Train foldwise RFs for each fold seperatly! For this loop over all folds
    #     [observed_folds contains all unique features observed for diff folds]!
    #     --> Results in a Forest of length 'lenght(observed_folds)' & each 
    #         entrance consits of 'num_trees' foldwise fitted trees!
    Forest <- list(); i_  <- 1; start_time <- Sys.time()
    
    for (fold_ in observed_folds) {
      
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
                          min_node_size     = min_node_size,
                          replace           = TRUE,  # always TRUE, as we need OOB!
                          splitrule         = NULL,  # always NULL!
                          unordered_factors = unorderd_factors)
      
      print(paste0("Fit FoldWise RF on current fold: ", i_))
      
      fold_RF <- mclapply(fold_RF, function(x) {
        x$grow(replace = TRUE)
        x
      }, mc.cores = 1)
      
      # 2-4-3-3 Check that all trees were grown correctly
      #         --> none w/o 'child_node_ID' after that!
      fold_RF <- all_trees_grown_correctly(fold_RF)
      
      # 2-4-4 Add the fitted tree to the Forest!
      # 2-4-4-1 Copy the tree first, so we do no overwriting!
      fold_RF_copy <- sapply(fold_RF, function(x) x$copy())
      
      # 2-4-4-2 Remove the fold_RF and add the copied fold_RF_copy to 'Forest'
      rm(fold_RF)
      Forest[[i_]] <- fold_RF_copy
      
      # 2-4-5 Count up i_ - used to fill Forest list!
      i_  = i_ + 1
    }
    
    # 2-5 Evaluate the RF on the different Testsets! Fromfull TestSet w/o any 
    #     missing blocks to TestSets w/ only one observed Block!
    # 2-5-1 Full TestSet!
    print("Evaluation full TestSet -------------------------------------------")
    curr_Forest <- copy_forrest(Forest)
    full[[i]]   <- do_evaluation(Forest = curr_Forest, testdata = test, 
                                 weighted = weighted, weight_metric = weight_metric)
    
    # 2-5-2 TestSet with 1 missing block!
    print("Evaluation TestSet w/ 1 missing omics block------------------------")
    curr_Forest  <- copy_forrest(Forest)
    miss1_A[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                  weight_metric = weight_metric,
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$A)])
    curr_Forest  <- copy_forrest(Forest)
    miss1_B[[i]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                  weight_metric = weight_metric,
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$B)])
    curr_Forest  <- copy_forrest(Forest)
    miss1_C[[i]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                  weight_metric = weight_metric,
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$C)])
    curr_Forest  <- copy_forrest(Forest)
    miss1_D[[i]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                  weight_metric = weight_metric,
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$D)])
    
    # 2-5-3 TestSet with 2 missing blocks!
    print("Evaluation TestSet w/ 2 missing omics blocks-----------------------")
    curr_Forest   <- copy_forrest(Forest)
    miss2_CD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                 curr_data$block_names$D))])
    curr_Forest   <- copy_forrest(Forest)
    miss2_BD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                 curr_data$block_names$D))])
    curr_Forest   <- copy_forrest(Forest)
    miss2_BC[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted, 
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                 curr_data$block_names$c))])
    curr_Forest   <- copy_forrest(Forest)
    miss2_AD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$D))])
    curr_Forest   <- copy_forrest(Forest)
    miss2_AC[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                 curr_data$block_names$A))])
    
    curr_Forest   <- copy_forrest(Forest)
    miss2_AB[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$B))])
    # 2-5-4 Testset with 3 missing blocks!
    print("Evaluation TestSet w/ 3 missing omics blocks-----------------------")
    curr_Forest    <- copy_forrest(Forest)
    miss3_ABC[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                    weight_metric = weight_metric,
                                    testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                  curr_data$block_names$B,
                                                                                  curr_data$block_names$C))])
    curr_Forest    <- copy_forrest(Forest)
    miss3_ACD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                    weight_metric = weight_metric,
                                    testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                  curr_data$block_names$C,
                                                                                  curr_data$block_names$D))])
    curr_Forest    <- copy_forrest(Forest)
    miss3_ABD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted, 
                                    weight_metric = weight_metric,
                                    testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                  curr_data$block_names$A,
                                                                                  curr_data$block_names$D))])
    
    curr_Forest    <- copy_forrest(Forest)
    miss3_BCD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                    weight_metric = weight_metric,
                                    testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                  curr_data$block_names$C,
                                                                                  curr_data$block_names$D))])
    # 2-5-5 Evaluation on single Block Testdata
    print("Evaluation TestSet w/ only 1 observed Block -----------------------")
    curr_Forest   <- copy_forrest(Forest)
    single_A[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                 curr_data$block_names$C,
                                                                                 curr_data$block_names$D,
                                                                                 curr_data$block_names$clin_block))])
    curr_Forest   <- copy_forrest(Forest)
    single_B[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$C,
                                                                                 curr_data$block_names$D,
                                                                                 curr_data$block_names$clin_block))])
    
    curr_Forest   <- copy_forrest(Forest)
    single_C[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$B,
                                                                                 curr_data$block_names$D,
                                                                                 curr_data$block_names$clin_block))])
    curr_Forest   <- copy_forrest(Forest)
    single_D[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$B,
                                                                                 curr_data$block_names$C,
                                                                                 curr_data$block_names$clin_block))])
    curr_Forest     <- copy_forrest(Forest)
    single_CL[[i]]  <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                     weight_metric = weight_metric,
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
  
  # 4-2 Collect the Settings, used to do the CV!
  settings <- list("data_path"     = path,
                   "num_folds"     = k_splits,
                   "response"      = response,
                   "weighted"      = weighted,
                   "weight_metric" = weight_metric,
                   "num_trees"     = num_trees,
                   "mtry"          = mtry, 
                   "min_node_size" = min_node_size,
                   "unorderd_factors" = unorderd_factors,
                   "time_for_CV"      = time_for_CV)
  
  # 4-3 Return both lists!
  return(list("res_all"  = res_all, 
              "settings" = settings))
}

# These Functions are not completly done yet                                 ----
# If all is correct in Situation1 we can complete these functions!!
do_CV_setting4                <- function(data_path = "data/external/Dr_Hornung/Data/ProcessedData_subsets/seed_1234/KIRC_Subset.RData",
                                          response = "gender", seed = 1312, 
                                          weighted = TRUE, weight_metric = NULL,
                                          num_trees = as.integer(10), mtry = NULL, 
                                          min_node_size = NULL, 
                                          unorderd_factors = "ignore") {
  " Function to evaluate RF Adaption on blockwise missing data!
  
    Data is split into test and train set [curently fixed to 5-fold], with the 
    little adjustment, the amount of traindata can be split into 2 folds w/o rest
      --> All train folds have same amount of observations!
      --> TestFold can be a bit smaller!
    
    Then each [equally sized] trainingsfold is censored to scenario 4: That is:
      1. Bind 2 random omics blocks to a single feature block 
         [e.g. 'cnv' & 'mirna' is one block + 'rna' & 'mutation' is one block]
      2. Then divide the testset in two equally sized folds, from which we censor 
         one of the omics feature blocks in each of the folds!

    Then we train a serperate RandomForest on the 2 different training folds 
    [where each fold has different observed omics features] and ensemble the 
    predicitons from these 2 different RFs to a single prediciton and rate these
    w/ Accuracy, Precision, Specifity, F1-Socre,....
    
    The TestingSituations are different, as we can test the models on fully 
    observed testdata, on testdata w/ 1 missing block etc. etc.
 
    Args:
      - data_path (char)    : Path to the data, we want to CV! This should lead 
                              to a file w/ multiple sub DFs 
                              [details see 'create_data()']
      - response (chr)      : The repsonse we want to model - 
                              MUST be in the 'clin'-block & MUST be binary!
      - seed (int)          : Needed to assign obs to fold & to assign which
                              blocks are melted togehter in a reproducible way
      - weighted (bool)     : When ensembling the prediciton from the single
                              RFs shall we weight the predictions by their 
                              'weight_metric' [e.g. Acc, F1]
                              [the higher, the higher the weight]
      - weight_metric (chr) : When assigning weights to the different predictions
                              which metric to use to calc the weight?!
                              [- must be 'Acc' or 'F1' 
                               - if weighted = FALSE, it will be ignored!]
      - num_trees (int)     : amount of trees, we shall grow on each[!] fold
      - mtry (int)          : amount of split-variables we try, when looking for 
                              a split variable! 
                              If 'NULL': mtry = sqrt(p)
      - min_node_size (int) : Amount of Observations a node must at least 
                              contain, so the model keeps on trying to split 
                              them!
                              If 'NULL: It is automatically set to 10 in 
                                        'simpleRF()'
      - unorderd_factors (chr) : How to handle non numeric features!
                                 --> must be in ['ignore', 'order_once', 
                                                 'order_split', 'partition']

    Return:
      - list filled w/:
        * 'res_all' [the CV Results on different testsets]:
            - full    : CV Results for each fold on the fully observed testdata!
            - miss1_A : CV Results for each fold on the testdata, w/ 1 missing 
                        omics-block [Block A]!
            - miss1_B : CV Results for each fold on the testdata, w/ 1 missing 
                        omics-block [Block B]!
            - single_CL: CV-Results for each fold on the testdata w/ only 
                         clinical features
            - single_A:  CV-Results for each fold on the testdata w/ only 
                         block A features
            - single_B:  CV-Results for each fold on the testdata w/ only 
                         block B features
                         
        * 'settings' [settings used to do the CV - all arguments!]
            - datapath, seed, response, mtry, .... 
  "
  # [0] Check Inputs -----------------------------------------------------------
  # 0-0 data_path, response are all checked within 'load_data_extract_block_names()'
  
  # 0-1 mtry, min_node_size & num_trees are all checked within simpleRF()
  
  # 0-2 weighted must be boolean & seed integer
  assert_logical(weighted, len = 1)
  assert_int(seed)
  
  # 0-3 unorderd factors must be a legit value
  if (!(unorderd_factors %in% c("ignore", "order_once", "order_split", "partition"))) {
    stop("Unknown value for argument 'unordered_factors'")
  }
  
  # [1] Get the data & dimensions of train & test folds! -----------------------
  # 1-1 Load the blockwise Omics-Data & create a single DF 
  data <- load_data_extract_block_names(path = data_path, response = response)

  # 1-2 Get amount of Obs. we need for equally sized train folds 
  #     double it, as it was originally for 4 not 2 folds!
  obs_per_fold <- get_obs_per_fold(data = data$data)
  obs_per_fold$amount_train_fold <- obs_per_fold$amount_train_fold * 2 
  
  # [2] Shuffle the IDs and create empty lists to store results ----------------
  # 2-1 Shuffle IDs from 'data' randomly, for splitting it to test & train!
  set.seed(seed)
  fold_ids <- sample(nrow(data$data), nrow(data$data), replace = FALSE)
  
  # 2-2 Create empty lists to store results in!
  # 2-2-1 Full TestSet
  full <- list()
  
  # 2-2-2 TestSet with 1 missing omics-block
  miss1_A <- list(); miss1_B <- list()
  
  # 2-2-5 Single BlockTestSet
  single_A <- list(); single_B <- list(); single_CL <- list()
  
  # [3] Assign which of the omics blocks should belong together! ---------------
  #     Randomly sample the letters "cnv", "rna", "mutation" & "mirna" & 
  #     bind the first 2 together and the last 2!
  set.seed(seed)
  blocks_together <- sample(c("cnv_block", "rna_block", "mutation_block", "mirna_block"),
                            size = 4, replace = F)
  
  # 3-1 Bind the colnams of the first 2 sampled blocks!
  block_A_names <- c(data$block_names[[which(names(data$block_names) == blocks_together[1])]],
                     data$block_names[[which(names(data$block_names) == blocks_together[2])]])
  
  # 3-2 Bind the colnames if the last 2 sampled blocks!
  block_B_names <- c(data$block_names[[which(names(data$block_names) == blocks_together[3])]],
                     data$block_names[[which(names(data$block_names) == blocks_together[4])]])

  # [4] Start the CV, split data to Test and Train and evaluate it! ------------
  for (i in 0:4) {
    
    print(paste0("FOLD: ", as.character(i + 1), "/5 -------------------------"))
    
    # [1] Get TestSet from 'data', by taking the first 'obs_per_fold$amount_test' 
    #     IDs in 'fold_ids'
    test_ids <- fold_ids[((i * obs_per_fold$amount_test) + 1):(((i + 1) * obs_per_fold$amount_test))]
    test_df  <- data$data[test_ids,]
    
    # [2] Get the TrainSet from 'data' [= IDs not in TestSet] & 
    #     induce blockwise missingness!
    train_ids <- fold_ids[-which(fold_ids %in% test_ids)]
    train_df  <- data$data[train_ids,]
    
    # [3] Induce blockwise missingness [SCENARIO_4] 
    # 3-1 Sample equally sized 'observed' blocks [according to SCENARIO_2]
    set.seed(seed + i)
    observed_blocks <- sample(c(rep("fold1", obs_per_fold$amount_train_fold), 
                                rep("fold2", obs_per_fold$amount_train_fold)),
                              obs_per_fold$amount_train, replace = FALSE)
    
    # 3-2 Split Traindata into observed blocks! The resulting blocks will only
    #     contain the features, in the blocks!
    # 3-2-1 Fold1
    block1 <- train_df[which(observed_blocks == "fold1"), 
                       c(response, data$block_names$clin_block, block_A_names)]
    
    # 3-2-2 Fold2
    block2 <- train_df[which(observed_blocks == "fold2"), 
                       c(response, data$block_names$clin_block, block_B_names)]
    
    # [4] Fit 'num_trees' decision trees on each block!
    # 4-1 Get the Formula we use to fit all DecisionTrees/ partial forrests!
    formula_all <- as.formula(paste(response, " ~ ."))
    
    # 4-2 BLOCK1 - grow the trees [as long, as all of them are grown correctly]
    trees1 <- simpleRF(formula = formula_all, data = block1, 
                       num_trees = num_trees, mtry = mtry, 
                       min_node_size = min_node_size, replace = TRUE,  
                       splitrule = NULL, unordered_factors = unorderd_factors)
    trees1 <- mclapply(trees1, function(x) {
      x$grow(replace = TRUE)
      x
    }, mc.cores = 1)
    
    # 4-3 BLOCK2 - grow the trees [as long, as all of them are grown correctly]
    trees2 <- simpleRF(formula = formula_all, data = block2, 
                       num_trees = num_trees, mtry = mtry, 
                       min_node_size = min_node_size, replace = TRUE,  
                       splitrule = NULL, unordered_factors = unorderd_factors)
    trees2 <- mclapply(trees2, function(x) {
      x$grow(replace = TRUE)
      x
    }, mc.cores = 1)
    
    # [5] Check, that all of the trees were grown correctly & create a forest from it!
    trees1 <- all_trees_grown_correctly(trees1)
    trees2 <- all_trees_grown_correctly(trees2)
    
    Forest <- list(trees1, trees2)
    
    # [6] Start Testing!
    # 6-1 FULL TESTSET - all blocks observed!
    #       --> copy the forrest, so we don't override original tree!
    print("Evaluation full TestSet -------------------------------------------")
    curr_Forest   <- copy_forrest(Forest)
    full[[i + 1]] <- do_evaluation(Forest = curr_Forest, testdata = test_df, 
                                   weighted = weighted, weight_metric = weight_metric)
    rm(curr_Forest); gc()
    
    # 6-2 TESTSET ONE OMICS BLOCK MISSING - one block is missing in TestData, 
    #     everytime before evaluation we need to copy the Forest, as the 
    #     evaluation can leed to pruned trees [also outside of the function] 
    print("Evaluation TestSet w/ 1 missing omics block------------------------")
    curr_Forest      <- copy_forrest(Forest)
    miss1_A[[i + 1]] <- do_evaluation(Forest   = curr_Forest, weighted = weighted,
                                      weight_metric = weight_metric,
                                      testdata = test_df[,-which(colnames(test_df) %in% block_A_names)])
    rm(curr_Forest); gc()
    
    curr_Forest      <- copy_forrest(Forest)
    miss1_B[[i + 1]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                      weight_metric = weight_metric,
                                      testdata = test_df[,-which(colnames(test_df) %in% block_B_names)])
    rm(curr_Forest); gc()
    
    # 6-3 TESTSET with one single block only!
    #     TestData, everytime before evaluation we need to copy the Forest, 
    #     as the  evaluation can leed to pruned trees [also outside of the function] 
    print("Evaluation TestSet w/ only 1 observed Block -----------------------")

    curr_Forest        <- copy_forrest(Forest)
    single_A[[i + 1]]  <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                        weight_metric = weight_metric,
                                        testdata = test_df[,-which(colnames(test_df) %in% c(data$block_names$clin_block,
                                                                                            block_B_names))])
    rm(curr_Forest); gc()
    
    curr_Forest        <- copy_forrest(Forest)
    single_B[[i + 1]]  <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                        weight_metric = weight_metric,
                                        testdata = test_df[,-which(colnames(test_df) %in% c(data$block_names$clin_block,
                                                                                            block_A_names))])
    rm(curr_Forest); gc()
    
    curr_Forest         <- copy_forrest(Forest)
    single_CL[[i + 1]]  <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                         weight_metric = weight_metric,
                                         testdata = test_df[,-which(colnames(test_df) %in% c(block_A_names,
                                                                                             block_B_names))])
    rm(curr_Forest); gc()
  }
  
  # [5] Return the metric & settings of the fitting! ---------------------------
  # 4-1 Collect all CV Results in a list!
  res_all <- list("full" = full,
                  "miss1_A" = miss1_A, "miss1_B" = miss1_B,
                  "single_A" = single_A, "single_B" = single_B,
                  "single_CL" = single_CL)
  
  # 4-2 Collect the Settings, used to do the CV!
  settings <- list("data_path"     = data_path,
                   "response"      = response, 
                   "seed"          = seed,
                   "weighted"      = weighted,
                   "weight_metric" = weight_metric,
                   "num_trees"     = num_trees,
                   "mtry"          = mtry, 
                   "min_node_size" = min_node_size,
                   "unorderd_factors" = unorderd_factors,
                   "blocks_together"  = blocks_together) # which blocks belonged together?!
  
  # 4-3 Return both lists!
  return(list("res_all"  = res_all, 
              "settings" = settings))
}