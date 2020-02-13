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
  "Load the subsetted, test-train splitted data!
  
  Args:
    path (str) : Path to the data we want for CV!
                 Path must point to a list, that consits of two more lists 
                 'data' & 'block_names'! 
                 These two lists are further checked in this function!
                 - $data (list) must contain n entrances of 'train' & 'test'
                                where each is a dataframe!
                 - $block_names (list) must contain all block names!
                 
  
  Return:
    A list, filled with 'data' & 'block_names' that can be used for CV!
  
  "
  # [0] Check Inputs  ----------------------------------------------------------
  # 0-1 'path' of type string 
  assert_string(path)
  
  # [1] Load the Data and check it  --------------------------------------------
  # 1-1 Load DF and and put data to 'data_'-variable
  data_ <- load(path)
  data_ <- eval(as.symbol(data_))
  
  # 1-2 Check that data_ has what we need!
  # 1-2-1 List with 2 entrances
  assert_list(data_, len = 2)
  
  # 1-2-2 'data' and 'block_names'
  if (!('data' %in% names(data_))) stop("'path' should lead to list w/ 'data' entry")
  if (!('block_names' %in% names(data_))) stop("'path' should lead to list w/ 'block_names' entry")
  
  # 1-2-3 'data' & 'block_names' should be lists aswell! 
  assert_list(data_$data)
  assert_list(data_$block_names)
  
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
  
  # 1-2-6 'block_names' must contain the correct block names!
  needed_names <- c("clin_block", "cnv_block", "rna_block", "mutation_block", "mirna_block")
  
  res <- sapply(needed_names, function(x) x %in% names(data_$block_names))
  
  if (any(!res)) stop("'block_names' do not have the correct names!")
  
  
  # [2] Return the data  -------------------------------------------------------
  return(data_)
}
all_trees_grown_correctly <- function(trees) {
  "Check, whether 'trees', were grown correctly & if not grow these trees again,
   as long, as they are grown correctly!
   --> Growning not correctly: >> No 'childNodeIDs' 
  
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
    if (trees[[x]]$child_nodeIDs[1] == "pruned") {
      return(NULL)
    } else {
      return(x)
    }
  })
  
  # 1-1-1 If all the trees were pruned in the first node, we can not do
  #       any OOB predictions --> OOB-Accuracy = 0
  if (length(usable_trees) < 1) return(0)
  
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
  " Get the aggregated predicition from all trees! 
    Evaluate the aggregated predicitons & return metrics!
  
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
    assert_character(weight_metric)
    
    # check that it is of length 1
    if (length(weight_metric) != 1) {
      stop("'weight_metric' has more than 1 element!")
    }
    
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
  #     Get a prediction for every observation in TestData from all the trees
  #     Obacht: the prediction is based on the amount of usable trees!          
  #             --> this might reduce the forrest even to a single tree!        
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
  if (weighted) {
    tree_weights <- c()
    
    # Really slow, think about parallel!
    for (l in 1:length(Forest)) {
      tree_weights[l] <- get_oob_weight_metric(trees = Forest[[l]], 
                                               weight_metric = weight_metric)
    }
    
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
  #     Select 'direction = "auto"' and it will choose the class to be positive
  #     automatically!
  roc <- pROC::auc(pROC::roc(testdata[,1], all_forrest_preds_probs_class_0, 
                             levels = levels(Forest[[1]][[1]]$data$data[,1]), 
                             direction = "auto"))
  
  # 5-3 MCC Matthews correlation coefficient [only for binary cases!]
  mcc <- mcc_metric(conf_matrix = confmat)
  
  # [6] Create a list to collect the results!
  res <- list("Accuracy"    = confmat$overall["Accuracy"],
              "Sensitifity" = confmat$byClass["Sensitivity"],
              "Specificity" = confmat$byClass["Specificity"],
              "Precision"   = confmat$byClass["Precision"],
              "Recall"      = confmat$byClass["Recall"],
              "F1"          = confmat$byClass["F1"],
              "Balance_Acc" = confmat$byClass["Balanced Accuracy"],
              "AUC"         = as.numeric(roc),
              "MCC"         = mcc)
  
  # 6-1 If the F1-Score/ Precision/ Recall is NA, then we set it to 0 [-1 for MCC]! 
  if (is.na(res$F1))        res$F1        <- 0
  if (is.na(res$Precision)) res$Precision <- 0
  if (is.na(res$Recall))    res$Recall    <- 0
  if (is.na(res$MCC))       res$MCC       <- -1
  
  return(as.vector(res))
}

path = "data/external/Dr_Hornung/subsetted_12345/missingness_1234/BLCA_1.RData"
weighted = TRUE
weight_metric = "Acc"
num_trees = 10
mtry = NULL
min_node_size = NULL
unorderd_factors = "ignore"

do_CV <- function(path = "data/external/Dr_Hornung/subsetted_12345/missingness_1234/BLCA_1.RData",
                  weighted = TRUE, weight_metric = "Acc", 
                  num_trees = 10, mtry = NULL, min_node_size = NULL,
                  unorderd_factors = "ignore") {
  "Function to do CV"
  # [0] Check Inputs  ----------------------------------------------------------
  
  # [1] Prepare CV  ------------------------------------------------------------
  # 1-1 Load CV-Data [already splitted - data checked in 'load_CV_data' itself]
  curr_data <- load_CV_data(path = path)
  
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
  i = 1
  for (i in seq_len(k_splits)) {
    
    # 2-1 Current Fold Status
    print(paste0("FOLD ", i, "/", k_splits, " -------------------------------"))
    
    # 2-2 Extract the test and train set from 'curr_data'
    train <- curr_data$data[[i]]$train
    test  <- curr_data$data[[i]]$test
    
    # 2-3 Get the Observations that belong to the same fold [same feature space]
    # 2-3-1 Get for each obs. the index of the observed feas
    obs_feas <- sapply(seq_len(nrow(train)), function(x) {
      paste0(which(!(is.na(train[x,]))), collapse = "_")
    })
    
    # 2-3-2 Keep the unique observed feas [equals the different folds]
    observed_folds <- unique(obs_feas)
    
    # 2-4 Train foldwise RFs for each fold seperatly!
    #     --> Forest is of length lenght(observed_folds), and each entrance & 
    #         each entrance consits of 'num_trees' trees!
    Forest <- list(); i_  = 1
    
    start_time = Sys.time()
    for (fold_ in observed_folds) {
      
      # 2-4-1 Get all Obs. with the feture space as in 'fold_'
      fold_obs_ <- which(obs_feas == fold_)
      
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
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$cnv_block)])
    curr_Forest  <- copy_forrest(Forest)
    miss1_B[[i]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                  weight_metric = weight_metric,
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$rna_block)])
    curr_Forest  <- copy_forrest(Forest)
    miss1_C[[i]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                  weight_metric = weight_metric,
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$mutation_block)])
    curr_Forest  <- copy_forrest(Forest)
    miss1_D[[i]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                  weight_metric = weight_metric,
                                  testdata = test[,-which(colnames(test) %in% curr_data$block_names$mirna_block)])
    
    # 2-5-3 TestSet with 2 missing blocks!
    print("Evaluation TestSet w/ 2 missing omics blocks-----------------------")
    curr_Forest   <- copy_forrest(Forest)
    miss2_CD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$mirna_block,
                                                                                 curr_data$block_names$mutation_block))])
    curr_Forest   <- copy_forrest(Forest)
    miss2_BD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$mirna_block,
                                                                                 curr_data$block_names$rna_block))])
    curr_Forest   <- copy_forrest(Forest)
    miss2_BC[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted, 
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$rna_block,
                                                                                 curr_data$block_names$mutation_block))])
    curr_Forest   <- copy_forrest(Forest)
    miss2_AD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$mirna_block,
                                                                                 curr_data$block_names$cnv_block))])
    curr_Forest   <- copy_forrest(Forest)
    miss2_AC[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$mutation_block,
                                                                                 curr_data$block_names$cnv_block))])
    
    curr_Forest   <- copy_forrest(Forest)
    miss2_AB[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$rna_block,
                                                                                 curr_data$block_names$cnv_block))])
    # 2-5-4 Testset with 3 missing blocks!
    print("Evaluation TestSet w/ 3 missing omics blocks-----------------------")
    curr_Forest    <- copy_forrest(Forest)
    miss3_ABC[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                    weight_metric = weight_metric,
                                    testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$rna_block,
                                                                                  curr_data$block_names$cnv_block,
                                                                                  curr_data$block_names$mutation_block))])
    curr_Forest    <- copy_forrest(Forest)
    miss3_ACD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                    weight_metric = weight_metric,
                                    testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$mirna_block,
                                                                                  curr_data$block_names$cnv_block,
                                                                                  curr_data$block_names$mutation_block))])
    curr_Forest    <- copy_forrest(Forest)
    miss3_ABD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted, 
                                    weight_metric = weight_metric,
                                    testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$mirna_block,
                                                                                  curr_data$block_names$cnv_block,
                                                                                  curr_data$block_names$rna_block))])
    
    curr_Forest    <- copy_forrest(Forest)
    miss3_BCD[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                    weight_metric = weight_metric,
                                    testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$mirna_block,
                                                                                  curr_data$block_names$mutation_block,
                                                                                  curr_data$block_names$rna_block))])
    # 2-5-5 Evaluation on single Block Testdata
    print("Evaluation TestSet w/ only 1 observed Block -----------------------")
    curr_Forest   <- copy_forrest(Forest)
    single_A[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$mirna_block,
                                                                                 curr_data$block_names$mutation_block,
                                                                                 curr_data$block_names$rna_block,
                                                                                 curr_data$block_names$clin_block))])
    curr_Forest   <- copy_forrest(Forest)
    single_B[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$mirna_block,
                                                                                 curr_data$block_names$mutation_block,
                                                                                 curr_data$block_names$cnv_block,
                                                                                 curr_data$block_names$clin_block))])
    
    curr_Forest   <- copy_forrest(Forest)
    single_C[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$mirna_block,
                                                                                 curr_data$block_names$cnv_block,
                                                                                 curr_data$block_names$rna_block,
                                                                                 curr_data$block_names$clin_block))])
    curr_Forest   <- copy_forrest(Forest)
    single_D[[i]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                   weight_metric = weight_metric,
                                   testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$cnv_block,
                                                                                 curr_data$block_names$mutation_block,
                                                                                 curr_data$block_names$rna_block,
                                                                                 curr_data$block_names$clin_block))])
    curr_Forest     <- copy_forrest(Forest)
    single_CL[[i]]  <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                     weight_metric = weight_metric,
                                     testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$mirna_block,
                                                                                   curr_data$block_names$mutation_block,
                                                                                   curr_data$block_names$rna_block,
                                                                                   curr_data$block_names$cnv_block))])
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

# Run a example and check the results!                                       ----
start_time <- Sys.time()
no_weight <- do_CV_setting1(num_trees = as.integer(250), weighted = FALSE,
                            data_path = "data/external/Dr_Hornung/Data/ProcessedData_subsets/seed_1238/LGG_subset.RData")
no_weight_time_all <- Sys.time() - start_time 

start_time <- Sys.time()
acc_weight <- do_CV_setting1(num_trees = as.integer(250),
                             weight_metric = "Acc", weighted = TRUE,
                             data_path = "data/external/Dr_Hornung/Data/ProcessedData_subsets/seed_1238/LGG_subset.RData")
acc_weight_time_all <- Sys.time() - start_time 

start_time <- Sys.time()
f1_weight <- do_CV_setting1(num_trees = as.integer(250),
                            weight_metric = "F1", weighted = TRUE,
                            data_path = "data/external/Dr_Hornung/Data/ProcessedData_subsets/seed_1238/LGG_subset.RData")
f1_weight_time_all <- Sys.time() - start_time 


all_res_roman <- list(no_weight, acc_weight, f1_weight)
save(all_res_roman, file = "Roman_seed1238_subset_LGG_diff_weight_approaches.RData")

# These Functions are not completly done yet                                 ----
# If all is correct in Situation1 we can complete these functions!!
do_CV_setting2                <- function(data_path = "data/external/Dr_Hornung/Data/ProcessedData_subsets/seed_1234/KIRC_Subset.RData",
                                          response = "gender", seed = 1312, 
                                          weighted = TRUE, weight_metric = NULL,
                                          num_trees = as.integer(10), mtry = NULL, 
                                          min_node_size = NULL, 
                                          unorderd_factors = "ignore") {
  " Function to evaluate RF Adaption on blockwise missing data!
  
    Data is split into test and train set [curently fixed to 5-fold], with the 
    little adjustment, the amount of traindata can be split into 4 folds w/o rest
      --> All train folds have same amount of observations!
      --> TestFold can be a bit smaller!
    
    Then each [equally sized] trainingsfold is censored to scenario 2, so that 
    each fold has an observed clinical block + an additional observed omics 
    block [1. fold has Clin + 4 omics blocks (fully observed),
           2. fold has Clin + 3 omics blocks, 
           3. fold has Clin + 2 omics blocks,
           4. fold has Clin + 1 Omics Block!]

    Then we train a serperate RandomForest on the 4 different training folds 
    [where each fold has different observed features] and ensemble the 
    predicitons from these 4 different RFs to a single prediciton and rate these
    w/ Accuracy, Precision, Specifity, F1-Socre,....
    
    The TestingSituations are different, as we can test the models on fully 
    observed testdata, on testdata w/ 1 missing block etc. etc.
 
    Args:
      - data_path (char)    : Path to the data, we want to CV! This should lead 
                              to a file w/ multiple sub DFs 
                              [details see 'create_data()']
      - response (chr)      : The repsonse we want to model - 
                              MUST be in the 'clin'-block & MUST be binary!
      - seed (int)          : Needed for assignin obs to certain folds &
                              to assign which block which feature block in 
                              a reproducable way!
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
                        omics-block [cnv]!
            - miss1_B : CV Results for each fold on the testdata, w/ 1 missing 
                        omics-block [rna]!
            - miss1_C : CV Results for each fold on the testdata, w/ 1 missing 
                        omics-block [mutation]!
            - miss1_D : CV Results for each fold on the testdata, w/ 1 missing 
                        omics-block [mirna]!
            - miss2_CD: CV Results for each fold on the testdata, w/ 2 missing 
                        omics-block [mutation & mirna]!
            - miss2_BD: CV Results for each fold on the testdata, w/ 2 missing 
                        omics-block [rna & mirna]!
            - miss2_BC: CV Results for each fold on the testdata, w/ 2 missing 
                        omics-block [rna & mutation]!
            - miss2_AD: CV Results for each fold on the testdata, w/ 2 missing 
                        omics-block [cnv & mirna]!
            - miss2_AC: CV Results for each fold on the testdata, w/ 2 missing 
                        omics-block [cnv & mutation]!
            - miss2_AB: CV Results for each fold on the testdata, w/ 2 missing 
                        omics-block [cnv & rna]!
            - miss3_ABC: CV Results for each fold on the testdata, w/ 3 missing 
                         omics-block [cnv & rna & mutation]!
            - miss3_ACD: CV Results for each fold on the testdata, w/ 3 missing 
                         omics-block [cnv & mutation & mirna]!
            - miss3_ABD: CV Results for each fold on the testdata, w/ 3 missing 
                         omics-block [cnv & rna & mirna]!
            - miss3_BCD: CV Results for each fold on the testdata w/ 3 missing 
                         omics -bloc [rna & mutation & mirna]
            - single_CL: CV-Results for each fold on the testdata w/ only 
                         clinical features
            - single_A:  CV-Results for each fold on the testdata w/ only 
                         block A features [CNV]
            - single_B:  CV-Results for each fold on the testdata w/ only 
                         block B features [RNA]
            - single_C:  CV-Results for each fold on the testdata w/ only 
                         block C features [Mutation]
            - single_D:  CV-Results for each fold on the testdata w/ only 
                         block D features [Mirna]
                         
        * 'settings' [settings used to do the CV - all arguments!]
            - datapath, seed, response, mtry, .... 
  "
  # [0] Check Inputs -----------------------------------------------------------
  # 0-0 data_path, response are all checked within 'create_data()'
  
  # 0-1 mtry, min_node_size & num_trees are all checked within load_data_extract_block_names()
  
  # 0-2 weighted must be boolean & seed an integer
  assert_logical(weighted, len = 1)
  assert_int(seed)
  
  # 0-3 unorderd factors must be a legit value
  if (!(unorderd_factors %in% c("ignore", "order_once", "order_split", "partition"))) {
    stop("Unknown value for argument 'unordered_factors'")
  }
  
  # 0-4 Check weight_metric to be meaningful [if weighted is true at all]!
  if (weighted) {
    # check that it is character
    assert_character(weight_metric)
    
    # check that it is of length 1
    if (length(weight_metric) != 1) {
      stop("'weight_metric' has more than 1 element!")
    }
    
    # check it has a valid value!
    if (!(weight_metric %in% c("Acc", "F1"))) {
      stop("'weight_metric' must be 'Acc' or 'F1'!")
    }
  }
  
  # [1] Get the data & dimensions of train & test folds! -----------------------
  # 1-1 Load the blockwise Omics-Data & create a single DF 
  data <- load_data_extract_block_names(path = data_path, response = response)
  
  # 1-2 Get amount of Obs. we need for equally sized train folds 
  obs_per_fold <- get_obs_per_fold(data = data$data)
  
  # [2] Shuffle the IDs and create empty lists to store results ----------------
  # 2-1 Shuffle IDs from 'data' randomly, for splitting it to test & train!
  set.seed(seed)
  fold_ids <- sample(nrow(data$data), nrow(data$data), replace = FALSE)
  
  # 2-2 Create empty lists to store results in!
  # 2-2-1 Full TestSet
  full <- list()
  
  # 2-2-2 TestSet with 1 missing omics-block
  miss1_A <- list(); miss1_B <- list(); miss1_C <- list(); miss1_D <- list()
  
  # 2-2-3 TestSet with 2 missing omics-blocks
  miss2_CD <- list(); miss2_BD <- list(); miss2_BC <- list(); miss2_AD <- list()
  miss2_AC <- list(); miss2_AB <- list()
  
  # 2-2-4 TestSet with 3 missing omics-blocks
  miss3_ABC <- list(); miss3_ABD <- list(); miss3_ACD <- list(); miss3_BCD <- list()
  
  # 2-2-5 Single BlockTestSet
  single_A <- list(); single_B <- list(); single_CL <- list() 
  single_C <- list(); single_D <- list()
  
  # 2-3 Randomly assign the omics blocks to the letters 'A', 'B', 'C', 'D', 
  #     as SCEANRIO2, highly depens on which block is where!
  #     ['A' only observed in 1.fold, whereas 'D' is observed in all folds]
  set.seed(seed)
  letter_feas <- sample(c("data$block_names$cnv_block", "data$block_names$rna_block", 
                          "data$block_names$mutation_block", "data$block_names$mirna_block"), 
                        4, replace = FALSE)
  names(letter_feas) <- c("A", "B", "C", "D")
  
  # [3] Start the CV, split data to Test and Train and evaluate it! ------------
  for (i in 0:4) {
    
    print(paste0("FOLD: ", as.character(i + 1), "/5 -------------------------"))
    
    # [1] Get TestSet from 'data', by taking the first 'obs_per_fold$amount_test' 
    #     IDs in 'fold_ids'
    test_ids <- fold_ids[((i * obs_per_fold$amount_test) + 1):(((i + 1) * obs_per_fold$amount_test))]
    test_df  <- data$data[test_ids,]
    
    # [2] Get the TrainSet from 'data' [= IDs not in TestSet]
    train_ids <- fold_ids[-which(fold_ids %in% test_ids)]
    train_df  <- data$data[train_ids,]
    
    # [3] Induce blockwise missingness [SCENARIO_2] 
    # 3-1 Sample equally sized 'observed' blocks [according to SCENARIO_2]
    set.seed(seed + i)
    observed_blocks <- sample(c(rep("Clin, A, B, C, D", obs_per_fold$amount_train_fold), 
                                rep("Clin, B, C, D", obs_per_fold$amount_train_fold),
                                rep("Clin, C, D", obs_per_fold$amount_train_fold), 
                                rep("Clin, D", obs_per_fold$amount_train_fold)),
                              obs_per_fold$amount_train, replace = FALSE)
    
    # 3-2 Split Traindata into observed blocks! The resulting blocks will only
    #     contain the features, in the blocks!
    # 3-2-1 Clin, A, B, C, D
    block1 <- train_df[which(observed_blocks == "Clin, A, B, C, D"), 
                       c(response, 
                         data$block_names$clin_block, 
                         eval(parse(text = letter_feas["A"])),
                         eval(parse(text = letter_feas["B"])),
                         eval(parse(text = letter_feas["C"])),
                         eval(parse(text = letter_feas["D"])))]
    
    # 3-2-2 Clin, B, C, D
    block2 <- train_df[which(observed_blocks == "Clin, B, C, D"), 
                       c(response, 
                         data$block_names$clin_block,
                         eval(parse(text = letter_feas["B"])),
                         eval(parse(text = letter_feas["C"])),
                         eval(parse(text = letter_feas["D"])))]
    
    # 3-2-3 Clin, C, D
    block3 <- train_df[which(observed_blocks == "Clin, C, D"),
                       c(response, 
                         data$block_names$clin_block,
                         eval(parse(text = letter_feas["C"])),
                         eval(parse(text = letter_feas["D"])))]
    
    # 3-2-4 Clin, D"
    block4 <- train_df[which(observed_blocks == "Clin, D"),
                       c(response, 
                         data$block_names$clin_block,
                         eval(parse(text = letter_feas["D"])))]
    
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
    
    # 4-4 BLOCK3 - grow the trees [as long, as all of them are grown correctly]
    trees3 <- simpleRF(formula = formula_all, data = block3, 
                       num_trees = num_trees, mtry = mtry, 
                       min_node_size = min_node_size, replace = TRUE,  
                       splitrule = NULL, unordered_factors = unorderd_factors)
    trees3 <- mclapply(trees3, function(x) {
      x$grow(replace = TRUE)
      x
    }, mc.cores = 1)
    
    # 4-5 BLOCK4 - grow the trees [as long, as all of them are grown correctly]
    trees4 <- simpleRF(formula = formula_all, data = block4, 
                       num_trees = num_trees, mtry = mtry, 
                       min_node_size = min_node_size, replace = TRUE,  
                       splitrule = NULL, unordered_factors = unorderd_factors)
    trees4 <- mclapply(trees4, function(x) {
      x$grow(replace = TRUE)
      x
    }, mc.cores = 1)
    
    # [5] Check, that all of the trees were grown correctly & create a forest from it!
    trees1 <- all_trees_grown_correctly(trees1)
    trees2 <- all_trees_grown_correctly(trees2)
    trees3 <- all_trees_grown_correctly(trees3)
    trees4 <- all_trees_grown_correctly(trees4)
    
    Forest <- list(trees1, trees2, trees3, trees4)
    
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
    miss1_A[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                      weight_metric = weight_metric,
                                      testdata = test_df[,-which(colnames(test_df) %in% eval(parse(text = letter_feas["A"])))])
    rm(curr_Forest); gc()
    
    curr_Forest      <- copy_forrest(Forest)
    miss1_B[[i + 1]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                      weight_metric = weight_metric,
                                      testdata = test_df[,-which(colnames(test_df) %in% eval(parse(text = letter_feas["B"])))])
    rm(curr_Forest); gc()
    
    curr_Forest      <- copy_forrest(Forest)
    miss1_C[[i + 1]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                      weight_metric = weight_metric,
                                      testdata = test_df[,-which(colnames(test_df) %in% eval(parse(text = letter_feas["C"])))])
    rm(curr_Forest); gc()
    
    curr_Forest      <- copy_forrest(Forest)
    miss1_D[[i + 1]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                      weight_metric = weight_metric,
                                      testdata = test_df[,-which(colnames(test_df) %in% eval(parse(text = letter_feas["D"])))])
    rm(curr_Forest); gc()
    
    # 6-3 TESTSET TWO OMICS BLOCKS MISSING - two ommic blocks are missing in
    #     TestData, everytime before evaluation we need to copy the Forest, 
    #     as the  evaluation can leed to pruned trees [also outside of the function] 
    print("Evaluation TestSet w/ 2 missing omics blocks-----------------------")
    
    curr_Forest       <- copy_forrest(Forest)
    miss2_CD[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                       weight_metric = weight_metric,
                                       testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["C"])),
                                                                                           eval(parse(text = letter_feas["D"]))))])
    rm(curr_Forest); gc()
    
    curr_Forest       <- copy_forrest(Forest)
    miss2_BD[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                       weight_metric = weight_metric,
                                       testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["B"])),
                                                                                           eval(parse(text = letter_feas["D"]))))])
    rm(curr_Forest); gc()
    
    curr_Forest       <- copy_forrest(Forest)
    miss2_BC[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted, 
                                       weight_metric = weight_metric,
                                       testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["B"])),
                                                                                           eval(parse(text = letter_feas["C"]))))])
    rm(curr_Forest); gc()
    
    curr_Forest       <- copy_forrest(Forest)
    miss2_AD[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted, 
                                       weight_metric = weight_metric,
                                       testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["A"])),
                                                                                           eval(parse(text = letter_feas["D"]))))])
    rm(curr_Forest); gc()
    
    curr_Forest       <- copy_forrest(Forest)
    miss2_AC[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                       weight_metric = weight_metric,
                                       testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["A"])),
                                                                                           eval(parse(text = letter_feas["C"]))))])
    rm(curr_Forest); gc()
    
    curr_Forest       <- copy_forrest(Forest)
    miss2_AB[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                       weight_metric = weight_metric,
                                       testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["A"])),
                                                                                           eval(parse(text = letter_feas["B"]))))])
    
    # 6-4 TESTSET THREE OMICS BLOCKS MISSING - three omic blocks are missing in
    #     TestData, everytime before evaluation we need to copy the Forest, 
    #     as the  evaluation can leed to pruned trees [also outside of the function] 
    print("Evaluation TestSet w/ 3 missing omics blocks-----------------------")
    
    curr_Forest        <- copy_forrest(Forest)
    miss3_ABC[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                        weight_metric = weight_metric,
                                        testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["A"])),
                                                                                            eval(parse(text = letter_feas["B"])),
                                                                                            eval(parse(text = letter_feas["C"]))))])
    rm(curr_Forest); gc()
    
    curr_Forest        <- copy_forrest(Forest)
    miss3_ACD[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted, 
                                        weight_metric = weight_metric,
                                        testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["A"])),
                                                                                            eval(parse(text = letter_feas["C"])),
                                                                                            eval(parse(text = letter_feas["D"]))))])
    rm(curr_Forest); gc()
    
    curr_Forest        <- copy_forrest(Forest)
    miss3_ABD[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                        weight_metric = weight_metric,
                                        testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["A"])),
                                                                                            eval(parse(text = letter_feas["B"])),
                                                                                            eval(parse(text = letter_feas["D"]))))])
    rm(curr_Forest); gc()
    
    curr_Forest        <- copy_forrest(Forest)
    miss3_BCD[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                        weight_metric = weight_metric,
                                        testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["B"])),
                                                                                            eval(parse(text = letter_feas["C"])),
                                                                                            eval(parse(text = letter_feas["D"]))))])
    rm(curr_Forest); gc()
    
    
    # 7-5 TESTSET THREE OMICS BLOCKS MISSING - three omic blocks are missing in
    #     TestData, everytime before evaluation we need to copy the Forest, 
    #     as the  evaluation can leed to pruned trees [also outside of the function] 
    print("Evaluation TestSet w/ only 1 observed Block -----------------------")
    
    curr_Forest        <- copy_forrest(Forest)
    single_A[[i + 1]]  <- do_evaluation(Forest = curr_Forest, weighted = weighted, 
                                        weight_metric = weight_metric,
                                        testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["B"])),
                                                                                            eval(parse(text = letter_feas["C"])),
                                                                                            eval(parse(text = letter_feas["D"])),
                                                                                            data$block_names$clin_block))])
    rm(curr_Forest); gc()
    
    curr_Forest        <- copy_forrest(Forest)
    single_B[[i + 1]]  <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                        weight_metric = weight_metric,
                                        testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["A"])),
                                                                                            eval(parse(text = letter_feas["C"])),
                                                                                            eval(parse(text = letter_feas["D"])),
                                                                                            data$block_names$clin_block))])
    rm(curr_Forest); gc()
    
    curr_Forest        <- copy_forrest(Forest)
    single_C[[i + 1]]  <- do_evaluation(Forest = curr_Forest, weighted = weighted, 
                                        weight_metric = weight_metric,
                                        testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["A"])),
                                                                                            eval(parse(text = letter_feas["B"])),
                                                                                            eval(parse(text = letter_feas["D"])),
                                                                                            data$block_names$clin_block))])
    rm(curr_Forest); gc()
    
    curr_Forest        <- copy_forrest(Forest)
    single_D[[i + 1]]  <- do_evaluation(Forest = curr_Forest, weighted = weighted, 
                                        weight_metric = weight_metric,
                                        testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["A"])),
                                                                                            eval(parse(text = letter_feas["B"])),
                                                                                            eval(parse(text = letter_feas["C"])),
                                                                                            data$block_names$clin_block))])
    rm(curr_Forest); gc()
    
    curr_Forest         <- copy_forrest(Forest)
    single_CL[[i + 1]]  <- do_evaluation(Forest = curr_Forest, weighted = weighted, 
                                         weight_metric = weight_metric,
                                         testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["A"])),
                                                                                             eval(parse(text = letter_feas["B"])),
                                                                                             eval(parse(text = letter_feas["C"])),
                                                                                             eval(parse(text = letter_feas["D"]))))])
    rm(curr_Forest); gc()
  }
  
  # [4] Return the metric & settings of the fitting! ---------------------------
  # 4-1 Collect all CV Results in a list!
  res_all <- list("full" = full,
                  "miss1_A" = miss1_A, "miss1_B" = miss1_B,
                  "miss1_C" = miss1_C, "miss1_D" = miss1_D,
                  "miss2_CD" = miss2_CD, "miss2_BD" = miss2_BD,
                  "miss2_BC" = miss2_BC, "miss2_AD" = miss2_AD,
                  "miss2_AC" = miss2_AC, "miss2_AB" = miss2_AB,
                  "miss3_ABC" = miss3_ABC, "miss3_ABD" = miss3_ABD,
                  "miss3_ACD" = miss3_ACD, "miss3_BCD" = miss3_BCD,
                  "single_A" = single_A, "single_B" = single_B,
                  "single_C" = single_C, "single_D" = single_D,
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
                   "block_letters"    = letter_feas) # which letter represented 
  # which blocks?!
  
  # 4-3 Return both lists!
  return(list("res_all"  = res_all, 
              "settings" = settings))
}
do_CV_setting3                <- function(data_path = "data/external/Dr_Hornung/Data/ProcessedData_subsets/seed_1234/KIRC_Subset.RData",
                                          response = "gender", seed = 1312, 
                                          weighted = TRUE, weight_metric = NULL,
                                          num_trees = as.integer(10), mtry = NULL, 
                                          min_node_size = NULL,
                                          unorderd_factors = "ignore") {
  " Function to evaluate RF Adaption on blockwise missing data!
  
    Data is split into test and train set [curently fixed to 5-fold], with the 
    little adjustment, the amount of traindata can be split into 4 folds w/o rest
      --> All train folds have same amount of observations!
      --> TestFold can be a bit smaller!
    
    Then each [equally sized] trainingsfold is censored to scenario 3: That is:
    For each fold [1-4] we randomly sample which blocks are observed for the
    different blocks --> each block is deleted w/ prob. of 1/3 & kept w/ pro. 2/3

    Then we train a serperate RandomForest on the 4 different training folds 
    [where each fold has different observed features] and ensemble the 
    predicitons from these 4 different RFs to a single prediciton and rate these
    w/ Accuracy, Precision, Specifity, F1-Socre,....
    
    The TestingSituations are different, as we can test the models on fully 
    observed testdata, on testdata w/ 1 missing block etc. etc.
 
    Args:
      - data_path (char)    : Path to the data, we want to CV! This should lead 
                              to a file w/ multiple sub DFs 
                              [details see 'create_data()']
      - response (chr)      : The repsonse we want to model - 
                              MUST be in the 'clin'-block & MUST be binary!
      - seed (int)          : Needed to assign the obs. to folds & to assign
                              the observed blocks to the folds in a reproducable
                              way!
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
                        omics-block [cnv]!
            - miss1_B : CV Results for each fold on the testdata, w/ 1 missing 
                        omics-block [rna]!
            - miss1_C : CV Results for each fold on the testdata, w/ 1 missing 
                        omics-block [mutation]!
            - miss1_D : CV Results for each fold on the testdata, w/ 1 missing 
                        omics-block [mirna]!
            - miss2_CD: CV Results for each fold on the testdata, w/ 2 missing 
                        omics-block [mutation & mirna]!
            - miss2_BD: CV Results for each fold on the testdata, w/ 2 missing 
                        omics-block [rna & mirna]!
            - miss2_BC: CV Results for each fold on the testdata, w/ 2 missing 
                        omics-block [rna & mutation]!
            - miss2_AD: CV Results for each fold on the testdata, w/ 2 missing 
                        omics-block [cnv & mirna]!
            - miss2_AC: CV Results for each fold on the testdata, w/ 2 missing 
                        omics-block [cnv & mutation]!
            - miss2_AB: CV Results for each fold on the testdata, w/ 2 missing 
                        omics-block [cnv & rna]!
            - miss3_ABC: CV Results for each fold on the testdata, w/ 3 missing 
                         omics-block [cnv & rna & mutation]!
            - miss3_ACD: CV Results for each fold on the testdata, w/ 3 missing 
                         omics-block [cnv & mutation & mirna]!
            - miss3_ABD: CV Results for each fold on the testdata, w/ 3 missing 
                         omics-block [cnv & rna & mirna]!
            - miss3_BCD: CV Results for each fold on the testdata w/ 3 missing 
                         omics -bloc [rna & mutation & mirna]
            - single_CL: CV-Results for each fold on the testdata w/ only 
                         clinical features
            - single_A:  CV-Results for each fold on the testdata w/ only 
                         block A features [CNV]
            - single_B:  CV-Results for each fold on the testdata w/ only 
                         block B features [RNA]
            - single_C:  CV-Results for each fold on the testdata w/ only 
                         block C features [Mutation]
            - single_D:  CV-Results for each fold on the testdata w/ only 
                         block D features [Mirna]
                         
        * 'settings' [settings used to do the CV - all arguments!]
            - datapath, seed, response, mtry, .... 
  "
  # [0] Check Inputs -----------------------------------------------------------
  # 0-0 data_path, response are all checked within 'load_data_extract_block_names()'
  
  # 0-1 mtry, min_node_size & num_trees are all checked within simpleRF()
  
  # 0-2 weighted must be boolean & seed an int
  assert_logical(weighted, len = 1)
  assert_int(seed)
  
  # 0-3 unorderd factors must be a legit value
  if (!(unorderd_factors %in% c("ignore", "order_once", "order_split", "partition"))) {
    stop("Unknown value for argument 'unordered_factors'")
  }
  
  # 0-4 Check weight_metric to be meaningful [if weighted is true at all]!
  if (weighted) {
    # check that it is character
    assert_character(weight_metric)
    
    # check that it is of length 1
    if (length(weight_metric) != 1) {
      stop("'weight_metric' has more than 1 element!")
    }
    
    # check it has a valid value!
    if (!(weight_metric %in% c("Acc", "F1"))) {
      stop("'weight_metric' must be 'Acc' or 'F1'!")
    }
  }
  
  # [1] Get the data & dimensions of train & test folds! -----------------------
  # 1-1 Load the blockwise Omics-Data & create a single DF 
  data <- load_data_extract_block_names(path = data_path, response = response)
  
  # 1-2 Get amount of Obs. we need for equally sized train folds 
  obs_per_fold <- get_obs_per_fold(data = data$data)
  
  # [2] Shuffle the IDs and create empty lists to store results ----------------
  # 2-1 Shuffle IDs from 'data' randomly, for splitting it to test & train!
  set.seed(seed)
  fold_ids <- sample(nrow(data$data), nrow(data$data), replace = FALSE)
  
  # 2-2 Create empty lists to store results in!
  # 2-2-1 Full TestSet
  full <- list()
  
  # 2-2-2 TestSet with 1 missing omics-block
  miss1_A <- list(); miss1_B <- list(); miss1_C <- list(); miss1_D <- list()
  
  # 2-2-3 TestSet with 2 missing omics-blocks
  miss2_CD <- list(); miss2_BD <- list(); miss2_BC <- list(); miss2_AD <- list()
  miss2_AC <- list(); miss2_AB <- list()
  
  # 2-2-4 TestSet with 3 missing omics-blocks
  miss3_ABC <- list(); miss3_ABD <- list(); miss3_ACD <- list(); miss3_BCD <- list()
  
  # 2-2-5 Single BlockTestSet
  single_A <- list(); single_B <- list(); single_CL <- list() 
  single_C <- list(); single_D <- list()
  
  # [3] Assign which blocks were observed for which folds ----------------------
  #     Randomly assign the 'observed' blocks to the 4 different folds!
  #     The Index of TRUE / FALSE are indicators, whether the blocks are observed
  #     [1] = clinical; [2] = CNV; [3] = RNA; [4] = Mutation; [5] = Mirna
  set.seed(seed)
  fold1_obs <- sample(c(TRUE, FALSE), size = 5, replace = T, prob = c(2/3, 1/3))
  names(fold1_obs) <- c("Clin", "A", "B", "C", "D")
  set.seed(seed + 1)
  fold2_obs <- sample(c(TRUE, FALSE), size = 5, replace = T, prob = c(2/3, 1/3))
  names(fold2_obs) <- c("Clin", "A", "B", "C", "D")
  set.seed(seed + 2)
  fold3_obs <- sample(c(TRUE, FALSE), size = 5, replace = T, prob = c(2/3, 1/3))
  names(fold3_obs) <- c("Clin", "A", "B", "C", "D")
  set.seed(seed + 3)
  fold4_obs <- sample(c(TRUE, FALSE), size = 5, replace = T, prob = c(2/3, 1/3))
  names(fold4_obs) <- c("Clin", "A", "B", "C", "D")
  
  # 3-1 Save observed blocks per fold in a overall list!
  all_folds <- list("fold1" = fold1_obs, "fold2" = fold2_obs,
                    "fold3" = fold3_obs, "fold4" = fold4_obs)
  
  
  # [4] Start the CV, split data to Test and Train and evaluate it! ------------
  for (i in 0:4) {
    
    print(paste0("FOLD: ", as.character(i + 1), "/5 -------------------------"))
    
    # [1] Get TestSet from 'data', by taking the first 'obs_per_fold$amount_test' 
    #     IDs in 'fold_ids'
    test_ids <- fold_ids[((i * obs_per_fold$amount_test) + 1):(((i + 1) * obs_per_fold$amount_test))]
    test_df  <- data$data[test_ids,]
    
    # [2] Get the TrainSet from 'data' [= IDs not in TestSet]
    train_ids <- fold_ids[-which(fold_ids %in% test_ids)]
    train_df  <- data$data[train_ids,]
    
    # [3] Induce blockwise missingness [SCENARIO_3] 
    # 3-1 Sample equally sized 'observed' blocks [according to SCENARIO_3]
    set.seed(seed + i)
    observed_blocks <- sample(c(rep("fold1", obs_per_fold$amount_train_fold), 
                                rep("fold2", obs_per_fold$amount_train_fold),
                                rep("fold3", obs_per_fold$amount_train_fold), 
                                rep("fold4", obs_per_fold$amount_train_fold)),
                              obs_per_fold$amount_train, replace = FALSE)
    
    # 3-2 Split Traindata into observed blocks! The resulting blocks will only
    #     contain the features, in the blocks!
    # 3-2-1 Fold1
    block1 <- train_df[which(observed_blocks == "fold1"), 
                       c(response, 
                         unlist(data$block_names[fold1_obs]))]
    
    # 3-2-2 Fold2
    block2 <- train_df[which(observed_blocks == "fold2"), 
                       c(response, 
                         unlist(data$block_names[fold2_obs]))]
    
    # 3-2-3 Fold3
    block3 <- train_df[which(observed_blocks == "fold3"), 
                       c(response, 
                         unlist(data$block_names[fold3_obs]))]
    
    # 3-2-4 Fold4
    block4 <- train_df[which(observed_blocks == "fold4"), 
                       c(response, 
                         unlist(data$block_names[fold4_obs]))]
    
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
    
    # 4-4 BLOCK3 - grow the trees [as long, as all of them are grown correctly]
    trees3 <- simpleRF(formula = formula_all, data = block3, 
                       num_trees = num_trees, mtry = mtry, 
                       min_node_size = min_node_size, replace = TRUE,  
                       splitrule = NULL, unordered_factors = unorderd_factors)
    trees3 <- mclapply(trees3, function(x) {
      x$grow(replace = TRUE)
      x
    }, mc.cores = 1)
    
    # 4-5 BLOCK4 - grow the trees [as long, as all of them are grown correctly]
    trees4 <- simpleRF(formula = formula_all, data = block4, 
                       num_trees = num_trees, mtry = mtry, 
                       min_node_size = min_node_size, replace = TRUE,  
                       splitrule = NULL, unordered_factors = unorderd_factors)
    trees4 <- mclapply(trees4, function(x) {
      x$grow(replace = TRUE)
      x
    }, mc.cores = 1)
    
    # [5] Check, that all of the trees were grown correctly & create a forest from it!
    trees1 <- all_trees_grown_correctly(trees1)
    trees2 <- all_trees_grown_correctly(trees2)
    trees3 <- all_trees_grown_correctly(trees3)
    trees4 <- all_trees_grown_correctly(trees4)
    
    Forest <- list(trees1, trees2, trees3, trees4)
    
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
                                      testdata = test_df[,-which(colnames(test_df) %in% data$block_names$cnv_block)])
    rm(curr_Forest); gc()
    
    curr_Forest      <- copy_forrest(Forest)
    miss1_B[[i + 1]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                      weight_metric = weight_metric,
                                      testdata = test_df[,-which(colnames(test_df) %in% data$block_name$rna_block)])
    rm(curr_Forest); gc()
    
    curr_Forest      <- copy_forrest(Forest)
    miss1_C[[i + 1]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                      weight_metric = weight_metric,
                                      testdata = test_df[,-which(colnames(test_df) %in% data$block_name$mutation_block)])
    rm(curr_Forest); gc()
    
    curr_Forest      <- copy_forrest(Forest)
    miss1_D[[i + 1]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                      weight_metric = weight_metric,
                                      testdata =  test_df[,-which(colnames(test_df) %in% data$block_name$mirna_block)])
    rm(curr_Forest); gc()
    
    # 6-3 TESTSET TWO OMICS BLOCKS MISSING - two ommic blocks are missing in
    #     TestData, everytime before evaluation we need to copy the Forest, 
    #     as the  evaluation can leed to pruned trees [also outside of the function] 
    print("Evaluation TestSet w/ 2 missing omics blocks-----------------------")
    
    curr_Forest       <- copy_forrest(Forest)
    miss2_CD[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                       weight_metric = weight_metric,
                                       testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$mutation_block,
                                                                                           data$block_name$mirna_block))])
    rm(curr_Forest); gc()
    
    curr_Forest       <- copy_forrest(Forest)
    miss2_BD[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                       weight_metric = weight_metric,
                                       testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$rna_block,
                                                                                           data$block_name$mirna_block))])
    rm(curr_Forest); gc()
    
    curr_Forest       <- copy_forrest(Forest)
    miss2_BC[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                       weight_metric = weight_metric,
                                       testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$rna_block,
                                                                                           data$block_name$mutation_block))])
    rm(curr_Forest); gc()
    
    curr_Forest       <- copy_forrest(Forest)
    miss2_AD[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted, 
                                       weight_metric = weight_metric,
                                       testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$cnv_block,
                                                                                           data$block_name$mirna_block))])
    rm(curr_Forest); gc()
    
    curr_Forest       <- copy_forrest(Forest)
    miss2_AC[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                       weight_metric = weight_metric,
                                       testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$cnv_block,
                                                                                           data$block_name$mutation_block))])
    rm(curr_Forest); gc()
    
    curr_Forest       <- copy_forrest(Forest)
    miss2_AB[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                       weight_metric = weight_metric,
                                       testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$cnv_block,
                                                                                           data$block_name$rna_block))])
    
    # 6-4 TESTSET THREE OMICS BLOCKS MISSING - three omic blocks are missing in
    #     TestData, everytime before evaluation we need to copy the Forest, 
    #     as the  evaluation can leed to pruned trees [also outside of the function] 
    print("Evaluation TestSet w/ 3 missing omics blocks-----------------------")
    
    curr_Forest        <- copy_forrest(Forest)
    miss3_ABC[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                        weight_metric = weight_metric,
                                        testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$cnv_block,
                                                                                            data$block_name$rna_block,
                                                                                            data$block_name$mutation_block))])
    rm(curr_Forest); gc()
    
    curr_Forest        <- copy_forrest(Forest)
    miss3_ACD[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted, 
                                        weight_metric = weight_metric,
                                        testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$cnv_block,
                                                                                            data$block_name$mutation_block,
                                                                                            data$block_name$mirna_block))])
    rm(curr_Forest); gc()
    
    curr_Forest        <- copy_forrest(Forest)
    miss3_ABD[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                        weight_metric = weight_metric,
                                        testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$cnv_block,
                                                                                            data$block_name$rna_block,
                                                                                            data$block_name$mirna_block))])
    rm(curr_Forest); gc()
    
    curr_Forest        <- copy_forrest(Forest)
    miss3_BCD[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                        weight_metric = weight_metric,
                                        testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$rna_block,
                                                                                            data$block_name$mutation_block,
                                                                                            data$block_name$mirna_block))])
    rm(curr_Forest); gc()
    
    
    # 7-5 TESTSET THREE OMICS BLOCKS MISSING - three omic blocks are missing in
    #     TestData, everytime before evaluation we need to copy the Forest, 
    #     as the  evaluation can leed to pruned trees [also outside of the function] 
    print("Evaluation TestSet w/ only 1 observed Block -----------------------")
    
    curr_Forest        <- copy_forrest(Forest)
    single_A[[i + 1]]  <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                        weight_metric = weight_metric,
                                        testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$rna_block,
                                                                                            data$block_name$mutation_block,
                                                                                            data$block_name$mirna_block,
                                                                                            data$block_names$clin_block))])
    rm(curr_Forest); gc()
    
    curr_Forest        <- copy_forrest(Forest)
    single_B[[i + 1]]  <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                        weight_metric = weight_metric,
                                        testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$cnv_block,
                                                                                            data$block_name$mutation_block,
                                                                                            data$block_name$mirna_block,
                                                                                            data$block_names$clin_block))])
    rm(curr_Forest); gc()
    
    curr_Forest        <- copy_forrest(Forest)
    single_C[[i + 1]]  <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                        weight_metric = weight_metric,
                                        testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$rna_block,
                                                                                            data$block_name$cnv_block,
                                                                                            data$block_name$mirna_block,
                                                                                            data$block_names$clin_block))])
    rm(curr_Forest); gc()
    
    curr_Forest        <- copy_forrest(Forest)
    single_D[[i + 1]]  <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                        weight_metric = weight_metric,
                                        testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$rna_block,
                                                                                            data$block_name$mutation_block,
                                                                                            data$block_name$cnv_block,
                                                                                            data$block_names$clin_block))])
    rm(curr_Forest); gc()
    
    curr_Forest         <- copy_forrest(Forest)
    single_CL[[i + 1]]  <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                         weight_metric = weight_metric,
                                         testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$rna_block,
                                                                                             data$block_name$mutation_block,
                                                                                             data$block_name$mirna_block,
                                                                                             data$block_names$cnv_block))])
    rm(curr_Forest); gc()
  }
  
  # [5] Return the metric & settings of the fitting! ---------------------------
  # 4-1 Collect all CV Results in a list!
  res_all <- list("full" = full,
                  "miss1_A" = miss1_A, "miss1_B" = miss1_B,
                  "miss1_C" = miss1_C, "miss1_D" = miss1_D,
                  "miss2_CD" = miss2_CD, "miss2_BD" = miss2_BD,
                  "miss2_BC" = miss2_BC, "miss2_AD" = miss2_AD,
                  "miss2_AC" = miss2_AC, "miss2_AB" = miss2_AB,
                  "miss3_ABC" = miss3_ABC, "miss3_ABD" = miss3_ABD,
                  "miss3_ACD" = miss3_ACD, "miss3_BCD" = miss3_BCD,
                  "single_A" = single_A, "single_B" = single_B,
                  "single_C" = single_C, "single_D" = single_D,
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
                   "observed_blocks_in_folds"    = all_folds) # which blocks were 
                                                              # observed in which fold!
  
  # 4-3 Return both lists!
  return(list("res_all"  = res_all, 
              "settings" = settings))
}
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

