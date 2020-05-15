" Script to crossValidate the adaption of N. Krautenbachers's block-wise RF Algorithm!
 
  This approach learns a seperate RF on each of the feature blocks in the train-
  set. For predicitions on testdata it only uses the RFs that use avaible features
  in testdata! 
  --> Then to obtain a final predicton we combine all created predicitons!
"
# SetWD and define/load functions
library(randomForestSRC)
library(checkmate)
library(caret)

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
do_evaluation_rfsrc   <- function(Forest, testdata, weights) {
  " For a collection of blockwise fitted RFs [each RF fitted on diff. block] 
    get the aggregated predicition from all trees!
    Evaluate the these on 'testdata' & return metrics!
  
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
      - list w/ metrics [accuracy, f1, mcc, roc, ...] 
        > Metrics w/ NA are replaced by their worst possible value!
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
  test_cols <- colnames(testdata)
  
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
  if (length(Forest) < 1) stop("Forest can not predicit on TestData, as all 
                                trees use split vars not avaible in 'testdata'")
  
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

  # [4] From Probabilities to classes! -----------------------------------------
  all_forrest_preds_class <- ifelse(prob_class0 >= 0.5, 0, 1)
  all_forrest_preds_class <- factor(all_forrest_preds_class, 
                                    levels = levels(Forest[[1]]$yvar))
  
  # [5] Get Metrics for the current setting!  ----------------------------------
  # 5-1 Confusion Matrix
  confmat <- caret::confusionMatrix(data      = all_forrest_preds_class, 
                                    reference = testdata[,1])
  
  # 5-2 Are under the ROC Curve
  roc1 <- tryCatch(pROC::auc(pROC::roc(as.numeric(testdata[,1]), 
                                       as.numeric(all_forrest_preds_class),
                                       levels = unique(as.numeric(testdata[,1])),
                                       direction = "<")),
                   error = function(e) "not defined!")
  if (is.numeric(roc1)) roc1 <- as.numeric(roc1)
  
  roc2 <- tryCatch(pROC::auc(pROC::roc(as.numeric(testdata[,1]), 
                                       as.numeric(all_forrest_preds_class),
                                       levels = unique(as.numeric(testdata[,1])),
                                       direction = ">")),
                   error = function(e) "not defined")
  if (is.numeric(roc2)) roc2 <- as.numeric(roc2)
  
  # 5-3 MCC Matthews correlation coefficient [only for binary cases!]
  mcc <- mcc_metric(conf_matrix = confmat)
  
  # 5-4 Get the Variables used as split variables!
  used_split_vars <- sapply(Forest, function(x) x$var.used)
  
  # [6] Create a list to collect the results!  ---------------------------------
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
  
  return(as.vector(res))
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

do_CV_NK_5_blocks     <- function(path = "data/processed/TCGA_subset_12345/missingness_1234/BLCA_1.RData",
                                  weighted = TRUE, weight_metric = NULL,
                                  num_trees = 300, mtry = NULL, 
                                  min_node_size = 5) {
  "CrossValidate the Approach when the Traindata has blockwise missingness
   according to scenario 1, 2 or 3!
   
   'path' must lead to a list with 2 entrances: 'data' & 'block_names'
     - 'data' is a list filled with 'k' test-train-splits
        --> k-fold-Validation on this test-train-splits!
     - 'block_names' is a list filled with the names of the single blocks 
        & must be ['A', 'B', 'C', 'D', 'clin_block']!
        (Attention: With Scenario2 the order is different, but this is wanted!)
      
   Based on the 'k' test-train-splits in 'data', we will fit blokwise RFs to the
   train data (that has blockwise missingness in it). Then we ensemble the 
   predicitons from blockwise fitted RFs to single predicitons & rate these with 
   Accuracy, Precision, Specifity, F1-Score,...
   
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
      - weighted (bool)     : When ensembling the prediciton from the blockwise
                              RFs shall we weight the predictions by the OOB 
                              performance of the foldwise fitted RFs? 
                              The better the OOB-Metric for a blockwise fitted RF, 
                              the higher its contribution to the final prediction!
      - weight_metric (chr) : Which metric to use to assigning weights to the 
                              different predictions? 
                                - must be 'Acc' or 'F1' 
                                - if weighted = FALSE, it will be ignored!
      - num_trees (int)     : Amount of trees, we shall grow on each block!
                                      > If NULL: 1000 is the default!
      - mtry (int)          : Amount of split-variables we try, when looking for 
                              a split variable! 
                                      > If 'NULL': mtry = sqrt(p)
      - min_node_size (int) : Amount of Observations a node must at least contain
                              so the model keeps on trying to split them!  
                                      > If NULL: min_node_size = 1 
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
            - datapath, response, mtry, time for CV
  "
  # [0] Check Inputs -----------------------------------------------------------
  # 0-1 path must be string and have '1.RData' | '2.RData' | '3.RData' in it!
  assert_string(path)
  if (!grepl("1.RData", path) & !grepl("2.RData", path) & !grepl("3.RData", path)) {
    stop("'path' must end in '1.RData' | '2.RData' | '3.RData'")
  }
  
  # 0-2 weighted must be a single boolean
  assert_logical(weighted, len = 1)
  
  # 0-3 'num_trees', 'min_node_size' & 'mtry' must be integer > 0 if NOT NULL
  if (!is.null(num_trees)) assert_int(num_trees, lower = 10)
  if (!is.null(mtry)) assert_int(mtry)
  if (!is.null(min_node_size)) assert_int(min_node_size, lower = 1)
  
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
  # 1-1 Load CV-Data [already splitted - data checked in 'load_CV_data' itself]
  curr_data <- load_CV_data(path = path)
  
  # 1-1-1 Must contain 'A', 'B', 'C' 'D' & 'clin_block' as block_names
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
  
  # 1-4 Start the Timer, so we know how long CV took!
  start_time <- Sys.time()
  
  # [2] Start the CV and loop over all test-train splits in data  --------------
  for (i in seq_len(k_splits)) {
    
    # 2-1 Current Fold Status
    print(paste0("FOLD ", i, "/", k_splits, " -------------------------------"))
    
    # 2-2 Extract the test and train set from 'curr_data'
    train <- curr_data$data[[i]]$train
    test  <- curr_data$data[[i]]$test
   
    # 2-3 Train blockwise RFs for each block seperatly! For this loop over all 
    #     blocks ['curr_data$block_names']
    Forest <- list(); i_  <- 1
    
    for (block_ in names(curr_data$block_names)) {
      
      # 2-3-1 Get the Observations that have the features of 'block_'
      observed <- sapply(seq_len(nrow(train)), function(j_) {
        sum(is.na(train[j_, curr_data$block_names[[block_]]])) == 0
      })
      
      # 2-3-2 Convert Boolean to indicies!
      observed <- which(observed)
      
      # 2-3-3 Fit a RF on this fully observed (fold-)subdata!
      # 2-3-3-1 Define formula
      response    <- colnames(train)[1]
      formula_all <- as.formula(paste(response, " ~ ."))

      # 2-3-3-2 Get all Obs. from the current block, its block features and the 
      #         corresponding response in a seperate DF!
      curr_fold_train_data <- train[observed, 
                                    c(response, curr_data$block_names[[block_]])]
      
      # 2-3-3-3 Fit a Tree on the block data
      blockwise_rf <- rfsrc(formula = formula_all, data = curr_fold_train_data, 
                            ntree = num_trees, mtry = mtry, nodesize = min_node_size, 
                            samptype = "swr", seed = 12345678, var.used = 'all.trees')
      
      print(paste0("Fit FoldWise RF on current block: '", block_, "'"))
      
      Forest[[i_]] <- blockwise_rf
      rm(blockwise_rf)
      
      # 2-3-5 Count up i_ - used to fill Forest list!
      i_  = i_ + 1
    }
    
    # 2-4 Get the OOB metrics of the blockwise fitted RFs
    # 2-4-1 If weighted, loop over the blockwise RFs and get the oob F1 / Acc
    #       of each blockwise fitted RF!
    if (weighted) {
      weights <- c()
      
      for (block_RF in Forest) {
        weights <- c(weights, 
                     get_oob_weight_metric(block_RF)[weight_metric])
      }
    } else {
      weights <- rep(1, times = length(Forest))
    }
    
    # 2-4-2 Norm the weights - if not all are 0 - if all are 0, 
    #       give each block same weights
    if (any(weights != 0)) weights <- weights / sum(weights)
    if (all(weights == 0)) weights <-  rep(1, times = length(Forest))
    
    # 2-5 Evaluate the RF on the different Testsets! Fromfull TestSet w/o any 
    #     missing blocks to TestSets w/ only one observed Block!
    # 2-5-1 Full TestSet!
    print("Evaluation full TestSet -------------------------------------------")
    full[[i]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                     testdata = test) 
    
    # 2-5-2 TestSet with 1 missing block!
    print("Evaluation TestSet w/ 1 missing omics block------------------------")
    miss1_A[[i]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                        testdata = test[,-which(colnames(test) %in% curr_data$block_names$A)])
    
    miss1_B[[i]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                        testdata = test[,-which(colnames(test) %in% curr_data$block_names$B)])
    
    miss1_C[[i]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                        testdata = test[,-which(colnames(test) %in% curr_data$block_names$C)])
    
    miss1_D[[i]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                        testdata = test[,-which(colnames(test) %in% curr_data$block_names$D)])
    
    # 2-5-3 TestSet with 2 missing blocks!
    print("Evaluation TestSet w/ 2 missing omics blocks-----------------------")
    miss2_CD[[i]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                         testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                       curr_data$block_names$D))])
    miss2_BD[[i]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                         testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                       curr_data$block_names$D))])
    miss2_BC[[i]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                         testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                       curr_data$block_names$B))])
    miss2_AD[[i]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                         testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                       curr_data$block_names$D))])
    miss2_AC[[i]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                         testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                       curr_data$block_names$A))])
    miss2_AB[[i]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                         testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                       curr_data$block_names$B))])
    # 2-5-4 Testset with 3 missing blocks!
    print("Evaluation TestSet w/ 3 missing omics blocks-----------------------")
    miss3_ABC[[i]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                         testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                       curr_data$block_names$A,
                                                                                       curr_data$block_names$B))])
    miss3_ACD[[i]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                          testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                        curr_data$block_names$A,
                                                                                        curr_data$block_names$D))])
    miss3_ABD[[i]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                          testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$D,
                                                                                        curr_data$block_names$A,
                                                                                        curr_data$block_names$B))])
    miss3_BCD[[i]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                          testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                        curr_data$block_names$D,
                                                                                        curr_data$block_names$B))])
    # 2-5-5 Evaluation on single Block Testdata
    print("Evaluation TestSet w/ only 1 observed Block -----------------------")
    single_A[[i]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                        testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                      curr_data$block_names$D,
                                                                                      curr_data$block_names$B,
                                                                                      curr_data$block_names$clin_block))])
    single_B[[i]] <-  do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                          testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                        curr_data$block_names$D,
                                                                                        curr_data$block_names$A,
                                                                                        curr_data$block_names$clin_block))])
    single_C[[i]] <-  do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                          testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                        curr_data$block_names$D,
                                                                                        curr_data$block_names$B,
                                                                                        curr_data$block_names$clin_block))])
    single_D[[i]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                         testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                       curr_data$block_names$A,
                                                                                       curr_data$block_names$B,
                                                                                       curr_data$block_names$clin_block))])
    single_CL[[i]] <-  do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                           testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                         curr_data$block_names$D,
                                                                                         curr_data$block_names$B,
                                                                                         curr_data$block_names$A))])
  }
  # 2-6 Stop the time and take the difference!
  time_for_CV <- Sys.time() - start_time
  
  # [3] Return the results & settings of parameters used to do CV! -------------
  # 3-1 Collect all CV Results in a list!
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
  
  # 3-2 Collect the Settings, used to do the CV!
  settings <- list("data_path"     = path,
                   "num_folds"     = k_splits,
                   "response"      = response,
                   "weighted"      = weighted,
                   "weight_metric" = weight_metric,
                   "num_trees"     = num_trees,
                   "mtry"          = mtry, 
                   "min_node_size" = min_node_size,
                   "time_for_CV"   = time_for_CV)
  
  # 3-3 Return both lists!
  return(list("res_all"  = res_all, 
              "settings" = settings))
}

do_CV_NK_3_blocks     <- function(path = "data/processed/TCGA_subset_12345/missingness_1234/BLCA_4.RData",
                                  weighted = TRUE, weight_metric = NULL,
                                  num_trees = as.integer(300), mtry = NULL, 
                                  min_node_size = 5) {
  "CrossValidate the Approach when the Traindata has blockwise missingness
   according to scenario 4!
   
   'path' must lead to a list with 2 entrances: 'data' & 'block_names'
       - 'data' is a list filled with 'k' test-train-splits
          --> k-fold-Validation on this test-train-splits!
       - 'block_names' is a list filled with the names of the single blocks 
          & must be ['A', 'B', 'clin_block']!
      
   Based on the 'k' test-train-splits in 'data', we will fit blokwise RFs to the
   train data (that has blockwise missingness in it). Then we ensemble the 
   predicitons from blockwise fitted RFs to single predicitons & rate these with 
   Accuracy, Precision, Specifity, F1-Socre,...
   
   The TestingSituations are different, as we can test the models on fully 
   observed testdata, on testdata w/ 1 missing block, etc...
    --> Results is list with all results from the k test-train splits for all 
        possible testsituations - 6 in total!
   
  Args:
      - path (char)         : path to data with blockwise missingness we want to CV
                              Must end in '4.RData'
                              --> List with 2 entrances: 'data' & 'block_names':
                                  - 'data' consitis of 'k' test-train-splits, 
                                     where train has missingness induced and the 
                                     test-set is fully observed!
                                  - 'block_names' contains all colnames of the 
                                     different blocks!
      - weighted (bool)     : When ensembling the prediciton from the blockwise
                              RFs shall we weight the predictions by the OOB 
                              performance of the foldwise fitted RFs? 
                              The better the OOB-Metric for a blockwise fitted RF, 
                              the higher its contribution to the prediction!
      - weight_metric (chr) : Which metric to use to assigning weights to the 
                              different predictions? ['Acc', 'F1']
                                - if weighted = FALSE, it will be ignored!
      - num_trees (int)     : Amount of trees, we shall grow on each block!
                                    > if NULL: 1000 is the default
      - mtry (int)          : Amount of split-variables we try, when looking for 
                              a split variable! 
                                    > If 'NULL': mtry = sqrt(p)
      - min_node_size (int) : Amount of Observations a node must at least contain,
                              so the model keeps on trying to split them!  
                                    > If NULL: 1 is the default!
  Return:
      - list filled w/:
        * 'res_all' [the CV Results on different testsets]:
           >> A = CNV-Block, B = RNA-Block, C = Mutation-Block, D = Mirna-Block
            - full    : CV Results for each fold on the fully observed testdata!
            - miss1_A : CV Results for each fold on the testdata, w/ missing A
                        block [actually 2 blocks randomly thrown together]
            - miss1_B : CV Results for each fold on the testdata, w/ missing B 
                        block  [actually 2 blocks randomly thrown together]
                .
                .
            - single_CL: CV-Results for each fold on the testdata w/ clinical 
                         data only!
                         
        * 'settings' [settings used to do the CV - all arguments!]
            - datapath, response, mtry, time for CV
  "
  # [0] Check Inputs -----------------------------------------------------------
  # 0-1 path must be numeric and have '4.RData' in it!
  assert_string(path, fixed = "4.RData")
  
  # 0-2 weighted must be a single boolean
  assert_logical(weighted, len = 1)
  
  # 0-3 'num_trees', 'min_node_size' & 'mtry' must be integer > 0 if NOT NULL
  if (!is.null(num_trees)) assert_int(num_trees, lower = 10)
  if (!is.null(mtry)) assert_int(mtry)
  if (!is.null(min_node_size)) assert_int(min_node_size, lower = 1)
  
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
  # 1-1 Load CV-Data [already splitted - data checked in 'load_CV_data' itself]
  curr_data <- load_CV_data(path = path)
  
  # 1-1-1 Must contain 'A', 'B' & 'clin_block' as block_names
  corr_block_names <- ("A" %in% names(curr_data$block_names) & 
                       "B" %in% names(curr_data$block_names) &
                       "clin_block" %in% names(curr_data$block_names))
  
  if (!corr_block_names) stop("'path' lead to a file without 'A', 'B', & 'clin_block' as blocknames!")
  
  # 1-2 Create empty lists to store results in!
  # 1-2-1 Full TestSet
  full <- list()
  # 1-2-2 TestSet with 1 missing omics-block
  miss1_A <- list(); miss1_B <- list()
  # 1-2-3 Single BlockTestSet [4 missing blocks!]
  single_A <- list(); single_B <- list(); single_CL <- list()
  
  # 1-3 Get the amount of test-train splits in data 
  k_splits <- length(curr_data$data)
  
  # 1-4 Start the Timer, so we know how long CV took!
  start_time <- Sys.time()
  
  # [2] Start the CV and loop over all test-train splits in data  --------------
  for (i in seq_len(k_splits)) {
    
    # 2-1 Current Fold Status
    print(paste0("FOLD ", i, "/", k_splits, " -------------------------------"))
    
    # 2-2 Extract the test and train set from 'curr_data'
    train <- curr_data$data[[i]]$train
    test  <- curr_data$data[[i]]$test
    
    # 2-3 Train blockwise RFs for each block seperatly! For this loop over all 
    #     blocks ['curr_data$block_names']
    Forest <- list(); i_  <- 1
    
    for (block_ in names(curr_data$block_names)) {
      
      # 2-3-1 Get the Observations that have the features of 'block_'
      observed <- sapply(seq_len(nrow(train)), function(j_) {
        sum(is.na(train[j_, curr_data$block_names[[block_]]])) == 0
      })
      
      # 2-3-2 Convert Boolean to indicies!
      observed <- which(observed)
      
      # 2-3-3 Fit a RF on this fully observed (fold-)subdata!
      # 2-3-3-1 Define formula
      response    <- colnames(train)[1]
      formula_all <- as.formula(paste(response, " ~ ."))
      
      # 2-3-3-2 Get all Obs. from the current block, its block features and the 
      #         corresponding response in a seperate DF!
      curr_fold_train_data <- train[observed, 
                                    c(response, curr_data$block_names[[block_]])]
      
      # 2-3-3-3 Fit a Tree on the block data
      blockwise_rf <- rfsrc(formula = formula_all, data = curr_fold_train_data, 
                            ntree = num_trees, mtry = mtry, nodesize = min_node_size, 
                            samptype = "swr", seed = 12345678, var.used = 'all.trees')
      
      print(paste0("Fit FoldWise RF on current block: '", block_, "'"))
      
      Forest[[i_]] <- blockwise_rf
      rm(blockwise_rf)
      
      # 2-3-5 Count up i_ - used to fill Forest list!
      i_  = i_ + 1
    }
    
    # 2-4 Get the OOB metrics of the blockwise fitted RFs
    # 2-4-1 If weighted, loop over the blockwise RFs and get the oob F1/ Acc
    #       of each blockwise fitted RF!
    if (weighted) {
      weights <- c()
      
      for (block_RF in Forest) {
        weights <- c(weights, 
                     get_oob_weight_metric(block_RF)[weight_metric])
      }
    } else {
      weights <- rep(1, times = length(Forest))
    }
    
    # 2-4-2 Norm the weights - if not all are 0 - if all are 0, 
    #       give each block same weights
    if (any(weights != 0)) weights <- weights / sum(weights)
    if (all(weights == 0)) weights <-  rep(1, times = length(Forest))
    
    # 2-5 Evaluate the RF on the different Testsets! Fromfull TestSet w/o any 
    #     missing blocks to TestSets w/ only one observed Block!
    # 2-5-1 Full TestSet!
    print("Evaluation full TestSet -------------------------------------------")
    full[[i]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                     testdata = test)
    
    # 2-5-2 TestSet with 1 missing block!
    print("Evaluation TestSet w/ 1 missing omics block------------------------")
    miss1_A[[i]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                        testdata = test[,-which(colnames(test) %in% curr_data$block_names$A)])
    
    miss1_B[[i]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                        testdata = test[,-which(colnames(test) %in% curr_data$block_names$B)])
    
    # 2-5-3 Evaluation on single Block Testdata
    print("Evaluation TestSet w/ only 1 observed Block -----------------------")
    single_A[[i]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                         testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                       curr_data$block_names$clin_block))])
    single_B[[i]] <-  do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                          testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                        curr_data$block_names$clin_block))])
    single_CL[[i]] <-  do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                           testdata = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                         curr_data$block_names$B))])
  }
  # 2-6 Stop the time and take the difference!
  time_for_CV <- Sys.time() - start_time
  
  # [3] Return the results & settings of parameters used to do CV! -------------
  # 3-1 Collect all CV Results in a list!
  res_all <- list("full" = full,
                  "miss1_A" = miss1_A, "miss1_B" = miss1_B,
                  "single_A" = single_A, "single_B" = single_B,
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
                   "time_for_CV"   = time_for_CV)
  
  # 4-3 Return both lists!
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
  curr_path <- paste0("data/processed/TCGA_subset_12345/missingness_1234/", DF, "_1.RData")
  
  print("Setting - 1/3")
  sit1_1 <- do_CV_NK_5_blocks(path = curr_path,
                              weighted = FALSE, weight_metric = NULL,
                              num_trees = 300, mtry = NULL, min_node_size = 5)
  save(sit1_1, file = paste0("./docs/CV_Res/TCGA/BlockWise_Approach/setting1/", DF, ".RData"))
  
  print("Setting - 2/3")
  sit1_2 <- do_CV_NK_5_blocks(path = curr_path,
                              weighted = TRUE, weight_metric = "F1",
                              num_trees = 300, mtry = NULL, min_node_size = 5)
  save(sit1_2, file = paste0("./docs/CV_Res/TCGA/BlockWise_Approach/setting1/", DF, "_f1.RData"))
  
  print("Setting - 3/3")
  sit1_3 <- do_CV_NK_5_blocks(path = curr_path,
                              weighted = TRUE, weight_metric = "Acc",
                              num_trees = 300, mtry = NULL, min_node_size = 5)
  save(sit1_3, file = paste0("./docs/CV_Res/TCGA/BlockWise_Approach/setting1/", DF, "_acc.RData"))
}

# ----- Situation 2
for (DF in DFs_w_gender) {
  
  print(paste0("----- Situation 2 for DF: '", DF, "' -----"))
  
  # Create the path for the current DF
  curr_path <- paste0("data/processed/TCGA_subset_12345/missingness_1234/", DF, "_2.RData")
  
  print("Setting - 1/3")
  sit2_1 <- do_CV_NK_5_blocks(path = curr_path,
                              weighted = FALSE, weight_metric = NULL,
                              num_trees = 300, mtry = NULL, min_node_size = 5)
  save(sit2_1, file = paste0("./docs/CV_Res/TCGA/BlockWise_Approach/setting2/", DF, ".RData"))
  
  print("Setting - 2/3")
  sit2_2 <- do_CV_NK_5_blocks(path = curr_path,
                              weighted = TRUE, weight_metric = "F1",
                              num_trees = 300, mtry = NULL, min_node_size = 5)
  save(sit2_2, file = paste0("./docs/CV_Res/TCGA/BlockWise_Approach/setting2/", DF, "_f1.RData"))
  
  print("Setting - 3/3")
  sit2_3 <- do_CV_NK_5_blocks(path = curr_path,
                              weighted = TRUE, weight_metric = "Acc",
                              num_trees = 300, mtry = NULL, min_node_size = 5)
  save(sit2_3, file = paste0("./docs/CV_Res/TCGA/BlockWise_Approach/setting2/", DF, "_acc.RData"))
}

# ----- Situation 3
for (DF in DFs_w_gender) {
  
  print(paste0("----- Situation 3 for DF: '", DF, "' -----"))
  
  # Create the path for the current DF
  curr_path <- paste0("data/processed/TCGA_subset_12345/missingness_1234/", DF, "_3.RData")
  
  print("Setting - 1/3")
  sit3_1 <- do_CV_NK_5_blocks(path = curr_path,
                              weighted = FALSE, weight_metric = NULL,
                              num_trees = 300, mtry = NULL, min_node_size = 5)
  save(sit3_1, file = paste0("./docs/CV_Res/TCGA/BlockWise_Approach/setting3/", DF, ".RData"))
  
  print("Setting - 2/3")
  sit3_2 <- do_CV_NK_5_blocks(path = curr_path,
                              weighted = TRUE, weight_metric = "F1",
                              num_trees = 300, mtry = NULL, min_node_size = 5)
  save(sit3_2, file = paste0("./docs/CV_Res/TCGA/BlockWise_Approach/setting3/", DF, "_f1.RData"))
  
  print("Setting - 3/3")
  sit3_3 <- do_CV_NK_5_blocks(path = curr_path,
                              weighted = TRUE, weight_metric = "Acc",
                              num_trees = 300, mtry = NULL, min_node_size = 5)
  save(sit3_3, file = paste0("./docs/CV_Res/TCGA/BlockWise_Approach/setting3/", DF, "_acc.RData"))
}

# ----- Situation 4
for (DF in DFs_w_gender) {
  
  print(paste0("----- Situation 4 for DF: '", DF, "' -----"))
  
  # Create the path for the current DF
  curr_path <- paste0("data/processed/TCGA_subset_12345/missingness_1234/", DF, "_4.RData")
  
  print("Setting - 1/3")
  sit4_1 <- do_CV_NK_3_blocks(path = curr_path,
                              weighted = FALSE, weight_metric = NULL,
                              num_trees = 300, mtry = NULL, min_node_size = 5)
  save(sit4_1, file =  paste0("./docs/CV_Res/TCGA/BlockWise_Approach/setting4/", DF, ".RData"))
  
  print("Setting - 2/3")
  sit4_2 <- do_CV_NK_3_blocks(path = curr_path,
                              weighted = TRUE, weight_metric = "F1",
                              num_trees = 300, mtry = NULL, min_node_size = 5)
  save(sit4_2, file = paste0("./docs/CV_Res/TCGA/BlockWise_Approach/setting4/", DF, "_f1.RData"))
  
  print("Setting - 3/3")
  sit4_3 <- do_CV_NK_3_blocks(path = curr_path,
                              weighted = TRUE, weight_metric = "Acc",
                              num_trees = 300, mtry = NULL, min_node_size = 5)
  save(sit4_3, file = paste0("./docs/CV_Res/TCGA/BlockWise_Approach/setting4/", DF, "_acc.RData"))
}