"
Script to analyse the TCGA multi-omics data set with the single-block approach
"
# [0] Set WD, load packages & define functions
setwd("C:/Users/kuche_000/Desktop/MS_omics/")
library(assertthat)
library(randomForestSRC)

mcc_metric        <- function(conf_matrix) {
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
load_CV_data      <- function(path) {
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
evaluate_RF       <- function(model, test_set) {
  "
  Evaluate 'model' on the 'test_set'. 
  
  Args:
    - model (rfsrc) : RF that shall predict on 'test_set' & be evaluated
    - test_set (DF) : DF to recieve predicitons on - must conatin a 'gender'
                      column that is binary!
    
  Return:
    - res (vector) : Vector filled with metrics [Accuracy, F-1 Score, ...]
  "
  # [0] Check Inputs  ----------------------------------------------------------
  # 0-1 'model' is rfsrc
  if (!('rfsrc' %in% class(model))) stop("Model not of type 'rfsrc'")
  
  # 0-2 test_set is DF w/ a binary 'gender' column!
  assert_data_frame(test_set)
  if (!("gender" %in% colnames(test_set))) stop("test_set has no 'gender' column")
  if (nlevels(test_set$gender) != 2) stop("gender in test_set is not binary!")
  
  # [1] Get Predicitions for the 'test_set'   ----------------------------------
  # 1-1 Define a Vector to save the predicted classes - fill it with the 
  #     opposite of the true class! Needed later, when predicition on the 
  #     test-set, not all test-obs. can be predicted - these are rated as
  #     wrongly classified for the calculation of the metrics!
  predicted <- sapply(test_set$gender, FUN = function(x) {
    base::ifelse(x == 1, yes = 0, no = 1)
  })
  
  # 1-2 Delete all features from 'test_set' that were not used in training the RF
  # 1-2-1 check whether any variable the RF has been trained with is in test?!
  # 1-2-1-1 If there is no overlap -> can not create predictions!
  #         --> all metrics have worst value then!
  if (!(any(model$xvar.names %in% colnames(test_set)))) {
    
    res <- list("Accuracy"    = 0,
                "Kappa"       = 0,
                "Sensitifity" = 0,
                "Specificity" = 0,
                "Precision"   = 0,
                "Recall"      = 0,
                "F1"          = 0,
                "Balance_Acc" = 0,
                "Pos_Pred_Value" =  0,
                "Neg_Pred_Value" =  0,
                "Prevalence"  = 0,      
                "AUC1"        = 0,
                "AUC2"        = 0,
                "MCC"         = -1)
    return(res)
  }
  
  # 1-2-1-2 If there is a overlap reduce the test_set, so it only contains
  #         Variables the model has been trained with
  test_set <- test_set[, c("gender", model$xvar.names)]
  
  # 1-3 Create predictions for all observations w/o missing features
  CC_test  <- which(complete.cases(test_set))
  
  predicitions       <- predict(RF, test_set)
  predicted[CC_test] <- as.character(predicitions$class)
  
  # [2] Evaluate the predicitions  ---------------------------------------------
  predicted <- as.factor(predicted)
  levels(predicted) <- c(0, 1)
  
  # 2-1 Confusion-Matrix
  confmat <- caret::confusionMatrix(as.factor(predicted), 
                                    as.factor(test_set$gender))
  
  # 2-2 Are under the ROC Curve
  if (length(unique(predicted)) == 1) {
    roc1 <- 0
    roc2 <- 0
  } else {
    roc1 <- pROC::auc(pROC::roc(as.numeric(as.character(predicted)), 
                                as.numeric(as.character(test_set$gender)),
                                levels = levels(as.factor(as.character(test_set$gender))),
                                direction = "<"))
    
    roc2 <- pROC::auc(pROC::roc(as.numeric(as.character(predicted)), 
                                as.numeric(as.character(test_set$gender)),
                                levels = levels(as.factor(as.character(test_set$gender))),
                                direction = ">"))
  }
  
  # 2-3 MCC Matthews correlation coefficient [only for binary cases!]
  mcc <- mcc_metric(conf_matrix = confmat)
  
  # 2-4 Collect all metrics in a list & replace the not defined values
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
  
  return(res)
}

do_CV_NK_5_blocks     <- function(path = "data/processed/RH_subsetted_12345/missingness_1234/BLCA_1.RData",
                                  num_trees = 300, mtry = NULL, min_node_size = 5) {
  "CrossValidate the single-block Approach when the Traindata has blockwise 
   missingness according to scenario 1, 2 or 3!
   
   'path' must lead to a list with 2 entrances: 'data' & 'block_names'
     - 'data' is a list filled with 'k' test-train-splits
        --> k-fold-Validation on this test-train-splits!
     - 'block_names' is a list filled with the names of the single blocks 
        & must be ['A', 'B', 'C', 'D', 'clin_block']!
        (Attention: With Scenario2 the order is different, but this is wanted!)
      
   Based on the 'k' test-train-splits in 'data', we will fit RFs to the single
   feature-blocks in the train data (that has blockwise missingness in it). 
   Then based on the RF that has been trained on a single feature-block, 
   predicitions for the testset are created! If a predicition is not possible
   for any observation, the prediction for these observations is labelld as wrong!
   
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
  
  # 0-2 'num_trees', 'min_node_size' & 'mtry' must be integer > 0 if NOT NULL
  if (!is.null(num_trees)) assert_int(num_trees, lower = 10)
  if (!is.null(mtry)) assert_int(mtry)
  if (!is.null(min_node_size)) assert_int(min_node_size, lower = 1)
  
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
  unique_fea_blocks <- length(names(curr_data$block_names))
  
  # 1-2-1 Full TestSet
  full        <- vector(mode = "list", length = unique_fea_blocks)
  names(full) <- names(curr_data$block_names)
  
  # 1-2-2 TestSet with 1 missing omics-block
  miss1_A        <- vector(mode = "list", length = unique_fea_blocks)
  names(miss1_A) <- names(curr_data$block_names)
  
  miss1_B        <- vector(mode = "list", length = unique_fea_blocks)
  names(miss1_B) <- names(curr_data$block_names)
  
  miss1_C        <- vector(mode = "list", length = unique_fea_blocks)
  names(miss1_C) <- names(curr_data$block_names)
  
  miss1_D        <- vector(mode = "list", length = unique_fea_blocks)
  names(miss1_D) <- names(curr_data$block_names)
  
  # 1-2-3 TestSet with 2 missing omics-blocks
  miss2_CD <- vector(mode = "list", length = unique_fea_blocks)
  names(miss2_CD) <- names(curr_data$block_names)
  
  miss2_BD <- vector(mode = "list", length = unique_fea_blocks)
  names(miss2_BD) <- names(curr_data$block_names)
  
  miss2_BC <- vector(mode = "list", length = unique_fea_blocks)
  names(miss2_BC) <- names(curr_data$block_names)
  
  miss2_AD <- vector(mode = "list", length = unique_fea_blocks)
  names(miss2_AD) <- names(curr_data$block_names)
  
  miss2_AC <- vector(mode = "list", length = unique_fea_blocks)
  names(miss2_AC) <- names(curr_data$block_names)
  
  miss2_AB <- vector(mode = "list", length = unique_fea_blocks)
  names(miss2_AB) <- names(curr_data$block_names)
  
  # 1-2-4 TestSet with 3 missing omics-blocks
  miss3_ABC <- vector(mode = "list", length = unique_fea_blocks)
  names(miss3_ABC) <- names(curr_data$block_names)
  
  miss3_ABD <- vector(mode = "list", length = unique_fea_blocks)
  names(miss3_ABD) <- names(curr_data$block_names)
  
  miss3_ACD <- vector(mode = "list", length = unique_fea_blocks)
  names(miss3_ACD) <- names(curr_data$block_names)
  
  miss3_BCD <- vector(mode = "list", length = unique_fea_blocks)
  names(miss3_BCD) <- names(curr_data$block_names)
  
  # 1-2-5 Single BlockTestSet [4 missing blocks!]
  single_A <- vector(mode = "list", length = unique_fea_blocks)
  names(single_A) <- names(curr_data$block_names)
  
  single_B <- vector(mode = "list", length = unique_fea_blocks)
  names(single_B) <- names(curr_data$block_names)
  
  single_C <- vector(mode = "list", length = unique_fea_blocks)
  names(single_C) <- names(curr_data$block_names)
  
  single_D <- vector(mode = "list", length = unique_fea_blocks)
  names(single_D) <- names(curr_data$block_names)
  
  single_CL <- vector(mode = "list", length = unique_fea_blocks)
  names(single_CL) <- names(curr_data$block_names)
  
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
    
    # 2-3 Loop over each feature-block in train & fit a RF and create preds for 
    #     the test-set based on the RF fitted on a single feature-block
    for (curr_block in names(curr_data$block_names)) {
      
      # --1 Get all variables from 'train' that belong to 'curr_block'
      curr_block_feas <- curr_data$block_names[[curr_block]]
      
      # --2 Only keep the observations w/o missing values in 'curr_block' to 
      #     train the model
      curr_train <- train[,c('gender', curr_block_feas)]
      curr_train <- curr_train[complete.cases(curr_train),]
      
      # --3 Train a RF on 'curr_train'
      # - define formula
      curr_formula <- as.formula('gender ~ .')
      
      # - fit RF on 'curr_train'
      RF <- rfsrc(formula = curr_formula, data = curr_train, 
                  ntree = num_trees, mtry = mtry, nodesize = min_node_size, 
                  samptype = "swr", seed = 12345678, var.used = 'all.trees')
      
      # --4 Evaluate the Model on the different Test.Sets
      # --4-1 FullTestSet
      print("Evaluation on full TestSet --------------------------------------")
      full[[curr_block]][[i]] <- evaluate_RF(model = RF, test_set = test)
      
      # --4-1 One Missing FeatureBlock in Test
      print("Evaluation TestSet w/ 2 missing Blocks --------------------------")
      miss1_A[[curr_block]][[i]] <- evaluate_RF(model = RF, 
                                                test_set = test[,-which(colnames(test) %in% curr_data$block_names$A)])
      
      miss1_B[[curr_block]][[i]] <- evaluate_RF(model = RF, 
                                                test_set = test[,-which(colnames(test) %in% curr_data$block_names$B)])
      
      miss1_C[[curr_block]][[i]] <- evaluate_RF(model = RF, 
                                                test_set = test[,-which(colnames(test) %in% curr_data$block_names$C)])
      
      miss1_D[[curr_block]][[i]] <- evaluate_RF(model = RF, 
                                                test_set = test[,-which(colnames(test) %in% curr_data$block_names$D)])
      
      # --4-2 Two Missing FeatureBlock in Test
      print("Evaluation TestSet w/ 2 missing Blocks --------------------------")
      miss2_CD[[curr_block]][[i]] <- evaluate_RF(model = RF, 
                                                 test_set = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                               curr_data$block_names$D))])
      
      miss2_BD[[curr_block]][[i]] <- evaluate_RF(model = RF, 
                                                 test_set = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                               curr_data$block_names$D))])
      
      miss2_BC[[curr_block]][[i]] <- evaluate_RF(model = RF, 
                                                 test_set = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                               curr_data$block_names$B))])
      
      miss2_AD[[curr_block]][[i]] <- evaluate_RF(model = RF, 
                                                 test_set = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                               curr_data$block_names$D))])
      
      miss2_AC[[curr_block]][[i]] <- evaluate_RF(model = RF, 
                                                 test_set = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                               curr_data$block_names$A))])
      
      miss2_AB[[curr_block]][[i]] <- evaluate_RF(model = RF, 
                                                 test_set = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                               curr_data$block_names$B))])
      
      # --4-3 Three Missing FeatureBlock in Test
      print("Evaluation TestSet w/ 3 missing Blocks --------------------------")
      miss3_ABC[[curr_block]][[i]] <- evaluate_RF(model = RF,
                                                  test_set = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                                curr_data$block_names$A,
                                                                                                curr_data$block_names$B))])
      
      miss3_ABD[[curr_block]][[i]] <- evaluate_RF(model = RF,
                                                  test_set = test[,-which(colnames(test) %in% c(curr_data$block_names$D,
                                                                                                curr_data$block_names$A,
                                                                                                curr_data$block_names$B))])
      
      miss3_ACD[[curr_block]][[i]] <- evaluate_RF(model = RF,
                                                  test_set = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                                curr_data$block_names$A,
                                                                                                curr_data$block_names$D))])
      
      miss3_BCD[[curr_block]][[i]] <- evaluate_RF(model = RF,
                                                  test_set = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                                curr_data$block_names$B,
                                                                                                curr_data$block_names$D))])
      
      # --4-4 Single FeatureBlock in Test
      print("Evaluation TestSet w/ only 1 observed Block -----------------------")
      single_A[[curr_block]][[i]] <- evaluate_RF(model = RF,
                                                 test_set = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                               curr_data$block_names$D,
                                                                                               curr_data$block_names$B,
                                                                                               curr_data$block_names$clin_block))])
      
      single_B[[curr_block]][[i]] <- evaluate_RF(model = RF,
                                                 test_set = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                               curr_data$block_names$D,
                                                                                               curr_data$block_names$A,
                                                                                               curr_data$block_names$clin_block))])
      
      single_C[[curr_block]][[i]] <- evaluate_RF(model = RF,
                                                 test_set = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                               curr_data$block_names$D,
                                                                                               curr_data$block_names$B,
                                                                                               curr_data$block_names$clin_block))])
      
      single_D[[curr_block]][[i]] <- evaluate_RF(model = RF,
                                                 test_set = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                               curr_data$block_names$A,
                                                                                               curr_data$block_names$B,
                                                                                               curr_data$block_names$clin_block))])
      
      single_CL[[curr_block]][[i]] <- evaluate_RF(model = RF,
                                                  test_set = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                                curr_data$block_names$A,
                                                                                                curr_data$block_names$B,
                                                                                                curr_data$block_names$D))])
    }
  }
  
  # 1-5 Stop the time and take the difference!
  time_for_CV <- Sys.time() - start_time
  
  # [2] Return the results & settings of parameters used to do CV! -------------
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
                   "response"      = "gender",
                   "num_trees"     = num_trees,
                   "mtry"          = mtry, 
                   "min_node_size" = min_node_size,
                   "time_for_CV"   = time_for_CV)
  
  # 3-3 Return both lists!
  return(list("res_all"  = res_all, 
              "settings" = settings))
}

# MAIN                                                                      ----
DFs_w_gender <- c("COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC","LGG", "BLCA",
                  "LUAD", "LUSC", "PAAD", "SARC", "SKCM", "STAD")

# ----- Situation 1
for (DF in DFs_w_gender) {
  
  print(paste0("----- Situation 1 for DF: '", DF, "' -----"))
  
  # Create the path for the current DF
  curr_path <- paste0("data/processed/RH_subsetted_12345/missingness_1234/", DF, "_1.RData")
  
  
  sit1 <- do_CV_NK_5_blocks(path = curr_path, num_trees = 300, mtry = NULL, 
                            min_node_size = 5)
  save(sit1, file = paste0("./docs/CV_Res/TCGA/SingleBlock_Approach/setting1/", DF, ".RData"))
}

# ----- Situation 2
for (DF in DFs_w_gender) {
  
  print(paste0("----- Situation 1 for DF: '", DF, "' -----"))
  
  # Create the path for the current DF
  curr_path <- paste0("data/processed/RH_subsetted_12345/missingness_1234/", DF, "_2.RData")
  
  
  sit2 <- do_CV_NK_5_blocks(path = curr_path, num_trees = 300, mtry = NULL, 
                            min_node_size = 5)
  save(sit2, file = paste0("./docs/CV_Res/TCGA/SingleBlock_Approach/setting2/", DF, ".RData"))
}

# ----- Situation 3
for (DF in DFs_w_gender) {
  
  print(paste0("----- Situation 1 for DF: '", DF, "' -----"))
  
  # Create the path for the current DF
  curr_path <- paste0("data/processed/RH_subsetted_12345/missingness_1234/", DF, "_3.RData")
  
  
  sit3 <- do_CV_NK_5_blocks(path = curr_path, num_trees = 300, mtry = NULL, 
                            min_node_size = 5)
  save(sit3, file = paste0("./docs/CV_Res/TCGA/SingleBlock_Approach/setting3/", DF, ".RData"))
}