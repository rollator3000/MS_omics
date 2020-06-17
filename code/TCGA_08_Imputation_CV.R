"Script to tackle the blockwise missingness with the missFores Imputation Approach!
  - Impute the missing values, in TrainData with blockwise missingness!
  - Train Classifier on the imputed DF - only use the blocks the test contains!
"
# SetWD and define/load functions
setwd("C:/Users/kuche_000/Desktop/MS-Thesis/")
library(missForest)
library(checkmate)
library(randomForestSRC)
library(tidyverse)
library(caret)
library(doParallel)

# Set Cores for parallel computaion!
detectCores()
registerDoParallel(cores = 6)

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
impute_train_data <- function(path, ntree_imp = 25, maxiter = 1, para = "forests") {
  " Impute the missing values in the trainsets of the test-train splits 'path'
    points to! Impute the missing values with the missForest approach!
    Then save the imputed DF, so we do not need to impute it the whole time and 
    can evaluate different approaches withit!
    
    Look into the 'train-data' of each split and get the variables we need to 
    impute! Then apply the missForest approach to impute these missing values!
     >> recode numeric columns to factors, if the numeric column has less than 5
        unique numeric values --> converted back to numeric before returning! <<
   
  Args:
    - path (str)       : Path to a list of Test-Train Splits, whereby the Train-
                         Set has block-wise missing values - must be '.RData'!
                         The data we load is checked within 'load_CV_data()'
    - ntree_imp (int)  : Amount of trees we use for the missForest imputation
    - maxiter (int)    : Max Number of iterations to be performed - as long 
                         as the stopping criterion is not fullfilled yet!
    - para (str)      :  Should 'missForest' be run parallel? If 'variables' the
                         data is split into pieces of the size equal to the
                         number of cores registered in the parallel backend. 
                         If 'forests' the total number of trees in each random
                         forests is split in the same way. 
                          > 'no'; 'forests' or 'variables'
                         
  Return:
    - data_all (list)  : Same structure of the loaded Test-Train Splits, but the 
                         Training set has no more missing values, as all were
                         imputed by missForest!
  "
  # [0] Check Inputs & load the data  ------------------------------------------
  # 0-1 Check for meaningful arguments
  assert_string(path, fixed = ".RData")
  assert_int(maxiter, lower = 1)
  assert_int(ntree_imp, lower = 10)
  if (!(para %in% c('no', 'forests', 'variables'))) {
    stop("'para' must be: 'no'; 'forests'; or 'variables'")
  }
  
  # 0-2 Load data and extract the data-splits itself! 
  curr_data <- load_CV_data(path = path)
  data_all  <- curr_data$data
  
  # [1] Run the impuation for each of the 5 Test-Train Splits  -----------------
  for (curr_ind in 1:length(data_all)) {
    
    print(paste("Imputation", curr_ind, "/", length(data_all)))
    
    # 1-1 Extract the current train data!
    curr_train <- data_all[[curr_ind]]$train
    
    # 1-2 Remove the response from the train data, as this shall not be used 
    #     for the imputation
    curr_resp      <- curr_train[,1]
    curr_resp_name <- colnames(curr_train)[1]
    curr_train[,1] <- NULL
    
    # 1-3 Recode numeic variables to factor if the numeric < 5 unique values!
    # 1-3-1 Count amount of factor levels in each variable
    unique_values_cols <- sapply(1:ncol(curr_train), 
                                 function(x) length(unique(curr_train[,x])))
    
    # 1-3-2 Get the cols with less than 5 unique values!
    cols_to_recode <- which(unique_values_cols < 5)
    
    # 1-3-3 Convert numeric variables to factor variables
    #       - if there even is a col to recode!
    if (length(cols_to_recode) > 1) {
      curr_train[cols_to_recode] <- lapply(curr_train[cols_to_recode], 
                                           function(x) as.factor(x))
    }
    
    # 1-4 Start the Imputation!
    data_imputed <- missForest(curr_train, verbose = TRUE, 
                               parallelize = para, maxiter = maxiter, 
                               ntree = ntree_imp)
    data_imputed <- data_imputed$ximp
    
    # 1-5 Check that the imputed data has no more missing values!
    if (any(sapply(1:nrow(data_imputed), function(x) sum(is.na(x)) > 1))) {
      stop("Imputed DF still has missing datapoints!")
    }
    
    # 1-6 Recode the variables back to numeric if we coded them to factors in 1-3
    if (length(cols_to_recode) > 1) {
      data_imputed[cols_to_recode] <- lapply(data_imputed[cols_to_recode], as.numeric)
    }
    
    # 1-7 Add the response again into the first column and assign correct name!
    data_imputed              <- cbind(curr_resp, data_imputed)
    colnames(data_imputed)[1] <- curr_resp_name
    
    # 1-8 Replace the train data with block-wise missingness with the imputed
    #     Training data
    data_all[[curr_ind]]$train <- data_imputed
    
    # 1-9 Save the prt that has been imputed!
    name_save <- strsplit(path, split = "1234/")[[1]][2]
    save_path <- "data/processed/RH_subsetted_12345/missingness_1234_imputed/"
    save_path <- paste0(save_path, curr_ind, "_", name_save)
    save(data_imputed, file = save_path)
  }
  
  # [2] Return the Test-Train Splits  ------------------------------------------
  #     Add the imputed test-train splits to curr_data [overwrite old data]!
  #     And return it! --> train data imputed all the time
  #                    --> 
  #     Test-Train Splits as they have been, but the Traindata does not contain 
  #     anymore missing data!
  curr_data$data <- data_all
  return(curr_data)
}
impute_train_data_foldwise <- function(path, ntree_imp = 25, maxiter = 1, 
                                       para = "forests", start_fold = 1) {
  " Impute the missing values in the trainsets of the test-train splits 'path'
    points to! Impute the missing values with the missForest approach!
    Then save the imputed DF, so we do not need to impute it the whole time and 
    can evaluate different approaches withit!
    
    Look into the 'train-data' of each split and get the variables we need to 
    impute! Then apply the missForest approach to impute these missing values!
     >> recode numeric columns to factors, if the numeric column has less than 5
        unique numeric values --> converted back to numeric before returning! <<
   
  Args:
    - path (str)       : Path to a list of Test-Train Splits, whereby the Train-
                         Set has block-wise missing values - must be '.RData'!
                         The data we load is checked within 'load_CV_data()'
    - ntree_imp (int)  : Amount of trees we use for the missForest imputation
    - maxiter (int)    : Max Number of iterations to be performed - as long 
                         as the stopping criterion is not fullfilled yet!
    - para (str)       : Should 'missForest' be run parallel? If 'variables' the
                         data is split into pieces of the size equal to the
                         number of cores registered in the parallel backend. 
                         If 'forests' the total number of trees in each random
                         forests is split in the same way. 
                          > 'no'; 'forests' or 'variables'
    - start_fold (int) : Which fold to start the imputation - all lower are 
                         ignored! 
  Return:                Save the imputed data of all folds in 
                         'data/processed/RH_subsetted_12345/missingness_1234_imputed/'
                         fold_data_imputed.R
  "
  # [0] Check Inputs & load the data  ------------------------------------------
  # 0-1 Check for meaningful arguments
  assert_string(path, fixed = ".RData")
  assert_int(maxiter, lower = 1)
  assert_int(ntree_imp, lower = 10)
  if (!(para %in% c('no', 'forests', 'variables'))) {
    stop("'para' must be: 'no'; 'forests'; or 'variables'")
  }
  assert_int(start_fold, lower = 1, upper = 5)
  
  # 0-2 Load data and extract the data-splits itself! 
  curr_data <- load_CV_data(path = path)
  data_all  <- curr_data$data
  
  # [1] Run the impuation for each of the 5 Test-Train Splits  -----------------
  for (curr_ind in start_fold:length(data_all)) {
    
    print(paste("Imputation", curr_ind, "/", length(data_all)))
    
    # 1-1 Extract the current train data!
    curr_train <- data_all[[curr_ind]]$train
    
    # 1-2 Remove the response from the train data, as this shall not be used 
    #     for the imputation
    curr_resp      <- curr_train[,1]
    curr_resp_name <- colnames(curr_train)[1]
    curr_train[,1] <- NULL
    
    # 1-3 Recode numeic variables to factor if the numeric < 5 unique values!
    # 1-3-1 Count amount of factor levels in each variable
    unique_values_cols <- sapply(1:ncol(curr_train), 
                                 function(x) length(unique(curr_train[,x])))
    
    # 1-3-2 Get the cols with less than 5 unique values!
    cols_to_recode <- which(unique_values_cols < 5)
    
    # 1-3-3 Convert numeric variables to factor variables
    #       - if there even is a col to recode!
    if (length(cols_to_recode) > 1) {
      curr_train[cols_to_recode] <- lapply(curr_train[cols_to_recode], 
                                           function(x) as.factor(x))
    }
    
    # 1-4 Start the Imputation!
    data_imputed <- missForest(curr_train, verbose = TRUE, 
                               parallelize = para, maxiter = maxiter, 
                               ntree = ntree_imp)
    data_imputed <- data_imputed$ximp
    
    # 1-5 Check that the imputed data has no more missing values!
    if (any(sapply(1:nrow(data_imputed), function(x) sum(is.na(x)) > 1))) {
      stop("Imputed DF still has missing datapoints!")
    }
    
    # 1-6 Recode the variables back to numeric if we coded them to factors in 1-3
    if (length(cols_to_recode) > 1) {
      data_imputed[cols_to_recode] <- lapply(data_imputed[cols_to_recode], as.numeric)
    }
    
    # 1-7 Add the response again into the first column and assign correct name!
    data_imputed              <- cbind(curr_resp, data_imputed)
    colnames(data_imputed)[1] <- curr_resp_name
    
    # 1-8 Replace the train data with block-wise missingness with the imputed
    #     Training data
    data_all[[curr_ind]]$train <- data_imputed
    
    # 1-9 Save the prt that has been imputed!
    name_save <- strsplit(path, split = "1234/")[[1]][2]
    save_path <- "data/processed/TCGA_subset_12345/missingness_1234_imputed/"
    save_path <- paste0(save_path, curr_ind, "_", name_save)
    save(data_imputed, file = save_path)
  }
}
do_evaluation_imputed       <- function(train, test, num_trees, min_node_size, mtry) {
  "Evaluate the Imputation Approach! MissForest was used to impute the missing 
   Values in the training data. Evaluate the Approach on 'test'. 
   For this we supply the features -only the ones avaible in 'test'- from 'train'
   to train a regular RF on it. Then we evaluate the predicitve performance of
   the RF trained on imputed DF on the testdata and take the Acc, F1, ...
   
   Args:
    - train (data.frame) : Dataframe that has no missing data, as all missing 
                           datapoints were imputed!
    - test (data.frame)  : Dataframe we use to test the impuation approach!
                           - Must not contain any variables not avaible in 
                             'train'.
    - num_trees (int) : Amount of trees to be fit on the imputed 'train'!
                        --> RF will be fit on all feas in Train, that are also
                            observed in 'test'
                        If NULL: 1000 is the default!
    - min_node_size(int) : Amount of Observations a node must at least contain,
                           so the model keeps on trying to split them!  
                        If NULL: 1 is the default!
    - mtry (int)         : Amount of split-variables we try, when looking for
                           a split variable! 
                        If 'NULL': mtry = sqrt(p)
                        
   Return:
    - list w/ metrics [accuracy, f1, mcc, roc, ...] 
        if F-1 Score is NA [as precision or recall = 0], we set it to 0 [not NA]
  "
  # [0] Check Inputs  ----------------------------------------------------------
  # 0-1 'train' & 'test' must be dataframes
  assert_data_frame(train)
  assert_data_frame(test)
  
  # 0-2 'test' must not contain any colnames not avaible in train
  if (any(!(colnames(test) %in% colnames(train)))) {
    stop("Test-Set has different features than the Train-Set!")
  }
  
  # 0-3 'num_trees', 'min_node_size' & 'mtry' must be integer > 0 if NOT NULL
  if (!is.null(num_trees)) {
    assert_int(num_trees, lower = 10)
  } else {
    num_trees = 1000
  }
  if (!is.null(mtry)) assert_int(mtry, lower = 5)
  if (!is.null(min_node_size)) assert_int(min_node_size, lower = 1)
  
  # [1] Train a RF  ------------------------------------------------------------
  #     Train a RF on 'train', but only use the variables, that are avaible in
  #     'test', as the RF must do predicitons on the 'test' set!
  # 1-1 Remove the Features from Train that are not avaible in Test
  test_col_names  <- colnames(test)
  curr_train_data <- train[,test_col_names]
  
  # 1-2 Train a RF on the 'curr_train_data'
  # 1-2-1 Extract Response and create the formula the RF to be fit on!
  response    <- colnames(train)[1]
  formula_all <- as.formula(paste(response, " ~ ."))
  
  # 1-2-2 Fit the actual RF
  RF <- rfsrc(formula = formula_all, data = curr_train_data, 
              ntree = num_trees, mtry = mtry, nodesize = min_node_size, 
              samptype = "swr", seed = 12345678, var.used = 'all.trees')
  
  # 1-3 Get Prediciton on the testset from the RF
  predicitons <- predict(RF, test)
  
  # [2] Get Metrics of predicitive power  --------------------------------------
  # 2-1 ConfusionMatrix
  confmat <- caret::confusionMatrix(predicitons$class, test[,1])
  
  # 2-2 Are under the ROC Curve
  roc1 <- pROC::auc(pROC::roc(as.numeric(test[,1]), 
                              as.numeric(predicitons$class),
                              levels = unique(as.numeric(test[,1])),
                              direction = "<"))
  
  roc2 <- pROC::auc(pROC::roc(as.numeric(test[,1]), 
                              as.numeric(predicitons$class),
                              levels = unique(as.numeric(test[,1])),
                              direction = ">"))
  
  # 2-3 MCC Matthews correlation coefficient [only for binary cases!]
  mcc <- mcc_metric(conf_matrix = confmat)
  
  # [3] Create a list to collect the results!  ---------------------------------
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
              "MCC"         = mcc,
              "Selected_Vars" = RF$var.used)
  
  # 3-1 If the F1-Score/ Precision/ Recall is NA, then we set it to 0 [-1 for MCC]! 
  if (is.na(res$F1))             res$F1             <- 0
  if (is.na(res$Precision))      res$Precision      <- 0
  if (is.na(res$Recall))         res$Recall         <- 0
  if (is.na(res$MCC))            res$MCC            <- -1
  if (is.na(res$Pos_Pred_Value)) res$Pos_Pred_Value <- 0
  if (is.na(res$Neg_Pred_Value)) res$Neg_Pred_Value <- 0
  
  return(as.vector(res))
}
do_CV_missforrest_5         <- function(path = "data/processed/TCGA_subset_12345/missingness_1234_imputed/BLCA_IMP_1.RData",
                                        num_trees = 300, min_node_size = 5, mtry = NULL) {
  "Evalute the Approach where the block-wise missingness in the DFs is removed by
   missforest impuation method!
   
   'path' must lead to a list with 2 entrances: 'data' & 'block_names'
   - 'data' is a list filled with 'k' test-train-splits
        --> the 'train' parts have been imputed with 'impute_train_data()'
        --> the 'test' parts are completly observed
   - 'block_names' is a list filled with the names of the single blocks 
      & must be ['A', 'B', 'C', 'D', 'clin_block']!
      Attention: With Scenario2 the order is different, but this is wanted!
      
   Based on the 'k' test-train-splits in 'data', we will evaluate the imputation
   approach. For this the block-wise missing datapoints have been imputed already
   and a random forest model is fit on these imputed DFs and creates predicitions
   on the test data then!

   Then for each testsituation (fully observed testset,.., single block testset) 
   we prune the imputed data, so that only features remain that are also in the 
   test set. On this data a RF is fitted and the performance is measured in the 
   testset & rated with Accuracy, Precision, Recall, F1-Socre,...
   
   Args:
    - path (str) : path to the test-train splits, whereby the train parts were 
                   already imputed! 
                   --> list with the 2 entrances 'data' & 'block_names'
                       - 'data' consitis of 'k' test-train-splits
                       - 'block_names' consits of the feature names for the 
                          different feature-blocks
                   --> must contain 'imputed' as the imputed DFs are saved in 
                       such a corresponding folder!
    - num_trees (int) : Amount of trees to be fit on the imputed data!
                        --> RF will be fit on all feas in Train, that are also
                            observed in 'Test'
    - min_node_size(int) : Amount of Observations a node must at least contain,
                           so the model keeps on trying to split them!  
    - mtry (int)         : Amount of split-variables we try, when looking for
                           a split variable! >> if 'NULL': mtry = sqrt(p)
   
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
            - single_D: CV-Results for each fold on the testdata w/ mirna block 
                        only! [or in scenario2, what was sampled to be block 'D']
                         
        * 'settings' [settings used to do the CV - all arguments!]
            - datapath, seed, response, mtry, time for 
  "
  # [0] Check Inputs  ----------------------------------------------------------
  # 0-0 mtry, min_node_size & num_trees are all checked within 'rfsrc()'
  # 0-1 path must be string - rest is checked in 'load_CV_data()'
  assert_string(path, pattern = "imputed")
  if (!grepl("1.RData", path) & !grepl("2.RData", path) & !grepl("3.RData", path)) {
    stop("'path' must end in '1.RData' | '2.RData' | '3.RData'")
  }
  
  # [1] Get the data and create list to save results  --------------------------
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
  miss2_CD <- list(); miss2_BD <- list(); miss2_BC <- list(); miss2_AD <- list(); miss2_AC <- list(); miss2_AB <- list()
  # 1-2-4 TestSet with 3 missing omics-blocks
  miss3_ABC <- list(); miss3_ABD <- list(); miss3_ACD <- list(); miss3_BCD <- list()
  # 1-2-5 Single BlockTestSet
  single_A <- list(); single_B <- list(); single_CL <- list(); single_C <- list(); single_D <- list()
  
  # 1-3 Get the amount of Splits
  k_splits <- length(curr_data$data)
  
  # 1-4 Start a timer, so we can calc. how long the CV took in total!
  start_time <- Sys.time()
  
  # [2] Start the CV  ----------------------------------------------------------
  #     Run over each Test-Train Split, impute the missing data from the Train-
  #     data, and fit a RF on it. Evaluate this Approach for different observed
  #     Testdata [e.g. fully observed <-> miss_AB]
  for (i in 1:k_splits) {
    
    # 2-1 Current Fold Status
    print(paste0("FOLD ", i, "/", k_splits, " -------------------------------"))
    
    # 2-2 Extract the test and train set from 'curr_data'
    train <- curr_data$data[[i]]$train
    test  <- curr_data$data[[i]]$test
    
    # 2-3 Evaluate the RF on the different Testsets! Fromfull TestSet w/o any 
    #     missing blocks to TestSets w/ only one observed Block!
    # 2-4-1 Full TestSet!
    print("Evaluation full TestSet -------------------------------------------")
    full[[i]] <- do_evaluation_imputed(train = train, test = test, 
                                       num_trees = num_trees, 
                                       min_node_size = min_node_size, 
                                       mtry = mtry)
    
    # 2-4-2 TestSet with 1 missing block!
    print("Evaluation TestSet w/ 1 missing omics block------------------------")
    miss1_A[[i]] <- do_evaluation_imputed(train = train, 
                                          test = test[,-which(colnames(test) %in% curr_data$block_names$A)], 
                                          num_trees = num_trees, 
                                          min_node_size = min_node_size, 
                                          mtry = mtry)
    
    miss1_B[[i]] <- do_evaluation_imputed(train = train, 
                                          test = test[,-which(colnames(test) %in% curr_data$block_names$B)], 
                                          num_trees = num_trees, 
                                          min_node_size = min_node_size, 
                                          mtry = mtry)
    
    miss1_C[[i]] <- do_evaluation_imputed(train = train, 
                                          test = test[,-which(colnames(test) %in% curr_data$block_names$C)], 
                                          num_trees = num_trees, 
                                          min_node_size = min_node_size, 
                                          mtry = mtry)
    
    miss1_D[[i]] <- do_evaluation_imputed(train = train, 
                                          test = test[,-which(colnames(test) %in% curr_data$block_names$D)], 
                                          num_trees = num_trees, 
                                          min_node_size = min_node_size, 
                                          mtry = mtry)
    
    # 2-4-3 TestSet with 2 missing blocks!
    print("Evaluation TestSet w/ 2 missing omics blocks-----------------------")
    miss2_CD[[i]] <- do_evaluation_imputed(train = train, 
                                           test = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                     curr_data$block_names$D))], 
                                           num_trees = num_trees, 
                                           min_node_size = min_node_size, 
                                           mtry = mtry)
    
    miss2_BD[[i]] <- do_evaluation_imputed(train = train, 
                                           test = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                     curr_data$block_names$D))], 
                                           num_trees = num_trees, 
                                           min_node_size = min_node_size, 
                                           mtry = mtry)
    
    miss2_BC[[i]] <- do_evaluation_imputed(train = train, 
                                           test = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                     curr_data$block_names$B))], 
                                           num_trees = num_trees, 
                                           min_node_size = min_node_size, 
                                           mtry = mtry)
    
    miss2_AB[[i]] <- do_evaluation_imputed(train = train, 
                                           test = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                     curr_data$block_names$B))], 
                                           num_trees = num_trees, 
                                           min_node_size = min_node_size, 
                                           mtry = mtry)
    
    miss2_AC[[i]] <- do_evaluation_imputed(train = train, 
                                           test = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                     curr_data$block_names$A))], 
                                           num_trees = num_trees, 
                                           min_node_size = min_node_size, 
                                           mtry = mtry)
    
    miss2_AD[[i]] <- do_evaluation_imputed(train = train, 
                                           test = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                     curr_data$block_names$D))], 
                                           num_trees = num_trees, 
                                           min_node_size = min_node_size, 
                                           mtry = mtry)
    
    # 2-4-4 Testset with 3 missing blocks!
    print("Evaluation TestSet w/ 3 missing omics blocks-----------------------")
    miss3_ABC[[i]] <- do_evaluation_imputed(train = train, 
                                            test = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                      curr_data$block_names$C,
                                                                                      curr_data$block_names$A))], 
                                            num_trees = num_trees, 
                                            min_node_size = min_node_size, 
                                            mtry = mtry)
    
    miss3_ACD[[i]] <- do_evaluation_imputed(train = train, 
                                            test = test[,-which(colnames(test) %in% c(curr_data$block_names$D,
                                                                                      curr_data$block_names$C,
                                                                                      curr_data$block_names$A))], 
                                            num_trees = num_trees, 
                                            min_node_size = min_node_size, 
                                            mtry = mtry)
    
    miss3_ABD[[i]] <- do_evaluation_imputed(train = train, 
                                            test = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                      curr_data$block_names$D,
                                                                                      curr_data$block_names$A))], 
                                            num_trees = num_trees, 
                                            min_node_size = min_node_size, 
                                            mtry = mtry)
    
    miss3_BCD[[i]] <- do_evaluation_imputed(train = train, 
                                            test = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                      curr_data$block_names$C,
                                                                                      curr_data$block_names$D))], 
                                            num_trees = num_trees, 
                                            min_node_size = min_node_size, 
                                            mtry = mtry)
    
    # 2-4-5 Evaluation on single Block Testdata
    print("Evaluation TestSet w/ only 1 observed Block -----------------------")
    single_A[[i]] <- do_evaluation_imputed(train = train, 
                                           test = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                     curr_data$block_names$C,
                                                                                     curr_data$block_names$D,
                                                                                     curr_data$block_names$clin_block))], 
                                           num_trees = num_trees, 
                                           min_node_size = min_node_size, 
                                           mtry = mtry)
    
    single_B[[i]] <-  do_evaluation_imputed(train = train, 
                                            test = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                      curr_data$block_names$C,
                                                                                      curr_data$block_names$D,
                                                                                      curr_data$block_names$clin_block))], 
                                            num_trees = num_trees, 
                                            min_node_size = min_node_size, 
                                            mtry = mtry)
    
    single_C[[i]] <-  do_evaluation_imputed(train = train, 
                                            test = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                      curr_data$block_names$A,
                                                                                      curr_data$block_names$D,
                                                                                      curr_data$block_names$clin_block))], 
                                            num_trees = num_trees, 
                                            min_node_size = min_node_size, 
                                            mtry = mtry)
    
    single_D[[i]] <- do_evaluation_imputed(train = train, 
                                           test = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                     curr_data$block_names$C,
                                                                                     curr_data$block_names$A,
                                                                                     curr_data$block_names$clin_block))], 
                                           num_trees = num_trees, 
                                           min_node_size = min_node_size, 
                                           mtry = mtry)
    
    single_CL[[i]] <- do_evaluation_imputed(train = train, 
                                            test = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                      curr_data$block_names$C,
                                                                                      curr_data$block_names$D,
                                                                                      curr_data$block_names$A))], 
                                            num_trees = num_trees, 
                                            min_node_size = min_node_size, 
                                            mtry = mtry)
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
                   "num_trees"     = num_trees,
                   "mtry"          = mtry, 
                   "min_node_size" = min_node_size,
                   "time_for_CV"   = time_for_CV)
  
  # 3-3 Return both lists!
  return(list("res_all"  = res_all, 
              "settings" = settings))
}

# TO DO TO DO TO DO TO DO TO DO TO DO TO DO TO DO TO DO TO DO TO DO TO DO TO DO  
do_CV_missforrest_3   <- function(path = "data/processed/TCGA_subset_12345/missingness_1234/BLCA_4.RData",
                                  num_trees = 300, min_node_size = 10, mtry = NULL,
                                  n_tree_impute = 100, maxiter_impute = 10) {
  "Evalute the Approach where we impute the block-wise missing values with the 
   missForest Approach - for the cases with 3 blocks! [Scenario4]
   
   'path' must lead to a list with 2 entrances: 'data' & 'block_names'
   - 'data' is a list filled with 'k' test-train-splits
      --> k-fold-Validation on this test-train-splits!
   - 'block_names' is a list filled with the names of the single blocks 
      & must be ['A', 'B', 'clin_block']!
      
   Based on the 'k' test-train-splits in 'data', we will evaluate the imputation
   approach. For this take the train data (that has blockwise missingness in it)
   and impute missing values with missForest.
   Then for each testsituation (fully observed testset,.., single block testset) 
   we prune the imputed data, so that only features that are also in test remain.
   On this data a RF is fitted and the performance is measured in the testset &
   rated with Accuracy, Precision, Specifity, F1-Socre,...
   
   Args:
    - path (str) : path to the data with blockwise missingness we want to CV
                   --> must lead to list with 2 entrances 'data' & 'block_names'
                       'data' consitis of 'k' test-train-splits, where train 
                        has missingness induced and test fully observed!
    - num_trees (int) : Amount of trees to be fit on the imputed data!
                        --> RF will be fit on all feas in Train, that are also
                            observed in 'Test'
                        If NULL: 1000 is the default!
    - min_node_size(int) : Amount of Observations a node must at least contain,
                           so the model keeps on trying to split them!  
                        If NULL: 1 is the default!
    - mtry (int)         : Amount of split-variables we try, when looking for
                           a split variable! 
                        If 'NULL': mtry = sqrt(p)
    - n_tree_impute (int)  : Amount of Trees we use for imputing the 
                             missing values in the 'data$train' section of the
                             data in 'path'
    - maxiter_impute (int) : Maximum amount of itterations for the Imputation 
                             with missForest!
   
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
            - datapath, seed, response, mtry, time for 
  "
  # [0] Check Inputs  ----------------------------------------------------------
  # 0-0 mtry, min_node_size & num_trees are all checked within 'rfsrc()'
  # 0-1 path must be string rest is checked in 'load_CV_data()'
  assert_string(path, fixed = "4.RData")
  
  # 0-2 'n_tree_impute' & 'maxiter_impute' must be numeric!
  assert_int(n_tree_impute)
  assert_int(maxiter_impute)
  
  # [1] Get the data and create list to save results  --------------------------
  # 1-1 Load CV-Data [already splitted - data checked in 'load_CV_data' itself]
  curr_data <- load_CV_data(path = path)
  
  # 1-1-1 Must contain 'A', 'B', 'C' 'D' & 'clin_block' as block_names
  corr_block_names <- ("A" %in% names(curr_data$block_names) & 
                         "B" %in% names(curr_data$block_names) &
                         "clin_block" %in% names(curr_data$block_names))
  
  if (!corr_block_names) stop("'path' lead to a file without 'A', 'B', 'C', 'D' & 'clin_block' as blocknames!")
  
  # 1-2 Create empty lists to store results in!
  # 1-2-1 Full TestSet
  full <- list()
  # 1-2-2 TestSet with 1 missing omics-block
  miss1_A <- list(); miss1_B <- list()
  # 1-2-3 Single BlockTestSet
  single_A <- list(); single_B <- list(); single_CL <- list()
  
  # 1-3 Get the amount of Splits
  k_splits <- length(curr_data$data)
  
  # 1-4 Start a timer, so we can calc. how long the CV took in total!
  start_time <- Sys.time()
  
  # [2] Start the CV  ----------------------------------------------------------
  #     Run over each Test-Train Split, impute the missing data from the Train-
  #     data, and fit a RF on it. Evaluate this Approach for different observed
  #     Testdata [e.g. fully observed <-> miss_AB]
  for (i in 1:k_splits) {
    
    # 2-1 Current Fold Status
    print(paste0("FOLD ", i, "/", k_splits, " -------------------------------"))
    
    # 2-2 Extract the test and train set from 'curr_data'
    train <- curr_data$data[[i]]$train
    test  <- curr_data$data[[i]]$test
    
    # 2-3 Call Imput Function and impute the missing values!
    train_imputed <- impute_missing_values(data = train, 
                                           ntree_imp = n_tree_impute, 
                                           maxiter = maxiter_impute)
    
    # 2-4 Evaluate the RF on the different Testsets! Fromfull TestSet w/o any 
    #     missing blocks to TestSets w/ only one observed Block!
    # 2-4-1 Full TestSet!
    print("Evaluation full TestSet -------------------------------------------")
    full[[i]] <- do_evaluation_imputed(train = train_imputed, test = test, 
                                       num_trees = num_trees, 
                                       min_node_size = min_node_size, 
                                       mtry = mtry)  
    
    # 2-4-2 TestSet with 1 missing block!
    print("Evaluation TestSet w/ 1 missing omics block------------------------")
    miss1_A[[i]] <- do_evaluation_imputed(train = train_imputed, 
                                          test = test[,-which(colnames(test) %in% curr_data$block_names$A)], 
                                          num_trees = num_trees, 
                                          min_node_size = min_node_size, 
                                          mtry = mtry)
    
    miss1_B[[i]] <- do_evaluation_imputed(train = train_imputed, 
                                          test = test[,-which(colnames(test) %in% curr_data$block_names$B)], 
                                          num_trees = num_trees, 
                                          min_node_size = min_node_size, 
                                          mtry = mtry)
    
    # 2-4-3 Evaluation on single Block Testdata
    print("Evaluation TestSet w/ only 1 observed Block -----------------------")
    single_A[[i]] <- do_evaluation_imputed(train = train_imputed, 
                                           test = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                     curr_data$block_names$clin_block))], 
                                           num_trees = num_trees, 
                                           min_node_size = min_node_size, 
                                           mtry = mtry)
    
    single_B[[i]] <-  do_evaluation_imputed(train = train_imputed, 
                                            test = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                      curr_data$block_names$clin_block))], 
                                            num_trees = num_trees, 
                                            min_node_size = min_node_size, 
                                            mtry = mtry)
    
    single_CL[[i]] <- do_evaluation_imputed(train = train_imputed, 
                                            test = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                      curr_data$block_names$A))], 
                                            num_trees = num_trees, 
                                            min_node_size = min_node_size, 
                                            mtry = mtry)
  }
  
  # 2-6 Stop the time and take the difference!
  time_for_CV <- Sys.time() - start_time
  
  # [3] Return the results & settings of parameters used to do CV! -------------
  # 3-1 Collect all CV Results in a list!
  res_all <- list("full" = full, "miss1_A" = miss1_A, 
                  "miss1_B" = miss1_B, "single_A" = single_A, 
                  "single_B" = single_B, "single_CL" = single_CL)
  
  # 3-2 Collect the Settings, used to do the CV!
  settings <- list("data_path"     = path,
                   "num_folds"     = k_splits,
                   "num_trees"     = num_trees,
                   "mtry"          = mtry, 
                   "min_node_size" = min_node_size,
                   "n_tree_impute" = n_tree_impute, 
                   "maxiter_impute" = maxiter_impute,
                   "time_for_CV"   = time_for_CV)
  
  # 3-3 Return both lists!
  return(list("res_all"  = res_all, 
              "settings" = settings))
} 
# MAIN -------------------------------------------------------------------------
# [1] Impute the missing values in Train for all of the DFs
# 1-1 Names of the usable DFs
DFs_w_gender <- c("COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC","LGG", "BLCA",
                  "LUAD", "LUSC", "PAAD", "SARC", "SKCM", "STAD")

# 1-2 Path to the folder with the subsetted DFs 
path         <- "data/processed/TCGA_subset_12345/missingness_1234/"

# 1-3 Path to the folder to save the imputed test-train-splits!
save_path    <- "data/processed/TCGA_subset_12345/missingness_1234_imputed/"

# 1-4 Which block-wise missingness setting [1, 2, 3 or 4]?
setting      <- "2"

for (curr_data in "ESCA") {
  
  # print current DF we do Imputation on!
  print(paste("Start Imputation for:", curr_data))
  
  # paste the single path elements and start imputation
  curr_path    <- paste0(path, curr_data, "_", setting, ".RData")
  
  # start the imputation and take the time
  start_time <- Sys.time()
  imputed_data <- impute_train_data_foldwise(path = curr_path, ntree_imp = 25, maxiter = 1,
                                             para = "forests", start_fold = 1)
  end_time <- Sys.time()
  
  # print the needed time
  print(paste("Needed Time:", end_time - start_time))
}

# [2] Evaluate approaches on the imputed data
#     Loop over all imputed datasets in 'DFs_w_gender'
DFs_w_gender <- c("COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC","LGG", "BLCA",
                  "LUAD", "LUSC", "PAAD", "SARC", "SKCM", "STAD")

# ----- Situation 1
for (curr_data in DFs_w_gender) {
  
  print(paste0("----- Situation 1 for DF: '", curr_data, "' -----------------"))
  
  # --1 Define current path
  curr_path <- paste0("data/processed/TCGA_subset_12345/missingness_1234_imputed/", 
                      curr_data, "_IMP_1.RData")
  
  # --2 Do the CV
  sit1 <- do_CV_missforrest_5(path = curr_path, num_trees = 300, 
                             min_node_size = 5, mtry = NULL)
  
  # --3 Save the results of the CV
  save(sit1, file = paste0("./docs/CV_Res/TCGA/Imputation_Approach/setting1/", curr_data, ".RData"))
}

# ----- Situation 2
for (curr_data in DFs_w_gender) {
  
  print(paste0("----- Situation 2 for DF: '", curr_data, "' -----------------"))
  
  # --1 Define current path
  curr_path <- paste0("data/processed/TCGA_subset_12345/missingness_1234_imputed/", 
                      curr_data, "_IMP_2.RData")
  
  # --2 Do the CV
  sit2 <- do_CV_missforrest_5(path = curr_path, num_trees = 300, 
                              min_node_size = 5, mtry = NULL)
  
  # --3 Save the results of the CV
  save(sit2, file = paste0("./docs/CV_Res/TCGA/Imputation_Approach/setting2/", curr_data, ".RData"))
}

# ----- Situation 3
for (curr_data in DFs_w_gender) {
  
  print(paste0("----- Situation 3 for DF: '", curr_data, "' -----------------"))
  
  # --1 Define current path
  curr_path <- paste0("data/processed/TCGA_subset_12345/missingness_1234_imputed/", 
                      curr_data, "_IMP_2.RData")
  
  # --2 Do the CV
  sit3 <- do_CV_missforrest_5(path = curr_path, num_trees = 300, 
                              min_node_size = 5, mtry = NULL)
  
  # --3 Save the results of the CV
  save(sit3, file = paste0("./docs/CV_Res/TCGA/Imputation_Approach/setting3/", curr_data, ".RData"))
}

# ----- Situation 4 TO DO TO DO TO DO TO DO TO DO TO DO TO DO TO DO TO DO TO DO 
'
for (DF in DFs_w_gender) {
  
  print(paste0("----- Situation 4 for DF: ", DF, "-----"))
  
  # Create the path for the current DF
  curr_path <- paste0("data/processed/TCGA_subset_12345/missingness_1234/", DF, "_4.RData")
  
  sit4 <- do_CV_5_blocks(path = curr_path, num_trees = 300, mtry = NULL, 
                         min_node_size = 5, unorderd_factors = "ignore")
  save(sit4, file = paste0("./docs/CV_Res/gender/Roman_final_subsets/setting4/", DF, ".RData"))
}
'