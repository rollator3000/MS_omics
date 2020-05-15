" Script for the complete cases approach in the setting of blockwise missingness!

  We have a training data set with blockwise missingness and a testset from 
  fully observed, to a testset consiting only one block!
  e.g. - 'fully': testdata completly obsered with all blocks 
       - 'miss1_A' testdata that misses block 'A', the rest is observed
       - 'miss2_BC' testdata that misses block 'B' & 'C' , the rest is observed
   
  As we can only fit a RF, when the data doesn't contain any missing data points,
  we only use the observations that are complete cases regarding the testset
  [only use obs. with the features in test observed to fit a model].
  If the train data doesn't contain any complete cases [regarding the testset]  
  then the complete cases approach can not be applied!
"
# SetWD, load packages, define Functions & set cores for parallel computing!
setwd("C:/Users/kuche_000/Desktop/MS-Thesis/")
library(checkmate)
library(caret)
library(randomForestSRC)
library(doParallel)

detectCores()
registerDoParallel(cores = 2)

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
do_evaluation_CC    <- function(train, test, num_trees, min_node_size, mtry, 
                                min_obs) {
  "Evaluate the complete cases approach for a given 'train' and 'test' set!
   Only the observations from 'train', that are completly observed in the 
   variables, that appear in 'test', are used for training a RF.
   If train consits of at least 'min_obs'-complete-observations, a RF can be
   trained on it! 
   The Evaluation on the testset results in Acc, F1, selected Vars, MMC, ....
   
   Args:
    - train (data.frame) : Dataframe that has blockwise missingness.
    - test  (data.frame) : Dataframe we want to do predicitons for!
    - num_trees (int)    : Amount of trees to be fit on the complete cases in 
                           'train' [refferting to the testset]!
                              > If NULL: 1000 is the default!
    - min_node_size(int) : Amount of Observations a node must at least contain,
                           so the model keeps on trying to split them!  
                              > If NULL: 1 is the default!
    - mtry (int)         : Amount of split-variables we try, when looking for
                           a split variable! 
                              > If 'NULL': mtry = sqrt(p)
    - min_obs (int)      : Amount of complete cases we need at least to fit a RF
                           Must be >= 'min_node_size', else the tree can not grow!
   Return:
    - list w/ metrics [accuracy, f1, mcc, roc, selected vars, ...] 
       > Metrics w/ NA are replaced by their worst possible value!
  "
  # [0] Check Inputs  ----------------------------------------------------------
  # 0-1 'train' & 'test' must be dataframes
  assert_data_frame(train)
  assert_data_frame(test)
  
  # 0-2 'num_trees', 'min_node_size' & 'mtry' must be integer > 0 if NOT NULL
  if (!is.null(num_trees)) assert_int(num_trees, lower = 10)
  if (!is.null(mtry)) assert_int(mtry, lower = 5)
  if (!is.null(min_node_size)) assert_int(min_node_size, lower = 1)
  
  # 0-3 'min_obs' must be integer > 'min_node_size'
  if (is.null(min_node_size)) assert_int(min_obs, lower = 5)
  if (!is.null(min_node_size)) assert_int(min_obs, lower = min_node_size)
  
  # [1] Subset the traindata so it only consits of complete cases!  ------------
  # 1-1 Extract the covariates from 'test'
  covariates_test <- colnames(test)
  
  # 1-2 Get the Observations from Train that do not miss any in 'covariates_test'
  miss_values <- foreach(x = 1:nrow(train)) %dopar% sum(is.na(train[x, covariates_test]))
  usable_obs  <- which(unlist(miss_values) == 0)
  
  # 1-3 Check whether we have found enough complete cases to fit a RF, if not 
  #     stop here and return the reason for abortion!
  if (length(usable_obs) < min_node_size) return("No Complete Cases!")
  complete_train <- train[usable_obs, covariates_test]
  
  # [2] Fit a RF - if enough complete cases were found in [1]  -----------------
  # 2-1 Fit a RF on the complete cases
  # 2-1-1 Extract Response and create the formula the RF to be fit on!
  response    <- colnames(train)[1]
  formula_all <- as.formula(paste(response, " ~ ."))
  
  # 2-1-2 Fit the actual RF
  RF <- rfsrc(formula = formula_all, data = complete_train, 
              ntree = num_trees, mtry = mtry, nodesize = min_node_size, 
              samptype = "swr", seed = 12345678, var.used = 'all.trees')
  
  # 2-2 Get Prediciton on the testset from the RF
  predicitons <- predict(RF, test)
  
  # [3] Get Metrics of predicitive power  --------------------------------------
  # 3-1 ConfusionMatrix
  confmat <- caret::confusionMatrix(predicitons$class, test[,1])
  
  # 3-2 Are under the ROC Curve
  roc1 <- tryCatch(pROC::auc(pROC::roc(as.numeric(test[,1]), 
                              as.numeric(predicitons$class),
                              levels = unique(as.numeric(test[,1])),
                              direction = "<")),
                   error = function(e) "not defined!")
  if (is.numeric(roc1)) roc1 <- as.numeric(roc1)
  
  roc2 <- tryCatch(pROC::auc(pROC::roc(as.numeric(test[,1]), 
                              as.numeric(predicitons$class),
                              levels = unique(as.numeric(test[,1])),
                              direction = ">")),
                   error = function(e) "not defined!")
  if (is.numeric(roc2)) roc2 <- as.numeric(roc2)
  
  # 3-3 MCC Matthews correlation coefficient [only for binary cases!]
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
              "AUC1"        = roc1,
              "AUC2"        = roc2,
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
do_CV_CC_5 <- function(path = "data/processed/TCGA_subset_12345/missingness_1234/BLCA_1.RData",
                       num_trees = 100, min_node_size = 5, mtry = NULL,
                       min_obs = 5) {
  "Evalute the Approach that uses only complete cases for the model fitting! 
   For the cases with 5 blocks! [Scenario1, 2 or 3]
   
   'path' must lead to a list with 2 entrances: 'data' & 'block_names'
    - 'data' is a list filled with 'k' test-train-splits
       --> k-fold-Validation on this test-train-splits!
    - 'block_names' is a list filled with the names of the single blocks 
       & must be ['A', 'B', 'C', 'D', 'clin_block']!
       (Attention: With Scenario2 the order is different, but this is wanted!)
      
   Based on the 'k' test-train-splits in 'data', we will evaluate the complete 
   cases approach. For this we remove all observations fromt the traindata that
   miss at least one covariate of the features in the testset!
   Then for each testsituation (fully observed testset,.., single block testset) 
   we remove unusable observations, and use the remaining ones to fit a RF.
   If the trainset has no observations, we can not train a RF with it & therefore
   this approach is not possible in this setting!
   If the trainset is not empty a RF is fitted and the performance is measured
   with the testset & rated with Accuracy, Precision, Specifity, F1-Socre,...
   
   Args:
    - path (str)         : path to the data w/ blockwise missingness for the CV.
                           Must end in '1.RData', '2.RData' or '3.RData'
                           --> List with 2 entrances: 'data' & 'block_names'
                                - 'data' consitis of 'k' test-train-splits, 
                                   where train has missingness induced and the 
                                   test-set is fully observed!
                                - 'block_names' contains all colnames of the 
                                   different blocks!
    - num_trees (int)    : Amount of trees to be fit on the complete cases
                             > If NULL: 1000 is the default!
    - min_node_size(int) : Amount of Observations a node must at least contain,
                           so the model keeps on trying to split them!  
                             > If NULL: 1 is the default!
    - mtry (int)         : Amount of split-variables we try, when looking for
                           a split variable! 
                             > If 'NULL': mtry = sqrt(p)
    - min_obs (int)      : Amout of observations the trainset has to contain at
                           least, so we fit a model on it!
                           Must be >= 'min_node_size', else a tree can not start
                           to grow!
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
  # 0-1 path must be string ending in 1.- 2.- or 3.RData! 
  assert_string(path)
  if (!grepl("1.RData", path) & !grepl("2.RData", path) & !grepl("3.RData", path)) {
    stop("'path' must end in '1.RData' | '2.RData' | '3.RData'")
  }
  
  # 0-2 'num_trees', 'min_node_size' & 'mtry' must be integer > 0 if NOT NULL
  if (!is.null(num_trees)) assert_int(num_trees, lower = 10)
  if (!is.null(mtry)) assert_int(mtry)
  if (!is.null(min_node_size)) assert_int(min_node_size, lower = 1)
  
  # 0-3 'min_obs' must be int > 'min_node_size' 
  if (is.null(min_node_size)) assert_int(min_obs, lower = 5)
  if (!is.null(min_node_size)) assert_int(min_obs, lower = min_node_size)
  
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
    
    # 2-1 Print the current Fold
    print(paste0("FOLD ", i, "/", k_splits, " -------------------------------"))
    
    # 2-2 Extract the test and train set from 'curr_data'
    train <- curr_data$data[[i]]$train
    test  <- curr_data$data[[i]]$test
    
    # 2-3 Evaluate the RF on the different Testsets! From a full TestSet w/o any 
    #     missing blocks to TestSets w/ only one observed Block!
    #     For this we look into the trainset, whether we have complete cases in
    #     reference to the testset - if yes, only use these observations to fit
    #                                        a model, which is evaluated on test 
    #                                        then!
    #                              - if no, a model can not be fitted with the
    #                                       complete cases approach & it will 
    #                                       return a message instead of metrics!
    # 2-3-1 Full TestSet!
    print("Evaluation full TestSet -------------------------------------------")
    full[[i]] <- do_evaluation_CC(train = train, test = test,
                                  num_trees = num_trees, 
                                  min_node_size = min_node_size, 
                                  mtry = mtry, min_obs = min_obs)  
    
    # 2-3-2 TestSet with 1 missing block!
    print("Evaluation TestSet w/ 1 missing omics block------------------------")
    miss1_A[[i]] <- do_evaluation_CC(train = train, 
                                     test = test[,-which(colnames(test) %in% curr_data$block_names$A)],
                                     num_trees = num_trees, min_node_size = min_node_size, 
                                     mtry = mtry, min_obs = min_obs)  
    
    miss1_B[[i]] <- do_evaluation_CC(train = train, 
                                     test = test[,-which(colnames(test) %in% curr_data$block_names$B)],
                                     num_trees = num_trees, min_node_size = min_node_size, 
                                     mtry = mtry, min_obs = min_obs)  
    
    miss1_C[[i]] <- do_evaluation_CC(train = train, 
                                     test = test[,-which(colnames(test) %in% curr_data$block_names$C)],
                                     num_trees = num_trees, min_node_size = min_node_size, 
                                     mtry = mtry, min_obs = min_obs)  
    
    miss1_D[[i]] <- do_evaluation_CC(train = train, 
                                     test = test[,-which(colnames(test) %in% curr_data$block_names$D)],
                                     num_trees = num_trees, min_node_size = min_node_size, 
                                     mtry = mtry, min_obs = min_obs)  
    
    # 2-3-3 TestSet with 2 missing blocks!
    print("Evaluation TestSet w/ 2 missing omics blocks-----------------------")
    miss2_CD[[i]] <- do_evaluation_CC(train = train, 
                                      test = test[,-which(colnames(test) %in% c(curr_data$block_names$D,
                                                                                curr_data$block_names$C))],
                                      num_trees = num_trees, min_node_size = min_node_size, 
                                      mtry = mtry, min_obs = min_obs)  
    
    miss2_BD[[i]] <- do_evaluation_CC(train = train, 
                                      test = test[,-which(colnames(test) %in% c(curr_data$block_names$D,
                                                                                curr_data$block_names$B))],
                                      num_trees = num_trees, min_node_size = min_node_size, 
                                      mtry = mtry, min_obs = min_obs)
    
    miss2_BC[[i]] <- do_evaluation_CC(train = train, 
                                      test = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                curr_data$block_names$C))],
                                      num_trees = num_trees, min_node_size = min_node_size, 
                                      mtry = mtry, min_obs = min_obs)
    
    miss2_AD[[i]] <- do_evaluation_CC(train = train, 
                                      test = test[,-which(colnames(test) %in% c(curr_data$block_names$D,
                                                                                curr_data$block_names$A))],
                                      num_trees = num_trees, min_node_size = min_node_size, 
                                      mtry = mtry, min_obs = min_obs)
    
    miss2_AC[[i]] <- do_evaluation_CC(train = train, 
                                      test = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                curr_data$block_names$C))],
                                      num_trees = num_trees, min_node_size = min_node_size, 
                                      mtry = mtry, min_obs = min_obs)
    
    miss2_AB[[i]] <- do_evaluation_CC(train = train, 
                                      test = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                curr_data$block_names$B))],
                                      num_trees = num_trees, min_node_size = min_node_size, 
                                      mtry = mtry, min_obs = min_obs)
    
    # 2-3-4 Testset with 3 missing blocks!
    print("Evaluation TestSet w/ 3 missing omics blocks-----------------------")
    miss3_ABC[[i]] <- do_evaluation_CC(train = train, 
                                       test = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$B,
                                                                                 curr_data$block_names$C))],
                                       num_trees = num_trees, min_node_size = min_node_size, 
                                       mtry = mtry, min_obs = min_obs)
    
    miss3_ACD[[i]] <- do_evaluation_CC(train = train, 
                                       test = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$D,
                                                                                 curr_data$block_names$C))],
                                       num_trees = num_trees, min_node_size = min_node_size, 
                                       mtry = mtry, min_obs = min_obs)
    
    miss3_ABD[[i]] <- do_evaluation_CC(train = train, 
                                       test = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$B,
                                                                                 curr_data$block_names$D))],
                                       num_trees = num_trees, min_node_size = min_node_size, 
                                       mtry = mtry, min_obs = min_obs)
    
    miss3_BCD[[i]] <- do_evaluation_CC(train = train, 
                                       test = test[,-which(colnames(test) %in% c(curr_data$block_names$D,
                                                                                 curr_data$block_names$B,
                                                                                 curr_data$block_names$C))],
                                       num_trees = num_trees, min_node_size = min_node_size, 
                                       mtry = mtry, min_obs = min_obs)
    
    # 2-3-5 Evaluation on single Block Testdata
    print("Evaluation TestSet w/ only 1 observed Block -----------------------")
    single_A[[i]] <- do_evaluation_CC(train = train, 
                                      test = test[,-which(colnames(test) %in% c(curr_data$block_names$D,
                                                                                curr_data$block_names$B,
                                                                                curr_data$block_names$C,
                                                                                curr_data$block_names$clin_block))],
                                      num_trees = num_trees, min_node_size = min_node_size, 
                                      mtry = mtry, min_obs = min_obs)
    
    single_B[[i]] <-  do_evaluation_CC(train = train, 
                                       test = test[,-which(colnames(test) %in% c(curr_data$block_names$D,
                                                                                 curr_data$block_names$A,
                                                                                 curr_data$block_names$C,
                                                                                 curr_data$block_names$clin_block))],
                                       num_trees = num_trees, min_node_size = min_node_size, 
                                       mtry = mtry, min_obs = min_obs)
    
    single_C[[i]] <-  do_evaluation_CC(train = train, 
                                       test = test[,-which(colnames(test) %in% c(curr_data$block_names$D,
                                                                                 curr_data$block_names$B,
                                                                                 curr_data$block_names$A,
                                                                                 curr_data$block_names$clin_block))],
                                       num_trees = num_trees, min_node_size = min_node_size, 
                                       mtry = mtry, min_obs = min_obs)
    
    single_D[[i]] <- do_evaluation_CC(train = train, 
                                      test = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                curr_data$block_names$B,
                                                                                curr_data$block_names$C,
                                                                                curr_data$block_names$clin_block))],
                                      num_trees = num_trees, min_node_size = min_node_size, 
                                      mtry = mtry, min_obs = min_obs)
    
    single_CL[[i]] <- do_evaluation_CC(train = train, 
                                       test = test[,-which(colnames(test) %in% c(curr_data$block_names$D,
                                                                                 curr_data$block_names$B,
                                                                                 curr_data$block_names$C,
                                                                                 curr_data$block_names$A))],
                                       num_trees = num_trees, min_node_size = min_node_size, 
                                       mtry = mtry, min_obs = min_obs)
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
                   "time_for_CV"   = time_for_CV,
                   "min_obs"       = min_obs)
  
  # 3-3 Return both lists!
  return(list("res_all"  = res_all, 
              "settings" = settings))
}
do_CV_CC_3 <- function(path = "data/processed/TCGA_subset_12345/missingness_1234/BLCA_4.RData",
                       num_trees = 100, min_node_size = 5, mtry = NULL,
                       min_obs = 5) {
  "Evalute the Approach that only uses complete cases for the model fitting! 
   For the cases with 3 blocks! [Scenario4]
   
   'path' must lead to a list with 2 entrances: 'data' & 'block_names'
    - 'data' is a list filled with 'k' test-train-splits
       --> k-fold-Validation on this test-train-splits!
    - 'block_names' is a list filled with the names of the single blocks 
       & must be ['A', 'B' & 'clin_block']!
       
   Based on the 'k' test-train-splits in 'data', we will evaluate the complete 
   cases approach. For this we remove all observations fromt the traindata that
   miss at least one covariate of the features in the testset!
   Then for each testsituation (fully observed testset,.., single block testset) 
   we remove unusable observations, and use the remaining ones to fit a RF.
   If the trainset has no observations, we can not train a RF with it & therefore
   this approach is not possible in this setting!
   If the trainset is not empty a RF is fitted and the performance is measured
   with the testset & rated with Accuracy, Precision, Specifity, F1-Socre,...
   
   Args:
    - path (str)         : path to the data w/ blockwise missingness for the CV.
                           Must end in '4.RData'
                           --> List with 2 entrances: 'data' & 'block_names'
                                - 'data' consitis of 'k' test-train-splits, 
                                   where train has missingness induced and the 
                                   test-set is fully observed!
                                - 'block_names' contains all colnames of the 
                                   different blocks!
    - num_trees (int)    : Amount of trees to be fit on the complete cases
                             > If NULL: 1000 is the default!
    - min_node_size(int) : Amount of Observations a node must at least contain,
                           so the model keeps on trying to split them!  
                             > If NULL: 1 is the default!
    - mtry (int)         : Amount of split-variables we try, when looking for
                           a split variable! 
                             > If 'NULL': mtry = sqrt(p)
    - min_obs (int)      : Amout of observations the trainset has to contain at
                           least, so we fit a model on it!
                           Must be >= 'min_node_size', else a tree can not start
                           to grow!
     Return:
      - list filled w/:
        * 'res_all' [the CV Results on different testsets]:
           >> A = CNV-Block, B = RNA-Block, C = Mutation-Block, D = Mirna-Block
            - full    : CV Results for each fold on the fully observed testdata!
            - miss1_A : CV Results for each fold on the testdata, w/ missing 
                        block 'A' [2 random omics blocks united]
            - miss1_B : CV Results for each fold on the testdata, w/ missing rna 
                        block 'B' [2 random omics blocks united]
                .
                .
            - single_B: CV-Results for each fold on the testdata w/ only block
                        'B' [2 random omics blocks united] as features in test!
                         
        * 'settings' [settings used to do the CV - all arguments!]
            - datapath, seed, response, mtry, time for Evaluation, ...
  "
  # [0] Check Inputs  ----------------------------------------------------------
  # 0-1 path must be string ending in 4.RData! 
  assert_string(path, fixed = "4.RData")
  
  # 0-2 'num_trees', 'min_node_size' & 'mtry' must be integer > 0 if NOT NULL
  if (!is.null(num_trees)) assert_int(num_trees, lower = 10)
  if (!is.null(mtry)) assert_int(mtry)
  if (!is.null(min_node_size)) assert_int(min_node_size, lower = 1)
  
  # 0-3 'min_obs' must be int > 'min_node_size' 
  if (is.null(min_node_size)) assert_int(min_obs, lower = 5)
  if (!is.null(min_node_size)) assert_int(min_obs, lower = min_node_size)
  
  # [1] Get the data and create list to save results  --------------------------
  # 1-1 Load CV-Data [already splitted - data checked in 'load_CV_data' itself]
  curr_data <- load_CV_data(path = path)
  
  # 1-1-1 Must contain 'A', 'B' & 'clin_block' as block_names
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
    
    # 2-3 Evaluate the RF on the different Testsets! Fromfull TestSet w/o any 
    #     missing blocks to TestSets w/ only one observed Block!
    #     For this we look into the trainset, whether we have complete cases 
    #     reffering to the testset - if yes, only use these observations to fit
    #                                        a model, which is evaluated on test 
    #                                        then!
    #                              - if no, a model can not be fitted with the
    #                                       complete cases approach!
    # 2-3-1 Full TestSet!
    print("Evaluation full TestSet -------------------------------------------")
    full[[i]] <- do_evaluation_CC(train = train, test = test,
                                  num_trees = num_trees, 
                                  min_node_size = min_node_size, 
                                  mtry = mtry, min_obs = min_obs)  
    
    # 2-3-2 TestSet with 1 missing block!
    print("Evaluation TestSet w/ 1 missing omics block------------------------")
    miss1_A[[i]] <- do_evaluation_CC(train = train, 
                                     test = test[,-which(colnames(test) %in% curr_data$block_names$A)],
                                     num_trees = num_trees, min_node_size = min_node_size, 
                                     mtry = mtry, min_obs = min_obs)  
    
    miss1_B[[i]] <- do_evaluation_CC(train = train, 
                                     test = test[,-which(colnames(test) %in% curr_data$block_names$B)],
                                     num_trees = num_trees, min_node_size = min_node_size, 
                                     mtry = mtry, min_obs = min_obs)  
    
    # 2-3-3 Evaluation on single Block Testdata
    print("Evaluation TestSet w/ only 1 observed Block -----------------------")
    single_A[[i]] <- do_evaluation_CC(train = train, 
                                      test = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                curr_data$block_names$clin_block))],
                                      num_trees = num_trees, min_node_size = min_node_size, 
                                      mtry = mtry, min_obs = min_obs)
    
    single_B[[i]] <-  do_evaluation_CC(train = train, 
                                       test = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                 curr_data$block_names$clin_block))],
                                       num_trees = num_trees, min_node_size = min_node_size, 
                                       mtry = mtry, min_obs = min_obs)
    
    single_CL[[i]] <- do_evaluation_CC(train = train, 
                                       test = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                 curr_data$block_names$A))],
                                       num_trees = num_trees, min_node_size = min_node_size, 
                                       mtry = mtry, min_obs = min_obs)
  }
  # 2-6 Stop the time and take the difference!
  time_for_CV <- Sys.time() - start_time
  
  # [3] Return the results & settings of parameters used to do CV! -------------
  # 3-1 Collect all CV Results in a list!
  res_all <- list("full" = full,  "miss1_A" = miss1_A, 
                  "miss1_B" = miss1_B, "single_A" = single_A, 
                  "single_B" = single_B, "single_CL" = single_CL)
  
  # 3-2 Collect the Settings, used to do the CV!
  settings <- list("data_path"     = path,
                   "num_folds"     = k_splits,
                   "num_trees"     = num_trees,
                   "mtry"          = mtry, 
                   "min_node_size" = min_node_size,
                   "time_for_CV"   = time_for_CV,
                   "min_obs"       = min_obs)
  
  # 3-3 Return both lists!
  return(list("res_all"  = res_all, 
              "settings" = settings))
}

# RunMain ----------------------------------------------------------------------
"Run the CV for all DFs from the missForest paper ------------------------------"
DFs_w_gender <- c("COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC","LGG", "BLCA",
                  "LUAD", "LUSC", "PAAD", "SARC", "SKCM", "STAD")

# ----- Situation 1
for (DF in DFs_w_gender) {
  
  print(paste0("----- Situation 1 for DF: '", DF, "' -----"))
  
  # Create the path for the current DF
  curr_path <- paste0("data/processed/TCGA_subset_12345/missingness_1234/", DF, "_1.RData")
  
  print("Setting - 1/1")
  sit1 <- do_CV_CC_5(path = curr_path,
                     num_trees = 300, min_node_size = 5, mtry = NULL, min_obs = 5)
  save(sit1, file = paste0("./docs/CV_Res/TCGA/CompleteCase_Approach/setting1/", DF, ".RData"))
}

# ----- Situation 2
for (DF in DFs_w_gender) {
  
  print(paste0("----- Situation 2 for DF: '", DF, "' -----"))
  
  # Create the path for the current DF
  curr_path <- paste0("data/processed/TCGA_subset_12345/missingness_1234/", DF, "_2.RData")
  
  print("Setting - 1/1")
  sit2_1 <- do_CV_CC_5(path = curr_path,
                       num_trees = 300, min_node_size = 5, mtry = NULL, min_obs = 5)
  save(sit2_1, file = paste0("./docs/CV_Res/TCGA/CompleteCase_Approach/setting2/", DF, ".RData"))
  
}

# ----- Situation 3
for (DF in DFs_w_gender) {
  
  print(paste0("----- Situation 3 for DF: '", DF, "' -----"))
  
  # Create the path for the current DF
  curr_path <- paste0("data/processed/TCGA_subset_12345/missingness_1234/", DF, "_3.RData")
  
  print("Setting - 1/1")
  sit3_1 <- do_CV_CC_5(path = curr_path,
                       num_trees = 300, min_node_size = 5, mtry = NULL, min_obs = 5)
  save(sit3_1, file = paste0("./docs/CV_Res/TCGA/CompleteCase_Approach/setting3/", DF, ".RData"))
  
}

# ----- Situation 4
for (DF in DFs_w_gender) {
  
  print(paste0("----- Situation 4 for DF: '", DF, "' -----"))
  
  # Create the path for the current DF
  curr_path <- paste0("data/processed/TCGA_subset_12345/missingness_1234/", DF, "_4.RData")
  
  print("Setting - 1/1")
  sit2_1 <- do_CV_CC_3(path = curr_path,
                       num_trees = 300, min_node_size = 5, mtry = NULL, min_obs = 5)
  save(sit2_1, file = paste0("./docs/CV_Res/TCGA/CompleteCase_Approach/setting4/", DF, ".RData"))
  
}