" Script for the complete cases approach:
  We have a training data set with blockwise missingness and testset in all ranges
  ['fully' testdata completly obsered with all blocks, 
   'miss1_A' testdata that misses block 'A', the rest is observed,
   'miss1_B' testdata that misses block 'B', the rest is observed, ...].
   
  As we can only fit a RF, when the data doesn't contain any missing data points,
  we only use the observations that are complete cases regarding the testset
  [only use obs. with the features in test observed to fit a model].
  
  If there are no complete cases [regarding the testset] in train, we can not 
  apply this approach!
  Aswell, we will only grow a Forest, when we have more than 5 complete cases!
"

# SetWD, load packages and define Functions
setwd("C:/Users/kuche_000/Desktop/MS-Thesis/")
library(checkmate)
library(caret)
library(randomForestSRC)

load_CV_data          <- function(path) {
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
mcc_metric            <- function(conf_matrix) {
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


path <- "data/processed/RH_subsetted_12345/missingness_1234/BLCA_1.RData"

do_CV_missforrest_5   <- function(path = "data/processed/RH_subsetted_12345/missingness_1234/BLCA_1.RData",
                                  num_trees = 100, min_node_size = 10, mtry = NULL,
                                  min_obs = 5) {
  "Evalute the Approach where we only use complete cases for the model fitting! 
   or the cases with 5 blocks! [Scenario1, 2 or 3]
   
   'path' must lead to a list with 2 entrances: 'data' & 'block_names'
   - 'data' is a list filled with 'k' test-train-splits
      --> k-fold-Validation on this test-train-splits!
   - 'block_names' is a list filled with the names of the single blocks 
      & must be ['A', 'B', 'C', 'D', 'clin_block']!
      (Attention: With Scenario2 the order is different, but this is wanted!)
      
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
  assert_string(path)
  if (!grepl("1.RData", path) & !grepl("2.RData", path) & !grepl("3.RData", path)) {
    stop("'path' must end in '1.RData' | '2.RData' | '3.RData'")
  }
  
  # 0-2 'n_tree_impute' & 'maxiter_impute' must be numeric!
  assert_int(n_tree_impute)
  assert_int(maxiter_impute)
  
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
    
    miss1_C[[i]] <- do_evaluation_imputed(train = train_imputed, 
                                          test = test[,-which(colnames(test) %in% curr_data$block_names$C)], 
                                          num_trees = num_trees, 
                                          min_node_size = min_node_size, 
                                          mtry = mtry)
    
    miss1_D[[i]] <- do_evaluation_imputed(train = train_imputed, 
                                          test = test[,-which(colnames(test) %in% curr_data$block_names$D)], 
                                          num_trees = num_trees, 
                                          min_node_size = min_node_size, 
                                          mtry = mtry)
    
    # 2-4-3 TestSet with 2 missing blocks!
    print("Evaluation TestSet w/ 2 missing omics blocks-----------------------")
    miss2_CD[[i]] <- do_evaluation_imputed(train = train_imputed, 
                                           test = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                     curr_data$block_names$D))], 
                                           num_trees = num_trees, 
                                           min_node_size = min_node_size, 
                                           mtry = mtry)
    
    miss2_BD[[i]] <- do_evaluation_imputed(train = train_imputed, 
                                           test = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                     curr_data$block_names$D))], 
                                           num_trees = num_trees, 
                                           min_node_size = min_node_size, 
                                           mtry = mtry)
    
    miss2_BC[[i]] <- do_evaluation_imputed(train = train_imputed, 
                                           test = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                     curr_data$block_names$B))], 
                                           num_trees = num_trees, 
                                           min_node_size = min_node_size, 
                                           mtry = mtry)
    
    miss2_BC[[i]] <- do_evaluation_imputed(train = train_imputed, 
                                           test = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                     curr_data$block_names$B))], 
                                           num_trees = num_trees, 
                                           min_node_size = min_node_size, 
                                           mtry = mtry)
    
    miss2_AC[[i]] <- do_evaluation_imputed(train = train_imputed, 
                                           test = test[,-which(colnames(test) %in% c(curr_data$block_names$C,
                                                                                     curr_data$block_names$A))], 
                                           num_trees = num_trees, 
                                           min_node_size = min_node_size, 
                                           mtry = mtry)
    
    miss2_BC[[i]] <- do_evaluation_imputed(train = train_imputed, 
                                           test = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                     curr_data$block_names$C))], 
                                           num_trees = num_trees, 
                                           min_node_size = min_node_size, 
                                           mtry = mtry)
    
    # 2-4-4 Testset with 3 missing blocks!
    print("Evaluation TestSet w/ 3 missing omics blocks-----------------------")
    miss3_ABC[[i]] <- do_evaluation_imputed(train = train_imputed, 
                                            test = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                      curr_data$block_names$C,
                                                                                      curr_data$block_names$A))], 
                                            num_trees = num_trees, 
                                            min_node_size = min_node_size, 
                                            mtry = mtry)
    
    miss3_ACD[[i]] <- do_evaluation_imputed(train = train_imputed, 
                                            test = test[,-which(colnames(test) %in% c(curr_data$block_names$D,
                                                                                      curr_data$block_names$C,
                                                                                      curr_data$block_names$A))], 
                                            num_trees = num_trees, 
                                            min_node_size = min_node_size, 
                                            mtry = mtry)
    
    miss3_ABD[[i]] <- do_evaluation_imputed(train = train_imputed, 
                                            test = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                      curr_data$block_names$D,
                                                                                      curr_data$block_names$A))], 
                                            num_trees = num_trees, 
                                            min_node_size = min_node_size, 
                                            mtry = mtry)
    
    miss3_BCD[[i]] <- do_evaluation_imputed(train = train_imputed, 
                                            test = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                      curr_data$block_names$C,
                                                                                      curr_data$block_names$D))], 
                                            num_trees = num_trees, 
                                            min_node_size = min_node_size, 
                                            mtry = mtry)
    
    # 2-4-5 Evaluation on single Block Testdata
    print("Evaluation TestSet w/ only 1 observed Block -----------------------")
    single_A[[i]] <- do_evaluation_imputed(train = train_imputed, 
                                           test = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                     curr_data$block_names$C,
                                                                                     curr_data$block_names$D,
                                                                                     curr_data$block_names$clin_block))], 
                                           num_trees = num_trees, 
                                           min_node_size = min_node_size, 
                                           mtry = mtry)
    
    single_B[[i]] <-  do_evaluation_imputed(train = train_imputed, 
                                            test = test[,-which(colnames(test) %in% c(curr_data$block_names$A,
                                                                                      curr_data$block_names$C,
                                                                                      curr_data$block_names$D,
                                                                                      curr_data$block_names$clin_block))], 
                                            num_trees = num_trees, 
                                            min_node_size = min_node_size, 
                                            mtry = mtry)
    
    single_C[[i]] <-  do_evaluation_imputed(train = train_imputed, 
                                            test = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                      curr_data$block_names$A,
                                                                                      curr_data$block_names$D,
                                                                                      curr_data$block_names$clin_block))], 
                                            num_trees = num_trees, 
                                            min_node_size = min_node_size, 
                                            mtry = mtry)
    
    single_D[[i]] <- do_evaluation_imputed(train = train_imputed, 
                                           test = test[,-which(colnames(test) %in% c(curr_data$block_names$B,
                                                                                     curr_data$block_names$C,
                                                                                     curr_data$block_names$A,
                                                                                     curr_data$block_names$clin_block))], 
                                           num_trees = num_trees, 
                                           min_node_size = min_node_size, 
                                           mtry = mtry)
    
    single_CL[[i]] <- do_evaluation_imputed(train = train_imputed, 
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
                   "n_tree_impute" = n_tree_impute, 
                   "maxiter_impute" = maxiter_impute,
                   "time_for_CV"   = time_for_CV)
  
  # 3-3 Return both lists!
  return(list("res_all"  = res_all, 
              "settings" = settings))
}