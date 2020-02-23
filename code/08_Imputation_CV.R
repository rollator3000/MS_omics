"Script to tackle the blockwise missingness with the any Imputation Approach!
  - missForest!
  - MICE!
  - mdd-sPLS
"
# SetWD and define/load functions
setwd("C:/Users/kuche_000/Desktop/MS-Thesis/")
library(missForest)
library(checkmate)
library(randomForestSRC)
library(tidyverse)
library(caret)
library(doParallel)
registerDoParallel(cores = 2)

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
impute_missing_values <- function(data, ntree = 100, maxiter = 10) {
  "Impute missin values in a dataframe with the missForest [RF for imputation]!
   For this we look into 'data' and get the values we need to impute! 
   Then we apply the missForest approach to impute these!
    (For this we need to recode numeric columns to factors,
     if the numeric column has less than 5 unique numeric values. 
     But these are converted back to numeric before returning!)
   
  Args:
    - data (dataframe) : Dataframe with missing values in any column!
    - ntree (int)      : Amount of trees we use for the imputation
    - maxiter (int)    : Maximum number of iterations to be performed, as long 
                         as the stopping criterion is not fullfilled yet!
                         
  Return:
    - data (dataframe) : Dataframe with no more missing values! All values
                         that were missing in 'data' were imputed!
                         Colnames, Coltypes, Dimension etc. stays as it was!
  "
  # [0] Check Inputs  ----------------------------------------------------------
  assert_data_frame(data, min.rows = 1, min.cols = 2)
  assert_int(ntree, lower = 10)
  assert_int(maxiter, lower = 1)
  
  # [1] Prepare the data for the missForest imputation  ------------------------
  # 1-1 Recode to factor if a missing numeric col has less than 5 unique values!
  # 1-1-1 Count how many unique values each column has:
  unique_values_cols <- sapply(1:ncol(data), function(x) length(unique(data[,x])))
  
  # 1-1-2 Get the cols with less than 5 unique values!
  cols_to_recode <- which(unique_values_cols < 5)
  
  # 1-1-3 Remove '1' if inside, as the response is always observed for all obs.!
  if (1 %in% cols_to_recode) {
    cols_to_recode <- cols_to_recode[-which(cols_to_recode == 1)]
  }
  
  # 1-1-4 Convert numeric to factor - if there even is a col to recode!
  if (length(cols_to_recode) > 1) {
    data[cols_to_recode] <- lapply(data[cols_to_recode], as.factor)
  }
  
  # [2] Impute the data with the missForet Approach  ---------------------------
  # 2-1 Impute the data
  data_imputed <- missForest(data, verbose = TRUE, parallelize = 'forests',
                             maxiter = maxiter, ntree = ntree)
  data_imputed <- data_imputed$ximp
  
  # 2-2 Recode the variables back to numeric!
  if (length(cols_to_recode) > 1) {
    data_imputed[cols_to_recode] <- lapply(data_imputed[cols_to_recode], as.numeric)
  }
  
  # [3] Return the imputed Dataframe  ------------------------------------------
  return(data_imputed)
}

path = "data/external/Dr_Hornung/subsetted_12345/missingness_1312/COAD_1.RData"
do_CV_missforrest     <- function(path = "data/external/Dr_Hornung/subsetted_12345/missingness_1312/COAD_1.RData",
                                  num_trees, min_node_size, mtry,
                                  n_tree_impute = 100, maxiter_impute = 10) {
  "Evalute the Approach where we impute the block-wise missing values with the 
   missForest Approach!
   
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
   
   Return:
  "
  # 2-2 Create empty lists to store results in!
  # 2-2-1 Full TestSet
  full <- list()
  # 2-2-2 TestSet with 1 missing omics-block
  miss1_A <- list(); miss1_B <- list(); miss1_C <- list(); miss1_D <- list()
  # 2-2-3 TestSet with 2 missing omics-blocks
  miss2_CD <- list(); miss2_BD <- list(); miss2_BC <- list(); miss2_AD <- list(); miss2_AC <- list(); miss2_AB <- list()
  # 2-2-4 TestSet with 3 missing omics-blocks
  miss3_ABC <- list(); miss3_ABD <- list(); miss3_ACD <- list(); miss3_BCD <- list()
  # 2-2-5 Single BlockTestSet
  single_A <- list(); single_B <- list(); single_CL <- list(); single_C <- list(); single_D <- list()
  
  # 2-3 Start a timer, so we can calc. how long the CV took in total!
  start_time <- Sys.time()
  
  # [3] Start the CV [5-fold per default!] -------------------------------------
  for (i in 0:4) {
    # Impute the missing values
    train_df_imputed <- missForest(train_df, verbose = TRUE)
    
    # Shouldn't contain any missing values!
    summary(apply(train_df_imputed$ximp, MARGIN = 1, function(x) sum(is.na(x))))
  }
}