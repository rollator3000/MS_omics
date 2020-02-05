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

load_data_extract_block_names <- function(path = "data/external/Dr_Hornung/Data/ProcessedData_subsets/seed_1234/KIRC_Subset.RData",
                                          response = 'gender') {
  "Function to load (the already subsetted) data & returning 1 big DataFrame, 
   of all single blocks incl. the colnames to each block!
   
   Args:
    - path (char)     : path to a DF w/ block wise structure! 
                        Shall contain 'rna_subset', 'cnv_subset', 'mirna_subset',
                        'clin' & 'mutation_subset' block!
    - response (char) : feature used as reponse class - must be in 'clin' block!
                        + MUST be binary - else it will throw an error!
   
   Return:
    list with 2 entrances:  
      1 - data: all blocks pasted together as single DF, the first col of it 
                equals the response!
      2 - block_names: names of the features in the different  blocks
  "
  # [0] Check Inputs -----------------------------------------------------------
  # 0-1 Load data from path & check whether it has all blocks
  #     if the path is not valid, load() will throw an error!
  load(path)
  if (any(!exists('clin_') & !exists('cnv_sub') & !exists('mirna_sub') & 
          !exists('mutation_sub') & !exists('rna_sub'))) {
    stop("'path' led to a DF with at least one missing block!")
  }
  
  # 0-2 Check that the response is in the clincial block!
  if (!(response %in% colnames(clin_))) stop("Clin Block has no 'response' feature")
  
  # 0-3 Check that response is binary!
  if (length(levels(as.factor(clin_[response][,1]))) != 2) {
    stop("The selected Response doesn't have 2 levels! 
         Please choose a binary response from the 'clin'-block!")
  }
  
  # [1] Create single DF -------------------------------------------------------
  # 1-1 Extract the response from the clinical block & rm tehe col from the df
  response_       <- clin_[response]
  clin_[response] <- NULL
  
  # 1-2 Bind the single blocks to a big DF!
  df <- cbind(response_, clin_, cnv_sub, rna_sub, mirna_sub, mutation_sub)
  
  # 1-3 Recode the response as factor!
  df[,colnames(df) == response] <- factor(df[,colnames(df) == response])
  
  # 1-4 Extract the colnames of the single blocks & save them in list:
  block_variables <- list("clin_block"     = colnames(clin_),
                          "cnv_block"      = colnames(cnv_sub),
                          "rna_block"      = colnames(rna_sub),
                          "mutation_block" = colnames(mutation_sub),
                          "mirna_block"    = colnames(mirna_sub))
  
  # [2] Return list with the df & the colnames of the single blocks ------------
  return(list("data" = df,
              "block_names" = block_variables))
}
get_obs_per_fold              <- function(data) {
  "Find the amount of observations needed in each test fold, so the 4 different 
   training folds are equally sized!
   --> This might lead to smaller testfolds than trainingfolds!
       
  Args:
    - data (data.frame) : dataframe on which we want to do CV on blockwise
                          missingness patterns!
  Return:
    - list filled with:
      - 'amount_train':      amount of Observations used for Training in total!
      - 'amount_train_fold': amount of Observations in each Trainfold!
      - 'amount_test':       amount of Observations in each Testfold
  "
  # [0] Check Inputs -----------------------------------------------------------
  assert_data_frame(data, min.rows = 1)
  
  # [1] Calculat the amount of Obs. in Test & Train Fold -----------------------
  #     Get the amount of Obs. in each Trainings- & Testfold, so the 
  #     Trainigfolds have the same size!
  #     Amount of trainig observations must be dividable by 4, so all folds
  #     are equally sized! If this is not true increase the amount of obs. in 
  #     trainingfolds by 1, as long as it is dividable by 4
  #     The remaining obs. are used as TestSet
  # 1-1 Split to Train & Test
  amount_train <- floor(4/5 * nrow(data)) 
  
  # 1-2 Count up 'amount_train' until it is dividable by 4!
  while ((amount_train %% 4) != 0) {
    amount_train <- amount_train + 1
  }
  
  # 1-3 Amount of Obs. in Test + in sach train fold!
  amount_train_fold <- amount_train / 4
  amount_test       <- nrow(data) - amount_train
  
  # 1-4 Print info
  writeLines(paste0("From ", nrow(data), " Observations, each TestFold will hold "
                    , amount_test,  " Test-Observations!\n", amount_train, 
                    " Observations will be used for training --> ", 
                    amount_train_fold, " Observations per Trainingsfold"))
  
  # [2] Return Train- / Test-fold sizes ----------------------------------------
  return(list("amount_train"      = amount_train,
              "amount_train_fold" = amount_train_fold,
              "amount_test"       = amount_test))
}
data_path = "data/external/Dr_Hornung/Data/ProcessedData_subsets/seed_1234/LGG_Subset.RData"
response = "gender"
seed = 1312
do_CV_missforrest             <- function(data_path = "data/external/Dr_Hornung/Data/ProcessedData_subsets/seed_1234/LGG_Subset.RData",
                                          response = "gender", seed = 1312,
                                          num_trees, min_node_size, mtry) {
  "
  Function to evaluate the missForest Approach in the settings of blockwise
  missingness in SCENARIO1!
  
  Data is split into test and train set [curently fixed to 5-fold], with the 
  little adjustment, the amount of traindata can be split into 4 folds w/o rest
    --> All train folds have same amount of observations!
    --> TestFold can be a bit smaller - not necessarily!
    
  Then each [equally sized] trainingsfold is censored to scenario 1, so that 
  each fold has an observed clinical block & exactly one observed omics block!
  Then we impute the missing values with the missForest Approach!
  As a result we'll have a data w/o any missing values - on this we can apply
  a regular RF! This RF will be evaluated w/ Accuracy, Precision, Specifity, 
  F1-Socre,....
  
  Args:
    - data_path (char)    : Path to the data, we want to CV! This should lead 
                            to a file w/ multiple sub DFs 
                            [details in 'create_data()']
    - response (chr)      : The repsonse we want to model - 
                            MUST be in the 'clin'-block & MUST be binary!
    - seed (int)          : Needed for reproducibility! Shuffle rows of DF & 
                            assign which blocks were observed for the obs.!
    - num_trees (int)     : Amount of trees we use to fit a RF on the imputed data!
  
  Return:
    - - list filled w/:
        * 'res_all' [the CV Results on different testsets]:
            - Results of the CV with the different 
                         
        * 'settings' [settings used to do the CV - all arguments!]
            - datapath, seed, response, mtry,time, .... 
  "
  # [0] Check Inputs -----------------------------------------------------------
  # 0-0 data_path, response are all checked within 'load_data_extract_block_names()'
  
  # 0-1 mtry, min_node_size & num_trees are all checked within 'rfsrc()'
  
  # 0-2 seed must be an integer
  assert_int(seed)
  
  # [1] Get the data & dimensions of train & test folds! -----------------------
  # 1-1 Load Data & the names of the single blocks in the data!
  data <- load_data_extract_block_names(path = data_path, response = response)
  
  # 1-2 Get Obs. per fold [Train & Test]
  obs_per_fold <- get_obs_per_fold(data = data$data)
  
  # [2] Shuffle the data & create lists to save the results in -----------------
  # 2-1 Shuffle IDs of data, we use to split data later on!
  set.seed(seed)
  fold_ids <- sample(nrow(data$data), nrow(data$data), replace = FALSE)
  
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
    
    print(paste0("FOLD: ", as.character(i + 1), "/5 -------------------------"))
    
    # 3-1 Get TestSet from 'data', by taking the first 'obs_per_fold$amount_test' 
    #     IDs in 'fold_ids'
    test_ids <- fold_ids[((i * obs_per_fold$amount_test) + 1):(((i + 1) * obs_per_fold$amount_test))]
    test_df  <- data$data[test_ids,]
    
    # 3-2 Get the TrainSet from 'data' [= IDs not in TestSet] 
    train_ids <- fold_ids[-which(fold_ids %in% test_ids)]
    train_df  <- data$data[train_ids,]
    
    # 3-3 Induce blockwise missingness [SCENARIO_1]
    # 3-3-1 Sample equally sized 'observed' blocks [SCENARIO_1]
    set.seed(seed + i)
    observed_blocks <- sample(c(rep("Clin, A", obs_per_fold$amount_train_fold), 
                                rep("Clin, B", obs_per_fold$amount_train_fold),
                                rep("Clin, C", obs_per_fold$amount_train_fold), 
                                rep("Clin, D", obs_per_fold$amount_train_fold)),
                              obs_per_fold$amount_train, replace = FALSE)
    
    # 3-4 Induce the blockwise missingness into the train DF
    train_df[which(observed_blocks == "Clin, A"),  # CNV + CLIN ONLY 
             c(data$block_names$mirna_block, data$block_names$mutation_block, 
               data$block_names$rna_block)] <- NA
    train_df[which(observed_blocks == "Clin, B"),  # RNA + CLIN 
             c(data$block_names$mirna_block, data$block_names$mutation_block, 
               data$block_names$cnv_block)] <- NA
    train_df[which(observed_blocks == "Clin, C"),  # MUTATION + CLIN
             c(data$block_names$mirna_block, data$block_names$rna_block, 
               data$block_names$cnv_block)] <- NA
    train_df[which(observed_blocks == "Clin, D"),  # MINRA # CLIN
             c(data$block_names$rna_block, data$block_names$mutation_block, 
               data$block_names$cnv_block)] <- NA
    
    after <- apply(train_df, MARGIN = 1, function(x) sum(is.na(x)))
    print("Inducing 'NAs' done! Summary to NAs per row!")
    print(summary(after))
    
    diff_values_per_var <- sapply(1:ncol(train_df), function(x) length(unique(train_df[,x])))
    table(diff_values_per_var) # --> Values with a less than 5 unique values should 
                               #     coded as factor before imputation!
    
    # Impute the missing values
    train_df_imputed <- missForest(train_df, verbose = TRUE)
    
    # Shouldn't contain any missing values!
    summary(apply(train_df_imputed$ximp, MARGIN = 1, function(x) sum(is.na(x))))
  }
}