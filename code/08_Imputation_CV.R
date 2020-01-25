"Script to tackle the blockwise missingness with the any Imputation Approach!
  - missForest!
  - MICE!
  - mdd-sPLS
"
library(missForest)
library(checkmate)

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

# MissForest Approach                                                       ----
# [1] Try Imputation for TrainSetting1                                      ----
# 1-1 Load the data and keep data as 'testdata'
data         <- load_data_extract_block_names()
obs_per_fold <- get_obs_per_fold(data = data$data)
testdata     <- data$data

# Induce NAs to one random block - count NAs in each row!
# Shuffle IDs of the DF
set.seed(1234)
fold_ids <- sample(nrow(testdata), nrow(testdata), replace = FALSE)

# Draw which rows belong to which fold!
set.seed(1234)
observed_blocks <- sample(c(rep("Clin, A", obs_per_fold$amount_train_fold), 
                            rep("Clin, B", obs_per_fold$amount_train_fold),
                            rep("Clin, C", obs_per_fold$amount_train_fold), 
                            rep("Clin, D", obs_per_fold$amount_train_fold)),
                          obs_per_fold$amount_train, replace = FALSE)

# Induce the blockwise missingness!
testdata[which(observed_blocks == "Clin, A"),  # CNV + CLIN ONLY 
         c(data$block_names$mirna_block, data$block_names$mutation_block, 
           data$block_names$rna_block)] <- NA
testdata[which(observed_blocks == "Clin, B"),  # RNA + CLIN 
         c(data$block_names$mirna_block, data$block_names$mutation_block, 
           data$block_names$cnv_block)] <- NA
testdata[which(observed_blocks == "Clin, C"),  # MUTATION + CLIN
         c(data$block_names$mirna_block, data$block_names$rna_block, 
           data$block_names$cnv_block)] <- NA
testdata[which(observed_blocks == "Clin, D"),  # MINRA # CLIN
         c(data$block_names$rna_block, data$block_names$mutation_block, 
           data$block_names$cnv_block)] <- NA


after <- apply(testdata, MARGIN = 1, function(x) sum(is.na(x)))
print("Inducing 'NAs' done! Summary to NAs per row!")
print(summary(after))

# Impute the missing values
testdata_imp <- missForest(testdata, verbose = TRUE)

# Shouldn't contain any missing values!
summary(apply(testdata_imp$ximp, MARGIN = 1, function(x) sum(is.na(x))))
