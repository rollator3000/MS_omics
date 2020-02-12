"Script to split the data into k different Train & TestSets. To each of the 'k' 
 TrainSets we will induce block-wise missingness of a certain scenario - totally 4.
"
# Load Funcitons, librarys & set the WD!
setwd("C:/Users/kuche_000/Desktop/MS-Thesis/")
library(checkmate)
library(caret)

load_data_extract_block_names <- function(path = "data/external/Dr_Hornung/subsetted_12345/LGG_subset.RData",
                                          response = 'gender') {
  "Function to load the - already subsetted - data!
   Return one single DF (where all blocks are joint together) & 
   the colnames of each block!
   
   Args:
    - path (char)     : path to a DF w/ block wise structure! 
                        Shall contain 'rna', 'cnv', 'mirna', 'clin' & 'mutation'
                        blocks!
    - response (char) : feature used as reponse - must be in 'clin' block & 
                        must be binary - have exact 2 factor levels
   
   Return:
    list with 2 entrances:  
      1 - data: dataframe, where all blocks are pasted together, and the first 
                column equals the response and is of class factor!!
      2 - block_names: list with all colnames of each block in data!
  "
  # [0] Check Inputs -----------------------------------------------------------
  # 0-1 Load data from path & check whether it has all blocks
  #     if the path is not valid, load() will throw an error!
  load(path)
  if (any(!exists('clin') & !exists('cnv') & !exists('mirna') & 
          !exists('rna')  & !exists('mutation'))) {
    stop("'path' led to DF not having 'rna'/ 'cnv'/ 'mirna'/ 'clin'/ 'mutation' as block!")
  }
  
  # 0-2 Check that the response is in the clincial block!
  if (!(response %in% colnames(clin))) stop("Clin Block has no 'response' feature")
  
  # 0-3 Check that response is binary!
  if (length(levels(as.factor(clin[response][,1]))) != 2) {
    stop("'response' doesn't have 2 levels! Choose binary response from 'clin'!")
  }
  
  # [1] Create single DF -------------------------------------------------------
  # 1-1 Extract the response from the clinical block & rm the col from the DF
  response_      <- clin[response]
  clin[response] <- NULL
  
  # 1-2 As 'cnv' and 'rna' blocks share some colnames, rename rna colnames
  colnames(rna) <- paste0(colnames(rna), "_rna")
  
  # 1-3 Bind the single blocks to a big DF! [1. column is response]
  df <- cbind(response_, clin, cnv, rna, mirna, mutation)
  
  # 1-4 Recode the response as factor!
  df[,1] <- factor(df[,1], levels = c(0, 1))
  
  # 1-5 Collect the colnames of the single blocks in list:
  block_variables <- list("clin_block"     = colnames(clin),
                          "cnv_block"      = colnames(cnv),
                          "rna_block"      = colnames(rna),
                          "mutation_block" = colnames(mutation),
                          "mirna_block"    = colnames(mirna))
  
  # [2] Return list with the df & the colnames of the single blocks ------------
  return(list("data"        = df,
              "block_names" = block_variables))
}

dataa <- load_data_extract_block_names()

split_data_to_k_sets <- function(data_and_names = dataa, k = 5, seed = 12345,
                                 scenario = 1) {
  " Function to split the data in 'data_and_names' into k Train- and Testsplits!
    
    Args:
      data_and_names (list) : list filled with 'data' and 'block_names', where 
                              data is a 'data.frame' & 'block_names' is a list 
                              w/ the names 'clin_block', 'cnv_block', 'rna_block'
                                           'mutation_block' & 'mirna_block'
      k    (int)            : Number of Splits we want to have! Must be int > 1!
      seed (int)            : Integer to keep results repoducible - needed when
                              splitting data into k test-train splits!
      scenario (int)        : Scenario that we use to induce blockwise 
                              missingness - must be in (1, 2, 3, 4)!
    
    Return:
      
  "
  # [0] Check Inputs  ----------------------------------------------------------
  # 0-1 Check 'data_and_names'
  # 0-1-1 Check 'data_and_names' for 2 entrances of class DF & list
  assert_list(data_and_names, len = 2)
  
  # 0-1-2 Check the names of the elements in 'data_and_names'
  if (!"data" %in% names(data_and_names)) stop("'data_and_names' needs a 'data' entrance")
  if (!"block_names" %in% names(data_and_names)) stop("'data_and_names' needs a 'block_names' entrance")
  
  # 0-1-3 Check the types of the 'data' & 'block_names' in 'data_and_names'
  assert_data_frame(data_and_names$data, min.rows = k)
  assert_list(data_and_names$block_names, len = 5)
  
  # 0-1-4 Check for correct names in 'block_names' list in 'data_and_names'
  block_names_     <- c('clin_block', 'cnv_block', 'rna_block', 'mutation_block', 'mirna_block')
  all_block_names_ <- sapply(block_names_, function(x) x %in% names(data_and_names$block_names)) 
  if (!all(all_block_names_)) stop("'block_names' block is not completly!")
  
  # 0-1-5 First Column in data must be of type factor with exactly 2 levels!
  assert_factor(data_and_names$data[,1], n.levels = 2, empty.levels.ok = F)
  
  # 0-2 Check 'k', 'scenario' & 'seed'
  assert_int(k, lower = 2, upper = nrow(data_and_names$data))
  assert_int(scenario, lower = 1, upper = 4)
  assert_int(seed)
  
  # [1] Split the data into k equally sized folds  -----------------------------
  # 1-1 Get the fold IDs for each of the k partioions
  set.seed(seed)
  fold_ids <- createFolds(data_and_names$data[,1], k = k)
  
  # 1-2 Split data into k distinct test- & train-sets! Put the test-train pairs 
  #    together into a list w/ train and test! --> Create k of these!
  splitted_dfs <- lapply(seq_len(k), function(x) {
    tmp_train_ids <- unlist(fold_ids[-x])
    tmp_test_ids  <- unlist(fold_ids[x])
    
    list("train" = data_and_names$data[tmp_train_ids,],
         "test"  = data_and_names$data[tmp_test_ids,])
  })
  
  # 1-3 Overwrite the original unsplit data in 'data_and_names$data' with the 
  #    'splitted_dfs' list! 
  data_and_names$data <- splitted_dfs
  
  # [2] To each pair of test and trainset, we induce blockwise- missingness to 
  #     the trainset according to the 'scenario' argument
  # --> Blockwise missingness only induced to the traindata!
  # --> Dimensions stay the same, but amount of NAs per row increases!
  
  # 2-1 Scenario 1
  if (scenario == 1) {
    data_and_names_w_blockwise <- induce_blockmiss_1(data_and_names, seed = seed)
    
    return(data_and_names_w_blockwise)
  }
  
}

induce_blockmiss_1 <- function(data_and_names, seed) {
  " Induce blockwise missingness according to scenario 1 to the train data!
    
    Args:
      data_and_names (list) : list filled with 'data' & 'block_names'! 
                              'data' itself is a list again with k entrances, 
                              where each of these entrances consits of 'train' &
                              'test' [both 'data.frames'].
      seed (int)            : Seed for reproducibility!
      
    Return:
      data_and_names (list): list with the exact same layout as the input, BUT 
                             the data entrance has been induced w/ blockwise 
                             missingness
  "
  # [0] Check Inputs  ----------------------------------------------------------
  # 0-1 Check 'data_and_names' 
  # 0-1-1 Must be a list filled with 2 more lists!
  assert_list(data_and_names, len = 2)
  assert_list(data_and_names$block_names, len = 5)
  assert_list(data_and_names$data, min.len = 2)
  
  # 0-1-2 'data_and_names$data' must be a list aswell, filled with dataframes
  all_DFs_ <- sapply(data_and_names$data, function(x) {
    test_data_frame(x$train) & test_data_frame(x$test)
  })
  
  if (!all(all_DFs_)) stop("Not all elements in 'data_and_names$data' are of type dataframe")
  
  # 0-2 Check seed
  assert_int(seed)
  
  # [1] Induce the blockwise missingness  --------------------------------------
  # 1-1 Loop over each train-DF in 'data_and_names$data' and assign the 
  #     observations to different folds [each fold has own observed blocks!]
  for (j in 1:length(data_and_names$data)) {
    
    # 1-1-1 We want ~equally sized folds, need to check whether the amount of 
    #       train obs. in current fold is dividable by 4, if not, we need to make 
    #       random blocks bigger by 1 observations, so that we can assign all 
    #       observations into different folds and have no remainers!
    amount_train_obs <- nrow(data_and_names$data[[j]]$train)
    obs_per_fold     <- rep(floor(amount_train_obs / 4), 4)      # rounded obs. per fold!
    
    # remaining obs, that need to be put in random folds!
    remaining_obs    <- amount_train_obs - (obs_per_fold[1] * 4) 
    
    # If we have remaining observations, we assign these randomly to any fold!
    if (remaining_obs > 0) {
      folds_w_extra_obs <- base::sample(x       = seq_len(length(obs_per_fold)), 
                                        size    = remaining_obs, 
                                        replace = FALSE)
      
      # Give the folds, that were sampled above, the extra observation so there is
      # no trainobservation that can not be assigned to any fold!
      obs_per_fold[folds_w_extra_obs] <- obs_per_fold[folds_w_extra_obs] + 1 
    }
    print("Amount of observations per fold:")
    print(obs_per_fold)
    
    # 1-2 Sample the observed blocks for each fold!
    set.seed(seed)
    observed_blocks <- sample(c(rep("Clin, cnv",      obs_per_fold[1]), 
                                rep("Clin, rna",      obs_per_fold[2]),
                                rep("Clin, mirna",    obs_per_fold[3]),
                                rep("Clin, mutation", obs_per_fold[4])),
                              amount_train_obs, replace = FALSE)
    
    # 1-3 Remove the observed blocks that were not observed for folds
    #     (e.g. all obs. w/ "Clin, cnv", will only have these, the rest is NA!)
    data_and_names$data[[j]]$train[which(observed_blocks == "Clin, cnv"), 
                                   c(data_and_names$block_names$mutation_block,
                                     data_and_names$block_names$rna_block,
                                     data_and_names$block_names$mirna_block)] <- NA
    
    data_and_names$data[[j]]$train[which(observed_blocks == "Clin, rna"), 
                                   c(data_and_names$block_names$mutation_block,
                                     data_and_names$block_names$cnv_block,
                                     data_and_names$block_names$mirna_block)] <- NA
    
    data_and_names$data[[j]]$train[which(observed_blocks == "Clin, mirna"), 
                                   c(data_and_names$block_names$mutation_block,
                                     data_and_names$block_names$rna_block,
                                     data_and_names$block_names$cnv_block)] <- NA
    
    data_and_names$data[[j]]$train[which(observed_blocks == "Clin, mutation"), 
                                   c(data_and_names$block_names$cnv_block,
                                     data_and_names$block_names$rna_block,
                                     data_and_names$block_names$mirna_block)] <- NA
  }
  
  # [2] Return the list, but with block wise missing data now!  ----------------
  return(data_and_names)
}
