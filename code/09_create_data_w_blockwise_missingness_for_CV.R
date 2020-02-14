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
                             missingness according to scenario 1 and the omics
                             block names are replaced by 'A', 'B', ...
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
  # 1-1 Loop over each of the k test-train-splits in 'data_and_names$data' and 
  #     assign the observations to different folds [each has own observed blocks!]
  for (j in 1:length(data_and_names$data)) {
    
    # 1-1-1 The 4 different folds need to have about the same amount of 
    #       observations. If the amount of train obs. in train-split is not
    #       dividable by 44, random folds get enlarged by 1 observation, so that
    #       all obs belong in a fold & all folds are about equally sized [+-1]!
    n_train    <- nrow(data_and_names$data[[j]]$train) # total train obs.
    n_per_fold <- rep(floor(n_train / 4), 4)           # rounded obs. per fold
    
    # Remaining Observations, that need to be put in random folds!
    remaining_obs <- n_train - sum(n_per_fold) 
    
    # If we have remaining observations, we assign these randomly to any fold!
    if (remaining_obs > 0) {
      
      # Sample randomly the folds, that get the extra observation
      set.seed(seed)
      folds_w_extra_obs <- base::sample(x       = seq_len(length(n_per_fold)), 
                                        size    = remaining_obs, 
                                        replace = FALSE)
      
      # Give the folds, the extra observation so there is no trainobservation 
      # that can not be assigned to any fold!
      n_per_fold[folds_w_extra_obs] <- n_per_fold[folds_w_extra_obs] + 1 
    }
    
    # 1-2 Print the amount of observations per fold
    print("Amount of observations per fold:")
    print(n_per_fold)
    
    # 1-3 Sample the observed blocks for each fold!
    set.seed(seed)
    observed_blocks <- sample(c(rep("Clin, cnv",      n_per_fold[1]), 
                                rep("Clin, rna",      n_per_fold[2]),
                                rep("Clin, mirna",    n_per_fold[3]),
                                rep("Clin, mutation", n_per_fold[4])),
                              n_train, replace = FALSE)
    
    # 1-4 Replace the observed values from the blocks that weren't sampled for
    #     certain folds with NA!
    #     (e.g. all obs. w/ "Clin, cnv", will only have these, the rest is NA!)
    
    # 1-4-1 Subset all Values with NA except from 'clin' & 'cnv' block!
    data_and_names$data[[j]]$train[which(observed_blocks == "Clin, cnv"), 
                                   c(data_and_names$block_names$mutation_block,
                                     data_and_names$block_names$rna_block,
                                     data_and_names$block_names$mirna_block)] <- NA
    
    # 1-4-2 Subset all Values with NA except from 'clin' & 'rna' block!
    data_and_names$data[[j]]$train[which(observed_blocks == "Clin, rna"), 
                                   c(data_and_names$block_names$mutation_block,
                                     data_and_names$block_names$cnv_block,
                                     data_and_names$block_names$mirna_block)] <- NA
    
    # 1-4-3 Subset all Values with NA except from 'clin' & 'mirna' block!
    data_and_names$data[[j]]$train[which(observed_blocks == "Clin, mirna"), 
                                   c(data_and_names$block_names$mutation_block,
                                     data_and_names$block_names$rna_block,
                                     data_and_names$block_names$cnv_block)] <- NA
    
    # 1-4-4 Subset all Values with NA except from 'clin' & 'mutation' block!
    data_and_names$data[[j]]$train[which(observed_blocks == "Clin, mutation"), 
                                   c(data_and_names$block_names$cnv_block,
                                     data_and_names$block_names$rna_block,
                                     data_and_names$block_names$mirna_block)] <- NA
  }
  
  # [2] Modify the 'block_names' of 'data_and_names'  --------------------------
  # 2-1 Instead of cnv', 'rna', ... subset the ommics block names by 'A', 'B'
  names(data_and_names$block_names)[2] <- "A"
  names(data_and_names$block_names)[3] <- "B"
  names(data_and_names$block_names)[4] <- "C"
  names(data_and_names$block_names)[5] <- "D"
  
  # [3] Return the list, but with block wise missing data now!  ----------------
  return(data_and_names)
}

induce_blockmiss_2 <- function(data_and_names, seed) {
  " Induce blockwise missingness according to scenario 2 to the train data!
    
    Args:
      data_and_names (list) : list filled with 'data' & 'block_names'! 
                              'data' itself is a list again with k entrances, 
                              where each of these entrances consits of 'train' &
                              'test' [both 'data.frames'].
      seed (int)            : Seed for reproducibility!
      
    Return:
      data_and_names (list): list with the exact same layout as the input, BUT 
                             the data entrance has been induced w/ blockwise 
                             missingness according to scenario 2 and the omics
                             block names are replaced by 'A', 'B', ...
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
  # 1-1 Loop over each of the k test-train-splits in 'data_and_names$data' and 
  #     assign the observations to different folds [each has own observed blocks!]
  for (j in 1:length(data_and_names$data)) {
    
    # 1-1-1 The 4 different folds need to have about the same amount of 
    #       observations. If the amount of train obs. in train-split is not
    #       dividable by 44, random folds get enlarged by 1 observation, so that
    #       all obs belong in a fold & all folds are about equally sized [+-1]!
    n_train    <- nrow(data_and_names$data[[j]]$train) # total train obs.
    n_per_fold <- rep(floor(n_train / 4), 4)           # rounded obs. per fold
    
    # Remaining Observations, that need to be put in random folds!
    remaining_obs <- n_train - sum(n_per_fold) 
    
    # If we have remaining observations, we assign these randomly to any fold!
    if (remaining_obs > 0) {
      
      # Sample randomly the folds, that get the extra observation
      set.seed(seed)
      folds_w_extra_obs <- base::sample(x       = seq_len(length(n_per_fold)), 
                                        size    = remaining_obs, 
                                        replace = FALSE)
      
      # Give the folds, the extra observation so there is no trainobservation 
      # that can not be assigned to any fold!
      n_per_fold[folds_w_extra_obs] <- n_per_fold[folds_w_extra_obs] + 1 
    }
    
    # 1-2 Print the amount of observations per fold
    print("Amount of observations per fold:")
    print(n_per_fold)
    
    # 1-3 Randomly assign the omics blocks to the letters 'A', 'B', 'C', 'D', 
    #     as SCEANRIO2, highly depens on which block is which letter!
    #     ['A' only observed in 1.fold, whereas 'D' is observed in all folds]
    set.seed(seed)
    letter_feas <- sample(c("cnv_block", "rna_block", "mutation_block", "mirna_block"), 
                          4, replace = FALSE)
    names(letter_feas) <- c("A", "B", "C", "D")
    
    # 1-4 Sample the observed blocks for each fold!
    set.seed(seed)
    observed_blocks <- sample(c(rep("A", n_per_fold[1]), 
                                rep("B", n_per_fold[2]),
                                rep("C", n_per_fold[3]),
                                rep("D", n_per_fold[4])),
                              n_train, replace = FALSE)
    
    # 1-5 Replace the observed values from the blocks that weren't sampled for
    #     certain folds with NA!
    #     (e.g. all obs. w/ "A", will only have the block A, the rest is NA!)
    
    # 1-4-1 Obs. of Block 'A': Subset all Values w/ NA except from 'clin' & Block A!
    data_and_names$data[[j]]$train[which(observed_blocks == "A"), 
                                   c(data_and_names$block_names[[letter_feas["B"]]],
                                     data_and_names$block_names[[letter_feas["C"]]],
                                     data_and_names$block_names[[letter_feas["D"]]])] <- NA
    
    # 1-4-2 Obs. of Block 'B': Subset all Values w/ NA except from 'clin' & Block A & B!
    data_and_names$data[[j]]$train[which(observed_blocks == "B"), 
                                   c(data_and_names$block_names[[letter_feas["C"]]],
                                     data_and_names$block_names[[letter_feas["D"]]])] <- NA
    
    ## 1-4-3 Obs. of Block 'B': Subset all Values w/ NA except from 'clin' & Block A & B & C!
    data_and_names$data[[j]]$train[which(observed_blocks == "C"), 
                                   c(data_and_names$block_names[[letter_feas["D"]]])] <- NA
    
    # 1-4-4 Obs. from Block 'D': Nothing has to be subsetted! 
  }
  
  # [2] Rename the 'block_names' in 'data_and_names'  --------------------------
  # 2-1 As we shuffled which letter, is which block, we need to adjust the block
  #     Names aswell!
  names(data_and_names$block_names)[2] <- names(letter_feas)[letter_feas == "cnv_block"]
  names(data_and_names$block_names)[3] <- names(letter_feas)[letter_feas == "rna_block"]
  names(data_and_names$block_names)[4] <- names(letter_feas)[letter_feas == "mutation_block"]
  names(data_and_names$block_names)[5] <- names(letter_feas)[letter_feas == "mirna_block"]
  
  # [3] Return the list, but with block wise missing data now!  ----------------
  return(data_and_names)
}

induce_blockmiss_3 <- function(data_and_names, seed) {
  " Induce blockwise missingness according to scenario 3 to the train data!
    
    Args:
      data_and_names (list) : list filled with 'data' & 'block_names'! 
                              'data' itself is a list again with k entrances, 
                              where each of these entrances consits of 'train' &
                              'test' [both 'data.frames'].
      seed (int)            : Seed for reproducibility!
      
    Return:
      data_and_names (list): list with the exact same layout as the input, BUT 
                             the data entrance has been induced w/ blockwise 
                             missingness according to scenario 3 and the omics
                             block names are replaced by 'A', 'B', ...
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
  # 1-1 Loop over each of the k test-train-splits in 'data_and_names$data' and 
  #     assign the observations to different folds [each has own observed blocks!]
  for (j in 1:length(data_and_names$data)) {
    
    # 1-1-1 The 4 different folds need to have about the same amount of 
    #       observations. If the amount of train obs. in train-split is not
    #       dividable by 44, random folds get enlarged by 1 observation, so that
    #       all obs belong in a fold & all folds are about equally sized [+-1]!
    n_train    <- nrow(data_and_names$data[[j]]$train) # total train obs.
    n_per_fold <- rep(floor(n_train / 4), 4)           # rounded obs. per fold
    
    # Remaining Observations, that need to be put in random folds!
    remaining_obs <- n_train - sum(n_per_fold) 
    
    # If we have remaining observations, we assign these randomly to any fold!
    if (remaining_obs > 0) {
      
      # Sample randomly the folds, that get the extra observation
      set.seed(seed)
      folds_w_extra_obs <- base::sample(x       = seq_len(length(n_per_fold)), 
                                        size    = remaining_obs, 
                                        replace = FALSE)
      
      # Give the folds, the extra observation so there is no trainobservation 
      # that can not be assigned to any fold!
      n_per_fold[folds_w_extra_obs] <- n_per_fold[folds_w_extra_obs] + 1 
    }
    
    # 1-2 Print the amount of observations per fold
    print("Amount of observations per fold:")
    print(n_per_fold)
    
    # 1-3 Randomly assign the 'observed' blocks to the 4 different folds!
    #     The Index of TRUE / FALSE are indicators, whether the blocks are observed
    #     [1] = clinical; [2] = CNV; [3] = RNA; [4] = Mutation; [5] = Mirna
    set.seed(seed)
    fold1_obs <- sample(c(TRUE, FALSE), size = 5, replace = T, prob = c(2/3, 1/3))
    names(fold1_obs) <- c("clin_block", "cnv_block", "rna_block", "mirna_block", "mutation_block")
    set.seed(seed + 1)
    fold2_obs <- sample(c(TRUE, FALSE), size = 5, replace = T, prob = c(2/3, 1/3))
    names(fold2_obs) <- c("clin_block", "cnv_block", "rna_block", "mirna_block", "mutation_block")
    set.seed(seed + 2)
    fold3_obs <- sample(c(TRUE, FALSE), size = 5, replace = T, prob = c(2/3, 1/3))
    names(fold3_obs) <- c("clin_block", "cnv_block", "rna_block", "mirna_block", "mutation_block")
    set.seed(seed + 3)
    fold4_obs <- sample(c(TRUE, FALSE), size = 5, replace = T, prob = c(2/3, 1/3))
    names(fold4_obs) <- c("clin_block", "cnv_block", "rna_block", "mirna_block", "mutation_block")
    
    # 1-4 Sample the observed blocks for each fold!
    set.seed(seed)
    observed_blocks <- sample(c(rep("fold1", n_per_fold[1]), 
                                rep("fold2", n_per_fold[2]),
                                rep("fold3", n_per_fold[3]),
                                rep("fold4", n_per_fold[4])),
                              n_train, replace = FALSE)
    
    # 1-5 Replace the observed values from the blocks that weren't sampled for
    #     certain folds with NA! (e.g. all obs. from  "fold1", will only 
    #     have the blocks that were sampled for "fold1" rest replaces by NA!)
    
    # 1-5-1 Obs. of 'fold1': Subset all Values w/ NA that were not sampled to be observed!
    cols_to_rm <- unlist(sapply(names(fold1_obs)[which(!fold1_obs)], function(x) {
      data_and_names$block_names[[x]]
    }))
    data_and_names$data[[j]]$train[which(observed_blocks == "fold1"), 
                                   c(cols_to_rm)] <- NA
    
    # 1-5-2 Obs. of 'fold2': Subset all Values w/ NA that were not sampled to be observed!
    cols_to_rm <- unlist(sapply(names(fold2_obs)[which(!fold2_obs)], function(x) {
      data_and_names$block_names[[x]]
    }))
    data_and_names$data[[j]]$train[which(observed_blocks == "fold2"), 
                                   c(cols_to_rm)] <- NA
    
    # 1-5-3 Obs. of 'fold3': Subset all Values w/ NA that were not sampled to be observed!
    cols_to_rm <- unlist(sapply(names(fold3_obs)[which(!fold3_obs)], function(x) {
      data_and_names$block_names[[x]]
    }))
    data_and_names$data[[j]]$train[which(observed_blocks == "fold3"), 
                                   c(cols_to_rm)] <- NA
    
    # 1-5-4 Obs. of 'fold4': Subset all Values w/ NA that were not sampled to be observed!
    cols_to_rm <- unlist(sapply(names(fold4_obs)[which(!fold4_obs)], function(x) {
      data_and_names$block_names[[x]]
    }))
    data_and_names$data[[j]]$train[which(observed_blocks == "fold4"), 
                                   c(cols_to_rm)] <- NA
  }
  
  # [2] Modify the 'block_names' of 'data_and_names'  --------------------------
  # 2-1 Instead of cnv', 'rna', ... subset the ommics block names by 'A', 'B'
  names(data_and_names$block_names)[2] <- "A"
  names(data_and_names$block_names)[3] <- "B"
  names(data_and_names$block_names)[4] <- "C"
  names(data_and_names$block_names)[5] <- "D"
  
  # [3] Return the list, but with block wise missing data now!  ----------------
  return(data_and_names)
}

induce_blockmiss_4 <- function(data_and_names, seed) {
  " Induce blockwise missingness according to scenario 4 to the train data!
    CAUTION - differs from the others as it holds 2 folds only!
    
    Args:
      data_and_names (list) : list filled with 'data' & 'block_names'! 
                              'data' itself is a list again with k entrances, 
                              where each of these entrances consits of 'train' &
                              'test' [both 'data.frames'].
      seed (int)            : Seed for reproducibility!
      
    Return:
      data_and_names (list): list with the exact same layout as the input, BUT 
                             the data entrance has been induced w/ blockwise 
                             missingness according to scenario 4 and the omics
                             block names are replaced by 'A', 'B', ...
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
  # 1-1 Loop over each of the k test-train-splits in 'data_and_names$data' and 
  #     assign the observations to different folds [each has own observed blocks!]
  for (j in 1:length(data_and_names$data)) {
    
    # 1-1-1 The 2 different folds need to have about the same amount of 
    #       observations. If the amount of train obs. in train-split is not
    #       dividable by 44, random folds get enlarged by 1 observation, so that
    #       all obs belong in a fold & all folds are about equally sized [+-1]!
    n_train    <- nrow(data_and_names$data[[j]]$train) # total train obs.
    n_per_fold <- rep(floor(n_train / 2), 2)           # rounded obs. per fold
    
    # Remaining Observations, that need to be put in random folds!
    remaining_obs <- n_train - sum(n_per_fold) 
    
    # If we have remaining observations, we assign these randomly to any fold!
    if (remaining_obs > 0) {
      
      # Sample randomly the folds, that get the extra observation
      set.seed(seed)
      folds_w_extra_obs <- base::sample(x       = seq_len(length(n_per_fold)), 
                                        size    = remaining_obs, 
                                        replace = FALSE)
      
      # Give the folds, the extra observation so there is no trainobservation 
      # that can not be assigned to any fold!
      n_per_fold[folds_w_extra_obs] <- n_per_fold[folds_w_extra_obs] + 1 
    }
    
    # 1-2 Print the amount of observations per fold
    print("Amount of observations per fold:")
    print(n_per_fold)
    
    # 1-3 Randomly paste 2 blocks together, so that there are only 2 blocks left!
    #     First 2 and last 2 of here will be sampled!
    set.seed(seed)
    blocks_together <- sample(c("cnv_block", "rna_block", "mutation_block", "mirna_block"),
                              size = 4, replace = F)
    
    # 1-3-1 Bind the colnams of the first 2 sampled blocks!
    block_A_names <- c(data_and_names$block_names[[blocks_together[1]]],
                       data_and_names$block_names[[blocks_together[2]]])
    
    # 1-3-2 Bind the colnames if the last 2 sampled blocks!
    block_B_names <- c(data_and_names$block_names[[blocks_together[3]]],
                       data_and_names$block_names[[blocks_together[4]]])
    
    # 1-4 Sample the observed blocks for each fold!
    set.seed(seed)
    observed_blocks <- sample(c(rep("fold1", n_per_fold[1]), 
                                rep("fold2", n_per_fold[2])),
                              n_train, replace = FALSE)
    
    # 1-5 Replace the observed values from the blocks that weren't sampled for
    #     certain folds with NA! (e.g. all obs. from  "fold1", will only 
    #     have the blocks that were sampled for "fold1" rest replaces by NA!)
    
    # 1-5-1 Obs. of 'fold1': Subset all Values w/ NA that were in block A!
    data_and_names$data[[j]]$train[which(observed_blocks == "fold1"), 
                                   block_B_names] <- NA
    
    # 1-5-2 Obs. of 'fold1': Subset all Values w/ NA that were in block B!
    data_and_names$data[[j]]$train[which(observed_blocks == "fold2"), 
                                   block_A_names] <- NA
  }
  
  # [2] Modify the 'block_names' of 'data_and_names'  --------------------------
  # 2-1 Instead of cnv', 'rna', ... subset the ommics block names by 'A', 'B' & 
  #     Put the blocknames we pasted together together!
  # 2-1-1 Remove the old entrances
  data_and_names$block_names$cnv_block      <- NULL
  data_and_names$block_names$rna_block      <- NULL
  data_and_names$block_names$mutation_block <- NULL
  data_and_names$block_names$mirna_block    <- NULL
  
  # 2-1-1 Add the new entrances!
  data_and_names$block_names$A <- block_A_names
  data_and_names$block_names$B <- block_B_names
  
  # [3] Return the list, but with block wise missing data now!  ----------------
  return(data_and_names)
}

split_data_to_k_sets <- function(data_and_names, scenario, k = 5, seed = 12345) {
  " Function to split the data in 'data_and_names' into k Train- and Testsplits!
    And induce the blockwise missingness into Train according to 'scenario'!
    
    Args:
      data_and_names (list) : list filled with 'data' and 'block_names', where 
                              data is a 'data.frame' & 'block_names' is a list 
                              w/ the names 'clin_block', 'cnv_block', 'rna_block'
                                           'mutation_block' & 'mirna_block'
      scenario (int)        : Scenario that we use to induce blockwise 
                              missingness - must be in (1, 2, 3, 4)!
      k    (int)            : Number of Splits we want to have! Must be int > 1!
      seed (int)            : Integer to keep results repoducible - needed when
                              splitting data into k test-train splits!
      
    Return:
      data_and_names (list) filled with 'data' & 'block_names'. 
      -'data' is a set of k-lists, where each list holds a specific test-train
        split, where the testdata is fully observed and the train data was induced 
        with a type of blockwise missingness!
      - 'block_names' has the colanmes of all blocks!
  "
  # [0] Check Inputs  ----------------------------------------------------------
  # 0-1 Check 'data_and_names'
  # 0-1-1 Check for being a list w/ 2 entrances
  assert_list(data_and_names, len = 2)
  
  # 0-1-2 Check names of the 2 entraces in 'data_and_names'
  if (!"data" %in% names(data_and_names)) stop("'data_and_names' needs a 'data' entrance")
  if (!"block_names" %in% names(data_and_names)) stop("'data_and_names' needs a 'block_names' entrance")
  
  # 0-1-3 Check the types of 'data' & 'block_names' in 'data_and_names'
  assert_data_frame(data_and_names$data, min.rows = k)
  assert_list(data_and_names$block_names, len = 5)
  
  # 0-1-4 Check for correct names in 'block_names' in 'data_and_names'
  block_names_     <- c('clin_block', 'cnv_block', 'rna_block', 'mutation_block', 'mirna_block')
  all_block_names_ <- sapply(block_names_, function(x) x %in% names(data_and_names$block_names)) 
  if (!all(all_block_names_)) {
    print(all_block_names_)
    stop("'block_names' block is not completly!")
  }
  
  # 0-1-5 First Column in data must be of type factor with exactly 2 levels!
  assert_factor(data_and_names$data[,1], n.levels = 2, empty.levels.ok = F)
  
  # 0-2 Check 'k', 'scenario' & 'seed' for bein integers
  assert_int(k, lower = 2, upper = nrow(data_and_names$data))
  assert_int(scenario, lower = 1, upper = 4)
  assert_int(seed)
  
  # [1] Split the data into k equally sized folds  -----------------------------
  # 1-1 Get the fold IDs for each of the k partioions
  set.seed(seed)
  fold_ids <- createFolds(data_and_names$data[,1], k = k)
  
  # 1-2 Split data into k distinct test- & train-sets! 
  #     Put the test-train pairs together into a list w/ 'train' and 'test' !
  #     --> k list entrances w/ k test-train-splits of data!
  splitted_dfs <- lapply(seq_len(k), function(x) {
    tmp_train_ids <- unlist(fold_ids[-x])
    tmp_test_ids  <- unlist(fold_ids[x])
    
    list("train" = data_and_names$data[tmp_train_ids,],
         "test"  = data_and_names$data[tmp_test_ids,])
  })
  
  # 1-3 Overwrite the original (not splitted) data in 'data_and_names$data' with
  #     the 'splitted_dfs' list (each entrance has its own test-train-split)! 
  data_and_names$data <- splitted_dfs
  
  # [2] Induce Blockwise Missingness  ------------------------------------------
  # To each pair of test and trainset, we induce blockwise- missingness to 
  # the trainset according to the 'scenario' argument.
  #     --> Blockwise missingness only induced to the traindata!
  #     --> Dimensions stay the same, but amount of NAs per row/ col increase!
  
  # 2-1 Scenario 1
  if (scenario == 1) {
    data_and_names_w_blockwise <- induce_blockmiss_1(data_and_names, seed = seed)
    return(data_and_names_w_blockwise)
  }
  
  # 2-2 Scenario 2
  if (scenario == 2) {
    data_and_names_w_blockwise <- induce_blockmiss_2(data_and_names, seed = seed)
    return(data_and_names_w_blockwise)
  }
  
  # 2-3 Scenario 3
  if (scenario == 3) {
    data_and_names_w_blockwise <- induce_blockmiss_3(data_and_names, seed = seed)
    return(data_and_names_w_blockwise)
  }
  
  # 2-4 Scenario 4
  if (scenario == 4) {
    data_and_names_w_blockwise <- induce_blockmiss_4(data_and_names, seed = seed)
    return(data_and_names_w_blockwise)
  }
}

# Run TestTrainSplitting - and induce blockwise missingness to TrainSets  ------
DFs_w_gender <- c("BLCA", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC","LGG", 
                  "LUAD", "LUSC", "PAAD", "SARC", "SKCM", "STAD")
data_path    <- "./data/external/Dr_Hornung/subsetted_12345/" # Path to the data
response_    <- "gender"                                      # response from 'clin' block
seed         <- 1234                                          # seed for reprducibility!

for (curr_df in "BLCA") {
  
  # [1] Create the path to the subsetted DF, we want to induce blockwise missingness!
  curr_path <- paste0(data_path, curr_df, "_subset.RData")
  
  # [2] Load the Data
  curr_data <- load_data_extract_block_names(path = curr_path, response = response_)
  print("Dimensions of loaded Data")
  print(dim(curr_data$data))
  
  # [3] Induce Blockwise missingness for each Scenario [1, 2, 3, 4]
  curr_data_1 <- split_data_to_k_sets(curr_data, scenario = 1, k = 5, seed = seed)
  curr_data_2 <- split_data_to_k_sets(curr_data, scenario = 2, k = 5, seed = seed)
  curr_data_3 <- split_data_to_k_sets(curr_data, scenario = 3, k = 5, seed = seed)
  curr_data_4 <- split_data_to_k_sets(curr_data, scenario = 4, k = 5, seed = seed)
  
  # [4] Save the Lists with the blockwise missingness!
  # 4-1 Define the basic path to save the results!
  path_to_save <- paste0(data_path, "missingness_", seed, "/", curr_df)
  
  save(curr_data_1, file = paste0(path_to_save, "_1.RData"))
  save(curr_data_1, file = paste0(path_to_save, "_2.RData"))
  save(curr_data_1, file = paste0(path_to_save, "_3.RData"))
  save(curr_data_1, file = paste0(path_to_save, "_4.RData"))
}
