" Investigate how good the predicitive peroformance of a model is on a single 
  block / in General! 
    --> full block with the whole features that are originally in the block!
    --> block with only 10% | 5% | 2.5%... of the original features!
    --> all blocks - w/ only 10% | 5% | 2.5%... of original features conecated!

  For each of these Scenarios, we create a DF, that tracks:
    -the dataframe  
    -the block  
    -Dim of Traindata  
    -OOB & TestSet Accuracy + Test F1 Score! 
    -Time in minutes
"
# Set Working Directory and load the needed packages!
setwd("C:/Users/kuche_000/Desktop/MS-Thesis/")
library(randomForest)
library(caret)
library(checkmate)
library(randomForestSRC)

# Names of the usable dataframes (w/ gender in 'clin'-block & 4 omics blocks!)
DFs_w_gender <- c("BLCA", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC","LGG", 
                  "LUAD", "LUSC", "PAAD", "SARC", "SKCM", "STAD")
data_path    <- "./data/external/Dr_Hornung/Data/ProcessedData/"

eval_single_block_subsets <- function(DFs_w_gender, seed_to_subset, fraction) {
  "
  Function to evaluate single block performance for different DFs ['DFs_w_gender']. 
  Each DF has 5 different blocks [clin, mirna, cnv, rna, mutation]! 
  For each block [except for clin] we subset the feature space by only keeping a
  'fraction' of the original feature space & get the performance, so that we can
  find out which subset to use in our final study!
  For this we do 5 fold CV to each single block and save the results, obtained
  when fitting a standard rfsrc model to it w/ 250 trees!
  
  Args:
    DFs_w_gender (vector) : Vector filled with strings of DFs that have 'gender'
                            within their clinical block! These DF names have to 
                            be in './data/external/Dr_Hornung/Data/ProcessedData/'
    seed_to_subset (int)  : Seed used, when subsetting the feature space, so it
                            is reproducible!
    fraction (double)     : How much of the original feature space shall be kept?
                            [fraction = 0.2 --> keep 20% and remove 80% of feas
                                                in each block except for clin!]
  Return:
    data.frame with: 'Data'      - data used
                     'Block'     - current block that was used to fit RF
                     'train_dim' - dimension of the data used to fit the RF
                     'OOB_Acc'   - Metrics obtained from OOB / CV Testset
                     'Test_Acc'  -               - ' - 
                     'Test_F1'   -               - ' -
                     'Fold'      - which fold of the 5 fold CV
                     'Time'      - how long did it take to fit a RF on the block
                     'Fraction'  - fraction used to subset the single blocks!
                     'subset_seed' -  seed used to create subset!
    will be saved to 'docs/CV/gender/explorative_subsets/' the name of the file
    itself is a mixture of seed and fraction!
    
  "
  # [0] Check Arguments
  assert_vector(DFs_w_gender, min.len = 1, unique = TRUE)
  assert_int(seed_to_subset)
  assert_double(fraction, lower = 0, upper = 1)
  
  # [1] Empty DF - for all results of the Evaluation:
  eval_res <- data.frame("Data"      = character(), "Block"    = numeric(), 
                         "train_dim" = numeric(),   "OOB_Acc"  = numeric(),   
                         "Test_Acc"  = numeric(),   "Test_F1"  = numeric(), 
                         "Fold"      = numeric(),   "Time"     = numeric(), # min
                         "Fraction"  = numeric(),   "subset_seed" = numeric(),
                         stringsAsFactors = F)
  
  # 1-1 Fixed Datapath, where we have all our raw Dataframes!
  data_path    <- "./data/external/Dr_Hornung/Data/ProcessedData/"
  
  # [2] Loop over all DFs w/ gender in their clinical variables!
  for (df in DFs_w_gender) {
    
    writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
    
    # 2-1 Load 'df' & only keep names of the omics blocks in 'df'!
    omics_blocks <- load(paste0(data_path, df, ".Rda"))
    omics_blocks <- omics_blocks[-which(omics_blocks %in% c("targetvar"))]
    
    # 2-1-1 Extract 'gender' as response and also remove it from the 'clin' block!
    resp <- clin["gender"]
    clin <- clin[-which(colnames(clin) == "gender")]
    
    # [2-2] Loop over all blocks in 'omics_blocks' of current 'df'
    for (curr_block in omics_blocks) {
      
      writeLines(paste0("Current OmicsBlock: \n", curr_block))
      
      # 2-2-1 If the block is not clinical only keep 'fraction' of the original 
      #       features!
      curr_data <- eval(as.symbol(curr_block))
      
      if (curr_block != "clin") {
        set.seed(seed_to_subset)
        sampled_cols <- sample(seq_len(ncol(curr_data)), 
                               round(ncol(curr_data) * fraction))
        assign(curr_block, curr_data[ ,sampled_cols])
      } else {
        sampled_cols <- c(1:ncol(curr_data))
        assign(curr_block, curr_data[ ,sampled_cols])
      }
      
      # 2-2-2 Bind gender from the clin block & 'curr_block' (=omics) 
      DF <- cbind(resp, eval(as.symbol(curr_block)))
      
      # 2-2-3 Convert to DF, and name the response variable
      DF <- as.data.frame(DF)
      colnames(DF)[1] <- "resp"
      DF[,1] <- factor(DF[,1], levels = c(0, 1))
      
      # 2-2-3 Do 5 fold CV
      set.seed(12345)
      fold_ids <- createFolds(DF[,1], k = 5)
      
      for (i in 1:5) {
        
        print(paste("Fold:",as.character(i), "-----------------------------------"))
        
        # 2-2-4 split to Test & Train
        test      <- DF[fold_ids[[i]],]
        train_obs <- which(!(seq_len(nrow(DF)) %in% fold_ids[[i]]))
        train     <- DF[train_obs,]
        
        # 2-2-5 Print dimension of TrainingDF
        dim_train <- paste(dim(train), collapse = " x ")
        print(paste("Dimension of Train-DF:", dim_train))
        
        # 2-2-6 Fit a model [standard settings] on the train DF & take the time:
        start        <- Sys.time()
        curr_forrest <- rfsrc(as.formula(resp ~ .), data = train, seed = 1312,
                              ntree = 250)
        end          <- Sys.time()
        time_diff    <- difftime(end, start, units = 'mins')
        
        # 2-2-7 Evaluate the Model:
        # OOB:      Get OOB Accuracy!
        curr_OOB_ACC <- 1 - curr_forrest$err.rate[length(curr_forrest$err.rate)]
        
        # TESTSET:  Get the ErrorRates w/ TestSet
        prob_preds   <- predict(curr_forrest, test)$predicted
        pred_class   <- ifelse(prob_preds[,1] > 0.5, 0, 1)
        truth        <- as.factor(test[,1])
        pred_class   <- factor(pred_class, levels = levels(truth))
        res          <- cbind(pred_class, truth)
        test_acc     <- sum(res[,1] == res[,2]) / nrow(res)
        test_F1      <- caret::confusionMatrix(data      = pred_class, 
                                               reference = truth)$byClass["F1"]
        if (is.na(test_F1)) test_F1 <- 0 
        
        # 2-2-8 Add results of this block to the Results-DF!
        eval_res[nrow(eval_res) + 1, ] <- c(df, curr_block, dim_train, 
                                            curr_OOB_ACC, test_acc, test_F1, 
                                            i, time_diff, fraction, seed_to_subset)
        
      }
    }
  }
  res_name <- paste0("RF_single_block_frac", fraction, "_seed", seed_to_subset)
  
  write.csv2(eval_res, 
             paste0("./docs/cv_Res/gender/explorative_subsets/", res_name, ".csv"),
             row.names = FALSE) 
}

eval_joint_block_subsets <- function(DFs_w_gender, seed_to_subset, fraction_cnv,
                                     fraction_rna, fraction_mirna,
                                     fraction_mutation, fraction_clin = 1) {
  "
  Function to evaluate joint block performance for different DFs ['DFs_w_gender']. 
  Each DF has 5 different blocks [clin, mirna, cnv, rna, mutation]! 
  For each block we can choose a subset we want to use from this block, then
  we subset each block according to ['fraction_cnv', 'fraction_rna', ...].
  Then we bind these subsetted blocks to a single big DF and use this whole DF
  as input for the RF! 
  With this data we do 5 fold CV to and save the results, obtained
  when fitting a standard rfsrc model to it w/ 250 trees!
  
  Args:
    DFs_w_gender (vector) : Vector filled with strings of DFs that have 'gender'
                            within their clinical block! These DF names have to 
                            be in './data/external/Dr_Hornung/Data/ProcessedData/'
    seed_to_subset (int)  : Seed used, when subsetting the feature space, so it
                            is reproducible!
    fraction_cnv (double) : How much of the original 'cnv' feature space shall 
                            be kept? [fraction = 0.2 --> keep 20% and remove 80%
                                                         of feas in this block]
    fraction_rna (double) : How much of the original 'rna' feature space shall 
                            be kept? [fraction = 0.2 --> keep 20% and remove 80%
                                                         of feas in this block]
    fraction_mirna        : How much of the original 'cnv' feature space shall 
       (double)             be kept? [fraction = 0.2 --> keep 20% and remove 80%
                                                         of feas in this block]
    fraction_mutation     : How much of the original 'cnv' feature space shall 
       (double)             be kept? [fraction = 0.2 --> keep 20% and remove 80%
                                                         of feas in this block]
    Return:
      data.frame with: 'Data'      - data used
                       'Block'     - current block that was used to fit RF
                       'train_dim' - dimension of the data used to fit the RF
                       'OOB_Acc'   - Metrics obtained from OOB / CV Testset
                       'Test_Acc'  -               - ' - 
                       'Test_F1'   -               - ' -
                       'Fold'      - which fold of the 5 fold CV
                       'Time'      - how long did it take to fit a RF on the block
                       'Fraction'  - fraction used to subset the single blocks!
      will be saved to 'docs/CV/gender/explorative_subsets/' the name of the file
      itself is a mixture of seed and fraction!
  "
  # [0] Check Arguments
  assert_vector(DFs_w_gender, min.len = 1, unique = TRUE)
  assert_int(seed_to_subset)
  assert_double(fraction_cnv, lower = 0, upper = 1)
  assert_double(fraction_rna, lower = 0, upper = 1)
  assert_double(fraction_mirna, lower = 0, upper = 1)
  assert_double(fraction_mutation, lower = 0, upper = 1)
  
  # [1] Empty DF - for all results of the Evaluation:
  eval_res <- data.frame("Data"    = character(), "train_dim" = numeric(),
                         "OOB_Acc" = numeric(),   "Test_Acc"  = numeric(),
                         "Test_F1" = numeric(),   "Fold"      = numeric(),
                         "Time"    = numeric(),   "Fraction"  = character(),
                         "subset_seed" = numeric(),
                         stringsAsFactors = F)
  
  # 1-1 Fixed Datapath, where we have all our raw Dataframes!
  data_path    <- "./data/external/Dr_Hornung/Data/ProcessedData/"
  
  # [2] Loop over all DFs w/ gender in their clinical variables!
  for (df in DFs_w_gender) {
    
    writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
    
    # 2-1 Load 'df' & only keep names of the omics blocks in 'df'!
    omics_blocks <- load(paste0(data_path, df, ".Rda"))
    omics_blocks <- omics_blocks[-which(omics_blocks %in% c("targetvar"))]
    
    # 2-1-1 Extract 'gender' as response and also remove it from the 'clin' block!
    resp <- clin["gender"]
    clin <- clin[-which(colnames(clin) == "gender")]
    
    # [2-2] Loop over all blocks in 'omics_blocks' of current 'df' and subset them
    for (curr_block in omics_blocks) {
      
      # If it is the 'clin' block keep only 'fraction_clin' features and assign
      # them again to the variable name it was originally saved in!
      # Analog for the rest!
      if (curr_block == 'clin') {
        curr_data <- eval(as.symbol(curr_block))
        set.seed(seed_to_subset)
        sampled_cols <- sample(seq_len(ncol(curr_data)), 
                               round(ncol(curr_data) * fraction_clin))
        assign(curr_block, curr_data[ ,sampled_cols])
      }
      
      if (curr_block == 'cnv') {
        curr_data <- eval(as.symbol(curr_block))
        set.seed(seed_to_subset)
        sampled_cols <- sample(seq_len(ncol(curr_data)), 
                               round(ncol(curr_data) * fraction_cnv))
        assign(curr_block, curr_data[ ,sampled_cols])
      }
      
      if (curr_block == 'rna') {
        curr_data <- eval(as.symbol(curr_block))
        set.seed(seed_to_subset)
        sampled_cols <- sample(seq_len(ncol(curr_data)), 
                               round(ncol(curr_data) * fraction_rna))
        assign(curr_block, curr_data[ ,sampled_cols])
        
        # add '_rna' to colnames, as cnv & rna have the same one else!
        colnames(rna) <- paste0(colnames(rna), "_rna")
      }
      
      if (curr_block == 'mutation') {
        curr_data <- eval(as.symbol(curr_block))
        set.seed(seed_to_subset)
        sampled_cols <- sample(seq_len(ncol(curr_data)), 
                               round(ncol(curr_data) * fraction_mutation))
        assign(curr_block, curr_data[ ,sampled_cols])
      }
      
      if (curr_block == 'mirna') {
        curr_data <- eval(as.symbol(curr_block))
        set.seed(seed_to_subset)
        sampled_cols <- sample(seq_len(ncol(curr_data)), 
                               round(ncol(curr_data) * fraction_mirna))
        assign(curr_block, curr_data[ ,sampled_cols])
      }
    }
    
    # 2-3 Bind all subsetted Blocks together!
    DF <- cbind(resp, 
                eval(as.symbol(omics_blocks[1])), 
                eval(as.symbol(omics_blocks[2])), 
                eval(as.symbol(omics_blocks[3])),
                eval(as.symbol(omics_blocks[4])),
                eval(as.symbol(omics_blocks[5])))
    
    # 2-3-1 Convert to DF, and name the response variable
    DF <- as.data.frame(DF)
    colnames(DF)[1] <- "resp"
    DF[,1] <- factor(DF[,1], levels = c(0, 1))
    
    # 2-4 Do 5 fold CV
    set.seed(12345)
    fold_ids <- createFolds(DF[,1], k = 5)
    
    for (i in 1:5) {
      
      print(paste("Fold:",as.character(i), "-----------------------------------"))
      
      # 2-2-4 split to Test & Train
      test      <- DF[fold_ids[[i]],]
      train_obs <- which(!(seq_len(nrow(DF)) %in% fold_ids[[i]]))
      train     <- DF[train_obs,]
      
      # 2-2-5 Print dimension of TrainingDF
      dim_train <- paste(dim(train), collapse = " x ")
      print(paste("Dimension of Train-DF:", dim_train))
      
      # 2-2-6 Fit a model [standard settings] on the train DF & take the time:
      start        <- Sys.time()
      curr_forrest <- rfsrc(as.formula(resp ~ .), data = train, seed = 1312,
                            ntree = 250)
      end          <- Sys.time()
      time_diff    <- difftime(end, start, units = 'mins')
      
      # 2-2-7 Evaluate the Model:
      # OOB:      Get OOB Accuracy!
      curr_OOB_ACC <- 1 - curr_forrest$err.rate[length(curr_forrest$err.rate)]
      
      # TESTSET:  Get the ErrorRates w/ TestSet
      prob_preds   <- predict(curr_forrest, test)$predicted
      pred_class   <- ifelse(prob_preds[,1] > 0.5, 0, 1)
      truth        <- as.factor(test[,1])
      pred_class   <- factor(pred_class, levels = levels(truth))
      res          <- cbind(pred_class, truth)
      test_acc     <- sum(res[,1] == res[,2]) / nrow(res)
      test_F1      <- caret::confusionMatrix(data      = pred_class, 
                                             reference = truth)$byClass["F1"]
      if (is.na(test_F1)) test_F1 <- 0 
      
      # 2-2-8 Add results of this block to the Results-DF!
      fraction_all <- paste("clin", fraction_clin,
                            "_cnv", fraction_cnv,
                            "_rna", fraction_rna,
                            "_mirna", fraction_mirna,
                            "_mutation", fraction_mutation, sep = "_")
      
      eval_res[nrow(eval_res) + 1, ] <- c(df, dim_train, curr_OOB_ACC, 
                                          test_acc, test_F1, i, 
                                          time_diff, fraction_all, seed_to_subset)
    }
  }
  res_name <- paste0("RF_joint_block_seed_", seed_to_subset,"__fraction_", fraction_all)
  
  write.csv2(eval_res, 
             paste0("./docs/cv_Res/gender/explorative_subsets/", res_name, ".csv"),
             row.names = FALSE) 
}

# Get Performances w/ single blocks as feature space!                        ----
# Single fully observed Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = 1234,
                          fraction = 1)

# Single 75% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = 1234,
                          fraction = 0.75)

# Single 50% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = 1234,
                          fraction = 0.5)
# Single 25% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = 1234,
                          fraction = 0.25)

# Single 15% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = 1234,
                          fraction = 0.15)

# Single 10% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = 1234,
                          fraction = 0.10)

# Single 5% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = 1234,
                          fraction = 0.05)

# Single 2.5% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = 1234,
                          fraction = 0.025)

# Single 1.25% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = 1234,
                          fraction = 0.0125)

# Single 0.5% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = 1234,
                          fraction = 0.05)

# ------------------------------------------ Different Seed!

# Single fully observed Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = 5678,
                          fraction = 1)

# Single 75% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = 5678,
                          fraction = 0.75)

# Single 50% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = 5678,
                          fraction = 0.5)
# Single 25% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = 5678,
                          fraction = 0.25)

# Single 15% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = 5678,
                          fraction = 0.15)

# Single 10% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = 5678,
                          fraction = 0.10)

# Single 5% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = 5678,
                          fraction = 0.05)

# Single 2.5% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = 5678,
                          fraction = 0.025)

# Single 1.25% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = 5678,
                          fraction = 0.0125)

# Single 0.5% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = 5678,
                          fraction = 0.05)

# Get Performances w/ joint subsetted blocks as feature spaces               ----
# eval_joint_block_subsets(DFs_w_gender = DFs_w_gender, seed_to_subset = 1234,
#                          fraction_cnv = 0.1, fraction_mirna = 0.5,
#                          fraction_mutation = 0.5, fraction_rna = 0.1,
#                          fraction_clin = 1)


