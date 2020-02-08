" Script to create DFs needed in the following Scripts!
  
  [1] Create an artififcal dataset w/ factors, characters & numeric features 
      - this was/ is needed for the adaption of the 'simpleRF' Package!
      
  [2] Create Subsets of the original datasets!
      --> After exploration of the different subsets in '02_...R', we create a 
          final subset we will use for the comparison study!
          For this we will subset each omics block:
            - 0.025 for cnv omics block!
            - 0.01  for mutationn omcis block!
            - 0.15  for rna omics block!
            - 0.5   for mirna omics block!
    
  [3] & [4] 
     Get Performance on the subsetted DFs, when a regular DF is fit 
     on 'Joint-' / 'Single-' Blocks!
        - [3] = JointBlock Performance
        - [4] = SingleBlock Performance
"
# Load Packages & define Functions!
setwd("C:/Users/kuche_000/Desktop/MS-Thesis")
library(randomForest)
library(randomForestSRC)
library(mlbench)
library(caret)

# [1] Create artifical DF for investigation                                  ----
# 1-1 Load raw 'iris' data
data("iris")

# 1-2 Add factor variable
# 1-2-1 Randomly assign factor variables
iris$factor <-  sample(c("good", "ok", "legend"), size = nrow(iris), 
                       replace = TRUE, prob = c(0.35, 0.4, 0.25))
#1-2-2 Set status to legend for a single species (--> high correlation)
iris$factor[iris$Species == "virginica"] <- "legend"

# 1-2-3 Conver to Factor
iris$factor <-  as.factor(iris$factor)

# 1-3 Add randomly a character Variable
iris$char <- sample(c("lol", "not_lol"), size = nrow(iris),
                    replace = T, prob = c(0.6, 0.4))
iris$char <- as.character(iris$char)

# 1-4 Save the DF for Investigation!   
write.csv2(iris, file = "./data/external/example_data/iris_example.csv",
           row.names = F)

# [2] Subset the single Blocks in the original DFs                           ----
# 2-1 Define needed Variables
# 2-2-1 Names of the usable dataframes (gender in 'clin'-block & 4 omics blocks!)
DFs_w_gender <- c("BLCA", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC","LGG", 
                  "LUAD", "LUSC", "PAAD", "SARC", "SKCM", "STAD")

# 2-2-2 path to the data
data_path <- "./data/external/Dr_Hornung/original_processed_data/"

# 2-2-3 seed for reproductibity [selected in 02_...R]
subset_seed <- 12345

# 2-2-4 set the fraction for the single blocks we want to keep [selected in 02_...R]
cnv_frac      <- 0.05
mutation_frac <- 0.1
rna_frac      <- 0.15
mirna_frac    <- 0.5

# 2-2 Loop over the DFs (w/ gender in 'clin') and create the subsets!
for (df in DFs_w_gender) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # 2-3 Load 'df' & only keep names of the omics blocks in 'df'!
  omics_blocks <- load(paste0(data_path, df, ".Rda"))
  omics_blocks <- omics_blocks[-which(omics_blocks %in% c("targetvar", "clin"))]
  
  # 2-4 Loop over all the blocks in 'omics_blocks'
  for (curr_block in omics_blocks) {
    
    writeLines(paste0("Current OmicsBlock: \n", curr_block))
    
    # 2-4-1 get the data to 'curr_data'
    curr_data <- eval(as.symbol(curr_block))
    
    # 2-5 subset the data based on the block type!
    # 2-5-1 MIRNA
    if (curr_block %in% c("mirna")) {
      set.seed(subset_seed)
      sampled_cols <- sample(seq_len(ncol(curr_data)), 
                             round(ncol(curr_data) * mirna_frac),
                             replace = F)
      
      assign(curr_block, curr_data[ ,sampled_cols])
    }
    
    # 2-5-2 MUTATION
    if (curr_block %in% c("mutation")) {
      set.seed(subset_seed)
      sampled_cols <- sample(seq_len(ncol(curr_data)), 
                             round(ncol(curr_data) * mutation_frac),
                             replace = F)
      
      assign(curr_block, curr_data[ ,sampled_cols])
    }
    
    # 2-5-3 CNV
    if (curr_block %in% c("cnv")) {
      set.seed(subset_seed)
      sampled_cols <- sample(seq_len(ncol(curr_data)), 
                             round(ncol(curr_data) * cnv_frac),
                             replace = F)
      
      assign(curr_block, curr_data[ ,sampled_cols])
    }
    
    # 2-5-5 RNA
    if (curr_block %in% c("rna")) {
      set.seed(subset_seed)
      sampled_cols <- sample(seq_len(ncol(curr_data)), 
                             round(ncol(curr_data) * rna_frac),
                             replace = F)
      
      assign(curr_block, curr_data[ ,sampled_cols])
    }
  }
  
  # 2-6 Save the subsetted DF
  save(mirna, mutation, rna, cnv, clin, 
       file = paste0("./data/external/Dr_Hornung/subsetted_", subset_seed, "/", 
                     df, "_subset.RData"))
    
}

# [3] Get Test& OOB Performance on the subsetted DFs - JOINT BLOCKS          ----
#     Define filenames!
DFs_w_gender        <- c("BLCA", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC",
                         "LGG", "LUAD", "LUSC", "PAAD", "SARC", "SKCM", "STAD")
DFs_w_gender_subset <- paste0(DFs_w_gender, "_subset.RData")

# DF to save results of single block performances!
eval_res <- data.frame("Data"     = character(), "fold"     = numeric(),
                       "Dim"      = numeric(),   "OOB_Acc"  = numeric(),   
                       "Test_Acc" = numeric(),   "Test_F1"  = numeric(), 
                       "Time"     = numeric(), # min
                       stringsAsFactors = F)

# Folder that holds all subsetted MUlti-Omics DFs!
subset_folder <- "./data/external/Dr_Hornung/subsetted_12345/"

# Loop over the subsetted DFs, fit a RF & get the OOB- & Test-Error
for (df in DFs_w_gender_subset) {  
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # Load current DF - all single blocks it consits of!
  omics_blocks <- load(paste0(subset_folder, df))
  
  # Extract response 'gender' and remove it from clinical block!
  resp        <- as.factor(clin$gender)
  clin$gender <- NULL
  
  # Bind all Blocks to a single DF, where the 1. column is the response!
  DF_all <- cbind(resp,
                  eval(as.symbol(omics_blocks[1])), 
                  eval(as.symbol(omics_blocks[2])), 
                  eval(as.symbol(omics_blocks[3])), 
                  eval(as.symbol(omics_blocks[4])), 
                  eval(as.symbol(omics_blocks[5])))
  
  # Get the IDs for the different folds for the 5 fold CV!
  set.seed(12345)
  fold_ids <- createFolds(DF_all[,1], k = 5)
  
  # 3-3-2-2 Do 5 fold CV
  for (i in 1:5) {
    
    print(paste("Fold:", as.character(i), "----------------------------------"))
    
    # 3-3-2-3 Split to Test & Train
    test      <- DF_all[fold_ids[[i]],]
    train_obs <- which(!(seq_len(nrow(DF_all)) %in% fold_ids[[i]]))
    train     <- DF_all[train_obs,]
    
    # 3-3-2-4 Print dimension of TrainingDF
    dim_train <- paste(dim(train), collapse = " x ")
    
    # 3-3-3 Fit a model [standard settings] on the train DF & take the time:
    start        <- Sys.time()
    curr_forrest <- rfsrc(as.formula(resp ~ .), data = train, seed = 12345678,
                          ntree = 250)
    end          <- Sys.time()
    time_diff    <- difftime(end, start, units = 'mins')
    
    # 3-3-4 Evaluate the Model:
    # OOB:      Get OOB Accuracy!
    curr_OOB_ACC <- 1 - curr_forrest$err.rate[nrow(curr_forrest$err.rate), 1]
    
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
    
    # 3-3-5 Add results of this block to the Results-DF!
    eval_res[nrow(eval_res) + 1, ] <- c(df, i, dim_train, curr_OOB_ACC, 
                                        test_acc, test_F1, time_diff)
    
  }
}

write.csv2(eval_res, row.names = FALSE, 
           "./docs/CV_Res/gender/performance_final_subsets/joint_blocks_DFseed_12345.csv")

# [4] Get Test& OOB Performance on the subsetted DFs - SINGLE BLOCKS         ----
#     Define files!
DFs_w_gender <- c("BLCA", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC","LGG", 
                  "LUAD", "LUSC", "PAAD", "SARC", "SKCM", "STAD")
DFs_w_gender_subset <- paste0(DFs_w_gender, "_subset.RData")

# Empty DF to store Results!
eval_res <- data.frame("Data"     = character(), 
                       "Block"    = numeric(),   "fold"     = numeric(),
                       "Dim"      = numeric(),   "OOB_Acc"  = numeric(),   
                       "Test_Acc" = numeric(),   "Test_F1"  = numeric(), 
                       "Time"     = numeric(), # min
                       stringsAsFactors = F)

# Set the folder we load the data from!
subset_folder <- "./data/external/Dr_Hornung/subsetted_12345/"

# Loop over all subsetted DFs, fit a RF & get the OOB- & Test-Error SINGLE BLOCK FEAS
for (df in DFs_w_gender_subset) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # Load Data - with all its 5 blocks!
  omics_blocks <- load(paste0(subset_folder, df))
  
  # Extract 'gender' response  & remove it from clincial block
  resp        <- as.factor(clin$gender)
  clin$gender <- NULL
  
  # Loop over all the omcis blocks of current 'df'
  for (curr_block in omics_blocks) {
    
    writeLines(paste0("Current OmicsBlock: \n", curr_block))
    
    # Bind response & the current omics block
    DF <- cbind(resp, eval(as.symbol(curr_block)))
    
    # Convert it to DF and the respone to factor!
    DF     <- as.data.frame(DF)
    DF[,1] <- as.factor(DF[,1])
    levels(DF[,1]) <- levels(resp)
    
    # Get the IDs for the 5 fold CV!
    set.seed(12345)
    fold_ids <- createFolds(DF[,1], k = 5)
    
    # Start CV
    for (i in 1:5) {
      
      print(paste("Fold:",as.character(i), "---------------------------------"))
      
      # Split to Test and Train
      test      <- DF[fold_ids[[i]],]
      train_obs <- which(!(seq_len(nrow(DF)) %in% fold_ids[[i]]))
      train     <- DF[train_obs,]
      
      # Dimension of TrainingDF
      dim_train <- paste(dim(train), collapse = " x ")
      
      # Fit a model [standard settings] on the train DF & take the time:
      start        <- Sys.time()
      curr_forrest <- rfsrc(as.formula(resp ~ .), data = train, seed = 12345678,
                            ntree = 250)
      end          <- Sys.time()
      time_diff    <- difftime(end, start, units = 'mins')
      
      # Evaluate the Model:
      # OOB:    Get OOB Accuracy of each of the 500 trees and average it!
      curr_OOB_ACC <- 1 - curr_forrest$err.rate[nrow(curr_forrest$err.rate), 1]
      
      # TESTSET:  Get the ErrorRates w/ TestSet
      prob_preds <- predict(curr_forrest, test)$predicted[,1]
      pred_class <- ifelse(prob_preds > 0.5, 0, 1)
      truth      <- as.factor(test$resp)
      pred_class <- factor(pred_class, levels = levels(truth))
      res        <- cbind(pred_class, truth)
      test_acc   <- sum(res[,1] == res[,2]) / nrow(res)
      test_F1    <- caret::confusionMatrix(data      = pred_class, 
                                           reference = truth)$byClass["F1"]
      if (is.na(test_F1)) test_F1 <- 0
      
      # Add results of this block to the Results-DF!
      eval_res[nrow(eval_res) + 1, ] <- c(df, curr_block, i, dim_train,
                                          curr_OOB_ACC, test_acc, test_F1, time_diff)
    }
  }
}

write.csv2(eval_res, row.names = FALSE, 
           "./docs/CV_Res/gender/performance_final_subsets/single_blocks_DFseed_12345.csv"))

