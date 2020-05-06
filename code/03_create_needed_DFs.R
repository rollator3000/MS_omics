" Script to create the subsets of the TCGA DFs needed for the evaluation!

  [1] Create Subsets of the original TCGA datasets!
      --> After exploration of the different subsets in 
          '02_explorative_performance.R', we create a final subset we will use 
          for the comparison study!
            For this we will subset each omics block:
              - 0.025 for cnv omics block!
              - 0.1   for mutationn omcis block!
              - 0.15  for rna omics block!
              - 0.5   for mirna omics block!
    
  [2] & [3] Get Performance on the subsetted DFs, when a regular DF is fit on 
            'Joint-' / 'Single-' Blocks!
          - [2] = JointBlock Performance
          - [3] = SingleBlock Performance
        
  [4] Get the dimensions of the subsetted data sets!
"
# Load Packages & define Functions!
setwd("C:/Users/kuche_000/Desktop/MS-Thesis")
library(randomForest)
library(randomForestSRC)
library(mlbench)
library(caret)

# [1] Subset the single Blocks in the original TCGA DFs                     ----
# 1-1 Define needed Variables
# 1-2-1 Names of the usable dataframes (gender in 'clin'-block & 4 omics blocks!)
DFs_w_gender <- c("BLCA", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC","LGG", 
                  "LUAD", "LUSC", "PAAD", "SARC", "SKCM", "STAD")

# 1-2-2 path to the data
data_path <- "./data/external/Dr_Hornung/"

# 1-2-3 seed for reproductibity [selected in 02_...R]
subset_seed <- 12345

# 1-2-4 set the fraction for the single blocks we want to keep [selected in 02_...R]
cnv_frac      <- 0.05
mutation_frac <- 0.1
rna_frac      <- 0.15
mirna_frac    <- 0.5

# 1-2 Loop over the DFs (w/ gender in 'clin') and create the subsets!
for (df in DFs_w_gender) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # 1-3 Load 'df' & only keep names of the omics blocks in 'df'!
  omics_blocks <- load(paste0(data_path, df, ".Rda"))
  omics_blocks <- omics_blocks[-which(omics_blocks %in% c("targetvar", "clin"))]
  
  # 1-4 Loop over all the blocks in 'omics_blocks'
  for (curr_block in omics_blocks) {
    
    writeLines(paste0("Current OmicsBlock: \n", curr_block))
    
    # 1-4-1 get the data to 'curr_data'
    curr_data <- eval(as.symbol(curr_block))
    
    # 1-5 subset the data based on the block type!
    # 1-5-1 MIRNA
    if (curr_block %in% c("mirna")) {
      set.seed(subset_seed)
      sampled_cols <- sample(seq_len(ncol(curr_data)), 
                             round(ncol(curr_data) * mirna_frac),
                             replace = F)
      
      assign(curr_block, curr_data[ ,sampled_cols])
    }
    
    # 1-5-2 MUTATION
    if (curr_block %in% c("mutation")) {
      set.seed(subset_seed)
      sampled_cols <- sample(seq_len(ncol(curr_data)), 
                             round(ncol(curr_data) * mutation_frac),
                             replace = F)
      
      assign(curr_block, curr_data[ ,sampled_cols])
    }
    
    # 1-5-3 CNV
    if (curr_block %in% c("cnv")) {
      set.seed(subset_seed)
      sampled_cols <- sample(seq_len(ncol(curr_data)), 
                             round(ncol(curr_data) * cnv_frac),
                             replace = F)
      
      assign(curr_block, curr_data[ ,sampled_cols])
    }
    
    # 1-5-5 RNA
    if (curr_block %in% c("rna")) {
      set.seed(subset_seed)
      sampled_cols <- sample(seq_len(ncol(curr_data)), 
                             round(ncol(curr_data) * rna_frac),
                             replace = F)
      
      assign(curr_block, curr_data[ ,sampled_cols])
    }
  }
  
  # 1-6 Save the subsetted DF
  save(mirna, mutation, rna, cnv, clin, 
       file = paste0("./data/processed/RH_subsetted_", subset_seed, "/", 
                     df, "_subset.RData"))
}

# [2] Get Test& OOB Performance on the subsetted DFs - JOINT BLOCKS          ----
# 2-1 Define file and foldernames!
DFs_w_gender        <- c("BLCA", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC",
                         "LGG", "LUAD", "LUSC", "PAAD", "SARC", "SKCM", "STAD")
DFs_w_gender_subset <- paste0(DFs_w_gender, "_subset.RData")
subset_folder       <- "./data/processed/RH_subsetted_12345/"

# 2-2 DF to save results of single block performances!
eval_res <- data.frame("Data"     = character(), "fold"     = numeric(),
                       "Dim"      = numeric(),   "OOB_Acc"  = numeric(),   
                       "Test_Acc" = numeric(),   "Test_F1"  = numeric(), 
                       "Time"     = numeric(), # min
                       stringsAsFactors = F)

# 2-3 Loop over the subsetted DFs, fit a RF & get the OOB- & Test-Error
for (df in DFs_w_gender_subset) {  
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # 2-3-1 Load current DF - all single blocks it consits of!
  omics_blocks <- load(paste0(subset_folder, df))
  
  # 2-3-2 Extract response 'gender' and remove it from clinical block!
  resp        <- as.factor(clin$gender)
  clin$gender <- NULL
  
  # 2-3-3 Bind all Blocks to a single DF, where the 1. column is the response!
  DF_all <- cbind(resp,
                  eval(as.symbol(omics_blocks[1])), 
                  eval(as.symbol(omics_blocks[2])), 
                  eval(as.symbol(omics_blocks[3])), 
                  eval(as.symbol(omics_blocks[4])), 
                  eval(as.symbol(omics_blocks[5])))
  
  # 2-3-4 Get the IDs for the different folds for the 5 fold CV!
  set.seed(12345)
  fold_ids <- createFolds(DF_all[,1], k = 5)
  
  # 2-3-5 Do 5 fold CV
  for (i in 1:5) {
    
    print(paste("Fold:", as.character(i), "----------------------------------"))
    
    # 2-3-5-1 Split to Test & Train
    test      <- DF_all[fold_ids[[i]],]
    train_obs <- which(!(seq_len(nrow(DF_all)) %in% fold_ids[[i]]))
    train     <- DF_all[train_obs,]
    
    # 2-3-5-2 Print dimension of TrainingDF
    dim_train <- paste(dim(train), collapse = " x ")
    
    # 2-3-5-3 Fit a model [standard settings] on the train DF & take the time:
    start        <- Sys.time()
    curr_forrest <- rfsrc(as.formula(resp ~ .), data = train, seed = 12345678,
                          ntree = 250)
    end          <- Sys.time()
    time_diff    <- difftime(end, start, units = 'mins')
    
    # 2-3-5-4 Evaluate the Model:
    # OOB:    Get OOB Accuracy!
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
    
    # 2-3-5-5 Add results of this block to the Results-DF!
    eval_res[nrow(eval_res) + 1, ] <- c(df, i, dim_train, curr_OOB_ACC, 
                                        test_acc, test_F1, time_diff)
    
  }
}

# 2-4 Save the result
write.csv2(eval_res, row.names = FALSE, 
           "./docs/CV_Res/gender/performance_final_subsets/joint_blocks_DFseed_12345.csv")

# [3] Get Test& OOB Performance on the subsetted DFs - SINGLE BLOCKS         ----
# 3-1 Define files and folder!
DFs_w_gender        <- c("BLCA", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC","LGG", 
                          "LUAD", "LUSC", "PAAD", "SARC", "SKCM", "STAD")
DFs_w_gender_subset <- paste0(DFs_w_gender, "_subset.RData")
subset_folder       <- "./data/processed/RH_subsetted_12345/"

# 3-2 Empty DF to store Results!
eval_res <- data.frame("Data"     = character(), 
                       "Block"    = numeric(),   "fold"     = numeric(),
                       "Dim"      = numeric(),   "OOB_Acc"  = numeric(),   
                       "Test_Acc" = numeric(),   "Test_F1"  = numeric(), 
                       "Time"     = numeric(), # min
                       stringsAsFactors = F)

# 3-3 Loop over all subsetted DFs, fit a RF & get the OOB- & Test-Error SINGLE BLOCK FEAS
for (df in DFs_w_gender_subset) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # 3-3-1 Load Data - with all its 5 blocks!
  omics_blocks <- load(paste0(subset_folder, df))
  
  # 3-3-2 Extract 'gender' response  & remove it from clincial block
  resp        <- as.factor(clin$gender)
  clin$gender <- NULL
  
  # 3-4 Loop over all the omcis blocks of current 'df'
  for (curr_block in omics_blocks) {
    
    writeLines(paste0("Current OmicsBlock: \n", curr_block))
    
    # 3-4-1 Bind response & the current omics block
    DF <- cbind(resp, eval(as.symbol(curr_block)))
    
    # 3-4-2 Convert it to DF and the respone to factor!
    DF     <- as.data.frame(DF)
    DF[,1] <- as.factor(DF[,1])
    levels(DF[,1]) <- levels(resp)
    
    # 3-4-3 Get the IDs for the 5 fold CV!
    set.seed(12345)
    fold_ids <- createFolds(DF[,1], k = 5)
    
    # 3-4-4 Start CV
    for (i in 1:5) {
      
      print(paste("Fold:",as.character(i), "---------------------------------"))
      
      # 3-4-4-1 Split to Test and Train
      test      <- DF[fold_ids[[i]],]
      train_obs <- which(!(seq_len(nrow(DF)) %in% fold_ids[[i]]))
      train     <- DF[train_obs,]
      
      # 3-4-4-2 Dimension of TrainingDF
      dim_train <- paste(dim(train), collapse = " x ")
      
      # 3-4-4-3 Fit a model [standard settings] on the train DF & take the time:
      start        <- Sys.time()
      curr_forrest <- rfsrc(as.formula(resp ~ .), data = train, seed = 12345678,
                            ntree = 250)
      end          <- Sys.time()
      time_diff    <- difftime(end, start, units = 'mins')
      
      # 3-4-4-4 Evaluate the Model:
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
      
      # 3-4-4-5 Add results of this block to the Results-DF!
      eval_res[nrow(eval_res) + 1, ] <- c(df, curr_block, i, dim_train,
                                          curr_OOB_ACC, test_acc, test_F1, time_diff)
    }
  }
}

# 3-5 Save the resuls
write.csv2(eval_res, row.names = FALSE, 
           "./docs/CV_Res/gender/performance_final_subsets/single_blocks_DFseed_12345.csv")
# [4] Get the amount of features in the reduced feature-blocks               ----
# 4-1 Define DFs we want to inspect!
DFs_w_gender <- c("BLCA", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC","LGG", 
                  "LUAD", "LUSC", "PAAD", "SARC", "SKCM", "STAD")

# 4-2 Loop over the possible feature-blocks in all different DFs, extract 
#     the dimension for each DF and print the summary of the dimensionality of
#     the current block over all 'DFs_w_gender'
clin_ <- c(); cnv_ <- c(); mirna_ <- c(); mutation_ <- c(); rna_ <- c()

for (df_ in DFs_w_gender) {
  
  # load 'df_' and its corresponding blocks
  df <- load(paste0("./data/processed/RH_subsetted_12345/", df_, "_subset.RData"))
  
  # Bind the amount of features to the corresponding vectors
  clin_     <- c(clin_, ncol(clin))
  cnv_      <- c(cnv_, ncol(cnv))
  mirna_    <- c(mirna_, ncol(mirna))
  mutation_ <- c(mutation_, ncol(mutation))
  rna_      <- c(rna_, ncol(rna))
}

# 4-3 Print the summaries for each
print(summary(clin_)) #  -1 as the response still inside!
print(summary(cnv_))
print(summary(mirna_))
print(summary(mutation_))
print(summary(rna_))