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
library(randomForestSRC)

# Names of the usable dataframes (w/ gender in 'clin'-block & 4 omics blocks!)
DFs_w_gender <- c("BLCA", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC","LGG", 
                  "LUAD", "LUSC", "PAAD", "SARC", "SKCM", "STAD")
data_path    <- "./data/external/Dr_Hornung/Data/ProcessedData/"

# Single fully observed Omics-Blocks as features to RF                      ----
" Use single omics Blocks and get a performance of a RF how it performs!
  Get performance, runtime, etc.
     --> local machine can maximally handle 12.5k feas for a tree
     --> not all blocks could be used as full block!
"
# [1] Empty DF - for all results of the Evaluation
eval_res <- data.frame("Data"     = character(), "Block"    = numeric(), 
                       "Dim"      = numeric(),   "OOB_Acc"  = numeric(),   
                       "Test_Acc" = numeric(),   "Test_F1"  = numeric(), 
                       "Time"     = numeric(), # min
                       stringsAsFactors = F)

# [2] Loop over all DFs w/ gender in their clinical variables!
#     For Blocks w/ more than 12.5k feas, we skip them, as it crashes!
for (df in DFs_w_gender) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # 2-1 Load 'df' & only keep names of the omics blocks in 'df'!
  omics_blocks <- load(paste0(data_path, df, ".Rda"))
  omics_blocks <- omics_blocks[-which(omics_blocks %in% c("targetvar"))]
  
  # 2-1-1 Extract 'gender' as response and also remove it from the 'clin' block!
  resp <- clin["gender"]
  clin <- clin[-which(colnames(clin) == "gender")]
  
  # [2-2] Loop over all blocks in 'omics_blocks' of current 'df'
  for (curr_block in  omics_blocks) {
    
    writeLines(paste0("Current OmicsBlock: \n", curr_block))
    
    # 2-2-1 Bind gender from the clin block & 'curr_block' (=omics) 
    DF <- cbind(resp, eval(as.symbol(curr_block)))
    
    # 2-2-2 Convert to DF, and name the response variable
    DF <- as.data.frame(DF); colnames(DF)[1] <- "gender"
    
    # 2-2-3 Split to train and test set [80:20] w/ seed for reproducibility!
    set.seed(12345)
    train_ind <- sample(x = seq_len(nrow(DF)), size = round(nrow(DF) * 0.8))
    test_ind  <- which(!(1:nrow(DF) %in% train_ind))
    train  <- DF[train_ind,]; test <- DF[test_ind,]
    
    # 2-2-3-1 Print dimension of TrainingDF
    dim_train <- paste(dim(train), collapse = " x ")
    print(paste("Dimension of Train-DF:", dim_train))
    
    # 2-2-4 If the block has more than 12500 feas, stop it, as overflow error,
    #       when fitting RF to it!
    if (dim(train)[2] > 12500) {
      print("too high dimensional! Will crash the computer!")
      next
    }
    
    # 2-2-5 Fit a model [standard settings] on the train DF & take the time:
    start        <- Sys.time()
    curr_forrest <- rfsrc(as.factor(gender) ~ ., data = train, seed = 123,
                          ntree = 250)
    end          <- Sys.time()
    time_diff    <- difftime(end, start, units = 'mins')
    
    # 2-2-6 Evaluate the Model:
    # OOB:      Get OOB Accuracy!
    curr_OOB_ACC <- 1 - curr_forrest$err.rate[length(curr_forrest$err.rate)]
    
    # TESTSET:  Get the ErrorRates w/ TestSet
    prob_preds   <- predict(curr_forrest, test)$predicted
    pred_class   <- ifelse(prob_preds > 0.5, 1, 0)
    truth        <- as.factor(test$gender)
    pred_class   <- factor(pred_class, levels = levels(truth))
    res          <- cbind(pred_class, truth)
    test_acc     <- sum(res[,1] == res[,2]) / nrow(res)
    test_F1      <- caret::confusionMatrix(data      = pred_class, 
                                           reference = truth)$byClass["F1"]
    if (is.na(test_F1)) test_F1 <- 0 
    
    # 2-2-7 Add results of this block to the Results-DF!
    eval_res[nrow(eval_res) + 1, ] <- c(df, curr_block, dim_train, curr_OOB_ACC, 
                                        test_acc, test_F1, time_diff)
    
  }
}

# Save the results!
write.csv2(eval_res, "./docs/CV_Res/gender/explorative_subsets/RF_whole_single_blocks.csv",
           row.names = FALSE)
# Single 10% subsetted Omics-Blocks as features to RF                       ----
"Use single Blocks w/ 10% of the original features & get performance, runtime...
 There shouldn't be any more computional issues now!
"
# [1] Empty DF - for all results of the Evaluation:
eval_res_sub <- data.frame("Data"     = character(), "Block"    = numeric(), 
                           "Dim"      = numeric(),   "OOB_Acc"  = numeric(),   
                           "Test_Acc" = numeric(),   "Test_F1"  = numeric(), 
                           "Time"     = numeric(), # min
                           stringsAsFactors = F)

# [2] Loop over all DFs w/ gender in their clinical variables!
for (df in DFs_w_gender) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # [2-1] Load 'df' & only keep names of the omics blocks in 'df'!
  omics_blocks <- load(paste0(data_path, df, ".Rda"))
  omics_blocks <- omics_blocks[-which(omics_blocks %in% c("targetvar"))]
  
  # 2-1-1 Extract 'gender' as response and also remove it from the 'clin' block!
  resp <- clin["gender"]
  clin <- clin[-which(colnames(clin) == "gender")]
  
  # [2-2] Loop over all blocks in 'omics_blocks' of current 'df'
  for (curr_block in omics_blocks) {
    
    writeLines(paste0("Current OmicsBlock: \n", curr_block))
    
    # 2-2-1 If the block is not clinical only keep 10% of the original features!
    curr_data <- eval(as.symbol(curr_block))
    
    if (curr_block != "clin") {
      set.seed(12345)
      sampled_cols <- sample(seq_len(ncol(curr_data)), round(ncol(curr_data) * 0.1))
    } else {
      sampled_cols <- c(1:ncol(curr_data))
    }
    
    # 2-2-2 Bind gender from the clin block & the (subsetted) curr block 
    DF <- cbind(resp, curr_data[ ,sampled_cols])
    
    # [2-2-2] Convert to DF, and name the response variable
    DF <- as.data.frame(DF); colnames(DF)[1] <- "gender"
    
    # [2-2-3] Split to train and test set [80:20]
    set.seed(12345)
    train_ind <- sample(x = seq_len(nrow(DF)), size = round(nrow(DF) * 0.8))
    test_ind  <- which(!(1:nrow(DF) %in% train_ind))
    train  <- DF[train_ind,]; test <- DF[test_ind,]
    
    # [2-2-3-1] Print dimension of TrainingDF
    dim_train <- paste(dim(train), collapse = " x ")
    print(paste("Dimension of Train-DF:", dim_train))
    
    # [2-2-4] If the block has more than 12500 feas, stop it, as overflow error,
    #         when fitting RF to it!
    if (dim(train)[2] > 12500) {
      print("too high dimensional! Will crash the computer!")
      next
    }
    
    # [2-2-5] Fit a model [standard settings] on the train DF & take the time:
    start        <- Sys.time()
    curr_forrest <- rfsrc(as.factor(gender) ~ ., data = train, seed = 123,
                          ntree = 250)
    end          <- Sys.time()
    time_diff    <- difftime(end, start, units = 'mins')
    
    # [2-2-6] Evaluate the Model:
    # OOB:    Get OOB Accuracy!
    curr_OOB_ACC <- 1 - curr_forrest$err.rate[length(curr_forrest$err.rate)]
    
    # TESTSET:  Get the ErrorRates w/ TestSet
    prob_preds   <- predict(curr_forrest, test)$predicted
    pred_class   <- ifelse(prob_preds > 0.5, 1, 0)
    truth        <- as.factor(test$gender)
    pred_class   <- factor(pred_class, levels = levels(truth))
    res          <- cbind(pred_class, truth)
    test_acc     <- sum(res[,1] == res[,2]) / nrow(res)
    test_F1      <- caret::confusionMatrix(data      = pred_class, 
                                           reference = truth)$byClass["F1"]
    if (is.na(test_F1)) test_F1 <- 0 
    
    # [2-2-7] Add results of this block to the Results-DF!
    eval_res_sub[nrow(eval_res_sub) + 1, ] <- c(df, curr_block, dim_train, 
                                                curr_OOB_ACC, test_acc, 
                                                test_F1, time_diff)
  }
}

# [3] Save the results!
write.csv2(eval_res_sub, "./docs/CV_Res/gender/explorative_subsets/RF_10%_of_single_blocks.csv",
           row.names = FALSE)

# 10% subsets of all blocks as feature to RF                                ----
"Remove 90% of all omics blocks features and paste them together to learn a 
 model on the joint data"

# [1] Empty DF - for all results of the Evaluation:
eval_res <- data.frame("Data"     = character(), "Dim"      = numeric(),   
                       "OOB_Acc"  = numeric(),   "Test_Acc" = numeric(),   
                       "Test_F1"  = numeric(),   "Time"     = numeric(), # min
                       stringsAsFactors = F)
# [2] Loop over all DFs w/ gender in their clinical variables!
for (df in DFs_w_gender) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # 2-1 Load 'df' & only keep names of the omics blocks in 'df'!
  omics_blocks <- load(paste0(data_path, df, ".Rda"))
  omics_blocks <- omics_blocks[-which(omics_blocks %in% c("targetvar"))]
  
  # 2-1-1 Extract 'gender' as response and also remove it from the 'clin' block!
  resp <- clin["gender"]
  clin <- clin[-which(colnames(clin) == "gender")]
  
  # [2-2] Loop over all blocks in 'omics_blocks' of current 'df' to prune them!
  for (curr_block in omics_blocks) {
    
    writeLines(paste0("Remove 90% of the features of OmicsBlock:\n", curr_block))
    
    # 2-2-1 Load data from current block and remove 90% feas, except for 'clin'
    curr_data <- eval(as.symbol(curr_block))
    
    if (curr_block != "clin") {
      set.seed(12345)
      sampled_cols <- sample(seq_len(ncol(curr_data)), round(ncol(curr_data) * 0.1))
    } else {
      sampled_cols <- c(1:ncol(curr_data))
    }
    
    # 2-2-2 Only keep the sampled columns in the block 
    #       [overwrite current variable w/ reduced dimension var!]
    assign(curr_block, curr_data[ ,sampled_cols])
  }
  
  # [2-3] Bind gender from the clin block & all blocks!
  # 2-3-1 As rna & cnv share colnames [which leads to an error when fitting RF]
  #       we rename the colnames of the rna DF
  colnames(rna) <- paste0(colnames(rna), "_rna")
  
  # [2-3-2] Paste the subsetted DFs together!
  DF <- cbind(resp, 
              eval(as.symbol(omics_blocks[1])), 
              eval(as.symbol(omics_blocks[2])), 
              eval(as.symbol(omics_blocks[3])),
              eval(as.symbol(omics_blocks[4])),
              eval(as.symbol(omics_blocks[5])))
  
  # [2-4] Convert to DF, and name the response variable
  DF <- as.data.frame(DF); colnames(DF)[1] <- "gender"
  
  # [2-5] Split to train and test set [80:20] w/ seed for reproducibility!
  set.seed(12345)
  train_ind <- sample(x = seq_len(nrow(DF)), size = round(nrow(DF) * 0.8))
  test_ind  <- which(!(1:nrow(DF) %in% train_ind))
  train  <- DF[train_ind,]; test <- DF[test_ind,]
  
  # 2-5-1 Print dimension of TrainingDF
  dim_train <- paste(dim(train), collapse = " x ")
  print(paste("Dimension of Train-DF:", dim_train))
  
  # 2-5-2 If the block has more than 12500 feas, stop it, as overflow error,
  #       when fitting RF to it!
  if (dim(train)[2] > 12500) {
    print("too high dimensional! Will crash the computer!")
    next
  }
  
  # [2-6] Fit a model [standard settings] on the train DF & take the time:
  start        <- Sys.time()
  curr_forrest <- rfsrc(as.factor(gender) ~ ., data = train, seed = 123,
                        ntree = 250)
  end          <- Sys.time()
  time_diff    <- difftime(end, start, units = 'mins')
  
  # [2-7] Evaluate the Model:
  # OOB:      Get OOB Accuracy!
  curr_OOB_ACC <- 1 - curr_forrest$err.rate[length(curr_forrest$err.rate)]
  
  # TESTSET:  Get the ErrorRates w/ TestSet
  prob_preds   <- predict(curr_forrest, test)$predicted
  pred_class   <- ifelse(prob_preds > 0.5, 1, 0)
  truth        <- as.factor(test$gender)
  pred_class   <- factor(pred_class, levels = levels(truth))
  res          <- cbind(pred_class, truth)
  test_acc     <- sum(res[,1] == res[,2]) / nrow(res)
  test_F1      <- caret::confusionMatrix(data      = pred_class, 
                                         reference = truth)$byClass["F1"]
  if (is.na(test_F1)) test_F1 <- 0 
  
  # [2-8] Add results of this block to the Results-DF!
  eval_res[nrow(eval_res) + 1, ] <- c(df, dim_train, curr_OOB_ACC, test_acc, 
                                      test_F1, time_diff)
}

# [3] Save the results!
write.csv2(eval_res, "./docs/cv_Res/gender/explorative_subsets/RF_10%_all_blocks_joint.csv",
           row.names = FALSE)

# 5% subsets of 'rna'/'cnv' as features to RF                               ----
" For the blocks 'rna' / 'cnv' we use 5% subsets and train 
  a RF on these [seperat] 
"
# [1] Empty DF - for all results of the Evaluation:
eval_res <- data.frame("Data"     = character(), "Block"    = numeric(), 
                       "Dim"      = numeric(),   "OOB_Acc"  = numeric(),   
                       "Test_Acc" = numeric(),   "Test_F1"  = numeric(), 
                       "Time"     = numeric(), # min
                       stringsAsFactors = F)

# [2] Loop over all DFs w/ gender in their clinical variables!
for (df in DFs_w_gender) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # [2-1] Load 'df' & only keep 'rna' & 'cnv' + the gender we use as response
  load(paste0(data_path, df, ".Rda"))
  omics_blocks <- c("rna", "cnv")
  resp         <- clin["gender"]
  
  # [2-2] Loop over all blocks in 'omics_blocks' of current 'df'
  for (curr_block in omics_blocks) {
    
    writeLines(paste0("Current OmicsBlock: \n", curr_block))
    
    # 2-2-1 Subset the current Block and  only keep 5% of the original feas
    curr_data    <- eval(as.symbol(curr_block))
    set.seed(12345)
    sampled_cols <- sample(ncol(curr_data), round(ncol(curr_data) * 0.05))
    
    # 2-2-2 Bind gender from the clin block & 'curr_block' (=omics), 
    #       convert it to a DF and name the response column
    DF <- cbind(resp, curr_data[ ,sampled_cols])
    DF <- as.data.frame(DF)
    colnames(DF)[1] <- "gender"
    
    # 2-2-3 Split to train and test set [80:20] w/ seed for reproducibility!
    set.seed(12345)
    train_ind <- sample(x = seq_len(nrow(DF)), size = round(nrow(DF) * 0.8))
    test_ind  <- which(!(1:nrow(DF) %in% train_ind))
    train  <- DF[train_ind,]; test <- DF[test_ind,]
    
    # 2-2-3-1 Print dimension of TrainingDF
    dim_train <- paste(dim(train), collapse = " x ")
    print(paste("Dimension of Train-DF:", dim_train))
    
    # 2-2-4 Fit a model [standard settings] on the train DF & take the time:
    start        <- Sys.time()
    curr_forrest <- rfsrc(as.factor(gender) ~ ., data = train, seed = 123,
                          ntree = 250)
    end          <- Sys.time()
    time_diff    <- difftime(end, start, units = 'mins')
    
    # 2-2-6 Evaluate the Model:
    # OOB:      Get OOB Accuracy of each of the 500 trees and average it!
    curr_OOB_ACC <- 1 - curr_forrest$err.rate[length(curr_forrest$err.rate)]
    
    # TESTSET:  Get the ErrorRates w/ TestSet
    prob_preds   <- predict(curr_forrest, test)$predicted
    pred_class   <- ifelse(prob_preds > 0.5, 1, 0)
    truth        <- as.factor(test$gender)
    pred_class   <- factor(pred_class, levels = levels(truth))
    res          <- cbind(pred_class, truth)
    test_acc     <- sum(res[,1] == res[,2]) / nrow(res)
    test_F1      <- caret::confusionMatrix(data      = pred_class, 
                                           reference = truth)$byClass["F1"]
    if (is.na(test_F1)) test_F1 <- 0 
    
    # [2-2-6] Add results of this block to the Results-DF!
    eval_res[nrow(eval_res) + 1, ] <- c(df, curr_block, dim_train, 
                                        curr_OOB_ACC, test_acc, test_F1, time_diff)
  }
}

# [3] Save the results!
write.csv2(eval_res, "./docs/CV_Res/gender/explorative_subsets/RF_5%_of_rna_&_cnv_single_block.csv",
           row.names = FALSE)

# Use 10% mirna & mutation + 5% cnv & rna all together                      ----
" Subset mirna & mutation blocks by 10% and use 5% subsets of the 'cnv' & 'rna'
  blocks, as featurespace to predcit the gender!
"
eval_res <- data.frame("Data"     = character(), "Dim"      = numeric(),   
                       "OOB_Acc"  = numeric(),   "Test_Acc" = numeric(),   
                       "Test_F1"  = numeric(),   "Time"     = numeric(), # min
                       stringsAsFactors = F)

# [2] Loop over all DFs w/ gender in their clinical variables!
for (df in DFs_w_gender) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # [2-1] Load 'df' & only keep names of the omics blocks in 'df'!
  omics_blocks <- load(paste0(data_path, df, ".Rda"))
  omics_10     <- c("mirna", "mutation")
  omics_5      <- c("rna", "cnv")
  omics_all    <- c("clin")
  
  resp <- clin$gender
  clin$gender <- NULL
  
  # [2-2] Loop over all blocks in 'omics_blocks' of current 'df' to prune them!
  # [2-2-1] Blocks that get reduced by 10%
  for (curr_block in omics_10) {
    
    writeLines(paste0("Remove 90% of the features of OmicsBlock:\n", curr_block))
    
    # 2-2-1-1 Subset the current Block and only keep 10% of the original feas
    curr_data <- eval(as.symbol(curr_block))
    set.seed(12345)
    sampled_cols <- sample(seq_len(ncol(curr_data)), round(ncol(curr_data) * 0.1))
    
    # 2-2-1-2 Only keep the sampled columns in the block 
    #        [overwrite current variable w/ reduced dimension var!]
    assign(curr_block, curr_data[ ,sampled_cols])
  }
  
  # [2-2-2] Blocks that get reduced by 5%
  for (curr_block in omics_5) {
    
    writeLines(paste0("Remove 95% of the features of OmicsBlock:\n", curr_block))
    
    # 2-2-1-1 Subset the current Block and only keep 10% of the original feas
    curr_data <- eval(as.symbol(curr_block))
    set.seed(12345)
    sampled_cols <- sample(seq_len(ncol(curr_data)), round(ncol(curr_data) * 0.05))
    
    # 2-2-1-2 Only keep the sampled columns in the block 
    #        [overwrite current variable w/ reduced dimension var!]
    assign(curr_block, curr_data[ ,sampled_cols])
  }
  
  # [2-3] Bind gender from the clin block & all blocks!
  # [2-3-1] As rna & cnv share colnames [which leads to an error when fitting RF]
  #         we rename the colnames of the rna DF
  colnames(rna) <- paste0(colnames(rna), "_rna")
  
  # [2-3-2] Paste the subsetted DFs together!
  DF <- cbind(resp, 
              eval(as.symbol(omics_10[1])), 
              eval(as.symbol(omics_10[2])), 
              eval(as.symbol(omics_5[1])),
              eval(as.symbol(omics_5[2])),
              clin)
  
  # [2-4] Convert to DF, and name the response variable
  DF <- as.data.frame(DF); colnames(DF)[1] <- "gender"
  
  # [2-5] Split to train and test set [80:20] w/ seed for reproducibility!
  set.seed(12345)
  train_ind <- sample(x = seq_len(nrow(DF)), size = round(nrow(DF) * 0.8))
  test_ind  <- which(!(1:nrow(DF) %in% train_ind))
  train  <- DF[train_ind,]; test <- DF[test_ind,]
  
  # [2-5-1] Print dimension of TrainingDF
  dim_train <- paste(dim(train), collapse = " x ")
  print(paste("Dimension of Train-DF:", dim_train))
  
  # [2-6] Fit a model [standard settings] on the train DF & take the time:
  start        <- Sys.time()
  curr_forrest <- rfsrc(as.factor(gender) ~ ., data = train, seed = 123,
                        ntree = 250)
  end          <- Sys.time()
  time_diff    <- difftime(end, start, units = 'mins')
  
  # [2-7] Evaluate the Model:
  # OOB:      Get OOB Accuracy of each of the 500 trees and average it!
  curr_OOB_ACC <- 1 - curr_forrest$err.rate[length(curr_forrest$err.rate)]
  
  # TESTSET:  Get the ErrorRates w/ TestSet
  prob_preds   <- predict(curr_forrest, test)$predicted
  pred_class   <- ifelse(prob_preds > 0.5, 1, 0)
  truth        <- as.factor(test$gender)
  pred_class   <- factor(pred_class, levels = levels(truth))
  res          <- cbind(pred_class, truth)
  test_acc     <- sum(res[,1] == res[,2]) / nrow(res)
  test_F1      <- caret::confusionMatrix(data      = pred_class, 
                                         reference = truth)$byClass["F1"]
  if (is.na(test_F1)) test_F1 <- 0 
  
  # [2-8] Add results of this block to the Results-DF!
  eval_res[nrow(eval_res) + 1, ] <- c(df, dim_train, 
                                      curr_OOB_ACC, test_acc, 
                                      test_F1, time_diff)
}

# [3] Save the results!
write.csv2(eval_res, "./docs/CV_Res/gender/explorative_subsets/RF_10%_mirna_&_mutation_5%_rna_&_cnv_all_blocks_joint.csv",
           row.names = FALSE)

# 2.5% subsets of 'rna'&'cnv' as features to RF                             ----
"Use single 'rna'/'cnv' Blocks w/ 2.5% of the original features 
 & get performance, runtime...
 There should be no issues with the computer now!
"
# [1] Empty DF - for all results of the Evaluation:
eval_res <- data.frame("Data"     = character(), "Block"    = numeric(), 
                       "Dim"      = numeric(),   "OOB_Acc"  = numeric(),   
                       "Test_Acc" = numeric(),   "Test_F1"  = numeric(), 
                       "Time"     = numeric(), # min
                       stringsAsFactors = F)

# [2] Loop over all DFs w/ gender in their clinical variables!
for (df in DFs_w_gender) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # [2-1] Load 'df' & only keep 'rna' & 'cnv' + the gender we use as response
  load(paste0(data_path, df, ".Rda"))
  omics_blocks <- c("rna", "cnv")
  resp         <- clin["gender"]
  
  # [2-2] Loop over all blocks in 'omics_blocks' of current 'df'
  for (curr_block in omics_blocks) {
    
    writeLines(paste0("Current OmicsBlock: \n", curr_block))
    
    # 2-2-1 Subset the current Block and  only keep 5% of the original feas
    curr_data    <- eval(as.symbol(curr_block))
    set.seed(12345)
    sampled_cols <- sample(ncol(curr_data), round(ncol(curr_data) * 0.025))
    
    # 2-2-2 Bind gender from the clin block & 'curr_block' (=omics), 
    #       convert it to a DF and name the response column
    DF <- cbind(resp, curr_data[ ,sampled_cols])
    DF <- as.data.frame(DF)
    colnames(DF)[1] <- "gender"
    
    # 2-2-3 Split to train and test set [80:20] w/ seed for reproducibility!
    set.seed(12345)
    train_ind <- sample(x = seq_len(nrow(DF)), size = round(nrow(DF) * 0.8))
    test_ind  <- which(!(1:nrow(DF) %in% train_ind))
    train  <- DF[train_ind,]; test <- DF[test_ind,]
    
    # 2-2-3-1 Print dimension of TrainingDF
    dim_train <- paste(dim(train), collapse = " x ")
    print(paste("Dimension of Train-DF:", dim_train))
    
    # 2-2-4 Fit a model [standard settings] on the train DF & take the time:
    start        <- Sys.time()
    curr_forrest <- rfsrc(as.factor(gender) ~ ., data = train, seed = 123,
                          ntree = 250)
    end          <- Sys.time()
    time_diff    <- difftime(end, start, units = 'mins')
    
    # 2-2-6 Evaluate the Model:
    # OOB:      Get OOB Accuracy of each of the 500 trees and average it!
    curr_OOB_ACC <- 1 - curr_forrest$err.rate[length(curr_forrest$err.rate)]
    
    # TESTSET:  Get the ErrorRates w/ TestSet
    prob_preds   <- predict(curr_forrest, test)$predicted
    pred_class   <- ifelse(prob_preds > 0.5, 1, 0)
    truth        <- as.factor(test$gender)
    pred_class   <- factor(pred_class, levels = levels(truth))
    res          <- cbind(pred_class, truth)
    test_acc     <- sum(res[,1] == res[,2]) / nrow(res)
    test_F1      <- caret::confusionMatrix(data      = pred_class, 
                                           reference = truth)$byClass["F1"]
    if (is.na(test_F1)) test_F1 <- 0 
    
    # [2-2-6] Add results of this block to the Results-DF!
    eval_res[nrow(eval_res) + 1, ] <- c(df, curr_block, dim_train, 
                                        curr_OOB_ACC, test_acc, test_F1, time_diff)
  }
}

# [3] Save the results!
write.csv2(eval_res, "./docs/CV_Res/gender/explorative_subsets/RF_2,5%_of_rna_&_cnv_single_block.csv",
           row.names = FALSE)

# 10% mirna & mutation + 2.5% cnv & rna as joint features to RF             ----
" Subset mirna & mutation blocks by 10% and use 2.5% subsets of the 'cnv' & 'rna'
  blocks, as featurespace to predcit the gender!
"
eval_res <- data.frame("Data"     = character(), "Dim"      = numeric(),   
                       "OOB_Acc"  = numeric(),   "Test_Acc" = numeric(),   
                       "Test_F1"  = numeric(),   "Time"     = numeric(), # min
                       stringsAsFactors = F)

# [2] Loop over all DFs w/ gender in their clinical variables!
for (df in DFs_w_gender) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # [2-1] Load 'df' & only keep names of the omics blocks in 'df'!
  omics_blocks <- load(paste0(data_path, df, ".Rda"))
  omics_10     <- c("mirna", "mutation")
  omics_5      <- c("rna", "cnv")
  omics_all    <- c("clin")
  
  resp <- clin$gender
  clin$gender <- NULL
  
  # [2-2] Loop over all blocks in 'omics_blocks' of current 'df' to prune them!
  # [2-2-1] Blocks that get reduced by 10%
  for (curr_block in omics_10) {
    
    writeLines(paste0("Remove 95% of the features of OmicsBlock:\n", curr_block))
    
    # 2-2-1-1 Subset the current Block and only keep 10% of the original feas
    curr_data <- eval(as.symbol(curr_block))
    set.seed(12345)
    sampled_cols <- sample(seq_len(ncol(curr_data)), round(ncol(curr_data) * 0.1))
    
    # 2-2-1-2 Only keep the sampled columns in the block 
    #        [overwrite current variable w/ reduced dimension var!]
    assign(curr_block, curr_data[ ,sampled_cols])
  }
  
  # [2-2-2] Blocks that get reduced by 5%
  for (curr_block in omics_5) {
    
    writeLines(paste0("Remove 97.5% of the features of OmicsBlock:\n", curr_block))
    
    # 2-2-1-1 Subset the current Block and only keep 10% of the original feas
    curr_data <- eval(as.symbol(curr_block))
    set.seed(12345)
    sampled_cols <- sample(seq_len(ncol(curr_data)), round(ncol(curr_data) * 0.025))
    
    # 2-2-1-2 Only keep the sampled columns in the block 
    #        [overwrite current variable w/ reduced dimension var!]
    assign(curr_block, curr_data[ ,sampled_cols])
  }
  
  # [2-3] Bind gender from the clin block & all blocks!
  # [2-3-1] As rna & cnv share colnames [which leads to an error when fitting RF]
  #         we rename the colnames of the rna DF
  colnames(rna) <- paste0(colnames(rna), "_rna")
  
  # [2-3-2] Paste the subsetted DFs together!
  DF <- cbind(resp, 
              eval(as.symbol(omics_10[1])), 
              eval(as.symbol(omics_10[2])), 
              eval(as.symbol(omics_5[1])),
              eval(as.symbol(omics_5[2])),
              clin)
  
  # [2-4] Convert to DF, and name the response variable
  DF <- as.data.frame(DF); colnames(DF)[1] <- "gender"
  
  # [2-5] Split to train and test set [80:20] w/ seed for reproducibility!
  set.seed(12345)
  train_ind <- sample(x = seq_len(nrow(DF)), size = round(nrow(DF) * 0.8))
  test_ind  <- which(!(1:nrow(DF) %in% train_ind))
  train  <- DF[train_ind,]; test <- DF[test_ind,]
  
  # [2-5-1] Print dimension of TrainingDF
  dim_train <- paste(dim(train), collapse = " x ")
  print(paste("Dimension of Train-DF:", dim_train))
  
  # [2-6] Fit a model [standard settings] on the train DF & take the time:
  start        <- Sys.time()
  curr_forrest <- rfsrc(as.factor(gender) ~ ., data = train, seed = 123,
                        ntree = 250)
  end          <- Sys.time()
  time_diff    <- difftime(end, start, units = 'mins')
  
  # [2-7] Evaluate the Model:
  # OOB:      Get OOB Accuracy of each of the 500 trees and average it!
  curr_OOB_ACC <- 1 - curr_forrest$err.rate[length(curr_forrest$err.rate)]
  
  # TESTSET:  Get the ErrorRates w/ TestSet
  prob_preds   <- predict(curr_forrest, test)$predicted
  pred_class   <- ifelse(prob_preds > 0.5, 1, 0)
  truth        <- as.factor(test$gender)
  pred_class   <- factor(pred_class, levels = levels(truth))
  res          <- cbind(pred_class, truth)
  test_acc     <- sum(res[,1] == res[,2]) / nrow(res)
  test_F1      <- caret::confusionMatrix(data      = pred_class, 
                                         reference = truth)$byClass["F1"]
  if (is.na(test_F1)) test_F1 <- 0 
  
  # [2-8] Add results of this block to the Results-DF!
  eval_res[nrow(eval_res) + 1, ] <- c(df, dim_train, 
                                      curr_OOB_ACC, test_acc, 
                                      test_F1, time_diff)
}

# [3] Save the results!
write.csv2(eval_res, "./docs/CV_Res/gender/explorative_subsets/RF_10%_mirna_&_mutation_2,5%_rna_&_cnv_all_blocks_joint.csv",
           row.names = FALSE)

# 1.25% subsets of cnv as features to RF                                    ----
"Use single 'cnv' Block w/ 1.25% of the original features & 
 get performance, runtime...
"

# [1] Empty DF - for all results of the Evaluation:
eval_res <- data.frame("Data"     = character(), "Block"    = numeric(), 
                       "Dim"      = numeric(),   "OOB_Acc"  = numeric(),   
                       "Test_Acc" = numeric(),   "Test_F1"  = numeric(), 
                       "Time"     = numeric(), # min
                       stringsAsFactors = F)

# [2] Loop over all DFs w/ gender in their clinical variables!
for (df in DFs_w_gender) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # [2-1] Load 'df' & only keep 'rna' & 'cnv' + the gender we use as response
  load(paste0(data_path, df, ".Rda"))
  omics_blocks <- c("cnv")
  resp         <- clin["gender"]
  
  # [2-2] Loop over all blocks in 'omics_blocks' of current 'df'
  for (curr_block in omics_blocks) {
    
    writeLines(paste0("Current OmicsBlock: \n", curr_block))
    
    # 2-2-1 Subset the current Block and  only keep 5% of the original feas
    curr_data    <- eval(as.symbol(curr_block))
    set.seed(12345)
    sampled_cols <- sample(ncol(curr_data), round(ncol(curr_data) * 0.0125))
    
    # 2-2-2 Bind gender from the clin block & 'curr_block' (=omics), 
    #       convert it to a DF and name the response column
    DF <- cbind(resp, curr_data[ ,sampled_cols])
    DF <- as.data.frame(DF)
    colnames(DF)[1] <- "gender"
    
    # 2-2-3 Split to train and test set [80:20] w/ seed for reproducibility!
    set.seed(12345)
    train_ind <- sample(x = seq_len(nrow(DF)), size = round(nrow(DF) * 0.8))
    test_ind  <- which(!(1:nrow(DF) %in% train_ind))
    train  <- DF[train_ind,]; test <- DF[test_ind,]
    
    # 2-2-3-1 Print dimension of TrainingDF
    dim_train <- paste(dim(train), collapse = " x ")
    print(paste("Dimension of Train-DF:", dim_train))
    
    # 2-2-4 Fit a model [standard settings] on the train DF & take the time:
    start        <- Sys.time()
    curr_forrest <- rfsrc(as.factor(gender) ~ ., data = train, seed = 123,
                          ntree = 250)
    end          <- Sys.time()
    time_diff    <- difftime(end, start, units = 'mins')
    
    # 2-2-6 Evaluate the Model:
    # OOB:      Get OOB Accuracy of each of the 500 trees and average it!
    curr_OOB_ACC <- 1 - curr_forrest$err.rate[length(curr_forrest$err.rate)]
    
    # TESTSET:  Get the ErrorRates w/ TestSet
    prob_preds   <- predict(curr_forrest, test)$predicted
    pred_class   <- ifelse(prob_preds > 0.5, 1, 0)
    truth        <- as.factor(test$gender)
    pred_class   <- factor(pred_class, levels = levels(truth))
    res          <- cbind(pred_class, truth)
    test_acc     <- sum(res[,1] == res[,2]) / nrow(res)
    test_F1      <- caret::confusionMatrix(data      = pred_class, 
                                           reference = truth)$byClass["F1"]
    if (is.na(test_F1)) test_F1 <- 0 
    
    # [2-2-6] Add results of this block to the Results-DF!
    eval_res[nrow(eval_res) + 1, ] <- c(df, curr_block, dim_train, 
                                        curr_OOB_ACC, test_acc, test_F1, time_diff)
  }
}

# [3] Save the results!
write.csv2(eval_res, "./docs/CV_Res/gender/explorative_subsets/RF_1,25%_of_cnv_single_block.csv",
           row.names = FALSE)

# 10% mirna & mutation + 2.5% rna + 1.25% cnv as joint features to RF       ----
" Subset mirna & mutation blocks by 10% and use 2.5% subsets of 'rna' & 1.25%
  subset of 'cnv' block, as featurespace to predcit the gender!
"
eval_res <- data.frame("Data"     = character(), "Dim"      = numeric(),   
                       "OOB_Acc"  = numeric(),   "Test_Acc" = numeric(),   
                       "Test_F1"  = numeric(),   "Time"     = numeric(), # min
                       stringsAsFactors = F)

# [2] Loop over all DFs w/ gender in their clinical variables!
for (df in DFs_w_gender) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # [2-1] Load 'df' & only keep names of the omics blocks in 'df'!
  omics_blocks <- load(paste0(data_path, df, ".Rda"))
  omics_10     <- c("mirna", "mutation")
  omics_25     <- c("rna")
  omics_125    <- c("cnv")
  omics_all    <- c("clin")
  
  resp <- clin$gender
  clin$gender <- NULL
  
  # [2-2] Loop over all blocks in 'omics_blocks' of current 'df' to prune them!
  # [2-2-1] Blocks that get reduced by 10%
  for (curr_block in omics_10) {
    
    writeLines(paste0("Remove 90% of the features of OmicsBlock:\n", curr_block))
    
    # 2-2-1-1 Subset the current Block and only keep 10% of the original feas
    curr_data <- eval(as.symbol(curr_block))
    set.seed(12345)
    sampled_cols <- sample(seq_len(ncol(curr_data)), round(ncol(curr_data) * 0.1))
    
    # 2-2-1-2 Only keep the sampled columns in the block 
    #        [overwrite current variable w/ reduced dimension var!]
    assign(curr_block, curr_data[ ,sampled_cols])
  }
  
  # [2-2-2] Blocks that get reduced by 5%
  for (curr_block in omics_25) {
    
    writeLines(paste0("Remove 97.5% of the features of OmicsBlock:\n", curr_block))
    
    # 2-2-1-1 Subset the current Block and only keep 10% of the original feas
    curr_data <- eval(as.symbol(curr_block))
    set.seed(12345)
    sampled_cols <- sample(seq_len(ncol(curr_data)), round(ncol(curr_data) * 0.025))
    
    # 2-2-1-2 Only keep the sampled columns in the block 
    #        [overwrite current variable w/ reduced dimension var!]
    assign(curr_block, curr_data[ ,sampled_cols])
  }
  
  # [2-2-3] Blocks that get reduced by 5%
  for (curr_block in omics_125) {
    
    writeLines(paste0("Remove 98.75% of the features of OmicsBlock:\n", curr_block))
    
    # 2-2-1-1 Subset the current Block and only keep 10% of the original feas
    curr_data <- eval(as.symbol(curr_block))
    set.seed(12345)
    sampled_cols <- sample(seq_len(ncol(curr_data)), round(ncol(curr_data) * 0.0125))
    
    # 2-2-1-2 Only keep the sampled columns in the block 
    #        [overwrite current variable w/ reduced dimension var!]
    assign(curr_block, curr_data[ ,sampled_cols])
  }
  
  # [2-3] Bind gender from the clin block & all blocks!
  # [2-3-1] As rna & cnv share colnames [which leads to an error when fitting RF]
  #         we rename the colnames of the rna DF
  colnames(rna) <- paste0(colnames(rna), "_rna")
  
  # [2-3-2] Paste the subsetted DFs together!
  DF <- cbind(resp, 
              eval(as.symbol(omics_10[1])), 
              eval(as.symbol(omics_10[2])), 
              eval(as.symbol(omics_25[1])),
              eval(as.symbol(omics_125[1])),
              clin)
  
  # [2-4] Convert to DF, and name the response variable
  DF <- as.data.frame(DF); colnames(DF)[1] <- "gender"
  
  # [2-5] Split to train and test set [80:20] w/ seed for reproducibility!
  set.seed(12345)
  train_ind <- sample(x = seq_len(nrow(DF)), size = round(nrow(DF) * 0.8))
  test_ind  <- which(!(1:nrow(DF) %in% train_ind))
  train  <- DF[train_ind,]; test <- DF[test_ind,]
  
  # [2-5-1] Print dimension of TrainingDF
  dim_train <- paste(dim(train), collapse = " x ")
  print(paste("Dimension of Train-DF:", dim_train))
  
  # [2-6] Fit a model [standard settings] on the train DF & take the time:
  start        <- Sys.time()
  curr_forrest <- rfsrc(as.factor(gender) ~ ., data = train, seed = 123,
                        ntree = 250)
  end          <- Sys.time()
  time_diff    <- difftime(end, start, units = 'mins')
  
  # [2-7] Evaluate the Model:
  # OOB:      Get OOB Accuracy of each of the 500 trees and average it!
  curr_OOB_ACC <- 1 - curr_forrest$err.rate[length(curr_forrest$err.rate)]
  
  # TESTSET:  Get the ErrorRates w/ TestSet
  prob_preds   <- predict(curr_forrest, test)$predicted
  pred_class   <- ifelse(prob_preds > 0.5, 1, 0)
  truth        <- as.factor(test$gender)
  pred_class   <- factor(pred_class, levels = levels(truth))
  res          <- cbind(pred_class, truth)
  test_acc     <- sum(res[,1] == res[,2]) / nrow(res)
  test_F1      <- caret::confusionMatrix(data      = pred_class, 
                                         reference = truth)$byClass["F1"]
  if (is.na(test_F1)) test_F1 <- 0 
  
  # [2-8] Add results of this block to the Results-DF!
  eval_res[nrow(eval_res) + 1, ] <- c(df, dim_train, 
                                      curr_OOB_ACC, test_acc, 
                                      test_F1, time_diff)
}

# [3] Save the results!
write.csv2(eval_res, "./docs/CV_Res/gender/explorative_subsets/RF_10%_mirna_&_mutation_2,5%_rna_1,25%_cnv_all_blocks_joint.csv",
           row.names = FALSE)