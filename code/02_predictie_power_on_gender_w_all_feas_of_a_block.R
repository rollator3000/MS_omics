" Investigate how good the predicitive peroformance of a model is on a single 
  block / in General! 
    --> full block with the whole features that are originally in the block!
    --> block with only 10% | 5% | 2.5%... of the original features!
    --> all blocks - w/ only 10% | 5% | 2.5%... of original features conecated!

  For each of these Scenarios, we create a DF, that tracks:
    -the dataframe  
    -the block  
    -Dim of Traindata  
    -OOB & TestSet Accuracy  
    -Time in minutes
"
# Set Working Directory and load the needed packages!
setwd("C:/Users/kuche_000/Desktop/MS-Thesis/")
library(randomForest)

# Names of the usable dataframes (w/ gender in 'clin'-block & 4 omics blocks!)
DFs_w_gender <- c("BLCA", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC","LGG", 
                  "LUAD", "LUSC", "PAAD", "SARC", "SKCM", "STAD")
data_path    <- "./data/external/Dr_Hornung/Data/ProcessedData/"

# Single fully observed Omics-Blocks as features to RF                       ----
" Use single omics Blocks and get a performance of a RF how it performs!
  Get performance, runtime, etc.
     --> local machine can maximally handle 12.5k feas for a tree
     --> not all blocks could be used as full block!
"
# [1] Empty DF - for all results of the Evaluation:
eval_res <- data.frame("Data"  = character(),  "Block" = character(),
                       "Dim"   = numeric(),    "OOB"   = numeric(),
                       "Test"  = numeric(),    "Time"  = numeric(), # in minutes
                       stringsAsFactors = F)

# [2] Loop over all DFs w/ gender in their clinical variables!
#     For Blocks w/ more than 12.5k feas, we skip them, as it would crash!
for (df in DFs_w_gender) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # [2-1] Load 'df' & only keep names of the omics blocks in 'df'!
  omics_blocks <- load(paste0(data_path, df, ".Rda"))
  omics_blocks <- omics_blocks[-which(omics_blocks %in% c("clin", "targetvar"))]
  
  # [2-2] Loop over all blocks in 'omics_blocks' of current 'df'
  for (curr_block in  omics_blocks) {
    
    writeLines(paste0("Current OmicsBlock: \n", curr_block))
    
    #[2-2-1] Bind gender from the clin block & 'curr_block' (=omics) 
    DF <- cbind(clin$gender, eval(as.symbol(curr_block)))
    
    # [2-2-2] Convert to DF, and name the response variable
    DF <- as.matrix(DF); colnames(DF)[1] <- "gender"
    
    # [2-2-3] Split to train and test set [80:20] 
    test  <- DF[(round((nrow(DF)*0.8 )) + 1):nrow(DF),]
    train <- DF[1:round(nrow(DF)*0.8),]
    
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
    curr_forrest <- randomForest(as.factor(gender) ~ ., data = train)
    end          <- Sys.time()
    time_diff    <- difftime(end, start, units = 'mins')
    
    # [2-2-6] Evaluate the Model:
    # OOB:    Get OOB Accuracy of each of the 500 trees and average it!
    curr_OOB_ACC <- mean(1 - curr_forrest$err.rate[,1])

    # TESTSET: Get the ErrorRates w/ TestSet
    prediction   <- predict(curr_forrest, test); truth <- as.factor(test[,1])
    res          <- cbind(prediction, truth)
    test_acc     <- sum(res[,1] == res[,2]) / nrow(res)
    
    # [2-2-7] Add results of this block to the Results-DF!
    eval_res[nrow(eval_res) + 1, ] <- c(df, curr_block, dim_train, curr_OOB_ACC, 
                                        test_acc, time_diff)
    
  }
}

# Save the results!
write.csv2(eval_res, "./docs/explorative_subsets/performance_RF_whole_single_blocks_gender_classif.csv",
           row.names = FALSE)
# Single 10% subsetted Omics-Blocks as features to RF                        ----
"Use single Blocks w/ 10% of the original features & get performance, runtime...
 There shouldn't be any more computional issues now!
"
# [1] Empty DF - for all results of the Evaluation:
eval_res_sub <- data.frame("Data"  = character(),  "Block" = character(),
                           "Dim"   = numeric(),    "OOB"   = numeric(),
                           "Test"  = numeric(),    "Time"  = numeric(), # in min
                           stringsAsFactors = F)

# [2] Loop over all DFs w/ gender in their clinical variables!
for (df in DFs_w_gender) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # [2-1] Load 'df' & only keep names of the omics blocks in 'df'!
  omics_blocks <- load(paste0(data_path, df, ".Rda"))
  omics_blocks <- omics_blocks[-which(omics_blocks %in% c("clin", "targetvar"))]
  
  # [2-2] Loop over all blocks in 'omics_blocks' of current 'df'
  for (curr_block in omics_blocks) {
    
    writeLines(paste0("Current OmicsBlock: \n", curr_block))
    
    # [2-2-1-1] Subset the current Block and  only keep 10% of the original feas
    curr_data    <- eval(as.symbol(curr_block))
    sampled_cols <- sample(ncol(curr_data), round(ncol(curr_data) * 0.1))
    
    # [2-2-1-2] Only keep the sampled columns in 'curr_block'
    assign(curr_block, curr_data[ ,sampled_cols])
    
    # [2-2-1-3] Bind gender from the clin block & 'curr_block' (=omics)
    DF <- cbind(clin$gender, eval(as.symbol(curr_block)))
    
    
    # [2-2-2] Convert to DF, and name the response variable
    DF <- as.matrix(DF); colnames(DF)[1] <- "gender"
    
    # [2-2-3] Split to train and test set [80:20] 
    test  <- DF[(round((nrow(DF)*0.8 )) + 1):nrow(DF),]
    train <- DF[1:round(nrow(DF)*0.8),]
    
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
    curr_forrest <- randomForest(as.factor(gender) ~ ., data = train)
    end          <- Sys.time()
    time_diff    <- difftime(end, start, units = 'mins')
    
    # [2-2-6] Evaluate the Model:
    # OOB:    Get OOB Accuracy of each of the 500 trees and average it!
    curr_OOB_ACC <- mean(1 - curr_forrest$err.rate[,1])
    
    # TESTSET: Get the ErrorRates w/ TestSet
    prediction   <- predict(curr_forrest, test); truth <- as.factor(test[,1])
    res          <- cbind(prediction, truth)
    test_acc     <- sum(res[,1] == res[,2]) / nrow(res)
    
    # [2-2-7] Add results of this block to the Results-DF!
    eval_res_sub[nrow(eval_res_sub) + 1, ] <- c(df, curr_block, dim_train, 
                                                curr_OOB_ACC, test_acc, time_diff)
  }
}

# [3] Save the results!
write.csv2(eval_res_sub, "./docs/explorative_subsets/performance_RF_10percent_single_blocks_gender_classif.csv",
           row.names = FALSE)

# [4] Litte analyse w/ results:
# [4-0] Load the results!
eval_res_sub = read.csv2("./docs/explorative_subsets/performance_RF_10percent_single_blocks_gender_classif.csv", 
                         stringsAsFactors = F)

# [4-1] convert to numeric:
cols_ <- c("OOB", "Test", "Time")
eval_res_sub[cols_] <- sapply(eval_res_sub[cols_], as.numeric)

# [4-2] Get the mean TestAccuracy / OOB-Accuracy of a certain block!
# [4-2-1] TEST
TEST <- sapply(unique(eval_res_sub$Block), 
       FUN = function(x) summary(eval_res_sub$Test[eval_res_sub$Block == x]))
# [4-2-2] OOB
OOB <- sapply(unique(eval_res_sub$Block), 
       FUN = function(x) summary(eval_res_sub$OOB[eval_res_sub$Block == x]))

# [4-4-3] Apply meaningful names and check the data
colnames(OOB)  <- paste0(colnames(OOB), "_oob")
colnames(TEST) <- paste0(colnames(TEST), "_test")

# [4-5] General Overview of Performance on the single blocks!
cbind(OOB, TEST)

# 10% subsets of all blocks as feature to RF                                 ----
"Remove 90% of all omics blocks features and paste them together to learn a 
 model on the joint data"

# [1] Empty DF - for all results of the Evaluation:
eval_res_sub_all <- data.frame("Data"  = character(), "Dim"   = numeric(),   
                               "OOB"   = numeric(),   "Test"  = numeric(),   
                               "Time"  = numeric(),    stringsAsFactors = F)

# [2] Loop over all DFs w/ gender in their clinical variables!
for (df in DFs_w_gender) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # [2-1] Load 'df' & only keep names of the omics blocks in 'df'!
  omics_blocks <- load(paste0(data_path, df, ".Rda"))
  omics_blocks <- omics_blocks[-which(omics_blocks %in% c("clin", "targetvar"))]
  
  # [2-2] Loop over all blocks in 'omics_blocks' of current 'df' to prune them!
  for (curr_block in omics_blocks) {
    
    writeLines(paste0("Remove 90% of the features of OmicsBlock:\n", curr_block))
    
    # [2-2-1] Subset the current Block and only keep 10% of the original feas
    curr_data    <- eval(as.symbol(curr_block))
    sampled_cols <- sample(ncol(curr_data), round(ncol(curr_data) * 0.1))
    
    # [2-2-2] Only keep the sampled columns in the block 
    #         [overwrite current variable w/ reduced dimension var!]
    assign(curr_block, curr_data[ ,sampled_cols])
  }
  
  # [2-3] Bind gender from the clin block & all blocks!
  # [2-3-1] As rna & cnv share colnames [which leads to an error when fitting RF]
  #         we rename the colnames of the rna DF
  colnames(rna) <- paste0(colnames(rna), "_rna")
  
  # [2-3-2] Paste the subsetted DFs together!
  DF <- cbind(clin$gender, 
              eval(as.symbol(omics_blocks[1])), 
              eval(as.symbol(omics_blocks[2])), 
              eval(as.symbol(omics_blocks[3])),
              eval(as.symbol(omics_blocks[4])))
  
  # [2-4] Convert to DF, and name the response variable
  DF <- as.matrix(DF); colnames(DF)[1] <- "gender"
  
  # [2-5] Split to train and test set [80:20] 
  test  <- DF[(round((nrow(DF)*0.8 )) + 1):nrow(DF),]
  train <- DF[1:round(nrow(DF)*0.8),]
  
  # [2-5-1] Print dimension of TrainingDF
  dim_train <- paste(dim(train), collapse = " x ")
  dim_test  <- paste(dim(test), collapse = " x ")
  print(paste("Dimension of Train-DF:", dim_train))
  print(paste("Dimension of Test-DF:", dim_test))
  
  # [2-5-2] If the block has more than 12500 feas, stop it, as overflow error,
  #         when fitting RF to it!
  if (dim(train)[2] > 12500) {
    print("too high dimensional! Will crash the computer!")
    next
  }
  
  # [2-6] Fit a model [standard settings] on the train DF & take the time:
  start        <- Sys.time()
  curr_forrest <- randomForest(as.factor(gender) ~ ., data = train)
  end          <- Sys.time()
  time_diff    <- difftime(end, start, units = 'mins')
  
  # [2-7] Evaluate the Model:
  # OOB:    Get OOB Accuracy of each of the 500 trees and average it!
  curr_OOB_ACC <- mean(1 - curr_forrest$err.rate[,1])
  
  # TESTSET: Get the ErrorRates w/ TestSet
  prediction   <- predict(curr_forrest, test); truth <- as.factor(test[,1])
  res          <- cbind(prediction, truth)
  test_acc     <- sum(res[,1] == res[,2]) / nrow(res)
  
  # [2-8] Add results of this block to the Results-DF!
  eval_res_sub_all[nrow(eval_res_sub_all) + 1, ] <- c(df, dim_train, 
                                                      curr_OOB_ACC, test_acc, 
                                                      time_diff)
}

# [3] Save the results!
write.csv2(eval_res_sub_all, "./docs/explorative_subsets/performance_RF_10percent_all_blocks_gender_classif.csv",
           row.names = FALSE)

# [4] Litte analyse w/ results:
# [4-0] Load the results!
eval_res_sub_all = read.csv2("./docs/explorative_subsets/performance_RF_10percent_all_blocks_gender_classif.csv", 
                             stringsAsFactors = F)

# [4-1] convert to correct types:
cols_ <- c("OOB", "Test", "Time")
eval_res_sub_all[cols_] <- sapply(eval_res_sub_all[cols_], as.numeric)
eval_res_sub_all$Data <- as.factor(eval_res_sub_all$Data)

# [4-2] Get average Test/OOB Accuracy  
# [4-2-1] TEST
TEST <- summary(eval_res_sub_all$Test)
# [4-2-2] OOB
OOB <- summary(eval_res_sub_all$OOB)

# [4-3] General Overview of Performance on the single blocks!
cbind(OOB, TEST)

# 5% subsets of 'rna'/'cnv' as features to RF                                ----
" For the blocks 'rna' / 'cnv' we use 5% subsets and train a RF on these [seperat] 
"

# [1] Empty DF - for all results of the Evaluation:
eval_res_sub_2 <- data.frame("Data"  = character(),  "Block" = character(),
                             "Dim"   = numeric(),    "OOB"   = numeric(),
                             "Test"  = numeric(),    "Time"  = numeric(), # in min
                             stringsAsFactors = F)

# [2] Loop over all DFs w/ gender in their clinical variables!
for (df in DFs_w_gender) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # [2-1] Load 'df' & only keep names of 'rna' & 'cnv'
  load(paste0(data_path, df, ".Rda"))
  omics_blocks <- c("rna", "cnv")
  
  # [2-2] Loop over all blocks in 'omics_blocks' of current 'df'
  for (curr_block in omics_blocks) {
    
    writeLines(paste0("Current OmicsBlock: \n", curr_block))
    
    # [2-2-1-1] Subset the current Block and  only keep 5% of the original feas
    curr_data    <- eval(as.symbol(curr_block))
    sampled_cols <- sample(ncol(curr_data), round(ncol(curr_data) * 0.05))
    
    # [2-2-1-2] Only keep the sampled columns in 'curr_block'
    assign(curr_block, curr_data[ ,sampled_cols])
    
    # [2-2-1-3] Bind gender from the clin block & 'curr_block' (=omics)
    DF <- cbind(clin$gender, eval(as.symbol(curr_block)))
    
    # [2-2-2] Convert to DF, and name the response variable
    DF <- as.matrix(DF); colnames(DF)[1] <- "gender"
    
    # [2-2-3] Split to train and test set [80:20] 
    test  <- DF[(round((nrow(DF)*0.8 )) + 1):nrow(DF),]
    train <- DF[1:round(nrow(DF)*0.8),]
    
    # [2-2-3-1] Print dimension of TrainingDF
    dim_train <- paste(dim(train), collapse = " x ")
    print(paste("Dimension of Train-DF:", dim_train))
    
    # [2-2-4] Fit a model [standard settings] on the train DF & take the time:
    start        <- Sys.time()
    curr_forrest <- randomForest(as.factor(gender) ~ ., data = train)
    end          <- Sys.time()
    time_diff    <- difftime(end, start, units = 'mins')
    
    # [2-2-5] Evaluate the Model:
    # OOB:    Get OOB Accuracy of each of the 500 trees and average it!
    curr_OOB_ACC <- mean(1 - curr_forrest$err.rate[,1])
    
    # TESTSET: Get the ErrorRates w/ TestSet
    prediction   <- predict(curr_forrest, test); truth <- as.factor(test[,1])
    res          <- cbind(prediction, truth)
    test_acc     <- sum(res[,1] == res[,2]) / nrow(res)
    
    # [2-2-6] Add results of this block to the Results-DF!
    eval_res_sub_2[nrow(eval_res_sub_2) + 1, ] <- c(df, curr_block, dim_train, 
                                                    curr_OOB_ACC, test_acc, time_diff)
  }
}

# [3] Save the results!
write.csv2(eval_res_sub_2, "./docs/explorative_subsets/performance_RF_5percent_rna_cnv_single_block_gender_classif.csv",
           row.names = FALSE)

# [4] Litte analyse w/ results:
# [4-0] Load the results!
eval_res_sub_2 = read.csv2("./docs/explorative_subsets/performance_RF_5percent_rna_cnv_single_block_gender_classif.csv",
                         stringsAsFactors = F)

# [4-1] convert to numeric:
cols_ <- c("OOB", "Test", "Time")
eval_res_sub_2[cols_] <- sapply(eval_res_sub_2[cols_], as.numeric)

# [4-2] Get the mean TestAccuracy / OOB-Accuracy of a certain block!
# [4-2-1] TEST
TEST <- sapply(unique(eval_res_sub_2$Block),
               FUN = function(x) summary(eval_res_sub_2$Test[eval_res_sub_2$Block == x]))
# [4-2-2] OOB
OOB <- sapply(unique(eval_res_sub_2$Block),
              FUN = function(x) summary(eval_res_sub_2$OOB[eval_res_sub_2$Block == x]))

# [4-4-3] Apply meaningful names and check the data
colnames(OOB)  <- paste0(colnames(OOB), "_oob")
colnames(TEST) <- paste0(colnames(TEST), "_test")

# [4-5] General Overview of Performance on the single blocks!
cbind(OOB, TEST)


# Use 10% mirna & mutation + 5% cnv & rna all together                       ----
" Subset mirna & mutation blocks by 10% and use 5% subsets of the 'cnv' & 'rna'
  blocks, as featurespace to predcit the gender!
"
eval_res_sub_all2 <- data.frame("Data"  = character(), "Dim"   = numeric(),   
                                "OOB"   = numeric(),   "Test"  = numeric(),   
                                "Time"  = numeric(),    stringsAsFactors = F)

# [2] Loop over all DFs w/ gender in their clinical variables!
for (df in DFs_w_gender) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # [2-1] Load 'df' & only keep names of the omics blocks in 'df'!
  omics_blocks <- load(paste0(data_path, df, ".Rda"))
  omics_10     <- c("mirna", "mutation")
  omics_5      <- c("rna", "cnv")
  
  # [2-2] Loop over all blocks in 'omics_blocks' of current 'df' to prune them!
  
  # [2-2-1] Blocks that get reduced by 10%
  for (curr_block in omics_10) {
    
    writeLines(paste0("Remove 90% of the features of OmicsBlock:\n", curr_block))
    
    # [2-2-1-1] Subset the current Block and only keep 10% of the original feas
    curr_data    <- eval(as.symbol(curr_block))
    sampled_cols <- sample(ncol(curr_data), round(ncol(curr_data) * 0.1))
    
    # [2-2-1-2] Only keep the sampled columns in the block 
    #          [overwrite current variable w/ reduced dimension var!]
    assign(curr_block, curr_data[ ,sampled_cols])
  }
  
  # [2-2-2] Blocks that get reduced by 5%
  for (curr_block in omics_5) {
    
    writeLines(paste0("Remove 95% of the features of OmicsBlock:\n", curr_block))
    
    # [2-2-2-1] Subset the current Block and only keep 10% of the original feas
    curr_data    <- eval(as.symbol(curr_block))
    sampled_cols <- sample(ncol(curr_data), round(ncol(curr_data) * 0.05))
    
    # [2-2-2-2] Only keep the sampled columns in the block 
    #          [overwrite current variable w/ reduced dimension var!]
    assign(curr_block, curr_data[ ,sampled_cols])
  }
  
  # [2-3] Bind gender from the clin block & all blocks!
  # [2-3-1] As rna & cnv share colnames [which leads to an error when fitting RF]
  #         we rename the colnames of the rna DF
  colnames(rna) <- paste0(colnames(rna), "_rna")
  
  # [2-3-2] Paste the subsetted DFs together!
  DF <- cbind(clin$gender, 
              eval(as.symbol(omics_10[1])), 
              eval(as.symbol(omics_10[2])), 
              eval(as.symbol(omics_5[1])),
              eval(as.symbol(omics_5[2])))
  
  # [2-4] Convert to DF, and name the response variable
  DF <- as.matrix(DF); colnames(DF)[1] <- "gender"
  
  # [2-5] Split to train and test set [80:20] 
  test  <- DF[(round((nrow(DF)*0.8 )) + 1):nrow(DF),]
  train <- DF[1:round(nrow(DF)*0.8),]
  
  # [2-5-1] Print dimension of TrainingDF
  dim_train <- paste(dim(train), collapse = " x ")
  dim_test  <- paste(dim(test), collapse = " x ")
  print(paste("Dimension of Train-DF:", dim_train))
  print(paste("Dimension of Test-DF:", dim_test))
  
  # [2-6] Fit a model [standard settings] on the train DF & take the time:
  start        <- Sys.time()
  curr_forrest <- randomForest(as.factor(gender) ~ ., data = train)
  end          <- Sys.time()
  time_diff    <- difftime(end, start, units = 'mins')
  
  # [2-7] Evaluate the Model:
  # OOB:    Get OOB Accuracy of each of the 500 trees and average it!
  curr_OOB_ACC <- mean(1 - curr_forrest$err.rate[,1])
  
  # TESTSET: Get the ErrorRates w/ TestSet
  prediction   <- predict(curr_forrest, test); truth <- as.factor(test[,1])
  res          <- cbind(prediction, truth)
  test_acc     <- sum(res[,1] == res[,2]) / nrow(res)
  
  # [2-8] Add results of this block to the Results-DF!
  eval_res_sub_all2[nrow(eval_res_sub_all2) + 1, ] <- c(df, dim_train, 
                                                        curr_OOB_ACC, test_acc, 
                                                        time_diff)
}

# [3] Save the results!
write.csv2(eval_res_sub_all2, "./docs/explorative_subsets/performance_RF_10percent_mirna_mutation_5percent_rna_cnv_all_blocks_gender_classif.csv",
           row.names = FALSE)

# [4] Litte analyse w/ results:
# [4-0] Load the results!
eval_res_sub_all = read.csv2("./docs/explorative_subsets/performance_RF_10percent_mirna_mutation_5percent_rna_cnv_all_blocks_gender_classif.csv", 
                             stringsAsFactors = F)

# [4-1] convert to correct types:
cols_ <- c("OOB", "Test", "Time")
eval_res_sub_all[cols_] <- sapply(eval_res_sub_all[cols_], as.numeric)
eval_res_sub_all$Data <- as.factor(eval_res_sub_all$Data)

# [4-2] Get average Test/OOB Accuracy  
# [4-2-1] TEST
TEST <- summary(eval_res_sub_all$Test)
# [4-2-2] OOB
OOB <- summary(eval_res_sub_all$OOB)

# [4-3] General Overview of Performance on the single blocks!
cbind(OOB, TEST)


# 2.5% subsets of 'rna'&'cnv' as features to RF                              ----
"Use single 'rna'/'cnv' Blocks w/ 2.5% of the original features 
 & get performance, runtime...
 There should be no issues with the computer now!
"

# [1] Empty DF - for all results of the Evaluation:
eval_res_sub_2 <- data.frame("Data"  = character(),  "Block" = character(),
                             "Dim"   = numeric(),    "OOB"   = numeric(),
                             "Test"  = numeric(),    "Time"  = numeric(), # in min
                             stringsAsFactors = F)

# [2] Loop over all DFs w/ gender in their clinical variables!
for (df in DFs_w_gender) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # [2-1] Load 'df' & only keep names of 'rna' & 'cnv'
  load(paste0(data_path, df, ".Rda"))
  omics_blocks <- c("rna", "cnv")
  
  # [2-2] Loop over all blocks in 'omics_blocks' of current 'df'
  for (curr_block in omics_blocks) {
    
    writeLines(paste0("Current OmicsBlock: \n", curr_block))
    
    # [2-2-1-1] Subset the current Block and  only keep 5% of the original feas
    curr_data    <- eval(as.symbol(curr_block))
    sampled_cols <- sample(ncol(curr_data), round(ncol(curr_data) * 0.025))
    
    # [2-2-1-2] Only keep the sampled columns in 'curr_block'
    assign(curr_block, curr_data[ ,sampled_cols])
    
    # [2-2-1-3] Bind gender from the clin block & 'curr_block' (=omics)
    DF <- cbind(clin$gender, eval(as.symbol(curr_block)))
    
    # [2-2-2] Convert to DF, and name the response variable
    DF <- as.matrix(DF); colnames(DF)[1] <- "gender"
    
    # [2-2-3] Split to train and test set [80:20] 
    test  <- DF[(round((nrow(DF)*0.8 )) + 1):nrow(DF),]
    train <- DF[1:round(nrow(DF)*0.8),]
    
    # [2-2-3-1] Print dimension of TrainingDF
    dim_train <- paste(dim(train), collapse = " x ")
    print(paste("Dimension of Train-DF:", dim_train))
    
    # [2-2-4] Fit a model [standard settings] on the train DF & take the time:
    start        <- Sys.time()
    curr_forrest <- randomForest(as.factor(gender) ~ ., data = train)
    end          <- Sys.time()
    time_diff    <- difftime(end, start, units = 'mins')
    
    # [2-2-5] Evaluate the Model:
    # OOB:    Get OOB Accuracy of each of the 500 trees and average it!
    curr_OOB_ACC <- mean(1 - curr_forrest$err.rate[,1])
    
    # TESTSET: Get the ErrorRates w/ TestSet
    prediction   <- predict(curr_forrest, test); truth <- as.factor(test[,1])
    res          <- cbind(prediction, truth)
    test_acc     <- sum(res[,1] == res[,2]) / nrow(res)
    
    # [2-2-6] Add results of this block to the Results-DF!
    eval_res_sub_2[nrow(eval_res_sub_2) + 1, ] <- c(df, curr_block, dim_train, 
                                                    curr_OOB_ACC, test_acc, time_diff)
  }
}

# [3] Save the results!
write.csv2(eval_res_sub_2, "./docs/explorative_subsets/performance_RF_2.5percent_rna_cnv_single_block_gender_classif.csv",
           row.names = FALSE)

# [4] Litte analyse w/ results:
# [4-0] Load the results!
eval_res_sub_2 = read.csv2("./docs/explorative_subsets/performance_RF_2.5percent_rna_cnv_single_block_gender_classif.csv",
                         stringsAsFactors = F)

# [4-1] convert to numeric:
cols_ <- c("OOB", "Test", "Time")
eval_res_sub_2[cols_] <- sapply(eval_res_sub_2[cols_], as.numeric)

# [4-2] Get the mean TestAccuracy / OOB-Accuracy of a certain block!
# [4-2-1] TEST
TEST <- sapply(unique(eval_res_sub_2$Block),
               FUN = function(x) summary(eval_res_sub_2$Test[eval_res_sub_2$Block == x]))
# [4-2-2] OOB
OOB <- sapply(unique(eval_res_sub_2$Block),
              FUN = function(x) summary(eval_res_sub_2$OOB[eval_res_sub_2$Block == x]))

# [4-4-3] Apply meaningful names and check the data
colnames(OOB)  <- paste0(colnames(OOB), "_oob")
colnames(TEST) <- paste0(colnames(TEST), "_test")

# [4-5] General Overview of Performance on the single blocks!
cbind(OOB, TEST)

# 10% mirna & mutation + 2.5% cnv & rna as joint features to RF              ----
" Subset mirna & mutation blocks by 10% and use 2.5% subsets of the 'cnv' & 'rna'
  blocks, as featurespace to predcit the gender!
"
eval_res_sub_all2 <- data.frame("Data"  = character(), "Dim"   = numeric(),   
                                "OOB"   = numeric(),   "Test"  = numeric(),   
                                "Time"  = numeric(),    stringsAsFactors = F)

# [2] Loop over all DFs w/ gender in their clinical variables!
for (df in DFs_w_gender) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # [2-1] Load 'df' & only keep names of the omics blocks in 'df'!
  omics_blocks <- load(paste0(data_path, df, ".Rda"))
  omics_10     <- c("mirna", "mutation")
  omics_5      <- c("rna", "cnv")
  
  # [2-2] Loop over all blocks in 'omics_blocks' of current 'df' to prune them!
  
  # [2-2-1] Blocks that get reduced by 10%
  for (curr_block in omics_10) {
    
    writeLines(paste0("Remove 90% of the features of OmicsBlock:\n", curr_block))
    
    # [2-2-1-1] Subset the current Block and only keep 10% of the original feas
    curr_data    <- eval(as.symbol(curr_block))
    sampled_cols <- sample(ncol(curr_data), round(ncol(curr_data) * 0.1))
    
    # [2-2-1-2] Only keep the sampled columns in the block 
    #          [overwrite current variable w/ reduced dimension var!]
    assign(curr_block, curr_data[ ,sampled_cols])
  }
  
  # [2-2-2] Blocks that get reduced by 5%
  for (curr_block in omics_5) {
    
    writeLines(paste0("Remove 97.5% of the features of OmicsBlock:\n", curr_block))
    
    # [2-2-2-1] Subset the current Block and only keep 10% of the original feas
    curr_data    <- eval(as.symbol(curr_block))
    sampled_cols <- sample(ncol(curr_data), round(ncol(curr_data) * 0.025))
    
    # [2-2-2-2] Only keep the sampled columns in the block 
    #          [overwrite current variable w/ reduced dimension var!]
    assign(curr_block, curr_data[ ,sampled_cols])
  }
  
  # [2-3] Bind gender from the clin block & all blocks!
  # [2-3-1] As rna & cnv share colnames [which leads to an error when fitting RF]
  #         we rename the colnames of the rna DF
  colnames(rna) <- paste0(colnames(rna), "_rna")
  
  # [2-3-2] Paste the subsetted DFs together!
  DF <- cbind(clin$gender, 
              eval(as.symbol(omics_10[1])), 
              eval(as.symbol(omics_10[2])), 
              eval(as.symbol(omics_5[1])),
              eval(as.symbol(omics_5[2])))
  
  # [2-4] Convert to DF, and name the response variable
  DF <- as.matrix(DF); colnames(DF)[1] <- "gender"
  
  # [2-5] Split to train and test set [80:20] 
  test  <- DF[(round((nrow(DF)*0.8 )) + 1):nrow(DF),]
  train <- DF[1:round(nrow(DF)*0.8),]
  
  # [2-5-1] Print dimension of TrainingDF
  dim_train <- paste(dim(train), collapse = " x ")
  dim_test  <- paste(dim(test), collapse = " x ")
  print(paste("Dimension of Train-DF:", dim_train))
  print(paste("Dimension of Test-DF:", dim_test))
  
  # [2-6] Fit a model [standard settings] on the train DF & take the time:
  start        <- Sys.time()
  curr_forrest <- randomForest(as.factor(gender) ~ ., data = train)
  end          <- Sys.time()
  time_diff    <- difftime(end, start, units = 'mins')
  
  # [2-7] Evaluate the Model:
  # OOB:    Get OOB Accuracy of each of the 500 trees and average it!
  curr_OOB_ACC <- mean(1 - curr_forrest$err.rate[,1])
  
  # TESTSET: Get the ErrorRates w/ TestSet
  prediction   <- predict(curr_forrest, test); truth <- as.factor(test[,1])
  res          <- cbind(prediction, truth)
  test_acc     <- sum(res[,1] == res[,2]) / nrow(res)
  
  # [2-8] Add results of this block to the Results-DF!
  eval_res_sub_all2[nrow(eval_res_sub_all2) + 1, ] <- c(df, dim_train, 
                                                        curr_OOB_ACC, test_acc, 
                                                        time_diff)
}

# [3] Save the results!
write.csv2(eval_res_sub_all2, "./docs/explorative_subsets/performance_RF_10percent_mirna_mutation_2.5percent_rna_cnv_all_blocks_gender_classif.csv",
           row.names = FALSE)


# [4] Litte analyse w/ results:
# [4-0] Load the results!
eval_res_sub_all = read.csv2("./docs/explorative_subsets/performance_RF_10percent_mirna_mutation_2.5percent_rna_cnv_all_blocks_gender_classif.csv", 
                             stringsAsFactors = F)

# [4-1] convert to correct types:
cols_ <- c("OOB", "Test", "Time")
eval_res_sub_all[cols_] <- sapply(eval_res_sub_all[cols_], as.numeric)
eval_res_sub_all$Data <- as.factor(eval_res_sub_all$Data)

# [4-2] Get average Test/OOB Accuracy  
# [4-2-1] TEST
TEST <- summary(eval_res_sub_all$Test)
# [4-2-2] OOB
OOB <- summary(eval_res_sub_all$OOB)

# [4-3] General Overview of Performance on the single blocks!
cbind(OOB, TEST)

# 1.25% subsets of cnv' as features to RF                                    ----
"Use single 'cnv' Block w/ 1.25%/ 1% / 0.5%/ 0.25%... [adjust withing the loop] 
 of the original features & get performance, runtime...
"

# [1] Empty DF - for all results of the Evaluation:
eval_res_sub_2 <- data.frame("Data"  = character(),  "Block" = character(),
                             "Dim"   = numeric(),    "OOB"   = numeric(),
                             "Test"  = numeric(),    "Time"  = numeric(), # in min
                             stringsAsFactors = F)

# [2] Loop over all DFs w/ gender in their clinical variables!
for (df in DFs_w_gender) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # [2-1] Load 'df' & only keep names of 'rna' & 'cnv'
  load(paste0(data_path, df, ".Rda"))
  omics_blocks <- c("cnv")
  
  # [2-2] Loop over all blocks in 'omics_blocks' of current 'df'
  for (curr_block in omics_blocks) {
    
    writeLines(paste0("Current OmicsBlock: \n", curr_block))
    
    # [2-2-1-1] Subset the current Block and  only keep 5% of the original feas
    curr_data    <- eval(as.symbol(curr_block))
    sampled_cols <- sample(ncol(curr_data), round(ncol(curr_data) * 0.001))
    
    # [2-2-1-2] Only keep the sampled columns in 'curr_block'
    assign(curr_block, curr_data[ ,sampled_cols])
    
    # [2-2-1-3] Bind gender from the clin block & 'curr_block' (=omics)
    DF <- cbind(clin$gender, eval(as.symbol(curr_block)))
    
    # [2-2-2] Convert to DF, and name the response variable
    DF <- as.matrix(DF); colnames(DF)[1] <- "gender"
    
    # [2-2-3] Split to train and test set [80:20] 
    test  <- DF[(round((nrow(DF)*0.8 )) + 1):nrow(DF),]
    train <- DF[1:round(nrow(DF)*0.8),]
    
    # [2-2-3-1] Print dimension of TrainingDF
    dim_train <- paste(dim(train), collapse = " x ")
    print(paste("Dimension of Train-DF:", dim_train))
    
    # [2-2-4] Fit a model [standard settings] on the train DF & take the time:
    start        <- Sys.time()
    curr_forrest <- randomForest(as.factor(gender) ~ ., data = train)
    end          <- Sys.time()
    time_diff    <- difftime(end, start, units = 'mins')
    
    # [2-2-5] Evaluate the Model:
    # OOB:    Get OOB Accuracy of each of the 500 trees and average it!
    curr_OOB_ACC <- mean(1 - curr_forrest$err.rate[,1])
    
    # TESTSET: Get the ErrorRates w/ TestSet
    prediction   <- predict(curr_forrest, test); truth <- as.factor(test[,1])
    res          <- cbind(prediction, truth)
    test_acc     <- sum(res[,1] == res[,2]) / nrow(res)
    
    # [2-2-6] Add results of this block to the Results-DF!
    eval_res_sub_2[nrow(eval_res_sub_2) + 1, ] <- c(df, curr_block, dim_train, 
                                                    curr_OOB_ACC, test_acc, time_diff)
  }
}

# [3] Save the results!
write.csv2(eval_res_sub_2, "./docs/explorative_subsets/performance_RF_0.1percent_cnv_single_block_gender_classif.csv",
           row.names = FALSE)

# [4] Litte analyse w/ results:
# [4-0] Load the results!
eval_res_sub_2 = read.csv2("./docs/explorative_subsets/performance_RF_0.1percent_cnv_single_block_gender_classif.csv",
                           stringsAsFactors = F)

# [4-1] convert to numeric:
cols_ <- c("OOB", "Test", "Time")
eval_res_sub_2[cols_] <- sapply(eval_res_sub_2[cols_], as.numeric)

# [4-2] Get the mean TestAccuracy / OOB-Accuracy of a certain block!
# [4-2-1] TEST
TEST <- sapply(unique(eval_res_sub_2$Block),
               FUN = function(x) summary(eval_res_sub_2$Test[eval_res_sub_2$Block == x]))
# [4-2-2] OOB
OOB <- sapply(unique(eval_res_sub_2$Block),
              FUN = function(x) summary(eval_res_sub_2$OOB[eval_res_sub_2$Block == x]))

# [4-4-3] Apply meaningful names and check the data
colnames(OOB)  <- paste0(colnames(OOB), "_oob")
colnames(TEST) <- paste0(colnames(TEST), "_test")

# [4-5] General Overview of Performance on the single blocks!
cbind(OOB, TEST)
# 10% mirna & mutation + 2.5% rna + 1.25% cnv as joint features to RF        ----
" Subset mirna & mutation blocks by 10% and use 2.5% subsets of 'rna' & 1%
  subset of 'cnv' block, as featurespace to predcit the gender!
"
eval_res_sub_all2 <- data.frame("Data"  = character(), "Dim"   = numeric(),   
                                "OOB"   = numeric(),   "Test"  = numeric(),   
                                "Time"  = numeric(),    stringsAsFactors = F)

# [2] Loop over all DFs w/ gender in their clinical variables!
for (df in DFs_w_gender) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # [2-1] Load 'df' & only keep names of the omics blocks in 'df'!
  omics_blocks <- load(paste0(data_path, df, ".Rda"))
  omics_10     <- c("mirna", "mutation")
  omics_025    <- c("rna")
  omics_0125   <- c("cnv")
  
  # [2-2] Loop over all blocks in 'omics_blocks' of current 'df' to prune them!
  
  # [2-2-1] Blocks that get reduced by 10%
  for (curr_block in omics_10) {
    
    writeLines(paste0("Remove 90% of the features of OmicsBlock:\n", curr_block))
    
    # [2-2-1-1] Subset the current Block and only keep 10% of the original feas
    curr_data    <- eval(as.symbol(curr_block))
    sampled_cols <- sample(ncol(curr_data), round(ncol(curr_data) * 0.1))
    
    # [2-2-1-2] Only keep the sampled columns in the block 
    #          [overwrite current variable w/ reduced dimension var!]
    assign(curr_block, curr_data[ ,sampled_cols])
  }
  
  # [2-2-2] Blocks that get reduced by 2.5%
  for (curr_block in omics_025) {
    
    writeLines(paste0("Remove 97.5% of the features of OmicsBlock:\n", curr_block))
    
    # [2-2-2-1] Subset the current Block and only keep 10% of the original feas
    curr_data    <- eval(as.symbol(curr_block))
    sampled_cols <- sample(ncol(curr_data), round(ncol(curr_data) * 0.025))
    
    # [2-2-2-2] Only keep the sampled columns in the block 
    #          [overwrite current variable w/ reduced dimension var!]
    assign(curr_block, curr_data[ ,sampled_cols])
  }
  
  # [2-2-3] Blocks that get reduced by 2.5%
  for (curr_block in omics_0125) {
    
    writeLines(paste0("Remove 99.9% of the features of OmicsBlock:\n", curr_block))
    
    # [2-2-3-1] Subset the current Block and only keep 10% of the original feas
    curr_data    <- eval(as.symbol(curr_block))
    sampled_cols <- sample(ncol(curr_data), round(ncol(curr_data) * 0.001))
    
    # [2-2-3-2] Only keep the sampled columns in the block 
    #          [overwrite current variable w/ reduced dimension var!]
    assign(curr_block, curr_data[ ,sampled_cols])
  }
  
  # [2-3] Bind gender from the clin block & all blocks!
  # [2-3-1] As rna & cnv share colnames [which leads to an error when fitting RF]
  #         we rename the colnames of the rna DF
  colnames(rna) <- paste0(colnames(rna), "_rna")
  
  # [2-3-2] Paste the subsetted DFs together!
  DF <- cbind(clin$gender, 
              eval(as.symbol(omics_10[1])), 
              eval(as.symbol(omics_10[2])), 
              eval(as.symbol(omics_025[1])),
              eval(as.symbol(omics_0125[1])))
  
  # [2-4] Convert to DF, and name the response variable
  DF <- as.matrix(DF); colnames(DF)[1] <- "gender"
  
  # [2-5] Split to train and test set [80:20] 
  test  <- DF[(round((nrow(DF)*0.8 )) + 1):nrow(DF),]
  train <- DF[1:round(nrow(DF)*0.8),]
  
  # [2-5-1] Print dimension of TrainingDF
  dim_train <- paste(dim(train), collapse = " x ")
  dim_test  <- paste(dim(test), collapse = " x ")
  print(paste("Dimension of Train-DF:", dim_train))
  print(paste("Dimension of Test-DF:", dim_test))
  
  # [2-6] Fit a model [standard settings] on the train DF & take the time:
  start        <- Sys.time()
  curr_forrest <- randomForest(as.factor(gender) ~ ., data = train)
  end          <- Sys.time()
  time_diff    <- difftime(end, start, units = 'mins')
  
  # [2-7] Evaluate the Model:
  # OOB:    Get OOB Accuracy of each of the 500 trees and average it!
  curr_OOB_ACC <- mean(1 - curr_forrest$err.rate[,1])
  
  # TESTSET: Get the ErrorRates w/ TestSet
  prediction   <- predict(curr_forrest, test); truth <- as.factor(test[,1])
  res          <- cbind(prediction, truth)
  test_acc     <- sum(res[,1] == res[,2]) / nrow(res)
  
  # [2-8] Add results of this block to the Results-DF!
  eval_res_sub_all2[nrow(eval_res_sub_all2) + 1, ] <- c(df, dim_train, 
                                                        curr_OOB_ACC, test_acc, 
                                                        time_diff)
}

# [3] Save the results!
write.csv2(eval_res_sub_all2, "./docs/explorative_subsets/performance_RF_10percent_mirna_mutation_2.5percent_rna_0.1percent_cnv_all_blocks_gender_classif.csv",
           row.names = FALSE)


# [4] Litte analyse w/ results:
# [4-0] Load the results!
eval_res_sub_all = read.csv2("./docs/explorative_subsets/performance_RF_10percent_mirna_mutation_2.5percent_rna_0.1percent_cnv_all_blocks_gender_classif.csv ", 
                             stringsAsFactors = F)

# [4-1] convert to correct types:
cols_ <- c("OOB", "Test", "Time")
eval_res_sub_all[cols_] <- sapply(eval_res_sub_all[cols_], as.numeric)
eval_res_sub_all$Data <- as.factor(eval_res_sub_all$Data)

# [4-2] Get average Test/OOB Accuracy  
# [4-2-1] TEST
TEST <- summary(eval_res_sub_all$Test)
# [4-2-2] OOB
OOB <- summary(eval_res_sub_all$OOB)

# [4-3] General Overview of Performance on the single blocks!
cbind(OOB, TEST)
