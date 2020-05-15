" Investigate how good the predicitive peroformance of a RF is on a single 
  block/ joint Blocks for different amount of subsetted features!
  This is nexessary to find out which feature-blocks need to trimmed!

  For each of these Scenarios, we create a DF, that tracks:
    - the dataframe  
    - the block  
    - dim of traindata  
    - OOB & TestSet Accuracy + Test F1 Score! 
    - time in minutes
"
# Set Working Directory and load the needed packages!
setwd("C:/Users/kuche_000/Desktop/MS-Thesis/")
library(randomForest)
library(randomForestSRC)
library(caret)
library(checkmate)
library(ggplot2)
require(gridExtra)
library(reshape2)
library(checkmate)

# Names of the usable dataframes (w/ gender in 'clin'-block & 4 omics blocks!)
DFs_w_gender <- c("BLCA", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC","LGG", 
                  "LUAD", "LUSC", "PAAD", "SARC", "SKCM", "STAD")
data_path    <- "./data/external/TCGA/"

eval_single_block_subsets <- function(DFs_w_gender, seed_to_subset, fraction) {
  "
  Function to evaluate single block performance for different DFs ['DFs_w_gender']. 
  Each DF has 5 different blocks [clin, mirna, cnv, rna & mutation]! 
  For each block [except clin] subset the feature space by only keeping a
  'fraction' of the original covariates & get the performance, to find out 
  which subset to use in our final study!
  For this: 5-fold-CV to each single block and save the results, obtained
  when fitting a standard rfsrc model to it w/ 250 trees!
  
  Args:
    DFs_w_gender (vector) : Vector filled with strings of DFs that have 'gender'
                            within their clinical block! These DF names have to 
                            be in './data/external/TCGA/'
    seed_to_subset (int)  : Seed used, when subsetting the feature space, so it
                            is reproducible!
    fraction (double)     : How much of the original feature space shall be kept?
                            [fraction = 0.2 --> keep 20% and remove 80% of feas
                                                in each block except for clin!]
  Return:
    data.frame with: 'Data'      - data used
                     'Block'     - current block that was used to fit RF
                     'train_dim' - dimension of the data used to fit the RF
                     'OOB_Acc'   - Metric obtained from OOB 
                     'Test_Acc'  - Metric obtained from CV-Testset
                     'Test_F1'   - Metric obtained from CV-Testset
                     'Fold'      - which fold of the 5 folds were used as hold-out
                     'Time'      - how long did it take to fit a RF on the block
                     'Fraction'  - fraction used to subset the single blocks!
                     'subset_seed' -  seed used to create subset!
    will be saved to 'docs/CV/TCGA/explorative_subsets/' the name of the file
    itself is a mixture of seed and fraction!
  "
  # [0] Check Arguments  -------------------------------------------------------
  assert_vector(DFs_w_gender, min.len = 1, unique = TRUE)
  assert_int(seed_to_subset)
  assert_double(fraction, lower = 0, upper = 1)
  
  # [1] Empty DF - for all results of the Evaluation:  -------------------------
  eval_res <- data.frame("Data"      = character(), "Block"    = numeric(), 
                         "train_dim" = numeric(),   "OOB_Acc"  = numeric(),   
                         "Test_Acc"  = numeric(),   "Test_F1"  = numeric(), 
                         "Fold"      = numeric(),   "Time"     = numeric(), # min
                         "Fraction"  = numeric(),   "subset_seed" = numeric(),
                         stringsAsFactors = F)
  
  # 1-1 Fixed Datapath, where we have all our raw Dataframes!
  data_path    <- "./data/external/TCGA/"
  
  # [2] Loop over all DFs w/ gender in their clinical variables!  --------------
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
             paste0("./docs/cv_Res/TCGA/explorative_subsets/", res_name, ".csv"),
             row.names = FALSE) 
}

eval_joint_block_subsets <- function(DFs_w_gender, seed_to_subset, fraction_cnv,
                                     fraction_rna, fraction_mirna,
                                     fraction_mutation, fraction_clin = 1) {
  "
  Function to evaluate joint block performance for different DFs ['DFs_w_gender']. 
  Each DF has 5 different blocks [clin, mirna, cnv, rna & mutation]! 
  For each block choose a subset of features to use from this block, then
  subset each block according to ['fraction_cnv', 'fraction_rna', ...].
  Then bind the subsetted blocks to a single big DF and use this whole DF
  as input for the RF! 
  With this data we do 5-fold-CV to and save the results, obtained
  when fitting a standard rfsrc model to it w/ 250 trees!
  
  Args:
    DFs_w_gender (vector) : Vector filled with strings of DFs that have 'gender'
                            within their clinical block! These DF names have to 
                            be in './data/external/TCGA/'
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
                       'OOB_Acc'   - Metric obtained from OOB
                       'Test_Acc'  - Metric obtained from test-set
                       'Test_F1'   - Metric obtained from test-set
                       'Fold'      - which fold was used as hold-iout fold   
                       'Time'      - how long did it take to fit a RF on the block
                       'Fraction'  - fraction used to subset the single blocks!
      will be saved to 'docs/CV/TCGA/explorative_subsets/' the name of the file
      itself is a mixture of seed and fraction!
  "
  # [0] Check Arguments  -------------------------------------------------------
  assert_vector(DFs_w_gender, min.len = 1, unique = TRUE)
  assert_int(seed_to_subset)
  assert_double(fraction_cnv, lower = 0, upper = 1)
  assert_double(fraction_rna, lower = 0, upper = 1)
  assert_double(fraction_mirna, lower = 0, upper = 1)
  assert_double(fraction_mutation, lower = 0, upper = 1)
  
  # [1] Empty DF - for all results of the Evaluation:  -------------------------
  eval_res <- data.frame("Data"    = character(), "train_dim" = numeric(),
                         "OOB_Acc" = numeric(),   "Test_Acc"  = numeric(),
                         "Test_F1" = numeric(),   "Fold"      = numeric(),
                         "Time"    = numeric(),   "Fraction"  = character(),
                         "subset_seed" = numeric(),
                         stringsAsFactors = F)
  
  # 1-1 Fixed Datapath, where we have all our raw Dataframes!
  data_path    <- "./data/external/TCGA/"
  
  # [2] Loop over all DFs w/ gender in their clinical variables!  --------------
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
             paste0("./docs/cv_Res/TCGA/explorative_subsets/", res_name, ".csv"),
             row.names = FALSE) 
}

# Get Performances w/ single blocks as feature space!                        ----
seed_to_subset = 12345

# Single fully observed Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = seed_to_subset,
                          fraction = 1)

# Single 75% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = seed_to_subset,
                          fraction = 0.75)

# Single 50% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = seed_to_subset,
                          fraction = 0.5)

# Single 25% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = seed_to_subset,
                          fraction = 0.25)

# Single 15% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = seed_to_subset,
                          fraction = 0.15)

# Single 10% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = seed_to_subset,
                          fraction = 0.10)

# Single 5% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = seed_to_subset,
                          fraction = 0.05)

# Single 2.5% subsetted Omics-Blocks as features to RF
eval_single_block_subsets(DFs_w_gender = DFs_w_gender,
                          seed_to_subset = seed_to_subset,
                          fraction = 0.025)

# Get Performances w/ joint subsetted blocks as feature spaces               ----
seed_to_subset = 12345

for (rna_sub in c(0.75, 0.5, 0.25, 0.15, 1)) {
  for (cnv_sub in c(0.005, 0.025, 0.0125, 1)) {
    for (mut_sub in c(1, 0.5, 0.1)) {
      eval_joint_block_subsets(DFs_w_gender = DFs_w_gender, 
                               seed_to_subset = seed_to_subset,
                               fraction_cnv = cnv_sub, 
                               fraction_mirna = 1,
                               fraction_mutation = mut_sub, 
                               fraction_rna = rna_sub,
                               fraction_clin = 1)
      
      eval_joint_block_subsets(DFs_w_gender = DFs_w_gender, 
                               seed_to_subset = seed_to_subset,
                               fraction_cnv = cnv_sub, 
                               fraction_mirna = 0.5,
                               fraction_mutation = mut_sub, 
                               fraction_rna = rna_sub,
                               fraction_clin = 1)
      
      eval_joint_block_subsets(DFs_w_gender = DFs_w_gender, 
                               seed_to_subset = seed_to_subset,
                               fraction_cnv = cnv_sub, 
                               fraction_mirna = 0.1,
                               fraction_mutation = mut_sub, 
                               fraction_rna = rna_sub,
                               fraction_clin = 1)
    }
  }
}
# PLOT the explorative singleblock single block Results                      ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/TCGA/explorative_subsets"

# [1] Load all explorative singleblock Results
# 1-1 get all files in the folder, that are singleblock performances
files <- list.files(data_path)

# 1-1-1 Only keep the explorative 'single-blocks' 
files <- files[grepl("single", files)]

# 1-2 Load the data and store it into a single DF!
DF_all <- data.frame()
for (curr_file in files) {
  DF_curr <-  read.csv2(paste0(data_path, "/", curr_file), stringsAsFactors = F)
  DF_all  <- rbind(DF_all, DF_curr)
}

# 1-3 Convert features to right data type!
str(DF_all)
num_cols <- c("OOB_Acc", "Test_Acc", "Test_F1", "Fraction")
DF_all[,num_cols] <- sapply(num_cols, function(x) as.numeric(DF_all[,x]))

# 1-4 Reshape the layout of data for the plot! 
plot_df <- melt(DF_all, id.vars = c("Block", "Fraction", "subset_seed", "Data"), 
                measure.vars = c("Test_Acc", "Test_F1"))

# 1-4-1 Remove the clincal block as it is not subsetted
plot_df       <- plot_df[-which(plot_df$Block == "clin"),]

# 1-4-2 Convert 'Fraction' to a factor variable and code it in '%'
plot_df$Fraction <- paste0(plot_df$Fraction * 100, "%")
plot_df$Fraction <- factor(plot_df$Fraction, 
                           levels = c("2.5%", "5%", "10%", "15%", "25%", 
                                      "50%", "75%", "100%"))

# 1-4-3 Adjust the names of the single Blocks, such that the spelling is the 
#       same as in the thesis
plot_df$Block[grep("cnv", plot_df$Block)]      <- "CNV"
plot_df$Block[grep("mirna", plot_df$Block)]    <- "miRNA"
plot_df$Block[grep("mutation", plot_df$Block)] <- "Mutation"
plot_df$Block[grep("rna", plot_df$Block)]      <- "RNA"


# 1-5 Plot the performance with the accuracy and the f1-score measure split by
#     the fraction we've used to subset the single blocks!
# 1-5-1 Accuracy
plot_acc <- plot_df[plot_df$variable %in% c("Test_Acc"),]

ggplot(data = plot_acc, aes(x = Fraction, y = value)) + 
  geom_boxplot(fill = "palegreen3") + 
  facet_grid(. ~ Block) +
  theme_bw() +
  ggtitle("Single Block Performance of a random forest over all 14 TCGA data sets",
          subtitle = "Split by the different feature-blocks - Evaluated with 5-fold CV") +
  xlab("Subset of used features") +
  ylab("Accuracy") +
  theme(axis.text.x = element_text(angle = 45 , hjust = 1),
        text = element_text(size = 28))

# 1-5-2 F-1-Score
plot_f1  <- plot_df[plot_df$variable %in% c("Test_F1"),]

ggplot(data = plot_f1, aes(x = Fraction, y = value)) + 
  geom_boxplot(fill = "palegreen3") + 
  facet_grid(. ~ Block) +
  theme_bw() +
  ggtitle("Single Block Performance of a random forest over all 14 TCGA data sets",
          subtitle = "Split by the different feature-blocks - Evaluated with 5-fold CV") +
  xlab("Subset of used features") +
  ylab("F-1-Score") +
  theme(axis.text.x = element_text(angle = 45 , hjust = 1),
        text = element_text(size = 28))

# PLOT the explorative joint block block Results                             ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/TCGA/explorative_subsets"

# [1] Load all explorative singleblock Results
# 1-1 get all files in the folder, that are singleblock performances
files <- list.files(data_path)
files <- files[grepl("joint", files)]

# 1-2 Add the amount of subsets to each block:
DF_all <- data.frame()
for (curr_file in files) {
  DF_curr <-  read.csv2(paste0(data_path, "/", curr_file), stringsAsFactors = F)
  DF_all  <- rbind(DF_all, DF_curr)
}

# 1-4 reshape DF_all for the plot!
plot_df <- melt(DF_all, id.vars = c("Data", "Fraction"), measure.vars = c("Test_Acc", "Test_F1", "Fold"))
plot_df <- plot_df[plot_df$variable %in% c("Test_Acc", "Test_F1"),]
plot_df$value <- as.numeric(plot_df$value)

# 1-5 Add the fraction used of the single blocks as meta data!
for (i in seq_len(nrow(plot_df))) {
  plot_df$mirna_subset[i]    <- strsplit(strsplit(plot_df$Fraction[i], split = "mirna_")[[1]][2], split = "__")[[1]][1]
  plot_df$mutation_subset[i] <- strsplit(plot_df$Fraction[i], split = "mutation_")[[1]][2]
  plot_df$rna_subset[i]      <- strsplit(strsplit(plot_df$Fraction[i], split = "_rna_")[[1]][2], split = "_")[[1]][1]
  plot_df$cnv_subset[i]      <- strsplit(strsplit(plot_df$Fraction[i], split = "_cnv_")[[1]][2], split = "_")[[1]][1]
  plot_df$Fraction_new[i]    <- strsplit(plot_df$Fraction[i], split = "__mirna_")[[1]][2]
  plot_df$Fraction_new[i]    <- paste0("mirna_", plot_df$Fraction_new[i])
}

plot_df$rna_subset_plot <- sapply(plot_df$rna_subset, function(x) paste0("rna_", x)) 
plot_df$cnv_subset_plot <- sapply(plot_df$cnv_subset, function(x) paste0("cnv_", x)) 

# [2] Do the plot, split by subsets 
ggplot(data = plot_df, aes(x = Fraction_new , y = value, fill = variable)) +
  geom_boxplot() + 
  theme_bw() +
  facet_grid(rna_subset_plot ~ cnv_subset_plot) +
  ggtitle("EXPLORATIVE - Joint Block Performance on all 14 DFs - w/ seed 12345", 
          subtitle = "Single Blocks were subsetted! y-axis RNA splits || x-axis CNV splits") +
  xlab("subsets for mirna & mutation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 17))