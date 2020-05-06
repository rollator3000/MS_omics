"
Script to extract the VariableImportance for the different approaches!
The importance of a variable is calculated according to the 'Permutation' of
the given variable of the out-of-bag examples! For this the RF has only to be
fitted on the training data for this!
"
# [0] Set working directory, load packages and define functions             ----
setwd("C:/Users/kuche_000/Desktop/MS-Thesis/")
library(randomForestSRC)
library(doParallel)

detectCores()
registerDoParallel(cores = 2)

# 0-1 Standard Values for the RandomForest Models
num_trees = 300
mtry = NULL
min_node_size = 5

# [1] Imputation Approach                                                   ----
# 1-1 Load an imputed dataset
load("C:/Users/kuche_000/Desktop/MS_omics/data/processed/RH_subsetted_12345/missingness_1234_imputed/BLCA_IMP_1.RData")

# Fit a random forest model on the imputed data and get the variable importance
s_ <- Sys.time()
formula_  = as.formula("gender ~ .")
fittedmodel <- rfsrc(formula_, data = imputed_data$data[[1]]$train)

importance <- vimp(fittedmodel)$importance
e_ <- Sys.time()
print(e_ - s_)

print(vimp(fittedmodel)$importance)


# [2] Fold-wise approach                                                    ----
# 1-1 Load an Dataset with block-wise missingness
load("C:/Users/kuche_000/Desktop/MS_omics/data/processed/RH_subsetted_12345/missingness_1234/BLCA_1.RData")

# 1-2 Create an empty list to save the results
var_importance_all <- list()

# 1-3 Loop over all Test-Train splits in the data 
for (curr_ind in 1:length(curr_data_1$data)) {
  
  # 1-3-1 Get the Train data
  train <- curr_data_1$data[[curr_ind]]$train
  
  # 1-3-2 Get the Observations that belong to the same fold [same feature space]
  # - Get for each obs. the index of the observed feas
  observed_feas <- foreach(x = seq_len(nrow(train)), .combine = 'c') %dopar% {
    paste0(which(!(is.na(train[x,]))), collapse = "_")
  }
  
  # - Keep the unique observed feas - equals the different folds
  observed_folds <- unique(observed_feas)
  print(paste0("Found ", length(observed_folds), " unique folds!"))
  
  # 1-3-3 Fit a random forest model on each fold & get the variable importance
  importance_list <- list()
  importance_list <- foreach(j_ = 1:length(observed_folds)) %do% {
    
    # - Get the observed features of the current fold
    fold_ = observed_folds[j_]
    
    # - Get all Obs. with the feture space as in 'fold_'
    fold_obs_ <- which(observed_feas == fold_)
    
    # - Get all the indeces of the columns that were observed with this fold!
    obs_columns_ <- as.numeric(unlist(strsplit(fold_, split = "_")))
    
    # - Get all Trainpoints from the obs. w/ same features +  only keep observed 
    # features of these --> fully observed subdata!
    curr_fold_train_data <- train[fold_obs_, obs_columns_]
    
    # 1-3-4 Fit a RF on this fully observed (fold-)subdata!
    # - Define formula
    response    <- colnames(train)[1]
    formula_all <- as.formula(paste(response, " ~ ."))
    
    # - Fit foldwise RF on the current fold!
    fold_RF <- rfsrc(formula_all, data = curr_fold_train_data, 
                     ntree = num_trees, mtry = mtry, nodesize = min_node_size, 
                     samptype = "swr", seed = 12345678, var.used = 'all.trees')
    
    # 1-3-5 Extract the importance of the variables!
    importance <- vimp(fold_RF)$importance
    importance <- as.data.frame(importance)
    importance$Variable <- row.names(importance) 
    
    return(importance)
  }
  
  # 1-4 Define DF to save the Importance for each Fold!
  var_import <- data.frame("Variable" = colnames(train)[2:ncol(train)],
                           "Fold1" = rep(NA, times = (ncol(train) - 1)),
                           "Fold2" = rep(NA, times = (ncol(train) - 1)),
                           "Fold3" = rep(NA, times = (ncol(train) - 1)),
                           "Fold4" = rep(NA, times = (ncol(train) - 1)),
                           stringsAsFactors = FALSE)
  
  # 1-4-1 Save the Variable from each fold - bind per VariableNames!
  for (itt_ in 1:length(importance_list)) {
    matching_rows <- match(importance_list[[itt_]]$Variable, var_import$Variable)
    var_import[matching_rows, (itt_ + 1)] <- importance_list[[itt_]]$all
  }
  
  average_all_folds <- sapply(1:nrow(var_import), 
                              FUN = function(x) mean(c(var_import$Fold1[x], var_import$Fold2[x], 
                                                       var_import$Fold3[x], var_import$Fold4[x]), 
                                                     na.rm = TRUE))
  names(average_all_folds) <- var_import$Variable
  
  var_importance_all[[curr_ind]] <- average_all_folds
}
# [3] Block-wise approach                                                   ----
# 1-1 Load an Dataset with block-wise missingness
load("C:/Users/kuche_000/Desktop/MS_omics/data/processed/RH_subsetted_12345/missingness_1234/BLCA_1.RData")

# 1-2 Create an empty list to save the results
var_importance_all <- list()

# 1-3 Loop over all Test-Train splits in the data 
for (curr_ind in 1:length(curr_data_1$data)) {
  
  # 1-3-1 Get the Train data & the names of the feature-blocks
  train      <- curr_data_1$data[[curr_ind]]$train
  block_vars <- curr_data_1$block_names
  
  # 1-4 Fit a RF on each feature-block
  v_imp <- list()
  
  for (block_ in 1:length(block_vars)) {
    
    # 1-4-1 Get the Observations that have the features of 'block_'
    observed <- sapply(seq_len(nrow(train)), function(j_) {
      sum(is.na(train[j_, block_vars[[block_]]])) == 0
    })
    
    # - Convert Boolean to indicies!
    observed <- which(observed)
    
    # 1-4-2 Fit a RF on this fully observed (fold-)subdata!
    # - Define formula
    response    <- colnames(train)[1]
    formula_all <- as.formula(paste(response, " ~ ."))
    
    # - Get all Obs. from the current block, its block features and the 
    #   corresponding response in a seperate DF!
    curr_fold_train_data <- train[observed, 
                                  c(response, block_vars[[block_]])]
    
    # - Fit a Tree on the block data
    blockwise_rf <- rfsrc(formula = formula_all, data = curr_fold_train_data, 
                          ntree = num_trees, mtry = mtry, nodesize = min_node_size, 
                          samptype = "swr", seed = 12345678, var.used = 'all.trees')
    
    # 1-4-3 Extract the importance of the variables!
    importance <- vimp(blockwise_rf)$importance
    importance <- as.data.frame(importance)
    importance$Variable <- row.names(importance) 
    
    v_imp[[block_]] <- importance
  }
  
  var_importance_all[[curr_ind]] <- v_imp
}