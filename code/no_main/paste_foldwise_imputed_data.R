"Script to paste the single fold-wise imputed data into the whole data"
# Define dataset and setting
df      <- "ESCA"
setting <- "2"

# Create name of the DF
df_name <- paste0(df, "_", setting)

# Load the 'df_name' with the block-wise missingness patterns
curr_DF <- load(paste0("./data/processed/TCGA_subset_12345/missingness_1234/", df_name, ".RData"))
curr_DF <- eval(as.symbol(curr_DF))

# Loop over the the different train folds with missing data and replace them by the imputed data
for (fold_ in 1:length(curr_DF$data)) {
  
  # Load the imputed data for the current fold!
  to_load <- paste0("./data/processed/TCGA_subset_12345/missingness_1234_imputed/", fold_, "_", df, "_", setting, ".RData")
  imputed <- load(to_load)
  imputed <- eval(as.symbol(imputed))
  
  curr_DF$data[[fold_]]$train <- imputed
}

path_to_save <- paste0("./data/processed/TCGA_subset_12345/missingness_1234_imputed/", df, "_IMP_", setting, ".RData") 
save(curr_DF, file = path_to_save)