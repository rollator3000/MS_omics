"Script to get a basic overview over the raw data used in this thesis!
 Data is already preprocessed by Dr. Hornung & was already used the 
 'block-forest' paper - where Dr. Hornung took part! "

# Set the WD, load Packages and define Functions
setwd("C:/Users/kuche_000/Desktop/MS-Thesis/")

get_cols_w_NAs     <- function(df) {
  "
  Check whether the passed dataframe (df) contains any columns w/ missing values
  
  Args:
    - df (data.frame) : Object of class data.frame, with n rows and p columns, 
                        we want to check for NAs in any of its p columns!
    
  Return:
    for each column a boolean, whether it cotains NA Values
  "
  
  results <- apply(df, 2, function(x) any(is.na(x)))
  
  if (length(which(results)) <= 0) {
    return("No missing Values in any Column")
  } else {
    return(which(results))
  }
}
get_rows_w_NAs     <- function(df) {
  "
  Check whether the passed dataframe (df) contains any rows w/ missing values
  
  Args:
    - df (data.frame) : Object of class data.frame, with n rows and p columns, 
                        we want to check for NAs in any of its n rows!
    
  Return:
    Proportion of rows, w/ at least 1 missing Value
  "
  
  results <- apply(df, 1, function(x) any(is.na(x)))
  
  if (length(which(results)) <= 0) {
    return("No missing Values in any of the Observations!")
  } else {
    return(paste0("Proportion of Observations w/ missing Values: ", 
                  length(which(results)) / nrow(df)))
  }
}
get_single_df_info <- function(df) {
  "
  Function to get basic describtive information on the passed data.frame!
  Which features contain missing data, which observations contain missing data,
  what dimension does the data have, what features does it contain and whi types
  are they?
  
  Args:
    - df (data.frame) :  Object of class data.frame  with n rows & p columns we 
                         want to get basic information from!
                         
  Return: 
    - list with:
      - the dimensions of the passed DF
      - the maximal first 10 feature names of the passed DF
      - the columns w/ at least 1 missing value
      - proportion of observations with at least 1 missing value
  "
  # Get dimensions of passed DF
  dim_df <- dim(df)
  
  # Get coltypes (incl. names) of the passed DF! If it has more than 10 columns
  # cut it, and only return coltypes (incl. names) of the first 10 columns!
  if (dim_df[2] > 10) {
    print(paste("Only first 10 cols of", dim_df[2], "in total, are printed"))
    cols_type <- apply(df[,1:10], MARGIN = 2, function(x) typeof(x))
  } else {
    cols_type <- apply(df, MARGIN = 2, function(x) typeof(x))
  }
  
  # Check for features w/ at least 1 NA
  na_cols <- get_cols_w_NAs(df)
  
  # Check for Rows w/ at least 1 NA
  na_rows <- get_rows_w_NAs(df)
  
  # List to save all results and retunrning it!
  result <- list("dimensions"  = dim_df,
                 "col_types"   = cols_type,
                 "Cols w/ NAs" = na_cols,
                 "Rows w/ NAs" = na_rows)
  
  return(result)
}
load_rda_get_infos <- function(rda_path) {
  "
  Function to load a '.Rda' datafile! 
  This can contain single or multiple objects, we will print the names of all
  projects connected with the '.Rda' file! For each 'data.frame' we will call 
  'get_df_info' to get the basic properties and will collect the properties for
  each data.frame/ matrix in a list! 
  The list only has named list elements according to the data.frame/ matrix!
  For none data.frame/ matrix objects we only return the type of the object and
  that we could not do a analysis!
  
  Args:
    - rda_path (string) : path [on top of WD] to a '.Rda' file, that can contain
                          different objects!
  
  Return:
    - list filled with as many named entrances (named by the data.frame/ matrix),
      as the '.Rda'-file contains data.frames / matrices! 
      Each of the lists in the list contains:
        - the dimensions of the passed DF
        - the maximal first 10 feature names of the passed DF
        - the columns w/ at least 1 missing value
        - proportion of observations with at least 1 missing value
  "
  
  # Load '.Rda" data and get all names of the objects with it!
  sub_dfs <- load(rda_path)
  
  print(paste0("The path '", rda_path ,"' lead to a file with >", length(sub_dfs),
               "< objects"))
  print(sub_dfs)
  
  # Define veector to save all describtive results for each DF!
  all_res  <- list()
  all_cols <- 0
  
  # For each of the objects we check whether it is 'data.frame' if they are, we 
  # call 'get_df_info' to get basic properties of the DF!
  for (curr_df in sub_dfs) {
    curr_df_loaded = get(curr_df)
    all_cols <- all_cols + ncol(curr_df_loaded)
    
    if (is.data.frame(curr_df_loaded) | is.matrix(curr_df_loaded)) {
      print(paste0("Analysis of the '", curr_df, "' data.frame"))
      all_res[[curr_df]] <- get_single_df_info(get(curr_df))
      
    } else {
      print(paste(curr_df, "is of type", class(curr_df_loaded)))
    }
  }
  
  print(paste("In Total we have", all_cols, "features in these DFs"))
  
  return(all_res)
}

# Datainspection -- needed fpr further analysis --> always run this         ----
# ----- BLCA DF
BLCA_Res <- load_rda_get_infos("./data/external/Dr_Hornung/original_processed_data/BLCA.Rda")
BLCA_Res
# ----- BRCA DF
BRCA_Res <- load_rda_get_infos("./data/external/Dr_Hornung/original_processed_data/BRCA.Rda")
BRCA_Res
# ----- CESC DF
CESC_Res <- load_rda_get_infos("./data/external/Dr_Hornung/original_processed_data/CESC.Rda")
CESC_Res
# ----- COAD DF
COAD_Res <- load_rda_get_infos("./data/external/Dr_Hornung/original_processed_data/COAD.Rda")
COAD_Res
# ----- ESCA DF
ESCA_Res <- load_rda_get_infos("./data/external/Dr_Hornung/original_processed_data/ESCA.Rda")
ESCA_Res
# ----- GBM DF
GBM_Res <- load_rda_get_infos("./data/external/Dr_Hornung/original_processed_data/GBM.Rda")
GBM_Res
# ----- HNSC DF
HNSC_Res <- load_rda_get_infos("./data/external/Dr_Hornung/original_processed_data/HNSC.Rda")
HNSC_Res
# ----- KIRC DF
KIRC_Res <- load_rda_get_infos("./data/external/Dr_Hornung/original_processed_data/KIRC.Rda")
KIRC_Res
# ----- KIRP DF
KIRP_Res <- load_rda_get_infos("./data/external/Dr_Hornung/original_processed_data/KIRP.Rda")
KIRP_Res
# ----- LGG DF
LGG_Res <- load_rda_get_infos("./data/external/Dr_Hornung/original_processed_data/LGG.Rda")
LGG_Res
# ----- LIHC DF
LIHC_Res <- load_rda_get_infos("./data/external/Dr_Hornung/original_processed_data/LIHC.Rda")
LIHC_Res
# ----- LUAD DF
LUAD_Res <- load_rda_get_infos("./data/external/Dr_Hornung/original_processed_data/LUAD.Rda")
LUAD_Res
# ----- LUSC DF
LUSC_Res <- load_rda_get_infos("./data/external/Dr_Hornung/original_processed_data/LUSC.Rda")
LUSC_Res
# ----- OV DF
OV_Res <- load_rda_get_infos("./data/external/Dr_Hornung/original_processed_data/OV.Rda")
OV_Res
# ----- PAAD DF
PAAD_Res <- load_rda_get_infos("./data/external/Dr_Hornung/original_processed_data/PAAD.Rda")
PAAD_Res
# ----- PRAD DF
PRAD_Res <- load_rda_get_infos("./data/external/Dr_Hornung/original_processed_data/PRAD.Rda")
PRAD_Res
# ----- READ DF
READ_Res <- load_rda_get_infos("./data/external/Dr_Hornung/original_processed_data/READ.Rda")
READ_Res
# ----- SARC DF
SARC_Res <- load_rda_get_infos("./data/external/Dr_Hornung/original_processed_data/SARC.Rda")
SARC_Res
# ----- SKCM DF
SKCM_Res <- load_rda_get_infos("./data/external/Dr_Hornung/original_processed_data/SKCM.Rda")
SKCM_Res
# ----- STAD DF
STAD_Res <- load_rda_get_infos("./data/external/Dr_Hornung/original_processed_data/STAD.Rda")
STAD_Res
# ----- UCEC DF
UCEC_Res <- load_rda_get_infos("./data/external/Dr_Hornung/original_processed_data/UCEC.Rda")
UCEC_Res

# Check DFs for 'gender' as clinical variable & get the distribution        ----
# ! Run the first part, where all 'DF_Res' are created ["Datainspection" Block] !
#   For each of the 21 DFs check whether a variable 'sex' or 'gender' is avaible
#   within the clinical features of the dataset!
BLCA_Res$clin  # --> has gender as clin Variable! --> 3 clinical feas remaining
COAD_Res$clin  # --> has gender as clin Variable! --> 4 clinical feas remaining
ESCA_Res$clin  # --> has gender as clin Variable! --> 3 clinical feas remaining
GBM_Res$clin   # --> has gender as clin Variable! --> 2 clinical feas remaining
HNSC_Res$clin  # --> has gender as clin Variable! --> 4 clinical feas remaining
KIRC_Res$clin  # --> has gender as clin Variable! --> 5 clinical feas remaining
KIRP_Res$clin  # --> has gender as clin Variable! --> 3 clinical feas remaining
LGG_Res$clin   # --> has gender as clin Variable! --> 2 clinical feas remaining
LIHC_Res$clin  # --> has gender as clin Variable! --> 3 clinical feas remaining
LUAD_Res$clin  # --> has gender as clin Variable! --> 5 clinical feas remaining
LUSC_Res$clin  # --> has gender as clin Variable! --> 6 clinical feas remaining
PAAD_Res$clin  # --> has gender as clin Variable! --> 3 clinical feas remaining
READ_Res$clin  # --> has gender as clin Variable! --> 4 clinical feas remaining
SARC_Res$clin  # --> has gender as clin Variable! --> 1 clinical fea remaining
SKCM_Res$clin  # --> has gender as clin Variable! --> 4 clinical feas remaining
STAD_Res$clin  # --> has gender as clin Variable! --> 5 clinical feas remaining

BRCA_Res$clin  # --> doesnt have gender in clin variables!
CESC_Res$clin  # --> doesnt have gender in clin variables!
OV_Res$clin    # --> doesnt have gender in clin variables!
PRAD_Res$clin  # --> doesnt have gender in clin variables!
UCEC_Res$clin  # --> has gender as clin Variable, but with 1 level only!

#   Save the names from the DFs that contain 'gender' in 'clin'-block!
DFs_w_gender <- c("BLCA", "COAD", "GBM", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC",
                  "LGG", "LUAD", "LUSC", "PAAD", "READ", "SARC", "SKCM", "STAD")

#   Check the distribution for all DFs w/ gender!
for (df_curr in DFs_w_gender) {
  load(paste0("./data/external/Dr_Hornung/original_processed_data/", df_curr, ".Rda"))
  print(paste0("Distribution of 'gender' for DF: '", df_curr, "'"))
  print(prop.table(table(factor(clin$gender))))
  print("-------------------------------------------------------")
}

# Get the distribution of the amount of features in each block over all dfs ----

# "DFs_w_gender" needed for the looop is defined in chapter above!
blocks <- c("cnv", "mirna", "mutation", "rna")
for (i in blocks) {
  feas <- c()
  for (df_curr in DFs_w_gender) { 
    feas <- c(feas, get(paste0(df_curr, "_Res"))[[i]]$dimensions[2])
  }
  print(paste0("For the Block: ", i, " -------------------------"))
  print(summary(feas))
}

# Do the same, but subsett the blocks before [same dim in CV later then!]
blocks <- c("cnv", "mirna", "mutation", "rna")
for (i in blocks) {
  feas <- c()
  for (df_curr in DFs_w_gender) { 
    
    # Select the proportion we subset [depending on the selected block!]
    if (i == "cnv")      prop = 0.025
    if (i == "rna")      prop = 0.15
    if (i == "mirna")    prop = 0.5
    if (i == "mutation") prop = 0.1
    
    feas <- c(feas, get(paste0(df_curr, "_Res"))[[i]]$dimensions[2] * prop)
  }
  print(paste0("For the Block: ", i, " -------------------------"))
  print(summary(feas))
}

# Get the type of the features used in all blocks!                          ----
# Names of the usable dataframes (w/ gender in 'clin'-block & 4 omics blocks!)
DFs_w_gender <- c("BLCA", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC","LGG", 
                  "LUAD", "LUSC", "PAAD", "SARC", "SKCM", "STAD")
data_path    <- "./data/external/Dr_Hornung/original_processed_data/"

res_all <- list()

for (df in DFs_w_gender) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # [0] empty list w/ name of the current DF
  assign(df, list())
  curr_res_list <- eval(as.symbol(df))
  
  # [1] Load 'df' & only keep names of the feature blocks!
  omics_blocks <- load(paste0(data_path, df, ".Rda"))
  omics_blocks <- omics_blocks[-which(omics_blocks %in% c("targetvar"))]
  
  # [2] Loop over all the blocks/subDFs
  for (block in omics_blocks) {
    
    writeLines(paste0("Block: --------------------------------------\n", block))
    
    DF    <- eval(as.symbol(block))
    types <- sapply(1:ncol(DF),  FUN = function(x) class(DF[,x]))
    curr_res_list[[length(curr_res_list) + 1]] <- table(types)
  }
  
  res_all[[df]] <- curr_res_list
}

# Get thy types of all DFs and all blocks:
sapply(names(res_all), FUN = function(x) unlist(res_all[[x]]))
lapply(names(res_all), FUN = function(x) unlist(res_all[[x]]))