" Script to create DFs needed in the following Scripts!
  
  [1] Create an artififcal dataset w/ factors, characters & numeric features 
      - this was/ is needed for the adaption of the 'simpleRF' Package!
      
  [2] Create Subsets of the original datasets!
      - 10% subset for the mirna & mutation blocks!
      - 5% subset for the rna block!
      - 0.25% ubset for the cnv block!
"
# Load Packages & define Functions!
setwd("C:/Users/kuche_000/Desktop/MS-Thesis")
library(randomForest)

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
data_path <- "./data/external/Dr_Hornung/Data/ProcessedData/"

# 2-2-3 seed for reproductibity
seed <- 1234 # OBACHT: there needs to be a folder 'seed_1234' [or whatever seed set]
             #          in '.data/external/Dr_Hornung/Data/ProcessedData_subsets'

# 2-2 Loop over the DFs and create the subsets!
for (DF in DFs_w_gender) {
  
  # 2-2-1 Load Data and the names of the single blocks 
  #       [except for 'target_var', as we only model classification problems!]
  curr_blocks <- load(paste0(data_path, DF, ".Rda"))
  curr_blocks <- curr_blocks[-which(curr_blocks %in% c("targetvar"))]
  
  # 2-2-2 Subset each of the blocks
  #       mirna - remove 90%
  mirna_ <- eval(as.symbol(curr_blocks[curr_blocks == "mirna"]))
  set.seed(seed)
  mirna_sub <- mirna_[, base::sample(ncol(mirna_), round(ncol(mirna_) * 0.1))]
  
  #       mutation - remove 90%  
  mutation_ <- eval(as.symbol(curr_blocks[curr_blocks == "mutation"]))
  set.seed(seed)
  mutation_sub <- mutation_[, base::sample(ncol(mutation_), round(ncol(mutation_) * 0.1))]
  
  #       rna - remove 95%  
  rna_ <- eval(as.symbol(curr_blocks[curr_blocks == "rna"]))
  set.seed(seed)
  rna_sub <- rna_[, base::sample(ncol(rna_), round(ncol(rna_) * 0.05))]
  colnames(rna_sub) <- paste0(colnames(rna_sub), "_rna") # add '_rna' to colnames
                                                         # else rna & cnv have 
                                                         # same colnames
  
  #       cnv - remove 99.75%  
  cnv_ <- eval(as.symbol(curr_blocks[curr_blocks == "cnv"]))
  set.seed(seed)
  cnv_sub <- cnv_[, base::sample(ncol(cnv_), round(ncol(cnv_) * 0.0025))]
  
  #       clin - no subset
  clin_ <- eval(as.symbol(curr_blocks[curr_blocks == "clin"])) 
  
  # 2-2-4 Print how far the DF was reduced!
  cols_org <- ncol(mirna_) + ncol(mutation_) + ncol(rna_) + ncol(cnv_) + ncol(clin_)
  cols_red <- ncol(cbind(clin_, mirna_sub, mutation_sub, rna_sub, cnv_sub))
  
  print(paste0("DF: '", DF, "'"))
  print(paste0("Reduced from orignally: ", cols_org, " features"))
  print(paste0("To: ", cols_red, " features"))
  
  # 2-2-5 Save the subsetted DF!
  #       'RData' --> load() will load all single blocks!
  save(mirna_sub, mutation_sub, rna_sub, cnv_sub, clin_, 
       file = paste0("./data/external/Dr_Hornung/Data/ProcessedData_subsets/seed_", seed, "/", 
                     DF, "_subset.RData"))
  
  # 2-2-6 Print info about sucessfully saving etc.
  print(paste0("Subset for ", DF, " sucessfully saved to ", 
               paste0("./data/external/Dr_Hornung/Data/ProcessedData_subsets/seed_", seed, "/", 
                      DF, "_subset.RData")))
}

# [3] Get Test& OOB Performance on the subsetted DFs - JOINT BLOCKS          ----
# 3-1 Table to save the results in
eval_res <- data.frame("Data"  = character(),  "Dim"   = numeric(),    
                       "OOB"   = numeric(),    "Test"  = numeric(),    
                       "Time"  = numeric(), # in minutes
                       stringsAsFactors = F)
# 3-2 Define File- & Foldernames!
DFs_w_gender <- c("BLCA", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC","LGG", 
                  "LUAD", "LUSC", "PAAD", "SARC", "SKCM", "STAD")
DFs_w_gender_subset <- paste0(DFs_w_gender, "_subset.RData")
subset_folder <- "./data/external/Dr_Hornung/Data/ProcessedData_subsets/seed_1234/"

# 3-3 Loop over the DFs, fit a RF & get the OOB- & Test-Error JOINT FEATURES
for (df in DFs_w_gender_subset) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # 3-3-1 Load 'df' & paste the single blocks to a single DF!
  omics_blocks <- load(paste0(subset_folder, df))
  DF_all <- cbind(eval(as.symbol(omics_blocks[1])), 
                  eval(as.symbol(omics_blocks[2])), 
                  eval(as.symbol(omics_blocks[3])), 
                  eval(as.symbol(omics_blocks[4])), 
                  eval(as.symbol(omics_blocks[5])))
  
  # 3-3-2 Use 10% as TestSubset!
  test  <- DF_all[(round((nrow(DF_all)*0.8)) + 1):nrow(DF_all),]
  train <- DF_all[1:round(nrow(DF_all)*0.8),]
  
  # 3-3-3 Fit a model [standard settings] on the train DF & take the time:
  start        <- Sys.time()
  curr_forrest <- randomForest(as.factor(gender) ~ ., data = train)
  end          <- Sys.time()
  time_diff    <- difftime(end, start, units = 'mins')
  
  # 3-3-4 Evaluate the Model:
  # OOB:      Get OOB Accuracy of each of the 500 trees and average it!
  curr_OOB_ACC <- mean(1 - curr_forrest$err.rate[,1])
  
  # TESTSET:  Get the ErrorRates w/ TestSet
  prediction   <- predict(curr_forrest, test)
  truth        <- as.factor(test$gender)
  res          <- cbind(prediction, truth)
  test_acc     <- sum(res[,1] == res[,2]) / nrow(res)
  
  # 3-3-5 Add results of this block to the Results-DF!
  dim_train <- paste(dim(DF_all), collapse = " x ")
  eval_res[nrow(eval_res) + 1, ] <- c(df, dim_train, curr_OOB_ACC, 
                                      test_acc, time_diff)
  
}

write.csv2(eval_res, "./docs/performance_final_subsets/joint_blocks_seed1234.csv",
           row.names = FALSE)

# 3-4 Analyse Results
eval_res <- read.csv2("./docs/performance_final_subsets/joint_blocks_seed1234.csv", 
                      stringsAsFactors = F)

# 3-4-1 Convert type of features!
cols_ <- c("OOB", "Test", "Time")
eval_res[cols_] <- sapply(eval_res[cols_], as.numeric)
eval_res$Data <- as.factor(eval_res$Data)

# 3-4-2 Get average Test/OOB Accuracy  
#       TEST
TEST <- summary(eval_res$Test)
#        OOB
OOB <- summary(eval_res$OOB)

# 3-4-3 General Overview of Performance on the single blocks!
cbind(OOB, TEST)

# [4] Get Test& OOB Performance on the subsetted DFs - SINGLE BLOCKS         ----
# 4-1 Table to save the results in
eval_res <- data.frame("Data"  = character(),  "Dim"   = numeric(),
                       "Block" = character(),  "OOB"   = numeric(),    
                       "Test"  = numeric(),    "Time"  = numeric(), # in minutes
                       stringsAsFactors = F)
# 4-2 Define File- & Foldernames!
DFs_w_gender <- c("BLCA", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC","LGG", 
                  "LUAD", "LUSC", "PAAD", "SARC", "SKCM", "STAD")
DFs_w_gender_subset <- paste0(DFs_w_gender, "_subset.RData")
subset_folder <- "./data/external/Dr_Hornung/Data/ProcessedData_subsets/seed_1234/"

# 4-3 Loop over the DFs, fit a RF & get the OOB- & Test-Error SINGLE BLOCK FEAS
for (df in DFs_w_gender_subset) {
  
  writeLines(paste0("Load Dataframe: ----------------------------------\n", df))
  
  # 4-3-1 Load 'df' & get names of the single blocks & extract response!
  omics_blocks <- load(paste0(subset_folder, df))
  gender <- eval(as.symbol("clin_"))
  gender <- gender$gender
  
  # 4-3-2 Loop over all blocks in 'omics_blocks' of current 'df'
  for (curr_block in  omics_blocks) {
    
    writeLines(paste0("Current OmicsBlock: \n", curr_block))
    
    # 4-3-3 Bind gender from the clin block & the current omics block! 
    DF <- cbind(gender, eval(as.symbol(curr_block)))
    
    # 4-3-4 Convert to DF, and name the response variable
    DF <- as.matrix(DF); colnames(DF)[1] <- "gender"
    
    # 4-3-5 Split to train and test set [80:20] 
    test  <- DF[(round((nrow(DF)*0.8)) + 1):nrow(DF),]
    train <- DF[1:round(nrow(DF)*0.8),]
    
    # 4-3-6 Print dimension of TrainingDF
    dim_train <- paste(dim(train), collapse = " x ")
    print(paste("Dimension of Train-DF:", dim_train))
    
    # 4-3-7 Fit a model [standard settings] on the train DF & take the time:
    start        <- Sys.time()
    curr_forrest <- randomForest(as.factor(gender) ~ ., data = train)
    end          <- Sys.time()
    time_diff    <- difftime(end, start, units = 'mins')
    
    # 4-3-8 Evaluate the Model:
    # OOB:      Get OOB Accuracy of each of the 500 trees and average it!
    curr_OOB_ACC <- mean(1 - curr_forrest$err.rate[,1])
    
    # TESTSET:  Get the ErrorRates w/ TestSet
    prediction   <- predict(curr_forrest, test)
    truth        <- as.factor(test[,1])
    res          <- cbind(prediction, truth)
    test_acc     <- sum(res[,1] == res[,2]) / nrow(res)
    
    # 4-3-9 Add results of this block to the Results-DF!
    eval_res[nrow(eval_res) + 1, ] <- c(df, dim_train, curr_block, curr_OOB_ACC, 
                                        test_acc, time_diff)
    
  }
  
  
}

write.csv2(eval_res, "./docs/performance_final_subsets/single_blocks_seed1234.csv",
           row.names = FALSE)

# 4-4 Analyse the Results
eval_res <- read.csv2("./docs/performance_final_subsets/single_blocks_seed1234.csv", 
                      stringsAsFactors = F)

# 4-4-1 Convert type of features!
cols_ <- c("OOB", "Test", "Time")
eval_res[cols_] <- sapply(eval_res[cols_], as.numeric)
eval_res$Data <- as.factor(eval_res$Data)

# 4-4-2 Get average Test/OOB Accuracy  
#       TEST
TEST <- sapply(unique(eval_res$Block),
               FUN = function(x) summary(eval_res$Test[eval_res$Block == x]))
#        OOB
OOB <- sapply(unique(eval_res$Block),
              FUN = function(x) summary(eval_res$OOB[eval_res$Block == x]))

# 4-4-3 Apply meaningful names and check the data
colnames(OOB)  <- paste0(colnames(OOB), "_oob")
colnames(TEST) <- paste0(colnames(TEST), "_test")

# 4-4-4 General Overview of Performance on the single blocks!
cbind(OOB, TEST)