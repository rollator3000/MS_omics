" Norberts Adjustment for the RF in blockwise missing data situations!"


# SetWD and define/load functions
setwd("C:/Users/kuche_000/Desktop/MS-Thesis/")

library(randomForestSRC)

create_data <- function(path = "./data/external/Dr_Hornung/Data/ProcessedData/KIRC.Rda",
                        seed = 1312, response = 'gender') {
  "Function to load data, subset the featurespace - each block will use only 5% 
   of their avaible feas - and returning 1 big DF, of all single (subsetted) 
   blocks!
   
   Args:
    - path (char)     : path to a DF w/ block structure! Shall contain
                        'rna', 'cnv', 'mirna', 'clin' & 'mutation' block!
    - seed (int)      : seed, so the subsetting of single blocks is reproducable
    - response (char) : feature used as reponse class - must be in 'clin' block!
   
   Return:
    list with 2 entrances:  
      1 - data: all blocks pasted together as single DF, the first col of it 
                equals the response!
      2 - block_names: names of the features in the different  blocks
  "
  # [0] Check Inputs
  # [0-1] Load data from path & check it 
  load(path)
  if (any(!exists('clin') & !exists('cnv') & !exists('mirna') & 
          !exists('mutation') & !exists('rna'))) {
    stop("'path' led to a DF with missing block!")
  }
  if (!(response %in% colnames(clin))) stop("Clin Block has no 'response' feature")
  
  # [0-2] Check the 'seed' argument!
  if (seed %% 1 != 0 ) stop("'seed' must be an integer!")
  
  # [1] Subsett the single blocks:
  # [1-1] Extract response from 'clin' block & all other clinical features
  resp     <- clin[which(colnames(clin) == response)]
  clin     <- clin[-which(colnames(clin) == response)]
  
  # [1-2] Subset the 'cnv' features
  set.seed(seed)
  cnv      <- cnv[,sample(ncol(cnv), round(ncol(cnv) * 0.05))]
  
  # [1-3] Subset 'rna' features + rename to avoid duplicates w/ 'cnv'
  set.seed(seed)
  rna           <- rna[,sample(ncol(rna), round(ncol(rna) * 0.05))]
  colnames(rna) <- paste0(colnames(rna), "_rna")
  
  # [1-4] Subset 'mirna' features
  set.seed(seed)
  mirna    <- mirna[,sample(ncol(mirna), round(ncol(mirna) * 0.05))]
  
  # [1-5] Subset 'mutation' features
  set.seed(seed)
  mutation <- mutation[,sample(ncol(mutation), round(ncol(mutation) * 0.05))]
  
  
  # [2] Create single DF
  # [2-1] Bind the subsetted blocks together to one big DF
  df <- cbind(resp, clin, cnv, rna, mirna, mutation)
  
  # [2-2] Recode the response as factor!
  df[,colnames(df) == response] <- factor(df[,colnames(df) == response])
  
  
  # [3] Extract the colnames of the single blocks & save them in list:
  block_variables <- list("clin_block"     = colnames(clin),
                          "cnv_block"      = colnames(cnv),
                          "mirna_block"    = colnames(mirna),
                          "mutation_block" = colnames(mutation),
                          "rna_block"      = colnames(rna))
  
  # [4] Return a list with a 'data' entrance & the colnames of the single blocks
  return(list("data" = df,
              "block_names" = block_variables))
}

get_obs_per_fold <- function(data) {
  "Find the amount of observations, we need in each test fold, so we have same
   sized training folds!
   --> This might lead to slightly smaller testfolds than trainingfolds, but the
       price are equally sized training folds
       
  Args:
    - data (data.frame) : dataframe on which we want to do CV oon blockwise
                          missingness patterns!
  Return:
    - list filled with:
      -'amount_train':       amount of Observations used for Training in total!
      - 'amount_train_fold': amount of Observations in each Trainfold!
      - 'amount_test':       amount of Observations in each Testfold
  "
  # [0] Check Inputs!
  assert_data_frame(data, min.rows = 1)
  
  # [1] Get the amount of Obs. in each Trainingsfold & Testfold, so the Trainig-
  #     folds have the same size!
  #     Amount of trainig observations we use must be dividable by 4! 
  #     Increase the amount of observations in trainingfolds by 1 & 
  #     Reduce the amount of observations in testfold by 1 until, amount of 
  #     observations in trainingfoldsis dividable by 4!
  amount_train <- floor(4/5 * nrow(data))
  amount_test  <- ceiling(1/5 * nrow(data)) 
  
  while ((amount_train %% 4) != 0) {
    amount_train      <- amount_train + 1
    amount_train_fold <- amount_train / 4
    amount_test       <- amount_test - 1
  }
  
  # print info
  writeLines(paste("Each TestFold will hold", amount_test,  "Test-Observations!"))
  writeLines(paste("The remaning", nrow(df) - amount_test, "Obs. are used for Training",
                   "-", amount_train_fold, "Observations per Trainingsfold"))
  
  # [2] Return the amount of observations needed in each fold, to create equally
  #     sized trainingfolds!
  return(list("amount_train"      = amount_train,
              "amount_train_fold" = amount_train_fold,
              "amount_test"       = amount_test))
}

do_evaluation_rfsrc <- function(Forest, testdata) {
  " Get the aggregated predicition from all trees!
    Evaluate the aggregated predicitons & return metrics!
  
     Args:
      - Forest (list)         : list filled with the objects of class 'rfsrc'
      - testdata (data.frame) : testdata we want to get predicitons for!
                                Must conatain the response we learned the Forest
                                on in the first column!
      
     Return:
      - list w/ metrics [accuracy, f1, ...]
  "
  # [0] Check Inputs:
  # [0-1] Reponse Class in the testdata
  assert_data_frame(testdata, min.rows = 1)
  if (!(Forest[[1]]$yvar.names %in% colnames(testdata))) stop("testdata is missing response column!")
  if (colnames(testdata)[1] != Forest[[1]]$yvar.names) stop("'testdata' needs the response column of the trees as first feature")
  
  # [0-2] All Elements of Forest should be of class "rfsrc"
  if (any(sapply(Forest, FUN = function(x) !("rfsrc" %in% class(x))))) {
    stop("not all elements in 'Forest' are of class 'rfsrc'")
  }

  
  # [1] Remove the Trees, that use split variables, not avaible in testdata!
  # [1-1] Get the feas of the testdata and remove all trees using any feature 
  #       not in the testdata [--> can't do predicitons then!]
  test_cols     <- colnames(testdata)
  forrest_to_rm <- sapply(Forest, FUN = function(x) !any(x$xvar.names %in% test_cols))
  if (any(forrest_to_rm)) Forest <- Forest[-c(which(forrest_to_rm))]
  
  
  # [1-2] Check whether there are any forrests left to do predicitons with 
  #       & if so, print the amount of usable trees!
  if (length(Forest) < 1) stop("Forest can not predicit on TestData, as all 
                                trees use split vars not avaible in 'testdata'")
  
  cat(paste(sum(forrest_to_rm), "RFs had to be removed from 'Forest', as they are using splitvariables, not existing in 'testdata'"))
  cat(paste(length(Forest), "RFs are used for prediciton now!"))
  
  
  # [2] Use the remaining trees to create a prediciton
  predicitions <- lapply(Forest, FUN = function(x) predict(x, testdata))
  
  # [3] Aggregate the Predicitons!
  prob_class0 <- c()
  for (j in 1:nrow(testdata)) {
    prob_class0 <- c(prob_class0,
                     mean(sapply(1:length(Forest), 
                               FUN = function(x) predicitions[[x]]$predicted[j, 1])))
  }

  # [4] From Probabilities to classes!
  all_forrest_preds_class <- ifelse(prob_class0 >= 0.5, 0, 1)
  all_forrest_preds_class <- factor(all_forrest_preds_class, 
                                    levels = levels(Forest[[1]]$yvar))
  
  # [5] Get Metrics for the current setting!
  #     Confusion Matrix
  confmat <- caret::confusionMatrix(data      = all_forrest_preds_class, 
                                    reference = testdata[,1])
  
  # [6] Create a list to collect the results!
  res = list("Accuracy"    = confmat$overall["Accuracy"],
             "Sensitifity" = confmat$byClass["Sensitivity"],
             "Specificity" = confmat$byClass["Specificity"],
             "Precision"   = confmat$byClass["Precision"],
             "Recall"      = confmat$byClass["Recall"],
             "F1"          = confmat$byClass["F1"],
             "Balance_Acc" = confmat$byClass["Balanced Accuracy"])
  
  return(as.vector(res))
}

######### IMPLEMENT THE CV FOR THE APPROACH OF NORBERT #########################
seed = 1312
response = "gender"

# [1] Create the data
data <- create_data()
dim(data$data)
length(data$block_names)
names(data$block_names)

# [2] Get Obs. per fold
obs_per_fold <- get_obs_per_fold(data = data$data)

# [3] Split to test and train
set.seed(seed)
fold_ids <- sample(nrow(data$data), nrow(data$data), replace = FALSE)

# [4] Create empty lists to store results in!
full <- list(); miss1_1 <- list(); miss1_2 <- list(); miss1_3 <- list()
miss1_4 <- list(); miss3_1 <- list(); miss3_2 <- list(); miss3_3 <- list()
miss3_4 <- list(); single1 <- list()

# [5] Start the CV, split data to Test and Train and evaluate it!
for (i in 0:4) {
  i = 0
  # [1] Get TestSet from 'data', by taking the first 'obs_per_fold$amount_test' 
  #     IDs in 'fold_ids'
  test_ids <- fold_ids[((i * obs_per_fold$amount_test) + 1):(((i + 1) * obs_per_fold$amount_test))]
  test_df  <- data$data[test_ids,]
  
  # [2] Get the TrainSet from 'data' [= IDs not in TestSet] & 
  #     induce blockwise missingness!
  train_ids <- fold_ids[-which(fold_ids %in% test_ids)]
  train_df  <- data$data[train_ids,]
  
  # [3] Print Infos to size of Test & Train
  print(paste("FOLD", as.character(i)))
  print(paste0("total TRAIN-obs. = ", as.character(nrow(train_df)), " Observations"))
  print(paste0("total TEST-obs. = ", as.character(nrow(test_df)), " Observations"))
  if (ncol(test_df) == ncol(train_df)) {
    print(paste0("total amount of Features = ", as.character(ncol(test_df))))
  } else {
    stop("Test and Trainfold have different amount of features!")
  }
  
  # [4] Induce blockwise missingness [SCENARIO_1]
  # [4-1] Sample equally sized 'observed' blocks [according to SCENARIO_1]
  set.seed(seed)
  observed_blocks <- sample(c(rep("Clin, A", obs_per_fold$amount_train_fold), 
                              rep("Clin, B", obs_per_fold$amount_train_fold),
                              rep("Clin, C", obs_per_fold$amount_train_fold), 
                              rep("Clin, D", obs_per_fold$amount_train_fold)),
                            obs_per_fold$amount_train, replace = FALSE)
  
  # [4-2] From the original TrainData censor the blocks, that were not observed!
  train_df[which(observed_blocks == "Clin, A"),  # CNV + CLIN ONLY 
           c(data$block_names$mirna_block, data$block_names$mutation_block, 
             data$block_names$rna_block)] <- NA
  train_df[which(observed_blocks == "Clin, B"),  # RNA + CLIN 
           c(data$block_names$mirna_block, data$block_names$mutation_block, 
             data$block_names$cnv_block)] <- NA
  train_df[which(observed_blocks == "Clin, C"),  # MUTATION + CLIN
           c(data$block_names$mirna_block, data$block_names$rna_block, 
             data$block_names$cnv_block)] <- NA
  train_df[which(observed_blocks == "Clin, D"),  # MINRA # CLIN
           c(data$block_names$rna_block, data$block_names$mutation_block, 
             data$block_names$cnv_block)] <- NA
  
  # [5] Fit a seperate RF on each of the blocks seperatly!
  # [5-1] Define the formulas we use to train the RFs on seperate blocks!
  formula1 <- paste(response, "~", 
                    paste(data$block_names$clin_block, collapse = " + "))
  formula2 <- paste(response, "~", 
                    paste(data$block_names$cnv_block, collapse = " + "))
  formula3 <- paste(response, "~", 
                    paste(data$block_names$mirna_block, collapse = " + "))
  formula4 <- paste(response, "~", 
                    paste(data$block_names$mutation_block, collapse = " + "))
  formula5 <- paste(response, "~", 
                    paste(data$block_names$rna_block, collapse = " + "))
  
  # [5-2] fit the RFs seperatly on each of the observed blocks!
  RF1 <- rfsrc(formula = as.formula(formula1),
               data = train_df, ntree = 15, mtry = 15, nodesize = 10,
               samptype = "swr")
  
  RF2 <- rfsrc(formula = as.formula(formula2),
               data = train_df, ntree = 15, mtry = 15, nodesize = 10,
               samptype = "swr")
  
  RF3 <- rfsrc(formula = as.formula(formula3),
               data = train_df, ntree = 15, mtry = 15, nodesize = 10,
               samptype = "swr")
  
  RF4 <- rfsrc(formula = as.formula(formula4),
               data = train_df, ntree = 15, mtry = 15, nodesize = 10,
               samptype = "swr")

  RF5 <- rfsrc(formula = as.formula(formula5),
               data = train_df, ntree = 15, mtry = 15, nodesize = 10,
               samptype = "swr")
  
  Forest <- list(RF1, RF2, RF3, RF4, RF5)
  saveRDS(Forest, "./data/interim/tmp_model/tmp_forrest_SRC.rds")
  
  # [6] Get (unweighted) predicitons on the different testsets!
  full[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, testdata = test_df) # FULL TESTSET
  
  miss1_1[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, # NO CNV
                                          testdata = test_df[,-which(colnames(test_df) %in% data$block_name$cnv_block)])
  
  miss1_2[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, # NO MINRA
                                          testdata = test_df[,-which(colnames(test_df) %in% data$block_name$mirna_block)])
  
  miss1_3[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, # NO MUTATION
                                          testdata = test_df[,-which(colnames(test_df) %in% data$block_name$mutation_block)])
  
  miss1_4[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, # NO RNA
                                          testdata = test_df[,-which(colnames(test_df) %in% data$block_name$rna_block)])
  
}

# Return the results!
res_all <- list("full" = full,
                "miss1_1" = miss1_1,
                "miss1_2" = miss1_2,
                "miss1_3" = miss1_3,
                "miss1_4" = miss1_4)

return(res_all)
  