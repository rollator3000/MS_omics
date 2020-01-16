" Script to crossValidate the adaption of Norbert Krautenbachers's RF Alogrithm!
  For this we learn a seperate RF-Model for each block of features!
  For predicitions on testdata we use the RFs, that use the features avaible in 
  testdata! --> Then to obtain a final predicton we combine all created predicitons!
"
# SetWD and define/load functions
setwd("C:/Users/kuche_000/Desktop/MS-Thesis/")
library(randomForestSRC)
library(checkmate)

load_data_extract_block_names <- function(path = "data/external/Dr_Hornung/Data/ProcessedData_subsets/seed_1234/KIRC_Subset.RData",
                                          seed = 1312, response = 'gender') {
  "Function to load (the already subsetted) data & returning 1 big DataFrame, 
   of all single (already subsetted) blocks!
   
   Args:
    - path (char)     : path to a DF w/ block structure! Shall contain
                        'rna_subset', 'cnv_subset', 'mirna_subset', 'clin' & 'mutation_subset' block!
    - seed (int)      : seed, so the subsetting of single blocks is reproducable
    - response (char) : feature used as reponse class - must be in 'clin' block!
                        and MUST be binary
   
   Return:
    list with 2 entrances:  
      1 - data: all blocks pasted together as single DF, the first col of it 
                equals the response!
      2 - block_names: names of the features in the different  blocks
  "
  # [0] Check Inputs -----------------------------------------------------------
  # 0-1 Load data from path & check whether it has all blocks
  #     if the path is not valid, load() will throw an error!
  load(path)
  if (any(!exists('clin_') & !exists('cnv_sub') & !exists('mirna_sub') & 
          !exists('mutation_sub') & !exists('rna_sub'))) {
    stop("'path' led to a DF with missing block!")
  }
  
  # 0-2 Check that the response is in the clincial block!
  if (!(response %in% colnames(clin_))) stop("Clin Block has no 'response' feature")
  
  # 0-3 Check the 'seed' argument!
  if (seed %% 1 != 0 ) stop("'seed' must be an integer!")
  
  # [1] Create single DF
  # 1-1 Bind the subsetted blocks together to one big DF
  response_       <- clin_[response]
  clin_[response] <- NULL
  df              <- cbind(response_, clin_, cnv_sub, rna_sub, mirna_sub, mutation_sub)
  
  # 1-2 Recode the response as factor!
  df[,colnames(df) == response] <- factor(df[,colnames(df) == response])
  
  # 1-3 Extract the colnames of the single blocks & save them in list:
  block_variables <- list("clin_block"     = colnames(clin_),
                          "cnv_block"      = colnames(cnv_sub),
                          "mirna_block"    = colnames(mirna_sub),
                          "mutation_block" = colnames(mutation_sub),
                          "rna_block"      = colnames(rna_sub))
  
  # [2] Return list with the df & the colnames of the single blocks ------------
  return(list("data" = df,
              "block_names" = block_variables))
}
get_obs_per_fold              <- function(data) {
  "Find the amount of observations needed in each test fold, so the 4 different 
   training folds are equally sized!
   --> This might lead to smaller testfolds than trainingfolds!
       
  Args:
    - data (data.frame) : dataframe on which we want to do CV oon blockwise
                          missingness patterns!
  Return:
    - list filled with:
      -'amount_train':       amount of Observations used for Training in total!
      - 'amount_train_fold': amount of Observations in each Trainfold!
      - 'amount_test':       amount of Observations in each Testfold
  "
  # [0] Check Inputs -----------------------------------------------------------
  assert_data_frame(data, min.rows = 1)
  
  # [1] Get the amount of Obs. in each Trainings- & Testfold, so the Trainig-
  #     folds have the same size!
  #     Amount of trainig observations we use must be dividable by 4! 
  #     Increase the amount of observations in trainingfolds by 1 & 
  #     Reduce the amount of observations in testfold by 1 until, amount of 
  #     observations in trainingfoldsis dividable by 4!
  # 1-1 Split to Train & Test
  amount_train <- floor(4/5 * nrow(data)) 
  
  # 1-2 Count up 'amount_train' until it is dividable by 4!
  while ((amount_train %% 4) != 0) {
    amount_train <- amount_train + 1
  }
  
  # 1-3 Amount of Obs. in Test + in sach train fold!
  amount_train_fold <- amount_train / 4
  amount_test       <- nrow(data) - amount_train
  
  # 1-4 Print info
  writeLines(paste0("From ", nrow(data), " Observations, each TestFold will hold "
                    , amount_test,  " Test-Observations!\n", amount_train, 
                    " Observations will be used for training --> ", 
                    amount_train_fold, " Observations per Trainingsfold"))
  
  # [2] Return the amount of observations needed in each fold, to create equally
  #     sized trainingfolds! ---------------------------------------------------
  return(list("amount_train"      = amount_train,
              "amount_train_fold" = amount_train_fold,
              "amount_test"       = amount_test))
}
mcc_metric                    <- function(conf_matrix) {
  "Function to calculate the MCC Metric
   [MCC = Matthews correlation coefficient]
   --> only for binary cases! If the Conf_Matrix has more than
       2 classes it will return NULL instead of the MCC!
       
    Args: 
      - conf_matrix (confusionMatrix) : Confusion Matrix created with the 
                                        'caret'-Package!
    Return: 
      Matthews correlation coefficient [1 is best, 0 is worst!]
  "
  if (nrow(conf_matrix$table) != 2) {
    warning("Can not calc the MCC-Metric! Return NULL")
    return(NULL)
  }
  
  
  TP <- conf_matrix$table[1,1]
  TN <- conf_matrix$table[2,2]
  FP <- conf_matrix$table[1,2]
  FN <- conf_matrix$table[2,1]
  
  mcc_num <- (TP*TN - FP*FN)
  mcc_den <- 
    as.double((TP+FP))*as.double((TP+FN))*as.double((TN+FP))*as.double((TN+FN))
  
  mcc_final <- mcc_num/sqrt(mcc_den)
  return(mcc_final)
}
do_evaluation_rfsrc           <- function(Forest, testdata, weighted) {
  " Get the aggregated predicition from all trees!
    Evaluate the aggregated predicitons & return metrics!
  
     Args:
      - Forest (list)         : list filled with the objects of class 'rfsrc'
      - testdata (data.frame) : testdata we want to get predicitons for!
                                Must conatain the response we learned the Forest
                                on in the first column!
      - weighted (bool)       : Shall the predicitons from the different RFs
                                be weighted according to their accuracy?!
                                [higher Acc --> higher weight!]
      
     Return:
      - list w/ metrics [accuracy, f1, ...]
  "
  # [0] Check Inputs:
  # 0-1 Reponse Class in the testdata
  assert_data_frame(testdata, min.rows = 1)
  if (!(Forest[[1]]$yvar.names %in% colnames(testdata))) stop("testdata is missing response column!")
  if (colnames(testdata)[1] != Forest[[1]]$yvar.names) stop("'testdata' needs the response column of the trees as first feature")
  
  # 0-2 All Elements of Forest should be of class "rfsrc"
  if (any(sapply(Forest, FUN = function(x) !("rfsrc" %in% class(x))))) {
    stop("not all elements in 'Forest' are of class 'rfsrc'")
  }

  # [1] Remove the Trees, that use split variables, not avaible in testdata!
  # 1-1 Get the feas of the testdata and remove all trees using any feature 
  #     not in the testdata [--> can't do predicitons then!]
  test_cols     <- colnames(testdata)
  forrest_to_rm <- sapply(Forest, FUN = function(x) !any(x$xvar.names %in% test_cols))
  if (any(forrest_to_rm)) Forest <- Forest[-c(which(forrest_to_rm))]
  
  
  # 1-2 Check whether there are any forrests left to do predicitons with 
  #     & if so, print the amount of usable trees!
  if (length(Forest) < 1) stop("Forest can not predicit on TestData, as all 
                                trees use split vars not avaible in 'testdata'")
  
  print(paste(sum(forrest_to_rm), "RFs had to be removed from 'Forest', as these use splitvars, not in 'testdata'"))
  
  # 1-3 Get the OOB-Accuracy of the remaining trees - if weighted is activated
  # 1-3-1 Initial the list!
  weights <- rep(1, times = length(Forest))
  if (weighted) {
    for (i in 1:length(Forest)) {
      weights[i] <- 1 - Forest[[i]]$err.rate[,1][!is.na(Forest[[i]]$err.rate[,1])]
    }
    # normalize the weights!
    weights <- weights / sum(weights)
  }
  
  # [2] Use the remaining trees to create a prediciton
  predicitions <- lapply(Forest, FUN = function(x) predict(x, testdata))
  
  # [3] Aggregate the Predicitons!
  prob_class0 <- c()
  for (j in 1:nrow(testdata)) {
    prob_class0 <- c(prob_class0,
                     weighted.mean(sapply(1:length(Forest), 
                                          FUN = function(x) predicitions[[x]]$predicted[j, 1]),
                                   w = weights, na.rm = TRUE))
  }

  # [4] From Probabilities to classes!
  all_forrest_preds_class <- ifelse(prob_class0 >= 0.5, 0, 1)
  all_forrest_preds_class <- factor(all_forrest_preds_class, 
                                    levels = levels(Forest[[1]]$yvar))
  
  # [5] Get Metrics for the current setting!
  # [5-1] Confusion Matrix
  confmat <- caret::confusionMatrix(data      = all_forrest_preds_class, 
                                    reference = testdata[,1])
  
  # [5-2] Are under the ROC Curve
  roc <- pROC::auc(pROC::roc(testdata[,1], prob_class0, 
                             levels = levels(testdata[,1])))
  
  # [5-3] MCC Matthews correlation coefficient [only for binary cases!]
  mcc <- mcc_metric(conf_matrix = confmat)
  
  
  # [6] Create a list to collect the results!
  res = list("Accuracy"    = confmat$overall["Accuracy"],
             "Sensitifity" = confmat$byClass["Sensitivity"],
             "Specificity" = confmat$byClass["Specificity"],
             "Precision"   = confmat$byClass["Precision"],
             "Recall"      = confmat$byClass["Recall"],
             "F1"          = confmat$byClass["F1"],
             "Balance_Acc" = confmat$byClass["Balanced Accuracy"],
             "AUC"         = as.numeric(roc),
             "MCC"         = mcc)
  
  # 6-1 If the F1-Score/ Precision/ Recall is NA, then we set it to 0 [-1 for MCC]! 
  if (is.na(res$F1))        res$F1        <- 0
  if (is.na(res$Precision)) res$Precision <- 0
  if (is.na(res$Recall))    res$Recall    <- 0
  if (is.na(res$MCC))       res$MCC       <- -1
  
  return(as.vector(res))
}

data_path = "data/external/Dr_Hornung/Data/ProcessedData_subsets/seed_1234/KIRC_Subset.RData"
response = "gender"
seed = 1312
weighted = TRUE
num_trees = as.integer(10)
mtry = NULL
min_node_size = 10

do_CV_NK_setting1             <- function(data_path = "data/external/Dr_Hornung/Data/ProcessedData_subsets/seed_1234/KIRC_Subset.RData",
                                          response = "gender", seed = 1312, weighted = TRUE,
                                          num_trees = as.integer(10), mtry = NULL, 
                                          min_node_size = 10) {
  
  " Function to evaluate RF Adaption on blockwise missing data from Norbert
    Krautenbacher [proposed in it's Ph.D.]
    Data is split into test and train set [5-fold], w/ little adjustment, that
    all train folds have same amount of obs., so the TestFold can be up to 3 obs.
    smaller than the testfold!
    Then each [equally sized] trainingsfold is censored to scenario 1, so that 
    each fold has an observed clinical block & an one observed omics block!
    Then we train a serperate RandomForest on  each of the feature blocks 
    [clin, omicsblock1,...] and ensemble the predicitons from these RFs to a 
    single prediciton and rate these w/ Accuracy, Precision, Specifity, F1-Socre,....
    The TestingSituations are different, as we can test the models on fully 
    observed testdata, on testdata w/ 1 missing block etc. etc.
    
    Args:
      - data_path (char)    : Path to the data, we want to CV! This should lead 
                              to a file w/ multiple sub DFs 
                              [details see 'create_data()']
      - response (char)     : The repsonse we want to model - 
                              MUST be in the 'clin'-block!
      - seed (int)          : Seed to keep results reproducible
      - weighted (bool)     : Shall the predicitons from the different blocks be
                              weighted by the oob accuracy of the blocks?!
      - num_trees (int)     : amount of trees, we shall grow on each[!] fold
      - mtry (int)          : amount of split-variables we try, when looking for 
                              a split variable!
                              If NULL: mtry = sqrt(p)
      - min_node_size (int) : Amount of Observations a node must at least 
                              contain, so the model keeps on trying to split 
                              them!
                              If NULL: min_node_size = 1 for classification...
                              [see ?randomForestSRC()]
    Return:
       - list filled w/:
        * 'res_all' [the CV Results on different testsets]:
            - full    : CV Results for each fold on the fully observed testdata!
            - miss1_A : CV Results for each fold on the testdata, w/ 1 missing 
                        omics-block [cnv]!
            - miss1_B : CV Results for each fold on the testdata, w/ 1 missing 
                        omics-block [rna]!
            - miss1_C : CV Results for each fold on the testdata, w/ 1 missing 
                        omics-block [mutation]!
            - miss1_D : CV Results for each fold on the testdata, w/ 1 missing 
                        omics-block [mirna]!
            - miss2_CD: CV Results for each fold on the testdata, w/ 2 missing 
                        omics-block [mutation & mirna]!
            - miss2_BD: CV Results for each fold on the testdata, w/ 2 missing 
                        omics-block [rna & mirna]!
            - miss2_BC: CV Results for each fold on the testdata, w/ 2 missing 
                        omics-block [rna & mutation]!
            - miss2_AD: CV Results for each fold on the testdata, w/ 2 missing 
                        omics-block [cnv & mirna]!
            - miss2_AC: CV Results for each fold on the testdata, w/ 2 missing 
                        omics-block [cnv & mutation]!
            - miss2_AB: CV Results for each fold on the testdata, w/ 2 missing 
                        omics-block [cnv & rna]!
        * 'settings' [settings used to do the CV - all arguments!]
            - datapath, seed, response, mtry,....   "
  # [0] Check Inputs
  # - data_path, seed, response are all checked within 'create_data()'
  # [0-1] mtry, min_node_size & num_trees are all checked within 'randomForestSRC()'
  
  # [0-2] weighted must be boolean
  assert_logical(weighted)
  
  # [1] Load data & get the names of the different blocks!
  data <- load_data_extract_block_names(path = data_path, seed = seed, 
                                        response = response)
  
  # [2] Get Obs. per fold & print Infos [done by function!]
  obs_per_fold <- get_obs_per_fold(data = data$data)
  
  # [3] Shuffle the data & create lists to save the results in
  # [3-1] Shuffle IDs of data, we use to split data later on!
  set.seed(seed)
  fold_ids <- sample(nrow(data$data), nrow(data$data), replace = FALSE)
  
  # [3-2] Empty lists, to store results of CV!
  full <- list()
  miss1_A <- list(); miss1_B <- list(); miss1_C <- list(); miss1_D <- list()
  miss2_CD <- list(); miss2_BD <- list(); miss2_BC <- list(); miss2_AD <- list()
  miss2_AC <- list(); miss2_AB <- list()
  
  # [4] Start the CV [5-fold per default!]
  for (i in 0:4) {
    
    print(paste0("FOLD: ", as.character(i + 1), "/5 -------------------------"))
    
    # [1] Get TestSet from 'data', by taking the first 'obs_per_fold$amount_test' 
    #     IDs in 'fold_ids'
    test_ids <- fold_ids[((i * obs_per_fold$amount_test) + 1):(((i + 1) * obs_per_fold$amount_test))]
    test_df  <- data$data[test_ids,]
    
    # [2] Get the TrainSet from 'data' [= IDs not in TestSet] & 
    #     induce blockwise missingness!
    train_ids <- fold_ids[-which(fold_ids %in% test_ids)]
    train_df  <- data$data[train_ids,]
    
    # [3] Induce blockwise missingness [SCENARIO_1]
    # [3-1] Sample equally sized 'observed' blocks [according to SCENARIO_1]
    set.seed(seed)
    observed_blocks <- sample(c(rep("Clin, A", obs_per_fold$amount_train_fold), 
                                rep("Clin, B", obs_per_fold$amount_train_fold),
                                rep("Clin, C", obs_per_fold$amount_train_fold), 
                                rep("Clin, D", obs_per_fold$amount_train_fold)),
                              obs_per_fold$amount_train, replace = FALSE)
    
    # [3-2] From the original TrainData censor the blocks, that were not observed!
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
    
    # [4] Fit a seperate RF on each of the blocks seperatly!
    # [4-1] Define the formulas we use to train the RFs on seperate blocks!
    formula1 <- paste(response, "~", 
                      paste(data$block_names$clin_block, collapse = " + "))
    formula2 <- paste(response, "~", 
                      paste(data$block_names$cnv_block, collapse = " + "))
    formula3 <- paste(response, "~", 
                      paste(data$block_names$rna_block, collapse = " + "))
    formula4 <- paste(response, "~", 
                      paste(data$block_names$mutation_block, collapse = " + "))
    formula5 <- paste(response, "~",
                      paste(data$block_names$mirna_block, collapse = " + "))
    
    # [4-2] fit the RFs seperatly on each of the observed blocks!
    RF1 <- rfsrc(formula = as.formula(formula1),
                 data = train_df, ntree = num_trees, mtry = mtry, 
                 nodesize = min_node_size, samptype = "swr")
    
    RF2 <- rfsrc(formula = as.formula(formula2),
                 data = train_df, ntree = num_trees, mtry = mtry, 
                 nodesize = min_node_size, samptype = "swr")
    
    RF3 <- rfsrc(formula = as.formula(formula3),
                 data = train_df, ntree = num_trees, mtry = mtry, 
                 nodesize = min_node_size, samptype = "swr")
    
    RF4 <- rfsrc(formula = as.formula(formula4),
                 data = train_df, ntree = num_trees, mtry = mtry, 
                 nodesize = min_node_size, samptype = "swr")
    
    RF5 <- rfsrc(formula = as.formula(formula5),
                 data = train_df, ntree = num_trees, mtry = mtry, 
                 nodesize = min_node_size, samptype = "swr")
    
    Forest <- list(RF1, RF2, RF3, RF4, RF5)
    
    # [6] Get (unweighted) predicitons on the different testsets!
    # [6-1] Full TestSet
    print("Evaluation full TestSet -------------------------------------------")
    full[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, testdata = test_df) # FULL TESTSET
    
    # [6-2] TestSet, where one of the omics blocks is missing!
    print("Evaluation TestSet w/ 1 missing omics block------------------------")
    miss1_A[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest,
                                            testdata = test_df[,-which(colnames(test_df) %in% data$block_name$cnv_block)])
    
    miss1_B[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, 
                                            testdata = test_df[,-which(colnames(test_df) %in% data$block_name$rna_block)])
    
    miss1_C[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, 
                                            testdata = test_df[,-which(colnames(test_df) %in% data$block_name$mutation_block)])
    
    miss1_D[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, 
                                            testdata = test_df[,-which(colnames(test_df) %in% data$block_name$mirna_block)])
    
    # [6-2] TestSet, where two of the omics blocks is missing!
    print("Evaluation TestSet w/ 2 missing omics blocks-----------------------")
    miss2_CD[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest,
                                             testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$mutation_block,
                                                                                                 data$block_name$mirna_block))])
    
    miss2_BD[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest,
                                             testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$mirna_block,
                                                                                                 data$block_names$rna_block))])
    
    miss2_BC[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, 
                                             testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$mutation_block,
                                                                                                 data$block_name$rna_block))])
    
    miss2_AD[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, 
                                             testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$cnv_block,
                                                                                                 data$block_name$mirna_block))])
    
    miss2_AC[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest,
                                             testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$cnv_block,
                                                                                                 data$block_name$mutation_block))])
    
    miss2_AB[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest,
                                             testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$cnv_block,
                                                                                                 data$block_name$rna_block))])
  }
  
  # [5] Return the results & settings of parameters used to do CV!
  res_all <- list("full" = full,
                  "miss1_A" = miss1_A,
                  "miss1_B" = miss1_B,
                  "miss1_C" = miss1_C,
                  "miss1_D" = miss1_D,
                  "miss2_CD" = miss2_CD,
                  "miss2_BD" = miss2_BD,
                  "miss2_BC" = miss2_BC,
                  "miss2_AD" = miss2_AD,
                  "miss2_AC" = miss2_AC,
                  "miss2_AB" = miss2_AB)
  
  # Collect the Settings, used to do the CV!
  settings <- list("data_path"     = data_path,
                   "response"      = response, 
                   "seed"          = seed, 
                   "num_trees"     = num_trees,
                   "mtry"          = mtry, 
                   "min_node_size" = min_node_size)
  
  # Return both lists!
  return(list("res_all" = res_all, 
              "settings" = settings))
}

# Run a example and check the results!                                       ----
start_time <- Sys.time()
a <- do_CV_NK_setting1(data_path = "./data/external/Dr_Hornung/Data/ProcessedData/KIRC.Rda",
                       response = "gender", seed = 1312, weighted = TRUE,
                       num_trees = as.integer(10), mtry = NULL, 
                       min_node_size = 10)
end_time <- Sys.time()
end_time - start_time # ~60 sek w/ 10 trees and mtry = NULL; min_node_size = 10

sapply(names(a$res_all), FUN = function(x) mean(c(a$res_all[[x]][[1]]$F1, 
                                                  a$res_all[[x]][[2]]$F1, 
                                                  a$res_all[[x]][[3]]$F1, 
                                                  a$res_all[[x]][[4]]$F1,
                                                  a$res_all[[x]][[5]]$F1),
                                                na.rm = TRUE))