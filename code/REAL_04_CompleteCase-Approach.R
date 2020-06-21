"
Script to analyse the real data with the Complete Case Approach
"
# [0] Set WD, load packages & define functions
setwd("C:/Users/kuche_000/Desktop/MS_omics/")
library(assertthat)
library(randomForestSRC)

detectCores()
registerDoParallel(cores = 2)

mcc_metric        <- function(conf_matrix) {
  "Calculate the MCC Metric [Matthews correlation coefficient]!
   Works only for binary cases! 
   If the Conf_Matrix has more than 2 classes it will return NA 
   instead of the MCC!
       
   Definition of the Metric:
      MCC takes into account true and false positives and negatives and is 
      generally regarded as a balanced measure which can be used even if the 
      classes are of very different sizes.
       
    Args: 
      - conf_matrix (confusionMatrix) : Confusion Matrix created with the 
                                        'caret'-Package!
    Return: 
      - Matthews correlation coefficient (numeric): '1' is best & '-1' is worst
  "
  # [0] Check Inputs -----------------------------------------------------------
  # 0-1 check amount of classes:
  if (nrow(conf_matrix$table) != 2) {
    warning("Can not calc the MCC-Metric! Return NULL")
    return(NA)
  }
  
  # 0-2 Check class of conf_matrix
  if (class(conf_matrix) != "confusionMatrix") {
    stop("conf_matrix not of class 'confusionMatrix'")
  }
  
  # [1] Calc the Score ---------------------------------------------------------
  # 1-1 Get the TruePositives, FalsePositives, FalseNegatives & TrueNegatives
  TP <- conf_matrix$table[1,1]
  TN <- conf_matrix$table[2,2]
  FP <- conf_matrix$table[1,2]
  FN <- conf_matrix$table[2,1]
  
  # 1-2 Calc single parts of MMC [bases on TP, TN, FP, FN]
  mcc_num <- (TP * TN - FP * FN)
  mcc_den <- as.double((TP + FP)) * as.double((TP + FN)) * as.double((TN + FP)) * as.double((TN + FN))
  
  # 1-3 Put together the single parts!
  mcc_final <- mcc_num/sqrt(mcc_den)
  
  # [2] Return  ----------------------------------------------------------------
  return(mcc_final)
}
get_predicition   <- function(train, test, ntree = 300, mtry = NULL, min_node_size = 5) {
  " Function to get predicitions with the imputation approach.
    Based on 'train' a RF is trained and used to predict on 'test'.
    Important:  All observations in 'test' need the same observed features! As
                we only keep the features in the train-set that are also available
                for the test-set. On this reduced train-set a RF is trained!
                
    Args:
      - train (data.frame) : DF with all covariates available that are in test!
                             + does not contain any NA only observations that
                             are completly observed regaring the variables in 
                             'test'!
      - test  (data.frame) : DF where all observations have the same observed
                             features! Must not contain features that are not 
                             available for train!
      - ntree (int)        : Amount of trees to be fit on 'train'!
                              > If NULL: 1000 is the default!
      - min_node_size (int): Amount of Observations a node must at least contain,
                             so the model keeps on trying to split them!  
                              > If NULL: 1 is the default!
      - mtry (int)         : Amount of split-variables we try, when looking for
                             a split variable! 
                              > If 'NULL': mtry = sqrt(p)
                              
    Return:
      - vector with predicted class for each observation in 'test'!
  "
  # [0] Check Inputs  ----------------------------------------------------------
  # 0-1 'train' & 'test' must be dataframes
  assert_data_frame(train, any.missing = FALSE)
  assert_data_frame(test, min.rows = 1)
  
  # 0-2 'test' must not contain any colnames not avaible in train
  if (!all((colnames(test) %in% colnames(train)))) {
    stop("Test-Set has different features than the Train-Set!")
  }
  
  # 0-3 'ntree', 'min_node_size' & 'mtry' must be integer > 0 if NOT NULL
  if (!is.null(ntree)) {
    assert_int(ntree, lower = 10)
  } else {
    ntree = 1000
  }
  if (!is.null(mtry)) assert_int(mtry, lower = 5)
  if (!is.null(min_node_size)) assert_int(min_node_size, lower = 1)
  
  # [1] Train a RF  ------------------------------------------------------------
  # 1-1 Train a RF on 'train'
  # 1-1-1 Extract Response and create the formula the RF to be fit on!
  response    <- colnames(train)[1]
  formula_all <- as.formula(paste(response, " ~ ."))
  
  # 1-1-2 Fit the actual RF
  RF <- rfsrc(formula = formula_all, data = train, 
              ntree = ntree, mtry = mtry, nodesize = min_node_size, 
              samptype = "swr", seed = 12345678, var.used = 'all.trees')
  
  # 1-3 Get Prediciton on the testset from the RF
  predicitons <- predict(RF, test)
  
  return(predicitons$class)
}

# [1] Read in the data and some preprocessing                               ----
# 1-1 Read in the data with Block-Wise missingness & only keep the merged DF!
load("./data/processed/real_data/data 05052020.RData")
rm(df1, df1_mod, df2, df2_mod, df3, df3_mod, 
   df4, df4_mod, df51, df51_mod, df53, df53_mod)

# 1-2 Create DF with information to the outcomes for the CV
#     Completlty copied from Hagenberg!
#         > needed for CV 
index_df <- data.frame(index = 1:521,
                       outcome = y$outcome,
                       group = 0)

#         > weight the obs. such that there is a 50/50 mix of 0 & 1
set.seed(8274)
index_outcome_0 <- sample(rep(1:5, each = 53))
index_outcome_1 <- sample(c(rep(1:5, each = 51), 1))
index_df[index_df$outcome == 0, "group"] <- index_outcome_0
index_df[index_df$outcome == 1, "group"] <- index_outcome_1

# 1-3 Read in the Block-Structure of the data. Tells us which observations
#     have which observed feature blocks
missing_str <- read.csv("./data/processed/real_data/Block Structure.csv", 
                        sep = ";")
missing_str$outcome <- NULL

# [2] Start the 5-fold CV - CompleteCase Approach                           ----
# 2-1 List for the results
CC_res_all <- list()

# 2-2 Impute the missing data for the different TrainSets in the CV
for (i in 1:5) {
  
  # --1 Print Current Fold Status
  print(paste0("FOLD ", i, "/ 5 -------------------------------"))
  
  # --2 Get the current Test- & Train-Set
  index_test  <- index_df[index_df$group == i, "index"]
  index_train <- index_df[index_df$group != i, "index"]
  train_x <- data_merged[index_train, ]
  train_y <- as.factor(index_df[index_train, "outcome"])
  test_x  <- data_merged[index_test, ]
  test_y  <- as.factor(index_df[index_test, "outcome"])
  
  test_data  <- cbind("outcome_Y" = test_y, test_x)
  train_data <- cbind("outcome_Y" = train_y, train_x)
  
  # --3 Get the observed feas per test-obs. - Important, as a modle needs to be
  #     trained on basis of the observations in train, that were observed in the 
  #     exact same features [CompleteCases]
  observed_feas_test <- foreach(x = seq_len(nrow(test_data)), .combine = 'c') %dopar% {
    paste0(which(!(is.na(test_data[x,]))), collapse = "_")
  }
  
  # --3-1 Keep the unique observed test feas [equals the different folds]
  observed_test_folds <- unique(observed_feas_test)
  print(paste0("Found ", length(observed_test_folds), " unique test folds!"))
  
  # --4 For each the 'observed_test_folds' train a RF on basis of the train obs.
  #     that were observed in the exact same covariates!
  # --4-1 Vector to save the results
  predicted    <- c()
  true_reponse <- c()
  for (curr_test in 1:length(observed_test_folds)) {
    
    # --0 Print progress
    print(paste("Evaluation:", curr_test, "/", length(observed_test_folds), "----"))
    
    # --1 Get the current TestFold
    fold_ = observed_test_folds[curr_test]
    
    # --2 Get all Obs. with the feture space as in 'fold_'
    fold_obs_ <- which(observed_feas_test == fold_)
    
    # --3 Get all the indeces of the columns that were observed with this fold!
    obs_columns_ <- as.numeric(unlist(strsplit(fold_, split = "_")))
    
    # --4 Get the current fold - all observations and corresponding features
    curr_fold_test_data <- test_data[fold_obs_, obs_columns_]
    
    # --5 Check for predicitions
    # --5-1 Only the features available for 'curr_fold_test_data'
    curr_train <- train_data[,colnames(curr_fold_test_data)]
    
    # --5-2 Only keep the observations that are completly observed - no NAs
    NAs_per_row <- sapply(1:nrow(curr_train), FUN = function(x) sum(is.na(curr_train[x,])))
    curr_train  <- curr_train[NAs_per_row == 0,]

    # --5-3 Train a RF with 'curr_train' and create predicitions
    #       If the trainset has only 0/ 1 obs. training a RF does not work...
    if (nrow(curr_train) <= 1) {
      curr_preds <- rep(NA, times = nrow(curr_fold_test_data))
    } else {
      curr_preds <- get_predicition(train = curr_train, test = curr_fold_test_data,
                                    ntree = 300, mtry = NULL, min_node_size = 5)
    }
    
    # --5-4 Bind it to the predicted vector
    predicted <- c(predicted, as.character(curr_preds))
    
    # --6 Extract the true response
    true_reponse <- c(true_reponse, 
                      as.character(curr_fold_test_data$outcome_Y))
  }
  
  print(paste(sum(is.na(predicted)), "/ of the test observations could not be predicted"))
  
  # --5 Calculate the metrics
  # --5-1 Get the confusion matrix
  confmat <- caret::confusionMatrix(as.factor(predicted), 
                                    as.factor(true_reponse))
  
  # --5-2 Are under the ROC Curve
  roc1 <- pROC::auc(pROC::roc(as.numeric(predicted), 
                              as.numeric(true_reponse),
                              levels = levels(as.factor(true_reponse)),
                              direction = "<"))
  
  roc2 <- pROC::auc(pROC::roc(as.numeric(predicted), 
                              as.numeric(true_reponse),
                              levels = levels(as.factor(true_reponse)),
                              direction = ">"))
  
  # --5-3 MCC Matthews correlation coefficient [only for binary cases!]
  mcc <- mcc_metric(conf_matrix = confmat)
  
  # --5-4 Collect all metrics in a list & replace the not defined values
  res <- list("Accuracy"    = confmat$overall["Accuracy"],
              "Kappa"       = confmat$overall["Kappa"],
              "Sensitifity" = confmat$byClass["Sensitivity"],
              "Specificity" = confmat$byClass["Specificity"],
              "Precision"   = confmat$byClass["Precision"],
              "Recall"      = confmat$byClass["Recall"],
              "F1"          = confmat$byClass["F1"],
              "Balance_Acc" = confmat$byClass["Balanced Accuracy"],
              "Pos_Pred_Value" =  confmat$byClass["Pos Pred Value"],
              "Neg_Pred_Value" =  confmat$byClass["Neg Pred Value"],
              "Prevalence"  = confmat$byClass["Prevalence"],      
              "AUC1"        = as.numeric(roc1),
              "AUC2"        = as.numeric(roc2),
              "MCC"         = mcc)
  
  #   If the F1-Score/ Precision/ Recall is NA, then we set it to 0 [-1 for MCC]! 
  if (is.na(res$F1))             res$F1             <- 0
  if (is.na(res$Precision))      res$Precision      <- 0
  if (is.na(res$Recall))         res$Recall         <- 0
  if (is.na(res$MCC))            res$MCC            <- -1
  if (is.na(res$Pos_Pred_Value)) res$Pos_Pred_Value <- 0
  if (is.na(res$Neg_Pred_Value)) res$Neg_Pred_Value <- 0
  
  CC_res_all[[i]] <- res
}
  
# 2-3 Save the results of the CV
save(CC_res_all, file = "./docs/CV_Res/REAL/CC_Approach.R")
