"
Script to analyse the real multi-omics data set with the single-block approach
"
# [0] Set WD, load packages & define functions
setwd("C:/Users/kuche_000/Desktop/MS_omics/")
library(assertthat)
library(randomForestSRC)

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

# [2] Start the 5-fold CV - Single-Block Approach                           ----
# 2-1 Create a list to save the metrics and define values for the RF when fitted!
single_block_res        <- vector(mode = "list", 
                                  length = length(colnames(missing_str)))
names(single_block_res) <- colnames(missing_str)      

num_trees         = 300
mtry              = NULL
min_node_size     = 5

# 2-2 Impute the missing data for the different TrainSets in the CV
for (i in 1:5) {
  
  # --1 Print Current Fold Status
  print(paste0("FOLD ", i, "/ 5 -------------------------------"))
  
  # --2 Get the current Test- & Train-Set
  index_test  <- index_df[index_df$group == i, "index"]
  index_train <- index_df[index_df$group != i, "index"]
  train_x <- data_merged[index_train, ]
  train_x <- cbind("outcome_Y" = as.factor(index_df[index_train, "outcome"]), 
                   train_x)
  test_x  <- data_merged[index_test, ]
  test_x  <- cbind("outcome_Y" = as.factor(index_df[index_test, "outcome"]), 
                   test_x)
  
  # --3 Loop over all the feature-blocks in 'train_x' and evalute the predicitions
  #     of the single blockwise fitted models!
  for (curr_block in colnames(missing_str)) {
    
    # --1 Get all variables from 'train_x' that belong to the block
    pattern_  <- paste0(curr_block, "_")
    used_feas <- colnames(train_x)[grepl(pattern_, colnames(train_x))]
    
    curr_train <- train_x[,c('outcome_Y', used_feas)]
    
    # --2 Only keep the observations that are observed in 'curr_block'
    curr_train <- curr_train[complete.cases(curr_train),]
    
    # --3 Get predicitons by training a RF on 'curr_train
    # --3-1 Define a Vector to save the predicted classes - fill it with the 
    #       opposite of the true class! Needed later, when predicition on the 
    #       test-set, not all test-obs. can be predicted - these are rated as
    #       wrongly classified for the calculation of the metrics!
    predicted <- sapply(test_x$outcome_Y, FUN = function(x) {
      base::ifelse(x == 1, yes = 0, no = 1)
    })
    
    # --3-2 Check if 'curr_train' has more than 2 obs. - else can not train RF
    #       and therefore not predict on the observations --> preds are NAs
    if (nrow(curr_train) >= 2) {
      
      # - define formula
      curr_formula <- as.formula('outcome_Y ~ .')
      
      # - fit RF on 'curr_train'
      RF <- rfsrc(formula = curr_formula, data = curr_train, 
                  ntree = num_trees, mtry = mtry, nodesize = min_node_size, 
                  samptype = "swr", seed = 12345678, var.used = 'all.trees')
      
      # - get predicitions on the testset - RF can only create preds for obs. 
      #   w/o any NAs in a covariate the model has orginally been trained with!
      #   --> only creates predicitons for CompleteCases! Observations with 
      #       missing values do not recieve a predicition 
      #       --> give them opposite label, so they are rated as wrongly classified
      curr_test   <- test_x[,c('outcome_Y', used_feas)]
      CC_test     <- which(complete.cases(curr_test))
      
      predicitions       <- predict(RF, curr_test)
      predicted[CC_test] <- as.character(predicitions$class)
    }
    
    # --5 Calculate the metrics
    # --5-1 Confusion-Matrix
    confmat <- caret::confusionMatrix(as.factor(predicted), 
                                      as.factor(test_x$outcome_Y))
    
    # --5-2 Are under the ROC Curve
    roc1 <- pROC::auc(pROC::roc(as.numeric(predicted), 
                                as.numeric(as.character(test_x$outcome_Y)),
                                levels = levels(as.factor(as.character(test_x$outcome_Y))),
                                direction = "<"))
    
    roc2 <- pROC::auc(pROC::roc(as.numeric(predicted), 
                                as.numeric(as.character(test_x$outcome_Y)),
                                levels = levels(as.factor(as.character(test_x$outcome_Y))),
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
    
    single_block_res[[curr_block]][[i]] <- res
  }
}
  
  