"Script to crossValidate the adaption of Roman's RF Alogrithm!
 Here for each fold [set of observations w/ the same observed feature space], we
 train a seperate RF [-> might result in e.g. 4 different RFs]
 For Testing we prune all the seperate RF [or at least the trees inside]
 Pruning: When a tree uses a split variable that is not avaible in the data we
          shall do predicitons with, we use the terminal node, before the tree 
          uses not avaible splitvar., as new terminal node! 
          If the 1. splitvar. is not known the tree can not be used at all!
  
 Then we combine the predicitons from the different RFs to obtain a final 
 prediciton! The Prediciton of a single tree in a RF is the terminal node the 
 test observation ends up in!
 
 > CV: Train RF-Adaption on TrainData w/ blockwise missingness!
       Then test it based on testobs. w/ all feas avaible/ 1 block missing / ...
"
# Load Funcitons, Classes, librarys & set the WD!
setwd("C:/Users/kuche_000/Desktop/MS-Thesis/")
library(pROC)
source("./code/04_simpleRF_adaption.R")

create_data               <- function(path = "./data/external/Dr_Hornung/Data/ProcessedData/KIRC.Rda",
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
get_obs_per_fold          <- function(data) {
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
  
  if ((amount_train %% 4) == 0) {
    amount_train      <- amount_train + 1
    amount_train_fold <- amount_train / 4
    amount_test       <- amount_test - 1
  }
  
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
all_trees_grown_correctly <- function(trees, replace_rf = replace_rf) {
  "Function to check, whether all trees, were grown correctly & if not, we grow
   these trees again, as long, as they are grown correctly!
  
   Args:
      trees (list) : list filled with object of the class 'Tree'! For each 
                     object in the 'trees' list we check, whether it was grown 
                     correctly! If it doesn't have any child nodes, it is grown 
                     again!
      replace_rf (bool) : When growing the trees, shall we draw from the 
                          observations with or without replacement!
                          
   Return: 
      list of trees, where all of these trees were grown correctly
  "
  # [1] Get the entrance of the objects, that miss child node IDs
  amount_trees <- length(trees)
  wrong_trees  <- unlist(lapply(1:amount_trees, FUN = function(x) if (length(trees[[x]]$child_nodeIDs) == 0) x))
  
  # [2] If there are any trees not grown correctly, grow them again until all
  #     of the trees were grown correctly!
  while (length(wrong_trees) > 0) {
    
    # grow the errours trees again
    trees[wrong_trees] <- mclapply(trees[wrong_trees], 
                                   function(x) {
      x$grow(replace = replace_rf)
      x
    }, mc.cores = 1)
    
    # check whether any of the trees is not grown correctly!
    wrong_trees <- unlist(lapply(1:amount_trees, FUN = function(x) if (length(trees[[x]]$child_nodeIDs) == 0) x))
  }
  
  return(trees)
}
copy_forrest              <- function(Forest) {
  "Funciton to copy all trees in forest!
   This is needed, as we change our trees when doing predicitons on data,
   where we need to prune our trees! 
   --> Need to copy the fitted forest, before using it to do predicitons 
       on different testsets!
       
    Args:
      - Forest (list) : list filled with objects of the class 'tree'
                        e.g. [tree1, tree2, tree3, tree4], where tree1, ..., tree4
                             are also a lists filled with objects of class 'tree' 
    Return:
      - Forest (list) : Same object as passed, but deeply copied!
  "
  # [0] Check Arguments
  for (i in 1:length(Forest)) {
   classes_in_list <- sapply(Forest[[i]], function(x) "Tree" %in% class(x))
   if (any(classes_in_list)) {
     msg <- paste("List", i, "in Forest contains at least one object not of class 'Tree'")
     stop(msg)
   }
  }
  
  # [1] Copy the trees!
  for (i in 1:length(Forest)) {
    assign(paste0("treescopy", as.character(i)), lapply(Forest[[i]], FUN = function(x) x$copy()))
  }
  
  Forest_copy <- list(treescopy1, treescopy2, treescopy3, treescopy4)
  
  return(Forest_copy)
}
get_oob_acc               <- function(trees) {
  " Calculate the OOB Accuracy of a list of trees (= partly RF)!
  
    Args: 
      - trees (list) : list filled w/ objects of class 'Tree'
                       For this list of trees we want to extrat the error rate
                       to e.g. weight the predicitons when using multiple
                               trees-lists to create an aggregated prediciton!
    
    Return:
      - average oob-error Accuracy over all trees in 'trees'-list!
  "
  # [0] Check Input:
  #     Make sure 'trees' is a list filled with 'Trees'
  assert_list(trees, min.len = 1, any.missing = FALSE)
  if (any(sapply(trees, function(x) 'Tree' %in% class(x)))) {
    stop("not all elements in 'trees' are of class 'Tree'")
  }
  
  # [1] Calc the OOB ErrorRate
  # [1-1] Get the OOB Predicitons - first as probability than convert to class!
  oob_preds_prob <- lapply(trees, function(x){
    
    # Get the Predicted Probabilites [from each tree]
    preds_prob <- x$predictOOB()
    sapply(1:ncol(preds_prob), function(y) {
      
      # Convert the predicited probabilities to classes
      # In case of ties between the classes - use the first class!
      rownames(preds_prob)[which(preds_prob[,y] == max(preds_prob[,y]))[1]]
    })
  })
  
  # [1-2] Compare predicted Classes with the true classes & get the error rate
  #       every single tree does --> average them so we know how good the 
  #       predicitive power of this 'subRF' is!
  Accuracy <- lapply(1:length(trees), function(x) {
    # get true response and the predicted response!
    true_resp <- trees[[x]]$data$subset(trees[[x]]$oob_sampleIDs  ,1)
    predicted <- factor(oob_preds_prob[[x]], levels = levels(true_resp))
    
    # use them to get a confusionmatrix and extract the Accuracy!
    confmat   <- caret::confusionMatrix(true_resp, predicted)
    confmat$overall["Accuracy"]
  })
  
  return(mean(unlist(Accuracy)))
}
mcc_metric                <- function(conf_matrix) {
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
do_evaluation             <- function(Forest, testdata, weighted) {
  " Get the aggregated predicition from all trees! 
    Evaluate the aggregated predicitons & return metrics!
  
     Args:
      - Forest (list)         : list filled with the objects of class 'Tree'
      - testdata (data.frame) : testdata we want to get predicitons for!
                                Must conatain the response we learned the Forest
                                on!
      - weighted (boolean)    : Shall the predicitons be weighted when aggregted
                                --> Blocks with higher Accuracy, recieve higher
                                    weight!
      
     Return:
      - list w/ metrics [accuracy, f1, mcc, roc, ...]
  "
  # [1] For each TreeSet in the Forest, we prepare the testdata
  #     Setting factors to same levels, ...
  testdata1 <- process_test_data(tree = Forest[[1]][[1]], test_data = testdata)
  testdata2 <- process_test_data(tree = Forest[[2]][[1]], test_data = testdata)
  testdata3 <- process_test_data(tree = Forest[[3]][[1]], test_data = testdata)
  testdata4 <- process_test_data(tree = Forest[[4]][[1]], test_data = testdata)
  
  # [2] Get a prediction for every observation in TestData from all the tree
  tree1_pred <- get_pruned_prediction(trees = Forest[[1]], test_set = testdata1)
  tree2_pred <- get_pruned_prediction(trees = Forest[[2]], test_set = testdata2)
  tree3_pred <- get_pruned_prediction(trees = Forest[[3]], test_set = testdata3)
  tree4_pred <- get_pruned_prediction(trees = Forest[[4]], test_set = testdata4)
  
  # [3] If we want to create weighted ensemble of the predicitons, we need to
  #     calc the OOB-Accuracy and per 'trees' and use these as weights!
  #     [lower ACC --> lower weight!]
  if (weighted) {
    trees1_acc <- round(get_oob_acc(Forest[[1]]), 2)
    trees2_acc <- round(get_oob_acc(Forest[[2]]), 2)
    trees3_acc <- round(get_oob_acc(Forest[[3]]), 2)
    trees4_acc <- round(get_oob_acc(Forest[[4]]), 2)
    tree_weights <- c(trees1_acc, trees2_acc, trees3_acc, trees4_acc) 
  } else {
    tree_weights <- c(1, 1, 1, 1)
  }
  
  # [4] Aggregate Predictions from the different trees!
  # [4-1] Get the probabilities of all test obs to be of class '0'
  all_forrest_preds_probs_class_0 <- sapply(1:nrow(testdata), FUN = function(x) {
    
    # Class1
    prob_class1 <- weighted.mean(c(tree1_pred$Probs[[x]][1], tree2_pred$Probs[[x]][1],
                                   tree3_pred$Probs[[x]][1], tree4_pred$Probs[[x]][1]), 
                                 w = tree_weights, na.rm = TRUE)
    prob_class1
  })
  
  # [4-2] Convert the probabilities to class predicitons and convert it to 'class'
  all_forrest_preds_class <- ifelse(all_forrest_preds_probs_class_0 >= 0.5, 0, 1)
  all_forrest_preds_class <- factor(all_forrest_preds_class, 
                                    levels = levels(Forest[[1]][[1]]$data$data[,1]))
  
  # [5] Get Metrics for the current setting!
  # [5-1] Confusion Matrix - from which we can calc/ extract most metrics
  confmat <- caret::confusionMatrix(data      = all_forrest_preds_class, 
                                    reference = testdata[,1])
  
  # [5-2] Are under the ROC Curve
  roc <- pROC::auc(pROC::roc(testdata[,1], all_forrest_preds_probs_class_0, 
                             levels = levels(Forest[[1]][[1]]$data$data[,1])))
  
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
  
  return(as.vector(res))
}
do_CV_setting1            <- function(data_path = "./data/external/Dr_Hornung/Data/ProcessedData/KIRC.Rda",
                                      response = "gender", seed = 1312, weighted = TRUE,
                                      num_trees = as.integer(10), mtry = NULL, 
                                      min_node_size = NULL, unorderd_factors = "ignore",
                                      replace_rf = TRUE) {
  " Function to evaluate RF Adaption on blockwise missing data!
    Data is split into test and train set [curently fixed to 5-fold], with the 
    little adjustment, the amount of traindata can be split into 4 folds w/o rest
    --> All train folds have same amount of observations
    --> TestFold can be a bit smaller as price for the equally sized training folds!
    
    Then each [equally sized] trainingsfold is censored to scenario 1, so that 
    each fold has an observed clinical block & an observed omics block!
    Then we train a serperate RandomForest on the 4 different training folds 
    [where each fold has different observed features] and ensemble the 
    predicitons from these 4 different RFs to a single prediciton and rate these
    w/ Accuracy, Precision, Specifity, F1-Socre,....
    The TestingSituations are different, as we can test the models on fully 
    observed testdata, on testdata w/ 1 missing block etc. etc.
    
    Args:
      - data_path (char)    : Path to the data, we want to CV! This should lead 
                              to a file w/ multiple sub DFs 
                              [details see 'create_data()']
      - response (chr)      : The repsonse we want to model - 
                              MUST be in the 'clin'-block!
      - seed (int)          : Seed to keep results reproducible
      - weighted (bool)     : Shall the different trees be weighted w/ 
                              OOB Accuracy of the trees
      - num_trees (int)     : amount of trees, we shall grow on each[!] fold
      - mtry (int)          : amount of split-variables we try, when looking for 
                              a split variable! 
                              If 'NULL': mtry = sqrt(p)
      - min_node_size (int) : Amount of Observations a node must at least 
                              contain, so the model keeps on trying to split 
                              them!
                              If 'NULL: Set automatically in 'simpleRF()'
      - unorderd_factors (chr) : How to handle non numeric features!
                                 --> must be in ['ignore', 'order_once', 
                                                 'order_split', 'partition']
      - replace_rf (bool)      : When growing the tree, shall we draw w/ or w/o 
                                 replacement!
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
            - datapath, seed, response, mtry,.... 
  "
  # [0] Check Inputs
  # [0-0] data_path, seed, response are all checked within 'create_data()'
  
  # [0-1] mtry, min_node_size & num_trees are all checked within simpleRF()
  
  # [0-2] replace_rf & weighted must be boolean
  assert_logical(replace_rf, len = 1)
  assert_logical(weighted, len = 1)
  
  # [0-3] unorderd factors must be a legit value
  if (!(unorderd_factors %in% c("ignore", "order_once", "order_split", "partition"))) {
    stop("Unknown value for argument 'unordered_factors'")
  }
  
  # [1] Get the data & dimensions of train & test folds!
  # [1-1] Load the blockwise Omics-Data & create a single DF 
  data <- create_data(path = data_path, seed = seed, response = response)
  
  # [1-2] Get amount of Obs. we need for equally sized train folds 
  obs_per_fold <- get_obs_per_fold(data = data$data)
  
  # [2] Split the data to Test & Train
  # [2-1] Shuffle IDs from 'data' randomly, for splitting it to test & train!
  set.seed(seed)
  fold_ids <- sample(nrow(data$data), nrow(data$data), replace = FALSE)
  
  # [2-2] Create empty lists to store results in!
  full <- list()
  miss1_A <- list(); miss1_B <- list(); miss1_C <- list(); miss1_D <- list()
  miss2_CD <- list(); miss2_BD <- list(); miss2_BC <- list(); miss2_AD <- list()
  miss2_AC <- list(); miss2_AB <- list()
 
  # [3] Start the CV, split data to Test and Train and evaluate it!
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
    
    # [4] Induce blockwise missingness [SCENARIO_1]
    # [4-1] Sample equally sized 'observed' blocks [according to SCENARIO_1]
    set.seed(seed)
    observed_blocks <- sample(c(rep("Clin, A", obs_per_fold$amount_train_fold), 
                                rep("Clin, B", obs_per_fold$amount_train_fold),
                                rep("Clin, C", obs_per_fold$amount_train_fold), 
                                rep("Clin, D", obs_per_fold$amount_train_fold)),
                              obs_per_fold$amount_train, replace = FALSE)
    
    # [4-2] Split Traindata into observed blocks! The resulting blocks will only
    #       contain the features, in the blocks!
    block1 <- train_df[which(observed_blocks == "Clin, A"), 
                       c(response, data$block_names$clin_block, data$block_names$cnv_block)]
    block2 <- train_df[which(observed_blocks == "Clin, B"), 
                       c(response, data$block_names$clin_block, data$block_names$rna_block)]
    block3 <- train_df[which(observed_blocks == "Clin, C"), 
                       c(response, data$block_names$clin_block, data$block_names$mutation_block)]
    block4 <- train_df[which(observed_blocks == "Clin, D"), 
                       c(response, data$block_names$clin_block, data$block_names$mirna_block)]
    
    # [5] Fit 'num_trees' decision trees on each block!
    # [5-1] Get the Formula we use to fit all DecisionTrees/ partial forrests!
    formula_all <- as.formula(paste(response, " ~ ."))
    
    # [5-2] BLOCK1 - grow the trees [as long, as all of them are grown correctly]
    trees1 <- simpleRF(formula = formula_all, data = block1, 
                       num_trees = num_trees, mtry = mtry, 
                       min_node_size = min_node_size, replace = replace_rf,  
                       splitrule = NULL, unordered_factors = unorderd_factors)
    trees1 <- mclapply(trees1, function(x) {
      x$grow(replace = replace_rf)
      x
    }, mc.cores = 1)
    
    # [5-3] BLOCK2 - grow the trees [as long, as all of them are grown correctly]
    trees2 <- simpleRF(formula = formula_all, data = block2, 
                       num_trees = num_trees, mtry = mtry, 
                       min_node_size = min_node_size, replace = replace_rf,  
                       splitrule = NULL, unordered_factors = unorderd_factors)
    trees2 <- mclapply(trees2, function(x) {
      x$grow(replace = replace_rf)
      x
    }, mc.cores = 1)

    # [5-4] BLOCK3 - grow the trees [as long, as all of them are grown correctly]
    trees3 <- simpleRF(formula = formula_all, data = block3, 
                       num_trees = num_trees, mtry = mtry, 
                       min_node_size = min_node_size, replace = replace_rf,  
                       splitrule = NULL, unordered_factors = unorderd_factors)
    trees3 <- mclapply(trees3, function(x) {
      x$grow(replace = replace_rf)
      x
    }, mc.cores = 1)
    
    # [5-5] BLOCK4 - grow the trees [as long, as all of them are grown correctly]
    trees4 <- simpleRF(formula = formula_all, data = block4, 
                       num_trees = num_trees, mtry = mtry, 
                       min_node_size = min_node_size, replace = replace_rf,  
                       splitrule = NULL, unordered_factors = unorderd_factors)
    trees4 <- mclapply(trees4, function(x) {
      x$grow(replace = replace_rf)
      x
    }, mc.cores = 1)
    
    # [6] Check, that all of the trees were grown correctly & create a forest from it!
    trees1 <- all_trees_grown_correctly(trees1, replace_rf = replace_rf)
    trees2 <- all_trees_grown_correctly(trees2, replace_rf = replace_rf)
    trees3 <- all_trees_grown_correctly(trees3, replace_rf = replace_rf)
    trees4 <- all_trees_grown_correctly(trees4, replace_rf = replace_rf)
    
    Forest <- list(trees1, trees2, trees3, trees4)
    
    # [7] Start Testing!
    # [7-1] FULL TESTSET - all blocks observed!
    #       --> copy the forrest, so we don't override original tree!
    print("Evaluation full TestSet -------------------------------------------")
    curr_Forest   <- copy_forrest(Forest)
    full[[i + 1]] <- do_evaluation(Forest = curr_Forest, testdata = test_df, 
                                   weighted = weighted)
    rm(curr_Forest); gc()
    
    # [7-2] TESTSET ONE OMICS BLOCK MISSING - one block is missing in TestData, 
    #       everytime before evaluation we need to copy the Forest, as the 
    #       evaluation can leed to pruned trees [also outside of the function] 
    print("Evaluation TestSet w/ 1 missing omics block------------------------")
    curr_Forest      <- copy_forrest(Forest)
    miss1_A[[i + 1]] <- do_evaluation(Forest   = curr_Forest, weighted = weighted,
                                      testdata = test_df[,-which(colnames(test_df) %in% data$block_name$cnv_block)])
    rm(curr_Forest); gc()
    
    curr_Forest      <- copy_forrest(Forest)
    miss1_B[[i + 1]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                      testdata = test_df[,-which(colnames(test_df) %in% data$block_name$rna_block)])
    rm(curr_Forest); gc()
    
    curr_Forest      <- copy_forrest(Forest)
    miss1_C[[i + 1]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                      testdata = test_df[,-which(colnames(test_df) %in% data$block_name$mutation_block)])
    rm(curr_Forest); gc()
    
    curr_Forest      <- copy_forrest(Forest)
    miss1_D[[i + 1]] <- do_evaluation(Forest = curr_Forest,  weighted = weighted,
                                      testdata =  test_df[,-which(colnames(test_df) %in% data$block_name$mirna_block)])
    rm(curr_Forest); gc()
    
    # [7-3] TESTSET TWO OMICS BLOCKS MISSING - two ommic blocks are missing in
    #       TestData, everytime before evaluation we need to copy the Forest, 
    #       as the  evaluation can leed to pruned trees [also outside of the function] 
    print("Evaluation TestSet w/ 2 missing omics blocks-----------------------")

    curr_Forest       <- copy_forrest(Forest)
    miss2_CD[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                       testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$mutation_block,
                                                                                           data$block_name$mirna_block))])
    rm(curr_Forest); gc()
    
    curr_Forest       <- copy_forrest(Forest)
    miss2_BD[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                       testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$mirna_block,
                                                                                           data$block_names$rna_block))])
    rm(curr_Forest); gc()
    
    curr_Forest       <- copy_forrest(Forest)
    miss2_BC[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted, 
                                       testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$mutation_block,
                                                                                           data$block_name$rna_block))])
    rm(curr_Forest); gc()

    curr_Forest       <- copy_forrest(Forest)
    miss2_AD[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted, 
                                       testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$cnv_block,
                                                                                           data$block_name$mirna_block))])
    rm(curr_Forest); gc()
    
    curr_Forest       <- copy_forrest(Forest)
    miss2_AC[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                       testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$cnv_block,
                                                                                          data$block_name$mutation_block))])
    rm(curr_Forest); gc()
    
    curr_Forest       <- copy_forrest(Forest)
    miss2_AB[[i + 1]] <- do_evaluation(Forest = curr_Forest, weighted = weighted,
                                       testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$cnv_block,
                                                                                           data$block_name$rna_block))])
  }
  
  # Collect all CV Results in a list!
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
                   "weighted"      = weighted,
                   "num_trees"     = num_trees,
                   "mtry"          = mtry, 
                   "min_node_size" = min_node_size,
                   "unorderd_factors" = unorderd_factors,
                   "replace_rf"   = replace_rf)
  
  # Return both lists!
  return(list("res_all"  = res_all, 
              "settings" = settings))
}

# Run a example and check the results!                                       ----
start_time <- Sys.time()
a <- do_CV_setting1(data_path = "./data/external/Dr_Hornung/Data/ProcessedData/KIRC.Rda",
                    response = "gender", seed = 1312, weighted = TRUE,
                    num_trees = as.integer(10), mtry = NULL, 
                    min_node_size = NULL, unorderd_factors = "ignore",
                    replace_rf = TRUE)
end_time <- Sys.time()
end_time - start_time # ~5.5min... w/ 15trees & 10 mtry!
                      # ~6.7min... w/ 10 trees & NULL mtry & KIRC.Rda 
                      # ~ 31mins... w/ 100trees and mtry = NULL

sapply(names(a$res_all), FUN = function(x) mean(c(a$res_all[[x]][[1]]$F1, 
                                                  a$res_all[[x]][[2]]$F1, 
                                                  a$res_all[[x]][[3]]$F1, 
                                                  a$res_all[[x]][[4]]$F1,
                                                  a$res_all[[x]][[5]]$F1),
                                                na.rm = TRUE))

sapply(names(a$res_all), FUN = function(x) mean(c(a$res_all[[x]][[1]]$Accuracy, 
                                                  a$res_all[[x]][[2]]$Accuracy, 
                                                  a$res_all[[x]][[3]]$Accuracy, 
                                                  a$res_all[[x]][[4]]$Accuracy,
                                                  a$res_all[[x]][[5]]$Accuracy),
                                                na.rm = TRUE))