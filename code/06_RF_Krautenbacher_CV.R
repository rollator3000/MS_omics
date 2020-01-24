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
                                          response = 'gender') {
  "Function to load (the already subsetted) data & returning 1 big DataFrame, 
   of all single blocks incl. the colnames to each block!
   
   Args:
    - path (char)     : path to a DF w/ block wise structure! 
                        Shall contain 'rna_subset', 'cnv_subset', 'mirna_subset',
                        'clin' & 'mutation_subset' block!
    - response (char) : feature used as reponse class - must be in 'clin' block!
                        + MUST be binary - else it will throw an error!
   
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
    stop("'path' led to a DF with at least one missing block!")
  }
  
  # 0-2 Check that the response is in the clincial block!
  if (!(response %in% colnames(clin_))) stop("Clin Block has no 'response' feature")
  
  # 0-3 Check that response is binary!
  if (length(levels(as.factor(clin_[response][,1]))) != 2) {
    stop("The selected Response doesn't have 2 levels! 
         Please choose a binary response from the 'clin'-block!")
  }
  
  # [1] Create single DF -------------------------------------------------------
  # 1-1 Extract the response from the clinical block & rm tehe col from the df
  response_       <- clin_[response]
  clin_[response] <- NULL
  
  # 1-2 Bind the single blocks to a big DF!
  df <- cbind(response_, clin_, cnv_sub, rna_sub, mirna_sub, mutation_sub)
  
  # 1-3 Recode the response as factor!
  df[,colnames(df) == response] <- factor(df[,colnames(df) == response])
  
  # 1-4 Extract the colnames of the single blocks & save them in list:
  block_variables <- list("clin_block"     = colnames(clin_),
                          "cnv_block"      = colnames(cnv_sub),
                          "rna_block"      = colnames(rna_sub),
                          "mutation_block" = colnames(mutation_sub),
                          "mirna_block"    = colnames(mirna_sub))
  
  # [2] Return list with the df & the colnames of the single blocks ------------
  return(list("data" = df,
              "block_names" = block_variables))
}
get_obs_per_fold              <- function(data) {
  "Find the amount of observations needed in each test fold, so the 4 different 
   training folds are equally sized!
   --> This might lead to smaller testfolds than trainingfolds!
       
  Args:
    - data (data.frame) : dataframe on which we want to do CV on blockwise
                          missingness patterns!
  Return:
    - list filled with:
      - 'amount_train':      amount of Observations used for Training in total!
      - 'amount_train_fold': amount of Observations in each Trainfold!
      - 'amount_test':       amount of Observations in each Testfold
  "
  # [0] Check Inputs -----------------------------------------------------------
  assert_data_frame(data, min.rows = 1)
  
  # [1] Calculat the amount of Obs. in Test & Train Fold -----------------------
  #     Get the amount of Obs. in each Trainings- & Testfold, so the 
  #     Trainigfolds have the same size!
  #     Amount of trainig observations must be dividable by 4, so all folds
  #     are equally sized! If this is not true increase the amount of obs. in 
  #     trainingfolds by 1, as long as it is dividable by 4
  #     The remaining obs. are used as TestSet
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
  
  # [2] Return Train- / Test-fold sizes ----------------------------------------
  return(list("amount_train"      = amount_train,
              "amount_train_fold" = amount_train_fold,
              "amount_test"       = amount_test))
}
mcc_metric                    <- function(conf_matrix) {
  "Function to calculate the MCC [Matthews correlation coefficient] Metric
   --> only for binary cases! If the Conf_Matrix has more than 2 classes 
       it will return NULL instead of the MCC!
       
    Definition of the Metric:
      MCC takes into account true and false positives and negatives and is 
      generally regarded as a balanced measure which can be used even if the 
      classes are of very different sizes.
       
    Args: 
      - conf_matrix (confusionMatrix) : Confusion Matrix created with the 
                                        'caret'-Package!
    Return: 
      Matthews correlation coefficient [1 is best, -1 is worst!]
  "
  # [0] Check Inputs -----------------------------------------------------------
  # 0-1 check amount of classes:
  if (nrow(conf_matrix$table) != 2) {
    warning("Can not calc the MCC-Metric! Return NULL")
    return(NULL)
  }
  
  # 0-2 Check class of conf_matrix
  if (class(conf_matrix) != "confusionMatrix") {
    stop("conf_matrix not of class 'confusionMatrix'")
  }
  
  # [1] Calc the Score ---------------------------------------------------------
  TP <- conf_matrix$table[1,1]
  TN <- conf_matrix$table[2,2]
  FP <- conf_matrix$table[1,2]
  FN <- conf_matrix$table[2,1]
  
  mcc_num <- (TP*TN - FP*FN)
  mcc_den <- 
    as.double((TP+FP))*as.double((TP+FN))*as.double((TN+FP))*as.double((TN+FN))
  
  mcc_final <- mcc_num/sqrt(mcc_den)
  
  # [2] Return it --------------------------------------------------------------
  return(mcc_final)
}
do_evaluation_rfsrc           <- function(Forest, testdata, weights) {
  " Get the aggregated predicition from all trees!
    Evaluate the aggregated predicitons & return metrics!
  
     Args:
      - Forest (list)         : list filled with the objects of class 'rfsrc'
      - testdata (data.frame) : testdata we want to get predicitons for!
                                Must conatain the response we learned the Forest
                                on in the first column!
      - weights (vector)      : The weights for the different BlockWise RFs!
                                As some of the blockwise RFs (in Forest) have a 
                                higher Accuracy/ F1/ ... than others, it means
                                that they have a higher predictive power!
                                  --> give weights to them!
                                  --> if we want equally weighted blocks, pass
                                      a vector filled w/ '1'
                                --> needs same length as Forest!
      
     Return:
      - list w/ metrics [accuracy, f1, ...]
  "
  # [0] Check Inputs -----------------------------------------------------------
  # 0-1 Reponse Class in the testdata
  assert_data_frame(testdata, min.rows = 1)
  if (!(Forest[[1]]$yvar.names %in% colnames(testdata))) stop("testdata is missing response column!")
  if (colnames(testdata)[1] != Forest[[1]]$yvar.names) stop("'testdata' needs the response column of the trees as first feature")
  
  # 0-2 All Elements of Forest should be of class "rfsrc"
  if (any(sapply(Forest, FUN = function(x) !("rfsrc" %in% class(x))))) {
    stop("not all elements in 'Forest' are of class 'rfsrc'")
  }
  
  # 0-3 weights should be vector of exactly the length of 'Forest'
  assert_vector(weights, len = length(Forest))
  if (any(sapply(weights, function(x) !is.numeric(x)))) {
    print("'weights:\n")
    print(weights)
    stop("'weights' is not only filled with numerics!")
  }

  # [1] Remove RFs, that use split variables not avaible in the testdata!  ----
  # 1-1 Get the feas of the testdata and remove all RFs using any feature 
  #     not in the testdata [--> can't do predicitons then!]
  #     ---> remove the corresponding weights to the RFs we removed
  test_cols     <- colnames(testdata)
  forrest_to_rm <- sapply(Forest, FUN = function(x) any(!(x$xvar.names %in% test_cols)))
  if (any(forrest_to_rm)) {
    Forest  <- Forest[-c(which(forrest_to_rm))]
    weights <- weights[-c(which(forrest_to_rm))]
  } 
  
  
  # 1-2 Check whether there are any forrests left to do predicitons with 
  #     & if so, print the amount of usable trees!
  if (length(Forest) < 1) stop("Forest can not predicit on TestData, as all 
                                trees use split vars not avaible in 'testdata'")
  
  print(paste0(sum(forrest_to_rm), "/5 blockwise RFs had to be removed from 'Forest', as these use splitvars, not in 'testdata'"))
  
  # [2] Use the remaining trees to create a prediciton -------------------------
  predicitions <- lapply(Forest, FUN = function(x) predict(x, testdata))
  
  # [3] Aggregate the Predicitons! ---------------------------------------------
  prob_class0 <- c()
  for (j in 1:nrow(testdata)) {
    prob_class0 <- c(prob_class0,
                     weighted.mean(sapply(1:length(Forest), 
                                          FUN = function(x) predicitions[[x]]$predicted[j, 1]),
                                   w = weights, na.rm = TRUE))
  }

  # [4] From Probabilities to classes! -----------------------------------------
  all_forrest_preds_class <- ifelse(prob_class0 >= 0.5, 0, 1)
  all_forrest_preds_class <- factor(all_forrest_preds_class, 
                                    levels = levels(Forest[[1]]$yvar))
  
  # [5] Get Metrics for the current setting!
  # [5-1] Confusion Matrix
  confmat <- caret::confusionMatrix(data      = all_forrest_preds_class, 
                                    reference = testdata[,1])
  
  # [5-2] Are under the ROC Curve
  #     Select 'direction = "auto"' and it will choose the class to be positive
  #     automatically!
  roc <- pROC::auc(pROC::roc(testdata[,1], prob_class0, 
                             levels = levels(testdata[,1]), 
                             direction = "auto"))
  
  # [5-3] MCC Matthews correlation coefficient [only for binary cases!]
  mcc <- mcc_metric(conf_matrix = confmat)
  
  # [6] Create a list to collect the results!
  res <- list("Accuracy"    = confmat$overall["Accuracy"],
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
do_CV_NK_setting1             <- function(data_path = "data/external/Dr_Hornung/Data/ProcessedData_subsets/seed_1234/KIRC_Subset.RData",
                                          response = "gender", seed = 1312, 
                                          weighted = TRUE, weight_metric = NULL,
                                          num_trees = as.integer(10), mtry = NULL, 
                                          min_node_size = 10) {
  
  " Function to evaluate RF Adaption on blockwise missing data from Norbert
    Krautenbacher [proposed in it's Ph.D.] Approach!
    
    Data is split into test and train set [5-fold], w/ little adjustment, that
    all train folds have same amount of obs., so the TestFold can be a bit
    smaller than the testfold!
    Then each [equally sized] trainingsfold is censored to scenario 1, so that 
    each fold has an observed clinical block & an one observed omics block!
    
    Then we train a serperate RandomForest on  each of the feature blocks 
    [clin, omicsblock1,...] and ensemble the predicitons from these RFs to a 
    single prediciton and rate these w/ Accuracy, Precision, F1-Socre,....
    The TestingSituations are different, as we can test the models on fully 
    observed testdata, on testdata w/ 1 missing block etc. etc.
    
    Args:
      - data_path (char)    : Path to the data, we want to CV! This should lead 
                              to a file w/ multiple sub DFs 
                              [details see 'create_data()']
      - response (char)     : The repsonse we want to model - 
                              MUST be in the 'clin'-block!
      - seed (int)          : Needed for shuffeling the rows of data before CV,
                              & the assignment of obs to folds[diff folds diff feas] 
                              in a reproducable way! Aswell needed so the RFs on
                              the different blocks to train in a reproducable way!
      - weighted (bool)     : Shall the predicitons from the different blocks be
                              weighted by the oob accuracy of the blocks?!
      - weight_metric (chr) : Which Metric shall be used to weight the different 
                              RFs - the better the single blocks according to the
                              metric, the higher their weight when ensembling!
                              Must be in c('Accuracy', 'F1')
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
            - miss3_ABC: CV Results for each fold on the testdata, w/ 3 missing 
                         omics-block [cnv & rna & mutation]!
            - miss3_ACD: CV Results for each fold on the testdata, w/ 3 missing 
                         omics-block [cnv & mutation & mirna]!
            - miss3_ABD: CV Results for each fold on the testdata, w/ 3 missing 
                         omics-block [cnv & rna & mirna]!
            - miss3_BCD: CV Results for each fold on the testdata w/ 3 missing 
                         omics -bloc [rna & mutation & mirna]
            - single_CL: CV-Results for each fold on the testdata w/ only 
                         clinical features
            - single_A:  CV-Results for each fold on the testdata w/ only 
                         block A features [CNV]
            - single_B:  CV-Results for each fold on the testdata w/ only 
                         block B features [RNA]
            - single_C:  CV-Results for each fold on the testdata w/ only 
                         block C features [Mutation]
            - single_D:  CV-Results for each fold on the testdata w/ only 
                         block D features [Mirna]
        * 'settings' [settings used to do the CV - all arguments!]
            - datapath, seed, response, mtry,....   "
  # [0] Check Inputs -----------------------------------------------------------
  # 0-0 data_path, response are all checked within 'load_data_extract_block_names()'
  # 0-1 mtry, min_node_size & num_trees are all checked within 'rfsrc()'
  
  # 0-2 weighted must be boolean & seed an integer
  assert_logical(weighted)
  assert_int(seed)
  
  # 0-3 Check weight_metric to be meaningful [if weighted is true at all]!
  if (weighted) {
    # check that it is character
    assert_character(weight_metric)
    
    # check that it is of length 1
    if (length(weight_metric) != 1) {
      stop("'weight_metric' has more than 1 element!")
    }
    
    # check it has a valid value!
    if (!(weight_metric %in% c("Accuracy", "F1"))) {
      stop("'weight_metric' must be 'Accuracy' or 'F1'!")
    }
  }
  
  # [1] Get the data & dimensions of train & test folds! -----------------------
  # 1-1 Load Data & the names of the single blocks in the data!
  data <- load_data_extract_block_names(path = data_path, response = response)
  
  # 1-2 Get Obs. per fold [Train & Test]
  obs_per_fold <- get_obs_per_fold(data = data$data)
  
  # [2] Shuffle the data & create lists to save the results in -----------------
  # 2-1 Shuffle IDs of data, we use to split data later on!
  set.seed(seed)
  fold_ids <- sample(nrow(data$data), nrow(data$data), replace = FALSE)
  
  # 2-2 Create empty lists to store results in!
  # 2-2-1 Full TestSet
  full <- list()
  
  # 2-2-2 TestSet with 1 missing omics-block
  miss1_A <- list(); miss1_B <- list(); miss1_C <- list(); miss1_D <- list()
  
  # 2-2-3 TestSet with 2 missing omics-blocks
  miss2_CD <- list(); miss2_BD <- list(); miss2_BC <- list(); miss2_AD <- list()
  miss2_AC <- list(); miss2_AB <- list()
  
  # 2-2-4 TestSet with 3 missing omics-blocks
  miss3_ABC <- list(); miss3_ABD <- list(); miss3_ACD <- list(); miss3_BCD <- list()
  
  # 2-2-5 Single BlockTestSet
  single_A <- list(); single_B <- list(); single_CL <- list() 
  single_C <- list(); single_D <- list()
  
  # [3] Start the CV [5-fold per default!] -------------------------------------
  for (i in 0:4) {
    
    print(paste0("FOLD: ", as.character(i + 1), "/5 -------------------------"))
    
    # [1] Get TestSet from 'data', by taking the first 'obs_per_fold$amount_test' 
    #     IDs in 'fold_ids'
    test_ids <- fold_ids[((i * obs_per_fold$amount_test) + 1):(((i + 1) * obs_per_fold$amount_test))]
    test_df  <- data$data[test_ids,]
    
    # [2] Get the TrainSet from 'data' [= IDs not in TestSet] 
    train_ids <- fold_ids[-which(fold_ids %in% test_ids)]
    train_df  <- data$data[train_ids,]
    
    # [3] Induce blockwise missingness [SCENARIO_1]
    # 3-1 Sample equally sized 'observed' blocks [according to SCENARIO_1]
    set.seed(seed + i)
    observed_blocks <- sample(c(rep("Clin, A", obs_per_fold$amount_train_fold), 
                                rep("Clin, B", obs_per_fold$amount_train_fold),
                                rep("Clin, C", obs_per_fold$amount_train_fold), 
                                rep("Clin, D", obs_per_fold$amount_train_fold)),
                              obs_per_fold$amount_train, replace = FALSE)
    
    # 3-2 Subset the observations in the single Blocks that have been observed
    #     Select the Observations that have been observed [s. 'observed_blocks']
    #     and only keep the variables in the certain block!
    clin_block      <- train_df[grep("Clin", observed_blocks), 
                                c(response, data$block_names$clin_block)]
    cnv_block       <- train_df[grep("A", observed_blocks), 
                                c(response, data$block_names$cnv_block)]
    rna_block       <- train_df[grep("B", observed_blocks), 
                                c(response, data$block_names$rna_block)]
    mutation_block  <- train_df[grep("C", observed_blocks), 
                                c(response, data$block_names$mutation_block)]
    mirna_block     <- train_df[grep("D", observed_blocks), 
                                c(response, data$block_names$mirna_block)]
    
    # [4] Fit RFs on the seperate blocks
    #     'Formula' for the fitting has form of 'response ~ .' 
    #      --> use all cols as features except for the 'response'
    # 4-1 Clincal Block
    RF_clin <- rfsrc(formula = as.formula(paste(eval(response), "~ .")),
                     data = clin_block, ntree = num_trees, mtry = mtry, 
                     nodesize = min_node_size, samptype = "swr",
                     seed = seed)
    
    # 4-1-2 On CNV Block
    RF_cnv <- rfsrc(formula = as.formula(paste(eval(response), "~ .")),
                    data = cnv_block, ntree = num_trees, mtry = mtry, 
                    nodesize = min_node_size, samptype = "swr",
                    seed = seed)
    
    # 4-1-3 On RNA Block only
    RF_rna <- rfsrc(formula = as.formula(paste(eval(response), "~ .")),
                    data = rna_block, ntree = num_trees, mtry = mtry, 
                    nodesize = min_node_size, samptype = "swr",
                    seed = seed)
    
    # 4-1-4 On Mutation Block only
    RF_mutation <- rfsrc(formula = as.formula(paste(eval(response), "~ .")),
                         data = mutation_block, ntree = num_trees, mtry = mtry, 
                         nodesize = min_node_size, samptype = "swr",
                         seed = seed)
    
    # 4-1-5 On Mirna Block only
    RF_mirna <- rfsrc(formula = as.formula(paste(eval(response), "~ .")),
                      data = mirna_block, ntree = num_trees, mtry = mtry, 
                      nodesize = min_node_size, samptype = "swr",
                      seed = seed)
    
    # 4-2 Collect all single RFs to the Forest!
    Forest <- list(RF_clin, RF_cnv, RF_rna, RF_mutation, RF_mirna)
    
    # [6] Get predicitons on the different testsets!
    # 6-1 TestSet, when there is only one observed block! 
    #     Also done to extract the predictive performance of the different RFs,
    #     that were fitted to the single blocks! 
    #     As testset has 1 block only, we don't need to weight the different 
    #     blockwise RFs, as only one single block is used to obtain predicitons!
    #     --> the resulting performance will be used as weight for the remaining
    #         testsettings!
    print("Evaluation TestSet w/ 1 observed block only------------------------")
    print("Obtain performance & corresponding weights on the single blocks!")
    
    single_A[[i + 1]]  <- do_evaluation_rfsrc(Forest = Forest, weights = c(1, 1, 1, 1, 1),
                                              testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$rna_block,
                                                                                                  data$block_name$mutation_block,
                                                                                                  data$block_name$mirna_block,
                                                                                                  data$block_names$clin_block))])
    
    single_B[[i + 1]]  <- do_evaluation_rfsrc(Forest = Forest, weights = c(1, 1, 1, 1, 1), 
                                              testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$cnv_block,
                                                                                                  data$block_name$mutation_block,
                                                                                                  data$block_name$mirna_block,
                                                                                                  data$block_names$clin_block))])
    
    single_C[[i + 1]]  <- do_evaluation_rfsrc(Forest = Forest, weights = c(1, 1, 1, 1, 1),
                                              testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$rna_block,
                                                                                                  data$block_name$cnv_block,
                                                                                                  data$block_name$mirna_block,
                                                                                                  data$block_names$clin_block))])
    
    single_D[[i + 1]]  <- do_evaluation_rfsrc(Forest = Forest, weights = c(1, 1, 1, 1, 1), 
                                              testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$rna_block,
                                                                                                  data$block_name$mutation_block,
                                                                                                  data$block_name$cnv_block,
                                                                                                  data$block_names$clin_block))])
    single_CL[[i + 1]]  <- do_evaluation_rfsrc(Forest = Forest, weights = c(1, 1, 1, 1, 1),
                                               testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$rna_block,
                                                                                                   data$block_name$mutation_block,
                                                                                                   data$block_name$mirna_block,
                                                                                                   data$block_names$cnv_block))])
    
    # 6-1-2 Extract the weights from the single block performance - weighted is true!
    weights <- c()
    if (weighted) {
      weights <- c(weights, single_CL[[i + 1]][[weight_metric]])
      weights <- c(weights, single_A[[i + 1]][[weight_metric]])
      weights <- c(weights, single_B[[i + 1]][[weight_metric]])
      weights <- c(weights, single_C[[i + 1]][[weight_metric]])
      weights <- c(weights, single_D[[i + 1]][[weight_metric]])
    } else {
      weights <- rep(1, times = 5)
    }
    
    # 6-1-3 Norm the weights
    weights <- weights / sum(weights)
    
    print("Extracted weights for different blocks [CL, A, B, C, D] & normalized them")
    print(weights)
    
    # 6-2 Full TestSet
    print("Evaluation full TestSet -------------------------------------------")
    full[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, testdata = test_df, 
                                         weights = weights) # FULL TESTSET
    
    # 6-3 TestSet, where one of the omics blocks is missing!
    print("Evaluation TestSet w/ 1 missing omics block------------------------")
    miss1_A[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                            testdata = test_df[,-which(colnames(test_df) %in% data$block_name$cnv_block)])
    
    miss1_B[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights, 
                                            testdata = test_df[,-which(colnames(test_df) %in% data$block_name$rna_block)])
    
    miss1_C[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                            testdata = test_df[,-which(colnames(test_df) %in% data$block_name$mutation_block)])
    
    miss1_D[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                            testdata = test_df[,-which(colnames(test_df) %in% data$block_name$mirna_block)])
    
    # 6-4 TestSet, where two of the omics blocks are missing!
    print("Evaluation TestSet w/ 2 missing omics blocks-----------------------")
    miss2_CD[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                             testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$mutation_block,
                                                                                                 data$block_name$mirna_block))])
    
    miss2_BD[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                             testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$mirna_block,
                                                                                                 data$block_names$rna_block))])
    
    miss2_BC[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                             testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$mutation_block,
                                                                                                 data$block_name$rna_block))])
    
    miss2_AD[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                             testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$cnv_block,
                                                                                                 data$block_name$mirna_block))])
    
    miss2_AC[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                             testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$cnv_block,
                                                                                                 data$block_name$mutation_block))])
    
    miss2_AB[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                             testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$cnv_block,
                                                                                                 data$block_name$rna_block))])
    # 6-5 TestSet, where three of the omics blocks are missing!
    print("Evaluation TestSet w/ 3 missing omics blocks-----------------------")
    miss3_ABC[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                              testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$cnv_block,
                                                                                                  data$block_name$rna_block,
                                                                                                  data$block_name$mutation_block))])
    miss3_ACD[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights, 
                                              testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$cnv_block,
                                                                                                  data$block_name$mutation_block,
                                                                                                  data$block_name$mirna_block))])
    
    miss3_ABD[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights, 
                                              testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$cnv_block,
                                                                                                  data$block_name$rna_block,
                                                                                                  data$block_name$mirna_block))])
    
    miss3_BCD[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights, 
                                        testdata = test_df[,-which(colnames(test_df) %in% c(data$block_name$rna_block,
                                                                                            data$block_name$mutation_block,
                                                                                            data$block_name$mirna_block))])
  }
  
  # [4] Return the results & settings of parameters used to do CV! -------------
  res_all <- list("full" = full,
                  "miss1_A" = miss1_A, "miss1_B" = miss1_B,
                  "miss1_C" = miss1_C, "miss1_D" = miss1_D,
                  "miss2_CD" = miss2_CD, "miss2_BD" = miss2_BD,
                  "miss2_BC" = miss2_BC, "miss2_AD" = miss2_AD,
                  "miss2_AC" = miss2_AC, "miss2_AB" = miss2_AB,
                  "miss3_ABC" = miss3_ABC, "miss3_ABD" = miss3_ABD,
                  "miss3_ACD" = miss3_ACD, "miss3_BCD" = miss3_BCD,
                  "single_A" = single_A, "single_B" = single_B,
                  "single_C" = single_C, "single_D" = single_D,
                  "single_CL" = single_CL)
  
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
do_CV_NK_setting2             <- function(data_path = "data/external/Dr_Hornung/Data/ProcessedData_subsets/seed_1234/KIRC_Subset.RData",
                                          response = "gender", seed = 1312,
                                          weighted = TRUE, weight_metric = NULL,
                                          num_trees = as.integer(10), mtry = NULL, 
                                          min_node_size = 10) {
  
  " Function to evaluate RF Adaption on blockwise missing data from Norbert
    Krautenbacher [proposed in it's Ph.D.] Approach!
    
    Data is split into test and train set [5-fold], w/ little adjustment, that
    all train folds have same amount of obs., so the TestFold can be a bit
    smaller than the testfold!
    Then each [equally sized] trainingsfold is censored to scenario 2, so that 
    each fold has an observed clinical block + an additional observed omics 
    block [1. fold has Clin + 4 omics blocks (fully observed),
           2. fold has Clin + 3 omics blocks, 
           3. fold has Clin + 2 omics blocks,
           4. fold has Clin + 1 Omics Block!]
    
    Then we train a serperate RandomForest on  each of the feature blocks 
    [clin, omicsblock1,...] and ensemble the predicitons from these RFs to a 
    single prediciton and rate these w/ Accuracy, Precision, F1-Socre,....
    The TestingSituations are different, as we can test the models on fully 
    observed testdata, on testdata w/ 1 missing block etc. etc.
    
    Args:
      - data_path (char)    : Path to the data, we want to CV! This should lead 
                              to a file w/ multiple sub DFs 
                              [details see 'create_data()']
      - response (char)     : The repsonse we want to model - 
                              MUST be in the 'clin'-block!
      - seed (int)          : Needed for assiging which blocks are which letters,
                              shuffeling of observations and assigning them in a 
                              reproducable way!
      - weighted (bool)     : Shall the predicitons from the different blocks be
                              weighted by the oob accuracy of the blocks?!
      - weight_metric (chr) : Which Metric shall be used to weight the different 
                              RFs - the better the single blocks according to the
                              metric, the higher their weight when ensembling!
                              Must be in c('Accuracy', 'F1')
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
            - miss3_ABC: CV Results for each fold on the testdata, w/ 3 missing 
                         omics-block [cnv & rna & mutation]!
            - miss3_ACD: CV Results for each fold on the testdata, w/ 3 missing 
                         omics-block [cnv & mutation & mirna]!
            - miss3_ABD: CV Results for each fold on the testdata, w/ 3 missing 
                         omics-block [cnv & rna & mirna]!
            - miss3_BCD: CV Results for each fold on the testdata w/ 3 missing 
                         omics -bloc [rna & mutation & mirna]
            - single_CL: CV-Results for each fold on the testdata w/ only 
                         clinical features
            - single_A:  CV-Results for each fold on the testdata w/ only 
                         block A features [CNV]
            - single_B:  CV-Results for each fold on the testdata w/ only 
                         block B features [RNA]
            - single_C:  CV-Results for each fold on the testdata w/ only 
                         block C features [Mutation]
            - single_D:  CV-Results for each fold on the testdata w/ only 
                         block D features [Mirna]
        * 'settings' [settings used to do the CV - all arguments!]
            - datapath, seed, response, mtry,....   "
  # [0] Check Inputs -----------------------------------------------------------
  # 0-0 data_path, response are all checked within 'load_data_extract_block_names()'
  # 0-1 mtry, min_node_size & num_trees are all checked within 'rfsrc()'
  
  # 0-2 weighted must be boolean & seed an integer
  assert_logical(weighted)
  assert_int(seed)
  
  # 0-3 Check weight_metric to be meaningful [if weighted is true at all]!
  if (weighted) {
    # check that it is character
    assert_character(weight_metric)
    
    # check that it is of length 1
    if (length(weight_metric) != 1) {
      stop("'weight_metric' has more than 1 element!")
    }
    
    # check it has a valid value!
    if (!(weight_metric %in% c("Accuracy", "F1"))) {
      stop("'weight_metric' must be 'Accuracy' or 'F1'!")
    }
  }
  
  # [1] Get the data & dimensions of train & test folds! -----------------------
  # 1-1 Load Data & the names of the single blocks in the data!
  data <- load_data_extract_block_names(path = data_path, response = response)
  
  # 1-2 Get Obs. per fold [Train & Test]
  obs_per_fold <- get_obs_per_fold(data = data$data)
  
  # [2] Shuffle the data & create lists to save the results in -----------------
  # 2-1 Shuffle IDs of data, we use to split data later on!
  set.seed(seed)
  fold_ids <- sample(nrow(data$data), nrow(data$data), replace = FALSE)
  
  # 2-2 Create empty lists to store results in!
  # 2-2-1 Full TestSet
  full <- list()
  
  # 2-2-2 TestSet with 1 missing omics-block
  miss1_A <- list(); miss1_B <- list(); miss1_C <- list(); miss1_D <- list()
  
  # 2-2-3 TestSet with 2 missing omics-blocks
  miss2_CD <- list(); miss2_BD <- list(); miss2_BC <- list(); miss2_AD <- list()
  miss2_AC <- list(); miss2_AB <- list()
  
  # 2-2-4 TestSet with 3 missing omics-blocks
  miss3_ABC <- list(); miss3_ABD <- list(); miss3_ACD <- list(); miss3_BCD <- list()
  
  # 2-2-5 Single BlockTestSet
  single_A <- list(); single_B <- list(); single_CL <- list() 
  single_C <- list(); single_D <- list()
  
  # 2-3 Randomly assign the omics blocks to the letters 'A', 'B', 'C', 'D', 
  #     as SCEANRIO2, highly depens on which block is where!
  #     ['A' only observed in 1.fold, whereas 'D' is observed in all folds]
  set.seed(seed)
  letter_feas <- sample(c("data$block_names$cnv_block", "data$block_names$rna_block", 
                          "data$block_names$mutation_block", "data$block_names$mirna_block"), 
                        4, replace = FALSE)
  names(letter_feas) <- c("A", "B", "C", "D")
  
  # [3] Start the CV [5-fold per default!] -------------------------------------
  for (i in 0:4) {
    
    print(paste0("FOLD: ", as.character(i + 1), "/5 -------------------------"))
    
    # [1] Get TestSet from 'data', by taking the first 'obs_per_fold$amount_test' 
    #     IDs in 'fold_ids'
    test_ids <- fold_ids[((i * obs_per_fold$amount_test) + 1):(((i + 1) * obs_per_fold$amount_test))]
    test_df  <- data$data[test_ids,]
    
    # [2] Get the TrainSet from 'data' [= IDs not in TestSet] 
    train_ids <- fold_ids[-which(fold_ids %in% test_ids)]
    train_df  <- data$data[train_ids,]
    
    # [3] Induce blockwise missingness [SCENARIO_1]
    # 3-1 Sample equally sized 'observed' blocks [according to SCENARIO_1]
    set.seed(seed + i)
    observed_blocks <- sample(c(rep("Clin, A, B, C, D", obs_per_fold$amount_train_fold), 
                                rep("Clin, B, C, D", obs_per_fold$amount_train_fold),
                                rep("Clin, C, D", obs_per_fold$amount_train_fold), 
                                rep("Clin, D", obs_per_fold$amount_train_fold)),
                              obs_per_fold$amount_train, replace = FALSE)
    
    # 3-2 Subset the observations in the single Blocks that have been observed
    #     Select the Observations that have been observed [s. 'observed_blocks']
    #     and only keep the variables in the certain block!
    clin_block  <- train_df[grep("Clin", observed_blocks), 
                            c(response, data$block_names$clin_block)]
    block_A     <- train_df[grep("A,", observed_blocks), 
                            c(response, eval(parse(text = letter_feas["A"])))]
    block_B     <- train_df[grep("B,", observed_blocks), 
                            c(response, eval(parse(text = letter_feas["B"])))]
    block_C     <- train_df[grep("C,", observed_blocks), 
                            c(response, eval(parse(text = letter_feas["C"])))]
    block_D     <- train_df[grep("D", observed_blocks), 
                            c(response, eval(parse(text = letter_feas["D"])))]
    
    # [4] Fit RFs on the seperate blocks
    #     'Formula' for the fitting has form of 'response ~ .' 
    #      --> use all cols as features except for the 'response'
    # 4-1 Clincal Block
    RF_clin <- rfsrc(formula = as.formula(paste(eval(response), "~ .")),
                     data = clin_block, ntree = num_trees, mtry = mtry, 
                     nodesize = min_node_size, samptype = "swr",
                     seed = seed)
    
    # 4-1-2 Block A
    RF_A <- rfsrc(formula = as.formula(paste(eval(response), "~ .")),
                  data = block_A, ntree = num_trees, mtry = mtry, 
                  nodesize = min_node_size, samptype = "swr",
                  seed = seed)
    
    # 4-1-3 Block B
    RF_B <- rfsrc(formula = as.formula(paste(eval(response), "~ .")),
                  data = block_B, ntree = num_trees, mtry = mtry, 
                  nodesize = min_node_size, samptype = "swr",
                  seed = seed)
    
    # 4-1-4 Block C
    RF_C <- rfsrc(formula = as.formula(paste(eval(response), "~ .")),
                  data = block_C, ntree = num_trees, mtry = mtry, 
                  nodesize = min_node_size, samptype = "swr",
                  seed = seed)
    
    # 4-1-5 Block D
    RF_D <- rfsrc(formula = as.formula(paste(eval(response), "~ .")),
                  data = block_D, ntree = num_trees, mtry = mtry, 
                  nodesize = min_node_size, samptype = "swr",
                  seed = seed)
    
    # 4-2 Collect all single RFs to the Forest!
    Forest <- list(RF_clin, RF_A, RF_B, RF_C, RF_D)
    
    # [6] Get (unweighted) predicitons on the different testsets!
    # 6-1 TestSet, when there is only one observed block! 
    #     Also done to extract the predictive performance of the different RFs,
    #     that were fitted to the single blocks! 
    #     As testset has 1 block only, we don't need to weight the different 
    #     blockwise RFs, as only one single block is used to obtain predicitons!
    #     --> the resulting performance will be used as weight for the remaining
    #         testsettings!
    print("Evaluation TestSet w/ 1 observed block only------------------------")
    print("Obtain performance & corresponding weights on the single blocks!")
    
    single_A[[i + 1]]  <- do_evaluation_rfsrc(Forest = Forest, weights = c(1, 1, 1, 1, 1),  
                                              testdata = test_df[,-which(colnames(test_df) %in% c(data$block_names$clin_block,
                                                                                                  eval(parse(text = letter_feas["B"])),
                                                                                                  eval(parse(text = letter_feas["C"])),
                                                                                                  eval(parse(text = letter_feas["D"]))))])
    single_B[[i + 1]]  <- do_evaluation_rfsrc(Forest = Forest, weights = c(1, 1, 1, 1, 1),  
                                              testdata = test_df[,-which(colnames(test_df) %in% c(data$block_names$clin_block,
                                                                                                  eval(parse(text = letter_feas["A"])),
                                                                                                  eval(parse(text = letter_feas["C"])),
                                                                                                  eval(parse(text = letter_feas["D"]))))])
    
    single_C[[i + 1]]  <- do_evaluation_rfsrc(Forest = Forest, weights = c(1, 1, 1, 1, 1),  
                                              testdata = test_df[,-which(colnames(test_df) %in% c(data$block_names$clin_block,
                                                                                                  eval(parse(text = letter_feas["B"])),
                                                                                                  eval(parse(text = letter_feas["A"])),
                                                                                                  eval(parse(text = letter_feas["D"]))))])
    single_D[[i + 1]]  <- do_evaluation_rfsrc(Forest = Forest, weights = c(1, 1, 1, 1, 1),  
                                              testdata = test_df[,-which(colnames(test_df) %in% c(data$block_names$clin_block,
                                                                                                  eval(parse(text = letter_feas["B"])),
                                                                                                  eval(parse(text = letter_feas["C"])),
                                                                                                  eval(parse(text = letter_feas["A"]))))])
    single_CL[[i + 1]]  <- do_evaluation_rfsrc(Forest = Forest, weights = c(1, 1, 1, 1, 1),  
                                               testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["A"])),
                                                                                                   eval(parse(text = letter_feas["B"])),
                                                                                                   eval(parse(text = letter_feas["C"])),
                                                                                                   eval(parse(text = letter_feas["D"]))))])
    
    # 6-1-2 Extract the weights from the single block performance - weighted is true!
    weights <- c()
    if (weighted) {
      weights <- c(weights, single_CL[[i + 1]][[weight_metric]])
      weights <- c(weights, single_A[[i + 1]][[weight_metric]])
      weights <- c(weights, single_B[[i + 1]][[weight_metric]])
      weights <- c(weights, single_C[[i + 1]][[weight_metric]])
      weights <- c(weights, single_D[[i + 1]][[weight_metric]])
    } else {
      weights <- rep(1, times = 5)
    }
    
    # 6-1-3 Norm the weights
    weights <- weights / sum(weights)
    
    print("Extracted weights for different blocks [CL, A, B, C, D] & normalized them")
    print(weights)
    
    # 6-2 Full TestSet
    print("Evaluation full TestSet -------------------------------------------")
    full[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, testdata = test_df, 
                                         weights = weights) # FULL TESTSET
    
    # 6-3 TestSet, where one of the omics blocks is missing!
    print("Evaluation TestSet w/ 1 missing omics block------------------------")
    miss1_A[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                            testdata = test_df[,-which(colnames(test_df) %in% eval(parse(text = letter_feas["A"])))])
    
    miss1_B[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights, 
                                            testdata = test_df[,-which(colnames(test_df) %in% eval(parse(text = letter_feas["B"])))])
    
    miss1_C[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                            testdata = test_df[,-which(colnames(test_df) %in% eval(parse(text = letter_feas["C"])))])
    
    miss1_D[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                            testdata = test_df[,-which(colnames(test_df) %in% eval(parse(text = letter_feas["D"])))])
    
    # 6-4 TestSet, where two of the omics blocks are missing!
    print("Evaluation TestSet w/ 2 missing omics blocks-----------------------")
    miss2_CD[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                             testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["C"])),
                                                                                                 eval(parse(text = letter_feas["D"]))))])
    
    miss2_BD[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                             testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["B"])),
                                                                                                 eval(parse(text = letter_feas["D"]))))])
    
    miss2_BC[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                             testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["C"])),
                                                                                                 eval(parse(text = letter_feas["B"]))))])
    
    miss2_AD[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                             testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["A"])),
                                                                                                 eval(parse(text = letter_feas["D"]))))])
    
    miss2_AC[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                             testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["C"])),
                                                                                                 eval(parse(text = letter_feas["A"]))))])
    
    miss2_AB[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                             testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["A"])),
                                                                                                 eval(parse(text = letter_feas["B"]))))])
    # 6-5 TestSet, where three of the omics blocks are missing!
    print("Evaluation TestSet w/ 3 missing omics blocks-----------------------")
    miss3_ABC[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                              testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["A"])),
                                                                                                  eval(parse(text = letter_feas["B"])),
                                                                                                  eval(parse(text = letter_feas["C"]))))])
    miss3_ACD[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights, 
                                              testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["A"])),
                                                                                                  eval(parse(text = letter_feas["D"])),
                                                                                                  eval(parse(text = letter_feas["C"]))))])
    
    miss3_ABD[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights, 
                                              testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["A"])),
                                                                                                  eval(parse(text = letter_feas["B"])),
                                                                                                  eval(parse(text = letter_feas["D"]))))])
    
    miss3_BCD[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights, 
                                              testdata = test_df[,-which(colnames(test_df) %in% c(eval(parse(text = letter_feas["D"])),
                                                                                                  eval(parse(text = letter_feas["B"])),
                                                                                                  eval(parse(text = letter_feas["C"]))))])
  }
  
  # [4] Return the results & settings of parameters used to do CV! -------------
  res_all <- list("full" = full,
                  "miss1_A" = miss1_A, "miss1_B" = miss1_B,
                  "miss1_C" = miss1_C, "miss1_D" = miss1_D,
                  "miss2_CD" = miss2_CD, "miss2_BD" = miss2_BD,
                  "miss2_BC" = miss2_BC, "miss2_AD" = miss2_AD,
                  "miss2_AC" = miss2_AC, "miss2_AB" = miss2_AB,
                  "miss3_ABC" = miss3_ABC, "miss3_ABD" = miss3_ABD,
                  "miss3_ACD" = miss3_ACD, "miss3_BCD" = miss3_BCD,
                  "single_A" = single_A, "single_B" = single_B,
                  "single_C" = single_C, "single_D" = single_D,
                  "single_CL" = single_CL)
  
  # Collect the Settings, used to do the CV!
  settings <- list("data_path"     = data_path,
                   "response"      = response, 
                   "seed"          = seed, 
                   "num_trees"     = num_trees,
                   "mtry"          = mtry, 
                   "min_node_size" = min_node_size,
                   "block_letters" = letter_feas)
  
  # Return both lists!
  return(list("res_all" = res_all, 
              "settings" = settings))
}
do_CV_NK_setting3             <- function(data_path = "data/external/Dr_Hornung/Data/ProcessedData_subsets/seed_1234/KIRC_Subset.RData",
                                          response = "gender", seed = 1312, 
                                          weighted = TRUE, weight_metric = NULL,
                                          num_trees = as.integer(10), mtry = NULL, 
                                          min_node_size = 10) {
  
  " Function to evaluate RF Adaption on blockwise missing data from Norbert
    Krautenbacher [proposed in it's Ph.D.] Approach!
    
    Data is split into test and train set [5-fold], w/ little adjustment, that
    all train folds have same amount of obs., so the TestFold can be a bit
    smaller than the testfold!
   
    Then each [equally sized] trainingsfold is censored to scenario 3: That is:
    For each fold [1-4] we randomly sample which blocks are observed for the
    different blocks --> each block is deleted w/ prob. of 1/3 & kept w/ pro. 2/3
    
    Then we train a serperate RandomForest on  each of the feature blocks 
    [clin, omicsblock1,...] and ensemble the predicitons from these RFs to a 
    single prediciton and rate these w/ Accuracy, Precision, F1-Socre,....
    The TestingSituations are different, as we can test the models on fully 
    observed testdata, on testdata w/ 1 missing block etc. etc.
    
    Args:
      - data_path (char)    : Path to the data, we want to CV! This should lead 
                              to a file w/ multiple sub DFs 
                              [details see 'create_data()']
      - response (char)     : The repsonse we want to model - 
                              MUST be in the 'clin'-block!
      - seed (int)          : Needed to assign Obs. to folds, plus which folds
                              have which observed feas in a reproducable way!
      - weighted (bool)     : Shall the predicitons from the different blocks be
                              weighted by the oob accuracy of the blocks?!
      - weight_metric (chr) : Which Metric shall be used to weight the different 
                              RFs - the better the single blocks according to the
                              metric, the higher their weight when ensembling!
                              Must be in c('Accuracy', 'F1')
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
            - miss3_ABC: CV Results for each fold on the testdata, w/ 3 missing 
                         omics-block [cnv & rna & mutation]!
            - miss3_ACD: CV Results for each fold on the testdata, w/ 3 missing 
                         omics-block [cnv & mutation & mirna]!
            - miss3_ABD: CV Results for each fold on the testdata, w/ 3 missing 
                         omics-block [cnv & rna & mirna]!
            - miss3_BCD: CV Results for each fold on the testdata w/ 3 missing 
                         omics -bloc [rna & mutation & mirna]
            - single_CL: CV-Results for each fold on the testdata w/ only 
                         clinical features
            - single_A:  CV-Results for each fold on the testdata w/ only 
                         block A features [CNV]
            - single_B:  CV-Results for each fold on the testdata w/ only 
                         block B features [RNA]
            - single_C:  CV-Results for each fold on the testdata w/ only 
                         block C features [Mutation]
            - single_D:  CV-Results for each fold on the testdata w/ only 
                         block D features [Mirna]
        * 'settings' [settings used to do the CV - all arguments!]
            - datapath, seed, response, mtry,....   "
  # [0] Check Inputs -----------------------------------------------------------
  # 0-0 data_path, response are all checked within 'load_data_extract_block_names()'
  # 0-1 mtry, min_node_size & num_trees are all checked within 'rfsrc()'
  
  # 0-2 weighted must be boolean & seed an integer
  assert_logical(weighted)
  assert_int(seed)
  
  # 0-3 Check weight_metric to be meaningful [if weighted is true at all]!
  if (weighted) {
    # check that it is character
    assert_character(weight_metric)
    
    # check that it is of length 1
    if (length(weight_metric) != 1) {
      stop("'weight_metric' has more than 1 element!")
    }
    
    # check it has a valid value!
    if (!(weight_metric %in% c("Accuracy", "F1"))) {
      stop("'weight_metric' must be 'Accuracy' or 'F1'!")
    }
  }
  
  # [1] Get the data & dimensions of train & test folds! -----------------------
  # 1-1 Load Data & the names of the single blocks in the data!
  data <- load_data_extract_block_names(path = data_path, response = response)
  
  # 1-2 Get Obs. per fold [Train & Test]
  obs_per_fold <- get_obs_per_fold(data = data$data)
  
  # [2] Shuffle the data & create lists to save the results in -----------------
  # 2-1 Shuffle IDs of data, we use to split data later on!
  set.seed(seed)
  fold_ids <- sample(nrow(data$data), nrow(data$data), replace = FALSE)
  
  # 2-2 Create empty lists to store results in!
  # 2-2-1 Full TestSet
  full <- list()
  
  # 2-2-2 TestSet with 1 missing omics-block
  miss1_A <- list(); miss1_B <- list(); miss1_C <- list(); miss1_D <- list()
  
  # 2-2-3 TestSet with 2 missing omics-blocks
  miss2_CD <- list(); miss2_BD <- list(); miss2_BC <- list(); miss2_AD <- list()
  miss2_AC <- list(); miss2_AB <- list()
  
  # 2-2-4 TestSet with 3 missing omics-blocks
  miss3_ABC <- list(); miss3_ABD <- list(); miss3_ACD <- list(); miss3_BCD <- list()
  
  # 2-2-5 Single BlockTestSet
  single_A <- list(); single_B <- list(); single_CL <- list() 
  single_C <- list(); single_D <- list()
  
  # [3] Assign which blocks were observed for which folds ----------------------
  #     Randomly assign the 'observed' blocks to the 4 different folds!
  #     The Index of TRUE / FALSE are indicators, whether the blocks are observed
  #     [1] = clinical; [2] = CNV; [3] = RNA; [4] = Mutation; [5] = Mirna
  set.seed(seed)
  fold1_obs <- sample(c(TRUE, FALSE), size = 5, replace = T, prob = c(2/3, 1/3))
  names(fold1_obs) <- c("Clin", "A", "B", "C", "D")
  set.seed(seed + 1)
  fold2_obs <- sample(c(TRUE, FALSE), size = 5, replace = T, prob = c(2/3, 1/3))
  names(fold2_obs) <- c("Clin", "A", "B", "C", "D")
  set.seed(seed + 2)
  fold3_obs <- sample(c(TRUE, FALSE), size = 5, replace = T, prob = c(2/3, 1/3))
  names(fold3_obs) <- c("Clin", "A", "B", "C", "D")
  set.seed(seed + 3)
  fold4_obs <- sample(c(TRUE, FALSE), size = 5, replace = T, prob = c(2/3, 1/3))
  names(fold4_obs) <- c("Clin", "A", "B", "C", "D")
  
  # 3-1 Save observed blocks per fold in a overall list!
  all_folds <- list("fold1" = fold1_obs, "fold2" = fold2_obs,
                    "fold3" = fold3_obs, "fold4" = fold4_obs)

  # [3] Start the CV [5-fold per default!] -------------------------------------
  for (i in 0:4) {
    
    print(paste0("FOLD: ", as.character(i + 1), "/5 -------------------------"))
    
    # [1] Get TestSet from 'data', by taking the first 'obs_per_fold$amount_test' 
    #     IDs in 'fold_ids'
    test_ids <- fold_ids[((i * obs_per_fold$amount_test) + 1):(((i + 1) * obs_per_fold$amount_test))]
    test_df  <- data$data[test_ids,]
    
    # [2] Get the TrainSet from 'data' [= IDs not in TestSet] 
    train_ids <- fold_ids[-which(fold_ids %in% test_ids)]
    train_df  <- data$data[train_ids,]
    
    # [3] Induce blockwise missingness [SCENARIO_3] 
    # 3-1 Sample equally sized 'observed' blocks [according to SCENARIO_3]
    set.seed(seed + i)
    observed_blocks <- sample(c(rep("fold1", obs_per_fold$amount_train_fold), 
                                rep("fold2", obs_per_fold$amount_train_fold),
                                rep("fold3", obs_per_fold$amount_train_fold), 
                                rep("fold4", obs_per_fold$amount_train_fold)),
                              obs_per_fold$amount_train, replace = FALSE)
    
    # 3-2 Now collect the different blocks from the different folds, where 
    #     the feas have been observed!
    # 3-2-1 Clin
    folds_w_clin <- which(sapply(all_folds, function(x) x["Clin"] == TRUE))
    folds_w_clin <- paste0("fold", folds_w_clin)
    rows_w_clin  <- which(observed_blocks %in% folds_w_clin)
    clin_block   <- train_df[rows_w_clin, c(response, data$block_names$clin_block)]
    
    # 3-2-2 CNV
    folds_w_cnv <- which(sapply(all_folds, function(x) x["A"] == TRUE))
    folds_w_cnv <- paste0("fold", folds_w_cnv)
    rows_w_cnv  <- which(observed_blocks %in% folds_w_cnv)
    cnv_block   <- train_df[rows_w_cnv, c(response, data$block_names$cnv_block)]
    
    # 3-2-3 RNA
    folds_w_rna <- which(sapply(all_folds, function(x) x["B"] == TRUE))
    folds_w_rna <- paste0("fold", folds_w_rna)
    rows_w_rna  <- which(observed_blocks %in% folds_w_rna)
    rna_block   <- train_df[rows_w_rna, c(response, data$block_names$rna_block)]
    
    # 3-2-4 MUTATION
    folds_w_mutation <- which(sapply(all_folds, function(x) x["C"] == TRUE))
    folds_w_mutation <- paste0("fold", folds_w_mutation)
    rows_w_mutation  <- which(observed_blocks %in% folds_w_mutation)
    mutation_block   <- train_df[rows_w_mutation, c(response, data$block_names$mutation_block)]
    
    # 3-2-5 MIRNA
    folds_w_mirna <- which(sapply(all_folds, function(x) x["D"] == TRUE))
    folds_w_mirna <- paste0("fold", folds_w_mirna)
    rows_w_mirna  <- which(observed_blocks %in% folds_w_mirna)
    mirna_block   <- train_df[rows_w_mirna, c(response, data$block_names$mirna_block)]
  
    
    # [4] Fit RFs on the seperate blocks
    #     'Formula' for the fitting has form of 'response ~ .' 
    #      --> use all cols as features except for the 'response'
    # 4-1 Clincal Block
    RF_clin <- rfsrc(formula = as.formula(paste(eval(response), "~ .")),
                     data = clin_block, ntree = num_trees, mtry = mtry, 
                     nodesize = min_node_size, samptype = "swr",
                     seed = seed)
    
    # 4-1-2 Block A
    RF_cnv <- rfsrc(formula = as.formula(paste(eval(response), "~ .")),
                    data = cnv_block, ntree = num_trees, mtry = mtry, 
                    nodesize = min_node_size, samptype = "swr",
                    seed = seed)
    
    # 4-1-3 Block B
    RF_rna <- rfsrc(formula = as.formula(paste(eval(response), "~ .")),
                    data = rna_block, ntree = num_trees, mtry = mtry, 
                    nodesize = min_node_size, samptype = "swr",
                    seed = seed)
    
    # 4-1-4 Block C
    RF_mutation <- rfsrc(formula = as.formula(paste(eval(response), "~ .")),
                         data = mutation_block, ntree = num_trees, mtry = mtry, 
                         nodesize = min_node_size, samptype = "swr",
                         seed = seed)
    
    # 4-1-5 Block D
    RF_mirna <- rfsrc(formula = as.formula(paste(eval(response), "~ .")),
                      data = mirna_block, ntree = num_trees, mtry = mtry, 
                      nodesize = min_node_size, samptype = "swr",
                      seed = seed)
    
    # 4-2 Collect all single RFs to the Forest!
    Forest <- list(RF_clin, RF_cnv, RF_rna, RF_mutation, RF_mirna)
    
    # [6] Get (unweighted) predicitons on the different testsets!
    # 6-1 TestSet, when there is only one observed block! 
    #     Also done to extract the predictive performance of the different RFs,
    #     that were fitted to the single blocks! 
    #     As testset has 1 block only, we don't need to weight the different 
    #     blockwise RFs, as only one single block is used to obtain predicitons!
    #     --> the resulting performance will be used as weight for the remaining
    #         testsettings!
    print("Evaluation TestSet w/ 1 observed block only------------------------")
    print("Obtain performance & corresponding weights on the single blocks!")
    
    single_A[[i + 1]]  <- do_evaluation_rfsrc(Forest = Forest, weights = c(1, 1, 1, 1, 1), 
                                              testdata = test_df[,-which(colnames(test_df) %in% c(data$block_names$mirna_block,
                                                                                                  data$block_names$rna_block,
                                                                                                  data$block_names$mutation_block,
                                                                                                  data$block_names$clin_block))])
    single_B[[i + 1]]  <- do_evaluation_rfsrc(Forest = Forest, weights = c(1, 1, 1, 1, 1), 
                                              testdata = test_df[,-which(colnames(test_df) %in% c(data$block_names$mirna_block,
                                                                                                  data$block_names$cnv_block,
                                                                                                  data$block_names$mutation_block,
                                                                                                  data$block_names$clin_block))])
    
    single_C[[i + 1]]  <- do_evaluation_rfsrc(Forest = Forest, weights = c(1, 1, 1, 1, 1), 
                                              testdata = test_df[,-which(colnames(test_df) %in% c(data$block_names$mirna_block,
                                                                                                  data$block_names$rna_block,
                                                                                                  data$block_names$cnv_block,
                                                                                                  data$block_names$clin_block))])
    single_D[[i + 1]]  <- do_evaluation_rfsrc(Forest = Forest, weights = c(1, 1, 1, 1, 1), 
                                              testdata = test_df[,-which(colnames(test_df) %in% c(data$block_names$cnv_block,
                                                                                                  data$block_names$rna_block,
                                                                                                  data$block_names$mutation_block,
                                                                                                  data$block_names$clin_block))])
    single_CL[[i + 1]]  <- do_evaluation_rfsrc(Forest = Forest, weights = c(1, 1, 1, 1, 1), 
                                               testdata = test_df[,-which(colnames(test_df) %in% c(data$block_names$mirna_block,
                                                                                                   data$block_names$rna_block,
                                                                                                   data$block_names$mutation_block,
                                                                                                   data$block_names$cnv_block))])
    
    # 6-1-2 Extract the weights from the single block performance - weighted is true!
    weights <- c()
    if (weighted) {
      weights <- c(weights, single_CL[[i + 1]][[weight_metric]])
      weights <- c(weights, single_A[[i + 1]][[weight_metric]])
      weights <- c(weights, single_B[[i + 1]][[weight_metric]])
      weights <- c(weights, single_C[[i + 1]][[weight_metric]])
      weights <- c(weights, single_D[[i + 1]][[weight_metric]])
    } else {
      weights <- rep(1, times = 5)
    }
    
    # 6-1-3 Norm the weights
    weights <- weights / sum(weights)
    
    print("Extracted weights for different blocks [CL, A, B, C, D] & normalized them")
    print(weights)
    
    # 6-2 Full TestSet
    print("Evaluation full TestSet -------------------------------------------")
    full[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, testdata = test_df, 
                                         weights = weights) # FULL TESTSET
    
    # 6-3 TestSet, where one of the omics blocks is missing!
    print("Evaluation TestSet w/ 1 missing omics block------------------------")
    miss1_A[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                            testdata = test_df[,-which(colnames(test_df) %in% data$block_names$cnv_block)])
    
    miss1_B[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights, 
                                            testdata = test_df[,-which(colnames(test_df) %in% data$block_names$rna_block)])
    
    miss1_C[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                            testdata = test_df[,-which(colnames(test_df) %in% data$block_names$mutation_block)])
    
    miss1_D[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                            testdata = test_df[,-which(colnames(test_df) %in% data$block_names$mirna_block)])
    
    # 6-4 TestSet, where two of the omics blocks are missing!
    print("Evaluation TestSet w/ 2 missing omics blocks-----------------------")
    miss2_CD[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest,weights = weights,
                                             testdata = test_df[,-which(colnames(test_df) %in% c(data$block_names$mutation_block,
                                                                                                 data$block_names$mirna_block))])
    
    miss2_BD[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                             testdata = test_df[,-which(colnames(test_df) %in% c(data$block_names$rna_block,
                                                                                                 data$block_names$mirna_block))])
    
    miss2_BC[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                             testdata = test_df[,-which(colnames(test_df) %in% c(data$block_names$mutation_block,
                                                                                                 data$block_names$rna_block))])
    
    miss2_AD[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                             testdata = test_df[,-which(colnames(test_df) %in% c(data$block_names$cnv_block,
                                                                                                 data$block_names$mirna_block))])
    
    miss2_AC[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                             testdata = test_df[,-which(colnames(test_df) %in% c(data$block_names$cnv_block,
                                                                                                 data$block_names$mutation_block))])
    
    miss2_AB[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                             testdata = test_df[,-which(colnames(test_df) %in% c(data$block_names$cnv_block,
                                                                                                 data$block_names$rna_block))])
    # 6-5 TestSet, where three of the omics blocks are missing!
    print("Evaluation TestSet w/ 3 missing omics blocks-----------------------")
    miss3_ABC[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                              testdata = test_df[,-which(colnames(test_df) %in% c(data$block_names$cnv_block,
                                                                                                  data$block_names$rna_block,
                                                                                                  data$block_names$mutation_block))])
    miss3_ACD[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights, 
                                              testdata = test_df[,-which(colnames(test_df) %in% c(data$block_names$cnv_block,
                                                                                                  data$block_names$mirna_block,
                                                                                                  data$block_names$mutation_block))])
    
    miss3_ABD[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest,weights = weights, 
                                              testdata = test_df[,-which(colnames(test_df) %in% c(data$block_names$cnv_block,
                                                                                                  data$block_names$rna_block,
                                                                                                  data$block_names$mirna_block))])
    
    miss3_BCD[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights, 
                                              testdata = test_df[,-which(colnames(test_df) %in% c(data$block_names$mirna_block,
                                                                                                  data$block_names$rna_block,
                                                                                                  data$block_names$mutation_block))])
    
    
  }
  
  # [4] Return the results & settings of parameters used to do CV! -------------
  res_all <- list("full" = full,
                  "miss1_A" = miss1_A, "miss1_B" = miss1_B,
                  "miss1_C" = miss1_C, "miss1_D" = miss1_D,
                  "miss2_CD" = miss2_CD, "miss2_BD" = miss2_BD,
                  "miss2_BC" = miss2_BC, "miss2_AD" = miss2_AD,
                  "miss2_AC" = miss2_AC, "miss2_AB" = miss2_AB,
                  "miss3_ABC" = miss3_ABC, "miss3_ABD" = miss3_ABD,
                  "miss3_ACD" = miss3_ACD, "miss3_BCD" = miss3_BCD,
                  "single_A" = single_A, "single_B" = single_B,
                  "single_C" = single_C, "single_D" = single_D,
                  "single_CL" = single_CL)
  
  # Collect the Settings, used to do the CV!
  settings <- list("data_path"     = data_path,
                   "response"      = response, 
                   "seed"          = seed, 
                   "num_trees"     = num_trees,
                   "mtry"          = mtry, 
                   "min_node_size" = min_node_size,
                   "observed_blocks_in_folds" = all_folds)
  
  # Return both lists!
  return(list("res_all" = res_all, 
              "settings" = settings))
}
do_CV_NK_setting4             <- function(data_path = "data/external/Dr_Hornung/Data/ProcessedData_subsets/seed_1234/KIRC_Subset.RData",
                                          response = "gender", seed = 1312, 
                                          weighted = TRUE, weight_metric = NULL,
                                          num_trees = as.integer(10), mtry = NULL, 
                                          min_node_size = 10) {
  
  " Function to evaluate RF Adaption on blockwise missing data from Norbert
    Krautenbacher [proposed in it's Ph.D.] Approach!
    
    Data is split into test and train set [5-fold], w/ little adjustment, that
    all train folds have same amount of obs., so the TestFold can be a bit
    smaller than the testfold!
   
    Then each [equally sized] trainingsfold is censored to scenario 4: That is:
      1. Bind 2 random omics blocks to a single feature block 
         [e.g. 'cnv' & 'mirna' is one block + 'rna' & 'mutation' is one block]
      2. Then divide the testset in two equally sized folds, from which we censor 
         one of the omics feature blocks in each of the folds!
    
    Then we train a serperate RandomForest on  each of the feature blocks 
    [clin, omicsblock1,...] and ensemble the predicitons from these RFs to a 
    single prediciton and rate these w/ Accuracy, Precision, F1-Socre,....
    The TestingSituations are different, as we can test the models on fully 
    observed testdata, on testdata w/ 1 missing block etc. etc.
    
    Args:
      - data_path (char)    : Path to the data, we want to CV! This should lead 
                              to a file w/ multiple sub DFs 
                              [details see 'create_data()']
      - response (char)     : The repsonse we want to model - 
                              MUST be in the 'clin'-block!
      - seed (int)          : Needed to assign 2 omics blocks to one, and assignig
                              the observations to the folds in a reproducable way!
      - weighted (bool)     : Shall the predicitons from the different blocks be
                              weighted by the oob accuracy of the blocks?!
      - weight_metric (chr) : Which Metric shall be used to weight the different 
                              RFs - the better the single blocks according to the
                              metric, the higher their weight when ensembling!
                              Must be in c('Accuracy', 'F1')
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
            - miss3_ABC: CV Results for each fold on the testdata, w/ 3 missing 
                         omics-block [cnv & rna & mutation]!
            - miss3_ACD: CV Results for each fold on the testdata, w/ 3 missing 
                         omics-block [cnv & mutation & mirna]!
            - miss3_ABD: CV Results for each fold on the testdata, w/ 3 missing 
                         omics-block [cnv & rna & mirna]!
            - miss3_BCD: CV Results for each fold on the testdata w/ 3 missing 
                         omics -bloc [rna & mutation & mirna]
            - single_CL: CV-Results for each fold on the testdata w/ only 
                         clinical features
            - single_A:  CV-Results for each fold on the testdata w/ only 
                         block A features [CNV]
            - single_B:  CV-Results for each fold on the testdata w/ only 
                         block B features [RNA]
            - single_C:  CV-Results for each fold on the testdata w/ only 
                         block C features [Mutation]
            - single_D:  CV-Results for each fold on the testdata w/ only 
                         block D features [Mirna]
        * 'settings' [settings used to do the CV - all arguments!]
            - datapath, seed, response, mtry,....   "
  # [0] Check Inputs -----------------------------------------------------------
  # 0-0 data_path, response are all checked within 'load_data_extract_block_names()'
  # 0-1 mtry, min_node_size & num_trees are all checked within 'rfsrc()'
  
  # 0-2 weighted must be boolean & seed an integer
  assert_logical(weighted)
  assert_int(seed)
  
  # 0-3 Check weight_metric to be meaningful [if weighted is true at all]!
  if (weighted) {
    # check that it is character
    assert_character(weight_metric)
    
    # check that it is of length 1
    if (length(weight_metric) != 1) {
      stop("'weight_metric' has more than 1 element!")
    }
    
    # check it has a valid value!
    if (!(weight_metric %in% c("Accuracy", "F1"))) {
      stop("'weight_metric' must be 'Accuracy' or 'F1'!")
    }
  }
  
  # [1] Get the data & dimensions of train & test folds! -----------------------
  # 1-1 Load Data & the names of the single blocks in the data!
  data <- load_data_extract_block_names(path = data_path, response = response)
  
  # 1-2 Get amount of Obs. we need for equally sized train folds 
  #     double it, as it was originally for 4 not 2 folds!
  obs_per_fold <- get_obs_per_fold(data = data$data)
  obs_per_fold$amount_train_fold <- obs_per_fold$amount_train_fold * 2 
  
  # [2] Shuffle the data & create lists to save the results in -----------------
  # 2-1 Shuffle IDs of data, we use to split data later on!
  set.seed(seed)
  fold_ids <- sample(nrow(data$data), nrow(data$data), replace = FALSE)
  
  # 2-2 Create empty lists to store results in!
  # 2-2-1 Full TestSet
  full <- list()
  
  # 2-2-2 TestSet with 1 missing omics-block
  miss1_A <- list(); miss1_B <- list()
  
  # 2-2-3 Single BlockTestSet
  single_A <- list(); single_B <- list(); single_CL <- list() 
  
  # [3] Assign which of the omics blocks should belong together! ---------------
  #     Randomly sample the letters "cnv", "rna", "mutation" & "mirna" & 
  #     bind the first 2 together and the last 2!
  set.seed(seed)
  blocks_together <- sample(c("cnv_block", "rna_block", "mutation_block", "mirna_block"),
                            size = 4, replace = F)
  
  # 3-1 Bind the colnams of the first 2 sampled blocks!
  block_A_names <- c(data$block_names[[which(names(data$block_names) == blocks_together[1])]],
                     data$block_names[[which(names(data$block_names) == blocks_together[2])]])
  
  # 3-2 Bind the colnames if the last 2 sampled blocks!
  block_B_names <- c(data$block_names[[which(names(data$block_names) == blocks_together[3])]],
                     data$block_names[[which(names(data$block_names) == blocks_together[4])]])
  
  # [3] Start the CV [5-fold per default!] -------------------------------------
  for (i in 0:4) {
    
    print(paste0("FOLD: ", as.character(i + 1), "/5 -------------------------"))
    
    # [1] Get TestSet from 'data', by taking the first 'obs_per_fold$amount_test' 
    #     IDs in 'fold_ids'
    test_ids <- fold_ids[((i * obs_per_fold$amount_test) + 1):(((i + 1) * obs_per_fold$amount_test))]
    test_df  <- data$data[test_ids,]
    
    # [2] Get the TrainSet from 'data' [= IDs not in TestSet] 
    train_ids <- fold_ids[-which(fold_ids %in% test_ids)]
    train_df  <- data$data[train_ids,]
    
    # [3] Induce blockwise missingness [SCENARIO_4] 
    # 3-1 Sample equally sized 'observed' blocks [according to SCENARIO_2]
    set.seed(seed + i)
    observed_blocks <- sample(c(rep("fold1", obs_per_fold$amount_train_fold), 
                                rep("fold2", obs_per_fold$amount_train_fold)),
                              obs_per_fold$amount_train, replace = FALSE)
    
    # 3-2 Now collect the different blocks from the different folds, where 
    #     the feas have been observed!
    # 3-2-1 Clin Block
    clin_block <- train_df[, c(response, data$block_names$clin_block)]
    
    # 3-2-2 Block A
    block_A <- train_df[which(observed_blocks == "fold1"), c(response, block_A_names)]
    
    # 3-2-3 Block B
    block_B <- train_df[which(observed_blocks == "fold2"), c(response, block_B_names)]
    
    # [4] Fit RFs on the seperate blocks
    #     'Formula' for the fitting has form of 'response ~ .' 
    #      --> use all cols as features except for the 'response'
    # 4-1 Clincal Block
    RF_clin <- rfsrc(formula = as.formula(paste(eval(response), "~ .")),
                     data = clin_block, ntree = num_trees, mtry = mtry, 
                     nodesize = min_node_size, samptype = "swr",
                     seed = seed)
    
    # 4-1-2 Block A
    RF_A <- rfsrc(formula = as.formula(paste(eval(response), "~ .")),
                  data = block_A, ntree = num_trees, mtry = mtry, 
                  nodesize = min_node_size, samptype = "swr",
                  seed = seed)
    
    # 4-1-3 Block B
    RF_B <- rfsrc(formula = as.formula(paste(eval(response), "~ .")),
                  data = block_B, ntree = num_trees, mtry = mtry, 
                  nodesize = min_node_size, samptype = "swr",
                  seed = seed)
    
    # 4-2 Collect all single RFs to the Forest!
    Forest <- list(RF_clin, RF_A, RF_B)
    
    # [5] Get (unweighted) predicitons on the different testsets!
    # 5-1 TestSet, when there is only one observed block! 
    #     Also done to extract the predictive performance of the different RFs,
    #     that were fitted to the single blocks! 
    #     As testset has 1 block only, we don't need to weight the different 
    #     blockwise RFs, as only one single block is used to obtain predicitons!
    #     --> the resulting performance will be used as weight for the remaining
    #         testsettings!
    print("Evaluation TestSet w/ 1 observed block only------------------------")
    print("Obtain performance & corresponding weights on the single blocks!")
    
    single_A[[i + 1]]  <- do_evaluation_rfsrc(Forest = Forest, weights = c(1, 1, 1), 
                                              testdata = test_df[,-which(colnames(test_df) %in% c(block_B_names,
                                                                                                  data$block_names$clin_block))])
    single_B[[i + 1]]  <- do_evaluation_rfsrc(Forest = Forest, weights = c(1, 1, 1),  
                                              testdata = test_df[,-which(colnames(test_df) %in% c(block_A_names,
                                                                                                  data$block_names$clin_block))])
    single_CL[[i + 1]]  <- do_evaluation_rfsrc(Forest = Forest, weights = c(1, 1, 1),  
                                               testdata = test_df[,-which(colnames(test_df) %in% c(block_B_names,
                                                                                                   block_A_names))])
    
    # 5-1-2 Extract the weights from the single block performance - weighted is true!
    weights <- c()
    if (weighted) {
      weights <- c(weights, single_CL[[i + 1]][[weight_metric]])
      weights <- c(weights, single_A[[i + 1]][[weight_metric]])
      weights <- c(weights, single_B[[i + 1]][[weight_metric]])
    } else {
      weights <- rep(1, times = 5)
    }
    
    # 5-1-3 Norm the weights
    weights <- weights / sum(weights)
    
    print("Extracted weights for different blocks [CL, A, B] & normalized them")
    print(weights)
    
    
    # 5-2 Full TestSet
    print("Evaluation full TestSet -------------------------------------------")
    full[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, testdata = test_df, 
                                         weights = weights) # FULL TESTSET
    
    # 5-3 TestSet, where one of the omics blocks is missing!
    print("Evaluation TestSet w/ 1 missing omics block------------------------")
    miss1_A[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights,
                                            testdata = test_df[,-which(colnames(test_df) %in% block_A_names)])
    
    miss1_B[[i + 1]] <- do_evaluation_rfsrc(Forest = Forest, weights = weights, 
                                            testdata = test_df[,-which(colnames(test_df) %in% block_B_names)])
    
  }
  
  # [4] Return the results & settings of parameters used to do CV! -------------
  res_all <- list("full" = full,
                  "miss1_A" = miss1_A, "miss1_B" = miss1_B,
                  "single_A" = single_A, "single_B" = single_B,
                  "single_CL" = single_CL)
  
  # Collect the Settings, used to do the CV!
  settings <- list("data_path"     = data_path,
                   "response"      = response, 
                   "seed"          = seed, 
                   "num_trees"     = num_trees,
                   "mtry"          = mtry, 
                   "min_node_size" = min_node_size,
                   "blocks_together"    = blocks_together) # the first and last two were merged to a single block!
  
  # Return both lists!
  return(list("res_all" = res_all, 
              "settings" = settings))
}

# Run a example and check the results!                                       ----
start_time <- Sys.time()
a1 <- do_CV_NK_setting1(num_trees = as.integer(250),
                        data_path = "data/external/Dr_Hornung/Data/ProcessedData_subsets/seed_1234/KIRC_Subset.RData",
                        weighted = TRUE, weight_metric = "F1")
end_time <- Sys.time()
a1_time <- end_time - start_time 

start_time <- Sys.time()
a2 <- do_CV_NK_setting2(num_trees = as.integer(250),
                        data_path = "data/external/Dr_Hornung/Data/ProcessedData_subsets/seed_1234/KIRC_Subset.RData",
                        weighted = TRUE, weight_metric = "F1")
end_time <- Sys.time()
a2_time <- end_time - start_time 

start_time <- Sys.time()
a3 <- do_CV_NK_setting3(num_trees = as.integer(250),
                        data_path = "data/external/Dr_Hornung/Data/ProcessedData_subsets/seed_1234/KIRC_Subset.RData",
                        weighted = TRUE, weight_metric = "F1")
end_time <- Sys.time()
a3_time <- end_time - start_time 

start_time <- Sys.time()
a4 <- do_CV_NK_setting4(num_trees = as.integer(250),
                        data_path = "data/external/Dr_Hornung/Data/ProcessedData_subsets/seed_1234/KIRC_Subset.RData",
                        weighted = TRUE, weight_metric = "F1")
end_time <- Sys.time()
a4_time <- end_time - start_time


res_all <- list(a1, a2, a3, a4)
save(res_all, file = "./docs/CV_Res/gender/Norbert_final_subsets/1312_subset.RData")
