"Script to implement crossValidation in our special setting!
 Explorative part and basis for the implementation!
"

# Load Funcitons, Classes, librarys & set the WD!
setwd("C:/Users/kuche_000/Desktop/Master_Omicsdaten/")
source("./code/04_simpleRF_adaption.R")

# [1] Load the blockwise Omics-Data & create a single DF                    ----
# [1-1] Raw Data
load("./data/external/Dr_Hornung/Data/ProcessedData/KIRC.Rda")

# [1-2] Subsett the single Blocks
resp     <- clin$gender
clin     <- clin[-which(colnames(clin) == "gender")]
cnv      <- cnv[,sample(ncol(cnv), round(ncol(cnv) * 0.05))]
rna      <- rna[,sample(ncol(rna), round(ncol(rna) * 0.05))]
mirna    <- mirna[,sample(ncol(mirna), round(ncol(mirna) * 0.05))]
mutation <- mutation[,sample(ncol(mutation), round(ncol(mutation) * 0.05))]

# [1-3] Create a single DF from all blocks [rename one block else same colnames]
#       & extract the colnames of all blocks seperatly!
colnames(rna) <- paste0(colnames(rna), "_rna");  rna_block <- colnames(rna)
clin_block    <- colnames(clin);                 cnv_block <- colnames(cnv)
mirna_block   <- colnames(mirna);           mutation_block <- colnames(mutation)

# [1-4] Bind all the blocks to a single DF & set response to factor! 
df <- cbind(resp, clin, cnv, rna, mirna, mutation); df$resp <- as.factor(df$resp)

# [2] Split the data into test and train data                               ----
# [2-1] Want same amount of obs. in all training folds! So the amount of trainig
#       observations we use must be dividable by 4! If 80% of nrow(df) is not 
#       dividable by 4, decrease amount_test by 1 & increase amount_train by 1!
#       --> each train fold same size & test fold can be used 
amount_train <- floor(4/5 * nrow(df)); amount_test <- ceiling(1/5 * nrow(df)) 
while ((amount_train %% 4) != 0) {
  amount_train      <- amount_train + 1
  amount_train_fold <- amount_train / 4
  amount_test       <- amount_test - 1
}

# [2-2] Print amount of observations we will use for the different folds!
writeLines(paste("Each TestFold will hold", amount_test,  "Test-Observations!"))
writeLines(paste("The remaning", nrow(df) - amount_test, "Obs. are used for Training",
                 "-", amount_train_fold, "Observations per Trainingsfold"))

# [2-3] 'Randomly' shuffle the order & assign Obs. into test & train!
set.seed(1312)
fold_ids <- sample(nrow(df), size = nrow(df), replace = FALSE)

# [3] Start the k-fold-CV, and use different folds!                         ----
i = 0
# for (i in 0:4) {

# [3-1] split data to test and train!
# [3-1-1] Test-Set w/ 'amount_test' observations + save it in 'data/temp'
test_ids <- fold_ids[((i * amount_test) + 1):(((i + 1) * amount_test))]
test_df  <- df[test_ids,] 
saveRDS(test_df, "./data/interim/tmp_data/tmp_test_set.rds")

# [3-1-2] Train-Set w/ 'amount_train' observations
train_ids <- fold_ids[-which(fold_ids %in% test_ids)]
train_df  <- df[train_ids,]

# [3-2] Induce blockwise missingness for train data!    
# [3-2-1] Example Case, where each obs, has 'clin' & one single block! 
#         Sample randomly 'amoun_train' of each of the observed blocks
observed_blocks <- sample(c(rep("Clin, A", amount_train_fold), 
                            rep("Clin, B", amount_train_fold),
                            rep("Clin, C", amount_train_fold), 
                            rep("Clin, D", amount_train_fold)),
                          amount_train, replace = FALSE)

# [3-2-2] Split the train DF into single folds - w/ different observed blocks!
block1 <- train_df[which(observed_blocks == "Clin, A"), 
                   c("resp", clin_block, cnv_block)]
block2 <- train_df[which(observed_blocks == "Clin, B"), 
                   c("resp", clin_block, rna_block)]
block3 <- train_df[which(observed_blocks == "Clin, C"), 
                   c("resp", clin_block, mutation_block)]
block4 <- train_df[which(observed_blocks == "Clin, D"), 
                   c("resp", clin_block, mirna_block)]

# [3-3] Based on each [fully observed block] we fit a RandomForest & check, 
#       whether each tree was grown correctly - if not, it will get buggy...
# [3-3-1] BLOCK1
TREES1 <- simpleRF(formula = resp ~ ., data = block1, 
                   num_trees = as.integer(10), mtry = as.integer(10), 
                   min_node_size = as.integer(5), replace = TRUE,  
                   splitrule = NULL, unordered_factors = "ignore")
TREES1 <- mclapply(TREES1, function(x) {
  x$grow(replace = TRUE)
  x
}, mc.cores = 1)
# Check, that the TREES were grown completly
if (any(unlist(lapply(TREES1, FUN = function(x) length(x$child_nodeIDs) == 0)))) {
  stop("TREES were not correctly grown")
}

# [3-3-2] BLOCK2
TREES2 <- simpleRF(formula = resp ~ ., data = block2, 
                   num_trees = as.integer(10), mtry = as.integer(10), 
                   min_node_size = as.integer(5), replace = TRUE,  
                   splitrule = NULL, unordered_factors = "ignore")
TREES2 <- mclapply(TREES2, function(x) {
  x$grow(replace = TRUE)
  x
}, mc.cores = 1)
# Check, that the TREES were grown completly
if (any(unlist(lapply(TREES2, FUN = function(x) length(x$child_nodeIDs)  == 0)))) {
  stop("TREES were not correctly grown")
}

# [3-3-3] BLOCK3
TREES3 <- simpleRF(formula = resp ~ ., data = block3, 
                   num_trees = as.integer(10), mtry = as.integer(10), 
                   min_node_size = as.integer(5), replace = TRUE,  
                   splitrule = NULL, unordered_factors = "ignore")
TREES3 <- mclapply(TREES3, function(x) {
  x$grow(replace = TRUE)
  x
}, mc.cores = 1)
# Check, that the TREES were grown completly
if (any(unlist(lapply(TREES3, FUN = function(x) length(x$child_nodeIDs)  == 0)))) {
  stop("TREES were not correctly grown")
}

# [3-3-4] BLOCK4
TREES4 <- simpleRF(formula = resp ~ ., data = block4, 
                   num_trees = as.integer(10), mtry = as.integer(10), 
                   min_node_size = as.integer(5), replace = TRUE,  
                   splitrule = NULL, unordered_factors = "ignore")
TREES4 <- mclapply(TREES4, function(x) {
  x$grow(replace = TRUE)
  x
}, mc.cores = 1)
# Check, that the TREES were grown completly
if (any(unlist(lapply(TREES4, FUN = function(x) length(x$child_nodeIDs)  == 0)))) {
  stop("TREES were not correctly grown")
}

# [3-4] Add all to the 'Forest' list!
Forest <- list(TREES1, TREES2, TREES3, TREES4)
saveRDS(Forest, "./data/interim/tmp_model/tmp_forrest.rds")

# [5] Start Testing! For each TestSituation, we need to reload the trees & data,
#     as pruning etc. will change the trees & for test situations we need to 
#     prune testdata!

# [5-1] Predictions - FULLY OBSERVED TESTSET W/ ALL BLOCKS!
# [5-1-1] Reload the data
Forest  <- readRDS("./data/interim/tmp_model/tmp_forrest.rds")
test_df <- readRDS( "./data/interim/tmp_data/tmp_test_set.rds")

# [5-1-2] Get Predicitons from the trees from the different blocks!
test_data   <- process_test_data(tree = Forest[[1]][[1]], test_data = test_df) 
TREES1_pred <- get_pruned_prediction(trees = Forest[[1]], test_set = test_data)

test_data   <- process_test_data(tree = Forest[[2]][[1]], test_data = test_df) 
TREES2_pred <- get_pruned_prediction(trees = Forest[[2]], test_set = test_data)

test_data   <- process_test_data(tree = Forest[[3]][[1]], test_data = test_df) 
TREES3_pred <- get_pruned_prediction(trees =  Forest[[3]], test_set = test_data)

test_data   <- process_test_data(tree = Forest[[4]][[1]], test_data = test_df) 
TREES4_pred <- get_pruned_prediction(trees = Forest[[4]], test_set = test_data)

# [5-1-3] Look into it exemplary!
TREES1_pred$Probs[[1]]; TREES2_pred$Probs[[1]]; TREES3_pred$Probs[[1]]; TREES4_pred$Probs[[1]]
test_df$resp[1]

# [5-1-4] Aggregate the predicitons from the different trees!
all_forrest_preds_probs_class_0 <- sapply(1:nrow(test_df), FUN = function(x) {
  # Class1
  prob_class1 <- mean(c(TREES1_pred$Probs[[x]][1], TREES2_pred$Probs[[x]][1],
                        TREES3_pred$Probs[[x]][1], TREES4_pred$Probs[[x]][1]))
  prob_class1
})

all_forrest_preds_class <- ifelse(all_forrest_preds_probs_class_0 >= 0.5, 0, 1)
all_forrest_preds_class <- factor(all_forrest_preds_class, 
                                  levels = levels(TREES1[[1]]$data$data$resp))

# [5-1-5] Get Metrics for the current setting!
#         Confusion Matrix
confmat <- caret::confusionMatrix(data      = all_forrest_preds_class, 
                                  reference = test_df$resp)

confmat$overall["Accuracy"]; confmat$byClass


# [5-2] Evaluation w/ a DF, that misses a block!
Forest  <- readRDS("./data/interim/tmp_model/tmp_forrest.rds")
test_df <- readRDS( "./data/interim/tmp_data/tmp_test_set.rds")

# Remove the 'rna' block from the test_df, to adjust data to our settings!
old     <- ncol(test_df)
test_df <- test_df[-which(colnames(test_df) %in% rna_block)]
print(paste("Removing of 'rna_block' --> ", old - ncol(test_df), "less features"))

# [5-2-2] Get Predicitons from the trees from the different blocks!
test_data   <- process_test_data(tree = Forest[[1]][[1]], test_data = test_df) 
TREES1_pred <- get_pruned_prediction(trees = Forest[[1]], test_set = test_data)

test_data   <- process_test_data(tree = Forest[[2]][[1]], test_data = test_df) 
TREES2_pred <- get_pruned_prediction(trees = Forest[[2]], test_set = test_data)

test_data   <- process_test_data(tree = Forest[[3]][[1]], test_data = test_df) 
TREES3_pred <- get_pruned_prediction(trees =  Forest[[3]], test_set = test_data)

test_data   <- process_test_data(tree = Forest[[4]][[1]], test_data = test_df) 
TREES4_pred <- get_pruned_prediction(trees = Forest[[4]], test_set = test_data)

# [5-1-3] Look into it exemplary!
TREES1_pred$Probs[[1]]; TREES2_pred$Probs[[1]]; TREES3_pred$Probs[[1]]; TREES4_pred$Probs[[1]]
test_df$resp[1]

# [5-1-4] Aggregate the predicitons from the different trees!
all_forrest_preds_probs_class_0 <- sapply(1:nrow(test_df), FUN = function(x) {
  # Class1
  prob_class1 <- mean(c(TREES1_pred$Probs[[x]][1], TREES2_pred$Probs[[x]][1],
                        TREES3_pred$Probs[[x]][1], TREES4_pred$Probs[[x]][1]), 
                      na.rm = TRUE)
  prob_class1
})

all_forrest_preds_class <- ifelse(all_forrest_preds_probs_class_0 >= 0.5, 0, 1)
all_forrest_preds_class <- factor(all_forrest_preds_class, 
                                  levels = levels(TREES1[[1]]$data$data$resp))

# [5-1-5] Get Metrics for the current setting!
#         Confusion Matrix
confmat <- caret::confusionMatrix(data      = all_forrest_preds_class, 
                                  reference = test_df$resp)

confmat$overall["Accuracy"]; confmat$byClass