"Script for the examplary usage of the pruning method!"

# Load all functions, objects etc. 
setwd("C:/Users/kuche_000/Desktop/Master_Omicsdaten/")
source("./code/04_simpleRF_adaption.R")
require(parallel)
library(checkmate)
library(ROCR)
library(caret)

# [1] Load Data
df <- read.csv2("./data/external/example_data/04_iris_example.csv", 
                stringsAsFactors = T)

# [2] Only run this code if you want to create new TREES and not use TREES grown
#     already, which will be loaded in the other chapters!
' DO NOT RUN LOAD THE TREES TRAINED ALREADY! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
# [2] Get 10 trees w/ same arguments
TREES <- simpleRF(formula = Species ~ ., data = df, num_trees = as.integer(10),
                  mtry = as.integer(3), min_node_size = as.integer(10), 
                  replace = TRUE,  splitrule = NULL, unordered_factors = "ignore")

# [3] Grow all of the trees & save them!
TREES <- mclapply(TREES, function(x) {
  x$grow(replace = TRUE)
  x
}, mc.cores = 1)
saveRDS(TREES, "./data/external/example_data/04_TREES.rds") # DO NOT OVERRIDE!
DO NOT RUN LOAD THE TREES TRAINED ALREADY! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! '

# No Pruning Needed                                                         ----
TREES <- readRDS("./data/external/example_data/04_TREES.rds")
TREES <- TREES[c(-4:-length(TREES))] # keep first 3 --> easier to understand!

# [N-4] Create a testset w/o missing features in correct layout etc.
test_data_full        <- df[c(1,75, 134),]
test_data_full_usable <- process_test_data(tree = TREES[[1]], test_data = test_data_full) 

# [N-5] Get Prediciton from the model
preds <- get_pruned_prediction(trees = TREES, test_set = test_data_full_usable)
preds$Class <- factor(preds$Class, levels = levels(TREES[[1]]$data$data$Species))

# [N-6] Evaluate the Model
true_classes <- test_data_full_usable$data$Species
confmat      <- caret::confusionMatrix(data = preds$Class,
                                       reference = true_classes)

# [N-7] Get the Metrics from the confusion matrix!
# [N-7-1] Extract´the metrics_we_care for all possible response classes
metrics_we_care <- c("Sensitivity", "Specificity", "F1", "Precision", "Recall", "Balanced Accuracy")
setosa <- confmat$byClass[grep("setosa", rownames(confmat$byClass)),][metrics_we_care]
versic <- confmat$byClass[grep("versicolor", rownames(confmat$byClass)),][metrics_we_care]
virgin <- confmat$byClass[grep("virginica", rownames(confmat$byClass)),][metrics_we_care]

# [N-7-2] Set Metrics w/ NA to 0 [na if the class wasn't predicted anytime!]
setosa[is.na(setosa)] <- 0
versic[is.na(versic)] <- 0
virgin[is.na(virgin)] <- 0

# [N-7-3] Compute the average metric over all classes!
metrics_overall <- sapply(metrics_we_care, 
                          FUN = function(x) mean(c(setosa[x], versic[x], virgin[x])))

# [8] Print Results:
confmat$overall["Accuracy"]
metrics_overall


# Pruning Needed                                                            ----
TREES <- readRDS("./data/external/example_data/04_TREES.rds")
TREES <- TREES[c(-4:-length(TREES))] # keep first 3 --> easier to understand!

# [W-4] Create the testdata w/ a missing feature variable!
test_data        <- df[c(1,75, 134), -which(colnames(df) == "Petal.Width")]
test_data_usable <- process_test_data(tree = TREES[[1]], test_data = test_data)

# [W-5] Get the pruned predicitons!
preds       <- get_pruned_prediction(trees = TREES, test_set = test_data_usable)
preds$Class <- factor(preds$Class, levels = levels(TREES[[1]]$data$data$Species))

# [7] Evaluate the Model on test_data
# [7-1] Create a ConfMatrix and compare Truth w/ Predicitons!
true_classes <- test_data_usable$data$Species
confmat      <- caret::confusionMatrix(data = preds$Class, reference = true_classes)

# [7-2] Get the Metrics from the confusion matrix!
metrics_we_care <- c("Sensitivity", "Specificity", "F1", "Precision", "Recall", "Balanced Accuracy")
setosa <- confmat$byClass[grep("setosa", rownames(confmat$byClass)),][metrics_we_care]
versic <- confmat$byClass[grep("versicolor", rownames(confmat$byClass)),][metrics_we_care]
virgin <- confmat$byClass[grep("virginica", rownames(confmat$byClass)),][metrics_we_care]

metrics_overall <- sapply(metrics_we_care, 
                          FUN = function(x) mean(c(setosa[x], versic[x], virgin[x])))

# [8] Print Results:
confmat$overall["Accuracy"]
metrics_overall

# --> Works as expected - manual checks confirmed!
# LITTLE HELPFUNCTION                                                       ----
t_ree = 3; no_de = 3
prop.table(table(TREES[[t_ree]]$data$subset(TREES[[t_ree]]$sampleIDs[[no_de]], 1)))
