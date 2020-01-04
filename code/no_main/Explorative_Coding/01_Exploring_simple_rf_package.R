"Research File - old scripts I used to understand and adjust code"
# Grow a tree, understand how it is built up and get predictons afer train  ----
# [0] Functions needed to create trees and get predicitons etc.
" ONLY DONE FOR CLASSIFICATION (PROBABILITY) not yet for SURVIVAL!!!
 orded_split & partition need to be inventigated as 'unordered_factors' argument!
"
# [1] Set Arguments - easier for investigations!
df <- read.csv2("./data/external/example_data/iris_example.csv", 
                stringsAsFactors = T)
formula           = Species ~ .
data              = df
num_trees         = as.integer(2)
mtry              = as.integer(2)
min_node_size     = as.integer(10)
replace           = TRUE
splitrule         = NULL
num_threads       = as.integer(1)
unordered_factors = "order_split"

# [2] Get 2 trees, with the same arguments!
TREES <- simpleRF(formula = Species ~ ., data = df, num_trees = as.integer(5),
                  mtry = as.integer(2), min_node_size = as.integer(10), 
                  replace = TRUE,  splitrule = NULL, num_threads = as.integer(1), 
                  unordered_factors = "order_once")

# [3] Inspect the trees!
TREES[[1]] # --> has all the information a tree needs!

TREES[[1]]$grow(replace = TRUE)
TREES[[2]]$grow(replace = TRUE)
# --> Sketch the trees based on "child_nodeIDs", "split_varIDs" & "split_values"
TREES[[1]]
TREES[[2]]

# [4] Get the Predicion of 3 different datapoints!
# [4-1] Get TestData and prepare it for prediciton [needs to do it on one tree!]
test_data <- df[c(49, 98, 147),]

test_data1 <- process_test_data(test_data = test_data, tree = TREES[[1]])

# [4-2] Get the true response % values of the testdata [check it manually!]
test_data1$data

# Get the Observations in Terminal Nodes
node_id = 15; TREE_ID = 1
TREES[[TREE_ID]]$data$subset(TREES[[TREE_ID]]$sampleIDs[[node_id]])
TREES[[TREE_ID]]$data$subset(TREES[[TREE_ID]]$sampleIDs[[node_id]], 1)

# [4-3] Get Prediciton from implemented method! Compare it with Truth + the 
#       manual prediction!
test_data1$data
TREES[[1]]$predict(test_data1) # --> ALL PREDICITONS OF TREE1 ARE CORRECT!

test_data1$data
TREES[[2]]$predict(test_data1) # --> ALL PREDICITONS ARE CORRECT!
# [4] PseudoFunction to prune a tree! --> needs 'tree' & 'data' only!       ----
# [1] Load the Trees from Chapter before
tree1 <- readRDS(file = "./data/external/example_data/example_tree_I.rds")
tree2 <- readRDS(file = "./data/external/example_data/example_tree_II.rds")

# [2] Create TestData from before and remove single varibales!
df         <- read.csv2("./data/external/example_data/iris_example.csv", 
                        stringsAsFactors = T)
# [2-1] DF w/ Sepal.Width missing!
df_test <- df[c(49, 98, 147), c("Sepal.Length","char", "Species",
                                "Petal.Width", "factor", "Petal.Length")]

# [2-2] Brings it to the layout needed and fills Variables w/ NA, if these
#       are not in the original DF!
test_data1 <- process_test_data(test_data = df_test, tree = tree1)
df_test; test_data1$data
tree = tree1
data = test_data1

# [3-1] Do we need to prune the tree?!
# [3-1-1] Get the Splitvariables the tree uses!
used_split_var       <- tree$split_varIDs[-which(is.na(tree$split_varIDs))]
used_split_var_names <- colnames(tree$data$data)[used_split_var]

# [3-1-2] Get the variables with missing data
missing_var       <- sapply(used_split_var_names, 
                            FUN = function(x) any(is.na(data$data[,x])))
missing_var_names <- used_split_var_names[missing_var]

# [3-1-3] Check whether we need to prune the tree
if (length(missing_var_names) < 1) {
  writeLines("Tree doesn't need to be pruned.\nNone of the used split variables is missing in test_data")
  stop()
}

# [3-1-4] Check whether we can prune the tree
if (any(missing_var_names %in% used_split_var_names[1])) {
  writeLines("Tree can not do predicitons, as the first splitVar the tree uses is missing in test_data")
  tree$child_nodeIDs[[1]] <- "pruned"
  stop()
}

# [3-2] PRUNE - If we have reached this point, we NEED TP PRUNE THE TREE
# [3-2-1] Find all Nodes, that use a split_point, not avaible in our data!
# [3-2-1-1] Get the names of all split_variables used in the tree [incl position]
split_var_names_orderd <- colnames(tree$data$data)[tree$split_varIDs]

# [3-2-1-2] Get all nodes, that use a missing variable for splitting!
nodes_to_prune <- c()
for (miss_var_name_curr in missing_var_names) {
  curr_ <- which(split_var_names_orderd == miss_var_name_curr)
  nodes_to_prune <- c(nodes_to_prune, curr_)
}

# [3-2-2] Set the 'child_nodeID' status of all the 'nodes_to_prune' to 'pruned'
for (curr_prune in nodes_to_prune) {
  tree$child_nodeIDs[[curr_prune]] <- NA
}

tree$predict(test_data1) 
# --> Works, needs to be implemented in the tree class!
#   --> OB8: tree was before the functionality was extended! 
#     --> works with NA & old Method
#       --> Proof of concept seems to fit!

# [5] Create the Trees and grow them                                        ----
#       Text extended functionality of the pruning
df <- read.csv2("./data/external/example_data/iris_example.csv", stringsAsFactors = T)

#       Create the trees
TREES <- simpleRF(formula = Species ~ ., data = df, num_trees = as.integer(5),
                  mtry = as.integer(2), min_node_size = as.integer(10), 
                  replace = TRUE,  splitrule = NULL, num_threads = as.integer(1), 
                  unordered_factors = "order_once")

#        Grow the trees
TREES[[1]]$grow(replace = TRUE)
TREES[[2]]$grow(replace = TRUE)

#        Save the trees
saveRDS(TREES[[1]], file = "./data/external/example_data/example_tree_I.rds")
saveRDS(TREES[[2]], file = "./data/external/example_data/example_tree_II.rds")

# [4-2] Check whether predicitons are correct!
#       Load the Trees
tree1 <- readRDS(file = "./data/external/example_data/example_tree_I.rds")
tree2 <- readRDS(file = "./data/external/example_data/example_tree_II.rds")

#       Get the Data
df             <- read.csv2("./data/external/example_data/iris_example.csv", 
                            stringsAsFactors = T)

# Check Predictions for fully observed!
test_data_full <- process_test_data(test_data = df[c(49, 98, 147),], tree = tree1)

test_data_full
tree1$predict(test_data_full) # all predicitons correct!
tree2$predict(test_data_full) # all predicitons correct!

# Check Predictions for data w/ missing values!
test_data_miss <- process_test_data(test_data = df[c(49, 98, 147), -which(colnames(df) %in% c("Petal.Width", 
                                                                                              "Sepal.Length"))], 
                                    tree = tree1)

tree1$prune(test_data = test_data_miss)
tree1$predict(test_data_miss)

tree2$prune(test_data = test_data_miss)
tree2$predict(test_data_miss)
# [6] Handle the trees as they were a forest!                               ----
"Pseudofunction from [4] already implemented!
 --> get aggregate predictions from an ensemble of trees!
"
# [0] Load Data
df <- read.csv2("./data/external/example_data/iris_example.csv", 
                stringsAsFactors = T)

# [2] Get 10 trees w/ same arguments
TREES <- simpleRF(formula = Species ~ ., data = df, num_trees = as.integer(10),
                  mtry = as.integer(5), min_node_size = as.integer(10), 
                  replace = TRUE,  splitrule = NULL, num_threads = as.integer(1), 
                  unordered_factors = "order_once")

# [3] Grow all of the trees
TREES <- mclapply(TREES, function(x) {
  x$grow(replace = TRUE)
  x
}, mc.cores = 1)

# [4] Get Predicitons of single trees & the aggregate of them!
# [4-1] Create TestData
test_data <- df[c(1:5, 50:55, 105:110), ]
test_data_usable <- process_test_data(tree = TREES[[1]], test_data = test_data)

# [4-2] Predicitions from all single trees!
predictions <- mclapply(TREES, function(x) {
  res <- x$predict(test_data_usable)
  res
}, mc.cores = num_threads)

# [4-3] From the single Trees, we form an aggregated prediction [w/ probs & class preds]
aggregated_predictions_prob <- mclapply(seq(nrow(test_data)), function(x) {
  'For each observation aggregate the prediciton from all single trees!'
  
  # Extract predciton for observation 'x' from all tree-predcitions!
  preds_obs_x <- sapply(predictions, FUN = function(y) y[,x])
  
  # Create probabilites from all the tree predicitons!
  preds_prob <- apply(preds_obs_x, MARGIN = 1, 
                      FUN = function(x) sum(x)/length(predictions))
  preds_prob
}, mc.cores = 1)

aggregated_predictions_class <- unlist(
  mclapply(aggregated_predictions_prob, FUN = function(x) {
    "From the class probabilites extract the class prediction itself"
    class_highest_prob <- names(x)[which(x == max(x))]
    class_highest_prob
  })
)

# [4-4] Get the true classes
true_classes <- test_data_usable$data$Species

# [5] Get the Metrics to rate the performance!
confmat <- caret::confusionMatrix(data = as.factor(aggregated_predictions_class), 
                                  reference =  test_data_usable$data$Species,
                                  positive = "setosa")

confmat$byClass[grep("setosa", rownames(confmat$byClass)),]
confmat$byClass[grep("versicolor", rownames(confmat$byClass)),]
confmat$byClass[grep("virginica", rownames(confmat$byClass)),]
