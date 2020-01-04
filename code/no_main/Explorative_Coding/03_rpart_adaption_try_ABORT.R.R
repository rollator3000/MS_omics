"
Script to play around and try to implement the RF Adaption of Roman by myself!

--------------------------------> OLD VERSION! <-------------------------------- 
--> DIDN'T WORK OUT IN THE END!
--> STILL IN REPO AS IT MIGHT BE USEFUL WHILE TRYING TO ADAPT IT IN AN OTHER WAY

Quick Explanation:
  - We want to adjust RF, so it can train & predit on block wise missing data
  - M different blocks of training data, where not all blocks are avaible for
    all observations! 
    [no NA's in the blocks themself, only whole blocks are missing]

--------------------------------> OLD VERSION! <-------------------------------- 
TRAINNIG:
  - for i in [1,...,ntree]:
    - draw a random block of the M blocks
      - draw a subsample from this block[sample Obs. & Feas.]
        - Grow a tree and consider sqrt(p_m) features 
          as split candidates at each splitpoint!
        - Grow the tree until terminal nodes!
    - save the decision tree and its splitpoints etc.

--> Grow ntree-trees that were grown on a subsetted observations & features of
    each Block

--------------------------------> OLD VERSION! <-------------------------------- 
PREDICTING:
  - TestData might contain features/blocks that were not avaible in the train
    data itself 
      --> remove all test-features/ blocks , that dom't appear in the trainset
  
  Predicton:
    - take the observation we want to predict the outcome for
    - remove all trees from the forrest, that use a feature as first split 
      criterion 
    - for all remaining trees, walk down their branches until a node is found,
      that is not existent in the observation 
        --> cut the tree at the these nodes and use the new terminal nodes for 
            prediction [e.g. majority vote, mean,...]
"
# Set WD, Load Packages & define Functions
# Set the WD to the main folder, that holds all things for "Master OmicsDaten"
setwd("C:/Users/kuche_000/Desktop/Master_Omicsdaten/")

library(rpart)
library(rpart.plot)
library(rpart.utils)
library(rattle)
library(RColorBrewer)
library(partykit)
library(dplyr)
library(checkmate)
library(Rcpp)
library(foreach)
library(doParallel)

# Define Functions (see doc_string in function itself for explanation!)
name_rules_df <- function(rules_df) {
  "Function to assign meaningful colnames to a set of rules from a decision tree!
  
   Args:
    - rules_df (rpart.rules, data.frame): Set of splitting rules obtained from a 
                                          decison tree! (see Example)
                                          
   Return:
    - vector of columnnames for the rules_df!
    
   Example:
    tree       <- rpart::rpart(Species ~ ., data = iris[1:100,])
    tree_rules <- rpart.plot::rpart.rules(tree, cover = T, nn = T, roundint = F)
    colnames(tree_rules) <- name_rules_df(tree_rules)
  "
  # [0] Check Inputs
  # [0-1] Tree from the right class
  if (!("rpart.rules" %in% class(rules_df) | "data.frame" %in% class(rules_df))) {
    stop("'rules_df' argument is not of tpyes: 'rpart.rules' | 'data.frame'")
  }
  
  # [1] Extract the names of the splitting variables!
  # [1-1] Create Colnames for the rules set!
  new_rules_colnames   <- c("nn", "response", "when")
  max_appliable_rules  <- round((ncol(rules_df) - 4) / 4)
  
  # Create meaningful colnames in the loop!
  for (i in seq_len(max_appliable_rules)) {
    if (i < max_appliable_rules) {
      var_names <- paste0(c("Var_", "Oper_", "Val_", "&_"), i)
      new_rules_colnames <- c(new_rules_colnames, var_names)
    }
    if (i == max_appliable_rules) {
      var_names <- paste0(c("Var_", "Oper_", "Val_"), i)
      new_rules_colnames <- c(new_rules_colnames, var_names)
    }
  }
  return(c(new_rules_colnames, "Cover"))
}

find_split_var <- function(rules_df, row) {
  "Function to extract the next split Variables from a certain row
   in the splitting rules of a decision tree 
   [see 'rpart.plot::rpart.rules(tree, cover = T, nn = T, roundint = F)',
    where 'tree' is a fitted decisiontree]
    
  Args:
    - rules_df (rpart.rules, data.frame): Set of splitting rules obtained from a 
                                          decison tree! 
                                          [decisiontree: rpart::rpart(),
                                          rules_set    : rpart.plot::rpart.rules()]
    - row (integer)                     : Indicating the row we want to get all 
                                          splitvariables from!
                                          
  Return:
    - Vector of all Variabels in the current row!
  
  Example:
    tree       <- rpart(Species ~ ., data = iris[1:100,])
    tree_rules <- rpart.plot::rpart.rules(tree, cover = T, nn = T, roundint = F)
    colnames(tree_rules) <- name_rules_df(tree_rules)
    find_split_var(rules_df = tree_rules, row = 1)
  "
  check_int(row, lower = 1)
  check_data_frame(rules_df, any.missing = F)
  
  rules_df <- rules_df[row,]
  
  var_names <- sapply(grep("Var", colnames(rules_df)), 
                      FUN = function(x) rules_df[,x])
  
  return(var_names[var_names != ""])
}

subset_intervall_operator <- function(tree_rules) {
  "Function, that changes the layout of the tree_rules, so that splits of the 
   form 'x in [4.6, 7.9]' can be presented in a single cell!
   --> Mainly needed to interpret meanigful!
   
   Args:
    - tree_rules (rpart.rules, data.frame): Set of splitting rules obtained 
                                            from a decison tree! 
                                            [decisiontree: rpart::rpart(),
                                             rules_set   : rpart.plot::rpart.rules()]
  Return:
    - tree rules with - for each intervall split - 2 columns less dimensional df!
  "
  # [0] Check Inputs
  # [0-1] Tree from the right class
  if (!("rpart.rules" %in% class(tree_rules))) stop("rules_df must be of class 'rpart.rules'")
  
  # Which column contains the 'to'
  cols_w_to <- which(unlist(lapply(1:ncol(tree_rules), FUN = function(x) "to" %in% tree_rules[,x])))
  
  for (i in cols_w_to) {
    
    # Check the rows that contain the 'to'
    rows_w_to <- which(unlist(lapply(1:nrow(tree_rules), FUN = function(x) "to" %in% tree_rules[x,i])))
    
    # Replace the operator 'is' by 'in_intervall' and put 2 limits in the same col!
    tree_rules[rows_w_to, i - 2] <- 'in_intervall'
    
    # Put the Limits from the columns around "to" into the same column
    tree_rules[, i - 1] <- paste(tree_rules[, i - 1], tree_rules[, i + 1]) 
    
    # Now Remove the column w/ 'to' and the column behind
    tree_rules[,i]     <- NA
    tree_rules[,i + 1] <- NA
  }
  
  # Now Remove the Columns, that just contain NA
  tree_rules <- tree_rules[,-which(colSums(is.na(tree_rules)) >= nrow(tree_rules))] 
  
  return(tree_rules)
}

get_prediciton_from_tree_rules <- function(rules, observation) {
  "Function to generate a prediciton from a set of rules, obtained by a single 
   decision tree from rpart::rpart()! (see Example) 
   This is especially useful, if we prune our trees [and therefore its rules]!
   Only made for binary classification cases, that have numeric or factor 
   features!
   
   Args:
    - rules (rpart.rules, data.frame): Set of rules obtained from a decison tree.
                                       Use this set of rules to generate a
                                       prediciton! Observation needs to have the
                                       1. split argument in rules, else a 
                                       prediciton can not be generated!
                                       ['rpart()' used to fit decisiontree]
                                        
    - observation (row of dataframe) : Observation we want to generate a 
                                       prediction for!
   Return:
    - numeric: Probabiltiy for a positive class in the prediciton!
    
    Example:
        1)  tree = rpart::rpart(Species ~ ., data = iris[1:100,])
        2)  tree_rules = rpart.plot::rpart.rules(tree, cover = T, nn = T, 
                                                 roundint = F)
        3)  get_prediciton_from_tree_rules(rules = tree_rules, 
                                           observation = iris[100, 1:4])
  "
  # [0] ----- Inut Checks
  # [0-1] Check for right datatypes of the arguments!
  assert_data_frame(observation, max.rows = 1, min.rows = 1)
  
  if (!("rpart.rules" %in% class(rules) | "data.frame" %in% class(rules))) {
    stop("'rules' argument is not of tpyes: 'rpart.rules' | 'data.frame'")
  }
  
  # [1] ----- Extract Informations needed to create a prediciton!
  # [1-0-1] Modify Rules DF, if it contains 'to' [intervals]
  if (sum(rules == "to") >= 1) rules <- subset_intervall_operator(rules)
  
  # [1-1] Assign meaningful colnames to the rules DF:
  colnames(rules) <- name_rules_df(rules)
  
  # [1-2] Extract split Variables, used in the rules [from rpart::rpart()]
  used_split_var <- unique(unlist(
    sapply(grep("Var", colnames(rules)), FUN = function(x) unique(rules[,x]))
  ))
  if ("" %in% used_split_var) used_split_var <- used_split_var[used_split_var != ""]
  
  
  # [2] ----- Generate a prediction for 'observation' based on 'rules'
  # Send the test_obs over all rows [rules] to get a prediciton [one has to fit!]
  for (row in seq_len(nrow(rules))) {
    
    print(paste("ROW:", row))
    
    # Get all variable names and Column-Index [of rules] that are considerd as 
    # split in the current row!
    row_split_var       <- find_split_var(rules, row)          
    row_split_var_index <- which(rules[row,] %in% row_split_var)
    last_split_var_row  <- row_split_var[length(row_split_var)]
    prediction          <- NA
    
    # Loop over these Variables and check whether our X fits the needs!
    for (possible_split_var in row_split_var) {
      
      # [1] Get the column, in which the 'possible_split_var' is, in 'rules'!
      current_split_var_index <- which(rules[row, ] == possible_split_var)
      
      # [2] Get Value of the current split-variable, split-value & opertor from tree!
      curr_obs_    <- observation[possible_split_var]
      split_value  <- rules[row, current_split_var_index + 2]
      log_operator <- rules[row, current_split_var_index + 1]
      
      # [3] Create a logical statement from X_Value, split_value & log_operator
      #     & evaluate it- based on logical operator needed!
      
      # [3-1] If we have a factor comparison, we need to transform "curr_obs_" 
      #       from factor to character & check whether this value is a split value!
      if (log_operator == "is") {
        # covert from factor to character!
        curr_obs_    <- levels(curr_obs_[[1]])[as.numeric(as.character(curr_obs_))]
        
        # multiple split values, paste them together so we can have easier comparison!
        # If Log_operator contains 'or' we do a factor check
        if (grepl("or", split_value)) split_value <- strsplit(split_value, " or ")[[1]]
          curr_check <- curr_obs_ %in% c(split_value)
        } else if (log_operator == "in_intervall") {
          
          # Get the lower and upper bounds of the intervall and only keep numerics!
          split_values <- as.numeric(strsplit(split_value, split = " ")[[1]])
          if (any(is.na(split_values))) split_values <- split_values[-which(is.na(split_values))]
          
          # check whether our obs. is within the intervall!
          curr_check   <- curr_obs_ >= min(split_values) & curr_obs_ <= max(split_values)
        } else {
        curr_check <- paste(curr_obs_, log_operator, split_value)
        curr_check <- eval(parse(text = curr_check))
      }
      
      # [4] If it's not true break and jump to next row, as the current row uses
      #     a split value, that is not present in 'observation' 
      #     --> no prediction possible!
      if (!curr_check) break
      
      # [5] If the current possible_split_var is the last one in the current row
      #     - so that it has fullfilled all requirements - we can do a prediciton!
      if (possible_split_var == last_split_var_row) {
        prediction = as.numeric(rules$response[row])
      }
    }
    if (!is.na(prediction)) {
      print("prediciton found!")
      print(paste("Prob. pos. class:", prediction))
      break
    }
  }
  return(prediction)
}

prune_tree <- function(tree_rules, observation) {
  "Function to prune all branches of a tree [= tree_rules], that use split 
   variables, not existent in 'observation' 
   --> we can create meaningful predictons for the observation, w/ the pruned 
       tree! Except, when the 1.split arguemnt is not existent in 'observation' 
             then w/ this tree no predicitons can be done!

   !!!'tree_rules' are the extracted rules from a decisiontree fitted w/ rpart!!!
      >>     tree = rpart::rpart(Species ~ ., data = iris[1:99,])       <<
      >> rpart.plot::rpart.rules(tree, cover = T, nn = T, roundint = F) <<
   
   Args:
    - tree_rules (rpart.rules, data.frame): Set of rules obtained from a decison 
                                            tree. If the 1. split argument in 
                                            tree_rules is not existent in 
                                            'observation' the tree_rules can not
                                            be pruned!
                                            ['rpart' used to fit decisiontree]
                                            
    - observation (data.frame)            : single row of a data.frame, we want
                                            to prune the tree for! 
                                            
   Return:
    - pruned tree_rules (rpart.rules, data.frame), that uses only split
      variables, that are existent in 'observation'!
  "
  # [0] ----- Inut Checks
  # [0-1] Check for right datatypes of the arguments!
  assert_data_frame(observation, max.rows = 1, min.rows = 1)
  
  if (!("rpart.rules" %in% class(tree_rules) & "data.frame" %in% class(tree_rules))) {
    stop("'rules' argument is not of tpyes: 'rpart.rules' & 'data.frame'")
  }
  
  # [1] ----- Get all variables that are used for splitting, but not existent in 
  #           'obsergation' -> so we can prune the tree in these variables!
  # [1-1] Prepare the 'tree_rules'
  # [1-1-1] Modify Rules DF, if it contains 'to' [= intervals]
  if (sum(tree_rules == "to") >= 1) tree_rules <- subset_intervall_operator(tree_rules)
  
  # [1-1-2] Assign Meaningful colnames
  colnames(tree_rules) <- name_rules_df(tree_rules)
  
  # [1-2] Extract all Split-Variables used in the tree/ tree_rules!
  split_variables <- unlist(sapply(grep("Var", colnames(tree_rules)), 
                                   FUN = function(x) unique(tree_rules[,x])))
  if ("" %in% split_variables) split_variables <- split_variables[split_variables != ""]
  
  # [1-3] Find Variables used in the tree, but not avaible in observation!
  missing_variables <- split_variables[which(!(split_variables %in% colnames(observation)))]
  
  # [1-4] If no splitting variable is missing in 'observation' we don't need to
  #       prune the tree!
  if (length(missing_variables) < 1) {
    print("Tree doesn't use variables, not appearing in model to predict")
    return(tree_rules)
  }
  
  # [2] If tree contains split variables not in observation, we need to prune it!
  # [2-1] Get indices & names in the 'tree_rules' of the variables missing in '
  #       observation!
  missing_var_index <- apply(tree_rules, MARGIN = 2, FUN = function(x) any(x %in% missing_variables))
  missing_var_names <- colnames(tree_rules[missing_var_index])
  
  # [2-2] If the 1. split Variable of the tree is not in the observation, we can
  #       not prune the tree --> return NA incl. a warning!
  if ("Var_1" %in% missing_var_names) {
    warning("Tree not usable! First split_variable not avaible as feature for 'observation'") 
    return(NA)
  }
  
  # [2-3] Prepare the 'tree_rules' for pruning, by swithcing columns, types etc.
  tree_rules$nn       <- tree_rules$Cover
  tree_rules$nn       <- sapply(tree_rules$nn, 
                                FUN = function(x) unlist(strsplit(x, split = "%")))
  tree_rules$nn       <- as.numeric(tree_rules$nn) 
  tree_rules$response <- as.numeric(tree_rules$response)
  
  
  # [2-4] Start pruning of the tree - start pruning at the terminal nodes and 
  #       go up & up!
  for (miss_variable in rev(missing_var_names)) {
    
    # Every Row representing a single track from the decision tree. Any row,
    # that contains a split_var that is not in 'observation' needs to be pruned
    rows_to_prune <- which(tree_rules[,miss_variable] != "")
    
    # Create new terminal nodes, by looking where the current node is originally
    # coming from!
    tree_rules_curr_row <- data.frame(tree_rules)[rows_to_prune, ]
    
    # Cut Everything, that comes behind the split! And get Histroy of the nodes!
    cols_to_prune <- which(colnames(tree_rules_curr_row) == miss_variable)
    history_nodes <- data.frame(tree_rules_curr_row)[,c(1:(cols_to_prune - 2))]
    
    # Check, that all came the same way to the current 'miss_variable' -->
    # the ones belonging together have the same history!
    diff_historys <- nrow(unique(history_nodes[,3:ncol(history_nodes)]))
    
    # Cummulate all nodes, that have the same history!
    for (i in seq(1:diff_historys)) {
      
      # Get the unique Histroy 
      curr_histroy <- unique(history_nodes[i, 3:ncol(history_nodes)])
      
      # Find all rows, that have the same history so we can combine them!
      identical_rows <- t(sapply(seq(1:nrow(tree_rules_curr_row)),
                                 FUN = function(x) 
                                   tree_rules_curr_row[x, 3:(cols_to_prune - 2)] == curr_histroy))
      rows_w_same_hist <- which(apply(identical_rows, MARGIN = 1, 
                                      FUN = function(x) all(x)))
      
      # Extract all rows with identical history and merge them to calculate the
      # new terminal node of the pruned tree!
      tree_rules_same_hist <- tree_rules_curr_row[rows_w_same_hist, ]
      
      pruned_node          <- tree_rules_same_hist[1,]
      pruned_node$nn       <- sum(tree_rules_same_hist$nn)
      pruned_node$response <- sum(tree_rules_same_hist$response * tree_rules_same_hist$nn) / pruned_node$nn
      
      # fill up the cut off rules w/ ""
      pruned_node[,(cols_to_prune - 1):ncol(pruned_node)] <- ""
      
      # replace the original rules with the pruned version, by replacing the rows
      # that use the missing split_var + same history!
      tree_rules[rows_to_prune[rows_w_same_hist], ] <- pruned_node
    }
    
    # remove duplicated rules now!
    tree_rules <- unique(tree_rules)
  }
  
  # [3] ------ Return the pruned tree rules!
  tree_rules <- unique(tree_rules)
  print("Pruning of the tree sucessfull")
  return(data.frame(tree_rules))
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Now the Example Code, that shows the usage of the functions defined above!

# [1] Create DataSet                                                        ----
# [1-1] Blockwise df from Roman -> paste 3-Blocks to a single df:
load("./data/external/Dr_Hornung/Data/ProcessedData/BLCA.Rda")
eg_data <- cbind(clin, mirna[, 1:100], mutation[,1:100])

# [1-2] Split to test and train!
train_df <- eg_data[1:300, ]
test_df  <- eg_data[301:310,]

# [2] Fit a standard tree on train data & extract its rules                 ----
# [2-1] Seed for reproducibility, when fitting the tree
set.seed(100)
tree  <- rpart(gender ~ ., data = train_df, method = "class", 
               control = rpart.control(minsplit = 50, maxdepth = 3))

# [2-2] Extract the rules of the fitted Trees!
tree_rules <- rpart.plot::rpart.rules(tree, cover = T, nn = T, roundint = F)

# [2-3] Create a Vizualization to manually check predicitons from function!
fancyRpartPlot(tree, caption = NULL) # ugly version via: 'drawTreeNodes(tree)' 


# [3] Use tree_rules to crear prediciton! Compare with true prediciton      ----
# Implemented Version vs Own Version
print("implemented prediction: "); predict(tree, test_df[4,])
print("Own Prediciton"); get_prediciton_from_tree_rules(rules = tree_rules, 
                                                        observation = test_df[4,]) 

# [4] Prune a tree and try to get a meaningful predicition                  ----
# [4-1] Remove a predictive variable from a testobservation, so we need to prune 
#       the tree for predicting on the test_variable! --> check tree to see what 
#                                                         to remove!
test_df_w_rm_col <- test_df[,-which(colnames(test_df) == "hsa.mir.132")]
tree_rules_old   <- tree_rules

pruend_tree_rules <- prune_tree(tree_rules, test_df_w_rm_col[7,])
print("old tree_rules ---- not pruned"); tree_rules_old
print("new pruned rules, after removing 'age'"); pruend_tree_rules


# [4-2] Remove 2 predictive variables!
test_df_w_rm_col <- test_df[,-which(colnames(test_df) %in% c("hsa.mir.132", "hsa.mir.130b"))]
tree_rules_old   <- tree_rules

pruend_tree_rules <- prune_tree(tree_rules, test_df_w_rm_col[7,])
print("old tree_rules ---- not pruned"); tree_rules_old
print("new pruned rules, after removing 'hsa.mir.132' & 'hsa.mir.130b'"); pruend_tree_rules

# [5] Tree that uses factor splits & check whether we can handle it         ----
data <- iris; data$Sepal.Width <- as.factor(data$Sepal.Width); str(iris)

# [5-1] Create a tree with a factor variable for prediciting:
set.seed(1312)
tree <- rpart(Species ~ . , data = data[1:99,c("Species", "Sepal.Width")])
fancyRpartPlot(tree)

# [5-2] Check the prediction with the vizualized tree!
data[100 ,c("Species", "Sepal.Width")]
get_prediciton_from_tree_rules(rules = rpart.plot::rpart.rules(tree,cover = T, 
                                                               nn = T, roundint = F), 
                               observation = data[100 ,c("Species", "Sepal.Width")])


# [6] Get a prediction from a ensemble of trees on a complex DF             ----
# [6-1] Create a slightly more complex dataset, by enlarging the feature space!
complex_train <- cbind(eg_data, mutation[,501:1500], mirna[, 501:750])

# Compare dimensions of featurespace in old & new data!
dim(complex_train); dim(eg_data) # --> 1200 more features!

# Register a parallel backend for faster computation!
registerDoParallel(cores = detectCores() - 2)

# [6-2] Fit 100 Trees, extract the data and take the time needed!
# [6-2-1] Here we do not have subsetting of the feature space only of Observations!
now       <- Sys.time()
bagging_trees <- foreach(i = 1:100) %dopar% {
  rpart.plot::rpart.rules(
    rpart::rpart(gender ~ ., data = complex_train, method = "class",
                 minsplit = 10, maxdepth = 5, 
                 subset = c(base::sample(1:nrow(complex_train), 
                                         size = nrow(complex_train), 
                                         replace = TRUE))),
    cover = T, nn = T, roundint = F)
} 
now_done <- Sys.time()
print(now_done - now) # --> ~ 3,5min

# [6-2-2] Fit Trees w/ bagging and subsetting the feature space!
# Get the index of the column w/ reponse!
resp <- which(colnames(complex_train) == "gender")

now       <- Sys.time() 
boosted_trees <- foreach(i = 1:100) %dopar% {
  
  # Sample Features used to fit the tree
  cols_used <- c(resp, 
                 base::sample(setdiff(1:ncol(complex_train), resp), 
                              sqrt(ncol(complex_train)),  # --> how many feas per tree?!
                              replace = FALSE))
  
  # Sample Rows used as observations for fitting the tree
  rows_used <- base::sample(1:nrow(complex_train), 
                            size = nrow(complex_train), 
                            replace = TRUE)
  
  # Get the Rows and Cols we use to fit a single Decision Tree!
  curr_data <- complex_train[rows_used, cols_used]
  
  rpart.plot::rpart.rules(
    rpart::rpart(gender ~ ., data = curr_data, method = "class",
                 minsplit = 10, maxdepth = 5),
    cover = T, nn = T, roundint = F)
}
now_done <- Sys.time()
print(now_done - now) # --> ~ 15sek

# [6-3] Based on a single Observations use all single trees to get predicitons
# [6-3-1] Parallel Compution
now       <- Sys.time() 
responses_parralel <- foreach(i = 1:length(bagging_trees)) %dopar% {
  library(checkmate)
  get_prediciton_from_tree_rules(rules = bagging_trees[[i]],
                                 observation = complex_train[13,])
}
now_done <- Sys.time()
print(now_done - now) # --> ~9sek

# [6-3-2] sapply
now  <- Sys.time()
responses <- sapply(1:length(bagging_trees), 
                    FUN = function(x) get_prediciton_from_tree_rules(rules = bagging_trees[[x]],
                                                                     observation = complex_train[13,]))
now_done <- Sys.time()
print(now_done - now) # --> ~14sek

# [6-3-3] Compare whether predicitons are equal!
all.equal(responses,
          unlist(responses_parralel))
mean(responses);
#fine!

# [7] Get a prediciton on a observation, that misses a feature              ----
# [7-1] Create Observation, that misses 1 split value in at least 1 decision tree!
boosted_trees[[2]] # --> rm one of the features, that were in the tree!
new_test <- complex_train[13, -which(colnames(complex_train) %in% c('hsa.mir.6511b.2'))]

# [7-2] Prune all Trees, so we can generate a prediciton for 'new_test'
pruned_trees <- lapply(1:length(boosted_trees),
                       FUN = function(x) prune_tree(tree_rules = boosted_trees[[x]],
                                                    observation = new_test))


# [7-2-1] Print Info about how many trees had to be removed!
print(paste("From originally", as.character(length(boosted_trees)), 
            "trees, we removed",  as.character(sum(is.na(pruned_trees))))) 

if (any(is.na(pruned_trees))) pruned_trees <- pruned_trees[-which(is.na(pruned_trees))]

# [7-3] Get Predicitonson pruned trees:
responses_pruned <- sapply(1:length(pruned_trees), 
                           FUN = function(x) get_prediciton_from_tree_rules(rules = pruned_trees[[x]],
                                                                            observation = new_test))
# pruned vs. not pruned, what's the difference?!
mean(responses_pruned)
mean(responses)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


# Try 2 fold CV w/ complex data!
train <- complex_train[1:290,]
test  <- complex_train[291:310,]


# Train the RF Model:
boosted_trees <- foreach(i = 1:100) %dopar% {
  
  # Sample Features used to fit the tree
  cols_used <- c(resp, 
                 base::sample(setdiff(1:ncol(complex_train), resp), 
                              sqrt(ncol(complex_train)),  # --> how many feas per tree?!
                              replace = FALSE))
  
  # Sample Rows used as observations for fitting the tree
  rows_used <- base::sample(1:nrow(complex_train), 
                            size = nrow(complex_train), 
                            replace = TRUE)
  
  # Get the Rows and Cols we use to fit a single Decision Tree!
  curr_data <- complex_train[rows_used, cols_used]
  
  rpart.plot::rpart.rules(
    rpart::rpart(response ~ ., data = curr_data, method = "class",
                 minsplit = 10, maxdepth = 5),
    cover = T, nn = T, roundint = F)
}










# AREA TO PLAY AROUND                                                       ----
# Inspect the randomforest function
library(randomForestSRC)
a <- rfsrc



# Get all split-points & Values!  --> BOTH are not really helpful...............
rpart.subrules.table(tree)  # OB8: --> logic not 100% clear...
rattle::asRules(tree)       # better suited?!
as.party(tree)              # Not really suited, too much info to extract nicely
rattle::asRules(tree, TRUE) # looks good, but can not be saved...


# Create 2 DFs: 1st with variables and applied rules... could actually work.....
#               2nd with
# 1st DF
rule_df <- rpart.rules.table(tree) %>%
  filter(Leaf == TRUE) %>%
  group_by(Rule) %>%
  summarise(Subrules = paste(Subrule, collapse = ","))

df <- train_df %>%
  mutate(Rule = row.names(tree$frame)[tree$where]) %>%
  left_join(rule_df, by = "Rule")

# 2nd DF
df; rpart.subrules.table(tree)
# --> based on both we know which observation is split by which rules!
# --> BUT we don't know the final nodes --> can not prune...



# RULES APPROACH BEST APPROACH TO NOW ------------------------------------------
# Extract rules --> do a predicton based on these!
rules <- rpart.plot::rpart.rules(tree, cover = T, nn = T)
dim(rules)
rules # "response" = Probability for positive class! 
      # "cover"    = fraction of Observations that are in the node!
      # "nn"       = the final node?!
" 'rules'might be most important thing! 
  - read out first split criterion  --> remove tree if it has the wrong one!
                                        - quite easy to do, first split equals 
                                          first column after 'when'
  - read oout second split criterion --> combine all node trees beyond it!
  -           .
  -           .
"




# Remove blank spaces in the columns with the logical operator!
cols_w_double_blank <- grep("Oper", colnames(rules))
for (i in cols_w_double_blank) gsub(" ", "", rules[, i])

# Convert the split values to numeric!
rules[grep("Val", colnames(rules))] <- sapply(
  grep("Val", colnames(rules)), FUN = function(x) as.numeric(rules[,x])
)



# Extract the split Variables, that were used in the tree!
split_var <- unique(unlist(
  sapply(grep("Var", colnames(rules)), FUN = function(x) unique(rules[,x]))
))
if ("" %in% split_var) split_var <- split_var[split_var != ""]

# Get a TestObservation, that only keeps the features, used in the tree!
x <- test_df[6, split_var] # True Response: >>Negative<<
