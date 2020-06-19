" Script to modify 'simpleRF' implementation for our use case!

  Create a RF-Class that has the option to dynamically prune the single decision
  trees of a random forest model. This is needed for the implementation of the 
  fold-wise appropach.
  Add the method 'prune'! This method prunes a tree at splitpoints, that 
  are not avaible in the test_set we want to get predicitons for!
  We do this by enlarging the functionality of the 'simpleRF' package!
"
# [0] Load packages
require(parallel)
library(checkmate)
library(ROCR)
library(caret)

# Define the Classes & Functions we need for the Random Forrest!
# [1] Data                                                                  ----
# Define the class Data, which makes usage of dataframe a bit easier!
Data <- setRefClass("Data", 
fields = list(
  data = "data.frame", 
  ncol = "integer", 
  nrow = "integer",
  names = "character"),

methods = list(
  initialize = function(...) {
    callSuper(...)
    ncol <<- ncol(data)
    nrow <<- nrow(data)
    names <<- colnames(data)
},

  column     = function(col) {
    return(data[, col])   
},

  subset     = function(row, col) {
    return(data[row, col])
})
)

# [2] Tree                                                                  ----
Tree <- setRefClass("Tree", 
fields = list(
  mtry              = "integer", 
  min_node_size     = "integer",
  splitrule         = "character",
  unordered_factors = "character",
  data              = "Data", 
  sampleIDs         = "list",    
  oob_sampleIDs     = "integer",
  child_nodeIDs     = "list",      
  split_varIDs      = "integer",    
  split_values      = "numeric", 
  split_levels_left = "list"), 

methods = list(
  grow = function(replace) {
    "Method to grow a tree based on data and formula
      - draw obs. w/ | w/o replacment for growing and OOB
      - call splitNode then on the 1. Node!
        [recursivly done then!] 
      
    "
    # Amount of Observations we draw to grow the tree
    num_samples <- data$nrow
    if (replace) {
      num_bootstrap_samples <- num_samples
    } else {
      num_bootstrap_samples <- num_samples * 0.6321
    }
    
    # Draw 'num_bootstrap_samples' Obs. from all Observations to grow the tree 
    # with! Unused Obs. will be used as OOB-Samples!
    bootstrap_sample <- sample(num_samples, num_bootstrap_samples, replace = replace)
    oob_sampleIDs    <<- (1:num_samples)[-bootstrap_sample]
    
    # Assign bootstrap samples to root node
    sampleIDs <<- list(bootstrap_sample)
    
    # Call splitting function on root node 
    # [function is recursive -> does this until terminal nodes are reached!]
    splitNode(1)
  }, 
  
  splitNode = function(nodeID) {
    "Function to split a node and assigning the Bootstrap IDs to the 
     corresponding child nodes and save split_value/variable etc.!
    "
    
    # Sample possible split variables [+1 & -1 to exclude response!]
    possible_split_varIDs <- sample(data$ncol - 1, mtry, replace = FALSE) + 1
    
    # Split the Node w/ one of the possible_split_varIDs! 
    # Return the best splitVariable & splitValue [if not terminal node]
    split <- splitNodeInternal(nodeID, possible_split_varIDs)
    
    # If we found a split, assign IDs to child nodes etc.
    if (!is.null(split)) {
      
      # Assign found split-VAriable & Value to the tree @ current node!
      split_varIDs[[nodeID]] <<- split$varID
      split_values[[nodeID]] <<- split$value
      
      # Create child nodes fo the tree
      left_child               <- length(sampleIDs) + 1
      right_child              <- length(sampleIDs) + 2
      child_nodeIDs[[nodeID]] <<- c(left_child, right_child)
      
      # ASsign Observations from current node into child nodes! This depends on
      # the type of split [numeric, factor]
      if (length(split_levels_left[[nodeID]]) == 0) {
        
        # NUMERIC_SPLITTING --> left node all boostrap obs. >= split value 
        idx <- data$subset(sampleIDs[[nodeID]], split$varID) <= split$value
      } else {
        
        # NON_NUMERIC_SPLITTING --> left node all boostrap obs. %in% split_levels_left
        idx <- data$subset(sampleIDs[[nodeID]], split$varID) %in% split_levels_left[[nodeID]]
      }
      sampleIDs[[left_child]]  <<- sampleIDs[[nodeID]][idx]
      sampleIDs[[right_child]] <<- sampleIDs[[nodeID]][!idx]                    
      
      # Recursivly call 'splitNode' until terminal node is reached!
      splitNode(left_child)
      splitNode(right_child)
    } else {
      
      # If no new split was found the current Node is the terminal node 
      # --> estimate probs.,.. for the terminal nodes! 
      split_values[[nodeID]] <<- estimate(nodeID)
      split_varIDs[[nodeID]] <<- NA
    }
  },
  
  splitNodeInternal = function(nodeID, possible_split_varIDs) { 
    ## Empty virtual function - called by child classes
  },
  
  estimate = function(nodeID) {
    ## Empty virtual function - called by child classes
  }, 
  
  getNodePrediction = function(nodeID) {
    ## Empty virtual function - called by child classes
  },
  
  predict = function(predict_data) {
    "Function to predicton on data! predict_data must be of type Data!
     OBACHT: data needs the same features as the data we trained the model with
             else this will raise an error!
    "
    # [0] Check that, predict data has the same layout as the data we trained 
    #     our tree with!
    if (class(predict_data) != "Data") stop("predict_data must be of class 'Data'")
    if (!all(data$names %in% predict_data$names)) {
      stop("'predict_data' needs the same columns as the data the trees were trained with")
    }
    same_order <- sapply(seq(data$ncol), function(x) data$names[x] != predict_data$names[x])
    if (any(same_order)) stop("tree$data & test_set do not have the same ordering!") 
    
    # [1] Initialize amount of predicitons + list to save predicitons!
    num_samples_predict <- predict_data$nrow
    predictions         <- list()
    
    # [2] For each sample start in rootnode and walk down the tree until we reach
    #     a terminal / pruned node!
    for (i in 1:num_samples_predict) {
      nodeID <- 1
      while (TRUE) {
        
        # [2-1] Stop if the current node is a terminal node or a pruned node!
        if (nodeID > length(child_nodeIDs) || is.null(child_nodeIDs[[nodeID]]) || grepl("pruned", child_nodeIDs[[nodeID]])) {
          break
        }
        
        # [2-2] Else move observations to the child node [based on split criterion]
        if (length(split_levels_left[[nodeID]]) == 0) { 
          
          # 2-2-1 ORDERD splitting & assigning to the child node!
          value <- as.numeric(predict_data$subset(i, split_varIDs[nodeID]))
          if (value <= split_values[nodeID]) {
            nodeID <- child_nodeIDs[[nodeID]][1]
          } else {
            nodeID <- child_nodeIDs[[nodeID]][2]
          }
        } else {
          
          # 2-2-2 UNORDERED SPLITTING & assiging to the child node!
          value <- predict_data$subset(i, split_varIDs[nodeID])
          if (value %in% split_levels_left[[nodeID]]) {
            nodeID <- child_nodeIDs[[nodeID]][1]
          } else {
            nodeID <- child_nodeIDs[[nodeID]][2]
          }
        }
      }
      
      # [2-3] Add to prediction, but make sure the tree was not pruned in the  
      #       first split!
      if (nodeID == 1) {
        predictions[[i]] <- NA
      } else {
        predictions[[i]] <- getNodePrediction(nodeID)
      }
    }
    return(simplify2array(predictions))
  }, 
  
  predictOOB = function() {
    "Get the Predicitions for the samples not used to grow the tree!
     Logic pretty much, as in 'predict' method - see there for comments!
    "
    
    ## Initialize
    num_samples_predict <- length(oob_sampleIDs)
    predictions         <- list()
    
    ## For each sample start in root and drop down tree
    for (i in 1:num_samples_predict) {
      nodeID <- 1
      while (TRUE) {
        ## Break if terminal or pruned node
        if (nodeID > length(child_nodeIDs) || is.null(child_nodeIDs[[nodeID]]) || grepl("pruned", child_nodeIDs[[nodeID]])) {
          break
        }
        
        ## Move to child
        if (length(split_levels_left[[nodeID]]) == 0) {
          ## Ordered splitting
          value <- as.numeric(data$subset(oob_sampleIDs[i], split_varIDs[nodeID]))
          if (value <= split_values[nodeID]) {
            nodeID <- child_nodeIDs[[nodeID]][1]
          } else {
            nodeID <- child_nodeIDs[[nodeID]][2]
          }
        } else {
          ## Unordered splitting
          value <- data$subset(oob_sampleIDs[i], split_varIDs[nodeID])
          if (value %in% split_levels_left[[nodeID]]) {
            nodeID <- child_nodeIDs[[nodeID]][1]
          } else {
            nodeID <- child_nodeIDs[[nodeID]][2]
          }
        }
        
      }
      
      # Grt Prediciton for the terminal node! If the tree was pruned 
      # in its frst split split var. we do no prediciton, as this 
      # would only represent the distribution in the first node!
      if (nodeID == 1) {
        predictions[[i]] <- NA
      } else {
        predictions[[i]] <- getNodePrediction(nodeID)
      }
    }
    
    return(simplify2array(predictions))
  },
  
  prune = function(test_data) {
    "Prune a tree, by cutting after nodes, that use unkown variables!
     data must be of type 'Data'
    "    
    # [0] Check that, test_data has the same layout as the data we trained 
    #     our tree with!
    if (class(test_data) != "Data") stop("predict_data must be of class 'Data'")
    if (!all(data$names %in% test_data$names)) {
      stop("'predict_data' needs the same columns as the data the trees were trained with")
    }
    same_order <- sapply(seq(data$ncol), function(x) data$names[x] != test_data$names[x])
    if (any(same_order)) stop("tree$data & test_data do not have the same ordering!") 
    
    # [1] Extract the splitvariables used in tree
    used_split_var       <- split_varIDs[-which(is.na(split_varIDs))]
    used_split_var_names <- colnames(data$data)[used_split_var]
    
    # [2] Get the variables in test_data with missing values!
    missing_var       <- sapply(used_split_var_names, 
                                FUN = function(x) any(is.na(test_data$data[,x])))
    missing_var_names <- used_split_var_names[missing_var]
    
    
    # [4] PRUNE - If we have reached this point, we NEED TP PRUNE THE TREE
    # [4-1] Find all Nodes, that use a split_point, not avaible in our data!
    # 4-1-1 Get the names of all split_variables used in the tree [incl position]
    split_var_names_orderd <- colnames(data$data)[split_varIDs]
    
    # 4-1-2 Get all nodes, that use a missing variable for splitting!
    nodes_to_prune <- c()
    for (miss_var_name_curr in missing_var_names) {
      curr_          <- which(split_var_names_orderd == miss_var_name_curr)
      nodes_to_prune <- c(nodes_to_prune, curr_)
    }
    
    # [4-2] Set the 'child_nodeID' status of all the 'nodes_to_prune' to 'pruned'
    for (curr_prune in nodes_to_prune) {
      child_nodeIDs[[curr_prune]] <<- "pruned"
    }
  })
)

# [3] ClassProbability Tree                                                 ----
TreeProbability <- setRefClass("TreeProbability", 
contains = "Tree",
fields = list(),
methods = list(
  splitNodeInternal = function(nodeID, possible_split_varIDs) {
   "Check whether the current NODE is further splitable and doesn't
    break any assumptions [min_node_size, is pure already, ...]
   "
   # Node size, stop if maximum reached
   if (length(sampleIDs[[nodeID]]) <= min_node_size) {
     return(NULL)
   }
   
   # Stop if node is pure      
   unique_response <- unique(data$subset(sampleIDs[[nodeID]], 1))
   if (length(unique_response) == 1) {
     return(NULL)
   }
   
   # If none of the stopping criteria is valid, we look for the best split!
   return(findBestSplit(nodeID, possible_split_varIDs))
}, 
  
  findBestSplit = function(nodeID, possible_split_varIDs) {
   " Find the best split variable & Value for the 'NodeID' and the
     'possible_split_varIDs'
   "
   # [0] Initalize best_split - tracks best split values!
   best_split <- NULL
   best_split$decrease <- -1
   best_split$varID <- -1
   best_split$value <- -1
   
   # [1] Get the response of all Observation in 'nodeID'
   response <- data$subset(sampleIDs[[nodeID]], 1)
   
   # [2] Loop over all possible split_variables and test all possible 
   #     split_values to find the best split for the current nodeID
   for (i in 1:length(possible_split_varIDs)) {
     
     # [2-1] Get current split_variable and possible its values!
     split_varID <- possible_split_varIDs[i]
     data_values <- data$subset(sampleIDs[[nodeID]], split_varID)
     
     # [2-2] Based on the class of the current 'split_varID' find the best split
     # 2-2-1 Non-Numeric unorderd values + 'order_split' as 'unordered_factors'
     if (!is.numeric(data_values) & !is.ordered(data_values) & unordered_factors == "order_split") {
      
       # Order the factor levels, by the mean value of the repsonse 
       # [for the certain factor level]
       num.response <- as.numeric(response)
       means <- sapply(levels(data_values), function(x) {
         mean(num.response[data_values == x])
       })
       levels.ordered <- as.character(levels(data_values)[order(means)])
       
       # Get all levels not in node
       levels.missing <- setdiff(levels(data_values), levels.ordered)
       levels.ordered <- c(levels.missing, levels.ordered)
       
       # Return reordered factor
       data_values <- factor(data_values, levels = levels.ordered, ordered = TRUE)
     }
     
     # 2-2-2 Non-Numeric unorderd values + NOT order_split' as 'unordered_factors'
     if (!is.numeric(data_values) & !is.ordered(data_values)) {
       
       # Find best split w/ extra method!
       best_split = findBestSplitValuePartition(split_varID, data_values, best_split, response)
       
       # Set split levels left [all w/ these lvls go to left child node!]
       if (best_split$varID == split_varID) {
         split_levels_left[[nodeID]] <<- best_split$values_left
       }
     } else {
       
       # 2-2-3 ORDERD features
       # Find best split value orderd (numeric, orderd factors etc.)
       best_split = findBestSplitValueOrdered(split_varID, data_values, best_split, response)
       
       # Set split levels left (empty if ordered splitting)
       if (unordered_factors == "order_split") {
         if (best_split$varID == split_varID) {
           split_levels_left[[nodeID]] <<- unique(data_values[data_values <= best_split$value])
           if (is.factor(data_values)) {
             ints <- as.integer(factor(split_levels_left[[nodeID]], levels = levels(data$subset(sampleIDs[[nodeID]], split_varID))))
             if (sum(2^(ints-1)) >= 2^(max(as.numeric(data$subset(sampleIDs[[nodeID]], split_varID))) - 1)) {
               split_levels_left[[nodeID]] <<- unique(data_values[data_values > best_split$value])
             }
           }
         }
       } else {
         if (best_split$varID == split_varID) {
           split_levels_left[[nodeID]] <<- list()
         }
       }
     }
   }
   
   # [3] If found a split_variable return best split else NULL
   if (best_split$varID < 0) {
     ## Stop if no good split found
     return(NULL)
   } else {
     ## Return best split
     result <- NULL
     result$varID <- as.integer(best_split$varID)
     result$value <- best_split$value
     return(result)
   }  
  }, 
  
  findBestSplitValueOrdered = function(split_varID, data_values, best_split, response) {
   "Find the best split value for odered featuers [numeric, orderd factors,...]"
   # [0] Get all possible splitvalues of the current split_Variable Canidate!
   possible_split_values <- unique(data_values)
   
   # [1] Loop over these possible splitvalues and record the gini indes [lower is better!]
   for (j in 1:length(possible_split_values)) {
     split_value <- possible_split_values[j]
     
     # Count classes in child nodes
     idx <- data_values <= split_value
     class_counts_left <- tabulate(response[idx])
     class_counts_right <- tabulate(response[!idx])
     
     # Skip if one child empty
     if (sum(class_counts_left) == 0 | sum(class_counts_right) == 0) {
       next
     }
     
     # Calculate the GiniScore
     if (splitrule == "Gini") {
       ## Decrease of impurity
       decrease <- sum(class_counts_left^2)/sum(class_counts_left) + 
         sum(class_counts_right^2)/sum(class_counts_right)
     } else {
       stop("Unknown splitrule.")
     }
     
     # Use this as split if better than the best split before
     if (decrease > best_split$decrease) {
       best_split$value <- split_value
       best_split$varID <- split_varID
       best_split$decrease <- decrease
     }
   }
   return(best_split)
  },
  
  findBestSplitValuePartition = function(split_varID, data_values, best_split, response) {
   "Find best SplitValue for features that are unorderd (factors, characters...)"
    
   # [0] Record all possible splits for current split_var
   possible_split_values <- sort(unique(data_values))
   
   # [1] Loop over all all '2^(n-1)-1' 2-partitions and record the gini score!
   num_partitions <- 2^(length(possible_split_values) - 1) - 1
   for (j in 1:num_partitions) {
     
     # Convert number to logic vector
     left_idx    <- as.bitvect(j, length = length(possible_split_values))
     values_left <- possible_split_values[left_idx]
     
     # Count classes in child nodes
     idx <- data_values %in% values_left
     class_counts_left  <- tabulate(response[idx])
     class_counts_right <- tabulate(response[!idx])
     
     # Skip if one child empty
     if (sum(class_counts_left) == 0 | sum(class_counts_right) == 0) {
       next
     }
     
     #Calculate the GiniScore for the current split variable/value!
     if (splitrule == "Gini") {
       ## Decrease of impurity
       decrease <- sum(class_counts_left^2)/sum(class_counts_left) + 
         sum(class_counts_right^2)/sum(class_counts_right)
     } else {
       stop("Unknown splitrule.")
     }
     
     # Use this split if better than best split before
     if (decrease > best_split$decrease) {
       best_split$values_left <- values_left
       best_split$varID <- split_varID
       best_split$decrease <- decrease
     }
   }
   return(best_split)
  },
  
  estimate = function(nodeID) {     
   ## Return only NA, value is not used later
   return(NA)
  }, 
  
  getNodePrediction = function(nodeID) {
   # Return class fractions in the current nodeID
   node_samples <- sampleIDs[[nodeID]]
   table(data$subset(node_samples, 1))/length(node_samples)
  })
)
# [4] Help-Functions from simpleRF-Package                                  ----
#     Not all of these are commented in detail!
which.max.random <- function(x) {
  if (all(is.na(x))) {
    return(NA)
  }
  which(rank(x, ties.method = "random", na.last = FALSE) == length(x))
}

which.min.random <- function(x) {
  if (all(is.na(x))) {
    return(NA)
  }
  which(rank(x, ties.method = "random", na.last = TRUE) == 1)
}

cor.order <- function(y, x) {
  ## Create contingency table of the nominal outcome with the nominal covariate
  tab <- table(droplevels(y), droplevels(x))
  
  ## Compute correlation matrix of the contingency table (correlation of the covariate levels w.r.t outcome)
  cr <- suppressWarnings(cor(tab))
  cr[is.na(cr)] <- 0                      
  diag(cr) <- NA
  
  ## Start with a random level and select as next level the level with highest correlation to the current level (excluding already selected levels)
  num_levels <- nlevels(droplevels(x))
  next_level <- sample(num_levels, 1)
  res <- c(next_level, rep(NA, num_levels - 1))
  for (i in 2:num_levels) {
    cr[, next_level] <- NA
    next_level <- which.max.random(cr[next_level, ])
    res[i] <- next_level
  }
  
  ## Return ordered factor levels
  as.character(levels(droplevels(x))[res])
}

pc.order <- function(y, x) {
  " Function to transform a non orderd factor/ character variable, into a orderd
    factor, depending on the distribution of the response!"
  
  # [0] If 'x' is character, transform it to factor
  if (class(x) == "character") x <- as.factor(x)
    
  # [1] If 'x' has less than 2 unique classes, transfer it to character!
  if (nlevels(droplevels(x)) < 2) {
    return(as.character(levels(droplevels(x))))
  }
  
  # [2] If 'x' has more than 2 classes, create contingency table of the nominal 
  #     outcome with the nominal covariate
  N <- table(droplevels(x), droplevels(y))
  
  # [3] PCA of weighted covariance matrix of class probabilites
  P     <- N/rowSums(N)                    # probs. for response given a factor 
                                           # lvl & vic versa
  S     <- cov.wt(P, wt = rowSums(N))$cov  # estimates of weighted cov-matrix
  pc1   <- prcomp(S, rank. = 1)$rotation   # 1. Principal component of the 
                                           # weighted cov-matrix
  score <- P %*% pc1                       # Main components for the factor levels!
  
  # [4] Return factor levels and order them by the score each factor levels has 
  #     scored!
  as.character(levels(droplevels(x))[order(score)])
}

reorder.factor.columns <- function(data) {
  
  
  # [1] Get the features, that are of class 'char' or  are unoreded factors
  #     These need to be recoded!
  character.idx <- sapply(data[, -1], is.character)
  ordered.idx   <- sapply(data[, -1], is.ordered)
  factor.idx    <- sapply(data[, -1], is.factor)
  recode.idx    <- character.idx | (factor.idx & !ordered.idx)
  
  # [2] Extract the response and transfer it to a numeric variable!
  response <- data[, 1]
  if (is.factor(response)) {
    num.response <- as.numeric(response)
  } else if ("Surv" %in% class(response)) { # ---> NEEDS TO BE CHECKED!
    num.response <- coin::logrank_trafo(response, ties.method = "Hothorn-Lausen")
  } else {
    num.response <- response
  }
  
  # [3] Recode the cols of class 'char' / unoreded factors - depending on type 
  #     of response!
  data[, -1][, recode.idx] <- lapply(data[, -1][, recode.idx, drop = FALSE], function(x) {
    
    if (is.factor(response) & nlevels(response) > 2) {
      
      # [3-1] Recoding for multiclass repsonse!
      levels.ordered <- pc.order(y = response, x = x)
    } else {
      
      # [3-2] Recoding for all other responses!
      means <- sapply(levels(x), function(y) {
        mean(num.response[x == y])
      })
      levels.ordered <- as.character(levels(x)[order(means)])
    }
    
    # [3-3] Return the reordered factors!
    factor(x, levels = levels.ordered, ordered = TRUE)
  })
  
  # [4] Return the data
  data
}

as.bitvect <- function(x, length = 32) {
  i <- 1
  string <- numeric(length)
  while(x > 0) {
    string[i] <- x %% 2
    x <- x %/% 2
    i <- i + 1 
  }
  as.logical(string)
}

# [5] Custom Help-Functions for Pruning, testdata-processing, ...           ----
process_test_data <- function(tree, test_data) {
  "Function to process the TestData in the same way, as our TrainData, we used
   to grow the tree!
   Important, that factor levels are orderd in the same way, as in train data +
   that we have the same structure as in the data used to grow the tree!
   [Exact same ordering, same feature variables & same level for the non numeric
    feature variables!]
    
    Args:
      - tree (Tree)            : Tree Object, that is already grown and we want 
                                 to use to make predicitons!
      - test_data (data.frame) : data we want to create predicitons for!
                                 In 'simpleRF' Implementation we have to take
                                 care, that the data we want to do predicitons
                                 for has exactly the same layout as our train 
                                 data used to grow the tree + the factor levels
                                 need the exact same ORDED levels!
    Return:
      - data (Data-Object), with the same features in the same layout as the
        tree$data$data dataframe! 
        Aswell same ordering of the factor levels!
  "
  # [0] Check Input
  if (!(grepl("Tree", x = class(tree)))) stop("tree must be of (sub)class 'Tree'")
  assert_data_frame(test_data, min.rows = 1)
  
  # [1] Make Sure, TestData and Train Data have exactly the same set of features!
  # [1-1] Get colnames of data, used to grow the tree + of testdata!
  train_cols_order <- colnames(tree$data$data)
  test_cols        <- colnames(test_data)
  
  # [1-2] Check, whether features are missing in testdata and add them w/ NA Value!
  missing_variables <- train_cols_order[which(!(train_cols_order %in% test_cols))]
  print(paste0("Amount missing Variables in TestData: ", length(missing_variables)))
  
  # [1-2-1] If there are missing variables, add them to test_data, but only w/ NA Values
  if (length(missing_variables) >= 1) {
    for (miss_var in missing_variables) {
      test_data[, miss_var] <-  NA
    }
  }
  
  # [2] Order the Testdata in the same way as tree$data$data
  test_data <- test_data[, train_cols_order]
  
  # [3] For unordered factors & character features, we check whether they have 
  #     the same levels & if this is fine, we set the ordering of the levels as 
  #     for the tree$data$data
  
  # [3-1] Get all non numeric feas [factors & characters columns!] (exclude response) 
  non_numeric_feas    <- sapply(1:tree$data$ncol, 
                                FUN = function(x) !is.numeric(tree$data$data[,x]))
  non_numeric_feas[1] <- FALSE    # exclude response!
  non_numeric_feas_names <- names(tree$data$data[,non_numeric_feas])
  
  # [3-2] Set the levels in test_data as they are in the data, used to train a tree
  test_data[,non_numeric_feas_names] <- lapply(non_numeric_feas_names, 
                                               FUN = function(x) {
                                                 factor(test_data[,x], ordered = TRUE, 
                                                        levels = levels(tree$data$data[,x]))
                                               })
  
  # [4] Return test_data as Data object!
  test_data = Data$new(data = test_data)
  return(test_data)
}

simpleRF <- function(formula, data, num_trees = 50, mtry = NULL, 
                     min_node_size = NULL, replace = TRUE, splitrule = NULL, 
                     unordered_factors = "ignore") {
  " Create a RandomForest, as a list of single TREEs!
  
    Function to check the passed arguments for validity & then preprocess the 
    arguments and set them to reasonable values!
    With these, we create 'num_trees' objects of the class 'Tree' 
    [type of tree we create [Classif, Survival,...] depends on the type of the 
    response column!
    
    Args:
      - formula (formula)        : Formula describing the relation between 
                                   response & features
      - data (data.frame)        : Data used to build the tree! Must contain all
                                   names defined in 'formula'
      - num_trees (int)          : Amount of trees [int > 1] to create!
      - mtry (int)               : Amount of possible split_candidates, we draw 
                                   at every split point! If NULL, we use 
                                   ceiling(sqrt(p))
      - min_node_size (int)      : Amount of Observations, that, at least, need 
                                   to be in a terminal node!
      - replace (boolean)        : When drawing the observations from data, we 
                                   use to build the tree, shall we draw w/ or
                                   w/o replacement
      - splitrule (char)         : Splitrule to use! - depends on type of 
                                   response - either NULL or it must be in
                                   'Gini', 'Variance', 'Logrank'
      - unordered_factors (char) : How to handle non numeric features, mut be in
                                   'ignore', 'order_once', 'order_split', 'partition'
    Return:
      - list filled with 'num_tree' objects of the class 'Tree'! All with the
        internal arguments, as they have been processed by this function!
  "
  # [0] Check Arguments - that need to be passed & can not set w/ 'NULL'
  assert_formula(formula)
  assert_data_frame(data, any.missing = FALSE)
  assert_int(num_trees, lower = 1)
  assert_logical(replace, len = 1)
  assert_character(unordered_factors, len = 1)
  
  # [0-1] All names from formula in data?
  formula_vars <- all.vars(formula)
  if (!(formula_vars[-which(formula_vars == ".")] %in% colnames(data))) {
    stop("not all variables, specified in 'formula' are feature names in data")
  }
  
  # [1] Based on the formula, order and prune our data + set type of response!
  #     [plain classification was excluded, as we can use probabilites w/o
  #     loosing any type of information!]
  model.data <- model.frame(formula, data)
  
  if (class(model.data[, 1]) == "factor") {
    treetype <- "Probability" 
  } else if (class(model.data[, 1]) == "numeric") {
    treetype <- "Regression"
  } else if (class(model.data[, 1]) == "Surv") {
    treetype <- "Survival"
  } else {
    stop("Unkown response type.")
  }
  
  # [1-1] Only support ClassifTreesuntil now --> stop, if the treetpye is 
  #       different from "Probability"
  if (treetype != "Probability") {
    stop(paste("SUPPORT FOR CLASSIFICATION ONLY! Your response yielded a tree of tpye: ", treetype))
  }
  
  # [2] Check Arguments, that can be passed w/ 'NULL'
  # [2-1] mtry - set to standard values from literature if NULL & < ncol(data)
  if (is.null(mtry)) {
    mtry <- as.integer(ceiling(sqrt(ncol(model.data) - 1))) # -1 to exclude response
  } else if (mtry > ncol(model.data) - 1) {
    stop("'mtry' cannot be larger than number of feature variables!")
  }
  
  # [2-2] min_node_size - choose depening on class of response 
  if (is.null(min_node_size)) {
    if (treetype == "Probability") {
      min_node_size <- as.integer(10)
    } else if (treetype == "Regression") {
      min_node_size <- as.integer(5)
    } else if (treetype == "Survival") {
      min_node_size <- as.integer(3)
    }
  } else {
    assert_int(min_node_size, lower = 1)
  }
  
  # [2-3] Splitrule - choose depening on class of response 
  if (is.null(splitrule)) {
    if (treetype == "Probability") {
      splitrule <- "Gini"
    } else if (treetype == "Regression") {
      splitrule <- "Variance"
    } else if (treetype == "Survival") {
      splitrule <- "Logrank"
    }
  } else {
    if (!(splitrule %in% c("Gini", "Variance", "Logrank"))) {
      stop("splitrule must be in 'Gini', 'Variance', 'Logrank'")
    }
  }
  
  # [2-4] Unordered factors - must have valid argument
  if (!(unordered_factors %in% c("ignore", "order_once", "order_split", "partition"))) {
    stop("Unknown value for unordered_factors.")
  }
  
  # [3] Reorder the character/ unorded factor features in 'data' / 'model_data'
  # [3-1] Create a list, where we save the ordering of the features etc.
  covariate_levels <- list()
  
  # [3-2] Reorder the levels [see 'reorder.factor.columns' for details]
  if (unordered_factors == "order_once") {
    model.data       <- reorder.factor.columns(model.data)
    covariate_levels <- lapply(model.data[, -1], levels)
  } 
  
  # [3-3] Keep the order as it is, BUT mark unorderd factor/ character features 
  #       as ordered and treat them as numerics, when considering splits!
  #       [as it is orderd we have classes > & < other classes!]
  if (unordered_factors == "ignore") {
    
    # [3-3-1] Get the Columns, that are of type character OR of type factor, but
    #         not orderd yet!
    character.idx <- sapply(model.data[, -1], is.character)
    ordered.idx   <- sapply(model.data[, -1], is.ordered)
    factor.idx    <- sapply(model.data[, -1], is.factor)
    recode.idx    <- character.idx | (factor.idx & !ordered.idx)
    
    # [3-3-2] Mark the columns [char or not orderd factors] as ordered
    #         Ordering & factors stay the, same jsut as orderd 
    #        (some level > other levels - rather arbitrary!)
    if (sum(recode.idx) == 0) {
      print("No factors/ characters to recode!")
    } else if (sum(recode.idx) == 1) {
      model.data[, -1][, recode.idx] <- as.ordered(model.data[, -1][, recode.idx])
    } else {
      model.data[, -1][, recode.idx] <- lapply(model.data[, -1][, recode.idx], as.ordered)
    }
    covariate_levels <- lapply(model.data[, -1], levels)
  }
  
  # [5] Create forest object, based on the passed/ corrected arguments!
  if (treetype == "Probability") {
    
    # Create 'num_trees' trees!
    trees    <- replicate(num_trees, TreeProbability$new())
    
    # Pass the necesary Information to the trees!
    temp <- lapply(trees, function(x) {
      x$mtry              <- mtry
      x$min_node_size     <- as.integer(min_node_size)
      x$splitrule         <- splitrule
      x$unordered_factors <- unordered_factors
      x$data              <- Data$new(data = model.data)
    })
    
    return(trees)
  } else {
    stop("Unkown tree type.")
  }
  
  # This is not done yet! Works for Classificaion only!
  # --> Can be adjusted! Most of the work has already been done!
    # else if (treetype == "Regression") { 
  #   # Create 'num_trees' trees!
  #   trees    <- replicate(num_trees, TreeRegression$new())
  #   
  #   # Pass the necesary Information to the trees!
  #   temp <- lapply(trees, function(x) {
  #     x$mtry              <- mtry
  #     x$min_node_size     <- min_node_size
  #     x$splitrule         <- splitrule
  #     x$unordered_factors <- unordered_factors
  #     x$data              <- Data$new(data = model.data)
  #   })
  # } else if (treetype == "Survival") { 
  #   idx.death <- model.data[, 1][, 2] == 1
  #   timepoints <- sort(unique(model.data[idx.death, 1][, 1]))
  #   forest <- ForestSurvival$new(num_trees = as.integer(num_trees), mtry = as.integer(mtry), 
  #                                min_node_size = as.integer(min_node_size), 
  #                                replace = replace, splitrule = splitrule,
  #                                data = Data$new(data = model.data), 
  #                                formula = formula, unordered_factors = unordered_factors, 
  #                                covariate_levels = covariate_levels,
  #                                timepoints = timepoints)
  # } else {
  #   stop("Unkown tree type.")
  # }
}

get_pruned_prediction <- function(trees, test_set) {
  " Recieve the Predicitive Performance of a RandomForest [list of Trees] on the 
    'test_set'. This 'test_set' can also miss variables, the original trees had 
    in their train_set! 
    Important: 'trees' is a list filled with objects of the class 'Tree' & 
               'test_set' is an object of class 'Data'!
                [use 'process_test_data()' to process the 'test_set']
               
    Args:
      - trees    (list) : list filled with Object of the Class 'Tree'
                          These trees were all grown on the same dataset w/ the
                          same trainset + same parameters!
      - test_set (Data) : Object of type 'Data' 
                          --> processed via 'process_test_data()'
    Return:
      - List w/ entrances for the aggregated predictions from all trees for all 
        observations!
         --> [[1]] For each observation the aggregated probabilities of all 
                   trees for all possbile classes!
         --> [[2]] For each observation the aggregated class prediciton of all
                   trees (~> class w/ highest probability!)
  "
  # [0] InputChecking ----------------------------------------------------------
  # [0-1] Check 'trees' argument - list filled w/ 'Tree' Objects!
  assert_list(trees, min.len = 2, any.missing = F)
  tree_obj <- sapply(trees, FUN = function(x) grepl("Tree", class(x)))
  if (any(!tree_obj)) stop("tree is not only filled with objects of class 'Tree'")
  
  # [0-2] Check 'test_set' argument
  # [0-2-1] Check for right class + at least 1 observation in the set
  if (class(test_set) != "Data") stop("'test_set' must be of type 'Data'")
  if (test_set$nrow < 1) stop("'test_set' needs at least one observation!")
  
  # [0-2-2] 'test_set' needs same colnames + order of colnams as the data the
  #         trees were trained with - else predictions are meaningless!
  if (!all(trees[[1]]$data$names %in% test_set$names)) {
    stop("'test_set' needs the same columns as the data the trees were trained with")
  }
  
  same_order <- sapply(seq(trees[[1]]$data$ncol), function(x) trees[[1]]$data$names[x] != test_set$names[x])
  if (any(same_order)) stop("tree$data & test_set do not have the same ordering!") 
  
  # [1] Pruning of the trees ---------------------------------------------------
  # [1-1] Prune the trees [child nodeIDs to 'pruned', if split_var not in test]
  tmp <- sapply(trees, FUN = function(x) x$prune(test_set))
  
  # [1-2] Remove the trees, of which the first split was pruned!
  trees_to_rm <- sapply(trees, FUN = function(x) x$child_nodeIDs[1] == 'pruned')
  if (sum(trees_to_rm) >= 1) trees <- trees[-c(which(trees_to_rm))]
  
  # [1-3] Print Info, about pruning process:
  print(paste("From originally", length(trees_to_rm), "Trees,",
              sum(trees_to_rm), "were removed in pruning"))
  print(paste("Forrest consiting of", length(trees), "Trees now!"))
  
  # [1-4] If all trees, were removed, we can not create any predicition
  if (length(trees) == 0) {
    print("Can not used for prediciton, as all trees were removed in the pruning")
    return(list("Probs" = rep(NA, test_set$nrow),
                "Class" = rep(NA, test_set$nrow)))
  }
  
  # [2] Get Prediciton of the trees --------------------------------------------
  # [2-1] Get the predicitons of the trees
  #       [each tree a prediciton for each observation in 'test_set']
  predictions <- lapply(trees, function(x) {
    res <- x$predict(test_set)
    res
  })
  
  # [2-2] Form an aggregated prediciton over all trees for all observations!
  # [2-2-1] Aggregated ClassProbabilities
  aggregated_predictions_prob <- lapply(seq(test_set$nrow), 
                                          function(x) {
                                            'For each observation aggregate the prediciton from all single trees!
                                             Trees, that can not be used for predicitons are not considered
                                            '
                                            
                                            # Extract predciton for observation 'x' from all tree-predcitions!
                                            preds_obs_x <- sapply(predictions, FUN = function(i) {
                                              as.data.frame(i)[,x]
                                            })
                                            
                                            # Create probabilites from all the tree predicitons!
                                            preds = apply(preds_obs_x, MARGIN = 1, 
                                                          FUN = function(x) sum(x)/length(predictions))
                                            preds
                                          })
  
  # [2-2-2] Get Class Predictions [= Class w/ highest probability]
  poss_classes <- as.character(levels(trees[[1]]$data$subset( ,1)))
  
  aggregated_predictions_class <- unlist(
    lapply(aggregated_predictions_prob, FUN = function(x) {
      "From the class probabilites extract the class prediction itself"
      class_highest_prob <- poss_classes[which(x == max(x))]
      class_highest_prob
    })
  )
  
  # [3] Return the aggregated predicitons & probabilites -----------------------
  results = list("Probs" = aggregated_predictions_prob,
                 "Class" = aggregated_predictions_class)
  return(results)
}