"
Script to visualize the different Results from the CrossValidations of the 
different random forest adaptions on the TCGA data. 
"
# [0] Load Packages and define functions
library(ggplot2)
require(gridExtra)
library(reshape2)
library(checkmate)

# Used for all Results except from 'FoldWise' & 'SingleBlock' approach
extract_metrics     <- function(x, metric, train_sit) {
  "
  Function to convert a list 'x' - w/ the 2 entrances 'res_all' & 'settings' - 
  to a DF, that we can use for plotting!
  'x' is the result of a CV of one trainsetting & all has results for all of 
  the 20 testsettings! From this list we extract the metric for all of the 
  different trainsettings!
  
   Args:
    x (list)           : list filled with the metrics for every single 
                         Itteration in CV + the settings it was created from!
                         Needs the names: 'res_all' & 'settings'
    metric (str)       : the metric we shall extraxt from the result list!
                         Needs to be in:
                         ['Accuracy', 'Sensitifity', 'Specificity', 'Precision', 
                          'Recall', 'F1', 'Balance_Acc', 'AUC', 'MCC']
    train_sit (int)    : Which TrainingSituation do the results come from!
                         [four different trainsettings in total s. MS Thesis!]

   Return:
    DF, suited for plotting, with all metric for all testsettings in 'x'
  "
  # [0] Check Inputs 
  # 0-1 'x' must be a list of length 2 w/ names 'res_all' & 'settings'
  assert_list(x, len = 2)
  if (!('res_all' %in% names(x)) | !('settings' %in% names(x))) {
    stop("'x' needs 2 entrances called 'res_all' & 'settings'")
  }
  
  # 0-2 'metric' must be in 'x'
  assert_string(metric)
  
  # extract all used  metrics and check whether 'metric' exist!
  met_meas <- unique(
    unlist(sapply(names(x$res_all), function(j) {
    names(x$res_all[[j]][[1]])
  })))
  
  if (!(metric %in% met_meas)) {
    stop("'metric' is not within these!")
  }
  
  # 0-3 'train_sit' must be integer [1; 4]
  assert_int(train_sit, lower = 1, upper = 4)
  
  # [1] Unpack the lists and put them into a single regular DF
  # 1-1 Seperate Results and Settings
  x_res      <- x$res_all
  x_settings <- x$settings
  
  # 1-2 Extract the metric from the wanted TrainSetting for all 20 TestSettings!
  #     ["full", "miss_1A", ...]
  res_curr_train_set <- sapply(names(x_res), function(x) c(x_res[[x]][[1]][metric],
                                                           x_res[[x]][[2]][metric],
                                                           x_res[[x]][[3]][metric],
                                                           x_res[[x]][[4]][metric],
                                                           x_res[[x]][[5]][metric]))
  
  # 1-3 Convert the results to a usable DF
  res_curr_train_set <- as.data.frame(res_curr_train_set)
  
  DF_final <- data.frame("Trainsituation" = numeric(),
                         "Testsituation"  = character(),
                         "Fold"           = numeric(),
                         "Metric"         = character())
  
  # 1-3-1 Loop over each column in 'res_curr_train_set' (each is 1 testsetting)
  #       & fill 'DF_final' with the additional information
  for (j in seq_len(ncol(res_curr_train_set))) {
    
    # 1-3-1-1 Extract TrainSetting, folds, test_situation and the metrics
    #         for the given column!
    train_sit_ <- rep(train_sit, times = nrow(res_curr_train_set))
    folds      <- seq_len(nrow(res_curr_train_set))
    test_sit   <- rep(colnames(res_curr_train_set)[j], times = nrow(res_curr_train_set))
    metric_c   <- unlist(res_curr_train_set[,j])
    
    # If metric_c contains any NA replace these by '0' or '-1' if the metric is MMC
    if (any(is.na(metric_c))) {
      to_replace <- which(is.na(metric_c))
      if (metric == "MCC") {
        replace_value <- -1
      } else {
        replace_value <- 0
      }
      
      metric_c[to_replace] <- replace_value
      metric_c             <- as.numeric(metric_c)
    }
    
    # 1-3-1-2 Bind it to the DF with all Results
    df_current <- data.frame("Trainsituation" = train_sit_,
                             "Testsituation"  = test_sit,
                             "Fold"           = folds,
                             "Metric"         = metric_c)
    DF_final <- rbind(DF_final, df_current)
  }
  
  # [2] Add the settings that have been used
  # 2-1 Add used DF and used Seed
  path_split       <- unlist(strsplit(x_settings$data_path, split = "/"))
  DF_final$DF      <- path_split[length(path_split)]
  DF_final$DF_seed <- path_split[length(path_split) - 1]
  
  # 2-2 Add the infos from "x_settings"
  DF_final$performance_metric <- metric
  DF_final$reponse            <- x_settings$response
  DF_final$model_seed         <- x_settings$seed
  DF_final$num_trees          <- x_settings$num_trees
  DF_final$mtry               <- x_settings$mtry
  DF_final$min_node_size      <- x_settings$min_node_size
  DF_final$weighted_ens       <- x_settings$weighted
  
  if (is.null(x_settings[["weighted"]]) | !x_settings[["weighted"]]) {
    DF_final$weight_metric <- NA
  } else {
    DF_final$weight_metric <- x_settings[["weighted"]]
  }
  
  # [3] Return 'DF_final'
  return(DF_final)
}
extract_avg_metrics <- function(x, metric, train_sit) {
  "
  Function to convert a list 'x' - w/ the 2 entrances 'res_all' & 'settings' - 
  to a DF, that we can use for plotting!
  'x' is the result of a CV of one trainsetting & all has results for all of 
  the 20 testsettings! From this list we extract the metric for all of the 
  different trainsettings!
    --> Only grab the average values for the metric of all folds!
  
   Args:
    x (list)           : list filled with the metrics for every single 
                         Itteration in CV + the settings it was created from!
                         Needs the names:
                         ['res_all' & 'settings']
    metric (str)       : the metric we shall extraxt from the result list!
                         Needs to be in:
                         ['Accuracy', 'Sensitifity', 'Specificity', 'Precision', 
                          'Recall', 'F1', 'Balance_Acc', 'AUC', 'MCC']
    train_sit (int)    : Which TrainingSituation do the results come from!
                         [four different trainsettings in total s. MS Thesis!]

   Return:
    DF, suited for plotting, with all metric for all testsettings in 'x'
  "
  # [0] Check Inputs 
  # 0-1 'x' must be a list of length 2 w/ names 'res_all' & 'settings'
  assert_list(x, len = 2)
  if (!('res_all' %in% names(x)) | !('settings' %in% names(x))) {
    stop("'x' needs 2 entrances called 'res_all' & 'settings'")
  }
  
  # 0-2 'metric' must be in 'x'
  assert_string(metric)
  
  # extract all used  metrics and check whether 'metric' exist!
  met_meas <- unique(
    unlist(sapply(names(x$res_all), function(j) {
      if (length(x$res_all[[j]]) != 0) {
        names(x$res_all[[j]][[1]])
      }
    })))
  
  if (!(metric %in% met_meas)) {
    stop("'metric' is not within these!")
  }
  
  # 0-3 'train_sit' must be integer [1; 4]
  assert_int(train_sit, lower = 1, upper = 4)
  
  # [1] Unpack the lists and put them into a single regular DF
  # 1-1 Seperate Results and Settings
  x_res      <- x$res_all
  x_settings <- x$settings
  
  # 1-2 Extract the metric from the wanted TrainSetting for all 20 TestSettings!
  #     ["full", "miss_1A", ...]
  res_curr_train_set <- sapply(names(x_res), function(x) {
    if (length(x_res[[x]]) == 0) {
      c(rep(NA, times = 5))
    }
    else {
      c(x_res[[x]][[1]][metric],
        x_res[[x]][[2]][metric],
        x_res[[x]][[3]][metric],
        x_res[[x]][[4]][metric],
        x_res[[x]][[5]][metric])
    }
  })
  
  # 1-3 Convert the results to a usable DF
  res_curr_train_set <- as.data.frame(res_curr_train_set)
  
  # 1-3-1 Average the values, so for each DF we only recieve the averge metric of
  #       all folds!
  res_curr_train_set_avg <- sapply(1:ncol(res_curr_train_set), 
                                   function(x) mean(as.numeric(unlist(res_curr_train_set[,x])), na.rm = T))
  
  # 1-3-2 Add the names of the testsituations
  names(res_curr_train_set_avg) <- colnames(res_curr_train_set)
  
  # 1-4 Create a DF holding all results and Metainfos!
  DF_final <- data.frame("Trainsituation" = numeric(),
                         "Testsituation"  = character(),
                         "Metric"         = character())
  
  # 1-4-1 Loop over each column in 'res_curr_train_set' (each is 1 testsetting)
  #       & fill 'DF_final'
  for (j in seq_len(length(res_curr_train_set_avg))) {
    
    # 1-4-1-1 Extract TrainSetting, folds, test_situation and the metrics
    #         for the given column!
    train_sit_ <- train_sit
    test_sit   <- names(res_curr_train_set_avg)[j]
    metric_c   <- res_curr_train_set_avg[j]
    
    # 1-4-1-2 Bind it to the DF with all Results
    df_current <- data.frame("Trainsituation" = train_sit_,
                             "Testsituation"  = test_sit,
                             "Metric"         = metric_c)
    DF_final <- rbind(DF_final, df_current)
  }
  
  # [2] Add the settings that have been used
  # 2-1 Add used DF and used Seed
  path_split       <- unlist(strsplit(x_settings$data_path, split = "/"))
  DF_final$DF      <- path_split[length(path_split)]
  DF_final$DF_seed <- path_split[length(path_split) - 1]
  
  # 2-2 Add the infos from "x_settings"
  DF_final$performance_metric <- metric
  DF_final$reponse            <- x_settings$response
  DF_final$model_seed         <- x_settings$seed
  DF_final$num_trees          <- x_settings$num_trees
  DF_final$mtry               <- x_settings$mtry
  DF_final$min_node_size      <- x_settings$min_node_size
  DF_final$weighted_ens       <- x_settings$weighted
  
  if ("weighted" %in% names(x_settings)) {
    if (is.null(x_settings[["weighted"]]) | !x_settings[["weighted"]]) {
      DF_final$weight_metric <- NA
    } else {
      DF_final$weight_metric <- x_settings[["weighted"]]
    }
  } else {
    DF_final$weight_metric <- NA
  }
  
  # [3] Return 'DF_final'
  return(DF_final)
}

# Used to extract Results from 'FoldWise_Approach'
# - these have a slightly different structure due to performance realted changes
extract_metrics_FW     <- function(x, metric, train_sit) {
  "
  Function to convert a list 'x' - w/ the 2 entrances 'res_all' & 'settings' - 
  to a DF, that we can use for plotting!
  'x' is the result of a CV of one trainsetting & all has results for all of 
  the 20 testsettings! From this list we extract the metric for all of the 
  different trainsettings!
    --> Only grab the average values for the metric of all folds!
  
   Args:
    x (list)           : list filled with the metrics for every single 
                         Itteration in CV + the settings it was created from!
                         Needs the names:
                         ['res_all' & 'settings']
    metric (str)       : the metric we shall extraxt from the result list!
                         Needs to be in:
                         ['Accuracy', 'Sensitifity', 'Specificity', 'Precision', 
                          'Recall', 'F1', 'Balance_Acc', 'AUC', 'MCC']
    train_sit (int)    : Which TrainingSituation do the results come from!
                         [four different trainsettings in total s. MS Thesis!]

   Return:
    DF, suited for plotting, with all metric for all testsettings in 'x'
  "
  # [0] Check Inputs  ---------------------------------------------------------- 
  # 0-1 'x' must be a list of length 2 w/ names 'res_all' & 'settings'
  assert_list(x, len = 2)
  if (!('res_all' %in% names(x)) | !('settings' %in% names(x))) {
    stop("'x' needs 2 entrances called 'res_all' & 'settings'")
  }
  
  # 0-2 'metric' must be a string, also existing in 'x'
  assert_string(metric)
  
  #      Get all exisiting metrics 'x' contains 
  #      - check for 'full' results, as aslways has results!
  exis_meas_no_weight <- names(x$res_all$full[[1]]$no_weighting)
  exis_meas_f1_weight <- names(x$res_all$full[[1]]$f1_weighting)
  exis_meas_acc_weight <- names(x$res_all$full[[1]]$acc_weighting)
  
  metric_exisits <- metric %in% exis_meas_no_weight && 
                    metric %in% exis_meas_f1_weight && 
                    metric %in% exis_meas_acc_weight
  
  if (!metric_exisits) stop("metric doesn't exist in 'x'")
  
  # 0-3 'train_sit' must be integer [1; 4]
  assert_int(train_sit, lower = 1, upper = 4)
  
  # [1] Extract the metric for all testsettings  -------------------------------
  # 1-1 Get the settings the model has been trained with
  x_settings <- x$settings
  
  # 1-2 Extract the metric for all 20 TestSettings ["full", "miss_1A", ...]!
  #     --> NA weights stand for trees that were pruned at the first splitvar.
  #         so that the trees could not be used for a predicition
  # 1-2-1 NO_WEIGHTING
  #       Loop over all possible TestSets
  res_no_weight <- sapply(names(x$res_all), function(setting) {
    #     For each fold check whether it contains only a single character 
    #     [-> only happens when the tree had to be pruned] - if not extract the 
    #     metric for the current 'Setting' and the current 'fold'
    sapply(1:5, function(fold) {
      if (is.character(x$res_all[[setting]][[fold]])) {
        NA
      }
      else {
        x$res_all[[setting]][[fold]][["no_weighting"]][metric]
      }
    })
  })
  
  # 1-2-2 ACC_WEIGHTING
  #       Loop over all possible TestSets
  res_acc_weight <- sapply(names(x$res_all), function(setting) {
    #     For each fold check whether it contains only a single character 
    #     [-> only happens when the tree had to be pruned] - if not extract the 
    #     metric for the current 'Setting' and the current 'fold'
    sapply(1:5, function(fold) {
      if (is.character(x$res_all[[setting]][[fold]])) {
        NA
      }
      else {
        x$res_all[[setting]][[fold]][["acc_weighting"]][metric]
      }
    })
  })
  
  # 1-2-3 F1_WEIGHTING
  #       Loop over all possible TestSets
  res_f1_weight <- sapply(names(x$res_all), function(setting) {
      #     For each fold check whether it contains only a single character 
      #     [-> only happens when the tree had to be pruned] - if not extract the 
      #     metric for the current 'Setting' and the current 'fold'
      sapply(1:5, function(fold) {
        if (is.character(x$res_all[[setting]][[fold]])) {
          NA
        }
        else {
          x$res_all[[setting]][[fold]][["f1_weighting"]][metric]
        }
      })
  })

  # 1-3 Convert the results to a usable DF
  DF_final <- data.frame("Trainsituation" = numeric(),
                         "Testsituation"  = character(),
                         "Fold"           = numeric(),
                         "Metric"         = character())
  
  # 1-3-1 Loop over each of the different weighted results
  for (res_curr_train_set in list(as.data.frame(res_no_weight), 
                                  as.data.frame(res_f1_weight),
                                  as.data.frame(res_acc_weight))) {
    
    # 1-3-2 Loop over each column in 'res_curr_train_set' (each is 1 testsetting)
    #       & fill 'DF_final'
    for (j in seq_len(ncol(res_curr_train_set))) {
      
      # 1-3-2-1 Extract TrainSetting, folds, test_situation and the metrics
      #         for the given column!
      train_sit_ <- rep(train_sit, times = nrow(res_curr_train_set))
      folds      <- seq_len(nrow(res_curr_train_set))
      test_sit   <- rep(colnames(res_curr_train_set)[j], times = nrow(res_curr_train_set))
      metric_c   <- unlist(res_curr_train_set[,j])
      
      # If metric_c contains any NA replace these by '0' or '-1' if the metric is MMC
      if (any(is.na(metric_c))) {
        to_replace <- which(is.na(metric_c))
        if (metric == "MCC") {
          replace_value <- -1
        } else {
          replace_value <- 0
        }
        
        metric_c[to_replace] <- replace_value
        metric_c             <- as.numeric(metric_c)
      }
      
      # 1-3-2-2 Bind it to the DF with all Results
      df_current <- data.frame("Trainsituation" = train_sit_,
                               "Testsituation"  = test_sit,
                               "Fold"           = folds,
                               "Metric"         = metric_c)
      DF_final <- rbind(DF_final, df_current)
    }
  }
  
  # 1-4 Add the Type of Weighting - always added to DF_final in same order
  # 1-4-1 Generate names of the weighting
  weighted_ens_names <- c(rep("no", times = length(res_no_weight)),
                          rep("F1", times = length(res_f1_weight)),
                          rep("Acc", times = length(res_acc_weight)))
  
  weighted_bool <- c(rep(FALSE, times = length(res_no_weight)),
                     rep(TRUE, times = length(res_f1_weight) + 
                                       length(res_acc_weight)))
  
  # 1-4-2 Add to DF_final
  DF_final$weighted_ens <- weighted_bool
  DF_final$weight_metric <- weighted_ens_names

  
  # [2] Add the settings that have been used
  # 2-1 Add used DF and used Seed
  path_split       <- unlist(strsplit(x_settings$data_path, split = "/"))
  DF_final$DF      <- path_split[length(path_split)]
  DF_final$DF_seed <- path_split[length(path_split) - 1]
  
  # 2-2 Add the infos from "x_settings"
  DF_final$performance_metric <- metric
  DF_final$reponse            <- x_settings$response
  DF_final$model_seed         <- x_settings$seed
  DF_final$num_trees          <- x_settings$num_trees
  DF_final$mtry               <- x_settings$mtry
  DF_final$min_node_size      <- x_settings$min_node_size
  
  # [3] Return 'DF_final'
  return(DF_final)
}
extract_avg_metrics_FW <- function(x, metric, train_sit) {
  "
  Function to convert a list 'x' - w/ the 2 entrances 'res_all' & 'settings' - 
  to a DF, that we can use for plotting!
  'x' is the result of a CV of one trainsetting & all has results for all of 
  the 20 testsettings! From this list we extract the metric for all of the 
  different trainsettings!
    --> Only grab the average values for the metric of all folds!
  
   Args:
    x (list)           : list filled with the metrics for every single 
                         Itteration in CV + the settings it was created from!
                         Needs the names:
                         ['res_all' & 'settings']
    metric (str)       : the metric we shall extraxt from the result list!
                         Needs to be in:
                         ['Accuracy', 'Sensitifity', 'Specificity', 'Precision', 
                          'Recall', 'F1', 'Balance_Acc', 'AUC', 'MCC']
    train_sit (int)    : Which TrainingSituation do the results come from!
                         [four different trainsettings in total s. MS Thesis!]

   Return:
    DF, suited for plotting, with all metric for all testsettings in 'x'
  "
  # [0] Check Inputs  ---------------------------------------------------------- 
  # 0-1 'x' must be a list of length 2 w/ names 'res_all' & 'settings'
  assert_list(x, len = 2)
  if (!('res_all' %in% names(x)) | !('settings' %in% names(x))) {
    stop("'x' needs 2 entrances called 'res_all' & 'settings'")
  }
  
  # 0-2 'metric' must be a string, also existing in 'x'
  assert_string(metric)
  
  #      Get all exisiting metrics 'x' contains
  exis_meas_no_weight <- names(x$res_all$full[[1]]$no_weighting)
  exis_meas_f1_weight <- names(x$res_all$full[[1]]$f1_weighting)
  exis_meas_acc_weight <- names(x$res_all$full[[1]]$acc_weighting)
  
  metric_exisits <- metric %in% exis_meas_no_weight && 
    metric %in% exis_meas_f1_weight && 
    metric %in% exis_meas_acc_weight
  
  if (!metric_exisits) stop("metric doesn't exist in 'x'")
  
  # 0-3 'train_sit' must be integer [1; 4]
  assert_int(train_sit, lower = 1, upper = 4)
  
  # [1] Extract the metric for all testsettings  -------------------------------
  # 1-1 Get the settings the model has been trained with
  x_settings <- x$settings
  
  # 1-2 Extract the metric for all 20 TestSettings ["full", "miss_1A", ...]!
  #     --> NA weights stand for trees that were pruned at the first splitvar.
  #         so that the trees could not be used for a predicition
  # 1-2-1 NO_WEIGHTING
  #       Loop over all possible TestSets
  res_no_weight <- sapply(names(x$res_all), function(setting) {
    #     For each fold check whether it contains only a single character 
    #     [-> only happens when the tree had to be pruned] - if not extract the 
    #     metric for the current 'Setting' and the current 'fold'
    sapply(1:5, function(fold) {
      if (is.character(x$res_all[[setting]][[fold]])) {
        NA
      }
      else {
        x$res_all[[setting]][[fold]][["no_weighting"]][metric]
      }
    })
  })
  
  # 1-2-2 ACC_WEIGHTING
  #       Loop over all possible TestSets
  res_acc_weight <- sapply(names(x$res_all), function(setting) {
    #     For each fold check whether it contains only a single character 
    #     [-> only happens when the tree had to be pruned] - if not extract the 
    #     metric for the current 'Setting' and the current 'fold'
    sapply(1:5, function(fold) {
      if (is.character(x$res_all[[setting]][[fold]])) {
        NA
      }
      else {
        x$res_all[[setting]][[fold]][["acc_weighting"]][metric]
      }
    })
  })
  
  # 1-2-3 F1_WEIGHTING
  #       Loop over all possible TestSets
  res_f1_weight <- sapply(names(x$res_all), function(setting) {
    #     For each fold check whether it contains only a single character 
    #     [-> only happens when the tree had to be pruned] - if not extract the 
    #     metric for the current 'Setting' and the current 'fold'
    sapply(1:5, function(fold) {
      if (is.character(x$res_all[[setting]][[fold]])) {
        NA
      }
      else {
        x$res_all[[setting]][[fold]][["f1_weighting"]][metric]
      }
    })
  })
  
  # 1-3 Convert the results to a usable DF
  DF_final <- data.frame("Trainsituation" = numeric(),
                         "Testsituation"  = character(),
                         "Metric"         = character())
  
  # 1-3-1 Loop over each of the different weighted results
  for (res_curr_train_set in list(as.data.frame(res_f1_weight), 
                                  as.data.frame(res_acc_weight),
                                  as.data.frame(res_no_weight))) {
    
    # 1-3-2 Loop over each column in 'res_curr_train_set' (each is 1 testsetting)
    #       & fill 'DF_final'
    for (j in seq_len(ncol(res_curr_train_set))) {
      
      # 1-3-2-1 Extract TrainSetting, folds, test_situation and the metrics
      #         for the given column!
      train_sit_ <- train_sit
      test_sit   <- colnames(res_curr_train_set)[j]
      metric_c   <- mean(unlist(res_curr_train_set[,j]), na.rm = T)
      
      # If metric_c contains any NA replace these by '0' or '-1' if the metric is MMC
      if (any(is.na(metric_c))) {
        to_replace <- which(is.na(metric_c))
        if (metric == "MCC") {
          replace_value <- -1
        } else {
          replace_value <- 0
        }
        
        metric_c[to_replace] <- replace_value
        metric_c             <- as.numeric(metric_c)
      }
      
      # 1-3-2-2 Bind it to the DF with all Results
      df_current <- data.frame("Trainsituation" = train_sit_,
                               "Testsituation"  = test_sit,
                               "Metric"         = metric_c)
      DF_final <- rbind(DF_final, df_current)
    }
  }
  
  # 1-4 Add the Type of Weighting - always added to DF_final in same order
  # 1-4-1 Generate names of the weighting
  weighted_ens_names <- c(rep("F1", times  = nrow(DF_final)/3),
                          rep("Acc", times = nrow(DF_final)/3),
                          rep("No", times  = nrow(DF_final)/3))
  
  weighted_bool <- c(rep(TRUE, times = nrow(DF_final)/3 + nrow(DF_final)/3),
                     rep(FALSE, times = nrow(DF_final)/3))
  
  # 1-4-2 Add to DF_final
  DF_final$weighted_ens  <- weighted_bool
  DF_final$weight_metric <- weighted_ens_names
  
  
  # [2] Add the settings that have been used
  # 2-1 Add used DF and used Seed
  path_split       <- unlist(strsplit(x_settings$data_path, split = "/"))
  DF_final$DF      <- path_split[length(path_split)]
  DF_final$DF_seed <- path_split[length(path_split) - 1]
  
  # 2-2 Add the infos from "x_settings"
  DF_final$performance_metric <- metric
  DF_final$reponse            <- x_settings$response
  DF_final$model_seed         <- x_settings$seed
  DF_final$num_trees          <- x_settings$num_trees
  DF_final$mtry               <- x_settings$mtry
  DF_final$min_node_size      <- x_settings$min_node_size
  
  # [3] Return 'DF_final'
  return(DF_final)
}

# Used for the 'SingleBlock' Approach
extract_avg_metrics_SB <- function(x, metric, train_sit) {
  "
  Function to convert a list 'x' - w/ the 2 entrances 'res_all' & 'settings' - 
  to a DF, that we can use for plotting!
  'x' is the result of a CV of one trainsetting & all has results for all of 
  the 20 testsettings! From this list we extract the metric for all of the 
  different trainsettings!
    --> Only grab the average values for the metric of all folds!
  
   Args:
    x (list)           : list filled with the metrics for every single 
                         Itteration in CV + the settings it was created from!
                         Needs the names:
                         ['res_all' & 'settings']
    metric (str)       : the metric we shall extraxt from the result list!
                         Needs to be in:
                         ['Accuracy', 'Sensitifity', 'Specificity', 'Precision', 
                          'Recall', 'F1', 'Balance_Acc', 'AUC', 'MCC']
    train_sit (int)    : Which TrainingSituation do the results come from!
                         [four different trainsettings in total s. MS Thesis!]

   Return:
    DF, suited for plotting, with all metric for all testsettings in 'x'
  "
  # [0] Check Inputs 
  # 0-1 'x' must be a list of length 2 w/ names 'res_all' & 'settings'
  assert_list(x, len = 2)
  if (!('res_all' %in% names(x)) | !('settings' %in% names(x))) {
    stop("'x' needs 2 entrances called 'res_all' & 'settings'")
  }
  
  # 0-2 'metric' must be in 'x'
  assert_string(metric)
  
  # extract all used  metrics and check whether 'metric' exist!
  met_meas <- unique(names(x$res_all[[1]]$clin_block[[1]]))
  
  if (!(metric %in% met_meas)) {
    stop("'metric' is not within these!")
  }
  
  # 0-3 'train_sit' must be integer [1; 4]
  assert_int(train_sit, lower = 1, upper = 4)
  
  # 0-4 There must be 5 single block results in 'x$res_all' [clin, A, B, C, D]
  res_learners <- sapply(c("clin_block", "A", "B", "C", "D"), 
                         function(block) block %in% names(x$res_all$full))
  
  if (!(all(res_learners))) {
    stop("x does not have 5 single_block results ")
  }
  
  # [1] Unpack the lists and put them into a single regular DF
  # 1-1 Seperate Results and Settings
  x_res      <- x$res_all
  x_settings <- x$settings
  
  # 1-2 Extract the metric from the wanted TrainSetting for all 20 TestSettings!
  #     ["full", "miss_1A", ...] & for all different single-block learners!
  # 1-2-1 Clinical RF - extract Results and average them then!
  clin_res <- as.data.frame(
    sapply(names(x_res), function(x) {
      c(x_res[[x]]$clin_block[[1]][metric],
        x_res[[x]]$clin_block[[2]][metric],
        x_res[[x]]$clin_block[[3]][metric],
        x_res[[x]]$clin_block[[4]][metric],
        x_res[[x]]$clin_block[[5]][metric])
    }))
  
  clin_res <- sapply(1:ncol(clin_res), 
                     function(col_) mean(as.numeric(unlist(clin_res[,col_])), na.rm = T))
  names(clin_res) <- names(x_res)
  
  # 1-2-2 Block 'A' - extract Results and average them then!
  b_A_res <- as.data.frame(
    sapply(names(x_res), function(x) {
      c(x_res[[x]]$A[[1]][metric],
        x_res[[x]]$A[[2]][metric],
        x_res[[x]]$A[[3]][metric],
        x_res[[x]]$A[[4]][metric],
        x_res[[x]]$A[[5]][metric])
    }))
  
  b_A_res <- sapply(1:ncol(b_A_res), 
                    function(col_) mean(as.numeric(unlist(b_A_res[,col_])), na.rm = T))
  names(b_A_res) <- names(x_res)
  
  # 1-2-3 Block 'B' - extract Results and average them then!
  b_B_res <- as.data.frame(
    sapply(names(x_res), function(x) {
      c(x_res[[x]]$B[[1]][metric],
        x_res[[x]]$B[[2]][metric],
        x_res[[x]]$B[[3]][metric],
        x_res[[x]]$B[[4]][metric],
        x_res[[x]]$B[[5]][metric])
    }))
  
  b_B_res <- sapply(1:ncol(b_B_res), 
                    function(col_) mean(as.numeric(unlist(b_B_res[,col_])), na.rm = T))
  names(b_B_res) <- names(x_res)
  
  # 1-2-4 Block 'C' - extract Results and average them then!
  b_C_res <- as.data.frame(
    sapply(names(x_res), function(x) {
      c(x_res[[x]]$C[[1]][metric],
        x_res[[x]]$C[[2]][metric],
        x_res[[x]]$C[[3]][metric],
        x_res[[x]]$C[[4]][metric],
        x_res[[x]]$C[[5]][metric])
    }))
  
  b_C_res <- sapply(1:ncol(b_C_res), 
                    function(col_) mean(as.numeric(unlist(b_C_res[,col_])), na.rm = T))
  names(b_C_res) <- names(x_res)
  
  # 1-2-5 Block 'D' - extract Results and average them then!
  b_D_res <- as.data.frame(
    sapply(names(x_res), function(x) {
      c(x_res[[x]]$D[[1]][metric],
        x_res[[x]]$D[[2]][metric],
        x_res[[x]]$D[[3]][metric],
        x_res[[x]]$D[[4]][metric],
        x_res[[x]]$D[[5]][metric])
    }))
  
  b_D_res <- sapply(1:ncol(b_D_res), 
                    function(col_) mean(as.numeric(unlist(b_D_res[,col_])), na.rm = T))
  names(b_D_res) <- names(x_res)
  
  # 1-3 Bind the results from the different single_blocks into a single DF
  DF_final <- data.frame()
  
  # 1-3-1 CLIN_RES
  DF_final <- rbind(DF_final,
                    data.frame("Trainsituation" = train_sit,
                               "Testsituation"  = names(clin_res),
                               "Metric"         = clin_res,
                               "Learn_Block"    = "Clinical"))
  
  # 1-3-2 A_RES
  DF_final <- rbind(DF_final,
                    data.frame("Trainsituation" = train_sit,
                               "Testsituation"  = names(b_A_res),
                               "Metric"         = b_A_res,
                               "Learn_Block"    = "A"))
  
  # 1-3-3 B_RES
  DF_final <- rbind(DF_final,
                    data.frame("Trainsituation" = train_sit,
                               "Testsituation"  = names(b_B_res),
                               "Metric"         = b_B_res,
                               "Learn_Block"    = "B"))
  
  # 1-3-4 C_RES
  DF_final <- rbind(DF_final,
                    data.frame("Trainsituation" = train_sit,
                               "Testsituation"  = names(b_C_res),
                               "Metric"         = b_C_res,
                               "Learn_Block"    = "C"))
  
  # 1-3-5 D_RES
  DF_final <- rbind(DF_final,
                    data.frame("Trainsituation" = train_sit,
                               "Testsituation"  = names(b_D_res),
                               "Metric"         = b_D_res,
                               "Learn_Block"    = "D"))
  
  # [2] Add the settings that have been used
  # 2-1 Add used DF and used Seed
  path_split       <- unlist(strsplit(x_settings$data_path, split = "/"))
  DF_final$DF      <- path_split[length(path_split)]
  DF_final$DF_seed <- path_split[length(path_split) - 1]
  
  # 2-2 Add the infos from "x_settings"
  DF_final$performance_metric <- metric
  DF_final$reponse            <- x_settings$response
  DF_final$model_seed         <- x_settings$seed
  DF_final$num_trees          <- x_settings$num_trees
  DF_final$mtry               <- x_settings$mtry
  DF_final$min_node_size      <- x_settings$min_node_size
  
  # [3] Return 'DF_final'
  return(DF_final)
}

# Analyse Results of the complete case Approach --- pattern 1               ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/TCGA/CompleteCase_Approach/setting1/"

# [1] Load all Results w/ Romans Approach
# 1-1 get all files in the folder, that are singleblock performances
files <- list.files(data_path)

# 1-2 Loop over all the files and extract the results
DF_all <- data.frame()
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(file_curr, metric = "F1", train_sit = 1)
  DF_all    <- rbind(DF_all, curr_df)
}

# [2] Plot the Results
# 2-1 The used metric for the comparison of the performance
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 2-2 The plot itself
ggplot(data = DF_all, aes(x = Testsituation, y = Metric)) +
  geom_boxplot(fill = 'darkolivegreen3') + 
  theme_bw() +
  ggtitle("Complete-Case Approach",
          subtitle = "TCGA - Pattern 1") +
  ylab(used_metric_) +
  xlab("Test-Situations") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24)) +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005)

# 2-3 Get summarys to the perfornance!
res_ <- sapply(unique(DF_all$Testsituation), FUN = function(x){
  summary(DF_all$Metric[DF_all$Testsituation == x])
})
names(res_) <- unique(DF_all$Testsituation)

# Analyse Results of the Single-Block Approach  --- pattern 1               ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/TCGA/SingleBlock_Approach/setting1/"

# [1] Load all Results w/ Romans Approach
# 1-1 get all files in the folder, that are singleblock performances
files <- list.files(data_path)

# 1-2 Loop over all the files and extract the results
DF_all <- data.frame()
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics_SB(x = file_curr, metric = "F1", train_sit = 1)
  DF_all    <- rbind(DF_all, curr_df)
}

# [2] Plot the Results
# 2-1-1 The used metric for the comparison of the performance
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 2-2 The plot itself
ggplot(data = DF_all, aes(x = Testsituation, y = Metric, fill = Learn_Block)) +
  geom_boxplot(position = position_dodge(preserve = "single")) + 
  theme_bw() +
  ylab(used_metric_) +
  xlab("Test-Situations") +
  ggtitle("Single-Block Approach",
          subtitle = "TCGA - Pattern 1") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24),
        legend.position = "top") +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005) +
  guides(fill = guide_legend(title = "Used feature-block: ")) +
  scale_fill_manual(values = c('darkolivegreen3', "darkorange3", "cyan4", 
                               "darkmagenta", "bisque3"))

# 2-3 Get summarys to the perfornance!
res_ <- sapply(unique(DF_all$Learn_Block), FUN = function(x){
  summary(DF_all$Metric[which(DF_all$Testsituation == "full" &
                                DF_all$Learn_Block == x)])
})
colnames(res_) <- unique(DF_all$Learn_Block)

# Analyse Results of the Imputation Approach    --- pattern 1               ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/TCGA/Imputation_Approach/setting1/"

# [1] Load all Results w/ Romans Approach
# 1-1 get all files in the folder, that are singleblock performances
files <- list.files(data_path)

# 1-2 Loop over all the files and extract the results
DF_all <- data.frame()
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(x = file_curr, metric = "F1", train_sit = 1)
  DF_all    <- rbind(DF_all, curr_df)
}

# [2] Plot the Results
# 2-1-1 The used metric for the comparison of the performance
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 2-2 The plot itself
ggplot(data = DF_all, aes(x = Testsituation, y = Metric)) +
  geom_boxplot(fill = 'darkolivegreen3') + 
  theme_bw() +
  ggtitle("Imputation Approach",
          subtitle = "TCGA - Pattern 1") +
  ylab(used_metric_) +
  xlab("Test-Situations") + 
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24)) +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005)

# 2-3 Get summarys to the perfornance!
res_ <- sapply(unique(DF_all$Testsituation), FUN = function(x){
  summary(DF_all$Metric[DF_all$Testsituation == x])
})
colnames(res_) <- unique(DF_all$Testsituation)

# Analyse Results of the BlockWise Approach     --- pattern 1               ----
data_path <- "./docs/CV_Res/TCGA/BlockWise_Approach/setting1"

# [1] Load all Results w/ Romans Approach
# 1-1 List all files from the 'data_path'
files <- list.files(data_path)

# 1-2 Loop over all the files and extract the results
DF_all <- data.frame()
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  # Extract the Metrics for the current file
  curr_df   <- extract_avg_metrics(file_curr, metric = "F1", train_sit = 1)
  
  if (grepl("acc", curr_file)) {
    curr_df$weight_metric <- "Accuracy"
  } else if (grepl("f1", curr_file)) {
    curr_df$weight_metric <- "F-1 Score"
  } else {
    curr_df$weight_metric <- "None"
  }
  
  # bind it to 'DF_all'
  DF_all <- rbind(DF_all, curr_df)
}

# [2] Plot the Results
# 2-1 The used metric for the comparison of the performance
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 2-2 Relevel the 'weight_metric' variable
DF_all$weight_metric <- factor(DF_all$weight_metric, 
                               levels = c("None", "Accuracy", "F-1 Score"))

# 2-2 Plot itself!
ggplot(data = DF_all, aes(x = Testsituation, y = Metric, fill = weight_metric)) +
  geom_boxplot() + 
  theme_bw() +
  ggtitle("Blockwise Approach",
          subtitle = "TCGA - Pattern 1") +
  ylab(used_metric_) +
  xlab("Test-Situations") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24),
        legend.position = "top") +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005) +
  guides(fill = guide_legend(title = "Weight Metric: ")) +
  scale_fill_manual(values = c( 'darkolivegreen3', "darkorange3", "cyan4"))

# 2-3 Count how often a approach has the highest median!
counter        <- c(0, 0, 0)
names(counter) <- c("No", "F1", "Acc")
for (curr_test in unique(DF_all$Testsituation)) {
  
  m_none <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$weight_metric == 'None'])
  
  m_acc <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                   DF_all$weight_metric == 'Accuracy'])
  
  m_f1s <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                  DF_all$weight_metric == 'F-1 Score'])
  
  all_res_  <- c(m_none, m_f1s, m_acc)
  max_index <- which(all_res_ == max(all_res_, na.rm = TRUE))
  
  if (length(max_index) == 1) counter[max_index] <- counter[max_index] + 1
}

# Analyse Results of the FoldWise_Approach      --- pattern 1               ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/TCGA/FoldWise_Approach/setting1"

# [1] Load all Results w/ Romans Approach
# 1-1 List all files from the 'data_path'
files <- list.files(data_path)

# 1-2 Loop over all the files and extract the results
DF_all <- data.frame()
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics_FW(file_curr, metric = "F1", train_sit = 1)
  DF_all    <- rbind(DF_all, curr_df)
}

# [2] Plot the Results
# 2-1 Extract the needed Information from 'DF_all'
# 2-1-1 The used metric for the comparison of the performance
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 2-1-2 Rename the Weight Metrics [--> clearer names]
DF_all$weight_metric[DF_all$weight_metric == "F1"]  <- "F-1 Score"
DF_all$weight_metric[DF_all$weight_metric == "Acc"] <- "Accuracy"
DF_all$weight_metric[DF_all$weight_metric == "No"]  <- "None"

DF_all$weight_metric <- factor(DF_all$weight_metric, 
                               levels = c("None", "Accuracy", "F-1 Score"))

# 2-2 Do the plot
ggplot(data = DF_all, aes(x = Testsituation, y = Metric, fill = weight_metric)) +
  geom_boxplot() + 
  theme_bw() +
  ggtitle("Foldwise Approach",
          subtitle = "TCGA - Pattern 1") +
  ylab(used_metric_) +
  xlab("Test-Situations") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24),
        legend.position = "top") +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005) +
  guides(fill = guide_legend(title = "Weight Metric: ")) +
  scale_fill_manual(values = c('darkolivegreen3', "darkorange3", "cyan4"))

# 2-3 Count how often a approach has the highest median!
counter        <- c(0, 0, 0)
names(counter) <- c("No", "F1", "Acc")
for (curr_test in unique(DF_all$Testsituation)) {
  
  m_none <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                   DF_all$weight_metric == 'None'])
  
  m_acc <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                  DF_all$weight_metric == 'Accuracy'])
  
  m_f1s <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                  DF_all$weight_metric == 'F-1 Score'])
  
  all_res_  <- c(m_none, m_f1s, m_acc)
  max_index <- which(all_res_ == max(all_res_, na.rm = TRUE))
  
  if (length(max_index) == 1) counter[max_index] <- counter[max_index] + 1
}


# Compare the different Approaches              --- pattern 1 [F1-Score]    ----
# [1] ----- CC Results
DF_CC <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/CompleteCase_Approach/setting1/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(file_curr, metric = "F1", train_sit = 1)
  DF_CC     <- rbind(DF_CC, curr_df)
}

# 1-3 Add the approach to the DF
DF_CC$Approach <- "Complete-Case"

# [2] ----- SB Results
DF_SB <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/SingleBlock_Approach/setting1/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics_SB(x = file_curr, metric = "F1", train_sit = 1)
  
  # Only keep the results of single 'A'
  curr_df   <- curr_df[curr_df$Learn_Block == 'A',]
  DF_SB     <- rbind(DF_SB, curr_df)
}

# 1-3 Add the approach to the DF
DF_SB$Approach <- "Single-Block"

# [3] ----- IMP Results
DF_IMP <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/Imputation_Approach//setting1/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(x = file_curr, metric = "F1", train_sit = 1)
  DF_IMP    <- rbind(DF_IMP, curr_df)
}

# 1-3 Add the approach to the DF
DF_IMP$Approach <- "Imputation"

# [4] ----- BW Results
DF_BW <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/BlockWise_Approach//setting1/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files[grep("f1", files)]) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(x = file_curr, metric = "F1", train_sit = 1)
  DF_BW     <- rbind(DF_BW, curr_df)
}

# 1-3 Add the approach to the DF
DF_BW$Approach <- "Block-wise"

# [5] ----- FW Results
DF_FW <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/FoldWise_Approach//setting1/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics_FW(file_curr, metric = "F1", train_sit = 1)
  DF_FW     <- rbind(DF_FW, curr_df)
}

# 1-3 Only keep Results with the F1-Score as weight metric
DF_FW <- DF_FW[DF_FW$weight_metric == "F1",]

# 1-4 Add the approach to the DF
DF_FW$Approach <- "Fold-wise"

# [6] ----- Bind the DFs and create a plot
# 6-1 Columns we need from every approach-DF
needed_cols <- c("Testsituation", "Metric", "DF", "performance_metric", "Approach")

# 6-2 Bind the DFs from the different approaches to a single DF for plots!
DF_all <- rbind(DF_CC[,needed_cols], DF_SB[,needed_cols], DF_IMP[,needed_cols],
                DF_BW[,needed_cols], DF_FW[,needed_cols])

# 6-2-1 Relevel the 'Approach' variable
DF_all$Approach <- factor(DF_all$Approach, levels = c("Fold-wise", "Block-wise",
                                                      "Imputation", "Single-Block",
                                                      "Complete-Case"), 
                          ordered = TRUE)

# 6-3 Extract the used metric
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 6-4 Do the plot
ggplot(data = DF_all, aes(x = Testsituation, y = Metric, fill = Approach)) +
  geom_boxplot(position = position_dodge(preserve = "single")) + 
  theme_bw() +
  ggtitle("Comparison of all Approaches",
          subtitle = "TCGA - Pattern 1") +
  ylab(used_metric_) +
  xlab("Test-Situations") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24),
        legend.position = "top") +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005) +
  guides(fill = guide_legend(title = "Approach: ")) +
  scale_fill_manual(values = c("darkolivegreen2", "darkorange1", 'darkorchid1', 'darkgoldenrod1', 'darkslategray3'))

# 6-5 Count how often a approach has the highest/ worst median!
best_counter        <- c(0, 0, 0, 0, 0)
names(best_counter) <- c("CC", "SB", "IMP", "BW", "FW")

worst_counter        <- c(0, 0, 0, 0, 0)
names(worst_counter) <- c("CC", "SB", "IMP", "BW", "FW")
for (curr_test in unique(DF_all$Testsituation)) {
  
  m_cc <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Complete-Case'])
  
  m_sb <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Single-Block'])
  
  m_imp <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                  DF_all$Approach == 'Imputation'])
  
  m_bw <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Block-wise'])
  
  m_fw <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Fold-wise'])
  
  all_res_  <- c(m_cc, m_sb, m_imp, m_bw, m_fw)
  
  max_index <- which(all_res_ == max(all_res_, na.rm = TRUE))
  min_index <- which(all_res_ == min(all_res_, na.rm = TRUE))
  
  if (length(max_index) == 1) best_counter[max_index]  <- best_counter[max_index] + 1
  if (length(min_index) == 1) worst_counter[min_index] <- worst_counter[min_index] + 1
}


# Compare the different Approaches              --- pattern 1 [Balanced Accuracy]   ----
# [1] ----- CC Results
DF_CC <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/CompleteCase_Approach/setting1/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(file_curr, metric = "Balance_Acc", train_sit = 1)
  DF_CC     <- rbind(DF_CC, curr_df)
}

# 1-3 Add the approach to the DF
DF_CC$Approach <- "Complete-Case"

# [2] ----- SB Results
DF_SB <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/SingleBlock_Approach/setting1/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics_SB(x = file_curr, metric = "Balance_Acc", train_sit = 1)
  
  # Only keep the results of single 'A'
  curr_df   <- curr_df[curr_df$Learn_Block == 'A',]
  DF_SB     <- rbind(DF_SB, curr_df)
}

# 1-3 Add the approach to the DF
DF_SB$Approach <- "Single-Block"

# [3] ----- IMP Results
DF_IMP <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/Imputation_Approach//setting1/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(x = file_curr, metric = "Balance_Acc", train_sit = 1)
  DF_IMP    <- rbind(DF_IMP, curr_df)
}

# 1-3 Add the approach to the DF
DF_IMP$Approach <- "Imputation"

# [4] ----- BW Results
DF_BW <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/BlockWise_Approach//setting1/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files[grep("f1", files)]) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(x = file_curr, metric = "Balance_Acc", train_sit = 1)
  DF_BW     <- rbind(DF_BW, curr_df)
}

# 1-3 Add the approach to the DF
DF_BW$Approach <- "Block-wise"

# [5] ----- FW Results
DF_FW <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/FoldWise_Approach//setting1/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics_FW(file_curr, metric = "Balance_Acc", train_sit = 1)
  DF_FW     <- rbind(DF_FW, curr_df)
}

# 1-3 Only keep Results with the F1-Score as weight metric
DF_FW <- DF_FW[DF_FW$weight_metric == "F1",]

# 1-4 Add the approach to the DF
DF_FW$Approach <- "Fold-wise"

# [6] ----- Bind the DFs and create a plot
# 6-1 Columns we need from every approach-DF
needed_cols <- c("Testsituation", "Metric", "DF", "performance_metric", "Approach")

# 6-2 Bind the DFs from the different approaches to a single DF for plots!
DF_all <- rbind(DF_CC[,needed_cols], DF_SB[,needed_cols], DF_IMP[,needed_cols],
                DF_BW[,needed_cols], DF_FW[,needed_cols])

# 6-2-1 Relevel the 'Approach' variable
DF_all$Approach <- factor(DF_all$Approach, levels = c("Fold-wise", "Block-wise",
                                                      "Imputation", "Single-Block",
                                                      "Complete-Case"), 
                          ordered = TRUE)

# 6-3 Extract the used metric
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else if (DF_all$performance_metric[1] == "Balance_Acc") {
  used_metric_ <- "Metric: Balanced Accuracy" 
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 6-4 Do the plot
ggplot(data = DF_all, aes(x = Testsituation, y = Metric, fill = Approach)) +
  geom_boxplot(position = position_dodge(preserve = "single")) + 
  theme_bw() +
  ggtitle("Comparison of all Approaches",
          subtitle = "TCGA - Pattern 1") +
  ylab(used_metric_) +
  xlab("Test-Situations") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24),
        legend.position = "top") +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005) +
  guides(fill = guide_legend(title = "Approach: ")) +
  scale_fill_manual(values = c("darkolivegreen2", "darkorange1", 'darkorchid1', 'darkgoldenrod1', 'darkslategray3'))

# 6-5 Count how often a approach has the highest median!
counter        <- c(0, 0, 0, 0, 0)
names(counter) <- c("CC", "SB", "IMP", "BW", "FW")
for (curr_test in unique(DF_all$Testsituation)) {
  
  m_cc <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Complete-Case'])
  
  m_sb <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Single-Block'])
  
  m_imp <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                  DF_all$Approach == 'Imputation'])
  
  m_bw <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Block-wise'])
  
  m_fw <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Fold-wise'])
  
  all_res_  <- c(m_cc, m_sb, m_imp, m_bw, m_fw)
  max_index <- which(all_res_ == max(all_res_, na.rm = TRUE))
  
  if (length(max_index) == 1) counter[max_index] <- counter[max_index] + 1
}
# Compare the different Approaches              --- pattern 1 [MCC]         ----
# [1] ----- CC Results
DF_CC <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/CompleteCase_Approach/setting1/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(file_curr, metric = "MCC", train_sit = 1)
  DF_CC     <- rbind(DF_CC, curr_df)
}

# 1-3 Add the approach to the DF
DF_CC$Approach <- "Complete-Case"

# [2] ----- SB Results
DF_SB <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/SingleBlock_Approach/setting1/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics_SB(x = file_curr, metric = "MCC", train_sit = 1)
  
  # Only keep the results of single 'A'
  curr_df   <- curr_df[curr_df$Learn_Block == 'A',]
  DF_SB     <- rbind(DF_SB, curr_df)
}

# 1-3 Add the approach to the DF
DF_SB$Approach <- "Single-Block"

# [3] ----- IMP Results
DF_IMP <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/Imputation_Approach//setting1/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(x = file_curr, metric = "MCC", train_sit = 1)
  DF_IMP    <- rbind(DF_IMP, curr_df)
}

# 1-3 Add the approach to the DF
DF_IMP$Approach <- "Imputation"

# [4] ----- BW Results
DF_BW <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/BlockWise_Approach//setting1/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files[grep("f1", files)]) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(x = file_curr, metric = "MCC", train_sit = 1)
  DF_BW     <- rbind(DF_BW, curr_df)
}

# 1-3 Add the approach to the DF
DF_BW$Approach <- "Block-wise"

# [5] ----- FW Results
DF_FW <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/FoldWise_Approach//setting1/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics_FW(file_curr, metric = "MCC", train_sit = 1)
  DF_FW     <- rbind(DF_FW, curr_df)
}

# 1-3 Only keep Results with the F1-Score as weight metric
DF_FW <- DF_FW[DF_FW$weight_metric == "F1",]

# 1-4 Add the approach to the DF
DF_FW$Approach <- "Fold-wise"

# [6] ----- Bind the DFs and create a plot
# 6-1 Columns we need from every approach-DF
needed_cols <- c("Testsituation", "Metric", "DF", "performance_metric", "Approach")

# 6-2 Bind the DFs from the different approaches to a single DF for plots!
DF_all <- rbind(DF_CC[,needed_cols], DF_SB[,needed_cols], DF_IMP[,needed_cols],
                DF_BW[,needed_cols], DF_FW[,needed_cols])

# 6-2-1 Relevel the 'Approach' variable
DF_all$Approach <- factor(DF_all$Approach, levels = c("Fold-wise", "Block-wise",
                                                      "Imputation", "Single-Block",
                                                      "Complete-Case"), 
                          ordered = TRUE)

# 6-3 Extract the used metric
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else if (DF_all$performance_metric[1] == "Balance_Acc") {
  used_metric_ <- "Metric: Balanced Accuracy" 
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 6-4 Do the plot
ggplot(data = DF_all, aes(x = Testsituation, y = Metric, fill = Approach)) +
  geom_boxplot(position = position_dodge(preserve = "single")) + 
  theme_bw() +
  ggtitle("Comparison of all Approaches",
          subtitle = "TCGA - Pattern 1") +
  ylab(used_metric_) +
  xlab("Test-Situations") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24),
        legend.position = "top") +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005) +
  guides(fill = guide_legend(title = "Approach: ")) +
  scale_fill_manual(values = c("darkolivegreen2", "darkorange1", 'darkorchid1', 'darkgoldenrod1', 'darkslategray3'))

# 6-5 Count how often a approach has the highest median!
counter        <- c(0, 0, 0, 0, 0)
names(counter) <- c("CC", "SB", "IMP", "BW", "FW")
for (curr_test in unique(DF_all$Testsituation)) {
  
  m_cc <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Complete-Case'])
  
  m_sb <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Single-Block'])
  
  m_imp <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                  DF_all$Approach == 'Imputation'])
  
  m_bw <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Block-wise'])
  
  m_fw <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Fold-wise'])
  
  all_res_  <- c(m_cc, m_sb, m_imp, m_bw, m_fw)
  max_index <- which(all_res_ == max(all_res_, na.rm = TRUE))
  
  if (length(max_index) == 1) counter[max_index] <- counter[max_index] + 1
} 

###################################################################### PATTERN 2
# Analyse Results of the complete case Approach --- pattern 2               ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/TCGA/CompleteCase_Approach/setting2/"

# [1] Load all Results w/ Romans Approach
# 1-1 get all files in the folder, that are singleblock performances
files <- list.files(data_path)

# 1-2 Loop over all the files and extract the results
DF_all <- data.frame()
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(file_curr, metric = "F1", train_sit = 1)
  DF_all    <- rbind(DF_all, curr_df)
}

# [2] Plot the Results
# 2-1-1 The used metric for the comparison of the performance
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 2-2 The plot itself
ggplot(data = DF_all, aes(x = Testsituation, y = Metric)) +
  geom_boxplot(fill = 'darkolivegreen3') + 
  theme_bw() +
  ggtitle("Complete-Case Approach",
          subtitle = "TCGA - Pattern 2") +
  ylab(used_metric_) +
  xlab("Test-Situations") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24)) +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005)

# 2-3 Get summarys to the perfornance!
res_           <- sapply(unique(DF_all$Testsituation), FUN = function(x){
  summary(DF_all$Metric[DF_all$Testsituation == x])
})
res_           <- as.data.frame(res_)
colnames(res_) <- unique(DF_all$Testsituation)

# Analyse Results of the Single-Block Approach  --- pattern 2               ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/TCGA/SingleBlock_Approach/setting2/"

# [1] Load all Results w/ Romans Approach
# 1-1 get all files in the folder, that are singleblock performances
files <- list.files(data_path)

# 1-2 Loop over all the files and extract the results
DF_all <- data.frame()
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics_SB(x = file_curr, metric = "F1", train_sit = 1)
  DF_all    <- rbind(DF_all, curr_df)
}

# [2] Plot the Results
# 2-1-1 The used metric for the comparison of the performance
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 2-2 The plot itself
# 2-2 The plot itself
ggplot(data = DF_all, aes(x = Testsituation, y = Metric, fill = Learn_Block)) +
  geom_boxplot(position = position_dodge(preserve = "single")) + 
  theme_bw() +
  ylab(used_metric_) +
  xlab("Test-Situations") +
  ggtitle("Single-Block Approach",
          subtitle = "TCGA - Pattern 2") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24),
        legend.position = "top") +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005) +
  guides(fill = guide_legend(title = "Used feature-block: ")) +
  scale_fill_manual(values = c('darkolivegreen3', "darkorange3", "cyan4", 
                               "darkmagenta", "bisque3"))


# 2-3 Get summarys to the perfornance!
res_ <- sapply(unique(DF_all$Learn_Block), FUN = function(x){
  summary(DF_all$Metric[which(DF_all$Testsituation == "full" &
                                DF_all$Learn_Block == x)])
})
colnames(res_) <- unique(DF_all$Learn_Block)

# Analyse Results of the Imputation Approach    --- pattern 2               ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/TCGA/Imputation_Approach/setting2/"

# [1] Load all Results w/ Romans Approach
# 1-1 get all files in the folder, that are singleblock performances
files <- list.files(data_path)

# 1-2 Loop over all the files and extract the results
DF_all <- data.frame()
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(x = file_curr, metric = "F1", train_sit = 1)
  DF_all    <- rbind(DF_all, curr_df)
}

# [2] Plot the Results
# 2-1-1 The used metric for the comparison of the performance
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 2-2 The plot itself
ggplot(data = DF_all, aes(x = Testsituation, y = Metric)) +
  geom_boxplot(fill = 'darkolivegreen3') + 
  theme_bw() +
  ggtitle("Imputation Approach",
          subtitle = "TCGA - Pattern 2") +
  ylab(used_metric_) +
  xlab("Test-Situations") + 
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24)) +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005)

# 2-3 Get summarys to the perfornance!
res_ <- sapply(unique(DF_all$Testsituation), FUN = function(x){
  summary(DF_all$Metric[DF_all$Testsituation == x])
})
colnames(res_) <- unique(DF_all$Testsituation)

# Analyse Results of the BlockWise Approach     --- pattern 2               ----
data_path <- "./docs/CV_Res/TCGA/BlockWise_Approach/setting2"

# [1] Load all Results w/ Romans Approach
# 1-1 List all files from the 'data_path'
files <- list.files(data_path)

# 1-2 Loop over all the files and extract the results
DF_all <- data.frame()
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  # Extract the Metrics for the current file
  curr_df   <- extract_avg_metrics(file_curr, metric = "F1", train_sit = 1)
  
  if (grepl("acc", curr_file)) {
    curr_df$weight_metric <- "Accuracy"
  } else if (grepl("f1", curr_file)) {
    curr_df$weight_metric <- "F-1 Score"
  } else {
    curr_df$weight_metric <- "None"
  }
  
  # bind it to 'DF_all'
  DF_all <- rbind(DF_all, curr_df)
}

# [2] Plot the Results
# 2-1 The used metric for the comparison of the performance
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 2-2 Relevel the 'weight_metric' variable
DF_all$weight_metric <- factor(DF_all$weight_metric, 
                               levels = c("None", "Accuracy", "F-1 Score"))

# 2-2 Plot itself!
ggplot(data = DF_all, aes(x = Testsituation, y = Metric, fill = weight_metric)) +
  geom_boxplot() + 
  theme_bw() +
  ggtitle("Blockwise Approach",
          subtitle = "TCGA - Pattern 2") +
  ylab(used_metric_) +
  xlab("Test-Situations") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24),
        legend.position = "top") +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005) +
  guides(fill = guide_legend(title = "Weight Metric: ")) +
  scale_fill_manual(values = c( 'darkolivegreen3', "darkorange3", "cyan4"))

# 2-3 Count how often a approach has the highest median!
counter        <- c(0, 0, 0)
names(counter) <- c("No", "F1", "Acc")
for (curr_test in unique(DF_all$Testsituation)) {
  
  m_none <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                   DF_all$weight_metric == 'None'])
  
  m_acc <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                  DF_all$weight_metric == 'Accuracy'])
  
  m_f1s <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                  DF_all$weight_metric == 'F-1 Score'])
  
  all_res_  <- c(m_none, m_f1s, m_acc)
  max_index <- which(all_res_ == max(all_res_, na.rm = TRUE))
  
  if (length(max_index) == 1) counter[max_index] <- counter[max_index] + 1
}

# Analyse Results of the FoldWise_Approach      --- pattern 2               ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/TCGA/FoldWise_Approach/setting2"

# [1] Load all Results w/ Romans Approach
# 1-1 List all files from the 'data_path'
files <- list.files(data_path)

# 1-2 Loop over all the files and extract the results
DF_all <- data.frame()
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics_FW(file_curr, metric = "F1", train_sit = 1)
  DF_all    <- rbind(DF_all, curr_df)
}

# [2] Plot the Results
# 2-1 Extract the needed Information from 'DF_all'
# 2-1-1 The used metric for the comparison of the performance
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 2-1-2 Rename the Weight Metrics [--> clearer names]
DF_all$weight_metric[DF_all$weight_metric == "F1"]  <- "F-1 Score"
DF_all$weight_metric[DF_all$weight_metric == "Acc"] <- "Accuracy"
DF_all$weight_metric[DF_all$weight_metric == "No"]  <- "None"

DF_all$weight_metric <- factor(DF_all$weight_metric, 
                               levels = c("None", "Accuracy", "F-1 Score"))

# 2-2 Do the plot
ggplot(data = DF_all, aes(x = Testsituation, y = Metric, fill = weight_metric)) +
  geom_boxplot() + 
  theme_bw() +
  ggtitle("Foldwise Approach",
          subtitle = "TCGA - Pattern 2") +
  ylab(used_metric_) +
  xlab("Test-Situations") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24),
        legend.position = "top") +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005) +
  guides(fill = guide_legend(title = "Weight Metric: ")) +
  scale_fill_manual(values = c('darkolivegreen3', "darkorange3", "cyan4"))

# 2-3 Count how often a approach has the highest median!
counter        <- c(0, 0, 0)
names(counter) <- c("No", "F1", "Acc")
for (curr_test in unique(DF_all$Testsituation)) {
  
  m_none <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                   DF_all$weight_metric == 'None'])
  
  m_acc <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                  DF_all$weight_metric == 'Accuracy'])
  
  m_f1s <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                  DF_all$weight_metric == 'F-1 Score'])
  
  all_res_  <- c(m_none, m_f1s, m_acc)
  max_index <- which(all_res_ == max(all_res_, na.rm = TRUE))
  
  if (length(max_index) == 1) counter[max_index] <- counter[max_index] + 1
}


# Compare the different Approaches              --- pattern 2 [F1-Score]    ----
# [1] ----- CC Results
DF_CC <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/CompleteCase_Approach/setting2/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(file_curr, metric = "F1", train_sit = 1)
  DF_CC     <- rbind(DF_CC, curr_df)
}

# 1-3 Add the approach to the DF
DF_CC$Approach <- "Complete-Case"

# [2] ----- SB Results
DF_SB <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/SingleBlock_Approach/setting2/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics_SB(x = file_curr, metric = "F1", train_sit = 1)
  
  # Only keep the results of single 'D'
  curr_df   <- curr_df[curr_df$Learn_Block == 'D',]
  DF_SB     <- rbind(DF_SB, curr_df)
}

# 1-3 Add the approach to the DF
DF_SB$Approach <- "Single-Block"

# [3] ----- IMP Results
DF_IMP <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/Imputation_Approach//setting2/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(x = file_curr, metric = "F1", train_sit = 1)
  DF_IMP    <- rbind(DF_IMP, curr_df)
}

# 1-3 Add the approach to the DF
DF_IMP$Approach <- "Imputation"

# [4] ----- BW Results
DF_BW <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/BlockWise_Approach//setting2/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files[grep("f1", files)]) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(x = file_curr, metric = "F1", train_sit = 1)
  DF_BW     <- rbind(DF_BW, curr_df)
}

# 1-3 Add the approach to the DF
DF_BW$Approach <- "Block-wise"

# [5] ----- FW Results
DF_FW <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/FoldWise_Approach//setting2/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics_FW(file_curr, metric = "F1", train_sit = 1)
  DF_FW     <- rbind(DF_FW, curr_df)
}

# 1-3 Only keep Results with the F1-Score as weight metric
DF_FW <- DF_FW[DF_FW$weight_metric == "F1",]

# 1-4 Add the approach to the DF
DF_FW$Approach <- "Fold-wise"

# [6] ----- Bind the DFs and create a plot
# 6-1 Columns we need from every approach-DF
needed_cols <- c("Testsituation", "Metric", "DF", "performance_metric", "Approach")

# 6-2 Bind the DFs from the different approaches to a single DF for plots!
DF_all <- rbind(DF_CC[,needed_cols], DF_SB[,needed_cols], DF_IMP[,needed_cols],
                DF_BW[,needed_cols], DF_FW[,needed_cols])

# 6-2-1 Relevel the 'Approach' variable
DF_all$Approach <- factor(DF_all$Approach, levels = c("Fold-wise", "Block-wise",
                                                      "Imputation", "Single-Block",
                                                      "Complete-Case"), 
                          ordered = TRUE)

# 6-3 Extract the used metric
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 6-4 Do the plot
ggplot(data = DF_all, aes(x = Testsituation, y = Metric, fill = Approach)) +
  geom_boxplot(position = position_dodge(preserve = "single")) + 
  theme_bw() +
  ggtitle("Comparison of all Approaches",
          subtitle = "TCGA - Pattern 2") +
  ylab(used_metric_) +
  xlab("Test-Situations") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24),
        legend.position = "top") +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005) +
  guides(fill = guide_legend(title = "Approach: ")) +
  scale_fill_manual(values = c("darkolivegreen2", "darkorange1", 'darkorchid1', 'darkgoldenrod1', 'darkslategray3'))

# 6-5 Count how often a approach has the highest/ worst median!
best_counter        <- c(0, 0, 0, 0, 0)
names(best_counter) <- c("CC", "SB", "IMP", "BW", "FW")

worst_counter        <- c(0, 0, 0, 0, 0)
names(worst_counter) <- c("CC", "SB", "IMP", "BW", "FW")
for (curr_test in unique(DF_all$Testsituation)) {
  
  m_cc <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Complete-Case'])
  
  m_sb <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Single-Block'])
  
  m_imp <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                  DF_all$Approach == 'Imputation'])
  
  m_bw <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Block-wise'])
  
  m_fw <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Fold-wise'])
  
  all_res_  <- c(m_cc, m_sb, m_imp, m_bw, m_fw)
  
  max_index <- which(all_res_ == max(all_res_, na.rm = TRUE))
  min_index <- which(all_res_ == min(all_res_, na.rm = TRUE))
  
  if (length(max_index) == 1) best_counter[max_index]  <- best_counter[max_index] + 1
  if (length(min_index) == 1) worst_counter[min_index] <- worst_counter[min_index] + 1
}

# Compare the different Approaches              --- pattern 2 [Balanced Accuracy]   ----
# [1] ----- CC Results
DF_CC <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/CompleteCase_Approach/setting2/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(file_curr, metric = "Balance_Acc", train_sit = 1)
  DF_CC     <- rbind(DF_CC, curr_df)
}

# 1-3 Add the approach to the DF
DF_CC$Approach <- "Complete-Case"

# [2] ----- SB Results
DF_SB <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/SingleBlock_Approach/setting2/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics_SB(x = file_curr, metric = "Balance_Acc", train_sit = 1)
  
  # Only keep the results of single 'D'
  curr_df   <- curr_df[curr_df$Learn_Block == 'D',]
  DF_SB     <- rbind(DF_SB, curr_df)
}

# 1-3 Add the approach to the DF
DF_SB$Approach <- "Single-Block"

# [3] ----- IMP Results
DF_IMP <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/Imputation_Approach//setting2/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(x = file_curr, metric = "Balance_Acc", train_sit = 1)
  DF_IMP    <- rbind(DF_IMP, curr_df)
}

# 1-3 Add the approach to the DF
DF_IMP$Approach <- "Imputation"

# [4] ----- BW Results
DF_BW <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/BlockWise_Approach//setting2/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files[grep("f1", files)]) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(x = file_curr, metric = "Balance_Acc", train_sit = 1)
  DF_BW     <- rbind(DF_BW, curr_df)
}

# 1-3 Add the approach to the DF
DF_BW$Approach <- "Block-wise"

# [5] ----- FW Results
DF_FW <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/FoldWise_Approach//setting2/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics_FW(file_curr, metric = "Balance_Acc", train_sit = 1)
  DF_FW     <- rbind(DF_FW, curr_df)
}

# 1-3 Only keep Results with the F1-Score as weight metric
DF_FW <- DF_FW[DF_FW$weight_metric == "F1",]

# 1-4 Add the approach to the DF
DF_FW$Approach <- "Fold-wise"

# [6] ----- Bind the DFs and create a plot
# 6-1 Columns we need from every approach-DF
needed_cols <- c("Testsituation", "Metric", "DF", "performance_metric", "Approach")

# 6-2 Bind the DFs from the different approaches to a single DF for plots!
DF_all <- rbind(DF_CC[,needed_cols], DF_SB[,needed_cols], DF_IMP[,needed_cols],
                DF_BW[,needed_cols], DF_FW[,needed_cols])

# 6-2-1 Relevel the 'Approach' variable
DF_all$Approach <- factor(DF_all$Approach, levels = c("Fold-wise", "Block-wise",
                                                      "Imputation", "Single-Block",
                                                      "Complete-Case"), 
                          ordered = TRUE)

# 6-3 Extract the used metric
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else if (DF_all$performance_metric[1] == "Balance_Acc") {
  used_metric_ <- "Metric: Balanced Accuracy" 
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 6-4 Do the plot
ggplot(data = DF_all, aes(x = Testsituation, y = Metric, fill = Approach)) +
  geom_boxplot(position = position_dodge(preserve = "single")) + 
  theme_bw() +
  ggtitle("Comparison of all Approaches",
          subtitle = "TCGA - Pattern 2") +
  ylab(used_metric_) +
  xlab("Test-Situations") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24),
        legend.position = "top") +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005) +
  guides(fill = guide_legend(title = "Approach: ")) +
  scale_fill_manual(values = c("darkolivegreen2", "darkorange1", 'darkorchid1', 'darkgoldenrod1', 'darkslategray3'))

# 6-5 Count how often a approach has the highest median!
counter        <- c(0, 0, 0, 0, 0)
names(counter) <- c("CC", "SB", "IMP", "BW", "FW")
for (curr_test in unique(DF_all$Testsituation)) {
  
  m_cc <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Complete-Case'])
  
  m_sb <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Single-Block'])
  
  m_imp <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                  DF_all$Approach == 'Imputation'])
  
  m_bw <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Block-wise'])
  
  m_fw <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Fold-wise'])
  
  all_res_  <- c(m_cc, m_sb, m_imp, m_bw, m_fw)
  max_index <- which(all_res_ == max(all_res_, na.rm = TRUE))
  
  if (length(max_index) == 1) counter[max_index] <- counter[max_index] + 1
}
# Compare the different Approaches              --- pattern 2 [MCC]         ----
# [1] ----- CC Results
DF_CC <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/CompleteCase_Approach/setting2/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(file_curr, metric = "MCC", train_sit = 1)
  DF_CC     <- rbind(DF_CC, curr_df)
}

# 1-3 Add the approach to the DF
DF_CC$Approach <- "Complete-Case"

# [2] ----- SB Results
DF_SB <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/SingleBlock_Approach/setting2/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics_SB(x = file_curr, metric = "MCC", train_sit = 1)
  
  # Only keep the results of single 'D'
  curr_df   <- curr_df[curr_df$Learn_Block == 'D',]
  DF_SB     <- rbind(DF_SB, curr_df)
}

# 1-3 Add the approach to the DF
DF_SB$Approach <- "Single-Block"

# [3] ----- IMP Results
DF_IMP <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/Imputation_Approach//setting2/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(x = file_curr, metric = "MCC", train_sit = 1)
  DF_IMP    <- rbind(DF_IMP, curr_df)
}

# 1-3 Add the approach to the DF
DF_IMP$Approach <- "Imputation"

# [4] ----- BW Results
DF_BW <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/BlockWise_Approach//setting2/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files[grep("f1", files)]) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(x = file_curr, metric = "MCC", train_sit = 1)
  DF_BW     <- rbind(DF_BW, curr_df)
}

# 1-3 Add the approach to the DF
DF_BW$Approach <- "Block-wise"

# [5] ----- FW Results
DF_FW <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/FoldWise_Approach//setting2/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics_FW(file_curr, metric = "MCC", train_sit = 1)
  DF_FW     <- rbind(DF_FW, curr_df)
}

# 1-3 Only keep Results with the F1-Score as weight metric
DF_FW <- DF_FW[DF_FW$weight_metric == "F1",]

# 1-4 Add the approach to the DF
DF_FW$Approach <- "Fold-wise"

# [6] ----- Bind the DFs and create a plot
# 6-1 Columns we need from every approach-DF
needed_cols <- c("Testsituation", "Metric", "DF", "performance_metric", "Approach")

# 6-2 Bind the DFs from the different approaches to a single DF for plots!
DF_all <- rbind(DF_CC[,needed_cols], DF_SB[,needed_cols], DF_IMP[,needed_cols],
                DF_BW[,needed_cols], DF_FW[,needed_cols])

# 6-2-1 Relevel the 'Approach' variable
DF_all$Approach <- factor(DF_all$Approach, levels = c("Fold-wise", "Block-wise",
                                                      "Imputation", "Single-Block",
                                                      "Complete-Case"), 
                          ordered = TRUE)

# 6-3 Extract the used metric
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else if (DF_all$performance_metric[1] == "Balance_Acc") {
  used_metric_ <- "Metric: Balanced Accuracy" 
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 6-4 Do the plot
ggplot(data = DF_all, aes(x = Testsituation, y = Metric, fill = Approach)) +
  geom_boxplot(position = position_dodge(preserve = "single")) + 
  theme_bw() +
  ggtitle("Comparison of all Approaches",
          subtitle = "TCGA - Pattern 2") +
  ylab(used_metric_) +
  xlab("Test-Situations") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24),
        legend.position = "top") +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005) +
  guides(fill = guide_legend(title = "Approach: ")) +
  scale_fill_manual(values = c("darkolivegreen2", "darkorange1", 'darkorchid1', 'darkgoldenrod1', 'darkslategray3'))

# 6-5 Count how often a approach has the highest median!
counter        <- c(0, 0, 0, 0, 0)
names(counter) <- c("CC", "SB", "IMP", "BW", "FW")
for (curr_test in unique(DF_all$Testsituation)) {
  
  m_cc <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Complete-Case'])
  
  m_sb <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Single-Block'])
  
  m_imp <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                  DF_all$Approach == 'Imputation'])
  
  m_bw <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Block-wise'])
  
  m_fw <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Fold-wise'])
  
  all_res_  <- c(m_cc, m_sb, m_imp, m_bw, m_fw)
  max_index <- which(all_res_ == max(all_res_, na.rm = TRUE))
  
  if (length(max_index) == 1) counter[max_index] <- counter[max_index] + 1
} 
# Analyse Results of the complete case Approach --- pattern 3               ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/TCGA/CompleteCase_Approach/setting3/"

# [1] Load all Results w/ Romans Approach
# 1-1 get all files in the folder, that are singleblock performances
files <- list.files(data_path)

# 1-2 Loop over all the files and extract the results
DF_all <- data.frame()
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(file_curr, metric = "F1", train_sit = 1)
  DF_all    <- rbind(DF_all, curr_df)
}

# [2] Plot the Results
# 2-1-1 The used metric for the comparison of the performance
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 2-2 The plot itself
ggplot(data = DF_all, aes(x = Testsituation, y = Metric)) +
  geom_boxplot(fill = 'darkolivegreen3') + 
  theme_bw() +
  ggtitle("Complete-Case Approach",
          subtitle = "TCGA - Pattern 3") +
  ylab(used_metric_) +
  xlab("Test-Situations") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24)) +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005)

# 2-3 Get summarys to the perfornance!
res_ <- sapply(unique(DF_all$Testsituation), FUN = function(x){
  summary(DF_all$Metric[DF_all$Testsituation == x])
})
names(res_) <- unique(DF_all$Testsituation)

# Analyse Results of the Single-Block Approach  --- pattern 3               ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/TCGA/SingleBlock_Approach/setting3/"

# [1] Load all Results w/ Romans Approach
# 1-1 get all files in the folder, that are singleblock performances
files <- list.files(data_path)

# 1-2 Loop over all the files and extract the results
DF_all <- data.frame()
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics_SB(x = file_curr, metric = "F1", train_sit = 1)
  DF_all    <- rbind(DF_all, curr_df)
}

# [2] Plot the Results
# 2-1-1 The used metric for the comparison of the performance
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 2-2 The plot itself
ggplot(data = DF_all, aes(x = Testsituation, y = Metric, fill = Learn_Block)) +
  geom_boxplot(position = position_dodge(preserve = "single")) + 
  theme_bw() +
  ylab(used_metric_) +
  xlab("Test-Situations") +
  ggtitle("Single-Block Approach",
          subtitle = "TCGA - Pattern 3") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24),
        legend.position = "top") +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005) +
  guides(fill = guide_legend(title = "Used feature-block: ")) +
  scale_fill_manual(values = c('darkolivegreen3', "darkorange3", "cyan4", 
                               "darkmagenta", "bisque3"))

# 2-3 Get summarys to the perfornance!
res_ <- sapply(unique(DF_all$Learn_Block), FUN = function(x){
  summary(DF_all$Metric[which(DF_all$Testsituation == "full" &
                                DF_all$Learn_Block == x)])
})
colnames(res_) <- unique(DF_all$Learn_Block)

# Analyse Results of the Imputation Approach    --- pattern 3               ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/TCGA/Imputation_Approach/setting3/"

# [1] Load all Results w/ Romans Approach
# 1-1 get all files in the folder, that are singleblock performances
files <- list.files(data_path)

# 1-2 Loop over all the files and extract the results
DF_all <- data.frame()
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(x = file_curr, metric = "F1", train_sit = 1)
  DF_all    <- rbind(DF_all, curr_df)
}

# [2] Plot the Results
# 2-1-1 The used metric for the comparison of the performance
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 2-2 The plot itself
ggplot(data = DF_all, aes(x = Testsituation, y = Metric)) +
  geom_boxplot(fill = 'darkolivegreen3') + 
  theme_bw() +
  ggtitle("Imputation Approach",
          subtitle = "TCGA - Pattern 3") +
  ylab(used_metric_) +
  xlab("Test-Situations") + 
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24)) +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005)

# 2-3 Get summarys to the perfornance!
res_ <- sapply(unique(DF_all$Testsituation), FUN = function(x){
  summary(DF_all$Metric[DF_all$Testsituation == x])
})
colnames(res_) <- unique(DF_all$Testsituation)

# Analyse Results of the BlockWise Approach     --- pattern 3               ----
data_path <- "./docs/CV_Res/TCGA/BlockWise_Approach/setting3"

# [1] Load all Results w/ Romans Approach
# 1-1 List all files from the 'data_path'
files <- list.files(data_path)

# 1-2 Loop over all the files and extract the results
DF_all <- data.frame()
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  # Extract the Metrics for the current file
  curr_df   <- extract_avg_metrics(file_curr, metric = "F1", train_sit = 1)
  
  if (grepl("acc", curr_file)) {
    curr_df$weight_metric <- "Accuracy"
  } else if (grepl("f1", curr_file)) {
    curr_df$weight_metric <- "F-1 Score"
  } else {
    curr_df$weight_metric <- "None"
  }
  
  # bind it to 'DF_all'
  DF_all <- rbind(DF_all, curr_df)
}

# [2] Plot the Results
# 2-1 The used metric for the comparison of the performance
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 2-2 Relevel the 'weight_metric' variable
DF_all$weight_metric <- factor(DF_all$weight_metric, 
                               levels = c("None", "Accuracy", "F-1 Score"))

# 2-2 Plot itself!
ggplot(data = DF_all, aes(x = Testsituation, y = Metric, fill = weight_metric)) +
  geom_boxplot() + 
  theme_bw() +
  ggtitle("Blockwise Approach",
          subtitle = "TCGA - Pattern 3") +
  ylab(used_metric_) +
  xlab("Test-Situations") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24),
        legend.position = "top") +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005) +
  guides(fill = guide_legend(title = "Weight Metric: ")) +
  scale_fill_manual(values = c( 'darkolivegreen3', "darkorange3", "cyan4"))

# 2-3 Count how often a approach has the highest median!
counter        <- c(0, 0, 0)
names(counter) <- c("No", "F1", "Acc")
for (curr_test in unique(DF_all$Testsituation)) {
  
  m_none <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                   DF_all$weight_metric == 'None'])
  
  m_acc <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                  DF_all$weight_metric == 'Accuracy'])
  
  m_f1s <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                  DF_all$weight_metric == 'F-1 Score'])
  
  all_res_  <- c(m_none, m_f1s, m_acc)
  max_index <- which(all_res_ == max(all_res_, na.rm = TRUE))
  
  if (length(max_index) == 1) counter[max_index] <- counter[max_index] + 1
}

# Analyse Results of the FoldWise_Approach      --- pattern 3               ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/TCGA/FoldWise_Approach/setting3"

# [1] Load all Results w/ Romans Approach
# 1-1 List all files from the 'data_path'
files <- list.files(data_path)

# 1-2 Loop over all the files and extract the results
DF_all <- data.frame()
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics_FW(file_curr, metric = "F1", train_sit = 1)
  DF_all    <- rbind(DF_all, curr_df)
}

# [2] Plot the Results
# 2-1 Extract the needed Information from 'DF_all'
# 2-1-1 The used metric for the comparison of the performance
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 2-1-2 Rename the Weight Metrics [--> clearer names]
DF_all$weight_metric[DF_all$weight_metric == "F1"]  <- "F-1 Score"
DF_all$weight_metric[DF_all$weight_metric == "Acc"] <- "Accuracy"
DF_all$weight_metric[DF_all$weight_metric == "No"]  <- "None"

DF_all$weight_metric <- factor(DF_all$weight_metric, 
                               levels = c("None", "Accuracy", "F-1 Score"))

# 2-2 Do the plot
ggplot(data = DF_all, aes(x = Testsituation, y = Metric, fill = weight_metric)) +
  geom_boxplot() + 
  theme_bw() +
  ggtitle("Foldwise Approach",
          subtitle = "TCGA - Pattern 3") +
  ylab(used_metric_) +
  xlab("Test-Situations") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24),
        legend.position = "top") +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005) +
  guides(fill = guide_legend(title = "Weight Metric: ")) +
  scale_fill_manual(values = c('darkolivegreen3', "darkorange3", "cyan4"))

# 2-3 Count how often a approach has the highest median!
counter        <- c(0, 0, 0)
names(counter) <- c("No", "F1", "Acc")
for (curr_test in unique(DF_all$Testsituation)) {
  
  m_none <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                   DF_all$weight_metric == 'None'])
  
  m_acc <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                  DF_all$weight_metric == 'Accuracy'])
  
  m_f1s <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                  DF_all$weight_metric == 'F-1 Score'])
  
  all_res_  <- c(m_none, m_f1s, m_acc)
  max_index <- which(all_res_ == max(all_res_, na.rm = TRUE))
  
  if (length(max_index) == 1) counter[max_index] <- counter[max_index] + 1
}


# Compare the different Approaches              --- pattern 3 [F1-Score]    ----
# [1] ----- CC Results
DF_CC <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/CompleteCase_Approach/setting3/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(file_curr, metric = "F1", train_sit = 1)
  DF_CC     <- rbind(DF_CC, curr_df)
}

# 1-3 Add the approach to the DF
DF_CC$Approach <- "Complete-Case"

# [2] ----- SB Results
DF_SB <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/SingleBlock_Approach/setting3/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics_SB(x = file_curr, metric = "F1", train_sit = 1)
  
  # Only keep the results of single 'A'
  curr_df   <- curr_df[curr_df$Learn_Block == 'A',]
  DF_SB     <- rbind(DF_SB, curr_df)
}

# 1-3 Add the approach to the DF
DF_SB$Approach <- "Single-Block"

# [3] ----- IMP Results
DF_IMP <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/Imputation_Approach//setting3/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(x = file_curr, metric = "F1", train_sit = 1)
  DF_IMP    <- rbind(DF_IMP, curr_df)
}

# 1-3 Add the approach to the DF
DF_IMP$Approach <- "Imputation"

# [4] ----- BW Results
DF_BW <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/BlockWise_Approach//setting3/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files[grep("f1", files)]) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(x = file_curr, metric = "F1", train_sit = 1)
  DF_BW     <- rbind(DF_BW, curr_df)
}

# 1-3 Add the approach to the DF
DF_BW$Approach <- "Block-wise"

# [5] ----- FW Results
DF_FW <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/FoldWise_Approach//setting3/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics_FW(file_curr, metric = "F1", train_sit = 1)
  DF_FW     <- rbind(DF_FW, curr_df)
}

# 1-3 Only keep Results with the F1-Score as weight metric
DF_FW <- DF_FW[DF_FW$weight_metric == "F1",]

# 1-4 Add the approach to the DF
DF_FW$Approach <- "Fold-wise"

# [6] ----- Bind the DFs and create a plot
# 6-1 Columns we need from every approach-DF
needed_cols <- c("Testsituation", "Metric", "DF", "performance_metric", "Approach")

# 6-2 Bind the DFs from the different approaches to a single DF for plots!
DF_all <- rbind(DF_CC[,needed_cols], DF_SB[,needed_cols], DF_IMP[,needed_cols],
                DF_BW[,needed_cols], DF_FW[,needed_cols])

# 6-2-1 Relevel the 'Approach' variable
DF_all$Approach <- factor(DF_all$Approach, levels = c("Fold-wise", "Block-wise",
                                                      "Imputation", "Single-Block",
                                                      "Complete-Case"), 
                          ordered = TRUE)

# 6-3 Extract the used metric
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 6-4 Do the plot
ggplot(data = DF_all, aes(x = Testsituation, y = Metric, fill = Approach)) +
  geom_boxplot(position = position_dodge(preserve = "single")) + 
  theme_bw() +
  ggtitle("Comparison of all Approaches",
          subtitle = "TCGA - Pattern 3") +
  ylab(used_metric_) +
  xlab("Test-Situations") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24),
        legend.position = "top") +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005) +
  guides(fill = guide_legend(title = "Approach: ")) +
  scale_fill_manual(values = c("darkolivegreen2", "darkorange1", 'darkorchid1', 'darkgoldenrod1', 'darkslategray3'))

# 6-5 Count how often a approach has the highest/ worst median!
best_counter        <- c(0, 0, 0, 0, 0)
names(best_counter) <- c("CC", "SB", "IMP", "BW", "FW")

worst_counter        <- c(0, 0, 0, 0, 0)
names(worst_counter) <- c("CC", "SB", "IMP", "BW", "FW")
for (curr_test in unique(DF_all$Testsituation)) {
  
  m_cc <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Complete-Case'])
  
  m_sb <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Single-Block'])
  
  m_imp <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                  DF_all$Approach == 'Imputation'])
  
  m_bw <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Block-wise'])
  
  m_fw <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Fold-wise'])
  
  all_res_  <- c(m_cc, m_sb, m_imp, m_bw, m_fw)
  
  max_index <- which(all_res_ == max(all_res_, na.rm = TRUE))
  min_index <- which(all_res_ == min(all_res_, na.rm = TRUE))
  
  if (length(max_index) == 1) best_counter[max_index]  <- best_counter[max_index] + 1
  if (length(min_index) == 1) worst_counter[min_index] <- worst_counter[min_index] + 1
}


# Compare the different Approaches              --- pattern 3 [Balanced Accuracy]   ----
# [1] ----- CC Results
DF_CC <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/CompleteCase_Approach/setting3/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(file_curr, metric = "Balance_Acc", train_sit = 1)
  DF_CC     <- rbind(DF_CC, curr_df)
}

# 1-3 Add the approach to the DF
DF_CC$Approach <- "Complete-Case"

# [2] ----- SB Results
DF_SB <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/SingleBlock_Approach/setting3/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics_SB(x = file_curr, metric = "Balance_Acc", train_sit = 1)
  
  # Only keep the results of single 'A'
  curr_df   <- curr_df[curr_df$Learn_Block == 'A',]
  DF_SB     <- rbind(DF_SB, curr_df)
}

# 1-3 Add the approach to the DF
DF_SB$Approach <- "Single-Block"

# [3] ----- IMP Results
DF_IMP <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/Imputation_Approach//setting3/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(x = file_curr, metric = "Balance_Acc", train_sit = 1)
  DF_IMP    <- rbind(DF_IMP, curr_df)
}

# 1-3 Add the approach to the DF
DF_IMP$Approach <- "Imputation"

# [4] ----- BW Results
DF_BW <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/BlockWise_Approach//setting3/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files[grep("f1", files)]) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(x = file_curr, metric = "Balance_Acc", train_sit = 1)
  DF_BW     <- rbind(DF_BW, curr_df)
}

# 1-3 Add the approach to the DF
DF_BW$Approach <- "Block-wise"

# [5] ----- FW Results
DF_FW <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/FoldWise_Approach//setting3/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics_FW(file_curr, metric = "Balance_Acc", train_sit = 1)
  DF_FW     <- rbind(DF_FW, curr_df)
}

# 1-3 Only keep Results with the F1-Score as weight metric
DF_FW <- DF_FW[DF_FW$weight_metric == "F1",]

# 1-4 Add the approach to the DF
DF_FW$Approach <- "Fold-wise"

# [6] ----- Bind the DFs and create a plot
# 6-1 Columns we need from every approach-DF
needed_cols <- c("Testsituation", "Metric", "DF", "performance_metric", "Approach")

# 6-2 Bind the DFs from the different approaches to a single DF for plots!
DF_all <- rbind(DF_CC[,needed_cols], DF_SB[,needed_cols], DF_IMP[,needed_cols],
                DF_BW[,needed_cols], DF_FW[,needed_cols])

# 6-2-1 Relevel the 'Approach' variable
DF_all$Approach <- factor(DF_all$Approach, levels = c("Fold-wise", "Block-wise",
                                                      "Imputation", "Single-Block",
                                                      "Complete-Case"), 
                          ordered = TRUE)

# 6-3 Extract the used metric
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else if (DF_all$performance_metric[1] == "Balance_Acc") {
  used_metric_ <- "Metric: Balanced Accuracy" 
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 6-4 Do the plot
ggplot(data = DF_all, aes(x = Testsituation, y = Metric, fill = Approach)) +
  geom_boxplot(position = position_dodge(preserve = "single")) + 
  theme_bw() +
  ggtitle("Comparison of all Approaches",
          subtitle = "TCGA - Pattern 3") +
  ylab(used_metric_) +
  xlab("Test-Situations") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24),
        legend.position = "top") +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005) +
  guides(fill = guide_legend(title = "Approach: ")) +
  scale_fill_manual(values = c("darkolivegreen2", "darkorange1", 'darkorchid1', 'darkgoldenrod1', 'darkslategray3'))

# 6-5 Count how often a approach has the highest median!
counter        <- c(0, 0, 0, 0, 0)
names(counter) <- c("CC", "SB", "IMP", "BW", "FW")
for (curr_test in unique(DF_all$Testsituation)) {
  
  m_cc <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Complete-Case'])
  
  m_sb <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Single-Block'])
  
  m_imp <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                  DF_all$Approach == 'Imputation'])
  
  m_bw <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Block-wise'])
  
  m_fw <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Fold-wise'])
  
  all_res_  <- c(m_cc, m_sb, m_imp, m_bw, m_fw)
  max_index <- which(all_res_ == max(all_res_, na.rm = TRUE))
  
  if (length(max_index) == 1) counter[max_index] <- counter[max_index] + 1
}
# Compare the different Approaches              --- pattern 3 [MCC]         ----
# [1] ----- CC Results
DF_CC <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/CompleteCase_Approach/setting3/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(file_curr, metric = "MCC", train_sit = 1)
  DF_CC     <- rbind(DF_CC, curr_df)
}

# 1-3 Add the approach to the DF
DF_CC$Approach <- "Complete-Case"

# [2] ----- SB Results
DF_SB <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/SingleBlock_Approach/setting3/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics_SB(x = file_curr, metric = "MCC", train_sit = 1)
  
  # Only keep the results of single 'A'
  curr_df   <- curr_df[curr_df$Learn_Block == 'A',]
  DF_SB     <- rbind(DF_SB, curr_df)
}

# 1-3 Add the approach to the DF
DF_SB$Approach <- "Single-Block"

# [3] ----- IMP Results
DF_IMP <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/Imputation_Approach//setting3/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(x = file_curr, metric = "MCC", train_sit = 1)
  DF_IMP    <- rbind(DF_IMP, curr_df)
}

# 1-3 Add the approach to the DF
DF_IMP$Approach <- "Imputation"

# [4] ----- BW Results
DF_BW <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/BlockWise_Approach//setting3/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files[grep("f1", files)]) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(x = file_curr, metric = "MCC", train_sit = 1)
  DF_BW     <- rbind(DF_BW, curr_df)
}

# 1-3 Add the approach to the DF
DF_BW$Approach <- "Block-wise"

# [5] ----- FW Results
DF_FW <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/TCGA/FoldWise_Approach//setting3/"
files     <- list.files(data_path)

# 1-2 Loop over the files and get the results
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics_FW(file_curr, metric = "MCC", train_sit = 1)
  DF_FW     <- rbind(DF_FW, curr_df)
}

# 1-3 Only keep Results with the F1-Score as weight metric
DF_FW <- DF_FW[DF_FW$weight_metric == "F1",]

# 1-4 Add the approach to the DF
DF_FW$Approach <- "Fold-wise"

# [6] ----- Bind the DFs and create a plot
# 6-1 Columns we need from every approach-DF
needed_cols <- c("Testsituation", "Metric", "DF", "performance_metric", "Approach")

# 6-2 Bind the DFs from the different approaches to a single DF for plots!
DF_all <- rbind(DF_CC[,needed_cols], DF_SB[,needed_cols], DF_IMP[,needed_cols],
                DF_BW[,needed_cols], DF_FW[,needed_cols])

# 6-2-1 Relevel the 'Approach' variable
DF_all$Approach <- factor(DF_all$Approach, levels = c("Fold-wise", "Block-wise",
                                                      "Imputation", "Single-Block",
                                                      "Complete-Case"), 
                          ordered = TRUE)

# 6-3 Extract the used metric
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else if (DF_all$performance_metric[1] == "Balance_Acc") {
  used_metric_ <- "Metric: Balanced Accuracy" 
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 6-4 Do the plot
ggplot(data = DF_all, aes(x = Testsituation, y = Metric, fill = Approach)) +
  geom_boxplot(position = position_dodge(preserve = "single")) + 
  theme_bw() +
  ggtitle("Comparison of all Approaches",
          subtitle = "TCGA - Pattern 1") +
  ylab(used_metric_) +
  xlab("Test-Situations") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24),
        legend.position = "top") +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005) +
  guides(fill = guide_legend(title = "Approach: ")) +
  scale_fill_manual(values = c("darkolivegreen2", "darkorange1", 'darkorchid1', 'darkgoldenrod1', 'darkslategray3'))

# 6-5 Count how often a approach has the highest median!
counter        <- c(0, 0, 0, 0, 0)
names(counter) <- c("CC", "SB", "IMP", "BW", "FW")
for (curr_test in unique(DF_all$Testsituation)) {
  
  m_cc <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Complete-Case'])
  
  m_sb <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Single-Block'])
  
  m_imp <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                  DF_all$Approach == 'Imputation'])
  
  m_bw <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Block-wise'])
  
  m_fw <- median(DF_all$Metric[DF_all$Testsituation == curr_test & 
                                 DF_all$Approach == 'Fold-wise'])
  
  all_res_  <- c(m_cc, m_sb, m_imp, m_bw, m_fw)
  max_index <- which(all_res_ == max(all_res_, na.rm = TRUE))
  
  if (length(max_index) == 1) counter[max_index] <- counter[max_index] + 1
} 



######################################################################### TEST ----
# Combine the results of all CC approaches [pattern1, .., 3] into a single plot!
#### CC
# [0] Define variables
# 0- 1 Paths to the results of the CV
data_path_1 <- "./docs/CV_Res/TCGA/CompleteCase_Approach/setting1/"
data_path_2 <- "./docs/CV_Res/TCGA/CompleteCase_Approach/setting2/"
data_path_3 <- "./docs/CV_Res/TCGA/CompleteCase_Approach/setting3/"

all_paths <- c(data_path_1, data_path_2, data_path_3)

# 0-2 empty DF to store the results
DF_CC <- data.frame()

# [1] Extract the results
# 1-1 Loop over all folders that contain the CV results for CC!
for (pattern_ in c(1, 2, 3)) {
  
  curr_path <- all_paths[pattern_]
  
  for (curr_file in list.files(curr_path)) {
    
    # Load the result and assign it to 'file_curr'
    file_curr <- load(paste0(curr_path, "/", curr_file))
    file_curr <- eval(as.symbol(file_curr))
    
    curr_df         <- extract_avg_metrics(file_curr, metric = "F1", train_sit = 1)
    curr_df$pattern <- pattern_
    DF_CC           <- rbind(DF_CC, curr_df)
  }
}

# [2] Do the plots and get informations
if (DF_CC$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_CC$performance_metric[1])
}

# 2-2 The plot itself
ggplot(data = DF_CC, aes(x = Testsituation, y = Metric)) +
  geom_boxplot(fill = 'darkolivegreen3') + 
  theme_bw() +
  ggtitle("Complete-Case Approach",
          subtitle = "TCGA - Patterns 1, 2 & 3") +
  ylab(used_metric_) +
  xlab("Test-Situations") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24)) +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005) +
  facet_grid(pattern ~ .)

# 2-3 Get summarys to the perfornance!
for (patt_ in c(1, 2, 3)) {
  
  print(paste0("Current Pattern: ", patt_, " --------------------------------"))
  
  res_ <- sapply(unique(DF_CC$Testsituation), FUN = function(x)  {
    
    DF_curr <- DF_CC[DF_CC$pattern == patt_,]
    DF_curr <- DF_curr[DF_curr$Testsituation == x,]
    
    if (all(is.na(DF_curr$Metric))) {
      c(0, 0, 0, 0, 0)
    } else {
      as.numeric(summary(DF_curr$Metric[DF_curr$Testsituation == x])) 
    }
  })
  
  
  names(res_) <- unique(DF_CC$Testsituation)

  print(res_)
}




# Analyse Results of the complete case Approach --- pattern 3               ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/TCGA/CompleteCase_Approach/setting3/"

# [1] Load all Results w/ Romans Approach
# 1-1 get all files in the folder, that are singleblock performances
files <- list.files(data_path)

# 1-2 Loop over all the files and extract the results
DF_all <- data.frame()
for (curr_file in files) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", curr_file))
  file_curr <- eval(as.symbol(file_curr))
  
  curr_df   <- extract_avg_metrics(file_curr, metric = "F1", train_sit = 1)
  DF_all    <- rbind(DF_all, curr_df)
}

# [2] Plot the Results
# 2-1-1 The used metric for the comparison of the performance
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 2-2 The plot itself
ggplot(data = DF_all, aes(x = Testsituation, y = Metric)) +
  geom_boxplot(fill = 'darkolivegreen3') + 
  theme_bw() +
  ggtitle("Complete-Case Approach",
          subtitle = "TCGA - Pattern 3") +
  ylab(used_metric_) +
  xlab("Test-Situations") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24)) +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005)

# 2-3 Get summarys to the perfornance!
res_ <- sapply(unique(DF_all$Testsituation), FUN = function(x){
  summary(DF_all$Metric[DF_all$Testsituation == x])
})
names(res_) <- unique(DF_all$Testsituation)
