"
Script to visualize the different Results from the CVs
"
library(ggplot2)
require(gridExtra)
library(reshape2)
library(checkmate)

# Used for all Results except from 'Romans Approach'
extract_metrics <- function(x, metric, train_sit) {
  "
  Function to convert a list 'x' - w/ the 2 entrances 'res_all' & 'settings' - 
  to a DF, that we can use for plotting!
  'x' is the result of a CV of one trainsetting & all has results for all of 
  the 20 testsettings! From this list we extract the metric for all of the 
  different trainsettings!
  
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
  #       & fill 'DF_final'
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
  
  if (is.null(x_settings[["weighted"]]) || !x_settings[["weighted"]]) {
    DF_final$weight_metric <- NA
  } else {
    ifelse(!x_settings$weighted,
           DF_final$weight_metric <- NA,
           DF_final$weight_metric <- x_settings$weight_metric)
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
  
  if (is.null(x_settings[["weighted"]] || !x_settings[["weighted"]])) {
    DF_final$weight_metric <- NA
  } else {
    ifelse(!x_settings$weighted,
           DF_final$weight_metric <- NA,
           DF_final$weight_metric <- x_settings$weight_metric)
  }
  
  # [3] Return 'DF_final'
  return(DF_final)
}


# Used to extract Results from Romans Approach 
# - as these have a slightly different structure due to 
#   performance realted changes
extract_metrics_FW <- function(x, metric, train_sit) {
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

# Analyse the explorative singleblock single block Results                  ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/gender/explorative_subsets"

# [1] Load all explorative singleblock Results
# 1-1 get all files in the folder, that are singleblock performances
files <- list.files(data_path)
files <- files[grepl("single", files)]

# 1-2 Load the data and store it into a single DF!
DF_all <- data.frame()
for (curr_file in files) {
  DF_curr <-  read.csv2(paste0(data_path, "/", curr_file), stringsAsFactors = F)
  DF_all  <- rbind(DF_all, DF_curr)
}

# 1-3 Convert features to right data type!
str(DF_all)
num_cols <- c("OOB_Acc", "Test_Acc", "Test_F1", "Fraction")
DF_all[,num_cols] <- sapply(num_cols, function(x) as.numeric(DF_all[,x]))

# 1-4 Reshape the layout of data for the plot! 
plot_df <- melt(DF_all, id.vars = c("Block", "Fraction", "subset_seed"), 
                measure.vars = c("Test_Acc", "Test_F1", "Block"))
plot_df <- plot_df[plot_df$variable %in% c("Test_Acc", "Test_F1"),]
plot_df$value <- as.numeric(plot_df$value)

# 1-5 Plot the performance 
#     Split by the fraction we've used to subset the single blocks!
ggplot(data = plot_df, aes(x = Block, y = value, fill = variable)) + 
  geom_boxplot() + 
  facet_grid(subset_seed ~ Fraction) +
  theme_bw() +
  ggtitle("Single Block Performance on all 14 DFs [w/ 5 fold CV]",
          subtitle = "split by the amount of features we kept") +
  xlab("Blocks used as Feature Space") +
  ylab("F1 Score // Accuracy") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 17))

# Analyse the explorative joint block block Results                         ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/gender/explorative_subsets"

# [1] Load all explorative singleblock Results
# 1-1 get all files in the folder, that are singleblock performances
files <- list.files(data_path)
files <- files[grepl("joint", files)]

# 1-2 Add the amount of subsets to each block:
DF_all <- data.frame()
for (curr_file in files) {
  DF_curr <-  read.csv2(paste0(data_path, "/", curr_file), stringsAsFactors = F)
  DF_all  <- rbind(DF_all, DF_curr)
}

# 1-4 reshape DF_all for the plot!
plot_df <- melt(DF_all, id.vars = c("Data", "Fraction"), measure.vars = c("Test_Acc", "Test_F1", "Fold"))
plot_df <- plot_df[plot_df$variable %in% c("Test_Acc", "Test_F1"),]
plot_df$value <- as.numeric(plot_df$value)

# 1-5 Add the fraction used of the single blocks as meta data!
for (i in seq_len(nrow(plot_df))) {
  plot_df$mirna_subset[i]    <- strsplit(strsplit(plot_df$Fraction[i], split = "mirna_")[[1]][2], split = "__")[[1]][1]
  plot_df$mutation_subset[i] <- strsplit(plot_df$Fraction[i], split = "mutation_")[[1]][2]
  plot_df$rna_subset[i]      <- strsplit(strsplit(plot_df$Fraction[i], split = "_rna_")[[1]][2], split = "_")[[1]][1]
  plot_df$cnv_subset[i]      <- strsplit(strsplit(plot_df$Fraction[i], split = "_cnv_")[[1]][2], split = "_")[[1]][1]
  plot_df$Fraction_new[i]    <- strsplit(plot_df$Fraction[i], split = "__mirna_")[[1]][2]
  plot_df$Fraction_new[i]    <- paste0("mirna_", plot_df$Fraction_new[i])
}

plot_df$rna_subset_plot <- sapply(plot_df$rna_subset, function(x) paste0("rna_", x)) 
plot_df$cnv_subset_plot <- sapply(plot_df$cnv_subset, function(x) paste0("cnv_", x)) 

# [2] Do the plot, split by subsets 
ggplot(data = plot_df, aes(x = Fraction_new , y = value, fill = variable)) +
  geom_boxplot() + 
  theme_bw() +
  facet_grid(rna_subset_plot ~ cnv_subset_plot) +
  ggtitle("EXPLORATIVE - Joint Block Performance on all 14 DFs - w/ seed 12345", 
          subtitle = "Single Blocks were subsetted! y-axis RNA splits || x-axis CNV splits") +
  xlab("subsets for mirna & mutation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 17))

# Analyse the fix subsettet DFs w/ joint blocks!                            ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/gender/performance_final_subsets"

# [1] Load all explorative singleblock Results
# 1-1 get all files in the folder, that are singleblock performances
files <- list.files(data_path)
files <- files[grepl("joint", files)]

# 1-2 Add the amount of subsets to each block:
DF <- read.csv2(paste0(data_path, "/", files[1]), stringsAsFactors = F)

# 1-3 Convert numeric columns to numeric
numeric_cols          <- c("OOB_Acc", "Test_Acc", "Test_F1")
DF[,numeric_cols] <- sapply(numeric_cols, 
                            function(x) as.numeric(DF[,x]))

# 1-4 reshape DF_all for the plot!
plot_df <- melt(DF, id.vars = c("Data", "fold"), measure.vars = c("Test_Acc", "Test_F1"))

# [2] Do the plots
# 2-1 Split by the different Seeds!
ggplot(data = plot_df, aes(x = Data, y = value, fill = variable)) +
  geom_boxplot() + 
  theme_bw() +
  ggtitle("Joint Block Performance on the 14 fixed DF w/ 5 fold CV",
          subtitle = "50% mirna, 10% mutation, 15% rna & 1.25% cnv") +
  xlab("Dataset") +
  ylab("Metric [Acc & F1]") +
  ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15))

# Analyse the fix subsettet DFs w/ single blocks!                           ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/gender/performance_final_subsets"

# [1] Load all explorative singleblock Results
# 1-1 get all files in the folder, that are singleblock performances
files <- list.files(data_path)
files <- files[grepl("single", files)]

# 1-2 Add the amount of subsets to each block:
DF <- read.csv2(paste0(data_path, "/", files[1]), stringsAsFactors = F)

# 1-3 Convert the columns to the right format
num_cols <- c("OOB_Acc", "Test_Acc", "Test_F1")
DF[,num_cols] <- sapply(num_cols, 
                        function(x) as.numeric(DF[,x]))

# 1-4 Make the name of the single DFs nicer
DF$Data <- sapply(DF$Data, function(x) strsplit(x, split = "_subset")[[1]][1])

# 1-4 reshape the layout of data for the plot! 
plot_df <- melt(DF, id.vars = c("Data", "Block"), 
                measure.vars = c("Test_Acc", "Test_F1"))

# [2] Do the plot
  ggplot(data = plot_df, aes(x = Data, y = value, fill = variable)) + 
    geom_boxplot() + 
    facet_grid(. ~ Block) +
    theme_bw() +
    ggtitle("Single Block Performance on the 14 fixed DFs w/ 5 fold CV",
            subtitle = "Final-Subset: 10% mirna & mutation, 5% rna & 2.5% cnv \n--- split by the seeds used to subset the featurespace") +
    xlab("Blocks used as Feature Space") +
    ylab("Metric [Acc & F1]") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 15))

# Analyse Results of Romans Approach on the fixed DFs                       ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/gender/Roman_final_subsets/setting2"

# [1] Load all Results w/ Romans Approach
# 1-1 List all files from the 'data_path'
files <- list.files(data_path)
files <- files[grepl("NEW", files)]

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
ggplot(data = DF_all, aes(x = Testsituation, y = Metric, fill = weight_metric)) +
  geom_boxplot() + 
  theme_bw() +
  ggtitle("Romans Method applied to a data subset",
          subtitle = "split by the weighting used for the single folds!") +
  xlab("TestSituations") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15)) +
  geom_vline(xintercept = seq(from = 1.5, to = 19.5, by = 1),
             col = "red", lty = 2) 

# Analyse Results of Norberts Approach on the fixed DFs                     ----
data_path <- "./docs/CV_Res/gender/Norbert_final_subsets/setting2"

# [1] Load all Results w/ Romans Approach
# 1-1 List all files from the 'data_path'
files <- list.files(data_path)

# 1-2 Assign Files to the different approaches
f1_weighted  <- files[c(grep("f1", files), grep("F1", files))]
acc_weighted <- files[grep("acc", files)]
no_weighted  <- files[!grepl("acc", files) & !grepl("f1", files)  & !grepl("F1", files)]

# 1-3 Extract the results Load the Results for all TrainSettings
# 1-3-1 F1-Weighting
f1_res <- data.frame()
for (file_ in f1_weighted) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", file_))
  file_curr <- eval(as.symbol(file_curr))
  
  # Extract the Metrics for the current file
  curr_df   <- extract_avg_metrics(file_curr, metric = "F1", train_sit = 1)
  
  # Merge curr_df into the f1_res
  f1_res <- rbind(f1_res, curr_df)
}

# 1-3-2 Acc-Weighting
acc_res <- data.frame()
for (file_ in acc_weighted) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", file_))
  file_curr <- eval(as.symbol(file_curr))
  
  # Extract the Metrics for the current file
  curr_df   <- extract_avg_metrics(file_curr, metric = "F1", train_sit = 1)
  
  # Merge curr_df into the f1_res
  acc_res <- rbind(acc_res, curr_df)
}

# 1-3-3 No weighting
no_res <- data.frame()
for (file_ in no_weighted) {
  
  # Load the result and assign it to 'file_curr'
  file_curr <- load(paste0(data_path, "/", file_))
  file_curr <- eval(as.symbol(file_curr))
  
  # Extract the Metrics for the current file
  curr_df   <- extract_avg_metrics(file_curr, metric = "F1", train_sit = 1)
  
  # Merge curr_df into the f1_res
  no_res <- rbind(no_res, curr_df)
}

# 1-4 Bind all DFs to a single DF
DF_all <- rbind(f1_res, acc_res, no_res)

# [2] Plot the Results
ggplot(data = DF_all, aes(x = Testsituation, y = Metric, fill = weight_metric)) +
  geom_boxplot() + 
  theme_bw() +
  ggtitle("Norberts Method applied to a data subset",
          subtitle = "split by the weighting used for the single folds!") +
  xlab("TestSituations") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15)) +
  geom_vline(xintercept = seq(from = 1.5, to = 19.5, by = 1),
             col = "red", lty = 2) 


# Analyse Results of the complete case Approach on the fixed DFs            ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/gender/complete_cases/setting2/"

# [1] Load all Results w/ Romans Approach
# 1-1 get all files in the folder, that are singleblock performances
files <- list.files(data_path)

# 1-2 Load the Results for all TrainSettings
BLCA <- load(paste0(data_path, "/", files[1]))
BLCA <- eval(as.symbol(BLCA))

# 1-3 Extract the Metrics and TestSetting as desired:
# BLCA
res_blca <- extract_metrics(x = BLCA, metric = "F1", train_sit = 1)

ggplot(data = res_blca, aes(x = Testsituation, y = Metric)) +
  geom_boxplot() + 
  theme_bw() +
  ggtitle("Complete Cases Approach applied to a data subset") +
  xlab("TestSituations") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15))