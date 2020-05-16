"
Script to visualize the different Results from the REAL CVs
"
# [0] Load Packages and define functions
library(ggplot2)
require(gridExtra)
library(reshape2)
library(checkmate)

extract_metrics       <- function(x, metric) {
  "
  Function to convert a list 'x' w/ the 5 entrances - one entrance full of 
  metrics for each fold of the CV - into a DF suitable for plotting!
  'x' is the result of a CV on the real data. 
  
   Args:
    x (list)           : list filled with the metrics for every single 
                         Itteration in CV + the settings it was created from!
                         Needs the names: 'res_all' & 'settings'
    metric (str)       : the metric we shall extraxt from the result list!
                         Needs to be in:
                         ['Accuracy', 'Sensitifity', 'Specificity', 'Precision', 
                          'Recall', 'F1', 'Balance_Acc', 'AUC', 'MCC']

   Return:
    DF, suited for plotting, with all metric for all testsettings in 'x'
  "
  # [0] Check Inputs 
  # 0-1 'x' must be a list of length 2 w/ names 'res_all' & 'settings'
  assert_list(x, len = 5)
  
  # 0-2 'metric' must be in 'x'
  assert_string(metric)
  
  # extract all used  metrics and check whether 'metric' exist!
  met_meas <- unique(
    unlist(sapply(names(x), function(j) {
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
extract_metrics_FW_BW <- function(x, metric) {
  "
  Function to convert a list 'x' w/ the 5 entrances - one entrance full of 
  metrics for each fold of the CV - into a DF suitable for plotting!
  'x' is the result of a CV on the real data, when applying the fold-wise/ 
  blockwise approach!!!
  
   Args:
    x (list)           : list filled with the metrics for every single 
                         Itteration in CV + the settings it was created from!
                         Needs the names: 'res_all' & 'settings'
    metric (str)       : the metric we shall extraxt from the result list!
                         Needs to be in:
                         ['Accuracy', 'Sensitifity', 'Specificity', 'Precision', 
                          'Recall', 'F1', 'Balance_Acc', 'AUC', 'MCC']

   Return:
    DF, suited for plotting, with all metric for all testsettings in 'x'
  "
  # [0] Check Inputs 
  # 0-1 'x' must be a list of length 2 w/ names 'res_all' & 'settings'
  assert_list(x, len = 5)
  
  # 0-2 'metric' must be in 'x'
  assert_string(metric)
  
  # extract all used  metrics and check whether 'metric' exist!
  met_meas <- unique(names(x[[1]]$no_weight))
  
  if (!(metric %in% met_meas)) {
    stop("'metric' is not within these!")
  }
  
  # [1] Unpack the lists and put them into a single regular DF
  # 1-1 Extract the metrics for each fold of the CV
  # 1-1-1 No weight when ensembling the block/fold-wise predictions 
  no_weight_res <- unlist(c(x[[1]][["no_weight"]][metric],
                            x[[2]][["no_weight"]][metric],
                            x[[3]][["no_weight"]][metric],
                            x[[4]][["no_weight"]][metric],
                            x[[5]][["no_weight"]][metric]))
  
  # 1-1-2 Acc weight when ensembling the block/fold-wise predictions 
  acc_weight_res <- unlist(c(x[[1]][["acc_weight"]][metric],
                             x[[2]][["acc_weight"]][metric],
                             x[[3]][["acc_weight"]][metric],
                             x[[4]][["acc_weight"]][metric],
                             x[[5]][["acc_weight"]][metric]))
  
  # 1-1-3 F1 weight when ensembling the block/fold-wise predictions 
  f1_weight_res <- unlist(c(x[[1]][["f1_weight"]][metric],
                            x[[2]][["f1_weight"]][metric],
                            x[[3]][["f1_weight"]][metric],
                            x[[4]][["f1_weight"]][metric],
                            x[[5]][["f1_weight"]][metric]))
  
  # 1-2
  DF_final <- data.frame("Fold"      = c(1:5, 1:5, 1:5),
                         "weighting" = c(rep("None", 5),
                                         rep("F1-Score", 5),
                                         rep("Accuracy", 5)),
                         "Weight_Metric"    = c(no_weight_res,
                                                f1_weight_res,
                                                acc_weight_res),
                         "Metric" = metric)
  
  # [2] Return 'DF_final'
  return(DF_final)
}

# # Analyse Results of the FoldWise_Approach                                  ----
# [0] Select the file with the results
curr_file <- "./docs/CV_Res/REAL/FoldWise_Approach.R"

# Load the result and assign it to 'file_curr'
file_curr <- load(curr_file)
file_curr <- eval(as.symbol(file_curr))

extract_metrics_FW(x  =file_curr, metric = "F1")



# [2] Plot the Results
# 2-1 Extract the needed Information from 'DF_all'
# 2-1-1 The used metric for the comparison of the performance
used_metric_ <- paste("Metric:", DF_all$performance_metric[1])

# 2-1-2 Rename the Weight Metrics [--> clearer names]
DF_all$weight_metric[DF_all$weight_metric == "F1"]  <- "F1-Score"
DF_all$weight_metric[DF_all$weight_metric == "Acc"] <- "Accuracy"
DF_all$weight_metric[DF_all$weight_metric == "No"]  <- "None"

# 2-2 Do the plot
ggplot(data = DF_all, aes(x = Testsituation, y = Metric, fill = weight_metric)) +
  geom_boxplot() + 
  theme_bw() +
  ylab(used_metric_) +
  xlab("Different Test-Situations") +
  ggtitle("Foldwise Approach - CV Results",
          subtitle = "split by the weighting of the foldwise predictions!") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 18)) +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 13.5,
                            14.5, 16.5, 17.5, 18.5, 19.5),
             col = "darkgreen", lty = 2) +
  geom_vline(xintercept = c(1.5, 5.5, 11.5, 15.5),
             col = "red", lty = 2, lwd = 1.005) +
  guides(fill = guide_legend(title = "Weight Metric"))

