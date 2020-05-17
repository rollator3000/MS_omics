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
    >> Only for complete cases & imputation approaches!
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
  met_meas <- unique(names(x[[1]]))
  
  if (!(metric %in% met_meas)) {
    stop("'metric' is not within these!")
  }
  
  # [1] Unpack the lists and put them into a single regular DF
  # 1-1 Extract the metric 
  curr_res <-  unlist(c(x[[1]][metric], x[[2]][metric], x[[3]][metric],
                        x[[4]][metric], x[[5]][metric]))
  
  
  DF_final <- data.frame("Fold"      = 1:5,
                         "Metric"    = curr_res,
                         "performance_metric" = metric)
  
  # [2] Return 'DF_final'
  return(DF_final)
}
extract_metrics_FW_BW <- function(x, metric) {
  "
    >> Only Block-/ FoldWise Approach
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
                         "weight_metric" = c(rep("None", 5),
                                         rep("F1-Score", 5),
                                         rep("Accuracy", 5)),
                         "Metric"    = c(no_weight_res,
                                         f1_weight_res,
                                         acc_weight_res),
                         "performance_metric" = metric)
  
  # [2] Return 'DF_final'
  return(DF_final)
}
extract_metrics_SBAPP <- function(x, metric) {
  "
    >> Only for single_block approach!
  Function to convert a list 'x' w/ the 6 entrances - one entrance for the results 
  of each feature-block in train in REAL.
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
  assert_list(x, len = 6)
  
  # 0-2 'metric' must be in 'x'
  assert_string(metric)
  
  # extract all used  metrics and check whether 'metric' exist!
  met_meas <- unique(names(x$df1[[1]]))
  
  if (!(metric %in% met_meas)) {
    stop("'metric' is not within these!")
  }
  
  # [1] Unpack the lists for each single block and put them into a single regular DF
  # 1-1 Extract the metric 'df1'
  DF1_res <-  unlist(c(x$df1[[1]][metric], x$df1[[2]][metric], x$df1[[3]][metric],
                       x$df1[[4]][metric], x$df1[[5]][metric]))
  
  # 1-2 Extract the metric 'df2'
  DF2_res <-  unlist(c(x$df2[[1]][metric], x$df2[[2]][metric], x$df2[[3]][metric],
                       x$df2[[4]][metric], x$df2[[5]][metric]))
  
  # 1-3 Extract the metric 'df3'
  DF3_res <-  unlist(c(x$df3[[1]][metric], x$df3[[2]][metric], x$df3[[3]][metric],
                       x$df3[[4]][metric], x$df3[[5]][metric]))
  
  # 1-4 Extract the metric 'df4'
  DF4_res <-  unlist(c(x$df4[[1]][metric], x$df4[[2]][metric], x$df4[[3]][metric],
                       x$df4[[4]][metric], x$df4[[5]][metric]))
  
  # 1-2 Extract the metric 'df2'
  DF51_res <-  unlist(c(x$df51[[1]][metric], x$df51[[2]][metric], x$df51[[3]][metric],
                        x$df51[[4]][metric], x$df51[[5]][metric]))
  
  # 1-2 Extract the metric 'df2'
  DF53_res <-  unlist(c(x$df53[[1]][metric], x$df53[[2]][metric], x$df53[[3]][metric],
                        x$df53[[4]][metric], x$df53[[5]][metric]))
  
  
  DF_final <- data.frame("Fold"      = c(1:5, 1:5, 1:5, 1:5, 1:5, 1:5),
                         "Metric"    = c(DF1_res, DF2_res, DF3_res, DF4_res,
                                         DF51_res, DF53_res),
                         "performance_metric" = metric,
                         "Feature_Block" = c(rep("df1", times = 5),
                                             rep("df2", times = 5),
                                             rep("df3", times = 5),
                                             rep("df4", times = 5),
                                             rep("df51", times = 5),
                                             rep("df53", times = 5)))
  
  # [2] Return 'DF_final'
  return(DF_final)
}

# Analyse Results of the FoldWise_Approach                                   ----
# [0] Select the file with the results & load it as 'file_curr'
curr_file <- "./docs/CV_Res/REAL/FoldWise_Approach.R"
file_curr <- load(curr_file)
file_curr <- eval(as.symbol(file_curr))

# [1] Extract the metrics of the CV
DF_all <- extract_metrics_FW_BW(x = file_curr, metric = "F1")

# [2] Plot the Results
# 2-1 Extract the needed Information from 'DF_all'
# 2-1-1 The used metric for the comparison of the performance
used_metric_ <- paste("Metric:", DF_all$performance_metric[1])

# 2-2 Do the plot
ggplot(data = DF_all, aes(x = weight_metric, y = Metric)) +
  geom_boxplot(fill = '#F8766D') + 
  theme_bw() +
  ylab(used_metric_) +
  xlab("Weighting of the foldwise predictions") +
  ggtitle("Foldwise Approach - CV Results") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 18)) 


# Analyse Results of the BlockWise_Approach                                  ----
# [0] Select the file with the results & load it as 'file_curr'
curr_file <- "./docs/CV_Res/REAL/BlockWise_Approach.R"
file_curr <- load(curr_file)
file_curr <- eval(as.symbol(file_curr))

# [1] Extract the metrics of the CV
DF_all <- extract_metrics_FW_BW(x = file_curr, metric = "F1")

# [2] Plot the Results
# 2-1 Extract the needed Information from 'DF_all'
# 2-1-1 The used metric for the comparison of the performance
used_metric_ <- paste("Metric:", DF_all$performance_metric[1])

# 2-2 Do the plot
ggplot(data = DF_all, aes(x = weight_metric, y = Metric)) +
  geom_boxplot(fill = '#F8766D') + 
  theme_bw() +
  ylab(used_metric_) +
  xlab("Weighting of the blockwise predictions") +
  ggtitle("Blockwise Approach - CV Results") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 18))


# Analyse Results of the Imputation                                          ----
# [0] Select the file with the results & load it as 'file_curr'
curr_file <- "./docs/CV_Res/REAL/Imputation_Approach.R"
file_curr <- load(curr_file)
file_curr <- eval(as.symbol(file_curr))

# [1] Extract the metrics of the CV
DF_all <- extract_metrics(x = file_curr, metric = "F1")

# [2] Plot the Results
# 2-1 Extract the needed Information from 'DF_all'
# 2-1-1 The used metric for the comparison of the performance
used_metric_ <- paste("Metric:", DF_all$performance_metric[1])

# 2-2 Do the plot
ggplot(data = DF_all, aes(x = 1, y = Metric)) +
  geom_boxplot(fill = '#F8766D') + 
  theme_bw() +
  ylab(used_metric_) +
  xlab("") +
  ggtitle("Imputation Approach - CV Results") +
  theme(text = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()) 


# Analyse Results of the CC Approach                                         ----
# [0] Select the file with the results & load it as 'file_curr'
curr_file <- "./docs/CV_Res/REAL/CC_Approach.R"
file_curr <- load(curr_file)
file_curr <- eval(as.symbol(file_curr))

# [1] Extract the metrics of the CV
DF_all <- extract_metrics(x = file_curr, metric = "F1")

# [2] Plot the Results
# 2-1 Extract the needed Information from 'DF_all'
# 2-1-1 The used metric for the comparison of the performance
used_metric_ <- paste("Metric:", DF_all$performance_metric[1])

# 2-2 Do the plot
ggplot(data = DF_all, aes(x = 1, y = Metric)) +
  geom_boxplot(fill = '#F8766D') + 
  theme_bw() +
  ylab(used_metric_) +
  xlab("") +
  ggtitle("Complete Cases Approach - CV Results") +
  theme(text = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())

# Analyse Results of the Single Block Approach                               ----
# [0] Select the file with the results & load it as 'file_curr'
curr_file <- "./docs/CV_Res/REAL/SingleBlockApproach.R"
file_curr <- load(curr_file)
file_curr <- eval(as.symbol (file_curr))

# [1] Extract the metrics of the CV
DF_all <- extract_metrics_SBAPP(x = file_curr, metric = "F1")

# [2] Plot the Results
# 2-1 Extract the needed Information from 'DF_all'
# 2-1-1 The used metric for the comparison of the performance
used_metric_ <- paste("Metric:", DF_all$performance_metric[1])

# 2-2 Do the plot
ggplot(data = DF_all, aes(x = Feature_Block, y = Metric)) +
  geom_boxplot(fill = '#F8766D') + 
  theme_bw() +
  ylab(used_metric_) +
  xlab("FeatureBlock for RF") +
  ggtitle("Single Block Approach - CV Results") +
  theme(text = element_text(size = 18))
