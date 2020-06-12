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
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 2-2 Do the plot
ggplot(data = DF_all, aes(x = 1, y = Metric)) +
  geom_boxplot(fill = 'darkolivegreen3') + 
  theme_bw() +
  ylab(used_metric_) +
  xlab("") +
  ggtitle("Complete-Case Approach",
          subtitle = "Clinical asthma data") +
  theme(text = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())

# 2-3 Get a summary of the resulting Metrics
summary(DF_all$Metric)


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
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 2-1-2 Rename the feature-blocks
DF_all$Feature_Block <- as.character(DF_all$Feature_Block)
DF_all$Feature_Block[DF_all$Feature_Block == 'df1']  <- "Questionaire"
DF_all$Feature_Block[DF_all$Feature_Block == 'df2']  <- "Clinical routine diagnostics"
DF_all$Feature_Block[DF_all$Feature_Block == 'df3']  <- "Allergen sensitization"
DF_all$Feature_Block[DF_all$Feature_Block == 'df4']  <- "Cytokine expression data"
DF_all$Feature_Block[DF_all$Feature_Block == 'df51'] <- "Gene expression data I"
DF_all$Feature_Block[DF_all$Feature_Block == 'df53'] <- "Gene expression data II"
DF_all$Feature_Block <- factor(DF_all$Feature_Block,
                               levels = c("Questionaire", "Clinical routine diagnostics", 
                                          "Allergen sensitization", "Cytokine expression data",
                                          "Gene expression data I", "Gene expression data II"))
                               
# 2-2 Do the plot
ggplot(data = DF_all, aes(x = Feature_Block, y = Metric)) +
  geom_boxplot(fill = 'darkolivegreen3') + 
  theme_bw() +
  ylab(used_metric_) +
  xlab("Feature-block") +
  ggtitle("Single-Block Approach",
          subtitle = "Clinical asthma data") +
  theme(text = element_text(size = 24),
        axis.text.x = element_text(angle = 25, hjust = 1))

# 2-3 Get summary for each single-block
res_curr <- sapply(unique(DF_all$Feature_Block), FUN = function(x_) {
  summary(DF_all$Metric[DF_all$Feature_Block == x_])
})
colnames(res_curr) <- unique(DF_all$Feature_Block)

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
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 2-2 Do the plot
ggplot(data = DF_all, aes(x = 1, y = Metric)) +
  geom_boxplot(fill = 'darkolivegreen3') + 
  theme_bw() +
  ylab(used_metric_) +
  xlab("") +
  ggtitle("Imputation Approach",
          subtitle = "Clinical asthma data") +
  theme(text = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())


# 2-3 Get a summary of the resulting Metrics
summary(DF_all$Metric)

# Analyse Results of the BlockWise_Approach                                  ----
# [0] Select the file with the results & load it as 'file_curr'
curr_file <- "./docs/CV_Res/REAL/BlockWise_Approach.R"
file_curr <- load(curr_file)
file_curr <- eval(as.symbol(file_curr))

# [1] Extract the metrics of the CV
DF_all <- extract_metrics_FW_BW(x = file_curr, metric = "F1")

# 1-1 Order and name factor levels correctly!
DF_all$weight_metric <- as.character(DF_all$weight_metric)
DF_all$weight_metric[DF_all$weight_metric == "F1-Score"] <- "F-1 Score"
DF_all$weight_metric <- factor(DF_all$weight_metric,
                               levels = c('None', 'Accuracy', 'F-1 Score'))

# [2] Plot the Results
# 2-1 Extract the needed Information from 'DF_all'
# 2-1-1 The used metric for the comparison of the performance
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 2-2 Do the plot
ggplot(data = DF_all, aes(x = weight_metric, y = Metric)) +
  geom_boxplot(fill = 'darkolivegreen3') + 
  theme_bw() +
  ylab(used_metric_) +
  xlab("Weight Metric for the block-wise predictions") +
  ggtitle("Block-wise Approach",
          subtitle = "Clinical asthma data") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        text = element_text(size = 24))

# 2-3 Get summary for each weighting metric used for the ensemble
res_curr <- sapply(unique(DF_all$weight_metric), FUN = function(x_) {
  summary(DF_all$Metric[DF_all$weight_metric == x_])
})
colnames(res_curr) <- unique(DF_all$weight_metric)


# Analyse Results of the FoldWise_Approach                                   ----
# [0] Select the file with the results & load it as 'file_curr'
curr_file <- "./docs/CV_Res/REAL/FoldWise_Approach.R"
file_curr <- load(curr_file)
file_curr <- eval(as.symbol(file_curr))

# [1] Extract the metrics of the CV
DF_all <- extract_metrics_FW_BW(x = file_curr, metric = "F1")

# 1-1 Order and name factor levels correctly!
DF_all$weight_metric <- as.character(DF_all$weight_metric)
DF_all$weight_metric[DF_all$weight_metric == "F1-Score"] <- "F-1 Score"
DF_all$weight_metric <- factor(DF_all$weight_metric,
                               levels = c('None', 'Accuracy', 'F-1 Score'))

# [2] Plot the Results
# 2-1 Extract the needed Information from 'DF_all'
# 2-1-1 The used metric for the comparison of the performance
if (DF_all$performance_metric[1] == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# 2-2 Do the plot
ggplot(data = DF_all, aes(x = weight_metric, y = Metric)) +
  geom_boxplot(fill = 'darkolivegreen3') + 
  theme_bw() +
  ylab(used_metric_) +
  xlab("Weight Metric for the fold-wise predictions") +
  ggtitle("Fold-wise Approach",
          subtitle = "Clinical asthma data") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        text = element_text(size = 24)) 

# 2-3 Get summary for each weighting metric used for the ensemble
res_curr <- sapply(unique(DF_all$weight_metric), FUN = function(x_) {
  summary(DF_all$Metric[DF_all$weight_metric == x_])
})
colnames(res_curr) <- unique(DF_all$weight_metric)
# Analyse the approaches of Hagenberg                                        ----
# All metrics were read out and calculated in 'REAL_07_get_metrics_for_mDD....R'
# Hagenberg --- Setting 5_3_1                                                ----
"Pririty of blocks is assigned by the amount of missing values! "

# [1] Load the Metrics of the CV w/ Hagenbergs Approaches
load("./docs/CV_Res/REAL/Hagenberg_5_3_1.RData")

# [2] Select a metric from:  ["Accuracy", "Kappa", "Sensitifity", "Specificity", 
#                             "Precision", "Recall", "F1", "Balance_Acc", 
#                             "Pos_Pred_Value", "Neg_Pred_Value", "Prevalence", 
#                             "AUC1", "AUC2", "MCC"]
metric__ = "F1"

if (metric__ == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# [3] Create a DF to plot the results!
plot_df <- data.frame("method" = c(rep('PL - ignore, zero', 5),
                                   rep('PL - ignore, intercept', 5),
                                   rep('PL - impute, maximise blocks', 5),
                                   rep('PL - impute, maximise n', 5),
                                   rep("mdd-sPLS", 5)),
                      "fold"   = rep(1:5, 5),
                      "Metric" = c(unlist(all_res$`ignore, zero`[metric__,]),
                                   unlist(all_res$`ignore, intercept`[metric__,]),
                                   unlist(all_res$`impute maximise blocks`[metric__,]),
                                   unlist(all_res$`impute, maximise n`[metric__,]),
                                   unlist(all_res$mdd_sPLS[metric__,])),
                      "used_Metric" = metric__)


# [4] Do the plot
ggplot(data = plot_df, aes(x = method, y = Metric)) +
  geom_boxplot(fill = 'darkolivegreen3') + 
  theme_bw() +
  ylab(used_metric_) +
  xlab("Different Approaches") +
  ggtitle("Priority-Lasso adaptions & mdd-sPLS method",
          subtitle = "Clinical asthma data") +
  theme(text = element_text(size = 24),
        axis.text.x = element_text(angle = 25, hjust = 1))

# 4-1 Get a summary for the different approaches
sapply(levels(plot_df$method), FUN = function(x) summary(plot_df$Metric[plot_df$method == x]))

# Hagenberg --- Setting 5_3_2                                                ----
"Different Block Combinations for the prediciton [Train on all blocks - but only
 use block YZ & QX for the prediciton, etc.
 Two further distinctions for PL approach - 'all_blocks' & 'pred_blocks' "

# [1] Load the Metrics of the CV w/ Hagenbergs Approaches
load("./docs/CV_Res/REAL/Hagenberg_5_3_2.RData")

# [2] Select a metric from:  ["Accuracy", "Kappa", "Sensitifity", "Specificity", 
#                             "Precision", "Recall", "F1", "Balance_Acc", 
#                             "Pos_Pred_Value", "Neg_Pred_Value", "Prevalence", 
#                             "AUC1", "AUC2", "MCC"]
metric__ = "F1"

if (metric__ == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# [3] Create data for the plots
# 3-1 `ignore, zero` & all blocks
df1_1 <- data.frame("method" = rep("ignore, zero", 30),
                    "Fold"   = rep(1:5),
                    "Block"  = rep("all"),
                    "used_block"  = rep(names(all_res_532$all_block$`ignore, zero`), 
                                        each = 5),
                    "used_metric" = used_metric_,
                    "metric" = c(unlist(all_res_532$all_block$`ignore, zero`$block_1[metric__,]),
                                 unlist(all_res_532$all_block$`ignore, zero`$block_1_2[metric__,]),
                                 unlist(all_res_532$all_block$`ignore, zero`$block_1_2_3[metric__,]),
                                 unlist(all_res_532$all_block$`ignore, zero`$block_1_2_3_4[metric__,]),
                                 unlist(all_res_532$all_block$`ignore, zero`$block_1_2_3_4_5[metric__,]),
                                 unlist(all_res_532$all_block$`ignore, zero`$block_1_2_3_4_5_6[metric__,])))

# 3-2 `ignore, zero` & pred_blocks
df1_2 <- data.frame("method" = rep("ignore, zero", 30),
                    "Fold"   = rep(1:5),
                    "Block"  = rep("prediction"),
                    "used_block"  = rep(names(all_res_532$pred_blocks$`ignore, zero`), 
                                        each = 5),
                    "used_metric" = used_metric_,
                    "metric" = c(unlist(all_res_532$pred_blocks$`ignore, zero`$block_1[metric__,]),
                                 unlist(all_res_532$pred_blocks$`ignore, zero`$block_1_2[metric__,]),
                                 unlist(all_res_532$pred_blocks$`ignore, zero`$block_1_2_3[metric__,]),
                                 unlist(all_res_532$pred_blocks$`ignore, zero`$block_1_2_3_4[metric__,]),
                                 unlist(all_res_532$pred_blocks$`ignore, zero`$block_1_2_3_4_5[metric__,]),
                                 unlist(all_res_532$pred_blocks$`ignore, zero`$block_1_2_3_4_5_6[metric__,])))


# 3-3 `ignore, intercept` & all blocks
df2_1 <- data.frame("method" = rep("ignore, intercept", 30),
                    "Fold"   = rep(1:5),
                    "Block"  = rep("all"),
                    "used_block"  = rep(names(all_res_532$all_block$`ignore, intercept`), 
                                        each = 5),
                    "used_metric" = used_metric_,
                    "metric" = c(unlist(all_res_532$all_block$`ignore, intercept`$block_1[metric__,]),
                                 unlist(all_res_532$all_block$`ignore, intercept`$block_1_2[metric__,]),
                                 unlist(all_res_532$all_block$`ignore, intercept`$block_1_2_3[metric__,]),
                                 unlist(all_res_532$all_block$`ignore, intercept`$block_1_2_3_4[metric__,]),
                                 unlist(all_res_532$all_block$`ignore, intercept`$block_1_2_3_4_5[metric__,]),
                                 unlist(all_res_532$all_block$`ignore, intercept`$block_1_2_3_4_5_6[metric__,])))

# 3-4 `ignore, intercept` & pred_blocks
df2_2 <- data.frame("method" = rep("ignore, intercept", 30),
                    "Fold"   = rep(1:5),
                    "Block"  = rep("prediction"),
                    "used_block"  = rep(names(all_res_532$pred_blocks$`ignore, intercept`), 
                                        each = 5),
                    "used_metric" = used_metric_,
                    "metric" = c(unlist(all_res_532$pred_blocks$`ignore, intercept`$block_1[metric__,]),
                                 unlist(all_res_532$pred_blocks$`ignore, intercept`$block_1_2[metric__,]),
                                 unlist(all_res_532$pred_blocks$`ignore, intercept`$block_1_2_3[metric__,]),
                                 unlist(all_res_532$pred_blocks$`ignore, intercept`$block_1_2_3_4[metric__,]),
                                 unlist(all_res_532$pred_blocks$`ignore, intercept`$block_1_2_3_4_5[metric__,]),
                                 unlist(all_res_532$pred_blocks$`ignore, intercept`$block_1_2_3_4_5_6[metric__,])))

# 3-5 `impute maximise blocks` & all blocks
df3_1 <- data.frame("method" = rep("impute maximise blocks", 30),
                    "Fold"   = rep(1:5),
                    "Block"  = rep("all"),
                    "used_block"  = rep(names(all_res_532$all_block$`impute maximise blocks`), 
                                        each = 5),
                    "used_metric" = used_metric_,
                    "metric" = c(unlist(all_res_532$all_block$`impute maximise blocks`$block_1[metric__,]),
                                 unlist(all_res_532$all_block$`impute maximise blocks`$block_1_2[metric__,]),
                                 unlist(all_res_532$all_block$`impute maximise blocks`$block_1_2_3[metric__,]),
                                 unlist(all_res_532$all_block$`impute maximise blocks`$block_1_2_3_4[metric__,]),
                                 unlist(all_res_532$all_block$`impute maximise blocks`$block_1_2_3_4_5[metric__,]),
                                 unlist(all_res_532$all_block$`impute maximise blocks`$block_1_2_3_4_5_6[metric__,])))

# 3-6 `impute maximise blocks` & pred_blocks
df3_2 <- data.frame("method" = rep("impute maximise blocks", 30),
                    "Fold"   = rep(1:5),
                    "Block"  = rep("prediction"),
                    "used_block"  = rep(names(all_res_532$pred_blocks$`impute maximise blocks`), 
                                        each = 5),
                    "used_metric" = used_metric_,
                    "metric" = c(unlist(all_res_532$pred_blocks$`impute maximise blocks`$block_1[metric__,]),
                                 unlist(all_res_532$pred_blocks$`impute maximise blocks`$block_1_2[metric__,]),
                                 unlist(all_res_532$pred_blocks$`impute maximise blocks`$block_1_2_3[metric__,]),
                                 unlist(all_res_532$pred_blocks$`impute maximise blocks`$block_1_2_3_4[metric__,]),
                                 unlist(all_res_532$pred_blocks$`impute maximise blocks`$block_1_2_3_4_5[metric__,]),
                                 unlist(all_res_532$pred_blocks$`impute maximise blocks`$block_1_2_3_4_5_6[metric__,])))

# 3-7 `impute, maximise n` & all blocks
df4_1 <- data.frame("method" = rep("impute, maximise n", 30),
                    "Fold"   = rep(1:5),
                    "Block"  = rep("all"),
                    "used_block"  = rep(names(all_res_532$all_block$`impute, maximise n`), 
                                        each = 5),
                    "used_metric" = used_metric_,
                    "metric" = c(unlist(all_res_532$all_block$`impute, maximise n`$block_1[metric__,]),
                                 unlist(all_res_532$all_block$`impute, maximise n`$block_1_2[metric__,]),
                                 unlist(all_res_532$all_block$`impute, maximise n`$block_1_2_3[metric__,]),
                                 unlist(all_res_532$all_block$`impute, maximise n`$block_1_2_3_4[metric__,]),
                                 unlist(all_res_532$all_block$`impute, maximise n`$block_1_2_3_4_5[metric__,]),
                                 unlist(all_res_532$all_block$`impute, maximise n`$block_1_2_3_4_5_6[metric__,])))

# 3-8 `impute maximise n` & pred_blocks
df4_2 <- data.frame("method" = rep("impute, maximise n", 30),
                    "Fold"   = rep(1:5),
                    "Block"  = rep("prediction"),
                    "used_block"  = rep(names(all_res_532$pred_blocks$`impute, maximise n`), 
                                        each = 5),
                    "used_metric" = used_metric_,
                    "metric" = c(unlist(all_res_532$pred_blocks$`impute, maximise n`$block_1[metric__,]),
                                 unlist(all_res_532$pred_blocks$`impute, maximise n`$block_1_2[metric__,]),
                                 unlist(all_res_532$pred_blocks$`impute, maximise n`$block_1_2_3[metric__,]),
                                 unlist(all_res_532$pred_blocks$`impute, maximise n`$block_1_2_3_4[metric__,]),
                                 unlist(all_res_532$pred_blocks$`impute, maximise n`$block_1_2_3_4_5[metric__,]),
                                 unlist(all_res_532$pred_blocks$`impute, maximise n`$block_1_2_3_4_5_6[metric__,])))

# 3-9 mdd_SPLS
df5_1 <- data.frame("method" = rep("mdd_SPLS", 5),
                    "Fold"   = rep(1:5),
                    "Block"  = "all",
                    "used_block"  = "block_1_2_3_4_5_6",
                    "used_metric" = used_metric_,
                    "metric" = c(unlist(all_res_532$mdd_splss["F1",])))

# [4] Create a DF to plot the results!
# 4-1 Bind the DFs
plot_df <- rbind(df1_1, df1_2, df2_1, df2_2, df3_1, df3_2, df4_1, df4_2)

# 4-2 Only keep the 'all' results, as the results are almost identical
plot_df <- plot_df[plot_df$Block == "all",]

# 4-3 Rename the used blocks [used for the prediction]
plot_df$used_block2 <- sapply(1:nrow(plot_df), FUN = function(x) {
  strsplit(as.character(plot_df$used_block[x]), split = "block_")[[1]][2]
})
plot_df$used_block2 <- gsub("_", ", ", plot_df$used_block2)


# 4-4 Rename the approaches so it is consitent!
levels(plot_df$method) <- c("PL - ignore, zero", "PL - ignore, intercept",
                            "PL - impute, maximise blocks", "PL - impute, maximise n")

plot_df$method <- factor(plot_df$method, levels = c( "PL - ignore, intercept", 
                                                     "PL - ignore, zero", 
                                                     "PL - impute, maximise blocks",
                                                     "PL - impute, maximise n"))

# 4-4 Do the plot
ggplot(data = plot_df, aes(x = method, y = metric, fill = used_block2)) +
  geom_boxplot() + 
  theme_bw() +
  ylab(used_metric_) +
  xlab("Different Approaches") +
  ggtitle("Priority-Lasso - split by the blocks used for the prediction",
          subtitle = "Clinical asthma data") +
  theme(text = element_text(size = 24),
        axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(fill = "Blocks used for \nthe predictions") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5))

# 4-5 Get a summary for the different approaches
for (method_ in levels(plot_df$method)) {
  print("METHOD ----------------------------")
  print(method_)
  
  print(sapply(levels(plot_df$used_block), FUN = function(blocks_) {
    summary(plot_df$metric[plot_df$method == method_ & plot_df$used_block == blocks_])
  }))
}

summary(df5_1$metric)


# Hagenberg --- Setting 5_3_4  --- 1 [4, 2, 1, 3, 5]                         ----
# [1] Load the Metrics of the CV w/ Hagenbergs Approaches
load("./docs/CV_Res/REAL/Hagenberg_5_3_4__Setting1.R")

# [2] Select a metric from:  ["Accuracy", "Kappa", "Sensitifity", "Specificity", 
#                             "Precision", "Recall", "F1", "Balance_Acc", 
#                             "Pos_Pred_Value", "Neg_Pred_Value", "Prevalence", 
#                             "AUC1", "AUC2", "MCC"]
metric__ = "F1"

if (metric__ == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# [3] Create a DF to plot the results!
# 3-0-1 Define a DF to store the results
df_all <- data.frame(approach = factor(),
                     fold     = integer(),
                     blocks   = factor(),
                     metrics  = numeric(),
                     used_metric = factor())

# 3-1 'ignore, zero' Approach
for (used_blocks in names(all_res$`ignore, zero`)) {
  
  # extract the metrics for 'used_blocks' 
  curr_metrics_ <- unlist(all_res$`ignore, zero`[[used_blocks]][metric__,])
  
  # bind it to 'df_all'
  df_all <- rbind(df_all, data.frame(approach = "ignore, zero",
                                     fold     = 1:5,
                                     blocks   = used_blocks,
                                     metrics  = curr_metrics_,
                                     used_metric = used_metric_))
}

# 3-2 'ignore, intercept' Approach
for (used_blocks in names(all_res$`ignore, intercept`)) {
  
  # extract the metrics for 'used_blocks' 
  curr_metrics_ <- unlist(all_res$`ignore, intercept`[[used_blocks]][metric__,])
  
  # bind it to 'df_all'
  df_all <- rbind(df_all, data.frame(approach = "ignore, intercept",
                                     fold     = 1:5,
                                     blocks   = used_blocks,
                                     metrics  = curr_metrics_,
                                     used_metric = used_metric_))
}

# 3-3 'impute, maximise blocks' Approach
for (used_blocks in names(all_res$`impute maximise blocks`)) {
  
  # extract the metrics for 'used_blocks' 
  curr_metrics_ <- unlist(all_res$`impute maximise blocks`[[used_blocks]][metric__,])
  
  # bind it to 'df_all'
  df_all <- rbind(df_all, data.frame(approach = "impute maximise blocks",
                                     fold     = 1:5,
                                     blocks   = used_blocks,
                                     metrics  = curr_metrics_,
                                     used_metric = used_metric_))
}

# 3-4 'impute, maximise n' Approach
for (used_blocks in names(all_res$`impute, maximise n`)) {
  
  # extract the metrics for 'used_blocks' 
  curr_metrics_ <- unlist(all_res$`impute, maximise n`[[used_blocks]][metric__,])
  
  # bind it to 'df_all'
  df_all <- rbind(df_all, data.frame(approach = "impute, maximise n",
                                     fold     = 1:5,
                                     blocks   = used_blocks,
                                     metrics  = curr_metrics_,
                                     used_metric = used_metric_))
}

# 3-5 Add the mdd-sPLS approach to df_all
df_all <- rbind(df_all, data.frame(approach = "mdd-sPLS",
                                   fold     = 1:5,
                                   blocks   = "block_4_2_1_3_5",
                                   metrics  = unlist(all_res$mdd_spls[metric__,]),
                                   used_metric = used_metric_))

# [4] Do the plot
# 4-1 Rename the used blocks [used for the prediction]
df_all$used_block2 <- sapply(1:nrow(df_all), FUN = function(x) {
  strsplit(as.character(df_all$blocks[x]), split = "block_")[[1]][2]
})
df_all$used_block2 <- gsub("_", ", ", df_all$used_block2)

# 4-2 Rename the approaches so it is consitent!
levels(df_all$approach) <- c("PL - ignore, zero", "PL - ignore, intercept",
                             "PL - impute, maximise blocks", "PL - impute, maximise n", "mdd-sPLS")

df_all$approach <- factor(df_all$approach, levels = c("mdd-sPLS",
                                                      "PL - ignore, intercept", 
                                                      "PL - ignore, zero", 
                                                      "PL - impute, maximise blocks",
                                                      "PL - impute, maximise n"))

# 4-3 Do the actual plot
ggplot(data = df_all, aes(x = approach, y = metrics, fill = used_block2)) +
  geom_boxplot() + 
  theme_bw() +
  ylab(used_metric_) +
  xlab("Different Approaches") +
  ggtitle("Priority-Lasso adaptions & mdd-sPLS approach",
          subtitle = "Clinical asthma data with block-priorities: 4, 2, 1, 3, 5") +
  theme(text = element_text(size = 24),
        axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(fill = "Blocks used for \nthe predictions") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5))

# 4-4 Get the summarys of the single approaches!
for (method_ in levels(df_all$approach)) {
  print("METHOD ----------------------------")
  print(method_)
  
  print(sapply(levels(df_all$blocks), FUN = function(blocks_) {
    summary(df_all$metric[df_all$approach == method_ & df_all$blocks == blocks_])
  }))
}

# Hagenberg --- Setting 5_3_4  --- 2 [4, 2, 1, 3, 6]                         ----
# [1] Load the Metrics of the CV w/ Hagenbergs Approaches
load("./docs/CV_Res/REAL/Hagenberg_5_3_4__Setting2.R")

# [2] Select a metric from:  ["Accuracy", "Kappa", "Sensitifity", "Specificity", 
#                             "Precision", "Recall", "F1", "Balance_Acc", 
#                             "Pos_Pred_Value", "Neg_Pred_Value", "Prevalence", 
#                             "AUC1", "AUC2", "MCC"]
metric__ = "F1"

if (metric__ == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# [3] Create a DF to plot the results!
# 3-0-1 Define a DF to store the results
df_all <- data.frame(approach = factor(),
                     fold     = integer(),
                     blocks   = factor(),
                     metrics  = numeric(),
                     used_metric = factor())

# 3-1 'ignore, zero' Approach
for (used_blocks in names(all_res$`ignore, zero`)) {
  
  # extract the metrics for 'used_blocks' 
  curr_metrics_ <- unlist(all_res$`ignore, zero`[[used_blocks]][metric__,])
  
  # bind it to 'df_all'
  df_all <- rbind(df_all, data.frame(approach = "ignore, zero",
                                     fold     = 1:5,
                                     blocks   = used_blocks,
                                     metrics  = curr_metrics_,
                                     used_metric = used_metric_))
}

# 3-2 'ignore, intercept' Approach
for (used_blocks in names(all_res$`ignore, intercept`)) {
  
  # extract the metrics for 'used_blocks' 
  curr_metrics_ <- unlist(all_res$`ignore, intercept`[[used_blocks]][metric__,])
  
  # bind it to 'df_all'
  df_all <- rbind(df_all, data.frame(approach = "ignore, intercept",
                                     fold     = 1:5,
                                     blocks   = used_blocks,
                                     metrics  = curr_metrics_,
                                     used_metric = used_metric_))
}

# 3-3 'impute, maximise blocks' Approach
for (used_blocks in names(all_res$`impute maximise blocks`)) {
  
  # extract the metrics for 'used_blocks' 
  curr_metrics_ <- unlist(all_res$`impute maximise blocks`[[used_blocks]][metric__,])
  
  # bind it to 'df_all'
  df_all <- rbind(df_all, data.frame(approach = "impute maximise blocks",
                                     fold     = 1:5,
                                     blocks   = used_blocks,
                                     metrics  = curr_metrics_,
                                     used_metric = used_metric_))
}

# 3-4 'impute, maximise n' Approach
for (used_blocks in names(all_res$`impute, maximise n`)) {
  
  # extract the metrics for 'used_blocks' 
  curr_metrics_ <- unlist(all_res$`impute, maximise n`[[used_blocks]][metric__,])
  
  # bind it to 'df_all'
  df_all <- rbind(df_all, data.frame(approach = "impute, maximise n",
                                     fold     = 1:5,
                                     blocks   = used_blocks,
                                     metrics  = curr_metrics_,
                                     used_metric = used_metric_))
}

# 3-5 Add the mdd-sPLS approach to df_all
df_all <- rbind(df_all, data.frame(approach = "mdd-sPLS",
                                   fold     = 1:5,
                                   blocks   = "block_4_2_1_3_6",
                                   metrics  = unlist(all_res$mdd_spls[metric__,]),
                                   used_metric = used_metric_))

# [4] Do the plot
# 4-1 Rename the used blocks [used for the prediction]
df_all$used_block2 <- sapply(1:nrow(df_all), FUN = function(x) {
  strsplit(as.character(df_all$blocks[x]), split = "block_")[[1]][2]
})
df_all$used_block2 <- gsub("_", ", ", df_all$used_block2)

# 4-2 Rename the approaches so it is consitent!
levels(df_all$approach) <- c("PL - ignore, zero", "PL - ignore, intercept",
                             "PL - impute, maximise blocks", "PL - impute, maximise n", 
                             "mdd-sPLS")

df_all$approach <- factor(df_all$approach, levels = c("mdd-sPLS",
                                                      "PL - ignore, intercept", 
                                                      "PL - ignore, zero", 
                                                      "PL - impute, maximise blocks",
                                                      "PL - impute, maximise n"))

# 4-3 Do the actual plot
ggplot(data = df_all, aes(x = approach, y = metrics, fill = used_block2)) +
  geom_boxplot() + 
  theme_bw() +
  ylab(used_metric_) +
  xlab("Different Approaches") +
  ggtitle("Priority-Lasso adaptions & mdd-sPLS approach",
          subtitle = "Clinical asthma data with block-priorities: 4, 2, 1, 3, 6") +
  theme(text = element_text(size = 24),
        axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(fill = "Blocks used for \nthe predictions") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5))

# 4-4 Get the summarys of the single approaches!
for (method_ in levels(df_all$approach)) {
  print("METHOD ----------------------------")
  print(method_)
  
  print(sapply(levels(df_all$blocks), FUN = function(blocks_) {
    summary(df_all$metric[df_all$approach == method_ & df_all$blocks == blocks_])
  }))
}

# 4-1 Rename the used blocks [used for the prediction]
df_all$used_block2 <- sapply(1:nrow(df_all), FUN = function(x) {
  strsplit(as.character(df_all$blocks[x]), split = "block_")[[1]][2]
})

# 4-2 Rename the approaches so it is consitent!
levels(df_all$approach) <- c("PL - ignore, zero", "PL - ignore, intercept",
                             "PL - impute, maximise blocks", "PL - impute, maximise n", "mdd-sPLS")

df_all$approach <- factor(df_all$approach, levels = c("mdd-sPLS",
                                                      "PL - ignore, intercept", 
                                                      "PL - ignore, zero", 
                                                      "PL - impute, maximise blocks",
                                                      "PL - impute, maximise n"))

# 4-3 Do the actual plot
ggplot(data = df_all, aes(x = approach, y = metrics, fill = blocks)) +
  geom_boxplot() + 
  theme_bw() +
  ylab(used_metric_) +
  xlab("Different Approaches") +
  ggtitle("Priority-Lasso with different block-priorities",
          subtitle = "Clinical asthma data") +
  theme(text = element_text(size = 24),
        axis.text.x = element_text(angle = 25, hjust = 1)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5))

# Hagenberg --- Setting 5_3_4  --- 3 [4, 2, 1, 5, 3]                         ----
# [1] Load the Metrics of the CV w/ Hagenbergs Approaches
load("./docs/CV_Res/REAL/Hagenberg_5_3_4__Setting3.R")

# [2] Select a metric from:  ["Accuracy", "Kappa", "Sensitifity", "Specificity", 
#                             "Precision", "Recall", "F1", "Balance_Acc", 
#                             "Pos_Pred_Value", "Neg_Pred_Value", "Prevalence", 
#                             "AUC1", "AUC2", "MCC"]
metric__ = "F1"

if (metric__ == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# [3] Create a DF to plot the results!
# 3-0-1 Define a DF to store the results
df_all <- data.frame(approach = factor(),
                     fold     = integer(),
                     blocks   = factor(),
                     metrics  = numeric(),
                     used_metric = factor())

# 3-1 'ignore, zero' Approach
for (used_blocks in names(all_res$`ignore, zero`)) {
  
  # extract the metrics for 'used_blocks' 
  curr_metrics_ <- unlist(all_res$`ignore, zero`[[used_blocks]][metric__,])
  
  # bind it to 'df_all'
  df_all <- rbind(df_all, data.frame(approach = "ignore, zero",
                                     fold     = 1:5,
                                     blocks   = used_blocks,
                                     metrics  = curr_metrics_,
                                     used_metric = used_metric_))
}

# 3-2 'ignore, intercept' Approach
for (used_blocks in names(all_res$`ignore, intercept`)) {
  
  # extract the metrics for 'used_blocks' 
  curr_metrics_ <- unlist(all_res$`ignore, intercept`[[used_blocks]][metric__,])
  
  # bind it to 'df_all'
  df_all <- rbind(df_all, data.frame(approach = "ignore, intercept",
                                     fold     = 1:5,
                                     blocks   = used_blocks,
                                     metrics  = curr_metrics_,
                                     used_metric = used_metric_))
}

# 3-3 'impute, maximise blocks' Approach
for (used_blocks in names(all_res$`impute maximise blocks`)) {
  
  # extract the metrics for 'used_blocks' 
  curr_metrics_ <- unlist(all_res$`impute maximise blocks`[[used_blocks]][metric__,])
  
  # bind it to 'df_all'
  df_all <- rbind(df_all, data.frame(approach = "impute maximise blocks",
                                     fold     = 1:5,
                                     blocks   = used_blocks,
                                     metrics  = curr_metrics_,
                                     used_metric = used_metric_))
}

# 3-4 'impute, maximise n' Approach
for (used_blocks in names(all_res$`impute, maximise n`)) {
  
  # extract the metrics for 'used_blocks' 
  curr_metrics_ <- unlist(all_res$`impute, maximise n`[[used_blocks]][metric__,])
  
  # bind it to 'df_all'
  df_all <- rbind(df_all, data.frame(approach = "impute, maximise n",
                                     fold     = 1:5,
                                     blocks   = used_blocks,
                                     metrics  = curr_metrics_,
                                     used_metric = used_metric_))
}

# [4] Do the plot
ggplot(data = df_all, aes(x = approach, y = metrics, fill = blocks)) +
  geom_boxplot() + 
  theme_bw() +
  ylab(used_metric_) +
  xlab("Different Approaches") +
  ggtitle("Priority-Lasso with different block-priorities",
          subtitle = "Clinical asthma data") +
  theme(text = element_text(size = 24),
        axis.text.x = element_text(angle = 25, hjust = 1)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5))

# Hagenberg --- Setting 5_3_4  --- 4 [4, 2, 1, 6, 3]                         ----
# [1] Load the Metrics of the CV w/ Hagenbergs Approaches
load("./docs/CV_Res/REAL/Hagenberg_5_3_4__Setting4.R")

# [2] Select a metric from:  ["Accuracy", "Kappa", "Sensitifity", "Specificity", 
#                             "Precision", "Recall", "F1", "Balance_Acc", 
#                             "Pos_Pred_Value", "Neg_Pred_Value", "Prevalence", 
#                             "AUC1", "AUC2", "MCC"]
metric__ = "F1"

if (metric__ == "F1") {
  used_metric_ <- "Metric: F-1 Score"
} else {
  used_metric_ <- paste("Metric:", DF_all$performance_metric[1])
}

# [3] Create a DF to plot the results!
# 3-0-1 Define a DF to store the results
df_all <- data.frame(approach = factor(),
                     fold     = integer(),
                     blocks   = factor(),
                     metrics  = numeric(),
                     used_metric = factor())

# 3-1 'ignore, zero' Approach
for (used_blocks in names(all_res$`ignore, zero`)) {
  
  # extract the metrics for 'used_blocks' 
  curr_metrics_ <- unlist(all_res$`ignore, zero`[[used_blocks]][metric__,])
  
  # bind it to 'df_all'
  df_all <- rbind(df_all, data.frame(approach = "ignore, zero",
                                     fold     = 1:5,
                                     blocks   = used_blocks,
                                     metrics  = curr_metrics_,
                                     used_metric = used_metric_))
}

# 3-2 'ignore, intercept' Approach
for (used_blocks in names(all_res$`ignore, intercept`)) {
  
  # extract the metrics for 'used_blocks' 
  curr_metrics_ <- unlist(all_res$`ignore, intercept`[[used_blocks]][metric__,])
  
  # bind it to 'df_all'
  df_all <- rbind(df_all, data.frame(approach = "ignore, intercept",
                                     fold     = 1:5,
                                     blocks   = used_blocks,
                                     metrics  = curr_metrics_,
                                     used_metric = used_metric_))
}

# 3-3 'impute, maximise blocks' Approach
for (used_blocks in names(all_res$`impute maximise blocks`)) {
  
  # extract the metrics for 'used_blocks' 
  curr_metrics_ <- unlist(all_res$`impute maximise blocks`[[used_blocks]][metric__,])
  
  # bind it to 'df_all'
  df_all <- rbind(df_all, data.frame(approach = "impute maximise blocks",
                                     fold     = 1:5,
                                     blocks   = used_blocks,
                                     metrics  = curr_metrics_,
                                     used_metric = used_metric_))
}

# 3-4 'impute, maximise n' Approach
for (used_blocks in names(all_res$`impute, maximise n`)) {
  
  # extract the metrics for 'used_blocks' 
  curr_metrics_ <- unlist(all_res$`impute, maximise n`[[used_blocks]][metric__,])
  
  # bind it to 'df_all'
  df_all <- rbind(df_all, data.frame(approach = "impute, maximise n",
                                     fold     = 1:5,
                                     blocks   = used_blocks,
                                     metrics  = curr_metrics_,
                                     used_metric = used_metric_))
}

# [4] Do the plot
ggplot(data = df_all, aes(x = approach, y = metrics, fill = blocks)) +
  geom_boxplot() + 
  theme_bw() +
  ylab(used_metric_) +
  xlab("Different Approaches") +
  ggtitle("Priority-Lasso with different block-priorities",
          subtitle = "Clinical asthma data") +
  theme(text = element_text(size = 24),
        axis.text.x = element_text(angle = 25, hjust = 1)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5))


# Comparison of the approaches                                               ----
# [1] ----- CC Results
DF_CC <- data.frame()

# 1-1 List the files in the path
data_path <- "./docs/CV_Res/REAL//CC_Approach.R"

# 1-2 Get results
file_curr <- load(data_path)
file_curr     <- eval(as.symbol(file_curr))

DF_CC <- extract_metrics(x = file_curr, metric = "F1")

# 1-3 Add the approach to the DF
DF_CC$Approach <- "Complete-Case"

# [2] ----- SB Results
DF_SB <- data.frame()

# 2-1 List the files in the path
data_path <- "./docs/CV_Res/REAL/SingleBlockApproach.R"

# 2-2 Get results
file_curr <- load(data_path)
file_curr <- eval(as.symbol(file_curr))

DF_SB <- extract_metrics_SBAPP(x = file_curr, metric = "F1")

# 2-3 Add the approach to the DF
DF_SB$Approach <- "Single-Block"

# 2-4 Only keep the block 'Questionaire' as the results were the best!
DF_SB <- DF_SB[DF_SB$Feature_Block == "df1",]

# [3] ----- IMP Results
DF_IMP <- data.frame()

# 3-1 List the files in the path
data_path <- "./docs/CV_Res/REAL/Imputation_Approach.R"

# 3-2 Get results
file_curr <- load(data_path)
file_curr <- eval(as.symbol(file_curr))

DF_IMP <- extract_metrics(x = file_curr, metric = "F1")

# 3-3 Add the approach to the DF
DF_IMP$Approach <- "Imputation"

# [4] ----- Block-Wise Results
DF_BW <- data.frame()

# 4-1 List the files in the path
data_path <- "./docs/CV_Res/REAL/BlockWise_Approach.R"

# 4-2 Get results
file_curr <- load(data_path)
file_curr <- eval(as.symbol(file_curr))

DF_BW <- extract_metrics_FW_BW(x = file_curr, metric = "F1")

# 4-3 Add the approach to the DF
DF_BW$Approach <- "Block-wise"

# 4-4 Only keep the ones with weight_metric == "F-1 Score" 
DF_BW <- DF_BW[DF_BW$weight_metric == 'F1-Score',]
DF_BW$weight_metric <- "F-1 Score"

# [5] ----- Fold-Wise Results
DF_FW <- data.frame()

# 5-1 List the files in the path
data_path <- "./docs/CV_Res/REAL/FoldWise_Approach.R"

# 5-2 Get results
file_curr <- load(data_path)
file_curr <- eval(as.symbol(file_curr))

DF_FW <- extract_metrics_FW_BW(x = file_curr, metric = "F1")

# 5-3 Add the approach to the DF
DF_FW$Approach <- "Fold-wise"

# 5-4 Only keep the ones with weight_metric == "F-1 Score" 
DF_FW <- DF_FW[DF_FW$weight_metric == 'F1-Score',]
DF_FW$weight_metric <- "F-1 Score"

# [6] ----- Bind the DFs and create a plot
# 6-1 Columns we need from every approach-DF
needed_cols <- c("Metric", "Approach")

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
ggplot(data = DF_all, aes(x = Approach, y = Metric)) +
  geom_boxplot(position = position_dodge(preserve = "single"),
               fill = c("darkolivegreen2", "darkorange1", 'darkorchid1', 
                        'darkgoldenrod1', 'darkslategray3')) + 
  theme_bw() +
  ggtitle("Comparison of all Approaches",
          subtitle = "clinical asthma data") +
  ylab(used_metric_) +
  xlab("Approach") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        text = element_text(size = 24),
        legend.position = "top") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5),
             col = "red", lty = 2, lwd = 1.005)