"
Script to visualize the different Results from the CVs
"
setwd("C:/Users/kuche_000/Desktop/MS-Thesis/")
library(ggplot2)
require(gridExtra)
library(reshape2)

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
  if (!(metric  %in% names(x$res_all[[1]][[1]]))) {
    print("There exist only the metrics:\n")
    print(names(x$res_all[[1]][[1]]))
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
  #       & 
  for (j in seq_len(ncol(res_curr_train_set))) {
    
    # 1-3-1-1 Extract TrainSetting, folds, test_situation and the metrics
    #         for the given column!
    train_sit_ <- rep(train_sit, times = nrow(res_curr_train_set))
    folds      <- seq_len(nrow(res_curr_train_set))
    test_sit   <- rep(colnames(res_curr_train_set)[j], times = nrow(res_curr_train_set))
    metric_c   <- unlist(res_curr_train_set[,j])
    
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
  ifelse(!x_settings$weighted,
         DF_final$weight_metric <- NA,
         DF_final$weight_metric <- x_settings$weight_metric)
  
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
DF1 <- read.csv2(paste0(data_path, "/", files[1]), stringsAsFactors = F)
DF1$seed <- "1234"

DF2 <- read.csv2(paste0(data_path, "/", files[2]), stringsAsFactors = F)
DF2$seed <- "123456789"

DF3 <- read.csv2(paste0(data_path, "/", files[3]), stringsAsFactors = F)
DF3$seed <- "1235"

DF4 <- read.csv2(paste0(data_path, "/", files[4]), stringsAsFactors = F)
DF4$seed <- "1236"

DF5 <- read.csv2(paste0(data_path, "/", files[5]), stringsAsFactors = F)
DF5$seed <- "1237"

DF6 <- read.csv2(paste0(data_path, "/", files[6]), stringsAsFactors = F)
DF6$seed <- "1238"

DF7 <- read.csv2(paste0(data_path, "/", files[7]), stringsAsFactors = F)
DF7$seed <- "1239"

DF8 <- read.csv2(paste0(data_path, "/", files[8]), stringsAsFactors = F)
DF8$seed <- "1240"

# 1-3 Bind them to a single DF & convert numeric cols to numeric!
DF_all <- rbind(DF1, DF2, DF3, DF4, DF5, DF6, DF7, DF8)

numeric_cols <- c("OOB_Acc", "Test_Acc", "Test_F1")
DF_all[,numeric_cols] <- sapply(numeric_cols, 
                                function(x) as.numeric(DF_all[,x]))

# 1-4 reshape DF_all for the plot!
plot_df <- melt(DF_all, id.vars = c("Data", "seed"), measure.vars = c("Test_Acc", "Test_F1"))

# [2] Do the plots
# 2-1 Split by the different Seeds!
ggplot(data = plot_df, aes(x = seed, y = value, fill = variable)) +
  geom_boxplot() + 
  theme_bw() +
  ggtitle("Joint Block Performance on the 14 fixed DF w/ 5 fold CV",
          subtitle = "Final-Subset: 10% mirna & mutation, 5% rna & 2.5% cnv") +
  xlab("Seed used to subset") +
  ylab("Metric [Acc & F1]") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15))

# 2-2 For a certain Dataframe!
df_tmp_sub = "PAAD_subset.RData"

df_temp <- plot_df[plot_df$Data == df_tmp_sub,]
ggplot(data = df_temp, aes(x = seed, y = value, fill = variable)) +
  geom_boxplot() + 
  theme_bw() +
  ggtitle(paste("Joint Block Performance on:", df_tmp_sub),
          subtitle = "finalsubset: 10% mirna & mutation, 5% rna & 2.5% cnv\n --- 5-fold- CV on DF") +
  xlab("Seed used to subset") +
  ylab("metric") +
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
# 1-2 Add the amount of subsets to each block:
DF1 <- read.csv2(paste0(data_path, "/", files[1]), stringsAsFactors = F)
DF1$seed <- "1234"

DF2 <- read.csv2(paste0(data_path, "/", files[2]), stringsAsFactors = F)
DF2$seed <- "123456789"

DF3 <- read.csv2(paste0(data_path, "/", files[3]), stringsAsFactors = F)
DF3$seed <- "1235"

DF4 <- read.csv2(paste0(data_path, "/", files[4]), stringsAsFactors = F)
DF4$seed <- "1236"

DF5 <- read.csv2(paste0(data_path, "/", files[5]), stringsAsFactors = F)
DF5$seed <- "1237"

DF6 <- read.csv2(paste0(data_path, "/", files[6]), stringsAsFactors = F)
DF6$seed <- "1238"

DF7 <- read.csv2(paste0(data_path, "/", files[7]), stringsAsFactors = F)
DF7$seed <- "1239"

DF8 <- read.csv2(paste0(data_path, "/", files[8]), stringsAsFactors = F)
DF8$seed <- "1240"

# 1-3 Bind them to a single DF & convert numeric cols to numeric!
DF_all <- rbind(DF1, DF2, DF3, DF4, DF5, DF6, DF7, DF8)

num_cols <- c("OOB_Acc", "Test_Acc", "Test_F1")
DF_all[,num_cols] <- sapply(num_cols, 
                            function(x) as.numeric(DF_all[,x]))

# 1-4 reshape the layout of data for the plot! 
plot_df <- melt(DF_all, id.vars = c("Data", "Block", "seed"), 
                measure.vars = c("Test_Acc", "Test_F1"))

# [2] Do the plot
ggplot(data = plot_df, aes(x = Block, y = value, fill = variable)) + 
  geom_boxplot() + 
  facet_grid(. ~ seed) +
  theme_bw() +
  ggtitle("Single Block Performance on the 14 fixed DFs w/ 5 fold CV",
          subtitle = "Final-Subset: 10% mirna & mutation, 5% rna & 2.5% cnv \n--- split by the seeds used to subset the featurespace") +
  xlab("Blocks used as Feature Space") +
  ylab("Metric [Acc & F1]") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15))

# Analyse Results of Romans Approach on the fixed DFs                       ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/gender/Roman_final_subsets/TrainSetting1"

# [1] Load all Results w/ Romans Approach
# 1-1 get all files in the folder, that are singleblock performances
files <- list.files(data_path)

# 1-2 Load the Results for all TrainSettings
DF1 <- load(paste0(data_path, "/", files[1]))
DF1 <- eval(as.symbol(DF1))
DF2 <- load(paste0(data_path, "/", files[2]))
DF2 <- eval(as.symbol(DF2))

# 1-3 Extract the Metrics and TestSetting as desired:
# DF1 with LGG DF 
DF1_1  <- extract_metrics(x = DF1[[1]], metric = "F1", train_sit = 1)
DF1_12 <- extract_metrics(x = DF1[[1]], metric = "Accuracy", train_sit = 1)
DF1_2  <- extract_metrics(x = DF1[[2]], metric = "F1", train_sit = 1)
DF1_22 <- extract_metrics(x = DF1[[2]], metric = "Accuracy", train_sit = 1)
DF1_3  <- extract_metrics(x = DF1[[3]], metric = "F1", train_sit = 1)
DF1_32 <- extract_metrics(x = DF1[[3]], metric = "Accuracy", train_sit = 1)

# DF2 with LIHC DF
DF2_1  <- extract_metrics(x = DF2[[1]], metric = "F1", train_sit = 1)
DF2_12 <- extract_metrics(x = DF2[[1]], metric = "Accuracy", train_sit = 1)
DF2_2  <- extract_metrics(x = DF2[[2]], metric = "F1", train_sit = 1)
DF2_22 <- extract_metrics(x = DF2[[2]], metric = "Accuracy", train_sit = 1)
DF2_3  <- extract_metrics(x = DF2[[3]], metric = "F1", train_sit = 1)
DF2_32 <- extract_metrics(x = DF2[[3]], metric = "Accuracy", train_sit = 1)

# 1-4 Bind the results to a single DF
DF1_all <- rbind(DF1_1, DF1_12, DF1_2, DF1_22, DF1_3, DF1_32)
DF2_all <- rbind(DF2_1, DF2_12, DF2_2, DF2_22, DF2_3, DF2_32)
DF_all  <- rbind(DF1_all, DF2_all)

ggplot(data = DF_all, aes(x = Testsituation, y = Metric, 
                          fill = performance_metric)) +
  geom_boxplot() + 
  theme_bw() +
  ggtitle("Romans Method applied to 2 fixed data subsets",
          subtitle = "split by the weighting used for ensembling!") +
  xlab("TestSituations") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15)) +
  facet_grid(. ~ weight_metric)


# Analyse Results of Norberts Approach on the fixed DFs                     ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/gender/Norbert_final_subsets/TrainSetting1"

# [1] Load all Results w/ Romans Approach
# 1-1 get all files in the folder, that are singleblock performances
files <- list.files(data_path)

# 1-2 Load the Results for all TrainSettings
DF1 <- load(paste0(data_path, "/", files[1]))
DF1 <- eval(as.symbol(DF1))
DF2 <- load(paste0(data_path, "/", files[2]))
DF2 <- eval(as.symbol(DF2))

# 1-3 Extract the Metrics and TestSetting as desired:
# DF1 with LGG DF 
DF1_1  <- extract_metrics(x = DF1[[1]], metric = "F1", train_sit = 1)
DF1_12 <- extract_metrics(x = DF1[[1]], metric = "Accuracy", train_sit = 1)
DF1_2  <- extract_metrics(x = DF1[[2]], metric = "F1", train_sit = 1)
DF1_22 <- extract_metrics(x = DF1[[2]], metric = "Accuracy", train_sit = 1)
DF1_3  <- extract_metrics(x = DF1[[3]], metric = "F1", train_sit = 1)
DF1_32 <- extract_metrics(x = DF1[[3]], metric = "Accuracy", train_sit = 1)

# DF2 with LIHC DF
DF2_1  <- extract_metrics(x = DF2[[1]], metric = "F1", train_sit = 1)
DF2_12 <- extract_metrics(x = DF2[[1]], metric = "Accuracy", train_sit = 1)
DF2_2  <- extract_metrics(x = DF2[[2]], metric = "F1", train_sit = 1)
DF2_22 <- extract_metrics(x = DF2[[2]], metric = "Accuracy", train_sit = 1)
DF2_3  <- extract_metrics(x = DF2[[3]], metric = "F1", train_sit = 1)
DF2_32 <- extract_metrics(x = DF2[[3]], metric = "Accuracy", train_sit = 1)

# 1-4 Bind the results to a single DF
DF1_all <- rbind(DF1_1, DF1_12, DF1_2, DF1_22, DF1_3, DF1_32)
DF2_all <- rbind(DF2_1, DF2_12, DF2_2, DF2_22, DF2_3, DF2_32)
DF_all  <- rbind(DF1_all, DF2_all)

ggplot(data = DF_all, aes(x = Testsituation, y = Metric, 
                          fill = performance_metric)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("Norberts Method applied to 2 fixed data subsets",
          subtitle = "split by the weighting used for ensembling!") +
  xlab("TestSituations") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15)) +
  facet_grid(. ~ weight_metric)
