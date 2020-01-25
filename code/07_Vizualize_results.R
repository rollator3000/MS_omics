"
Script to visualize the different Results from the CVs
"
setwd("C:/Users/kuche_000/Desktop/MS-Thesis/")
library(ggplot2)
require(gridExtra)
library(reshape2)

unpack_list_res <- function(x, metric, trainsetting) {
  "unpack the results from a list (x) and store it in regular DFs!
   x is the list of results we get when doing CV w/ Roman or Norberts Method!
   Each list entrance is an own Setting!
   
   Args:
    x (list)           : list filled with the metrics for every single 
                         Itteration in CV
    metric (str)       : the metric we shall extraxt from the result list!
                          --> must be in 'x'!
    trainsetting (int) : the trainsetting we want to extract the results for
                         [4 different settings in total:
                          first element  = 1. Trainsetting, ..., 
                          fourth element = 4. Trainsetting]
    
   Return:
    Dataframe suited for plots for the trainsetting 'i' and the 'metric'
  "
  # [0] Check Inputs 
  # 0-1 'x' must be a list of length 4 [1 for each TrainSetting]
  assert_list(x, len = 4)
  
  # 0-2 'trainsetting' must be integer within [1; 4]
  assert_int(trainsetting, lower = 0, upper = 4)
  
  # 0-3 Make sure, that 'metric' is within x[['trainsetting']]
  assert_string(metric)
  if (!(metric  %in% names(x[[trainsetting]]$res_all[["full"]][[1]]))) {
    print("There exist only the metrics:\n")
    print(names(x[[trainsetting]]$res_all[["full"]][[1]]))
    stop("'metric' is not within these!")
  }
  
  # [1] Unpack the lists and put them into a single regular DF
  # 1-1 Extract the needed TestSetting & get the settings it was trained on!
  x_res      <- x[[trainsetting]]$res_all
  x_settings <- x[[trainsetting]]$settings
  
  # 1-2 Extract the metric from the wanted TrainSetting for all TestSettings!
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
  
  # 1-3-1 Loop over each column in 'res_curr_train_set' [= 1 TestSetting each]
  #       and fill the lists with results!
  for (j in seq_len(ncol(res_curr_train_set))) {
    
    # 1-3-1-1 Extract TrainSetting, folds, test_situation and the metrics
    #         for the given column!
    train_sit <- rep(colnames(res_curr_train_set)[j], times = nrow(res_curr_train_set))
    folds     <- seq_len(nrow(res_curr_train_set))
    test_sit  <- rep(trainsetting, times = nrow(res_curr_train_set))
    metric_c  <- unlist(res_curr_train_set[,j])
    
    # 1-3-1-2 Bind it to the DF with all Results
    df_current <- data.frame("Trainsituation" = train_sit,
                             "Testsituation"  = test_sit,
                             "Fold"           = folds,
                             "Metric"         = metric_c)
    DF_final <- rbind(DF_final, df_current)
  }
  
  # [2] Return 'DF_final'
  # 2-1 Adjust the colname for metric, so we know which metric we extracted!
  colnames(DF_final)[which(colnames(DF_final) == "Metric")] <- metric
  return(DF_final)
}

# Analyse the explorative singleblock single block Results                  ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/gender/explorative_subsets"

# [1] Load all explorative singleblock Results
# 1-1 get all files in the folder, that are singleblock performances
files <- list.files(data_path)
files <- files[grepl("single", files)]

# 1-2 Add the amount of subsets to each block:
DF1 <- read.csv2(paste0(data_path, "/", files[1]), stringsAsFactors = F)
DF1$subset <- 1.25

DF2 <- read.csv2(paste0(data_path, "/", files[2]), stringsAsFactors = F)
DF2$subset                         <- 100
DF2[DF2$Block != "clin", "subset"] <- 10

DF3 <- read.csv2(paste0(data_path, "/", files[3]), stringsAsFactors = F)
DF3$subset <- 2.5

DF4 <- read.csv2(paste0(data_path, "/", files[4]), stringsAsFactors = F)
DF4$subset <- 5

DF5 <- read.csv2(paste0(data_path, "/", files[5]), stringsAsFactors = F)
DF5$subset <- 100

# 1-3 Bind them to a single DF
DF_all <- rbind(DF1, DF2, DF3, DF4, DF5)

# 1-4 reshape the layout of data for the plot! 
plot_df <- melt(DF_all, id.vars = c("Block", "subset"), measure.vars = c("Test_Acc", "Test_F1", "subset", "Block"))
plot_df <- plot_df[plot_df$variable %in% c("Test_Acc", "Test_F1"),]
plot_df$value <- as.numeric(plot_df$value)

# [2] Do the plot
ggplot(data = plot_df, aes(x = Block, y = value, fill = variable)) + 
  geom_boxplot() + 
  facet_grid(. ~ subset) +
  theme_bw() +
  ggtitle("EXPLORATIVE - Single Block Performance on all 14 DFs",
          subtitle = "split by the amount of subset we used as feas") +
  xlab("Blocks used as Feature Space") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15))

# Analyse the explorative joint block block Results                         ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/gender/explorative_subsets"

# [1] Load all explorative singleblock Results
# 1-1 get all files in the folder, that are singleblock performances
files <- list.files(data_path)
files <- files[grepl("joint", files)]

# 1-2 Add the amount of subsets to each block:
DF1 <- read.csv2(paste0(data_path, "/", files[1]), stringsAsFactors = F)
DF1$subset <- "10%_all"

DF2 <- read.csv2(paste0(data_path, "/", files[2]), stringsAsFactors = F)
DF2$subset <- "10%_Mirna&Mutation_2.5%CNV&RNA"

DF3 <- read.csv2(paste0(data_path, "/", files[3]), stringsAsFactors = F)
DF3$subset <- "10%_Mirna&Mutation_2.5%RNA_1.25%CNV"

DF4 <- read.csv2(paste0(data_path, "/", files[4]), stringsAsFactors = F)
DF4$subset <- "10%_Mirna&Mutation_5%RNA&CNV"

# 1-3 Bind them to a single DF
DF_all <- rbind(DF1, DF2, DF4, DF3)

# 1-4 reshape DF_all for the plot!
plot_df <- melt(DF_all, id.vars = c("Data", "subset"), measure.vars = c("Test_Acc", "Test_F1", "subset", "Data"))
plot_df <- plot_df[plot_df$variable %in% c("Test_Acc", "Test_F1"),]
plot_df$value <- as.numeric(plot_df$value)

# [2] Do the plot, split by subsets 
ggplot(data = plot_df, aes(x = subset, y = value, fill = variable)) +
  geom_boxplot() + 
  theme_bw() +
  ggtitle("EXPLORATIVE - Joint Block Performance on all 14 DFs") +
  xlab("Subsetting of single blocks") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15))

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

# 1-3 Bind them to a single DF & convert numeric cols to numeric!
DF_all <- rbind(DF1, DF2, DF3, DF4, DF5, DF6)

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

# 1-3 Bind them to a single DF & convert numeric cols to numeric!
DF_all <- rbind(DF1, DF2, DF3, DF4, DF5, DF6)

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
data_path <- "./docs/CV_Res/gender/Roman_final_subsets"

# [1] Load all Results w/ Romans Approach
# 1-1 get all files in the folder, that are singleblock performances
files <- list.files(data_path)

# 1-2 Load the Results for all TrainSettings
DF1 <- load(paste0(data_path, "/", files[1]))
DF1 <- eval(as.symbol(DF1))
DF2 <- load(paste0(data_path, "/", files[2]))
DF2 <- eval(as.symbol(DF2))

# 1-3 Extract the Metrics and TestSetting as desired:
# 1-3-1 ACCURACY
curr_res      <- unpack_list_res(x = DF1, metric = "F1", trainsetting = 1)
curr_res$seed <- "1234"

curr_res2      <- unpack_list_res(x = DF2, metric = "F1", trainsetting = 1)
curr_res2$seed <- "1312"


# 1-4 Bind the results [different seeds when substting data] to a single DF
DF_all <- rbind(curr_res, curr_res2)

DF_all[,4] <- as.numeric(DF_all[,4])

ggplot(data = DF_all, aes(x = Trainsituation, y = F1, fill = seed)) +
  geom_boxplot() + 
  theme_bw() +
  ggtitle("Romans Method applied to fixed data subsets",
          subtitle = "split by the amount of subset we used as feas") +
  xlab("TestSituations") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15))

# Analyse Results of Norberts Approach on the fixed DFs                     ----
# [0] Define needed Variables
data_path <- "./docs/CV_Res/gender/Norbert_final_subsets"

# [1] Load all Results w/ Romans Approach
# 1-1 get all files in the folder, that are singleblock performances
files <- list.files(data_path)

# 1-2 Load the Results for all TrainSettings
DF1 <- load(paste0(data_path, "/", files[1]))
DF1 <- eval(as.symbol(DF1))
DF2 <- load(paste0(data_path, "/", files[2]))
DF2 <- eval(as.symbol(DF2))

# 1-3 Extract the Metrics and TestSetting as desired:
# 1-3-1 ACCURACY
curr_res      <- unpack_list_res(x = DF1, metric = "F1", trainsetting = 2)
curr_res$seed <- "1234"

curr_res2      <- unpack_list_res(x = DF2, metric = "F1", trainsetting = 2)
curr_res2$seed <- "1312"


# 1-4 Bind the results [different seeds when substting data] to a single DF
DF_all <- rbind(curr_res, curr_res2)

DF_all[,4] <- as.numeric(DF_all[,4])

ggplot(data = DF_all, aes(x = Trainsituation, y = F1, fill = seed)) +
  geom_boxplot() + 
  theme_bw() +
  ggtitle("Norberts Method applied to fixed data subsets",
          subtitle = "split by the amount of subset we used as feas") +
  xlab("TestSituations") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15))
