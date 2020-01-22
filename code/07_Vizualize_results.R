"
Script to visualize the different Results from the CVs
"
setwd("C:/Users/kuche_000/Desktop/MS-Thesis/")
library(ggplot2)
require(gridExtra)
library(reshape2)

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

# Analyse the 