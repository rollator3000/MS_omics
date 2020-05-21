"
Script to vizualize the block-wise missingness of REAL
"
# [0] Load needed packages
library(reshape2)
library(ggplot2)

# [1] Get the observed feature_blocks for each observation
# 1-1 Load Data 
load("./data/processed/real_data/data 05052020.RData")
rm(df1, df1_mod, df2, df2_mod, df3, df3_mod, 
   df4, df4_mod, df51, df51_mod, df53, df53_mod)


# 1-2 Get info whether a observation is observed in a given block
# 1-2-1 DF1
df1_res <- sapply(1:nrow(data_merged), function(x) {
  any(!(is.na(data_merged[x, grep("df1", colnames(data_merged))])))
})

# 1-2-2 DF1
df2_res <- sapply(1:nrow(data_merged), function(x) {
  any(!(is.na(data_merged[x, grep("df2", colnames(data_merged))])))
})

# 1-2-3 DF3
df3_res <- sapply(1:nrow(data_merged), function(x) {
  any(!(is.na(data_merged[x, grep("df3", colnames(data_merged))])))
})

# 1-2-4 DF4
df4_res <- sapply(1:nrow(data_merged), function(x) {
  any(!(is.na(data_merged[x, grep("df4", colnames(data_merged))])))
})

# 1-2-5 DF51
df51_res <- sapply(1:nrow(data_merged), function(x) {
  any(!(is.na(data_merged[x, grep("df51", colnames(data_merged))])))
})

# 1-2-6 DF53
df53_res <- sapply(1:nrow(data_merged), function(x) {
  any(!(is.na(data_merged[x, grep("df53", colnames(data_merged))])))
})

# 1-3 Merge all together to a a single matrix
all_res <- matrix(data = c(df1_res, df2_res, df3_res, df4_res, df51_res, df53_res),
                  nrow = nrow(data_merged),
                  ncol = 6,
                  byrow = FALSE)

colnames(all_res) <- c("Questionaire", "Clinical routine diagnostics", 
                       "Allergen sensitization", "Cytokine expression data", 
                       "Gene expression data I", "Gene expression data II")

# 1-4 Create a plot from 'all_res'
# 1-4-1 Change Layout of the melted data
melted <- reshape2::melt(all_res)

# 1-4-2 Create the plot
ggplot(melted, aes(x = Var2, y = Var1, fill = value)) + 
  geom_tile() +
  scale_fill_manual(values = c("coral2", "darkolivegreen2"),
                    labels = c("NO", "YES")) +
  labs(fill = "Observed") +
  theme_bw() +
  ylab("Observation ID") + ylim(0, 521) +
  xlab("Data Source") +
  theme(text = element_text(size = 22),
        axis.text.x = element_text(angle = 35, hjust = 1)) +
  ggtitle("Structure of the block-wise missingness for the clinical asthma data") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5), lty = 2)

