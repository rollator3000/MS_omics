# load the needed libraries
devtools::load_all("../../prioritylasso/")
library(pROC)
library(dplyr)
library(ddsPLS)

# read in the data
load("../data/data 12032020.RData")

# include the NAs of the missing persons
df1_mod <- full_join(df1, y, by = "id") %>%
  dplyr::select(-outcome) %>%
  arrange(id) %>%
  dplyr::select(-id)
df2_mod <- full_join(df2, y, by = "id") %>%
  dplyr::select(-outcome) %>%
  arrange(id) %>%
  dplyr::select(-id)
df3_mod <- full_join(df3, y, by = "id") %>%
  dplyr::select(-outcome) %>%
  arrange(id) %>%
  dplyr::select(-id)
df4_mod <- full_join(df4, y, by = "id") %>%
  dplyr::select(-outcome) %>%
  arrange(id) %>%
  dplyr::select(-id)
df51_mod <- full_join(df51, y, by = "id") %>%
  dplyr::select(-outcome) %>%
  arrange(id) %>%
  dplyr::select(-id)
df53_mod <- full_join(df53, y, by = "id") %>%
  dplyr::select(-outcome) %>%
  arrange(id) %>%
  dplyr::select(-id)

data_merged <- cbind(df1_mod,
                     df2_mod,
                     df3_mod,
                     df4_mod,
                     df51_mod,
                     df53_mod)

# delete factor variables with more than 2 levels
data_merged[, c("df1_6", "df1_42", "df1_43", "df1_44", "df1_45")] <- NULL
# delete the variable df1_7 -> in the questionnaire, this is: does your child was
# diagnosed with asthma?
data_merged[, c("df1_7")] <- NULL

# convert factors to numeric
index_factors <- unlist(lapply(1:ncol(data_merged), function(i) {
  is.factor(data_merged[, i])
}))

data_merged <- as.matrix(data_merged)
data_merged <- apply(data_merged, 2, as.numeric)

# recode 1/2 to 0/1
data_merged[, index_factors] <- data_merged[, index_factors] - 1

# block information
block_information <- list(1:44, 45:60, 61:79, 80:108, 109:190, 191:274)

# information about the outcome
index_df <- data.frame(index = 1:521,
                       outcome = missing_str$outcome,
                       group = 0)

# make a 5 fold CV, weight the observations so that there is a 50/50 mix of 0
# and 1 outcomes
# 8274
set.seed(8274)
index_outcome_0 <- sample(rep(1:5, each = 53))
index_outcome_1 <- sample(c(rep(1:5, each = 51), 1))
index_df[index_df$outcome == 0, "group"] <- index_outcome_0
index_df[index_df$outcome == 1, "group"] <- index_outcome_1

# for every group, make a prediction on the other data and use the data of this
# group as the test data
results <- lapply(1:5, function(i) {
  index_test <- index_df[index_df$group == i, "index"]
  index_train <- index_df[index_df$group != i, "index"]
  train_x <- data_merged[index_train, ]
  train_y <- index_df[index_train, "outcome"]
  test_x <- data_merged[index_test, ]
  test_y <- index_df[index_test, "outcome"]
  
  # in the current CV, there is one observation in the test set of block 1
  # for which prioritylasso can't make a prediction
  # -> therefore exclude this observation
  if (i == 1) {
    test_x <- test_x[-c(15, 68), ]
    test_y <- test_y[-c(15, 68)]
  }
  
  # make lists for ddsPLS
  train_x_list <- lapply(1:length(block_information), function(i) {
    train_x[, block_information[[i]]]
  })
  test_x_list <- lapply(1:length(block_information), function(i) {
    test_x[, block_information[[i]]]
  })
  
  result_cv <- perf_mddsPLS(Xs = train_x_list,
                            Y = as.factor(train_y),
                            n_lambda = 10,
                            R = 1,
                            NCORES = 1,
                            mode = "logit",
                            plot_result = F,
                            kfolds = 10,
                            weight = TRUE)
  model <- mddsPLS(Xs = train_x_list,
                   Y = as.factor(train_y),
                   lambda = result_cv$Optim$optim_para_all$Lambdas[1],
                   R = result_cv$Optim$optim_para_all$R[1],
                   mode = "logit",
                   weight = TRUE)
  print(paste0("fold ", i))
  
  
  # prediction for mdd-sPLS (all blocks)
  pred_ddspls <- predict(object = model,
                         newdata = test_x_list)$probY
  
  
  
  list(model = model,
       test_y = test_y,
       pred_value_ddspls = pred_ddspls)
})


# save the results
saveRDS(results, "../data/results_different_blocks_ddsPLS_weight_2020_05_20.Rds")
