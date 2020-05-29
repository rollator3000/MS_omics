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

data_merged <- cbind(df4_mod,
                     df2_mod,
                     df1_mod,
                     df3_mod,
                     df51_mod)

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
block_information <- list(1:29, 30:45, 46:89, 90:108, 109:190)

# information about the outcome
index_df <- data.frame(index = 1:521,
                       outcome = missing_str$outcome,
                       group = 0)

# make a 5 fold CV, weight the observations so that there is a 50/50 mix of 0
# and 1 outcomes
set.seed(827)
index_outcome_0 <- sample(rep(1:5, each = 53))
index_outcome_1 <- sample(c(rep(1:5, each = 51), 1))
index_df[index_df$outcome == 0, "group"] <- index_outcome_0
index_df[index_df$outcome == 1, "group"] <- index_outcome_1

# for every group, make a prediction on the other data and use the data of this
# group as the test data
results <- lapply(c(1:5), function(i) {
  index_test <- index_df[index_df$group == i, "index"]
  index_train <- index_df[index_df$group != i, "index"]
  train_x <- data_merged[index_train, ]
  train_y <- index_df[index_train, "outcome"]
  test_x <- data_merged[index_test, ]
  test_y <- index_df[index_test, "outcome"]
  
  list_models <- list()
  print(paste0("start fold ", i))
  list_models[[1]] <- prioritylasso(X = train_x,
                                    Y = train_y,
                                    family = "binomial",
                                    type.measure = "auc",
                                    blocks = block_information,
                                    standardize = TRUE,
                                    mcontrol = missing.control(handle.missingdata = "ignore",
                                                               offset.firstblock = "zero"))
  print("done with model 1")
  list_models[[2]] <- prioritylasso(X = train_x,
                                    Y = train_y,
                                    family = "binomial",
                                    type.measure = "auc",
                                    blocks = block_information,
                                    standardize = TRUE,
                                    mcontrol = missing.control(handle.missingdata = "ignore",
                                                               offset.firstblock = "intercept"))
  print("done with model 2")
  list_models[[3]] <- prioritylasso(X = train_x,
                                    Y = train_y,
                                    family = "binomial",
                                    type.measure = "auc",
                                    blocks = block_information,
                                    standardize = TRUE,
                                    mcontrol = missing.control(handle.missingdata = "impute.offset",
                                                               impute.offset.cases = "available.cases",
                                                               select.available.cases = "maximise.blocks"))

  print("done with model 3")
  list_models[[4]] <- prioritylasso(X = train_x,
                                    Y = train_y,
                                    family = "binomial",
                                    type.measure = "auc",
                                    blocks = block_information,
                                    standardize = TRUE,
                                    mcontrol = missing.control(handle.missingdata = "impute.offset",
                                                               impute.offset.cases = "available.cases",
                                                               select.available.cases = "max"))
  print("done with model 4")
  
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
  
  list_models[[5]] <- mddsPLS(Xs = train_x_list,
                   Y = as.factor(train_y),
                   lambda = result_cv$Optim$optim_para_all$Lambdas[1],
                   R = result_cv$Optim$optim_para_all$R[1],
                   mode = "logit",
                   weight = TRUE)
  
  print("done with model 5")
  
  pred_ddspls <- predict(object = list_models[[5]],
                         newdata = test_x_list)$probY
  
  # in the current CV, there is one observation in the test set of block 1
  # for which prioritylasso can't make a prediction
  # -> therefore exclude this observation
  # if (i == 1) {
  #   test_x <- test_x[-c(15, 68), ]
  #   test_y <- test_y[-c(15, 68)]
  # }
  
  pred_value_list <- list()
  # blocks 1 # index 1-4
  pred_value_list[[1]] <- tryCatch({
    predict(object = list_models[[1]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "set.zero",
            include.allintercepts = TRUE,
            use.blocks = c(1))
  }, error = function(e) {
    print(e)
  })
  pred_value_list[[2]] <- tryCatch({
    predict(object = list_models[[2]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "set.zero",
            include.allintercepts = TRUE,
            use.blocks = c(1))
  }, error = function(e) {
    print(e)
  })
  pred_value_list[[3]] <- tryCatch({
    predict(object = list_models[[3]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "impute.block",
            include.allintercepts = TRUE,
            use.blocks = c(1))
  }, error = function(e) {
    print(e)
  })
  pred_value_list[[4]] <- tryCatch({
    predict(object = list_models[[4]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "impute.block",
            include.allintercepts = TRUE,
            use.blocks = c(1))
  }, error = function(e) {
    print(e)
  })
  print("done with block 1")
  # blocks 1,2 # index 5-8
  pred_value_list[[5]] <- tryCatch({
    predict(object = list_models[[1]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "set.zero",
            include.allintercepts = TRUE,
            use.blocks = c(1, 2))
  }, error = function(e) {
    print(e)
  })
  pred_value_list[[6]] <- tryCatch({
    predict(object = list_models[[2]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "set.zero",
            include.allintercepts = TRUE,
            use.blocks = c(1, 2))
  }, error = function(e) {
    print(e)
  })
  pred_value_list[[7]] <- tryCatch({
    predict(object = list_models[[3]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "impute.block",
            include.allintercepts = TRUE,
            use.blocks = c(1, 2))
  }, error = function(e) {
    print(e)
  })
  pred_value_list[[8]] <- tryCatch({
    predict(object = list_models[[4]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "impute.block",
            include.allintercepts = TRUE,
            use.blocks = c(1, 2))
  }, error = function(e) {
    print(e)
  })
  print("done with blocks 1,2")
  # blocks 1,2,3 # index 9-12
  pred_value_list[[9]] <- tryCatch({
    predict(object = list_models[[1]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "set.zero",
            include.allintercepts = TRUE,
            use.blocks = c(1, 2, 3))
  }, error = function(e) {
    print(e)
  })
  pred_value_list[[10]] <- tryCatch({
    predict(object = list_models[[2]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "set.zero",
            include.allintercepts = TRUE,
            use.blocks = c(1, 2, 3))
  }, error = function(e) {
    print(e)
  })
  pred_value_list[[11]] <- tryCatch({
    predict(object = list_models[[3]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "impute.block",
            include.allintercepts = TRUE,
            use.blocks = c(1, 2, 3))
  }, error = function(e) {
    print(e)
  })
  pred_value_list[[12]] <- tryCatch({
    predict(object = list_models[[4]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "impute.block",
            include.allintercepts = TRUE,
            use.blocks = c(1, 2, 3))
  }, error = function(e) {
    print(e)
  })
  print("done with blocks 1,2,3")
  # blocks 1,2,3,4 # index 13-16
  pred_value_list[[13]] <- tryCatch({
    predict(object = list_models[[1]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "set.zero",
            include.allintercepts = TRUE,
            use.blocks = c(1, 2, 3, 4))
  }, error = function(e) {
    print(e)
  })
  pred_value_list[[14]] <- tryCatch({
    predict(object = list_models[[2]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "set.zero",
            include.allintercepts = TRUE,
            use.blocks = c(1, 2, 3, 4))
  }, error = function(e) {
    print(e)
  })
  pred_value_list[[15]] <- tryCatch({
    predict(object = list_models[[3]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "impute.block",
            include.allintercepts = TRUE,
            use.blocks = c(1, 2, 3, 4))
  }, error = function(e) {
    print(e)
  })
  pred_value_list[[16]] <- tryCatch({
    predict(object = list_models[[4]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "impute.block",
            include.allintercepts = TRUE,
            use.blocks = c(1, 2, 3, 4))
  }, error = function(e) {
    print(e)
  })
  print("done with blocks 1,2,3,4")
  # blocks 1,2,3,4,5 # index 17-20
  pred_value_list[[17]] <- tryCatch({
    predict(object = list_models[[1]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "set.zero",
            include.allintercepts = TRUE,
            use.blocks = c(1, 2, 3, 4, 5))
  }, error = function(e) {
    print(e)
  })
  pred_value_list[[18]] <- tryCatch({
    predict(object = list_models[[2]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "set.zero",
            include.allintercepts = TRUE,
            use.blocks = c(1, 2, 3, 4, 5))
  }, error = function(e) {
    print(e)
  })
  pred_value_list[[19]] <- tryCatch({
    predict(object = list_models[[3]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "impute.block",
            include.allintercepts = TRUE,
            use.blocks = c(1, 2, 3, 4, 5))
  }, error = function(e) {
    print(e)
  })
  pred_value_list[[20]] <- tryCatch({
    predict(object = list_models[[4]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "impute.block",
            include.allintercepts = TRUE,
            use.blocks = c(1, 2, 3, 4, 5))
  }, error = function(e) {
    print(e)
  })
  print("done with blocks 1,2,3,4,5")
  
  # blocks 1,2,5 # index 21-24
  pred_value_list[[21]] <- tryCatch({
    predict(object = list_models[[1]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "set.zero",
            include.allintercepts = TRUE,
            use.blocks = c(1, 2, 5))
  }, error = function(e) {
    print(e)
  })
  pred_value_list[[22]] <- tryCatch({
    predict(object = list_models[[2]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "set.zero",
            include.allintercepts = TRUE,
            use.blocks = c(1, 2, 5))
  }, error = function(e) {
    print(e)
  })
  pred_value_list[[23]] <- tryCatch({
    predict(object = list_models[[3]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "impute.block",
            include.allintercepts = TRUE,
            use.blocks = c(1, 2, 5))
  }, error = function(e) {
    print(e)
  })
  pred_value_list[[24]] <- tryCatch({
    predict(object = list_models[[4]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "impute.block",
            include.allintercepts = TRUE,
            use.blocks = c(1, 2, 5))
  }, error = function(e) {
    print(e)
  })
  print("done with blocks 1,2,5")
  
  
  
  list(list_model = list_models,
       pred_value_list = pred_value_list,
       test_y = test_y,
       pred_ddspls = pred_ddspls)
})


# save the results
saveRDS(results, "../data/results_permute_blocks_small_1_2020_05_20.Rds")
