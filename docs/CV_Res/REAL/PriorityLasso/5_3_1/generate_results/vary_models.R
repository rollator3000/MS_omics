# load the needed libraries
devtools::load_all("../../prioritylasso/")
library(pROC)
library(dplyr)


# read in the data
load("../data/data 12032020.RData")

# include the NAs of the missing persons
df1_mod <- full_join(df1, y, by = "id") %>%
  select(-outcome) %>%
  arrange(id) %>%
  select(-id)
df2_mod <- full_join(df2, y, by = "id") %>%
  select(-outcome) %>%
  arrange(id) %>%
  select(-id)
df3_mod <- full_join(df3, y, by = "id") %>%
  select(-outcome) %>%
  arrange(id) %>%
  select(-id)
df4_mod <- full_join(df4, y, by = "id") %>%
  select(-outcome) %>%
  arrange(id) %>%
  select(-id)
df51_mod <- full_join(df51, y, by = "id") %>%
  select(-outcome) %>%
  arrange(id) %>%
  select(-id)
df53_mod <- full_join(df53, y, by = "id") %>%
  select(-outcome) %>%
  arrange(id) %>%
  select(-id)

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
  
  list_models <- list()
  
  list_models[[1]] <- prioritylasso(X = train_x,
                           Y = train_y,
                           family = "binomial",
                           type.measure = "auc",
                           blocks = block_information,
                           standardize = TRUE,
                           mcontrol = missing.control(handle.missingdata = "ignore",
                                                      offset.firstblock = "zero"))
  list_models[[2]] <- prioritylasso(X = train_x,
                           Y = train_y,
                           family = "binomial",
                           type.measure = "auc",
                           blocks = block_information,
                           standardize = TRUE,
                           mcontrol = missing.control(handle.missingdata = "ignore",
                                                      offset.firstblock = "intercept"))
  list_models[[3]] <- prioritylasso(X = train_x,
                           Y = train_y,
                           family = "binomial",
                           type.measure = "auc",
                           blocks = block_information,
                           standardize = TRUE,
                           mcontrol = missing.control(handle.missingdata = "impute.offset",
                                                      impute.offset.cases = "available.cases",
                                                      select.available.cases = "maximise.blocks"))
  list_models[[4]] <- prioritylasso(X = train_x,
                           Y = train_y,
                           family = "binomial",
                           type.measure = "auc",
                           blocks = block_information,
                           standardize = TRUE,
                           mcontrol = missing.control(handle.missingdata = "impute.offset",
                                                      impute.offset.cases = "available.cases",
                                                      select.available.cases = "max"))
  print(paste0("fold ", i))
  # in the current CV, there is one observation in the test set of block 1
  # for which prioritylasso can't make a prediction
  # -> therefore exclude this observation
  if (i == 1) {
    test_x <- test_x[-c(15, 68), ]
    test_y <- test_y[-c(15, 68)]
  }
  
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
  
  
  # blocks 1,2,3,4,5,6 # index 21-24
  pred_value_list[[21]] <- tryCatch({
    predict(object = list_models[[1]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "set.zero",
            include.allintercepts = TRUE)
  }, error = function(e) {
    print(e)
  })
  pred_value_list[[22]] <- tryCatch({
    predict(object = list_models[[2]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "set.zero",
            include.allintercepts = TRUE)
  }, error = function(e) {
    print(e)
  })
  pred_value_list[[23]] <- tryCatch({
    predict(object = list_models[[3]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "impute.block",
            include.allintercepts = TRUE)
  }, error = function(e) {
    print(e)
  })
  pred_value_list[[24]] <- tryCatch({
    predict(object = list_models[[4]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "impute.block",
            include.allintercepts = TRUE)
  }, error = function(e) {
    print(e)
  })
  
  # blocks 1,2,6 # index 25-28
  pred_value_list[[25]] <- tryCatch({
    predict(object = list_models[[1]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "set.zero",
            include.allintercepts = TRUE,
            use.blocks = c(1, 2, 6))
  }, error = function(e) {
    print(e)
  })
  pred_value_list[[26]] <- tryCatch({
    predict(object = list_models[[2]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "set.zero",
            include.allintercepts = TRUE,
            use.blocks = c(1, 2, 6))
  }, error = function(e) {
    print(e)
  })
  pred_value_list[[27]] <- tryCatch({
    predict(object = list_models[[3]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "impute.block",
            include.allintercepts = TRUE,
            use.blocks = c(1, 2, 6))
  }, error = function(e) {
    print(e)
  })
  pred_value_list[[28]] <- tryCatch({
    predict(object = list_models[[4]],
            newdata = test_x,
            type = "response",
            handle.missingtestdata = "impute.block",
            include.allintercepts = TRUE,
            use.blocks = c(1, 2, 6))
  }, error = function(e) {
    print(e)
  })
  
  ##############################################################################
  list_models_single <- list()
  
  pred_value_single <- list()
  # blocks 1 # index 1-4
  list_models_single[[1]] <- prioritylasso(X = train_x[, unlist(block_information[c(1)])],
                                    Y = train_y,
                                    family = "binomial",
                                    type.measure = "auc",
                                    blocks = block_information[c(1)],
                                    standardize = TRUE,
                                    mcontrol = missing.control(handle.missingdata = "ignore",
                                                               offset.firstblock = "zero"))
  list_models_single[[2]] <- prioritylasso(X = train_x[, unlist(block_information[c(1)])],
                                    Y = train_y,
                                    family = "binomial",
                                    type.measure = "auc",
                                    blocks = block_information[c(1)],
                                    standardize = TRUE,
                                    mcontrol = missing.control(handle.missingdata = "ignore",
                                                               offset.firstblock = "intercept"))
  list_models_single[[3]] <- prioritylasso(X = train_x[, unlist(block_information[c(1)])],
                                    Y = train_y,
                                    family = "binomial",
                                    type.measure = "auc",
                                    blocks = block_information[c(1)],
                                    standardize = TRUE,
                                    mcontrol = missing.control(handle.missingdata = "impute.offset",
                                                               impute.offset.cases = "available.cases",
                                                               select.available.cases = "maximise.blocks"))
  list_models_single[[4]] <- prioritylasso(X = train_x[, unlist(block_information[c(1)])],
                                    Y = train_y,
                                    family = "binomial",
                                    type.measure = "auc",
                                    blocks = block_information[c(1)],
                                    standardize = TRUE,
                                    mcontrol = missing.control(handle.missingdata = "impute.offset",
                                                               impute.offset.cases = "available.cases",
                                                               select.available.cases = "max"))
  pred_value_single[[1]] <- tryCatch({
    predict(object = list_models_single[[1]],
            newdata = test_x[, unlist(block_information[c(1)])],
            type = "response",
            handle.missingtestdata = "none",
            include.allintercepts = TRUE)
  }, error = function(e) {
    print(e)
  })
  pred_value_single[[2]] <- tryCatch({
    predict(object = list_models_single[[2]],
            newdata = test_x[, unlist(block_information[c(1)])],
            type = "response",
            handle.missingtestdata = "none",
            include.allintercepts = TRUE)
  }, error = function(e) {
    print(e)
  })
  pred_value_single[[3]] <- tryCatch({
    predict(object = list_models_single[[3]],
            newdata = test_x[, unlist(block_information[c(1)])],
            type = "response",
            handle.missingtestdata = "none",
            include.allintercepts = TRUE)
  }, error = function(e) {
    print(e)
  })
  pred_value_single[[4]] <- tryCatch({
    predict(object = list_models_single[[4]],
            newdata = test_x[, unlist(block_information[c(1)])],
            type = "response",
            handle.missingtestdata = "none",
            include.allintercepts = TRUE)
  }, error = function(e) {
    print(e)
  })
  
  # blocks 1,2 # index 5-8
  list_models_single[[5]] <- prioritylasso(X = train_x[, unlist(block_information[c(1, 2)])],
                                           Y = train_y,
                                           family = "binomial",
                                           type.measure = "auc",
                                           blocks = block_information[c(1, 2)],
                                           standardize = TRUE,
                                           mcontrol = missing.control(handle.missingdata = "ignore",
                                                                      offset.firstblock = "zero"))
  list_models_single[[6]] <- prioritylasso(X = train_x[, unlist(block_information[c(1, 2)])],
                                           Y = train_y,
                                           family = "binomial",
                                           type.measure = "auc",
                                           blocks = block_information[c(1, 2)],
                                           standardize = TRUE,
                                           mcontrol = missing.control(handle.missingdata = "ignore",
                                                                      offset.firstblock = "intercept"))
  list_models_single[[7]] <- prioritylasso(X = train_x[, unlist(block_information[c(1, 2)])],
                                           Y = train_y,
                                           family = "binomial",
                                           type.measure = "auc",
                                           blocks = block_information[c(1, 2)],
                                           standardize = TRUE,
                                           mcontrol = missing.control(handle.missingdata = "impute.offset",
                                                                      impute.offset.cases = "available.cases",
                                                                      select.available.cases = "maximise.blocks"))
  list_models_single[[8]] <- prioritylasso(X = train_x[, unlist(block_information[c(1, 2)])],
                                           Y = train_y,
                                           family = "binomial",
                                           type.measure = "auc",
                                           blocks = block_information[c(1, 2)],
                                           standardize = TRUE,
                                           mcontrol = missing.control(handle.missingdata = "impute.offset",
                                                                      impute.offset.cases = "available.cases",
                                                                      select.available.cases = "max"))
  # in the 5th fold, there are no missing values in the test data
  if (i == 5) {
    pred_value_single[[5]] <- tryCatch({
      predict(object = list_models_single[[5]],
              newdata = test_x[, unlist(block_information[c(1, 2)])],
              type = "response",
              handle.missingtestdata = "none",
              include.allintercepts = TRUE)
    }, error = function(e) {
      print(e)
    })
    pred_value_single[[6]] <- tryCatch({
      predict(object = list_models_single[[6]],
              newdata = test_x[, unlist(block_information[c(1, 2)])],
              type = "response",
              handle.missingtestdata = "none",
              include.allintercepts = TRUE)
    }, error = function(e) {
      print(e)
    })
    pred_value_single[[7]] <- tryCatch({
      predict(object = list_models_single[[7]],
              newdata = test_x[, unlist(block_information[c(1, 2)])],
              type = "response",
              handle.missingtestdata = "none",
              include.allintercepts = TRUE)
    }, error = function(e) {
      print(e)
    })
    pred_value_single[[8]] <- tryCatch({
      predict(object = list_models_single[[8]],
              newdata = test_x[, unlist(block_information[c(1, 2)])],
              type = "response",
              handle.missingtestdata = "none",
              include.allintercepts = TRUE)
    }, error = function(e) {
      print(e)
    })
  } else {
    pred_value_single[[5]] <- tryCatch({
      predict(object = list_models_single[[5]],
              newdata = test_x[, unlist(block_information[c(1, 2)])],
              type = "response",
              handle.missingtestdata = "set.zero",
              include.allintercepts = TRUE)
    }, error = function(e) {
      print(e)
    })
    pred_value_single[[6]] <- tryCatch({
      predict(object = list_models_single[[6]],
              newdata = test_x[, unlist(block_information[c(1, 2)])],
              type = "response",
              handle.missingtestdata = "set.zero",
              include.allintercepts = TRUE)
    }, error = function(e) {
      print(e)
    })
    pred_value_single[[7]] <- tryCatch({
      predict(object = list_models_single[[7]],
              newdata = test_x[, unlist(block_information[c(1, 2)])],
              type = "response",
              handle.missingtestdata = "impute.block",
              include.allintercepts = TRUE)
    }, error = function(e) {
      print(e)
    })
    pred_value_single[[8]] <- tryCatch({
      predict(object = list_models_single[[8]],
              newdata = test_x[, unlist(block_information[c(1, 2)])],
              type = "response",
              handle.missingtestdata = "impute.block",
              include.allintercepts = TRUE)
    }, error = function(e) {
      print(e)
    })
  }

  
  # blocks 1,2,3 # index 9-12
  list_models_single[[9]] <- prioritylasso(X = train_x[, unlist(block_information[c(1, 2, 3)])],
                                           Y = train_y,
                                           family = "binomial",
                                           type.measure = "auc",
                                           blocks = block_information[c(1, 2, 3)],
                                           standardize = TRUE,
                                           mcontrol = missing.control(handle.missingdata = "ignore",
                                                                      offset.firstblock = "zero"))
  list_models_single[[10]] <- prioritylasso(X = train_x[, unlist(block_information[c(1, 2, 3)])],
                                           Y = train_y,
                                           family = "binomial",
                                           type.measure = "auc",
                                           blocks = block_information[c(1, 2, 3)],
                                           standardize = TRUE,
                                           mcontrol = missing.control(handle.missingdata = "ignore",
                                                                      offset.firstblock = "intercept"))
  list_models_single[[11]] <- prioritylasso(X = train_x[, unlist(block_information[c(1, 2, 3)])],
                                           Y = train_y,
                                           family = "binomial",
                                           type.measure = "auc",
                                           blocks = block_information[c(1, 2, 3)],
                                           standardize = TRUE,
                                           mcontrol = missing.control(handle.missingdata = "impute.offset",
                                                                      impute.offset.cases = "available.cases",
                                                                      select.available.cases = "maximise.blocks"))
  list_models_single[[12]] <- prioritylasso(X = train_x[, unlist(block_information[c(1, 2, 3)])],
                                           Y = train_y,
                                           family = "binomial",
                                           type.measure = "auc",
                                           blocks = block_information[c(1, 2, 3)],
                                           standardize = TRUE,
                                           mcontrol = missing.control(handle.missingdata = "impute.offset",
                                                                      impute.offset.cases = "available.cases",
                                                                      select.available.cases = "max"))
  pred_value_single[[9]] <- tryCatch({
    predict(object = list_models_single[[9]],
            newdata = test_x[, unlist(block_information[c(1, 2, 3)])],
            type = "response",
            handle.missingtestdata = "set.zero",
            include.allintercepts = TRUE)
  }, error = function(e) {
    print(e)
  })
  pred_value_single[[10]] <- tryCatch({
    predict(object = list_models_single[[10]],
            newdata = test_x[, unlist(block_information[c(1, 2, 3)])],
            type = "response",
            handle.missingtestdata = "set.zero",
            include.allintercepts = TRUE)
  }, error = function(e) {
    print(e)
  })
  pred_value_single[[11]] <- tryCatch({
    predict(object = list_models_single[[11]],
            newdata = test_x[, unlist(block_information[c(1, 2, 3)])],
            type = "response",
            handle.missingtestdata = "impute.block",
            include.allintercepts = TRUE)
  }, error = function(e) {
    print(e)
  })
  pred_value_single[[12]] <- tryCatch({
    predict(object = list_models_single[[12]],
            newdata = test_x[, unlist(block_information[c(1, 2, 3)])],
            type = "response",
            handle.missingtestdata = "impute.block",
            include.allintercepts = TRUE)
  }, error = function(e) {
    print(e)
  })
  
  # blocks 1,2,3,4 # index 13-16
  list_models_single[[13]] <- prioritylasso(X = train_x[, unlist(block_information[c(1, 2, 3, 4)])],
                                           Y = train_y,
                                           family = "binomial",
                                           type.measure = "auc",
                                           blocks = block_information[c(1, 2, 3, 4)],
                                           standardize = TRUE,
                                           mcontrol = missing.control(handle.missingdata = "ignore",
                                                                      offset.firstblock = "zero"))
  list_models_single[[14]] <- prioritylasso(X = train_x[, unlist(block_information[c(1, 2, 3, 4)])],
                                           Y = train_y,
                                           family = "binomial",
                                           type.measure = "auc",
                                           blocks = block_information[c(1, 2, 3, 4)],
                                           standardize = TRUE,
                                           mcontrol = missing.control(handle.missingdata = "ignore",
                                                                      offset.firstblock = "intercept"))
  list_models_single[[15]] <- prioritylasso(X = train_x[, unlist(block_information[c(1, 2, 3, 4)])],
                                           Y = train_y,
                                           family = "binomial",
                                           type.measure = "auc",
                                           blocks = block_information[c(1, 2, 3, 4)],
                                           standardize = TRUE,
                                           mcontrol = missing.control(handle.missingdata = "impute.offset",
                                                                      impute.offset.cases = "available.cases",
                                                                      select.available.cases = "maximise.blocks"))
  list_models_single[[16]] <- prioritylasso(X = train_x[, unlist(block_information[c(1, 2, 3, 4)])],
                                           Y = train_y,
                                           family = "binomial",
                                           type.measure = "auc",
                                           blocks = block_information[c(1, 2, 3, 4)],
                                           standardize = TRUE,
                                           mcontrol = missing.control(handle.missingdata = "impute.offset",
                                                                      impute.offset.cases = "available.cases",
                                                                      select.available.cases = "max"))
  pred_value_single[[13]] <- tryCatch({
    predict(object = list_models_single[[13]],
            newdata = test_x[, unlist(block_information[c(1, 2, 3, 4)])],
            type = "response",
            handle.missingtestdata = "set.zero",
            include.allintercepts = TRUE)
  }, error = function(e) {
    print(e)
  })
  pred_value_single[[14]] <- tryCatch({
    predict(object = list_models_single[[14]],
            newdata = test_x[, unlist(block_information[c(1, 2, 3, 4)])],
            type = "response",
            handle.missingtestdata = "set.zero",
            include.allintercepts = TRUE)
  }, error = function(e) {
    print(e)
  })
  pred_value_single[[15]] <- tryCatch({
    predict(object = list_models_single[[15]],
            newdata = test_x[, unlist(block_information[c(1, 2, 3, 4)])],
            type = "response",
            handle.missingtestdata = "impute.block",
            include.allintercepts = TRUE)
  }, error = function(e) {
    print(e)
  })
  pred_value_single[[16]] <- tryCatch({
    predict(object = list_models_single[[16]],
            newdata = test_x[, unlist(block_information[c(1, 2, 3, 4)])],
            type = "response",
            handle.missingtestdata = "impute.block",
            include.allintercepts = TRUE)
  }, error = function(e) {
    print(e)
  })
  
  # blocks 1,2,3,4,5 # index 17-20
  list_models_single[[17]] <- prioritylasso(X = train_x[, unlist(block_information[c(1, 2, 3, 4, 5)])],
                                           Y = train_y,
                                           family = "binomial",
                                           type.measure = "auc",
                                           blocks = block_information[c(1, 2, 3, 4, 5)],
                                           standardize = TRUE,
                                           mcontrol = missing.control(handle.missingdata = "ignore",
                                                                      offset.firstblock = "zero"))
  list_models_single[[18]] <- prioritylasso(X = train_x[, unlist(block_information[c(1, 2, 3, 4, 5)])],
                                           Y = train_y,
                                           family = "binomial",
                                           type.measure = "auc",
                                           blocks = block_information[c(1, 2, 3, 4, 5)],
                                           standardize = TRUE,
                                           mcontrol = missing.control(handle.missingdata = "ignore",
                                                                      offset.firstblock = "intercept"))
  list_models_single[[19]] <- prioritylasso(X = train_x[, unlist(block_information[c(1, 2, 3, 4, 5)])],
                                           Y = train_y,
                                           family = "binomial",
                                           type.measure = "auc",
                                           blocks = block_information[c(1, 2, 3, 4, 5)],
                                           standardize = TRUE,
                                           mcontrol = missing.control(handle.missingdata = "impute.offset",
                                                                      impute.offset.cases = "available.cases",
                                                                      select.available.cases = "maximise.blocks"))
  list_models_single[[20]] <- prioritylasso(X = train_x[, unlist(block_information[c(1, 2, 3, 4, 5)])],
                                           Y = train_y,
                                           family = "binomial",
                                           type.measure = "auc",
                                           blocks = block_information[c(1, 2, 3, 4, 5)],
                                           standardize = TRUE,
                                           mcontrol = missing.control(handle.missingdata = "impute.offset",
                                                                      impute.offset.cases = "available.cases",
                                                                      select.available.cases = "max"))
  pred_value_single[[17]] <- tryCatch({
    predict(object = list_models_single[[17]],
            newdata = test_x[, unlist(block_information[c(1, 2, 3, 4, 5)])],
            type = "response",
            handle.missingtestdata = "set.zero",
            include.allintercepts = TRUE)
  }, error = function(e) {
    print(e)
  })
  pred_value_single[[18]] <- tryCatch({
    predict(object = list_models_single[[18]],
            newdata = test_x[, unlist(block_information[c(1, 2, 3, 4, 5)])],
            type = "response",
            handle.missingtestdata = "set.zero",
            include.allintercepts = TRUE)
  }, error = function(e) {
    print(e)
  })
  pred_value_single[[19]] <- tryCatch({
    predict(object = list_models_single[[19]],
            newdata = test_x[, unlist(block_information[c(1, 2, 3, 4, 5)])],
            type = "response",
            handle.missingtestdata = "impute.block",
            include.allintercepts = TRUE)
  }, error = function(e) {
    print(e)
  })
  pred_value_single[[20]] <- tryCatch({
    predict(object = list_models_single[[20]],
            newdata = test_x[, unlist(block_information[c(1, 2, 3, 4, 5)])],
            type = "response",
            handle.missingtestdata = "impute.block",
            include.allintercepts = TRUE)
  }, error = function(e) {
    print(e)
  })
  
  
  # blocks 1,2,3,4,5,6 # index 21-24
  # doesn't need to be redone, can use the predictions from above
  pred_value_single[[21]] <- NULL
  pred_value_single[[22]] <- NULL
  pred_value_single[[23]] <- NULL
  pred_value_single[[24]] <- NULL
  list_models_single[[21]] <- NULL
  list_models_single[[22]] <- NULL
  list_models_single[[23]] <- NULL
  list_models_single[[24]] <- NULL
  
  # blocks 1,2,6 # index 25-28
  # adjust the block_information because the data (X) has now less columns,
  # therefore the original block information are invald
  # block 2 ends at 60, and block 6 has 84 entries, therefore the new boundaries
  # are:
  # 61-144
  list_models_single[[25]] <- prioritylasso(X = train_x[, unlist(block_information[c(1, 2, 6)])],
                                           Y = train_y,
                                           family = "binomial",
                                           type.measure = "auc",
                                           blocks = c(block_information[c(1, 2)], list(61:144)),
                                           standardize = TRUE,
                                           mcontrol = missing.control(handle.missingdata = "ignore",
                                                                      offset.firstblock = "zero"))
  list_models_single[[26]] <- prioritylasso(X = train_x[, unlist(block_information[c(1, 2, 6)])],
                                           Y = train_y,
                                           family = "binomial",
                                           type.measure = "auc",
                                           blocks = c(block_information[c(1, 2)], list(61:144)),
                                           standardize = TRUE,
                                           mcontrol = missing.control(handle.missingdata = "ignore",
                                                                      offset.firstblock = "intercept"))
  list_models_single[[27]] <- prioritylasso(X = train_x[, unlist(block_information[c(1, 2, 6)])],
                                           Y = train_y,
                                           family = "binomial",
                                           type.measure = "auc",
                                           blocks = c(block_information[c(1, 2)], list(61:144)),
                                           standardize = TRUE,
                                           mcontrol = missing.control(handle.missingdata = "impute.offset",
                                                                      impute.offset.cases = "available.cases",
                                                                      select.available.cases = "maximise.blocks"))
  list_models_single[[28]] <- prioritylasso(X = train_x[, unlist(block_information[c(1, 2, 6)])],
                                           Y = train_y,
                                           family = "binomial",
                                           type.measure = "auc",
                                           blocks = c(block_information[c(1, 2)], list(61:144)),
                                           standardize = TRUE,
                                           mcontrol = missing.control(handle.missingdata = "impute.offset",
                                                                      impute.offset.cases = "available.cases",
                                                                      select.available.cases = "max"))
  pred_value_single[[25]] <- tryCatch({
    predict(object = list_models_single[[25]],
            newdata = test_x[, unlist(block_information[c(1, 2, 6)])],
            type = "response",
            handle.missingtestdata = "set.zero",
            include.allintercepts = TRUE)
  }, error = function(e) {
    print(e)
  })
  pred_value_single[[26]] <- tryCatch({
    predict(object = list_models_single[[26]],
            newdata = test_x[, unlist(block_information[c(1, 2, 6)])],
            type = "response",
            handle.missingtestdata = "set.zero",
            include.allintercepts = TRUE)
  }, error = function(e) {
    print(e)
  })
  pred_value_single[[27]] <- tryCatch({
    predict(object = list_models_single[[27]],
            newdata = test_x[, unlist(block_information[c(1, 2, 6)])],
            type = "response",
            handle.missingtestdata = "impute.block",
            include.allintercepts = TRUE)
  }, error = function(e) {
    print(e)
  })
  pred_value_single[[28]] <- tryCatch({
    predict(object = list_models_single[[28]],
            newdata = test_x[, unlist(block_information[c(1, 2, 6)])],
            type = "response",
            handle.missingtestdata = "impute.block",
            include.allintercepts = TRUE)
  }, error = function(e) {
    print(e)
  })
  
  
  list(list_model = list_models,
       pred_value_list = pred_value_list,
       list_models_single = list_models_single,
       pred_value_single = pred_value_single,
       test_y = test_y)
})


# save the results
saveRDS(results, "../data/results_different_blocks_2020_05_06.Rds")
