library(pROC)

# read in the results from the different block combinations
comb_1 <- readRDS("../data/results_permute_blocks_small_1_2020_05_20.Rds")
comb_2 <- readRDS("../data/results_permute_blocks_small_2_2020_05_20.Rds")
comb_3 <- readRDS("../data/results_permute_blocks_small_3_2020_05_20.Rds")
comb_4 <- readRDS("../data/results_permute_blocks_small_4_2020_05_20.Rds")


generate_auc_method_comp <- function(object) {

  true_values_list <- lapply(object, function(x) x$test_y)
  true_values <- do.call("c", true_values_list)
  
  # get the indices (all combinations from block 1-5)
  indices <- 1:20
  
  # calculate the AUC for ddsPLS
  predicted_values_ddspls_list <- lapply(object, function(x) {
    as.vector(x$pred_ddspls[, "1"])
  })
  
  predicted_values_ddspls <- do.call("c", predicted_values_ddspls_list)
  
  auc_ddspls <- roc(true_values, predicted_values_ddspls)$auc
  
  # initialise data.frame
  results_auc <- data.frame(ign_zero = rep(NA, 5),
                            ign_incp = rep(NA, 5),
                            imp_maxb = rep(NA, 5),
                            imp_maxn = rep(NA, 5))

  
  # calculate the AUC for the different methods and blocks for predictions
  for (i in indices) {
    predicted_values_list <- lapply(object, function(x) {
      x[["pred_value_list"]][[i]]
      
    })
    predicted_values <- do.call("c", predicted_values_list)
    roc_curve <- roc(true_values, predicted_values)
    
    column_index <- ifelse(i %% 4 == 0, 4, i %% 4)
    row_index <- ceiling(i / 4)
    results_auc[row_index, column_index] <- roc_curve$auc
  }

  
  # return the AUCs
  list(results_auc = results_auc,
       auc_ddspls = auc_ddspls)
}

comb_1_auc <- generate_auc_method_comp(comb_1)
print(xtable(t(comb_1_auc$results_auc)))
comb_1_auc$auc_ddspls
# ddspls 0.8653

comb_2_auc <- generate_auc_method_comp(comb_2)
print(xtable(t(comb_2_auc$results_auc)))
comb_2_auc$auc_ddspls
# ddspls 0.8969

comb_3_auc <- generate_auc_method_comp(comb_3)
print(xtable(t(comb_3_auc$results_auc)))
comb_3_auc$auc_ddspls
# ddspls 0.882

comb_4_auc <- generate_auc_method_comp(comb_4)
print(xtable(t(comb_4_auc$results_auc)))
comb_4_auc$auc_ddspls
# ddspls 0.8925