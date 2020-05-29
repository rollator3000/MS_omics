library(pROC)
library(RColorBrewer)

# load the results
results <- readRDS("../data/results_different_blocks_2020_05_06.Rds")
results_ddspls <- readRDS("../data/results_different_blocks_ddsPLS_weight_2020_05_20.Rds")

true_values_list <- lapply(results, function(x) x$test_y)
true_values <- do.call("c", true_values_list)

generate_roc_method_comp <- function(object = results,
                                     object_2 = results_ddspls,
                                     type = c("pred_value_list")) {
  type <- match.arg(type)
  
  # get the indices depending on the method
  indices <- 21:24
  
  # record the AUCs
  results_auc <- rep(NA, 5)
  
  # plot the first ROC curve
  # define a colour palette
  col_pal <- brewer.pal(n = 5, name = "Set1")
  # prioritylasso
  predicted_values_list <- lapply(object, function(x) {
    x[[type]][[indices[1]]]
  })
  predicted_values <- do.call("c", predicted_values_list)
  
  # ddsPLS
  predicted_values_ddspls_list <- lapply(object_2, function(x) {
    as.vector(x$pred_value_ddspls[, "1"])
  })
  
  predicted_values_ddspls <- do.call("c", predicted_values_ddspls_list)
  
  roc_curve <- roc(true_values, predicted_values_ddspls)
  results_auc[1] <- roc_curve$auc
  plot(roc_curve, col = col_pal[1])
  
  # plot the other curves into this plot
  for (i in 1:length(indices)) {
    predicted_values_list <- lapply(object, function(x) {
      x[[type]][[indices[i]]]
      
    })
    predicted_values <- do.call("c", predicted_values_list)
    roc_curve <- roc(true_values, predicted_values)
    results_auc[i + 1] <- roc_curve$auc
    lines(roc_curve, col = col_pal[i + 1])
  }
  # add a legend
  legend("bottomright",
         legend = c("mdd-sPLS",
                    "ignore, zero", "ignore, intercept",
                    "impute, maximise blocks", "impute, max. n"),
         col = col_pal, lwd = 2)
  
  # return the AUCs
  results_auc
}

# generate the plot
pdf("../plots/compare_5_methods.pdf")
generate_roc_method_comp()
dev.off()
# AUCs: 0.8780642 0.6193731 0.6469098 0.8756824 0.8742980