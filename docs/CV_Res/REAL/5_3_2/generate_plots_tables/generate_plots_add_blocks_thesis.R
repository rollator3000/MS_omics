library(pROC)
library(RColorBrewer)

# load the results
results <- readRDS("../data/results_different_blocks_2020_05_06.Rds")

true_values_list <- lapply(results, function(x) x$test_y)
true_values <- do.call("c", true_values_list)

generate_plot_add_blocks <- function(object = results,
                                     type = c("pred_value_list", "pred_value_single"),
                                     method) {
  type <- match.arg(type)
  
  # get the indices depending on the method
  indices <- c(4, 8, 20, 24) + method
  
  # record the AUCs
  results_auc <- rep(NA, 4)
  
  # plot the first ROC curve
  # define a colour palette
  col_pal <- brewer.pal(n = 4, name = "Set2")
  predicted_values_list <- lapply(object, function(x) {
    x[[type]][[indices[1]]]
  })
  predicted_values <- do.call("c", predicted_values_list)
  roc_curve <- roc(true_values, predicted_values)
  results_auc[1] <- roc_curve$auc
  plot(roc_curve, col = col_pal[1])
  
  # plot the other curves into this plot
  for (i in 2:length(indices)) {
    predicted_values_list <- lapply(object, function(x) {
      if (i != 3) {
        x[[type]][[indices[i]]]
      } else {
        # for all blocks, "pred_value_single" doesn't have its own model
        x[["pred_value_list"]][[indices[i]]]
      }
    })
    predicted_values <- do.call("c", predicted_values_list)
    roc_curve <- roc(true_values, predicted_values)
    results_auc[i] <- roc_curve$auc
    lines(roc_curve, col = col_pal[i])
  }
  # add a legend
  legend("bottomright",
         legend = c("blocks 1, 2", "blocks 1, 2, 3",
                    "blocks 1, 2, 3, 4, 5, 6", "blocks 1, 2, 6"),
         col = col_pal, lwd = 2)
  
  # return the AUCs
  results_auc
}

# this one is in the thesis
# impute, max observations
pdf("../plots/model_all_blocks_add_blocks_impute_max_observations.pdf")
generate_plot_add_blocks(method = 4)
dev.off()
# AUCs
# 0.9032239 0.9108453 0.8742980 0.8488932