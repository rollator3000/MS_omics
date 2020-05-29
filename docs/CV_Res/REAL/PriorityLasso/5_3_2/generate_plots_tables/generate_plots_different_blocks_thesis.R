library(pROC)
library(RColorBrewer)

# load the results
results <- readRDS("../data/results_different_blocks_2020_05_06.Rds")
results_ddspls <- readRDS("../data/results_different_blocks_ddsPLS_weight_2020_05_20.Rds")

true_values_list <- lapply(results, function(x) x$test_y)
true_values <- do.call("c", true_values_list)

generate_roc_ddspls <- function(object = results_ddspls) {
  predicted_values_list <- lapply(object, function(x) {
    x$pred_value_ddspls[, "1"]
  })
  predicted_values <- do.call("c", predicted_values_list)
  roc(true_values, predicted_values)$auc
}

generate_plot_one_method <- function(object = results,
                                     type = c("pred_value_list", "pred_value_single"),
                                     method,
                                     show_legend = TRUE,
                                     plot = TRUE) {
  type <- match.arg(type)
  
  # get the indices depending on the method
  indices <- seq(from = method, to = 20 + method, by = 4)
  
  # record the AUCs
  results_auc <- rep(NA, 6)
  
  # plot the first ROC curve
  # define a colour palette
  col_pal <- brewer.pal(n = 6, name = "Set2")
  predicted_values_list <- lapply(object, function(x) {
    x[[type]][[indices[1]]]
  })
  predicted_values <- do.call("c", predicted_values_list)
  roc_curve <- roc(true_values, predicted_values)
  results_auc[1] <- roc_curve$auc
  if (plot) {
  plot(roc_curve, col = col_pal[1])
  }
  
  # plot the other curves into this plot
  for (i in 2:length(indices)) {
    predicted_values_list <- lapply(object, function(x) {
      if (i != length(indices)) {
      x[[type]][[indices[i]]]
      } else {
        # for all blocks, "pred_value_single" doesn't have its own model
        x[["pred_value_list"]][[indices[i]]]
      }
    })
    predicted_values <- do.call("c", predicted_values_list)
    roc_curve <- roc(true_values, predicted_values)
    results_auc[i] <- roc_curve$auc
    if (plot) {
    lines(roc_curve, col = col_pal[i])
    }
  }
  # add a legend
  if (show_legend && plot) {
  legend("bottomright",
         legend = c("block 1", "block 1, 2", "blocks 1, 2, 3",
                    "blocks 1, 2, 3, 4", "blocks 1, 2, 3, 4, 5",
                    "blocks 1, 2, 3, 4, 5, 6"),
         col = col_pal, lwd = 2, bg = "white")
  }
  
  # return the AUCs
  results_auc
}

# generate the plot
# combine these plots to one plot
op <- par(no.readonly = TRUE)
par(mfrow=c(2,2))
generate_plot_one_method(method = 1, show_legend = FALSE)
generate_plot_one_method(method = 2, show_legend = FALSE)
generate_plot_one_method(method = 3, show_legend = FALSE)
generate_plot_one_method(method = 4, show_legend = TRUE)
par(op)

################################################################################
# generate the AUC for the tables

# plots when all blocks are used
# ignore, zero
generate_plot_one_method(method = 1, plot = FALSE)
# AUCs
# 0.8896152 0.9053930 0.9161937 0.9165354 0.6494874 0.6193731

# ignore, intercept
generate_plot_one_method(method = 2, plot = FALSE)
# AUCs
# 0.8919923 0.9073095 0.9162532 0.9156737 0.6459367 0.6469098

# impute, maximise blocks
generate_plot_one_method(method = 3, plot = FALSE)
# AUCs
# 0.8883970 0.9017234 0.9167137 0.9156143 0.8939394 0.8756824

# impute, max observations
generate_plot_one_method(method = 4, plot = FALSE)
# AUCs
# 0.8884415 0.9032239 0.9108453 0.9061209 0.8856188 0.8742980



# plots when only the limited number of blocks are used
# ignore, zero
generate_plot_one_method(method = 1, type = "pred_value_single", plot = FALSE)
# AUCs
# 0.8920963 0.9073986 0.9212450 0.9201159 0.6448076 0.6193731

# ignore, intercept
generate_plot_one_method(method = 2, type = "pred_value_single", plot = FALSE)
# AUCs
# 0.8951567 0.9106373 0.9164314 0.9230426 0.6587431 0.6469098

# impute, maximise blocks
generate_plot_one_method(method = 3, type = "pred_value_single", plot = FALSE)
# AUCs
# 0.8927797 0.9003268 0.9162383 0.9136830 0.8861090 0.8756824

# impute, max observations
generate_plot_one_method(method = 4, type = "pred_value_single", plot = FALSE)
# AUCs
# 0.8927797 0.9093597 0.9139355 0.9050346 0.8835004 0.8742980

################################################################################
# mdd-sPLS (uses all blocks)
generate_roc_ddspls()
# 0.8781