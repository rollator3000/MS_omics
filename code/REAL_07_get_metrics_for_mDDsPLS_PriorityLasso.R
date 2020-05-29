"Script to calculate the Metrics of the Priority-Lasso / mdd-sPLS
  - only recived the predictions of the models!
"



# [1] 5_3_1 Results
# 1-1 Read in the results
# 1-1-1 Priority Lasso
res_PL <- readRDS("./docs/CV_Res/REAL/PriorityLasso/5_3_1/data/results_different_blocks_2020_05_06.Rds")

# res_PL[[2]]$list_model           --> info zum gefitteten Model
# res_PL[[2]]$pred_value_list      --> kp
# res_PL[[2]]$list_models_single   --> Info zu fitted models
res_PL[[2]]$pred_value_single    # --> predicted values
res_PL[[2]]$test_y               # --> true response values

# 1-1-2 MD-sPLS
res_MD <- readRDS("./docs/CV_Res/REAL/PriorityLasso/5_3_1/data/results_different_blocks_ddsPLS_weight_2020_05_20.Rds")

# res_MD[[2]]$model             --> Info to fitted model
res_MD[[2]]$test_y            # --> true response values
res_MD[[2]]$pred_value_ddspls # --> predicted probabilities 



