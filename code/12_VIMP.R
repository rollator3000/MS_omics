"
Script to extract the VariableImportance for the different approaches!
The importance of a variable is calculated according to the 'Permutation' of
the given variable of the out-of-bag examples!
"

library(randomForestSRC)


formula_  = as.formula("gender ~ .")

fittedmodel <- rfsrc(formula_, data = imputed_data$data[[1]]$train)


print(vimp(fittedmodel)$importance)
