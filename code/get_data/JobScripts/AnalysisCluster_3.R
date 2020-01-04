# Set working directory:

setwd("~/")


# Make table of settings:

method <- c("BlockVarSel", "SplitWeights", "BlockForest", 
"RandomBlock", "VarProb", "randomsurvivalforest")
dat <- c("BRCA.Rda", "GBM.Rda", "KIRC.Rda", "KIRP.Rda", "LIHC.Rda", "STAD.Rda")
cvind <- 1:5

cvfoldind <- 1:5

scenariogrid <- expand.grid(cvind=cvind, dat=dat, cvfoldind=cvfoldind, method=method, stringsAsFactors = FALSE)
scenariogrid <- scenariogrid[,ncol(scenariogrid):1]

set.seed(1234)
seeds <- sample(1000:10000000, size=length(dat)*length(cvind))

scenariogrid$seed <- rep(seeds, times=length(method)*length(cvfoldind))


# Randomly permute rows of the table containing the settings.
# This is performed to ensure a comparable computational burden for
# the jobs to be performed in parallel:

set.seed(1234)
reorderind <- sample(1:nrow(scenariogrid))
scenariogrid <- scenariogrid[reorderind,]


# Save scenariogrid, needed in evaluation of the results:

save(scenariogrid, file="./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Results/scenariogrid3.Rda")


# Source the functions that are used in performing the calculations 
# on the cluster:

source("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Functions/Functions_AnalysisCluster.R")


# Start the cluster:

# NOTE: This syntax requires the use of the RMPISNOW script, see the README file
# contained in the root folder "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION".

cl <- makeCluster()


# Export the objects in the workspace to the
# parallel jobs:

clusterExport(cl, list=ls())


# Perform the calculations:

Results <- parLapply(cl, 1:nrow(scenariogrid), function(z)
  try({evaluatesetting(z)}))

  
# Stop the cluster:

stopCluster(cl)
