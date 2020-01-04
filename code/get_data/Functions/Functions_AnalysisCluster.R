# This function performs one iteration of the five times repeated 
# 5-fold cross-validation on a specific data set using a specific
# method for the analysis of the multi-blocks case.
#
# It takes the whole number 'iter', which corresponds to the iter-th line 
# of 'scenariogrid', which contains the necessary information
# on the iter-th setting.

evaluatesetting <- function(iter) {
  
  
  # Initiate lists in which the results
  # will be stored:
  
  ytrue <- list()
  ytruestatus <- list()
  riskpreds <- list()
  paramvalues <- list()
  mtrys <- list()
  
  
  # Obtain information for the iter-th setting:
  
  dat <- scenariogrid$dat[iter]
  seed <- scenariogrid$seed[iter]
  method <- scenariogrid$method[iter]
  
  cvind <- scenariogrid$cvind[iter]
  cvfoldind <- scenariogrid$cvfoldind[iter]
  
  
  # Load data set:
  
  load(paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/ProcessedData/", dat, sep=""))
  
  
  # Make covariate matrix and target variable:
  
  blocknames <- c("clin", 
                  "mirna",
                  "mutation",
                  "cnv",
                  "rna")
  
  blocknamesav <-  c("")
  
  if(class(try(ncol(clin)))!="try-error")
    blocknamesav <- c(blocknamesav, "clin")
  if(class(try(ncol(mirna)))!="try-error")
    blocknamesav <- c(blocknamesav, "mirna")
  if(class(try(ncol(mutation)))!="try-error")
    blocknamesav <- c(blocknamesav, "mutation")
  if(class(try(ncol(cnv)))!="try-error")
    blocknamesav <- c(blocknamesav, "cnv")
  if(class(try(ncol(rna)))!="try-error")
    blocknamesav <- c(blocknamesav, "rna")  
  
  blocknamesav <- blocknamesav[-1]
  
  eval(parse(text=paste("X <- cbind(", paste(blocknamesav, collapse=", "), ")", sep="")))
  
  eval(parse(text=paste("block <- rep(1:length(blocknamesav), times=c(", paste(paste("ncol(", blocknamesav, ")", sep=""), collapse=", "), "))", sep="")))
  
  block <- lapply(1:length(blocknamesav), function(x) which(block==x))
  
  y <- cbind(targetvar$time, targetvar$status)
    
  
  # Divide data set into training and test data:

  ncv <- 5
  set.seed(seed)

  cvdiv <- makeCVdiv(n=nrow(X), ncv=ncv)
  
  Xtrain <- X[cvdiv!=cvfoldind,]
  ytrain <- y[cvdiv!=cvfoldind,]
  
  Xtest <- X[cvdiv==cvfoldind,]
  ytest <- y[cvdiv==cvfoldind,]
  
  
  # Random survival forest:
  
  if(method=="randomsurvivalforest") {
    predobj <- randomsurvivalforestwrap(Xtrain=Xtrain, ytrain=ytrain, Xtest=Xtest, ytest=ytest, block=block)
    riskpreds[[1]] <- predobj$riskpred
    paramvalues[[1]] <- predobj$paramvalues
    mtrys[[1]] <- predobj$mtry
  }
  
  # Random survival forest variant:
  if(method %in% c("BlockVarSel", "SplitWeights", "BlockForest", "RandomBlock", "VarProb")) {
    predobj <- blockforestwrap(Xtrain=Xtrain, ytrain=ytrain, Xtest=Xtest, ytest=ytest, block=block, block.method=method)
    riskpreds[[1]] <- predobj$riskpred
    paramvalues[[1]] <- predobj$paramvalues
    mtrys[[1]] <- predobj$mtry
  }
  
  
  # Test data: true survival times and values of the status variable:
  
  ytrue[[1]] <- ytest[,1]
  ytruestatus[[1]] <- ytest[,2]
  
  
  # Combine results in list:
  
  res <- list(ytrue=ytrue, ytruestatus=ytruestatus, riskpreds=riskpreds, paramvalues=paramvalues, mtrys=mtrys, settingind=iter)
  
  
  # Save results in different folder, depending on which data set
  # was used (different data sets were considered in "AnalysisCluster_1.R",
  # "AnalysisCluster_2.R", and "AnalysisCluster_3.R"):
  
  if(dat %in% c("COAD.Rda", "LGG.Rda", "OV.Rda", "PAAD.Rda", "SKCM.Rda", "UCEC.Rda"))
    save(res, file=paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Results/Results1/res", iter, ".Rda", sep=""))
  
  if(dat %in% c("BLCA.Rda", "CESC.Rda", "ESCA.Rda", "HNSC.Rda", "LUAD.Rda", 
                "LUSC.Rda", "PRAD.Rda", "READ.Rda", "SARC.Rda"))
    save(res, file=paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Results/Results2/res", iter, ".Rda", sep=""))
  
  if(dat %in% c("BRCA.Rda", "GBM.Rda", "KIRC.Rda", "KIRP.Rda", "LIHC.Rda", "STAD.Rda"))
    save(res, file=paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Results/Results3/res", iter, ".Rda", sep=""))
  
}






# This function performs one iteration of the five times repeated 
# 5-fold cross-validation on a specific data set using a specific
# method for the analysis of the two-blocks case.
#
# It takes the whole number 'iter', which corresponds to the iter-th line 
# of 'scenariogrid', which contains the necessary information
# on the iter-th setting.


evaluatesettingtwoblocks <- function(iter) {
    
  
  # Initiate lists in which the results
  # will be stored:
  
  ytrue <- list()
  ytruestatus <- list()
  riskpreds <- list()
  paramvalues <- list()
  mtrys <- list()
  
  
  # Obtain information for the iter-th setting:
  
  dat <- scenariogrid$dat[iter]
  seed <- scenariogrid$seed[iter]
  method <- scenariogrid$method[iter]
  
  cvind <- scenariogrid$cvind[iter]
  cvfoldind <- scenariogrid$cvfoldind[iter]
  
  
  # Load data set:
  
  load(paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/ProcessedData/", dat, sep=""))
  
  
  # Make covariate matrix and target variable:
  
  X <- cbind(clin, rna)
  
  block <- rep(1:2, times=c(ncol(clin), ncol(rna)))
  
  block <- lapply(1:2, function(x) which(block==x))
  
  y <- cbind(targetvar$time, targetvar$status)
    
  
  # Divide data set into training and test data:

  ncv <- 5
  set.seed(seed)

  cvdiv <- makeCVdiv(n=nrow(X), ncv=ncv)

  Xtrain <- X[cvdiv!=cvfoldind,]
  ytrain <- y[cvdiv!=cvfoldind,]
  
  Xtest <- X[cvdiv==cvfoldind,]
  ytest <- y[cvdiv==cvfoldind,]
  
  
  # Random survival forest:
  
  if(method=="randomsurvivalforest") {
    predobj <- randomsurvivalforestwrap(Xtrain=Xtrain, ytrain=ytrain, Xtest=Xtest, ytest=ytest, block=block)
    riskpreds[[1]] <- predobj$riskpred
    paramvalues[[1]] <- predobj$paramvalues
    mtrys[[1]] <- predobj$mtry
  }
  
  # Random survival forest variant:
  
  if(method %in% c("BlockVarSel", "SplitWeights", "BlockForest", "RandomBlock", "VarProb")) {
    predobj <- blockforestwrap(Xtrain=Xtrain, ytrain=ytrain, Xtest=Xtest, ytest=ytest, block=block, block.method=method)
    riskpreds[[1]] <- predobj$riskpred
    paramvalues[[1]] <- predobj$paramvalues
    mtrys[[1]] <- predobj$mtry
  }
  
  
  # Test data: true survival times and values of the status variable:
  
  ytrue[[1]] <- ytest[,1]
  ytruestatus[[1]] <- ytest[,2]
  
  
  # Combine results in list:
  
  res <- list(ytrue=ytrue, ytruestatus=ytruestatus, riskpreds=riskpreds, paramvalues=paramvalues, mtrys=mtrys, settingind=iter)
  
  
  # Save results:
  
  save(res, file=paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Results/ResultsTwoBlocks/res", iter, ".Rda", sep=""))

}






# Function to generate the splittings for stratified cross-validation:

# Input parameters:

# n    - number of observations in the data set
# ncv  - number of folds to use

makeCVdiv <- function(n, ncv) {
  
  nperfold <- ceiling(n/ncv)
  partition <- sample(rep(1:ncv, nperfold), n)
  
  partition
  
}






# Wrapper function that performs supervised variable selection 
# using univariate Cox regressions on training data, subsequently
# trains a random survival forest variant on the training data and, lastly,
# obtains predictions on test data using the random survival forest
# variant trained on the training data.

# Input parameters:

# Xtrain        - covariate matrix of the training data set
# ytrain        - target variable of the training data set. Matrix with two 
#                 columns, where the first column contains the vector of 
#                 survival/censoring times (one for each observation) and the 
#                 second column contains the status variable, that has the 
#                 value '1' if the corresponding time is a survival time and 
#                 '0' if that time is a censoring time.
# Xtest         - covariate matrix of the test data set
# ytest         - target variable in the test data set. Matrix with two
#                 columns, where the first column contains the vector of 
#                 survival/censoring times (one for each observation) and the 
#                 second column contains the status variable, that has the 
#                 value '1' if the corresponding time is a survival time and 
#                 '0' if that time is a censoring time.
# block         - A list of length equal to the number of blocks considered. 
#                 Each entry contains the vector of column indices in 'Xtrain' 
#                 of the covariates in one of the blocks.
# block.method  - Forest variant to use. One of the following: "BlockForest", 
#                 "RandomBlock", "BlockVarSel", "VarProb", "SplitWeights".

blockforestwrap <- function(Xtrain, ytrain, Xtest, ytest, block, block.method) {
  
  library("survival")
  library("blockForest")
  
  # Perform supervised variable selection per block using univariate
  # Cox regressions:
  
  indsel <- c()
  blocksub <- c()
  
  for(i in seq(along=block)) {
    if(length(block[[i]]) <= 2500) {
      indsel <- c(indsel, block[[i]])
      blocksub <- c(blocksub, rep(i, length(block[[i]])))
    }
    else {
      varinds <- block[[i]]
      pvals <- apply(Xtrain[,varinds], 2, function(y2) summary(coxph(Surv(ytrain[,1], ytrain[,2]) ~ y2))$coef[,"Pr(>|z|)"][1])
      indsel <- c(indsel, varinds[order(pvals)[1:2500]])
      blocksub <- c(blocksub, rep(i, 2500))
    }
  }
  
  blocksub <- blocksub[order(indsel)]
  indsel <- indsel[order(indsel)]
  
  Xtrain <- Xtrain[,indsel]
  Xtest <- Xtest[,indsel]
  
  blocksub <- lapply(1:length(block), function(x) which(blocksub==x))
  
  
  # Train random survival forest variant:
  
  blockforobj <- blockfor(Xtrain, ytrain, num.trees = 2000, replace = FALSE, probability = FALSE, blocks=blocksub,
                          nsets = 300, num.trees.pre = 1500, splitrule="extratrees", 
                          block.method = block.method, num.threads=1)						  
  
  paramvalues <- blockforobj$paramvalues
  
  
  # Obtain predictions on test data:
  
  riskpred <- rowSums(predict(blockforobj$forest, data=Xtest, block.method=block.method, num.threads=1)$chf)
  
  
  # Return results:
  
  return(list(riskpred=riskpred, paramvalues=paramvalues, mtry=NULL))
  
}






# Wrapper function that performs supervised variable selection 
# using univariate Cox regressions on training data, subsequently
# trains a random survival forest on the training data and, lastly,
# obtains predictions on test data using the random survival forest
# trained on the training data.

# Input parameters:

# Xtrain        - covariate matrix of the training data set
# ytrain        - target variable of the training data set. Matrix with two 
#                 columns, where the first column contains the vector of 
#                 survival/censoring times (one for each observation) and the 
#                 second column contains the status variable, that has the 
#                 value '1' if the corresponding time is a survival time and 
#                 '0' if that time is a censoring time.
# Xtest         - covariate matrix of the test data set
# ytest         - target variable in the test data set. Matrix with two
#                 columns, where the first column contains the vector of 
#                 survival/censoring times (one for each observation) and the 
#                 second column contains the status variable, that has the 
#                 value '1' if the corresponding time is a survival time and 
#                 '0' if that time is a censoring time.
# block         - A list of length equal to the number of blocks considered. 
#                 Each entry contains the vector of column indices in 'Xtrain' 
#                 of the covariates in one of the blocks.

randomsurvivalforestwrap <- function(Xtrain, ytrain, Xtest, ytest, block) {
  
  library("survival")
  library("ranger")
  
  # Perform supervised variable selection per block using univariate
  # Cox regressions:
  
  indsel <- c()
  blocksub <- c()
  
  for(i in seq(along=block)) {
    if(length(block[[i]]) <= 2500) {
      indsel <- c(indsel, block[[i]])
      blocksub <- c(blocksub, rep(i, length(block[[i]])))
    }
    else {
      varinds <- block[[i]]
      pvals <- apply(Xtrain[,varinds], 2, function(y2) summary(coxph(Surv(ytrain[,1], ytrain[,2]) ~ y2))$coef[,"Pr(>|z|)"][1])
      indsel <- c(indsel, varinds[order(pvals)[1:2500]])
      blocksub <- c(blocksub, rep(i, 2500))
    }
  }
  
  blocksub <- blocksub[order(indsel)]
  indsel <- indsel[order(indsel)]
  
  Xtrain <- Xtrain[,indsel]
  Xtest <- Xtest[,indsel]
  
  blocksub <- lapply(1:length(block), function(x) which(blocksub==x))
  
  
  # Train random survival forest variant:
  
  datatemp <- data.frame(ytrain[,1], ytrain[,2], Xtrain)
  names(datatemp)[1:2] <- c("time", "status")
  
  mtrygrid <- ceiling(c(0.1, 0.25, 0.5, 1, 2)*sqrt(ncol(Xtrain)))
  errmtry <- 0
  for(i in seq(along=mtrygrid)) {
    errmtry[i] <- ranger(Surv(time, status) ~ ., data=datatemp, mtry=mtrygrid[i], num.trees = 1500, splitrule="extratrees", replace = FALSE, num.threads=1)$prediction.error
  }
  if(sum(errmtry==min(errmtry))>1)
    mtryopt <- mtrygrid[sample(which(errmtry==min(errmtry)), size=1)]
  else
    mtryopt <- mtrygrid[which.min(errmtry)]
  
  rangerobj <- ranger(Surv(time, status) ~ ., data=datatemp, mtry=mtryopt, num.trees = 2000, splitrule="extratrees", replace = FALSE, num.threads=1)
  
  
  # Obtain predictions on test data:
  
  riskpred <- rowSums(predict(rangerobj, data=Xtest, num.threads=1)$chf)
  
  
  # Return results:
  
  return(list(riskpred=riskpred, paramvalues=NULL, mtry=mtryopt))
  
}
