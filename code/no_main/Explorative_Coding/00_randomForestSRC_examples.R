setwd("C:/Users/kuche_000/Desktop/Master_Omicsdaten/")
library(randomForestSRC)
library(survival)
library(pec)
library(microbenchmark)

# Example for regular survival analysis w/ RF-SRC Algorithm! -------------------
# Get the data
data(veteran, package = "randomForestSRC")
head(veteran)

# Randomized trial of two treatment regimens for lung cancer
v.obj <- rfsrc(Surv(time, status) ~ ., data = veteran,
               ntree = 100, block.size = 1)

# Print and plot the grow object
print(v.obj)
plot(v.obj)

# Plot survival curves for first 10 individuals -- direct way
matplot(v.obj$time.interest, 100 * t(v.obj$survival.oob[1:10, ]),
        xlab = "Time", ylab = "Survival", type = "l", lty = 1)

# Plot survival curves for first 10 individuals using function "plot.survival"
plot.survival(v.obj, subset = 1:10)


# Fast nodesize optimization for veteran data optimal nodesize in survival is 
# larger than other families see the function "tune" for more examples
tune.nodesize(Surv(time,status) ~ ., veteran)


# Example of imputation in survival analysis -----------------------------------
data(pbc, package = "randomForestSRC")
head(pbc)

# Fit RF without imputing missing data
pbc.obj2 <- rfsrc(Surv(days, status) ~ ., pbc, nsplit = 10, 
                  na.action = "na.impute")


# Same as above but we iterate the missing data algorithm
pbc.obj3 <- rfsrc(Surv(days, status) ~ ., pbc,
                  na.action = "na.impute", nimpute = 3)


# Fast way to impute the data (no inference is done)
pbc.imp <- impute(Surv(days, status) ~ ., pbc, splitrule = "random")

# see impute for more details!
# Compare RF-SRC to Cox regression ---------------------------------------------
# C-index and Brier score measures of performance
predictSurvProb.rfsrc <- function(object, newdata, times, ...){
  ptemp <- predict(object,newdata=newdata,...)$survival
  pos <- sindex(jump.times = object$time.interest, eval.times = times)
  p <- cbind(1,ptemp)[, pos + 1]
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
    stop("Prediction failed")
  p
}

# Data, formula specifications
data(pbc, package = "randomForestSRC")
head(pbc)
pbc.na <- na.omit(pbc) # remove NA's

# Create formulas we use to fit the model!
surv.f <- as.formula(Surv(days, status) ~ .)
pec.f  <- as.formula(Hist(days,status) ~ 1)

# Run cox/rfsrc models [use smaller number of trees for illustration]
cox.obj   <- coxph(surv.f, data = pbc.na, x = TRUE)
rfsrc.obj <- rfsrc(surv.f, pbc.na, ntree = 150)

# Compute bootstrap cross-validation estimate of expected Brier score
#     ---> Not working probably yet <---
# [see Mogensen, Ishwaran & Gerds (2012) - Journal of Statistical Software]
set.seed(17743)
prederror.pbc <- pec(list(cox.obj,rfsrc.obj), data = pbc.na, formula = pec.f,
                     splitMethod = "bootcv", B = 50)

print(prederror.pbc)
plot(prederror.pbc)

# Compute out-of-bag C-index for cox regression and compare to rfsrc
rfsrc.obj <- rfsrc(surv.f, pbc.na)
cat("out-of-bag Cox Analysis ...", "\n")
cox.err <- sapply(1:100, function(b) {
  if (b %% 10 == 0) cat("cox bootstrap:", b, "\n")
  train <- sample(1:nrow(pbc.na), nrow(pbc.na), replace = TRUE)
  cox.obj <- tryCatch({coxph(surv.f, pbc.na[train, ])}, error = function(ex){NULL})
  if (!is.null(cox.obj)) {
    get.cindex(pbc.na$days[-train], pbc.na$status[-train], predict(cox.obj, pbc.na[-train, ]))
  } else NA
})
cat("\n\tOOB error rates\n\n")
cat("\tRSF : ", rfsrc.obj$err.rate[rfsrc.obj$ntree], "\n")
cat("\tCox regression : ", mean(cox.err, na.rm = TRUE), "\n")

# competing risks --------------------------------------------------------------
# WIHS analysis - cumulative incidence function (CIF) for HAART and AIDS 
#                 stratified by IDU
data(wihs, package = "randomForestSRC")
head(wihs)

# Fit a RF on the data
wihs.obj <- rfsrc(Surv(time, status) ~ ., wihs, nsplit = 3, ntree = 100)

# plot the competing risks
plot.competing.risk(wihs.obj)

# Plot them for different settings
cif <- wihs.obj$cif.oob
Time <- wihs.obj$time.interest
idu <- wihs$idu
cif.haart <- cbind(apply(cif[,,1][idu == 0,], 2, mean),
                   apply(cif[,,1][idu == 1,], 2, mean))
cif.aids <- cbind(apply(cif[,,2][idu == 0,], 2, mean),
                  apply(cif[,,2][idu == 1,], 2, mean))
matplot(Time, cbind(cif.haart, cif.aids), type = "l",
        lty = c(1,2,1,2), col = c(4, 4, 2, 2), lwd = 3,
        ylab = "Cumulative Incidence")
legend("topleft",
       legend = c("HAART (Non-IDU)", "HAART (IDU)", "AIDS (Non-IDU)", "AIDS (IDU)"),
       lty = c(1,2,1,2), col = c(4, 4, 2, 2), lwd = 3, cex = 1.5)

# illustrates the various splitting rules
# illustrates event specific and non-event specific variable selection
if (library("survival", logical.return = TRUE)) {
  ## use the pbc data from the survival package
  ## events are transplant (1) and death (2)
  data(pbc, package = "survival")
  pbc$id <- NULL
  ## modified Gray's weighted log-rank splitting
  pbc.cr <- rfsrc(Surv(time, status) ~ ., pbc)
  ## log-rank event-one specific splitting
  pbc.log1 <- rfsrc(Surv(time, status) ~ ., pbc,
                    splitrule = "logrank", cause = c(1,0), importance = TRUE)
  ## log-rank event-two specific splitting
  pbc.log2 <- rfsrc(Surv(time, status) ~ ., pbc,
                    splitrule = "logrank", cause = c(0,1), importance = TRUE)
  ## extract VIMP from the log-rank forests: event-specific
  ## extract minimal depth from the Gray log-rank forest: non-event specific
  var.perf <- data.frame(md = max.subtree(pbc.cr)$order[, 1],
                         vimp1 = 100 * pbc.log1$importance[ ,1],
                         vimp2 = 100 * pbc.log2$importance[ ,2])
  print(var.perf[order(var.perf$md), ])
}
# Applied to OMICs-Data --------------------------------------------------------
setwd("C:/Users/kuche_000/Desktop/Master_Omicsdaten/")
load("./data/external/Dr_Hornung/Data/ProcessedData/UCEC.Rda")

# Use MIRNA data, with ~866 features
dim(mutation)
dim(targetvar)

# bind them together to one DF
data_combined <- cbind(targetvar, mutation[,1:2500])

dim(data_combined)
colnames(data_combined)[1:10]

# Remove entrances with negative timeentries!
data_combined <- data_combined[data_combined$time >= 0,]

# Randomized trial of two treatment regimens for lung cancer
start_100 <- Sys.time()
rf_surv_100 <- rfsrc(Surv(time, status) ~ ., data = data_combined,
                 ntree = 100, block.size = 1, nsplit = 1000)
end_100 <- Sys.time()

print("Time needed 100trees, 250 feas:")
print(end_100 - start_100)

start_1000 <- Sys.time()
rf_surv_1000 <- rfsrc(Surv(time, status) ~ ., data = data_combined, mtry = 1000,
                 ntree = 1000, block.size = 1, nsplit = 1000)
end_1000 <- Sys.time()

print("Time needed 1000trees, 250 feas:")
print(end_1000 - start_1000)

# Print and plot the grow object
print(rf_surv_1000)
plot(rf_surv_1000)

# Plot survival curves for first 10 individuals -- direct way
matplot(rf_surv_1000$time.interest, 100 * t(rf_surv_1000$survival.oob[1:10, ]),
        xlab = "Time", ylab = "Survival", type = "l", lty = 1)

# Plot survival curves for first 10 individuals using function "plot.survival"
plot.survival(rf_surv_1000, subset = 1:10)


# Fit Survival Trees on missing data -------------------------------------------
# ----- Data
data(veteran, package = "randomForestSRC")

# check basic properties and check for missing values!
dim(veteran); head(veteran)
any(apply(veteran, 1, function(x) any(is.na(x))))

# Add missing Values in block-wise style!
veteran$age[1:68]     <- NA
veteran$karno[69:137] <- NA

# check data, each row should contain a NA
head(veteran); tail(veteran)
all(apply(veteran, 1, function(x) any(is.na(x))))


# ----- Now fit the RFforSV with the 2 implemented Methods!

# [1] throw away Obs. that contain at least 1 missing value
#     --> Error in blockwise settings!
RF_no_imputation <- rfsrc(Surv(time, status) ~ ., data = veteran,
                          ntree = 100, block.size = 1, na.action = "na.omit")

# [2] Impute missing Values
#     --> This works directly! 
#     --> don't know imputation technique!
RF_w_imputation <- rfsrc(Surv(time, status) ~ ., data = veteran,
                         ntree = 100, block.size = 1, na.action = "na.impute",
                         nimpute = 3)

#   Higher Amount of 'nimpute' makes it better
RF_w_imputation2 <- rfsrc(Surv(time, status) ~ ., data = veteran,
                          ntree = 100, block.size = 1, na.action = "na.impute",
                          nimpute = 10)

#   Lower Amount of 'nimpute' makes it worse
RF_w_imputation3 <- rfsrc(Surv(time, status) ~ ., data = veteran,
                          ntree = 100, block.size = 1, na.action = "na.impute")

# In overall, the only reasonable, directly integrated method is the imputation
# eventhough this is not really appliable in our setting!




# Own Usages -------------------------------------------------------------------


# Fit SurvTree based on OMICS Features and check OOB Rate ----------------------
load("./data/external/Dr_Hornung/Data/ProcessedData/BLCA.Rda")

data <- cbind(clin, targetvar)

microbenchmark(rfsrc(Surv(time, status) ~ ., data = data,
                     ntree = 100, block.size = 1))

rfsrc(Surv(time, status) ~ ., data = data,
      ntree = 100, block.size = 1)


data <- cbind(rna, targetvar)

rfsrc(Surv(time, status) ~ ., data = data,
      ntree = 100, block.size = 1)


data <- cbind(rna, clin, targetvar)

rfsrc(Surv(time, status) ~ ., data = data,
      ntree = 100, block.size = 1)


data <- cbind(rna, clin, mutation, targetvar)

rfsrc(Surv(time, status) ~ ., data = data,
      ntree = 1500, block.size = 1)


data <- cbind(rna, clin, targetvar)
data[1:100, 1:100]         <- NA
data[200:300, 10000:10100] <- NA

rfsrc(Surv(time, status) ~ ., data = data,
      ntree = 100, block.size = 1, na.action = 'na.impute', nimpute = 25)


# Roman's Method Adaption w randomForestSRC ------------------------------------
# ----- Data
data(veteran, package = "randomForestSRC")

# check basic properties and check for missing values!
dim(veteran); head(veteran)
any(apply(veteran, 1, function(x) any(is.na(x))))

# Add missing Values in block-wise style!
veteran$age[1:68]     <- NA
veteran$karno[69:137] <- NA

# check data, each row should contain a NA
head(veteran); tail(veteran)
all(apply(veteran, 1, function(x) any(is.na(x))))


# ----- Adaption of the Method
# [1] Slice the data into blocks w/o any missing values!
get_blocks_no_NAs <- function(df){
  "
  Function to slice a df into X sub-DF, where none of the subdf contains any NAs
  "
  rows_w_NAs <- apply(df, 1, function(x) any(is.na(x)))
  print(paste0("Amount of Rows w/ missing data: ", sum(rows_w_NAs), 
               " corresping to: ", sum(rows_w_NAs)/nrow(df) * 100, "% of rows"))
  
  cols_w_NAs <- colnames(df)[apply(df, 2, function(x) any(is.na(x)))]
  print("Cols with missing Values:"); print(cols_w_NAs)
  
  # Now do the following for all columns with NAs:
  # Remove the column with at least one NA from the original DF [copying] 
  # rows that have a NA in the rows in one of the cols of 'cols_w_NAs'
  # then check whether 
  df_curr       <- data.frame(df[,!(colnames(df) %in% cols_w_NAs[1])])
  df_curr_no_na <- na.omit(df_curr)
  
  if (nrow(df_curr_no_na) < 3) {
    break
  } else {
    print(paste("removing Variable", cols_w_NAs[1], 
                "led to a DF without any missing data of the dimensions", 
                nrow(df), "x", ncol(df)))
    tree <- rfsrc(Surv(time, status) ~ ., data = df_curr_no_na, ntree = 1,
                  statistics = TRUE, forest = F)
    
  }
  
  
  df[,cols_w_NAs]
  
}


# Fit a RF with 'ntree = 1' multiple times on the subsets!
rfsrc(Surv(time, status) ~ ., data = veteran, ntree = 1)

# Roman's Mehtod Adaption w/ rangers -------------------------------------------
library(ranger)
library(survival)
rf <- ranger(Surv(time, status) ~ ., data = veteran, num.trees = 1, mtry = 3)
rf$num.trees
rf$forest$split.varIDs


rf$forest$split.varIDs


# Roman's Mehtod Adaption w/ LTRCtrees -----------------------------------------
## Adjust data & clean data
library(LTRCtrees)
library(survival)
set.seed(0)
## Since LTRCART uses cross-validation to prune the tree, specifying the seed 
## guarantees that the results given here will be duplicated in other analyses
Data <- flchain
Data <- Data[!is.na(Data$creatinine),]
Data$End <- Data$age + Data$futime/365
DATA <- Data[Data$End > Data$age,]
names(DATA)[6] <- "FLC"

## Setup training set and test set
Train = DATA[1:500,]
Test = DATA[1000:1020,]


LTRCART.obj <- LTRCART(Surv(age, End, death) ~ sex + FLC + creatinine, Train)
LTRCIT.obj <- LTRCIT(Surv(age, End, death) ~ sex + FLC + creatinine, Train)

predict(LTRCART.obj, newdata=Test, type = "response")
