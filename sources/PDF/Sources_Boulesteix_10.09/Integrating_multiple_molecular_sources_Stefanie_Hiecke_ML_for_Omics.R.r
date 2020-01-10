###########################################################################################################################################################################################################
#                                                                                     
#   Project     :   R code for manuscript "Intergrating multiple molecular sources into a clinical risk prediction signature by extracting complementary information"                                                                    
#   Author      :   Stefanie Hieke                                                                    
#   Date        :   26.05.2016
#																					                                                                                                                                                                                                                                               
###########################################################################################################################################################################################################

###########################################################################################################################################################################################################
##   Required librarys        ###
library(CoxBoost)    # Componentwise likelihood-based boosting adapted for Cox proportional hazards models
library(GAMBoost)    # Componentwise likelihood-based boosting adapted for continuous endpoints
library(peperr)      # Generation of training and test set indices for resampling procedures (e.g. bootstrap with and without replacement)
library(penalized)   # L1 (lasso) penalized estimation in Cox proportional hazards models

###########################################################################################################################################################################################################

###########################################################################################################################################################################################################
## 1. SNP risk prediction model (1) estimated by componentwise likelihood-based boosting

### load(additional_file_1.Rdata)

### current.status: event indicator (SNP data)
### current.time  : event time (SNP data)
### xmat          : covariate (SNP) matrix in addition to the five clinical covariates
### SNP_con_GEP   : Overlap sample ID (SNP and GEP data)

penalty <- 19*sum(current.status==1)
boost.stepno <- 500 
K <- 10

## Main tuning parameter (boosting steps) determined by cross-validation

set.seed(22)
cv.res.full.SNP <- cv.CoxBoost(current.time,current.status,xmat,
                               unpen.index=1:5,
                               standardize=FALSE,
                               maxstepno=boost.stepno,
                               penalty=(1-1/K)*penalty,
                               x.is.01=TRUE,
                               K=K,
                               multicore=10)

res.full.SNP <- CoxBoost(current.time,current.status,xmat,
                         unpen.index=1:5,
                         standardize=FALSE,
                         stepno=cv.res.full.SNP$optimal.step,
                         penalty=penalty,
                         x.is.01=TRUE)

colnames(res.full.SNP$coefficients) <- res.full.SNP$xnames

final.coef.SNP <- res.full.SNP$coefficients[nrow(res.full.SNP$coefficient),]
betas <- final.coef.SNP[-c(1:5)]

lp.SNP <- xmat[,-c(1:5)] %*% betas                                             # SNP signature
lp.SNP.con <- lp.SNP[SNP_con_GEP[,"SNP"]]                            # SNP signature conditioning on the overlap samples (SNP ID)
names(lp.SNP.con) <- SNP_con_GEP[,"GEP"]                             # SNP signature conditioning on the overlap samples (GEP ID)

###########################################################################################################################################################################################################

###########################################################################################################################################################################################################

## 2. Linking model (3) estimated by componentwise likelihood-based boosting

### load(additional_file_2.Rdata)

### totalmat : covariate (GEP) matrix in addition to the five clinical covariates

lp.response <- lp.SNP.con                                                      # SNP signature conditioning on the overlap samples
totalmat_con<- totalmat[names(lp.response),]                                   # totalmat conditioning on the biological samples with measurements from both sources in parallel (overlap samples)
totalmat_pred <- totalmat[-which(rownames(totalmat) %in% names(lp.response)),] # totalmat conditioning on biological samples with measurements from GEP only (non-overlap samples)

lamb <- (dim(totalmat_con)[1])*((1/0.005)-1)
expr.index <- dim(totalmat_con)[2]-5

## Main tuning parameter (boosting steps) determined by cross-validation

set.seed(5544)
res.GEP <- cv.GAMBoost(y=lp.response,
                       x.linear=totalmat_con,
                       K=10,
                       family=gaussian(),
                       standardize.linear=TRUE,
                       penalty.linear=c(rep(0,length=5),
                           rep(lamb,length=expr.index)),
                       criterion="score",
                       multicore=10)

pred.lp <- predict(object=res.GEP, newdata.linear=totalmat_pred,
                   type="response")

lp <- c(lp.response,pred.lp)                                                   # SNP signature including known SNP signatur (overlap samples) and predicted SNP signature (non-overlap samples)

###########################################################################################################################################################################################################
## 2.1. Prediction performance: Linking model (3)

###########################################################################################################################################################################################################
## 2.1.1. bootstrap .632+ prediction error estimate

### res.GEP: linking model (3) fitted on full overlap samples

K <- 10
boost.stepno <- 500
boot.samples <- 1:100
do.full <- TRUE
do.bootstrap <- TRUE

if (do.full) {

    null.apparent <- sum((mean(lp.response)-lp.response)^2)/length(lp.response)               # benchmark null model (determined by the apparent error)

    resp.mat <- matrix(rep(as.vector(lp.response), 26), nrow=26, byrow=T)
    
    pred <- predict(res.GEP, newdata.linear=totalmat_con,
                    type="response")                                            
 
    pred.full <- sum((pred-lp.response)^2)/length(lp.response)                                # apparent error

    pred.noinf <- mean((resp.mat-pred)^2)                                                     # no-informative error
}

if (do.bootstrap) {
    boot.n <- 100
    set.seed(2002121)
    boot.indices <- resample.indices(n = nrow(totalmat_con),
                                     method = "sub632",sample.n=boot.n)                       # without replacement
    
    for (actual.boot in boot.samples){
    
        avail <- (1:nrow(totalmat_con) %in% boot.indices$sample.index[[actual.boot]])
        test.index <- (1:nrow(totalmat_con) %in% boot.indices$not.in.sample[[actual.boot]])
        
        penalty <- nrow(totalmat_con[avail,])*((1/0.005)-1)
        expr.index <- dim(totalmat_con)[2]-5
        
        res.boot <- cv.GAMBoost(y=lp.response[avail],
                                x.linear=totalmat_con[avail,],
                                K=K,family=gaussian(),
                                standardize.linear=TRUE,
                                penalty.linear=c(rep(0,length=5),
                                    rep(penalty,length=expr.index)),
                                criterion="score",
                                multicore=10)
        
        pred.test <- predict(res.boot,newdata.linear=totalmat_con[test.index,],
                                   type="response")                                        
        
        pred.boot <- sum((pred.boost.test-lp.response[test.index])^2)/
            length(lp.response[test.index])                                                    # MSE in the bth bootstrap
        
        save(pred.boot,file=paste("pred_result",actual.boot,".RData",sep=""))
    }
}

boot.boost <- NULL
boot.null <- NULL

for (i in c(1:100)) {
    load(paste("pred_result",i,".RData",sep=""))
    boot.boost <- c(boot.boost,pred.boot)
}

mean.boot <- mean(boot.boost)                                                                  # bootstrap error

relative.overfit <- ifelse(pred.noinf > pred.full,
                           (ifelse(mean.boot < pred.noinf,mean.boot,pred.noinf) - pred.full)/
                           (pred.noinf - pred.full),0)                                         # relative overfitting rate

weights <- .632/(1-.368*relative.overfit) 

boot632 <- (1-weights)*pred.full + weights*ifelse(mean.boot < pred.noinf,mean.boot.,prd.noinf) # .632+ estimate

###########################################################################################################################################################################################################
## 2.1.2. bootstrap .632+ prediction error estimate individually for each bootstrap sample

modboot <- unlist(lapply(seq_along(1:100), function(i) {
    temp <- boot.boost[i]
    temp1 <- ifelse(pred.noinf > pred.full,
                    (ifelse(temp < pred.noinf,temp,pred.noinf) - pred.full)/
                    (pred.noinf - pred.full.boost),0)
    weights <- .632/(1-.368*temp1)
    temp2 <- (1-weights)*pred.full + weights*ifelse(temp < pred.noinf,temp,pred.noinf)
    temp2
}))


###########################################################################################################################################################################################################

###########################################################################################################################################################################################################
## 3. Sequential complementary strategy fitted on the original data:
##    risk prediction model (2) for the GEP data estimated by componentwise likelihood-based boosting incorporating the SNP signature as a fixed offset

### current.status: event indicator (GEP data)
### current.time  : event time (GEP data)
### totalmat      : covariate (GEP) matrix in addition to the five clinical covariates
### clon.gene     : genes linked to cDNA IMAGE clones

lamb <- (1/0.02-1)*sum(current.status)
boost.stepno <- 500 
K <- 10

lp <- lp[rownames(totalmat)]
totalmatlp <- cbind(lp,totalmat)

## Main tuning parameter (boosting steps) determined by cross-validation

set.seed(297113)
cv.seq.GEP <- cv.CoxBoost(current.time,current.status,
                          totalmatlp,unpen.index=1:6,
                          standardize=TRUE,
                          maxstepno=boost.stepno,penalty=lamb,
                          criterion="pscore",x.is.01=FALSE,
                          K=K,trace=TRUE,multicore=10)

seq.GEP <- CoxBoost(current.time,current.status,
                    totalmatlp,unpen.index=1:6,
                    standardize=TRUE,
                    stepno=cv.seq.GEP$optimal.step,
                    penalty=lamb,criterion="pscore",
                    x.is.01=FALSE)

colnames(seq.GEP$coefficients) <- seq.GEP$xnames
final.coef.seq.GEP <- seq.GEP$coefficients[nrow(seq.GEP$coefficient),] # coefficients in the optimal boosting step
chos.index <- which(final.coef.seq.GEP !=0)
beta.chos <- final.coef.cDNA.cv[chos.index]
seq.GEP.beta <- beta.chos[-c(1:6)]                                     # GEP features which have a coefficient !=0
gene.found <- clon.gene[names(seq.GEP.beta),]
seq.frame <- data.frame(ID=names(seq.GEP.beta) ,coef=round(seq.GEP.beta,5), gene=as.factor(gene.found$Name))

###########################################################################################################################################################################################################

###########################################################################################################################################################################################################
## 4. Reference approach fitted on the original data:
##    risk prediction model (2) estimated by componentwise likelihood-based boosting without adjusting for the effects of the SNP data:

### current.status: event indicator (GEP data)
### current.time  : event time (GEP data)
### totalmat      : covariate (GEP) matrix in addition to the five clinical covariates
### clon.gene     : genes linked to cDNA IMAGE clones

lamb <- (1/0.02-1)*sum(current.status)
boost.stepno <- 500 
K <- 10

## Main tuning parameter (boosting steps) determined by cross-validation

set.seed(11)
cv.ref.GEP <- cv.CoxBoost(current.time,current.status,
                          totalmat,unpen.index=1:5,
                          standardize=TRUE,
                          maxstepno=boost.stepno,penalty=lamb,
                          criterion="pscore",x.is.01=FALSE,
                          K=K,trace=FALSE,multicore=10)

ref.GEP <- CoxBoost(current.time,current.status,
                    totalmat,unpen.index=1:5,
                    standardize=TRUE,
                    stepno=cv.ref.GEP$optimal.step,
                    penalty=lamb,criterion="pscore",
                    x.is.01=FALSE,trace=FALSE)

colnames(ref.GEP$coefficients) <- ref.GEP$xnames
final.coef.ref.GEP <- ref.GEP$coefficients[nrow(ref.GEP$coefficient),] # coefficients in the optimal boosting step
chos.index <- which(final.coef.ref.GEP !=0)
beta.chos <- final.coef.cDNA.cv[chos.index]
ref.GEP.beta <- beta.chos[-c(1:5)]                                     # GEP features which have a coefficient !=0
gene.found <- clon.gene[names(ref.GEP.beta),]
ref.frame <- data.frame(ID=names(ref.GEP.beta) ,coef=round(ref.GEP.beta,5), gene=as.factor(gene.found$Name))

###########################################################################################################################################################################################################

###########################################################################################################################################################################################################
## 5. Model stability

### current.status: event indicator (GEP data)
### current.time  : event time (GEP data)
### clon.gene     : genes linked to cDNA IMAGE clones

set.seed(227111)
resample.index <- resample.indices(n=dim(totalmat)[1],
                                   sample.n=100, method="sub632") # without replacement
sample.indices <- resample.index$sample.index

###########################################################################################################################################################################################################
## 5.1. Reference approach and risk prediction model (2) within the sequential complementarys strategy are estimated by boosting

boost.stepno <- 500 
K <- 10

###########################################################################################################################################################################################################
## 5.1.1. Sequential complementary strategy

### totalmatlp    : covariate (GEP) matrix in addition to the five clinical covariates and linera predictor including SNP information
### seq.frame     : gene features selected by the sequential complementary strategy

seq.inclu_freq <- lapply(sample.indices, function(x) {
  temp.data <- totalmatlp[x,]
  temp.time <- current.time[x]
  temp.status <- current.status[x]
  lamb <- (1/0.02-1)*sum(current.status[x]==1)

  cv.Boost <- cv.CoxBoost(temp.time,temp.status,temp.data,
                          unpen.index=1:6,
                          standardize=TRUE,
                          maxstepno=boost.stepno, penalty=lamb,
                          criterion="pscore", x.is.01=FALSE,
                          K=K, multicore=10)

  if(length(cv.Boost$optimal.step)==0 | cv.Boost$optimal.step==0){ # stop boosting if optimal boosting number equals to zero
      gene.frame <- 0
      optno <- 0
  }else{
      optno <- cv.Boost$optimal.step
      
      fit <- CoxBoost(temp.time,temp.status,temp.data,
                      unpen.index=c(1:5,dim(temp.data)[2]), standardize=TRUE,
                      stepno=cv.Boost$optimal.step, penalty=lamb,
                      criterion="pscore", x.is.01=FALSE)
      
      colnames(fit$coefficients) <- colnames(temp.data)
      fit.betas <- fit$coefficients[,-c(1:6), drop=FALSE]
      
      if(sum(fit.betas[nrow(fit.betas),]==0)==dim(fit.betas)[[2]]){
          gene.frame <- NULL
      }else{
          chosen.index <- which(fit.betas[nrow(fit.betas),]!=0)
          chosen.betas <- fit.betas[nrow(fit.betas),chosen.index]
          genes <- clon.gene[names(chosen.betas),]
          gene.frame <- data.frame(ID=names(chosen.betas), coef=round(chosen.betas,5),
                                   genes=as.factor(genes$Name))
          rownames(gene.frame) <- NULL
      }
  }
  res <- list(stepnrs=optno, genes=gene.frame)
  res
})

sampled_genes_seq<- table(unlist(lapply(seq.inclu_freq, function(x) {
    genes <- unique(x$genes$genes)
    genes
})))

###########################################################################################################################################################################################################
## 5.1.2. Reference approach

### totalmat      : covariate (GEP) matrix in addition to the five clinical covariates
### ref.frame     : gene features selected by the reference approach

ref.inclu_freq <- lapply(sample.indices, function(x) {
    temp.data <- totalmat[x,]
    temp.time <- current.time[x]
    temp.status <- current.status[x]
    
    cv.Boost <- cv.CoxBoost(temp.time,temp.status,temp.data,
                            unpen.index=c(1:5), standardize=TRUE,
                            maxstepno=boost.stepno, penalty=lamb,
                            criterion="pscore", x.is.01=FALSE,
                            K=K, multicore=10)
    
    if(length(cv.Boost$optimal.step)==0 | cv.Boost$optimal.step==0){
        gene.frame <- 0
        optno <- 0
    }else{
        optno <- cv.Boost$optimal.step
        
        fit <- CoxBoost(temp.time,temp.status,temp.data,
                        unpen.index=c(1:5), standardize=TRUE,
                        stepno=cv.Boost$optimal.step, penalty=lamb,
                        criterion="pscore", x.is.01=FALSE)
        
        colnames(fit$coefficients) <- colnames(temp.data)
        fit.betas <- fit$coefficients[,-c(1:5) , drop=FALSE]
        
        if(sum(fit.betas[nrow(fit.betas),]==0)==dim(fit.betas)[[2]]){
            gene.frame <- NULL
        }else{
            chosen.index <- which(fit.betas[nrow(fit.betas),]!=0)
            chosen.betas <- fit.betas[nrow(fit.betas),chosen.index]
            genes <- clon.gene[names(chosen.betas),]
            gene.frame <- data.frame(ID=names(chosen.betas), coef=round(chosen.betas,5),
                                     genes=as.factor(genes$Name))
            rownames(gene.frame) <- NULL            
        }
    }
    res <- list(stepnrs=optno, genes=gene.frame)
    res
})

sampled_genes_ref <- table(unlist(lapply(ref.inclu_freq, function(x) {
    genes <- unique(x$genes$genes)
    genes
})))

###########################################################################################################################################################################################################
## 5.1.3. Compare inclusion frequencies between sequential complementary strategy and reference approach

if_ref_genes <- data.frame(genes=names(sampled_genes_ref[as.character(ref.frame$gene)]),
                           count.ref=sampled_genes_ref[as.character(ref.frame$gene)],
                           count.seq=sampled_genes_seq[as.character(ref.frame$gene)])

if_seq_genes <- data.frame(genes=names(sampled_genes_seq[as.character(seq.frame$gene)]),
                           count.ref=sampled_genes_ref[as.character(seq.frame$gene)],
                           count.seq=sampled_genes_seq[as.character(seq.frame$gene)])

###########################################################################################################################################################################################################
## 5.2. Reference approach and risk prediction model (2) within the sequential complementarys strategy are estimated by lasso

###########################################################################################################################################################################################################
## 5.2.1. Sequential complementary strategy

### totalmatlp    : covariate (GEP) matrix in addition to the five clinical covariates and linera predictor including SNP information
### seq.frame     : gene features selected by the sequential complementary strategy

lassomatlp <- as.data.frame(totalmatlp)

## Sequential complementary strategy (risk prediction model (2)) fitted on the original data and estimated by lasso
set.seed(27021430)
proflambda <- profL1(response=Surv(current.time,current.status),
                     penalized=lassomatlp[,7:dim(lassomatlp)[2]],
                     unpenalized=~lp+age+FLT3+NPM1+WBC+low.risk,
                     data=lassomatlp,
                     model="cox",
                     standardize=TRUE,
                     fold=10,
                     lambda2=0,
                     trace=FALSE, plot=TRUE)

lambda <- proflambda$lambda[which.max(proflambda$cvl)]

seq.fit <- penalized(response=Surv(current.time,current.status),
                     penalized=lassomatlp[,7:dim(lassomatlp)[2]],
                     unpenalized=~lp+age+FLT3+NPM1+WBC+low.risk,
                     data=lassomatlp,
                     model="cox",
                     standardize=TRUE,
                     lambda1=lambda,
                     lambda2=0,
                     steps=100)

res.seq <- coefficients(seq.fit[[100]])[-c(1:6)]
gene.found.seq <- clon.gene[names(res.seq),]

### Inclusion frequencies
set.seed(272201440)
seq.inclu_freq.lasso <- mclapply(sample.indices, function(x) {
    temp.data <- lassomatlp[x,]
    temp.time <- current.time[x]
    temp.status <- current.status[x]
    
    proflambda <- profL1(response=Surv(temp.time,temp.status),
                         penalized=temp.data[,7:dim(temp.data)[2]],
                         unpenalized=~lp+age+FLT3+NPM1+WBC+low.risk,
                         data=temp.data,
                         model="cox",
                         standardize=TRUE,
                         fold=10,
                         lambda2=0,
                         trace=FALSE)
    
    optlam <- proflambda$lambda[which.max(proflambda$cvl)]
    
    fit <- penalized(response=Surv(temp.time,temp.status),
                     penalized=temp.data[,7:dim(temp.data)[2]],
                     unpenalized=~lp+age+FLT3+NPM1+WBC+low.risk,
                     data=temp.data,
                     model="cox",
                     standardize=TRUE,
                     lambda1=optlam,
                     lambda2=0)
    
    fit.betas <- coefficients(fit)
    
    if(length(fit.betas)==5){
        gene.frame <- NULL
    }else{
        chosen.betas <- fit.betas[-c(1:6)]
        genes <- clon.gene[names(chosen.betas),]
        gene.frame <- data.frame(ID=names(chosen.betas),
                                 genes=as.factor(genes$Name))
        rownames(gene.frame) <- NULL
    }
    res <- list(optlam=optlam, genes=gene.frame)
    res
}, mc.cores=20, mc.preschedule=FALSE)

sampled_genes_seq.lasso <- table(unlist(lapply(seq.inclu_freq.lasso, function(x) {
  genes <- unique(x$genes$genes)
  genes
})))

###########################################################################################################################################################################################################
## 5.2.2. Reference approach

### totalmat      : covariate (GEP) matrix in addition to the five clinical covariates
### ref.frame     : gene features selected by the reference approach

lassomat <- as.data.frame(totalmat)

## Reference approach (risk prediction model (2)) fitted on the original data and estimated by lasso
set.seed(27022014)
proflambda <- profL1(response=Surv(current.time,current.status),
                     penalized=lassomat[,6:dim(lassomat)[2]],
                     unpenalized=~age+FLT3+NPM1+WBC+low.risk,
                     data=lassomat,
                     model="cox",
                     standardize=TRUE,
                     fold=10,
                     lambda2=0,
                     trace=FALSE, plot=TRUE)

lambda <- proflambda$lambda[which.max(proflambda$cvl)]

ref.fit <- penalized(response=Surv(current.time,current.status),
                     penalized=lassomat[,6:dim(lassomat)[2]],
                     unpenalized=~age+FLT3+NPM1+WBC+low.risk,
                     data=lassomat,
                     model="cox",
                     standardize=TRUE,
                     lambda1=lambda,
                     lambda2=0,
                     steps=100)

res.ref <- coefficients(ref.fit[[100]])[-c(1:5)]
gene.found.ref <- clon.gene[names(res.ref),]

### Inclusion frequencies
set.seed(272201410)
ref.inclu_freq.lasso <- mclapply(sample.indices, function(x) {
    temp.data <- lassomat[x,]
    temp.time <- current.time[x]
    temp.status <- current.status[x]
    
    proflambda <- profL1(response=Surv(temp.time,temp.status),
                         penalized=temp.data[,6:dim(temp.data)[2]],
                         unpenalized=~age+FLT3+NPM1+WBC+low.risk,
                         data=temp.data,
                         model="cox",
                         standardize=TRUE,
                         fold=10,
                         lambda2=0,
                         trace=FALSE)
    
    optlam <- proflambda$lambda[which.max(proflambda$cvl)]
    
    fit <- penalized(response=Surv(temp.time,temp.status),
                     penalized=temp.data[,6:dim(temp.data)[2]],
                     unpenalized=~age+FLT3+NPM1+WBC+low.risk,
                     data=temp.data,
                     model="cox",
                     standardize=TRUE,
                     lambda1=optlam,
                     lambda2=0)
    
    fit.betas <- coefficients(fit)
    
    if(length(fit.betas)==5){
        gene.frame <- NULL
    }else{
        chosen.betas <- fit.betas[-c(1:5)]
        genes <- clon.gene[names(chosen.betas),]
        gene.frame <- data.frame(ID=names(chosen.betas),
                                 genes=as.factor(genes$Name))
        rownames(gene.frame) <- NULL
    }
    res <- list(optlam=optlam, genes=gene.frame)
    res
}, mc.cores=20, mc.preschedule=FALSE)

sampled_genes_ref.lasso <- table(unlist(lapply(ref.inclu_freq.lasso, function(x) {
  genes <- unique(x$genes$genes)
  genes
})))

###########################################################################################################################################################################################################
## 5.2.3. Compare inclusion frequencies between sequential complementary strategy and reference approach

if_ref_genes.lasso <- data.frame(genes=names(sampled_genes_ref.lasso[as.character(gene.found.ref$Name)]),
                                 count.ref=sampled_genes_ref.lasso[as.character(gene.found.ref$Name)],
                                 count.seq=sampled_genes_seq.lasso[as.character(gene.found.ref$Name)])

if_seq_genes.lasso <- data.frame(genes=names(sampled_genes_seq.lasso[as.character(gene.found.seq$Name)]),
                                 count.ref=sampled_genes_ref.lasso[as.character(gene.found.seq$Name)],
                                 count.seq=sampled_genes_seq.lasso[as.character(gene.found.seq$Name)])

###########################################################################################################################################################################################################

###########################################################################################################################################################################################################
## 6. Prediction performance:
##    Sequential complementary strategy and reference approach (risk prediction model (2)) estimated by boosting

### current.status: event indicator (GEP data)
### current.time  : event time (GEP data)
### totalmat      : covariate (GEP) matrix in addition to the five clinical covariates
### cv.ref.GEP    : Main tuning parameter for the reference approach (risk prediction model (2)) fitted on the original data
### ref.GEP       : reference approach (risk prediction model (2)) fitted on the original data
### totalmatlp    : covariate (GEP) matrix in addition to the five clinical covariates and linera predictor including SNP information
### cv.seq.GEP    : Main tuning parameter for the sequential complementary strategy (risk prediction model (2)) fitted on the original data
### seq.GEP       : sequential complementary strategy (risk prediction model (2)) fitted on the original data
### lp.response   : SNP signature conditioning on the overlap samples

K <- 10
boost.stepno <- 500
eval.times <- seq(from=0,to=1000,length=100)

do.full <- TRUE
do.bootstrap <- TRUE
boot.samples <- 1:100

###########################################################################################################################################################################################################
## 6.1. bootstrap .632+ prediction error estimates

if (do.full) {
    km <- survfit(Surv(time,status) ~ 1, data=data.frame(time=current.time,status=current.status)) 
    
    pec.km.full <- pmpec(km,response=Surv(current.time,current.status),          # Kaplan-Meier benchmark (determined by the apparent error)
                         x=totalmatlp[,1:6],times=eval.times)

    clin.frame <- data.frame(time=current.time,status=current.status,
                             age=totalmatlp[,"age"],
                             low.risk=totalmatlp[,"low.risk"],
                             FLT3=totalmatlp[,"FLT3"],
                             NPM1=totalmatlp[,"NPM1"],
                             WBC=totalmatlp[,"WBC"])
  
    cp1.full <- coxph(Surv(time,status) ~ age + low.risk + FLT3 + NPM1 + WBC,
                      data=clin.frame)
 
    pec.clin.full <- pmpec(cp1.full,response=Surv(current.time,current.status),  # apparent error Cox model
                           x=clin.frame,times=eval.times)
  
    pec.clin.noinf <- pmpec(cp1.full,response=Surv(current.time,current.status), # no-informative error Cox model
                            x=clin.frame,times=eval.times,type="NoInf")
    
    crr.response <- cbind(current.time,current.status)
    colnames(crr.response) <- c("time","status")

    pec.boost.full <- pmpec(seq.GEP,response=crr.response,                       # apparent error reference approach
                            x=totalmat,
                            times=eval.times,
                            model.args=list(complexity=cv.seq.GEP$stepno))

    pec.boost.noinf <- pmpec(seq.GEP,response=crr.response,                      # no-informtive error reference approacj
                             x=totalmat, 
                             times=eval.times,
                             model.args=list(complexity=cv.seq.GEP$stepno),      
                             type="NoInf")
    
    pec.boost.full.lp <- pmpec(seq.GEP,response=crr.response,                    # apparent error sequential complementary strategy
                               x=totalmatlp,
                               times=eval.times,
                               model.args=list(complexity=cv.seq.GEP$stepno))

    pec.boost.noinf.lp <- pmpec(seq.GEP,response=crr.response,                   # no-informtive error sequential complementary strategy
                                x=totalmatlp,
                                times=eval.times,
                                model.args=list(complexity=cv.seq.GEP$stepno),   
                                type="NoInf")
}

if (do.bootstrap) {
    boot.n <- 100
    set.seed(42)
    boot.indices <- resample.indices(n = nrow(totalmat),
                                     method = "sub632",sample.n=boot.n)          # without replacement
    lp_overlap <- lp.response
    crr.response <- cbind(current.time,current.status)
    colnames(crr.response) <- c("time","status")
    for (actual.boot in boot.samples){
        
        avail <- (1:nrow(totalmat) %in% boot.indices$sample.index[[actual.boot]])
        test.index <- (1:nrow(totalmat) %in% boot.indices$not.in.sample[[actual.boot]])
        data_avail <- totalmat[avail,]
        data_test <- totalmat[test.index,]

        ## Cox proportional hazards model
        clin.frame <- data.frame(time=current.time[avail],status=current.status[avail],
                                 age=data_avail[,"age"],low.risk=data_avail[,"low.risk"],
                                 FLT3=data_avail[,"FLT3"],NPM1=data_avail[,"NPM1"],
                                 WBC=data_avail[,"WBC"])

        clin.frame.test <- data.frame(time=current.time[test.index],status=current.status[test.index],
                                      age=data_test[,"age"],low.risk=data_test[,"low.risk"],
                                      FLT3=data_test[,"FLT3"],NPM1=data_test[,"NPM1"],
                                      WBC=data_test[,"WBC"])
        
        cp1.boot <- coxph(Surv(time,status) ~ age + low.risk + FLT3 + NPM1 + WBC,
                          data=clin.frame)

        pec.clin.boot <- pmpec(cp1.boot,response=Surv(current.time[test.index], #  MSE Cox proportional hazards model in the bth bootstrap
                                            current.status[test.index]),
                               x=clin.frame.test,times=eval.times)

        ## Reference approach
        penalty <- 19*sum(current.status[avail]==1)

        set.seed(22)
        cv.res.boot <- cv.CoxBoost(current.time[avail],current.status[avail],
                                   data_avail,unpen.index=1:5,
                                   standardize=TRUE,maxstepno=boost.stepno,
                                   penalty=(1-1/K)*penalty,x.is.01=FALSE,
                                   K=K,multicore=10)

        res.boot <- CoxBoost(current.time[avail],current.status[avail],
                             data_avail,unpen.index=1:5,
                             standardize=TRUE,stepno=cv.res.boot$optimal.step,
                             penalty=penalty,x.is.01=FALSE,
                             return.score=FALSE)

        pec.boost.boot <- pmpec(res.boot,response=crr.response[test.index,],    #  MSE reference approach in the bth bootstrap
                                x=data_test,times=eval.times,
                                model.args=list(complexity=cv.res.boot$optimal.step))
        
        ## Sequential complementary strategy
        overlap_avail.lp <- (names(lp_overlap) %in% rownames(data_avail))
        lp_overlap_avail <- lp_overlap[overlap_avail.lp]
        data.avail.overlap <- data_avail[names(lp_overlap_avail),]
        data.avail.non.overlap <- data_avail[-which(rownames(data_avail) %in% rownames(data.avail.overlap)),]
        
        overlap_test.lp <- (names(lp_overlap) %in% rownames(data_test))
        lp_overlap_test <- lp_overlap[overlap_test.lp]
        data.test.overlap <- data_test[names(lp_overlap_test),]
        data.test.non.overlap <- data_test[-which(rownames(data_test) %in% rownames(data.test.overlap)),]
        
        data.to.pred <- rbind(data.avail.non.overlap,data.test.non.overlap)
        
        lamb <- (nrow(data.avail.overlap))*((1/0.005)-1)
        expr.index <- ncol(data.avail.overlap)-5

        ## Linking model
        set.seed(666)
        fit_lp <- cv.GAMBoost(y=lp_overlap_avail,
                              x.linear=data.avail.overlap,
                              K=K,family=gaussian(),
                              standardize.linear=TRUE,
                              penalty.linear=c(rep(0,length=5),
                                  rep(lamb,length=expr.index)),
                              criterion="score",multicore=10)
                
        pred.lp.non.overlap <- predict(object=fit_lp, newdata.linear=data.to.pred,type="response")
        
        lp.all <- c(lp_overlap_avail,lp_overlap_test,pred.lp.non.overlap)
        lp.full <- lp.all[rownames(totalmat)]
        gamboost.cv.fit <- list(res.subsamples=fit_lp)
        
        lp.train <- lp.full[avail]
        lp.test <- lp.full[test.index]

        
        data_avail.lp <- cbind(lp.train,data_avail)
        colnames(data_avail.lp)[1] <- "lp"

        ## Risk prediction model (2)
        penalty <- 19*sum(current.status[avail]==1)
        
        set.seed(33)
        cv.res.boot.lp <- cv.CoxBoost(current.time[avail],current.status[avail],
                                      data_avail.lp,unpen.index=1:6,
                                      standardize=TRUE,maxstepno=boost.stepno,
                                      penalty=(1-1/K)*penalty,x.is.01=FALSE,
                                      K=K,multicore=10)
        
        res.boot.lp <- CoxBoost(current.time[avail],current.status[avail],
                                data_avail.lp,unpen.index=1:6,
                                standardize=TRUE,stepno=cv.res.boot.lp$optimal.step,
                                penalty=penalty,x.is.01=FALSE,return.score=FALSE)
        
        data_test.lp <- cbind(lp.test,data_test)
        colnames(data_test.lp)[1] <- "lp"
        
        pec.boost.boot.lp <- pmpec(res.boot.lp,response=crr.response[test.index,],   # MSE sequential complementary strategy in the bth bootstrap
                                   x=data_test.lp,times=eval.times,
                                   model.args=list(complexity=cv.res.boot.lp$optimal.step))
        
        boot.res <- list(res.subsample=res.boot,res.subsample.lp=res.boot.lp,GAMB.fit=gamboost.cv.fit)

        save(boot.res,pec.clin.boot,pec.boost.boot,pec.boost.boot.lp,file=paste("pred_resultII",actual.boot,".RData",sep=""))
    }
}

boot.boost <- NULL
boot.boost.lp <- NULL
boot.clin <- NULL

for (i in c(1:100)) {
   load(paste("pred_resultII",i,".RData",sep=""))
   boot.clin <- rbind(boot.clin,pec.clin.boot)
   boot.boost <- rbind(boot.boost,pec.boost.boot)
   boot.boost.lp <- rbind(boot.boost.lp,pec.boost.boot.lp)
}

mean.boot.clin <- colMeans(boot.clin)                                                # bootstrap error Cox proportional hazards model
mean.boot.boost <- colMeans(boot.boost)                                              # bootstrap error reference approach
mean.boot.boost.lp <- colMeans(boot.boost.lp)                                        # bootstrap error sequential complementary strategy

relative.overfit.clin <- ifelse(pec.noinf.clin > pec.full.clin,
                                (ifelse(mean.boot.clin < pec.noinf.clin,mean.boot.clin,pec.noinf.clin) - pec.full.clin)/
                                (pec.noinf.clin - pec.full.clin),0)
weights <- .632/(1-.368*relative.overfit.clin)
boot632p.clin <- (1-weights)*pec.full.clin + weights*ifelse(mean.boot.clin < pec.noinf.clin,mean.boot.clin,pec.noinf.clin)                         # .632+ estimate Cox proportional hazards model

relative.overfit.boost <- ifelse(pec.noinf.boost > pec.full.boost,
                                 (ifelse(mean.boot.boost < pec.noinf.boost,mean.boot.boost,pec.noinf.boost) - pec.full.boost)/
                                 (pec.noinf.boost - pec.full.boost),0)
weights <- .632/(1-.368*relative.overfit.boost)
boot632p.boost <- (1-weights)*pec.full.boost + weights*ifelse(mean.boot.boost < pec.noinf.boost,mean.boot.boost,pec.noinf.boost)                   # .632+ estimate reference approach

relative.overfit.boost.lp <- ifelse(pec.noinf.boost.lp > pec.full.boost.lp,
                                    (ifelse(mean.boot.boost.lp < pec.noinf.boost.lp,mean.boot.boost.lp,pec.noinf.boost.lp) - pec.full.boost.lp)/
                                    (pec.noinf.boost.lp - pec.full.boost.lp),0)
weights <- .632/(1-.368*relative.overfit.boost.lp)
boot632p.boost.lp <- (1-weights)*pec.full.boost.lp + weights*ifelse(mean.boot.boost.lp < pec.noinf.boost.lp,mean.boot.boost.lp,pec.noinf.boost.lp) # .632+ estimate sequential complementary strategy


###########################################################################################################################################################################################################
## 6.2. bootstrap .632+ prediction error estimates individually for each bootstrap sample

modsubclin <- apply(boot.clin, 1, function(x) {
    temp1 <- ifelse(pec.noinf.clin > pec.full.clin,
                    (ifelse(x < pec.noinf.clin,x,pec.noinf.clin) - pec.full.clin)/
                    (pec.noinf.clin - pec.full.clin),0)
    weights <- .632/(1-.368*temp1)
    temp2 <- (1-weights)*pec.full.clin + weights*ifelse(x < pec.noinf.clin,x,pec.noinf.clin)
    temp2
})

modsubref <- apply(boot.boost, 1, function(x) {
    temp1 <- ifelse(pec.noinf.boost > pec.full.boost,
                    (ifelse(x < pec.noinf.boost,x,pec.noinf.boost) - pec.full.boost)/
                    (pec.noinf.boost - pec.full.boost),0)
    weights <- .632/(1-.368*temp1)
    temp2 <- (1-weights)*pec.full.boost + weights*ifelse(x < pec.noinf.boost,x,pec.noinf.boost)
    temp2
})

modsubseq <- apply(boot.boost.lp, 1, function(x) {
    temp1 <- ifelse(pec.noinf.boost.lp > pec.full.boost.lp,
                    (ifelse(x < pec.noinf.boost.lp,x,pec.noinf.boost.lp) - pec.full.boost.lp)/
                    (pec.noinf.boost.lp - pec.full.boost.lp),0)
    weights <- .632/(1-.368*temp1)
    temp2 <- (1-weights)*pec.full.boost.lp + weights*ifelse(x < pec.noinf.boost.lp,x,pec.noinf.boost.lp)
    temp2
})

ripec <- function(actual.error){
    y.ipec <- apply(cbind(actual.error[2:length(actual.error)],
                          actual.error[1:(length(actual.error)-1)]),1,
                    function(arg)mean(arg,na.rm = T))
    x.ipec <- eval.times[2:length(eval.times)]-eval.times[1:(length(eval.times)-1)]
    ripec <- x.ipec%*%y.ipec
    return(ripec)
}

km_ipec <- ripec(pec.full.km)
clin_ipec <- apply(modsubclin, 2, ripec)
ref_ipec <- apply(modsubref, 2, ripec)
seq_ipec <- apply(modsubseq, 2, ripec)

###########################################################################################################################################################################################################

###########################################################################################################################################################################################################
## 7. Effect of the size of overlap

### lp.response: SNP signature conditioning on the overlap samples

set.seed(109144)
overlap.15 <- sample(lp.response, size=15) # decreased overlap size of 15 biological samples

set.seed(1010143)
overlap.10 <- sample(lp.response, size=10) # decreased overlap size of 10 biological samples

###########################################################################################################################################################################################################
## 7.1. Effect of the size of overlap on the performance of the sequential complementary stratey

### Repeat for decreased size of overlap from 26 to 10 biological samples:

#### 1. SNP risk prediction model (1):
####    biological samples which are nor randomly drawn are treated as if they would have measurements from GEP data only

#### 2. Linking model (3):
####    impute the information from the SNP data for the individuals with measurements from the GEP data only

#### 3. Sequential complementary strategy:
####    risk prediction model (2) for the GEP data estimated by componentwise likelihood-based boosting incorporating the SNP signature as a fixed offset

#### 6. Prediction performance:
####    Sequential complementary strategy estimated by componentwise boosting

###########################################################################################################################################################################################################
## 7.2. Effect of the size of overlap on the performance of the linking model (2)

### Repeat for decreased size of overlap from 26 to 10 biological samples:

#### 2.1. Prediction performance: Linking model (3)

###########################################################################################################################################################################################################

###########################################################################################################################################################################################################
## 8. Swapped order of the sequential complementary strategy

###########################################################################################################################################################################################################
## 8.1. GEP risk prediction model (1) estimated by componentwise likelihood-based boosting

### ref.GEP: results from step 4, reference approach fitted on the GEP data 

lp.GEP <- ref.GEP$linear.predictor[nrow(ref.GEP$linear.predictor),]   # GEP signature
lp.GEP.con <- lp.GEP[SNP_con_GEP[,"GEP"]]                             # GEP signature conditioning on the overlap samples
names(lp.GEP.con) <- SNP_con_GEP[,"SNP"]                              # GEP signature conditioning on the overlap samples (SNP ID)

###########################################################################################################################################################################################################
## 8.2. Linking model (3) estimated by componentwise likelihood-based boosting

### xmat: covariate (SNP) matrix in addition to the five clinical covariates

lp.response <- lp.GEP.con                                                  # GEP signature conditioning on the overlap samples
totalmat_con<- xmat[names(lp.response),]                                   # totalmat conditioning on the biological samples with measurements from both sources in parallel (overlap samples)
totalmat_pred <- xmat[-which(rownames(totalmat) %in% names(lp.response)),] # totalmat conditioning on biological samples with measurements from GEP only (non-overlap samples)

lamb <- (dim(totalmat_con)[1])*((1/0.05)-1)
expr.index <- dim(totalmat_con)[2]-5

## Main tuning parameter (boosting steps) determined by cross-validation

set.seed(6320141)
res.SNP <- cv.GAMBoost(y=lp.response,
                       x.linear=totalmat_con,
                       K=10,
                       family=gaussian(),
                       standardize.linear=TRUE,
                       penalty.linear=c(rep(0,length=5),
                           rep(lamb,length=expr.index)),
                       criterion="score",
                       multicore=10)

pred.lp <- predict(object=res.SNP, newdata.linear=totalmat_pred,
                   type="response")

lp <- c(lp.response,pred.lp)                                                 # GEP signature including known GEP signatur (overlap samples) and predicted GEP signature (non-overlap samples)

###########################################################################################################################################################################################################
## 8.3. Sequential complementary strategy fitted on the original data:
##      risk prediction model (2) for the SNP data estimated by componentwise likelihood-based boosting incorporating the GEP signature as a fixed offset

### current.status: event indicator (SNP data)
### current.time  : event time (SNP data)
### xmat      : covariate (SNP) matrix in addition to the five clinical covariates

lamb <- 19*sum(current.status==1)
boost.stepno <- 500 
K <- 10

lp <- lp[rownames(xmat)]
totalmatlp <- cbind(lp,xmat)

## Main tuning parameter (boosting steps) determined by cross-validation

set.seed(70320141)
cv.seq.SNP <- cv.CoxBoost(current.time,current.status,
                          totalmatlp,unpen.index=1:6,
                          standardize=TRUE,
                          maxstepno=boost.stepno,penalty=lamb,
                          criterion="pscore",x.is.01=FALSE,
                          K=K,trace=TRUE,multicore=10)

seq.SNP <- CoxBoost(current.time,current.status,
                    totalmatlp,unpen.index=1:6,
                    standardize=TRUE,
                    stepno=cv.seq.SNP$optimal.step,
                    penalty=lamb,criterion="pscore",
                    x.is.01=FALSE)


###########################################################################################################################################################################################################
## 8.4. Reference approach fitted on the SNP data:
##      to risk prediction model (2) estimated by componentwise likelihood-based boosting without adjusting for the effects of the GEP data:

### res.full.SNP: results from step 1, SNP risk prediction model (1) fitted on the SNP data 

ref.SNP <- res.full.SNP

###########################################################################################################################################################################################################
## 8.5. Prediction performance:
##      Sequential complementary strategy (swapped order) and reference approach (risk prediction model (2)) estimated by boosting

### current.status : event indicator (SNP data)
### current.time   : event time (SNP data)
### xmat           : covariate (SNP) matrix in addition to the five clinical covariates
### cv.res.full.SNP: results from step 1, main tuning parameter for the reference approach (risk prediction model (1)) fitted on the original SNP data
### ref.SNP        : reference approach, results from step 1, SNP risk prediction model fitted on the SNP data 
### totalmatlp     : covariate (SNP) matrix in addition to the five clinical covariates and linera predictor including SNP information
### cv.seq.SNP     : Main tuning parameter for the sequential complementary strategy (risk prediction model (2)) fitted on the original SNP data
### seq.SNP        : sequential complementary strategy (risk prediction model (2)) fitted on the original SNP data
### lp.response    : GEP signature conditioning on the overlap samples

cv.ref.SNP <- cv.res.full.SNP

K <- 10
boost.stepno <- 500
eval.times <- seq(from=0,to=1000,length=100)

do.full <- TRUE
do.bootstrap <- TRUE
boot.samples <- 1:100

###########################################################################################################################################################################################################
## 8.5.1. bootstrap .632+ prediction error estimates

if (do.full) {
    km <- survfit(Surv(time,status) ~ 1, data=data.frame(time=current.time,status=current.status)) 
    
    pec.km.full <- pmpec(km,response=Surv(current.time,current.status),          # Kaplan-Meier benchmark (determined by the apparent error)
                         x=totalmatlp[,1:6],times=eval.times)

    clin.frame <- data.frame(time=current.time,status=current.status,
                             age=totalmatlp[,"age"],
                             low.risk=totalmatlp[,"low.risk"],
                             FLT3=totalmatlp[,"FLT3"],
                             NPM1=totalmatlp[,"NPM1"],
                             WBC=totalmatlp[,"WBC"])
  
    cp1.full <- coxph(Surv(time,status) ~ age + low.risk + FLT3 + NPM1 + WBC,
                      data=clin.frame)
 
    pec.clin.full <- pmpec(cp1.full,response=Surv(current.time,current.status),  # apparent error Cox model
                           x=clin.frame,times=eval.times)
  
    pec.clin.noinf <- pmpec(cp1.full,response=Surv(current.time,current.status), # no-informative error Cox model
                            x=clin.frame,times=eval.times,type="NoInf")
    
    crr.response <- cbind(current.time,current.status)
    colnames(crr.response) <- c("time","status")

    pec.boost.full <- pmpec(seq.SNP,response=crr.response,                       # apparent error reference approach
                            x=xmat,
                            times=eval.times,
                            model.args=list(complexity=cv.seq.GEP$stepno))

    pec.boost.noinf <- pmpec(seq.SNP,response=crr.response,                      # no-informtive error reference approacj
                             x=xmat, 
                             times=eval.times,
                             model.args=list(complexity=cv.seq.GEP$stepno),      
                             type="NoInf")
    
    pec.boost.full.lp <- pmpec(seq.SNP,response=crr.response,                    # apparent error sequential complementary strategy
                               x=totalmatlp,
                               times=eval.times,
                               model.args=list(complexity=cv.seq.GEP$stepno))

    pec.boost.noinf.lp <- pmpec(seq.SNP,response=crr.response,                   # no-informtive error sequential complementary strategy
                                x=totalmatlp,
                                times=eval.times,
                                model.args=list(complexity=cv.seq.GEP$stepno),   
                                type="NoInf")
}

if (do.bootstrap) {
    boot.n <- 100
    set.seed(42)
    boot.indices <- resample.indices(n = nrow(xmat),
                                     method = "sub632",sample.n=boot.n)          # without replacement
    lp_overlap <- lp.response
    crr.response <- cbind(current.time,current.status)
    colnames(crr.response) <- c("time","status")
    for (actual.boot in boot.samples){
        
        avail <- (1:nrow(xmat) %in% boot.indices$sample.index[[actual.boot]])
        test.index <- (1:nrow(xmat) %in% boot.indices$not.in.sample[[actual.boot]])
        data_avail <- xmat[avail,]
        data_test <- xmat[test.index,]

        ## Cox proportional hazards model
        clin.frame <- data.frame(time=current.time[avail],status=current.status[avail],
                                 age=data_avail[,"age"],low.risk=data_avail[,"low.risk"],
                                 FLT3=data_avail[,"FLT3"],NPM1=data_avail[,"NPM1"],
                                 WBC=data_avail[,"WBC"])

        clin.frame.test <- data.frame(time=current.time[test.index],status=current.status[test.index],
                                      age=data_test[,"age"],low.risk=data_test[,"low.risk"],
                                      FLT3=data_test[,"FLT3"],NPM1=data_test[,"NPM1"],
                                      WBC=data_test[,"WBC"])
        
        cp1.boot <- coxph(Surv(time,status) ~ age + low.risk + FLT3 + NPM1 + WBC,
                          data=clin.frame)

        pec.clin.boot <- pmpec(cp1.boot,response=Surv(current.time[test.index], #  MSE Cox proportional hazards model in the bth bootstrap
                                            current.status[test.index]),
                               x=clin.frame.test,times=eval.times)

        ## Reference approach
        penalty <- 19*sum(current.status[avail]==1)

        set.seed(22)
        cv.res.boot <- cv.CoxBoost(current.time[avail],current.status[avail],
                                   data_avail,unpen.index=1:5,
                                   standardize=TRUE,maxstepno=boost.stepno,
                                   penalty=(1-1/K)*penalty,x.is.01=FALSE,
                                   K=K,multicore=10)

        res.boot <- CoxBoost(current.time[avail],current.status[avail],
                             data_avail,unpen.index=1:5,
                             standardize=TRUE,stepno=cv.res.boot$optimal.step,
                             penalty=penalty,x.is.01=FALSE,
                             return.score=FALSE)

        pec.boost.boot <- pmpec(res.boot,response=crr.response[test.index,],    #  MSE reference approach in the bth bootstrap
                                x=data_test,times=eval.times,
                                model.args=list(complexity=cv.res.boot$optimal.step))
        
        ## Sequential complementary strategy
        overlap_avail.lp <- (names(lp_overlap) %in% rownames(data_avail))
        lp_overlap_avail <- lp_overlap[overlap_avail.lp]
        data.avail.overlap <- data_avail[names(lp_overlap_avail),]
        data.avail.non.overlap <- data_avail[-which(rownames(data_avail) %in% rownames(data.avail.overlap)),]
        
        overlap_test.lp <- (names(lp_overlap) %in% rownames(data_test))
        lp_overlap_test <- lp_overlap[overlap_test.lp]
        data.test.overlap <- data_test[names(lp_overlap_test),]
        data.test.non.overlap <- data_test[-which(rownames(data_test) %in% rownames(data.test.overlap)),]
        
        data.to.pred <- rbind(data.avail.non.overlap,data.test.non.overlap)
        
        lamb <- (nrow(data.avail.overlap))*((1/0.005)-1)
        expr.index <- ncol(data.avail.overlap)-5

        ## Linking model
        set.seed(70320145)
        fit_lp <- cv.GAMBoost(y=lp_overlap_avail,
                              x.linear=data.avail.overlap,
                              K=K,family=gaussian(),
                              standardize.linear=TRUE,
                              penalty.linear=c(rep(0,length=5),
                                  rep(lamb,length=expr.index)),
                              criterion="score",multicore=10)
                
        pred.lp.non.overlap <- predict(object=fit_lp, newdata.linear=data.to.pred,type="response")
        
        lp.all <- c(lp_overlap_avail,lp_overlap_test,pred.lp.non.overlap)
        lp.full <- lp.all[rownames(totalmat)]
        gamboost.cv.fit <- list(res.subsamples=fit_lp)
        
        lp.train <- lp.full[avail]
        lp.test <- lp.full[test.index]

        
        data_avail.lp <- cbind(lp.train,data_avail)
        colnames(data_avail.lp)[1] <- "lp"

        ## Risk prediction model (2)
        penalty <- 19*sum(current.status[avail]==1)
        
        set.seed(70320146)
        cv.res.boot.lp <- cv.CoxBoost(current.time[avail],current.status[avail],
                                      data_avail.lp,unpen.index=1:6,
                                      standardize=TRUE,maxstepno=boost.stepno,
                                      penalty=(1-1/K)*penalty,x.is.01=FALSE,
                                      K=K,multicore=10)
        
        res.boot.lp <- CoxBoost(current.time[avail],current.status[avail],
                                data_avail.lp,unpen.index=1:6,
                                standardize=TRUE,stepno=cv.res.boot.lp$optimal.step,
                                penalty=penalty,x.is.01=FALSE,return.score=FALSE)
        
        data_test.lp <- cbind(lp.test,data_test)
        colnames(data_test.lp)[1] <- "lp"
        
        pec.boost.boot.lp <- pmpec(res.boot.lp,response=crr.response[test.index,],   # MSE sequential complementary strategy in the bth bootstrap
                                   x=data_test.lp,times=eval.times,
                                   model.args=list(complexity=cv.res.boot.lp$optimal.step))
        
        boot.res <- list(res.subsample=res.boot,res.subsample.lp=res.boot.lp,GAMB.fit=gamboost.cv.fit)

        save(boot.res,pec.clin.boot,pec.boost.boot,pec.boost.boot.lp,file=paste("pred_SNP_resultII",actual.boot,".RData",sep=""))
    }
}

boot.boost <- NULL
boot.boost.lp <- NULL
boot.clin <- NULL

for (i in c(1:100)) {
   load(paste("pred_SNP_resultII",i,".RData",sep=""))
   boot.clin <- rbind(boot.clin,pec.clin.boot)
   boot.boost <- rbind(boot.boost,pec.boost.boot)
   boot.boost.lp <- rbind(boot.boost.lp,pec.boost.boot.lp)
}

mean.boot.clin <- colMeans(boot.clin)                                                # bootstrap error Cox proportional hazards model
mean.boot.boost <- colMeans(boot.boost)                                              # bootstrap error reference approach
mean.boot.boost.lp <- colMeans(boot.boost.lp)                                        # bootstrap error sequential complementary strategy

relative.overfit.clin <- ifelse(pec.noinf.clin > pec.full.clin,
                                (ifelse(mean.boot.clin < pec.noinf.clin,mean.boot.clin,pec.noinf.clin) - pec.full.clin)/
                                (pec.noinf.clin - pec.full.clin),0)
weights <- .632/(1-.368*relative.overfit.clin)
boot632p.clin <- (1-weights)*pec.full.clin + weights*ifelse(mean.boot.clin < pec.noinf.clin,mean.boot.clin,pec.noinf.clin)                         # .632+ estimate Cox proportional hazards model

relative.overfit.boost <- ifelse(pec.noinf.boost > pec.full.boost,
                                 (ifelse(mean.boot.boost < pec.noinf.boost,mean.boot.boost,pec.noinf.boost) - pec.full.boost)/
                                 (pec.noinf.boost - pec.full.boost),0)
weights <- .632/(1-.368*relative.overfit.boost)
boot632p.boost <- (1-weights)*pec.full.boost + weights*ifelse(mean.boot.boost < pec.noinf.boost,mean.boot.boost,pec.noinf.boost)                   # .632+ estimate reference approach

relative.overfit.boost.lp <- ifelse(pec.noinf.boost.lp > pec.full.boost.lp,
                                    (ifelse(mean.boot.boost.lp < pec.noinf.boost.lp,mean.boot.boost.lp,pec.noinf.boost.lp) - pec.full.boost.lp)/
                                    (pec.noinf.boost.lp - pec.full.boost.lp),0)
weights <- .632/(1-.368*relative.overfit.boost.lp)
boot632p.boost.lp <- (1-weights)*pec.full.boost.lp + weights*ifelse(mean.boot.boost.lp < pec.noinf.boost.lp,mean.boot.boost.lp,pec.noinf.boost.lp) # .632+ estimate sequential complementary strategy


#######################################################################################################################################################################################################
## 8.5.2 bootstrap .632+ prediction error estimates individually for each bootstrap sample

modsubclin <- apply(boot.clin, 1, function(x) {
    temp1 <- ifelse(pec.noinf.clin > pec.full.clin,
                    (ifelse(x < pec.noinf.clin,x,pec.noinf.clin) - pec.full.clin)/
                    (pec.noinf.clin - pec.full.clin),0)
    weights <- .632/(1-.368*temp1)
    temp2 <- (1-weights)*pec.full.clin + weights*ifelse(x < pec.noinf.clin,x,pec.noinf.clin)
    temp2
})

modsubref <- apply(boot.boost, 1, function(x) {
    temp1 <- ifelse(pec.noinf.boost > pec.full.boost,
                    (ifelse(x < pec.noinf.boost,x,pec.noinf.boost) - pec.full.boost)/
                    (pec.noinf.boost - pec.full.boost),0)
    weights <- .632/(1-.368*temp1)
    temp2 <- (1-weights)*pec.full.boost + weights*ifelse(x < pec.noinf.boost,x,pec.noinf.boost)
    temp2
})

modsubseq <- apply(boot.boost.lp, 1, function(x) {
    temp1 <- ifelse(pec.noinf.boost.lp > pec.full.boost.lp,
                    (ifelse(x < pec.noinf.boost.lp,x,pec.noinf.boost.lp) - pec.full.boost.lp)/
                    (pec.noinf.boost.lp - pec.full.boost.lp),0)
    weights <- .632/(1-.368*temp1)
    temp2 <- (1-weights)*pec.full.boost.lp + weights*ifelse(x < pec.noinf.boost.lp,x,pec.noinf.boost.lp)
    temp2
})

ripec <- function(actual.error){
    y.ipec <- apply(cbind(actual.error[2:length(actual.error)],
                          actual.error[1:(length(actual.error)-1)]),1,
                    function(arg)mean(arg,na.rm = T))
    x.ipec <- eval.times[2:length(eval.times)]-eval.times[1:(length(eval.times)-1)]
    ripec <- x.ipec%*%y.ipec
    return(ripec)
}

km_ipec <- ripec(pec.full.km)
clin_ipec <- apply(modsubclin, 2, ripec)
ref_ipec <- apply(modsubref, 2, ripec)
seq_ipec <- apply(modsubseq, 2, ripec)

#######################################################################################################################################################################################################
