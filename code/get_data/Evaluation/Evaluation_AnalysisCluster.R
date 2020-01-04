########################################################

# NOTE: Before the code can be excecuted, the R working directory *MUST* 
# be set to the directory that contains the folder 'Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION':

# Remove '#' from the two lines below and replace 'Z:/here/is/my/path/' by the path
# to the directory that contains 'Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION' ( the path 'Z:/here/is/my/path/'
# *MUST* end with a slash / ):

# basedir <- "Z:/here/is/my/path/"
# setwd(basedir)

########################################################


# Set working directory:

setwd("Z:/Projects/BlockForests/")



# Load in and prepare the results:


load("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Results/scenariogrid1.Rda")
scenariogrid1 <- scenariogrid

allinds <- 1:nrow(scenariogrid1)

allresults <- list.files("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Results/Results1")
indres <- gsub("res", "", allresults)
indres1 <- sort(as.numeric(gsub(".Rda", "", indres)))

(anymissing <- length(setdiff(allinds, indres1)) > 0)

# --> No missing iterations for this part.



load("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Results/scenariogrid2.Rda")
scenariogrid2 <- scenariogrid

allinds <- 1:nrow(scenariogrid2)

allresults <- list.files("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Results/Results2")
indres <- gsub("res", "", allresults)
indres2 <- sort(as.numeric(gsub(".Rda", "", indres)))

(anymissing <- length(setdiff(allinds, indres2)) > 0)

(scenmiss <- scenariogrid2[setdiff(allinds, indres2),])

# --> Missing results only for the data set PRAD that was
#     excluded from the evaluation (see paper for details):

scenariogrid2 <- scenariogrid2[indres2,]



load("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Results/scenariogrid3.Rda")
scenariogrid3 <- scenariogrid

allinds <- 1:nrow(scenariogrid3)

allresults <- list.files("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Results/Results3")
indres <- gsub("res", "", allresults)
indres3 <- sort(as.numeric(gsub(".Rda", "", indres)))

(anymissing <- length(setdiff(allinds, indres3)) > 0)





scenariogrid <- rbind(scenariogrid1, scenariogrid2, scenariogrid3)


Results <- list()

count <- 1

counttemp <- 1
for(i in 1:nrow(scenariogrid1)) {
  load(paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Results/Results1/res", indres1[counttemp], ".Rda", sep=""))
  Results[[count]] <- res
  count <- count + 1
  counttemp <- counttemp + 1
}

counttemp <- 1
for(i in 1:nrow(scenariogrid2)) {
  load(paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Results/Results2/res", indres2[counttemp], ".Rda", sep=""))
  Results[[count]] <- res
  count <- count + 1
  counttemp <- counttemp + 1
}

counttemp <- 1
for(i in 1:nrow(scenariogrid3)) {
  load(paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Results/Results3/res", indres3[counttemp], ".Rda", sep=""))
  Results[[count]] <- res
  count <- count + 1
  counttemp <- counttemp + 1
}


scenariogrid$method[scenariogrid$method=="randomsurvivalforest"] <- "RSF"



ytrue <- lapply(rapply(sapply(Results, function(x) x$ytrue), enquote, how="unlist"), eval)
ytruestatus <- lapply(rapply(sapply(Results, function(x) x$ytruestatus), enquote, how="unlist"), eval)

riskpreds <- lapply(rapply(sapply(Results, function(x) x$riskpreds), enquote, how="unlist"), eval)


paramvalues <- sapply(Results, function(x) x$paramvalues)
zeroparamvalues <- sapply(paramvalues, function(x) length(x)==0)
paramvalues[zeroparamvalues] <- as.list(rep(NA, sum(zeroparamvalues)))
paramvalues <- lapply(rapply(paramvalues, enquote, how="unlist"), eval)


mtrys <- sapply(Results, function(x) x$mtrys)
zeromtrys <- sapply(mtrys, function(x) length(x)==0)
mtrys[zeromtrys] <- as.list(rep(NA, sum(zeromtrys)))
mtrys <- lapply(rapply(mtrys, enquote, how="unlist"), eval)






# Compute the C index values from the raw results:


# First, check the cases for which the C index is not computable
# or result in 'NA':

options(warn=2)

cs <- list()

library("survcomp")
for(i in seq(along=ytrue)) {
  
  cs[[i]] <- try(concordance.index(x=riskpreds[[i]], surv.time=ytrue[[i]], surv.event=ytruestatus[[i]])$c.index)
  
}

options(warn=0)


table(sapply(cs, class))
table(sapply(cs, is.na))


all(scenariogrid[sapply(cs, class)=="try-error",]$dat=="PRAD.Rda")

# --> All iterations for which the C index could not be computed
# were obtained with the data set PRAD that was excluded from the analysis.



scenariogrid[sapply(cs, is.na) & scenariogrid$dat!="PRAD.Rda",]

nrow(scenariogrid[sapply(cs, is.na) & scenariogrid$dat!="PRAD.Rda",])
nrow(scenariogrid[scenariogrid$dat!="PRAD.Rda",])

nrow(scenariogrid[sapply(cs, is.na) & scenariogrid$dat!="PRAD.Rda",])/
  nrow(scenariogrid[scenariogrid$dat!="PRAD.Rda",])

# --> The cases for which the computed C index was 'NA' (and not associated
# with the data set PRDA) were very rare.




# --> Calculate the C index values:

library("survcomp")
results <- scenariogrid
results$cindex <- mapply(function(x, y, z) concordance.index(x=x, surv.time=y, surv.event=z)$c.index, riskpreds, ytrue, ytruestatus)
results$mtry <- unlist(mtrys)

results$seed <- NULL



# Calculate mean C index values per data set:

library("plyr")
resultsum <- ddply(results, .variables=c("method", "dat", "cvind"), .fun=summarise, cindex=mean(cindex, na.rm=TRUE))

resultsum$method <- factor(resultsum$method, levels=c("RSF", 
                                                      "VarProb",
                                                      "SplitWeights",
                                                      "BlockVarSel",
                                                      "RandomBlock",
                                                      "BlockForest"))

datasets <- sort(unique(resultsum$dat))






# Calculate the mean C index values across data sets:

resultsumsum <- ddply(resultsum, .variables=c("method", "dat"), .fun=summarise, cindex=mean(cindex))

resultswide <- reshape(resultsumsum, idvar = c("dat"), timevar = "method", direction = "wide")
resultswide

resultswide2 <- resultswide[,-1]




# Result obtained with the data set PRAD that was excluded from
# the results:

resultswide[resultswide$dat=="PRAD.Rda",]



# Exclude data set PRAD:

results_withPRAD <- results
resultsum_withPRAD <- resultsum
resultsumsum_withPRAD <- resultsumsum
resultswide_withPRAD <- resultswide
resultswide2_withPRAD <- resultswide2
paramvalues_withPRAD <- paramvalues
mtrys_withPRAD <- mtrys

results <- results[results$dat!="PRAD.Rda",]
resultsum <- resultsum[resultsum$dat!="PRAD.Rda",]
resultsumsum <- resultsumsum[resultsumsum$dat!="PRAD.Rda",]
resultsum <- resultsum[resultsum$dat!="PRAD.Rda",]
resultswide <- resultswide[resultswide$dat!="PRAD.Rda",]
resultswide2 <- resultswide[,-1]
paramvalues <- paramvalues[results_withPRAD$dat!="PRAD.Rda"]
mtrys <- mtrys[results_withPRAD$dat!="PRAD.Rda"]





# Calculate ranks of the methods:


resultranks <- t(apply(resultswide2, 1, function(x) rank(-x)))

colnames(resultranks) <- gsub("cindex.", "", colnames(resultranks))



resultranks2 <- reshape(data.frame(resultranks), varying=colnames(resultranks), 
                        v.names="rank", 
                        timevar="method", times=colnames(resultranks),
                        direction="long")
resultranks2$method <- factor(resultranks2$method, levels=c("RSF", 
                                                            "VarProb",
                                                            "SplitWeights",
                                                            "BlockVarSel",
                                                            "RandomBlock",
                                                            "BlockForest"))

datasets <- datasets[datasets!="PRAD.Rda"]







# Figure: Multi-omics data: C index values obtained for the individual repetitions 
# of the cross-validation separately for each data set and method:


resultsumtemp <- resultsum
datasetstemp <- datasets

resultsumtemp$dat <- gsub(".Rda", "", resultsumtemp$dat)
datasetstemp <- gsub(".Rda", "", datasetstemp)

library("ggplot2")
p <- ggplot(data=resultsumtemp[resultsumtemp$dat %in% datasetstemp[1:10],], aes(x=method, y=cindex)) + facet_wrap(~dat, ncol=2, scales = 'free_y') + 
  geom_boxplot() + theme_bw() + theme(strip.text.x=element_text(size=14)) + labs(y="C index values") + 
  theme(axis.text.x=element_text(size=14, angle=45, vjust=1, hjust=1, colour="black"),
        axis.text.y=element_text(size=11, colour="black"), axis.title.x=element_blank(),
        axis.title.y=element_text(size=14))
p

ggsave("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Figures/Results_AnalysisCluster_1.pdf", width=6.5*1.35, height=8*1.35)




library("ggplot2")
p <- ggplot(data=resultsumtemp[resultsumtemp$dat %in% datasetstemp[11:20],], aes(x=method, y=cindex)) + facet_wrap(~dat, ncol=2, scales = 'free_y') + 
  geom_boxplot() + theme_bw() + theme(strip.text.x=element_text(size=14)) + labs(y="C index values") + 
  theme(axis.text.x=element_text(size=14, angle=45, vjust=1, hjust=1, colour="black"),
        axis.text.y=element_text(size=11, colour="black"), axis.title.x=element_blank(),
        axis.title.y=element_text(size=14))
p

ggsave("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Figures/Results_AnalysisCluster_2.pdf", width=6.5*1.35, height=8*1.35)





# Mean ranks:


round(sort(colMeans(resultranks)), 2)




# one-sided paired Student's t-tests to compare the performance each variant against RSF:


variants <- c("cindex.VarProb", "cindex.SplitWeights", "cindex.BlockVarSel", "cindex.RandomBlock", "cindex.BlockForest")

ps <- sapply(variants, function(x) t.test(resultswide[, x], resultswide[,"cindex.RSF"], paired=TRUE, alternative="greater")$p.value)

# Adjust the p values:
p.adjust(p=ps, method = "holm")




# Both, BlockForest and RandomBlock perform better than RSF for
# the same 14 data sets:


table(resultswide2[,"cindex.RandomBlock"] > resultswide2[,"cindex.RSF"],
      resultswide2[,"cindex.BlockForest"] > resultswide2[,"cindex.RSF"])




# Figure: Multi-omics data: Performances of all six considered methods.


library("RColorBrewer")

resultsumsumtemp <- resultsumsum
resultsumsumtemp$dat

datasetlevels <- resultsumsumtemp[resultsumsumtemp$method=="RSF",]$dat[order(resultsumsumtemp[resultsumsumtemp$method=="RSF",]$cindex)]

resultsumsumtemp$dat <- factor(resultsumsumtemp$dat, levels=datasetlevels)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

p1 <- ggplot(data=resultsumsumtemp, aes(x=method, y=cindex)) + theme_bw() + 
  geom_line(aes(group=dat, color=dat, linetype=dat)) + ylab("Mean C index values") + 
  scale_linetype_manual(values = c(rep("12345678", 5), rep("solid", 5), rep("dotted", 5), rep("dashed", 5))) + # scale_linetype_manual(values = c(rep("solid", 5), rep("longdash", 5), rep("dotted", 5), rep("dotdash", 5))) +
  scale_color_manual(values = rep(gg_color_hue(5), times=4)) +
  geom_boxplot(fill="transparent") + ggtitle("    a") +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=14, angle=45, vjust=1, hjust=1, colour="black"),
        axis.text.y=element_text(size=11, colour="black"), axis.title=element_text(size=14), legend.position="none",
        plot.title = element_text(size = 25, face = "bold"))
p1

resultsumsumtemp2 <- resultsumsumtemp 
datstemp <- unique(resultsumsumtemp$dat)
for(i in seq(along=datstemp)) {
  rsfc <- resultsumsumtemp2$cindex[resultsumsumtemp2$method=="RSF" & resultsumsumtemp2$dat==datstemp[i]]
  resultsumsumtemp2$cindex[resultsumsumtemp2$dat==datstemp[i]] <- resultsumsumtemp2$cindex[resultsumsumtemp2$dat==datstemp[i]] - rsfc
}

resultsumsumtemp2$dat <- factor(resultsumsumtemp2$dat, levels=datasetlevels)

p2 <- ggplot(data=resultsumsumtemp2, aes(x=method, y=cindex)) + theme_bw() + 
  geom_line(aes(group=dat, color=dat, linetype=dat)) + ylab("Differences between the mean C index values obtained\nfor the variants and the reference method RSF") + 
  geom_hline(yintercept=0, col="darkgrey", size=1) +
  scale_linetype_manual(values = c(rep("12345678", 5), rep("solid", 5), rep("dotted", 5), rep("dashed", 5))) + # scale_linetype_manual(values = c(rep("solid", 5), rep("longdash", 5), rep("dotted", 5), rep("dotdash", 5))) +
  scale_color_manual(values = rep(gg_color_hue(5), times=4)) +
  geom_boxplot(fill="transparent") + ggtitle("b") +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=14, angle=45, vjust=1, hjust=1, colour="black"),
        axis.text.y=element_text(size=11, colour="black"), axis.title=element_text(size=14), legend.position="none",
        plot.title = element_text(size = 25, face = "bold"))
p2

set.seed(1234)
p3 <- ggplot(data=resultranks2, aes(x=method, y=rank)) + theme_bw() + geom_boxplot() + 
  ylab("Data set specific ranks") + ggtitle("    c") +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=14, angle=45, vjust=1, hjust=1, colour="black"),
        axis.text.y=element_text(size=11, colour="black"), axis.title=element_text(size=14),
        plot.title = element_text(size = 25, face = "bold"))
p3

library("grid")
library("gridExtra")
pall <- grid.arrange(p1, p2, p3, ncol = 1)

grid.draw(pall)

ggsave(file="./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Figures/Figure1.pdf", plot=pall, width=7.5, height=14)






# Figure: Multi-omics data: Performances of BlockForest and RandomBlock relative to that of RSF.


dataordered <- resultswide$dat[order(resultswide[,"cindex.BlockForest"] - resultswide[,"cindex.RSF"])]
dataordered <- gsub(".Rda", "", dataordered)

dataordered <- factor(dataordered, levels=dataordered)

datatemp <- data.frame(x=dataordered, y1=sort(resultswide[,"cindex.BlockForest"] - resultswide[,"cindex.RSF"]),
                       y2=sapply(paste(as.character(dataordered), ".Rda", sep=""), function(x) resultswide$cindex.RandomBlock[resultswide$dat==x] -
                                   resultswide$cindex.RSF[resultswide$dat==x]))

datatemp <- reshape(datatemp, varying=c("y1", "y2"), 
                    v.names="difference", 
                    timevar="Method", times=c("BlockForest", "RandomBlock"),
                    direction="long")


p <- ggplot(data=datatemp, aes(x=x, y=difference)) + theme_bw() + geom_hline(yintercept=0, col="grey", size=1) + 
  geom_point(aes(color=Method)) + labs(x="data set", y="Differences between mean C index values\nobtained for BlockForest / RandomBlock and RSF") +  
  theme(axis.text.x=element_text(size=12, angle=45, vjust=1, hjust=1, colour="black"),
        axis.text.y=element_text(size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=14),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14))
p

ggsave("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Figures/Figure2.pdf", width=12*0.85, height=6.5*0.85)










# Figures: Multi-omics data: tuning parameter values optimized for the different variants:


resultsparamvalues <- results

cvalbool <- sapply(paramvalues, function(x) all(!is.na(unlist(x))))
resultsparamvalues <- results[cvalbool,]
paramvaluesnoNA <- paramvalues[cvalbool]



alldata <- sort(unique(resultsparamvalues$dat))


# The following code is commented out, because for
# this part the processed data files would be necessary,
# which are, however, not available with the Electronic
# Appendix (available upon request).
# However, the result of this part is stored as an
# Rda files in the subfolder "Results" of the Electronic
# Appendix, which is why the code following this out-commented
# part is still executable.


# blocknamesavall <- list()
# 
# try(rm(clin))
# try(rm(mirna))
# try(rm(mutation))
# try(rm(cnv))
# try(rm(rna))
# gc()
# 
# for(i in seq(along=alldata)) {
#   
#   load(paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/ProcessedData/", alldata[i], sep=""))
#   
#   blocknames <- c("clin", 
#                   "mirna",
#                   "mutation",
#                   "cnv",
#                   "rna")
#   
#   blocknamesav <-  c("")
#   
#   if(class(try(ncol(clin)))!="try-error")
#     blocknamesav <- c(blocknamesav, "clin")
#   if(class(try(ncol(mirna)))!="try-error")
#     blocknamesav <- c(blocknamesav, "mirna")
#   if(class(try(ncol(mutation)))!="try-error")
#     blocknamesav <- c(blocknamesav, "mutation")
#   if(class(try(ncol(cnv)))!="try-error")
#     blocknamesav <- c(blocknamesav, "cnv")
#   if(class(try(ncol(rna)))!="try-error")
#     blocknamesav <- c(blocknamesav, "rna")  
#   
#   blocknamesav <- blocknamesav[-1]
#   
#   blocknamesavall[[i]] <- blocknamesav
# 
#   try(rm(clin))
#   try(rm(mirna))
#   try(rm(mutation))
#   try(rm(cnv))
#   try(rm(rna))
#   gc()
#   
#   cat(i, "\n")
#   
# }
# 
# save(blocknamesavall, file="./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Results/blocknamesavall.Rda")


# Load the result of the out-commented part above:

load("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Results/blocknamesavall.Rda")


availmat <- t(sapply(blocknamesavall, function(x) blocknamesavall[[1]] %in% x))
colnames(availmat) <- blocknamesavall[[1]]


resultsparamvalues$cval_clin <- NA
resultsparamvalues$cval_mirna <- NA
resultsparamvalues$cval_mutation <- NA 
resultsparamvalues$cval_cnv <- NA
resultsparamvalues$cval_rna <- NA

for(i in seq(along=alldata)) {
  
  dataav <- availmat[i,]
  
  resultsparamvalues[resultsparamvalues$dat==alldata[i], c("cval_clin", "cval_mirna", 
                                                           "cval_mutation", "cval_cnv", 
                                                           "cval_rna")][,dataav] <- 
    matrix(nrow=length(paramvaluesnoNA[resultsparamvalues$dat==alldata[i]]), 
           data=unlist(paramvaluesnoNA[resultsparamvalues$dat==alldata[i]]), byrow=TRUE)
  
}





load("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/numbercovariates.Rda")

numbercovariates


numbercovariates <- numbercovariates[rownames(numbercovariates)!="PRAD.Rda",]

for(i in 1:nrow(numbercovariates)) {
  ncovs <- numbercovariates[i,]
  restemp <- resultsparamvalues[resultsparamvalues$dat==rownames(numbercovariates)[i] & 
                                  resultsparamvalues$method=="VarProb",][,c("cval_clin", "cval_mirna", "cval_mutation", "cval_cnv", "cval_rna")]
  resultsparamvalues[resultsparamvalues$dat==rownames(numbercovariates)[i] &
                       resultsparamvalues$method=="VarProb",][,c("cval_clin", "cval_mirna", "cval_mutation", "cval_cnv", "cval_rna")] <-
    t(apply(restemp, 1, function(x) x/sum(ncovs*x, na.rm=TRUE)))
  
}



methods <- c("VarProb",
             "SplitWeights", 
             "BlockVarSel",
             "RandomBlock", 
             "BlockForest")


alldatatemp <- gsub(".Rda", "", alldata) 



i <- 1

resultsparamvaluestemp <- resultsparamvalues[resultsparamvalues$method==methods[i],]
resultsparamvaluestempsafe <- resultsparamvaluestemp

resultsparamvaluestemp <- reshape(resultsparamvaluestemp, varying=c("cval_clin", "cval_mirna", "cval_mutation", "cval_cnv", "cval_rna"), 
                                  v.names="cval", 
                                  timevar="block", times=c("clin", "mirna", "mutation", "cnv", "rna"),
                                  direction="long")

resultsparamvaluestemp$dat <- gsub(".Rda", "", resultsparamvaluestemp$dat)

library("ggplot2")
p <- ggplot(data=resultsparamvaluestemp[resultsparamvaluestemp$dat %in% alldatatemp[1:10],], aes(x=block, y=cval)) + facet_wrap(~dat, ncol=2, scales = 'free_y') + 
  geom_boxplot() + theme_bw() + theme(strip.text.x=element_text(size=14)) + labs(x="Block", y="Tuning parameter values") + 
  theme(axis.text.x=element_text(size=14, angle=45, vjust=1, hjust=1, colour="black"),
        axis.text.y=element_text(size=11, colour="black"), axis.title=element_text(size=14))
p

ggsave(paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Figures/TunParam_", methods[i], "_1.pdf", sep=""), width=6*1.4, height=8*1.4)


p <- ggplot(data=resultsparamvaluestemp[resultsparamvaluestemp$dat %in% alldatatemp[11:20],], aes(x=block, y=cval)) + facet_wrap(~dat, ncol=2, scales = 'free_y') + 
  geom_boxplot() + theme_bw() + theme(strip.text.x=element_text(size=14)) + labs(x="Block", y="Tuning parameter values") + 
  theme(axis.text.x=element_text(size=14, angle=45, vjust=1, hjust=1, colour="black"),
        axis.text.y=element_text(size=11, colour="black"), axis.title=element_text(size=14))
p

ggsave(paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Figures/TunParam_", methods[i], "_2.pdf", sep=""), width=6*1.4, height=7*1.4)




i <- 2

resultsparamvaluestemp <- resultsparamvalues[resultsparamvalues$method==methods[i],]
resultsparamvaluestempsafe <- resultsparamvaluestemp

resultsparamvaluestemp <- reshape(resultsparamvaluestemp, varying=c("cval_clin", "cval_mirna", "cval_mutation", "cval_cnv", "cval_rna"), 
                                  v.names="cval", 
                                  timevar="block", times=c("clin", "mirna", "mutation", "cnv", "rna"),
                                  direction="long")

resultsparamvaluestemp$dat <- gsub(".Rda", "", resultsparamvaluestemp$dat)

library("ggplot2")
p <- ggplot(data=resultsparamvaluestemp[resultsparamvaluestemp$dat %in% alldatatemp[1:10],], aes(x=block, y=cval)) + facet_wrap(~dat, ncol=2, scales = 'free_y') + 
  geom_boxplot() + theme_bw() + theme(strip.text.x=element_text(size=14)) + labs(x="Block", y="Tuning parameter values") + 
  theme(axis.text.x=element_text(size=14, angle=45, vjust=1, hjust=1, colour="black"),
        axis.text.y=element_text(size=11, colour="black"), axis.title=element_text(size=14))
p

ggsave(paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Figures/TunParam_", methods[i], "_1.pdf", sep=""), width=6*1.4, height=8*1.4)


p <- ggplot(data=resultsparamvaluestemp[resultsparamvaluestemp$dat %in% alldatatemp[11:20],], aes(x=block, y=cval)) + facet_wrap(~dat, ncol=2, scales = 'free_y') + 
  geom_boxplot() + theme_bw() + theme(strip.text.x=element_text(size=14)) + labs(x="Block", y="Tuning parameter values") + 
  theme(axis.text.x=element_text(size=14, angle=45, vjust=1, hjust=1, colour="black"),
        axis.text.y=element_text(size=11, colour="black"), axis.title=element_text(size=14))
p

ggsave(paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Figures/TunParam_", methods[i], "_2.pdf", sep=""), width=6*1.4, height=7*1.4)




i <- 3

resultsparamvaluestemp <- resultsparamvalues[resultsparamvalues$method==methods[i],]
resultsparamvaluestempsafe <- resultsparamvaluestemp

resultsparamvaluestemp <- reshape(resultsparamvaluestemp, varying=c("cval_clin", "cval_mirna", "cval_mutation", "cval_cnv", "cval_rna"), 
                                  v.names="cval", 
                                  timevar="block", times=c("clin", "mirna", "mutation", "cnv", "rna"),
                                  direction="long")

resultsparamvaluestemp$dat <- gsub(".Rda", "", resultsparamvaluestemp$dat)

library("ggplot2")
p <- ggplot(data=resultsparamvaluestemp[resultsparamvaluestemp$dat %in% alldatatemp[1:10],], aes(x=block, y=cval)) + facet_wrap(~dat, ncol=2, scales = 'free_y') + 
  geom_boxplot() + theme_bw() + theme(strip.text.x=element_text(size=14)) + labs(x="Block", y="Tuning parameter values") + 
  theme(axis.text.x=element_text(size=14, angle=45, vjust=1, hjust=1, colour="black"),
        axis.text.y=element_text(size=11, colour="black"), axis.title=element_text(size=14))
p

ggsave(paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Figures/TunParam_", methods[i], "_1.pdf", sep=""), width=6*1.4, height=8*1.4)


p <- ggplot(data=resultsparamvaluestemp[resultsparamvaluestemp$dat %in% alldatatemp[11:20],], aes(x=block, y=cval)) + facet_wrap(~dat, ncol=2, scales = 'free_y') + 
  geom_boxplot() + theme_bw() + theme(strip.text.x=element_text(size=14)) + labs(x="Block", y="Tuning parameter values") + 
  theme(axis.text.x=element_text(size=14, angle=45, vjust=1, hjust=1, colour="black"),
        axis.text.y=element_text(size=11, colour="black"), axis.title=element_text(size=14))
p

ggsave(paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Figures/TunParam_", methods[i], "_2.pdf", sep=""), width=6*1.4, height=7*1.4)





i <- 4

resultsparamvaluestemp <- resultsparamvalues[resultsparamvalues$method==methods[i],]
resultsparamvaluestempsafe <- resultsparamvaluestemp

resultsparamvaluestemp <- reshape(resultsparamvaluestemp, varying=c("cval_clin", "cval_mirna", "cval_mutation", "cval_cnv", "cval_rna"), 
                                  v.names="cval", 
                                  timevar="block", times=c("clin", "mirna", "mutation", "cnv", "rna"),
                                  direction="long")

resultsparamvaluestemp$dat <- gsub(".Rda", "", resultsparamvaluestemp$dat)

library("ggplot2")
p <- ggplot(data=resultsparamvaluestemp[resultsparamvaluestemp$dat %in% alldatatemp[1:10],], aes(x=block, y=cval)) + facet_wrap(~dat, ncol=2, scales = 'free_y') + 
  geom_boxplot() + theme_bw() + theme(strip.text.x=element_text(size=14)) + labs(x="Block", y="Tuning parameter values") + 
  theme(axis.text.x=element_text(size=14, angle=45, vjust=1, hjust=1, colour="black"),
        axis.text.y=element_text(size=11, colour="black"), axis.title=element_text(size=14))
p

ggsave(paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Figures/TunParam_", methods[i], "_1.pdf", sep=""), width=6*1.4, height=8*1.4)


p <- ggplot(data=resultsparamvaluestemp[resultsparamvaluestemp$dat %in% alldatatemp[11:20],], aes(x=block, y=cval)) + facet_wrap(~dat, ncol=2, scales = 'free_y') + 
  geom_boxplot() + theme_bw() + theme(strip.text.x=element_text(size=14)) + labs(x="Block", y="Tuning parameter values") + 
  theme(axis.text.x=element_text(size=14, angle=45, vjust=1, hjust=1, colour="black"),
        axis.text.y=element_text(size=11, colour="black"), axis.title=element_text(size=14))
p

ggsave(paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Figures/TunParam_", methods[i], "_2.pdf", sep=""), width=6*1.4, height=7*1.4)



# Additional analysis of the block selection probabilites optimized using RandomBlock:

library("plyr")
resultsparamvaluestempsum <- ddply(resultsparamvaluestemp, .variables=c("method", "dat", "block"), .fun=summarise, cval=mean(cval, na.rm=TRUE))

p <- ggplot(data=resultsparamvaluestempsum, aes(x=block, y=cval)) + theme_bw() + geom_line(aes(group=dat))
p

# --> Either the optimized block selection probability of the mutation block is large or that
# of the RNA block. --> Strong negative correlation between these two block selection
# probabilities.



# Mean block selection probabilities:

resultsparamvalueswide <- reshape(resultsparamvaluestempsum[,-1], idvar = c("dat"), timevar = "block", direction = "wide")

round(sort(colMeans(resultsparamvalueswide[,-1], na.rm=TRUE), decreasing=TRUE), 2)



# Correlations between the mean block selection probabilities calculated
# per data set:

cors <- apply(combn(1:ncol(resultsparamvalueswide[,-1]), 2), 2, function(x) cor(resultsparamvalueswide[,-1][,x][complete.cases(resultsparamvalueswide[,-1][,x]),])[1,2])
names(cors) <- apply(combn(1:ncol(resultsparamvalueswide[,-1]), 2), 2, function(x) {
  so <- names(resultsparamvalueswide[,-1])[x]
  so <- gsub("cval.", "", so)
  paste(so, collapse="_")
})

sort(cors)





i <- 5

resultsparamvaluestemp <- resultsparamvalues[resultsparamvalues$method==methods[i],]
resultsparamvaluestempsafe <- resultsparamvaluestemp

resultsparamvaluestemp <- reshape(resultsparamvaluestemp, varying=c("cval_clin", "cval_mirna", "cval_mutation", "cval_cnv", "cval_rna"), 
                                  v.names="cval", 
                                  timevar="block", times=c("clin", "mirna", "mutation", "cnv", "rna"),
                                  direction="long")

resultsparamvaluestemp$dat <- gsub(".Rda", "", resultsparamvaluestemp$dat)

library("ggplot2")
p <- ggplot(data=resultsparamvaluestemp[resultsparamvaluestemp$dat %in% alldatatemp[1:10],], aes(x=block, y=cval)) + facet_wrap(~dat, ncol=2, scales = 'free_y') + 
  geom_boxplot() + theme_bw() + theme(strip.text.x=element_text(size=14)) + labs(x="Block", y="Tuning parameter values") + 
  theme(axis.text.x=element_text(size=14, angle=45, vjust=1, hjust=1, colour="black"),
        axis.text.y=element_text(size=11, colour="black"), axis.title=element_text(size=14))
p

ggsave(paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Figures/TunParam_", methods[i], "_1.pdf", sep=""), width=6*1.4, height=8*1.4)


p <- ggplot(data=resultsparamvaluestemp[resultsparamvaluestemp$dat %in% alldatatemp[11:20],], aes(x=block, y=cval)) + facet_wrap(~dat, ncol=2, scales = 'free_y') + 
  geom_boxplot() + theme_bw() + theme(strip.text.x=element_text(size=14)) + labs(x="Block", y="Tuning parameter values") + 
  theme(axis.text.x=element_text(size=14, angle=45, vjust=1, hjust=1, colour="black"),
        axis.text.y=element_text(size=11, colour="black"), axis.title=element_text(size=14))
p

ggsave(paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Figures/TunParam_", methods[i], "_2.pdf", sep=""), width=6*1.4, height=7*1.4)









# Figure: Multi-omics data: Differences between the mean C index values obtained using BlockForest
# and that obtained using RSF plotted against the values of 'n' (left panel), 
# 'oneblockimp' (middle panel), and 'signal' (right panel).


alldata <- sort(unique(resultsumsum$dat))


# The following code is commented out, because for
# this part the processed data files would be necessary,
# which are, however, not available with the Electronic
# Appendix (available upon request).
# However, the result of this part is stored as an
# Rda files in the subfolder "Results" of the Electronic
# Appendix, which is why the code following this out-commented
# part is still executable.


# datasetlist <- list()
# blocknamesavall <- list()
# 
# try(rm(clin))
# try(rm(mirna))
# try(rm(mutation))
# try(rm(cnv))
# try(rm(rna))
# gc()
# 
# for(i in seq(along=alldata)) {
#   
#   
#   load(paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/ProcessedData/", alldata[i], sep=""))
#   
#   blocknames <- c("clin", 
#                   "mirna",
#                   "mutation",
#                   "cnv",
#                   "rna")
#   
#   blocknamesav <-  c("")
#   
#   if(class(try(ncol(clin)))!="try-error")
#     blocknamesav <- c(blocknamesav, "clin")
#   if(class(try(ncol(mirna)))!="try-error")
#     blocknamesav <- c(blocknamesav, "mirna")
#   if(class(try(ncol(mutation)))!="try-error")
#     blocknamesav <- c(blocknamesav, "mutation")
#   if(class(try(ncol(cnv)))!="try-error")
#     blocknamesav <- c(blocknamesav, "cnv")
#   if(class(try(ncol(rna)))!="try-error")
#     blocknamesav <- c(blocknamesav, "rna")  
#   
#   blocknamesav <- blocknamesav[-1]
#   
#   blocknamesavall[[i]] <- blocknamesav
#   
#   eval(parse(text=paste("datasetlist[[i]] <- list(targetvar, ", paste(blocknamesav, collapse=", "), ")", sep="")))
#   eval(parse(text=paste("names(datasetlist[[i]]) <- c(\"targetvar\", \"", paste(blocknamesav, collapse="\", \""), "\")", sep="")))
#   
#   
#   try(rm(clin))
#   try(rm(mirna))
#   try(rm(mutation))
#   try(rm(cnv))
#   try(rm(rna))
#   gc()
#   
#   cat(i, "\n")
#   
# }
# 
# names(datasetlist) <- alldata
# 
# 
# 
# 
# resultstemp <- resultsparamvalues[resultsparamvalues$method=="RandomBlock",]
# 
# resultstemp$method <- resultstemp$cvfoldind <- resultstemp$cvind <- resultstemp$cindex <- resultstemp$mtry <- NULL
# 
# library("plyr")
# resultstempsum <- ddply(resultstemp, .variables=c("dat"), .fun=summarise, cval_clin=mean(cval_clin, na.rm=TRUE), 
#                         cval_mirna=mean(cval_mirna, na.rm=TRUE), cval_mutation=mean(cval_mutation, na.rm=TRUE),
#                         cval_cnv=mean(cval_cnv, na.rm=TRUE),
#                         cval_rna=mean(cval_rna, na.rm=TRUE))
# 
# resultstempsum$maxprob <- apply(resultstempsum[,2:6], 1, max, na.rm=TRUE)
# 
# resultstempsum$cindex.BlockForest <- resultswide$cindex.BlockForest[sapply(resultstempsum$dat, function(x) which(resultswide$dat==x))]
# resultstempsum$cindex.RandomBlock <- resultswide$cindex.RandomBlock[sapply(resultstempsum$dat, function(x) which(resultswide$dat==x))]
# resultstempsum$cindex.RSF <- resultswide$cindex.RSF[sapply(resultstempsum$dat, function(x) which(resultswide$dat==x))]
# 
# resultstempsum$n <- sapply(resultstempsum$dat, function(x) nrow(datasetlist[names(datasetlist)==x][[1]]$targetvar))
# 
# resultstempsum$cindex.MeanRSFRandomBlock <- apply(resultswide[,c("cindex.RSF", "cindex.RandomBlock")], 1, mean, na.rm=TRUE)[sapply(resultstempsum$dat, function(x) which(resultswide$dat==x))]
# resultstempsum$cindex.MeanRSFBlockForest <- apply(resultswide[,c("cindex.RSF", "cindex.BlockForest")], 1, mean, na.rm=TRUE)[sapply(resultstempsum$dat, function(x) which(resultswide$dat==x))]
# 
# save(resultstempsum, file="./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Results/resultstempsum_multiomics.Rda")


# Load the result of the out-commented part above:

load("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Results/resultstempsum_multiomics.Rda")


x <- resultstempsum$maxprob
y <- resultstempsum$cindex.BlockForest - resultstempsum$cindex.RSF

pred <- loess.smooth(x, y)
dataloess <- data.frame(x=pred$x, y=pred$y)

p2 <- ggplot(resultstempsum, aes(x=maxprob, y=cindex.BlockForest - cindex.RSF)) + theme_bw() + geom_point() + 
  geom_hline(yintercept=0, col="grey") +
  geom_line(data=dataloess, aes(x=x, y=y), col="blue") + ggtitle("Influence of 'oneblockimp'") + labs(x="oneblockimp", y="diffC") +
  theme(plot.title = element_text(size=16), axis.title=element_text(size=15), axis.text=element_text(size=11))

p2


x <- resultstempsum$n
y <- resultstempsum$cindex.BlockForest - resultstempsum$cindex.RSF

pred <- loess.smooth(x, y)
dataloess <- data.frame(x=pred$x, y=pred$y)

p1 <- ggplot(resultstempsum, aes(x=n, y=cindex.BlockForest - cindex.RSF)) + theme_bw() + geom_point() + 
  geom_hline(yintercept=0, col="grey") +
  geom_line(data=dataloess, aes(x=x, y=y), col="blue") + ggtitle("Influence of 'n'") + ylab("diffC") +
  theme(plot.title = element_text(size=16), axis.title=element_text(size=15), axis.text=element_text(size=11))

p1


x <- resultstempsum$cindex.MeanRSFBlockForest
y <- resultstempsum$cindex.BlockForest - resultstempsum$cindex.RSF

pred <- loess.smooth(x, y)
dataloess <- data.frame(x=pred$x, y=pred$y)

p3 <- ggplot(resultstempsum, aes(x=cindex.MeanRSFBlockForest, y=cindex.BlockForest - cindex.RSF)) + theme_bw() + geom_point() + 
  geom_hline(yintercept=0, col="grey") +
  geom_line(data=dataloess, aes(x=x, y=y), col="blue") + ggtitle("Influence of 'signal'") + labs(x="signal", y="diffC") +
  theme(plot.title = element_text(size=16), axis.title=element_text(size=15), axis.text=element_text(size=11))

p3

library("grid")
library("gridExtra")
pall <- grid.arrange(p1, p2, p3, ncol = 3)

grid.draw(pall)

ggsave(file="./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Figures/FactorsDifferenceInPerformance.pdf", plot=pall, width=14*0.75, height=5*0.75)
