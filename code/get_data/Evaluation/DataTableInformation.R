########################################################

# NOTE: This code requires the preprocessed data sets to
# be contained in the folder "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/ProcessedData".
# These preprocessed data sets are, however, not contained
# there due to their large sizes. They are, however, available
# from the corresponding author upon request.
# 
# Provided that the preprocessed data sets are available, in order
# for the code to be executable, the R working directory *MUST* 
# be set to the directory that contains the folder 'Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION':
#
# Remove '#' from the two lines below and replace 'Z:/here/is/my/path/' by the path
# to the directory that contains 'Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION' ( the path 'Z:/here/is/my/path/'
# *MUST* end with a slash / ):

# basedir <- "Z:/here/is/my/path/"
# setwd(basedir)

########################################################


# Set working directory:

setwd("Z:/Projects/BlockForests/PaperAnalysis")



# Read in data sets and extract necessary information:


datasets <- c("BLCA.Rda", "BRCA.Rda", "CESC.Rda", "COAD.Rda", 
              "ESCA.Rda", "GBM.Rda", "HNSC.Rda", "KIRC.Rda", 
              "KIRP.Rda", "LGG.Rda", "LIHC.Rda", "LUAD.Rda", 
              "LUSC.Rda", "OV.Rda", "PAAD.Rda", "PRAD.Rda", 
              "READ.Rda", "SARC.Rda", "SKCM.Rda", "STAD.Rda", 
              "UCEC.Rda")

ns <- 0
percobserved <- ""
clinpres <- cnvpres <- mirnapres <- mutationpres <- rnapres <- rep(FALSE, length(datasets))
nclin <- ncnv <- nmirna <- nmutation <- nrna <- rep(NA, length(datasets))

omicssttring <- ""


for(i in seq(along=datasets)) {
  
  load(paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/ProcessedData/", datasets[i], sep=""))

  ns[i] <- nrow(targetvar)
  percobserved[i] <- paste(100*round(mean(targetvar$status), 2), "%")
  
  omicssttringtemp <- ""
  count <- 1
  
  if(exists("clin")) {
    clinpres[i] <- TRUE
    nclin[i] <- ncol(clin)
    omicssttringtemp[count] <- paste("clinical (p = ", nclin[i], ")", sep="")
    rm(clin)
    count <- count + 1
  }
  
  if(exists("cnv")) {
    cnvpres[i] <- TRUE
    ncnv[i] <- ncol(cnv)
    omicssttringtemp[count] <- paste("CNV (p = ", ncnv[i], ")", sep="")
    rm(cnv)
    count <- count + 1
  }
  
  if(exists("mirna")) {
    mirnapres[i] <- TRUE
    nmirna[i] <- ncol(mirna)
    omicssttringtemp[count] <- paste("miRNA (p = ", nmirna[i], ")", sep="")
    rm(mirna)
    count <- count + 1
  }
  
  if(exists("mutation")) {
    mutationpres[i] <- TRUE
    nmutation[i] <- ncol(mutation)
    omicssttringtemp[count] <- paste("mutation (p = ", nmutation[i], ")", sep="")
    rm(mutation)
    count <- count + 1
  }
  
  if(exists("rna")) {
    rnapres[i] <- TRUE
    nrna[i] <- ncol(rna)
    omicssttringtemp[count] <- paste("rna (p = ", nrna[i], ")", sep="")
    rm(rna)
    count <- count + 1
  }
  
  omicssttring[i] <- paste(omicssttringtemp, collapse="\", \"")
  
  cat(i, "\n")
  
}



# Numbers of covariates:


numbercovariates <- matrix(nrow=length(datasets), ncol=5)
colnames(numbercovariates) <- c("clin", "mirna", "mutation", "cnv", "rna")

rownames(numbercovariates) <- datasets

numbercovariates[,"clin"] <- nclin
numbercovariates[,"mirna"] <- nmirna
numbercovariates[,"mutation"] <- nmutation
numbercovariates[,"cnv"] <- ncnv
numbercovariates[,"rna"] <- nrna

numbercovariates

save(numbercovariates, file="./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/numbercovariates.Rda")




# Mean numbers of covariates per block:


meannums <- c(mean(nclin, na.rm=TRUE), mean(ncnv, na.rm=TRUE), mean(nmirna, na.rm=TRUE), 
              mean(nmutation, na.rm=TRUE), mean(nrna, na.rm=TRUE))
round(meannums, 1)




# Extract cancer names from folder names containing the different
# corresponding data sets:


foldernames <- list.dirs("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing", full.names=FALSE, recursive=FALSE)

datasets <- c("BLCA.Rda", "BRCA.Rda", "CESC.Rda", "COAD.Rda", 
              "ESCA.Rda", "GBM.Rda", "HNSC.Rda", "KIRC.Rda", 
              "KIRP.Rda", "LGG.Rda", "LIHC.Rda", "LUAD.Rda", 
              "LUSC.Rda", "OV.Rda", "PAAD.Rda", "PRAD.Rda", 
              "READ.Rda", "SARC.Rda", "SKCM.Rda", "STAD.Rda", 
              "UCEC.Rda")

datasets2 <- gsub(".Rda", "", datasets)


datanames <- foldernames[sapply(datasets2, function(x) grep(x, foldernames))]
datanames

datanames <- sapply(datanames, function(x) {
  y <- strsplit(x, split="_")[[1]]
  y <- y[-length(y)]
  paste(y, collapse=" ")
})
names(datanames) <- NULL

datanames



# Other information included in the tables in the main paper and
# in Additional file 1:


datasets2
ns
percobserved
