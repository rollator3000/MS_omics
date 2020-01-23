########################################################

# NOTE: Before the code can be excecuted, the R working directory *MUST* 
# be set to the directory that contains the folder 'Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION':

# Remove '#' from the two lines below and replace 'Z:/here/is/my/path/' by the path
# to the directory that contains 'Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION' ( the path 'Z:/here/is/my/path/'
# *MUST* end with a slash / ):

# basedir <- "Z:/here/is/my/path/"
# setwd(basedir)

########################################################



setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Bladder_Urothelial_Carcinoma_BLCA", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for BLCA data

load("./Raw/CNV_BLCA.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_BLCA.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Brain_Lower_Grade_Glioma_LGG", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for LGG data

load("./Raw/CNV_LGG.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_LGG.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Breast_Invasive_Carcinoma_BRCA", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for BRCA data

load("./Raw/CNV_BRCA.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_BRCA.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Cervical_Squamous_Cell_Carcinoma_and_Endocervical_Adenocarcinoma_CESC", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for CESC data

load("./Raw/CNV_CESC.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_CESC.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Colon_Adenocarcinoma_COAD", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for COAD data

load("./Raw/CNV_COAD.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_COAD.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Esophageal_Carcinoma_ESCA", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for ESCA data

load("./Raw/CNV_ESCA.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_ESCA.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Glioblastoma_Multiforme_GBM", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for GBM data

load("./Raw/CNV_GBM.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_GBM.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Head_and_Neck_Squamous_Cell_Carcinoma_HNSC", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for HNSC data

load("./Raw/CNV_HNSC.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_HNSC.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Kidney_Renal_Clear_Cell_Carcinoma_KIRC", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for KIRC data

load("./Raw/CNV_KIRC.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_KIRC.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Kidney_Renal_Papillary_Cell_Carcinoma_KIRP", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for KIRP data

load("./Raw/CNV_KIRP.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_KIRP.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Liver_Hepatocellular_Carcinoma_LIHC", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for LIHC data

load("./Raw/CNV_LIHC.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_LIHC.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Lung_Adenocarcinoma_LUAD", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for LUAD data

load("./Raw/CNV_LUAD.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_LUAD.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Lung_Squamous_Cell_Carcinoma_LUSC", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for LUSC data

load("./Raw/CNV_LUSC.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_LUSC.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Ovarian_Serous_Cystadenocarcinoma_OV", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for OV data

load("./Raw/CNV_OV.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_OV.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Pancreatic_Adenocarcinoma_PAAD", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for PAAD data

load("./Raw/CNV_PAAD.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_PAAD.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Pheochromocytoma_and_Paraganglioma_PCPG", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for PCPG data

load("./Raw/CNV_PCPG.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_PCPG.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Prostate_Adenocarcinoma_PRAD", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for PRAD data

load("./Raw/CNV_PRAD.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_PRAD.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Rectum_Adenocarcinoma_READ", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for READ data

load("./Raw/CNV_READ.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_READ.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Sarcoma_SARC", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for SARC data

load("./Raw/CNV_SARC.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_SARC.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Skin_Cutaneous_Melanoma_SKCM", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for SKCM data

load("./Raw/CNV_SKCM.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_SKCM.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Stomach_Adenocarcinoma_STAD", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for STAD data

load("./Raw/CNV_STAD.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_STAD.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Testicular_Germ_Cell_Tumors_TGCT", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for TGCT data

load("./Raw/CNV_TGCT.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_TGCT.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Thymoma_THYM", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for THYM data

load("./Raw/CNV_THYM.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_THYM.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Thyroid_Carcinoma_THCA", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for THCA data

load("./Raw/CNV_THCA.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_THCA.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Uterine_Corpus_Endometrial_Carcinoma_UCEC", sep=""))

load(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/loci_all.Rdata", sep=""))
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

###################################################
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  
  return(res)
}
###################################################

# build CNV matrix for UCEC data

load("./Raw/CNV_UCEC.Rdata")
names.cnv <- substr(cnv$Sample, 1, 12)

cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cnv$Sample <- substr(cnv$Sample, 1, 12)
cols <- unique(cnv$Sample)

mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
ind <- which(is.na(mat.cnv))         # set missing values to zero
mat.cnv[ind] <- 0                                 

ind.in <- which(rowSums(abs(mat.cnv))!=0) # remove rows with no detected CNV
mat.cnv <- mat.cnv[ind.in, ]
save(mat.cnv, file="./Processed/CNV_Matrix_UCEC.Rdata")
