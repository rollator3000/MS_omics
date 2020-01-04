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

library(limma); library(edgeR)

# RNA matrix for BLCA data

load("./Raw/RNASeq_Counts_BLCA.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_BLCA.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Brain_Lower_Grade_Glioma_LGG", sep=""))

library(limma); library(edgeR)

# RNA matrix for LGG data

load("./Raw/RNASeq_Counts_LGG.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_LGG.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Breast_Invasive_Carcinoma_BRCA", sep=""))

library(limma); library(edgeR)

# RNA matrix for BRCA data

load("./Raw/RNASeq_Counts_BRCA.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_BRCA.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Cervical_Squamous_Cell_Carcinoma_and_Endocervical_Adenocarcinoma_CESC", sep=""))

library(limma); library(edgeR)

# RNA matrix for CESC data

load("./Raw/RNASeq_Counts_CESC.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_CESC.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Colon_Adenocarcinoma_COAD", sep=""))

library(limma); library(edgeR)

# RNA matrix for COAD data

load("./Raw/RNASeq_Counts_COAD.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_COAD.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Esophageal_Carcinoma_ESCA", sep=""))

library(limma); library(edgeR)

# RNA matrix for ESCA data

load("./Raw/RNASeq_Counts_ESCA.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_ESCA.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Glioblastoma_Multiforme_GBM", sep=""))

library(limma); library(edgeR)

# RNA matrix for GBM data

load("./Raw/RNASeq_Counts_GBM.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_GBM.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Head_and_Neck_Squamous_Cell_Carcinoma_HNSC", sep=""))

library(limma); library(edgeR)

# RNA matrix for HNSC data

load("./Raw/RNASeq_Counts_HNSC.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_HNSC.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Kidney_Renal_Clear_Cell_Carcinoma_KIRC", sep=""))

library(limma); library(edgeR)

# RNA matrix for KIRC data

load("./Raw/RNASeq_Counts_KIRC.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_KIRC.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Kidney_Renal_Papillary_Cell_Carcinoma_KIRP", sep=""))

library(limma); library(edgeR)

# RNA matrix for KIRP data

load("./Raw/RNASeq_Counts_KIRP.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_KIRP.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Liver_Hepatocellular_Carcinoma_LIHC", sep=""))

library(limma); library(edgeR)

# RNA matrix for LIHC data

load("./Raw/RNASeq_Counts_LIHC.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_LIHC.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Lung_Adenocarcinoma_LUAD", sep=""))

library(limma); library(edgeR)

# RNA matrix for LUAD data

load("./Raw/RNASeq_Counts_LUAD.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_LUAD.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Lung_Squamous_Cell_Carcinoma_LUSC", sep=""))

library(limma); library(edgeR)

# RNA matrix for LUSC data

load("./Raw/RNASeq_Counts_LUSC.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_LUSC.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Ovarian_Serous_Cystadenocarcinoma_OV", sep=""))

library(limma); library(edgeR)

# RNA matrix for OV data

load("./Raw/RNASeq_Counts_OV.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_OV.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Pancreatic_Adenocarcinoma_PAAD", sep=""))

library(limma); library(edgeR)

# RNA matrix for PAAD data

load("./Raw/RNASeq_Counts_PAAD.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_PAAD.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Pheochromocytoma_and_Paraganglioma_PCPG", sep=""))

library(limma); library(edgeR)

# RNA matrix for PCPG data

load("./Raw/RNASeq_Counts_PCPG.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_PCPG.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Prostate_Adenocarcinoma_PRAD", sep=""))

library(limma); library(edgeR)

# RNA matrix for PRAD data

load("./Raw/RNASeq_Counts_PRAD.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_PRAD.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Rectum_Adenocarcinoma_READ", sep=""))

library(limma); library(edgeR)

# RNA matrix for READ data

load("./Raw/RNASeq_Counts_READ.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_READ.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Sarcoma_SARC", sep=""))

library(limma); library(edgeR)

# RNA matrix for SARC data

load("./Raw/RNASeq_Counts_SARC.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_SARC.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Skin_Cutaneous_Melanoma_SKCM", sep=""))

library(limma); library(edgeR)

# RNA matrix for SKCM data

load("./Raw/RNASeq_Counts_SKCM.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_SKCM.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Stomach_Adenocarcinoma_STAD", sep=""))

library(limma); library(edgeR)

# RNA matrix for STAD data

load("./Raw/RNASeq_Counts_STAD.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_STAD.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Testicular_Germ_Cell_Tumors_TGCT", sep=""))

library(limma); library(edgeR)

# RNA matrix for TGCT data

load("./Raw/RNASeq_Counts_TGCT.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_TGCT.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Thymoma_THYM", sep=""))

library(limma); library(edgeR)

# RNA matrix for THYM data

load("./Raw/RNASeq_Counts_THYM.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_THYM.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Thyroid_Carcinoma_THCA", sep=""))

library(limma); library(edgeR)

# RNA matrix for THCA data

load("./Raw/RNASeq_Counts_THCA.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_THCA.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Uterine_Corpus_Endometrial_Carcinoma_UCEC", sep=""))

library(limma); library(edgeR)

# RNA matrix for UCEC data

load("./Raw/RNASeq_Counts_UCEC.Rdata")
ind.in <- which(!duplicated(substr(colnames(rna), 1, 12)))  # remove duplicate samples
mat.rna <- rna[ ,ind.in]
dge <- DGEList(counts=mat.rna)
cc <- cpm(dge)
keep <- which(rowSums(cc>=0.5)>=0.05*ncol(cc))             # remove underexpressed RNAs 
rna <- mat.rna[keep, ]                                     # for each gene, at least 5% of patients should have values over 0.5
  
dge <- DGEList(counts=rna)
dge <- calcNormFactors(dge)
v <- voom(dge)
rna.matrix <- v$E
colnames(rna.matrix) <- substr(colnames(rna.matrix), 1, 12)
  
save(rna.matrix, file="./Processed/RNA_Matrix_UCEC.Rdata")
