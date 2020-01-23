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

library("limma")

# miRNA matrix for BLCA data

load("./Raw/miRNA_BLCA.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_BLCA.Rdata")  








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Brain_Lower_Grade_Glioma_LGG", sep=""))

library("limma")

# miRNA matrix for LGG data

load("./Raw/miRNA_LGG.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_LGG.Rdata")  








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Breast_Invasive_Carcinoma_BRCA", sep=""))

library("limma")

# miRNA matrix for BRCA data

load("./Raw/miRNA_BRCA.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_BRCA.Rdata")  








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Cervical_Squamous_Cell_Carcinoma_and_Endocervical_Adenocarcinoma_CESC", sep=""))

library("limma")

# miRNA matrix for CESC data

load("./Raw/miRNA_CESC.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_CESC.Rdata")  








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Colon_Adenocarcinoma_COAD", sep=""))

library("limma")

# miRNA matrix for COAD data

load("./Raw/miRNA_COAD.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_COAD.Rdata")  








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Esophageal_Carcinoma_ESCA", sep=""))

library("limma")

# miRNA matrix for ESCA data

load("./Raw/miRNA_ESCA.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_ESCA.Rdata")  








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Glioblastoma_Multiforme_GBM", sep=""))

library("limma")

# miRNA matrix for GBM data

load("./Raw/miRNA_GBM.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_GBM.Rdata")  








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Head_and_Neck_Squamous_Cell_Carcinoma_HNSC", sep=""))

library("limma")

# miRNA matrix for HNSC data

load("./Raw/miRNA_HNSC.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_HNSC.Rdata")  








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Kidney_Renal_Clear_Cell_Carcinoma_KIRC", sep=""))

library("limma")

# miRNA matrix for KIRC data

load("./Raw/miRNA_KIRC.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_KIRC.Rdata")  








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Kidney_Renal_Papillary_Cell_Carcinoma_KIRP", sep=""))

library("limma")

# miRNA matrix for KIRP data

load("./Raw/miRNA_KIRP.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_KIRP.Rdata")  








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Liver_Hepatocellular_Carcinoma_LIHC", sep=""))

library("limma")

# miRNA matrix for LIHC data

load("./Raw/miRNA_LIHC.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_LIHC.Rdata")  








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Lung_Adenocarcinoma_LUAD", sep=""))

library("limma")

# miRNA matrix for LUAD data

load("./Raw/miRNA_LUAD.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_LUAD.Rdata")  








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Lung_Squamous_Cell_Carcinoma_LUSC", sep=""))

library("limma")

# miRNA matrix for LUSC data

load("./Raw/miRNA_LUSC.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_LUSC.Rdata")  








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Ovarian_Serous_Cystadenocarcinoma_OV", sep=""))

library("limma")

# miRNA matrix for OV data

load("./Raw/miRNA_OV.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_OV.Rdata")  








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Pancreatic_Adenocarcinoma_PAAD", sep=""))

library("limma")

# miRNA matrix for PAAD data

load("./Raw/miRNA_PAAD.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_PAAD.Rdata")  








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Pheochromocytoma_and_Paraganglioma_PCPG", sep=""))

library("limma")

# miRNA matrix for PCPG data

load("./Raw/miRNA_PCPG.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_PCPG.Rdata")  








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Prostate_Adenocarcinoma_PRAD", sep=""))

library("limma")

# miRNA matrix for PRAD data

load("./Raw/miRNA_PRAD.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_PRAD.Rdata")  








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Rectum_Adenocarcinoma_READ", sep=""))

library("limma")

# miRNA matrix for READ data

load("./Raw/miRNA_READ.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_READ.Rdata")  








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Sarcoma_SARC", sep=""))

library("limma")

# miRNA matrix for SARC data

load("./Raw/miRNA_SARC.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_SARC.Rdata")  








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Skin_Cutaneous_Melanoma_SKCM", sep=""))

library("limma")

# miRNA matrix for SKCM data

load("./Raw/miRNA_SKCM.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_SKCM.Rdata")  








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Stomach_Adenocarcinoma_STAD", sep=""))

library("limma")

# miRNA matrix for STAD data

load("./Raw/miRNA_STAD.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_STAD.Rdata")  








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Testicular_Germ_Cell_Tumors_TGCT", sep=""))

library("limma")

# miRNA matrix for TGCT data

load("./Raw/miRNA_TGCT.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_TGCT.Rdata")  








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Thymoma_THYM", sep=""))

library("limma")

# miRNA matrix for THYM data

load("./Raw/miRNA_THYM.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_THYM.Rdata")  








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Thyroid_Carcinoma_THCA", sep=""))

library("limma")

# miRNA matrix for THCA data

load("./Raw/miRNA_THCA.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_THCA.Rdata")  








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Uterine_Corpus_Endometrial_Carcinoma_UCEC", sep=""))

library("limma")

# miRNA matrix for UCEC data

load("./Raw/miRNA_UCEC.Rdata")
ind <- grep("reads_per_million_miRNA_mapped", colnames(mirna))
mirna.matrix <- as.matrix(mirna[ ,ind])
rownames(mirna.matrix) <- mirna$miRNA_ID
cols <- colnames(mirna.matrix)
cols <- gsub("reads_per_million_miRNA_mapped_", "", cols)
cols <- substr(cols, 1, 12)
colnames(mirna.matrix) <- cols
ind.in <- which(!duplicated(cols))                                       # remove duplicate samples
mirna.matrix <- mirna.matrix[ ,ind.in]

keep <- which(rowSums(mirna.matrix>=0.5)>=0.05*ncol(mirna.matrix))       # remove underexpressed miRNAs 
mirna.matrix <- mirna.matrix[keep, ]

v <- voom(mirna.matrix)
mirna.matrix <- v$E

save(mirna.matrix, file="./Processed/miRNA_Matrix_UCEC.Rdata")  








