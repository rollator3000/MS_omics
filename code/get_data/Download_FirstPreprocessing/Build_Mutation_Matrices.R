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

# Mutation matrix for BLCA data

load("./Raw/Mutations_BLCA.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_BLCA.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Brain_Lower_Grade_Glioma_LGG", sep=""))

# Mutation matrix for LGG data

load("./Raw/Mutations_LGG.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_LGG.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Breast_Invasive_Carcinoma_BRCA", sep=""))

# Mutation matrix for BRCA data

load("./Raw/Mutations_BRCA.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_BRCA.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Cervical_Squamous_Cell_Carcinoma_and_Endocervical_Adenocarcinoma_CESC", sep=""))

# Mutation matrix for CESC data

load("./Raw/Mutations_CESC.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_CESC.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Colon_Adenocarcinoma_COAD", sep=""))

# Mutation matrix for COAD data

load("./Raw/Mutations_COAD.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_COAD.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Esophageal_Carcinoma_ESCA", sep=""))

# Mutation matrix for ESCA data

load("./Raw/Mutations_ESCA.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_ESCA.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Glioblastoma_Multiforme_GBM", sep=""))

# Mutation matrix for GBM data

load("./Raw/Mutations_GBM.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_GBM.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Head_and_Neck_Squamous_Cell_Carcinoma_HNSC", sep=""))

# Mutation matrix for HNSC data

load("./Raw/Mutations_HNSC.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_HNSC.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Kidney_Renal_Clear_Cell_Carcinoma_KIRC", sep=""))

# Mutation matrix for KIRC data

load("./Raw/Mutations_KIRC.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_KIRC.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Kidney_Renal_Papillary_Cell_Carcinoma_KIRP", sep=""))

# Mutation matrix for KIRP data

load("./Raw/Mutations_KIRP.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_KIRP.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Liver_Hepatocellular_Carcinoma_LIHC", sep=""))

# Mutation matrix for LIHC data

load("./Raw/Mutations_LIHC.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_LIHC.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Lung_Adenocarcinoma_LUAD", sep=""))

# Mutation matrix for LUAD data

load("./Raw/Mutations_LUAD.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_LUAD.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Lung_Squamous_Cell_Carcinoma_LUSC", sep=""))

# Mutation matrix for LUSC data

load("./Raw/Mutations_LUSC.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_LUSC.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Ovarian_Serous_Cystadenocarcinoma_OV", sep=""))

# Mutation matrix for OV data

load("./Raw/Mutations_OV.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_OV.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Pancreatic_Adenocarcinoma_PAAD", sep=""))

# Mutation matrix for PAAD data

load("./Raw/Mutations_PAAD.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_PAAD.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Pheochromocytoma_and_Paraganglioma_PCPG", sep=""))

# Mutation matrix for PCPG data

load("./Raw/Mutations_PCPG.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_PCPG.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Prostate_Adenocarcinoma_PRAD", sep=""))

# Mutation matrix for PRAD data

load("./Raw/Mutations_PRAD.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_PRAD.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Rectum_Adenocarcinoma_READ", sep=""))

# Mutation matrix for READ data

load("./Raw/Mutations_READ.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_READ.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Sarcoma_SARC", sep=""))

# Mutation matrix for SARC data

load("./Raw/Mutations_SARC.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_SARC.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Skin_Cutaneous_Melanoma_SKCM", sep=""))

# Mutation matrix for SKCM data

load("./Raw/Mutations_SKCM.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_SKCM.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Stomach_Adenocarcinoma_STAD", sep=""))

# Mutation matrix for STAD data

load("./Raw/Mutations_STAD.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_STAD.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Testicular_Germ_Cell_Tumors_TGCT", sep=""))

# Mutation matrix for TGCT data

load("./Raw/Mutations_TGCT.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_TGCT.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Thymoma_THYM", sep=""))

# Mutation matrix for THYM data

load("./Raw/Mutations_THYM.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_THYM.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Thyroid_Carcinoma_THCA", sep=""))

# Mutation matrix for THCA data

load("./Raw/Mutations_THCA.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_THCA.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Uterine_Corpus_Endometrial_Carcinoma_UCEC", sep=""))

# Mutation matrix for UCEC data

load("./Raw/Mutations_UCEC.Rdata")

genes <- unique(mutations$Hugo_Symbol)
sample <- unique(substr(mutations$Tumor_Sample_Barcode, 1, 12))
mut.sample <- substr(mutations$Tumor_Sample_Barcode, 1, 12)
  
mutation.matrix <- matrix(0, nr=length(genes), ncol=length(sample))
rownames(mutation.matrix) <- genes
colnames(mutation.matrix) <- sample
  
for(j in 1:length(genes)) {
  ind <- which(mutations$Hugo_Symbol==genes[j])
  mutation.matrix[j, mut.sample[ind]] <- 1
}

save(mutation.matrix, file="./Processed/Mutation_Matrix_UCEC.Rdata")
