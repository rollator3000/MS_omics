########################################################

# NOTE: Before the code can be excecuted, the R working directory *MUST* 
# be set to the directory that contains the folder 'Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION':

# Remove '#' from the two lines below and replace 'Z:/here/is/my/path/' by the path
# to the directory that contains 'Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION' ( the path 'Z:/here/is/my/path/'
# *MUST* end with a slash / ):

# basedir <- "Z:/here/is/my/path/"
# setwd(basedir)

########################################################



setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Bladder_Urothelial_Carcinoma_BLCA/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-BLCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_BLCA.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-BLCA", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_BLCA.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_BLCA.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_BLCA.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-BLCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_BLCA.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-BLCA",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_BLCA.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("BLCA", pipelines = "mutect2")
save(mutations, file="Mutations_BLCA.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Brain_Lower_Grade_Glioma_LGG/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-LGG",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_LGG.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-LGG", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_LGG.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_LGG.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_LGG.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-LGG",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_LGG.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-LGG",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_LGG.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("LGG", pipelines = "mutect2")
save(mutations, file="Mutations_LGG.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Breast_Invasive_Carcinoma_BRCA/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_BRCA.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_BRCA.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_BRCA.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_BRCA.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_BRCA.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_BRCA.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("BRCA", pipelines = "mutect2")
save(mutations, file="Mutations_BRCA.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Cervical_Squamous_Cell_Carcinoma_and_Endocervical_Adenocarcinoma_CESC/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-CESC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_CESC.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-CESC", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_CESC.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_CESC.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_CESC.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-CESC",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_CESC.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-CESC",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_CESC.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("CESC", pipelines = "mutect2")
save(mutations, file="Mutations_CESC.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Colon_Adenocarcinoma_COAD/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-COAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_COAD.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-COAD", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_COAD.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_COAD.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_COAD.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-COAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_COAD.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-COAD",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_COAD.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("COAD", pipelines = "mutect2")
save(mutations, file="Mutations_COAD.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Esophageal_Carcinoma_ESCA/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-ESCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_ESCA.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-ESCA", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_ESCA.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_ESCA.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_ESCA.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-ESCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_ESCA.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-ESCA",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_ESCA.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("ESCA", pipelines = "mutect2")
save(mutations, file="Mutations_ESCA.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Glioblastoma_Multiforme_GBM/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-GBM",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_GBM.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-GBM", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_GBM.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_GBM.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_GBM.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-GBM",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_GBM.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-GBM",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_GBM.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("GBM", pipelines = "mutect2")
save(mutations, file="Mutations_GBM.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Head_and_Neck_Squamous_Cell_Carcinoma_HNSC/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-HNSC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_HNSC.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-HNSC", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_HNSC.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_HNSC.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_HNSC.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-HNSC",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_HNSC.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-HNSC",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_HNSC.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("HNSC", pipelines = "mutect2")
save(mutations, file="Mutations_HNSC.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Kidney_Renal_Clear_Cell_Carcinoma_KIRC/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-KIRC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_KIRC.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-KIRC", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_KIRC.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_KIRC.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_KIRC.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-KIRC",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_KIRC.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-KIRC",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_KIRC.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("KIRC", pipelines = "mutect2")
save(mutations, file="Mutations_KIRC.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Kidney_Renal_Papillary_Cell_Carcinoma_KIRP/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-KIRP",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_KIRP.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-KIRP", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_KIRP.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_KIRP.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_KIRP.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-KIRP",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_KIRP.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-KIRP",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_KIRP.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("KIRP", pipelines = "mutect2")
save(mutations, file="Mutations_KIRP.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Liver_Hepatocellular_Carcinoma_LIHC/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-LIHC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_LIHC.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-LIHC", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_LIHC.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_LIHC.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_LIHC.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-LIHC",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_LIHC.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-LIHC",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_LIHC.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("LIHC", pipelines = "mutect2")
save(mutations, file="Mutations_LIHC.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Lung_Adenocarcinoma_LUAD/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_LUAD.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-LUAD", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_LUAD.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_LUAD.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_LUAD.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_LUAD.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_LUAD.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("LUAD", pipelines = "mutect2")
save(mutations, file="Mutations_LUAD.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Lung_Squamous_Cell_Carcinoma_LUSC/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-LUSC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_LUSC.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-LUSC", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_LUSC.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_LUSC.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_LUSC.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-LUSC",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_LUSC.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-LUSC",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_LUSC.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("LUSC", pipelines = "mutect2")
save(mutations, file="Mutations_LUSC.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Ovarian_Serous_Cystadenocarcinoma_OV/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-OV",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_OV.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-OV", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_OV.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_OV.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_OV.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-OV",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_OV.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-OV",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_OV.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("OV", pipelines = "mutect2")
save(mutations, file="Mutations_OV.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Pancreatic_Adenocarcinoma_PAAD/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-PAAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_PAAD.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-PAAD", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_PAAD.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_PAAD.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_PAAD.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-PAAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_PAAD.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-PAAD",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_PAAD.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("PAAD", pipelines = "mutect2")
save(mutations, file="Mutations_PAAD.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Pheochromocytoma_and_Paraganglioma_PCPG/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-PCPG",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_PCPG.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-PCPG", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_PCPG.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_PCPG.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_PCPG.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-PCPG",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_PCPG.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-PCPG",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_PCPG.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("PCPG", pipelines = "mutect2")
save(mutations, file="Mutations_PCPG.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Prostate_Adenocarcinoma_PRAD/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-PRAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_PRAD.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-PRAD", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_PRAD.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_PRAD.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_PRAD.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-PRAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_PRAD.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-PRAD",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_PRAD.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("PRAD", pipelines = "mutect2")
save(mutations, file="Mutations_PRAD.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Rectum_Adenocarcinoma_READ/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-READ",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_READ.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-READ", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_READ.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_READ.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_READ.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-READ",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_READ.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-READ",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_READ.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("READ", pipelines = "mutect2")
save(mutations, file="Mutations_READ.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Sarcoma_SARC/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-SARC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_SARC.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-SARC", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_SARC.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_SARC.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_SARC.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-SARC",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_SARC.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-SARC",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_SARC.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("SARC", pipelines = "mutect2")
save(mutations, file="Mutations_SARC.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Skin_Cutaneous_Melanoma_SKCM/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-SKCM",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_SKCM.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-SKCM", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_SKCM.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_SKCM.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_SKCM.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-SKCM",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_SKCM.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-SKCM",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_SKCM.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("SKCM", pipelines = "mutect2")
save(mutations, file="Mutations_SKCM.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Stomach_Adenocarcinoma_STAD/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-STAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_STAD.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-STAD", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_STAD.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_STAD.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_STAD.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-STAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_STAD.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-STAD",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_STAD.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("STAD", pipelines = "mutect2")
save(mutations, file="Mutations_STAD.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Testicular_Germ_Cell_Tumors_TGCT/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-TGCT",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_TGCT.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-TGCT", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_TGCT.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_TGCT.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_TGCT.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-TGCT",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_TGCT.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-TGCT",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_TGCT.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("TGCT", pipelines = "mutect2")
save(mutations, file="Mutations_TGCT.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Thymoma_THYM/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-THYM",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_THYM.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-THYM", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_THYM.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_THYM.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_THYM.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-THYM",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_THYM.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-THYM",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_THYM.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("THYM", pipelines = "mutect2")
save(mutations, file="Mutations_THYM.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Thyroid_Carcinoma_THCA/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-THCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_THCA.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-THCA", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_THCA.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_THCA.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_THCA.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-THCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_THCA.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-THCA",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_THCA.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("THCA", pipelines = "mutect2")
save(mutations, file="Mutations_THCA.Rdata")








setwd(paste(basedir, "Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/Uterine_Corpus_Endometrial_Carcinoma_UCEC/Raw", sep=""))

library(TCGAbiolinks)
library(SummarizedExperiment)

#----- RNA expression, measured by RNAseq, raw counts -----#
query <- GDCquery(project = "TCGA-UCEC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk=20)
data <- GDCprepare(query)
rna <- assay(data) 
save(rna, file="RNASeq_Counts_UCEC.Rdata")     # matrix 


#----- clinical data -----#
clin.query <- GDCquery(project = "TCGA-UCEC", data.category = "Clinical")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical <- GDCprepare_clinic(clin.query, clinical.info = "patient") 
save(clinical, file="Clinical_UCEC.Rdata")     # data frame

#----- clinical follow-up data -----#
clinical.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up") 
save(clinical.followup, file="Clinical_Followup_UCEC.Rdata")     # data frame

#----- clinical time to relapse -----#
clinical.relapse <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
save(clinical.relapse, file="Clinical_Relapse_UCEC.Rdata")     # data frame


#----- miRNA expression, measured by RNAseq  -----#
query <- GDCquery(project = "TCGA-UCEC",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query, files.per.chunk=20)
mirna <- GDCprepare(query) 
save(mirna, file="miRNA_UCEC.Rdata")     # data frame


#----- CNV -----#
query <- GDCquery(project = "TCGA-UCEC",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)
save(cnv, file="CNV_UCEC.Rdata")


#----- mutation data -----#
mutations <- GDCquery_Maf("UCEC", pipelines = "mutect2")
save(mutations, file="Mutations_UCEC.Rdata")
