########################################################

# NOTE: Before the code can be excecuted, the R working directory *MUST* 
# be set to the directory that contains the folder 'Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION':

# Remove '#' from the two lines below and replace 'Z:/here/is/my/path/' by the path
# to the directory that contains 'Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION' ( the path 'Z:/here/is/my/path/'
# *MUST* end with a slash / ):

# basedir <- "Z:/here/is/my/path/"
# setwd(basedir)

########################################################



cancertypeshort <- "ESCA"
cancertypelong <- "Esophageal_Carcinoma_ESCA"


# Load clinical data and follow-up data:

load(paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/", cancertypelong, "/Raw/Clinical_", cancertypeshort, ".Rdata", sep=""))
load(paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/", cancertypelong, "/Raw/Clinical_Followup_", cancertypeshort, ".Rdata", sep=""))



# Check whether each patient is only represented once in the clinical data:

length(clinical$bcr_patient_barcode)==length(unique(clinical$bcr_patient_barcode))

# --> OK.



# In the follow-up data there are multiple measurements per patient and we want
# to use only the last follow-up time points.

# Make a vector 'inclind' that contains the indices of measurements in 
# 'clinical.followup' to use:

# We loop through all patients:
#
# - If a patient features only one measurement in 'clinical.followup,
#   use that measurement.
#
# - If a patient features more than one measurement: If there is a death
#   time for this patient, use the sample with the death time in
#   'clinical.followup' and if there is no death time, use the sample
#   with the longest follow-up time.

set.seed(1234)

library("nnet")

inclind <- c()

allpatients <- unique(clinical.followup$bcr_patient_barcode)
for(i in 1:length(allpatients)) {
  
  if(sum(clinical.followup$bcr_patient_barcode==allpatients[i])==1)
    inclind[i] <- which(clinical.followup$bcr_patient_barcode==allpatients[i])
  else {
    
    vectemp <- clinical.followup$days_to_death[clinical.followup$bcr_patient_barcode==allpatients[i]]
    if(!all(is.na(vectemp))) {
      vectemp[is.na(vectemp)] <- min(vectemp, na.rm=TRUE)
      
      inclind[i] <- which(clinical.followup$bcr_patient_barcode==allpatients[i])[which.is.max(vectemp)]
    }
    else {
      vectemp <- clinical.followup$days_to_last_followup[clinical.followup$bcr_patient_barcode==allpatients[i]]
      vectemp[is.na(vectemp)] <- min(vectemp, na.rm=TRUE)
      
      inclind[i] <- which(clinical.followup$bcr_patient_barcode==allpatients[i])[which.is.max(vectemp)]
    }   
    
  }
  
}


# Subset 'clinical.followup' correspondingly:

clinical.followup2 <- clinical.followup[inclind,]

dim(clinical.followup)
dim(clinical.followup2)


# Subset 'clinical' to feature only patients present
# in 'clinical.followup2':

inclind <- clinical$bcr_patient_barcode %in% clinical.followup2$bcr_patient_barcode
sum(inclind)

clinical2 <- clinical[inclind,]



# Check whether the patients in clinical2 and clinical.followup2 have 
# the same ordering and if not rearrange ordering of clinical.followup2:

all(clinical2$bcr_patient_barcode==clinical.followup2$bcr_patient_barcode)

all(clinical2$bcr_patient_barcode %in% clinical.followup2$bcr_patient_barcode)
all(clinical.followup2$bcr_patient_barcode %in% clinical2$bcr_patient_barcode)

clinical.followup2 <- clinical.followup2[as.numeric(factor(clinical2$bcr_patient_barcode, levels=clinical.followup2$bcr_patient_barcode)),]

all(clinical2$bcr_patient_barcode==clinical.followup2$bcr_patient_barcode)


# For some patients there exists neither a follow-up time
# nor a follow-up time:

table(data.frame(x1=is.na(clinical.followup2$days_to_last_followup),
                 x2=is.na(clinical.followup2$days_to_death)))

# --> Exclude these patients.

inclind <- !(is.na(clinical.followup2$days_to_last_followup) & is.na(clinical.followup2$days_to_death))
clinical2 <- clinical2[inclind,]
clinical.followup2 <- clinical.followup2[inclind,]

table(data.frame(x1=is.na(clinical.followup2$days_to_last_followup),
                 x2=is.na(clinical.followup2$days_to_death)))

# --> Some patients have both a follow-up time and a death time.
# Check the vital status of these patients:

table(clinical.followup2$vital_status[!is.na(clinical.followup2$days_to_last_followup) & !is.na(clinical.followup2$days_to_death)])

# --> These patients are all dead.
# --> Simple delete the follow-up times of these patients:

clinical.followup2$days_to_last_followup[!is.na(clinical.followup2$days_to_last_followup) & !is.na(clinical.followup2$days_to_death)] <- NA



table(data.frame(x1=is.na(clinical.followup2$days_to_last_followup),
                 x2=is.na(clinical.followup2$days_to_death)))

# --> OK.



# if(!is.null(clinical2$tobacco_smoking_history)) {
#   
#   # Smoking status will be used as a clinical covariate:
#   
#   # Categories:
#   
#   # Lifelong Non-smoker (less than 100 cigarettes smoked in Lifetime) = 1
#   # Current smoker (includes daily smokers and non-daily smokers or occasional smokers) = 2
#   # Current reformed smoker for > 15 years (greater than 15 years) = 3
#   # Current reformed smoker for <=15 years (less than or equal to 15 years) = 4
#   # Current reformed smoker, duration not specified = 5
#   
#   table(clinical2$tobacco_smoking_history)
#   
#   # --> Unfortunately, categoires 1 and 5 are too small.
#   # --> Exclude these (as in Zhao et al. (2013)).
#   
#   
#   inclind <- !(clinical2$tobacco_smoking_history %in% c(1, 5))
#   
#   clinical2 <- clinical2[inclind,]
#   clinical.followup2 <- clinical.followup2[inclind,]
#   
#   table(clinical2$tobacco_smoking_history)
#   table(clinical2$tobacco_smoking_history, useNA="ifany")
#   
#   # --> Also include NAs as some of these patients might also
#   # be in categories 1 and 5.
#   
#   inclind <- !is.na(clinical2$tobacco_smoking_history)
#   
#   clinical2 <- clinical2[inclind,]
#   clinical.followup2 <- clinical.followup2[inclind,]
#   
#   table(clinical2$tobacco_smoking_history, useNA="ifany")
#   
# }



# Make clinical variables data frame:

clin <- data.frame(status=as.numeric(!is.na(clinical.followup2$days_to_death)))

clin$time <- ifelse(clin$status==1, clinical.followup2$days_to_death, clinical.followup2$days_to_last_followup)


# Variable 'age':

# clinical2$age_at_initial_pathologic_diagnosis
# clin$age <- as.numeric(as.character(clinical2$age_at_initial_pathologic_diagnosis))


# Variable 'race':

table(clinical2$race_list)

clin$race <- 0
clin$race <- as.numeric(clinical2$race_list=="WHITE")
clin$race[clinical2$race_list==""] <- NA


# Variable 'gender':

table(clinical2$gender, useNA="ifany")

clin$gender <- 0
clin$gender <- as.numeric(clinical2$gender=="MALE")


# if(!is.null(clinical2$tobacco_smoking_history)) {
#   
#   # Variable 'smoking':
#   
#   # Lifelong Non-smoker (less than 100 cigarettes smoked in Lifetime) = 1
#   # Current smoker (includes daily smokers and non-daily smokers or occasional smokers) = 2
#   # Current reformed smoker for > 15 years (greater than 15 years) = 3
#   # Current reformed smoker for <= 15 years (less than or equal to 15 years) = 4
#   # Current reformed smoker, duration not specified = 5
#   
#   table(clinical2$tobacco_smoking_history)
#   
#   clin$smoking_longer15 <- NA
#   clin$smoking_longer15[clinical2$tobacco_smoking_history==3] <- 1
#   clin$smoking_longer15[clinical2$tobacco_smoking_history %in% c(2,4)] <- 0
#   table(clin$smoking_longer15, useNA="ifany")
#   
#   clin$smoking_shorter15 <- NA
#   clin$smoking_shorter15[clinical2$tobacco_smoking_history==4] <- 1
#   clin$smoking_shorter15[clinical2$tobacco_smoking_history %in% c(2,3)] <- 0
#   table(clin$smoking_shorter15, useNA="ifany")
#   
# }


# # Variable 'ER status':
#
# table(clinical2$breast_carcinoma_estrogen_receptor_status)
# 
# clin$er_status <- 0
# clin$er_status[clinical2$breast_carcinoma_estrogen_receptor_status=="Positive"] <- 1
# clin$er_status[clinical2$breast_carcinoma_estrogen_receptor_status==""] <- NA


# # Variable 'PR status':
# 
# clin$pr_status <- 0
# clin$pr_status[clinical2$breast_carcinoma_progesterone_receptor_status=="Positive"] <- 1
# clin$pr_status[clinical2$breast_carcinoma_progesterone_receptor_status==""] <- NA


# # Variable 'HER2 final status', dummy-coded with reference category 'Negative':
# 
# table(clinical2$lab_proc_her2_neu_immunohistochemistry_receptor_status)
#
# clin$her2_positive <- NA
# clin$her2_positive[clinical2$lab_proc_her2_neu_immunohistochemistry_receptor_status=="Positive"] <- 1
# clin$her2_positive[clinical2$lab_proc_her2_neu_immunohistochemistry_receptor_status %in% c("Equivocal", "Negative")] <- 0
# table(clin$her2_positive, useNA="ifany")
# 
# clin$her2_equivocal <- NA
# clin$her2_equivocal[clinical2$lab_proc_her2_neu_immunohistochemistry_receptor_status=="Equivocal"] <- 1
# clin$her2_equivocal[clinical2$lab_proc_her2_neu_immunohistochemistry_receptor_status %in% c("Positive", "Negative")] <- 0
# table(clin$her2_equivocal, useNA="ifany")
# 
# # --> Include HER2 final status.



# Variable 'T':
# Tumor stage code (T1 versusT_other)

table(clinical2$stage_event_tnm_categories)

clin$T_positive <- NA
clin$T_positive <- 0
tempind <- grep("T1", clinical2$stage_event_tnm_categories)
clin$T_positive[tempind] <- 1

table(clin$T_positive)

# --> Too few positive observations.
# --> Exclude variable again:

# clin$T_positive <- NULL



# Variable 'N':
# Lymph node stage (positive versus negative)

table(clinical2$stage_event_tnm_categories)

clin$N_positive <- NA
tempind <- c(grep("N1", clinical2$stage_event_tnm_categories),
             grep("N2", clinical2$stage_event_tnm_categories),
             grep("N3", clinical2$stage_event_tnm_categories))
clin$N_positive[tempind] <- 1
tempind <- grep("N0", clinical2$stage_event_tnm_categories)
clin$N_positive[tempind] <- 0

table(clin$N_positive)



# Variable 'M':
# Metastasis stage code (positive versus negative)

table(clinical2$stage_event_tnm_categories)

clin$M_positive <- NA
clin$M_positive <- 0
tempind <- grep("M1", clinical2$stage_event_tnm_categories)
clin$M_positive[tempind] <- 1

table(clin$M_positive)

# # --> Too few positive observations.
# # --> Exclude variable again:
# 
clin$M_positive <- NULL




# Include patient ID:

clin$pat_id <- clinical2$bcr_patient_barcode

head(clin)



# Exclude patients with more than four missing values
# in the clinical covariates:

table(apply(clin, 1, function(x) sum(is.na(x))))
nrow(clin)
clin <- clin[apply(clin, 1, function(x) sum(is.na(x))) < 5,]
nrow(clin)




# Make a new data.frame for the clinical covariates only that
# excludes the survival information and the variables for observation
# and sample ID:

clinsafe <- clin

clin <- clinsafe[,!(names(clinsafe) %in% c("status", "time", "pat_id"))]


head(clin)



# Impute the clinical covariates using kNN imputation:

apply(clin, 2, function(x) sum(is.na(x)))
apply(clin, 2, function(x) mean(is.na(x)))


isfactor <- apply(clin, 2, function(x) all(x[!is.na(x)] %in% c(0,1)))
factorvars <- names(clin)[isfactor]

for(i in seq(along=factorvars))
  eval(parse(text=paste("clin$", factorvars[i], " <- factor(clin$", factorvars[i], ", levels=c(1,0))", sep="")))

library("DMwR")
clin <- knnImputation(clin, k = 10)

for(i in seq(along=factorvars))
  eval(parse(text=paste("clin$", factorvars[i], " <- as.numeric(as.character(clin$", factorvars[i], "))", sep="")))

head(clin)

# As rownames use the patient ID:

rownames(clin) <- clinsafe$pat_id



# Make own data.frame for target variable:

targetvar <- clinsafe[, names(clinsafe) %in% c("time", "status")]

rownames(targetvar) <- clinsafe$pat_id




# Prepare omics data types:

# miRNA data:

obj <- ls()

load(paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/", cancertypelong, "/Processed/miRNA_Matrix_", cancertypeshort, ".Rdata", sep=""))

setdiff(ls(), obj)

dim(mirna.matrix)

mirna.matrix[1:4, 1:2]

class(mirna.matrix)

any(is.na(mirna.matrix))

mirna <- t(mirna.matrix)



# mutation data:

obj <- ls()

load(paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/", cancertypelong, "/Processed/Mutation_Matrix_", cancertypeshort, ".Rdata", sep=""))

setdiff(ls(), obj)

dim(mutation.matrix)

mutation.matrix[1:4, 1:2]

class(mutation.matrix)

any(is.na(mutation.matrix))

mutation <- t(mutation.matrix)



# CNV data:

obj <- ls()

load(paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/", cancertypelong, "/Processed/CNV_Matrix_", cancertypeshort, ".Rdata", sep=""))

setdiff(ls(), obj)

dim(mat.cnv)

mat.cnv[1:4, 1:2]

class(mat.cnv)

any(is.na(mat.cnv))

cnv <- t(mat.cnv)



# RNA data:

obj <- ls()

load(paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/Download_FirstPreprocessing/", cancertypelong, "/Processed/RNA_Matrix_", cancertypeshort, ".Rdata", sep=""))

setdiff(ls(), obj)

dim(rna.matrix)

rna.matrix[1:4, 1:2]

class(rna.matrix)

any(is.na(rna.matrix))

rna <- t(rna.matrix)




# Determine set of patients that have measurements of all data types:

c(nrow(clin), nrow(mirna), nrow(mutation), nrow(cnv), nrow(rna))

comids <- intersect(intersect(intersect(intersect(rownames(clin), rownames(mirna)), rownames(mutation)), rownames(cnv)), rownames(rna))
length(comids)




# Subset the different data types to feature all the same patients:

clin <- clin[rownames(clin) %in% comids,]
targetvar <- targetvar[rownames(targetvar) %in% comids,]

mirna <- mirna[rownames(mirna) %in% comids,]
mutation <- mutation[rownames(mutation) %in% comids,]
cnv <- cnv[rownames(cnv) %in% comids,]
rna <- rna[rownames(rna) %in% comids,]



# Re-order the observations in the different data types in order for them to have
# the same ordering:

all(rownames(clin)==rownames(targetvar))
all(rownames(clin)==rownames(mirna))
all(rownames(clin)==rownames(mutation))
all(rownames(clin)==rownames(cnv))
all(rownames(clin)==rownames(rna))

mirna <- mirna[as.numeric(factor(rownames(clin), levels=rownames(mirna))),]
mutation <- mutation[as.numeric(factor(rownames(clin), levels=rownames(mutation))),]
cnv <- cnv[as.numeric(factor(rownames(clin), levels=rownames(cnv))),]
rna <- rna[as.numeric(factor(rownames(clin), levels=rownames(rna))),]

all(rownames(clin)==rownames(targetvar))
all(rownames(clin)==rownames(mirna))
all(rownames(clin)==rownames(mutation))
all(rownames(clin)==rownames(cnv))
all(rownames(clin)==rownames(rna))



# predict.ranger does not accept "-" in variable names. Therefore, replace
# all occurences of "-" by "." in the variable names:

names(clin) <- gsub("-", ".", names(clin))
colnames(mirna) <- gsub("-", ".", colnames(mirna))
colnames(mutation) <- gsub("-", ".", colnames(mutation))
colnames(cnv) <- gsub("-", ".", colnames(cnv))
colnames(rna) <- gsub("-", ".", colnames(rna))



# Save processed data:

save(clin, targetvar, mirna, mutation, cnv, rna, file=paste("./Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION/Data/ProcessedData/", cancertypeshort, ".Rda", sep=""))


