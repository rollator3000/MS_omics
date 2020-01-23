########################################################

# NOTE: Before the code can be excecuted, the R working directory *MUST* 
# be set to the directory that contains the folder 'Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION':

# Remove '#' from the two lines below and replace 'Z:/here/is/my/path/' by the path
# to the directory that contains 'Additional_file_2_HornungWright_WITH_DATA_EXECUTABLE_VERSION' ( the path 'Z:/here/is/my/path/'
# *MUST* end with a slash / ):

# basedir <- "Z:/here/is/my/path/"
# setwd(basedir)

########################################################



cancertypeshort <- "TGCT"
cancertypelong <- "Testicular_Germ_Cell_Tumors_TGCT"


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

# --> Only one non-censored observation.



# --> Data set cannot be used for analysis.
############################################