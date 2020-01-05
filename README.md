# A comparison study of prediction approaches for multiple training data sets & test data with block-wise missing values
This is the README to the repository of Frederik Ludwigs' Master-Thesis, supervised by: <br>
***Dr. rer. nat. Roman Hornung - Ludwig-Maximilians University - IBE***

## Project description
This project compares different methods to deal with blockwise missing data! In a dataset with blockwise missingness, for different folds *(collection of observations)* different features were collected! So basically, multiple training-sets *(not necessarily sharing features)* with the same response! In these settings it could be benefical for prediciton models to incorporate the different features from the different folds/ training-sets ! <p/>
**Example:**
Different hospitals do reseach regarding the same response *(e.g. Breast Cancer)*, but the different hospitals *(different folds)* do not collect the same features:
``` 
- Hospital_1: Clinical + RNA Data
- Hospital_2: Clinical + CopyNumberVariation Data
- Hospital_3: Clinical + CopyNumberVariation + Mutation Data
```
In this setting we compare different adaptions of the RandomForest algorithm to deal with these kind of settings without imputation *(often unreliable thesesettings)* <p/>
For Details regarinding methods, data situations etc. please have a look at: <br> 
MasterThesis / code in repository

## Code
Short describtion of the scripts in './code'!
``` 
- 01_data_overview: 
    Get a overview of the different DFs & extract the ones usable for our project + some extra information to the different DFs

- 02_predictie_power_on_gender_w_all_feas_of_a_block:
    Check how good the predicitons of a RF are, when it is trained on single blocks, on all blocks joint, single subsetted blocks,  all subsetted blocks joint....

- 03_create_artifical_DF:   
    Create an artifical DF, used for the implementation of the pruning!

- 04_simpleRF_adaption:
    Implementation of the RF-Adjustment + additional functions needed for evaluation, creating trees, ....

- 05_RF_Adaption_CV:
    Code to CrossValidate the adjusted RF-Algorithm from Dr. Hornung

- 06_RF_Krautenbacher_CV:
    Code to CrossValidate the adjusted RF-Algorithm from Dr. Krautenbacher
```

## Project Organization
------------
    ├── LICENSE
    ├── .git               <- Folder for VersionControl [created by GIT]
    ├── README.md          <- Top-level 'README' for developers working with
    │                         this repository
    │
    ├── code               <- Folder w/ all the R-Code needed to run this 
    │   │                     repository! For an overview & details see 'Code'
    │   │                     section in this README
    │   │
    │   ├── get_data      <- Code from Dr. Hornung to preprocess and download 
    │   │                    multi-omics data used in this project!
    │   │     
    │   └── no_main        <- Subfolder, with all the code, that was used as 
    │                         template/ Inspiration for final code in './code'
    │  
    ├── data               <- All data used in this project!
    │   │   
    │   ├── raw            <- The original, immutable data dump - empty
    │   │
    │   ├── external       <- Data from third party sources
    │   │   | 
    │   │   ├── example_data <- data used for investigation & adjusting code!
    │   │   | 
    │   │   └─ Dr_Hornung <- Clinical Omics-Data in preprocessed form, used
    │   │      │             in this project! - Code to download and preprocess
    │   │      │             by Dr.Hornung [also in './code/get_data'] 
    │   │      │
    │   │      └─ Data/ProcessedData  <- processed data used in this project!
    │   │
    │   ├── processed      <- Processed Data that can be directly used!
    │   │
    │   └─ interim         <-  internal data - e.g. temporarily saved data
    │ 
    ├── sources            <- Folder, that holds all sources, used in/ for
    │   │                     this project!
    │   │   
    │   ├── online         <- Online sources from the internet 
    │   │                     - could not be downloaded!
    │   │                     - downloaded whole page with 'Google Chrome'
    │   │
    │   └── PDF            <- All sources that could be saved/ dowloaded as PDF
    │
    └── docs               <- Documents used within this repository! 
                              Explainatory/ Describtive material, Notes etc.


--------
