# A comparison study of prediction approaches for multiple training data sets & test data with block-wise missing values
This is the README to the repository of Frederik Ludwigs' Master-Thesis: <br>
***A comparison study of prediction approaches for multiple training data sets & test data with block-wise missing values*** <br> 
supervised by: <br>
***Dr. rer. nat. Roman Hornung - Ludwig-Maximilians University - IBE***

## Project description
This project compares different adaptions of the RandomForest-Method & other approaches to deal with blockwise missing data! A dataset with blockwise missingness, consists of different folds, where each fold has different observed features - basically, it's like having different training-sets for the same response! *(different observed features for each training-set, where each trainingset has the same response value)*  <br>
As in these settings imputation techniques are only reliable if the different datasets are not too heterogenous & the set of covariates highly overlap, we try a new approach by trying to incorporate the different features from the different folds/ training-sets into a single RandomForest / Approach!

**Example:**
Different hospitals do reseach regarding the same response *(e.g. Breast Cancer)*, but the different hospitals do not necessarily collect the same features - in this setting the data from the different hospitals can be seen as the data from different folds - where none of the folds have the exact same feature space!
``` 
- Hospital_1: Clinical + RNA Data
- Hospital_2: Clinical + CopyNumberVariation Data
- Hospital_3: Clinical + CopyNumberVariation + Mutation Data
```
For Details regarinding the methods, data situations etc. please have a look at: <br> 
     MasterThesis / code in repository

## Code
Short describtion of the scripts in './code'!
``` 
- 01_data_overview: 
    Get a overview of the different DFs & extract the ones usable for our project + 
    some extra information to the different DFs

- 02_predictie_power_on_gender_w_all_feas_of_a_block:
    Check how good the predicitons of a RF are, when it is trained on single blocks, 
    on all blocks joint, single subsetted blocks,  all subsetted blocks joint....

- 03_create_artifical_DF:   
    1. Create an artifical DF, used for the implementation of the pruning!
    2. Do final Subsets of the processed omics data
    3. Get the Performance of the final subsets when fitting a model on
        - all omics blocks joint together
        - each single omics block! 

- 04_simpleRF_adaption:
    Implementation of the RF-Adjustment + additional functions needed for evaluation,
    creating trees, ....

- 05_RF_Adaption_CV:
    Code to CrossValidate the adjusted RF-Algorithm from Dr. Hornung for different
    Training Settings & all possible combinations of blockwise missingness in 
    the TestSet

- 06_RF_Krautenbacher_CV:
    Code to CrossValidate the adjusted RF-Algorithm from Dr. Krautenbacher 
    for different Training Settings & all possible combinations of blockwise 
    missingness in the TestSet
```

## Project Organization
------------
    ├── LICENSE
    ├── .git               <- Folder for VersionControl [created by GIT]
    ├── README.md          <- Top-level 'README' for developers working with
    │                         this repository
    │
    ├── code               <- Folder w/ all the R-Code needed to run this 
    │   │                     repository! For an overview see:
    │   │                     'Code' section in this README
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
    │   │      │             in this project!
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
        │                      Explainatory/ Describtive material, Notes etc.
        │ 
        ├─ explorative_subsets <- Performance of a RF with differnt percent 
        │                         of subsets! 
        │                         --> Explorative Results! To find out which
        │                             blocks do need which subsets for our 
        │                             project!
        │
        └── performance_final_subsets <- RF Performance on the final subsetted DFs
                                         -> once for joint blocks
                                            [all blocks as features to RF]
                                         -> once for single blocks
                                            [fit a seperate RF on each of the blocks!] 
--------
