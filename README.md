# A comparison study of prediction approaches for multiple training data sets & test data with block-wise missing values
This is the README to the repository of Frederik Ludwigs' Master-Thesis: <br>
***A comparison study of prediction approaches for multiple training data sets & test data with block-wise missing values*** <br> 
supervised by: <br>
***Dr. rer. nat. Roman Hornung - Ludwig-Maximilians University - IBE***

## Project description
This project compares different approaches capable to deal with blockwise missingness in Multi-Omics data.  
Blockwise missingness is a special type of missingness that appears frequently in the context of Multi-Omics data.  
Data with blockwise missingness always consits of different **folds** and **blocks**.
  - A block describes a set of covariatescontaining all features collected on the basis of a characteristic. Basically all covariates that are related in content (e.g. physical properties: height, weight, skin).  
  - A fold represents a set of observations with the same observed blocks. Basically all observations with the same observed features - each fold is unique and every obserbation in the data belongs to exactly one of them.

Blockwise Missingness can for example arise when concatenating multiple training-sets for the same target variable. And could have the following form:
| ID  | Weight  | Height  | Income  | Education   | g1      | ...   | g100    | Y   |
|---- |-------- |-------- |-------- |-----------  |-------  |-----  |-------  |---  |
| 1   | 65.4    | 187     | 2.536   | Upper       |         | ...   |         | 1   |
| 2   | 83.9    | 192     | 1.342   | Lower       |         | ...   |         | 0   |
| 3   | 67.4    | 167     | 5.332   | Upper       |         | ...   |         | 1   |
| 4   |         |         | 743     | Lower       | -0.42   | ...   | 1.43    | 1   |
| 5   |         |         | 2.125   | Lower       | 0.52    | ...   | -1.37   | 0   |
| 6   | 105.2   | 175     |         |             | -1.53   | ...   | 2.01    | 0   |
| 7   | 71.5    | 173     |         |             | 0.93    | ...   | 0.53    | 0   |
| 8   | 73.0    | 169     |         |             | 0.31    | ...   | -0.07   | 1   |

Regular model fitting on these kind of data is for most methods not directly possible, so that either the method needs to be adjusted or the data processed! 
As the testdata can also consist of block-wise missingness - e.g. all blocks are avaible | just 1 block avaible | Combination of 2 blocks avaible - the approaches must be able to deal with this aswell <br>

Multiple Approaches will be compared:
-  Baseline Approach 1: Only use complete cases (regarding the testset) to fit a RF on - if the testset consists of 'Weight', 'Height', 
'g1',...,'g100' then only Observations with these observed features will be used for the trainig of the model!
-  Baseline Approach 2: Only use complete cases from a single feature block to fit a RF on - if the testset consists of 'Weight', 'Height', 
'g1',...,'g100' then fit a model, that uses all observations with the features 'Weight' & 'Height' and do predicitons on the test data only using these variables. Analog fit a model, that uses all observations with the features 'g1',...,'g100' and do predicitons on the test data only using these variables. 
-  Imputation Approach: Use the 'missForest' approach to impute the missing values, so that we can regualry fit a RF on it - Remove all features that are not part of the testst and then Impute the missing values in the train set and fit a regular RF on the imputed DF.
-  Adaption Approach: Fit multiple RandomForests blockwise on the concatenated data - fit a RF on each feature block of the data and aggregatre the predicitons of these when prediciton on testdata  
-  Adaption Approach: Fit multiple RandomForests foldwise on the concatenated data. - fit a RF on each fold of the data and aggregate the predicitons of these when prediciton on testdata 

***Closer Information to the approaches, aswell as to the results are in the MS-Thesis itself!***

<br>
-------------------- STOPPED HERE
## Code
Short describtion of the scripts in './code'!
``` 
- 01_data_overview: 
    Get a overview of the different DFs.
    Get the DFs we can use for our work, and get some 
    additional Information!

- 02_explorative_performance:
    Check how good the predicitons of a RF are, when it 
    is trained on single blocks, on all blocks joint, 
    single subsetted blocks, all subsetted blocks joint....

- 03_create_needed_DFs:
    1. Do final Subsets of the processed omics data
    2. Get the Performance of the final subsets when 
       fitting a model on:
           - all omics blocks joint together
           - each single omics block! 

- 04_simpleRF_adaption:
    Implementation of the RF-Adjustment + additional functions 
    needed for evaluation, creating trees, ....

- 05_RF_Adaption_CV:
    Code to CrossValidate the adjusted RF-Algorithm from Dr. 
    Hornung (foldwise) for different Training Settings & all  
    possible combinations of blockwise missingness in the TestSet

- 06_RF_Krautenbacher_CV:
    Code to CrossValidate the adjusted RF-Algorithm from Dr. 
    Krautenbacher (blockwise) for different Training Settings &
    all possible combinations of blockwise missingness in the TestSet

- 07_Vizualize_results:
    Code to plot the results of the different approaches!

- 08_Imputatiojn_CV:
    Code to CrossValidate the (missforest-) imputation approach! 

- 09_create_data_w_blockwise_missingness_for_CV:
    Script to create the final subsets of the original data,
    split the data to test and train & inducde the blockwise
    missingness into the trainsets of the splitted DFs!
```

## Data
The original data is not part of this Repo, only the subsetted data!
If interested in the original data, please contact me on github or send 
an E-Mail to 'f.ludwigs@yahoo.de'.

## Project Organization
------------
    ├── LICENSE
    ├── .git               <- Folder for VersionControl [created by GIT]
    ├── README.md          <- Top-level 'README' for developers working with
    │                         this repository
    │
    ├── code               <- All the R-Code needed to run this repository! 
    │   │                     For an overview of the code in here see:
    │   │                     'Code' section in README
    │   │     
    │   └── no_main        <- All the code, that was used as template/ 
    │                         inspiration for final scripts in './code'
    │  
    ├── data               <- All data used in this project!
    │   │   
    │   ├── raw            <- The original, immutable data dump - empty!
    │   │
    │   ├── external       <- External Data from not public sources
    │   │   │  
    │   │   └─ Dr_Hornung  <- Clinical preprocessed & fully observed Omics-Data 
    │   │                     from Dr. Hornung - used it in 'BlockForest' 
    │   │                     originally - original preprocessed data!
    │   │    
    │   ├── processed      <- Processed Data that can be directly used!
    │   │    │
    │   │    └─RH_subsetted_12345 <- Data from Dr. Hornung where the single 
    │   │        │                   blocks were subsetted already! Seed 
    │   │        │                   used for this was the '12345'   
    │   │        │
    │   │        └─ missingness_WXYZ <- TestTrainSplits of the subsetted data
    │   │                               with the seed 'WXYZ', where each DF is       
    │   │                               saved 4 times - once for each scenario
    │   │                               we use to induce blockwise missingness! 
    │   │
    │   └─ interim         <- internal data - e.g. temporarily saved data
    │         │
    │         └─ example_data <- files and DFs used to implement the foldwise RF
    │                            was used for investigation etc.
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
        │                     Explainatory/ Describtive material, Notes etc.
        │
        ├─ Graphics        <- Folder for each section in the MS-Thesis with all 
        │                     graphics that have been used in this section or
        │                     belog to this secction! 
        │ 
        └─ CV_Res          <- Results of different CVs/ explorative subsetting 
            │
            └── gender    <- All results where gender was the response for the data
                  │          for the data from Dr_Hornung! 
                  │
                  ├── explorative  <- results from the explorative subsetting
                  │                   done to find out which blocks predict how
                  │                   strong!
                  │                    - Joint / Single Blocks as feature space!
                  │
                  ├── final_subsets <- RF Performance on the final subsetted DFs
                  │                    [want to use them to CV new approaches]
                  │                     - Joint / Single Blocks as feature space!
                  │
                  ├── Norbert_final_subsets <- Results w/ the approach from Norberts
                  │                            RF Adjustment - BlockWise Fitting
                  │
                  └── Roman_final_subsets   <- Results w/ the approach from Romans
                                               RF Adjustment - <FoldWise Fitting 
-------- 
