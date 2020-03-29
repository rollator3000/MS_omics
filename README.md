# A comparison study of prediction approaches for multiple training data sets & test data with block-wise missing values
This is the README to the repository of Frederik Ludwigs' Master-Thesis: <br>
***A comparison study of prediction approaches for multiple training data sets & test data with block-wise missing values*** <br> 
supervised by: <br>
***Dr. rer. nat. Roman Hornung - Ludwig-Maximilians University - IBE***

---

## Project description
This project compares different approaches capable to deal with block-wise missingness in Multi-Omics data.  
Block-wise missingness is a special type of missingness that appears frequently in the context of Multi-Omics data.  
Data with blockwise missingness always consits of different **folds** and **blocks**.
  - A **block** describes a set of covariates containing all features collected on the basis of a characteristic.  
    Basically all covariates that are related in content (e.g. physical properties: Height & Weight | educational properties: Income & Education').  
  - A **fold** represents a set of observations with the same observed blocks.  
    Basically all observations with the same observed features.  
    Each fold is unique and every obserbation belongs to exactly one of them.

#### A dataset with blockwise missingness could habe the following form:  
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

  - Consits of three blocks:
     - Physical properties:   Weight, Height
     - Education properties:  Income, Education
     - Biological properties: g1, ..., g100
  - Consists of three folds:
     - Fold1: All observations with observed Physical & Education properties
     - Fold2: All observations with observed Education & Biological properties
     - Fold3: All observations with observed Physical & Biological properties

In the context of Multi-Omics Data block-wise missingness is a common problem, where reliable analysis strategies are urgently needed!  
Block-wise Missingness can for example arise when concatenating multiple clinical studies for the same target variable. Eventhough all different datasets do have the same target variable they still can differ in the collected features, such that the concatination of these results in a DF with block-wise missingness!  

Regular model fitting on data with block-wise missingness is for most statistical appropaches not directly possible, so that either the method needs to be adjusted or the data processed! As the testdata can also consist of block-wise missingness - e.g. 1 block avaible; Combination of 2 blocks avaible; ... - the approaches must be able to deal with block-wise missing data in the test data as well <br>

#### Approaches:
-  **Baseline Approach 1:** Only use complete cases (regarding the testset) to fit a RF
   - if the testset consists of 'Weight', 'Height' & 'g1',...,'g100' then only Observations with these observed features will be used for the trainig of the model!
-  **Baseline Approach 2:** Only use the features from a single feature block to fit a RF
   - if the testset consists of 'Weight', 'Height' & 'g1',...,'g100' then fit a model on all observations with the features 'Weight' & 'Height' **OR** on all observations with the features 'g1',...,'g100'. 
   - For predicitons only use the features the model has been trained with and discard all other variables in test!
   - If a observation misses a feature the RF has been trained with predictions are not possible
-  **Imputation Approach:** Use the 'missForest' approach to impute the missing values
   - Impute the missing data in the TrainingSet with the missForest Approach 
   - For predicition on testset, remove all features from the (imputed) trainset that are not part of the testst
   - On this pruned (imputed) trainset fit a RF and generate predicitions for testset
-  **Block-Wise Approach:** Fit a seperate RF on each feature block and create a final prediciton by combining the different block-wise predicitons
   - On each feature block of the data fit a seperate RF -1 RF on the Physical properties, 1 RF on the Education properties, ...
   - To predict on a TestSet, each block-wise RF generates a seperate prediciton that are combined for a final prediciton 
-  **Fold-Wise Approach:** Fit a seperate RF on each fold and create a final prediciton by combining the different fold-wise predicitons
   - On each fold of the data fit a seperate RF -1 RF on Fold1, 1 RF on Fold2, ...
   - To predict on a TestSet, each fold-wise RF generates a seperate prediciton - for these predicitons it might be that the single decision trees the RF consists of need to pruned - that are combined for a final prediciton 
   - Pruning: If a decision tree uses a split variable that is not avaible in the testset, cut the decision tree and use the node before that split as new terminal node


***! ! ! Closer Information to the approaches, aswell as to the results are in the MS-Thesis itself! ! !***

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
