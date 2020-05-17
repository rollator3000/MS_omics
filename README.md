# A comparison study of prediction approaches for multiple training data sets & test data with block-wise missing values
This is the README to the repository of Frederik Ludwigs' Master-Thesis: <br>
***A comparison study of prediction approaches for multiple training data sets & test data with block-wise missing values*** <br> 
supervised by: <br>
***Dr. rer. nat. Roman Hornung - Ludwig-Maximilians University - IBE***  
<br> 
Block-wise missingness is a common problem in the context of Multi-Omics Data. To my knowledge there are no standard approaches, nor compariosons studies for this type of missingness yet. In 2018 Norbert Krautenbacher has already stated that a reliable analysis strategy for multi-omics data with block-wise missingness are urgently needed! This thesis aims to provide such a comparison study and shall help finding a reliable analysis strategy.

---

## Project description
This project compares different approaches capable to deal with block-wise missingness in Multi-Omics data.  
Block-wise missingness is a special type of missingness that appears frequently in the context of Multi-Omics data. 
It can for example arise when concatenating multiple clinical studies with the same target variable. Eventhough the datasets from the different studies have the same target variable, the observed features can still can differ! The concatination of such datasets results then in a DF with block-wise missingness!  

Data with blockwise missingness always consits of different **folds** and **blocks**.
  - A **block** describes a set of covariates containing all features collected on the basis of a characteristic.  
    Basically all covariates that are related in content (e.g. *physical properties*: Height & Weight | *educational properties*: Income & Education').  
  - A **fold** represents a set of observations with the same observed blocks.  
    Basically all observations with the same observed features. Each fold is unique and every obserbation belongs to exactly one of them.

### Example for a dataset with blockwise missingness:  
| ID  | Weight  | Height  | Income  | Education   | g1      | ...   | g100    | Y   |
|---- |-------- |-------- |-------- |-----------  |-------  |-----  |-------  |---  |
| 1   | 65.4    | 187     | 2.536   | Upper       |         |       |         | 1   |
| 2   | 83.9    | 192     | 1.342   | Lower       |         |       |         | 0   |
| 3   | 67.4    | 167     | 5.332   | Upper       |         |       |         | 1   |
| 4   |         |         | 743     | Lower       | -0.42   | ...   | 1.43    | 1   |
| 5   |         |         | 2.125   | Lower       | 0.52    | ...   | -1.37   | 0   |
| 6   | 105.2   | 175     |         |             | -1.53   | ...   | 2.01    | 0   |
| 7   | 71.5    | 173     |         |             | 0.93    | ...   | 0.53    | 0   |
| 8   | 73.0    | 169     |         |             | 0.31    | ...   | -0.07   | 1   |

  - The data consits of three feature-blocks:
     - **Physical properties:**     Weight, Height
     - **Educational properties:**  Income, Education
     - **Biological properties:**   g1, ..., g100
  -  The data consits of three folds:
     - **Fold1:** All observations with observed Physical & Educational properties
     - **Fold2:** All observations with observed Educational & Biological properties
     - **Fold3:** All observations with observed Physical & Biological properties
  

Regular model fitting on data with block-wise missingness is for most statistical appropaches not directly possible, so that either the method needs to be adjusted or the data processed! Besides the data we use to fit a model, the testdata can also consist of block-wise missingness - e.g. 1 block avaible; Combination of 2 blocks avaible; .... Therefore the approaches must be able to deal with block-wise missing data in the test data as well <br>

### Approaches:
-  **Complete Case Approach:** Only use complete cases - regarding the testset - to fit a RF
   - if the testset consists of 'Weight', 'Height' & 'g1',...,'g100' then only the observations with these observed features are used for the trainig of the model!  
-  **Single Block Approach:** Only use the features from a single feature block to fit a RF
   - if the testset consists of 'Weight', 'Height' & 'g1',...,'g100' then fit a model on all observations with the features 'Weight' & 'Height' **OR** on all observations with the features 'g1',...,'g100'. 
   - For predicitons on the test data only use the features the model has been trained with and discard all other variables in test!
   - If a test observation misses a feature the RF has been trained with predictions are not possible
-  **Imputation Approach:** Use the 'missForest' approach to impute the missing values and fit a RF on this fully observed data then
   - Impute the missing data in the TrainingSet with the missForest Approach 
   - For predicition on testset, remove all features from the (imputed) trainset that are not part of the testst
   - On this pruned (imputed) trainset fit a RF and generate predicitions for the testset then
-  **Block-Wise Approach:** Fit a seperate RF on each feature block and create a final prediciton by combining the different block-wise predicitons
   - On each feature block of the data, fit a seperate RF *- one RF on the Physical properties, one RF on the Educational properties, ...*
   - For a prediction, each block-wise RF is asked for a predicition - only the RFs that were trained on a feature-block that is available for the test-set can return a predicition 
   - Average the seperate block-wise predicitons for a final prediciton - weighted/ unweighted
-  **Fold-Wise Approach:** Fit a seperate RF on each fold and create a final prediciton by combining the different fold-wise predicitons
   - On each fold of the data fit a seperate RF *- one RF on Fold1, one RF on Fold2, ...*
   - For a prediction, each fold-wise RF generates a seperate prediciton *- for these predicitons it might be that the single decision trees the RF consists of need to pruned.* 
   - **Pruning:** If a decision tree uses a split variable that is not avaible in the testset, cut the decision tree and use the node before that split as new terminal node
   - Average the seperate fold-wise predicitons for a final prediciton - weighted/ unweighted

### ! ! ! Closer Information to the approaches, aswell as to the results are in the MS-Thesis itself! ! !

---

## Data
#### TCGA
The original TCGA data is not part of this Repo. If interested in the original data send an E-Mail to 'f.ludwigs@yahoo.de'.  
Only the subsetted TCGA data can be found in the repository under: "data/processed/TCGA_subset_12345"   
With the script 'code/TCGA_03_subset_DFs' the data is subsetted & the script 'code/TCGA_04_TestTrain_splits_on_subsetted_DFs.R' splits the subsetted DFs to test- and train-set, whereby the training set is induced with block-wise missingness! Based on these test-train splits the predictive performance of the different approaches are investigated.

#### Real Data
The 'real' dataset used in the thesis comes from a coperation with the 'Hospital of the University of Munich'. For data protection reasons, the data can not be stored in the repository nor be shared.

---

## Code
Short describtion of the scripts in './code'!  
The code scripts either refer to the 'TCGA' data, the 'real' data or is 'general' script needed for both data sources!

#### General
``` 
- GENERAL_DecisionTreeExample:
    Get an example for the featurespace splitting of a single decision tree

- GENERAL_simpleRF_adaption:
    Implement RF function where the single trees can be pruned.
      --> builds up on 'github.com/mnwright/simpleRF'
``` 

#### TCGA
``` 
- TCGA_01_data_overview: 
    Get a overview of the different DFs in TCGA & get some additional Information!

- TCGA_02_explorative_performance:
    Check the performance of a RF when trained on a single feature-block/ on all blocks   
    joint feature blocks. Then create different subsets of the original data & get the 
    the predictive performance on these subsets [incl. plots]

- TCGA_03_subset_DFs:
    Create final subsets of the original TCGA datasets & get the predictive performance
    on the final subsets (single-block & joint-block). Plot the results and get dimensions
    of the subsetted DFs

- TCGA_04_TestTrain_splits_on_subsetted_DFs:
    Split the subsetted DFs into Test-Train splits for the CV. The training data is induced
    with block-wise missingness. The resulting Test-Train-Splits can then be used for the CV 
    of the different approaches.  

- TCGA_05_Foldwise_RF_CV:
    CrossValidate the foldwise Approach for all different settings & all possible combinations 
    of blockwise missingness in the TestSet

- TCGA_06_BlockWise_RF_CV:
    CrossValidate the blockwise Approach for all different settings & all possible combinations 
    of blockwise missingness in the TestSet

- TCGA_07_CompleteCases_CV:
    CrossValidate the complete case Approach for all different settings & all possible combinations 
    of blockwise missingness in the TestSet

- TCGA_08_Imputation_CV:
    CrossValidate the Imputation Approach for all different settings & all possible combinations 
    of blockwise missingness in the TestSet

- TCGA_09_SingleBlock_CV:
   CrossValidate the SingleBlock Approach for all different settings & all possible combinations 
    of blockwise missingness in the TestSet

- TCGA_10_Plot_CV_Results:
    Vizualize the Results from the CV of the TCGA data for all the different approaches!
``` 

#### REAL
``` 
- REAL_01_Fold-Wise-Approach: 
    CrossValidate the foldwise Approach on the real dataset 

- REAL_02_Block-Wise-Approach:
    CrossValidate the blockwise Approach on the real dataset 

- REAL_03_Imputation-Approach:
   CrossValidate the Imputation Approach on the real dataset 

- REAL_04_CompleteCase-Approach:
    CrossValidate the CompleteCase Approach on the real dataset

- REAL_05_Single-Block-Approach:
    CrossValidate the SingleBlock Approach on the real dataset

- REAL_06_Plot_CV_Results:
    Vizualize the Results from the CV of the real data for all the different approaches!
``` 
---

## Project Organization
```
    ├── LICENSE
    ├── .git               <- Folder for VersionControl [created by GIT]
    ├── README.md          <- Top-level 'README' for developers
    │
    ├── code               <- All the R-Code needed to run this repository! 
    │   │                     Overview of the code: 'Code' section in README
    │   │     
    │   └─ no_main         <- Code that was used as template/ inspiration 
    │  
    ├── data               <- All data used in this project!
    │   │
    │   ├── external       <- External raw Data from not public sources!
    │   │   │  
    │   │   └─ TCGA        <- Processed & fully observed TCGA data from
    │   │                     Dr. Hornung!
    │   │    
    │   ├── processed      <- Processed Data that can be directly used!
    │   │    │
    │   │    ├─ real_data  <- real clincal data
    │   │    │
    │   │    └─ TCGA_subset_12345    <- Subsetted TCGA Data - single blocks were subsetted. 
    │   │        │                      Basically original data with lower dimensions!   
    │   │        │
    │   │        ├─ missingness_1234 <- Subsetted data split to test and train, whereby
    │   │        │                      the trainsplit contains blockwise missingness 
    │   │        │
    │   │        └─ missingness_1234_imputed  <- TestTrainSplits of the subsetted data,
    │   │                                        whereby the blockwise missingness in 
    │   │                                        train was removed by imputation!
    │   │
    │   └─ interim         <- internal data - e.g. dummy data for implementation
    │ 
    ├── sources            <- All sources used in this project
    │   │ 
    │   ├── GENERAL        <- General Sources - not related to Block-Wise Missingness 
    │   │
    │   ├── INTRO          <- Sources needed for the introduction
    │   │
    │   ├── Metrics        <- Sources to the used metrics
    │   │
    │   ├── Omics_Methods  <- Sources related to Block-Wise missingness 
    │   │                     in Multi-Omics data
    │   │
    │   ├── Online Sources <- Online sources from the internet 
    │   │                    - downloaded whole web-pages
    │   │
    │   └── revieved_MS    <- reviewed versions of the MS thesis
    │
    └── docs               <- Documents used within this repository! 
        │                     Explainatory/ Describtive material, Notes etc.
        │
        ├─ Graphics        <- Folder for each section in the MS-Thesis - all 
        │                     graphics that have been used in the corrsponding 
        │                     section or it belongs to the secction! 
        │ 
        └─ CV_Res          <- Results of CrossValidations
            │
            ├─ REAL        <- CV-Results for the different approaches on the 'REAL' data
            │
            └── TCGA       <- CV-Results for the different approaches on the 'TCGA' data
                 │
                 ├── BlockWise_Approach    <- CV Results of the BlockWise approach for the 
                 │                            different pattern of missingness in TCGA
                 │
                 ├── CompleteCase_Approach <- CV Results of the CC approach for the 
                 │                            different pattern of missingness in TCGA
                 │
                 ├── explorative_subsets   <- CV Results on joint-/ single blocks - for 
                 │                            different fractions of subsets for the blocks
                 |
                 ├── final_subsets         <- CV Results on joint-/ single blocks for the 
                 |                            the final subsets of TCGA used for the whole CV!
                 │
                 ├── FoldWise_Approach     <- CV Results of the FoldWise approach for the 
                 │                            different pattern of missingness in TCGA
                 │
                 ├── Imputation_Approach   <- CV Results of the Imputation approach for the 
                 │                            different pattern of missingness in TCGA
                 │
                 └── SingleBlock_Approach  <- CV Results of the SingleBlock approach for  
                                              the different pattern of missingness in TCGA
``` 