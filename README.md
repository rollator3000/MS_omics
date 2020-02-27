# A comparison study of prediction approaches for multiple training data sets & test data with block-wise missing values
This is the README to the repository of Frederik Ludwigs' Master-Thesis: <br>
***A comparison study of prediction approaches for multiple training data sets & test data with block-wise missing values*** <br> 
supervised by: <br>
***Dr. rer. nat. Roman Hornung - Ludwig-Maximilians University - IBE***

## Project description
This project compares different approaches to deal with blockwise missingness in training- & testsets! <br>
When having multiple training-sets for the same target variable usually all these training sets recored different feature variables. These Trainingsets usually consit of differnet blocks. A block is a set of features from a Topic (e.g. DNA, Clinical,...) and when combining datasets with different feature blocks the concatenation results in a dataframe with blockwise missing data. For each observation in the concatenated dataframe a feature block is either fully observed or not all. As a regular model fitting on these data is not possible either the method needs to be adjusted or the data processed! <br> 
Multiple Approaches will be compared:
-  Baseline Approach: Only use complete cases (regarding the testset) to fit a RF on *(e.g. testset constits of Clinical Block and 2 Omics Block -> only use train observations, that were fully observed in these 3 blocks)*. <br>
-  Baseline Approach: Only use complete cases from 1 single feature block to fit a RF on *(e.g. only observations, that were fully observed in the 'clin'-block)* <br>
-  Imputation Approach: Use the 'missForest' approach to impute the missing values, so that we can regualry fit a RF on it. <br>
-  Adaption Approach: Fit multiple RandomForests blockwise on the concatenated data <br>
-  Adaption Approach: Fit multiple RandomForests foldwise on the concatenated data <br>

Closer Information to the approaches, aswell as to the results are in the MS-Thesis itself!

**Example**
Different hospitals do research regarding the same response *(e.g. different BreastCancer)*, but these different hospitals do not necessarily collect the same feature for this - e.g. Hospital_1 collects clinical & RNA data whereas Hospital_2 collects  clinical & CNV data, etc.. <br>
``` 
- Hospital_1: Clinical + RNA Data
- Hospital_2: Clinical + CopyNumberVariation Data
- Hospital_3: CopyNumberVariation + Mutation Data
```
The data from the different hospitals can be seen as different folds - where none of the folds have the exact same feature space. In these cases it could be benefical to have approaches, that can learn with all the avaible features at once. So that we fit a RandomForest on the data, eventhough not all features were observed for all observations:
```
RandomForest(formula = Cancertype ~ Clinical + RNA + CopyNumberVariation + Mutation, 
             data    = rbind(Hospital_1, Hospital_2, Hospital_3)) 
```
<br>
- The two baseline approaches try to process the data in a way that the data only consits of complete cases. For this variables and/ or observations can be removed from the original concatenated data. On the resulting dataframe a regular RF can be trained with.
- The imputation approach imputes the missing values in the concatenated data, so that the resulting dataframe has no missing values anymore. On this imputed dataframe a regular RF can be trained.
- I the first RF-adaption (blockwise-adaption) a RF is fit on each feature block of the concatenated trainingdata seperatly. For the fitting of the blockwise-RFs it only uses complete cases from each block. When doing predicitons on a testset it creates predicitons based on each block in the testset & aggregates these to a final prediciton - e.g. TestSet consits of two blocks. Block-Wise RF creates a prediciton for each block, and then aggregates the seperate predicitons to a final one!
- In the second RF-adaption (foldwise-adaption) a RF is fit on each fold of the the concatenated trainingdata seperatly. So a RF is fit on each fold [fold is set of observations with the same features] - e.g. fit a RF to all data from Hospital_1, fit a RF to all data from  Hospital_2, etc.. For the predicitons on testdata, we need to aggregate the predicitons from the different foldwise fitted RFs. If the testset is missing any feature, that is used as splitvariable in any of the foldwise fitted RFs, the RFs have to be pruned (details see MS-Thesis). 

As the testdata can also consist of block-wise missingness - e.g. all blocks avaible; just 1 block avaible; Combination of 2 blocks avaibale - the approaches must be able to deal with this aswell <br>
For Details regarding the methods, the data, the metrics and the results,  please have a look at the **MS-Thesis**

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
