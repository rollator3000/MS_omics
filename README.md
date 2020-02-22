# A comparison study of prediction approaches for multiple training data sets & test data with block-wise missing values
This is the README to the repository of Frederik Ludwigs' Master-Thesis: <br>
***A comparison study of prediction approaches for multiple training data sets & test data with block-wise missing values*** <br> 
supervised by: <br>
***Dr. rer. nat. Roman Hornung - Ludwig-Maximilians University - IBE***

## Project description
This project compares different approaches to deal with blockwise missing data in training- & testsets! <br>
When having multiple different training-sets for the same target variable, the concatenation of these training-sets results in a dataframe with missing data. For each observation in the concatenated dataframe each feature block is either fully observed or not all. In these settings either the method needs to be adjusted or the data processed, before a model can be fit on it! <br> 
Multiple Approaches will be compared:
- 1. Baseline Approach: Only use complete cases (regarding the testset) to fit a RF on *(e.g. testset constits of Clinical Block and 2 Omics Block -> only use train observations, that were fully observed in these 3 blocks)*. <br>
- 2. Baseline Approach: Only use complete cases from 1 single feature block to fit a RF on *(e.g. only observations, that were fully observed in the 'clin'-block)* <br>
- 3. Imputation Approach: Use the 'missForest' approach to impute the missing values, so that we can regualry fit a RF on it. <br>
- 4. Adaption Approach: Fit multiple RandomForests blockwise on the concatenated data <br>
- 5. Adaption Approach: Fit multiple RandomForests foldwise on the concatenated data <br>
Approaches 4. /5. are trying to incorporate the features from the different folds/ blocks into a single Approach! <br>
<br>
Closer Information to the approaches, aswell as the results, can be seen in the MS-Thesis itself!

**Example**
Different hospitals do research regarding the same response *(e.g. different Cancertypes)*, but these hospitals do not necessarily collect the same feature blocks *(e.g. omics data)* - the data from the different hospitals can be seen as different folds - where none of the folds have the exact same feature space!
``` 
- Hospital_1: Clinical + RNA Data
- Hospital_2: Clinical + CopyNumberVariation Data
- Hospital_3: CopyNumberVariation + Mutation Data
```
Now it would benefical to have a model, that can learn with all the avaible features for the different folds!
```
RandomForest(Cancertype ~ Clinical + RNA + CopyNumberVariation + Mutation) 
```
The first set of approaches *(1., 2., 3.)* try to process the data in a way, that the block-wise missingness disappears from the concatenated dataframe.<br>
The 4. & 5. Approach try to take all features into account - eventhough some were not observed for all. <br>

The models also need to be able to predict on test-obs., with different combinations of observedblocks!
All blocks avaible, or just 1 block, or a combination of 2 blocks etc.. <br>
For Details regarding the methods, data situations etc. please have a look at the **MS-Thesis**

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
    │   │   └─ Dr_Hornung <- Clinical Omics-Data from Dr. Hornung
    │   │      │
    │   │      ├─ original_processed_data  <- original preprocessed data
    │   │      │ 
    │   │      └─ subsetted_12345 <-  subsetted data used for CV etc.
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
        └─ CV_Res          <- Results of different CVs/ explorative subsetting 
            │
            └── gender    <- all results where gender was the response
                  │
                  ├── explorative  <- results from the subsetting, we did 
                  │                   to find out which blocks predict how
                  │                   strong!
                  │
                  ├── final_subsets <- RF Performance on the final subsetted DFs
                  │                    [want to use them to CV new approaches]
                  │                     -> once for joint blocks
                  │                        [all blocks as features to RF]
                  │                     -> once for single blocks
                  │                        [fit a seperate RF on each of the blocks!] 
                  │
                  ├── Norbert_final_subsets <- Results w/ the approach from Norberts
                  │                            RF Adjustment 
                  │
                  └── Roman_final_subsets   <- Results w/ the approach from Romans
                                               RF Adjustment
-------- 
