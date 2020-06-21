# A comparison study of prediction approaches for multiple training data sets & test data with block-wise missing values
This is the README to the repository of Frederik Ludwigs' Master-Thesis:
A comparison study of prediction approaches for multiple training data sets & test data with block-wise missing values
supervised by:
Dr. rer. nat. Roman Hornung - Ludwig-Maximilians University - IBE

Block-wise missingness describes a special type of missingness that is common in the context of Multi-Omics data. To my knowledge, there are no standard approaches, nor comparison studies for this type of missingness yet. In 2018 Norbert Krautenbacher has already stated that a reliable analysis strategy for multi-omics data with block-wise missingness is urgently needed! This thesis aims to provide such a comparison study and shall help to find a reliable analysis strategy for data with block-wise missingness.

---

## Project description
This project compares different approaches capable to deal with block-wise missingness in Multi-Omics data. For this different random forest based adaptions that are capable to deal with block-wise missingness are introduced. Penalised regression adaptions from Hagenberg's thesis are also compared to the random forest based adaptions, even though these are not introduced theoretically.

### Block-wise missingness:
Block-wise missingness is a special type of missingness that appears frequently in the context of Multi-Omics data. 
It can, for example, arise when concatenating multiple clinical studies with the same target variable. Even though the datasets from the different studies have the same target variable, the observed features can still differ! The concatenation of such datasets results then in a DF with block-wise missingness!  

Data with blockwise missingness always consists of different **folds** and **blocks**.
  - A **block** describes a set of covariates containing all features collected based on a characteristic.  
    Basically all covariates that are related in content (e.g. *physical properties*: Height & Weight | *educational properties*: Income & Education').  
  - A **fold** represents a set of observations with the same observed blocks.  
    All observations with the same observed features. Each fold is unique and every observation belongs to exactly one of them.

#### Example for data with blockwise missingness:  
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
   
Regular model fitting on data with block-wise missingness is for most statistical approaches not directly possible so that either the method needs to be adjusted or the data processed! Besides the training, the test data can also consist of block-wise missingness. Therefore the approaches must be able to deal with block-wise missing data in the test data as well as in the train data. <br>

### Approaches:
The different random forest - RF - adaptions/ approaches are listed below and briefly explained using the example data above.

-  **Complete Case Approach:** Only use complete cases - regarding the testset - to fit a RF
   - If the test-set consists of the features 'Weight', 'Height' (Physical properties) & 'g1',...,'g100' (Biological properties) then only the observations with these observed features are used for the trainig of the model (Fold3)
   - The fitted model can then predict on the test-set regularly
-  **Single Block Approach:** Only use the features from a single feature block to fit a RF
   - If the test-set consists of the features 'Weight', 'Height' (Physical properties) & 'g1',...,'g100' (Biological properties) then a model can either be fit on all observations with the features 'Weight' & 'Height' (Physical properties) **OR** on all observations with the features 'g1',...,'g100' (Biological properties)
   - For predicitons on the test-set only use the features the model has been trained with and discard all other variables from the test-set
-  **Imputation Approach:** Impute the missing values with the 'missForest' approach and fit a RF - regarding the testset - on this fully observed data
   - Impute the missing data in the train-set with the missForest Approach
   - For predicition on test-set, remove all features from the fully observed train-set that are not part of the test-set
   - On this pruned (imputed) train-set fit a RF and predicit on the test-set then
-  **Block-Wise Approach:** Fit a seperate RF on each feature-block and create a final prediciton by combining the different block-wise predicitons
   - Fit a seperate RF on each feature block of the data *- one RF on the Physical properties, one RF on the Educational properties, ...*
   - For a prediction, each block-wise RF is asked for a predicition - only the RFs that were trained on a feature-block that is available for the test-set can predict on the test-set  
   - The seperate block-wise predicitons are averaged in a weighted/ unweighted way for the final predicitons
-  **Fold-Wise Approach:** Fit a seperate RF on each fold and create a final prediciton by combining the different fold-wise predicitons
   - Fit a seperate RF on each fold of the data  *- one RF on Fold1, one RF on Fold2, ...*
   - For a prediction, each fold-wise RF is asked for a predicition - only the RFs that were trained on a fold with at least one feature-block that is also in the test-set can can try it  
   - A RF that was trained on a fold that has a feature-block that is not available for the test-set might use split variables that are not available for the test-set can try to create a prediciton
   - The RFs might need to be pruned before they can genterate a prediciton
   - **Pruning:** If a decision tree of a RF uses a split variable that is not avaible in the test-set, this split needs to be 'cut off', such that the node before that split is a new terminal node then
   - The seperate fold-wise predicitons are averaged in a weighted/ unweighted way for the final predicitons

#### ! ! ! Closer Information to approaches above are in the MS-Thesis itself! ! !  


The approaches from Hagenberg's master thesis are listed below - closer Information to these approaches in the MS-Thesis of Hagenberg.
- **mdd-sPLS:**
   - Method from Lorenzo et al. that can directly deal with block-wise missingness
- **piority-Lasso**
   - Hagenberg's adaption of the original priority-Lasso method to deal with block-wise missingness

---

## Data
Two different data sources are used for the comparison of the differnt approaches for data with block-wise missingness.

### TCGA
14 real multi-omics data sets, where each of these data sets contains the measurements of patients with a certain cancer type. The data was not directly accessed via TCGA, but provided by R. Hornung who has used the data in one of his articles already.  

The original TCGA data is not part of this Repo. If interested in the original data send an E-Mail to 'f.ludwigs@yahoo.de'.  
Only the subsetted TCGA data can be found in the repository under: "data/processed/TCGA_subset_12345" 

### Clinical asthma data
It is a real world data set with block-wise missingness that was provided by the group of Prof. Dr. med. Bianca Schaub at the 'paediatric clinic Dr. von Haunersches Kinderspital'. The data was collected as part of a clinical case-control study in the field of asthma research, whereby the target variable of the data is binary and defined as the presence of asthma.  

For data protection reasons, the data can not be stored in the repository nor be shared.  


#### ! ! ! Closer Information to the data above are in the MS-Thesis itself! ! !  

---

## Code
Short describtion of the scripts in './code'!  
The code scripts either refer to the 'TCGA' data, the 'real' data or is 'general' script needed for both data sources!

#### General
``` 
- GENERAL_DecisionTreeExample:
    Get an example figure for the splitting of the feature-space of a single decision tree

- GENERAL_simpleRF_adaption:
    Implemention of a random forest class that has the option to dynamically prune the single     
    decision trees of a RF. This is needed for the implementation of the 'fold-wise' approach.  
    Whole code builds up on 'github.com/mnwright/simpleRF'
``` 

#### TCGA
``` 
- TCGA_01_data_overview: 
    Get a overview of the different DFs in TCGA & get some additional Information!

- TCGA_02_explorative_performance:
    Check the performance of a RF when trained on a single feature-block/ on all   
    joint feature blocks. Then create different subsets of the original data & get the 
    the predictive performance on these subsets

- TCGA_03_subset_DFs:
    Create final subsets of the original TCGA datasets & get the predictive performance
    on the final subsets (single-block & joint-block). Plot the results and get dimensions
    of the subsetted DFs

- TCGA_04_TestTrain_splits_on_subsetted_DFs:
    Split the subsetted DFs into Test-Train splits for the 5-fold CV. The training data is
    induced with different patterns of block-wise missingness. The resulting Test-Train splits 
    can then be used for the CV of the different approaches.  

- TCGA_05_Foldwise_RF_CV:
    CrossValidate the fold-wise Approach for all different settings & all possible combinations 
    of block-wise missingness in the test-set

- TCGA_06_BlockWise_RF_CV:
    CrossValidate the block-wise Approach for all different settings & all possible combinations 
    of block-wise missingness in the test-set

- TCGA_07_CompleteCases_CV:
    CrossValidate the complete-case Approach for all different settings & all possible combinations 
    of block-wise missingness in the test-set

- TCGA_08_Imputation_CV:
    CrossValidate the Imputation Approach for all different settings & all possible combinations 
    of block-wise missingness in the test-set

- TCGA_09_SingleBlock_CV:
   CrossValidate the Single-Block Approach for all different settings & all possible combinations 
    of block-wise missingness in the test-set

- TCGA_10_Plot_CV_Results:
    Vizualize the Results from the CV of the TCGA data for all the different approaches!
``` 

#### Clinical Asthma Data
``` 
- REAL_01_Fold-Wise-Approach: 
    CrossValidate the fold-wise Approach on the clinical asthma data 

- REAL_02_Block-Wise-Approach:
    CrossValidate the block-wise Approach on the clinical asthma data   

- REAL_03_Imputation-Approach:
   CrossValidate the Imputation Approach on the clinical asthma data 

- REAL_04_CompleteCase-Approach:
    CrossValidate the Complete-Case Approach on the clinical asthma data 

- REAL_05_Single-Block-Approach:
    CrossValidate the Single-Block Approach on the clinical asthma data 

- REAL_06_Plot_CV_Results:
    Vizualize the Results from the CV of the clinical asthma data for all the different approaches!
``` 


All scripts were written in the 'R' language - version 3.6.2.  
The used packages and their corresponding versions are listed below:  

- Attached base packages:  
``` 
[1] grid      parallel  stats     graphics  grDevices utils     datasets 
[8] methods   base     
``` 

- Other attached packages:  
``` 
 [1] rattle_5.3.0          rpart.plot_3.0.8      rpart_4.1-15         
 [4] missForest_1.4        itertools_0.1-3       ROCR_1.0-7           
 [7] gplots_3.0.1.1        e1071_1.7-3           doParallel_1.0.15    
[10] iterators_1.0.12      foreach_1.4.7         assertthat_0.2.1     
[13] pROC_1.15.3           mlbench_2.1-1         reshape2_1.4.3       
[16] gridExtra_2.3         checkmate_1.9.4       caret_6.0-84         
[19] ggplot2_3.2.1         lattice_0.20-38       randomForestSRC_2.9.2
[22] randomForest_4.6-14  
``` 
- Loaded via a namespace (and not attached):  
``` 
 [1] Rcpp_1.0.3         lubridate_1.7.4    class_7.3-15       gtools_3.8.1      
 [5] ipred_0.9-9        R6_2.4.1           plyr_1.8.5         backports_1.1.5   
 [9] stats4_3.6.2       pillar_1.4.3       rlang_0.4.2        lazyeval_0.2.2    
[13] rstudioapi_0.10    data.table_1.12.8  gdata_2.18.0       Matrix_1.2-18     
[17] splines_3.6.2      gower_0.2.1        stringr_1.4.0      munsell_0.5.0     
[21] compiler_3.6.2     pkgconfig_2.0.3    nnet_7.3-12        tidyselect_0.2.5  
[25] tibble_2.1.3       prodlim_2019.11.13 codetools_0.2-16   crayon_1.3.4      
[29] dplyr_0.8.5        withr_2.1.2        MASS_7.3-51.4      bitops_1.0-6      
[33] recipes_0.1.8      ModelMetrics_1.2.2 nlme_3.1-142       gtable_0.3.0      
[37] lifecycle_0.1.0    magrittr_1.5       scales_1.1.0       KernSmooth_2.23-16
[41] stringi_1.4.3      timeDate_3043.102  generics_0.0.2     lava_1.6.6        
[45] tools_3.6.2        glue_1.4.0         purrr_0.3.3        survival_3.1-8    
[49] colorspace_1.4-1   caTools_1.17.1.3  
``` 

---

## Project Organization
```
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
    │   │        ├─ missingness_1234 <- Subsetted data split to test and train set, whereby
    │   │        │                      the train-set contains blockwise missingness 
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
    │   ├── Hagenberg      <- Part of the thesis' of Hagenberg + describtion
    │   │                     of clinical asthma data from Hagenberg
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