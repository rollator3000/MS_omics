# A comparison study of prediction approaches for multiple training data sets & test data with block-wise missing values

This is the README to the repository of Frederik Ludwigs' Master-Thesis, supervised by:
*Dr. Hornung @ Ludwig-Maximilians University - IBE*

## Project description
In this project different methods to deal with blockwise missing data are compared! Blockwise missing data describes data, where for a set of observations whole feature-sets are missing! So the overlap of measurements between different observations is limited! In these settings incorporating data structures, where different features are measured are observed for different observations, could be beneficial!
**Example:**
Different hospitals do reseach regarding the same response e.g. 'BreastCancer'. For this the different hospitals collect data, but the features, that were collected, do not necessarily have to be the same!
``` 
- Hospital_1: Clinical + RNA Data
- Hospital_2: Clinical + CopyNumberVariation Data
- Hospital_3: Clinical + CopyNumberVariation + Mutation Data
```

For Details regarinding methods, data situations etc. please have a look at: 
    --> the MasterThesis itself
    --> the Repository

## Project Organization
------------
    ├── LICENSE
    ├── .git               <- Folder for VersionControl [created by GIT]
    ├── README.md          <- Top-level 'README' for developers working with
    │                         this repository
    │
    ├── code               <- Folder w/ all the R-Code needed to run this 
    │   │                     repository!
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