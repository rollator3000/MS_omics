Adapting the RandomForrest-Algorithm for the case with blockwise missing data!
==============================
This is the README for the Masterthesis of Frederik Ludwigs, supervised by: </p>
Dr. Hornung @ Ludwig-Maximilians University

## Project description
Topic of the Masterthesis is the comparison of methods that can deal with blockwise missing data!
The Algorithms will be applied to clinical data, with clinical features and 4-5 blocks of OmicsData
The OmicsBlocks are extremly highdimensional.
We apply different Methods, to deal with the blockwise missingness. 
For further details on the methods, the results and the data itself please have a look at the thesis!

## Project Organization
------------
    ├── LICENSE
    ├── .git               <- Folder for VersionControl [created by GIT]
    ├── README.md          <- Top-level 'README' for developers working with
    │                         this repository
    │
    ├── code               <- Folder that holds all the code, that was used in
    │   |                     is project!
    │   │     
    │   └── examples       <- Code that is not directly used for the project,
    │                         but served as template, for different situations
    │  
    ├── data               <- All data used in this project
    │   │   
    │   ├── raw            <- The original, immutable data dump.
    │   │
    │   ├── external       <- Data from third party sources
    │   │    | 
    │   │    └─ Dr_Hornung <- Clinical Omics-Data from Dr. Hornung in raw & 
    │   │                     preprocessed form! Also holds the code to down-
    │   │                     load the raw data & preprocess it!
    │   │
    │   └── processed      <- Final processed Data, used for modelling, 
    │                         performance Evaluation,...
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
    └── docs               <- Documents used within this repository! Explainatory material
                              Describtive material, Notes taken, etc.


--------