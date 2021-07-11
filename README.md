# Testis_Slide-seq

### This repository contains custom code for analyzing data reported in the following paper:
     
     Dissecting Mammalian Spermatogenesis Using Spatial Transcriptomics.
     Haiqi Chen, Evan Murray, Anisha Laumas, Jilong Li, Xichen Nie, Jim Hotaling, Jingtao Guo, Bradley R. Cairns, Evan Z. Macosko, C. Yan Cheng, Fei Chen
     bioRxiv 2020.10.17.343335; doi: https://doi.org/10.1101/2020.10.17.343335

For the processed wild type and diabetic mouse Slide-seq datasets, please go to https://www.dropbox.com/s/ygzpj0d0oh67br0/Testis_Slideseq_Data.zip?dl=0.
   
   There are three files for each dataset:
    
    1. MappedDGE is the digital gene expression matrix;
    
    2. Beadlocations is the bead location matrix;
    
    3. bead_maxct_df is the cell type assignment by NMFreg.
        
        In this file, each cell type corresponds to a number:
          
          {‘1’:”Elongating/elongated Spermatid”, 
           
           ‘2’:”Round Spermatid”,
           
           ‘3’:”Myoid Cell”,
           
           ‘4’:”Spermatocyte”, 
           
           ‘5’:”Spermatogonium”,
           
           ‘6’:”Sertoli cell”, 
           
           ‘7’:”Leydig cell”, 
           
           ‘8:”Endothelial cell”, 
           
           ‘9’:”Macrophage"}

### For information on NMFreg, please go to https://github.com/tudaga/NMFreg_tutorial.

### Seminiferous Tubule Assignment Workflow.ipynb 
   
    Convert Slide-seq data to image data and segment individual seminiferous tubules.
   
### SPG_Compartment_Analysis.ipynb

    Identify beads belonging to either the undifferentiated or differentiating SPG neighorbood.
    
### Purity Score Calculation.ipynb

    Calculate the ES purity score for wild type and diabetic seminiferous tubules. 
    
### Pairwise Spatial Contact Frequency Analysis.m 

    Calculate the pairwise spatial contact frequency for wild type and diabetic seminiferous tubules. 

### 3D_Segmentation.cpproj
    
    Perform 3D segmentation of the DAPI image stacks from the targeted in situ RNA sequencing experiment. 
    
### Bleedthrough_subtraction.m

    Subtract the bleedthrough from the Cy3 channel in the Texas Red channel of the image stacks from the targeted in situ RNA sequencing experiment.

### Targeted_inSitu_RNA_Seq_Pipeline_Testis.m

    Perform image registration, peak calling, gene assignment and counting on image stacks from the targeted in situ RNA sequencing experiment. 
    The Scripts folder contain subsidiary code needed to run this pipeline. 
    A bioformats java package (bioformats_package.jar) is also needed, and can be downloaded here:       https://www.dropbox.com/s/656b06xvar6n531/bioformats_package.jar?dl=0. 
    After downloading, place the java package in the bfmatlab subfoler under the Scripts/helpers folder. 
    All code in the Scripts folder need to be added to the same directory as the pipeline code. 
