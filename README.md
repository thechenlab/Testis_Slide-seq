# Testis_Slide-seq

This repositaory contains custom code for analyzing testis Slide-seq data.

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

For information on NMFreg, please go to https://github.com/tudaga/NMFreg_tutorial.
