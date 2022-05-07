# IRRG:a method for constructing cellular communication networks by integrating receptor-regulated gene expression information
*by Shuqi Guo (guoshuqi9805@gmail.com),Shaowu Zhang\* (zhangsw@nwpu.edu.cn),Yan Li (linay0124@outlook.com), Shihua Zhang (zsh@amss.ac.cn)*

## Introduction
In this repository, We provide an R package for integrating receptor-regulated gene expression information on the basis of ligand-receptor interactions to construct cellular communication networks between cell typese.

## Directory Tree
```
├─Data                  [The database used by IRRG and the single cell data from the case study]
│  │  LRdb.rda                  [Ligand-Receptor Database]
│  │  mm2Hs.rda                 [List of mouse and human homologous gene interconversions]
│  │  PwC_ReactomeKEGG.rda      [Signaling Pathway Database]
│  │  
│  ├─Mouse_IFE                  [Mouse interfollicular epidermal scRNA-seq dataset]
│  │      cell_ann.csv              [Cell type annotations of IFE]
│  │      IFE_express.csv           [Gene expression matrix of IFE]
│  │      
│  └─RCC                        [Renal cell carcinoma scRNA-seq dataset]
│         SI_22368_tumor.csv        [Gene expression matrix of tumor tissues]
│         SI_22369_normal.csv       [Gene expression matrix of normal tissues]
│              
└─Program               [The code of IRRG framework]
        cal_Rscore.cpp          [The function uses a random walk algorithm to iteratively update the receptor gene regulatory network]
        data_process.R          [Script for simple pre-processing of gene expression matrix]
        IFE_Demo.R              [Example of building a cellular communication network for IFE]   
        intra_network.R         [The function used to build the receptor gene regulatory network] 
        LRscore.R               [The function used to calculate the ligand-receptor interaction score]
        net2adj.R               [This function converts the network into an adjacency matrix]
        R_info.R                [The function used to calculate the receptor score]
        separate_data.R         [This function groups gene expression matrices according to whether the receptor is expressed or not]
        simplify_interactions.R [This function simplifies the types of gene interactions in signaling pathways]
```

## Computational flow of IRRG algorithm
We provide an example of IRRG applied to the IFE dataset [`Program/IFE_Demo.R`](https://github.com/NWPU-903PR/IRRG/tree/master/Program/IFE_Demo.R), which consists of the following steps.

**1.Data pre-processing**

Here we first use a simple R script[`Program/data_process.R`](https://github.com/NWPU-903PR/IRRG/tree/master/Program/data_process.R) to normalize the gene expression matrix for pre-processing, the rows of the expression matrix represent genes and the columns represent cells.

**2.Calculation of receptor scores**

The receptor scores in each cell type can be calculated using the following function，and the results will be saved in the *"R_info"* folder under the working path.
```R
R_info(data = data,cluster = cluster,c.names = c.names,cell.prop = 0.2,
       LRdb=LRdb,mm2Hs = mm2Hs,PWC = PwC_ReactomeKEGG,species = "mus musculus")
```
This function consists of two main steps: 1) call function [`Program/intra_network.R`](https://github.com/NWPU-903PR/IRRG/tree/master/Program/intra_network.R) to build a gene regulatory network for each receptor in each cell type; 2) use function [`Program/cal_Rscore.cpp`] to calculate the stability value of the receptor nodes in the network in order to normalize for the receptor score.

**3.Calculation of ligand-receptor interaction scores**

The ligand-receptor co-expression was multiplied by the receptor score using the following function[`Program/LRscore.R`](https://github.com/NWPU-903PR/IRRG/tree/master/Program/LRscore.R) to obtain all ligand-receptor interaction scores between each cell type.
```R
LRscore(data = data,c.names = c.names,species = "mus musculus",celltype = celltype,Permutation.test = TRUE)
```
The results will be saved in the *"LR_score_notest"* folder in the current path.

By setting the parameter "Permutation.test = TRUE", we can screen for significantly specific ligand-receptor pairs between cell types via the permutation test,and the results will be saved in the *"LR_score_test"* folder of the current path.

**4.Cellular communication network construction**

After the calculation in step 3, we will find a summary of the intercellular communication strength and the number of ligand-receptor pairs in the *"Celltype_communication_summary"* folder of the current pathway