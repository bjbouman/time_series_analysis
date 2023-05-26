# Analysis of single-cell time series

Welcome! This Github repository is complementary to our paper "_Single-cell time series analysis reveals the dynamics of HSPCs response to inflammation_", which is out as [preprint on Biorxiv](https://www.biorxiv.org/content/10.1101/2023.03.09.531881v1.abstract) All of the figures in the manuscript and the extended figures can be reproduced using the scripts in this folder.

If you have any questions, comments, improvements, feel free to reach out by opening an issue. 

### Overview repository
This repository contains four folders:
- **code**: R scripts, bash scripts and Jupyter Notebooks that were used to analyze the data.
- **data**: all data that has been used as input for the different types of analysis.
- **figures**: this folder is empty, but if you decide to run the scripts the folder will be used to store figures.
- **results**: all non-figure results (like tables, lists of genes, etc.) are stored in this folder. If you don't want to run any of the scripts all results can be downloaded directly from here.

### Description of code

0. **producing count matrices** XXX
1. **processing and data description** The scripts 1.1-1.9 are used to: process the dataset (1.1, 1.3 and 1.7), find the cell types assignments for each cluster in our dataset (1.2.1-1.2.5b) and describe certain properties of the dataset (1.4, 1.5, 1.6, 1.8 and 1.9). The scripts produce the (sub)figures in Figure 1 and Extended Figure 1 of our manuscript.
2. XXX
3. XXX
4. XXX
5. XXX
6. XXX

A detailed description of each individual script can be found in the header of the script. 

### Packages used
This work was build on the shoulders of giants and I would like to give credit where credit is due. So here a list of R/Python packages that have been used in this project and you should definitely check out:

- [Scanpy](https://scanpy.readthedocs.io/en/stable/): a scalable toolkit for analyzing single-cell gene expression data.
- [Seurat](https://satijalab.org/seurat/): an R package designed for QC, analysis, and exploration of single-cell RNA-seq data.
- [Scanorama](https://github.com/brianhie/scanorama): enables batch-correction and integration of heterogeneous scRNA-seq datasets.
- [Libra](https://github.com/neurorestore/Libra): an R package to perform differential expression on single-cell data.
- 

### Running our scripts
XXX 
