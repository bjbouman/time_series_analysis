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
2. **change scores** Script 2.1 is used to find the response genes. In script 2.2. the expression of some of those response genes over time is visualized. In 2.3 the change score is calculated and each response genes is grouped based on the change score. In 2.4 and 2.5 the gene ontology (GO) terms for each of the change score groups are found and visualized.
3. **pseudotime** Script 3.1 is used to calculate response pseudotime. In 3.2. the pseudotime and pseudotemporal order are visualized in a 3D plot. In 3.3. we explore the dynamics of response genes in pseudotime and categorize the response genes based on their dynamics. In 3.4. and 3.5. the gene ontology (GO) terms for each of the categories of responses are found and visualized. Scripts 3.6., 3.7. and 3.8 explores the weights of the response genes and especially which genes get the highest weights. 
4. **downstream analyses** Script 4.1 creates a heatmap that combines the change score groups and the pseudotime patterns of the response genes. Script 4.3 plots the expression of a few interesting genes in pseudotime.
5. **abundance analysis** Script 5.1 plots the abundance of each cluster as a percentage. Script 5.2 uses a more advanced approach for abundance analysis, relying on neighborhoods. Script 5.3 is used to plot the expression of Ly6a/Sca-1 in the different experimental timepoints.
6. XXX

A detailed description of each individual script can be found in the header of the script. 

### Packages used
This work was build on the shoulders of giants and I would like to give credit where credit is due. So here a list of R/Python packages that have been used in this project and you should definitely check out:

- [Scanpy](https://scanpy.readthedocs.io/en/stable/): a scalable toolkit for analyzing single-cell gene expression data.
- [Seurat](https://satijalab.org/seurat/): an R package designed for QC, analysis, and exploration of single-cell RNA-seq data.
- [Scanorama](https://github.com/brianhie/scanorama): enables batch-correction and integration of heterogeneous scRNA-seq datasets.
- [Libra](https://github.com/neurorestore/Libra): an R package to perform differential expression on single-cell data.
- [topGO](https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf): a package designed to facilitate semi-automated enrichment analysis for Gene Ontology (GO) terms.
- [Milo](https://marionilab.github.io/miloR/articles/milo_demo.html): a tool for analysis of complex single cell datasets generated from replicated multi-condition experiments, which detects changes in composition between conditions.

### Running our scripts
XXX 
