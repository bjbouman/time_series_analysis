###############################################################
# 1.2.4a load Baccin dataset
# 
# In this R script we will extract the cell types and t-SNE coordinates
# from the Seurat object that is supplementary to the Baccin et al. 2019
# publication (https://doi.org/10.1038/s41556-019-0439-6). In the 
# Supplementary info you can download the Supplementary Data. This folder
# contains the "NicheData10X.rda" file. From this file we will extract the
# cell types and t-SNE coordinates, so we can use those in a Jupyter Notebook.

###############################################################

### set working directory to "code" folder
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

### load libraries
library(Seurat)

### load file
data_path <- "../data/public_datasets/Baccin_2019/NicheData10x.rda"
load(data_path)

###############################################################

### get and save t-SNE coordinates
tsne_coordinates <- NicheData10x@reductions$tsne@cell.embeddings
save_path <- "../data/public_datasets/Baccin_2019/tsne_coordinates.csv"
write.csv(tsne_coordinates, save_path, row.names=T)

### get and save clusters
clusters <- data.frame(NicheData10x@active.ident)
colnames(clusters) <- c("clusters")
save_path <- "../data/public_datasets/Baccin_2019/clusters.csv"
write.csv(clusters, save_path, row.names=T)

