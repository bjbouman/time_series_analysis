###############################################################
# 3.1 calculate pseudotime
# 
# This R script is used to calculate the pseudotemporal ordering for each
# cell in the HSC dataset. In order to calculate the pseudotime we will use
# the response genes that we found in R script 2.1. In case you don't want to
# run that script, you can also find the result on the Github repository 
# ("results/response_genes_all.csv"). To calculate pseudotime we will use a
# least square solution to find the transformation matrix W, that best 
# recovers the experimental timepoints. 

# Note: in this script we are using a L2 normalized and scaled version of
# the time series dataset. This object can be made using notebook 1.7. 

# Note 2: this script will only work if the number of response genes is
# lower than the number of cells.

###############################################################

### set working directory to "code" folder
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

### load libraries
library(Seurat)
library(SeuratDisk)
library(ggplot2)

### load dataset
adata_path <- "../data/count_matrices/all_filtered_L2-normalized_scaled.h5seurat"
adata <- LoadH5Seurat(adata_path)

### load response genes
response_genes_path <- "../results/response_genes/response_genes_all.csv"
response_genes <- read.csv(response_genes_path, header=TRUE)
response_genes <- unique(response_genes$gene)
adata <- adata[response_genes,]

### make subsets for each timepoint
HSCcontrol <- subset(x = adata, subset = (time=="control"))
HSC3h <- subset(x = adata, subset = (time=="3h"))
HSC24h <- subset(x = adata, subset = (time=="24h"))
HSC72h <- subset(x = adata, subset = (time=="72h"))

###############################################################
### calculate pseudotime

X <- adata@assays[["RNA"]]@counts
X <- as.matrix(X)

# get experimental timepoint for each cell
Tvect=c(rep(0,ncol(HSCcontrol)),rep(1,ncol(HSC3h)),rep(2,ncol(HSC24h)),rep(3,ncol(HSC72h)))

# get covariance matrix
XXt = X %*% t(X)

# find the inverse of XXt
XXtinv <- solve(XXt) 

# find the transformation (W) in genes space
W = XXtinv %*% X %*% Tvect 

# derive pseudotime (PT) from transformation
PT= t(W) %*% X

remove(HSCcontrol, HSC3h, HSC24h, HSC72h)
remove(XXt, X, Tvect)

###############################################################
### save calculated matrices as .csv files

# save transformation (W)
transformation_save_path <- "../results/pseudotime/transformation_response_genes.csv"
write.csv(W, transformation_save_path, row.names = T)

# save pseudotime
pseudotime_save_path <- "../results/pseudotime/pseudotime_all_cells.csv"
PT <- t(PT) #cell IDs as row names
write.csv(PT, pseudotime_save_path, row.names = T)
