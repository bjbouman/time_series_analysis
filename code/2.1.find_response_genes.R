###############################################################
# 2.1 find response genes
# 
# This R script uses the Libra package to find the differentially
# expressed genes (DEGs) in the HSPC dataset. For each cluster
# we find the DEGs between PBS and one other timepoint. In the
# end we take the DEGs from all clusters with the highest
# significance to use for downstream analyses. We will refer to
# to those genes as the response genes.
# 
# In order to run this code you will need the filtered dataset with 
# all timepoint, which we created in Notebook 1.3.

###############################################################

### set working directory to "code" folder
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

### load libraries
library(limma)
library(edgeR)
library(Seurat)
library(SeuratDisk)
library(Libra)
library(ggplot2)

### load dataset
adata_path <- "../data/count_matrices/all_filtered.h5seurat"
adata <- LoadH5Seurat(adata_path)

###############################################################
### get DEGs for each timepoint
# Here we are selecting DEGs between PBS and one other timepoint. 
# We do this using the edgeR method, which has been included in the
# Libra package. Libra contains many other DEG-finding methods, so
# feel free to play around with different methods too. The four 
# hashtags in each timepoint are used to separate the replicates.
# Therefore, we use the hashtags as entry for "replicate_col".

# (edgeR: https://doi.org/10.1093/bioinformatics/btp616)
# (Libra: https://doi.org/10.1038/s41467-021-25960-2)
# (Libra tutorial: https://github.com/neurorestore/Libra)

### get timepoints
timepoints <- levels(adata@meta.data[["time"]])

all_degs <- data.frame() #create empty dataframe

for (timepoint in timepoints[2:4])
{
  ### subset adata to contain only the entries from PBS and one other timepoint
  adata_subset <- subset(x = adata, subset = ((time==timepoints[1])|(time==timepoint)))
  adata_subset@meta.data[["time"]] <- droplevels(adata_subset@meta.data[["time"]])
  
  ### run DEG analysis
  degs_timepoint <- run_de(adata_subset, meta = adata_subset@meta.data,
                           replicate_col = "hashtags",
                           cell_type_col = "clusters",
                           label_col = "time", de_method = "edgeR", de_type = "LRT")
  
  ### add timepoint to DEGs
  degs_timepoint$time <- timepoint
  
  ### keep DEGs above log fold change and adjusted p-value threshold
  degs_timepoint <- degs_timepoint[(abs(degs_timepoint$avg_logFC)>=1),]
  degs_timepoint <- degs_timepoint[(degs_timepoint$p_val_adj <= 0.05),]
  
  ### add found DEGs for this timepoint to dataframe with all DEGs
  all_degs <- rbind(all_degs, degs_timepoint)
}

###############################################################
### save files with response DEGs

### print total number of unique DEGs found 
# (some DEGs will be found in multiple clusters or timepoints)
print(length(unique(all_degs$gene)))

### sort DEGs by adjusted p-value (lowest p-value first)
all_degs <- all_degs[order(all_degs$p_val_adj),]

### save complete DEG table
all_degs_save_path <- "../results/response_genes/response_genes_all.csv"
write.csv(all_degs, all_degs_save_path, row.names = F)

### select top x DEGs with the highest adjusted p-value
top_degs <- all_degs[!duplicated(all_degs[ , c("gene")]), ] #keep only one entry of each DEG 
top_x <- 500
top_degs <- top_degs$gene[1:top_x]
top_degs <- data.frame(top_degs)

### save list of top x DEGs
top_degs_save_path <- paste0("../results/response_genes/response_genes_top", top_x, ".csv")
write.csv(top_degs, top_degs_save_path, row.names = F)
