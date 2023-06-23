###############################################################
# 5.2. abundance analysis
# 
# In this R script we will use the package Milo to perform an abundance analysis 
# on the HSPC dataset. Milo is a package published by Dann et al. (Nature 
# Biotechnology, 2021) A more detailed description of Milo can be found here:
# https://github.com/MarioniLab/miloR (Github page) and here: 
# https://doi.org/10.1038/s41587-021-01033-z (publication). In this script
# we will mostly follow the steps as in the milo tutorial on the mouse
# gastrulation dataset.
# 
# In order to use this script you will need the filtered AND L2-normalized
# version of the  HSPC dataset that can be produced using Notebook 1.3 in 
# combination with Notebook 1.7. Additionaly, you will need the Scanorama-reduced 
# dimension, which were calculated in Notebook 1.3 too. These can also be found
# on the Github repository under "results/processing/scanorama_reduced_dimensions.csv".

###############################################################

### set working directory to "code" folder
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

### load libraries
library(Seurat) 
library(SeuratDisk)
library(scater)
library(SingleCellExperiment)
library(gridExtra) #needed for plotting plots next to each other
library(cowplot) #also used for arranging plots
library(biomaRt) #used for extracting ribo genes
library(miloR)
library(scran)
library(dplyr) 
library(patchwork)
library(batchelor)

### load dataset
adata_path <- "../data/count_matrices/all_filtered_L2-normalized.h5seurat"
adata <- LoadH5Seurat(adata_path)

### load HVGs (all datasets)
hvgs_path <- "../results/cluster_annotation/HVGs_top2000_all_timepoints_merged.csv"
#hvgs_path <- "../results/processing/HVGs_top500_HSC-PBS.csv"
hvgs <- read.csv(hvgs_path, header=TRUE)
hvgs <- hvgs$hvgs #make vector from dataframe

### select two timepoints and subset seurat object
seurat_pbs <- subset(x = adata, subset = (time == "PBS"))
seurat_3h <- subset(x = adata, subset = (time == "3h"))
seurat_24h <- subset(x = adata, subset = (time == "24h"))
seurat_72h <- subset(x = adata, subset = (time == "72h"))

### convert seurat subsets to sce (milo works on sce objects)
sce_pbs <- as.SingleCellExperiment(seurat_pbs)
sce_3h <- as.SingleCellExperiment(seurat_3h)
sce_24h <- as.SingleCellExperiment(seurat_24h)
sce_72h <- as.SingleCellExperiment(seurat_72h)

### save order of cell types
order_clusters <- levels(adata@meta.data[["clusters"]])

### remove unwanted object
remove(seurat_pbs, seurat_3h, seurat_24h, seurat_72h, adata_path, hvgs_path)

####################################################################################
### batch correction using fastMNN 
# In order to use milo on our dataset we will need our data in a batch-correct,
# dimension-reduced space. We will use the Scanorama-reduced dimensions that were
# calculated earlier in Notebook 1.3. 

### create single-cell experiment object with complete dataset
sce_object <- as.SingleCellExperiment(adata)

### add dimensions to single-cell experiment object
scanorama_dim_path <- "../results/processing/scanorama_reduced_dimensions.csv"
scanorama_dim <- read.csv(scanorama_dim_path, sep = ",", header = FALSE, stringsAsFactors = FALSE, row.names = 1)
scanorama_dim <- scanorama_dim[sce_object@colData@rownames,]
scanorama_dim <- data.matrix(scanorama_dim)
reducedDim(sce_object, "pca.corrected") <- scanorama_dim

### add sample (combination of hashtag and pbs) to sce object
sce_object@colData@listData[["sample"]] <- paste0(sce_object@colData@listData[["time"]], "_", sce_object@colData@listData[["hashtags"]])

### remove unwanted objects
remove(f.out, sce_3h, sce_pbs, sce_24h, sce_72h, adata, hvgs)

####################################################################################
### Use milo to perform abundance analysis
# Here, we use milo to perform an abundance analysis between PBS and every other
# timepoint. In order to do so we first calculate the neighbourhoods using all
# timepoints. Then we focus on PBS and one other timepoint. This is slightly
# different from the milo tutorial, where they are working with only two datasets.

### create milo object
milo_object <- Milo(sce_object)

### parameters for NN graphs
k_select = 30
d_select = 30
prop_select = 0.2

#k_select = 30
#d_select = 10
#prop_select = 0.3

### create neighbourhoods
milo_object <- buildGraph(milo_object, k = k_select, d = d_select, reduced.dim = "pca.corrected") #create kNN graph
milo_object <- makeNhoods(milo_object, prop = prop_select, k = k_select, d=d_select, refined = TRUE, reduced_dims = "pca.corrected") #assign the neighbourhoods
plotNhoodSizeHist(milo_object) #create plot of neighbourhood sizes
milo_object <- countCells(milo_object, meta.data = as.data.frame(colData(milo_object)), sample="sample")

### perform abundance analysis (on neighhbourhoods) per timepoint
timepoints <- c("3h", "24h", "72h")

for (timepoint in timepoints)
{
  print(paste0("Calculating DA for PBS vs ", timepoint))
  
  del_timepoints <- timepoints[timepoints != timepoint]
  
  # delete two timepoints to keep just PBS and one other timepoint 
  milo_object@colData@listData[["time_print"]] <- milo_object@colData@listData[["time"]]
  milo_object@colData@listData[["time_print"]][milo_object@colData@listData[["time_print"]] == del_timepoints[1]] <- NA
  milo_object@colData@listData[["time_print"]][milo_object@colData@listData[["time_print"]] == del_timepoints[2]] <- NA
  milo_object@colData@listData[["time_print"]] <- droplevels(milo_object@colData@listData[["time_print"]])
  
  # run the analysis
  study_design <- data.frame(colData(milo_object))[,c("sample", "time_print")] #create table with study design
  study_design <- distinct(study_design)
  study_design <- na.omit(study_design)
  rownames(study_design) <- study_design$sample
  milo_object <- calcNhoodDistance(milo_object, d=d_select, reduced.dim = "pca.corrected")
  da_results <- testNhoods(milo_object, design = ~ time_print, design.df = study_design) #da analysis on neighbourhoods
  
  # name object for use outside loop
  object_name <- paste0("da_results_pbs",timepoint)
  assign(object_name, da_results) #assign name to the object that is created (so it will be a separate variable outside the for-loop)
}

### remove unwanted objects
remove(da_results, study_design)

####################################################################################
### Plotting the results of DA analysis

### Create combined UMAP plots with neighbourhoods
milo_object <- buildNhoodGraph(milo_object)

vmax=6
vmin=-6

plot1 <- plotNhoodGraphDA(milo_object, da_results_pbs3h, layout="UMAP",alpha=0.1) +
  scale_size_continuous(range = c(1, 4)) +
  scale_fill_gradientn(limits = c(vmin,vmax), colours=c("blue","white","red")) +
  #ggtitle("control vs. 3h") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") 
plot2 <- plotNhoodGraphDA(milo_object, da_results_pbs24h, layout="UMAP",alpha=0.1) +
  scale_size_continuous(range = c(1, 4)) +
  scale_fill_gradientn(limits = c(vmin,vmax), colours=c("blue","white","red")) +
  #ggtitle("control vs. 24h") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
plot3 <- plotNhoodGraphDA(milo_object, da_results_pbs72h, layout="UMAP",alpha=0.1) +
  scale_size_continuous(range = c(1, 4)) +
  scale_fill_gradientn(limits = c(vmin,vmax), colours=c("blue","white","red")) +
  #ggtitle("control vs. 72h") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

all_plots <- plot_grid(plot1, plot2, plot3, hjust = -1, nrow = 1)
legend <- get_legend(plot1 + theme(legend.position="bottom"))
plot_grid(all_plots, legend, ncol = 1, rel_heights = c(1, .2))

# save all 3 plots together
#save_path = "../figures/2.2.milo_abundances_UMAPs.png"
#png(filename=save_path, width = 1000, height = 400)
#plot_grid(all_plots, legend, ncol = 1, rel_heights = c(1, .2))
#dev.off()

### save plots individually
ggsave("5.2.milo_UMAP_ctrl-3h.png", plot1, "png", path="../figures/", width=5, height=5, dpi=300)
ggsave("5.2.milo_UMAP_ctrl-24h.png", plot2, "png", path="../figures/", width=5, height=5, dpi=300)
ggsave("5.2.milo_UMAP_ctrl-72h.png", plot3, "png", path="../figures/", width=5, height=5, dpi=300)

### save legend 
legend <- get_legend(plot1 + theme(legend.position='right'))
ggsave("5.2.milo_legend.pdf", legend, "pdf", path="../figures/", width=2, height=6, dpi=300)


### plot beeswarmplots together
da_results_pbs3h <- annotateNhoods(milo_object, da_results_pbs3h, coldata_col = "clusters")
da_results_pbs24h <- annotateNhoods(milo_object, da_results_pbs24h, coldata_col = "clusters")
da_results_pbs72h <- annotateNhoods(milo_object, da_results_pbs72h, coldata_col = "clusters")

da_results_pbs3h$new_clusters <- ifelse(da_results_pbs3h$clusters_fraction < 0.7, "Mixed", da_results_pbs3h$clusters)
da_results_pbs24h$new_clusters <- ifelse(da_results_pbs24h$clusters_fraction < 0.7, "Mixed", da_results_pbs24h$clusters)
da_results_pbs72h$new_clusters <- ifelse(da_results_pbs72h$clusters_fraction < 0.7, "Mixed", da_results_pbs72h$clusters)

# remove mixed neighbourhoods
subset_pbs3h <- da_results_pbs3h[!is.na(da_results_pbs3h$new_clusters),]
subset_pbs24h <- da_results_pbs24h[!is.na(da_results_pbs24h$new_clusters),]
subset_pbs72h <- da_results_pbs72h[!is.na(da_results_pbs72h$new_clusters),]

subset_pbs3h$new_clusters <- factor(subset_pbs3h$new_clusters, levels = rev(order_clusters)) 
subset_pbs24h$new_clusters <- factor(subset_pbs24h$new_clusters, levels = rev(order_clusters))
subset_pbs72h$new_clusters <- factor(subset_pbs72h$new_clusters, levels = rev(order_clusters))

plot1 <- plotDAbeeswarm(subset_pbs3h, group.by = "new_clusters") +
  scale_color_gradientn(limits = c(vmin,vmax), colours=c("blue","white","red")) +
  scale_y_continuous(limits=c(vmin, vmax))
plot2 <- plotDAbeeswarm(subset_pbs24h, group.by = "new_clusters") +
  scale_color_gradientn(limits = c(vmin,vmax), colours=c("blue","white","red"))+
  scale_y_continuous(limits=c(vmin, vmax))
plot3 <- plotDAbeeswarm(subset_pbs72h, group.by = "new_clusters") +
  scale_color_gradientn(limits = c(vmin,vmax), colours=c("blue","white","red"))+
  scale_y_continuous(limits=c(vmin,vmax))

patchwork = plot1 + plot2 + plot3

# Remove title from second subplot
patchwork[[2]] = patchwork[[2]] + theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank() )

# Remove title from third subplot
patchwork[[3]] = patchwork[[3]] + theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank() )
patchwork

# save figure beeswarmplot
ggsave("5.2.milo_beeswarmplots.pdf", patchwork, "pdf", path="../figures/", width=13.5, height=9, dpi=300)

