###############################################################
# 3.7. weighted genes - GO analysis
# 
# In this R script we will analyze the GO terms from the top weighted 
# genes in the W matrix that we used to find pseudotime (R script 3.1).
# We will analyze both the top 50 most positively weighted genes and 
# the top 50 most negatively weighted genes. The package that
# we will use for GO analysis is called topGO and detailed instructions can be 
# found here: https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf.
# 
# Instead of showing the GO terms for all of the patterns in RStudio, the GO terms are directly 
# saved into an excel file in the results folder. If you don't want to run this code, but
# just have a look at the calculated GO terms, you can find the files in the Github repository.

###############################################################
### load libraries and data

# set working directory to "code" folder
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# load libraries
library(topGO)
library("xlsx")

# load patterns
gene_weights_path <- "../results/pseudotime/transformation_response_genes.csv"
gene_weights <- read.csv(gene_weights_path, sep = ",", header = TRUE, stringsAsFactors = FALSE)
names(gene_weights) <- c("gene", "weight")

# get vector with all genes
all_genes <- gene_weights$gene

###############################################################
### calculate and save GO terms for top 50 most positive genes

# get top 50 most positive genes
W_genes_positive <- gene_weights[order(gene_weights$weight, decreasing = TRUE),]$gene[1:50]

# specify path were results will be stored
save_path <- "../results/pseudotime/pseudotime_top50_positive_weighted_genes_GO_terms.xlsx"

# make dummy table of genes
test_genes <- as.factor(as.integer(all_genes %in% W_genes_positive))
names(test_genes) <- all_genes
  
# run GO analysis
go_data <- new("topGOdata", ontology = "BP", allGenes = test_genes, 
               annotationFun = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol")
  
resultFisher <- runTest(go_data, algorithm = "classic", statistic = "fisher")
goEnrichment <- GenTable(go_data, KS=resultFisher, orderBy="KS", topNodes=40, numChar=1000)
goEnrichment$KS <- as.numeric(goEnrichment$KS)
goEnrichment <- goEnrichment[goEnrichment$KS<0.05,]

# get genes significant in each GO term
goEnrichment$genes <- sapply(goEnrichment$GO.ID, function(x)
{
  genes <- genesInTerm(go_data, x) 
  genes <- genes[[1]][genes[[1]] %in% W_genes_positive]
  genes <- paste(genes, collapse=",")
})

# save GO terms in excel file
write.xlsx(goEnrichment, file = save_path, sheetName = "positive_weighted_genes", append = FALSE)

###############################################################
### calculate and save GO terms for top 50 most negative genes

# get top 50 most negative genes
W_genes_negative <- gene_weights[order(gene_weights$weight),]$gene[1:50]

# specify path were results will be stored
save_path <- "../results/pseudotime/pseudotime_top50_negative_weighted_genes_GO_terms.xlsx"

# make dummy table of genes
test_genes <- as.factor(as.integer(all_genes %in% W_genes_negative))
names(test_genes) <- all_genes

# run GO analysis
go_data <- new("topGOdata", ontology = "BP", allGenes = test_genes, 
               annotationFun = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol")

resultFisher <- runTest(go_data, algorithm = "classic", statistic = "fisher")
goEnrichment <- GenTable(go_data, KS=resultFisher, orderBy="KS", topNodes=40, numChar=1000)
goEnrichment$KS <- as.numeric(goEnrichment$KS)
goEnrichment <- goEnrichment[goEnrichment$KS<0.05,]

# get genes significant in each GO term
goEnrichment$genes <- sapply(goEnrichment$GO.ID, function(x)
{
  genes <- genesInTerm(go_data, x) 
  genes <- genes[[1]][genes[[1]] %in% W_genes_negative]
  genes <- paste(genes, collapse=",")
})

# save GO terms in excel file
write.xlsx(goEnrichment, file = save_path, sheetName = "negative_weighted_genes", append = FALSE)

