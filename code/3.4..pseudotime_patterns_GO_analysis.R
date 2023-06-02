###############################################################
# 3.4. pseudotime patterns - GO analysis
# 
# This R script takes the response genes in each pseudotime pattern (Jupyter Notebook 3.3)
# and analyses the GO terms that are associated with those genes. The package that
# we will use for GO analysis is called topGO and detailed instructions can be 
# found here: https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf.
# 
# Instead of showing the GO terms for all of the patterns in RStudio, the GO terms are directly 
# saved into an excel file in the results folder. If you don't want to run this code, but
# just have a look at the calculated GO terms, you can find the file in the Github repository
# (results/pseudotime/pseudotime_patterns_GO_terms.xlsx).

###############################################################
### load libraries and data

# set working directory to "code" folder
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# load libraries
library(topGO)
library("xlsx")

# load patterns
patterns_path <- "../results/pseudotime/pseudotime_pattern_clusters_response_genes.csv"
patterns <- read.csv(patterns_path, sep = ",", header = TRUE, stringsAsFactors = FALSE)
patterns[patterns==""] <- NA #remove empty strings

# get vector with all genes
all_genes <- unique(unlist(patterns))
all_genes <- na.omit(all_genes)

###############################################################
### get GO terms for each pattern

# create a dummy table with all genes
dummy_table <- apply(patterns, 2, function(x){as.integer(all_genes %in% x)})
dummy_table <- as.data.frame(dummy_table)

# specify path were results will be stored
save_path <- "../results/pseudotime/pseudotime_patterns_GO_terms.xlsx"

# calculate and save GO terms per pattern
for (i in 1:length(dummy_table))
{
  test_genes <- as.factor(dummy_table[ , i])
  names(test_genes) <- all_genes
  
  go_data <- new("topGOdata", 
                 ontology = "BP", 
                 allGenes = test_genes, 
                 annotationFun = annFUN.org,
                 mapping = "org.Mm.eg.db",
                 ID = "symbol")
  
  resultFisher <- runTest(go_data, algorithm = "classic", statistic = "fisher")
  print(i)
  goEnrichment <- GenTable(go_data, KS=resultFisher, orderBy="KS", topNodes=40, numChar=1000)
  goEnrichment$KS <- as.numeric(goEnrichment$KS)
  goEnrichment <- goEnrichment[goEnrichment$KS<0.05,]
  
  # get genes significant in each GO term
  goEnrichment$genes <- sapply(goEnrichment$GO.ID, function(x)
  {
    genes <- genesInTerm(go_data, x) 
    genes <- genes[[1]][genes[[1]] %in% patterns[i][[1]]]
    genes <- paste(genes, collapse=",")
  })
  
  if(dim(goEnrichment)[1] != 0)
    if(i == 1)
    {
      write.xlsx(goEnrichment, file = save_path, sheetName = paste0("pattern_", i), append = FALSE)
    }
    else
    {
      write.xlsx(goEnrichment, file = save_path, sheetName = paste0("pattern_", i), append = TRUE)
    }
  
  remove(test_genes, go_data, resultFisher, goEnrichment)
  
}


