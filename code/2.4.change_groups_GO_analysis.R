###############################################################
# 2.4 change groups - GO analysis
# 
# This R script takes the groups of response genes that were formed in Jupyter Notebook 2.3
# (based on their change scores in different clusters) and analyses the GO terms that are 
# associated with those DEGs. The package that we will use for GO analysis is called topGO 
# and detailed instructions can be found here: 
# https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf.
# 
# Instead of showing the GO terms for all of the groups in RStudio, the GO terms are directly 
# saved into an excel file in the results folder. If you don't want to run this code, but
# just have a look at the calculated GO terms, you can find the file in the Github repository
# (results/change_score/change_score_groups_GO_terms.xlsx).

###############################################################
### load libraries and data

# set working directory to "code" folder
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# load libraries
library(topGO)
library("xlsx")
new('topGOdata')

# load groups
groups_path <- "../results/change_scores/change_groups_response_genes.csv"
groups <- read.csv(groups_path, sep = ",", header = TRUE, stringsAsFactors = FALSE)
groups[groups==""] <- NA #remove empty strings

# get vector with all genes
all_genes <- unique(unlist(groups))
all_genes <- na.omit(all_genes)

###############################################################
### get DEGs for each timepoint

# create a dummy table with all genes
dummy_table <- apply(groups, 2, function(x){as.integer(all_genes %in% x)})
dummy_table <- as.data.frame(dummy_table)

# specify path were results will be stored
save_path <- "../results/change_scores/change_groups_GO_terms.xlsx"

# calculate and save GO terms per DEG group
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
  goEnrichment <- GenTable(go_data, KS=resultFisher, orderBy="KS", topNodes=40, numChar=1000)
  goEnrichment$KS <- as.numeric(goEnrichment$KS)
  goEnrichment <- goEnrichment[goEnrichment$KS<0.05,]
  
  # get genes significant in each GO term
  goEnrichment$genes <- sapply(goEnrichment$GO.ID, function(x)
  {
    genes <- genesInTerm(go_data, x) 
    genes <- genes[[1]][genes[[1]] %in% groups[i][[1]]]
    genes <- paste(genes, collapse=",")
  })
  
  if(i == 1)
  {
    write.xlsx(goEnrichment, file = save_path, sheetName = paste0("group_", i), append = FALSE)
  }
  else
  {
    write.xlsx(goEnrichment, file = save_path, sheetName = paste0("group_", i), append = TRUE)
  }
}
