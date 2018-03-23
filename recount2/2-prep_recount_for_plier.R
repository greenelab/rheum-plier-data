# Qiwen Hu & Jaclyn N. Taroni 2017-2018
# Data prep for PLIER model for recount datasets:
#   * Annotation with HGNC gene symbol
#   * z-score (PLIER::rowNorm)
#   * Determine k (with PLIER::num.pcs) for use with PLIER::PLIER
# USAGE: Rscript recount2/2-prep_recount_for_plier.R
#
# Output: 
# recount2/recount_data_prep_PLIER.RDS -- a list with the following elements:
# row-normalized RPKM from recount (genes common with pathways, only), 
# prior information matrix (contains bloodCellMarkersIRISDMAP, svmMarkers,
# canonicalPathways and only genes in common with recount RPKM), k the number
# of significant PCs

`%>%` <- dplyr::`%>%`
library(recount)
library(PLIER)
library(biomaRt)

# prior information (pathways) from PLIER
data(bloodCellMarkersIRISDMAP)
data(svmMarkers)
data(canonicalPathways)

# normalized recount2 data
rpkm.df <- readRDS(file.path("recount2", "recount_rpkm.RDS"))

# set seed for reproducibility
set.seed(12345)

# Transform ensembl id to genesymbol
mart <- biomaRt::useDataset("hsapiens_gene_ensembl", 
                            biomaRt::useMart("ensembl"))
genes <- unlist(lapply(strsplit(rpkm.df$ENSG, "[.]"), `[[`, 1))
rpkm.df$ensembl_gene_id <- unlist(lapply(strsplit(rpkm.df$ENSG, "[.]"), 
                                         `[[`, 1))
gene.df <- biomaRt::getBM(filters = "ensembl_gene_id",
                          attributes = c("ensembl_gene_id", "hgnc_symbol"),
                          values = genes, 
                          mart = mart)
# filter to remove genes without a gene symbol
gene.df <- gene.df %>% dplyr::filter(complete.cases(.))
# add gene symbols to expression df
rpkm.df <- dplyr::inner_join(gene.df, rpkm.df, 
                             by = "ensembl_gene_id")
# set symbols as rownames (req'd for PLIER)
rownames(rpkm.df) <- make.names(rpkm.df$hgnc_symbol, unique = TRUE)
# remove gene identifier columns
rpkm.df <- rpkm.df %>% dplyr::select(-c(ensembl_gene_id:ENSG))

# PLIER prior information (pathways)
allPaths <- PLIER::combinePaths(bloodCellMarkersIRISDMAP, svmMarkers,
                                canonicalPathways)
cm.genes <- PLIER::commonRows(allPaths, rpkm.df)

# filter to common genes before row normalization to save on computation
rpkm.cm <- rpkm.df[cm.genes, ]

# remove objects that are not needed
rm(mart,
   gene.list,
   genes,
   rpkm,
   bloodCellMarkersIRISDMAP,
   svmMarkers,
   canonicalPathways)

# row-normalize (z-score)
rpkm.cm <- PLIER::rowNorm(rpkm.cm)

# determine number of significant PCs
k <- PLIER::num.pc(rpkm.cm)

# save z-scored expression data, the prior information matrix to be supplied
# to PLIER::PLIER and the number of PCs
plier.data.list <- list("rpkm.cm" = rpkm.cm,
                        "all.paths.cm" = allPaths[cm.genes, ],
                        "k" = k)
saveRDS(plier.data.list, file = file.path("recount2",
                                          "recount_data_prep_PLIER.RDS"))
