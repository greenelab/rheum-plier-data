# Qiwen Hu - 2017
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

library(recount)
library(PLIER)
library(biomaRt)

data(bloodCellMarkersIRISDMAP)
data(svmMarkers)
data(canonicalPathways)
rpkm <- readRDS(file.path("recount2", "recount_rpkm.RDS"))
set.seed(12345)

# Transform ensemble id to genesymbol
mart <- biomaRt::useDataset("hsapiens_gene_ensembl", 
                            biomaRt::useMart("ensembl"))
genes <- unlist(lapply(strsplit(rpkm$ENSG, "[.]"), `[[`, 1))
rpkm$ensembl_gene_id <- unlist(lapply(strsplit(rpkm$ENSG, "[.]"), `[[`, 1))
gene.list <- biomaRt::getBM(filters = "ensembl_gene_id",
		attributes = c("ensembl_gene_id", "hgnc_symbol"),
		values = genes, mart = mart)
gene.list <- gene.list[gene.list$hgnc_symbol != "", ]

rpkm <- merge(rpkm, gene.list, by = c("ensembl_gene_id"))
rownames(rpkm) <- make.names(rpkm$hgnc_symbol, unique = TRUE)

# Remove redundant hgnc_symbol here
rpkm <- rpkm[, -ncol(rpkm)]

# Remove redundant ensembl_gene_id and gene_id and last column "by"
rpkm <- rpkm[, -1*c(1, 2, ncol(rpkm))]

# PLIER prior information (pathways)
allPaths <- combinePaths(bloodCellMarkersIRISDMAP, svmMarkers,
                         canonicalPathways)
cm.genes <- commonRows(allPaths, rpkm)

# filter to common genes before row normalization to save on computation
rpkm.cm <- rpkm[cm.genes, ]

# remove objects that are not needed
rm(mart,
   gene.list,
   genes,
   rpkm,
   bloodCellMarkersIRISDMAP,
   svmMarkers,
   canonicalPathways)

# row-normalize (z-score)
rpkm.cm <- rowNorm(rpkm.cm)

# determine number of significant PCs
k <- num.pc(rpkm.cm)

# save z-scored expression data, the prior information matrix to be supplied
# to PLIER::PLIER and the number of PCs
plier.data.list <- list("rpkm.cm" = rpkm.cm,
                        "all.paths.cm" = allPaths[cm.genes, ],
                        "k" = k)
saveRDS(plier.data.list, file = file.path("recount2",
                                          "recount_data_prep_PLIER.RDS"))
