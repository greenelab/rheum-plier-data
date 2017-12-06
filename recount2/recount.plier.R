# Qiwen Hu - 2017
# PLIER model for recount datasets
# Run PLIER for all recount datasets, genes are normalized by RPKM
# It can be run from the command line using Rscript recount.plier.R
#
# Output: 
# PLIER model for all the samples

library(recount)
library(PLIER)
library('biomaRt')

data(bloodCellMarkersIRISDMAP)
data(svmMarkers)
data(canonicalPathways)
load("recount.rpkm.RData")
set.seed(12345)

# Transform ensemble id to genesymbol
mart <- biomaRt::useDataset("hsapiens_gene_ensembl", biomaRt::useMart("ensembl"))
genes <- unlist(lapply(strsplit(rpkm$id, "[.]"), `[[`, 1))
rpkm$ensembl_gene_id <- unlist(lapply(strsplit(rpkm$id, "[.]"), `[[`, 1))
gene.list <- biomaRt::getBM(filters = "ensembl_gene_id",
		attributes = c("ensembl_gene_id", "hgnc_symbol"),
		values = genes, mart = mart)
gene.list <- gene.list[gene.list$hgnc_symbol != "", ]

rpkm <- merge(rpkm, gene.list, by = c("ensembl_gene_id"))
rownames(rpkm) <- make.names(rpkm$hgnc_symbol, unique = TRUE)

# Remove redundant hgnc_symbol here
rpkm <- rpkm[, -ncol(rpkm)]

# Remove redundant ensembl_gene_id and gene_id
rpkm <- rpkm[, -1*c(1, 2)]

# Run PLIER
allPaths <- combinePaths(bloodCellMarkersIRISDMAP, svmMarkers,canonicalPathways)
cm.genes <- commonRows(allPaths, rpkm)
rpkm <- rowNorm(rpkm)
k <- num.pc(rpkm[cm.genes, ])

rpkm <- as.matrix(as.data.frame(rpkm))
plierResult <- PLIER(rpkm[cm.genes, ], allPaths[cm.genes, ],
		     k = round((k + k*0.3), 0), trace = TRUE)
save(plierResult, file = "recount.all.data.plier.RData")
