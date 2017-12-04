library(recount)
library(PLIER)
library('biomaRt')

data(bloodCellMarkersIRISDMAP)
data(svmMarkers)
data(canonicalPathways)
load("recount.rpkm.RData")

#transform ensemble id to genesymbol
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- unlist(lapply(strsplit(rpkm$id, "[.]"), `[[`, 1))
rpkm$ensembl_gene_id <- unlist(lapply(strsplit(rpkm$id, "[.]"), `[[`, 1))
G_list <- getBM(filters= "ensembl_gene_id",
		attributes= c("ensembl_gene_id", "hgnc_symbol"),
		values=genes, mart= mart)
G_list <- G_list[G_list$hgnc_symbol != "",]

rpkm <- merge(rpkm, G_list, by=c("ensembl_gene_id"))
rownames(rpkm) <- make.names(rpkm$hgnc_symbol, unique=TRUE)

#remove redundant information
rpkm <- rpkm[, -1*ncol(rpkm)]
rpkm <- rpkm[, -1*c(1,2)]

#run PLIER###
allPaths=combinePaths(bloodCellMarkersIRISDMAP, svmMarkers,canonicalPathways)
cm.genes=commonRows(allPaths, rpkm)
rpkm <- rowNorm(rpkm)
k <- num.pc(rpkm[cm.genes, ])

rpkm <- as.matrix(as.data.frame(rpkm))
plierResult <- PLIER(rpkm[cm.genes, ], allPaths[cm.genes, ],
		     k=round((k + k*0.3),0), trace=T)
save(plierResult, file = "recount.all.data.plier.RData")

