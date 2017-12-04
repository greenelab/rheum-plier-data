library(recount)

#get RPKM value for each gene - function modified from recount package
getRPKM <- function(rse, length_var = 'bp_length', mapped_var = NULL) { 
  mapped <- if(!is.null(mapped_var)) colData(rse)[, mapped_var] 
		else colSums(assays(rse)$counts) 
  bg <- matrix(mapped, ncol = ncol(rse), nrow = nrow(rse), byrow = TRUE) 
  len <- if(!is.null(length_var)) rowData(rse)[, length_var] 
		else width(rowRanges(rse)) 
  wid <- matrix(len, nrow = nrow(rse), ncol = ncol(rse), byrow = FALSE) 
  assays(rse)$counts / (wid/1000) / (bg/1e6) 
} 


dir <- '/home/qhu/Recount'

##get all samples from recount database
metasample.sra <- all_metadata(subset = "sra", verbose = TRUE)
metasample.sra <- as.data.frame(metasample.sra)
metadata.nonempty <- metasample.sra[!is.na(metasample.sra$characteristics), ] 
sample_list <- unique(metadata.nonempty$project)

rpkm <- data.frame()

for(i in 1:length(sample_list)) {
  rse_gene <- get(load(file.path(dir, sample_list[i], 'rse_gene.Rdata')))
  if(i == 1) {
    rpkm <- as.data.frame(getRPKM(rse_gene))
    rpkm$id <- rownames(rpkm)
    } else {
    rpkm.tmp <- as.data.frame(getRPKM(rse_gene))
    rpkm.tmp$id <- rownames(rpkm.tmp)
    rpkm <- merge(rpkm, rpkm.tmp, by=c("id"))
    }
}
save(rpkm,file = "recount.rpkm.RData")
