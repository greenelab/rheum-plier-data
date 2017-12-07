# Qiwen Hu - 2017
# Processing all recounts datasets
# genes expression are normalized by RPKM
# It can be run from the command line using Rscript get.all.recount.dataset.R
#
# Output: 
# Normalized gene expression for each sample

library(recount)

# Get RPKM value for each gene - adapted from recount package
getRPKM <- function(rse, length_var = 'bp_length', mapped_var = NULL) { 
  # Computes the RPKM value for each gene in the sample.
  #
  # Args: 
  #  rse: A RangedSummarizedExperiment-class object in recount package
  #  length_var: A length 1 character vector with the column name from rowData(rse) that has
  #              the coding length. For gene level objects from recount this is bp_length. If
  #              NULL, then it will use width(rowRanges(rse)) which should be used for exon RSEs.
  #  mapped_var: A length 1 character vector with the column name from colData(rse) that has
  #              the number of reads mapped. If NULL (default) then it will use the column 
  #              sums of the counts matrix
  # Returns:
  #   RPKM value for each sample
  if(!is.null(mapped_var)){
    mapped <- colData(rse)[, mapped_var] 
  } else {
    mapped <- colSums(assays(rse)$counts) 
  } 
  bg <- matrix(mapped, ncol = ncol(rse), nrow = nrow(rse), byrow = TRUE) 
  if(!is.null(length_var)){
    len <- rowData(rse)[, length_var] 
  } else {
    len <- width(rowRanges(rse))
  }
  wid <- matrix(len, nrow = nrow(rse), ncol = ncol(rse), byrow = FALSE) 
  rpkm <- assays(rse)$counts / (wid/1000) / (bg/1e6) 
  return(rpkm)
} 

data.dir <- file.path("recount2", "data")

# Get all samples from recount database
metasample.sra <- all_metadata(subset = "sra", verbose = TRUE)
metasample.sra <- as.data.frame(metasample.sra)

# Remove samples without description
metadata.nonempty <- metasample.sra[!is.na(metasample.sra$characteristics), ] 
included.sample.list <- unique(metadata.nonempty$project)

# Download all recount2 samples in included.sample.list
lapply(included.sample.list, 
       function(x) download_study(x, type = "rse-gene", outdir = file.path(data.dir, x)))

rpkm <- data.frame()

for(i in 1:length(included.sample.list)) {
  rse.gene <- get(load(file.path(data.dir, included.sample.list[i], 'rse_gene.Rdata')))
  if(i == 1) {
    rpkm <- as.data.frame(getRPKM(rse.gene))
    rpkm$id <- rownames(rpkm)
  } else {
    rpkm.tmp <- as.data.frame(getRPKM(rse.gene))
    rpkm.tmp$id <- rownames(rpkm.tmp)
    rpkm <- merge(rpkm, rpkm.tmp, by=c("id"))
  }
}
save(rpkm, file = "recount.rpkm.RData")
