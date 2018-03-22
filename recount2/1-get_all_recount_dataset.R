# Qiwen Hu - 2017
# Processing all recounts datasets
# genes expression are normalized by RPKM
# It can be run from the command line using 
# Rscript recount2/1-get_all_recount_dataset.R
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
dir.create(data.dir, recursive = TRUE, showWarnings = FALSE)

# Get all samples from recount database
metasample.sra <- all_metadata(subset = "sra", verbose = TRUE)
metasample.sra <- as.data.frame(metasample.sra)

# Remove samples without description
metadata.nonempty <- metasample.sra[!is.na(metasample.sra$characteristics), ] 
included.sample.list <- unique(metadata.nonempty$project)

# Download all recount2 samples in included.sample.list
lapply(included.sample.list, 
       function(x) download_study(x, type = "rse-gene", 
                                  outdir = file.path(data.dir, x)))

# get RPKM for each experiment and add to list
rpkm.list <- list()
for(experiment in included.sample.list) {
  load(file.path(data.dir, experiment, 'rse_gene.Rdata'))
  rpkm <- as.data.frame(getRPKM(rse_gene))
  rpkm$id <- rownames(rpkm)
  rpkm.list[[experiment]] <- rpkm
}

# combine experiments -- this is the most memory efficient way to go about this
# that I've found -- will need to drop extraneous gene id columns
rpkm.df <- do.call(base::cbind, c(rpkm.list, by = "id"))
id.cols <- grep("id", colnames(rpkm.df))
rpkm.df <- rpkm.df[, -id.cols[2:length(id.cols)]]
rpkm.df <- rpkm.df[, c(id.cols[1], 1:(id.cols[1] - 1),
                       (id.cols[1] + 1):ncol(rpkm.df))]
colnames(rpkm.df)[1] <- "ENSG"

# save to file
saveRDS(rpkm.df, file = file.path("recount2", "recount_rpkm.RDS"))
