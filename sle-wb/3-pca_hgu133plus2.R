# J. Taroni 2017
# The purpose of this analysis is to explore the existence of "batch effects" --
# here we are referring to differences between experiments (e.g., GEO accessions
# ) -- without additional processing beyond normalization or following some 
# transformation s.t. the expression values are between zero and one. We are
# examining RMA normalization, SCANfast normalization, and UPC (SCAN.UPC) v. 
# zero to one transformation (gene-level feature scaling).
# USAGE: Rscript sle-wb/3-pca_hgu133plus2.R

# req'd for TDM
library(data.table)

set.seed(5689)

source(file.path("util", "viz_functions.R"))
source(file.path("util", "aggregate_norm_other_platforms.R"))

# path for plotting data
plot.dir <- file.path("sle-wb", "plots", "PCA")
dir.create(plot.dir, recursive = TRUE)

# path to processed data 
processed.dir <- file.path("sle-wb", "processed", "single_accession")

#### main ----------------------------------------------------------------------

plus2.files <- list.files(processed.dir, 
                          pattern = "hgu133plus2", 
                          full.names = TRUE)

# RMA normalized data
RMA.files <- plus2.files[grep("_RMA.pcl", plus2.files)]
PCAWrapper(list.of.pcl = RMA.files,
           UPC.arg = FALSE,
           png.file.lead = file.path(plot.dir, "HGU133PLUS2_RMA_PC1-5_pairs"))

# SCAN normalized data
SCAN.files <- plus2.files[grep("_SCAN.pcl", plus2.files)]
PCAWrapper(list.of.pcl = SCAN.files,
           UPC.arg = FALSE,
           png.file.lead = file.path(plot.dir, "HGU133PLUS2_SCAN_PC1-5_pairs"))

# SCANfast normalized data
SCANfast.files <- plus2.files[grep("_SCANfast.pcl", plus2.files)]
PCAWrapper(list.of.pcl = SCANfast.files,
           UPC.arg = FALSE,
           png.file.lead = file.path(plot.dir, 
                                     "HGU133PLUS2_SCANfast_PC1-5_pairs"))
