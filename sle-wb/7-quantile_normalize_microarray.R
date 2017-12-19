# J. Taroni 2017
# 
# This processing script quantile normalizes submitter processed microarray 
# data (E-GEOD-49454, E-GEOD-65391, and E-GEOD-78193) using SCANfast processed 
# Affymetrix data as a reference.
# 
# USAGE: Rscript sle-wb/7-quantile_normalize_microarray.R
# 

source(file.path("util", "aggregate_norm_other_platforms.R"))
# use aggregated SCANfast affy data
reference.file <- file.path("sle-wb", "processed", "aggregated_data", 
                            "SLE_WB_aggregated_affy_data.pcl")
single.acc.dir <- file.path("sle-wb", "processed", "single_accession")

#### main ----------------------------------------------------------------------

# read in reference data 
reference.df <- ReadInPCL(reference.file)
# get rid of "_at" at the end of Entrez IDs from brainarray
reference.df$Gene <- gsub("_at$", "", reference.df$Gene)

# quantile normalize submitter processed data
# E-GEOD-49454 
e.49454.file <- file.path(single.acc.dir, 
                        "E-GEOD-49454_submitter_processed_entrez_mean.pcl")
QNfromPCL(ref.df = reference.df,
          pcl.filename = e.49454.file)

# E-GEOD-65391 
e.65391.file <- file.path(single.acc.dir, 
                          "E-GEOD-65391_submitter_processed_entrez_mean.pcl")
QNfromPCL(ref.df = reference.df,
          pcl.filename = e.65391.file)

# E-GEOD-78193 
e.78193.file <- file.path(single.acc.dir, 
                          "E-GEOD-78193_submitter_processed_entrez_mean.pcl")
QNfromPCL(ref.df = reference.df,
          pcl.filename = e.78193.file)
