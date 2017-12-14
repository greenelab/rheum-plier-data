# J. Taroni 2017
# The purpose of this script is to aggregate all Affymetrix data in this 
# compendium. This "unscaled" version is suitable for use as a reference 
# distribution for quantile normalization of other platforms (specifically those
# without raw data available on ArrayExpress). We'll be using SCANfast 
# processed Affy data.
# 
# USAGE: Rscript sle-wb/5-aggregate_all_affy_data.R

source(file.path("util", "aggregate_norm_other_platforms.R"))

# aggregated data directory
agg.dir <- file.path("sle-wb", "processed", "aggregated_data")

#### read in data --------------------------------------------------------------

# identify SCANfast processed HGU133PLUS2 pcl files
plus2.files <- list.files(path = file.path("sle-wb", "processed", 
                                           "single_accession"),
                          pattern = "hgu133plus2_SCANfast",
                          full.names = TRUE)

# E-GEOD-11907 file -- has been aggregated in 4-aggregate_E-GEOD-11907.R
agg.ab.file <- file.path(agg.dir,
                        "E-GEOD-11907_aggregated_SLE_HC_only.pcl")

# initialize list to expression data from each accession/experiment/PCL file
pcl.list <- list()
# for each file -- read in data
for (fl in c(plus2.files, agg.ab.file)) {
  pcl.list[[length(pcl.list) + 1]] <- ReadInPCL(fl)
}

#### aggregation ---------------------------------------------------------------

# combine HGU133plus2 data
plus2.df <- GetCombinedDataset(pcl.list[1:3], 
                               return.class = "data.frame",
                               join.type = "inner")

# join -- only want genes represented on concatenated HGU133A 
# + HGU133B and HGU133plus2 -- only 6 genes from A+B are lost this way currently
all.agg.df <- dplyr::inner_join(pcl.list[[4]], plus2.df, by = "Gene")

#### write to file -------------------------------------------------------------

# aggregated, no scaling
no.transform.file <- file.path(agg.dir,
                               "SLE_WB_aggregated_affy_data.pcl")
readr::write_tsv(all.agg.df, path = no.transform.file)
