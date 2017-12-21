# J. Taroni 2017
# 
# The purpose of the analysis, much like the one performed in 
# 3-pca_hgu133plus2.R, is to visualize differences between experiments using
# Principal Components Analysis. This analysis is designed to determine if
# quantile normalization of submitter processed data using 
# SCANfast processed Affymetrix data is appropriate and to explore how
# [0, 1] scaling affects these differences between experiments.
# 
# USAGE: sle-wb/8-pca_visualization_all_microarray.R
# 

library(magrittr)

source(file.path("util", "viz_functions.R"))
source(file.path("util", "aggregate_norm_other_platforms.R"))
single.acc.dir <- file.path("sle-wb", "processed", "single_accession")
plot.dir <- file.path("sle-wb", "plots", "PCA")

plot.color.pal <- c("#54FF9F", "#43CD80", "#2E8B57", "#006400", "#FF8C00",
                    "#8B4500", "#000080")

#### with QN -------------------------------------------------------------------

qn.pcl.files <- 
  c(file.path(single.acc.dir, 
              "E-GEOD-39088_hgu133plus2_SCANfast.pcl"),
    file.path(single.acc.dir,
              "E-GEOD-61635_hgu133plus2_SCANfast.pcl"),
    file.path(single.acc.dir, 
              "E-GEOD-72747_hgu133plus2_SCANfast.pcl"),
    file.path("sle-wb", "processed", "aggregated_data", 
              "E-GEOD-11907_aggregated_SLE_HC_only.pcl"),
    file.path(single.acc.dir, 
              "E-GEOD-49454_submitter_processed_entrez_mean_QN.pcl"),
    file.path(single.acc.dir, 
              "E-GEOD-65391_submitter_processed_entrez_mean_QN.pcl"),
    file.path(single.acc.dir, 
              "E-GEOD-78193_submitter_processed_entrez_mean_QN.pcl"))

dataset.names <- paste0("E-GEOD-", 
                        sub(".*E-GEOD- *(.*?) *_.*", "\\1", qn.pcl.files))

PCAWrapper(list.of.pcl = qn.pcl.files,
           UPC.arg = FALSE,
           png.file.lead = file.path(plot.dir, 
                                     "SLE_WB_all_microarray_QN_PC1-5"),
           dataset.labels = dataset.names,
           color.pal = plot.color.pal)

#### without QN ----------------------------------------------------------------

pcl.files <- 
  c(file.path(single.acc.dir, 
              "E-GEOD-39088_hgu133plus2_SCANfast.pcl"),
    file.path(single.acc.dir,
              "E-GEOD-61635_hgu133plus2_SCANfast.pcl"),
    file.path(single.acc.dir, 
              "E-GEOD-72747_hgu133plus2_SCANfast.pcl"),
    file.path("sle-wb", "processed", "aggregated_data", 
              "E-GEOD-11907_aggregated_SLE_HC_only.pcl"),
    file.path(single.acc.dir, 
              "E-GEOD-49454_submitter_processed_entrez_mean.pcl"),
    file.path(single.acc.dir, 
              "E-GEOD-65391_submitter_processed_entrez_mean.pcl"),
    file.path(single.acc.dir, 
              "E-GEOD-78193_submitter_processed_entrez_mean.pcl"))

dataset.names <- paste0("E-GEOD-", 
                        sub(".*E-GEOD- *(.*?) *_.*", "\\1", pcl.files))

PCAWrapper(list.of.pcl = pcl.files,
           UPC.arg = FALSE,
           png.file.lead = file.path(plot.dir, 
                                     "SLE_WB_all_microarray_without_QN_PC1-5"),
           dataset.labels = dataset.names,
           color.pal = plot.color.pal)

#### make own color key --------------------------------------------------------

# get matrix to use as palette
palette.mat <- matrix(rep(1:length(dataset.names), 2), 
                      nrow = length(dataset.names), 
                      ncol = 2)
rownames(palette.mat) <- dataset.names

# use heatmap.2 to make key
png(file.path(plot.dir, "viz_all_microarray_colorkey.png"))
gplots::heatmap.2(palette.mat, Rowv = NA, Colv = NA, col = plot.color.pal,
                  key = FALSE, trace = "none", margins = c(1, 10))
dev.off()
