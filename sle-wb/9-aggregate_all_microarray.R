# J. Taroni 2017
# 
# This processing script aggregates all microarray data in this compendium -- 
# includes [0, 1] scaling. Outputs two files -- one where the [0, 1] scaling
# follows concatenation and one where the scaling precedes concatenation.
# We're using the versions of submitter processed data that have been quantile
# normalized using the Affymetrix SCANfast data as a 
# reference, based on the PCA results (see 
# sle-wb/8-pca_visualization_all_microarray.R).
# 
# USAGE: Rscript sle-wb/9-aggregate_all_microarray.R

source(file.path("util", "aggregate_norm_other_platforms.R"))
source(file.path("util", "viz_functions.R"))
single.acc.dir <- file.path("sle-wb", "processed", "single_accession")
agg.dir <- file.path("sle-wb", "processed", "aggregated_data")

#### functions -----------------------------------------------------------------

WriteExprsMat2PCL <- function(exprs.mat, pcl.filename) {
  # This function takes a matrix of expression data where rows are genes,
  # columns are sample, and the rownames are gene ids and writes this to
  # a PCL file where the first column contains the gene identifiers.
  # 
  # Args:
  #   exprs.mat: a matrix of expression data where rows are genes,
  #              columns are sample, and the rownames are gene ids
  #   pcl.filename: output filename (including path info where necessary)
  #   
  # Returns:
  #   NULL - writes a PCL file
  #
  
  # convert to data.frame
  exprs.df <- as.data.frame(exprs.mat)
  # add gene identifiers as first column
  exprs.df <- cbind(rownames(exprs.df), exprs.df)
  colnames(exprs.df)[1] <- "Gene"
  # write to file
  write.table(exprs.df, file = pcl.filename, quote = FALSE, row.names = FALSE,
              sep = "\t")
}


#### main ----------------------------------------------------------------------

# pcl files
qn.pcl.files <- 
  c(file.path(single.acc.dir, 
              "E-GEOD-39088_hgu133plus2_SCANfast.pcl"),
    file.path(single.acc.dir,
              "E-GEOD-61635_hgu133plus2_SCANfast.pcl"),
    file.path(single.acc.dir, 
              "E-GEOD-72747_hgu133plus2_SCANfast.pcl"),
    file.path(agg.dir, 
              "E-GEOD-11907_aggregated_SLE_HC_only.pcl"),
    file.path(single.acc.dir, 
              "E-GEOD-49454_submitter_processed_entrez_mean_QN.pcl"),
    file.path(single.acc.dir, 
              "E-GEOD-65391_submitter_processed_entrez_mean_QN.pcl"),
    file.path(single.acc.dir, 
              "E-GEOD-78193_submitter_processed_entrez_mean_QN.pcl"))

# combine datasets
combined.datasets <- CombineDatasets(list.of.pcl.files = qn.pcl.files,
                                     UPC = FALSE)

# write [0,1] scaling after concatenation expression matrix to file
WriteExprsMat2PCL(exprs.mat = combined.datasets$expression.matrices$zto.after,
                  pcl.filename = file.path(agg.dir,
                                           paste0("SLE_WB_all_microarray_QN", 
                                           "_zto_after.pcl")))

# write [0,1] scaling before concatenation expression matrix to file
# write [0,1] scaling after concatenation expression matrix to file
WriteExprsMat2PCL(exprs.mat = combined.datasets$expression.matrices$zto.before,
                  pcl.filename = file.path(agg.dir,
                                           paste0("SLE_WB_all_microarray_QN", 
                                                  "_zto_before.pcl")))
