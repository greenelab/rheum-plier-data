# J. Taroni 2017
# The purpose of this analysis is to examine how well RMA (a heavily-used
# multi-sample normalization method) values correlate with SCAN.UPC::SCAN,
# and SCAN.UPC::SCANfast values; it calculates the Spearman correlation 
# coefficient (for each array/sample) and outputs density plots of the 
# correlation values.
# USAGE: Rscript sle-wb/2-calculate_spearman_affy_norm.R

source(file.path("util", "aggregate_norm_other_platforms.R"))

# directory where processed PCL
pcl.dir <- file.path("sle-wb", "processed", "single_accession")
png.output.path <- file.path("sle-wb", "plots", "norm_method_correlation")
dir.create(png.output.path, recursive = TRUE)

#### functions -----------------------------------------------------------------

CompareSampleCorrelation <- function(first.data.frame, second.data.frame) {
  # This function calculates the Spearman correlation coefficient on a per
  # sample (per array) basis between the two data.frames supplied as arguments.
  # 
  # Args:
  #   first.data.frame: a data.frame of processed/normalized expression data
  #                     where genes are rows and columns are samples; the first
  #                     column should contain gene identifiers
  #   second.data.frame: a data.frame of processed/normalized expression data
  #                      where genes are rows and columns are samples; the first
  #                      column should contain gene identifiers
  #
  # Returns:
  #   correlation.vector: a vector of Spearman correlation coefficients, for
  #                       each sample
  
  ## Error-handling ##
  # data.frames should be of equal dimensions
  if (!all.equal(dim(first.data.frame), dim(second.data.frame))) {
    stop("first.data.frame and second.data.frame must have equal dimensions")
  }
  # colnames (sample names) must be the same
  if (!all.equal(colnames(first.data.frame), colnames(second.data.frame))) {
    stop("colnames of first.data.frame and second.data.frame should
         all be the same")
  }
  # matching gene identifiers
  if (!all.equal(first.data.frame[[1]], second.data.frame[[1]])) {
    stop("the first column of both data.frames should contain matching
         gene identifiers (probes)")
  }
  
  # initialize vector to hold Spearman
  correlation.vector <- vector()
  # for all sample columns (excluding the first column which contains gene
  # identifiers)
  for (col.iter in 2:ncol(first.data.frame)) {
    correlation.vector[(col.iter - 1)] <- 
      cor(first.data.frame[, col.iter], 
          second.data.frame[, col.iter],
          method = "spearman")
  }
  
  # return Spearman values
  return(correlation.vector)
  
}

CompareSamplesPlotWrapper <- function(processed.dir, png.lead) {
  # This function is a wrapper for calculating Spearman correlations between
  # gene expression data processed with different normalization methods; also 
  # outputs .pngs of density plots of the Spearman coefficients
  # 
  # Args:
  #   processed.dir: path to processed PCL files directory (1 each: RMA, SCAN, 
  #                  SCANfast)  
  #   png.lead: png file "lead" (including path), platform and comparison
  #             of normalization will be appended to the file name
  #             example: "path/to/plot/E-GEOD-123XX" ->
  #                      "path/to/plot/E-GEOD-123XX_hgu133a_SCAN_v_SCANfast.png"
  #   
  # Returns:
  #   NULL - outputs png files (number of unique platforms) * 2 comparisons
  
  # identify all PCL files in the process directory supplied by the user
  pcl.files <- list.files(processed.dir, full.names = TRUE)
  pcl.files <- pcl.files[grep(".pcl", pcl.files)]
  
  # identify the unique platforms in the directory, will be indicated in the 
  # filename if processed with 1-process_affy_data.R
  platforms <- unique(sapply(pcl.files,
                             function(f) sub(".*\\_", "", 
                                             sub("_([^_]*)$", "", f))))
  
  # for each platform
  for (plt in platforms) {
    
    # initialize list to hold the pcl.files
    pcl.list <- list()
    
    # list the files in the directory that are from the platform under 
    # consideration
    plt.files.in.dir <- pcl.files[grep(plt, pcl.files)]
    # RMA processed file
    rma.file <- plt.files.in.dir[grep("_RMA.pcl", plt.files.in.dir)]
    # SCAN processed file
    scan.file <- plt.files.in.dir[grep("_SCAN.pcl", plt.files.in.dir)]
    # SCANfast processed file
    scan.fast.file <- plt.files.in.dir[grep("_SCANfast.pcl", plt.files.in.dir)]
    

    if (length(rma.file) > 1) {  
      pcl.list$rma <- GetCombinedDataset(lapply(rma.file, ReadInPCL), 
                                         return.class = "data.frame",
                                         join.type = "full")
      pcl.list$scan <- GetCombinedDataset(lapply(scan.file, ReadInPCL),
                                          return.class = "data.frame",
                                          join.type = "full")
      pcl.list$scan.fast <- GetCombinedDataset(lapply(scan.fast.file, 
                                                      ReadInPCL),
                                               return.class = "data.frame",
                                               join.type = "full")
    } else {  # if only one file per platform exists
      # read in all 3 files from different normalization methods
      pcl.list$rma <- ReadInPCL(rma.file)
      pcl.list$scan <- ReadInPCL(scan.file)
      pcl.list$scan.fast <- ReadInPCL(scan.fast.file)
    }
    
    # number of assays 
    n.assays <- dim(dplyr::select(pcl.list$rma, -Gene))[2]
    
    # initialize list to hold comparisons (Spearman)
    compare.list <- list()
    
    ### comparisons (Spearman) ###
    # RMA v. SCAN
    compare.list$rma.v.scan <- 
      as.data.frame(cbind(CompareSampleCorrelation(pcl.list$rma, 
                                                   pcl.list$scan), 
            rep("RMA v. SCAN", n.assays)))
    # RMA v. SCANfast
    compare.list$rma.v.scanfast <- 
      as.data.frame(cbind(CompareSampleCorrelation(pcl.list$rma, 
                                                   pcl.list$scan.fast), 
                    rep("RMA v. SCANfast", n.assays)))
    # SCAN v. SCANfast
    compare.list$scan.v.scanfast <- 
      as.data.frame(cbind(CompareSampleCorrelation(pcl.list$scan, 
                                                   pcl.list$scan.fast), 
                    rep("SCAN v. SCANfast", n.assays)))

    # rbind list of data.frames -- to be used to plot
    compare.df <- data.table::rbindlist(compare.list)
    colnames(compare.df) <- c("Spearman.coef", "Norm.compare")
      
    # get spearman coefficients as numeric
    if (class(compare.df$Spearman.coef) != "numeric") {
      compare.df$Spearman.coef <- 
        as.numeric(as.character(compare.df$Spearman.coef))
    }
    
    ### plotting ###
    # will plot SCANfast v. SCAN separately, as they are expected to be
    # nearly perfectly correlated and thus on a different scale than RMA 
    # comparison
    
    # RMA v. SCAN values
    rma.png.file <- paste0(png.lead, "_", plt, "_RMA_v_SCAN.png")
    ggplot2::ggplot(dplyr::filter(compare.df, 
                                  Norm.compare != "SCAN v. SCANfast"), 
                    ggplot2::aes(x = Spearman.coef)) + 
      ggplot2::geom_density() +  
      ggplot2::facet_grid(~ Norm.compare) + 
      ggplot2::theme_bw() +
      ggplot2::labs(x = "Sample Spearman correlation")
    ggplot2::ggsave(filename = rma.png.file,
                    plot = ggplot2::last_plot(),
                    width = 7,
                    height = 5)
    
    # SCAN v. SCANfast values
    scan.png.file <- paste0(png.lead, "_", plt, "_SCAN_v_SCANfast.png")
    ggplot2::ggplot(dplyr::filter(compare.df, 
                                  Norm.compare == "SCAN v. SCANfast"), 
                    ggplot2::aes(x = Spearman.coef)) + 
      ggplot2::geom_density() +  
      ggplot2::facet_grid(~ Norm.compare) + 
      ggplot2::theme_bw() +
      ggplot2::labs(x = "Sample Spearman correlation")
    ggplot2::ggsave(filename = scan.png.file,
                    plot = ggplot2::last_plot(),
                    width = 3.5,
                    height = 5)
    
  }
  
}

#### main ----------------------------------------------------------------------

CompareSamplesPlotWrapper(processed.dir = pcl.dir,
						  png.lead = file.path(png.output.path, 
						  					   "SLE-WB_affy_norm_correlation"))
