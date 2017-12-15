# J. Taroni 2017

PlotPCA <- function(exprs.mat, group.vector, png.filename, 
                    legend.labels, pcs.to.viz = 1:5) {
  # This function takes an expression matrix, information regarding the
  # dataset of origin for each array (group.vector), and a .png filename and
  # outputs a "pairs" plot of Principal Components specified with
  # pcs.to.viz and a .tsv that contains information about the cumulative 
  # variance 
  # 
  # Args:
  #   exprs.mat: a matrix of expression data (gene ids as rownames, not 
  #              first column), where genes are rows and samples are columns
  #   group.vector: a vector of integers that indicate what data set each array/
  #                 sample came from
  #   png.filename: filename for .png output (including path information)
  #   legend.labels: a character vector of the labels to be used in the pairs
  #                  plot legend
  #   pcs.to.viz: a vector of integers specifying the PCs that should
  #               be visualized (1:5 by default)
  #   
  # Returns:
  #   NULL - outputs a .png (pairs plot) and .tsv (cumulative variance 
  #          explained)
  #
  
  # colorblind friendly palette
  cb.palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                  "#0072B2", "#D55E00", "#CC79A7")
  
  # perform PCA
  pc <- prcomp(t(exprs.mat))
  
  # save PNG of "paired" plot, specified PCs
  png(file = png.filename, width = 600, height = 480, units = "px")
  pairs(pc$x[, pcs.to.viz], col = cb.palette[group.vector], 
        oma = c(4, 4, 6, 15))
  par(xpd = TRUE)
  legend(x = 0.825, y = 0.7, legend = legend.labels, 
         fill = cb.palette[1:max(group.vector)], bty = "n")
  dev.off()
  
  # what is the cumulative variance explained by PCs?
  cumvar.df <- as.data.frame(cumsum(pc$sdev^2 / sum(pc$sdev^2))[pcs.to.viz])
  cumvar.df <- cbind(colnames(pc$x)[pcs.to.viz], cumvar.df)
  colnames(cumvar.df) <- c("Principal Component", 
                           "cumulative variance explained")
  # write cum variance info to file
  tbl.name <- sub(".png", ".tsv", png.filename)
  write.table(cumvar.df, file = tbl.name, quote = FALSE, row.names = FALSE, 
              sep = "\t")
  
}

PCAWrapper <- function(list.of.pcl, UPC.arg, png.file.lead,
                       dataset.labels, pcs.to.plot = 1:5) {
  # A wrapper function for combining datasets and ultimately PCA plots -- 
  # motivated by the fact that multiple normalization methods will be under
  # consideration.
  # 
  # Args:
  #   list.of.pcl: a list of pcl files (should be from one normalization method)
  #   UPC.arg: logical - should UPC generic processing be performed?
  #   png.file.lead: .png "file lead" for PC pairs plot (includes path info)
  #   dataset.labels: a character vector of the dataset names to be used in 
  #                   the pairs plot legend
  #   pcs.to.plot: a vector of integers specifying the PCs that should
  #                be visualized (1:5 by default)
  #   
  # Returns:
  #   NULL - outputs a .png (pairs plot) and .tsv (cumulative variance 
  #          explained)
  
  # get combined expression data list and group information
  combined.list <- CombineDatasets(list.of.pcl.files = list.of.pcl,
                                   UPC = UPC.arg)
  
  # for each processed expression matrix in combined list
  for (exprs.iter in seq_along(combined.list$expression.matrices)) {
    
    # append png to include transformation info
    png.transform.file <- 
      paste0(png.file.lead, "_", 
             names(combined.list$expression.matrices)[exprs.iter], ".png")
    
    # plot PC1-5 "paired"
    PlotPCA(exprs.mat = combined.list$expression.matrices[[exprs.iter]],
            group.vector = combined.list$groups,
            png.filename = png.transform.file,
            legend.labels = dataset.labels,
            pcs.to.viz = pcs.to.plot)
    
  }
  
}

