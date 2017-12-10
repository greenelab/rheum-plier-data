# J. Taroni 2017

# set missing values to zero
NAToZero <- function(dt, un = 0) suppressWarnings(gdata::NAToUnknown(dt, un))

GetCombinedMatrix <- function(list.of.df) {
  # Concatenate a list of data.frames, return a matrix of expression data
  # without the first column that contains the gene identifiers
  # 
  # Args:
  #   list.of.df: a list of data.frames of gene expression data, rows are
  #               genes, columns are samples/arrays; first column ("Gene"),
  #               contains gene identifiers
  #               
  # Returns:
  #   a combined expression matrix, with gene identifiers as rownames, rather
  #    than the first column -- suitable for principal components
  
  combined.df <- plyr::join_all(list.of.df, by = "Gene", type = "inner")
  rownames(combined.df) <- combined.df$Gene
  combined.df <- dplyr::select(combined.df, -Gene)
  
  # any missing values? set to zero
  if (any(is.na(combined.df))) {
    combined.df <- NAToZero(combined.df)
  }
  
  return(as.matrix(combined.df))
  
}

ZTOProcessing <- function(list.of.df, before = TRUE) {
  # This function is for the [0, 1] scaling of a list of expression data.frames.
  # The transformation can either be performed prior to combining the datasets
  # (before = TRUE), or afterwards (before = FALSE).
  # 
  # Args:
  #   list.of.df: a list of data.frames that contain expression data, genes are
  #               rows, samples/arrays are columns, gene identifiers are in the
  #               first column ("Gene")
  #   before: logical - should [0, 1] scaling be performed before concatenation?
  #           default is TRUE
  #           
  # Returns:
  #   exprs.mat: an expression matrix that has been rescaled, genes are rows,
  #              samples/arrays are columns, gene identifiers are rownames;
  #              suitable for PCA/plotting
  
  require(data.table)
  
  # if before = TRUE
  if (before) {
    
    # convert to data.table for use with TDM::zero_to_one_transform()
    list.of.dt <- lapply(list.of.df,
                         function(x) data.table(x))
    
    # any NAs? will cause errors with scaling, so set to zero
    any.na.list <- lapply(list.of.dt, function(x) any(is.na(x)))
    
    # for data.tables that have NAs, set to zero
    for (na.iter in which(unlist(any.na.list))) {
      list.of.dt[[na.iter]] <- NAToZero(list.of.dt[[na.iter]])
    }
    
    # zero to one transform (each gene)
    list.of.zto <- lapply(list.of.dt, 
                          function(x) TDM::zero_to_one_transform(x))
    
    # convert back to data.frame for use with GetCombinedMatrix
    list.of.zto <- lapply(list.of.zto,
                          function(x) as.data.frame(x))
    
    # get concatenated (combined) expression matrix
    exprs.mat <- GetCombinedMatrix(list.of.zto)
    
  } else {
    
    # concatenate
    combined.mat <- GetCombinedMatrix(list.of.df = list.of.df)
    # convert to data.table
    combined.dt <- data.table::data.table(rownames(combined.mat), 
                                          combined.mat)
    
    # if any missing values exist - set to zero
    if (any(is.na(combined.dt))) {
      combined.dt <- NAToZero(combined.dt)
    }
    
    # zero to one transform
    zto.dt <- TDM::zero_to_one_transform(combined.dt)
    
    # convert to matrix
    exprs.mat <- as.matrix(zto.dt[, 2:ncol(zto.dt), with = FALSE])
    rownames(exprs.mat) <- zto.dt[[1]]
    
  }
  
  return(exprs.mat)
  
}

CombineDatasets <- function(list.of.pcl.files, UPC = FALSE) {
  # Combines a list of pcl files in a number of ways -- no [0,1] scaling
  # ("no.transform"), scaling before concatenation ("zto.before"), scaling
  # after concatenation ("zto.after"), and if UPC == TRUE, using 
  # transformation with SCAN.UPC::UPCGeneric ("upc.generic"). (The list of pcl 
  # files should be normalized using one method [e.g., RMA].)
  #
  # Args:
  #   list.of.pcl.files: a list of processed pcl files
  #   UPC: logical - should UPC generic processing be performed? default is 
  #        FALSE
  
  # internal function
  UPCGenericProcessing <- function(list.of.data.frames) {
    # UPC Generic processing of combined data
    combined.mat <- GetCombinedMatrix(list.of.df = list.of.data.frames)
    
    upc.mat <- apply(combined.mat, 2, function(x) SCAN.UPC::UPC_Generic(x))
    rownames(upc.mat) <- rownames(combined.mat)
    
    return(upc.mat)
    
  }
  
  pcl.list <- lapply(list.of.pcl.files, ReadInPCL)
  
  # initialize master list to hold all transformed data
  master.list <- list()
  # without transformation s.t. values fall [0, 1]
  master.list[["no.transform"]] <- GetCombinedMatrix(list.of.df = pcl.list)
  
  # zero to one transform before combining data
  master.list[["zto.before"]] <- ZTOProcessing(list.of.df = pcl.list,
                                               before = TRUE)
  # zero to one transform after combining the data
  master.list[["zto.after"]] <- ZTOProcessing(list.of.df = pcl.list,
                                              before = FALSE)
  # if UPC = TRUE (would be the case for RMA processed data, typically)
  if (UPC) {  # UPC generic processing
    master.list[["upc.generic"]] <- 
      UPCGenericProcessing(list.of.data.frames = pcl.list)
  }
  
  # how many arrays are in each pcl? need to exclude
  array.counts <- lapply(pcl.list, 
                         function(x) dim(dplyr::select(x, -Gene))[2])
  # get number of arrays infomration in a form suitable for visualization
  # with pairs (base graphics)
  group.vector <- vector()
  for (grp.iter in seq_along(array.counts)) {
    group.vector <- base::append(group.vector, rep(grp.iter, 
                                                   array.counts[[grp.iter]]))
  }
  
  # return a list of expression matrices and the group vector
  return(list("expression.matrices" = master.list, "groups" = group.vector))
  
}

PlotPCA <- function(exprs.mat, group.vector, png.filename) {
  # This function takes an expression matrix, information regarding the
  # dataset of origin for each array (group.vector), and a .png filename and
  # outputs a "pairs" plot of Principal Components 1-5 and a .tsv that contains
  # information about the cumulative variance explained for PC1-5
  # 
  # Args:
  #   exprs.mat: a matrix of expression data (gene ids as rownames, not 
  #              first column), where genes are rows and samples are columns
  #   group.vector: a vector of integers that indicate what data set each array/
  #                 sample came from
  #   png.filename: filename for .png output (including path information)
  #   
  # Returns:
  #   NULL - outputs a .png (pairs plot) and .tsv (cumulative variance 
  #          explained)
  
  # colorblind friendly palette
  cb.palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                  "#0072B2", "#D55E00", "#CC79A7")
  
  # perform PCA
  pc <- prcomp(t(exprs.mat))
  
  # save PNG of "paired" plot, PC1-5
  png(file = png.filename)
  pairs(pc$x[, 1:5], col = cb.palette[group.vector])
  dev.off()
  
  # what is the cumulative variance explained by PC1-5?
  cumvar.df <- as.data.frame(cumsum(pc$sdev^2 / sum(pc$sdev^2))[1:5])
  cumvar.df <- cbind(colnames(pc$x)[1:5], cumvar.df)
  colnames(cumvar.df) <- c("Principal Component", 
                           "cumulative variance explained")
  # write cum variance info to file
  tbl.name <- sub(".png", ".tsv", png.filename)
  write.table(cumvar.df, file = tbl.name, quote = FALSE, row.names = FALSE, 
              sep = "\t")
  
}

PCAWrapper <- function(list.of.pcl, UPC.arg, png.file.lead) {
  # A wrapper function for combining datasets and ultimately PCA plots -- 
  # motivated by the fact that multiple normalization methods will be under
  # consideration.
  # 
  # Args:
  #   list.of.pcl: a list of pcl files (should be from one normalization method)
  #   UPC.arg: logical - should UPC generic processing be performed?
  #   png.file.lead: .png "file lead" for PC1-5 pairs plot (includes path info)
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
            png.filename = png.transform.file)
    
  }
  
}

