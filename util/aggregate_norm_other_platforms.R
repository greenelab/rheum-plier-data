# J. Taroni 2017

#### general purpose -----------------------------------------------------------

ReadInPCL <- function(pcl.file) {
  # Reads in PCL file and sets first colname to "Gene"
  # 
  # Args:
  #   pcl.file: pcl filenmae
  #   
  # Returns:
  #   pcl: a data.frame of expression data, genes are rows, columns are samples
  #        or arrays -- first column contains gene identifiers and is named
  #        "Gene"
  #        
  pcl <- data.table::fread(pcl.file, data.table = FALSE)
  colnames(pcl)[1] <- "Gene"
  
  # brainarray adds "_at" to the end of Entrez IDs, check if that's the case 
  # here, and if so remove
  if (any(grepl("_at$", pcl$Gene))) {
    pcl$Gene <- gsub("_at$", "", pcl$Gene)
  }
  
  # as.character
  pcl$Gene <- as.character(pcl$Gene)
  
  # return pcl
  return(pcl)
}

PrepExpressionDF <- function(exprs){
  # Takes a data.frame where the first columns contains gene identifiers.
  # Returns matrix of expression data, collapsed to gene level (mean), where
  # the gene identifiers are rownames.
  # 
  # Args:
  #   exprs: A data.frame of (normalized) gene expression data, where the
  #          rows are genes and the samples are columns. The first column
  #          should contain gene identifiers and have the column name "Gene".
  #          
  # Returns: 
  #   exprs.agg: a data.frame of aggregated expression data
  
  # error handling
  if (colnames(exprs)[1] != "Gene") {
    stop("The first column name of exprs must be named 'Gene'.")
  }
  
  require(dplyr)
  
  # for duplicate gene identifiers, take the average
  exprs.agg <- exprs %>%
    group_by(Gene) %>%
    summarise_each(funs(mean(., na.rm = TRUE)))
  
  # return aggregated expression data.frame
  return(exprs.agg)
  
}

GetCombinedDF <- function(list.of.df) {
  # Concatenate a list of data.frames, return a data of expression data
  # where the first column that contains the gene identifiers
  # 
  # Args:
  #   list.of.df: a list of data.frames of gene expression data, rows are
  #               genes, columns are samples/arrays; first column ("Gene"),
  #               contains gene identifiers
  #               
  # Returns:
  #   a combined expression data.frame, with gene identifiers as first column
  #   
  
  combined.df <- do.call(base::cbind, c(list.of.df, by = "Gene"))
  cols.to.rm <- c(which(colnames(combined.df) == 
                          "Gene")[2:length(list.of.df)], 
                  ncol(combined.df))
  combined.df <- combined.df[, -cols.to.rm]  
  return(combined.df)
  
}

#### submitter processed data --------------------------------------------------

ReformatSubmitterProcessedData <- function(processed.dir, output.file) {
  # This function reformats submitter processed data from ArrayExpress --
  # stored as individual sample tables -- into PCL files. All files in the 
  # user supplied directory that contain "_sample_table.txt" will be included
  # in the output PCL file.
  #
  # Args:
  #   processed.dir: directory that contains submitter processed data in
  #                  the form of sample tables
  #   output.file: filename for PCL output (include path info where 
  #                necessary)
  #
  # Returns:
  #   NULL - outputs PCL file
  #   
  
  # list sample table files in supplied directory
  processed.files <- list.files(path = processed.dir, 
                                pattern = "sample_table.txt", 
                                full.names = TRUE)
  
  # initialize list to hold processed data
  processed.list <- list()
  
  # for each file (sample table)
  for (fl in processed.files) {
    
    # read in sample table
    processed.list[[fl]] <- data.table::fread(fl, data.table = FALSE)

  }
  
  join.by.id <- colnames(processed.list[[1]])[1]
  processed.df <- plyr::join_all(processed.list, by = join.by.id,
                                 type = "inner")
  
  # use GEO sample accession -- get from file names
  smpl.names <- gsub(".*/", "", names(processed.list))
  smpl.names <- gsub("_.*", "", smpl.names)
  colnames(processed.df)[2:ncol(processed.df)] <- smpl.names
  # remove list
  rm(processed.list)
  
  # write to file
  write.table(processed.df, file = output.file, quote = FALSE, sep = "\t",
              row.names = FALSE)
  
}

AnnotateCollapseExprs <- function(exprs.df, platform) {
  # This function takes a reformatted expression data.frame and aggregates it 
  # by: 1) mapping the probe identifiers to Entrez ID and 2) for duplicate
  # Entrez ID, collapse to gene mean
  # 
  # Args:
  #   exprs.df: a data.frame of gene expression data, genes are rows, samples/
  #             arrays are columns, probe IDs are in the first column
  #   platform: What platform is the exprs.df from? in the format s.t.
  #             annotation package from bioconductor can be loaded 
  #             (e.g. "hgug4112a")
  #
  # Returns:
  #   agg.exprs.df: gene expression data, collapsed to Entrez ID mean
  #
  
  # load library for annotation
  library(paste0(platform, ".db"), character.only = TRUE)
  # get entrez id object from annotation package
  entrez.obj <- get(paste0(platform, "ENTREZID"))
  
  # get a list of entrez ids with mapped probes
  entrez.map.df <- reshape2::melt(as.list(entrez.obj[mappedkeys(entrez.obj)]))
  colnames(entrez.map.df) <- c("Gene", "ProbeID")
  
  # rename probe identifier column in expression data.frame
  colnames(exprs.df)[1] <- "ProbeID"
  
  # match probe ids to entrez identifiers, dropping any probes that do not 
  # have a corresponding entrez id by doing a left_join here
  filt.exprs.df <- dplyr::left_join(entrez.map.df, exprs.df, by = "ProbeID")
  
  # get rid of probe identifiers
  entrez.exprs.df <- dplyr::select(filt.exprs.df, -ProbeID)
  
  # average over duplicate entrez ids
  agg.exprs.df <- PrepExpressionDF(entrez.exprs.df)
  
}

#### quantile normalization ----------------------------------------------------

QNwithRef <- function(ref.df, target.df) {
  # This function quantile normalizes the target.df using the ref.df as a 
  # target. These two data.frames of gene expression data do not yet have to
  # be collapsed to the overlapping set of genes; this function performs that
  # task first.
  # 
  # Args:
  #   ref.df: a data.frame of gene expression data, rows are genes, samples/
  #           arrays are columns, the first column contains gene ids; 
  #           to be used as target
  #   target.df: a data.frame of gene expression data, rows are genes, samples/
  #              arrays are columns, the first column contains gene ids; 
  #              to be quantile normalized
  #              
  # Returns:
  #   qn.targ.df: the quantile normalized target.df gene expression data
  #   
  
  # error-handling
  if ((base::colnames(ref.df)[1] != "Gene") | 
      (base::colnames(target.df)[1] != "Gene")) {
    stop("The first column of ref.df and target.df should be called 'Gene' 
         and contain gene identifiers")
  }
  
  # find overlapping genes
  overlap.genes <- base::intersect(ref.df$Gene, target.df$Gene)
  
  # if there is no overlap, return an error
  if (length(overlap.genes) == 0) {
    stop("No overlapping genes found between ref.df and target.df")
  }
  # output number of overlapping genes -- would not want to use this approach
  # if the number of genes is (subjectively) low
  cat(paste("\nNumber of overlapping genes:", length(overlap.genes)))
  
  # filter expression data.frames to just overlapping genes 
  ref.df <- dplyr::filter(ref.df, Gene %in% overlap.genes)
  target.df <- dplyr::filter(target.df, Gene %in% overlap.genes)
 
  # make both Gene columns character
  ref.df$Gene <- as.character(ref.df$Gene)
  target.df$Gene <- as.character(target.df$Gene)
  
  # order by gene
  ref.df <- ref.df[order(ref.df$Gene), ]
  target.df <- target.df[order(target.df$Gene), ]
  
  if (!all.equal(ref.df$Gene, target.df$Gene)) {
    stop("Something went wrong: ref.df$Gene, target.df$Gene are not equal")
  }
  
  ## quantile normalization ##
  # determine target
  qn.ref <- 
    preprocessCore::normalize.quantiles.determine.target(
      data.matrix(ref.df[, 2:ncol(ref.df)]), 
      target.length = nrow(ref.df)
    )
  
  # quantile normalize using target
  qn.target <- 
    preprocessCore::normalize.quantiles.use.target(
      data.matrix(target.df[, 2:ncol(target.df)]), qn.ref, 
      copy = FALSE
    )
  
  # get into data.frame format with first column containing entrez gene ids
  qn.targ.df <- as.data.frame(cbind(target.df$Gene, qn.target))
  colnames(qn.targ.df)[1] <- "Gene"
  
  # return quantile normalized expression data
  return(qn.targ.df)
  
}

QNfromPCL <- function(ref.df, pcl.filename) {
  # This function is a wrapper for target quantile normalization. It takes
  # the reference gene expression data.frame and the PCL filename of the data
  # to be quantile normalized, appends the PCL filename to reflect that it has
  # been quantile normalized and writes to that new PCL filename.
  #
  # Args:
  #   ref.df: a data.frame of gene expression data, rows are genes, samples/
  #           arrays are columns, the first column contains gene ids; 
  #           to be used as target
  #   pcl.filename: the PCL (with path info if necessary) of the data to be
  #                 quantile normalized; used as target.df in QNwithRef
  #                 
  # Returns:
  #   NULL - outputs quantile normalized PCL 
  #   example: "path/to/pcl/target.pcl" input -> "path/to/pcl/target_QN.pcl"
  #                 
  # read in target PCL file
  target.df <- ReadInPCL(pcl.filename)
  # perform QN
  qn.target.df <- QNwithRef(ref.df = ref.df,
                            target.df = target.df)
  # append file name
  new.pcl.name <- sub(".pcl", "_QN.pcl", pcl.filename)
  # write to file
  write.table(qn.target.df, file = new.pcl.name, quote = FALSE, 
              row.names = FALSE, sep = "\t")
}
