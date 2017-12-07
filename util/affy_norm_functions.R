# J. Taroni 2017

AffyMultiNormWrapper <- function(cel.dir, output.file.lead,
                                 norm.methods = "all", fast = TRUE) {
  # This function is a wrapper for multiple normalization methods to be used 
  # to process raw Affymetrix data (CEL files). It takes as arguments the 
  # the directory that contains the CEL files and a "file lead" (including path
  # information) for the output of the function. The platform and normalization 
  # method will be appended to the "file lead." All CEL files from a single 
  # platform within a directory will be processed together. (This is more 
  # pertinent for RMA, a multi-array normalization method.) 
  # SCANfast is a fast version of SCAN that uses a smaller number of probes. See
  # the SCAN.UPC primer for more information.
  #
  # Args:
  #   cel.dir: path to directory that contains CEL files
  #   output.file.lead: file "leader" (including path information), platform
  #                     and normalization method will be appended to filename
  #                     example: "path/to/processed/E-GEOD-123XX" ->
  #                     "path/to/processed/E-GEOD-123XX_hgu133plus2_SCAN.pcl"
  #   norm.methods: string that indcates which normalization methods should be 
  #                 used.
  #                 options are "RMA", "SCAN", or "all" which run RMA only, 
  #                 SCAN only (SCAN vs. SCANfast depends on fast arg), or all
  #                 methods (RMA, SCAN, SCANfast), respectively
  #   fast: if SCAN is going to be used for normalization, should 
  #         SCAN.UPC::SCANfast be used? default = TRUE; only relevant if 
  #         methods = "SCAN"
  #
  # Returns: 
  #   NULL - outputs (number of unique platforms) * no. of normalization 
  #          methods that are specified
  #          files into the directory specified as part of arg output.file.lead
  #

  # internal functions
  IdentifyAffyPlatforms <- function(raw.dir) {
    # identify the platform for each of the files in raw.dir
    
    # a list of all CEL files in the directory to be processed
    cel.files <- affy::list.celfiles(raw.dir)
  
    # platform.type is now a list with the type of each array
    platform.type <-
      sapply(cel.files,
             function (f) affyio::read.celfile.header(file.path(raw.dir, f))[1])
    
    # all lowercase, no dashes (or other punctuation)
    platform.type <- lapply(platform.type, 
                            function(x) tolower(gsub("[[:punct:]]", "", x)))
    
    return(platform.type)
  }
  
  RMAWrapper <- function (list.of.cel.files, platform, output.file) {
    # perform RMA normalization (affy package)
    pkg.name <- paste0(platform, "hsentrezgcdf")
    library(pkg.name, character.only = TRUE)
    
    cdf.name <- paste0(platform, "hsentrezg")
    
    data <- affy::ReadAffy(cdfname = cdf.name, filenames = list.of.cel.files)
    express <- affy::rma(data)
    
    Biobase::write.exprs(express, file = output.file)
      
  }
  
  SCANWrapper <- function(list.of.cel.files, 
                          platform, output.file, fast = FALSE) {
    # perform SCAN normalization, either SCAN.UPC::SCAN() or 
    # SCAN.UPC::SCANfast() if fast = TRUE
    library(foreach)
    # start parallel backend
    # parallel backend
    cl <- parallel::makeCluster(2)
    doParallel::registerDoParallel(cl)
    
    # load the correct package
    pkg.name <- paste0(platform, "hsentrezgprobe")
    library(pkg.name, character.only = TRUE)
    
    # initialize list to hold scan normalized data
    norm.list <- list()
    for (cel.iter in seq_along(list.of.cel.files)) {
      
      cel.file <- list.of.cel.files[[cel.iter]]
      cel.name <- gsub('.*\\/', '', cel.file)
      
      if (fast) {
        scan.norm <- SCAN.UPC::SCANfast(cel.file, 
                                        probeSummaryPackage = pkg.name)
      } else {
        scan.norm <- SCAN.UPC::SCAN(cel.file, 
                                    probeSummaryPackage = pkg.name)
      }
    
      norm.list[[cel.name]] <- scan.norm
      
    }
    
    # stop parallel backend
    parallel::stopCluster(cl)
    
    master.eset <- norm.list[[1]]
    for (lst.iter in 2:length(norm.list)) {
      master.eset <- Biobase::combine(master.eset, norm.list[[lst.iter]])
    }
    
    Biobase::write.exprs(master.eset, file = output.file)

  }
  
  ## Normalization methods 
  # sort out methods to be used based on norm.methods arg
  all.methods <- norm.methods == "all"
  rma.only <- norm.methods == "RMA"
  scan.only <- norm.methods == "SCAN"

  # error-handling based on norm.methods argument
  if (!any(c(all.methods, rma.only, scan.only))) {
    stop("Accepted norm.methods: 'RMA', 'SCAN', or 'all'")
  }

  # read the platform information from the .CEL files
  platform.info <- IdentifyAffyPlatforms(cel.dir)
  # get the unique list of platforms in the directory, they will be processed
  # independently
  unique.platforms <- unique(unlist(platform.info))

  # for each platform represented in the directory
  for (plt in unique.platforms) {
    message(paste("Normalizing data from", plt, "..."))
    
    # list the cel files on that platform
    platform.cel.files <- gsub(".cdfName", "", 
                               names(which(unlist(platform.info) == plt)))
    
    platform.cel.files <- file.path(cel.dir, platform.cel.files)
    
    # if the platform ends in 'v1' as will be the case for the Affy 
    # hugene arrays -- remove 'v1' from the platform name for the purpose
    # of loading the correct BrainArray package
    if (substr(plt, (nchar(plt) - 1), nchar(plt)) == "v1") {
      plt <- substr(plt, 1, (nchar(plt) - 2))
    }

    # if norm.methods = "RMA" or norm.methods = "all"
    if (rma.only | all.methods) {
    # run RMA
      rma.output.file <- paste(output.file.lead, plt, "RMA.pcl", sep = "_")
      message("  Running RMA...")
      RMAWrapper(list.of.cel.files = platform.cel.files,
                 platform = plt,
                 output.file = rma.output.file)
    }

    # if norm.methods = "SCAN" or norm.methods = "all"
    if (scan.only | all.methods) {
      # logic for SCAN v. SCANfast
      if (fast | all.methods) {
      # run SCANfast
        scan.fast.output.file <- paste(output.file.lead, plt, "SCANfast.pcl",
                                       sep = "_")
        message("  Running SCANfast...")
        SCANWrapper(list.of.cel.files = platform.cel.files,
                    platform = plt,
                    output.file = scan.fast.output.file,
                    fast = TRUE)
      } 
      
      if (!(fast) | all.methods) {
        # run SCAN
        scan.output.file <- paste(output.file.lead, plt, "SCAN.pcl", sep = "_")
        message("  Running SCAN...")
        SCANWrapper(list.of.cel.files = platform.cel.files,
                    platform = plt,
                    output.file = scan.output.file,
                    fast = FALSE)
      }
    }
  }
  
}  
