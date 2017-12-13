# J. Taroni 2017
# The purpose of this script is to process Affymetrix data from raw CEL files
# using multiple normalization methods (affy::RMA, SCAN.UPC::SCAN, 
# SCAN.UPC::SCANfast). This is currently written to specifically process data
# for use in the SLE whole blood compendium (directories and output are
# hardcoded in this file).
# USAGE: Rscript sle-wb/1-process_affy_data.R
#

# load wrapper function for multiple normalization methods
source(file.path("util", "affy_norm_functions.R"))

# raw data resides in arrayexpress
ae.dir <- file.path("sle-wb", "arrayexpress")
# processed data dir
process.dir <- file.path("sle-wb", "processed", "single_accession")
dir.create(process.dir, recursive = TRUE)

#### main ----------------------------------------------------------------------

# list of pairs of cel.dir and output.file.lead -- the arguments to 
# AffyMultiNormWrapper
pairs.dir.list <- list(c(file.path(ae.dir, "E-GEOD-11907", "raw"),
                         file.path(process.dir, "E-GEOD-11907")), 
                       c(file.path(ae.dir, "E-GEOD-39088", "raw"),
                         file.path(process.dir, "E-GEOD-39088")),
                       c(file.path(ae.dir, "E-GEOD-61635", "raw"),
                         file.path(process.dir, "E-GEOD-61635")),
                       c(file.path(ae.dir, "E-GEOD-72747", "raw"),
                         file.path(process.dir, "E-GEOD-72747")))

# for each accession, perform multi-method normalization
lapply(pairs.dir.list,
       function(x) AffyMultiNormWrapper(cel.dir = x[1],
                                        output.file.lead = x[2],
                                        norm.methods = "all"))
