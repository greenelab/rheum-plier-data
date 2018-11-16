# J. Taroni 2018
# This script processes E-GEOD-26576 data from raw with SCANfast
# 
# USAGE: Rscript DIPG/E-GEOD-26576/process_DIPG.R
# 

# functions for Affymetrix processing
source(file.path("util", "affy_norm_functions.R"))

# directory to store the PCL
processed.dir <- file.path("DIPG", "E-GEOD-26576", "processed")
if (!dir.exists(processed.dir)) {
  dir.create(processed.dir, recursive = TRUE)
}

# normalization itself
AffyMultiNormWrapper(cel.dir = file.path("DIPG", "E-GEOD-26576", "raw"), 
                     output.file.lead = file.path(processed.dir, 
                                                  "DIPG_E-GEOD-26576"),
                     norm.methods = "SCAN", 
                     fast = TRUE)
