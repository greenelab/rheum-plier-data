# J. Taroni 2017
# Install Brainarray packages (v21.0.0) for use with SCAN.UPC and affy packages.
# The platforms we know will be used for this compendium are hardcoded into this 
# script; it does not accept command line arguments.
# USAGE (from top directory): Rscript annotation/install_brainarray.R

#### functions -----------------------------------------------------------------

InstallBrainarray <- function(platform, org.code, ba.version) {
	# This function makes use of devtools::install_url to install Brainarray
	# packages for the annotation of Affymetrix data. Specifically, the packages
	# required for use with SCAN.UPC and affy (RMA) are installed.
	#
	# Args:
	#   platform: The Affymetrix platform for which brainarray
	#             packages are to be installed (e.g., "hgu133plus2")
	#   org.code: Two letter organism code -- human would be "hs"
	#		ba.version: What version of brainarray should be used? (e.g., "21.0.0")
	#
	#	Returns:
	#	  NULL - this package completes installation of these packages and does not
	#					 return any values 

	# make sure platform and org.code are all lowercase and lack punctuation
	platform <- tolower(gsub("[[:punct:]]", "", platform))
	org.code <- tolower(gsub("[[:punct:]]", "", org.code))

	# probe version for use with SCAN.UPC
	probe.pkg.name <- paste0(platform, org.code, "entrezgprobe_",
                         ba.version, ".tar.gz")
	probe.url <- paste0("http://mbni.org/customcdf/", ba.version, 
	                    "/entrezg.download/", probe.pkg.name)
	devtools::install_url(probe.url)

	# cdf version for use with affy::RMA
	cdf.pkg.name <- paste0(platform, org.code, "entrezgcdf_",
						      	     ba.version, ".tar.gz")
	cdf.url <- paste0("http://mbni.org/customcdf/", ba.version, 
										"/entrezg.download/", cdf.pkg.name)
	devtools::install_url(cdf.url)

}

#### install brainarray main ---------------------------------------------------

# HGU133Plus2
InstallBrainarray(platform = "hgu133plus2",
                  org.code = "hs",
                  ba.version = "21.0.0")

# HGU133A
InstallBrainarray(platform = "hgu133a",
                  org.code = "hs",
                  ba.version = "21.0.0")

# HGU133B
InstallBrainarray(platform = "hgu133b",
                  org.code = "hs",
                  ba.version = "21.0.0")

# HGU95Av2
InstallBrainarray(platform = "hgu95av2",
                  org.code = "hs",
                  ba.version = "21.0.0")

# hugene11st
InstallBrainarray(platform = "hugene11st",
                  org.code = "hs",
                  ba.version = "21.0.0")

# hugene10st
InstallBrainarray(platform = "hugene10st",
                  org.code = "hs",
                  ba.version = "21.0.0")
