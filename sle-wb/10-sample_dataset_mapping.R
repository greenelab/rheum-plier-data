# J. Taroni 2018 
# This script is for obtaining a sample id to dataset mapping for the SLE
# whole blood data.
# 
# Usage: Rscript sle-wb/10-sample_dataset_mapping.R
# 

# list files in the single_accession directory
pcl.list <- list.files(file.path("sle-wb", "processed", "single_accession"), 
                       full.names = TRUE)

# get ArrayExpress accessions from the list of pcl files
accession.list <- lapply(pcl.list, 
                         function(x) gsub("[_].*", "", gsub(".*/", "", x)))

# get vector of unique accession/experiment numbers
accessions <- unique(unlist(accession.list))

# remove E-GEOD-11907 -- this was hgu133a + b, which was concatenated -- the
# concatenated file will have the correct sample names for the purposes of
# this processing
accessions <- accessions[which(accessions != "E-GEOD-11907")]

# initialize list to hold sample names for each experiment
sample.list <- list()
for (acc in accessions) {
  current.file <- pcl.list[[grep(acc, accession.list)[1]]]
  # read in file
  current.pcl <- data.table::fread(current.file, data.table = FALSE)
  # get sample names from columns
  sample.names <- colnames(current.pcl)[2:ncol(current.pcl)]
  # put in data.frame form
  current.df <- as.data.frame(cbind(sample.names, 
                                    rep(acc, length(sample.names))))
  colnames(current.df) <- c("SampleID", "Dataset")
  # add data.frame to sample.list
  sample.list[[acc]] <- current.df
}

# now E-GEOD-11907 sample information
file.11907 <- file.path("sle-wb", "processed", "aggregated_data", 
                        "E-GEOD-11907_aggregated_SLE_HC_only.pcl")
# read in file
pcl.11907 <- data.table::fread(file.11907, data.table = FALSE)
# data.frame 
sample.list[["E-GEOD-11907"]] <- 
  as.data.frame(cbind(colnames(pcl.11907)[2:ncol(pcl.11907)], 
                      rep("E-GEOD-11907", ncol(pcl.11907) - 1)))
colnames(sample.list[["E-GEOD-11907"]]) <- c("SampleID", "Dataset")

# bind rows to get "master" df
sample.dataset.df <- plyr::rbind.fill(sample.list)
# "master" filename
sdf.file <- file.path("sle-wb", "processed", 
                      "sle-wb_sample_dataset_mapping.tsv")
# write to file
readr::write_tsv(sample.dataset.df, path = sdf.file)
