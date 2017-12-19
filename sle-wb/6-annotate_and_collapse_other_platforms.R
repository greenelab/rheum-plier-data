# J. Taroni 2017
# 
# The purpose of this script is to annotate submitter processed data --
# probe IDs are mapped to Entrez gene identifiers, and then collapse to gene
# mean.
# 
# USAGE: Rscript sle-wb/6-annotation_and_collapse_other_platforms.R

source(file.path("util", "aggregate_norm_other_platforms.R"))
ae.dir <- file.path("sle-wb", "arrayexpress")
process.dir <- file.path("sle-wb", "processed", "single_accession")

#### reformat submitter processed data -----------------------------------------

# list of processed data directory, output file pairs
dir.output.pair.list <- 
  list(c(file.path(ae.dir, "E-GEOD-49454", "processed"),
         file.path(process.dir,
                   "E-GEOD-49454_submitter_processed_data.pcl")), 
       c(file.path(ae.dir, "E-GEOD-65391", "processed"),
         file.path(process.dir,
                   "E-GEOD-65391_submitter_processed_data.pcl")),
       c(file.path(ae.dir, "E-GEOD-78193", "processed"), 
         file.path(process.dir,
                   "E-GEOD-78193_submitter_processed_data.pcl")))

lapply(dir.output.pair.list, 
       function(x) ReformatSubmitterProcessedData(processed.dir = x[1], 
                                                  output.file = x[2]))

#### annotate and collapse data ------------------------------------------------

# submitter processed files
sub.process.list <- lapply(dir.output.pair.list, function(x) x[2])

# for each file, annotate and collapse
for (fl in sub.process.list) {
  
  # read in reformatted, submitter processed data
  exprs.df <- data.table::fread(fl, data.table = FALSE)
    
  # E-GEOD-78193 is Agilent 4x44K data
  if (grepl("E-GEOD-78193", fl)) {
    plt <- "hgug4112a"
  } else {  # the other two experiments are Illumina HT-12 v4.0
    plt <- "illuminaHumanv4"
  }

  # annotate and collapse
  ann.agg.df <- AnnotateCollapseExprs(exprs.df = exprs.df, 
                                      platform = plt)
  
  # alter submitter processed file name to reflect annotation and collapsing
  collapse.file <- gsub("_data.pcl", "_entrez_mean.pcl", fl)
  # write to file
  write.table(ann.agg.df, file = collapse.file, row.names = FALSE, 
              quote = FALSE, sep = "\t")

}
