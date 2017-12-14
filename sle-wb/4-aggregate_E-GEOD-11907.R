# J. Taroni 2017
# 
# One of the accessions (experiments) under consideration in this project
# was run on the Affymetrix HGU133A and HGU133B -- which can be considered an
# earlier version of the Affymetrix HGU133plus2 platform, which is the most 
# common platform in this small compendium (and most likely, the most commonly 
# used platform in all of ArrayExpress). This experiment is from the original
# blood modular framework paper (Chaussabel, et al. Immunity. 2008) and contains
# samples from patients that do not have SLE (including those with acute 
# infections, type I diabetes, etc.).
# 
# The purpose of this script is to aggregate/concatenate "matched" samples run
# on A and B platforms (such that they are more comparable to HGU133plus2
# data) and to only include samples from patients with SLE. We will be using 
# SCANfast processed data. We will be using only sample metadata to "match" A
# and B -- future analyses will use the transcriptomic data itself / imputation
# strategies.
# 
# USAGE: Rscript sle-wb/4-aggregate_E-GEOD-11907.R

# magrittr pipe
`%>%` <- dplyr::`%>%`
# aggregate functions
source(file.path("util", "aggregate_norm_other_platforms.R"))

# aggregated data directory
agg.dir <- file.path("sle-wb", "processed", "aggregated_data")
dir.create(agg.dir, recursive = TRUE)

#### functions -----------------------------------------------------------------

# this function is specific to the E-GEOD-11907
JoinAandBDF <- function(a.colname, a.df, b.colname, b.df, sample.id) {
  # This function is for binding rows from HGU133A and HGU133B matched samples
  # together.
  # 
  # Args:
  #   a.colname: the column name of the sample on HGU133A
  #   a.df: a data.frame of gene expression data, rows are genes, samples are
  #         columns, the first column contains gene identifiers and is named
  #         "Gene"
  #   b.colname: the column name of the sample on HGU133B
  #   b.df: a data.frame of gene expression data, rows are genes, samples are
  #         columns, the first column contains gene identifiers and is named
  #         "Gene"
  #   
  #
  #
  # get A sample, include gene identifiers
  a.single.df <- dplyr::select_(a.df, .dots = c("Gene", a.colname))
  colnames(a.single.df)[2] <- "Sample"
  # get B sample, include gene identifiers
  b.single.df <- dplyr::select_(b.df, .dots = c("Gene", b.colname))
  colnames(b.single.df)[2] <- "Sample"
  
  # bind rows together
  agg.df <- dplyr::bind_rows(a.single.df, b.single.df)
  # rename "Sample" colname to user supplied arg sample.id
  colnames(agg.df)[2] <- sample.id

  # return "aggregated" data.frame
  return(agg.df)
    
}

#### sample and data relationship data -----------------------------------------

sdrf.file <- file.path("sle-wb", "arrayexpress", "E-GEOD-11907", 
                       "E-GEOD-11907.sdrf.txt")
sample.rel.df <- data.table::fread(sdrf.file , data.table = FALSE)
colnames(sample.rel.df) <- make.unique(make.names(colnames(sample.rel.df)))

# select a subset of the metadata to be altered
select.cols <- c("Source.Name", "Hybridization.Name", "Array.Design.REF" ,
                 "Array.Data.File")
sample.meta.df <- dplyr::select_(sample.rel.df, .dots = select.cols)

# add a column that holds disease state information -- we only want healthy 
# and SLE
sample.meta.df$Disease.state <- rep(NA, nrow(sample.meta.df))
sample.meta.df$Disease.state[grep("SLE", sample.meta.df$Hybridization.Name)] <-
  "SLE"

# healthy identifiers
healthy.ids <- c("H 27", "H 28", "H 30", "H 36", "H 37", "H 45", "H 46")
sample.meta.df$Disease.state[grep(paste(healthy.ids, collapse = "|"), 
                                  sample.meta.df$Hybridization.Name)] <-
  "Control"
sample.meta.df$Disease.state[grep(paste(sub(" ", "", healthy.ids), 
                                        collapse = "|"), 
                                  sample.meta.df$Hybridization.Name)] <-
  "Control"

# filter to only the relevant disease states
sle.meta.df <- dplyr::filter(sample.meta.df, !is.na(Disease.state))

#### HGU133A + B processed -----------------------------------------------------

a.file <- file.path("sle-wb", "processed", "single_accession", 
                    "E-GEOD-11907_hgu133a_SCANfast.pcl")
b.file <- file.path("sle-wb", "processed", "single_accession", 
                    "E-GEOD-11907_hgu133b_SCANfast.pcl")

hgu133a.full <- ReadInPCL(a.file)
hgu133b.full <- ReadInPCL(b.file)

#### matching samples ----------------------------------------------------------

# each sample that ends in "B" -- only the SLE, not the controls marked this way
b.sle.samples <- 
  dplyr::filter(sle.meta.df, 
                grepl("B$", Hybridization.Name))$Hybridization.Name

# initialize list to hold matched samples
matched.ab.list <- list()

# for each sample marked B
for (smpl.iter in seq_along(b.sle.samples)) {
  
  # get the .CEL file name for HGU133B (now the colname of the expression data)
  matched.smpl.B <- 
    sle.meta.df$Array.Data.File[which(sle.meta.df$Hybridization.Name == 
                                        b.sle.samples[smpl.iter])]
  
  # either substitute A for B, or leave off the letter -- this is the pattern
  # used by the submitters that can be discerned
  match.a.possibilities <- c(gsub("B$", "A", b.sle.samples[smpl.iter]),
                             gsub("B$", "", b.sle.samples[smpl.iter]))
  
  # any matches?
  any.a.match <- any(match.a.possibilities %in% sle.meta.df$Hybridization.Name)
  
  # if a HGU133A match exists
  if (any.a.match) {
    
    # use %in% because we want an exact match
    match.a <- match.a.possibilities[which(match.a.possibilities %in%
                                             sle.meta.df$Hybridization.Name)]
    
    # get the .CEL file name for HGU133A
    matched.smpl.A <- 
      sle.meta.df$Array.Data.File[which(sle.meta.df$Hybridization.Name == 
                                          match.a)]
    
    # add to list
    matched.ab.list[[b.sle.samples[smpl.iter]]] <- c(matched.smpl.A, 
                                                     matched.smpl.B)
  }
  
}

### matching for healthy controls

for (h.id in healthy.ids) {
  
  # can we match the healthy id in the hybridization name?
  any.match <- any(grepl(h.id, sle.meta.df$Hybridization.Name))
  # is that match on the correct platform (A-AFFY-33, "A")
  correct.platform <- 
    (sle.meta.df$Array.Design.REF[grep(h.id, 
                                      sle.meta.df$Hybridization.Name)] ==
       "A-AFFY-33")
  
  if (correct.platform & any.match) {
    
    # sample file name from A
    match.hid.a <- 
      sle.meta.df$Array.Data.File[grep(h.id, sle.meta.df$Hybridization.Name)]
    
    # sample file name from B
    match.hid.b <- 
      sle.meta.df$Array.Data.File[grep(sub(" ", "", h.id), 
                                       sle.meta.df$Hybridization.Name)] 
    
    # add to list of matched samples pairs
     matched.ab.list[[h.id]] <- c(match.hid.a, match.hid.b)

  } 
  
}


#### aggregate matched samples -------------------------------------------------

# initialize list to hold aggregated data.frames
agg.df.list <- list()

# for each pair of matched samples
for (match.iter in seq_along(matched.ab.list)) {

  # error handling 
  # is a.colname in a.df?
  in.a <- matched.ab.list[[match.iter]][1] %in% colnames(hgu133a.full)
  # is b.colname in b.df?
  in.b <- matched.ab.list[[match.iter]][2] %in% colnames(hgu133b.full)
  if (all(c(in.a, in.b))) {
    agg.df.list[[length(agg.df.list) + 1]] <- 
      JoinAandBDF(a.colname = matched.ab.list[[match.iter]][1],
                  a.df = hgu133a.full,
                  b.colname = matched.ab.list[[match.iter]][2],
                  b.df = hgu133b.full,
                  sample.id = names(matched.ab.list)[match.iter])
  }
  
} 

# this is the most memory efficient way to do this, albeit not the most
# pretty 
agg.df <- do.call(base::cbind, c(agg.df.list, by = "Gene"))
cols.to.rm <- c(which(colnames(agg.df) == 
                        "Gene")[2:length(agg.df.list)], 
                ncol(agg.df))
agg.df <- agg.df[, -cols.to.rm]  

# for duplicate gene identifiers, take the gene mean
agg.df <- PrepExpressionDF(agg.df)

#### write processed data to file ----------------------------------------------

# expression data
exprs.output.file <- file.path(agg.dir,
                               "E-GEOD-11907_aggregated_SLE_HC_only.pcl")
write.table(agg.df, file = exprs.output.file, quote = FALSE, row.names = FALSE,
            sep = "\t")

# altered metadata
sdrf.output.file <- file.path("sle-wb", "processed", 
                              "E-GEOD-11907_SLE_only.sdrf.tsv")
write.table(sle.meta.df, file = sdrf.output.file, quote = FALSE, 
            row.names = FALSE, sep = "\t")
