# J. Taroni 2017
# This script processes NARES data from raw using SCANfast and also performs
# batch correction with ComBat. There is some data cleaning along the way --
# specifically we must ensure the demographic data IDs match the array headers.
# 
# USAGE: Rscript NARES/process_NARES.R

# magrittr pipe + affy normalization wrapper
`%>%` <- dplyr::`%>%`
source(file.path("util", "affy_norm_functions.R"))

# processed NARES data directiory
processed.dir <- file.path("NARES", "processed")
dir.create(processed.dir, recursive = TRUE)

# set seed
set.seed(3487)

#### SCANfast ------------------------------------------------------------------

AffyMultiNormWrapper(cel.dir = file.path("NARES", "raw"), 
                     output.file.lead = file.path("NARES", "processed", 
                                                  "NARES"),
                     norm.methods = "SCAN", 
                     fast = TRUE)

#### Sample ID cleaning --------------------------------------------------------

# read in data
nares.data <- read.delim(file.path(processed.dir,
                                   "NARES_hugene10st_SCANfast.pcl"))
rownames(nares.data) <- nares.data[, 1]

# read in clinical data, output as .tsv
demo.data <- readxl::read_xlsx(file.path("NARES",
                                         "Demographic info for NARES DNA.xlsx"))
colnames(demo.data)[1] <- "Sample"
readr::write_tsv(demo.data, 
                 path = file.path("NARES", "NARES_demographic_data.tsv"))


# expression data matrix
exprs.mat <- as.matrix(dplyr::select(nares.data, -X))

# need to get the names to match the demographic info file -- in format
# 'NXXXX'

# identify characters between "_N" and "_" in colnames(exprs.mat) 
# and paste an "N" at the beginning -- this will return the correct ids for
# all but one array (N1093) that does not follow the "NXXXX" 
# format (lacks the N)
colnames(exprs.mat) <- paste0("N", 
                              sub(".*_N *(.*?) *_.*", "\\1", 
                                  colnames(exprs.mat)))
# correct id for N1093
colnames(exprs.mat)[grep("1093", colnames(exprs.mat))] <- "N1093"

#### PCA -----------------------------------------------------------------------

# PCA to look at batch effect
pca.results <- prcomp(t(exprs.mat))

# get PC1-2 as data.frame for use for ggplot2 -> scatterplot
pca.df <- cbind(rownames(pca.results$x),
                dplyr::select(as.data.frame(pca.results$x), PC1:PC2))
colnames(pca.df) <- c("Sample", "PC1", "PC2")

# extract relevant demographic data
rel.demo.df <- demo.data[, c("Sample", "Disease", "Batch", "Any_Immune")]

# join
pca.demo.df <- dplyr::full_join(rel.demo.df, pca.df, by = "Sample")
pca.demo.df$Batch <- as.factor(pca.demo.df$Batch)

# scatterplot
no.adjust.scatter <-
  ggplot2::ggplot(pca.demo.df, ggplot2::aes(x = PC1, y = PC2, colour = Batch,
                                            shape = Disease)) +
  ggplot2::geom_point(size = 3.5) +
  ggplot2::theme_bw() +
  ggplot2::annotate("text", x = pca.demo.df$PC1, y = pca.demo.df$PC2,
                    label = pca.demo.df$Any_Immune, colour = "#000000") +
  ggplot2::scale_color_manual(values = c("#FFB90F", "#1E90FF", "#00CD66")) +
  ggplot2::ggtitle("No Adjustment")

#### combat --------------------------------------------------------------------

assay.data <- exprs.mat[, which(colnames(exprs.mat) %in% rel.demo.df$Sample)]

pheno.data <- as.data.frame(rel.demo.df)
rownames(pheno.data) <- rel.demo.df$Sample
pheno.data <- dplyr::select(pheno.data, -Sample)
pheno.data <- Biobase::AnnotatedDataFrame(data = pheno.data)

eset <- Biobase::ExpressionSet(assayData = assay.data,
                               phenoData = pheno.data)

combat.output <- SCAN.UPC::BatchAdjust(eset, "Batch", c("Disease"))

combat.exprs <- Biobase::exprs(combat.output)

# redo PCA
cb.pca.results <- prcomp(t(combat.exprs))

# get PC1-2 as data.frame for use for ggplot2 -> scatterplot
cb.pca.df <- cbind(rownames(cb.pca.results$x),
                   dplyr::select(as.data.frame(cb.pca.results$x), PC1:PC2))
colnames(cb.pca.df) <- c("Sample", "PC1", "PC2")

# join
cb.pca.demo.df <- dplyr::full_join(rel.demo.df, cb.pca.df, by = "Sample")
cb.pca.demo.df$Batch <- as.factor(cb.pca.demo.df$Batch)

# scatterplot
combat.scatter <-
  ggplot2::ggplot(cb.pca.demo.df, ggplot2::aes(x = PC1, y = PC2, colour = Batch,
                                               shape = Disease)) +
  ggplot2::geom_point(size = 3.5) +
  ggplot2::theme_bw() +
  ggplot2::annotate("text", x = cb.pca.demo.df$PC1, y = cb.pca.demo.df$PC2,
                    label = cb.pca.demo.df$Any_Immune, colour = "#000000") +
  ggplot2::scale_color_manual(values = c("#FFB90F", "#1E90FF", "#00CD66")) +
  ggplot2::ggtitle("ComBat Adjustment", subtitle = "Disease as covariate")

# get scatterplots on same plot
pdf(file.path("NARES", "NARES_PC_scatter_with_batch_correct.pdf"),
    width = 11, height = 4.25)
gridExtra::grid.arrange(no.adjust.scatter, combat.scatter, ncol = 2)
dev.off()

# save combat adjusted data ----------------------------------------------------

# first remove trailing "_at" on Entrez ID from BrainArray & set colname
# to "Gene
combat.df <- as.data.frame(cbind(gsub("_at", "", rownames(combat.exprs)),
                                 combat.exprs))
colnames(combat.df)[1] <- "Gene"

# write to file
readr::write_tsv(combat.df, file.path(processed.dir,
                                      "NARES_SCANfast_ComBat.pcl"))
