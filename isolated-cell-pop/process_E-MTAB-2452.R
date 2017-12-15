# J. Taroni 2017
# This script processes E-MTAB-2452 raw data using SCANfast. It also outputs
# a PCA plot, where samples/arrays are colored by the cell type assayed. 
# Separation by cell type (CD4 [T cells], CD14 [Monocytes], CD16 [Neutrophils])
# should be evident in this plot.
# 
# USAGE: Rscript isolated-cell-pop/process_E-MTAB-2452.R

library(ggplot2)
`%>%` <- dplyr::`%>%`

source(file.path("util", "affy_norm_functions.R"))

processed.dir <- file.path("isolated-cell-pop", "processed")
dir.create(processed.dir, recursive = TRUE, showWarnings = FALSE)
  
#### process from raw ----------------------------------------------------------

AffyMultiNormWrapper(cel.dir = file.path("isolated-cell-pop", "raw"),
                     output.file.lead = file.path(processed.dir, 
                                                  "E-MTAB-2452"),
                     norm.methods = "SCAN", fast = TRUE)

#### read in processed expression data -----------------------------------------
# read in E-MTAB-2452 expression data that has been processed with SCANfast
exprs.df <- readr::read_tsv(file.path(processed.dir,
                                      "E-MTAB-2452_hugene11st_SCANfast.pcl"))
exprs.df <- as.data.frame(exprs.df)
rownames(exprs.df) <- exprs.df[, 1]
exprs.mat <- as.matrix(exprs.df[, -1])

#### PCA plot examining cell types ---------------------------------------------

# PCA to look at batch effect vs. cell type
pca.results <- prcomp(t(exprs.mat))

# get PC1-2 as data.frame for use for ggplot2 -> scatterplot
pca.df <- cbind(rownames(pca.results$x),
                dplyr::select(as.data.frame(pca.results$x), PC1:PC2))
colnames(pca.df) <- c("Sample", "PC1", "PC2")
pca.df <- pca.df %>%
  dplyr::mutate("Cell Type" = stringr::str_extract(Sample, "[^_]+"))

# scatter plot
ggplot(pca.df, aes(x = PC1, y = PC2, colour = `Cell Type`, 
                   shape = `Cell Type`)) +
  geom_point(size = 3, alpha = 0.6) +
  theme_bw() +
  labs(title = "E-MTAB-2452 PCA, Cell type") +
  theme(plot.title = element_text(hjust = 0.5, lineheight = 18,
                                  face = "bold")) +
  scale_colour_manual(values = c("#000000", "#32CD32", "#1E90EF"))
ggsave(file.path("isolated-cell-pop", "E-MTAB-2452_PCA_scatter_cell_type.pdf"))
