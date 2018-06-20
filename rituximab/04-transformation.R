# J. Taroni 2018
# Transform RNA-seq data processed with Salmon + tximport with DESeq2
# USAGE (from rituximab/): Rscript 04-transformation.R
#

# tximport output
tximport.object <- "tximport.RDS"
gene.summary <- readRDS(tximport.object)

# read in select covariates 
covariates.df <- data.frame(readr::read_tsv("select_covariates.tsv"))

# sample names as rownames
fastq.sample.names <- readr::read_tsv("sample_list.txt",
                                      col_names = FALSE)
rownames(covariates.df) <- fastq.sample.names[[1]]
covariates.df$mainclass[which(covariates.df$comb_class == "ITN_required")] <-
  "ITN_required"

# get DESeqDataSet - w/ and w/o experimental design info
dds <- DESeq2::DESeqDataSetFromTximport(txi = gene.summary,
                                        colData = covariates.df,
                                        design = ~ procbatch + timepoint + 
                                          mainclass)
dds.blind <- DESeq2::DESeqDataSetFromTximport(txi = gene.summary,
                                              colData = covariates.df,
                                              design = ~ 1)

# variance stabilizing transformation - w/ and w/o experimental design taken
# into account
vst.design <- DESeq2::vst(object = dds, blind = FALSE)
vst.blind <- DESeq2::vst(object = dds.blind, blind = TRUE)

# extract and write transformed data to file
# with experimental design information
vst.design.assay <- SummarizedExperiment::assay(vst.design)
# we'd like to keep this as a matrix in order to keep the column names as is
vst.design.assay <- cbind(rownames(vst.design.assay), vst.design.assay)
colnames(vst.design.assay)[1] <- "Gene"
write.table(vst.design.assay, "VST_design.pcl", row.names = FALSE, 
            quote = FALSE, sep = "\t")

# blinded to experimental design
vst.blind.assay <- SummarizedExperiment::assay(vst.blind)
vst.blind.assay <- cbind(rownames(vst.blind.assay), vst.blind.assay)
colnames(vst.blind.assay)[1] <- "Gene"
write.table(vst.blind.assay, "VST_blind.pcl", row.names = FALSE, 
            quote = FALSE, sep = "\t")
