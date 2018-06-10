# J. Taroni 2018
# Summarize data to the gene-level with tximport

library(ensembldb)

StripTrailing <- function(filename) {
  # in order to get the tx2gene data.frame ENST ids to match, strip ".1", etc.
  df <- readr::read_tsv(filename)
  df$Name <- sub("\\..*", "", df$Name)
  readr::write_tsv(df, path = filename)
}

# derive transcript to gene mapping from gtf
ensembldb::ensDbFromGtf("Homo_sapiens.GRCh38.92.gtf")
edb <- ensembldb::EnsDb("Homo_sapiens.GRCh38.92.sqlite")
tx <- ensembldb::transcriptsBy(edb, by = "gene")
tx.df <- as.data.frame(tx@unlistData)
tx2gene <- tx.df[, c("tx_name", "gene_id")]

# vector of quant.sf files
sf.files <- list.files("quants", recursive = TRUE, full.names = TRUE,
                       pattern = "quant.sf")
sapply(sf.files, StripTrailing)

# main tximport function
txi <- tximport::tximport(sf.files, type = "salmon", tx2gene = tx2gene)
saveRDS(txi, file = "tximport.RDS")
