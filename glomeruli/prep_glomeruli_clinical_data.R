# J. Taroni 2018
# Cleaning ERCB clinical data for use in the multi-plier project
#
# USAGE: Rscript glomeruli/prep_glomeruli_clinical_data.R

# magrittr pipe req'd for data wrangling
`%>%` <- dplyr::`%>%`

# read in clinical file to be cleaned
clinical.file <- 
  file.path("glomeruli", 
            "Neptune_ERCB_GE_Clinical_Data_2016-08-24_15-11-58-1.txt")
clinical.df <- readr::read_tsv(clinical.file)

# read in expression data because we want to use this to filter the clinical
# data
ercb.file <- file.path("glomeruli", "ERCB_Glom_CustCDF19_forVCRC.txt")
exprs.df <- readr::read_tsv(ercb.file)

# only retain samples that are in the expression data -- the first two columns
# are gene identifiers
microarray.samples <- colnames(exprs.df)[3:ncol(exprs.df)]
if (!all(microarray.samples %in% clinical.df$`Microarray Sample ID`)) {
  stop("Not all samples are in expression data and clinical data")
}

# we only want the info for the samples in the ERCB glomeruli data -- that's
# the microarray data we're looking at
clinical.df <- clinical.df %>% 
  dplyr::filter(`Microarray Sample ID` %in% microarray.samples,
                `Tissue Source` == "Glomeruli")

# diagnosis information only - this will serve as group labels
# recode diagnosis such that nephrotic syndrome diagnoses are grouped
diagnosis.df <- clinical.df %>%
  dplyr::select(c(`Microarray Sample ID`, Diagnosis)) %>%
  dplyr::mutate(Diagnosis = 
                  dplyr::recode(Diagnosis,
                                MCD = "Nephrotic syndrome",
                                MN = "Nephrotic syndrome",
                                FSGS = "Nephrotic syndrome",
                                `FSGS/MCD` = "Nephrotic syndrome"),
                Sample = `Microarray Sample ID`) %>%
  dplyr::select(c(Sample, Diagnosis)) %>%
  readr::write_tsv(path = file.path("glomeruli", "ERCB_glom_diagnosis.tsv"))
