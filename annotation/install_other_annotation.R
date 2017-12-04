# J. Taroni 2017
# Install other annotation packages from bioconductor for use with the rheum
# PLIER datasets.
# USAGE: Rscript annotation/install_other_annotation.R

#### main ----------------------------------------------------------------------

# bioconductor installer
source("https://bioconductor.org/biocLite.R")

# Illumina HT-12 V4.0 BeadChip (E-GEOD-49454, E-GEOD-65391)
biocLite("illuminaHumanv4.db", suppressUpdates = TRUE)

# Agilent 4x44K G4112
biocLite("hgug4112a.db", suppressUpdates = TRUE)

# platform design information for Affymetrix
biocLite(c("affxparser", "pd.hg.u133.plus.2", "pd.hg.u133a",
           "pd.hg.u133b", "pd.hg.u95a", "pd.hugene.1.0.st.v1",
           "pd.hugene.1.1.st.v1"), suppressUpdates = TRUE)

lapply(list("affxparser", "pd.hg.u133.plus.2", "pd.hg.u133a",
            "pd.hg.u133b", "pd.hg.u95a", "pd.hugene.1.0.st.v1",
            "pd.hugene.1.1.st.v1", "illuminaHumanv4.db",
            "hgug4112a.db"), 
       function(x) library(x, character.only = TRUE))

sink(file.path("annotation", "other_annotation_sessionInfo.txt"))
sessionInfo()
sink()
