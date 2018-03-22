# Qiwen Hu & Jaclyn N. Taroni 2018
# Run main PLIER function on recount2 data prepared in 
# 2-prep_recount_for_plier.R. 
# Outputs a PLIER model: recount2/recount_PLIER_model.RDS
# Usage: Rscript recount2/3-run_recount_plier.R

library(PLIER)
set.seed(12345)

# read in data
plier.data.list <- readRDS(file = file.path("recount2",
                                           	"recount_data_prep_PLIER.RDS"))
# run PLIER
plierResult <- PLIER(as.matrix(plier.data.list$rpkm.cm), 
                     plier.data.list$all.paths.cm,
                     k = round((plier.data.list$k + plier.data.list$k * 0.3), 
                               0), 
                     trace = TRUE)
saveRDS(plierResult, file = file.path("recount2", "recount_PLIER_model.RDS"))
