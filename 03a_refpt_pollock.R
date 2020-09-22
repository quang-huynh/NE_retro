library(MSEtool)
library(mfSCA)

# Calc ref points (M = 0.2)
SRA <- lapply(1:3, function(i) readRDS(paste0("pollock/SRA_NR", i, ".rds")))
med_rec <- vapply(SRA, function(i) quantile(i@mean_fit$report$R, 0.5), numeric(1))

######## Reference points
ref_pt <- lapply(1:3, function(i) calc_refpt(SRA[[i]], med_rec = med_rec[i]))
saveRDS(ref_pt, file = "pollock/pollock_ref_pt.rds")

