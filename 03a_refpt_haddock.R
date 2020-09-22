library(MSEtool)
library(mfSCA)

# Calc OM ref points with original M or new M
SRA <- lapply(1:3, function(i) readRDS(paste0("GoM_haddock/SRA_NR", i, ".rds")))
med_rec <- vapply(SRA, function(i) quantile(i@mean_fit$report$R, 0.5), numeric(1))

######## Reference points
ref_pt <- list()
ref_pt[[1]] <- lapply(1:3, function(i) calc_refpt(SRA[[i]], med_rec = med_rec[i]))
ref_pt[[2]] <- lapply(1:3, function(i) calc_refpt(SRA[[i]], M = "true", med_rec = med_rec[i]))
saveRDS(ref_pt, file = "GoM_haddock/haddock_ref_pt.rds")

