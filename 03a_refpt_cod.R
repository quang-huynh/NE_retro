library(MSEtool)
library(mfSCA)

# Calc OM ref points with M = 0.20 or new M
OM <- c("MC", "MCIM", "IM")

SRA <- lapply(c("MC", "MCIM", "IM"), function(i) readRDS(paste0("GoM_cod/SRA_", i, ".rds")))
med_rec <- vapply(SRA, function(i) quantile(i@mean_fit$report$R, 0.75), numeric(1))

######## Reference points
ref_pt <- list()
ref_pt[[1]] <- lapply(1:3, function(i) calc_refpt(SRA[[i]], med_rec = med_rec[i]))
ref_pt[[2]] <- lapply(1:3, function(i) calc_refpt(SRA[[i]], M = "true", med_rec = med_rec[i]))
saveRDS(ref_pt, file = "GoM_cod/cod_ref_pt.rds")

##### Compensation ratio - original
lapply(1:3, function(i) SRA[[i]]@mean_fit$report$Arec * mean(SRA[[i]]@mean_fit$report$EPR0[1:3]))
lapply(1:3, function(i) SRA[[i]]@mean_fit$report$Arec * mean(SRA[[i]]@mean_fit$report$EPR0[37]))


# Calculate F40% for assessment models
SRA <- lapply(c("M02", "MRAMP"), function(i) readRDS(paste0("GoM_cod/SRA_cod_", i, ".rds")))
med_rec <- vapply(SRA, function(i) quantile(i@mean_fit$report$R, 0.75), numeric(1))
AM_ref <- list()
AM_ref[[1]] <- lapply(1:2, function(i) calc_refpt(SRA[[i]], med_rec = med_rec[i]))
AM_ref[[2]] <- lapply(1:2, function(i) calc_refpt(SRA[[i]], M = "true", med_rec = med_rec[i]))
