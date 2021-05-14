library(MSEtool)
library(mfSCA)

# Calc OM ref points with M = 0.20 or new M
SRA <- lapply(1:3, function(i) readRDS(paste0("GoM_cod/SRA_NR", i, ".rds")))
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

lapply(AM_ref[[2]], function(x) x[1, ])




i = 1
RR <- SRA[[i]]@mean_fit$report$R[2:38]
SS <- SRA[[i]]@mean_fit$report$E[1:37]
Rp <- SRA[[i]]@mean_fit$report$Arec * SS / (1 + SRA[[i]]@mean_fit$report$Brec * SS)

plot(SS, RR, xlim = c(0, 60000),  ylim = c(0, 40000))
lines(SS, Rp)
text(SS, RR, labels = 2019 - 36:0, pos = 3)
abline(h = 0, v = c(0, 0.2 * SRA[[i]]@mean_fit$report$E0_SR), col = 'grey')


plot(SS, c(SRA[[i]]@mean_fit$report$log_rec_dev[-1], 0), typ = 'o')
abline(h = 0, lty = 2)


