
library(MSEtool)
library(mfSCA)
library(dplyr)

setup(12)
sfLibrary(mfSCA)

SRA_M02 <- readRDS("GoM_cod/SRA_cod_M02.rds")
SRA_MRAMP <- readRDS("GoM_cod/SRA_cod_MRAMP.rds")
Cbias <- c(2.25, 1.25, 1)

OM <- c("MC", "MCIM", "IM")

for(i in 1:3) {
  SRA_OM <- readRDS(paste0("GoM_cod/SRA_", OM[i], ".rds")) %>% 
    generate_OM_args(AddInd_from_residuals = FALSE, AC_from_residuals = TRUE, Cobs = 0.05,
                     Cbias = Cbias[i])
  SRA_OM@OM@qcv <- SRA_OM@OM@qinc <- c(0, 0)
  
  args_M02 <- get_cod_M02(nsim = SRA_OM@OM@nsim) %>% 
    generate_EM_args(SRA_EM = SRA_M02, Cbias = Cbias[i])
  args_MRAMP <- get_cod_MRAMP(nsim = SRA_OM@OM@nsim) %>% 
    generate_EM_args(SRA_EM = SRA_MRAMP, Cbias = Cbias[i])
  
  # Current models with and without rho adjust
  M02_nra <- create_MP(args_M02)
  M02_ra <- create_MP(args_M02, rho_adjust = TRUE)
  MRAMP_nra <- create_MP(args_MRAMP)
  MRAMP_ra <- create_MP(args_MRAMP, rho_adjust = TRUE)
  
  # Model average
  MA <- create_MP(args_M02, args_MRAMP)
  
  # Reference MP
  FMSY <- readRDS("GoM_cod/cod_ref_pt.rds")[[2]][[i]][1, 1] %>% as.numeric()
  `75%FMSY` <- generate_Eff_MP_with_EM(0.75, FMSY = FMSY, terminalF = SRA_OM@OM@cpars$Find[1, SRA_OM@OM@nyears], 
                                       args = args_M02, args2 = args_MRAMP)
  
  # Export MPs to cores
  sfExport(list = c("M02_nra", "M02_ra", "MRAMP_nra", "MRAMP_ra", "MA", "75%FMSY"))
  
  message("Running OM ", i)
  
  # Batch 1
  MSE <- runMSE(SRA_OM@OM, MPs = c("M02_nra", "M02_ra", "MRAMP_nra", "MRAMP_ra"), parallel = TRUE)
  
  message("Finished with OM ", i, " Batch 1")
  saveRDS(MSE, file = paste0("GoM_cod/MSE_cod_t1_", OM[i], ".rds"))
  rm(MSE)
  
  # Batch 2
  MSE <- runMSE(SRA_OM@OM, MPs = c("MA", "NFref", "FMSYref", "75%FMSY"), parallel = TRUE,
                save_name = paste0("prelim_t2_NR", i))
  
  message("Finished with OM ", i, " Batch 2")
  saveRDS(MSE, file = paste0("GoM_cod/MSE_cod_t2_", OM[i], ".rds"))
  
}
sfStop()

##### Merge MSE
out <- list()
for(i in 1:3) {
  out[[1]] <- readRDS(paste0("GoM_cod/MSE_cod_t1_", OM[i], ".rds"))
  out[[2]] <- readRDS(paste0("GoM_cod/MSE_cod_t2_", OM[i], ".rds"))
  out[[2]]@Obs <- out[[1]]@Obs
  res <- do.call(merge_MSE, out)
  saveRDS(res, file = paste0("GoM_cod/MSE_cod_", OM[i], ".rds"))
}


