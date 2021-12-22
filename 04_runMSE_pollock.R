
library(MSEtool)
library(mfSCA)
library(dplyr)

setup(12)
sfLibrary(mfSCA)

SRA_base <- readRDS("pollock/SRA_pollock_base.rds")
SRA_flatsel <- readRDS("pollock/SRA_pollock_flatsel.rds")

OM <- c("SS", "SWB", "SWF")

for(i in 1:3) {
  SRA_OM <- readRDS(paste0("pollock/SRA_", OM[i], ".rds")) %>% 
    generate_OM_args(AddInd_from_residuals = FALSE, AC_from_residuals = TRUE, Cobs = 0.05)
  SRA_OM@OM@qcv <- SRA_OM@OM@qinc <- c(0, 0)
  
  args_base <- get_pollock_base(nsim = SRA_OM@OM@nsim) %>% 
    generate_EM_args(SRA_EM = SRA_base, SRA_OM = SRA_OM)
  args_flatsel <- get_pollock_flatsel(nsim = SRA_OM@OM@nsim) %>% 
    generate_EM_args(SRA_EM = SRA_flatsel, SRA_OM = SRA_OM)
  
  # Current models with and without rho adjust
  Base_nra <- create_MP(args_base)
  Base_ra <- create_MP(args_base, rho_adjust = TRUE)
  Flatsel_nra <- create_MP(args_flatsel)
  Flatsel_ra <- create_MP(args_flatsel, rho_adjust = TRUE)
  
  # Model average
  MA <- create_MP(args_base, args_flatsel)
  
  # Reference MP
  FMSY <- readRDS("pollock/pollock_ref_pt.rds")[[i]][1, 1] %>% as.numeric()
  `75%FMSY` <- generate_Eff_MP_with_EM(0.75, FMSY = FMSY, terminalF = SRA_OM@OM@cpars$Find[1, SRA_OM@OM@nyears], 
                                       args = args_base, args2 = args_flatsel)
  
  # Export MPs to cores
  sfExport(list = c("Base_nra", "Base_ra", "Flatsel_nra", "Flatsel_ra", "MA", "75%FMSY"))
  
  message("Running OM ", i)
  
  # Batch 1
  MSE <- runMSE(SRA_OM@OM, MPs = c("Base_nra", "Base_ra", "Flatsel", "Flatsel_ra"), parallel = TRUE)
  
  message("Finished with NR ", i, " Batch 1")
  
  saveRDS(MSE, file = paste0("pollock/MSE_pollock_t1_", OM[i], ".rds"))
  rm(MSE)
  
  # Batch 2
  MSE <- runMSE(SRA_OM@OM, MPs = c("MA", "NFref", "75%FMSY"), parallel = TRUE)
  
  message("Finished with NR ", i, " Batch 2")
  
  saveRDS(MSE, file = paste0("pollock/MSE_pollock_t2_", OM[i], ".rds"))
  rm(MSE) 
  
}
sfStop()

##### Merge MSE
out <- list()
for(i in 1:3) {
  out[[1]] <- readRDS(paste0("pollock/MSE_pollock_t1_", OM[i], ".rds"))
  out[[2]] <- readRDS(paste0("pollock/MSE_pollock_t2_", OM[i], ".rds"))
  out[[2]]@Obs <- out[[1]]@Obs
  res <- do.call(merge_MSE, out)
  saveRDS(res, file = paste0("pollock/MSE_pollock_", OM[i], ".rds"))
}


