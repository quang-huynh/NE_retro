
library(MSEtool)
library(mfSCA)
library(dplyr)

setup(12)
sfLibrary(mfSCA)

### Status quo - model average M02 and MRAMP
### Test 1: no/rho adjust current EMs (cod M02 and MRAMP) 

SRA_base <- readRDS("GoM_haddock/SRA_haddock_sel9.rds")

for(i in 1:3) {
  SRA_OM <- readRDS(paste0("GoM_haddock/SRA_NR", i, ".rds")) %>% 
    generate_OM_args(AddInd_from_residuals = FALSE, AC_from_residuals = TRUE, Cobs = 0.05)
  SRA_OM@OM@cpars$Data@AddInd <- SRA_base@OM@cpars$Data@AddInd
  SRA_OM@OM@cpars$Data@CV_AddInd <- SRA_base@OM@cpars$Data@CV_AddInd
  
  args_base <- get_haddock(nsim = SRA_OM@OM@nsim) %>% generate_EM_args(SRA_EM = SRA_base)
  
  # Current models with and without rho adjust
  base <- create_MP(args_base)
  base_ra <- create_MP(args_base, rho_adjust = TRUE)
  
  sfExport(list = c("base", "base_ra"))
  
  message("Running NR ", i)
  
  MSE <- runMSE(SRA_OM@OM, MPs = c("base", "base_ra", "NFref", "FMSYref"), parallel = TRUE,
                save_name = paste0("GoM_haddock/prelim_NR", i))
  
  message("Finished with NR ", i)
  
  saveRDS(MSE, file = paste0("GoM_haddock/MSE_haddock_NR", i, ".rds"))
  rm(MSE) 
  
}
sfStop()

##### Merge MSE
#out <- list()
#for(i in 1:3) {
#  out[[1]] <- readRDS(paste0("GoM_cod/MSE_cod_t1_NR", i, ".rds"))
#  out[[2]] <- readRDS(paste0("GoM_cod/MSE_cod_t2_NR", i, ".rds"))
#  res <- do.call(merge_MSE, out)
#  saveRDS(res, file = paste0("GoM_cod/MSE_cod_NR", i, ".rds"))
#}
#