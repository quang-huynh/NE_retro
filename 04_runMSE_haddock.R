
library(MSEtool)
library(mfSCA)
library(dplyr)

setup(12)
sfLibrary(mfSCA)

### Status quo
### Test 1: no/rho adjust current EM
### Tets 5: tuning MP
run_test <- 4

SRA_base <- readRDS("GoM_haddock/SRA_haddock_sel9.rds")

for(i in 1:3) {
  SRA_OM <- readRDS(paste0("GoM_haddock/SRA_NR", i, ".rds")) 
  if(i == 1) SRA_OM@OM@cpars$Data@AddInd <- SRA_base@OM@cpars$Data@AddInd
  if(i == 2) SRA_OM@OM@cpars$Data@CV_AddInd <- SRA_base@OM@cpars$Data@CV_AddInd
  SRA_OM <- generate_OM_args(SRA_OM, AddInd_from_residuals = FALSE, AC_from_residuals = TRUE, Cobs = 0.05)
  
  args_base <- get_haddock(nsim = SRA_OM@OM@nsim) %>% generate_EM_args(SRA_EM = SRA_base)
  
  # Current models with and without rho adjust
  base <- create_MP(args_base)
  base_ra <- create_MP(args_base, rho_adjust = TRUE)
  
  sfExport(list = c("base", "base_ra"))
  
  message("Running NR ", i)
  
  if(1 %in% run_test) {
    MSE <- runMSE(SRA_OM@OM, MPs = c("base", "base_ra", "NFref", "FMSYref"), parallel = TRUE,
                  save_name = paste0("GoM_haddock/prelim_NR", i))
    
    message("Finished with NR ", i)
    
    saveRDS(MSE, file = paste0("GoM_haddock/MSE_haddock_NR", i, ".rds"))
    rm(MSE)
  }
  
  
  if(4 %in% run_test) {
    SRA_OM@OM@qcv <- SRA_OM@OM@qinc <- c(0, 0)
    
    FMSY <- readRDS("GoM_haddock/haddock_ref_pt.rds")[[2]][[i]][1, 1] %>% as.numeric()
    `75%FMSY` <- generate_Eff_MP_with_EM(0.75, FMSY = FMSY, terminalF = SRA_OM@OM@cpars$Find[1, SRA_OM@OM@nyears], 
                                         args = args_base, assess = FALSE)
    
    sfExport(list = "75%FMSY")
    MSE <- runMSE(SRA_OM@OM, MPs = c("75%FMSY", "FMSYref75"), parallel = FALSE)
    message("Finished with NR ", i, " Test 4")
    
    saveRDS(MSE, file = paste0("GoM_haddock/MSE_haddock_t4_NR", i, ".rds"))
    rm(MSE)
  }
  
  if(5 %in% run_test) {
    SRA_OM@OM@qcv <- SRA_OM@OM@qinc <- c(0, 0)
    
    FMSY <- readRDS("GoM_haddock/haddock_ref_pt.rds")[[2]][[i]][1, 1] %>% as.numeric()
    
    set.seed(98273)
    relF <- runif(SRA_OM@OM@nsim, 0.25, 1.75)
    tuning_MP <- generate_Eff_MP_with_EM(relF, FMSY = FMSY, terminalF = SRA_OM@OM@cpars$Find[1, SRA_OM@OM@nyears], 
                                         args = args_base, assess = FALSE)
    
    sfExport(list = "tuning_MP")
    MSE <- runMSE(SRA_OM@OM, MPs = "tuning_MP", parallel = FALSE)
    message("Finished with NR ", i, " Test 5")
    
    saveRDS(MSE, file = paste0("GoM_haddock/MSE_haddock_tuningMP_NR", i, ".rds"))
    rm(MSE)
  }
  
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