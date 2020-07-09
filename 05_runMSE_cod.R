
library(MSEtool)
library(mfSCA)
library(dplyr)

setup(12)
sfLibrary(mfSCA)

### Status quo - model average M02 and MRAMP
### Test 1: no/rho adjust current EMs (cod M02 and MRAMP) 
### Test 2: model averaging them
### Test 3: incorrect EM 
run_test <- c(1, 2)

SRA_M02 <- readRDS("GoM_cod/SRA_cod_M02.rds")
SRA_MRAMP <- readRDS("GoM_cod/SRA_cod_MRAMP.rds")
Cbias <- c(2.25, 1.25, 1)

for(i in 1:3) {
  SRA_OM <- readRDS(paste0("GoM_cod/SRA_NR", i, ".rds")) %>% 
    generate_OM_args(AddInd_from_residuals = FALSE, AC_from_residuals = TRUE, Cobs = 0.05,
                     Cbias = Cbias[i])
  
  args_M02 <- get_cod_M02(nsim = SRA_OM@OM@nsim) %>% 
    generate_EM_args(SRA_EM = SRA_M02, Cbias = Cbias[i])
  args_MRAMP <- get_cod_MRAMP(nsim = SRA_OM@OM@nsim) %>% 
    generate_EM_args(SRA_EM = SRA_MRAMP, Cbias = Cbias[i])
  
  # Current models with and without rho adjust
  M02 <- create_MP(args_M02)
  M02_ra <- create_MP(args_M02, rho_adjust = TRUE)
  MRAMP <- create_MP(args_MRAMP)
  MRAMP_ra <- create_MP(args_MRAMP, rho_adjust = TRUE)
  
  # Model average in lieu of rho adjust
  ma <- create_MP(args_M02, args_MRAMP)
  
  sfExport(list = c("M02", "M02_ra", "MRAMP", "MRAMP_ra", "ma"))
  
  message("Running NR ", i)
  
  if(1 %in% run_test) {
    MSE <- runMSE(SRA_OM@OM, MPs = c("M02", "M02_ra", "MRAMP", "MRAMP_ra"), parallel = TRUE)
    
    message("Finished with NR ", i, " Test 1")
    
    saveRDS(MSE, file = paste0("GoM_cod/MSE_cod_t1_NR", i, ".rds"), save_name = paste0("prelim_t1_NR", i))
    rm(MSE)
  }
  
  if(2 %in% run_test) {
    MSE <- runMSE(SRA_OM@OM, MPs = c("ma", "NFref", "FMSYref"), parallel = TRUE)
    
    message("Finished with NR ", i, " Test 2")
    
    saveRDS(MSE, file = paste0("GoM_cod/MSE_cod_t2_NR", i, ".rds"), save_name = paste0("prelim_t2_NR", i))
    rm(MSE) 
  }
  
}
sfStop()

