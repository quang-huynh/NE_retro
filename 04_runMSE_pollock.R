
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

SRA_base <- readRDS("pollock/SRA_pollock_base.rds")
SRA_flatsel <- readRDS("pollock/SRA_pollock_flatsel.rds")

for(i in 1:3) {
  SRA_OM <- readRDS(paste0("pollock/SRA_NR", i, ".rds")) %>% 
    generate_OM_args(AddInd_from_residuals = FALSE, AC_from_residuals = TRUE, Cobs = 0.05)
  
  args_base <- get_pollock_base(nsim = SRA_OM@OM@nsim) %>% 
    generate_EM_args(SRA_EM = SRA_base, SRA_OM = SRA_OM)
  args_flatsel <- get_pollock_flatsel(nsim = SRA_OM@OM@nsim) %>% 
    generate_EM_args(SRA_EM = SRA_flatsel, SRA_OM = SRA_OM)
  
  # Current models with and without rho adjust
  base <- create_MP(args_base)
  base_ra <- create_MP(args_base, rho_adjust = TRUE)
  flatsel <- create_MP(args_flatsel)
  flatsel_ra <- create_MP(args_flatsel, rho_adjust = TRUE)
  
  # Model average in lieu of rho adjust
  ma <- create_MP(args_base, args_flatsel)
  
  sfExport(list = c("base", "base_ra", "flatsel", "flatsel_ra", "ma"))
  
  message("Running NR ", i)
  
  if(1 %in% run_test) {
    MSE <- runMSE(SRA_OM@OM, MPs = c("base", "base_ra", "flatsel", "flatsel_ra"), parallel = TRUE,
                  save_name = paste0("pollock/prelim_t1_NR", i))
    
    message("Finished with NR ", i, " Test 1")
    
    saveRDS(MSE, file = paste0("pollock/MSE_pollock_t1_NR", i, ".rds"))
    rm(MSE)
  }
  
  if(2 %in% run_test) {
    MSE <- runMSE(SRA_OM@OM, MPs = c("ma", "NFref", "FMSYref"), parallel = TRUE, 
                  save_name = paste0("pollock/prelim_t2_NR", i))
    
    message("Finished with NR ", i, " Test 2")
    
    saveRDS(MSE, file = paste0("pollock/MSE_pollock_t2_NR", i, ".rds"))
    rm(MSE) 
  }
  
}
sfStop()

##### Merge MSE
out <- list()
for(i in 1:3) {
  out[[1]] <- readRDS(paste0("pollock/MSE_pollock_t1_NR", i, ".rds"))
  out[[2]] <- readRDS(paste0("pollock/MSE_pollock_t2_NR", i, ".rds"))
  res <- do.call(merge_MSE, out)
  saveRDS(res, file = paste0("pollock/MSE_pollock_NR", i, ".rds"))
}


