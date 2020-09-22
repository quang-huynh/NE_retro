
library(MSEtool)
library(mfSCA)
library(dplyr)

setup(12)
sfLibrary(mfSCA)

### Status quo - model average base and flatsel
### Test 1: no/rho adjust current EMs
### Test 2: model averaging them
### Test 3: ref MP with no indicators
### Test 4: ref MP with indicators
### Test 5: tuning MP
run_test <- 5

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
  
  if(3 %in% run_test) {
    SRA_OM@OM@qcv <- SRA_OM@OM@qinc <- c(0, 0)
    Eff <- seq(0.5, 3, 0.25)
    EffMP <- paste0("FMSY", Eff)
    
    FMSY <- readRDS("pollock/pollock_ref_pt.rds")[[i]][1, 1] %>% as.numeric()
    Map(function(x, y) assign(y, generate_Eff_MP(x, FMSY = FMSY, terminalF = SRA_OM@OM@cpars$Find[1, SRA_OM@OM@nyears]), 
                              envir = globalenv()), x = Eff, y = EffMP)
    sfExport(list = EffMP)
    
    MSE <- runMSE(SRA_OM@OM, MPs = EffMP, parallel = TRUE)
    
    message("Finished with NR ", i, " Test 3")
    
    saveRDS(MSE, file = paste0("pollock/MSE_pollock_t3_NR", i, ".rds"))
    rm(MSE)
  }
  
  if(4 %in% run_test) {
    SRA_OM@OM@qcv <- SRA_OM@OM@qinc <- c(0, 0)
    
    FMSY <- readRDS("pollock/pollock_ref_pt.rds")[[i]][1, 1] %>% as.numeric()
    `75%FMSY` <- generate_Eff_MP_with_EM(0.75, FMSY = FMSY, terminalF = SRA_OM@OM@cpars$Find[1, SRA_OM@OM@nyears], 
                                         args = args_base, args2 = args_flatsel)
    
    sfExport(list = "75%FMSY")
    MSE <- runMSE(SRA_OM@OM, MPs = c("75%FMSY", "FMSYref75"), parallel = TRUE)
    message("Finished with NR ", i, " Test 4")
    
    saveRDS(MSE, file = paste0("pollock/MSE_pollock_t4_NR", i, ".rds"))
    rm(MSE)
  }
  
  if(5 %in% run_test) {
    SRA_OM@OM@qcv <- SRA_OM@OM@qinc <- c(0, 0)
    
    FMSY <- readRDS("pollock/pollock_ref_pt.rds")[[i]][1, 1] %>% as.numeric()
    
    set.seed(315)
    relF <- runif(SRA_OM@OM@nsim, 0.25, 3)
    tuning_MP <- generate_Eff_MP_with_EM(relF, FMSY = FMSY, terminalF = SRA_OM@OM@cpars$Find[1, SRA_OM@OM@nyears], 
                                         args = args_base, assess = FALSE)
    
    sfExport(list = "tuning_MP")
    MSE <- runMSE(SRA_OM@OM, MPs = "tuning_MP", parallel = TRUE)
    message("Finished with NR ", i, " Test 5")
    
    saveRDS(MSE, file = paste0("pollock/MSE_pollock_tuningMP_NR", i, ".rds"))
    rm(MSE)
  }
  
}
sfStop()

##### Merge MSE
#out <- list()
#for(i in 1:3) {
#  out[[1]] <- readRDS(paste0("pollock/MSE_pollock_NR", i, ".rds"))
#  out[[2]] <- readRDS(paste0("pollock/MSE_pollock_t4_NR", i, ".rds"))
#  out[[1]]@OM$qcv <- out[[1]]@OM$qinc <- 0
#  out[[2]]@Obs <- out[[1]]@Obs
#  res <- do.call(merge_MSE, out)
#  saveRDS(res, file = paste0("pollock/MSE_pollock_NR", i, ".rds"))
#}
#

