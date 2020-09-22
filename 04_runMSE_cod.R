
library(MSEtool)
library(mfSCA)
library(dplyr)

setup(12)
sfLibrary(mfSCA)

### Status quo - model average M02 and MRAMP
### Test 1: no/rho adjust current EMs (cod M02 and MRAMP) 
### Test 2: model averaging them
### Test 3: ref MP with no indicators
### Test 4: ref MP with indicators
### Test 5: tuning MP
run_test <- 5

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
    MSE <- runMSE(SRA_OM@OM, MPs = c("M02", "M02_ra", "MRAMP", "MRAMP_ra"), parallel = TRUE,
                  save_name = paste0("prelim_t1_NR", i))
    
    message("Finished with NR ", i, " Test 1")
    
    saveRDS(MSE, file = paste0("GoM_cod/MSE_cod_t1_NR", i, ".rds"))
    rm(MSE)
  }
  
  if(2 %in% run_test) {
    MSE <- runMSE(SRA_OM@OM, MPs = c("ma", "NFref", "FMSYref"), parallel = TRUE,
                  save_name = paste0("prelim_t2_NR", i))
    
    message("Finished with NR ", i, " Test 2")
    
    saveRDS(MSE, file = paste0("GoM_cod/MSE_cod_t2_NR", i, ".rds"))
    rm(MSE) 
  }
  
  if(3 %in% run_test) {
    SRA_OM@OM@qcv <- SRA_OM@OM@qinc <- c(0, 0)
    Eff <- seq(0.5, 3, 0.25)
    EffMP <- paste0("FMSY", Eff)
    
    FMSY <- readRDS("GoM_cod/cod_ref_pt.rds")[[2]][[i]][1, 1] %>% as.numeric()
    Map(function(x, y) assign(y, generate_Eff_MP(x, FMSY = FMSY, terminalF = SRA_OM@OM@cpars$Find[1, SRA_OM@OM@nyears]), 
                              envir = globalenv()), x = Eff, y = EffMP)
    sfExport(list = EffMP)
    
    MSE <- runMSE(SRA_OM@OM, MPs = EffMP, parallel = TRUE)
    
    message("Finished with NR ", i, " Test 3")
    
    saveRDS(MSE, file = paste0("GoM_cod/MSE_cod_t3_NR", i, ".rds"))
    rm(MSE)
  }
  
  if(4 %in% run_test) {
    SRA_OM@OM@qcv <- SRA_OM@OM@qinc <- c(0, 0)
    
    FMSY <- readRDS("GoM_cod/cod_ref_pt.rds")[[2]][[i]][1, 1] %>% as.numeric()
    `75%FMSY` <- generate_Eff_MP_with_EM(0.75, FMSY = FMSY, terminalF = SRA_OM@OM@cpars$Find[1, SRA_OM@OM@nyears], 
                                         args = args_M02, args2 = args_MRAMP)
    
    sfExport(list = "75%FMSY")
    MSE <- runMSE(SRA_OM@OM, MPs = c("75%FMSY", "FMSYref75"), parallel = TRUE)
    message("Finished with NR ", i, " Test 4")
    
    saveRDS(MSE, file = paste0("GoM_cod/MSE_cod_t4_NR", i, ".rds"))
    rm(MSE)
  }
  
  if(5 %in% run_test) {
    SRA_OM@OM@qcv <- SRA_OM@OM@qinc <- c(0, 0)
    
    FMSY <- readRDS("GoM_cod/cod_ref_pt.rds")[[2]][[i]][1, 1] %>% as.numeric()
    
    set.seed(421)
    relF <- runif(SRA_OM@OM@nsim, 0.25, 1.75)
    tuning_MP <- generate_Eff_MP_with_EM(relF, FMSY = FMSY, terminalF = SRA_OM@OM@cpars$Find[1, SRA_OM@OM@nyears], 
                                         args = args_M02, assess = FALSE)
    
    sfExport(list = "tuning_MP")
    MSE <- runMSE(SRA_OM@OM, MPs = "tuning_MP", parallel = FALSE)
    message("Finished with NR ", i, " Test 5")
    
    saveRDS(MSE, file = paste0("GoM_cod/MSE_cod_tuningMP_NR", i, ".rds"))
    rm(MSE)
  }
  
}
sfStop()

##### Merge MSE
#out <- list()
#for(i in 1:3) {
#  out[[1]] <- readRDS(paste0("GoM_cod/MSE_cod_NR", i, ".rds"))
#  out[[2]] <- readRDS(paste0("GoM_cod/MSE_cod_t4_NR", i, ".rds"))
#  out[[1]]@OM$qinc <- out[[1]]@OM$qcv <- 0
#  out[[2]]@Obs <- out[[1]]@Obs
#  res <- do.call(merge_MSE, out)
#  saveRDS(res, file = paste0("GoM_cod/MSE_cod_NR", i, ".rds"))
#}
#