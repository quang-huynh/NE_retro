
library(MSEtool)
library(mfSCA)
library(dplyr)


SRA_M02 <- readRDS("GoM_cod/SRA_cod_M02.rds")
SRA_MRAMP <- readRDS("GoM_cod/SRA_cod_MRAMP.rds")
Cbias <- c(2.25, 1.25, 1)
create_target_MP <- function(mod, pstar = 0.5, delta = c(0.8, 1.20)) {
  
  force(mod)
  force(pstar)
  
  get_index_target <- function(yy, mod, pstar = 0.5) {
    solver_fn <- function(x) {
      predict(mod, newdata = list(Year = factor(yy), Ind_1_mu = x), type = "response") - pstar
    }
    sol <- uniroot(solver_fn, interval = c(1e-8, 15))
    log_Itarget <- sol$root
    return(log_Itarget)
  }
  
  MP <- function(x, Data, reps, AddInd = 1) {
    Irecent <- Data@AddInd[x, AddInd, (length(Data@Year)-5):length(Data@Year)] %>% log() %>% mean() %>% exp()
    Itarget <- try(get_index_target(yy = factor(2018 + length(Data@Year) - Data@LHYear + 3), 
                                    mod = mod, pstar = pstar) %>% exp(), silent = TRUE)
    
    Rec <- new("Rec")
    if(is.character(Itarget)) {
      TAC <- NA_real_
    } else {
      Cat <- Data@Cat[x, length(Data@Year)]
      
      TAC <- Cat * Irecent/Itarget
      if(TAC < delta[1] * Cat) TAC <- delta[1] * Cat
      if(TAC > delta[2] * Cat) TAC <- delta[2] * Cat
    }
    Rec@TAC <- TAC
    return(Rec)
  }
  class(MP) <- "MP"
  return(MP)
}

mod <- readRDS("GoM_cod/tuneMP_model.rds")
tuned_MP <- create_target_MP(mod)
debug(tuned_MP)
for(i in 1:3) {
  SRA_OM <- readRDS(paste0("GoM_cod/SRA_NR", i, ".rds")) %>% 
    generate_OM_args(AddInd_from_residuals = FALSE, AC_from_residuals = TRUE, Cobs = 0.05,
                     Cbias = Cbias[i])
  
  MSE <- runMSE(SRA_OM@OM, MPs = "tuned_MP")
  
  message("Finished with NR ", i)
  
  saveRDS(MSE, file = paste0("GoM_cod/MSE_cod_tuned_MP_NR", i, ".rds"))
  rm(MSE)
  
}


