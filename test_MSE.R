
library(MSEtool)
library(dplyr)
library(mfSCA)
### Cod

#Example with OM = M02, EM = MRAMP

Cbias <- 2

SRA_M02 <- do.call(SRA_scope, get_cod_M02())
SRA_MRAMP <- do.call(SRA_scope, get_cod_MRAMP())

#args_M02 <- get_cod_M02() %>% generate_EM_args(SRA_EM = SRA_M02, Cbias = Cbias[i])
args_MRAMP <- get_cod_MRAMP() %>% generate_EM_args(SRA_EM = SRA_MRAMP, Cbias = Cbias)

SRA_OM <- SRA_M02 %>% generate_OM_args(AddInd_from_residuals = FALSE, AC_from_residuals = TRUE, Cobs = 0.05)

MP_MRAMP <- create_MP(args_MRAMP, rho_adjust = TRUE)
undebug(MP_MRAMP)

MSE <- runMSE(SRA_OM@OM, MPs = c("MP_MRAMP"))


# Haddock
# Example with OM = EM
SRA_base <- do.call(SRA_scope, get_haddock())
args_base <- get_haddock() %>% generate_EM_args(SRA_EM = SRA_base)

SRA_OM <- SRA_base %>% 
    generate_OM_args(AddInd_from_residuals = FALSE, AC_from_residuals = TRUE, Cobs = 0.05)
  
# Current models with and without rho adjust
base <- create_MP(args_base)
base_ra <- create_MP(args_base, rho_adjust = TRUE)
  
  
MSE <- runMSE(SRA_OM@OM, MPs = c("base"))

res <- mf_SCA(Data = MSE@Misc$Data[[1]], args = args_base)

# Pollock
# Example with OM = base, EM = flatsel

SRA_base <- do.call(SRA_scope, get_pollock_base())
SRA_flatsel <- do.call(SRA_scope, get_pollock_flatsel())
args_flatsel <- get_pollock_flatsel() %>% generate_EM_args(SRA_EM = SRA_flatsel, SRA_OM = SRA_base)

SRA_OM <- SRA_base %>% generate_OM_args(AddInd_from_residuals = FALSE, AC_from_residuals = TRUE, Cobs = 0.05)

flatsel <- create_MP(args_flatsel, calc_rho = FALSE)
#debug(mf_SCA)

MSE <- runMSE(SRA_OM@OM, MPs = "flatsel")
