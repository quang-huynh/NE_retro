
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
