
library(MSEtool)
library(mfSCA)
library(dplyr)

# Cod with nfleet = 1, nsurvey = 5

############## nsim = 100 for closed-loop
# M02
args <- get_cod_M02(100)
saveRDS(args, file = "GoM_cod/args_cod_M02.rds")

res <- do.call(SRA_scope, args)
ret <- retrospective(res, 7)
saveRDS(res, file = 'GoM_cod/SRA_cod_M02.rds')

plot(res, retro = ret, file = "report/GoM_cod/cod_M02", title = "ASAP3 import for GM cod M02", 
     dir = getwd(), s_name = c("NEFSCspring", "NEFSCfall", "MAspring"), open_file = FALSE)

# MRAMP
args <- get_cod_MRAMP(100)
saveRDS(args, file = "GoM_cod/args_cod_MRAMP.rds")

res <- do.call(SRA_scope, args)
ret <- retrospective(res, 7)
saveRDS(res, file = 'GoM_cod/SRA_cod_MRAMP.rds')

plot(res, retro = ret, file = "report/GoM_cod/cod_MRAMP", title = "ASAP3 import for GM cod MRAMP", 
     dir = getwd(), s_name = c("NEFSCspring", "NEFSCfall", "MAspring"), open_file = FALSE)
