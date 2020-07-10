library(MSEtool)
library(mfSCA)
library(dplyr)

# Base
args <- get_pollock_base(100)

res <- do.call(SRA_scope, args)
ret <- retrospective(res, 7)

saveRDS(res, file = "pollock/SRA_pollock_base.rds")

plot(res, retro = ret, file = "pollock_base", title = "ASAP3 import for pollock", dir = "report/pollock", 
     s_name = c("NEFSCspring", "NEFSCfall"), f_name = c("Commercial", "Recreational"), open_file = FALSE)

# Flat sel
args <- get_pollock_flatsel(100)

res <- do.call(SRA_scope, args)
ret <- retrospective(res, 7)

saveRDS(res, file = "pollock/SRA_pollock_flatsel.rds")

plot(res, retro = ret, file = "pollock_flatsel", title = "ASAP3 import for pollock", dir = "report/pollock", 
     s_name = c("NEFSCspring", "NEFSCfall"), f_name = c("Commercial", "Recreational"), open_file = FALSE)

