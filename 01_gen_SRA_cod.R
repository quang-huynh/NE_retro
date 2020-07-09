
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


# Wt age
#matplot(t(args$OM@cpars$Wt_age[1, , 1:args$OM@nyears]), lty = 1, type = "o", pch = 16)

# Mat age
#matplot(t(args$OM@cpars$Mat_age[1, , 1:args$OM@nyears]), lty = 1, type = "o", pch = 16)

#compare_SRA(res, res2, filename = "compare_cod", dir = getwd(), title = "ASAP3 GM cod comparison",
#            s_name = c("NEFSCspring", "NEFSCfall", "MAspring"), 
#            scenario = list(names = c("M02", "MRAMP")),
#            open_file = FALSE)
#
#
## Profile over steepness
#h_prof <- seq(0.6, 0.99, 0.01)
#setup(12)
#sfExport(list = c('args'))
#
#M02_prof <- sfClusterApplyLB(h_prof, function(x, args) {
#  args$OM@h <- rep(x, 2)
#  do.call(SRA_scope, args)
#}, args = args)
#
#profLL <- vapply(M02_prof, function(x) ifelse(all(x@conv), x@mean_fit$opt$objective, NA), numeric(1))
#plot(h_prof, profLL)
#which.min(profLL)
#
#plot(M02_prof[[31]])
#compare_SRA(res, M02_prof[[31]])
