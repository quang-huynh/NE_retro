library(mfSCA)

# Cod with nfleet = 1, nsurvey = 2
args <- get_haddock(nsim = 100)

res <- do.call(SRA_scope, args)
ret <- retrospective(res, 7)

saveRDS(res, file = "GoM_haddock/SRA_haddock_sel9.rds")

plot(res, retro = ret, file = "report/GoM_haddock/haddock_sel9", title = "ASAP3 import for GM haddock", 
     dir = getwd(), s_name = c("NEFSCspring", "NEFSCfall"), open_file = FALSE)


# Wt age
matplot(t(args$OM@cpars$Wt_age[1, , 1:args$OM@nyears]), lty = 1, type = "o", pch = 16)

# Mat age
matplot(t(args$OM@cpars$Mat_age[1, , 1:args$OM@nyears]), lty = 1, type = "o", pch = 16)
