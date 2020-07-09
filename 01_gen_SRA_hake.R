library(mfSCA)

# Cod with nfleet = 1, nsurvey = 2
get_hake <- function(nsim = 2, h = 0.99) {
  args <- ASAP2SRA("ASAP/2019_HKW_UNIT_MODEL_BASE.dat", nsim = nsim)
  args$OM@h <- rep(h, 2)
  
  args$vul_par <- matrix(c(0.003, 0.3, 0.8, rep(0.999, 6)), 9, 2)
  
  map_vul_par <- matrix(NA, 9, 2)
  map_vul_par[1:5, ] <- 1:10
  args$map_vul_par <- map_vul_par
  
  sel_par <- matrix(c(0.01, 0.5, rep(0.999999, 7)), 9, 2)
  
  m_sel_par <- matrix(NA, 9, 2)
  m_sel_par[-3, ] <- 1:16
  
  args$s_vul_par <- sel_par
  args$map_s_vul_par <- m_sel_par
  
  return(args)
}

args <- get_hake(100)

res <- do.call(SRA_scope, args)
ret <- retrospective(res, 7)

saveRDS(res, file = "hake/SRA_hake_base.rds")

plot(res, retro = ret, file = "report/hake/hake_base", title = "ASAP3 import for white hake", 
     dir = getwd(), s_name = c("NEFSCspring", "NEFSCfall"), open_file = FALSE)


# Wt age
matplot(t(args$OM@cpars$Wt_age[1, , 1:args$OM@nyears]), lty = 1, type = "o", pch = 16)

# Mat age
matplot(t(args$OM@cpars$Mat_age[1, , 1:args$OM@nyears]), lty = 1, type = "o", pch = 16)
