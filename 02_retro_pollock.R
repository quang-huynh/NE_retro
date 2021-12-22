
library(MSEtool)
library(dplyr)
library(mfSCA)

# Commerical sel 2-4 at maxage to 0.6
args <- get_pollock_base(h = 0.81)

hist_com_sel <- function(x, sel = 2:4, args) {
  args$vul_par[9, sel] <- x
  args$map_vul_par[9, sel] <- NA
  
  res <- do.call(SRA_scope, args)
  ret <- retrospective(res, 7)
  return(list(res, ret))
}

setup(6)
sfExport(list = "args")

com_sel <- seq(0.2, 1, 0.1)
com_sel_prof <- sfClusterApplyLB(com_sel, hist_com_sel, sel = 2:4, args = args)

rho <- vapply(com_sel_prof, function(xx) summary(xx[[2]])[3, 1], numeric(1))
plot(com_sel, rho)
plot(com_sel_prof[[5]][[2]])

lapply(com_sel_prof, function(x) MSEtool:::ilogit(x[[1]]@mean_fit$report$vul_par[, 1:5]))
lapply(com_sel_prof, function(x) x[[1]]@mean_fit$report$vul[, , 1])

#saveRDS(args, file = "pollock/args_flatcomsel.rds")
#saveRDS(res, file = "pollock/pollock_flatcomsel.rds")
#saveRDS(ret, file = "pollock/ret_pollock_flatcomsel.rds")
#
#plot(res, retro = ret, file = "pollock_flatcomsel", title = "ASAP3 import for pollock", dir = "pollock", 
#     s_name = c("NEFSCspring", "NEFSCfall"), f_name = c("Commercial", "Recreational"))
#

## Eval retro
eval_retro <- function(x, type = c("survey_sel", "catch_com", "catch_rec", 
                                   "I_sd", "CAA_mult", "s_CAA_mult", "I_mult", "M_mult"), 
                       change_com_sel = TRUE, M_mult = 1, years = NULL, 
                       est_s_sel = FALSE, com_sel = 0.5, args) {
  require(dplyr)
  type <- match.arg(type, several.ok = TRUE)
  if(est_s_sel) {
    args$map_s_vul_par <- matrix(c(1:5, NA, NA, 6:12, NA, NA, 13:14), 9, 2)
  }
  
  if(!is.na(com_sel)) {
    args$vul_par[9, 2:4] <- com_sel
    args$map_vul_par[9, 2:4] <- NA
  }
  
  for(i in 1:length(type)) {
    if (type[i] == "survey_sel") {
      args$s_vul_par[9, ] <- x[i] 
    }
    if (type[i] == "catch_com") {
      args$data$Chist[years[[i]], 1] <- x[i] * args$data$Chist[years[[i]], 1]
    }
    if (type[i] == "catch_rec") {
      args$data$Chist[years[[i]], 2] <- x[i] * args$data$Chist[years[[i]], 2]
    }
    if (type[i] == "I_mult") {
      args$data$Index[years[[i]], ] <- x[i] * args$data$Index[years[[i]], ]
    }
    if (type[i] == "I_sd") {
      #args$data$I_sd <- x * args$data$I_sd 
      args$LWT <- list(Index = x[i])
    }
    if (type[i] == "CAA_mult") {
      args$LWT <- list(CAA = x[i]) 
    }
    if (type[i] == "s_CAA_mult") {
      args$LWT <- list(s_CAA = x[i]) 
    }
    if (type[i] == "M_mult") {
      args$OM@cpars$M_ageArray[, , years[[i]]] <- x[i] * args$OM@cpars$M_ageArray[, , years[[i]]]
    }
  }
  args$OM@cpars$M_ageArray <- M_mult * args$OM@cpars$M_ageArray
  
  if (!change_com_sel) {
    args$selectivity <- args$selectivity[-5]
    args$data$nsel_block <- args$data$nsel_block - 1
    args$data$sel_block[args$data$sel_block >= 5] <- args$data$sel_block[args$data$sel_block >= 5] - 1
    args$vul_par <- args$vul_par[, -5]
    args$map_vul_par <- args$map_vul_par[, -5]
    args$map_vul_par[!is.na(args$map_vul_par)] <- 1:sum(!is.na(args$map_vul_par))
  }
  
  res <- do.call(SRA_scope, args)
  
  if(all(!res@conv) && max(abs(res@mean_fit$SD$gradient.fixed)) < 0.1) {
    args$vul_par[, 1:5] <- MSEtool:::ilogit(res@mean_fit$report$vul_par[, 1:5])
    args$s_vul_par <- MSEtool:::ilogit(res@mean_fit$report$s_vul_par)
    
    args$map_vul_par[, 1:5][abs(res@mean_fit$report$vul_par[, 1:5]) >= 8] <- NA
    args$map_s_vul_par[abs(res@mean_fit$report$s_vul_par) >= 8] <- NA
    
    res <- do.call(SRA_scope, args)
  }
  ret <- retrospective(res, 7)
  
  return(c(res, ret))
}


check_conv <- function(x) vapply(x, function(xx) unique(xx[[1]]@conv), logical(1))
calc_rho <- function(x) vapply(x, function(xx) ifelse(unique(xx[[1]]@conv), summary(xx[[2]])[3, 1], NA), numeric(1))


############### Generate pollock OM SS
############### Profile index selectivity at maxage
args <- get_pollock_base(h = 0.81)
setup(6)
sfExport(list = c("args", "eval_retro"))

sel_maxage <- seq(0.05, 0.95, 0.05)
sel_profile <- sfClusterApplyLB(sel_maxage, eval_retro, type = "survey_sel", args = args)
check_conv(sel_profile)
rho <- vapply(sel_profile, function(xx) summary(xx[[2]])[3, 1], numeric(1))

png("report/pollock/retro_survey_sel.png", width = 4, height = 4, units = "in", res = 400)
plot(sel_maxage, rho, typ = "o", xlab = "Survey selectivity at max age", ylab = "SSB Mohn's rho")
dev.off()

plot(sel_profile[[1]][[1]], retro = sel_profile[[1]][[2]], file = "report/pollock/NR1", 
     title = "Pollock Base, ComSel 0.5, Survey MaxAgeSel 0.05", dir = getwd(), 
     s_name = c("NEFSCspring", "NEFSCfall"), f_name = c("Commercial", "Recreational"), open_file = FALSE)

# Generate the OM
args <- get_pollock_base(100, h = 0.81)
NR1 <- lapply(0.05, eval_retro, type = "survey_sel", args = args)
saveRDS(NR1[[1]][[1]], file = "pollock/SRA_SS.rds")



############## Generate pollock OM SWB
############## Profile by lambda Index
args <- get_pollock_flatsel(h = 0.81)
setup(6)
sfExport(list = c("args", "eval_retro"))
Isd_mult <- seq(0.1, 0.5, 0.05) 
Isd_mult <- seq(1, 15, 1)
Isd_profile <- sfClusterApplyLB(Isd_mult, eval_retro, type = "I_sd", args = args)
#Isd_profile2 <- sfClusterApplyLB(Isd_mult, eval_retro, type = "I_sd", args = args, est_s_sel = TRUE)
conv <- vapply(Isd_profile, function(xx) unique(xx[[1]]@conv), logical(1))
rho <- vapply(Isd_profile, function(xx) summary(xx[[2]])[3, 1], numeric(1))

png("report/pollock/retro_index_lambda.png", width = 4, height = 4, units = "in", res = 400)
plot(Isd_mult, rho, typ = "o", xlab = "Index lambda", ylab = "SSB Mohn's rho")
dev.off()

plot(Isd_profile[[8]][[1]], retro = Isd_profile[[8]][[2]], file = "report/pollock/NR2", 
     title = "Pollock Base, ComSel 0.5, Survey Lambda 8x", dir = getwd(), 
     s_name = c("NEFSCspring", "NEFSCfall"), f_name = c("Commercial", "Recreational"), open_file = FALSE)

# Generate the OM
args <- get_pollock_base(100, h = 0.81)
NR2 <- lapply(8, eval_retro, type = "I_sd", args = args)
saveRDS(NR2[[1]][[1]], file = "pollock/SRA_SWB.rds")




############## Generate pollock OM SWF
############## Profile by lambda Index
args <- get_pollock_flatsel(h = 0.81)
setup(6)
sfExport(list = c("args", "eval_retro"))
Isd_mult <- seq(0.1, 0.5, 0.05) 
Isd_mult <- seq(1, 15, 1)
Isd_profile <- sfClusterApplyLB(Isd_mult, eval_retro, type = "I_sd", args = args)
#Isd_profile2 <- sfClusterApplyLB(Isd_mult, eval_retro, type = "I_sd", args = args, est_s_sel = TRUE)
conv <- vapply(Isd_profile, function(xx) unique(xx[[1]]@conv), logical(1))
rho <- vapply(Isd_profile, function(xx) summary(xx[[2]])[3, 1], numeric(1))

png("report/pollock/retro_index_lambda.png", width = 4, height = 4, units = "in", res = 400)
plot(Isd_mult, rho, typ = "o", xlab = "Index lambda", ylab = "SSB Mohn's rho")
dev.off()

plot(Isd_profile[[10]][[1]], retro = Isd_profile[[10]][[2]], file = "report/pollock/NR3", 
     title = "Pollock FlatSel, ComSel 0.5, Survey Lambda 10x", dir = getwd(), 
     s_name = c("NEFSCspring", "NEFSCfall"), f_name = c("Commercial", "Recreational"), open_file = FALSE)

# Generate the OM
args <- get_pollock_flatsel(100, h = 0.81)
NR3 <- lapply(10, eval_retro, type = "I_sd", args = args)
saveRDS(NR3[[1]][[1]], file = "pollock/SRA_NR3.rds")


