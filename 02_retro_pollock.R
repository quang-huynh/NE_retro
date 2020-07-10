
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



############### Profile index selectivity
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
saveRDS(NR1[[1]][[1]], file = "pollock/SRA_NR1.rds")




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
saveRDS(NR2[[1]][[1]], file = "pollock/SRA_NR2.rds")

##### Flat com sel Flatsel
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



############## Profile by lambda fishery CAA
CAA_mult <- seq(0.25, 5, 0.25)
CAA_profile <- sfClusterApplyLB(CAA_mult, eval_retro, type = "CAA_mult", args = args)
rho <- calc_rho(CAA_profile)
png("pollock/retro_CAA_lambda.png", width = 4, height = 4, units = "in", res = 400)
plot(CAA_mult, rho, typ = "o", xlab = "Fishery age comps lambda", ylab = "SSB Mohn's rho")
dev.off()


# Catch multiplier - commercial
#C_mult <- seq(0.25, 5, 0.25)
#C_mult_5 <- sfClusterApplyLB(C_mult, eval_retro, type = "catch_com", args = args, years = 45:49)
#C_mult_10 <- sfClusterApplyLB(C_mult, eval_retro, type = "catch_com", args = args, years = 39:49)
#C_mult_15 <- sfClusterApplyLB(C_mult, eval_retro, type = "catch_com", args = args, years = 35:49)
#
#com_base <- lapply(list(C_mult_5, C_mult_10, C_mult_15), calc_rho)
#plot(C_mult_15[[20]][[1]], retro = C_mult_15[[20]][[2]], compare = FALSE)



#C_mult_5 <- sfClusterApplyLB(C_mult, eval_retro, type = "catch_com", args = args, years = 1:5)
#C_mult_10 <- sfClusterApplyLB(C_mult, eval_retro, type = "catch_com", args = args, years = 1:10)
#C_mult_15 <- sfClusterApplyLB(C_mult, eval_retro, type = "catch_com", args = args, years = 1:15)
#
#com_base <- lapply(list(C_mult_5, C_mult_10, C_mult_15), calc_rho)


#C_mult_5 <- sfClusterApplyLB(C_mult, eval_retro, type = "catch_rec", args = args, years = 45:49)
#C_mult_10 <- sfClusterApplyLB(C_mult, eval_retro, type = "catch_rec", args = args, years = 39:49)
#C_mult_15 <- sfClusterApplyLB(C_mult, eval_retro, type = "catch_rec", args = args, years = 35:49)

# The right number of positive and negative peels to get mean rho = 0
#rec_base <- lapply(list(C_mult_5, C_mult_10, C_mult_15), calc_rho)
#plot(C_mult_10[[17]][[1]], retro = C_mult_10[[17]][[2]], compare = FALSE)


# I multiplier
I_mult <- seq(0.5, 5, 0.25)
I_mult_5 <- sfClusterApplyLB(I_mult, eval_retro, type = "I_mult", args = args, years = 45:49)
I_mult_10 <- sfClusterApplyLB(I_mult, eval_retro, type = "I_mult", args = args, years = 39:49)
I_mult_15 <- sfClusterApplyLB(I_mult, eval_retro, type = "I_mult", args = args, years = 35:49)

late_base <- lapply(list(I_mult_5, I_mult_10, I_mult_15), calc_rho)
plot(I_mult_5[[11]][[1]], retro = I_mult_5[[11]][[2]], file = "pollock_flatcomsel_noret", 
     title = "ASAP3, Imult 3 2014-2018", dir = "pollock", 
     s_name = c("NEFSCspring", "NEFSCfall"), f_name = c("Commercial", "Recreational"))
saveRDS(I_mult_5[[11]][[1]], file = "pollock/pollock_flatcomsel_noret.rds")

png("pollock/pollock_ret.png", height = 6, width = 5, units = 'in', res = 400)
matplot(I_mult, do.call(cbind, late_base), typ = "o", pch = 16, xlab = "Index multiplier", 
        ylab = "SSB Mohn's rho")
abline(h = 0, lty = 3)
legend("topright", c("2014-2018", "2009-2018", "2004-2018"), col = 1:3, lty = 1:3, pch = 16)
dev.off()


M_mult <- seq(0.25, 4, 0.25)
M_mult_5 <- sfClusterApplyLB(M_mult, eval_retro, type = "M_mult", args = args, years = 45:49)
M_mult_10 <- sfClusterApplyLB(M_mult, eval_retro, type = "M_mult", args = args, years = 39:49)

M_base <- lapply(list(M_mult_5, M_mult_10), calc_rho)
plot(M_mult_10[[6]][[1]], retro = M_mult_10[[6]][[2]], file = "pollock_flatcomsel_noret2", 
     title = "ASAP3 M = 0.3 2009-2018", dir = "pollock", 
     s_name = c("NEFSCspring", "NEFSCfall"), f_name = c("Commercial", "Recreational"))
saveRDS(M_mult_10[[6]][[1]], file = "pollock/pollock_flatcomsel_noret2.rds")

png("pollock/pollock_ret2.png", height = 6, width = 5, units = 'in', res = 400)
matplot(M_mult * 0.2, do.call(cbind, M_base), typ = "o", pch = 16, xlab = "M", 
        ylab = "SSB Mohn's rho")
abline(h = 0, lty = 3)
legend("topright", c("2014-2018", "2009-2018"), col = 1:3, lty = 1:3, pch = 16)
dev.off()




#### I mult
png("pollock/pollock_ret3.png", height = 6, width = 5, units = 'in', res = 400)
I_mult_5 <- sfClusterApplyLB(I_mult, eval_retro, type = "I_mult", args = args, years = 45:49)
I_mult_10 <- sfClusterApplyLB(I_mult, eval_retro, type = "I_mult", args = args, years = 39:49)
I_mult_15 <- sfClusterApplyLB(I_mult, eval_retro, type = "I_mult", args = args, years = 35:49)

late_base <- lapply(list(I_mult_5, I_mult_10, I_mult_15), calc_rho)
plot(I_mult_5[[15]][[1]], retro = I_mult_5[[15]][[2]], file = "pollock_flatcomsel_noret3", 
     title = "ASAP3 flatsel Imult 4 2014-2018", dir = "pollock", 
     s_name = c("NEFSCspring", "NEFSCfall"), f_name = c("Commercial", "Recreational"))
saveRDS(I_mult_5[[1]][[1]], file = "pollock/pollock_flatcomsel_noret3.rds")

png("pollock/pollock_ret3.png", height = 6, width = 5, units = 'in', res = 400)
matplot(I_mult, do.call(cbind, late_base), typ = "o", pch = 16, xlab = "Index multiplier", 
        ylab = "SSB Mohn's rho")
abline(h = 0, lty = 3)
legend("topright", c("2014-2018", "2009-2018", "2004-2018"), col = 1:3, lty = 1:3, pch = 16)
dev.off()


# M mult
M_mult <- seq(0.25, 4, 0.25)
M_mult_5 <- sfClusterApplyLB(M_mult, eval_retro, type = "M_mult", args = args, years = 45:49)
M_mult_10 <- sfClusterApplyLB(M_mult, eval_retro, type = "M_mult", args = args, years = 39:49)

M_base <- lapply(list(M_mult_5, M_mult_10), calc_rho)
png("pollock/pollock_ret4.png", height = 6, width = 5, units = 'in', res = 400)
matplot(M_mult * 0.2, do.call(cbind, M_base), typ = "o", pch = 16, xlab = "M", 
        ylab = "SSB Mohn's rho")
abline(h = 0, lty = 3)
legend("topright", c("2014-2018", "2009-2018"), col = 1:3, lty = 1:3, pch = 16)
dev.off()

M_noret <- eval_retro(0.33/0.2, type = "M_mult", args = args, years = 39:49)

plot(M_noret[[1]], retro = M_noret[[2]], file = "pollock_flatcomsel_noret4", 
     title = "ASAP3 M = 0.33 2009-2018", dir = "pollock", 
     s_name = c("NEFSCspring", "NEFSCfall"), f_name = c("Commercial", "Recreational"))
saveRDS(M_noret[[1]], file = "pollock/pollock_flatcomsel_noret4.rds")

# Base with combination of M mult and I mult
I_mult <- seq(0.5, 5, 0.25)
combo <- lapply(I_mult, function(x) c(x, 1.25) %>% structure(names = c("I_mult", "M_mult")))

combo_mult <- sfClusterApplyLB(combo, eval_retro, type = c("I_mult", "M_mult"), 
                               args = args, years = list(45:49, 39:49))

plot(I_mult, calc_rho(combo_mult))
abline(h = 0, lty = 3)

plot(combo_mult[[6]][[1]], retro = combo_mult[[6]][[2]], file = "pollock_flatcomsel_noret5", 
     title = "ASAP3, Imult 1.75 2014-2018, M = 0.25 2009-2018", dir = "pollock", 
     s_name = c("NEFSCspring", "NEFSCfall"), f_name = c("Commercial", "Recreational"))
saveRDS(combo_mult[[6]][[1]], file = "pollock/pollock_flatcomsel_noret5.rds")

png("pollock/pollock_ret5.png", height = 6, width = 5, units = 'in', res = 400)
plot(I_mult, calc_rho(combo_mult), typ = "o", pch = 16, xlab = "Index multiplier", 
        ylab = "SSB Mohn's rho")
abline(h = 0, lty = 3)
dev.off()



# Compare
res <- readRDS("pollock/pollock_base.rds")
res2 <- readRDS("pollock/pollock_flatsel.rds")
res3 <- readRDS("pollock/pollock_flatcomsel_noret.rds")
res4 <- readRDS("pollock/pollock_flatcomsel_noret2.rds")
res5 <- readRDS("pollock/pollock_flatcomsel_noret3.rds")
res6 <- readRDS("pollock/pollock_flatcomsel_noret4.rds")
res7 <- readRDS("pollock/pollock_flatcomsel_noret5.rds")


compare_SRA(res, res2, res3, res4, res5, res6, res7,
            filename = "compare_pollock_ret", dir = "pollock", title = "ASAP3 pollock comparison",
            s_name = c("NEFSCspring", "NEFSCfall"), f_name = c("Commercial", "Recreational"),
            scenario = list(names = c("Base", "FlatSel", 
                                      "Base_FlatComSel_Imult", "Base_FlatComSel_Mmult",
                                      "FlatSel_FlatComSel_Imult", "FlatSel_FlatComSel_Mmult",
                                      "Base_FlatComSel_Combo")))



SSB <- cbind %>% do.call(lapply(list(res, res2, res3, res4, res5), function(x) x@SSB[1, ]))

png("pollock/SSB.png", height = 4, width = 5, units = "in", res = 400)
par(mar = c(5, 4, 1, 1))
matplot(1970:2019, SSB, xlab = "Year", type = "l", lty = 1, lwd = 2, ylim = c(0, 6e5))
abline(h = 0, col = "grey")
legend("top", c("Base", "FlatSel", "NR1", "NR2", "NR3"), pch = 16, col = 1:5)
dev.off()

RR <- cbind %>% do.call(lapply(list(res, res2, res3, res4, res5), function(x) x@Misc[[1]]$R))

png("pollock/Recruits.png", height = 4, width = 5, units = "in", res = 400)
par(mar = c(5, 4, 1, 1))
matplot(1970:2019, RR, ylab = "Recruitment", xlab = "Year", type = "l", lty = 1, lwd = 2, ylim = c(0, 1.1e5))
abline(h = 0, col = "grey")
legend("topleft", c("Base", "FlatSel", "NR1", "NR2", "NR3"), pch = 16, col = 1:5)
dev.off()

setup(3)
rr <- sfLapply(list(res3, res4, res5), retrospective, nyr = 7, figure = FALSE)

png("pollock/rho.png", height = 4, width = 6.5, units = "in", res = 400)
par(mfcol = c(2, 3), mar = c(4, 4, 1, 1))
for(i in 1:3) {
  matplot(rr[[i]]@TS[, , 1] %>% t(), type = 'l', col = gplots::rich.colors(7), lty = 1, ylim = c(0, 1), 
          xlab = "Year", ylab = "Commerical F")
  legend("topleft", c(paste0("NR", i), latex2exp::TeX(paste0("$\\rho = ", summary(rr[[i]])[1,1], "$"))), bty = "n")
  matplot(rr[[i]]@TS[, , 3] %>% t(), type = 'l', col = gplots::rich.colors(7), lty = 1, ylim = c(0, 6e5), 
          xlab = "Year", ylab = "SSB")
  legend("topleft", c(paste0("NR", i), latex2exp::TeX(paste0("$\\rho = ", summary(rr[[i]])[2,1], "$"))), bty = "n")
}
dev.off()


setup(3)
rr <- sfLapply(list(res, res2), retrospective, nyr = 7, figure = FALSE)
EM <- c("Base", "FlatSel")

png("pollock/rho_EM.png", height = 4, width = 2 * 6.5 / 3, units = "in", res = 400)
par(mfcol = c(2, 2), mar = c(4, 4, 1, 1))
for(i in 1:2) {
  matplot(rr[[i]]@TS[, , 1] %>% t(), type = 'l', col = gplots::rich.colors(7), lty = 1, ylim = c(0, 1), 
          xlab = "Year", ylab = "Commerical F")
  legend("topleft", c(EM[i], latex2exp::TeX(paste0("$\\rho = ", summary(rr[[i]])[1,1], "$"))), bty = "n")
  matplot(rr[[i]]@TS[, , 3] %>% t(), type = 'l', col = gplots::rich.colors(7), lty = 1, ylim = c(0, 6e5), 
          xlab = "Year", ylab = "SSB")
  legend("topleft", c(EM[i], latex2exp::TeX(paste0("$\\rho = ", summary(rr[[i]])[2,1], "$"))), bty = "n")
}
dev.off()
