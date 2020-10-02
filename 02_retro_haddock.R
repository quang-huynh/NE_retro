
library(MSEtool)
library(dplyr)
library(mfSCA)

eval_retro <- function(x, y, type = c("M_inc", "index_mult", "index_omit", "index_error", "M_age", 
                                      "catch_mult", "M_inc2", "sigmaR"), years_I, years = NULL, args) {
  require(dplyr)
  type <- match.arg(type, several.ok = TRUE)
  
  if ("M_age" %in% type) {
    args$OM@cpars$M_ageArray[, y, years] <- x
  }
  
  if ("sigmaR" %in% type) {
    args$OM@Perr <- rep(x, 2)
  }
  
  if ("catch_mult" %in% type) {
    args$data$Chist[years, 1] <- x * args$data$Chist[years, 1]
  }
  
  if ("index_error" %in% type) {
    args$data$I_sd[years_I, ] <- x * args$data$I_sd[years_I, ]
  }
  
  if ("index_mult" %in% type) {
    if (missing(y)) y <- x
    if (missing(years_I)) years_I <- years
    args$data$Index[years_I, ] <- y * args$data$Index[years_I, ]
  }
  
  if ("index_omit" %in% type) {
    args$data$Index[(years - x + 1):years, ] <- NA
  }
  
  if ("M_inc" %in% type) { # M2 continues ramp (x = nyears)
    
    dM <- dim(args$OM@cpars$M_ageArray)
    
    nyears <- args$OM@nyears
    
    new_M <- round(0.2 + c(1:length(years)) * x, 2)
    new_M <- c(new_M, rep(new_M[length(new_M)], dM[3] - max(years)))
    
    new_Marr <- array(new_M, c(length(new_M), dM[1], dM[2])) %>% aperm(c(2, 3, 1))
    args$OM@cpars$M_ageArray[, , min(years):dM[3]] <- new_Marr

    # plot(args$OM@cpars$M_ageArray[1,4, ])
  }
  
  if ("M_inc2" %in% type) { # M2 continues ramp (x = nyears)
    dM <- dim(args$OM@cpars$M_ageArray)
    
    nyears <- args$OM@nyears
    
    new_M <- round(0.2 - c(1:length(years)) * x, 2) %>% rev()
    new_Marr <- array(new_M, c(length(years), args$OM@maxage, args$OM@nsim)) %>% aperm(c(3, 2, 1))
    args$OM@cpars$M_ageArray[, , years] <- new_Marr
    args$OM@cpars$M_ageArray[, , 1:(length(years) - 1)] <- max(new_M)
    
  }
  
  res <- do.call(SRA_scope, args)
  ret <- retrospective(res, 7, figure = FALSE)
  
  return(c(res, ret))
}

# Start here
args <- get_haddock(h = 0.74)
setup(6)
sfExport(list = "args")

#### Index omit - omit 3 years works!
Iomit <- 1:10
ret_cat_mult <- list()
ret_Iomit <- sfClusterApplyLB(Iomit, eval_retro, years = 42, type = "index_omit", args = args)  # Multiply for most recent ten years

rho <- vapply(ret_Iomit, function(xx) summary(xx[[2]])[2, 1], numeric(1))

plot(ret_Iomit[[3]][[1]], retro = ret_Iomit[[3]][[2]], file = "report/GoM_haddock/NR1", 
     title = "Haddock omit 3 recent yr surveys", dir = getwd(), 
     s_name = c("NEFSCspring", "NEFSCfall"), open_file = FALSE)

# Generate OM
args <- get_haddock(100, h = 0.74)
sfExport(list = "args")
NR1 <- lapply(3, eval_retro, years = 42, type = "index_omit", args = args)  # Multiply for most recent ten years
saveRDS(NR1[[1]][[1]], file = "GoM_haddock/SRA_NR1.rds")


#### Index error (better fit to index)
Isd_mult <- seq(0.05, 0.95, 0.05)
ret_Imult <- list()
ret_Imult[[1]] <- sfClusterApplyLB(Isd_mult, eval_retro, type = "index_error", years_I = 38:42, args = args) # 5
ret_Imult[[2]] <- sfClusterApplyLB(Isd_mult, eval_retro, type = "index_error", years_I = 33:42, args = args) # 10
ret_Imult[[3]] <- sfClusterApplyLB(Isd_mult, eval_retro, type = "index_error", years_I = 28:42, args = args) # 15

rho <- cbind %>% do.call(lapply(ret_Imult, function(x) vapply(x, function(xx) summary(xx[[2]])[2, 1], numeric(1))))
matplot(Isd_mult, rho)

plot(ret_Imult[[2]][[1]][[1]], retro = ret_Imult[[2]][[1]][[2]], file = "report/GoM_haddock/NR2", 
     title = "Haddock precise surveys 2009-2018", dir = getwd(), 
     s_name = c("NEFSCspring", "NEFSCfall"), open_file = FALSE)

# Generate OM
args <- get_haddock(100, h = 0.74)
NR2 <- lapply(0.05, eval_retro, years_I = 33:42, type = "index_error", args = args)
saveRDS(NR2[[1]][[1]], file = "GoM_haddock/SRA_NR2.rds")

# M prof
Mprof <- seq(0.3, 0.5, 0.01)
ret_Mprof <- list()
ret_Mprof[[1]] <- sfClusterApplyLB(Mprof, eval_retro, type = "M_age", y = 1:9, args = args, years = 1:28)
ret_Mprof[[2]] <- sfClusterApplyLB(Mprof, eval_retro, type = "M_age", y = 1:9, args = args, years = 1:33)
ret_Mprof[[3]] <- sfClusterApplyLB(Mprof, eval_retro, type = "M_age", y = 1:9, args = args, years = 1:38)

rho <- cbind %>% do.call(lapply(ret_Mprof, function(x) vapply(x, function(xx) summary(xx[[2]])[2, 1], numeric(1))))
matplot(Mprof, rho, typ = 'o', pch = 16, ylim = c(-0.5, 0.1))
abline(h = 0)

#plot(ret_Mprof[[1]][[11]][[1]], retro = ret_Mprof[[1]][[11]][[2]])
plot(ret_Mprof[[1]][[21]][[1]], retro = ret_Mprof[[1]][[21]][[2]], file = "report/GoM_haddock/NR3", 
     title = "Haddock M 0.5 -> 0.2 in 2003", dir = getwd(), 
     s_name = c("NEFSCspring", "NEFSCfall"), open_file = FALSE)


# Generate OM
args <- get_haddock(100, h = 0.74)
NR3 <- lapply(.50, eval_retro, years = 1:28, type = "M_age", y = 1:9, args = args)
saveRDS(NR3[[1]][[1]], file = "GoM_haddock/SRA_NR3.rds")



# I mult
Imult <- seq(0.1, 0.95, 0.05)
ret_Imult <- list()
ret_Imult[[1]] <- sfClusterApplyLB(Imult, eval_retro, type = "index_mult", years = 39:42, args = args) # Multiply for most recent 4 years
ret_Imult[[2]] <- sfClusterApplyLB(Imult, eval_retro, type = "index_mult", years = 38:42, args = args) # 5
ret_Imult[[3]] <- sfClusterApplyLB(Imult, eval_retro, type = "index_mult", years = 37:42, args = args) # 6

rho <- cbind %>% do.call(lapply(ret_Imult, function(x) vapply(x, function(xx) summary(xx[[2]])[2, 1], numeric(1))))

png("report/GoM_haddock/haddock_Imult_ret.png", height = 6, width = 5, units = 'in', res = 400)
par(mfrow = c(2, 1), mar = c(5, 4, 1, 1))
matplot(1977:2018, args$data$Index, col = "black", typ = "o", lty = 1:2, pch = c(16, 1), xlab = "Year", ylab = "Index")
matlines(2015:2018, 0.25 * args$data$Index[39:42, ], col = "red", typ = "o", lty = 1:2, pch = c(16, 1))
legend("topleft", c("NEFSCspring", "NEFSCfall"), lty = 1:2, pch = c(16, 1))

matplot(Imult, rho, typ = 'o', pch = 16, xlab = "Index multiplier", ylab = "SSB Mohn's rho", col = 2:4)
legend("topright", c('2015-2018', '2014-2018', '2013-2018'), col = 2:4, pch = 16)
abline(h = 0)

dev.off()

# Generate OM
args <- get_haddock(100, h = 0.74)
NR4 <- lapply(0.3, eval_retro, type = "index_mult", years = 39:42, args = args) 
saveRDS(NR4[[1]][[1]], file = "GoM_haddock/SRA_NR4.rds")



# Both Imult = 0.4 and profile for Cmult
c_mult <- seq(0.25, 0.75, 0.05)
ret_M_mult <- list()
ret_M_mult[[1]] <- sfClusterApplyLB(c_mult, eval_retro, type = c("catch_mult", "index_mult"), y = 0.4, years_I = 39:42, years = 26:42, args = args)  # Multiply for mo
ret_M_mult[[2]] <- sfClusterApplyLB(c_mult, eval_retro, type = c("catch_mult", "index_mult"), y = 0.4, years_I = 39:42, years = 31:42, args = args) # 5
ret_M_mult[[3]] <- sfClusterApplyLB(c_mult, eval_retro, type = c("catch_mult", "index_mult"), y = 0.4, years_I = 39:42, years = 36:42, args = args) # 15


png("report/GoM_haddock/haddock_ImultCmult_ret.png", height = 6, width = 5, units = 'in', res = 400)

par(mfrow = c(2, 1), mar = c(5, 4, 1, 1))
plot(1977:2018, ret_M_mult[[1]][[10]][[1]]@data$Chist[, 1], col = "red", typ = "o", pch = 16, xlab = "Year", ylab = "Catch")
lines(1977:2018, args$data$Chist, col = "black", typ = "o", pch = 16)
legend("topleft", c("Base", "Cmult = 0.7"), col = c("black", "red"), pch = 16)

rho <- cbind %>% do.call(lapply(ret_M_mult, function(x) vapply(x, function(xx) summary(xx[[2]])[2, 1], numeric(1))))
matplot(c_mult, rho, typ = 'o', pch = 16, xlab = "Catch multiplier", ylab = "SSB Mohn's rho", col = 2:4)
abline(h = 0)
legend("topright", c('2002-2018', '2007-2018', '2012-2018'), col = 2:4, pch = 16)
dev.off()

# Generate OM
args <- get_haddock(100, h = 0.74)
args$resample <- TRUE
NR5 <- lapply(0.7, eval_retro, type = c("catch_mult", "index_mult"), y = 0.4, years_I = 39:42, years = 26:42, args = args)
saveRDS(NR5[[1]][[1]], file = "GoM_haddock/SRA_NR5.rds")

















#Mprof <- seq(0.1, 0.4, 0.01)
#ret_Mprof <- list()
#ret_Mprof[[1]] <- sfClusterApplyLB(Mprof, eval_retro, type = "M_age", y = 1:9, args = args, years = 28:92)
#ret_Mprof[[2]] <- sfClusterApplyLB(Mprof, eval_retro, type = "M_age", y = 1:9, args = args, years = 33:92)
#ret_Mprof[[3]] <- sfClusterApplyLB(Mprof, eval_retro, type = "M_age", y = 1:9, args = args, years = 38:92)
#
#rho <- vapply(ret_Mprof, function(xx) summary(xx[[2]])[2, 1], numeric(1))
#plot(Mprof, rho)
#
#rho <- cbind %>% do.call(lapply(ret_Mprof, function(x) vapply(x, function(xx) summary(xx[[2]])[2, 1], numeric(1))))
#matplot(Mprof, rho, typ = 'o', pch = 16)
#abline(h = 0)

#### Index multiplier (gear catchability)
Imult <- seq(0.1, 0.95, 0.05)
ret_Imult <- list()
ret_Imult[[1]] <- sfClusterApplyLB(Imult, eval_retro, type = "index_mult", years = 39:42, args = args) # Multiply for most recent 4 years
ret_Imult[[2]] <- sfClusterApplyLB(Imult, eval_retro, type = "index_mult", years = 38:42, args = args) # 5
ret_Imult[[3]] <- sfClusterApplyLB(Imult, eval_retro, type = "index_mult", years = 37:42, args = args) # 6

rho <- cbind %>% do.call(lapply(ret_Imult, function(x) vapply(x, function(xx) summary(xx[[2]])[2, 1], numeric(1))))

saveRDS(ret_Imult[[1]][[4]], file = "GoM_haddock/haddock_no_ret.rds")

png("GoM_haddock/GoM_haddock_ret.png", height = 6, width = 5, units = 'in', res = 400)
par(mfrow = c(2, 1), mar = c(5, 4, 1, 1))
matplot(1977:2018, args$data$Index, col = "black", typ = "o", lty = 1:2, pch = c(16, 1), xlab = "Year", ylab = "Index")
matlines(2015:2018, 0.25 * args$data$Index[39:42, ], col = "red", typ = "o", lty = 1:2, pch = c(16, 1))
legend("topleft", c("NEFSCspring", "NEFSCfall"), lty = 1:2, pch = c(16, 1))

matplot(Imult, rho, typ = 'o', pch = 16, xlab = "Index multiplier", ylab = "SSB Mohn's rho", col = 2:4)
legend("topright", c('2015-2018', '2014-2018', '2013-2018'), col = 2:4, pch = 16)
abline(h = 0)

dev.off()

# M_inc
Minc <- seq(-0.01, -0.05, -0.01)
ret_M_mult <- list()
ret_M_mult[[1]] <- sfClusterApplyLB(Minc, eval_retro, type = c("M_inc", "index_mult"), y = 1, years = 33:42, args = args)  # Multiply for mo
ret_M_mult[[2]] <- sfClusterApplyLB(Minc, eval_retro, type = c("M_inc", "index_mult"), y = 1, years = 38:42, args = args) # 5
ret_M_mult[[3]] <- sfClusterApplyLB(Minc, eval_retro, type = c("M_inc", "index_mult"), y = 1, years = 28:42, args = args) # 15

# Min if decrease by 0.015 for 15 yrs
rho <- cbind %>% do.call(lapply(ret_M_mult, function(x) vapply(x, function(xx) summary(xx[[2]])[2, 1], numeric(1))))
matplot(Minc, rho, typ = 'o', pch = 16)
abline(h = 0)


Minc <- seq(-0.005, -0.01, -0.001)
ret_M_mult <- list()
ret_M_mult[[1]] <- sfClusterApplyLB(Minc, eval_retro, type = c("M_inc", "index_mult"), y = 0.4, years_I = 39:42, years = 31:35, args = args)  # Multiply for mo
ret_M_mult[[2]] <- sfClusterApplyLB(Minc, eval_retro, type = c("M_inc", "index_mult"), y = 0.4, years_I = 39:42, years = 36:40, args = args) # 5
ret_M_mult[[3]] <- sfClusterApplyLB(Minc, eval_retro, type = c("M_inc", "index_mult"), y = 0.4, years_I = 39:42, years = 38:42, args = args) # 15


png("GoM_haddock/GoM_haddock_ret2.png", height = 6, width = 5, units = 'in', res = 400)
par(mfrow = c(2, 1), mar = c(5, 4, 1, 1))
plot(1977:2018, rep(0.2, length(1977:2018)), col = "black", typ = "o", pch = 16, xlab = "Year", ylab = "M")
lines(1977:2018, ret_M_mult[[1]][[3]][[1]]@OM@cpars$M_ageArray[1, 1, 1:args$OM@nyears], col = "red", typ = "o", pch = 16)
legend("topleft", c("Base M", "M decrease"), col = c("black", "red"), pch = 16)

rho <- cbind %>% do.call(lapply(ret_M_mult, function(x) vapply(x, function(xx) summary(xx[[2]])[2, 1], numeric(1))))
matplot(Minc, rho, typ = 'o', pch = 16, xlab = "Annual decrease in M", ylab = "SSB Mohn's rho", col = 2:4)
abline(h = 0)
legend("topright", c('2005-2010', '2011-2016', '2014-2018'), col = 2:4, pch = 16)

dev.off()

plot(ret_M_mult[[1]][[3]][[1]], retro = ret_M_mult[[1]][[3]][[2]], compare = FALSE)


saveRDS(ret_M_mult[[1]][[3]], file = "GoM_haddock/haddock_no_ret2.rds")


# Index mult = 0.4 and catch mult
c_mult <- seq(0.25, 0.75, 0.05)
ret_M_mult <- list()
ret_M_mult[[1]] <- sfClusterApplyLB(c_mult, eval_retro, type = c("catch_mult", "index_mult"), y = 0.4, years_I = 39:42, years = 26:42, args = args)  # Multiply for mo
ret_M_mult[[2]] <- sfClusterApplyLB(c_mult, eval_retro, type = c("catch_mult", "index_mult"), y = 0.4, years_I = 39:42, years = 31:42, args = args) # 5
ret_M_mult[[3]] <- sfClusterApplyLB(c_mult, eval_retro, type = c("catch_mult", "index_mult"), y = 0.4, years_I = 39:42, years = 36:42, args = args) # 15


png("GoM_haddock/GoM_haddock_ret3.png", height = 6, width = 5, units = 'in', res = 400)

par(mfrow = c(2, 1), mar = c(5, 4, 1, 1))
plot(1977:2018, ret_M_mult[[1]][[10]][[1]]@data$Chist[, 1], col = "red", typ = "o", pch = 16, xlab = "Year", ylab = "Catch")
lines(1977:2018, args$data$Chist, col = "black", typ = "o", pch = 16)
legend("topleft", c("Base M", "M decrease"), col = c("black", "red"), pch = 16)

rho <- cbind %>% do.call(lapply(ret_M_mult, function(x) vapply(x, function(xx) summary(xx[[2]])[2, 1], numeric(1))))
matplot(c_mult, rho, typ = 'o', pch = 16, xlab = "Catch multiplier", ylab = "SSB Mohn's rho", col = 2:4)
abline(h = 0)
legend("topright", c('2002-2018', '2007-2018', '2012-2018'), col = 2:4, pch = 16)
dev.off()

saveRDS(ret_M_mult[[1]][[10]][[1]], file = "GoM_haddock/haddock_no_ret3.rds")

# Cmult only
c_mult <- seq(0.25, 0.75, 0.05)
ret_C_mult <- list()
ret_C_mult[[1]] <- sfClusterApplyLB(c_mult, eval_retro, type = "catch_mult", years = 26:42, args = args)  # Multiply for mo
ret_C_mult[[2]] <- sfClusterApplyLB(c_mult, eval_retro, type = "catch_mult", years = 31:42, args = args) # 5
ret_C_mult[[3]] <- sfClusterApplyLB(c_mult, eval_retro, type = "catch_mult", years = 36:42, args = args) # 15

#### SigmaR
sigmaR <- seq(0.85, 1.4, 0.05)

sR <- sfClusterApplyLB(sigmaR, eval_retro, type = "sigmaR", args = args)
rho <- vapply(sR, function(x) summary(x[[2]])[2, 1], numeric(1))
plot(sigmaR, rho)

#### M-age
M_age <- seq(0.1, 0.5, 0.05)
#M_age <- seq(0.25, 0.5, 0.05)

Mprof <- sfClusterApplyLB(M_age, eval_retro, type = "M_age", y = 3:9, years = 30:42, args = args)
Mprof <- sfClusterApplyLB(M_age, eval_retro, type = "M_age", y = 1:3, years = 30:42, args = args)

rho <- cbind %>% do.call(lapply(Mprof, function(x) summary(x[[2]])[2, 1]))
matplot(M_age, t(rho))
