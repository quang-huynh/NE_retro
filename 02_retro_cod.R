
library(MSEtool)
library(dplyr)

eval_retro <- function(x, type = c("catch_mult", "M1", "M2", "oldMinc", "youngMinc", "maxRAMP", "highM"), maxM = NULL, 
                       age = NULL, years = NULL, args) {
  require(dplyr)
  type <- match.arg(type, several.ok = TRUE)
  
  if ("highM" %in% type) {
    args$OM@cpars$M_ageArray <- NULL
    args$OM@M <- rep(x, 2)
  }
  if ("oldMinc" %in% type) { # In MRAMP, increase only for certain ages
    args$OM@cpars$M_ageArray[, 1:(age-1), ] <- 0.2
  }
  if ("youngMinc" %in% type) { # In MRAMP, increase only for certain ages
    args$OM@cpars$M_ageArray[, -c(1:(age-1)), ] <- 0.2
  }
  if ("maxRAMP" %in% type) { # In MRAMP, increase only for certain ages
    args$OM@cpars$M_ageArray[args$OM@cpars$M_ageArray >= maxM] <- maxM
  }
  
  if ("catch_mult" %in% type) {
    args$data$Chist[years, 1] <- x * args$data$Chist[years, 1]
  }
  if ("M2" %in% type) { # M2 continues ramp (x = nyears)
    
    dM <- dim(args$OM@cpars$M_ageArray)
    
    slope <- (0.4 - 0.2)/(2003 - 1988)
    new_M <- round(0.4 + c(1:x) * slope, 2)
    
    new_M <- c(new_M, rep(max(new_M), dM[3] - 23 - x + 1))
    
    new_Marr <- array(new_M, c(length(new_M), dM[1], dM[2])) %>% aperm(c(2, 3, 1))
    
    args$OM@cpars$M_ageArray[, , 23:dM[3]] <- new_Marr
    # plot(1982:2018, args$OM@cpars$M_ageArray[1,1, 1:length(1982:2018)])
  }
  if ("M1" %in% type) { # M1 increases to plateau in 2003 (x = max M),
    dM <- dim(args$OM@cpars$M_ageArray)
    
    slope <- (x - 0.2)/(2003 - 1988)
    new_M <- round(0.2 + c(1:length(8:21)) * slope, 2)
    
    new_M <- c(new_M, rep(x, dM[3] - 23))
    
    new_Marr <- array(new_M, c(length(8:dM[3]), dM[1], dM[2])) %>% aperm(c(2, 3, 1))
    
    args$OM@cpars$M_ageArray[, , 8:dM[3]] <- new_Marr
    #plot(1982:2018, args$OM@cpars$M_ageArray[1,1, 1:length(1982:2018)])
  }
  
  res <- do.call(SRA_scope, args)
  ret <- retrospective(res, 7)
  
  return(c(res, ret))
}

############### Catch multiplier for M 0.2
args <- get_cod_M02(0.81)
setup(6)
sfExport(list = "args")

c_mult <- seq(1.25, 3, 0.25)
ret_cat_mult <- list()
ret_cat_mult[[1]] <- sfClusterApplyLB(c_mult, eval_retro, type = "catch_mult", years = 28:37, args = args)  # Multiply for most recent ten years
ret_cat_mult[[2]] <- sfClusterApplyLB(c_mult, eval_retro, type = "catch_mult", years = 33:37, args = args) # 5
ret_cat_mult[[3]] <- sfClusterApplyLB(c_mult, eval_retro, type = "catch_mult", years = 23:37, args = args) # 15

rho <- cbind %>% do.call(lapply(ret_cat_mult, function(x) vapply(x, function(xx) summary(xx[[2]])[2, 1], numeric(1))))

png("report/GoM_cod/cod_M02_ret.png", height = 6, width = 5, units = 'in', res = 400)
par(mfrow = c(2, 1), mar = c(5, 4, 1, 1))
plot(1982:2018, args$data$Chist[, 1], typ = "o", pch = 16, xlab = "Year", ylab = "Catch", ylim = c(0, 25000))
lines(2008:2018, c(1, rep(2.25, 10)) * args$data$Chist[27:37, 1], col = 2, typ = "o", pch = 16)

matplot(c_mult, rho, typ = 'o', pch = 16, xlab = "Catch multiplier", ylab = "SSB Mohn's rho", col = 2:4)
legend("topright", c('2009-2018', '2014-2018', '2004-2018'), col = 2:4, pch = 16)
abline(h = 0)
dev.off()

plot(ret_cat_mult[[1]][[5]][[1]], retro = ret_cat_mult[[1]][[5]][[2]], file = "report/GoM_cod/cod_NR1", 
     title = "GM cod M02, Cmult 2.25 2009-2018", 
     dir = getwd(), s_name = c("NEFSCspring", "NEFSCfall", "MAspring"), open_file = FALSE)

# Generate the OM
args <- get_cod_M02(100)
args$OM@h <- rep(0.81, 2)
NR1 <- lapply(2.25, eval_retro, type = "catch_mult", years = 28:37, args = args)
saveRDS(NR1[[1]][[1]], file = "GoM_cod/SRA_NR1.rds")


############### Catch multiplier for MRAMP
args <- get_cod_MRAMP(0.81)
setup(6)
sfExport(list = "args")

c_mult <- seq(1, 1.5, 0.05)
ret_cat_mult <- list()
ret_cat_mult[[1]] <- sfClusterApplyLB(c_mult, eval_retro, type = "catch_mult", years = 28:37, args = args) # Multiply for most recent ten years
ret_cat_mult[[2]] <- sfClusterApplyLB(c_mult, eval_retro, type = "catch_mult", years = 33:37, args = args) # 5
ret_cat_mult[[3]] <- sfClusterApplyLB(c_mult, eval_retro, type = "catch_mult", years = 23:37, args = args) # 15

rho <- cbind %>% do.call(lapply(ret_cat_mult, function(x) vapply(x, function(xx) summary(xx[[2]])[2, 1], numeric(1))))

png("report/GoM_cod/cod_MRAMP_ret.png", height = 6, width = 5, units = 'in', res = 400)
par(mfrow = c(2, 1), mar = c(5, 4, 1, 1))
plot(1982:2018, args$data$Chist[, 1], typ = "o", pch = 16, xlab = "Year", ylab = "Catch", ylim = c(0, 25000))
lines(2008:2018, c(1, rep(1.25, 10)) * args$data$Chist[27:37, 1], col = 2, typ = "o", pch = 16)

matplot(c_mult, rho, typ = 'o', pch = 16, xlab = "Catch multiplier", ylab = "SSB Mohn's rho", col = 2:4)
legend("topright", c('2009-2018', '2014-2018', '2004-2018'), col = 2:4, pch = 16)
abline(h = 0)
dev.off()

plot(ret_cat_mult[[1]][[6]][[1]], retro = ret_cat_mult[[1]][[6]][[2]], file = "report/GoM_cod/NR2", 
     title = "GM cod MRAMP, Cmult 2 2009-2018", 
     dir = getwd(), s_name = c("NEFSCspring", "NEFSCfall", "MAspring"), open_file = FALSE)

# Generate the OM
args <- get_cod_MRAMP(h = 0.81)
NR2 <- lapply(1.25, eval_retro, type = "catch_mult", years = 28:37, args = args)
saveRDS(NR2[[1]][[1]], file = "GoM_cod/SRA_NR2.rds")


############### MRAMP only for older ages - doesn't work
# args <- get_cod_MRAMP(h = 0.81)
# 
# eval2 <- function(x, args) eval_retro(1, type = "oldMinc", args = args, age = x) 
# setup(6)
# sfExport(list = c("args", "eval_retro"))
# age_inc <- 2:9
# ret_Mage <- sfClusterApplyLB(age_inc, eval2, args = args)
# 
# rho <- vapply(ret_Mage, function(xx) summary(xx[[2]])[2, 1], numeric(1))
# plot(age_inc, rho)

############### MRAMP only for younger ages
args <- get_cod_MRAMP(h = 0.81)
eval2 <- function(x, args) eval_retro(1, type = "youngMinc", args = args, age = x) 
setup(6)
sfExport(list = c("args", "eval_retro"))
age_inc <- 2:9
ret_Mage <- sfClusterApplyLB(age_inc, eval2, args = args)

rho <- vapply(ret_Mage, function(xx) summary(xx[[2]])[2, 1], numeric(1))
plot(age_inc, rho)

# It's promising but let's now try increase in juvenile M with catch underreporting
c_mult <- seq(1, 3, 0.25)
ret_cm <- sfClusterApplyLB(c_mult, eval_retro, type = c("catch_mult", "youngMinc"), 
                           years = 28:37, args = args, age = 4)

rho <- vapply(ret_cm, function(xx) summary(xx[[2]])[2, 1], numeric(1))
plot(c_mult, rho)
abline(h = 0)

plot(ret_cm[[3]][[1]], retro = ret_cm[[3]][[2]])


########### Catch multiplier with MRAMP to only 0.3
args <- get_cod_MRAMP(h = 0.81)

setup(6)
sfExport(list = "args")

c_mult <- seq(1.25, 3, 0.25)
ret_M3 <- sfClusterApplyLB(c_mult, eval_retro, type = c("maxRAMP", "catch_mult"), maxM = 0.3,
                           years = 28:37, args = args)

rho <- vapply(ret_M3, function(xx) summary(xx[[2]])[2, 1], numeric(1))
plot(c_mult, rho)
abline(h = 0)
plot(ret_M3[[5]][[1]], retro = ret_M3[[5]][[2]])


############# M multiplier for MRAMP
#max_Mramp <- seq(0.65, 0.75, 0.01)
#ret_M <- sfClusterApplyLB(max_Mramp, eval_retro, type = "M1", args = args)
#eval_retro(5, type = "M2", args = args)
#
#rho <- vapply(ret_M, function(xx) summary(xx[[2]])[2, 1], numeric(1))
#plot(1982:2018, args$OM@cpars$M_ageArray[1,1, 1:37], xlab = "Year", ylab = "Natural mortality", type = "o", ylim = c(0, 0.5))

############# M continue ramp
args <- get_cod_MRAMP(h = 0.81)

setup(6)
sfExport(list = "args")
max_Mramp <- 1:15
ret_M <- sfClusterApplyLB(max_Mramp, eval_retro, type = "M2", args = args)

rho <- vapply(ret_M, function(xx) summary(xx[[2]])[2, 1], numeric(1))

Mramp <- cbind %>% do.call(lapply(ret_M, function(x) x[[1]]@OM@cpars$M_ageArray[1, 1, 1:37]))

png("report/GoM_cod/GoM_cod_MRAMP_ret2.png", height = 6, width = 5, units = 'in', res = 400)
par(mfrow = c(2, 1), mar = c(5, 4, 1, 1))
matplot(1982:2018, Mramp[, c(5, 1)], col = c("red", "black"), lty = 1, 
        type = "o", pch = 16, xlab = "Year", ylab = "Natural mortality", ylim = c(0, 0.6))
legend("topleft", c("MRAMP", "max M = 0.45"), col = c("black", "red"), lty = 1, pch = 16)

plot(apply(Mramp, 2, max), rho, xlab = "maximum M", ylab = "SSB Mohn's rho", pch = 16, type = "o")
abline(h = 0, lty = 2)
dev.off()

plot(ret_M[[4]][[1]], retro = ret_M[[5]][[2]], file = "report/GoM_cod/NR3", 
     title = "GM cod MRAMP to 0.45", 
     dir = getwd(), s_name = c("NEFSCspring", "NEFSCfall", "MAspring"), open_file = FALSE)

##### Generate OM
args <- get_cod_MRAMP(nsim = 100, h = 0.81)
sfExport(list = "args")
NR3 <- sfClusterApplyLB(4, eval_retro, type = "M2", args = args)
saveRDS(NR3[[1]][[1]], file = "GoM_cod/SRA_NR3.rds")


