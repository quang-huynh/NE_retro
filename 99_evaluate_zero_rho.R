
## View why there are no retros in NR1
library(MSEtool)
library(mfSCA)

## Get MSE object
MSE <- readRDS("GoM_cod/MSE_cod_NR1.rds")
MSE@MPs

# Perfect implementation = MP #8
Data <- MSE@Misc$Data[[8]]

# M02 = MP #1
#Data <- MSE@Misc$Data[[1]]

# Get retro matrix
rho_fn <- function(x, y) {
  if(is.matrix(x$rho)) {
    rho <- x$rho[y, ]
  } else {
    rho <- x$rho
  }
  structure(rho, names = sapply(x$diagnostic, getElement, "Year"))
}
rho <- sapply(Data@Misc, rho_fn, y = 1)

apply(rho, 1, median)
matplot(rho, typ = 'o')


rho[15, ]

# Fit M02 in year 79, sim 1
SRA_M02 <- readRDS("GoM_cod/SRA_cod_M02.rds")
SRA_MRAMP <- readRDS("GoM_cod/SRA_cod_MRAMP.rds")
Cbias <- c(2.25, 1.25, 1)

SRA_OM <- readRDS(paste0("GoM_cod/SRA_NR1.rds")) %>% 
  generate_OM_args(AddInd_from_residuals = FALSE, AC_from_residuals = TRUE, Cobs = 0.05,
                   Cbias = 2.25)
args_M02 <- get_cod_M02(nsim = SRA_OM@OM@nsim) %>% 
  generate_EM_args(SRA_EM = SRA_M02, Cbias = 2.25)
M02 <- create_MP(args_M02)
debug(M02)

#M02(1, Data = Data, reps = 1)
#do_Assessment <- mf_SCA(x, Data, args = args, ymax = 79)
#saveRDS(do_Assessment, file = 'cod_no_rho_sim1.rds')


M02(36, Data = Data, reps = 1)
do_Assessment <- mf_SCA(x, Data, args = args, ymax = 79)
saveRDS(do_Assessment, file = 'cod_no_rho_sim36.rds')

do_assess <- readRDS("cod_no_rho_sim36.rds")
ret <- mfSCA:::retrospective_mf_SCA(do_assess, 7)
#plot(ret)

png("cod_rho_0.01_sim36.png", width = 4, height = 3, res = 400, units = 'in')
par(mar = c(5, 4, 1, 1))
matplot(Data@Year[1:80] + 2018-37, ret@TS[,,2] %>% t(), typ = 'l', lty = 1, xlab = "Year", ylab = "SSB")
#matplot(Data@Year[1:80] + 2018-37, ret@TS[,,2] %>% t(), typ = 'l', lty = 1, xlab = "Year", ylab = "SSB",
#        ylim = c(40000, 100000))
lines(Data@Year[1:80] + 2018-37, ret@TS[1,,2], lwd = 2)
#legend("topleft", as.character(0:7), lwd = 1, col = 1:5, bty = "n")
abline(v = 2018, lty = 2)
abline(h = 0, col = "grey")
dev.off()


do_assess <- readRDS("cod_no_rho_sim1.rds")
ret <- mfSCA:::retrospective_mf_SCA(do_assess, 7)
#plot(ret)

png("cod_rho_0.06_sim1.png", width = 4, height = 3, res = 400, units = 'in')
par(mar = c(5, 4, 1, 1))
matplot(Data@Year[1:80] + 2018-37, ret@TS[,,2] %>% t(), typ = 'l', lty = 1, xlab = "Year", ylab = "SSB")
#matplot(Data@Year[1:80] + 2018-37, ret@TS[,,2] %>% t(), typ = 'l', lty = 1, xlab = "Year", ylab = "SSB",
#        ylim = c(40000, 100000))
lines(Data@Year[1:80] + 2018-37, ret@TS[1,,2], lwd = 2)
#legend("topleft", as.character(0:7), lwd = 1, col = 1:5, bty = "n")
abline(v = 2018, lty = 2)
abline(h = 0, col = "grey")
dev.off()








############# Example of negative rho in cod
library(MSEtool)
library(mfSCA)

## Get MSE object
MSE <- readRDS("GoM_cod/MSE_cod_NR1.rds")
MSE@MPs

# MRAMP MP = 3
Data <- MSE@Misc$Data[[3]]

# Get retro matrix
rho_fn <- function(x, y) {
  if(is.matrix(x$rho)) {
    rho <- x$rho[y, ]
  } else {
    rho <- x$rho
  }
  structure(rho, names = sapply(x$diagnostic, getElement, "Year"))
}
rho <- sapply(Data@Misc, rho_fn, y = 1)

apply(rho, 1, median)
matplot(rho, typ = 'o')
abline(h = 0, lty = 3)

rho[15, ]

# Fit M02 in year 79, sim 1
SRA_M02 <- readRDS("GoM_cod/SRA_cod_M02.rds")
SRA_MRAMP <- readRDS("GoM_cod/SRA_cod_MRAMP.rds")
Cbias <- c(2.25, 1.25, 1)

SRA_OM <- readRDS(paste0("GoM_cod/SRA_NR1.rds")) %>% 
  generate_OM_args(AddInd_from_residuals = FALSE, AC_from_residuals = TRUE, Cobs = 0.05,
                   Cbias = 2.25)
args_MRAMP <- get_cod_MRAMP(nsim = SRA_OM@OM@nsim) %>% 
  generate_EM_args(SRA_EM = SRA_MRAMP, Cbias = 2.25)
MRAMP <- create_MP(args_MRAMP)
debug(MRAMP)

#M02(1, Data = Data, reps = 1)
#do_Assessment <- mf_SCA(x, Data, args = args, ymax = 79)
#saveRDS(do_Assessment, file = 'cod_no_rho_sim1.rds')


MRAMP(1, Data = Data, reps = 1)
do_Assessment <- mf_SCA(x, Data, args = args, ymax = 79)
saveRDS(do_Assessment, file = 'cod_negative_rho_sim1.rds')

do_assess <- readRDS("cod_negative_rho_sim1.rds")
ret <- mfSCA:::retrospective_mf_SCA(do_assess, 7)
ret
plot(ret)

# Simulated Index
plot(Data@AddInd[1, 1, ], typ = 'o')

png("cod_negative_rho_sim1.png", width = 4, height = 3, res = 400, units = 'in')
par(mar = c(5, 4, 1, 1))
matplot(Data@Year[1:80] + 2018-37, ret@TS[,,2] %>% t(), typ = 'l', lty = 1, xlab = "Year", ylab = "SSB")
#matplot(Data@Year[1:80] + 2018-37, ret@TS[,,2] %>% t(), typ = 'l', lty = 1, xlab = "Year", ylab = "SSB",
#        ylim = c(40000, 100000))
lines(Data@Year[1:80] + 2018-37, ret@TS[1,,2], lwd = 2)
#legend("topleft", as.character(0:7), lwd = 1, col = 1:5, bty = "n")
abline(v = 2018, lty = 2)
abline(h = 0, col = "grey")
dev.off()


MRAMP(10, Data = Data, reps = 1)
do_Assessment <- mf_SCA(x, Data, args = args, ymax = 67)
saveRDS(do_Assessment, file = 'cod_negative_rho_sim1_early.rds')

do_Assessment <- readRDS("cod_negative_rho_sim1_early.rds")
ret <- mfSCA:::retrospective_mf_SCA(do_Assessment, 7)
ret
plot(ret)

png("cod_negative_rho_sim1_early.png", width = 4, height = 3, res = 400, units = 'in')
par(mar = c(5, 4, 1, 1))
matplot(Data@Year[1:68] + 2018-37, ret@TS[,,2] %>% t(), typ = 'l', lty = 1, xlab = "Year", ylab = "SSB")
#matplot(Data@Year[1:80] + 2018-37, ret@TS[,,2] %>% t(), typ = 'l', lty = 1, xlab = "Year", ylab = "SSB",
#        ylim = c(40000, 100000))
lines(Data@Year[1:68] + 2018-37, ret@TS[1,,2], lwd = 2)
#legend("topleft", as.character(0:7), lwd = 1, col = 1:5, bty = "n")
abline(v = 2018, lty = 2)
abline(h = 0, col = "grey")
dev.off()