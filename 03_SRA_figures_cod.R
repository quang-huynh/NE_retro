


################# Compare cod report
res1 <- readRDS("GoM_cod/SRA_cod_M02.rds")
res2 <- readRDS("GoM_cod/SRA_cod_MRAMP.rds")
res3 <- readRDS("GoM_cod/SRA_NR1.rds")
res4 <- readRDS("GoM_cod/SRA_NR2.rds")
res5 <- readRDS("GoM_cod/SRA_NR3.rds")


compare_SRA(res1, res2, res3, res4, res5,
            filename = "report/GoM_cod/compare_cod", dir = getwd(), title = "ASAP3 GM cod comparison",
            s_name = c("NEFSCspring", "NEFSCfall", "MAspring"), 
            scenario = list(names = c("M02", "MRAMP", "M02 Cmult", "MRAMP Cmult", "MRAMP Minc")),
            open_file = FALSE)


################# Compare SSB and R

SSB <- cbind %>% do.call(lapply(list(res1, res2, res3, res4, res5), function(x) x@SSB[1, ]))
RR <- cbind %>% do.call(lapply(list(res1, res2, res3, res4, res5), function(x) x@Misc[[1]]$R))

png("report/GoM_cod/SSB_R.png", height = 6, width = 5, units = "in", res = 400)
par(mar = c(2, 4, 1, 1), oma = c(2, 0, 0, 0), mfrow = c(2, 1))

matplot(1982:2019, RR, ylab = "Recruitment", xlab = "Year", type = "l", lty = 1, lwd = 2, ylim = c(0, 50000))
abline(h = 0, col = "grey")
legend("topright", c("M02", "MRAMP", "NR1", "NR2", "NR3"), pch = 16, col = 1:5)

matplot(1982:2019, SSB, xlab = "Year", type = "l", lty = 1, lwd = 2, ylim = c(0, 60000))
abline(h = 0, col = "grey")

mtext("Year", side = 1, outer = TRUE, line = 1)
dev.off()



################# OM rho
setup(3)
rr <- sfLapply(list(res3, res4, res5), retrospective, nyr = 7, figure = FALSE)

png("report/GoM_cod/cod_OM_rho.png", height = 4, width = 6.5, units = "in", res = 400)
par(mfcol = c(2, 3), mar = c(4, 4, 1, 1))
for(i in 1:3) {
  matplot(1982:2019, rr[[i]]@TS[, , 1] %>% t(), type = 'l', col = gplots::rich.colors(7), lty = 1, ylim = c(0, 3), 
          xlab = "Year", ylab = "F")
  abline(h = 0, col = "grey")
  legend("topleft", c(paste0("NR", i), latex2exp::TeX(paste0("$\\rho = ", summary(rr[[i]])[1,1], "$"))), bty = "n")
  matplot(1982:2019, rr[[i]]@TS[, , 2] %>% t(), type = 'l', col = gplots::rich.colors(7), lty = 1, ylim = c(0, 6e4), 
          xlab = "Year", ylab = "SSB")
  abline(h = 0, col = "grey")
  legend("topleft", c(paste0("NR", i), latex2exp::TeX(paste0("$\\rho = ", summary(rr[[i]])[2,1], "$"))), bty = "n")
}
dev.off()


################# EM rho
setup(3)
rr <- sfLapply(list(res1, res2), retrospective, nyr = 7, figure = FALSE)
EM <- c("M02", "MRAMP")

png("report/GoM_cod/cod_EM_rho.png", height = 4, width = 2 * 6.5 / 3, units = "in", res = 400)
par(mfcol = c(2, 2), mar = c(4, 4, 1, 1))
for(i in 1:2) {
  matplot(1982:2019, rr[[i]]@TS[, , 1] %>% t(), type = 'l', col = gplots::rich.colors(7), lty = 1, ylim = c(0, 3), 
          xlab = "Year", ylab = "F")
  abline(h = 0, col = "grey")
  legend("topleft", c(EM[i], latex2exp::TeX(paste0("$\\rho = ", summary(rr[[i]])[1,1], "$"))), bty = "n")
  matplot(1982:2019, rr[[i]]@TS[, , 2] %>% t(), type = 'l', col = gplots::rich.colors(7), lty = 1, ylim = c(0, 6e4), 
          xlab = "Year", ylab = "SSB")
  abline(h = 0, col = "grey")
  legend("topleft", c(EM[i], latex2exp::TeX(paste0("$\\rho = ", summary(rr[[i]])[2,1], "$"))), bty = "n")
}
dev.off()
