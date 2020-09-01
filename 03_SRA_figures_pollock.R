

################# Compare
res1 <- readRDS("pollock/SRA_pollock_base.rds")
res2 <- readRDS("pollock/SRA_pollock_flatsel.rds")
res3 <- readRDS("pollock/SRA_NR1.rds")
res4 <- readRDS("pollock/SRA_NR2.rds")
res5 <- readRDS("pollock/SRA_NR3.rds")


compare_SRA(res1, res2, res3, res4, res5,
            filename = "report/pollock/compare_pollock_ret", dir = getwd(),
            title = "ASAP3 pollock comparison",
            s_name = c("NEFSCspring", "NEFSCfall"), f_name = c("Commercial", "Recreational"),
            scenario = list(names = c("Base", "FlatSel", 
                                      "Base_LowSurvSel", "Base_UpSurv", "FlatSel_UpSurv")))



################# Compare SSB and R

SSB <- cbind %>% do.call(lapply(list(res1, res2, res3, res4, res5), function(x) x@SSB[1, ]))
RR <- cbind %>% do.call(lapply(list(res1, res2, res3, res4, res5), function(x) x@Misc[[1]]$R))

png("report/pollock/SSB_R.png", height = 6, width = 5, units = "in", res = 400)
par(mar = c(2, 4, 1, 1), oma = c(2, 0, 0, 0), mfrow = c(2, 1))

matplot(1970:2019, SSB, xlab = "Year", type = "l", col = c(1, 3, 4, 5, 2), lty = c(3, 3, 1, 1, 1), lwd = 2, ylim = c(0, 1e6))
abline(h = 0, col = "grey")
legend("top", c("Base", "FlatSel", "NR1", "NR2", "NR3"), col = c(1, 3, 4, 5, 2), lty = c(3, 3, 1, 1, 1), lwd = 2)

matplot(1970:2019, RR, ylab = "Recruitment", xlab = "Year", type = "l", lty = c(3, 3, 1, 1, 1), lwd = 2, 
        col = c(1, 3, 4, 5, 2), ylim = c(0, 1.1 * max(RR)))
abline(h = 0, col = "grey")

mtext("Year", outer = TRUE, side = 1, line = 1)
dev.off()


################# OM rho
setup(3)
rr <- sfLapply(list(res3, res4, res5), retrospective, nyr = 7, figure = FALSE)

png("report/pollock/pollock_OM_rho.png", height = 4, width = 6.5, units = "in", res = 400)
par(mfcol = c(2, 3), mar = c(4, 4, 1, 1))
for(i in 1:3) {
  matplot(1970:2019, rr[[i]]@TS[, , 1] %>% t(), type = 'l', col = gplots::rich.colors(7), lty = 1, ylim = c(0, 1), 
          xlab = "Year", ylab = "Commerical F")
  abline(h = 0, col = "grey")
  legend("topleft", c(paste0("NR", i), latex2exp::TeX(paste0("$\\rho = ", summary(rr[[i]])[1,1], "$"))), bty = "n")
  matplot(1970:2019, rr[[i]]@TS[, , 3] %>% t(), type = 'l', col = gplots::rich.colors(7), lty = 1, ylim = c(0, 6e5), 
          xlab = "Year", ylab = "SSB")
  abline(h = 0, col = "grey")
  legend("topleft", c(paste0("NR", i), latex2exp::TeX(paste0("$\\rho = ", summary(rr[[i]])[2,1], "$"))), bty = "n")
}
dev.off()


setup(3)
rr <- sfLapply(list(res1, res2), retrospective, nyr = 7, figure = FALSE)
EM <- c("Base", "FlatSel")

png("report/pollock/pollock_EM_rho.png", height = 4, width = 2 * 6.5 / 3, units = "in", res = 400)
par(mfcol = c(2, 2), mar = c(4, 4, 1, 1))
for(i in 1:2) {
  matplot(1970:2019, rr[[i]]@TS[, , 1] %>% t(), type = 'l', col = gplots::rich.colors(7), lty = 1, ylim = c(0, 1), 
          xlab = "Year", ylab = "Commerical F")
  abline(h = 0, col = "grey")
  legend("topleft", c(EM[i], latex2exp::TeX(paste0("$\\rho = ", summary(rr[[i]])[1,1], "$"))), bty = "n")
  matplot(1970:2019, rr[[i]]@TS[, , 3] %>% t(), type = 'l', col = gplots::rich.colors(7), lty = 1, ylim = c(0, 6e5), 
          xlab = "Year", ylab = "SSB")
  abline(h = 0, col = "grey")
  legend("topleft", c(EM[i], latex2exp::TeX(paste0("$\\rho = ", summary(rr[[i]])[2,1], "$"))), bty = "n")
}
dev.off()


