

# Compare haddock
res1 <- readRDS(file = "GoM_haddock/SRA_haddock_sel9.rds")
res2 <- readRDS(file = "GoM_haddock/SRA_NR1.rds")
res3 <- readRDS(file = "GoM_haddock/SRA_NR2.rds")
res4 <- readRDS(file = "GoM_haddock/SRA_NR3.rds")


compare_SRA(res1, res2, res3, res4, dir = getwd(),
            filename = "report/GoM_haddock/compare_haddock", title = "ASAP3 haddock comparison",
            s_name = c("NEFSCspring", "NEFSCfall"), open_file = FALSE,
            scenario = list(names = c("Base", "Index omit", "Index precise", "M decrease")))


res5 <- readRDS(file = "GoM_haddock/SRA_NR4.rds")
res6 <- readRDS(file = "GoM_haddock/SRA_NR5.rds")

compare_SRA(res1, res2, res3, res4, res5, res6, dir = getwd(),
            filename = "report/GoM_haddock/compare_haddock2", title = "ASAP3 haddock comparison",
            s_name = c("NEFSCspring", "NEFSCfall"), open_file = FALSE,
            scenario = list(names = c("Base", "Index omit", "Index precise", "M decrease", "Index q", "Index q High catch")))



################# Compare SSB and R

SSB <- cbind %>% do.call(lapply(list(res1, res2, res3, res4), function(x) x@SSB[1, ]))
RR <- cbind %>% do.call(lapply(list(res1, res2, res3, res4), function(x) x@Misc[[1]]$R))

png("report/GoM_haddock/SSB_R.png", height = 6, width = 5, units = "in", res = 400)
par(mar = c(2, 4, 1, 1), oma = c(2, 0, 0, 0), mfrow = c(2, 1))

matplot(1977:2019, RR, ylab = "Recruitment", xlab = "Year", type = "l", lty = 1, lwd = 2, ylim = c(0, 1.1 * max(RR)))
abline(h = 0, col = "grey")
legend("topleft", c("Base", "NR1", "NR2", "NR3"), pch = 16, col = 1:4)

matplot(1977:2019, SSB, xlab = "Year", type = "l", lty = 1, lwd = 2, ylim = c(0, 1.1 * max(SSB)))
abline(h = 0, col = "grey")

mtext("Year", side = 1, outer = TRUE, line = 1)
dev.off()


setup(4)
rr <- sfLapply(list(res1, res2, res3, res4), retrospective, nyr = 7, figure = FALSE)

png("report/GoM_haddock/rho.png", height = 4, width = 6.5, units = "in", res = 400)
par(mfcol = c(2, 4), mar = c(4, 4, 1, 1))
mod_vec <- c("Base", "NR1", "NR2", "NR3")
for(i in 1:4) {
  matplot(1977:2019, rr[[i]]@TS[, , 1] %>% t(), type = 'l', col = gplots::rich.colors(7), lty = 1, ylim = c(0, 1.5), 
          xlab = "Year", ylab = "F")
  legend("topleft", c(mod_vec[i], latex2exp::TeX(paste0("$\\rho = ", summary(rr[[i]])[1,1], "$"))), bty = "n")
  matplot(1977:2019, rr[[i]]@TS[, , 2] %>% t(), type = 'l', col = gplots::rich.colors(7), lty = 1, ylim = c(0, 4e5), 
          xlab = "Year", ylab = "SSB")
  legend("topleft", c(mod_vec[i], latex2exp::TeX(paste0("$\\rho = ", summary(rr[[i]])[2,1], "$"))), bty = "n")
}
dev.off()




