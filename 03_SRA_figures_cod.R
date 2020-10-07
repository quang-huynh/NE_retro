


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

SSB <- cbind %>% do.call(lapply(list(res3, res5, res4, res1, res2), function(x) x@SSB[1, ]))
RR <- cbind %>% do.call(lapply(list(res3, res5, res4, res1, res2), function(x) x@Misc[[1]]$R))

png("report/GoM_cod/SSB_R.png", height = 5, width = 4, units = "in", res = 400)
par(mar = c(2, 4, 1, 1), oma = c(2, 0, 0, 0), mfrow = c(2, 1))

matplot(1982:2019, RR, ylab = "Recruitment", xlab = "Year", xlim = c(1980, 2020), type = "l", lty = c(1, 1, 1, 2, 2), 
        col = c(4, 2, 5, 1, 3), lwd = 2, ylim = c(0, 50000))
abline(h = 0, col = "grey")
legend("topright", c("MC", "IM", "MCIM", "M02", "MRAMP"), col = c(4, 2, 5, 1, 3), lwd = 2, lty = c(1, 1, 1, 2, 2), ncol = 2, bty = "n")

matplot(1982:2019, SSB, xlab = "Year", xlim = c(1980, 2020), type = "l", lty = c(1, 1, 1, 2, 2), 
        col = c(4, 2, 5, 1, 3), lwd = 2, ylim = c(0, 60000))
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
          xlim = c(1980, 2020), xlab = "Year", ylab = "F")
  abline(h = 0, col = "grey")
  legend("topleft", c(paste0("NR", i), latex2exp::TeX(paste0("$\\rho = ", summary(rr[[i]])[1,1], "$"))), bty = "n")
  matplot(1982:2019, rr[[i]]@TS[, , 2] %>% t(), type = 'l', col = gplots::rich.colors(7), lty = 1, ylim = c(0, 6e4), 
          xlim = c(1980, 2020), xlab = "Year", ylab = "SSB")
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
          xlim = c(1980, 2020), xlab = "Year", ylab = "F")
  abline(h = 0, col = "grey")
  legend("topleft", c(EM[i], latex2exp::TeX(paste0("$\\rho = ", summary(rr[[i]])[1,1], "$"))), bty = "n")
  matplot(1982:2019, rr[[i]]@TS[, , 2] %>% t(), type = 'l', col = gplots::rich.colors(7), lty = 1, ylim = c(0, 6e4), 
          xlim = c(1980, 2020), xlab = "Year", ylab = "SSB")
  abline(h = 0, col = "grey")
  legend("topleft", c(EM[i], latex2exp::TeX(paste0("$\\rho = ", summary(rr[[i]])[2,1], "$"))), bty = "n")
}
dev.off()





# Model peels
setup(3)
rr_OM <- sfLapply(list(res3, res5, res4), retrospective, nyr = 7, figure = FALSE)
rr_EM <- sfLapply(list(res1, res2), retrospective, nyr = 7, figure = FALSE)

model_name <- c("MC", "IM", "MCIM", "M02", "MRAMP")
rr_out <- lapply(1:5, function(i, x) x[[i]]@TS[, , 2] %>% reshape2::melt(value.name = "SSB") %>% mutate(Model = model_name[i]), 
                 x = c(rr_OM, rr_EM))
rr_out <- do.call(rbind, rr_out) %>% mutate(Model = factor(Model, levels = model_name))

rho_label <- data.frame(Model = model_name %>% factor(),
                        rho = sapply(c(rr_OM, rr_EM), function(x) summary(x)[2, 1] %>% format(trim = TRUE) %>% paste0("rho==", .)))

no_panel_gap <- theme(panel.spacing = unit(0, "in"))
legend_bottom <- theme(legend.position = "bottom")

ggplot(rr_out %>% filter(Peel > 0) %>% mutate(Peel = factor(Peel)), aes(Year, SSB)) + 
  geom_line(aes(group = Peel, colour = Peel)) + facet_wrap(~ Model) + 
  geom_text(data = rho_label, x = 1982, y = 0, vjust = "inward", hjust = "inward", aes(label = rho), parse = TRUE) +
  geom_line(data = filter(rr_out, Peel == 0), size = 0.5, colour = "black") + 
  coord_cartesian(ylim = c(0, 6e4)) + #scale_x_continuous(breaks = c(1970, 1990, 2010)) + 
  gfplot::theme_pbs() + no_panel_gap + legend_bottom
ggsave("report/GoM_cod/cod_cond_rho.png", width = 6.5, height = 5)



ggplot(rr_out %>% filter(Peel > 0) %>% mutate(Peel = factor(Peel)), aes(Year, SSB)) + 
  geom_line(aes(group = Peel, colour = Peel)) + facet_wrap(~ Model) + 
  geom_text(data = rho_label, x = 1998, y = 0, vjust = "inward", hjust = "inward", aes(label = rho), parse = TRUE) +
  geom_line(data = filter(rr_out, Peel == 0), size = 0.5, colour = "black") + 
  coord_cartesian(ylim = c(0, 4e4), xlim = c(1998, 2020)) + #scale_x_continuous(breaks = c(1970, 1990, 2010)) + 
  gfplot::theme_pbs() + no_panel_gap + legend_bottom
ggsave("report/GoM_cod/cod_cond_rho2.png", width = 6.5, height = 5)









