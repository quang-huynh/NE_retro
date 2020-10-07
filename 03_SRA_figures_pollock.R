

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

SSB <- cbind %>% do.call(lapply(list(res3, res4, res5, res1, res2), function(x) x@SSB[1, ]))
RR <- cbind %>% do.call(lapply(list(res3, res4, res5, res1, res2), function(x) x@Misc[[1]]$R))

png("report/pollock/SSB_R.png", height = 5, width = 4, units = "in", res = 400)
par(mar = c(2, 4, 1, 1), oma = c(2, 0, 0, 0), mfrow = c(2, 1))

matplot(1970:2019, SSB, xlab = "Year", type = "l", col = c(4, 5, 2, 1, 3), lty = c(1, 1, 1, 3, 3), lwd = 2, ylim = c(0, 1e6))
abline(h = 0, col = "grey")
legend("topright", c("SS", "SWB", "SWF", "Base", "FlatSel"), col = c(4, 5, 2, 1, 3), lty = c(1, 1, 1, 3, 3), lwd = 2, ncol = 2, bty = "n")

matplot(1970:2019, RR, ylab = "Recruitment", xlab = "Year", type = "l", lty = c(1, 1, 1, 3, 3), lwd = 2, 
        col = c(4, 5, 2, 1, 3), ylim = c(0, 1.1 * max(RR)))
abline(h = 0, col = "grey")

mtext("Year", outer = TRUE, side = 1, line = 1)
dev.off()

################# OM rho
setup(3)
rr_OM <- sfLapply(list(res3, res4, res5), retrospective, nyr = 7, figure = FALSE)


png("report/pollock/pollock_OM_rho.png", height = 4, width = 6.5, units = "in", res = 400)
par(mfcol = c(2, 3), mar = c(4, 4, 1, 1))
for(i in 1:3) {
  matplot(1970:2019, rr_OM[[i]]@TS[, , 1] %>% t(), type = 'l', col = gplots::rich.colors(7), lty = 1, ylim = c(0, 1), 
          xlab = "Year", ylab = "Commerical F")
  abline(h = 0, col = "grey")
  legend("topleft", c(paste0("NR", i), latex2exp::TeX(paste0("$\\rho = ", summary(rr_OM[[i]])[1,1], "$"))), bty = "n")
  matplot(1970:2019, rr_OM[[i]]@TS[, , 3] %>% t(), type = 'l', col = gplots::rich.colors(7), lty = 1, ylim = c(0, 6e5), 
          xlab = "Year", ylab = "SSB")
  abline(h = 0, col = "grey")
  legend("topleft", c(paste0("NR", i), latex2exp::TeX(paste0("$\\rho = ", summary(rr_OM[[i]])[2,1], "$"))), bty = "n")
}
dev.off()


setup(3)
rr_EM <- sfLapply(list(res1, res2), retrospective, nyr = 7, figure = FALSE)
EM <- c("Base", "FlatSel")

png("report/pollock/pollock_EM_rho.png", height = 4, width = 2 * 6.5 / 3, units = "in", res = 400)
par(mfcol = c(2, 2), mar = c(4, 4, 1, 1))
for(i in 1:2) {
  matplot(1970:2019, rr_EM[[i]]@TS[, , 1] %>% t(), type = 'l', col = gplots::rich.colors(7), lty = 1, ylim = c(0, 1), 
          xlab = "Year", ylab = "Commerical F")
  abline(h = 0, col = "grey")
  legend("topleft", c(EM[i], latex2exp::TeX(paste0("$\\rho = ", summary(rr_EM[[i]])[1,1], "$"))), bty = "n")
  matplot(1970:2019, rr_EM[[i]]@TS[, , 3] %>% t(), type = 'l', col = gplots::rich.colors(7), lty = 1, ylim = c(0, 6e5), 
          xlab = "Year", ylab = "SSB")
  abline(h = 0, col = "grey")
  legend("topleft", c(EM[i], latex2exp::TeX(paste0("$\\rho = ", summary(rr_EM[[i]])[2,1], "$"))), bty = "n")
}
dev.off()





no_panel_gap <- theme(panel.spacing = unit(0, "in"))
legend_bottom <- theme(legend.position = "bottom")


# SSB and R
Model <- data.frame(name = c("SS", "SWB", "SWF", "Base", "FlatSel"), type = c("OM", "OM", "OM", "AM", "AM"))
SSB <- lapply(list(res3, res4, res5, res1, res2), function(x) x@SSB[1, ]) %>% do.call(cbind, .) %>% 
  structure(dimnames = list(Year = 1970:2019, name = Model$name)) %>% reshape2::melt() %>% mutate(y = "SSB") %>%
  left_join(Model, by = "name")

RR <- lapply(list(res3, res4, res5, res1, res2), function(x) x@Misc[[1]]$R) %>% do.call(cbind, .) %>% 
  structure(dimnames = list(Year = 1970:2019, name = Model$name)) %>% reshape2::melt() %>% mutate(y = "Recruitment") %>%
  left_join(Model, by = "name")

ggplot(rbind(SSB, RR) %>% mutate(y = factor(y, levels = c("SSB", "Recruitment")), 
                                 type = factor(type, levels = c("OM", "AM")),
                                 name = factor(name, levels = Model$name)) %>% filter(value <= 1e6), 
       aes(x = Year, y = value, linetype = type, colour = name)) + geom_line() +
  geom_hline(yintercept = 0, colour = "white") +
  facet_grid(y~., scales = "free_y", switch = "y") + ylab(NULL) +
  scale_y_continuous(labels = scales::comma, n.breaks = 5) + guides(linetype = FALSE) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom + 
  theme(legend.title = element_blank(), strip.placement = "outside")
ggsave("report/pollock/model_compare.png", height = 4, width = 4)

# Model peels
setup(3)
rr_OM <- sfLapply(list(res3, res4, res5), retrospective, nyr = 7, figure = FALSE)
rr_EM <- sfLapply(list(res1, res2), retrospective, nyr = 7, figure = FALSE)

model_name <- c("SS", "SWB", "SWF", "Base", "FlatSel")
rr_out <- lapply(1:5, function(i, x) x[[i]]@TS[, , 3] %>% reshape2::melt(value.name = "SSB") %>% mutate(Model = model_name[i]), 
                 x = c(rr_OM, rr_EM))
rr_out <- do.call(rbind, rr_out) %>% mutate(Model = factor(Model, levels = model_name))

rho_label <- data.frame(Model = model_name %>% factor(),
                        rho = sapply(c(rr_OM, rr_EM), function(x) summary(x)[3, 1] %>% round(2) %>% paste0("rho==", .)))


ggplot(rr_out %>% filter(Peel > 0) %>% mutate(Peel = factor(Peel)), aes(Year, SSB)) + 
  geom_line(aes(group = Peel, colour = Peel)) + facet_wrap(~ Model) + 
  geom_text(data = rho_label, x = 1970, y = 0, vjust = "inward", hjust = "inward", aes(label = rho), parse = TRUE) +
  geom_line(data = filter(rr_out, Peel == 0), size = 0.75, colour = "black") + 
  coord_cartesian(ylim = c(0, 6e5)) + scale_x_continuous(breaks = c(1970, 1990, 2010)) + 
  gfplot::theme_pbs() + no_panel_gap + legend_bottom
ggsave("report/pollock/pollock_cond_rho.png", width = 6.5, height = 5)



ggplot(rr_out %>% filter(Peel > 0) %>% mutate(Peel = factor(Peel)), aes(Year, SSB)) + 
  geom_line(aes(group = Peel, colour = Peel)) + facet_wrap(~ Model) + 
  geom_text(data = rho_label, x = 1998, y = 0, vjust = "inward", hjust = "inward", aes(label = rho), parse = TRUE) +
  geom_line(data = filter(rr_out, Peel == 0), size = 0.75, colour = "black") + 
  coord_cartesian(ylim = c(0, 6e5), xlim = c(1998, 2020)) + #scale_x_continuous(breaks = c(1970, 1990, 2010)) + 
  gfplot::theme_pbs() + no_panel_gap + legend_bottom
ggsave("report/pollock/pollock_cond_rho2.png", width = 6.5, height = 5)





