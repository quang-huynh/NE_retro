library(MSEtool)
library(dplyr)
library(mfSCA)
library(ggplot2)

# Calc ref points
SRA <- lapply(1:3, function(i) readRDS(paste0("pollock/SRA_NR", i, ".rds")))
med_rec <- vapply(SRA, function(i) quantile(i@mean_fit$report$R, 0.5), numeric(1))

######## Reference points
ref_pt <- lapply(1:3, function(i) calc_refpt(SRA[[i]], med_rec = med_rec[i]))
#saveRDS(ref_pt, file = "pollock/pollock_ref_pt.rds")

ref_pt_plot <- do.call(rbind, lapply(ref_pt, function(i) i[1, ])) %>% as.data.frame()
ref_pt_plot$OM <- paste0("NR", 1:3)

######## Load MSE
MSE <- lapply(1:3, function(i) readRDS(paste0("pollock/MSE_pollock_NR", i, ".rds")))
MSE <- lapply(1:3, function(x) {
  y <- MSE[[x]]
  y@OM$SSBMSY <- ref_pt_plot$SSBMSY[x]
  y@B_BMSY <- y@SSB/y@OM$SSBMSY
  
  y@OM$FMSY <- ref_pt_plot$FMSY[x]
  y@F_FMSY <- y@FM/y@OM$FMSY 
  
  y@OM$RefY <- ref_pt_plot$OY[x]
  
  return(y)
})


######## Plot OM
resOM <- plot_OM(MSE, MPs = c("base", "base_ra", "flatsel", "flatsel_ra", "ma", "75%FMSY"), 0.10, 0.90, yfilter = 10)

no_legend <- theme(legend.position = "none")
no_panel_gap <- theme(panel.spacing = unit(0, "in"))
legend_bottom <- theme(legend.position = "bottom")

ggplot(resOM[[1]], aes(Year, y = med, ymin = low, ymax = high)) + facet_grid(OM ~ MP, scales = "free_y") + 
  #geom_hline(yintercept = 0, col = "grey") + 
  geom_vline(xintercept = MSE[[1]]@OM$CurrentYr[1], linetype = 3) +
  geom_ribbon(fill = "grey80") + geom_line(size = 0.9) + 
  geom_hline(data = ref_pt_plot, aes(yintercept = FMSY), linetype = 2, size = 0.5) +
  ylab("Fishing mortality") + scale_y_continuous(n.breaks = 4, expand = c(0.01, 0.01)) +
  gfplot::theme_pbs() + no_panel_gap
ggsave("report/pollock/MSE_OM_F.png", height = 3.5, width = 8.5)

ggplot(resOM[[2]], aes(Year, y = med, ymin = low, ymax = high)) + facet_grid(OM ~ MP, scales = "free_y") + 
  #geom_hline(yintercept = 0, col = "grey") + 
  geom_vline(xintercept = MSE[[1]]@OM$CurrentYr[1], linetype = 3) +
  geom_ribbon(fill = "grey80") + geom_line(size = 0.9) + 
  geom_hline(data = ref_pt_plot, aes(yintercept = SSBMSY), linetype = 2) +
  geom_hline(data = ref_pt_plot, aes(yintercept = 0.5 * SSBMSY), linetype = 2) +
  geom_hline(yintercept = 0, alpha = 0) +
  ylab("SSB") + scale_y_continuous(labels = scales::comma, expand = c(0.01, 0.01)) +
  gfplot::theme_pbs() + no_panel_gap
ggsave("report/pollock/MSE_OM_SSB.png", height = 3.5, width = 8.5)

ggplot(resOM[[3]], aes(Year, y = med, ymin = low, ymax = high)) + facet_grid(OM ~ MP, scales = "free_y") + 
  #geom_hline(yintercept = 0, col = "grey") + 
  geom_vline(xintercept = MSE[[1]]@OM$CurrentYr[1], linetype = 3) +
  geom_ribbon(fill = "grey80") + geom_line(size = 0.9) + 
  geom_hline(data = ref_pt_plot, aes(yintercept = OY), linetype = 3) +
  ylab("Catch") + scale_y_continuous(labels = scales::comma, n.breaks = 4, expand = c(0.01, 0.01)) +
  gfplot::theme_pbs() + no_panel_gap
ggsave("report/pollock/MSE_OM_C.png", height = 3.5, width = 8.5)


####### Plot rho 
res <- list()
res[[1]] <- plot_EM(MSE, MPs = c("base", "base_ra", "flatsel", "flatsel_ra")) %>% 
  lapply(function(x) {
    x$EM <- ifelse(substr(x$MP, nchar(x$MP)-1, nchar(x$MP)) == "ra", 
                   substr(x$MP, 1, nchar(x$MP)-3), x$MP)
    return(x)
  })

res[[2]] <- plot_EM(MSE, MPs = "ma") %>% lapply(function(x) cbind(x, EM = "base"))
res[[3]] <- plot_EM(MSE, MPs = "ma", model_index = 2) %>% 
  lapply(function(x) cbind(x, EM = "flatsel"))

res[[4]] <- plot_EM(MSE, MPs = "75%FMSY")
res[[4]][[3]] <- cbind(res[[4]][[3]], EM = "base")

res[[5]] <- plot_EM(MSE, MPs = "75%FMSY", model_index = 2)
res[[5]][[3]] <- cbind(res[[5]][[3]], EM = "flatsel")

out <- lapply(c("F_FMSY", "B_BMSY", "rho"), function(x) {
  y <- do.call(rbind, lapply(res, getElement, x))
  y$MP <- factor(y$MP, levels = c("base", "base_ra", "flatsel", "flatsel_ra", "ma", "75%FMSY"))
  return(y)
})

ggplot(out[[3]] %>% filter(Year > 2020), aes(Year, y = `50%`, ymin = `25%`, ymax = `75%`, shape = EM)) + facet_grid(OM ~ MP, scales = "free_y") + 
  geom_hline(yintercept = 0, linetype = 3) + 
  geom_point(size = 1.25) + geom_linerange(size = 0.5) + 
  ylab(expression(rho[SSB])) + theme(legend.position = "bottom") +
  scale_x_continuous(breaks = c(2020, 2040, 2060)) + scale_shape_manual(values = c(16, 1)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom
ggsave("report/pollock/MSE_EM_rho.png", height = 4, width = 6.5)


######## Plot Index
Ind <- plot_Index(MSE, MPs = c("base", "base_ra", "flatsel", "flatsel_ra", "ma"), 
                  ind.names = c("Spring", "Fall"), yfilter = 10)
head(Ind[[1]])


ggplot(Ind[[1]], aes(Year, y = `50%`, ymin = `25%`, ymax = `75%`)) + facet_grid(OM ~ MP, scales = "free_y") + 
  #geom_hline(yintercept = 0, col = "grey") + 
  geom_vline(xintercept = MSE[[1]]@OM$CurrentYr[1], linetype = 3) +
  geom_ribbon(fill = "grey80") + geom_line(size = 0.9) +
  gfplot::theme_pbs() + theme_extra + ylab("NEFSC Fall survey")

ggplot(Ind[[2]], aes(Year, y = `50%`, ymin = `25%`, ymax = `75%`)) + facet_grid(OM ~ MP, scales = "free_y") + 
  #geom_hline(yintercept = 0, col = "grey") + 
  geom_vline(xintercept = MSE[[1]]@OM$CurrentYr[1], linetype = 3) +
  geom_ribbon(fill = "grey80") + geom_line(size = 0.9) +
  gfplot::theme_pbs() + theme_extra + ylab("NEFSC Spring survey")







######## Performance metric
PM_fn <- function(x, y) {
  out <- data.frame(MP = x@MPs, 
                    PNOF = PNOF(x)@Mean, 
                    PB50 = P50(x)@Mean, 
                    POY = LTY(x, Ref = 1, Yrs = c(1, 50))@Mean,
                    STOY = LTY(x, Ref = 1, Yrs = c(1, 20))@Mean,
                    stringsAsFactors = FALSE)  %>% 
    filter(MP != "FMSYref" & MP != "FMSYref75" & MP != "NFref")
  out$OM <- paste0("NR", y)
  
  out %>% reshape2::melt(id.vars = c("MP", "OM"), variable.name = "PM")
}

pm <- do.call(rbind, Map(PM_fn, x = MSE, y = 1:3)) %>% filter(PM != "STOY") %>% 
  mutate(MP = factor(MP, levels = c("base", "base_ra", "flatsel", "flatsel_ra", "ma", "75%FMSY")),
         label = ifelse(value > 0.99, ">99", ifelse(value < 0.01, "<1", round(100 * value, 0))) %>% paste0("%"))

ggplot(pm %>% filter(PM != "STOY"), aes(PM, value, shape = PM)) + facet_grid(OM ~ MP) + 
  geom_text(size = 2.5, nudge_y = 0.1, aes(label = label)) + 
  geom_point(size = 1.5) + geom_linerange(size = 0.25, ymin = 0, aes(ymax = value)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom + theme(axis.text.x = element_blank()) + 
  coord_cartesian(ylim = c(0, 1.1)) + scale_y_continuous(breaks = seq(0, 1, by = 0.25))
ggsave("report/pollock/PM.png", height = 4, width = 6.5)




######## Indicators - MPs
s_CAA_hist <- get_pollock_base()$data$s_CAA
indicators_raw <- lapply(1:length(MSE), get_indicators, ind_interval = 6, MSE = MSE, 
                         MPs = c("base", "base_ra", "flatsel", "flatsel_ra", "ma", "75%FMSY"), s_CAA_hist = s_CAA_hist,
                         mah_ind = c("SSB_rho", "Cat_slp", "Cat_mu", "Ind_1_slp", "MAge_1_slp"))
indicators <- do.call(rbind, lapply(indicators_raw, getElement, 1))

######## Indicators - refMPs
Year_vec <- vapply(MSE[[1]]@Misc$Data[[1]]@Misc[[1]]$diagnostic, getElement, numeric(1), "Year")
indicators_raw_ref <- lapply(1:length(MSE_ref), get_indicators, ind_interval = 6, MSE = MSE_ref, 
                         MPs = MSE_ref[[1]]@MPs, s_CAA_hist = s_CAA_hist,
                         mah_ind = c("SSB_rho", "Cat_slp", "Cat_mu", "Ind_1_slp", "MAge_1_slp"), Year_vec = Year_vec)
indicators_ref <- do.call(rbind, lapply(indicators_raw_ref, getElement, 1))
indicators_ts_ref <- summarise(group_by(indicators_ref, Ind, Year, MP, OM),
                           y = median(value, na.rm = TRUE),
                           ymin = quantile(value, 0.25, na.rm = TRUE), ymax = quantile(value, 0.75, na.rm = TRUE)) %>% filter(Year > 2020)

# Mahalanobis
#mah <- do.call(rbind, lapply(indicators_raw, getElement, 2))
#
#mah_gr <- expand.grid(Year = unique(mah$Year), MP = unique(mah$MP))
#mah_gr$pval <- vapply(1:nrow(mah_gr), function(x) {
#  xx <- mah %>% filter(MP == mah_gr$MP[x] & Year == mah_gr$Year[x])
#  pval <- aov(value ~ OM, data = xx) %>% summary() %>% getElement(1) %>% getElement("Pr(>F)") %>% getElement(1)
#  return(pval)
#}, numeric(1))
#mah_gr$power <- vapply(1:nrow(mah_gr), function(x) {
#  xx <- mah %>% filter(MP == mah_gr$MP[x] & Year == mah_gr$Year[x])
#  xx_aov <- aov(value ~ OM, data = xx)
#  power.anova.test(groups = xx_aov$rank, n = nrow(xx)/xx_aov$rank, between.var = summary(xx_aov)[[1]]$`Mean Sq`[1], 
#                   within.var = summary(xx_aov)[[1]]$`Mean Sq`[2], sig.level = 0.1)$power
#}, numeric(1))
#
##ggplot(mah_gr, aes(x = Year, y = power)) + geom_point() + geom_line() + facet_wrap(~MP)
#
#mah_summary <- mah %>% group_by(MP, OM, Year) %>% 
#  summarise(y = median(value, na.rm = TRUE), ymin = quantile(value, 0.25, na.rm = TRUE), ymax = quantile(value, 0.75, na.rm = TRUE))
#ggplot(mah_summary, aes(x = Year, y = y)) + geom_point(aes(colour = OM, shape = OM)) + geom_line(aes(colour = OM, shape = OM)) + 
#  facet_wrap(~ factor(MP, levels = c("base", "base_ra", "flatsel", "flatsel_ra", "ma")), scales = "free") + 
#  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = OM, colour = OM), alpha = 0.1) + 
#  ylab("Mahalanobis distance") +
#  geom_point(data = mah_gr %>% filter(pval <= 0.05), y = 0, shape = 3) + 
#  scale_x_continuous(breaks = c(2020, 2040, 2060)) +
#  gfplot::theme_pbs() + legend_bottom
#ggsave("report/pollock/mahalanobis.png", height = 4.5, width = 7)


# Time series plots
indicators_ts <- summarise(group_by(indicators, Ind, Year, MP, OM),
                           y = median(value, na.rm = TRUE),
                           ymin = quantile(value, 0.25, na.rm = TRUE), ymax = quantile(value, 0.75, na.rm = TRUE)) %>% filter(Year > 2020)

# rho, Cat, Ind
ggplot(filter(indicators_ts, !grepl("MAge", Ind)), aes(x = Year, shape = OM, colour = OM)) + facet_grid(Ind ~ MP, scales = "free_y") + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = OM), alpha = 0.1) + 
  geom_point(aes(y = y)) + geom_line(aes(y = y)) + 
  xlab("Year") + ylab("Indicator Value") + 
  gfplot::theme_pbs() + no_panel_gap + legend_bottom

# Mean Age
ggplot(filter(indicators_ts, grepl("MAge", Ind)), aes(x = Year, shape = OM, colour = OM)) + facet_grid(Ind ~ MP, scales = "free_y") + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = OM), alpha = 0.1) + 
  geom_point(aes(y = y)) + geom_line(aes(y = y)) + 
  xlab("Year") + ylab("Indicator Value") +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom


ggplot(filter(indicators_ts, grepl("1", Ind) | grepl("rho", Ind) | grepl("Cat", Ind)), 
       aes(x = Year, shape = OM, colour = OM)) + 
  facet_grid(Ind ~ factor(MP, levels = c("base", "base_ra", "flatsel", "flatsel_ra", "ma")), scales = "free_y") + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = OM), alpha = 0.1) + 
  geom_point(aes(y = y)) + geom_line(aes(y = y)) + 
  xlab("Year") + ylab("Indicator Value") + scale_x_continuous(breaks = c(2020, 2040, 2060)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom

# TS to plot
indicators_plot <- indicators_ts %>%
  filter(match(Ind, c("SSB_rho", "Cat_slp", "Cat_mu", "Ind_1_slp", "MAge_1_slp"), nomatch = 0) %>% as.logical()) %>%
  mutate(MP = factor(MP, levels = c("base", "base_ra", "flatsel", "flatsel_ra", "ma", "75%FMSY"))) %>%
  mutate(Ind = factor(Ind, levels = c("SSB_rho", "Cat_slp", "Cat_mu", "Ind_1_slp", "MAge_1_slp"))) #%>% filter(Year >= 2024)
ggplot(indicators_plot, aes(x = Year, shape = OM, colour = OM)) + 
  facet_grid(Ind ~ MP, scales = "free_y") + 
  geom_hline(data = data.frame(Ind = factor("SSB_rho", levels = c("SSB_rho", "Cat_slp", "Cat_mu", "Ind_1_slp", "MAge_1_slp")), yy = 0), aes(yintercept = yy), linetype = 2) + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = OM), alpha = 0.1) + 
  geom_point(aes(y = y)) + geom_line(aes(y = y)) + 
  xlab("Year") + ylab("Indicator Value") + scale_x_continuous(breaks = c(2020, 2040, 2060)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom
ggsave("report/pollock/indicator_ts.png", height = 7, width = 7)


##### TS of indicators relative to FMSY75

relative_ind <- filter(indicators, MP == "75%FMSY") %>% mutate(value2 = value, value = NULL, MP = NULL) %>% 
  right_join(filter(indicators, MP != "75%FMSY")) %>% filter(grepl("1", Ind) | grepl("Cat", Ind) | grepl("SSB", Ind)) %>%
  mutate(diff = value - value2, value = NULL, value2 = NULL) %>% filter(!is.na(diff)) %>% group_by(Ind, Year, MP, OM) %>% 
  summarise(y = median(diff, na.rm = TRUE), ymin = quantile(diff, 0.25, na.rm = TRUE), ymax = quantile(diff, 0.75, na.rm = TRUE)) %>% 
  filter(Year > 2020) %>% mutate(MP = factor(MP, levels = c("base", "base_ra", "flatsel", "flatsel_ra", "ma")))
ggplot(relative_ind, aes(x = Year, shape = OM, colour = OM)) + 
  geom_hline(yintercept = 0, linetype = 2) + 
  geom_point(aes(y = y)) + geom_line(aes(y = y)) + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = OM), alpha = 0.1) +
  facet_grid(Ind ~ MP, scales = "free_y") +  
  xlab("Year") + ylab(latex2exp::TeX("$\\Delta$ Indicator")) + scale_x_continuous(breaks = c(2020, 2040, 2060)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom
ggsave("report/pollock/indicator_delta.png", height = 7, width = 7)


# Scatter plots 
indicators_scatter <- reshape2::dcast(indicators, formula = Year + Sim + MP + OM ~ Ind, value.var = "value")

x <- MSE[[1]]
MPs <- x@MPs[1:2]
x@F_FMSY[, 1:2, ] %>% structure(dimnames = list(1:x@nsim, MPs, 1:x@proyears)) %>% 
  reshape2::melt(varnames = c("sim", "MP", "time"), value.name = expression(F/F[MSY])) %>% 
  mutate(OM = paste0("NR", i))

F_FMSY <- lapply(1:length(MSE), function(x, MPs) MSE[[x]]@F_FMSY[, match(MPs, MSE[[x]]@MPs), , drop = FALSE] %>% 
                   structure(dimnames = list(1:MSE[[x]]@nsim, MPs, MSE[[x]]@nyears + 1:MSE[[x]]@proyears)) %>% 
                   reshape2::melt(varnames = c("Sim", "MP", "time"), value.name = "F.FMSY") %>% 
                   mutate(Year = time + MSE[[x]]@OM$CurrentYr[1] - MSE[[x]]@nyears, time = NULL) %>%
                   mutate(OM = paste0("NR", x)), MPs = "base")
F_FMSY <- do.call(rbind, F_FMSY)

ggplot(filter(indicators_scatter, MP == "base") %>% left_join(F_FMSY) %>% filter(!is.na(F.FMSY)), 
       aes(x = SSB_rho, y = F.FMSY, colour = OM , shape = OM)) + facet_wrap(~ Year) + geom_point(alpha = 0.9) + 
  geom_hline(yintercept = 1, linetype = 2) + geom_vline(xintercept = 0, linetype = 2) + labs(x = expression(rho[SSB]), y = expression(F/F[MSY])) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom
ggsave("report/pollock/indicator_scatter.png", height = 6, width = 6)

#ggplot(filter(indicators_scatter, MP == "base") %>% left_join(F_FMSY) %>% filter(!is.na(F.FMSY)), 
#       aes(x = SSB_rho, y = F.FMSY, colour = OM)) + facet_wrap(~ Year) + geom_point() + 
#  geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2) +
#  gfplot::theme_pbs() + no_panel_gap + legend_bottom


unique(indicators_scatter$Year)
ggplot(filter(indicators_scatter, Year == 2027), aes(x = SSB_rho, y = Cat_slp, shape = OM, colour = OM)) + geom_point() + 
  facet_wrap(~ factor(MP, levels = c("base", "base_ra", "flatsel", "flastel_ra", "ma")), scales = "free") + ggtitle("Year 2027") +
  gfplot::theme_pbs() + legend_bottom
ggsave("report/pollock/indicator_scatter.png", height = 4.5, width = 7)

#filter(indicators_scatter, Year == 2024) %>% ggplot(aes(x = Cat_mu, y = Cat_slp, shape = OM, colour = OM)) + geom_point() + 
#  facet_wrap(~MP) + gfplot::theme_pbs()

#filter(indicators_scatter, Year == 2024) %>% ggplot(aes(x = Cat_mu, y = Ind_1_mu, shape = OM, colour = OM)) + geom_point() + 
#  facet_wrap(~MP) + gfplot::theme_pbs()


##### Indicators of ref MPs
ggplot(indicators_ts_ref, aes(x = Year, shape = OM, colour = OM)) + facet_grid(Ind ~ MP, scales = "free_y") + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = OM), alpha = 0.1) + 
  geom_point(aes(y = y)) + geom_line(aes(y = y)) + 
  xlab("Year") + ylab("Indicator Value") + 
  gfplot::theme_pbs() + no_panel_gap + legend_bottom
