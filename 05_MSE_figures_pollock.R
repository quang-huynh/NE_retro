library(MSEtool)
library(dplyr)
library(mfSCA)
library(ggplot2)

######## OM order
OM_names <- c("SS", "SWB", "SWF")

######## MP names
MPs <- c("Base", "Base_ra", "FlatSel", "FlatSel_ra", "MA", "75%FMSY")

######## Reference points
ref_pt <- readRDS(file = "pollock/pollock_ref_pt.rds")

ref_pt_plot <- do.call(rbind, lapply(ref_pt, function(i) i[1, ])) %>% as.data.frame()
ref_pt_plot$OM <- OM_names

######## Load MSE
MSE <- lapply(1:3, function(x) {
  y <- readRDS(paste0("pollock/MSE_pollock_NR", x, ".rds"))
  
  y@MPs[1:5] <- c("Base", "Base_ra", "FlatSel", "FlatSel_ra", "MA")
  
  y@OM$SSBMSY <- ref_pt_plot$SSBMSY[x]
  y@B_BMSY <- y@SSB/y@OM$SSBMSY
  
  y@OM$FMSY <- ref_pt_plot$FMSY[x]
  y@F_FMSY <- y@FM/y@OM$FMSY 
  
  y@OM$RefY <- ref_pt_plot$OY[x]
  
  return(y)
})


######## Plot OM
resOM <- plot_OM(MSE, MPs = MPs, 0.10, 0.90, yfilter = 10, OM_names = OM_names)

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
res[[1]] <- plot_EM(MSE, MPs = MPs[1:4], OM_names = OM_names) %>% 
  lapply(function(x) {
    x$EM <- ifelse(substr(x$MP, nchar(x$MP)-1, nchar(x$MP)) == "ra", 
                   substr(x$MP, 1, nchar(x$MP)-3), x$MP)
    return(x)
  })

res[[2]] <- plot_EM(MSE, MPs = "MA", OM_names = OM_names) %>% lapply(function(x) cbind(x, EM = "Base"))
res[[3]] <- plot_EM(MSE, MPs = "MA", OM_names = OM_names, model_index = 2) %>% 
  lapply(function(x) cbind(x, EM = "FlatSel"))

res[[4]] <- plot_EM(MSE, MPs = "75%FMSY", OM_names = OM_names)
res[[4]][[3]] <- cbind(res[[4]][[3]], EM = "Base")

res[[5]] <- plot_EM(MSE, MPs = "75%FMSY", OM_names = OM_names, model_index = 2)
res[[5]][[3]] <- cbind(res[[5]][[3]], EM = "FlatSel")

out <- lapply(c("F_FMSY", "B_BMSY", "rho"), function(x) {
  lapply(res, getElement, x) %>% do.call(rbind, .) %>% mutate(MP = factor(MP, levels = MPs))
})

ggplot(out[[3]] %>% filter(Year > 2020), aes(Year, y = `50%`, ymin = `25%`, ymax = `75%`, shape = EM)) + facet_grid(OM ~ MP, scales = "free_y") + 
  geom_hline(yintercept = 0, linetype = 3) + 
  geom_point(size = 1.25) + geom_linerange(size = 0.5) + 
  ylab(expression(rho[SSB])) + theme(legend.position = "bottom") +
  scale_x_continuous(breaks = c(2020, 2040, 2060)) + scale_shape_manual(values = c(16, 1)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("report/pollock/MSE_EM_rho.png", height = 4, width = 6.5)


######## Plot Index
#Ind <- plot_Index(MSE, MPs = MPs[1:5], ind.names = c("Spring", "Fall"), yfilter = 10)
#head(Ind[[1]])
#
#
#ggplot(Ind[[1]], aes(Year, y = `50%`, ymin = `25%`, ymax = `75%`)) + facet_grid(OM ~ MP, scales = "free_y") + 
#  #geom_hline(yintercept = 0, col = "grey") + 
#  geom_vline(xintercept = MSE[[1]]@OM$CurrentYr[1], linetype = 3) +
#  geom_ribbon(fill = "grey80") + geom_line(size = 0.9) +
#  gfplot::theme_pbs() + theme_extra + ylab("NEFSC Fall survey")
#
#ggplot(Ind[[2]], aes(Year, y = `50%`, ymin = `25%`, ymax = `75%`)) + facet_grid(OM ~ MP, scales = "free_y") + 
#  #geom_hline(yintercept = 0, col = "grey") + 
#  geom_vline(xintercept = MSE[[1]]@OM$CurrentYr[1], linetype = 3) +
#  geom_ribbon(fill = "grey80") + geom_line(size = 0.9) +
#  gfplot::theme_pbs() + theme_extra + ylab("NEFSC Spring survey")







######## Performance metric
PM_fn <- function(x, y) {
  data.frame(MP = x@MPs, 
             PNOF = PNOF(x)@Mean, 
             PB50 = P50(x)@Mean, 
             POY = LTY(x, Ref = 1, Yrs = c(1, 50))@Mean,
             STOY = LTY(x, Ref = 1, Yrs = c(1, 20))@Mean,
             stringsAsFactors = FALSE)  %>% 
    filter(MP != "FMSYref" & MP != "FMSYref75" & MP != "NFref") %>%
    mutate(OM = y) %>% 
    reshape2::melt(id.vars = c("MP", "OM"), variable.name = "PM")
}

pm <- Map(PM_fn, x = MSE, y = OM_names) %>% do.call(rbind, .) %>% filter(PM != "STOY") %>% 
  mutate(MP = factor(MP, levels = MPs),
         label = ifelse(value > 0.99, ">99", ifelse(value < 0.01, "<1", round(100 * value, 0))) %>% paste0("%"))

ggplot(pm %>% filter(PM != "STOY"), aes(PM, value, shape = PM)) + facet_grid(OM ~ MP) + 
  geom_text(size = 2.5, nudge_y = 0.1, aes(label = label)) + 
  geom_point(size = 1.5) + geom_linerange(size = 0.25, ymin = 0, aes(ymax = value)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom + theme(axis.text.x = element_blank()) + 
  coord_cartesian(ylim = c(0, 1.1)) + scale_y_continuous(breaks = seq(0, 1, by = 0.25))
ggsave("report/pollock/PM.png", height = 4, width = 6.5)

# Trade - off
pm2 <- reshape2::dcast(pm, MP + OM ~ PM)
ggplot(pm2, aes(PB50, POY, label = MP)) + facet_wrap(~OM) + geom_point() + ggrepel::geom_text_repel()

Yield_fn <- function(MSE, Yrs = c(1, 10), MPs, OM_names = "OM", Cbias = 1) {
  geo_mean <- function(x) (sum(log(x))/length(x)) %>% exp()
  CC <- (MSE@C[, match(MPs, MSE@MPs), Yrs[1]:Yrs[2], drop = FALSE]/Cbias) %>% apply(2, geo_mean)
  data.frame(MP = MPs, OM = OM_names, MeanC = CC)
}

color_fn <- function(x) {
  y <- rep(0, length(x))
  y[grep("Base", x)] <- 1
  y[grep("FlatSel", x)] <- 2
  return(y)
}

ypm <- Map(Yield_fn, MSE = MSE, OM_names = OM_names, MoreArgs = list(MPs = MPs)) %>% do.call(rbind, .) %>% 
  left_join(pm2, by = c("MP", "OM")) %>% 
  mutate(OM = factor(OM, levels = OM_names), col = color_fn(MP), 
         label = ifelse(PB50 > 0.99, ">99", ifelse(PB50 < 0.01, "<1", round(100 * PB50, 0))) %>% paste0(MP, " (", ., "%)"),
         PNOF = 100 * PNOF)

ggplot(ypm, aes(PNOF, MeanC, label = label, colour = factor(col), shape = factor(col))) + facet_grid(OM ~ ., scales = "free_y") + 
  geom_point(size = 2) + ggrepel::geom_text_repel(size = 2.5) + 
  scale_shape_manual(values = c(8, 1, 16)) + coord_cartesian(xlim = c(50, 100), ylim = c(10e3, 30e3)) + 
  xlab("PNOF (%)") + ylab("Observed short-term catch") + 
  scale_colour_manual(values = c("black", "blue", "red")) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom + no_legend
ggsave("report/pollock/pm_tradeoff.png", width = 3, height = 5)



#ypm <- Map(Yield_fn, MSE = MSE, OM_names = OM_names) %>% do.call(rbind, .) %>% 
#  left_join(pm2, by = c("MP", "OM")) %>% 
#  mutate(OM = factor(OM, levels = OM_names), col = color_fn(MP), 
#         label = ifelse(PNOF > 0.99, ">99", ifelse(PNOF < 0.01, "<1", round(100 * PNOF, 0))) %>% paste0(MP, " (", ., "%)"),
#         PB50 = 100 * PB50)
#
#ggplot(ypm, aes(PB50, MeanC, label = label, colour = factor(col), shape = factor(col))) + facet_grid(OM~., scales = "free_y") + 
#  geom_point(size = 2) + ggrepel::geom_text_repel(size = 2.5) + 
#  scale_shape_manual(values = c(8, 1, 16)) + coord_cartesian(xlim = c(0, 100), ylim = c(15e3, 30e3)) + 
#  xlab("PB50 (%)") + ylab("Observed short-term catch") + 
#  gfplot::theme_pbs() + no_panel_gap + legend_bottom + no_legend
#ggsave("report/pollock/pm_tradeoff.png", width = 3, height = 5)



#ypm <- Map(Yield_fn, MSE = MSE, OM_names = OM_names) %>% do.call(rbind, .) %>% 
#  left_join(pm, by = c("MP", "OM")) %>% mutate(col = color_fn(MP))
#ggplot(ypm %>% filter(PM == "PNOF" | PM == "PB50"), aes(value, MeanC, label = MP, colour = factor(col), shape = factor(col))) + 
#  facet_grid(OM ~ PM) + geom_point() + 
#  ggrepel::geom_text_repel(size = 2.5) + 
#  scale_shape_manual(values = c(8, 1, 16)) + coord_cartesian(ylim = c(15e3, 30e3)) + xlab("Performance metric") + 
#  ylab("Observed short-term catch") + 
#  gfplot::theme_pbs() + no_panel_gap + legend_bottom + no_legend
#ggsave("report/pollock/pm_tradeoff.png", width = 5, height = 5)



######## Indicators - MPs
s_CAA_hist <- get_pollock_base()$data$s_CAA
indicators_raw <- lapply(1:length(MSE), get_indicators, ind_interval = 3, MSE = MSE, 
                         MPs = MPs, s_CAA_hist = s_CAA_hist,
                         OM_names = OM_names,
                         mat_age = c(0.88, 0.294, 0.643, 0.887, 0.971, 0.993, 0.998, 1, 1))
indicators <- do.call(rbind, lapply(indicators_raw, getElement, 1)) #%>% mutate(value = ifelse(grepl("mu", Ind), exp(value), value))

indicators_cast <- indicators %>% reshape2::dcast(Year + MP + OM + Sim ~ Ind, value.var = "value") 

######### LDA
rr <- lda(OM ~ Ind_1_mu * Year + SSB_rho * Year + PMat_1_mu * Year + MAge_1_mu * Year, 
          data = filter(indicators_cast, Year > 2018 & MP == "Base_ra") %>% mutate(Year = factor(Year)))


#### All MPs for rho vs. Prop mature for 3 years
Yind <- c(2024, 2030, 2036)
rr_pred <- rr <- rr_CV <- list()
data_list <- lapply(1:6, function(i) filter(indicators_cast, Year > 2018 & MP == MPs[i]) %>% mutate(Year = factor(Year)))
for(i in 1:length(MPs[1:6])) {
  rr[[i]] <- MASS::lda(OM ~ PMat_1_mu * Year + SSB_rho * Year, data = data_list[[i]], method = "mle")
  rr_CV[[i]] <- MASS::lda(OM ~ PMat_1_mu * Year + SSB_rho * Year, data = data_list[[i]], method = "mle", CV = TRUE)
  rr_grid <- expand.grid(Year = table(indicators_cast$Year)[-1] %>% names(),
                         PMat_1_mu = seq(-0.4, 0, 0.01), SSB_rho = seq(-1, 1, 0.1)) 
  
  rr_pred[[i]] <- predict(rr[[i]], newdata = rr_grid) %>% getElement("posterior") %>% as.data.frame() %>% cbind(rr_grid) %>%
    mutate(MP = MPs[i])
}
rr_pred2 <- do.call(rbind, rr_pred) %>% mutate(MP = factor(MP, levels = MPs)) %>% 
  filter(match(Year, Yind, nomatch = 0) %>% as.logical())

# SVD
lapply(rr, getElement, 'svd')
lapply(rr, function(x) x$svd^2/sum(x$svd^2))

misclass <- Map(function(x, y, MP) {
  yy <- data.frame(misclass = abs(as.numeric(x$class) - as.numeric(y$OM)) %>% as.logical(), Year = y$Year)
  summarise(group_by(yy, Year), class_correct = round(1 - sum(misclass)/length(misclass), 2) %>% format()) %>% mutate(MP = MP)
}, x = rr_CV, y = data_list, MP = MPs) %>% do.call(rbind, .) %>% filter(match(Year, Yind, nomatch = 0) %>% as.logical()) %>%
  mutate(MP = factor(MP, MPs))

ggplot(indicators_cast %>% filter(match(Year, Yind, nomatch = 0) %>% as.logical()) %>%
         mutate(MP = factor(MP, levels = MPs)), aes(PMat_1_mu, SSB_rho)) + 
  facet_grid(Year ~ MP) + geom_hline(yintercept = 0, linetype = 3) + geom_point(aes(colour = OM), alpha = 0.6) +
  geom_contour(data = rr_pred2, breaks = 0.5, colour = "black", aes(z = SWF)) + 
  #metR::geom_label_contour(data = rr_pred2, breaks = 0.5, colour = "black", aes(z = SWF)) +
  geom_text(data = misclass, aes(label = class_correct), x = -0.05, y = 1, vjust = "inward", hjust = "inward") +
  coord_cartesian(xlim = c(-0.3, -0.05), ylim = c(-0.25, 1)) + labs(x = expression(P[mat]), y = expression(rho[SSB])) +
  scale_x_continuous(breaks = c(-0.3, -0.2, -0.1)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom
ggsave("report/pollock/indicators_LDA2_3y.png", width = 6.5, height = 4)


# Time series plots
indicators_ts <- summarise(group_by(indicators, Ind, Year, MP, OM),
                           y = median(value, na.rm = TRUE),
                           ymin = quantile(value, 0.25, na.rm = TRUE), ymax = quantile(value, 0.75, na.rm = TRUE)) %>% filter(Year > 2020)

# rho, Cat, Ind
#ggplot(filter(indicators_ts, !grepl("MAge", Ind)), aes(x = Year, shape = OM, colour = OM)) + facet_grid(Ind ~ MP, scales = "free_y") + 
#  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = OM), alpha = 0.1) + 
#  geom_point(aes(y = y)) + geom_line(aes(y = y)) + 
#  xlab("Year") + ylab("Indicator Value") + 
#  gfplot::theme_pbs() + no_panel_gap + legend_bottom
#
## Mean Age
#ggplot(filter(indicators_ts, grepl("MAge", Ind)), aes(x = Year, shape = OM, colour = OM)) + facet_grid(Ind ~ MP, scales = "free_y") + 
#  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = OM), alpha = 0.1) + 
#  geom_point(aes(y = y)) + geom_line(aes(y = y)) + 
#  xlab("Year") + ylab("Indicator Value") +
#  gfplot::theme_pbs() + no_panel_gap + legend_bottom

# TS to plot
Ind_order <- c("SSB_rho", "Cat_mu", "Ind_1_mu", "MAge_1_mu", "PMat_1_mu")
Ind_order2 <- c("rho", "Catch", "Index", "Mean Age", "Prop. Mature")
indicators_plot <- filter(indicators_ts, match(Ind, Ind_order, nomatch = 0) %>% as.logical()) %>%
  mutate(MP = factor(MP, levels = c("base", "base_ra", "flatsel", "flatsel_ra", "ma", "75%FMSY"))) %>%
  left_join(data.frame(Ind = Ind_order, Ind2 = factor(Ind_order2, levels = Ind_order2))) %>%
  mutate(Ind = factor(Ind, levels = c("SSB_rho", "Cat_slp", "Cat_mu", "Ind_1_slp", "MAge_1_slp"))) #%>% filter(Year >= 2024)
ggplot(indicators_plot, aes(x = Year, shape = OM)) + 
  geom_hline(data = data.frame(Ind2 = factor("rho", levels = Ind_order2), yy = 0), aes(yintercept = yy), linetype = 2) + 
  geom_point(aes(colour = OM, y = y)) + geom_line(aes(colour = OM, y = y)) + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = OM), alpha = 0.3) +
  facet_grid(Ind2 ~ MP, scales = "free_y") +  
  xlab("Year") + ylab("Predicted indicator value") + scale_x_continuous(breaks = c(2020, 2040, 2060)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom
ggsave("report/pollock/indicator_ts.png", height = 7, width = 7)




#### Plot recruitment
SRA <- lapply(1:3, function(i) readRDS(paste0("pollock/SRA_NR", i, ".rds")))

plot_recruitment <- function(MSE, SRA, OM_names, MPs) {
  #Rhist <- lapply(SRA@Misc, getElement, "R") %>% do.call(rbind, .)
  SSBpred <- MSE@SSB[, match(MPs, MSE@MPs), , drop = FALSE] %>% aperm(c(1, 3, 2))
  Rdevpred <- SRA@OM@cpars$Perr_y[, -c(1:(SRA@OM@maxage - 1 + MSE@nyears))] %>% as.vector()
  Arec <- vapply(SRA@Misc, getElement, numeric(1), "Arec")
  Brec <- vapply(SRA@Misc, getElement, numeric(1), "Brec")
  Rpred <- (Arec * SSBpred / (1 + Brec * SSBpred) * Rdevpred) %>% 
    structure(dimnames = list(Sim = 1:MSE@nsim, Time = 1:MSE@proyears, MP = MPs)) %>%
    reshape2::melt(value.name = "R") %>% mutate(Year = Time + MSE@OM$CurrentYr, Time = NULL, OM = OM_names)
  Rhist <- data.frame(R = SRA@mean_fit$report$R, Year = MSE@OM$CurrentYr[1] + 1 - (MSE@nyears):0, OM = OM_names)
  return(list(Rpred = Rpred, Rhist = Rhist))
}

rec <- Map(plot_recruitment, MSE = MSE, SRA = SRA, OM_names = OM_names,
           MoreArgs = list(MPs = c("base", "base_ra", "flatsel", "flatsel_ra", "ma", "75%FMSYref")))

recpred <- lapply(rec, getElement, "Rpred") %>% do.call(rbind, .) %>% group_by(MP, OM, Year) %>% summarise(R = mean(R)) %>%
  mutate(OM = factor(OM, levels = OM_names))
rhist <- lapply(rec, getElement, "Rhist") %>% do.call(rbind, .) %>% mutate(OM = factor(OM, levels = OM_names)) %>%
  filter(Year >= 2009)

ggplot(recpred, aes(Year, R)) + geom_line(aes(colour = MP)) + geom_point(aes(colour = MP, shape = MP)) + 
  facet_grid(OM ~ .) + geom_line(data = rhist) +
  geom_vline(xintercept = 2018, linetype = 2) + labs(y = "Recruitment") + 
  gfplot::theme_pbs() + no_panel_gap + legend_bottom
ggsave("report/pollock/mean_recruitment.png", width = 3, height = 4.5)  
