library(MSEtool)
library(dplyr)
library(mfSCA)
library(ggplot2)

######## OM order
OM_order <- c(1, 3, 2)
OM_names <- c("MC", "MCIM", "IM")

######## Reference points
ref_pt <- readRDS("GoM_cod/cod_ref_pt.rds")

ref_pt_plot <- lapply(ref_pt, function(i) lapply(i, function(ii) ii[1, ]) %>% do.call(rbind, .)) %>% do.call(rbind, .) %>%
  as.data.frame()
ref_pt_plot$OM <- rep(OM_names, 2)
ref_pt_plot$type <- rep(c("old", "new"), each = 3)
#ref_pt_plot <- ref_pt_plot[-4, ]

######## Load MSE
MSE <- lapply(1:3, function(i) readRDS(paste0("GoM_cod/MSE_cod_NR", i, ".rds")))
MSE <- lapply(1:3, function(x) {
  y <- MSE[[x]]
  y@OM$SSBMSY <- filter(ref_pt_plot, type == 'new')$SSBMSY[x]
  y@B_BMSY <- y@SSB/y@OM$SSBMSY
  
  y@OM$FMSY <- filter(ref_pt_plot, type == 'new')$FMSY[x]
  y@F_FMSY <- y@FM/y@OM$FMSY 
  
  y@OM$RefY <- filter(ref_pt_plot, type == 'new')$OY[x]
  
  return(y)
})

######## Plot OM
resOM <- plot_OM(MSE[OM_order], MPs = c("M02", "M02_ra", "MRAMP", "MRAMP_ra", "ma", "75%FMSY"), 0.10, 0.90, yfilter = 10, 
                 OM_names = OM_names[OM_order])

no_legend <- theme(legend.position = "none")
no_panel_gap <- theme(panel.spacing = unit(0, "in"))
legend_bottom <- theme(legend.position = "bottom")

ggplot(resOM[[1]], aes(Year, y = med, ymin = low, ymax = high)) + facet_grid(factor(OM) ~ MP, scales = "free_y") + 
  #geom_hline(yintercept = 0, col = "grey") + 
  geom_vline(xintercept = MSE[[1]]@OM$CurrentYr[1], linetype = 3) +
  geom_ribbon(fill = "grey80") + geom_line(size = 0.9) + 
  geom_hline(data = ref_pt_plot[4:6, ], aes(yintercept = FMSY), linetype = 2, size = 0.5) +
  ylab("Fishing mortality") + scale_y_continuous(expand = c(0.01, 0.01)) + 
  gfplot::theme_pbs() + no_panel_gap
ggsave("report/GoM_cod/MSE_OM_F.png", height = 3.5, width = 8.5)

ggplot(resOM[[2]], aes(Year, y = med, ymin = low, ymax = high)) + facet_grid(factor(OM) ~ MP, scales = "free_y") +  
  #geom_hline(yintercept = 0, col = "grey") + 
  geom_vline(xintercept = MSE[[1]]@OM$CurrentYr[1], linetype = 3) +
  geom_ribbon(fill = "grey80") + geom_line(size = 0.9) + 
  geom_hline(data = ref_pt_plot[4:6, ], aes(yintercept = SSBMSY), linetype = 2) +
  geom_hline(data = ref_pt_plot[4:6, ], aes(yintercept = 0.5 * SSBMSY), linetype = 2) +
  ylab("SSB") + scale_y_continuous(labels = scales::comma, expand = c(0.01, 0.01)) +
  gfplot::theme_pbs() + no_panel_gap
ggsave("report/GoM_cod/MSE_OM_SSB.png", height = 3.5, width = 8.5)

ggplot(resOM[[3]], aes(Year, y = med, ymin = low, ymax = high)) + facet_grid(factor(OM) ~ MP, scales = "free_y") + 
  #geom_hline(yintercept = 0, col = "grey") + 
  geom_vline(xintercept = MSE[[1]]@OM$CurrentYr[1], linetype = 3) +
  geom_ribbon(fill = "grey80") + geom_line(size = 0.9) + 
  geom_hline(data = ref_pt_plot[4:6, ], aes(yintercept = OY), linetype = 3) +
  ylab("Catch") + scale_y_continuous(labels = scales::comma, expand = c(0.01, 0.01)) +
  gfplot::theme_pbs() + no_panel_gap
ggsave("report/GoM_cod/MSE_OM_C.png", height = 3.5, width = 8.5)


####### Plot rho 
res <- list()
res[[1]] <- plot_EM(MSE[OM_order], MPs = c("M02", "M02_ra", "MRAMP", "MRAMP_ra"), OM_names = OM_names[OM_order]) %>% 
  lapply(function(x) {
    x$EM <- ifelse(substr(x$MP, nchar(x$MP)-1, nchar(x$MP)) == "ra", 
                   substr(x$MP, 1, nchar(x$MP)-3), x$MP)
    return(x)
  })

res[[2]] <- plot_EM(MSE[OM_order], MPs = "ma", OM_names = OM_names[OM_order]) %>% lapply(function(x) cbind(x, EM = "M02"))
res[[3]] <- plot_EM(MSE[OM_order], MPs = "ma", model_index = 2, OM_names = OM_names[OM_order]) %>% lapply(function(x) cbind(x, EM = "MRAMP"))

res[[4]] <- plot_EM(MSE[OM_order], MPs = "75%FMSY", model_index = 1, OM_names = OM_names[OM_order])
res[[4]][[3]] <- cbind(res[[4]][[3]], EM = "M02")

res[[5]] <- plot_EM(MSE[OM_order], MPs = "75%FMSY", model_index = 2, OM_names = OM_names[OM_order])
res[[5]][[3]] <- cbind(res[[5]][[3]], EM = "MRAMP")



out <- lapply(c("F_FMSY", "B_BMSY", "rho"), function(x) {
  lapply(res, getElement, x) %>% do.call(rbind, .) %>% 
    mutate(MP = factor(MP, levels = c("M02", "M02_ra", "MRAMP", "MRAMP_ra", "ma", "75%FMSY")))
})

ggplot(out[[3]] %>% filter(Year > 2020), aes(x = as.numeric(Year), y = `50%`, ymin = `25%`, ymax = `75%`, shape = EM)) + 
  facet_grid(factor(OM) ~ MP, scales = "free_y") + 
  geom_hline(yintercept = 0, linetype = 3) + 
  geom_point(size = 1.25) + geom_linerange(size = 0.5) + 
  xlab("Year") + ylab(expression(rho[SSB])) + theme(legend.position = "bottom") +
  scale_x_continuous(breaks = c(2020, 2040, 2060)) + scale_shape_manual(values = c(16, 1)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("report/GoM_cod/MSE_EM_rho.png", height = 4, width = 6.5)


######## Plot Index
#Ind <- plot_Index(MSE, MPs = c("M02", "M02_ra", "MRAMP", "MRAMP_ra", "ma"), 
#                  ind.names = c("Spring", "Fall", "MA"), yfilter = 10)
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
#  gfplot::theme_pbs() + theme_extra + ylab("MA survey")
#
#ggplot(Ind[[3]], aes(Year, y = `50%`, ymin = `25%`, ymax = `75%`)) + facet_grid(OM ~ MP, scales = "free_y") + 
#  #geom_hline(yintercept = 0, col = "grey") + 
#  geom_vline(xintercept = MSE[[1]]@OM$CurrentYr[1], linetype = 3) +
#  geom_ribbon(fill = "grey80") + geom_line(size = 0.9) +
#  gfplot::theme_pbs() + theme_extra + ylab("NEFSC Spring survey")









########  Performance metric
PM_fn <- function(x, y) {
  data.frame(MP = x@MPs, 
             PNOF = PNOF(x)@Mean, 
             PB50 = P50(x)@Mean, 
             POY = LTY(x, Ref = 1, Yrs = c(1, 50))@Mean,
             STOY = LTY(x, Ref = 1, Yrs = c(1, 20))@Mean,
             stringsAsFactors = FALSE) %>%
    filter(MP != "FMSYref" & MP != "FMSYref75" & MP != "NFref") %>% 
    mutate(OM = y) %>%
    reshape2::melt(id.vars = c("MP", "OM"), variable.name = "PM")
}
pm <- Map(PM_fn, x = MSE, y = OM_names) %>% do.call(rbind, .) %>% filter(PM != "STOY") %>% 
  mutate(MP = factor(MP, levels = c("M02", "M02_ra", "MRAMP", "MRAMP_ra", "ma", "75%FMSY")),
         label = ifelse(value > 0.99, ">99", ifelse(value < 0.01, "<1", round(100 * value, 0))) %>% paste0("%"))

ggplot(pm, aes(PM, value, shape = PM)) + facet_grid(OM ~ MP) + 
  geom_text(size = 2.5, nudge_y = 0.1, aes(label = label)) + 
  geom_point(size = 1.5) + geom_linerange(size = 0.25, ymin = 0, aes(ymax = value)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom + theme(axis.text.x = element_blank()) + 
  coord_cartesian(ylim = c(0, 1.1)) + scale_y_continuous(breaks = seq(0, 1, by = 0.25))
ggsave("report/GoM_cod/PM.png", height = 4, width = 6.5)

# Trade - off
pm2 <- reshape2::dcast(pm, MP + OM ~ PM)
ggplot(pm2, aes(PB50, POY, label = MP)) + facet_wrap(~OM) + geom_point() + ggrepel::geom_text_repel()

Yield_fn <- function(MSE, Yrs = c(1, 10), MPs = c("M02", "M02_ra", "MRAMP", "MRAMP_ra", "ma", "75%FMSY"),
                     OM_names = "OM", Cbias = 1) {
  geo_mean <- function(x) (sum(log(x))/length(x)) %>% exp()
  CC <- (MSE@C[, match(MPs, MSE@MPs), Yrs[1]:Yrs[2], drop = FALSE]/Cbias) %>% apply(2, geo_mean)
  data.frame(MP = MPs, OM = OM_names, MeanC = CC)
}

color_fn <- function(x) {
  y <- rep(0, length(x))
  y[grep("M02", x)] <- 1
  y[grep("MRAMP", x)] <- 2
  return(y)
}


ypm <- Map(Yield_fn, MSE = MSE[OM_order], OM_names = OM_names[OM_order], Cbias = c(2.25, 1.25, 1)) %>% do.call(rbind, .) %>% 
  left_join(pm2, by = c("MP", "OM")) %>% 
  mutate(OM = factor(OM, levels = OM_names[OM_order]), col = color_fn(MP), 
         label = ifelse(PB50 > 0.99, ">99", ifelse(PB50 < 0.01, "<1", round(100 * PB50, 0))) %>% paste0(MP, " (", ., "%)"),
         PNOF = 100 * PNOF)

#ypm <- Map(Yield_fn, MSE = MSE[OM_order], OM_names = OM_names[OM_order]) %>% do.call(rbind, .) %>% 
#  left_join(pm2, by = c("MP", "OM")) %>% 
#  mutate(OM = factor(OM, levels = OM_names[OM_order]), col = color_fn(MP), 
#         label = ifelse(PNOF > 0.99, ">99", ifelse(PNOF < 0.01, "<1", round(100 * PNOF, 0))) %>% paste0(MP, " (", ., "%)"),
#         PB50 = 100 * PB50)

ggplot(ypm, aes(PNOF, MeanC, label = label, colour = factor(col), shape = factor(col))) + facet_grid(OM~., scales = "free_y") + 
  geom_point(size = 2) + ggrepel::geom_text_repel(size = 2.5) + 
  scale_shape_manual(values = c(8, 1, 16)) + coord_cartesian(ylim = c(0, 1.1e3)) + xlab("PNOF (%)") + ylab("Observed short-term catch") +
  scale_y_continuous(breaks = seq(0, 1000, 250)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom + no_legend
ggsave("report/GoM_cod/pm_tradeoff.png", width = 3, height = 5)


#ypm <- Map(Yield_fn, MSE = MSE[OM_order], OM_names = OM_names[OM_order], Cbias = c(2.25, 1.25, 1)) %>% do.call(rbind, .) %>% 
#  left_join(pm2, by = c("MP", "OM")) %>% 
#  mutate(OM = factor(OM, levels = OM_names[OM_order]), col = color_fn(MP), 
#         label = ifelse(PNOF > 0.99, ">99", ifelse(PNOF < 0.01, "<1", round(100 * PNOF, 0))) %>% paste0(MP, " (", ., "%)"),
#         PB50 = 100 * PB50)

#ypm <- Map(Yield_fn, MSE = MSE[OM_order], OM_names = OM_names[OM_order]) %>% do.call(rbind, .) %>% 
#  left_join(pm2, by = c("MP", "OM")) %>% 
#  mutate(OM = factor(OM, levels = OM_names[OM_order]), col = color_fn(MP), 
#         label = ifelse(PNOF > 0.99, ">99", ifelse(PNOF < 0.01, "<1", round(100 * PNOF, 0))) %>% paste0(MP, " (", ., "%)"),
#         PB50 = 100 * PB50)

#ggplot(ypm, aes(PB50, MeanC, label = label, colour = factor(col), shape = factor(col))) + facet_grid(OM~., scales = "free_y") + 
#  geom_point(size = 2) + ggrepel::geom_text_repel(size = 2.5) + 
#  scale_shape_manual(values = c(8, 1, 16)) + coord_cartesian(ylim = c(0, 1e3)) + xlab("PB50 (%)") + ylab("Observed short-term catch") + 
#  gfplot::theme_pbs() + no_panel_gap + legend_bottom + no_legend
#ggsave("report/GoM_cod/pm_tradeoff.png", width = 3, height = 5)



#ypm <- Map(Yield_fn, MSE = MSE[OM_order], OM_names = OM_names[OM_order], Cbias = c(2.25, 1.25, 1)) %>% do.call(rbind, .) %>% 
#  left_join(pm, by = c("MP", "OM")) %>% mutate(col = color_fn(MP))
#ggplot(ypm %>% filter(PM == "PNOF" | PM == "PB50"), aes(value, MeanC, label = MP, colour = factor(col), shape = factor(col))) + 
#  facet_grid(OM ~ PM) + geom_point() + 
#  ggrepel::geom_text_repel(size = 2.5) + 
#  scale_shape_manual(values = c(8, 1, 16)) + coord_cartesian(ylim = c(0, 1e3)) + xlab("Performance metric") + 
#  ylab("Observed short-term catch") + 
#  gfplot::theme_pbs() + no_panel_gap + legend_bottom + no_legend
#ggsave("report/GoM_cod/pm_tradeoff.png", width = 5, height = 5)



########  Indicators
s_CAA_hist <- get_cod_M02()$data$s_CAA
indicators_raw <- lapply(1:length(MSE), get_indicators, ind_interval = 6, MSE = MSE[OM_order], 
                         MPs = c("M02", "M02_ra", "MRAMP", "MRAMP_ra", "ma", "75%FMSY"), Cbias = c(2.25, 1.25, 1), s_CAA_hist = s_CAA_hist,
                         OM_names = OM_names[OM_order],
                         mat_age = c(0.09, 0.32, 0.70, 0.92, 0.98, 1.00, 1.00, 1.00, 1.00))
indicators <- do.call(rbind, lapply(indicators_raw, getElement, 1)) #%>% mutate(value = ifelse(grepl("mu", Ind), exp(value), value))

indicators_cast <- indicators %>% reshape2::dcast(Year + MP + OM + Sim ~ Ind, value.var = "value") 

######### LDA
rr <- lda(OM ~ Ind_1_mu * Year + SSB_rho * Year + PMat_1_mu * Year + MAge_1_mu * Year, 
          data = filter(indicators_cast, Year > 2018 & MP == "M02_ra") %>% mutate(Year = factor(Year)))

# Index vs. rho
rr <- lda(OM ~ Ind_1_mu * Year + SSB_rho * Year, 
          data = filter(indicators_cast, Year > 2018 & MP == "M02_ra") %>% mutate(Year = factor(Year)))
rr_grid <- expand.grid(Year = table(indicators_cast$Year)[-1] %>% names(),
                       Ind_1_mu = seq(-2, 10, 0.25), SSB_rho = seq(-1, 3, 0.25)) 
rr_pred <- predict(rr, newdata = rr_grid) %>% getElement("posterior") %>% as.data.frame() %>% cbind(rr_grid)

ggplot(rr_pred, aes(Ind_1_mu, SSB_rho)) + facet_wrap(~ Year) + geom_raster(aes(fill = MC)) + 
  gfplot::theme_pbs() + no_panel_gap
ggplot(indicators_cast %>% filter(MP == "M02_ra" & Year > 2018), aes(Ind_1_mu, SSB_rho, colour = OM)) + facet_wrap(~ Year) + geom_point() +
  coord_cartesian(xlim = c(-2, 10), ylim = c(-1, 3))


# Mean age
rr <- lda(OM ~ MAge_1_mu * Year + SSB_rho * Year,
          data = filter(indicators_cast, Year > 2018 & MP == "M02") %>% mutate(Year = factor(Year)))
rr_grid <- expand.grid(Year = table(indicators_cast$Year)[-1] %>% names(),
                       MAge_1_mu = seq(1.2, 1.8, 0.025), SSB_rho = seq(-1, 3, 0.25)) 
rr_pred <- predict(rr, newdata = rr_grid) %>% getElement("posterior") %>% as.data.frame() %>% cbind(rr_grid)

ggplot(rr_pred, aes(MAge_1_mu, SSB_rho)) + facet_wrap(~ Year) + geom_raster(aes(fill = MC)) + 
  gfplot::theme_pbs() + no_panel_gap
ggplot(indicators_cast %>% filter(MP == "M02" & Year > 2018), aes(MAge_1_mu, SSB_rho, colour = OM)) + facet_wrap(~ Year) + geom_point() +
  coord_cartesian(ylim = c(-1, 3))


# Prop. mature age
rr <- lda(OM ~ PMat_1_mu * Year + SSB_rho * Year, 
          data = filter(indicators_cast, Year > 2018 & MP == "M02") %>% mutate(Year = factor(Year)))
rr_grid <- expand.grid(Year = table(indicators_cast$Year)[-1] %>% names(),
                       PMat_1_mu = seq(-0.5, -0.1, 0.01), SSB_rho = seq(-1, 3, 0.25)) 
rr_pred <- predict(rr, newdata = rr_grid) %>% getElement("posterior") %>% as.data.frame() %>% cbind(rr_grid)

#ggplot(rr_pred, aes(PMat_1_mu, SSB_rho)) + facet_wrap(~ Year) + geom_raster(aes(fill = MC)) + 
#  gfplot::theme_pbs() + no_panel_gap
ggplot(indicators_cast %>% filter(MP == "M02" & Year > 2018), aes(PMat_1_mu, SSB_rho)) + facet_wrap(~ Year) + 
  geom_point(aes(colour = OM), alpha = 0.8) +
  geom_contour(data = rr_pred, aes(z = MC), colour = "grey70") + 
  metR::geom_text_contour(data = rr_pred, skip = 1, rotate = TRUE, check_overlap = TRUE, colour = "black", aes(z = MC)) +
  coord_cartesian(xlim = c(-0.5, -0.1), ylim = c(-1, 3)) + labs(x = "Prop. mature", y = "rho") +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom
ggsave("report/GoM_cod/indicators_LDA.png", width = 7, height = 7)


#### All MPs for rho vs. Prop mature for 3 years
MPs <- c("M02", "M02_ra", "MRAMP", "MRAMP_ra", "ma")
Yind <- c(2024, 2030, 2036)
rr_pred <- list()
for(i in 1:length(MPs)) {
  rr <- MASS::lda(OM ~ PMat_1_mu * Year + SSB_rho * Year, 
                  data = filter(indicators_cast, Year > 2018 & MP == MPs[i]) %>% mutate(Year = factor(Year)))
  rr_grid <- expand.grid(Year = table(indicators_cast$Year)[-1] %>% names(),
                         PMat_1_mu = seq(-0.5, -0.1, 0.01), SSB_rho = seq(-1, 3, 0.25)) 
  
  rr_pred[[i]] <- predict(rr, newdata = rr_grid) %>% getElement("posterior") %>% as.data.frame() %>% cbind(rr_grid) %>%
    mutate(MP = MPs[i])
}
rr_pred2 <- do.call(rbind, rr_pred) %>% mutate(MP = factor(MP, levels = MPs)) %>% 
  filter(match(Year, Yind, nomatch = 0) %>% as.logical())

ggplot(indicators_cast %>% filter(MP != "75%FMSY" & match(Year, Yind, nomatch = 0) %>% as.logical()) %>%
  mutate(MP = factor(MP, levels = MPs)), aes(PMat_1_mu, SSB_rho)) + 
  facet_grid(Year ~ MP) + geom_hline(yintercept = 0, linetype = 3) + geom_point(aes(colour = OM), alpha = 0.6) +
  geom_contour(data = rr_pred2, aes(z = MC), colour = "grey70") + 
  metR::geom_text_contour(data = rr_pred2, skip = 1, rotate = TRUE, check_overlap = TRUE, colour = "black", aes(z = MC)) +
  coord_cartesian(xlim = c(-0.5, -0), ylim = c(-1, 3)) + labs(x = "Prop. mature (log)", y = expression(rho[SSB])) +
  scale_x_continuous(breaks = c(-0.4, -0.2, 0)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom
ggsave("report/GoM_cod/indicators_LDA2.png", width = 6, height = 4)




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
#  facet_wrap(~ factor(MP, levels = c("M02", "M02_ra", "MRAMP", "MRAMP_ra", "ma")), scales = "free") + 
#  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = OM, colour = OM), alpha = 0.1) + 
#  ylab("Mahalanobis distance") +
#  geom_point(data = mah_gr %>% filter(pval <= 0.05), y = 0, shape = 3) + 
#  scale_x_continuous(breaks = c(2020, 2040, 2060)) +
#  gfplot::theme_pbs() + legend_bottom
#ggsave("report/GoM_cod/mahalanobis.png", height = 4.5, width = 7)




# Time series plots
indicators_ts <- summarise(group_by(indicators, Ind, Year, MP, OM),
                           y = median(value, na.rm = TRUE),
                           ymin = quantile(value, 0.25, na.rm = TRUE), ymax = quantile(value, 0.75, na.rm = TRUE))

# rho, Cat, Ind
ggplot(indicators_ts, aes(x = Year, shape = OM, colour = OM)) + facet_grid(Ind ~ MP, scales = "free_y") + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = OM), alpha = 0.1) + 
  geom_point(aes(y = y)) + geom_line(aes(y = y)) + xlab("Year") + ylab("Indicator Value") +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom


# Mean Age
ggplot(filter(indicators_ts, grepl("PMat", Ind)), aes(x = Year, shape = OM, colour = OM)) + facet_grid(Ind ~ MP, scales = "free_y") + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = OM), alpha = 0.1) + 
  geom_point(aes(y = y)) + geom_line(aes(y = y)) + 
  xlab("Year") + ylab("Indicator Value") +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom


indicators_plot <- filter(indicators_ts, match(Ind, c("SSB_rho", "Cat_slp", "Cat_mu", "Ind_1_slp", "MAge_1_slp"), nomatch = 0) %>% as.logical()) %>%
  mutate(MP = factor(MP, levels = c("M02", "M02_ra", "MRAMP", "MRAMP_ra", "ma", "75%FMSY"))) %>%
  mutate(Ind = factor(Ind, levels = c("SSB_rho", "Cat_slp", "Cat_mu", "Ind_1_slp", "MAge_1_slp"))) #%>% filter(Year >= 2024)
ggplot(indicators_plot, aes(x = Year, shape = OM)) + 
  geom_hline(data = data.frame(Ind = factor(c("SSB_rho", "Cat_slp", "Ind_1_slp", "MAge_1_slp"), 
                                            levels = c("SSB_rho", "Cat_slp", "Cat_mu", "Ind_1_slp", "MAge_1_slp")), yy = 0), 
             aes(yintercept = yy), linetype = 2) + 
  geom_point(aes(colour = OM, y = y)) + geom_line(aes(colour = OM, y = y)) + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = OM), alpha = 0.3) +
  facet_grid(Ind ~ MP, scales = "free_y") +  
  xlab("Year") + ylab("Indicator Value") + scale_x_continuous(breaks = c(2020, 2040, 2060)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom


Ind_order <- c("SSB_rho", "Cat_mu", "Ind_1_mu", "MAge_1_mu", "PMat_1_mu")
Ind_order2 <- c("rho", "Catch", "Index", "Mean Age", "Prop. Mature")
indicators_plot <- filter(indicators_ts, match(Ind, Ind_order, nomatch = 0) %>% as.logical()) %>%
  mutate(MP = factor(MP, levels = c("M02", "M02_ra", "MRAMP", "MRAMP_ra", "ma", "75%FMSY"))) %>%
  left_join(data.frame(Ind = Ind_order, Ind2 = factor(Ind_order2, levels = Ind_order2)))

ggplot(indicators_plot, aes(x = Year, shape = OM)) + 
  geom_hline(data = data.frame(Ind2 = factor("rho", levels = Ind_order2), yy = 0), aes(yintercept = yy), linetype = 2) + 
  geom_point(aes(colour = OM, y = y)) + geom_line(aes(colour = OM, y = y)) + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = OM), alpha = 0.3) +
  facet_grid(Ind2 ~ MP, scales = "free_y") +  
  xlab("Year") + ylab("Predicted indicator value") + scale_x_continuous(breaks = c(2020, 2040, 2060)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("report/GoM_cod/indicator_ts_log.png", height = 7, width = 7)


##### TS of indicators relative to FMSY75
relative_ind <- filter(indicators, MP == "75%FMSY") %>% mutate(value2 = value, value = NULL, MP = NULL) %>% 
  right_join(filter(indicators, MP != "75%FMSY")) %>% filter(grepl("1", Ind) | grepl("Cat", Ind) | grepl("SSB", Ind)) %>%
  mutate(diff = value - value2, value = NULL, value2 = NULL) %>% filter(!is.na(diff)) %>% group_by(Ind, Year, MP, OM) %>% 
  summarise(y = median(diff, na.rm = TRUE), ymin = quantile(diff, 0.25, na.rm = TRUE), ymax = quantile(diff, 0.75, na.rm = TRUE)) %>% 
  filter(Year > 2020) %>% mutate(MP = factor(MP, levels = c("M02", "M02_ra", "MRAMP", "MRAMP_ra", "ma")))
ggplot(relative_ind, aes(x = Year, shape = OM, colour = OM)) + 
  geom_hline(yintercept = 0, linetype = 2) + 
  geom_point(aes(y = y)) + geom_line(aes(y = y)) + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = OM), alpha = 0.1) +
  facet_grid(Ind ~ MP, scales = "free_y") +  
  xlab("Year") + ylab(latex2exp::TeX("$\\Delta$ Indicator")) + scale_x_continuous(breaks = c(2020, 2040, 2060)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom
ggsave("report/GoM_cod/indicator_delta.png", height = 7, width = 7)





# Scatter plots M02 only
indicators_scatter <- reshape2::dcast(indicators, formula = Year + Sim + MP + OM ~ Ind, value.var = "value") 

F_FMSY <- lapply(1:length(MSE), function(x, MPs) MSE[[x]]@F_FMSY[, match(MPs, MSE[[x]]@MPs), , drop = FALSE] %>% 
                   structure(dimnames = list(1:MSE[[x]]@nsim, MPs, MSE[[x]]@nyears + 1:MSE[[x]]@proyears)) %>% 
                   reshape2::melt(varnames = c("Sim", "MP", "time"), value.name = "F.FMSY") %>% 
                   mutate(Year = time + MSE[[x]]@OM$CurrentYr[1] - MSE[[x]]@nyears, time = NULL) %>%
                   mutate(OM = paste0("NR", x)), MPs = c("M02", "M02_ra", "MRAMP", "MRAMP_ra", "ma"))
F_FMSY <- do.call(rbind, F_FMSY)


indicators_regress <- left_join(indicators_scatter, F_FMSY) %>% filter(!is.na(F.FMSY))

ggplot(filter(indicators_regress, MP == "M02"), aes(x = SSB_rho, y = F.FMSY)) + 
  facet_wrap(~ Year) + geom_point(alpha = 0.9, aes(colour = OM, shape = OM)) + 
  geom_hline(yintercept = 1, linetype = 2) + geom_vline(xintercept = 0, linetype = 2) + labs(x = expression(rho[SSB]), y = expression(F/F[MSY])) +
  coord_cartesian(xlim = c(-1, 5), ylim = c(0, 15)) +
  #geom_smooth(data = indicators_scatter, method = lm, se = FALSE, colour = "black", formula = y ~ x) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom
ggsave("report/GoM_cod/indicator_scatter.png", height = 6, width = 6)

out <- indicators_regress %>% group_by(Year, MP) %>% summarise(slp = lm(F.FMSY ~ SSB_rho) %>% coef() %>% getElement(2))
ggplot(out, aes(Year, slp, shape = MP)) + geom_line() + geom_point() +
  geom_hline(yintercept = 0, linetype = 2) + labs(y = "Regression slope") + 
  gfplot::theme_pbs() + no_panel_gap + legend_bottom
ggsave("report/GoM_cod/indicator_slope.png", height = 3, width = 4)

#### Plot recruitment
SRA <- lapply(1:3, function(i) readRDS(paste0("GoM_cod/SRA_NR", i, ".rds")))

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

rec <- Map(plot_recruitment, MSE = MSE[OM_order], SRA = SRA[OM_order], OM_names = OM_names[OM_order],
           MoreArgs = list(MPs = c("M02", "M02_ra", "MRAMP", "MRAMP_ra", "ma", "75%FMSY")))

recpred <- lapply(rec, getElement, "Rpred") %>% do.call(rbind, .) %>% group_by(MP, OM, Year) %>% summarise(R = mean(R)) %>%
  mutate(OM = factor(OM, levels = OM_names[OM_order]))
rhist <- lapply(rec, getElement, "Rhist") %>% do.call(rbind, .) %>% mutate(OM = factor(OM, levels = OM_names[OM_order])) %>%
  filter(Year >= 2009)

ggplot(recpred, aes(Year, R)) + geom_line(aes(colour = MP)) + geom_point(aes(colour = MP, shape = MP)) + 
  facet_grid(OM ~ .) + geom_line(data = rhist) +
  geom_vline(xintercept = 2018, linetype = 2) + labs(y = "Recruitment") + 
  gfplot::theme_pbs() + no_panel_gap + legend_bottom
ggsave("report/GoM_cod/mean_recruitment.png", width = 3, height = 4.5)  
