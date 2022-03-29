library(MSEtool)
library(dplyr)
library(mfSCA)
library(ggplot2)

######## OM order
OM_order <- c(1, 3, 2)
OM_names <- c("MC", "MCIM", "IM")

######## MP names
MPs <- c("M02_nra", "M02_ra", "MRAMP_nra", "MRAMP_ra", "MA", "75%FMSY")

######## Reference points
ref_pt <- readRDS("GoM_cod/cod_ref_pt.rds")

ref_pt_plot <- lapply(ref_pt, function(i) lapply(i, function(ii) ii[1, ]) %>% do.call(rbind, .)) %>% do.call(rbind, .) %>%
  as.data.frame()
ref_pt_plot$OM <- rep(OM_names, 2)
ref_pt_plot$type <- rep(c("old", "new"), each = 3)

######## Load MSE and update reference points
MSE <- lapply(1:3, function(x) {
  y <- readRDS(paste0("GoM_cod/MSE_cod_", OM_names[x], ".rds"))
  
  y@OM$SSBMSY <- filter(ref_pt_plot, type == 'new')$SSBMSY[x]
  y@B_BMSY <- y@SSB/y@OM$SSBMSY
  
  y@OM$FMSY <- filter(ref_pt_plot, type == 'new')$FMSY[x]
  y@F_FMSY <- y@FM/y@OM$FMSY 
  
  y@OM$RefY <- filter(ref_pt_plot, type == 'new')$OY[x]
  
  return(y)
})

######## Plot OM trajectories (SSB, F, catch)
resOM <- plot_OM(MSE[OM_order], MPs = MPs, 0.10, 0.90, yfilter = 10, OM_names = OM_names[OM_order])

no_legend <- theme(legend.position = "none")
no_panel_gap <- theme(panel.spacing = unit(0, "in"))
legend_bottom <- theme(legend.position = "bottom")

OM_letters <- expand.grid(OM = OM_names[OM_order], MP = MPs)
OM_letters$label <- paste0("(", letters[1:nrow(OM_letters)], ")")

# F
ggplot(resOM[[1]], aes(Year, y = med, ymin = low, ymax = high)) + facet_grid(factor(OM) ~ MP, scales = "free_y") + 
  geom_hline(yintercept = 0, col = "white") + 
  geom_vline(xintercept = MSE[[1]]@OM$CurrentYr[1], linetype = 3) +
  geom_ribbon(fill = "grey80") + geom_line(size = 0.9) + 
  geom_label(data = OM_letters, inherit.aes = FALSE, x = -Inf, y = Inf, colour = NA, hjust = -0.1, vjust = 1, 
             aes(label = label)) +
  geom_text(data = OM_letters, inherit.aes = FALSE, x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, 
            aes(label = label)) +
  geom_hline(data = ref_pt_plot[4:6, ], aes(yintercept = FMSY), linetype = 2, size = 0.5) +
  ylab("Fishing mortality") + scale_y_continuous(expand = c(0.01, 0.01)) + 
  gfplot::theme_pbs() + no_panel_gap
ggsave("report/GoM_cod/MSE_OM_F.png", height = 3.5, width = 8.5)

# SSB
resOM_label <- plot_OM(MSE[OM_order], MPs = MPs, 0.10, 0.90, yfilter = 10, OM_names = OM_names[OM_order],
                       aggregate_across_years = TRUE, rel = TRUE)
SSB_label <- resOM_label[[2]] %>% 
  mutate(label = round(med, 2) %>% format())

g <- ggplot(resOM[[2]], aes(Year, y = med, ymin = low, ymax = high)) + facet_grid(factor(OM) ~ MP, scales = "free_y") +  
  geom_hline(yintercept = 0, col = "white") + 
  geom_vline(xintercept = MSE[[1]]@OM$CurrentYr[1], linetype = 3) +
  geom_ribbon(fill = "grey80") + geom_line(size = 0.9) + 
  geom_text(data = SSB_label, inherit.aes = FALSE, x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,
            aes(label = label)) +
  geom_label(data = OM_letters, inherit.aes = FALSE, x = -Inf, y = Inf, colour = NA, hjust = -0.1, vjust = 1, 
             aes(label = label)) +
  geom_text(data = OM_letters, inherit.aes = FALSE, x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, 
            aes(label = label)) +
  geom_hline(data = ref_pt_plot[4:6, ], aes(yintercept = SSBMSY), linetype = 2) +
  geom_hline(data = ref_pt_plot[4:6, ], aes(yintercept = 0.5 * SSBMSY), linetype = 2) +
  ylab("SSB") + scale_y_continuous(labels = scales::comma, expand = c(0.01, 0.01)) +
  gfplot::theme_pbs() + no_panel_gap + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
ggsave("report/GoM_cod/MSE_OM_SSB.png", height = 3.5, width = 8.5)
ggsave("report/GoM_cod/MSE_OM_SSB.png", g, height = 80, width = 170, dpi = 600, units = "mm")

# Catch
ggplot(resOM[[3]], aes(Year, y = med, ymin = low, ymax = high)) + facet_grid(factor(OM) ~ MP, scales = "free_y") + 
  geom_hline(yintercept = 0, col = "white") + 
  geom_vline(xintercept = MSE[[1]]@OM$CurrentYr[1], linetype = 3) +
  geom_ribbon(fill = "grey80") + geom_line(size = 0.9) + 
  geom_label(data = OM_letters, inherit.aes = FALSE, x = -Inf, y = Inf, colour = NA, hjust = -0.1, vjust = 1, 
             aes(label = label)) +
  geom_text(data = OM_letters, inherit.aes = FALSE, x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, 
            aes(label = label)) +
  geom_hline(data = ref_pt_plot[4:6, ], aes(yintercept = OY), linetype = 3) +
  ylab("Catch") + scale_y_continuous(labels = scales::comma, expand = c(0.01, 0.01)) +
  gfplot::theme_pbs() + no_panel_gap
ggsave("report/GoM_cod/MSE_OM_C.png", height = 3.5, width = 8.5)


####### Plot SSB rho
rho <- local({
  
  res <- list()
  res[[1]] <- plot_EM(MSE[OM_order], MPs = MPs[1:4], OM_names = OM_names[OM_order]) %>% 
    lapply(function(x) {
      x$EM <- vapply(x$MP, function(y) strsplit(y, "_")[[1]][1], character(1))
      return(x)
    })
  
  res[[2]] <- plot_EM(MSE[OM_order], MPs = "MA", OM_names = OM_names[OM_order]) %>% 
    lapply(function(x) cbind(x, EM = "M02"))
  res[[3]] <- plot_EM(MSE[OM_order], MPs = "MA", model_index = 2, OM_names = OM_names[OM_order]) %>% 
    lapply(function(x) cbind(x, EM = "MRAMP"))
  
  res[[4]] <- plot_EM(MSE[OM_order], MPs = "75%FMSY", model_index = 1, OM_names = OM_names[OM_order])
  res[[4]][[3]] <- cbind(res[[4]][[3]], EM = "M02")
  
  res[[5]] <- plot_EM(MSE[OM_order], MPs = "75%FMSY", model_index = 2, OM_names = OM_names[OM_order])
  res[[5]][[3]] <- cbind(res[[5]][[3]], EM = "MRAMP")
  
  lapply(c("F_FMSY", "B_BMSY", "rho"), function(x) {
    lapply(res, getElement, x) %>% do.call(rbind, .) %>% mutate(MP = factor(MP, levels = MPs))
  })
})


ggplot(rho[[3]] %>% filter(Year > 2020), aes(x = as.numeric(Year), y = `50%`, ymin = `25%`, ymax = `75%`, shape = EM)) + 
  facet_grid(factor(OM) ~ MP, scales = "free_y") + 
  geom_hline(yintercept = 0, linetype = 3) + 
  geom_point(size = 1.25) + geom_linerange(size = 0.5) + 
  xlab("Year") + ylab(expression(rho[SSB])) + theme(legend.position = "bottom") +
  geom_text(data = OM_letters, inherit.aes = FALSE, x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, 
            aes(label = label)) +
  scale_x_continuous(breaks = c(2020, 2040, 2060)) + scale_shape_manual(values = c(16, 1)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(shape = "")
ggsave("report/GoM_cod/MSE_EM_rho.png", height = 4, width = 6.5)



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
  mutate(MP = factor(MP, levels = MPs), OM = factor(OM, levels = OM_names[OM_order]),
         label = ifelse(value > 0.99, ">99", ifelse(value < 0.01, "<1", round(100 * value, 0))) %>% paste0("%"))

# Trade - off
pm2 <- reshape2::dcast(pm, MP + OM ~ PM)

# Calculate short term yield
Yield_fn <- function(MSE, Yrs = c(1, 10), MPs, OM_names = "OM", Cbias = 1) {
  geo_mean <- function(x) (sum(log(x))/length(x)) %>% exp()
  CC <- (MSE@C[, match(MPs, MSE@MPs), Yrs[1]:Yrs[2], drop = FALSE]/Cbias) %>% apply(2, geo_mean)
  data.frame(MP = MPs, OM = OM_names, MeanC = CC)
}

AM_fn <- function(x) {
  y <- sapply(as.character(x), function(xx) strsplit(xx, "_")[[1]][1])
  y[grep("MA", y)] <- "Both"
  y[grep("75%FMSY", y)] <- "None"
  factor(y, levels = c("M02", "MRAMP", "Both", "None"))
}

ypm <- Map(Yield_fn, MSE = MSE[OM_order], OM_names = OM_names[OM_order], Cbias = c(2.25, 1.25, 1),
           MoreArgs = list(MPs = MPs)) %>% do.call(rbind, .) %>% 
  left_join(pm2, by = c("MP", "OM")) %>% 
  mutate(OM = factor(OM, levels = OM_names[OM_order]), AM = AM_fn(MP), 
         label = ifelse(PB50 > 0.99, ">99", ifelse(PB50 < 0.01, "<1", round(100 * PB50, 0))) %>% paste0(MP, " (", ., "%)"),
         PNOF = 100 * PNOF)

g <- ggplot(ypm, aes(PNOF, MeanC, label = label, colour = AM)) + facet_grid(OM~., scales = "free_y") + 
  geom_point(aes(shape = AM), size = 2) + 
  ggrepel::geom_text_repel(aes(label = label), size = 2.5) + 
  geom_text(data = data.frame(OM = OM_names[OM_order], txt = c("(a)", "(b)", "(c)")),
            aes(label = txt), x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,
            inherit.aes = FALSE) +
  scale_shape_manual(values = c(8, 1, 16, 2)) + 
  scale_colour_manual(values = c("blue", "red", "#13F240FF", "black")) + 
  coord_cartesian(ylim = c(0, 1.1e3)) + xlab("PNOF (%)") + ylab("Observed short-term catch") +
  scale_y_continuous(breaks = seq(0, 1000, 250)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom +
  guides(label = "none")
ggsave("report/GoM_cod/pm_tradeoff.png", width = 3, height = 5.5)
ggsave("report/GoM_cod/pm_tradeoff.png", g, width = 85, height = 140, dpi = 600, units = "mm")









########  Indicators
s_CAA_hist <- get_cod_M02()$data$s_CAA
indicators_raw <- lapply(1:length(MSE), get_indicators, ind_interval = 3, MSE = MSE[OM_order], 
                         MPs = MPs, Cbias = c(2.25, 1.25, 1), s_CAA_hist = s_CAA_hist,
                         OM_names = OM_names[OM_order],
                         mat_age = c(0.09, 0.32, 0.70, 0.92, 0.98, 1.00, 1.00, 1.00, 1.00))
indicators <- do.call(rbind, lapply(indicators_raw, getElement, 1))
indicators_cast <- indicators %>% reshape2::dcast(Year + MP + OM + Sim ~ Ind, value.var = "value") 

######### LDA
#### All MPs for rho vs. Prop mature
Yind <- c(2024, 2030, 2036)
rr_pred <- rr <- rr_CV <- list()
data_list <- lapply(1:6, function(i) dplyr::filter(indicators_cast, Year > 2018 & MP == MPs[i]) %>% mutate(Year = factor(Year)))

# For each MP, do the LDA on the training set (rr), perform the cross_validation (rr_CV), and calculate the posterior surface (rr_pred)
for(i in 1:length(MPs)) {
  rr[[i]] <- MASS::lda(OM ~ PMat_1_mu * Year + SSB_rho * Year, data = data_list[[i]], method = "mle")
  
  rr_CV[[i]] <- MASS::lda(OM ~ PMat_1_mu * Year + SSB_rho * Year, data = data_list[[i]], method = "mle", CV = TRUE)
  
  rr_grid <- expand.grid(Year = table(indicators_cast$Year)[-1] %>% names(),
                         PMat_1_mu = seq(-0.5, -0.1, 0.01), SSB_rho = seq(-1, 3, 0.25)) 
  
  rr_pred[[i]] <- predict(rr[[i]], newdata = rr_grid) %>% getElement("posterior") %>% as.data.frame() %>% cbind(rr_grid) %>%
    mutate(MP = MPs[i])
}

# Subset for Years in Yind
rr_pred2 <- do.call(rbind, rr_pred) %>% mutate(MP = factor(MP, levels = MPs)) %>% 
  filter(match(Year, Yind, nomatch = 0) %>% as.logical())
  
# Correct classification rate (cross validation)
class <- Map(function(x, y, MP) {
  yy <- data.frame(misclass = abs(as.numeric(x$class) - as.numeric(y$OM)) %>% as.logical(), Year = y$Year)
  summarise(group_by(yy, Year), class_correct = round(1 - sum(misclass)/length(misclass), 2) %>% format()) %>% mutate(MP = MP)
}, x = rr_CV, y = data_list, MP = MPs) %>% do.call(rbind, .) %>% filter(match(Year, Yind, nomatch = 0) %>% as.logical()) %>%
  mutate(MP = factor(MP, MPs))


LDA_letters <- expand.grid(Year = Yind, MP = MPs)
LDA_letters$label <- paste0("(", letters[1:nrow(LDA_letters)], ")")

g <- ggplot(indicators_cast %>% filter(match(Year, Yind, nomatch = 0) %>% as.logical()) %>%
  mutate(MP = factor(MP, levels = MPs)), aes(PMat_1_mu, SSB_rho)) + 
  facet_grid(Year ~ MP) + geom_hline(yintercept = 0, linetype = 3) + geom_point(aes(colour = OM, shape = OM), alpha = 0.6) +
  geom_contour(data = rr_pred2, breaks = 0.5, colour = "black", aes(z = MC)) + 
  geom_text(data = class, aes(label = class_correct), x = Inf, y = Inf, vjust = 1.1, hjust = 1.1) +
  geom_text(data = LDA_letters, inherit.aes = FALSE, x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, 
            aes(label = label)) +
  coord_cartesian(xlim = c(-0.5, -0), ylim = c(-0.5, 3.5)) + labs(x = expression(P[mat]), y = expression(rho[SSB])) +
  scale_x_continuous(breaks = c(-0.4, -0.2, 0)) +
  scale_shape_manual(values = c("MC" = 16, "IM" = 1, "MCIM" = 4)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom
ggsave("report/GoM_cod/indicators_LDA2_3y.png", width = 6.5, height = 4)
ggsave("report/GoM_cod/indicators_LDA2_3y.png", g, height = 100, width = 170, dpi = 600, units = "mm")



#### Plot mean recruitment (supplemental figure)
SRA <- lapply(1:3, function(i) readRDS(paste0("GoM_cod/SRA_", OM_names[i], ".rds")))

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

rec <- Map(plot_recruitment, MSE = MSE[OM_order], SRA = SRA[OM_order], OM_names = OM_names[OM_order], MoreArgs = list(MPs = MPs))

recpred <- lapply(rec, getElement, "Rpred") %>% do.call(rbind, .) %>% group_by(MP, OM, Year) %>% summarise(R = mean(R)) %>%
  mutate(OM = factor(OM, levels = OM_names[OM_order]))
recpred$OM <- factor(recpred$OM, levels = OM_names[OM_order])
rhist <- lapply(rec, getElement, "Rhist") %>% do.call(rbind, .) %>% mutate(OM = factor(OM, levels = OM_names[OM_order])) %>%
  filter(Year >= 2009)

ggplot(recpred, aes(Year, R)) + geom_line(aes(colour = MP)) + geom_point(aes(colour = MP, shape = MP)) + 
  facet_grid(OM ~ .) + geom_line(data = rhist) +
  geom_text(data = data.frame(OM = OM_names[OM_order], label = c("(a)", "(b)", "(c)")), aes(label = label),
            x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1,
            inherit.aes = FALSE) +
  geom_vline(xintercept = 2018, linetype = 2) + labs(y = "Recruitment") + 
  scale_y_continuous(labels = scales::comma) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom
ggsave("report/GoM_cod/mean_recruitment.png", width = 3, height = 4.5)  
