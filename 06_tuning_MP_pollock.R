
library(MSEtool)
library(mfSCA)
library(dplyr)
library(ggplot2)


######## Reference points
ref_pt <- readRDS(file = "pollock/pollock_ref_pt.rds")
ref_pt_plot <- do.call(rbind, lapply(ref_pt, function(i) i[1, ])) %>% as.data.frame()
ref_pt_plot$OM <- paste0("NR", 1:3)

######## Load MSE
MSE <- lapply(1:3, function(x) {
  y <- readRDS(paste0("pollock/MSE_pollock_tuningMP_NR", x, ".rds"))
  y@OM$SSBMSY <- ref_pt_plot$SSBMSY[x]
  y@B_BMSY <- y@SSB/y@OM$SSBMSY
  
  y@OM$FMSY <- ref_pt_plot$FMSY[x]
  y@F_FMSY <- y@FM/y@OM$FMSY 
  
  y@OM$RefY <- ref_pt_plot$OY[x]
  
  return(y)
})

######## FM df
FM <- lapply(1:3, function(x) MSE[[x]]@FM[, 1, ] %>% 
               structure(dimnames = list(Sim = 1:MSE[[x]]@nsim, Year = MSE[[x]]@OM$CurrentYr[1] + 1:MSE[[x]]@proyears)) %>%
               reshape2::melt(value.name = "FM") %>% mutate(F_FMSY = FM/MSE[[x]]@OM$FMSY[1], OM = paste0("NR", x)))
FM <- do.call(rbind, FM)

######## Generate indicators
s_CAA_hist <- get_pollock_base()$data$s_CAA
indicators_raw <- lapply(1:length(MSE), get_indicators, ind_interval = 6, MSE = MSE, 
                         MPs = "tuning_MP",  s_CAA_hist = s_CAA_hist,
                         mah_ind = c("MAge_1_slp", "MAge_1_mu", "Cat_slp", "Cat_mu", "Ind_1_slp", "MAge_1_slp"),
                         Year_vec = seq(MSE[[1]]@nyears, MSE[[1]]@nyears + MSE[[1]]@proyears, 3))
indicators <- do.call(rbind, lapply(indicators_raw, getElement, 1)) %>% mutate(MP = NULL) %>% left_join(FM) %>% 
  filter(!is.na(value) & Year > 2020 & (grepl("Cat", Ind) | grepl("1", Ind)))
indicators_cast <- reshape2::dcast(indicators, Sim + OM + Year + FM + F_FMSY ~ Ind, value.var = "value") %>% mutate(Year = factor(Year))


######## Plot OM
no_legend <- theme(legend.position = "none")
no_panel_gap <- theme(panel.spacing = unit(0, "in"))
legend_bottom <- theme(legend.position = "bottom")

ggplot(indicators, aes(value, F_FMSY)) + geom_point(aes(colour = OM, shape = OM)) + facet_grid(Year ~ Ind, scales = "free") +
  geom_hline(yintercept = 1, linetype = 2) + 
  geom_smooth(data = indicators %>% mutate(OM = NULL), aes(y = ifelse(F_FMSY > 1, 1, 0)), colour = "black",
              method = "glm", method.args = list(family = "binomial")) + 
  gfplot::theme_pbs() + no_panel_gap + legend_bottom


######### Plot mah
mah <- do.call(rbind, lapply(indicators_raw, getElement, 2)) %>% left_join(FM)

ggplot(mah, aes(value, F_FMSY)) + facet_wrap(~ Year, scales = "free_x") + geom_point(aes(colour = OM, shape = OM)) +
  coord_cartesian(ylim = c(0, 2)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom

ggplot(mah, aes(value)) + facet_wrap(~ Year, scales = "free") + geom_density(aes(fill = OM, colour = OM)) +
  scale_colour_grey() + scale_fill_grey() + 
  gfplot::theme_pbs() + no_panel_gap + legend_bottom








rr <- glm(ifelse(F_FMSY > 1, 1, 0) ~ Cat_mu + Year, data = indicators_cast, family = "binomial")
rr2 <- glm(ifelse(F_FMSY > 1, 1, 0) ~ Cat_mu * Year + Year * Ind_1_mu, data = indicators_cast, family = "binomial")
rr3 <- glm(ifelse(F_FMSY > 1, 1, 0) ~ Ind_1_mu * Year, data = indicators_cast, family = "binomial")
rr4 <- glm(ifelse(F_FMSY > 1, 1, 0) ~ Cat_mu * Year, data = indicators_cast, family = "binomial")
AIC(rr, rr2, rr3, rr4)
