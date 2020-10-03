
library(MSEtool)
library(mfSCA)
library(dplyr)
library(ggplot2)


######## Reference points
ref_pt <- readRDS(file = "GoM_cod/cod_ref_pt.rds")

ref_pt_plot <- do.call(rbind, lapply(ref_pt, function(i) do.call(rbind, lapply(i, function(ii) ii[1, ])))) %>% 
  as.data.frame()
ref_pt_plot$OM <- rep(paste0("NR", 1:3), 2)
ref_pt_plot$type <- rep(c("old", "new"), each = 3)


######## Load MSE
MSE <- lapply(1:3, function(x) {
  y <- readRDS(paste0("GoM_cod/MSE_cod_tuningMP_NR", x, ".rds"))
  y@OM$SSBMSY <- filter(ref_pt_plot, type == 'new')$SSBMSY[x]
  y@B_BMSY <- y@SSB/y@OM$SSBMSY
  
  y@OM$FMSY <- filter(ref_pt_plot, type == 'new')$FMSY[x]
  y@F_FMSY <- y@FM/y@OM$FMSY 
  
  y@OM$RefY <- filter(ref_pt_plot, type == 'new')$OY[x]
  
  return(y)
})

######## FM df
FM <- lapply(1:3, function(x) MSE[[x]]@FM[, 1, ] %>% 
               structure(dimnames = list(Sim = 1:MSE[[x]]@nsim, Year = MSE[[x]]@OM$CurrentYr[1] + 1:MSE[[x]]@proyears)) %>%
               reshape2::melt(value.name = "FM") %>% mutate(F_FMSY = FM/MSE[[x]]@OM$FMSY[1], OM = paste0("NR", x))) %>% do.call(rbind, .)

######## Generate indicators
s_CAA_hist <- get_cod_M02()$data$s_CAA
indicators_raw <- lapply(1:length(MSE), get_indicators, ind_interval = 6, MSE = MSE, 
                         MPs = c("tuning_MP"), Cbias = c(2.25, 1.25, 1), s_CAA_hist = s_CAA_hist,
                         mah_ind = c("MAge_1_slp", "MAge_1_mu", "Cat_slp", "Cat_mu", "Ind_1_slp", "MAge_1_slp"),
                         Year_vec = seq(MSE[[1]]@nyears, MSE[[1]]@nyears + MSE[[1]]@proyears, 3))
indicators <- lapply(indicators_raw, getElement, 1) %>% do.call(rbind, .) %>% mutate(MP = NULL) %>% left_join(FM) %>% 
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
  coord_cartesian(ylim = c(0, 2)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom

ggplot(indicators %>% filter(grepl("Ind_1", Ind)), aes(value, F_FMSY)) + geom_point(aes(colour = OM, shape = OM)) + facet_grid(Year ~ Ind, scales = "free") +
  geom_hline(yintercept = 1, linetype = 2) + 
  geom_smooth(data = indicators %>% mutate(OM = NULL), aes(y = ifelse(F_FMSY > 1, 1, 0)), colour = "black",
              method = "glm", method.args = list(family = "binomial")) + 
  coord_cartesian(ylim = c(0, 2)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom



######### Plot mah
mah <- lapply(indicators_raw, getElement, 2) %>% do.call(rbind, .) %>% left_join(FM)

ggplot(mah, aes(value, F_FMSY)) + facet_wrap(~ Year) + geom_point(aes(colour = OM, shape = OM)) +
  coord_cartesian(ylim = c(0, 2)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom

ggplot(mah, aes(value)) + facet_wrap(~ Year, scales = "free") + geom_density(aes(fill = OM, colour = OM)) +
  scale_colour_grey() + scale_fill_grey() + 
  gfplot::theme_pbs() + no_panel_gap + legend_bottom





rr <- glm(ifelse(F_FMSY > 1, 1, 0) ~ Cat_mu + Year, data = indicators_cast, family = "binomial")
rr2 <- glm(ifelse(F_FMSY > 1, 1, 0) ~ Cat_mu * Year + Year * Ind_1_mu, data = indicators_cast, family = "binomial")
rr3 <- glm(ifelse(F_FMSY > 1, 1, 0) ~ Ind_1_mu * Year, data = indicators_cast, family = "binomial")
rr4 <- glm(ifelse(F_FMSY > 1, 1, 0) ~ Cat_mu * Year, data = indicators_cast, family = "binomial")
rr5 <- glm(ifelse(F_FMSY > 1, 1, 0) ~ Ind_1_slp * Year, data = indicators_cast, family = "binomial")
rr6 <- glm(ifelse(F_FMSY > 1, 1, 0) ~ Ind_1_slp * Year + Ind_1_mu * Year, data = indicators_cast, family = "binomial")
AIC(rr, rr2, rr3, rr4, rr5, rr6)

predict(rr6, newdata = list(Year = factor(2057), Ind_1_slp = -0.1, Ind_1_mu = 5), type = "response", se.fit = TRUE)
predict(rr5, newdata = list(Year = factor(2057), Ind_1_slp = -1), type = "response", se.fit = TRUE)

plot(jitter(POF_slp_mu) ~ Year, ind_summary, ylim = c(0, 1), typ = 'o')


lines(POF_max ~ Year, ind_summary, typ = "o")
abline(h = 0.5)

out <- indicators_cast %>% mutate(POF = predict(rr5, newdata = list(Year = Year, Ind_1_slp = Ind_1_slp), type = "response"))
ggplot(out, aes(Ind_1_slp, F_FMSY)) + geom_point(aes(colour = OM, shape = OM)) + facet_wrap(~ Year, scales = "free") +
  geom_hline(yintercept = c(0.5, 1), linetype = 2) + 
  geom_line(aes(y = POF)) + coord_cartesian(ylim = c(0, 2)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom



ind_POF <- indicators_cast %>% mutate(POF_slp = predict(rr5, type = "response"), POF_slp_mu = predict(rr6, type = 'response')) %>%
  reshape2::melt(value.name = "POF", id.vars = c("Ind_1_slp", "Ind_1_mu", "Year"), measure.vars = c("POF_slp", "POF_slp_mu")) %>%
  reshape2::melt(id.vars = c("variable", "POF", "Year"), variable.name = "Ind", value.name = "value")

ind_summary <- ind_POF %>% group_by(Year, variable) %>% summarise(POF_min = min(POF), POF_max = max(POF))

plot(POF_max ~ Year, data = filter(ind_summary, variable == "POF_slp"), typ = 'o', ylim = c(0, 1))
lines(POF_min ~ Year, data = filter(ind_summary, variable == "POF_slp"), typ = 'o')
abline(h = 0.5, lty = 3)

plot(POF_max ~ Year, data = filter(ind_summary, variable == "POF_slp_mu"), typ = 'o', ylim = c(0, 1))
lines(POF_min ~ Year, data = filter(ind_summary, variable == "POF_slp_mu"), typ = 'o')
abline(h = 0.5, lty = 3)

ggplot(indicators %>% filter(grepl("Ind", Ind)), aes(value, F_FMSY)) + geom_point(aes(colour = OM, shape = OM)) + facet_grid(Year ~ Ind, scales = "free") +
  geom_hline(yintercept = 1, linetype = 2) + 
  geom_line(data = ind_POF, aes(y = POF, group = variable, linetype = variable)) +
  coord_cartesian(ylim = c(0, 2)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom

ggplot(ind_POF %>% filter(variable == "POF_slp_mu"), aes(x = value, y = POF, group = variable, colour = variable)) + facet_grid(Year ~ Ind, scales = "free") +
  coord_cartesian(ylim = c(0, 1)) + geom_hline(yintercept = 0.5) + geom_line() + 
  gfplot::theme_pbs() + no_panel_gap + legend_bottom




ind_POF <- indicators_cast %>% mutate(POF_slp = predict(rr6, newdata = list(Year = Year, Ind_1_slp = Ind_1_slp, Ind_1_mu = mean(Ind_1_mu) %>% rep(length(Year))), 
                                                        type = "response"),
                                      POF_mu = predict(rr6, newdata = list(Year = Year, Ind_1_slp = mean(Ind_1_slp) %>% rep(length(Year)), Ind_1_mu = Ind_1_mu), 
                                                       type = "response")) %>%
  reshape2::melt(value.name = "POF", id.vars = c("Ind_1_slp", "Ind_1_mu", "Year"), measure.vars = c("POF_slp", "POF_slp_mu")) %>%
  reshape2::melt(id.vars = c("variable", "POF", "Year", "POF_slp", "POF_mu"), variable.name = "Ind", value.name = "value")

ind_summary <- ind_POF %>% group_by(Year, variable) %>% summarise(POF_min = min(POF), POF_max = max(POF))

# Index target
uniroot(function(x) predict(rr6, newdata = list(Year = factor(2066), Ind_1_slp = 0, Ind_1_mu = x), type = "response") - 0.5, interval = c(0, 10))