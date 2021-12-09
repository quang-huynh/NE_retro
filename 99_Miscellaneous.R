
library(dplyr)
library(ggplot2)

no_legend <- theme(legend.position = "none")
no_panel_gap <- theme(panel.spacing = unit(0, "in"))
legend_bottom <- theme(legend.position = "bottom")


halibut <- read.csv('halibut_retro.csv') %>% mutate(facet = ifelse(Model <= 1994, 1, 2))

TAC <- halibut %>% group_by(Model) %>% summarise(TAC = 0.1 * VB[which.max(Year)]) %>% 
  mutate(facet = ifelse(Model <= 1994, 1, 2), TAC2 = TAC/ifelse(facet == 1, 23.3, 33.6))
terminal <- halibut %>% group_by(Model) %>% summarise(VB = VB[which.max(Year)]) %>% 
  mutate(facet = ifelse(Model <= 1994, 1, 2), Year = Model)
  

a <- ggplot(halibut, aes(x = Year, y = VB, group = Model, colour = factor(Model))) + 
  facet_wrap(~facet, ncol = 1, scales = "free") + geom_line() + geom_point(data = terminal) + 
  geom_hline(yintercept = 0, col = "white") + 
  geom_vline(data = data.frame(facet = c(1, 2), yy = c(1996, 1995)), aes(xintercept = yy), col = "white") +
  ylab("Biomass estimate") + 
  geom_text(data = data.frame(x = Inf, y = Inf, facet = 1:2, txt = c("(a)", "(c)")), 
            vjust = 1.1, hjust = 1.1,
            inherit.aes = FALSE, aes(x = x, y = y, label = txt)) + 
  gfplot::theme_pbs() + no_legend + theme(strip.text.x = element_blank(), strip.text.y = element_blank())
b <- ggplot(TAC, aes(x = Model, y = TAC2)) + facet_wrap(~facet, ncol = 1, scales = "free") +
  geom_line() + geom_point(aes(colour = factor(Model))) + 
  geom_hline(yintercept = 0, col = "white") + 
  geom_text(data = data.frame(x = Inf, y = Inf, facet = 1:2, txt = c("(b)", "(d)")), 
            vjust = 1.1, hjust = 1.1,
            inherit.aes = FALSE, aes(x = x, y = y, label = txt)) + 
  labs(x = "Year", y = "Relative catch advice") + 
  coord_cartesian(ylim = c(0, 1.5)) + 
  gfplot::theme_pbs() + no_legend + theme(strip.text.x = element_blank(), strip.text.y = element_blank())

cowplot::plot_grid(a, b, nrow = 1)
ggsave("report/intro_figure.png", height = 3, width = 5)


# Well - behaved assessments

make_well_behaved <- function(yy, halibut) {
  #browser()
  x <- x_out <- halibut %>% filter(Model == yy[1]-1)
  for(i in 1:length(yy)) {
    y <- halibut %>% filter(Model == yy[i])
    
    new_years <- y$VB[y$Year > yy[i]-1]/y$VB[y$Year == yy[1]-1]
    x_terminal <- x$VB[x$Year == yy[1]-1]
    
    x <- y %>% mutate(VB = c(x$VB, new_years * x_terminal))
    x_out <- rbind(x_out, x)
  }
  return(x_out)
}


out <- make_well_behaved(1990:1994, halibut = halibut %>% filter(facet == 1)) %>% 
  mutate(Model = factor(Model, 1994:1989))
terminal <- out %>% group_by(Model) %>% summarise(VB = VB[which.max(Year)]) %>% mutate(Year = 1994:1989)


ggplot(out, aes(x = Year, y = VB, group = Model)) + 
  geom_line(aes(colour = factor(Model))) + geom_point(data = terminal, aes(colour = factor(Model))) + 
  geom_text(data = terminal, nudge_x = -0.5, hjust = "right", aes(label = Year)) + 
  geom_hline(yintercept = 0, col = "white") + 
  ylab("Biomass estimate") + coord_cartesian(xlim = c(1970, 1996)) +
  gfplot::theme_pbs() + no_legend + theme(strip.text.x = element_blank(), strip.text.y = element_blank())
ggsave("report/ppt/well_behaved.png", width = 4, height = 3.25)


terminal <- halibut %>% group_by(Model) %>% summarise(VB = VB[which.max(Year)]) %>% 
  mutate(facet = ifelse(Model <= 1994, 1, 2), Year = Model)

ggplot(halibut %>% filter(facet == 1) %>% mutate(Model = factor(Model, 1994:1989)), 
       aes(x = Year, y = VB, group = Model)) + 
  geom_line(aes(colour = factor(Model))) + geom_point(data = terminal %>% filter(Model <= 1995), aes(colour = factor(Model))) + 
  #geom_text(data = terminal, nudge_x = -0.5, hjust = "right", aes(label = Year)) + 
  geom_hline(yintercept = 0, col = "white") + 
  ylab("Biomass estimate") + coord_cartesian(xlim = c(1970, 1996)) +
  gfplot::theme_pbs() + no_legend + theme(strip.text.x = element_blank(), strip.text.y = element_blank())
ggsave("report/ppt/retro.png", width = 4, height = 3.25)


ggplot(halibut %>% filter(facet == 2) %>% mutate(Model = factor(Model, 2007:2012)), 
       aes(x = Year, y = VB, group = Model)) + 
  geom_line(aes(colour = factor(Model))) + geom_point(data = terminal %>% filter(Model > 1995), aes(colour = factor(Model))) + 
  #geom_text(data = terminal, nudge_x = -0.5, hjust = "right", aes(label = Year)) + 
  geom_hline(yintercept = 0, col = "white") + 
  ylab("Biomass estimate") + #coord_cartesian(xlim = c(1970, 1996)) +
  gfplot::theme_pbs() + no_legend + theme(strip.text.x = element_blank(), strip.text.y = element_blank())
ggsave("report/ppt/retro2.png", width = 4, height = 3.25)


###### Single-panel indicator
ggplot(indicators_cast %>% filter(MP == "M02_ra" & Year == 2030), aes(PMat_1_mu, SSB_rho)) + 
  facet_grid(Year ~ MP) + geom_hline(yintercept = 0, linetype = 3) + geom_point(aes(colour = OM), alpha = 0.6) +
  geom_contour(data = rr_pred2 %>% filter(MP == "M02_ra" & Year == 2030), breaks = 0.5, colour = "black", aes(z = MC)) + 
  #metR::geom_label_contour(data = rr_pred2, breaks = 0.5, colour = "black", aes(z = MC)) +
  geom_text(data = misclass %>% filter(MP == "M02_ra" & Year == 2030), aes(label = class_correct), x = 0, y = 3, vjust = "inward", hjust = "inward") +
  coord_cartesian(xlim = c(-0.5, -0), ylim = c(-0.5, 3)) + labs(x = "Prop. mature (log)", y = expression(rho[SSB])) +
  scale_x_continuous(breaks = c(-0.4, -0.2, 0)) +
  gfplot::theme_pbs() + no_panel_gap + legend_bottom
ggsave("report/ppt/LDA_singlepanel.png", height = 3, width = 3)



###### Combine comparison plot with both cod and pollock
no_panel_gap <- theme(panel.spacing = unit(0, "in"))
legend_bottom <- theme(legend.position = "bottom")


cod <- local({
  
  res1 <- readRDS("GoM_cod/SRA_cod_M02.rds")
  res2 <- readRDS("GoM_cod/SRA_cod_MRAMP.rds")
  res3 <- readRDS("GoM_cod/SRA_NR1.rds")
  res4 <- readRDS("GoM_cod/SRA_NR2.rds")
  res5 <- readRDS("GoM_cod/SRA_NR3.rds")
  
  Model <- data.frame(name = c("MC", "IM", "MCIM", "M02", "MRAMP"), type = c("OM", "OM", "OM", "AM", "AM"))
  SSB <- lapply(list(res3, res4, res5, res1, res2), function(x) x@SSB[1, ]) %>% do.call(cbind, .) %>% 
    structure(dimnames = list(Year = 1982:2019, name = Model$name)) %>% reshape2::melt() %>% mutate(y = "SSB") %>%
    left_join(Model, by = "name")
  
  RR <- lapply(list(res3, res4, res5, res1, res2), function(x) x@Misc[[1]]$R) %>% do.call(cbind, .) %>% 
    structure(dimnames = list(Year = 1982:2019, name = Model$name)) %>% reshape2::melt() %>% mutate(y = "Recruitment") %>%
    left_join(Model, by = "name")
  
  rbind(SSB, RR) %>% 
    mutate(y = factor(y, levels = c("SSB", "Recruitment")), 
           name = factor(name, levels = Model$name), Stock = "Cod") %>% 
    filter(value <= 1e6) %>%
    ggplot(aes(x = Year, y = value, linetype = name, colour = name)) + geom_line() +
    geom_hline(yintercept = 0, colour = "white") +
    facet_grid(y ~ Stock, scales = "free_y", switch = "y") + ylab(NULL) + 
    geom_text(data = data.frame(x = Inf, yy = Inf, y = c("SSB", "Recruitment"), txt = c("(a)", "(b)")),
              aes(x, yy, label = txt), hjust = 1.1, vjust = 1.1,
              inherit.aes = FALSE) + 
    scale_y_continuous(labels = scales::comma, n.breaks = 5) +
    scale_linetype_manual(values = c("MC" = 1, "IM" = 1, "MCIM" = 1, "M02" = 2, "MRAMP" = 2)) +
    gfplot::theme_pbs() + no_panel_gap + legend_bottom + 
    theme(legend.title = element_blank(), strip.placement = "outside")
})


pollock <- local({
  
  res1 <- readRDS("pollock/SRA_pollock_base.rds")
  res2 <- readRDS("pollock/SRA_pollock_flatsel.rds")
  res3 <- readRDS("pollock/SRA_NR1.rds")
  res4 <- readRDS("pollock/SRA_NR2.rds")
  res5 <- readRDS("pollock/SRA_NR3.rds")
  
  # SSB and R
  Model <- data.frame(name = c("SS", "SWB", "SWF", "Base", "FlatSel"), type = c("OM", "OM", "OM", "AM", "AM"))
  SSB <- lapply(list(res3, res4, res5, res1, res2), function(x) x@SSB[1, ]) %>% do.call(cbind, .) %>% 
    structure(dimnames = list(Year = 1970:2019, name = Model$name)) %>% reshape2::melt() %>% mutate(y = "SSB") %>%
    left_join(Model, by = "name")
  
  RR <- lapply(list(res3, res4, res5, res1, res2), function(x) x@Misc[[1]]$R) %>% do.call(cbind, .) %>% 
    structure(dimnames = list(Year = 1970:2019, name = Model$name)) %>% reshape2::melt() %>% mutate(y = "Recruitment") %>%
    left_join(Model, by = "name")
  
  rbind(SSB, RR) %>% 
    mutate(y = factor(y, levels = c("SSB", "Recruitment")), 
           name = factor(name, levels = Model$name), Stock = "Pollock") %>% 
    filter(value <= 1e6) %>%
    ggplot(aes(x = Year, y = value, linetype = name, colour = name)) + geom_line() +
    geom_hline(yintercept = 0, colour = "white") +
    geom_text(data = data.frame(x = Inf, yy = Inf, y = c("SSB", "Recruitment"), txt = c("(c)", "(d)")),
              aes(x, yy, label = txt), hjust = 1.1, vjust = 1.1,
              inherit.aes = FALSE) + 
    facet_grid(y ~ Stock, scales = "free_y", switch = "y") + ylab(NULL) +
    scale_y_continuous(labels = scales::comma, n.breaks = 5) +
    scale_linetype_manual(values = c("SS" = 1, "SWB" = 1, "SWF" = 1, "Base" = 2, "FlatSel" = 2)) +
    gfplot::theme_pbs() + no_panel_gap + legend_bottom + 
    theme(legend.title = element_blank(), strip.placement = "outside")
})


ggpubr::ggarrange(cod, pollock, ncol = 2)
ggsave("report/model_compare.png", height = 4, width = 8)

