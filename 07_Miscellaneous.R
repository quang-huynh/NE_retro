
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
  gfplot::theme_pbs() + no_legend + theme(strip.text.x = element_blank(), strip.text.y = element_blank())
b <- ggplot(TAC, aes(x = Model, y = TAC2)) + facet_wrap(~facet, ncol = 1, scales = "free") +
  geom_line() + geom_point(aes(colour = factor(Model))) + 
  geom_hline(yintercept = 0, col = "white") + labs(x = "Year", y = "Relative catch advice") + 
  gfplot::theme_pbs() + no_legend + theme(strip.text.x = element_blank(), strip.text.y = element_blank())

plot_grid(a, b, nrow = 1)
ggsave("report/intro_figure.png", height = 3, width = 5)
