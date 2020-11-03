
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
