

################# Compare
res1 <- readRDS("pollock/SRA_pollock_base.rds")
res2 <- readRDS("pollock/SRA_pollock_flatsel.rds")
res3 <- readRDS("pollock/SRA_SS.rds")
res4 <- readRDS("pollock/SRA_SWB.rds")
res5 <- readRDS("pollock/SRA_SWF.rds")


compare_SRA(res1, res2, res3, res4, res5,
            filename = "report/pollock/compare_pollock_ret", dir = getwd(),
            title = "ASAP3 pollock comparison",
            s_name = c("NEFSCspring", "NEFSCfall"), f_name = c("Commercial", "Recreational"),
            scenario = list(names = c("Base", "FlatSel", 
                                      "Base_LowSurvSel", "Base_UpSurv", "FlatSel_UpSurv")))


no_panel_gap <- theme(panel.spacing = unit(0, "in"))
legend_bottom <- theme(legend.position = "bottom")


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
  geom_text(data = rho_label, x = -Inf, y = 0, vjust = "inward", hjust = -0.1, aes(label = rho), parse = TRUE) +
  geom_line(data = filter(rr_out, Peel == 0), size = 0.75, colour = "black") + 
  geom_text(data = data.frame(Model = model_name, label = paste0("(", letters[1:5], ")")),
            aes(label = label),
            inherit.aes = FALSE, x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1
  ) + 
  coord_cartesian(ylim = c(0, 6e5)) + 
  scale_y_continuous(labels = scales::comma) + scale_x_continuous(breaks = c(1970, 1990, 2010)) + 
  gfplot::theme_pbs() + no_panel_gap + legend_bottom
ggsave("report/pollock/pollock_cond_rho.png", width = 6.5, height = 5)


