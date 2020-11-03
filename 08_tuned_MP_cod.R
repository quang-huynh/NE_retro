
library(MSEtool)
library(mfSCA)
library(dplyr)
library(ggplot2)

######## OM order
OM_order <- c(1, 3, 2)
OM_names <- c("MC", "MCIM", "IM")

######## Reference points
ref_pt <- readRDS(file = "GoM_cod/cod_ref_pt.rds")

ref_pt_plot <- do.call(rbind, lapply(ref_pt, function(i) do.call(rbind, lapply(i, function(ii) ii[1, ])))) %>% 
  as.data.frame()
ref_pt_plot$OM <- rep(OM_names, 2)
ref_pt_plot$type <- rep(c("old", "new"), each = 3)


######## Load MSE
MSE <- lapply(1:3, function(x) {
  y <- readRDS(paste0("GoM_cod/MSE_cod_tuned_MP_NR", x, ".rds"))
  y@OM$SSBMSY <- filter(ref_pt_plot, type == 'new')$SSBMSY[x]
  y@B_BMSY <- y@SSB/y@OM$SSBMSY
  
  y@OM$FMSY <- filter(ref_pt_plot, type == 'new')$FMSY[x]
  y@F_FMSY <- y@FM/y@OM$FMSY 
  
  y@OM$RefY <- filter(ref_pt_plot, type == 'new')$OY[x]
  
  return(y)
})


######## Plot OM
resOM <- plot_OM(MSE[OM_order], MPs = c("tMP_20cap", "tMP_nocap"), 0.10, 0.90, yfilter = 10, OM_names = OM_names[OM_order])

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
#ggsave("report/GoM_cod/MSE_OM_F.png", height = 3.5, width = 8.5)

ggplot(resOM[[2]], aes(Year, y = med, ymin = low, ymax = high)) + facet_grid(factor(OM) ~ MP, scales = "free_y") +  
  #geom_hline(yintercept = 0, col = "grey") + 
  geom_vline(xintercept = MSE[[1]]@OM$CurrentYr[1], linetype = 3) +
  geom_ribbon(fill = "grey80") + geom_line(size = 0.9) + 
  geom_hline(data = ref_pt_plot[4:6, ], aes(yintercept = SSBMSY), linetype = 2) +
  geom_hline(data = ref_pt_plot[4:6, ], aes(yintercept = 0.5 * SSBMSY), linetype = 2) +
  ylab("SSB") + scale_y_continuous(labels = scales::comma, expand = c(0.01, 0.01)) +
  gfplot::theme_pbs() + no_panel_gap
#ggsave("report/GoM_cod/MSE_OM_SSB.png", height = 3.5, width = 8.5)

ggplot(resOM[[3]], aes(Year, y = med, ymin = low, ymax = high)) + facet_grid(factor(OM) ~ MP, scales = "free_y") + 
  #geom_hline(yintercept = 0, col = "grey") + 
  geom_vline(xintercept = MSE[[1]]@OM$CurrentYr[1], linetype = 3) +
  geom_ribbon(fill = "grey80") + geom_line(size = 0.9) + 
  geom_hline(data = ref_pt_plot[4:6, ], aes(yintercept = OY), linetype = 3) +
  ylab("Catch") + scale_y_continuous(labels = scales::comma, expand = c(0.01, 0.01)) +
  gfplot::theme_pbs() + no_panel_gap
#ggsave("report/GoM_cod/MSE_OM_C.png", height = 3.5, width = 8.5)

