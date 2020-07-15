library(MSEtool)
library(dplyr)
library(mfSCA)
library(ggplot2)

# Calc OM ref points with M = 0.02
SRA <- lapply(1:3, function(i) readRDS(paste0("GoM_cod/SRA_NR", i, ".rds")))
med_rec <- vapply(SRA, function(i) quantile(i@mean_fit$report$R, 0.75), numeric(1))

######## Reference points
ref_pt <- list()
ref_pt[[1]] <- lapply(1:3, function(i) calc_refpt(SRA[[i]], med_rec = med_rec[i]))
ref_pt[[2]] <- lapply(1:3, function(i) calc_refpt(SRA[[i]], M = "true", med_rec = med_rec[i]))
saveRDS(ref_pt, file = "GoM_cod/cod_ref_pt.rds")

ref_pt_plot <- do.call(rbind, lapply(ref_pt, function(i) do.call(rbind, lapply(i, function(ii) ii[1, ])))) %>% 
  as.data.frame()
ref_pt_plot$OM <- rep(paste0("NR", 1:3), 2)
ref_pt_plot$type <- rep(c("old", "new"), each = 3)
#ref_pt_plot <- ref_pt_plot[-4, ]

######## Load MSE
MSE <- lapply(1:3, function(i) readRDS(paste0("GoM_cod/MSE_cod_NR", i, ".rds")))

######## Plot OM
resOM <- plot_OM(MSE, MPs = c("M02", "M02_ra", "MRAMP", "MRAMP_ra", "ma", "NFref"), 0.10, 0.90, yfilter = 10)

theme_extra <- theme(panel.spacing = unit(0, "in"), legend.position="none")

ggplot(resOM[[1]], aes(Year, y = med, ymin = low, ymax = high)) + facet_grid(OM ~ MP, scales = "free_y") + 
  #geom_hline(yintercept = 0, col = "grey") + 
  geom_vline(xintercept = MSE[[1]]@OM$CurrentYr[1], linetype = 3) +
  geom_ribbon(fill = "grey80") + geom_line(size = 0.9) + 
  geom_hline(data = ref_pt_plot[4:6, ], aes(yintercept = FMSY), linetype = 2, size = 0.5) +
  gfplot::theme_pbs() + theme_extra + ylab("Fishing mortality")
ggsave("report/GoM_cod/MSE_OM_F.png", height = 3.5, width = 8.5)

ggplot(resOM[[2]], aes(Year, y = med, ymin = low, ymax = high)) + facet_grid(OM ~ MP, scales = "free_y") + 
  #geom_hline(yintercept = 0, col = "grey") + 
  geom_vline(xintercept = MSE[[1]]@OM$CurrentYr[1], linetype = 3) +
  geom_ribbon(fill = "grey80") + geom_line(size = 0.9) + 
  geom_hline(data = ref_pt_plot[4:6, ], aes(yintercept = SSBMSY), linetype = 2) +
  gfplot::theme_pbs() + theme_extra + ylab("SSB")
ggsave("report/GoM_cod/MSE_OM_SSB.png", height = 3.5, width = 8.5)

ggplot(resOM[[3]], aes(Year, y = med, ymin = low, ymax = high)) + facet_grid(OM ~ MP, scales = "free_y") + 
  #geom_hline(yintercept = 0, col = "grey") + 
  geom_vline(xintercept = MSE[[1]]@OM$CurrentYr[1], linetype = 3) +
  geom_ribbon(fill = "grey80") + geom_line(size = 0.9) + 
  geom_hline(data = ref_pt_plot[4:6, ], aes(yintercept = MSY), linetype = 3) +
  gfplot::theme_pbs() + theme_extra + ylab("Catch")
ggsave("report/GoM_cod/MSE_OM_C.png", height = 3.5, width = 8.5)


####### Plot rho 

res <- list()
res[[1]] <- plot_EM(MSE, MPs = c("M02", "M02_ra", "MRAMP", "MRAMP_ra")) %>% 
  lapply(function(x) {
    x$EM <- ifelse(substr(x$MP, nchar(x$MP)-1, nchar(x$MP)) == "ra", 
                   substr(x$MP, 1, nchar(x$MP)-3), x$MP)
    return(x)
  })

res[[2]] <- plot_EM(MSE, MPs = "ma") %>% lapply(function(x) cbind(x, EM = "M02"))
res[[3]] <- plot_EM(MSE, MPs = "ma", model_index = 2) %>% 
  lapply(function(x) cbind(x, EM = "MRAMP"))

out <- lapply(c("F_FMSY", "B_BMSY", "rho"), function(x) {
  y <- do.call(rbind, lapply(res, getElement, x))
  y$MP <- factor(y$MP, levels = c("M02", "M02_ra", "MRAMP", "MRAMP_ra", "ma"))
  return(y)
})

ggplot(out[[3]], aes(Year, y = `50%`, ymin = `25%`, ymax = `75%`, colour = EM)) + facet_grid(OM ~ MP, scales = "free_y") + 
  geom_hline(yintercept = 0, linetype = 3) + 
  geom_point(size = 0.75) + geom_linerange(size = 0.5) + 
  gfplot::theme_pbs() + ylab(expression(rho[SSB])) + theme(legend.position = "bottom")
ggsave("report/GoM_cod/cod_EM_rho.png", height = 4, width = 6.5)


######## Plot Index
Ind <- plot_Index(MSE, MPs = c("M02", "M02_ra", "MRAMP", "MRAMP_ra", "ma"), 
                  ind.names = c("Spring", "Fall", "MA"), yfilter = 10)
head(Ind[[1]])


ggplot(Ind[[1]], aes(Year, y = `50%`, ymin = `25%`, ymax = `75%`)) + facet_grid(OM ~ MP, scales = "free_y") + 
  #geom_hline(yintercept = 0, col = "grey") + 
  geom_vline(xintercept = MSE[[1]]@OM$CurrentYr[1], linetype = 3) +
  geom_ribbon(fill = "grey80") + geom_line(size = 0.9) +
  gfplot::theme_pbs() + theme_extra + ylab("NEFSC Fall survey")

ggplot(Ind[[2]], aes(Year, y = `50%`, ymin = `25%`, ymax = `75%`)) + facet_grid(OM ~ MP, scales = "free_y") + 
  #geom_hline(yintercept = 0, col = "grey") + 
  geom_vline(xintercept = MSE[[1]]@OM$CurrentYr[1], linetype = 3) +
  geom_ribbon(fill = "grey80") + geom_line(size = 0.9) +
  gfplot::theme_pbs() + theme_extra + ylab("MA survey")

ggplot(Ind[[3]], aes(Year, y = `50%`, ymin = `25%`, ymax = `75%`)) + facet_grid(OM ~ MP, scales = "free_y") + 
  #geom_hline(yintercept = 0, col = "grey") + 
  geom_vline(xintercept = MSE[[1]]@OM$CurrentYr[1], linetype = 3) +
  geom_ribbon(fill = "grey80") + geom_line(size = 0.9) +
  gfplot::theme_pbs() + theme_extra + ylab("NEFSC Spring survey")








# Performance metric
PM_fn <- function(x, y) {
  out <- data.frame(MP = x@MPs, PNOF = PNOF(x)@Mean, PB50 = P50(x)@Mean, 
                    LTOY = LTY(x, Ref = 0.75, Yrs = c(20, 50))@Mean,
                    STOY = LTY(x, Ref = 0.75, Yrs = c(1, 20))@Mean,
                    stringsAsFactors = FALSE)  %>% 
    filter(MP != "FMSYref" & MP != "NFref")
  out$OM <- paste0("NR", y)
  
  out %>% reshape2::melt(id.vars = c("MP", "OM"), variable.name = "PM")
}
pm <- do.call(rbind, Map(PM_fn, x = MSE, y = 1:3))
pm$MP <- pm$MP %>% factor(levels = c("M02", "M02_ra", "MRAMP", "MRAMP_ra", "ma"))

ggplot(pm, aes(PM, value, colour = PM)) + geom_point() + facet_grid(OM ~ MP) + gfplot::theme_pbs() + 
  geom_linerange(ymin = 0, aes(ymax = value)) + theme(legend.position = "bottom", axis.text.x = element_blank())
ggsave("MSE/cod/PM.png", height = 4, width = 6.5)

#ggplot(out[[1]], aes(Year, y = `50%`, ymin = `25%`, ymax = `75%`, colour = EM)) + facet_grid(OM ~ MP, scales = "free_y") + 
#  geom_hline(yintercept = 1, linetype = 3) + geom_point(size = 0.75) + geom_linerange(size = 0.5) +
#  gfplot::theme_pbs() + ylab(expression(paste("Estimated Terminal ", F/F[MSY]))) + theme(legend.position = "bottom")
#ggsave("MSE/cod/EM_F.png", height = 4, width = 6.5)
#
#ggplot(out[[2]], aes(Year, y = `50%`, ymin = `25%`, ymax = `75%`, colour = EM)) + facet_grid(OM ~ MP, scales = "free_y") + 
#  geom_hline(yintercept = 1, linetype = 3) + geom_point(size = 0.75) + geom_linerange(size = 0.5) + 
#  gfplot::theme_pbs() + ylab(expression(paste("Estimated Terminal ", SSB/SSB[MSY]))) + theme(legend.position = "bottom")
#ggsave("MSE/cod/EM_SSB.png", height = 4, width = 6.5)
