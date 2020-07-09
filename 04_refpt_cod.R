library(MSEtool)
library(dplyr)
library(mfSCA)
library(ggplot2)

# Calc OM ref points with M = 0.02
M02_refpt <- function(sra, M = c("M02", "true"), med_rec = 1) {
  M <- match.arg(M)
  if(M == "M02") {
    M_ageArray <- array(0.2, dim(sra@OM@cpars$Wt_age))
  } else {
    M_ageArray <- sra@OM@cpars$M_ageArray
  }
  
  ref_out <- lapply(1:sra@OM@nsim, mfSCA:::optMSY, 
                    M_ageArray = M_ageArray, Wt_age = sra@OM@cpars$Wt_age, Mat_age = sra@OM@cpars$Mat_age,
                    V = sra@OM@cpars$V, maxage = sra@OM@maxage, 
                    R0 = sra@OM@cpars$R0, SRrel = rep(sra@OM@SRrel, sra@OM@nsim), 
                    SSBPR0 = vapply(sra@Misc, getElement, numeric(1), "EPR0_SR"),
                    hs = sra@OM@cpars$h, yr.ind = sra@OM@nyears, plusgroup = 1)
  ref_pt <- lapply(c("F", "SB", "Yield"), function(x) vapply(ref_out, getElement, numeric(1), x))
  
  OY <- mfSCA:::MSYCalcs(log(0.75 * ref_pt[[1]][1]), M_at_Age = M_ageArray[1, , sra@OM@nyears], 
                         Wt_at_Age = sra@OM@cpars$Wt_age[1, , sra@OM@nyears], 
                         Mat_at_Age = sra@OM@cpars$Mat_age[1, , sra@OM@nyears], 
                         V_at_Age = sra@OM@cpars$V[1, , sra@OM@nyears],
                         maxage = sra@OM@maxage, 
                         R0x = sra@OM@cpars$R0[1], SRrelx = sra@OM@SRrel, Egg0 = sra@Misc[[1]]$EPR0_SR,
                         hx = sra@OM@cpars$h[1], opt = 2, plusgroup = 1)
  
  eq_fn <- function(logF) {
    Fout <- mfSCA:::MSYCalcs(logF, M_at_Age = M_ageArray[1, , sra@OM@nyears], 
                             Wt_at_Age = sra@OM@cpars$Wt_age[1, , sra@OM@nyears], 
                             Mat_at_Age = sra@OM@cpars$Mat_age[1, , sra@OM@nyears], 
                             V_at_Age = sra@OM@cpars$V[1, , sra@OM@nyears],
                             maxage = sra@OM@maxage, 
                             R0x = sra@OM@cpars$R0[1], SRrelx = sra@OM@SRrel, Egg0 = sra@Misc[[1]]$EPR0_SR,
                             hx = sra@OM@cpars$h[1], opt = 2, plusgroup = 1)
    return(Fout)
  }
  FF <- seq(0.001, 3, 0.001)
  eq_curve <- lapply(log(FF), eq_fn)
  RR <- sapply(eq_curve, function(x) x["RelRec"])
  Fcrash <- FF[which(RR <= 0)[1]]
  
  # F40 %
  ref_pt$Fproxy <- vapply(1:sra@OM@nsim, function(x) {
    uniroot(mfSCA:::SPR_opt, log(c(1e-8, 3)), SPR_ratio = 0.4, wt = sra@OM@cpars$Wt_age[x, , sra@OM@nyears], 
            mat = sra@OM@cpars$Mat_age[x, , sra@OM@nyears], M = M_ageArray[x, , sra@OM@nyears], 
            V = sra@OM@cpars$V[x, , sra@OM@nyears], 
            max_age = sra@OM@maxage, plusgroup = TRUE, R = 1, opt = TRUE)[[1]] %>% exp()
  }, numeric(1))
  
  # SSBPR40%
  ref_pt$SSBPRproxy <- vapply(1:sra@OM@nsim, function(x) {
    mfSCA:::SPR_opt(log(ref_pt$Fproxy[x]), wt = sra@OM@cpars$Wt_age[x, , sra@OM@nyears], 
                    mat = sra@OM@cpars$Mat_age[x, , sra@OM@nyears], M = M_ageArray[x, , sra@OM@nyears], 
                    V = sra@OM@cpars$V[x, , sra@OM@nyears], 
                    max_age = sra@OM@maxage, plusgroup = TRUE, R = med_rec, opt = FALSE)[[2]]
  }, numeric(1))
  
  ref_pt$YPRproxy <- vapply(1:sra@OM@nsim, function(x) {
    mfSCA:::YPR_opt(log(ref_pt$Fproxy[x]), wt = sra@OM@cpars$Wt_age[x, , sra@OM@nyears], 
                    M = M_ageArray[x, , sra@OM@nyears], V = sra@OM@cpars$V[x, , sra@OM@nyears], 
                    max_age = sra@OM@maxage, plusgroup = TRUE, R = med_rec, opt = FALSE)
  }, numeric(1))
  
  ref_pt$YPR75proxy <- vapply(1:sra@OM@nsim, function(x) {
    mfSCA:::YPR_opt(log(0.75 * ref_pt$Fproxy[x]), wt = sra@OM@cpars$Wt_age[x, , sra@OM@nyears], 
                    M = M_ageArray[x, , sra@OM@nyears], V = sra@OM@cpars$V[x, , sra@OM@nyears], 
                    max_age = sra@OM@maxage, plusgroup = TRUE, R = med_rec, opt = FALSE)
  }, numeric(1))
  
  out <- do.call(cbind, ref_pt) %>% cbind(rep(OY["Yield"], length(ref_pt[[1]]))) %>%
    cbind(rep(Fcrash, length(ref_pt[[1]]))) %>%
    structure(dimnames = list(NULL, c("FMSY", "SSBMSY", "MSY", "Fproxy", "SSBproxy", "MSYproxy", 
                                      "OYproxy", "OY", "Fcrash")))
  
  return(out)
}

SRA <- lapply(1:3, function(i) readRDS(paste0("GoM_cod/SRA_NR", i, ".rds")))
med_rec <- vapply(SRA, function(i) quantile(i@mean_fit$report$R, 0.75), numeric(1))

#ref_pt <- lapply(1:3, function(i) M02_refpt(SRA[[i]], med_rec = med_rec[i]))
ref_pt_M02 <- lapply(1:3, function(i) M02_refpt(SRA[[i]], med_rec = med_rec[i]))
ref_pt_newM <- lapply(1:3, function(i) M02_refpt(SRA[[i]], M = "true", med_rec = med_rec[i]))




######## Plot EM
#MSE <- lapply(1:3, function(i) {
#  m1 <- readRDS(paste0("MSE/MSE_cod_test1_NR", i, "_t1.rds"))
#  m2 <- readRDS(paste0("MSE/MSE_cod_test1_NR", i, "_t2.rds"))
#  out <- merge_MSE(m1, m2)
#  saveRDS(out, file = paste0("MSE/MSE_cod_NR", i, ".rds"))
#})
MSE <- lapply(1:3, function(i) readRDS(paste0("MSE/MSE_cod_NR", i, ".rds")))

### update ref pt
MSE <- Map(function(x, y, ind = c(1, 2, 5)) {
  x@F_FMSY <- x@FM / y[, 1] #FMSY
  x@B_BMSY <- x@SSB / y[, 2] #$SSBMSY
  x@OM$RefY <- x@OM$MSY <- y[, 5] #OY
  return(x)
}, x = MSE, y = ref_pt)

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




ggplot(out[[1]], aes(Year, y = `50%`, ymin = `25%`, ymax = `75%`, colour = EM)) + facet_grid(OM ~ MP, scales = "free_y") + 
  geom_hline(yintercept = 1, linetype = 3) + geom_point(size = 0.75) + geom_linerange(size = 0.5) +
  gfplot::theme_pbs() + ylab(expression(paste("Estimated Terminal ", F/F[MSY]))) + theme(legend.position = "bottom")
ggsave("MSE/cod/EM_F.png", height = 4, width = 6.5)

ggplot(out[[2]], aes(Year, y = `50%`, ymin = `25%`, ymax = `75%`, colour = EM)) + facet_grid(OM ~ MP, scales = "free_y") + 
  geom_hline(yintercept = 1, linetype = 3) + geom_point(size = 0.75) + geom_linerange(size = 0.5) + 
  gfplot::theme_pbs() + ylab(expression(paste("Estimated Terminal ", SSB/SSB[MSY]))) + theme(legend.position = "bottom")
ggsave("MSE/cod/EM_SSB.png", height = 4, width = 6.5)

ggplot(out[[3]], aes(Year, y = `50%`, ymin = `25%`, ymax = `75%`, colour = EM)) + facet_grid(OM ~ MP, scales = "free_y") + 
  geom_hline(yintercept = 0, linetype = 3) + geom_point(size = 0.75) + geom_linerange(size = 0.5) + 
  gfplot::theme_pbs() + ylab(expression(rho[SSB])) + theme(legend.position = "bottom")
ggsave("MSE/cod/EM_rho.png", height = 4, width = 6.5)



res <- plot_OM(MSE, MPs = c("M02", "M02_ra", "MRAMP", "MRAMP_ra", "ma"), ref = ref_pt_newM)

ggplot(res[[1]], aes(Year, y = med, ymin = low, ymax = high)) + facet_grid(OM ~ MP, scales = "free_y") + 
  geom_point(size = 0.5) + geom_linerange(size = 0.25) + gfplot::theme_pbs() + ylab(expression(F/F[MSY])) +
  geom_hline(yintercept = 1, linetype = 2)
ggsave("MSE/cod/OM_relF.png", height = 3.5, width = 6.5) 

ggplot(res[[2]], aes(Year, y = med, ymin = low, ymax = high)) + facet_grid(OM ~ MP, scales = "free_y") + 
  geom_point(size = 0.5) + geom_linerange(size = 0.25) + gfplot::theme_pbs() + ylab(expression(SSB/SSB[MSY]))
ggsave("MSE/cod/OM_relB.png", height = 3.5, width = 6.5)

ggplot(res[[3]], aes(Year, y = med, ymin = low, ymax = high)) + facet_grid(OM ~ MP, scales = "free_y") + 
  geom_point(size = 0.5) + geom_linerange(size = 0.25) + gfplot::theme_pbs() + ylab("Catch/OY") 
ggsave("MSE/cod/OM_relC.png", height = 3.5, width = 6.5)



res <- plot_OM(MSE, MPs = c("M02", "M02_ra", "MRAMP", "MRAMP_ra", "ma", "NFref"), ref_pt = "abs")

ggplot(res[[1]], aes(Year, y = med, ymin = low, ymax = high)) + facet_grid(OM ~ MP, scales = "free_y") + 
  geom_point(size = 0.5) + geom_line() + geom_linerange(size = 0.25) + gfplot::theme_pbs() + ylab("Fishing mortality") +
  geom_vline(xintercept = MSE[[1]]@nyears, linetype = 3)
ggsave("MSE/cod/OM_F.png", height = 3.5, width = 6.5)

ggplot(res[[2]], aes(Year, y = med, ymin = low, ymax = high)) + facet_grid(OM ~ MP, scales = "free_y") + 
  geom_point(size = 0.5) + geom_line() + geom_linerange(size = 0.25) + gfplot::theme_pbs() + ylab("SSB") +
  geom_vline(xintercept = MSE[[1]]@nyears, linetype = 3)
ggsave("MSE/cod/OM_B.png", height = 3.5, width = 6.5)

ggplot(res[[3]], aes(Year, y = med, ymin = low, ymax = high)) + facet_grid(OM ~ MP, scales = "free_y") + 
  geom_point(size = 0.5) + geom_line() + geom_linerange(size = 0.25) + gfplot::theme_pbs() + ylab("Catch") +
  geom_vline(xintercept = MSE[[1]]@nyears, linetype = 3)
  
ggsave("MSE/cod/OM_C.png", height = 3.5, width = 6.5)


res <- plot_Index(MSE, MPs = c("M02", "ma"), ind.names = c("Spring", "Fall", "MA"))

summary(MSE[[1]])



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
