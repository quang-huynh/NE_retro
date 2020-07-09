
#' @export
profile_fn <- function(x, args, .slot) {
  slot(args$OM, .slot) <- rep(x, 2)
  do.call(SRA_scope, args)
}

#' @export
get_profileLL <- function(x) {
  vapply(x, function(xx) ifelse(all(xx@conv), xx@mean_fit$opt$objective, NA), numeric(1))
}


fit_SR <- function(logR0, h, SSBPR0, sigmaR, SSB, Robs, opt = TRUE) {
  R0 <- exp(logR0)
  SSB0 <- R0 * SSBPR0
  Arec <- 4*h/(1-h)/SSBPR0
  Brec <- (5*h - 1)/(1-h)/SSB0
  
  Rpred <- Arec * SSB / (1 + Brec * SSB)
  
  if(opt) {
    return(sum(-dnorm(log(Robs/Rpred), -0.5 * sigmaR^2, sigmaR, log = TRUE)))
  } else {
    return(Rpred)
  }
}

get_Rexp <- function(x, SRA_h99, h) {
  nyears <- SRA_h99@OM@nyears
  SSB <- SRA_h99@SSB[x, 1:(nyears-1)]
  SSBPR0 <- SRA_h99@Misc[[1]]$EPR0[1:SRA_h99@mean_fit$obj$env$data$ageM] %>% mean()
  Robs <- SRA_h99@Misc[[x]]$R[2:nyears]
  opt <- optimize(fit_SR, log(c(1e-8, 2 * max(Robs))), h = h, SSBPR0 = SSBPR0, sigmaR = SRA_h99@OM@cpars$Perr[x], 
                  SSB = SSB, Robs = Robs)
  R0 <- exp(opt[[1]])
  Rexp <- fit_SR(opt[[1]], h = h, SSBPR0 = SSBPR0, sigmaR = SRA_h99@OM@cpars$Perr[x],
                 SSB = SSB, Robs = Robs, opt = FALSE)
  list(R0 = R0, SSB0 = R0 * SSBPR0, Rexp = Rexp)
}

# Condition OM and retrospectives with h = 0.99. But ensure that there is 
# a stock-recruit relationship in the projection period. Update h, historical Perr_y, R0, and D
#' @export
update_SRA_h <- function(SRA_h99, h, figure = TRUE) {
  maxage <- SRA_h99@OM@maxage
  nyears <- SRA_h99@OM@nyears
  
  OM <- SRA_h99@OM
  OM@h <- rep(h, 2)
  OM@cpars$h <- rep(h, OM@nsim)
  h99 <- runMSE(SRA_h99@OM, Hist = TRUE, silent = TRUE, parallel = SRA_h99@OM@nsim >= 48 & snowfall::sfIsRunning())
  
  SR_pars <- lapply(1:OM@nsim, get_Rexp, SRA_h99 = SRA_h99, h = h)
  R0 <- vapply(SR_pars, getElement, numeric(1), "R0")
  SSB0 <- vapply(SR_pars, getElement, numeric(1), "SSB0")
  Rexp <- do.call(rbind, lapply(SR_pars, getElement, "Rexp"))
  
  OM@cpars$R0 <- vapply(SR_pars, getElement, numeric(1), "R0")
  OM@cpars$D <- SRA_h99@SSB[, nyears]/SSB0
  
  Perr_y <- h99@AtAge$Nage[, 1, 2:nyears]/Rexp
  early_Perr_y <- SRA_h99@OM@cpars$Perr_y[, 1:maxage] * SRA_h99@OM@cpars$R0/R0
  
  OM@cpars$Perr_y[, 1:(maxage + nyears - 1)] <- cbind(early_Perr_y, Perr_y)
  
  if(figure) {
    SSB <- SRA_h99@SSB[1, 1:(ncol(SRA_h99@SSB) - 2)]
    Robs <- SRA_h99@Misc[[1]]$R[2:(ncol(SRA_h99@SSB)-1)]
    plot(SSB, Rexp[1, ], ylab = "Recruitment", xlim = c(0, 1.1 * max(c(SSB0, SSB))), typ = "l", col = "red",
         ylim = c(0, 1.1 * max(c(R0, Robs))))
    points(SSB, Robs)
    points(SSB0, R0, pch = 16, col = "red")
  }
  
  return(OM)
}

#' @export
generate_OM_args <- function(SRA, seed = 4, interval = 3, AddInd_from_residuals = TRUE, 
                             AC_from_residuals = TRUE, Cobs = 0.1, Cbias = 1) {
  SRA@OM@interval <- 3
  SRA@OM@Prob_staying <- SRA@OM@Frac_area_1 <- SRA@OM@Size_area_1 <- rep(0.5, 2)
  
  ###### Override indices
  n_ind <- SRA@OM@cpars$Data@AddInd %>% dim() %>% getElement(2)
  SRA@OM@cpars$AddIbeta <- matrix(1, SRA@OM@nsim, n_ind)
  
  Iobs <- SRA@OM@cpars$Data@AddInd[1, , ] %>% t()
  Ipred <- lapply(SRA@Misc, getElement, "Ipred")
  
  if(AddInd_from_residuals) {
    Isd <- rbind %>% do.call(lapply(Ipred, function(x, y) apply(log(y/x), 2, function(xx) sd(xx, na.rm = TRUE)), y = Iobs))
  } else {
    Isd <- sdconv(1, SRA@OM@cpars$Data@CV_AddInd[,,SRA@OM@nyears])
  }
  
  if(AC_from_residuals) {
    IAC <- rbind %>% do.call(lapply(Ipred, function(x, y) apply(log(y/x), 2, function(xx) acf(xx[!is.na(xx)], lag.max = 1, plot = FALSE)$acf[2, 1, 1]),
                                    y = Iobs))
  } else {
    IAC <- matrix(0, SRA@OM@nsim, ncol(Iobs))
  }
  
  set.seed(seed)
  I_dev_mu <- -0.5 * Isd^2 * (1 - IAC)/sqrt(1 - IAC^2)
  I_devs <- rnorm(SRA@OM@nsim * ncol(Iobs) * SRA@OM@proyears, I_dev_mu, Isd) %>%
    array(c(SRA@OM@nsim, ncol(Iobs), SRA@OM@proyears))
  
  hist_dev <- lapply(Ipred, function(x, y) y/x, y = Iobs) %>% unlist() %>% 
    array(c(SRA@OM@nyears, ncol(Iobs), SRA@OM@nsim)) %>% aperm(3:1)
  
  I_devs[,,1] <- IAC * log(hist_dev[, , SRA@OM@nyears]) + I_devs[,,1] * sqrt(1 - IAC^2)
  for(i in 2:dim(I_devs)[3]) I_devs[,,i] <- IAC * I_devs[,,i-1] + I_devs[,,i] * sqrt(1 - IAC^2)
  
  SRA@OM@cpars$AddIerr <- abind::abind(hist_dev, exp(I_devs), along = 3)
  
  ###### Remove catch from cpars$Data
  SRA@OM@cpars$Data@Cat <- matrix(NA, 1, 1)
  SRA@OM@Cobs <- rep(Cobs, 2)
  SRA@OM@TACFrac <- rep(Cbias, 2)
  if(ncol(SRA@data$Chist) == 1) {
    SRA@OM@CAA_ESS <- SRA@OM@CAA_nsamp <- SRA@data$CAA[SRA@OM@nyears, , 1] %>% sum() %>% rep(2)
  }
  SRA@OM@Cbiascv <- 1e-8
  
  return(SRA)
}

#' @export
generate_EM_args <- function(args, SRA_EM, SRA_OM, Ibias = 1, Cbias = 1) {
  #Hist <- runMSE(SRA_EM@OM, Hist = TRUE, parallel = snowfall::sfIsRunning())
  
  StockNames <- c("M_ageArray", "Mat_age", "Wt_age", "Len_age", "R0", "LenCV")
  StockPars <- lapply(StockNames, function(x) getElement(SRA_EM@OM@cpars, x)) %>% structure(names = StockNames)
  
  StockPars$R0 <- rep(args$OM@R0, args$OM@nsim)
  StockPars$procsd <- SRA_EM@OM@cpars$Perr
  StockPars$hs <- SRA_EM@OM@cpars$h
  
  StockPars$ageM <- matrix(SRA_EM@mean_fit$obj$env$data$ageM, args$OM@nsim, 1)
  StockPars$Linf <- rep(SRA_EM@mean_fit$obj$env$data$Linf, args$OM@nsim)
  
  FleetNames <- c("L5", "LFS", "Vmaxlen")
  FleetPars <- lapply(FleetNames, function(x) {
    slot(SRA_EM@OM, x)[1] %>% matrix(SRA_EM@OM@nyears + SRA_EM@OM@proyears, SRA_EM@OM@nsim)
  }) %>% structure(names = FleetNames)
  
  if(!is.null(args$vul_par) || any(args$selectivity == "free")) {
    vul_par <- args$vul_par
  } else {
    vul_par <- rbind(SRA_EM@mean_fit$report$LFS, SRA_EM@mean_fit$report$L5, SRA_EM@mean_fit$report$Vmaxlen)
  }
  dots <- list(map_log_early_rec_dev = args$map_log_early_rec_dev, 
               vul_par = vul_par, map_vul_par = args$map_vul_par,
               s_vul_par = args$s_vul_par, map_s_vul_par = args$map_s_vul_par)
  
  data_ind <- c("Chist", "CAA", "Index", "I_sd", "s_CAA", "sel_block", "nsel_block", "length_bin", "nfleet", "Ehist",
                "condition", "nsurvey", "CAL", "MS", "MS_cv", "MS_type", "C_eq", "E_eq", "s_CAL", "I_units", "abs_I", "age_error")
  args2 <- list(data = SRA_EM@data[match(data_ind, names(SRA_EM@data))], 
                selectivity = MSEtool:::int_sel(args$selectivity), 
                s_selectivity = MSEtool:::int_sel(args$s_selectivity), ESS = args$ESS,
                StockPars = StockPars, FleetPars = FleetPars, LWT = SRA_EM@data$LWT, dots = dots)
  
  args2$AddIndV <- SRA_EM@OM@cpars$Data@AddIndV[1, , ]
  
  # The following are the only things from the OM that is provided to the EM
  if(ncol(args2$data$Chist) > 1) { 
    args2$relF <- SRA_OM@mean_fit$report$F[nrow(SRA_OM@mean_fit$report$F), ]
    args2$V <- SRA_OM@mean_fit$report$vul[dim(SRA_OM@mean_fit$report$vul)[1], , ]
  }
  args2$Ibias <- Ibias
  args2$Cbias <- Cbias
  
  return(args2)
}


