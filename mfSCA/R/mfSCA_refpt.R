
SPR_opt <- function(logF, SPR_ratio, wt, mat, M, V, max_age, plusgroup, R = 1, opt = TRUE) {
  Z <- V * exp(logF) + M
  NPR_F <- c(1, exp(-cumsum(Z[1:(max_age-1)])))
  if(plusgroup) NPR_F[max_age] <- NPR_F[max_age]/(1 - exp(-Z[max_age]))
  SPR_F <- sum(NPR_F * wt * mat)
  
  NPR_unfished <- c(1, exp(-cumsum(M[1:(max_age-1)])))
  if(plusgroup) NPR_unfished[max_age] <- NPR_unfished[max_age]/(1 - exp(-M[max_age]))
  SPR0 <- sum(NPR_unfished * wt * mat)
  
  if(opt) {
    return(SPR_F/SPR0 - SPR_ratio)
  } else {
    VPR_F <- sum(NPR_F * wt * V)
    BPR_F <- sum(NPR_F * wt) 
    return(list(F = V * exp(logF), E = SPR_F * R, VB = VPR_F * R, B = BPR_F * R, SPR0 = SPR0))
  }
}

YPR_opt <- function(logF, wt, M, V, max_age, plusgroup, R = 1, opt = TRUE) {
  FF <- V * exp(logF)
  Z <- FF + M
  NPR_F <- c(1, exp(-cumsum(Z[1:(max_age-1)])))
  if(plusgroup) NPR_F[max_age] <- NPR_F[max_age]/(1 - exp(-Z[max_age]))
  CPR_F <- NPR_F * FF / Z * (1 - exp(-Z))
  Yield <- R * sum(CPR_F * wt)
  
  if(opt) {
    return(-1 * Yield)
  } else {
    return(Yield)
  }
}

F01_calc <- function(logF, wt, M, V, max_age, plusgroup) {
  # Note the grad fn returns dY/dlogF, but dY/dF = dY/dlogF * dlogF/dF where dlogF/dF = 1/F
  deriv_F0 <- numDeriv::grad(YPR_opt, x = log(1e-8), wt = wt, M = M, V = V, max_age = max_age, 
                             plusgroup = plusgroup, opt = FALSE) / 1e-8
  deriv_F <- numDeriv::grad(YPR_opt, x = logF, wt = wt, M = M, V = V, max_age = max_age, 
                            plusgroup = plusgroup, opt = FALSE) / exp(logF)
  return(deriv_F - 0.1 * deriv_F0)
}

mf_SCA_ref_pt <- function(SRA_out, R_proxy, F_proxy = c("SPR", "F01", "Fmax"), SPR,
                          M_basis = c("init", "terminal")) {
  F_proxy <- match.arg(F_proxy)
  M_basis <- match.arg(M_basis)
  if(F_proxy == "SPR" && missing(SPR)) SPR <- 0.4
  
  nyears <- SRA_out$obj$env$data$n_y
  max_age <- SRA_out$obj$env$data$max_age
  
  if(M_basis == "init") M_y <- 1 else M_y <- nyears
  
  plusgroup <- SRA_out$obj$env$data$plusgroup
  M <- SRA_out$obj$env$data$M[M_y, ]
  wt <- SRA_out$obj$env$data$wt[nyears, ]
  mat <- SRA_out$obj$env$data$mat[nyears, ]
  V <- SRA_out$report$F_at_age[nyears, ]/max(SRA_out$report$F_at_age[nyears, ])
  
  if(F_proxy == "SPR") {
    opt <- uniroot(SPR_opt, log(c(1e-8, 3)), SPR_ratio = SPR, wt = wt, mat = mat, M = M, V = V,
                   max_age = max_age, plusgroup = plusgroup)
  } else if(F_proxy == "Fmax") {
    opt <- optimize(YPR_opt, log(c(1e-8, 3)), wt = wt, M = M, V = V, max_age = max_age, plusgroup = plusgroup)
  } else {
    opt <- uniroot(F01_calc, log(c(1e-6, 3)), wt = wt, M = M, V = V, max_age = max_age, plusgroup = plusgroup)
  }
  output <- SPR_opt(opt[[1]], wt = wt, mat = mat, M = M, V = V, 
                    max_age = max_age, plusgroup = plusgroup, R = R_proxy, opt = FALSE)
  output$Y <- YPR_opt(opt[[1]], wt = wt, M = M, V = V, max_age = max_age, plusgroup = plusgroup, 
                      R = R_proxy, opt = FALSE)
  return(output)
}

MSYCalcs <- function(logF, M_at_Age, Wt_at_Age, Mat_at_Age, V_at_Age, maxage, 
                     R0x, SRrelx, Egg0, hx, opt = 1, plusgroup = 0) {
  FF <- exp(logF)
  lx <- rep(1, maxage)
  surv <- exp(-M_at_Age - FF * V_at_Age)
  for (a in 2:maxage) {
    lx[a] <- lx[a - 1] * surv[a - 1]
  }
  if (plusgroup == 1) {
    lx[length(lx)] <- lx[length(lx)]/(1 - surv[length(lx)])
  }
  EggF <- SBF <- sum(lx * Wt_at_Age * Mat_at_Age)
  vBF <- sum(lx * Wt_at_Age * V_at_Age)
  BF <- sum(lx * Wt_at_Age)
  hx[hx > 0.999] <- 0.999
  recK <- (4 * hx)/(1 - hx)
  reca <- recK/Egg0
  SPR <- EggF/Egg0
  if (SRrelx == 1) {
    recb <- (reca * Egg0 - 1)/(R0x * Egg0)
    RelRec <- (reca * EggF - 1)/(recb * EggF)
  }
  if (SRrelx == 2) {
    bR <- (log(5 * hx)/(0.8 * Egg0))
    aR <- exp(bR * Egg0)/(Egg0/R0x)
    RelRec <- (log(aR * EggF/R0x))/(bR * EggF/R0x)
  }
  RelRec[RelRec < 0] <- 0
  Z_at_Age <- FF * V_at_Age + M_at_Age
  YPR <- sum(lx * Wt_at_Age * FF * V_at_Age * (1 - exp(-Z_at_Age))/Z_at_Age)
  Yield <- YPR * RelRec
  if (opt == 1) 
    return(-Yield)
  if (opt == 2) {
    out <- c(Yield = Yield, F = FF, SB = SBF * RelRec, SB_SB0 = (SBF * RelRec)/(Egg0 * R0x), 
             #B_B0 = (BF * RelRec)/(Egg0 * R0x), 
             B = BF * RelRec, VB = vBF * RelRec, 
             #VB_VB0 = (vBF * RelRec)/(vB0 * R0x), 
             RelRec = RelRec, SB0 = Egg0 * R0x, SPR = SPR) #B0 = B0 * R0x)
    return(out)
  }
}

optMSY <- function(x, M_ageArray, Wt_age, Mat_age, V, maxage, R0, SRrel, 
                   hs, SSBPR0, yr.ind = 1, plusgroup = 0) {
  if (length(yr.ind) == 1) {
    M_at_Age <- M_ageArray[x, , yr.ind]
    Wt_at_Age <- Wt_age[x, , yr.ind]
    Mat_at_Age <- Mat_age[x, , yr.ind]
    V_at_Age <- V[x, , yr.ind]
  }
  else {
    M_at_Age <- apply(M_ageArray[x, , yr.ind], 1, mean)
    Wt_at_Age <- apply(Wt_age[x, , yr.ind], 1, mean)
    Mat_at_Age <- apply(Mat_age[x, , yr.ind], 1, mean)
    V_at_Age <- apply(V[x, , yr.ind], 1, mean)
  }
  boundsF <- c(1e-08, 3)
  doopt <- optimise(MSYCalcs, log(boundsF), M_at_Age, Wt_at_Age, 
                    Mat_at_Age, V_at_Age, maxage, R0x = R0[x], SRrelx = SRrel[x], Egg0 = SSBPR0[x],
                    hx = hs[x], opt = 1, plusgroup = plusgroup)
  MSYs <- MSYCalcs(doopt$minimum, M_at_Age, Wt_at_Age, Mat_at_Age, 
                   V_at_Age, maxage, R0x = R0[x], SRrelx = SRrel[x], Egg0 = SSBPR0[x], hx = hs[x], 
                   opt = 2, plusgroup = plusgroup)
  return(MSYs)
}
