
#' @export
create_MP <- function(args, args2 = NULL, rho_adjust = FALSE, calc_rho = TRUE) {
  force(args)
  force(args2)
  force(rho_adjust)
  force(calc_rho)
  
  HCR_MSY_mf_SCA <- function(Assessment, reps = 1, MSY_frac = 0.75, rho) {
    if(missing(rho)) rho <- 0
    Rec <- new("Rec")
    if(Assessment@conv) {
      y <- Assessment@obj$env$data$n_y + 1
      
      N_at_age <- Assessment@N_at_age[y, ]
      V <- Assessment@Selectivity[y-1, ]
      wt <- Assessment@obj$env$data$wt[y, ]
      M <- Assessment@obj$env$data$M[y-1, ]
      
      FF <- MSY_frac * Assessment@FMSY * V
      Z <- FF + M
      VB <- sum(V/Z * N_at_age * (1 - exp(-Z)) * wt)
      if(!is.na(rho)) VB <- VB/(1 + rho)
      ABC <- rep(MSY_frac * Assessment@FMSY * VB, reps)
    } else {
      ABC <- rep(NA, reps)
    }
    Rec@TAC <- TACfilter(ABC)
    return(Rec)
  }
  
  if(!is.null(args2)) { # Model average in lieu of rho adjustment
    
    myMP <- function(x, Data, reps) {
      models <- lapply(list(args, args2), function(xx) mf_SCA(x, Data, args = xx))
      conv <- vapply(models, getElement, logical(1), "conv")
      
      Rec <- new("Rec")
      if(all(conv)) { # Model average only if both models converge
        ABC <- lapply(models, function(xx) HCR_MSY_mf_SCA(xx, reps, rho = NA)@TAC)
        Rec@TAC <- do.call(c, ABC) %>% mean(na.rm = TRUE) %>% rep(reps) %>% TACfilter()
      } else {
        Rec@TAC <- NA_real_ %>% rep(reps)
      }
      
      Assess_diag <- MSEtool:::Assess_diagnostic(x, Data, models[[1]], include_assessment = FALSE)
      Assess_diag2 <- MSEtool:::Assess_diagnostic(x, Data, models[[2]], include_assessment = FALSE)
      
      Rec@Misc <- c(Assess_diag, models[[1]]@info) # Note that identical(models[[1]]@info, model[[2]]@info) = TRUE
      Rec@Misc$diagnostic2 <- Assess_diag2
      
      Rec@Misc$F_FMSY <- cbind(Data@Misc[[x]]$F_FMSY, 
                               vapply(models, function(xx) ifelse(xx@conv, xx@F_FMSY[length(xx@F_FMSY)], NA_real_), numeric(1)))
      Rec@Misc$B_BMSY <- cbind(Data@Misc[[x]]$B_BMSY,
                               vapply(models, function(xx) ifelse(xx@conv, xx@SSB_SSBMSY[length(xx@SSB_SSBMSY)], NA_real_), numeric(1)))
      
      if(calc_rho) {
        ret <- lapply(models, retrospective_mf_SCA, nyr = 7)
        rho <- vapply(ret, function(xx) summary(xx)[, 1]["Spawning biomass"], numeric(1))
        Rec@Misc$rho <- cbind(Data@Misc[[x]]$rho, rho)
      }
      
      return(Rec)
    }
    
  } else { # Assessment models with or without rho adjustment
    
    myMP <- function(x, Data, reps) {
      do_Assessment <- mf_SCA(x, Data, args = args)
      
      if(rho_adjust || calc_rho) {
        ret <- retrospective_mf_SCA(do_Assessment, nyr = 7)
        rho <- summary(ret)[, 1]["Spawning biomass"]
      } else {
        rho <- NA
      }
      
      Rec <- HCR_MSY_mf_SCA(do_Assessment, reps, rho = ifelse(rho_adjust, rho, NA))
      Assess_diag <- MSEtool:::Assess_diagnostic(x, Data, do_Assessment, include_assessment = FALSE)
      Rec@Misc <- c(Assess_diag, do_Assessment@info)
      
      Rec@Misc$F_FMSY <- c(Data@Misc[[x]]$F_FMSY, 
                           ifelse(do_Assessment@conv, do_Assessment@F_FMSY[length(do_Assessment@F_FMSY)], NA_real_))
      Rec@Misc$B_BMSY <- c(Data@Misc[[x]]$B_BMSY,
                           ifelse(do_Assessment@conv, do_Assessment@SSB_SSBMSY[length(do_Assessment@SSB_SSBMSY)], NA_real_))
      Rec@Misc$rho <- c(Data@Misc[[x]]$rho, rho)
      return(Rec)
    }
    
  }
  class(myMP) <- "MP"
  return(myMP)
}


#' @export
generate_Eff_MP <- function(relF, FMSY, terminalF) {
  force(relF)
  force(FMSY)
  force(terminalF)
  MP <- function(x, Data, reps) {
    Rec <- new("Rec")
    Ftarget <- relF * FMSY
    Rec@Effort <- rep(Ftarget/terminalF, reps)
    
    return(Rec)
  }
  class(MP) <- "MP"
  return(MP)
}



#' @export
generate_Eff_MP_with_EM <- function(relF, FMSY, terminalF, args, args2 = NULL) {
  force(args)
  force(args2)
  force(relF)
  force(FMSY)
  force(terminalF)
  MP <- function(x, Data, reps) {
    Rec <- new("Rec")
    Ftarget <- relF * FMSY
    Rec@Effort <- rep(Ftarget/terminalF, reps)
    
    if(!is.null(args2)) {
      models <- lapply(list(args, args2), function(xx) mf_SCA(x, Data, args = xx))
      conv <- vapply(models, getElement, logical(1), "conv")
      
      Assess_diag <- MSEtool:::Assess_diagnostic(x, Data, models[[1]], include_assessment = FALSE)
      Assess_diag2 <- MSEtool:::Assess_diagnostic(x, Data, models[[2]], include_assessment = FALSE)
      
      Rec@Misc <- c(Assess_diag, models[[1]]@info) # Note that identical(models[[1]]@info, model[[2]]@info) = TRUE
      Rec@Misc$diagnostic2 <- Assess_diag2
      
      ret <- lapply(models, retrospective_mf_SCA, nyr = 7)
      rho <- vapply(ret, function(xx) summary(xx)[, 1]["Spawning biomass"], numeric(1))
      Rec@Misc$rho <- cbind(Data@Misc[[x]]$rho, rho)
    } else {
      do_Assessment <- mf_SCA(x, Data, args = args)
      
      ret <- retrospective_mf_SCA(do_Assessment, nyr = 7)
      rho <- summary(ret)[, 1]["Spawning biomass"]
      
      Assess_diag <- MSEtool:::Assess_diagnostic(x, Data, do_Assessment, include_assessment = FALSE)
      Rec@Misc <- c(Assess_diag, do_Assessment@info)
      Rec@Misc$rho <- c(Data@Misc[[x]]$rho, rho)
    }
    return(Rec)
  }
  class(MP) <- "MP"
  return(MP)
}

