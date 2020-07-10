
#' @import MSEtool dplyr

# Multi-fleet SCA similar to ASAP
#' @export
#' @importFrom abind abind
mf_SCA <- function(x = 1, Data, args = list(), R_proxy = expression(median(SRA_out$report$R)),
                   ymax = expression(ncol(Data@Cat))) {  
  data <- args$data # arguments from SRA scoping (historical data only)
  ymax <- eval(ymax)
  data$nyears <- ymax
  set.seed(x * ymax + data$nsurvey)
  
  scalls <- lapply(sys.calls(), as.character) %>% unlist()
  inside_MSE <- length(grep("runMSE_int", scalls)) > 0
  
  if(data$nyears > Data@LHYear) {
    data$CAL <- array(0, c(data$nyears, Data@MaxAge, data$nfleet))
    data$s_CAL <- array(0, c(data$nyears, Data@MaxAge, data$nsurvey))
    data$MS <- data$Ehist <- matrix(NA, data$nyears, data$nfleet)
    
    update_Index <- matrix(Data@AddInd[x, , (Data@LHYear+1):data$nyears]/args$Ibias, 
                           data$nyears - Data@LHYear, ncol(args$data$Index), byrow = TRUE)
    data$Index <- rbind(args$data$Index, update_Index)
    data$I_sd <- rbind(args$data$I_sd, 
                       sdconv(1, Data@CV_AddInd[x, , (Data@LHYear+1):data$nyears]) %>%
                         matrix(data$nyears - Data@LHYear, ncol(args$data$Index), byrow = TRUE))
    
    data$sel_block <- rbind(args$data$sel_block, 
                             matrix(args$data$sel_block[nrow(args$data$sel_block), , drop = FALSE], 
                                    data$nyears - Data@LHYear, data$nfleet, byrow = TRUE))
    
    # Sample age comps from survey
    if(inside_MSE) {
      if(is.null(Data@Misc[[x]]$Nind)) {
        Nsearch <- lapply(1:length(sys.calls()), function(xx) try(get("N_P", envir = sys.frames()[[xx]], inherits = FALSE)[x, , , ] %>%
                                                                    apply(1:2, sum) %>% t(), silent = TRUE))
        Nind <- vapply(Nsearch, function(xx) !is.character(xx), logical(1)) %>% which()
        N <- Nsearch[[Nind]]
      } else {
        Nind <- Data@Misc[[x]]$Nind
        N <- get("N_P", envir = sys.frames()[[Nind]], inherits = FALSE)[x, , , ] %>% apply(1:2, sum) %>% t()
      }
    }
    if(inside_MSE && is.null(Data@Misc[[x]]$s_CAA)) { # Pass out s_CAA through Assessment@info
      s_CAA <- array(NA_real_, c(dim(N), data$nsurvey))
    } else {
      s_CAA <- Data@Misc[[x]]$s_CAA
    }
    s_CAA_n <- colSums(args$data$s_CAA[Data@LHYear, , ])
    
    styear <- which(apply(s_CAA, 1, sum, na.rm = TRUE) == 0)[1]
    endyear <- data$nyears - Data@LHYear
    
    if(inside_MSE && data$nyears >= styear + Data@LHYear) {
      for(y in styear:endyear) {
        s_CAA[y, , ] <- cbind %>% 
          do.call(lapply(1:data$nsurvey, function(i) rmultinom(1, s_CAA_n[i], args$AddIndV[i, ] * N[y, ])))
      }
    }
    data$s_CAA <- abind::abind(args$data$s_CAA, s_CAA[1:endyear, , , drop = FALSE], along = 1)
    
    if(data$nfleet > 1) { # Sample fleet-specific fishery age comps and catches
      if(inside_MSE && is.null(Data@Misc[[x]]$CAA)) {
        CAA <- array(NA_real_, c(dim(N), data$nfleet))
      } else {
        CAA <- Data@Misc[[x]]$CAA
      }
      if(inside_MSE && is.null(Data@Misc[[x]]$Chist)) {
        Chist <- array(NA_real_, c(dim(N)[1], data$nfleet))
      } else {
        Chist <- Data@Misc[[x]]$Chist
      }
      
      if(inside_MSE && data$nyears >= styear + Data@LHYear) { # age comps
        FM <- get("FM_Pret", envir = sys.frames()[[Nind]], inherits = FALSE)[x, , , 1] %>% t() # Assume equal size areas
        Fapical <- apply(FM[1:endyear, , drop = FALSE], 1, max) # Annual apical F
        relV <- args$relF/sum(args$relF) * t(args$V)
        F_fleet <- lapply(1:data$nfleet, function(xx) outer(Fapical, relV[xx, ], )) %>% unlist() %>% 
          array(c(endyear, Data@MaxAge, data$nfleet))
        M <- get("StockPars", envir = sys.frames()[[Nind]], inherits = FALSE)$M_ageArray[x, , Data@LHYear + 1:endyear]
        Z <- apply(F_fleet, 1:2, sum) + t(M)
        CAA_n <- colSums(args$data$CAA[Data@LHYear, , ])
        
        CB <- get("CB_Pret", envir = sys.frames()[[Nind]])[x, , , ] %>% apply(1:2, sum) %>% t()
        Cerr <- get("ErrList", envir = sys.frames()[[Nind]])$Cerr[x, ]
        for(y in styear:endyear) {
          Chist[y, ] <- vapply(1:data$nfleet, function(i) sum(CB[y, ] * relV[i, ]/colSums(relV)), numeric(1))/args$Cbias * Cerr[y]
          CAA[y, , ] <- cbind %>%
            do.call(lapply(1:data$nfleet, function(i) rmultinom(1, CAA_n[i], F_fleet[y, , i] * N[y, ]/Z[y, ] * (1 - exp(-Z[y, ])))))
        }
      }
      data$Chist <- rbind(args$data$Chist, Chist[1:endyear, , drop = FALSE])
      data$CAA <- abind::abind(args$data$CAA, CAA[1:endyear, , , drop = FALSE], along = 1)
      
    } else {
      data$Chist <- rbind(args$data$Chist, matrix(Data@Cat[x, (Data@LHYear+1):data$nyears], ncol = 1)/args$Cbias)
      data$CAA <- args$data$CAA[,,1] %>% rbind(Data@CAA[x, (Data@LHYear+1):data$nyears, ]) %>% 
        array(c(data$nyears, Data@MaxAge, 1))
    }
  }
  SRA_out <- MSEtool:::SRA_scope_est(x = x, data = data, selectivity = args$selectivity, s_selectivity = args$s_selectivity, 
                                     ESS = args$ESS, LWT = args$LWT, StockPars = args$StockPars, FleetPars = args$FleetPars,
                                     ObsPars = list(Isd = rep(0.1, nrow(data$Chist))), dots = args$dots)
  
  if(!SRA_out$SD$pdHess && (any(args$selectivity == -2) | any(args$s_selectivity == -2))) {
    if(any(args$selectivity == -2)) {
      args$dots$vul_par[, args$selectivity == -2] <- MSEtool:::ilogit(SRA_out$report$vul_par[, args$selectivity == -2])
      
      map2 <- args$dots$map_vul_par[, args$selectivity == -2]
      map2[abs(SRA_out$report$vul_par[, args$selectivity == -2]) >= 7] <- NA
      args$dots$map_vul_par[, args$selectivity == -2] <- map2
    }
    if(any(args$s_selectivity == -2)) {
      args$dots$s_vul_par[, args$s_selectivity == -2] <- MSEtool:::ilogit(SRA_out$report$s_vul_par[, args$s_selectivity == -2])
      
      map2 <- args$dots$map_s_vul_par[, args$s_selectivity == -2]
      map2[abs(SRA_out$report$s_vul_par[, args$s_selectivity == -2]) >= 7] <- NA
      args$dots$map_s_vul_par[, args$s_selectivity == -2] <- map2
    }
    SRA_out <- MSEtool:::SRA_scope_est(x = x, data = data, selectivity = args$selectivity, s_selectivity = args$s_selectivity, 
                                       ESS = args$ESS, LWT = args$LWT, StockPars = args$StockPars, FleetPars = args$FleetPars,
                                       ObsPars = list(Isd = rep(0.1, nrow(data$Chist))), dots = args$dots)
  }
  
  R_proxy <- eval(R_proxy)
  ref_pt <- mf_SCA_ref_pt(SRA_out, R_proxy = R_proxy) %>% # By default, use F40% and median recruitment
    structure(names = c("FMSY", "EMSY", "VBMSY", "BMSY", "MSY"))
  
  F_report <- which.max(ref_pt$FMSY)[1] # This must be the apical F in the terminal year
  F_out <- SRA_out$report$F_at_age[, F_report]
  
  Year <- Data@Year[1:ymax]
  Yearplusone <- c(Year, max(Year) + 1)
  YearEarly <- (Year[1] - Data@MaxAge + 1):(Year[1] - 1)
  YearDev <- c(YearEarly, Year)
  YearR <- c(YearDev, max(YearDev) + 1)
  R <- c(rev(SRA_out$report$R_early), SRA_out$report$R)
  
  Dev <- structure(c(rev(SRA_out$report$log_early_rec_dev), SRA_out$report$log_rec_dev), names = YearDev)
  
  Assessment <- new("Assessment", Model = "mf_SCA", Name = Data@Name, 
                    conv = !is.character(SRA_out$SD) && SRA_out$SD$pdHess,
                    h = SRA_out$report$h, FMort = structure(F_out, names = Year),
                    B = structure(SRA_out$report$B, names = Yearplusone),
                    SSB = structure(SRA_out$report$E, names = Yearplusone),
                    R = structure(R, names = YearR),
                    N = structure(rowSums(SRA_out$report$N), names = Yearplusone),
                    N_at_age = SRA_out$report$N,
                    Selectivity = SRA_out$report$F_at_age/apply(SRA_out$report$F_at_age, 1, max),
                    Dev = Dev, Dev_type = "log-Recruitment deviations",
                    NLL = ifelse(is.character(SRA_out$opt), NA_real_, SRA_out$opt$objective),
                    obj = SRA_out$obj, opt = SRA_out$opt, SD = SRA_out$SD, TMB_report = c(SRA_out$report, ref_pt))
  
  Assessment@FMSY <- mean(ref_pt$FMSY[F_report])
  Assessment@MSY <- ref_pt$MSY
  Assessment@BMSY <- ref_pt$BMSY
  Assessment@SSBMSY <- ref_pt$EMSY
  Assessment@VBMSY <- ref_pt$VBMSY
  Assessment@F_FMSY <- structure(Assessment@FMort/Assessment@FMSY, names = Year)
  Assessment@B_BMSY <- structure(Assessment@B/Assessment@BMSY, names = Yearplusone)
  Assessment@SSB_SSBMSY <- structure(Assessment@SSB/Assessment@SSBMSY, names = Yearplusone)
  Assessment@VB <- rowSums(rbind(Assessment@Selectivity, Assessment@Selectivity[data$nyears, ]) * 
                             Assessment@N_at_age * SRA_out$obj$env$data$wt)
  Assessment@VB_VBMSY <- structure(Assessment@VB/Assessment@VBMSY, names = Yearplusone)
  
  Assessment@info <- list()
  if(exists("Nind", inherits = FALSE)) Assessment@info$Nind <- Nind
  if(exists("s_CAA", inherits = FALSE)) Assessment@info$s_CAA <- s_CAA
  if(data$nfleet > 1) {
    if(exists("CAA", inherits = FALSE)) Assessment@info$CAA <- CAA
    if(exists("Chist", inherits = FALSE)) Assessment@info$Chist <- Chist
  }
  
  if(Assessment@conv) {
    SE_Early <- ifelse(SRA_out$obj$env$data$est_early_rec_dev, 
                       sqrt(diag(SRA_out$SD$cov.fixed)[names(SRA_out$SD$par.fixed) == "log_early_rec_dev"]), NA)
    SE_Main <- ifelse(SRA_out$obj$env$data$est_rec_dev, 
                      sqrt(diag(SRA_out$SD$cov.fixed)[names(SRA_out$SD$par.fixed) == "log_rec_dev"]), NA)
    SE_Dev <- structure(c(rev(SE_Early), SE_Main), names = YearDev)
    
    first_non_zero <- which(!is.na(SE_Dev))[1]
    if(!is.na(first_non_zero) && first_non_zero > 1) {
      Dev <- Dev[-c(1:(first_non_zero - 1))]
      SE_Dev <- SE_Dev[-c(1:(first_non_zero - 1))]
      SE_Dev[is.na(SE_Dev)] <- 0
    }
    Assessment@Dev <- Dev
    Assessment@SE_Dev <- SE_Dev
  }
  return(Assessment)
}
class(mf_SCA) <- "Assess"

retrospective_mf_SCA <- function(Assessment, nyr = 5) {
  OM <- DLMtool::testOM
  OM@CurrentYr <- names(Assessment@FMort) %>% as.numeric() %>% max()
  SRA <- new("SRA", OM = testOM, mean_fit = list(obj = Assessment@obj, opt = Assessment@opt, 
                                                 SD = Assessment@SD, report = Assessment@TMB_report))
  MSEtool:::SRA_retro(SRA, nyr)
}
