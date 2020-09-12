
#' @export
calc_refpt <- function(sra, M = c("firstyear", "true"), med_rec = 1) {
  M <- match.arg(M)
  if(M == "firstyear") {
    M_ageArray <- array(sra@OM@cpars$M_ageArray[1,1,1], dim(sra@OM@cpars$Wt_age))
  } else {
    M_ageArray <- sra@OM@cpars$M_ageArray
  }
  
  ref_out <- lapply(1:sra@OM@nsim, optMSY, 
                    M_ageArray = M_ageArray, Wt_age = sra@OM@cpars$Wt_age, Mat_age = sra@OM@cpars$Mat_age,
                    V = sra@OM@cpars$V, maxage = sra@OM@maxage, 
                    R0 = sra@OM@cpars$R0, SRrel = rep(sra@OM@SRrel, sra@OM@nsim), 
                    SSBPR0 = vapply(sra@Misc, getElement, numeric(1), "EPR0_SR"),
                    hs = sra@OM@cpars$h, yr.ind = sra@OM@nyears, plusgroup = 1)
  ref_pt <- lapply(c("F", "SB", "Yield"), function(x) vapply(ref_out, getElement, numeric(1), x))
  
  OY <- MSYCalcs(log(0.75 * ref_pt[[1]][1]), M_at_Age = M_ageArray[1, , sra@OM@nyears], 
                 Wt_at_Age = sra@OM@cpars$Wt_age[1, , sra@OM@nyears], 
                 Mat_at_Age = sra@OM@cpars$Mat_age[1, , sra@OM@nyears], 
                 V_at_Age = sra@OM@cpars$V[1, , sra@OM@nyears],
                 maxage = sra@OM@maxage, 
                 R0x = sra@OM@cpars$R0[1], SRrelx = sra@OM@SRrel, Egg0 = sra@Misc[[1]]$EPR0_SR,
                 hx = sra@OM@cpars$h[1], opt = 2, plusgroup = 1)
  
  eq_fn <- function(logF) {
    Fout <- MSYCalcs(logF, M_at_Age = M_ageArray[1, , sra@OM@nyears], 
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
    uniroot(SPR_opt, log(c(1e-8, 3)), SPR_ratio = 0.4, wt = sra@OM@cpars$Wt_age[x, , sra@OM@nyears], 
            mat = sra@OM@cpars$Mat_age[x, , sra@OM@nyears], M = M_ageArray[x, , sra@OM@nyears], 
            V = sra@OM@cpars$V[x, , sra@OM@nyears], 
            max_age = sra@OM@maxage, plusgroup = TRUE, R = 1, opt = TRUE)[[1]] %>% exp()
  }, numeric(1))
  
  # SSBPR40%
  ref_pt$SSBPRproxy <- vapply(1:sra@OM@nsim, function(x) {
    SPR_opt(log(ref_pt$Fproxy[x]), wt = sra@OM@cpars$Wt_age[x, , sra@OM@nyears], 
            mat = sra@OM@cpars$Mat_age[x, , sra@OM@nyears], M = M_ageArray[x, , sra@OM@nyears], 
            V = sra@OM@cpars$V[x, , sra@OM@nyears], 
            max_age = sra@OM@maxage, plusgroup = TRUE, R = med_rec, opt = FALSE)[[2]]
  }, numeric(1))
  
  ref_pt$YPRproxy <- vapply(1:sra@OM@nsim, function(x) {
    YPR_opt(log(ref_pt$Fproxy[x]), wt = sra@OM@cpars$Wt_age[x, , sra@OM@nyears], 
            M = M_ageArray[x, , sra@OM@nyears], V = sra@OM@cpars$V[x, , sra@OM@nyears], 
            max_age = sra@OM@maxage, plusgroup = TRUE, R = med_rec, opt = FALSE)
  }, numeric(1))
  
  ref_pt$YPR75proxy <- vapply(1:sra@OM@nsim, function(x) {
    YPR_opt(log(0.75 * ref_pt$Fproxy[x]), wt = sra@OM@cpars$Wt_age[x, , sra@OM@nyears], 
            M = M_ageArray[x, , sra@OM@nyears], V = sra@OM@cpars$V[x, , sra@OM@nyears], 
            max_age = sra@OM@maxage, plusgroup = TRUE, R = med_rec, opt = FALSE)
  }, numeric(1))
  
  out <- do.call(cbind, ref_pt) %>% cbind(rep(OY["Yield"], length(ref_pt[[1]]))) %>%
    cbind(rep(Fcrash, length(ref_pt[[1]]))) %>%
    structure(dimnames = list(NULL, c("FMSY", "SSBMSY", "MSY", "Fproxy", "SSBproxy", "MSYproxy", 
                                      "OYproxy", "OY", "Fcrash")))
  
  return(out)
}

#' @export
EM_stock_status <- function(MSE, MP = NULL, model_index = 1) {
  ind <- match(MP, MSE@MPs)
  Misc <- MSE@Misc$Data[[ind]]@Misc
  
  out_names <- c("F_FMSY", "B_BMSY", "rho")
  Yr <- vapply(Misc[[1]]$diagnostic, getElement, numeric(1), "Year")
  
  get_n_models <- function(xx, y) { ## xx is entry in Misc list
    metric <- getElement(xx, y)
    if(is.matrix(metric) && nrow(metric) > 1) metric <- metric[model_index, ]
    return(metric)
  }
  
  out <- lapply(out_names, function(x) {
    metric <- lapply(Misc, get_n_models, y = x)
    if(!is.null(metric[[1]])) {
      do.call(rbind, metric) %>% structure(dimnames = list(1:length(Misc), Yr))
    } else return(NULL)
  }) %>% structure(names = out_names)
  
  return(out)
}

#' @export
plot_EM <- function(MSE, MPs = c("cod_M02", "cod_M02_ra", "cod_MRAMP", "cod_MRAMP_ra"),
                    model_index = 1, probs = c(0.25, 0.5, 0.75)) {
  
  Map_fn <- function(mse, i) {
    res <- lapply(MPs, function(x) {
      output <- mfSCA::EM_stock_status(mse, x, model_index = model_index) %>% 
        lapply(function(xx) {
          if(is.null(xx)) {
            return(NULL)
          } else {
            apply(xx, 2, quantile, probs = probs, na.rm = TRUE) %>% 
              t() %>% as.data.frame() %>% mutate(Year = as.numeric(rownames(.)) + mse@OM$CurrentYr[1] - mse@nyears, MP = x)
          }
        })
      return(output)
    })
    
    res2 <- lapply(c("F_FMSY", "B_BMSY", "rho"), function(x) {
      out <- do.call(rbind, lapply(res, getElement, x))
      if(is.null(out)) {
        return(NULL)
      } else {
        return(mutate(out, OM = paste0("NR", i)))
      }
    }) %>% structure(names = c("F_FMSY", "B_BMSY", "rho"))
    
    return(res2)
  }
  out <- Map(Map_fn, mse = MSE, i = 1:length(MSE))
  
  out2 <- lapply(c("F_FMSY", "B_BMSY", "rho"), function(x) {
    do.call(rbind, lapply(out, getElement, x))
  }) %>% structure(names = c("F_FMSY", "B_BMSY", "rho"))
  return(out2)
}

#' @export
plot_OM <- function(MSE, MPs, qlow = 0.25, qhigh = 0.75, yfilter = NULL) {
  
  names_list <- c("F", "SSB", "Catch")
  
  Map_fn <- function(mse, i, ref, MPs) {
    MP_ind <- match(MPs, mse@MPs)
    
    SSB_hist <- apply(mse@SSB_hist, c(1, 3), sum) %>% array(dim = c(mse@nsim, mse@nyears, length(MP_ind))) %>%
      aperm(c(1, 3, 2))
    FM_hist <- apply(mse@FM_hist, c(1, 3), max) %>% array(dim = c(mse@nsim, mse@nyears, length(MP_ind))) %>%
      aperm(c(1, 3, 2))
    CB_hist <- apply(mse@CB_hist, c(1, 3), sum) %>% array(dim = c(mse@nsim, mse@nyears, length(MP_ind))) %>%
      aperm(c(1, 3, 2))
    
    MSE_out <- list(abind::abind(FM_hist, mse@FM[, MP_ind, , drop = FALSE], along = 3), 
                    abind::abind(SSB_hist, mse@SSB[, MP_ind, , drop = FALSE], along = 3),
                    abind::abind(CB_hist, mse@C[, MP_ind, , drop = FALSE], along = 3))
    
    quants <- lapply(MSE_out,
                     function(x, mse) {
                       df_out <- 
                         structure(x, dimnames = list(c(1:mse@nsim), mse@MPs[MP_ind], 
                                                      c(1:(mse@nyears + mse@proyears)))) %>%
                         reshape2::melt(varnames = c("quant", "MP", "Year")) %>%
                         group_by(MP, Year) %>% 
                         summarise(med = median(value), low = quantile(value, qlow), 
                                   high = quantile(value, qhigh)) %>% mutate(OM = paste0("NR", i))
                       df_out$Year <- df_out$Year + mse@OM$CurrentYr[1] - mse@nyears
                       if(!is.null(yfilter)) df_out <- filter(df_out, Year >= mse@OM$CurrentYr[1] - yfilter)
                       return(df_out)
                     }, mse = mse) %>%
      structure(names = names_list)
    return(quants)
  }
  
  res <- Map(Map_fn, mse = MSE, i = 1:length(MSE), MoreArgs = list(MPs = MPs))
  
  #if(length(ref) > 0) {
  #  res <- Map(Map_fn, mse = MSE, i = 1:length(MSE), ref = ref, MoreArgs = list(MPs = MPs))
  #} else {
  #  res <- Map(Map_fn, mse = MSE, i = 1:length(MSE), MoreArgs = list(MPs = MPs))
  #}
  res2 <- lapply(names_list, function(x) {
    do.call(rbind, lapply(res, getElement, x))
  }) %>% structure(names = names_list)
  return(res2)
}

#' @export
plot_Index <- function(MSE, ind.names, MPs, probs = c(0.25, 0.5, 0.75), yfilter = NULL) {
  
  Map_fn <- function(mse, i, MPs) {
    MP_ind <- match(MPs, mse@MPs)
    DataList <- mse@Misc$Data[MP_ind]
    n_index <- DataList[[1]]@AddInd %>% dim() %>% getElement(2)
    
    out <- lapply(1:length(DataList), function(x) {
      rr <- lapply(1:n_index, function(xx) {
        index_quantiles <- apply(DataList[[x]]@AddInd[, xx, ], 2, quantile, na.rm = TRUE, probs = probs) %>% 
          t() %>% as.data.frame()
        index_quantiles$Year <- 1:nrow(index_quantiles) + mse@OM$CurrentYr[1] - mse@nyears
        index_quantiles$Index <- ind.names[xx]
        if(!is.null(yfilter)) index_quantiles <- filter(index_quantiles, Year >= mse@OM$CurrentYr[1] - yfilter)
        return(index_quantiles)
      })
      rr <- do.call(rbind, rr)
      rr$MP <- MPs[x]
      return(rr)
    })
    out <- do.call(rbind, out)
    out$OM <- paste0("NR", i)
    return(out)
  }
  
  res <- do.call(rbind, Map(Map_fn, mse = MSE, i = 1:length(MSE), MoreArgs = list(MPs = MPs)))
  res <- split(res, f = res$Index)
  return(res)
}
