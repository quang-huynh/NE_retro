
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
    rbind %>% do.call(lapply(Misc, get_n_models, y = x)) %>% structure(dimnames = list(1:length(Misc), Yr))
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
          apply(xx, 2, quantile, probs = probs, na.rm = TRUE) %>% 
            t() %>% as.data.frame() %>% mutate(Year = as.numeric(rownames(.)), MP = x)
        })
      return(output)
    })
    
    res2 <- lapply(c("F_FMSY", "B_BMSY", "rho"), function(x) {
      do.call(rbind, lapply(res, getElement, x)) %>% mutate(OM = paste0("NR", i))
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
plot_OM <- function(MSE, ref_pt = c("rel", "abs"), ref = NULL,
                    MPs = c("cod_M02", "cod_M02_ra", "cod_MRAMP", "cod_MRAMP_ra", "cod_ma")) {
  
  ref_pt <- match.arg(ref_pt)
  if(ref_pt == "rel") {
    names_list <- c("F_FMSY", "B_BMSY", "Relative Catch")
  } else {
    names_list <- c("F", "SSB", "Catch")
  }
  
  Map_fn <- function(mse, i, ref, MPs) {
    MP_ind <- match(MPs, mse@MPs)
    
    SSB_hist <- apply(mse@SSB_hist, c(1, 3), sum) %>% array(dim = c(mse@nsim, mse@nyears, length(MP_ind))) %>%
      aperm(c(1, 3, 2))
    FM_hist <- apply(mse@FM_hist, c(1, 3), max) %>% array(dim = c(mse@nsim, mse@nyears, length(MP_ind))) %>%
      aperm(c(1, 3, 2))
    CB_hist <- apply(mse@CB_hist, c(1, 3), sum) %>% array(dim = c(mse@nsim, mse@nyears, length(MP_ind))) %>%
      aperm(c(1, 3, 2))
    if(ref_pt == "rel") {
      if(missing(ref)) {
        MSE_out <- list(mse@F_FMSY[, MP_ind, , drop = FALSE], 
                        mse@B_BMSY[, MP_ind, , drop = FALSE], 
                        mse@C[, MP_ind, , drop = FALSE]/mse@OM$RefY)
      } else {
        MSE_out <- list(mse@FM[, MP_ind, , drop = FALSE]/ref[, 1], 
                        mse@SSB[, MP_ind, , drop = FALSE]/ref[, 2], 
                        mse@C[, MP_ind, , drop = FALSE]/ref[, 5])
      }
        
    } else {
      MSE_out <- list(abind::abind(FM_hist, mse@FM[, MP_ind, , drop = FALSE], along = 3), 
                      abind::abind(SSB_hist, mse@SSB[, MP_ind, , drop = FALSE], along = 3),
                      abind::abind(CB_hist, mse@C[, MP_ind, , drop = FALSE], along = 3))
    }
    
    quants <- lapply(MSE_out,
                     function(x) {
                       structure(x, dimnames = list(c(1:mse@nsim), mse@MPs[MP_ind], 
                                                    c(1:(mse@nyears + mse@proyears)))) %>%
                         reshape2::melt(varnames = c("quant", "MP", "Year")) %>%
                         group_by(MP, Year) %>% 
                         summarise(med = median(value), low = quantile(value, 0.25), 
                                   high = quantile(value, 0.75)) %>% mutate(OM = paste0("NR", i))
                     }) %>%
      structure(names = names_list)
    return(quants)
  }
  
  if(length(ref) > 0) {
    res <- Map(Map_fn, mse = MSE, i = 1:length(MSE), ref = ref, MoreArgs = list(MPs = MPs))
  } else {
    res <- Map(Map_fn, mse = MSE, i = 1:length(MSE), MoreArgs = list(MPs = MPs))
  }
  res2 <- lapply(names_list, function(x) {
    do.call(rbind, lapply(res, getElement, x))
  }) %>% structure(names = names_list)
  return(res2)
}

#' @export
plot_Index <- function(MSE, ind.names, MPs, probs = c(0.25, 0.5, 0.75)) {
  
  Map_fn <- function(mse, i, MPs) {
    MP_ind <- match(MPs, mse@MPs)
    
    DataList <- mse@Misc$Data[MP_ind]
    n_index <- DataList[[1]]@AddInd %>% dim() %>% getElement(2)
    MSE_out <- lapply(1:n_index, function(x) {
      rbind %>% do.call(lapply(MP_ind, function(xx) {
        apply(mse@Misc$Data[[xx]]@AddInd[, x, -c(1:mse@nyears)], 2, quantile, na.rm = TRUE, probs = probs) %>% 
          structure(dimnames = list(rownames(.), mse@nyears + 1:ncol(.))) %>% reshape2::melt(varnames = c("quant", "Year")) %>%
          cbind(Index = ind.names[x], MP = mse@MPs[xx])
      }))
    })
    return(do.call(rbind, MSE_out) %>% cbind(OM = paste0("NR", i)))
  }
  browser()
  
  res <- do.call(rbind, Map(Map_fn, mse = MSE, i = 1:length(MSE), MoreArgs = list(MPs = MPs))) %>%
    reshape2::dcast(list("quant", "value"))
  res2 <- lapply(names_list, function(x) {
    do.call(rbind, lapply(res, getElement, x))
  }) %>% structure(names = names_list)
  return(res2)
}
