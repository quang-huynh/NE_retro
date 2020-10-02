


#' @export
get_indicators <- function(i, MSE_list, MPs = NULL, time_stanza = NULL, ind_type = NULL, ind_interval = 5, Cbias = rep(1, 3),
                           s_CAA_hist = NULL, mah_ind = NULL, Year_vec = NULL, OM_names = NULL, mat_age = NULL) {
  MSE <- MSE_list[[i]]
  Cbias_i <- Cbias[i]
  if(is.null(MPs)) MPs <- MSE@MPs
  out <- lapply(match(MPs, MSE@MPs), get_indicators_fn, MSE = MSE, ind_interval = ind_interval, time_stanza = time_stanza, ind_type = ind_type,
                Cbias = Cbias_i, s_CAA_hist = s_CAA_hist, mah_ind = mah_ind, Year_vec = Year_vec, mat_age = mat_age)
  
  ind <- lapply(out, function(x) x[[1]] %>% dplyr::mutate(OM = OM_names[i], Year = Time + MSE@OM$CurrentYr[1] - MSE@nyears, Time = NULL)) %>%
    do.call(rbind, .) %>% mutate(OM = factor(OM, levels = OM_names))
  if(!is.null(mah_ind)) {
    mah <- do.call(rbind, lapply(out, function(x) x[[2]] %>% dplyr::mutate(OM = OM_names[i], Year = Time + MSE@OM$CurrentYr[1] - MSE@nyears, Time = NULL)))
  } else {
    mah <- NULL
  }
  return(list(ind = ind, mah = mah))
}

get_indicators_fn <- function(ii, MSE, ind_interval = 5, time_stanza = NULL, ind_type = NULL, Cbias = 1, s_CAA_hist, mah_ind, Year_vec = NULL, mat_age = mat_age) {
  Data <- MSE@Misc$Data[[ii]]
  if(is.null(Year_vec)) {
    Year_vec <- vapply(Data@Misc[[1]]$diagnostic, getElement, numeric(1), "Year")
    if(length(Year_vec) == 0) Year_vec <- vapply(MSE@Misc$Data[[1]]@Misc[[1]]$diagnostic, getElement, numeric(1), "Year")
  }
  
  Cout <- getinds(Data, Year_vec = Year_vec, res = ind_interval, tsd = c("Cat", "Cat"), stat = c("slp", "mu"), Cbias = Cbias)
  nindex <- Data@AddInd %>% dim() %>% getElement(2)
  
  Iout <- list()
  for(i in 1:nindex) {
    Data@Ind <- Data@AddInd[, i, ]
    Iout[[i]] <- getinds(Data, Year_vec = Year_vec, res = ind_interval, tsd = c("Ind", "Ind"), stat = c("slp", "mu"))
    dimnames(Iout[[i]])[[1]] <- paste0("Ind_", i, c("_slp", "_mu"))
  }
  
  rho <- sapply(Data@Misc, function(x) if(is.matrix(x$rho)) x$rho[1, ] else x$rho)
  if(is.null(rho[[1]])) {
    rho <- NULL
  } else {
    rho <- array(rho, c(dim(rho), 1)) %>% aperm(c(3, 1, 2)) %>% structure(dimnames = list("SSB_rho", Year_vec, 1:MSE@nsim))
  }
  s_CAA <- lapply(Data@Misc, getElement, "s_CAA")
  
  if(is.null(s_CAA[[1]])) {
    MAout <- NULL 
  } else {
    MAout <- PMout <- array(NA, c(dim(s_CAA[[1]])[3]*2, length(Year_vec), length(s_CAA)))
    
    for(i in 1:dim(s_CAA[[1]])[3]) { # loop over index
      mean_age <- apply(s_CAA_hist[, , i], 1, function(xx) weighted.mean(1:length(xx), xx, na.rm = TRUE)) %>% 
        matrix(length(s_CAA), MSE@nyears, byrow = TRUE)
      mean_age2 <- sapply(s_CAA, function(x) apply(x[, , i], 1, function(xx) weighted.mean(1:length(xx), xx, na.rm = TRUE))) %>% t()
      mean_age <- abind::abind(mean_age, mean_age2, along = 2)
      
      Pmat <- apply(s_CAA_hist[, , i], 1, function(xx) weighted.mean(mat_age, xx, na.rm = TRUE)) %>% 
        matrix(length(s_CAA), MSE@nyears, byrow = TRUE)
      Pmat2 <- sapply(s_CAA, function(x) apply(x[, , i], 1, function(xx) weighted.mean(mat_age, xx, na.rm = TRUE))) %>% t()
      Pmat <- abind::abind(Pmat, Pmat2, along = 2)
      
      for(j in 1:length(Year_vec)) {
        MAout[2*(i-1)+1, j, ] <- vapply(1:nrow(mean_age), MSEtool:::slp, numeric(1), mat = mean_age, ind = seq(Year_vec[j] - ind_interval + 1, Year_vec[j]))
        MAout[2*(i-1)+2, j, ] <- vapply(1:nrow(mean_age), MSEtool:::mu, numeric(1), mat = mean_age, ind = seq(Year_vec[j] - ind_interval + 1, Year_vec[j]))
        PMout[2*(i-1)+1, j, ] <- vapply(1:nrow(Pmat), MSEtool:::slp, numeric(1), mat = Pmat, ind = seq(Year_vec[j] - ind_interval + 1, Year_vec[j]))
        PMout[2*(i-1)+2, j, ] <- vapply(1:nrow(Pmat), MSEtool:::mu, numeric(1), mat = Pmat, ind = seq(Year_vec[j] - ind_interval + 1, Year_vec[j]))
      }
    }
    dimnames(MAout) <- list(paste0("MAge_", rep(1:dim(s_CAA[[1]])[3], each = 2), "_", rep(c("slp", "mu"), nindex)), Year_vec, 1:length(s_CAA))
    dimnames(PMout) <- list(paste0("PMat_", rep(1:dim(s_CAA[[1]])[3], each = 2), "_", rep(c("slp", "mu"), nindex)), Year_vec, 1:length(s_CAA))
    
  }
  
  args <- c(list(rho, Cout, MAout, PMout), Iout, list(along = 1))
  out <- do.call(abind::abind, args)
  
  if(is.null(time_stanza)) time_stanza <- seq_len(dim(out)[2])
  if(is.null(ind_type)) ind_type <- seq_len(dim(out)[1])
  out <- out[ind_type, time_stanza, , drop = FALSE]
  ind_df <- out %>% reshape2::melt(varnames = c("Ind", "Time", "Sim")) %>% dplyr::mutate(MP = MSE@MPs[ii])
  
  if(!is.null(mah_ind)) {
    out <- out[match(mah_ind, dimnames(out)[[1]]), -1, ]
    out[is.na(out)] <- 0.01
    
    mah <- MSEtool::Probs(out, out)[[1]][1, , ] %>% structure(dimnames = dimnames(out)[-1]) %>% 
      reshape2::melt(varnames = c("Time", "Sim")) %>% dplyr::mutate(MP = MSE@MPs[ii])
  } else {
    mah <- NULL
  }
  #mcont <- mah_contribution(out) %>% structure(dimnames = dimnames(out)[-1]) %>% 
  #  reshape2::melt(varnames = c("Ind", "Time", "Sim")) %>% dplyr::mutate(MP = MSE@MPs[ii])
  #return(list(ind_df = ind_df, mah = mah, mcont = mcont))
  return(list(ind_df = ind_df, mah = mah))
}

getinds <- function(Data, Year_vec, res, tsd, stat, Cbias) {
  out <- array(NA, c(length(tsd), length(Year_vec), nrow(Data@Cat)))
  for(i in 1:length(tsd)) {
    dat <- slot(Data, tsd[i])
    if(tsd[i] == "Cat") dat <- dat/Cbias
    for(j in 1:length(Year_vec)) {
      if(stat[i] == "slp") out[i, j, ] <- vapply(1:nrow(dat), MSEtool:::slp, numeric(1), mat = dat, ind = seq(Year_vec[j] - res + 1, Year_vec[j]))
      if(stat[i] == "mu") out[i, j, ] <- vapply(1:nrow(dat), MSEtool:::mu, numeric(1), mat = dat, ind = seq(Year_vec[j] - res + 1, Year_vec[j]))
    }
  }
  dimnames(out) <- list(paste(tsd, stat, sep = "_"), Year_vec, 1:nrow(dat))
  return(out)
}


#mah_contribution <- function(indPPD){
#  ntsd<-dim(indPPD)[1]
#  np<-dim(indPPD)[2]
#  nsim<-dim(indPPD)[3]
#  PRB<-array(NA,c(2,np))  # False Positive, True Positive
#  mah<-array(NA,c(2,np,nsim))
#  
#  mcont <- array(NA, dim(indPPD)) %>% structure(dimnames = dimnames(indPPD))
#  
#  for(pp in 1:np){
#    
#    keep<-array(TRUE,c(ntsd,pp))
#    ni<-sum(keep)
#    keepind<-as.matrix(expand.grid(1:ntsd,1:pp,1:nsim))[rep(as.vector(keep),nsim),]
#    ind3PPD<-t(matrix(indPPD[keepind],nrow=ni))
#    
#    covr <- cov(ind3PPD)
#    mu<-apply(ind3PPD,2,median)
#    
#    mcont[, pp, ] <- MSEtool:::mahalanobis_contribution(ind3PPD,mu,covr) %>% t()
#  }
#  return(mcont)
#}