
#' @export
ASAP2SRA <- function(file = "GoM_cod/GOM_COD_2019_UPDATE_M02.dat", 
                     proyears = 50, nsim = 2, interval = 3, WAA_pointer = 2, discards = FALSE,
                     condition = "catch2") {
  out <- ASAPplots::ReadASAP3DatFile(file)
  dat <- out$dat
  
  ###### OM
  OM <- new("OM", Stock = Blue_shark, Fleet = Generic_Fleet, Obs = Precise_Unbiased, Imp = Perfect_Imp)
  OM@interval <- interval
  OM@Name <- "ASAP3"
  OM@Source <- "ASAP3"
  OM@M <- OM@AC <- OM@L50 <- OM@L50_95 <- OM@D <- OM@L5 <- OM@LFS <- rep(0, 2)
  OM@proyears <- proyears
  OM@nsim <- nsim
  
  OM@nyears <- dat$n_years
  OM@CurrentYr <- dat$year1 + OM@nyears - 1
  OM@maxage <- maxage <- dat$n_ages
  
  OM@h <- dat$steepness_ini %>% rep(2)
  OM@SRrel <- 1L
  
  OM@Linf <- rep(100, 2)
  
  OM@Perr <- sqrt(log(mean(dat$recruit_cv)^2 + 1)) %>% rep(2)
  
  Year <- dat$year1:OM@CurrentYr
  
  # Natural mortality
  M_hist <- array(dat$M, dim = c(OM@nyears, OM@maxage, nsim)) %>% aperm(3:1)
  M_proj <- array(dat$M[OM@nyears, ], c(OM@maxage, OM@proyears, nsim)) %>% aperm(c(3, 1, 2))
  OM@cpars$M_ageArray <- abind::abind(M_hist, M_proj, along = 3)
  
  # Maturity
  Mat_hist <- array(dat$maturity, dim = c(OM@nyears, OM@maxage, nsim)) %>% aperm(3:1)
  Mat_proj <- array(dat$maturity[OM@nyears, ], c(OM@maxage, OM@proyears, nsim)) %>% aperm(c(3, 1, 2))
  OM@cpars$Mat_age <- abind::abind(Mat_hist, Mat_proj, along = 3)
  
  # WAA
  Wt_hist <- array(dat$WAA_mats[[WAA_pointer]], dim = c(OM@nyears, OM@maxage, nsim)) %>% aperm(3:1)
  Wt_proj <- array(dat$WAA_mats[[WAA_pointer]][OM@nyears, ], c(OM@maxage, OM@proyears, nsim)) %>% aperm(c(3, 1, 2))
  OM@cpars$Wt_age <- abind::abind(Wt_hist, Wt_proj, along = 3)
  
  ###### SRA_scope args
  nfleet <- dat$n_fleets
  SRA_args <- list()
  SRA_args$max_F <- OM@maxF <- dat$Fmax
  
  # Selectivity
  ASAP_sel <- data.frame(code = 1:3, tx = c("free", "logistic", "dome"), stringsAsFactors = FALSE)
  SRA_args$selectivity <- ASAP_sel$tx[match(dat$sel_block_option, ASAP_sel$code)]
  #SRA_args$s_selectivity <- ASAP_sel$tx[match(dat$sel_block_option, ASAP_sel$code)]
  
  ###### SRA_scope data
  SRA_data <- list()
  
  # Catch
  SRA_data$Chist <- do.call(cbind, lapply(dat$CAA_mats, function(x) x[, OM@maxage + 1]))
  
  # Discards
  if(discards) {
    Dis = do.call(cbind, lapply(dat$DAA_mats, function(x) x[, OM@maxage + 1]))
    ndis_fleet <- ncol(Dis)
    SRA_data$Chist <- cbind(SRA_data$Chist, Dis)
  }
  
  # CAA
  SRA_data$CAA <- lapply(dat$CAA_mats, function(x) x[, 1:OM@maxage]) %>% unlist() %>% array(c(OM@nyears, OM@maxage, nfleet))
  for(ff in 1:nfleet) {
    SRA_data$CAA[,,ff] <- SRA_data$CAA[,,ff] * dat$catch_Neff[, ff] / rowSums(SRA_data$CAA[,,ff])
  }
  
  # CAA discards
  if(discards) {
    DAA <- lapply(dat$DAA_mats, function(x) x[, 1:OM@maxage]) %>% unlist() %>% array(c(OM@nyears, OM@maxage, ndis_fleet))
    for(ff in 1:ndis_fleet) DAA[,,ff] <- DAA[,,ff] * dat$discard_Neff[, ff] / rowSums(DAA[,,ff])
    SRA_data$CAA <- abind::abind(SRA_data$CAA, DAA, along = 3)
  }
  
  # Index
  Index_fn <- function(x) {
    out <- rep(NA, length(Year))
    out[match(x[, 1], Year)] <- x[, 2]
    return(out)
  }
  SRA_data$Index <- do.call(cbind, lapply(dat$IAA_mats[as.logical(dat$use_index)], Index_fn))
  SRA_data$Index[SRA_data$Index <= 0] <- NA
  
  
  # I_sd - convert from CV to sd
  Isd_fn <- function(x) {
    out <- rep(NA, length(Year))
    out[match(x[, 1], Year)] <- sqrt(log(x[, 3]^2 + 1))
    return(out)
  }
  SRA_data$I_sd <- do.call(cbind, lapply(dat$IAA_mats[as.logical(dat$use_index)], Isd_fn))
  
  # s_CAA
  s_CAA_fn <- function(x) {
    out <- matrix(NA, length(Year), OM@maxage)
    out[match(x[, 1], Year), ] <- x[, -c(1:3, ncol(x))]
    out <- out * x[, ncol(x)] / rowSums(out)
    return(out)
  }
  SRA_data$s_CAA <- do.call(cbind, lapply(dat$IAA_mats[as.logical(dat$use_index_acomp)], s_CAA_fn)) %>% 
    unlist() %>% array(c(OM@nyears, OM@maxage, sum(dat$use_index_acomp)))
  
  # Index units. ASAP (1 = biomass, 2 = numbers). SRA(1 = biomass, 0 = numbers)
  SRA_data$I_units <- dat$index_units[as.logical(dat$use_index)]
  SRA_data$I_units[SRA_data$I_units == 2] <- 0
  
  #SRA_data$I_type <- rep("est", length(SRA_data$I_units))
  
  # sel_block
  SRA_data$sel_block <- matrix(do.call(cbind, dat$sel_block_assign), OM@nyears, nfleet)
  if(discards) SRA_data$sel_block <- cbind(SRA_data$sel_block, SRA_data$sel_block)
  SRA_data$nsel_block <- SRA_data$sel_block %>% as.vector() %>% unique() %>% length()
  
  # Fill in OM
  OM@h <- rep(0.99, 2)
  OM@R0 <- 1e5
  OM@AC <- rep(0, 2)
  OM@L5 <- rep(2, 2)
  OM@LFS <- rep(4, 2)
  OM@Vmaxlen <- rep(0.4, 2)
  OM@isRel <- "FALSE"
  
  OM@Linf <- rep(10, 2)
  OM@K <- rep(0.2, 2)
  OM@t0 <- rep(0, 2)
  
  OM@cpars$Len_age <- c(1:maxage) %>% array(c(OM@maxage, OM@nsim, OM@proyears + OM@nyears)) %>%
    aperm(c(2, 1, 3))
  
  # Fill in args
  SRA_args$data <- SRA_data
  SRA_args$OM <- OM
  MLAA <- SRA_args$data$length_bin <- 1:maxage
  
  SRA_args$condition <- condition
  SRA_args$ESS <- rep(1e5, 2)
  SRA_args$mean_fit <- TRUE
  SRA_args$map_log_early_rec_dev <- 1:(maxage-1) 
  SRA_args$s_selectivity <- rep("free", length(SRA_data$I_units))
  
  if(discards) {
    return(list(SRA_args = SRA_args, prop_rel = dat$prop_rel_mats))
  } else {
    return(SRA_args)
  }
}


#' @export
get_cod_M02 <- function(nsim = 2, h = 0.99) {
  args <- ASAP2SRA("ASAP/GOM_COD_2019_UPDATE_M02.dat", condition = "catch", nsim = nsim)
  args$OM@h <- rep(h, 2)
  
  sel_par <- list(c(0.05, 0.2, 0.4, 0.79, rep(0.999999, 5)),
                  c(0.05, 0.2, 0.4, rep(0.999999, 6)),
                  c(0.999999, 0.8, 0.6, 0.5, 0.4, 0.2, rep(0.0001, 3))) %>% unlist() %>% matrix(9, 3)
  m_sel_par <- matrix(NA, 9, 3)
  m_sel_par[1:4, 1] <- 1:4
  m_sel_par[1:3, 2] <- 5:7
  m_sel_par[2:6, 3] <- 8:12
  
  args$s_vul_par <- sel_par
  args$map_s_vul_par <- m_sel_par
  return(args)
}

#' @export
get_cod_MRAMP <- function(nsim = 2, h = 0.99) {
  args <- ASAP2SRA("ASAP/GOM_COD_2019_UPDATE_MRAMP.dat", condition = "catch", nsim = nsim)
  args$OM@h <- rep(h, 2)
  
  sel_par <- list(c(0.05, 0.2, 0.4, 0.79, 0.9, rep(0.999999, 4)),
                  c(0.05, 0.2, 0.4, 0.79, 0.9, rep(0.999999, 4)),
                  c(0.999999, 0.8, 0.6, 0.5, 0.4, 0.2, rep(0.0001, 3))) %>% unlist() %>% matrix(9, 3)
  m_sel_par <- matrix(NA, 9, 3)
  m_sel_par[1:5, 1] <- 1:5
  m_sel_par[1:5, 2] <- 6:10
  m_sel_par[2:6, 3] <- 11:15
  
  args$s_vul_par <- sel_par
  args$map_s_vul_par <- m_sel_par
  return(args)
}


#' @export
get_haddock <- function(nsim = 2, h = 0.99) {
  args <- ASAP2SRA("ASAP/GOM_HADDOCK_ASAP_2019_BASE_SEL9.dat", nsim = nsim)
  args$OM@h <- rep(h, 2)
  
  args$vul_par <- matrix(c(0.01, 0.1, 0.3, 0.5, 0.8, 0.9, rep(0.99999, 3)), 9, 3)
  
  map_vul_par <- matrix(NA, 9, 3)
  map_vul_par[1:6, ] <- 1:18
  args$map_vul_par <- map_vul_par
  
  sel_par <- list(c(0.2, 0.4, 0.8, rep(0.999999, 6)),
                  c(0.2, 0.4, 0.8, rep(0.999999, 6))) %>% unlist() %>% matrix(9, 2)
  
  m_sel_par <- matrix(NA, 9, 2)
  m_sel_par[1:3, 1] <- 1:3
  m_sel_par[1:3, 2] <- 4:6
  
  args$s_vul_par <- sel_par
  args$map_s_vul_par <- m_sel_par
  return(args)
}
