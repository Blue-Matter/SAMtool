
RCM_est <- function(x = 1, RCMdata, selectivity, s_selectivity, LWT = list(),
                    comp_like = c("multinomial", "lognormal", "mvlogistic", "dirmult1", "dirmult2"), prior = list(),
                    max_F = 3, integrate = FALSE, StockPars, FleetPars, mean_fit = FALSE,
                    control = list(iter.max = 2e+05, eval.max = 4e+05), inner.control = list(maxit = 1e3), 
                    start = list(), map = list(), dots = list()) {
  
  comp_like <- match.arg(comp_like)
  
  if (is.null(dots$rescale)) {
    if (length(RCMdata@Chist) && any(RCMdata@Chist > 0, na.rm = TRUE)) {
      dots$rescale <- 1/mean(RCMdata@Chist, na.rm = TRUE)
    } else {
      dots$rescale <- 1
    }
  }
  
  TMB_params <- RCM_est_params(x, RCMdata, selectivity, s_selectivity, prior, LWT, comp_like, StockPars, FleetPars, 
                               start, map, dots)
  TMB_data <- RCM_est_data(x, RCMdata, selectivity, s_selectivity, LWT, comp_like, prior, max_F, 
                           StockPars, FleetPars, mean_fit, TMB_params$map, dots)
  
  if (integrate) {
    random <- c("log_early_rec_dev", "log_rec_dev") 
  } else {
    random <- NULL
  }
  
  obj <- MakeADFun(data = TMB_data, parameters = TMB_params$params, map = TMB_params$map, random = random,
                   inner.control = inner.control, DLL = "SAMtool", silent = TRUE)
  
  if (any(RCMdata@Misc$condition == "catch2")) {
    if (any(is.na(obj$report(obj$par)$F)) || any(is.infinite(obj$report(obj$par)$F))) {
      for(i in 1:10) {
        obj$par["R0x"] <- 0.5 + obj$par["R0x"]
        if (all(!is.na(obj$report(obj$par)$F)) && all(!is.infinite(obj$report(obj$par)$F))) break
      }
    }
  }
  
  mod <- optimize_TMB_model(obj, control)
  opt <- mod[[1]]
  SD <- mod[[2]]
  report <- obj$report(obj$env$last.par.best) %>% RCM_posthoc_adjust(obj)
  report$conv <- !is.character(opt) && SD$pdHess
  
  return(list(obj = obj, opt = opt, SD = SD, report = report))
}


RCM_est_data <- function(x, RCMdata, selectivity, s_selectivity, LWT = list(), comp_like, prior, max_F,
                         StockPars, FleetPars, mean_fit = FALSE, map, dots = list()) {
  
  SR_type <- switch(StockPars$SRrel[x],
                    "1" = "BH",
                    "2" = "Ricker",
                    "3" = "Mesnil-Rochet")
  if (is.null(SR_type)) stop("Can not identify stock-recruit function in OM@SRrel")
  
  selectivity <- as.integer(selectivity)
  s_selectivity <- as.integer(s_selectivity)
  
  nyears <- RCMdata@Misc$nyears
  nfleet <- RCMdata@Misc$nfleet
  n_age <- dim(RCMdata@CAA)[2]
  nsurvey <- ncol(RCMdata@Index)
  
  # Convert to proportions - add tiny to zeros
  RCMdata@CAA <- apply(RCMdata@CAA, c(1, 3), tiny_comp) %>% aperm(c(2, 1, 3))
  RCMdata@CAL <- apply(RCMdata@CAL, c(1, 3), tiny_comp) %>% aperm(c(2, 1, 3))
  
  RCMdata@IAA <- apply(RCMdata@IAA, c(1, 3), tiny_comp) %>% aperm(c(2, 1, 3))
  RCMdata@IAL <- apply(RCMdata@IAL, c(1, 3), tiny_comp) %>% aperm(c(2, 1, 3))
  
  if (mean_fit) {
    # Average across simulations for arrays: M_ageArray, Len_age, Mat_age, Fec_Age (mean across index 1)
    # Matrices: L5, LFS, Vmaxlen (mean across rows)
    # Vectors: Isd, Len_CV, hs, R0 (effort only)
    mean_vector <- function(x) rep(mean(x), length(x))
    mean_matrix <- function(x) matrix(apply(x, 1, mean), nrow(x), ncol(x))
    mean_array <- function(x) {
      means <- apply(x, c(2, 3), mean)
      new_output <- array(means, dim(x)[c(2,3,1)])
      return(aperm(new_output, c(3,1,2)))
    }
    
    StockPars_ind <- match(c("M_ageArray", "Len_age", "Mat_age", "Wt_age", "Fec_Age"), names(StockPars))
    StockPars[StockPars_ind] <- lapply(StockPars[StockPars_ind], mean_array)
    
    StockPars_ind <- match(c("hs", "LenCV", "procsd"), names(StockPars))
    StockPars[StockPars_ind] <- lapply(StockPars[StockPars_ind], mean_vector)
    
    StockPars_ind <- match("ageMarray", names(StockPars))
    StockPars[StockPars_ind] <- lapply(StockPars[StockPars_ind], mean_matrix)
    
    if (all(RCMdata@Misc$condition == "effort")) {
      StockPars$R0 <- mean_vector(StockPars$R0)
      if (!is.null(dots$OMeff) && dots$OMeff) {
        FleetPars$Find <- apply(FleetPars$Find, 2, mean) %>% matrix(length(StockPars$R0), nyears, byrow = TRUE)
      }
    }
    
    FleetPars_ind <- match(c("L5_y", "LFS_y", "Vmaxlen_y"), names(FleetPars))
    FleetPars[FleetPars_ind] <- lapply(FleetPars[FleetPars_ind], mean_matrix)
    
  }
  
  if (length(RCMdata@Chist) && any(RCMdata@Chist > 0, na.rm = TRUE)) {
    C_hist <- RCMdata@Chist
  } else {
    C_hist <- matrix(0, nyears, nfleet)
  }
  
  if (all(RCMdata@Misc$condition == "effort") && !is.null(dots$OMeff) && dots$OMeff) {
    RCMdata@Ehist <- matrix(FleetPars$Find[x, ], nyears, nfleet)
  }
  
  if (length(RCMdata@Ehist) && any(RCMdata@Ehist > 0, na.rm = TRUE)) {
    rescale_effort <- 1/mean(RCMdata@Ehist, na.rm = TRUE)
    E_hist <- RCMdata@Ehist * rescale_effort
    E_eq <- RCMdata@E_eq * rescale_effort
  } else {
    rescale_effort <- 1
    E_hist <- matrix(0, nyears, nfleet)
    E_eq <- rep(0, nfleet)
  }
  
  # Likelihood weights
  LWT_fleet <- cbind(LWT$Chist, LWT$C_eq, LWT$CAA, LWT$CAL, LWT$MS)
  LWT_index <- cbind(LWT$Index, LWT$IAA, LWT$IAL)
  
  if (is.null(dots$n_itF)) n_itF <- 3L else n_itF <- dots$n_itF
  if (is.null(dots$plusgroup)) plusgroup <- 1L else plusgroup <- as.integer(dots$plusgroup)
  
  TMB_data <- list(model = "RCM", C_hist = C_hist, C_eq = RCMdata@C_eq, 
                   sigma_C = RCMdata@C_sd, sigma_Ceq = RCMdata@C_eq_sd, 
                   E_hist = E_hist, E_eq = E_eq,
                   condition = sapply(RCMdata@Misc$condition, switch, "catch" = 0, "catch2" = 1, "effort" = 2) %>% as.integer(), 
                   nll_C = as.integer(any(RCMdata@Misc$condition != "catch2")),
                   I_hist = RCMdata@Index, sigma_I = RCMdata@I_sd, 
                   CAA_hist = RCMdata@CAA, CAA_n = RCMdata@CAA_ESS,
                   CAL_hist = RCMdata@CAL, CAL_n = RCMdata@CAL_ESS,
                   IAA_hist = RCMdata@IAA, IAA_n = RCMdata@IAA_ESS,
                   IAL_hist = RCMdata@IAL, IAL_n = RCMdata@IAL_ESS,
                   lbin = RCMdata@Misc$lbin, lbinmid = RCMdata@Misc$lbinmid, 
                   msize = RCMdata@MS, msize_type = RCMdata@MS_type,
                   sel_block = rbind(RCMdata@sel_block, RCMdata@sel_block[nyears, ]), 
                   nsel_block = RCMdata@Misc$nsel_block,
                   n_y = nyears, n_age = n_age, nfleet = nfleet, nsurvey = nsurvey,
                   M_data = if (prior$use_prior[3]) matrix(1, 1, 1) else t(StockPars$M_ageArray[x, , 1:nyears]), 
                   len_age = t(StockPars$Len_age[x, , 1:(nyears+1)]),
                   Linf = StockPars$Linf[x],
                   SD_LAA = t(StockPars$LatASD[x, , 1:(nyears+1)]), 
                   wt = t(StockPars$Wt_age[x, , 1:(nyears+1)]),
                   mat = t(StockPars$Mat_age[x, , 1:(nyears+1)]),
                   #mat = if (any(s_selectivity == -3L)) t(StockPars$Mat_age[x, , 1:(nyears+1)]) else matrix(1, 1, 1),
                   fec = if (is.null(StockPars$Fec_Age)) {
                     t(StockPars$Wt_age[x, , 1:(nyears+1)] * StockPars$Mat_age[x, , 1:(nyears+1)])
                   } else {
                     t(StockPars$Fec_Age[x, , 1:(nyears+1)])
                   },
                   vul_type = selectivity,
                   ivul_type = s_selectivity, 
                   abs_I = RCMdata@abs_I, 
                   I_units = as.integer(RCMdata@I_units), 
                   I_delta = if (.hasSlot(RCMdata, "I_delta")) RCMdata@I_delta else rep(0, nsurvey), 
                   age_error = RCMdata@age_error,
                   SR_type = SR_type, 
                   LWT_fleet = LWT_fleet, LWT_index = LWT_index, 
                   comp_like = comp_like,
                   max_F = max_F, 
                   rescale = dots$rescale, 
                   ageM = min(nyears, ceiling(StockPars$ageMarray[x, 1])),
                   yind_F = as.integer(rep(0.5 * nyears, nfleet)), 
                   n_itF = n_itF, plusgroup = plusgroup,
                   use_prior = prior$use_prior, 
                   prior_dist = prior$pr_matrix, 
                   nll_gr = 0L,
                   sim_process_error = 0L,
                   spawn_time_frac = ifelse(is.null(StockPars$spawn_time_frac), 0, StockPars$spawn_time_frac[x]),
                   est_q = ifelse(is.na(map$log_q), 0L, 1L),
                   pbc_recdev = if (is.null(dots$pbc_recdev)) rep(1, nyears) else dots$pbc_recdev,
                   pbc_early_recdev = if (is.null(dots$pbc_early_recdev)) rep(1, nyears) else dots$pbc_early_recdev)
  
  TMB_data$est_vul <- ifelse(is.na(map$vul_par) | duplicated(map$vul_par), 0, 1) %>%
    matrix(length(map$vul_par)/RCMdata@Misc$nsel_block, RCMdata@Misc$nsel_block)
  TMB_data$est_ivul <- ifelse(is.na(map$ivul_par) | duplicated(map$ivul_par), 0, 1) %>%
    matrix(length(map$ivul_par)/nsurvey, nsurvey)
  
  # 0 = fixed parameter: no bias correction, no contribution to likelihood
  # 1 = estimated parameter, not previous duplicated: bias correction, contribute to likelihood
  # 2 = estimated parameter, but duplicated: bias correction, no contribution to likelihood
  TMB_data$est_early_rec_dev <- ifelse(is.na(map$log_early_rec_dev), 0, 
                                       ifelse(duplicated(map$log_early_rec_dev), 2, 1))
  TMB_data$est_rec_dev <- ifelse(is.na(map$log_rec_dev), 0, 
                                 ifelse(duplicated(map$log_rec_dev), 2, 1))
  
  return(TMB_data)
}

RCM_est_params <- function(x, RCMdata, selectivity, s_selectivity, prior = list(), LWT = list(), comp_like,
                           StockPars, FleetPars, start = list(), map = list(), dots = list()) {
  
  ## Reference for parameters that can passed through start and map ----
  start_check <- c("vul_par", "ivul_par", "log_early_rec_dev", "log_rec_dev", "q", "MR_SRR")
  map_check <- start_check
  
  ## Logical to indicate whether to transfrom selectivity parameters to estimation space
  vul_transform <- TRUE
  if (!is.null(dots$vul_transform)) vul_transform <- dots$vul_transform
  
  ## Check for obsolete arguments ----
  if (!is.null(dots$vul_par)) stop("Specify vul_par via start$vul_par")
  if (!is.null(dots$ivul_par)) stop("Specify ivul_par via start$vul_par")
  if (!is.null(dots$log_rec_dev)) stop("Specify log_rec_dev via start$log_rec_dev")
  if (!is.null(dots$log_early_rec_dev)) stop("Specify log_early_rec_dev via start$log_early_rec_dev")
  
  if (!is.null(dots$map_vul_par)) stop("Specify map_vul_par via map$vul_par")
  if (!is.null(dots$map_ivul_par)) stop("Specify map_ivul_par via map$ivul_par")
  if (!is.null(dots$map_log_rec_dev)) stop("Specify map_log_rec_dev via map$log_rec_dev")
  if (!is.null(dots$map_log_early_rec_dev)) stop("Specify map_log_early_rec_dev via map$log_early_rec_dev")
  
  if (!is.null(dots$map_s_vul_par)) stop("Pass map argument of index selectivity parameters via map$ivul_par")
  if (!is.null(dots$s_vul_par)) stop("Pass index selectivity parameters via start$ivul_par")
  
  nyears <- RCMdata@Misc$nyears
  nfleet <- RCMdata@Misc$nfleet
  n_age <- dim(RCMdata@CAA)[2]
  maxage <- n_age - 1
  nsurvey <- ncol(RCMdata@Index)
  
  ## Fishery selectivity ----
  LFS <- FleetPars$LFS_y[x, nyears]
  L5 <- FleetPars$L5_y[x, nyears]
  Vmaxlen <- FleetPars$Vmaxlen_y[x, nyears]
  Linf <- StockPars$Linf[x]
  
  # Check for functional selectivity functions (dome or logistic)
  sel_check <- selectivity %in% c(0, -1, -5, -6)
  sel_check_age <- selectivity %in% c(-5, -6)
  sel_check_len <- sel_check & !sel_check_age
  
  if (is.null(start$vul_par)) {
    if (any(selectivity == -2)) {
      stop("Some fleet selectivity specified to be free parameters. Provide start$vul_par matrix to RCM.")
    }
    start$vul_par <- matrix(c(LFS, L5, Vmaxlen), 3, RCMdata@Misc$nsel_block)
    
    start$vul_par[1, sel_check_age] <- 0.5 * maxage
    start$vul_par[2, sel_check_age] <- 0.25 * maxage
  }
  if (ncol(start$vul_par) != RCMdata@Misc$nsel_block) {
    stop("start$vul_par needs to be a matrix with ", RCMdata@Misc$nsel_block, " columns")
  }
  if (any(selectivity == -2) && nrow(start$vul_par) != n_age) {
    stop("start$vul_par needs to be a matrix with ", n_age, " (maxage + 1) rows")
  }
  
  if (vul_transform && any(sel_check)) {  # Parametric sel
    start$vul_par[2, sel_check] <- log(start$vul_par[1, sel_check] - start$vul_par[2, sel_check])
    
    start$vul_par[1, sel_check_len] <- logit(pmin(start$vul_par[1, sel_check_len]/Linf, 0.99))
    start$vul_par[1, sel_check_age] <- logit(pmin(start$vul_par[1, sel_check_age]/maxage, 0.99))
    
    start$vul_par[3, sel_check] <- logit(pmin(start$vul_par[3, sel_check], 0.99))
  }
  
  if (is.null(map$vul_par)) {
    if (any(selectivity == -2)) {
      stop("Some (fleet) selectivity specified to be free parameters. Provide map$vul_par matrix to RCM.")
    }
    map$vul_par <- matrix(0, 3, RCMdata@Misc$nsel_block)
    map$vul_par[3, selectivity %in% c(-1, -6)] <- NA # Fix third parameter for logistic sel
    if (!is.null(dots$fix_dome) && dots$fix_dome) { # Obsolete
      warning("fix_dome is obsolete. Recommend using the map argument in RCM")
      map$vul_par[3, selectivity %in% c(0, -5)] <- NA # Fix dome
    }
    
    for(ff in 1:nfleet) {
      if (all(RCMdata@CAA[,,ff] <= 0, na.rm = TRUE) && all(RCMdata@CAL[,,ff] <= 0, na.rm = TRUE)) {
        map$vul_par[, unique(RCMdata@sel_block[, ff])] <- NA # Fix sel if no comp data
      }
    }
    if (any(!is.na(map$vul_par))) {
      map$vul_par[!is.na(map$vul_par)] <- 1:sum(!is.na(map$vul_par))
    }
  }
  if (ncol(map$vul_par) != RCMdata@Misc$nsel_block) {
    stop("map$vul_par needs to be a matrix with ", RCMdata@Misc$nsel_block, " columns")
  }
  
  if (any(selectivity == -2) ) { # Free parameters, convert start matrix to logit space
    if (nrow(map$vul_par) != n_age) stop("map$vul_par needs to be a matrix with ", n_age, " (maxage + 1) rows")
    
    if (vul_transform) {
      test <- start$vul_par[, selectivity == -2] >= 0 & start$vul_par[, selectivity == -2] <= 1
      if (any(!test)) stop("There are free selectivity parameters in start$vul_par that are less than zero or greater than one.")
      test_est <- !is.na(map$vul_par[, selectivity == -2])
      test_fix <- test & !test_est
      start$vul_par[, selectivity == -2][test_fix] <- logit(start$vul_par[, selectivity == -2][test_fix], soft_bounds = FALSE)
      start$vul_par[, selectivity == -2][test_est] <- logit(start$vul_par[, selectivity == -2][test_est], soft_bounds = TRUE)
    }
  }
  
  ## Index selectivity ----
  # Check for function selectivity (dome/logistic)
  sel_check <- s_selectivity %in% c(0, -1, -5, -6)
  sel_check_age <- s_selectivity %in% c(-5, -6)
  sel_check_len <- sel_check & !sel_check_age
  
  if (is.null(start$ivul_par)) {
    if (any(s_selectivity == -2)) {
      stop("Some s_selectivity specified to be free parameters. Provide start$ivul_par matrix to RCM.")
    }
    start$ivul_par <- matrix(c(LFS, L5, Vmaxlen), 3, nsurvey)
    
    start$ivul_par[1, sel_check_age] <- 0.5 * maxage
    start$ivul_par[2, sel_check_age] <- 0.25 * maxage
  }
  if (ncol(start$ivul_par) != nsurvey) {
    stop("start$ivul_par needs to be a matrix with ", nsurvey, " columns.")
  }
  if (any(s_selectivity == -2) && nrow(start$ivul_par) != n_age) {
    stop("start$ivul_par needs to be a matrix with ", n_age, " (maxage + 1) rows.")
  }
  
  if (vul_transform && any(sel_check)) { # Parametric sel
    start$ivul_par[2, sel_check] <- log(start$ivul_par[1, sel_check] - start$ivul_par[2, sel_check])
    
    start$ivul_par[1, sel_check_len] <- logit(pmin(start$ivul_par[1, sel_check_len]/Linf, 0.99))
    start$ivul_par[1, sel_check_age] <- logit(pmin(start$ivul_par[1, sel_check_age]/maxage, 0.99))
    
    start$ivul_par[3, sel_check] <- logit(pmin(start$ivul_par[3, sel_check], 0.99))
  }
  
  if (is.null(map$ivul_par)) {
    if (any(s_selectivity == -2)) {
      stop("Some s_selectivity specified to be free parameters. Provide map$ivul_par matrix to RCM.")
    }
    map$ivul_par <- matrix(0, 3, nsurvey)
    map$ivul_par[3, s_selectivity < 0] <- NA # if logistic
    for(sur in 1:nsurvey) {
      if (s_selectivity[sur] %in% c(-3, 4) || s_selectivity[sur] > 0 ||
         (all(RCMdata@IAA[,,sur] <= 0, na.rm = TRUE) & all(RCMdata@IAL[,,sur] <= 0, na.rm = TRUE))) {
        map$ivul_par[, sur] <- NA
      }
    }
    if (any(!is.na(map$ivul_par))) {
      map$ivul_par[!is.na(map$ivul_par)] <- 1:sum(!is.na(map$ivul_par))
    }
  }
  if (ncol(map$ivul_par) != nsurvey) {
    stop("map$ivul_par needs to be a matrix with ", nsurvey, " columns.")
  }
  
  if (any(s_selectivity == -2)) {  # Free parameters, convert start matrix to logit space
    if (nrow(map$ivul_par) != n_age) stop("map$ivul_par needs to be a matrix with ", n_age, " (maxage + 1) rows.")
    
    if (vul_transform) {
      test <- start$ivul_par[, s_selectivity == -2] >= 0 & start$ivul_par[, s_selectivity == -2] <= 1
      
      if (any(!test)) stop("There are free selectivity parameters in start$ivul_par that are less than zero or greater than one.")
      
      test_est <- !is.na(map$ivul_par[, s_selectivity == -2])
      test_fix <- test & !test_est
      start$ivul_par[, s_selectivity == -2][test_fix] <- logit(start$ivul_par[, s_selectivity == -2][test_fix], soft_bounds = FALSE)
      start$ivul_par[, s_selectivity == -2][test_est] <- logit(start$ivul_par[, s_selectivity == -2][test_est], soft_bounds = TRUE)
    }
  }
  
  if (is.null(start$log_early_rec_dev)) start$log_early_rec_dev <- rep(0, n_age - 1)
  if (length(start$log_early_rec_dev) != n_age - 1) stop("start$log_early_rec_dev is not a vector of length n_age - 1")
  
  if (is.null(start$log_rec_dev)) start$log_rec_dev <- rep(0, nyears)
  if (length(start$log_rec_dev) != nyears) stop("start$log_rec_dev is not a vector of length nyears")
  
  ## Mesnil Rochet parameters ----
  if (is.null(start$MR_SRR)) start$MR_SRR <- c(0.8, 0.001)
  if (length(start$MR_SRR) != 2) stop("start$MR_SRR needs to be a length 2 vector.")
  
  ## Index q ----
  # can be either estimated by TMB or analytically calculated - see map argument later
  if (is.null(start$q)) start$q <- rep(1, nsurvey)
  if (length(start$q) != nsurvey) stop("start$q needs to be a length ", nsurvey, " vector.")
  
  ## SRR and steepness ----
  SR_type <- switch(StockPars$SRrel[x],
                    "1" = "BH",
                    "2" = "Ricker",
                    "3" = "Mesnil-Rochet")
  if (is.null(SR_type)) stop("Can not identify stock-recruit function in OM@SRrel")
  
  transformed_h <- switch(SR_type,
                          "BH" = logit((StockPars$hs[x] - 0.2)/0.8),
                          "Ricker" = log(StockPars$hs[x] - 0.2),
                          "Mesnil-Rochet" = 0)
  
  ## Full parameter list ----
  TMB_params <- list(R0x = ifelse(!is.na(StockPars$R0[x]), log(StockPars$R0[x] * dots$rescale), 0),
                     transformed_h = transformed_h, 
                     MR_SRR = c(logit(start$MR_SRR[1]), log(start$MR_SRR[2])), # Shinge/S0, gamma 
                     log_M = log(mean(StockPars$M_ageArray[x, , nyears])),
                     vul_par = start$vul_par, 
                     ivul_par = start$ivul_par,
                     log_q_effort = rep(log(0.1), nfleet), 
                     log_F_dev = matrix(0, nyears, nfleet),
                     log_F_equilibrium = rep(log(0.05), nfleet),
                     log_CV_msize = log(RCMdata@MS_cv), 
                     log_tau = log(StockPars$procsd[x]),
                     log_early_rec_dev = start$log_early_rec_dev, 
                     log_rec_dev = start$log_rec_dev,
                     log_compf = matrix(0, nfleet, 2),
                     log_compi = matrix(0, nsurvey, 2),
                     log_q = log(start$q))
  
  if (any(RCMdata@Misc$condition == "catch")) {
    TMB_params$log_F_dev[as.integer(0.5 * nyears) + 1, 
                         RCMdata@Misc$condition == "catch"] <- log(0.5 * mean(StockPars$M_ageArray[x, , nyears]))
  }
  
  ## Map list (to fix parameters) ----
  map_out <- list()
  
  if (all(RCMdata@Misc$condition == "effort") && !sum(RCMdata@Chist, na.rm = TRUE) && !prior$use_prior[1]) {
    map_out$R0x <- factor(NA) # Fix if condition on effort, no catches, and no prior on R0
  }
  if (SR_type == "Mesnil-Rochet" || !prior$use_prior[2]) map_out$transformed_h <- factor(NA)
  
  if (SR_type == "Mesnil-Rochet") {
    if (is.null(map$MR_SRR)) map$MR_SRR <- c(1, NA)  # Estimates Ehinge/E0, fixes gamma = 0.001 --> hockey-stick
    map_out$MR_SRR <- factor(map$MR_SRR)
  } else {
    map_out$MR_SRR <- factor(c(NA, NA))
  }
  
  if (!prior$use_prior[3]) map_out$log_M <- factor(NA)
  
  map_out$log_tau <- factor(NA)
  map_out$vul_par <- factor(map$vul_par)
  map_out$ivul_par <- factor(map$ivul_par)
  
  map_out$log_F_equilibrium <- local({
    m <- rep(FALSE, nfleet)
    m[grepl("catch", RCMdata@Misc$condition) & RCMdata@C_eq > 0] <- TRUE
    
    m[m] <- 1:sum(m, na.rm = TRUE)
    m[!m] <- NA
    factor(m)
  })
  
  map_out$log_q_effort <- local({
    m <- grepl("effort", RCMdata@Misc$condition)
    m[m] <- 1:sum(m)
    m[!m] <- NA
    factor(m)
  })
  
  map_out$log_F_dev <- local({
    m <- matrix(FALSE, nyears, nfleet)
    for(ff in 1:nfleet) {
      if (RCMdata@Misc$condition[ff] == "catch") m[, ff] <- TRUE
    }
    m[m] <- 1:sum(m)
    m[!m] <- NA
    factor(m)
  })
  
  map_out$log_CV_msize <- factor(rep(NA, nfleet))
  
  if (is.null(map$log_early_rec_dev)) {
    map_out$log_early_rec_dev <- factor(rep(NA, n_age - 1))
  } else {
    map_out$log_early_rec_dev <- factor(map$log_early_rec_dev)
    if (length(map_out$log_early_rec_dev) != n_age - 1) {
      stop("map$log_early_rec_dev needs to be a vector of length ", n_age - 1)
    }
  }
  
  if (is.null(map$log_rec_dev)) {
    map_out$log_rec_dev <- factor(1:nyears)
  } else {
    map_out$log_rec_dev <- factor(map$log_rec_dev)
    if (length(map_out$log_rec_dev) != nyears) {
      stop("map$log_rec_dev needs to be a vector of length ", nyears)
    }
  }
  
  map_out$log_compf <- local({
    mapf <- matrix(NA, nfleet, 2)
    
    if (comp_like %in% c("dirmult1", "dirmult2")) {
      CAAv <- sapply(1:nfleet, function(ff) LWT$CAA[ff] > 0 && sum(RCMdata@CAA_ESS[, ff]) > 0)
      CALv <- sapply(1:nfleet, function(ff) LWT$CAL[ff] > 0 && sum(RCMdata@CAL_ESS[, ff]) > 0)
      
      est <- cbind(CAAv, CALv)
      mapf[est] <- 1:sum(est)
    }
    
    factor(mapf)
  })
  
  map_out$log_compi <- local({
    mapi <- matrix(NA, nsurvey, 2)
    
    if (comp_like %in% c("dirmult1", "dirmult2")) {
      IAAv <- sapply(1:nsurvey, function(sur) LWT$IAA[sur] > 0 && sum(RCMdata@IAA_ESS[, sur]) > 0)
      IALv <- sapply(1:nsurvey, function(sur) LWT$IAL[sur] > 0 && sum(RCMdata@IAL_ESS[, sur]) > 0)
      
      est <- cbind(IAAv, IALv)
      mapi[est] <- 1:sum(est)
    }
    factor(mapi)
  })
  
  # If q is not estimated, it is solved analytically - can still fix it to 1 in RCMdata@abs_I
  if (is.null(map$q)) map$q <- rep(NA, nsurvey)
  if (length(map$q) != nsurvey) stop("map$q must be length ", nsurvey)
  map_out$log_q <- factor(map$q)
  
  list(params = TMB_params, map = map_out)
}

par_identical_sims_fn <- function(StockPars, FleetPars, RCMdata, dots) {
  vector_fn <- function(x) sum(mean(x) - x) == 0
  array_fn <- function(x) {
    x_mean <- apply(x, 2:length(dim(x)), mean)
    all(apply(x, 1, identical, x_mean))
  }
  run_test <- function(x) if (is.null(dim(x))) vector_fn(x) else array_fn(x)
  
  StockPars_subset <- StockPars[c("hs", "procsd", "ageMarray", "M_ageArray", "Linf", "Len_age", "Wt_age", "Mat_age", "Fec_Age")]
  if (any(RCMdata@CAL > 0, na.rm = TRUE) || any(RCMdata@IAL > 0, na.rm = TRUE) || 
     (RCMdata@MS_type == "length" & any(RCMdata@MS > 0, na.rm = TRUE))) {
    StockPars_subset <- c(StockPars_subset, StockPars["LenCV"])
  }
  S_test <- vapply(StockPars_subset, run_test, logical(1))
  
  if (RCMdata@Misc$nfleet == 1 && !any(RCMdata@CAL > 0, na.rm = TRUE) && !any(RCMdata@CAA > 0, na.rm = TRUE)) {
    FleetPars_subset <- FleetPars[c("L5_y", "LFS_y", "Vmaxlen_y")]
    FleetPars_subset <- lapply(FleetPars_subset, function(x) x[, RCMdata@Misc$nyears])
    F_test_sel <- vapply(FleetPars_subset, run_test, logical(1))
  } else F_test_sel <- NULL
  
  if (all(RCMdata@Misc$condition == "effort") && !is.null(dots$OMeff) && dots$OMeff) {
    F_test_Find <- run_test(FleetPars$Find) %>% structure(names = "Find")
  } else F_test_Find <- NULL
  
  F_test <- c(F_test_sel, F_test_Find)
  
  test_all <- c(S_test, F_test)
  return(test_all)
}


RCM_dynamic_SSB0 <- function(obj, par = obj$env$last.par.best) {
  
  if (any(obj$env$data$condition == 1L)) { # all will be catch2
    
    new_data <- obj$env$data
    new_data$C_hist[] <- 1e-8
    
    obj2 <- MakeADFun(
      data = new_data, parameters = clean_tmb_parameters(obj), 
      map = obj$env$map, random = obj$env$random,
      DLL = "SAMtool", silent = TRUE
    )
    
  } else {
    obj2 <- obj
  }
  
  par[names(par) == "log_F_dev" | names(par) == "log_F_equilibrium"] <- log(1e-8)
  par[names(par) == "log_q_effort"] <- log(1e-8)
  
  out <- obj2$report(par)$E
  return(out)
}

# Calculate annual values of:
# F-at-age, dynamic and static SPR, dynamic SSB0, compensation ratio,
# steepness, unfished reference points (remove from output unfished numbers per recruit),
# selectivity at length
RCM_posthoc_adjust <- function(report, obj, par = obj$env$last.par.best, dynamic_SSB0 = TRUE) {
  data <- obj$env$data
  if (data$use_prior[3]) {
    M <- matrix(report$Mest, data$n_y, data$n_age)
  } else {
    M <- data$M_data
  }
  report$F_at_age <- report$Z - M
  report$NPR_unfished <- do.call(rbind, report$NPR_unfished)
  report$SPR_eq <- RCM_SPR(F_at_age = report$F_at_age, M = M, fec = data$fec, 
                           plusgroup = data$plusgroup, spawn_time_frac = data$spawn_time_frac)
  report$SPR_dyn <- RCM_SPR(F_at_age = report$F_at_age, M = M, fec = data$fec,
                            N_at_age = report$N, R = report$R, R_early = report$R_early, equilibrium = FALSE, 
                            plusgroup = data$plusgroup, spawn_time_frac = data$spawn_time_frac)
  
  report$ageM <- data$ageM
  if (data$SR_type == "BH") {
    report$E0 <- pmax((report$Arec * report$EPR0 - 1)/report$Brec, 0)
    report$h_annual <- local({
      h <- report$Arec * report$EPR0/(4 + report$Arec * report$EPR0)
      ifelse(h < 0.2, NA_real_, h)
    })
    
    report$CR <- report$Arec * report$EPR0 # Annual compensation ratio recalculated from annual EPR0
  } else if (data$SR_type == "Ricker") {
    report$E0 <- pmax(log(report$Arec * report$EPR0)/report$Brec, 0)
    report$h_annual <- local({
      h <- 0.2 * (report$Arec * report$EPR0)^0.8
      ifelse(h < 0.2, NA_real_, h)
    })
    
    report$CR <- report$Arec * report$EPR0 # Annual compensation ratio recalculated from annual EPR0
  } else if (data$SR_type == "Mesnil-Rochet") {
    
    MR_K <- sqrt(report$MRhinge^2 + 0.25 * report$MRgamma^2)
    MR_beta <- report$MRRmax/(report$MRhinge + MR_K)
    
    report$E0 <- local({
      num <- 2 * MR_K / report$EPR0 / MR_beta - 2 * (report$MRhinge + MR_K)
      den <- 1/report$EPR0/report$EPR0/MR_beta/MR_beta - 2/report$EPR0/MR_beta
      
      ifelse(1/report$EPR0 > 2 * MR_beta, 0, num/den)
    })
    
    report$CR <- 2 * MR_beta * report$EPR0 # Annual compensation ratio recalculated from annual EPR0
  }
  report$R0_annual <- report$E0/report$EPR0
  
  report$N0 <- apply(report$NPR_unfished * report$R0_annual, 1, sum)
  report$B0 <- apply(report$NPR_unfished * report$R0_annual * data$wt[1:data$n_y, ], 1, sum)

  lmid <- obj$env$data$lbinmid
  nlbin <- length(lmid)
  
  spawn_time_frac <- data$spawn_time_frac
  if (spawn_time_frac > 0) report$E[length(report$E)] <- NA
  if (dynamic_SSB0) {
    report$dynamic_SSB0 <- RCM_dynamic_SSB0(obj, par)
    if (spawn_time_frac > 0) report$dynamic_SSB0[length(report$dynamic_SSB0)] <- NA
  }
  
  if (data$comp_like == "mvlogistic") {
    
    CAAmv <- calc_mvlogistic_loglike(obs = data$CAA_hist, pred = report$CAApred, 
                                     LWT = data$LWT_fleet[, 3], nllv = report$nll_fleet[1, , 3])
    CALmv <- calc_mvlogistic_loglike(obs = data$CAL_hist, pred = report$CALpred, 
                                     LWT = data$LWT_fleet[, 4], nllv = report$nll_fleet[1, , 4])
    
    IAAmv <- calc_mvlogistic_loglike(obs = data$IAA_hist, pred = report$IAApred, 
                                     LWT = data$LWT_index[, 2], nllv = report$nll_index[1, , 2])
    
    if (is.null(report$IALpred)) {
      IALmv <- list(tau = rep(NaN, data$nsurvey), nll = matrix(0, data$n_y, data$nsurvey))
    } else {
      IALmv <- calc_mvlogistic_loglike(obs = data$IAL_hist, pred = report$IALpred, 
                                       LWT = data$LWT_index[, 3], nllv = report$nll_index[1, , 3])
    }
    
    report$nll_fleet[, 1:data$nfleet, 3] <- CAAmv$nll
    report$nll_fleet[, 1:data$nfleet, 4] <- CALmv$nll
    
    report$nll_index[, 1:data$nsurvey, 2] <- IAAmv$nll
    report$nll_index[, 1:data$nsurvey, 3] <- IALmv$nll
    
    report$log_compf <- log(cbind(CAAmv$tau, CALmv$tau))
    report$log_compi <- log(cbind(IAAmv$tau, IALmv$tau))
  }
  return(report)
}

get_vul_len <- function(report, selectivity, lmid, Linf) {
  vul <- matrix(NA_real_, length(lmid), length(selectivity))
  sel_ind  <- selectivity == 0 | selectivity == -1
  LFS <- report$LFS[sel_ind]
  L5 <- report$L5[sel_ind]
  Vmaxlen <- report$Vmaxlen[sel_ind]
  
  sls <- (LFS - L5)/sqrt(-log(0.05, 2))
  srs <- (Linf - LFS)/sqrt(-log(Vmaxlen, 2))
  asc <- Map(function(x, y) 2^-((lmid - y)/x * (lmid - y)/x), x = sls, y = LFS)
  dsc <- Map(function(x, y, z) ifelse(z > rep(0.99, length(lmid)), 1, 2^-((lmid - y)/x * (lmid - y)/x)),
             x = srs, y = LFS, z = Vmaxlen)
  vul_out <- Map(function(x, y, z) ifelse(lmid > x, y, z), x = LFS, y = dsc, z = asc)
  vul[, sel_ind] <- do.call(cbind, vul_out)
  return(vul)
}

get_ivul_len <- function(report, s_selectivity, lmid, Linf) {
  ivul_len <- matrix(NA_real_, length(lmid), length(s_selectivity)) # length-based: matrix of dimension nlbin, nsurvey
  for(i in 1:ncol(ivul_len)) {
    if (s_selectivity[i] == -1 || s_selectivity[i] == 0) {
      sls <- (report$iLFS[i] - report$iL5[i])/sqrt(-log(0.05, 2))
      srs <- (Linf - report$iLFS[i])/sqrt(-log(report$iVmaxlen[i], 2))
      
      asc <- 2^-((lmid - report$iLFS[i])/sls * (lmid - report$iLFS[i])/sls)
      dsc <- ifelse(report$iVmaxlen[i] > rep(0.99, length(lmid)), 1,
                    2^-((lmid - report$iLFS[i])/srs * (lmid - report$iLFS[i])/srs))
      ivul_len[, i] <- ifelse(lmid > report$iLFS[i], dsc, asc)
    } else if (s_selectivity[i] > 0) {
      ivul_len[, i] <- report$vul_len[, s_selectivity[i]]
    }
  }
  return(ivul_len)
}

RCM_SPR <- function(F_at_age, M, fec, N_at_age, R, R_early, equilibrium = TRUE, 
                    plusgroup = TRUE, spawn_time_frac = 0) {
  n_y <- nrow(F_at_age)
  n_age <- ncol(F_at_age)
  Z <- F_at_age + M
  
  if (equilibrium) {
    SSPR_F <- vapply(1:n_y, function(y) {
      NPR <- calc_NPR(exp(-Z[y, ]), n_age, plusgroup)
      sum(NPR * exp(-spawn_time_frac * Z[y, ]) * fec[y, ])
    }, numeric(1))
    
    SSPR_0 <- vapply(1:n_y, function(y) {
      NPR <- calc_NPR(exp(-M[y, ]), n_age, plusgroup)
      sum(NPR * exp(-spawn_time_frac * M[y, ]) * fec[y, ])
    }, numeric(1))
    
  } else {
    NPR_M <- NPR_F <- matrix(1, n_y, n_age)
    
    RR <- R_early %>% rev() %>% c(R[1])
    NPR_F[1, ] <- N_at_age[1, ]/rev(RR)
    
    NPR_M[1, -1] <- exp(-cumsum(M[1, -n_age]))
    if (plusgroup) NPR_M[1, n_age] <- NPR_M[1, n_age]/(1 - exp(-M[1, n_age]))
    for(y in 2:n_y) {
      for(a in 2:n_age) {
        NPR_M[y, a] <- NPR_M[y-1, a-1] * exp(-M[y-1, a-1])
        NPR_F[y, a] <- NPR_F[y-1, a-1] * exp(-Z[y-1, a-1])
      }
      if (plusgroup) {
        NPR_M[y, n_age] <- NPR_M[y, n_age] + NPR_M[y-1, n_age] * exp(-M[y-1, n_age])
        NPR_F[y, n_age] <- NPR_F[y, n_age] + NPR_F[y-1, n_age] * exp(-Z[y-1, n_age])
      }
    }
    SSPR_F <- rowSums(NPR_F * exp(-spawn_time_frac * Z) * fec[1:n_y, ])
    SSPR_0 <- rowSums(NPR_M * exp(-spawn_time_frac * M) * fec[1:n_y, ])
  }
  SPR <- SSPR_F/SSPR_0
  return(SPR)
}


dmvlogistic <- function(x, p, sd, xmin = 1e-8, log = FALSE) {
  resid <- log(x[x > xmin]) - log(p[x > xmin])
  if (!length(resid)) stop("Density function can not be calculated from x")
  accum <- log(sum(x[x <= xmin])) - log(sum(p[x <= xmin]))
  
  eta <- c(resid, accum) - mean(c(resid, accum))
  
  A <- length(eta)
  
  # No normalizing constants!
  log_like <- -(A-1) * log(sd) - 0.5 * sum(eta^2)/sd/sd
  if (log) {
    log_like
  } else {
    exp(log_like)
  }
}


calc_mvlogistic_loglike <- function(obs, pred, nllv, LWT, obsmin = 1e-8) {
  nf <- dim(obs)[3]
  n_y <- dim(obs)[1]
  
  vars <- lapply(1:nf, function(ff) {
    sum_count <- sapply(1:n_y, function(y) {
      if (sum(obs[y, , ff]) > 0 && LWT[ff] > 0) {
        A <- sum(obs[y, , ff] > obsmin)
        if (any(obs[y, , ff] <= obsmin)) A <- A + 1
        return(A-1)
      } else {
        return(0)
      }
    }) %>% sum()
    
    #tau2 <- exp((-nllv[ff] + 0.5 * sum_count)/(-0.5 * sum_count))
    tau2 <- exp(nllv[ff]/(0.5 * sum_count) - 1)
    tau <- sqrt(tau2)
    
    log_like <- sapply(1:n_y, function(y) {
      if (sum(obs[y, , ff]) > 0 && LWT[ff] > 0) {
        log_like <- dmvlogistic(x = obs[y, , ff], p = pred[y, , ff]/sum(pred[y, , ff]), 
                                sd = tau, xmin = obsmin, log = TRUE)
      } else {
        log_like <- 0
      }
      return(log_like)
    })
    
    list(tau = tau, log_like = log_like)
  })
  
  list(tau = sapply(vars, getElement, "tau"),
       nll = -1 * sapply(vars, getElement, "log_like"))
}

