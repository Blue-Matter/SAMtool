
RCM_est <- function(x = 1, RCMdata, selectivity, s_selectivity, LWT = list(),
                    comp_like = c("multinomial", "lognormal", "mvlogistic", "dirmult1", "dirmult2"), prior = list(),
                    max_F = 3, integrate = FALSE, StockPars, ObsPars, FleetPars, mean_fit = FALSE,
                    control = list(iter.max = 2e+05, eval.max = 4e+05), inner.control = list(maxit = 1e3), dots = list()) {
  
  comp_like <- match.arg(comp_like)
  
  if(is.null(dots$rescale)) {
    if(length(RCMdata@Chist) && any(RCMdata@Chist > 0, na.rm = TRUE)) {
      dots$rescale <- 1/mean(RCMdata@Chist, na.rm = TRUE)
    } else {
      dots$rescale <- 1
    }
  }
  
  TMB_params <- RCM_est_params(x, RCMdata, selectivity, s_selectivity, prior, LWT, comp_like, StockPars, FleetPars, dots)
  TMB_data <- RCM_est_data(x, RCMdata, selectivity, s_selectivity, LWT, comp_like, prior, max_F, 
                           StockPars, ObsPars, FleetPars, mean_fit, TMB_params$map, dots)
  
  if(integrate) {
    random <- c("log_early_rec_dev", "log_rec_dev") 
  } else {
    random <- NULL
  }
  
  obj <- MakeADFun(data = TMB_data, parameters = TMB_params$params, map = TMB_params$map, random = random,
                   inner.control = inner.control, DLL = "SAMtool", silent = TRUE)
  
  if(RCMdata@Misc$condition == "catch2") {
    if(any(is.na(obj$report(obj$par)$F)) || any(is.infinite(obj$report(obj$par)$F))) {
      for(i in 1:10) {
        obj$par["R0x"] <- 0.5 + obj$par["R0x"]
        if(all(!is.na(obj$report(obj$par)$F)) && all(!is.infinite(obj$report(obj$par)$F))) break
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
                         StockPars, ObsPars, FleetPars, mean_fit = FALSE, map, dots = list()) {
  
  SR_type <- ifelse(StockPars$SRrel[x] == 1, "BH", "Ricker")
  
  nyears <- RCMdata@Misc$nyears
  nfleet <- RCMdata@Misc$nfleet
  n_age <- dim(RCMdata@CAA)[2]
  nsurvey <- ncol(RCMdata@Index)
  
  # Convert to proportions - add tiny to zeros
  RCMdata@CAA <- apply(RCMdata@CAA, c(1, 3), tiny_comp) %>% aperm(c(2, 1, 3))
  RCMdata@CAL <- apply(RCMdata@CAL, c(1, 3), tiny_comp) %>% aperm(c(2, 1, 3))
  
  RCMdata@IAA <- apply(RCMdata@IAA, c(1, 3), tiny_comp) %>% aperm(c(2, 1, 3))
  RCMdata@IAL <- apply(RCMdata@IAL, c(1, 3), tiny_comp) %>% aperm(c(2, 1, 3))
  
  if(mean_fit) {
    # Average across simulations for arrays: M_ageArray, Len_age, Mat_age (mean across index 1)
    # Matrices: L5, LFS, Vmaxlen (mean across rows)
    # Vectors: Isd, Len_CV, hs, R0 (effort only)
    mean_vector <- function(x) rep(mean(x), length(x))
    mean_matrix <- function(x) matrix(apply(x, 1, mean), nrow(x), ncol(x))
    mean_array <- function(x) {
      means <- apply(x, c(2, 3), mean)
      new_output <- array(means, dim(x)[c(2,3,1)])
      return(aperm(new_output, c(3,1,2)))
    }
    
    StockPars_ind <- match(c("M_ageArray", "Len_age", "Mat_age"), names(StockPars))
    StockPars[StockPars_ind] <- lapply(StockPars[StockPars_ind], mean_array)
    
    StockPars_ind <- match(c("hs", "LenCV", "procsd"), names(StockPars))
    StockPars[StockPars_ind] <- lapply(StockPars[StockPars_ind], mean_vector)
    
    StockPars_ind <- match("ageMarray", names(StockPars))
    StockPars[StockPars_ind] <- lapply(StockPars[StockPars_ind], mean_matrix)
    
    if(RCMdata@Misc$condition == "effort") {
      StockPars$R0 <- mean_vector(StockPars$R0)
      if(!is.null(dots$OMeff) && dots$OMeff) {
        FleetPars$Find <- apply(FleetPars$Find, 2, mean) %>% matrix(length(StockPars$R0), nyears, byrow = TRUE)
      }
    }
    
    FleetPars_ind <- match(c("L5_y", "LFS_y", "Vmaxlen_y"), names(FleetPars))
    FleetPars[FleetPars_ind] <- lapply(FleetPars[FleetPars_ind], mean_matrix)
    
    ObsPars$Isd <- mean_vector(ObsPars$Isd)
  }
  
  if(!length(RCMdata@I_sd) || !any(RCMdata@I_sd > 0, na.rm = TRUE)) {
    RCMdata@I_sd <- matrix(sdconv(1, ObsPars$Isd[x]), nyears, nsurvey)
  }
  
  if(length(RCMdata@Chist) && any(RCMdata@Chist > 0, na.rm = TRUE)) {
    C_hist <- RCMdata@Chist
  } else {
    C_hist <- matrix(0, nyears, nfleet)
  }
  
  if(RCMdata@Misc$condition == "effort" && !is.null(dots$OMeff) && dots$OMeff) {
    RCMdata@Ehist <- matrix(FleetPars$Find[x, ], nyears, nfleet)
  }
  
  if(length(RCMdata@Ehist) && any(RCMdata@Ehist > 0, na.rm = TRUE)) {
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
  
  if(is.null(dots$n_itF)) n_itF <- 3L else n_itF <- dots$n_itF
  if(is.null(dots$plusgroup)) plusgroup <- 1L else plusgroup <- as.integer(dots$plusgroup)
  age_only_model <- StockPars$Len_age[x, , 1:(nyears+1)] %>%
    apply(2, function(xx) all(xx == 1:n_age)) %>% 
    all()
  
  TMB_data <- list(model = "RCM", C_hist = C_hist, C_eq = RCMdata@C_eq, 
                   sigma_C = RCMdata@C_sd, sigma_Ceq = RCMdata@C_eq_sd, 
                   E_hist = E_hist, E_eq = E_eq,
                   condition = RCMdata@Misc$condition, nll_C = as.integer(RCMdata@Misc$condition != "catch2"),
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
                   M_data = if(prior$use_prior[3]) matrix(1, 1, 1) else t(StockPars$M_ageArray[x, , 1:nyears]), 
                   len_age = t(StockPars$Len_age[x, , 1:(nyears+1)]),
                   Linf = ifelse(age_only_model, n_age, StockPars$Linf[x]),
                   SD_LAA = t(StockPars$LatASD[x, , 1:nyears]), 
                   wt = t(StockPars$Wt_age[x, , 1:(nyears+1)]),
                   mat = t(StockPars$Mat_age[x, , 1:(nyears+1)]),
                   vul_type = as.integer(selectivity),
                   ivul_type = as.integer(s_selectivity), 
                   abs_I = RCMdata@abs_I, 
                   I_units = as.integer(RCMdata@I_units), 
                   age_error = RCMdata@age_error,
                   SR_type = SR_type, 
                   LWT_fleet = LWT_fleet, LWT_index = LWT_index, 
                   comp_like = comp_like,
                   max_F = max_F, 
                   rescale = dots$rescale, 
                   ageM = min(nyears, ceiling(StockPars$ageM[x, 1])),
                   yind_F = as.integer(rep(0.5 * nyears, nfleet)), 
                   n_itF = n_itF, plusgroup = plusgroup,
                   use_prior = prior$use_prior, 
                   prior_dist = prior$pr_matrix, 
                   nll_gr = 0L)
                   
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
                           StockPars, FleetPars, dots = list()) {
  
  SR_type <- ifelse(StockPars$SRrel[x] == 1, "BH", "Ricker")
  
  nyears <- RCMdata@Misc$nyears
  nfleet <- RCMdata@Misc$nfleet
  n_age <- dim(RCMdata@CAA)[2]
  nsurvey <- ncol(RCMdata@Index)
  
  if(SR_type == "BH") {
    transformed_h <- logit((StockPars$hs[x] - 0.2)/0.8)
  } else transformed_h <- log(StockPars$hs[x] - 0.2)
  
  LFS <- FleetPars$LFS_y[x, nyears]
  L5 <- FleetPars$L5_y[x, nyears]
  Vmaxlen <- FleetPars$Vmaxlen_y[x, nyears]
  
  if(is.null(dots$vul_par)) {
    if(any(selectivity == -2)) {
      stop("Some fleet selectivity specified to be free parameters. Provide vul_par matrix to RCM.")
    }
    vul_par <- matrix(c(LFS, L5, Vmaxlen), 3, RCMdata@Misc$nsel_block)
  } else {
    vul_par <- dots$vul_par
    
    if(ncol(vul_par) != RCMdata@Misc$nsel_block) {
      stop("vul_par needs to be a matrix with ", RCMdata@Misc$nsel_block, " columns")
    }
    if(any(selectivity == -2) && nrow(vul_par) != n_age) {
      stop("vul_par needs to be a matrix with ", n_age, " (maxage + 1) rows")
    }
    
  }
  
  # Check for functional selectivity functions (dome or logistic)
  age_only_model <- StockPars$Len_age[x, , 1:(nyears+1)] %>%
    apply(2, function(xx) all(xx == 1:n_age)) %>% 
    all()
  Linf <- ifelse(age_only_model, n_age, StockPars$Linf[x])
  
  sel_check <- selectivity == -1 | selectivity == 0
  vul_par[2, sel_check] <- log(vul_par[1, sel_check] - vul_par[2, sel_check])
  vul_par[1, sel_check] <- logit(pmin(vul_par[1, sel_check]/Linf/0.99, 0.99))
  vul_par[3, sel_check] <- logit(pmin(vul_par[3, sel_check], 0.99))
  
  if(any(selectivity == -2)) {
    vul_par[, selectivity == -2] <- logit(vul_par[, selectivity == -2], soft_bounds = TRUE)
  }
  
  if(is.null(dots$map_vul_par)) {
    if(any(selectivity == -2)) {
      stop("Some (fleet) selectivity specified to be free parameters. Provide map_vul_par matrix to RCM.")
    }
    map_vul_par <- matrix(0, 3, RCMdata@Misc$nsel_block)
    map_vul_par[3, selectivity == -1] <- NA # Fix third parameter for logistic sel
    if(!is.null(dots$fix_dome) && dots$fix_dome) { # Obsolete
      map_vul_par[3, selectivity == 0] <- NA # Fix dome
    }
    
    for(ff in 1:nfleet) {
      if(all(RCMdata@CAA[,,ff] <= 0, na.rm = TRUE) && all(RCMdata@CAL[,,ff] <= 0, na.rm = TRUE)) {
        map_vul_par[, unique(RCMdata@sel_block[, ff])] <- NA # Fix sel if no comp data
      }
    }
    if(any(!is.na(map_vul_par))) {
      map_vul_par[!is.na(map_vul_par)] <- 1:sum(!is.na(map_vul_par))
    }
  } else {
    map_vul_par <- dots$map_vul_par
    
    if(ncol(map_vul_par) != RCMdata@Misc$nsel_block) {
      stop("map_vul_par needs to be a matrix with ", RCMdata@Misc$nsel_block, " columns")
    }
    if(any(selectivity == -2) && nrow(map_vul_par) != n_age) {
      stop("map_vul_par needs to be a matrix with ", n_age, " (maxage + 1) rows")
    }
  }
  
  if(any(selectivity == -2)) {
    test <- dots$vul_par[, selectivity == -2] %in% c(0, 1) & is.na(map_vul_par[, selectivity == -2])
    vul_par[, selectivity == -2][test] <- logit(dots$vul_par[, selectivity == -2][test], soft_bounds = FALSE)
  }
  
  # Index selectivity (ivul_par)
  if(!is.null(dots$s_vul_par)) dots$ivul_par <- dots$s_vul_par # Backwards compatibility
  if(is.null(dots$ivul_par)) {
    if(any(s_selectivity == -2)) {
      stop("Some s_selectivity specified to be free parameters. Provide ivul_par matrix to RCM.")
    }
    ivul_par <- matrix(c(LFS, L5, Vmaxlen), 3, nsurvey)
  } else {
    ivul_par <- dots$ivul_par
    
    if(ncol(ivul_par) != nsurvey) {
      stop("ivul_par needs to be a matrix with ", nsurvey, " columns.")
    }
    if(any(s_selectivity == -2) && nrow(ivul_par) != n_age) {
      stop("ivul_par needs to be a matrix with ", n_age, " (maxage + 1) rows.")
    }
  }
  
  # Check for function selectivity (dome/logistic)
  parametric_sel <- s_selectivity == -1 | s_selectivity == 0
  ivul_par[2, parametric_sel] <- log(ivul_par[1, parametric_sel] - ivul_par[2, parametric_sel])
  ivul_par[1, parametric_sel] <- logit(ivul_par[1, parametric_sel]/Linf/0.99)
  ivul_par[3, parametric_sel] <- logit(ivul_par[3, parametric_sel])
  
  if(any(s_selectivity == -2)) {
    ivul_par[, s_selectivity == -2] <- logit(ivul_par[, s_selectivity == -2], soft_bounds = TRUE)
  }
  if(!is.null(dots$map_s_vul_par)) dots$map_ivul_par <- dots$map_s_vul_par # Backwards compatibility
  
  if(is.null(dots$map_ivul_par)) {
    if(any(s_selectivity == -2)) {
      stop("Some s_selectivity specified to be free parameters. Provide map_ivul_par matrix to RCM.")
    }
    map_ivul_par <- matrix(0, 3, nsurvey)
    map_ivul_par[3, s_selectivity < 0] <- NA # if logistic
    for(sur in 1:nsurvey) {
      if(s_selectivity[sur] < -2 || s_selectivity[sur] > 0 ||
         (all(RCMdata@IAA[,,sur] <= 0, na.rm = TRUE) & all(RCMdata@IAL[,,sur] <= 0, na.rm = TRUE))) {
        map_ivul_par[, sur] <- NA
      }
    }
    if(any(!is.na(map_ivul_par))) {
      map_ivul_par[!is.na(map_ivul_par)] <- 1:sum(!is.na(map_ivul_par))
    }
  } else {
    map_ivul_par <- dots$map_ivul_par
    
    if(ncol(map_ivul_par) != nsurvey) {
      stop("ivul_par needs to be a matrix with ", nsurvey, " columns.")
    }
    if(any(s_selectivity == -2) && nrow(map_ivul_par) != n_age) {
      stop("ivul_par needs to be a matrix with ", n_age, " (maxage + 1) rows.")
    }
  }
  
  if(any(s_selectivity == -2)) {
    test <- dots$ivul_par[, s_selectivity == -2] %in% c(0, 1) & is.na(map_ivul_par[, s_selectivity == -2])
    ivul_par[, s_selectivity == -2][test] <- logit(dots$ivul_par[, s_selectivity == -2][test], soft_bounds = FALSE)
  }
  
  if(!is.null(dots$log_early_rec_dev)) {
    if(length(dots$log_early_rec_dev) != n_age - 1) stop("early_rec_dev is not a vector of length n_age - 1")
    log_early_rec_dev <- dots$log_early_rec_dev
  } else {
    log_early_rec_dev <- rep(0, n_age - 1)
  }
  
  if(!is.null(dots$log_rec_dev)) {
    if(length(dots$log_rec_dev) != nyears) stop("log_rec_dev is not a vector of length nyears")
    log_rec_dev <- dots$log_rec_dev
  } else {
    log_rec_dev <- rep(0, nyears)
  }
  
  TMB_params <- list(R0x = ifelse(!is.na(StockPars$R0[x]), log(StockPars$R0[x] * dots$rescale), 0),
                     transformed_h = transformed_h, log_M = log(mean(StockPars$M_ageArray[x, , nyears])),
                     vul_par = vul_par, ivul_par = ivul_par,
                     log_q_effort = rep(log(0.1), nfleet), log_F_dev = matrix(0, nyears, nfleet),
                     log_F_equilibrium = rep(log(0.05), nfleet),
                     log_CV_msize = log(RCMdata@MS_cv), log_tau = log(StockPars$procsd[x]),
                     log_early_rec_dev = log_early_rec_dev, log_rec_dev = log_rec_dev,
                     log_compf = matrix(0, nfleet, 2), log_compi = matrix(0, nsurvey, 2))
  
  if(RCMdata@Misc$condition == "catch") {
    TMB_params$log_F_dev[as.integer(0.5 * nyears) + 1, ] <- log(0.5 * mean(StockPars$M_ageArray[x, , nyears]))
  }
  
  # Map list (to fix parameters)
  map <- list()
  
  if(RCMdata@Misc$condition == "effort" && !sum(RCMdata@Chist, na.rm = TRUE) && !prior$use_prior[1]) {
    map$R0x <- factor(NA) # Fix if condition on effort, no catches, and no prior on R0
  }
  if(!prior$use_prior[2]) map$transformed_h <- factor(NA)
  if(!prior$use_prior[3]) map$log_M <- factor(NA)
  
  map$log_tau <- factor(NA)
  map$vul_par <- factor(map_vul_par)
  map$ivul_par <- factor(map_ivul_par)
  
  if(RCMdata@Misc$condition == "effort") {
    map$log_F_equilibrium <- factor(rep(NA, nfleet))
  } else {
    map$log_q_effort <- factor(rep(NA, nfleet))
    
    if(any(RCMdata@C_eq == 0)) {
      map$log_F_equilibrium <- local({
        m <- rep(NA, nfleet)
        m[RCMdata@C_eq > 0] <- 1:sum(RCMdata@C_eq > 0)
        factor(m)
      })
    }
  }
  
  if(RCMdata@Misc$condition != "catch") {
    map$log_F_dev <- factor(matrix(NA, nyears, nfleet))
  }
  map$log_CV_msize <- factor(rep(NA, nfleet))
  
  if(is.null(dots$map_log_early_rec_dev)) {
    map$log_early_rec_dev <- factor(rep(NA, n_age - 1))
  } else {
    map$log_early_rec_dev <- factor(dots$map_log_early_rec_dev)
    if(length(map$log_early_rec_dev) != n_age - 1) {
      stop("map_log_early_rec_dev needs to be a vector of length ", n_age - 1)
    }
  }
  
  if(is.null(dots$map_log_rec_dev)) {
    map$log_rec_dev <- factor(1:nyears)
  } else {
    map$log_rec_dev <- factor(dots$map_log_rec_dev)
    if(length(map$log_rec_dev) != nyears) {
      stop("map_log_rec_dev needs to be a vector of length ", nyears)
    }
  }
  
  map$log_compf <- local({
    mapf <- matrix(NA, nfleet, 2)
    
    if(comp_like %in% c("dirmult1", "dirmult2")) {
      CAAv <- sapply(1:nfleet, function(ff) LWT$CAA[ff] > 0 && sum(RCMdata@CAA_ESS[, ff]) > 0)
      CALv <- sapply(1:nfleet, function(ff) LWT$CAL[ff] > 0 && sum(RCMdata@CAL_ESS[, ff]) > 0)
      
      est <- cbind(CAAv, CALv)
      mapf[est] <- 1:sum(est)
    }
    
    factor(mapf)
  })
  
  map$log_compi <- local({
    mapi <- matrix(NA, nsurvey, 2)
    
    if(comp_like %in% c("dirmult1", "dirmult2")) {
      IAAv <- sapply(1:nsurvey, function(sur) LWT$IAA[sur] > 0 && sum(RCMdata@IAA_ESS[, sur]) > 0)
      IALv <- sapply(1:nsurvey, function(sur) LWT$IAL[sur] > 0 && sum(RCMdata@IAL_ESS[, sur]) > 0)
      
      est <- cbind(IAAv, IALv)
      mapi[est] <- 1:sum(est)
    }
    factor(mapi)
  })
  
  list(params = TMB_params, map = map)
}

par_identical_sims_fn <- function(StockPars, FleetPars, ObsPars, RCMdata, dots) {
  vector_fn <- function(x) sum(mean(x) - x) == 0
  array_fn <- function(x) {
    x_mean <- apply(x, 2:length(dim(x)), mean)
    all(apply(x, 1, identical, x_mean))
  }
  run_test <- function(x) if(is.null(dim(x))) vector_fn(x) else array_fn(x)
  
  StockPars_subset <- StockPars[c("hs", "procsd", "ageMarray", "M_ageArray", "Linf", "Len_age", "Wt_age", "Mat_age")]
  if(any(RCMdata@CAL > 0, na.rm = TRUE) || any(RCMdata@IAL > 0, na.rm = TRUE) || 
     (RCMdata@MS_type == "length" & any(RCMdata@MS > 0, na.rm = TRUE))) {
    StockPars_subset <- c(StockPars_subset, StockPars["LenCV"])
  }
  S_test <- vapply(StockPars_subset, run_test, logical(1))
  
  if(RCMdata@Misc$nfleet == 1 && !any(RCMdata@CAL > 0, na.rm = TRUE) && !any(RCMdata@CAA > 0, na.rm = TRUE)) {
    FleetPars_subset <- FleetPars[c("L5_y", "LFS_y", "Vmaxlen_y")]
    FleetPars_subset <- lapply(FleetPars_subset, function(x) x[, RCMdata@Misc$nyears])
    F_test_sel <- vapply(FleetPars_subset, run_test, logical(1))
  } else F_test_sel <- NULL
  
  if(RCMdata@Misc$condition == "effort" && !is.null(dots$OMeff) && dots$OMeff) {
    F_test_Find <- run_test(FleetPars$Find) %>% structure(names = "Find")
  } else F_test_Find <- NULL
  
  F_test <- c(F_test_sel, F_test_Find)
  
  if(RCMdata@Misc$nsurvey > 0 && !any(RCMdata@I_sd > 0, na.rm = TRUE)) {
    O_test <- run_test(ObsPars$Isd) %>% structure(names = "Isd")
  } else {
    O_test <- NULL
  }
  
  test_all <- c(S_test, F_test, O_test)
  return(test_all)
}


RCM_dynamic_SSB0 <- function(obj, par = obj$env$last.par.best) {
  
  if(obj$env$data$condition == "catch") {
    
    par[names(par) == "log_F_dev" | names(par) == "log_F_equilibrium"] <- log(1e-8)
    out <- obj$report(par)$E
    
  } else if(obj$env$data$condition == "catch2") {
    
    new_data <- obj$env$data
    new_data$C_hist[] <- 1e-8
    
    new_params <- lapply(obj$env$parameters, function(x) if(!is.null(attr(x, "map"))) attr(x, "shape") else x)
    
    obj2 <- MakeADFun(data = new_data, parameters = new_params, 
                      map = obj$env$map, random = obj$env$random,
                      DLL = "SAMtool", silent = TRUE)
    out <- obj2$report(par)$E
    
  } else {
    
    par[names(par) == "log_q_effort"] <- log(1e-8)
    out <- obj$report(par)$E
    
  }
  
  return(out)
}

# Calculate annual values of:
# F-at-age, dynamic and static SPR, dynamic SSB0, compensation ratio,
# steepness, unfished reference points (remove from output unfished numbers per recruit),
# selectivity at length
RCM_posthoc_adjust <- function(report, obj, par = obj$env$last.par.best, dynamic_SSB0 = TRUE) {
  data <- obj$env$data
  if(data$use_prior[3]) {
    M <- matrix(report$Mest, data$n_y, data$n_age)
  } else {
    M <- data$M_data
  }
  report$F_at_age <- report$Z - M
  report$NPR_unfished <- do.call(rbind, report$NPR_unfished)
  report$SPR_eq <- RCM_SPR(F_at_age = report$F_at_age, M = M, mat = data$mat, wt = data$wt)
  report$SPR_dyn <- RCM_SPR(F_at_age = report$F_at_age, M = M, mat = data$mat, wt = data$wt, 
                            N_at_age = report$N, R = report$R, R_early = report$R_early,
                            equilibrium = FALSE)
  
  report$ageM <- data$ageM
  if(data$SR_type == "BH") {
    report$E0 <- pmax((report$Arec * report$EPR0 - 1)/report$Brec, 0)
    report$h_annual <- local({
      h <- report$Arec * report$EPR0/(4 + report$Arec * report$EPR0)
      ifelse(h < 0.2, NA_real_, h)
    })
  } else {
    report$E0 <- pmax(log(report$Arec * report$EPR0)/report$Brec, 0)
    report$h_annual <- local({
      h <- 0.2 * (report$Arec * report$EPR0)^0.8
      ifelse(h < 0.2, NA_real_, h)
    })
  }
  report$R0_annual <- report$E0/report$EPR0
  report$N0 <- apply(report$NPR_unfished * report$R0_annual, 1, sum)
  report$B0 <- apply(report$NPR_unfished * report$R0_annual * data$wt[1:data$n_y, ], 1, sum)

  lmid <- obj$env$data$lbinmid
  nlbin <- length(lmid)
  
  report$CR <- report$Arec * report$EPR0 # Annual compensation ratio recalculated from annual EPR0
  
  age_only_model <- data$len_age %>%
    apply(1, function(x) all(x == 1:data$n_age)) %>% 
    all()
  if(age_only_model) {
    report$vul_len <- matrix(NA_real_, nlbin, data$nsel_block)
    report$ivul_len <- matrix(NA_real_, nlbin, dim(report$ivul)[3])
    
    report$MLpred <- array(NA_real_, dim(report$F))
    report$CALpred <- array(NA_real_, dim(report$CALpred))
    report$IALpred <- array(NA_real_, dim(report$IALpred))
  } else {
    report$vul_len <- get_vul_len(report, data$vul_type, lmid, data$Linf)
    report$ivul_len <- get_ivul_len(report, data$ivul_type, lmid, data$Linf)
  }
  if(dynamic_SSB0) report$dynamic_SSB0 <- RCM_dynamic_SSB0(obj, par)
  
  if(data$comp_like == "mvlogistic") {
    
    CAAmv <- calc_mvlogistic_loglike(obs = data$CAA_hist, pred = report$CAApred, 
                                     LWT = data$LWT_fleet[, 3], nllv = report$nll_fleet[1, , 3])
    CALmv <- calc_mvlogistic_loglike(obs = data$CAL_hist, pred = report$CALpred, 
                                     LWT = data$LWT_fleet[, 4], nllv = report$nll_fleet[1, , 4])
    
    IAAmv <- calc_mvlogistic_loglike(obs = data$IAA_hist, pred = report$IAApred, 
                                     LWT = data$LWT_index[, 2], nllv = report$nll_index[1, , 2])
    
    if(is.null(report$IALpred)) {
      IALmv <- list(tau = rep(NaN, data$nsurvey), nll = matrix(0, data$n_y, data$nsurvey))
    } else {
      IALmv <- calc_mvlogistic_loglike(obs = data$IAL_hist, pred = report$IALpred, 
                                       LWT = data$LWT_index[, 3], nllv = report$nll_index[1, , 3])
    }
    
    report$nll_fleet[, 1:data$nfleet, 3] <- CAAmv$nll
    report$nll_fleet[, 1:data$nfleet, 4] <- CALmv$nll
    
    report$nll_index[, 1:data$nsurvey, 2] <- IAAmv$nll
    report$nll_index[, 1:data$nsurvey, 3] <- IALmv$nll
    
    report$compf <- cbind(CAAmv$tau, CALmv$tau)
    report$compi <- cbind(IAAmv$tau, IALmv$tau)
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
    if(s_selectivity[i] == -1 || s_selectivity[i] == 0) {
      sls <- (report$iLFS[i] - report$iL5[i])/sqrt(-log(0.05, 2))
      srs <- (Linf - report$iLFS[i])/sqrt(-log(report$iVmaxlen[i], 2))
      
      asc <- 2^-((lmid - report$iLFS[i])/sls * (lmid - report$iLFS[i])/sls)
      dsc <- ifelse(report$iVmaxlen[i] > rep(0.99, length(lmid)), 1,
                    2^-((lmid - report$iLFS[i])/srs * (lmid - report$iLFS[i])/srs))
      ivul_len[, i] <- ifelse(lmid > report$iLFS[i], dsc, asc)
    } else if(s_selectivity[i] > 0) {
      ivul_len[, i] <- report$vul_len[, s_selectivity[i]]
    }
  }
  return(ivul_len)
}

RCM_SPR <- function(F_at_age, M, mat, wt, N_at_age, R, R_early, equilibrium = TRUE) {
  n_y <- nrow(F_at_age)
  n_age <- ncol(F_at_age)
  
  if(equilibrium) { # Plusgroup always on
    SSPR_F <- vapply(1:n_y, function(y) {
      yield_fn_SCA_int(max(F_at_age[y, ]), M = M[y, ], mat = mat[y, ], weight = wt[y, ], 
                       vul = F_at_age[y, ]/max(F_at_age[y, ]), Arec = 1, Brec = 1, opt = FALSE)["EPR"]
    }, numeric(1))
    
    SSPR_0 <- vapply(1:n_y, function(y) {
      yield_fn_SCA_int(0, M = M[y, ], mat = mat[y, ], weight = wt[y, ], 
                       vul = rep(1, n_age), Arec = 1, Brec = 1, opt = FALSE)["EPR"]
    }, numeric(1))
  } else {
    NPR_M <- NPR_F <- matrix(1, n_y, n_age)
    
    RR <- R_early %>% rev() %>% c(R[1])
    NPR_F[1, ] <- N_at_age[1, ]/rev(RR)
    
    NPR_M[1, -1] <- exp(-cumsum(M[1, -n_age]))
    NPR_M[1, n_age] <- NPR_M[1, n_age]/(1 - exp(-M[1, n_age]))
    for(y in 2:n_y) {
      for(a in 2:n_age) {
        NPR_M[y, a] <- NPR_M[y-1, a-1] * exp(-M[y-1, a-1])
        NPR_F[y, a] <- NPR_F[y-1, a-1] * exp(-F_at_age[y-1, a-1] - M[y-1, a-1])
      }
      NPR_M[y, n_age] <- NPR_M[y, n_age] + NPR_M[y-1, n_age] * exp(-M[y-1, n_age])
      NPR_F[y, n_age] <- NPR_F[y, n_age] + NPR_F[y-1, n_age] * exp(-F_at_age[y-1, n_age] - M[y-1, n_age])
    }
    SSPR_F <- rowSums(NPR_F * wt[1:n_y, ] * mat[1:n_y, ])
    SSPR_0 <- rowSums(NPR_M * wt[1:n_y, ] * mat[1:n_y, ])
  }
  SPR <- SSPR_F/SSPR_0
  return(SPR)
}


dmvlogistic <- function(x, p, sd, xmin = 1e-8, log = FALSE) {
  resid <- log(x[x > xmin]) - log(p[x > xmin])
  if(!length(resid)) stop("Density function can not be calculated from x")
  accum <- log(sum(x[x <= xmin])) - log(sum(p[x <= xmin]))
  
  eta <- c(resid, accum) - mean(c(resid, accum))
  
  A <- length(eta)
  
  # No normalizing constants!
  log_like <- -(A-1) * log(sd) - 0.5 * sum(eta^2)/sd/sd
  if(log) {
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
      if(sum(obs[y, , ff]) > 0 && LWT[ff] > 0) {
        A <- sum(obs[y, , ff] > obsmin)
        if(any(obs[y, , ff] <= obsmin)) A <- A + 1
        return(A-1)
      } else {
        return(0)
      }
    }) %>% sum()
    
    #tau2 <- exp((-nllv[ff] + 0.5 * sum_count)/(-0.5 * sum_count))
    tau2 <- exp(nllv[ff]/(0.5 * sum_count) - 1)
    tau <- sqrt(tau2)
    
    log_like <- sapply(1:n_y, function(y) {
      if(sum(obs[y, , ff]) > 0 && LWT[ff] > 0) {
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

