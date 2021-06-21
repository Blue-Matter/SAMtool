
RCM_est <- function(x = 1, RCMdata, selectivity, s_selectivity, LWT = list(),
                    comp_like = c("multinomial", "lognormal"), prior = list(),
                    max_F = 3, integrate = FALSE, StockPars, ObsPars, FleetPars, mean_fit = FALSE,
                    control = list(iter.max = 2e+05, eval.max = 4e+05), inner.control = list(maxit = 1e3), dots = list()) {
  
  comp_like <- match.arg(comp_like)
  
  SR_type <- ifelse(StockPars$SRrel[x] == 1, "BH", "Ricker")
  
  nyears <- RCMdata@Misc$nyears
  nfleet <- RCMdata@Misc$nfleet
  n_age <- dim(RCMdata@CAA)[2]
  nsurvey <- ncol(RCMdata@Index)
  
  # Convert to proportions?
  RCMdata@CAA <- apply(RCMdata@CAA, c(1, 3), tiny_comp) %>% aperm(c(2, 1, 3))
  RCMdata@CAL <- apply(RCMdata@CAL, c(1, 3), tiny_comp) %>% aperm(c(2, 1, 3))
  
  RCMdata@IAA <- apply(RCMdata@IAA, c(1, 3), tiny_comp) %>% aperm(c(2, 1, 3))
  RCMdata@IAL <- apply(RCMdata@IAL, c(1, 3), tiny_comp) %>% aperm(c(2, 1, 3))
  
  #Backwards compatibility
  if(!is.null(dots$ESS) && length(dots$ESS) == 2) {
    RCMdata@CAA_ESS <- pmin(RCMdata@CAA_ESS, dots$ESS[1])
    RCMdata@CAL_ESS <- pmin(RCMdata@CAL_ESS, dots$ESS[2])
    
    RCMdata@IAA_ESS <- pmin(RCMdata@IAA_ESS, dots$ESS[1])
    RCMdata@IAL_ESS <- pmin(RCMdata@IAL_ESS, dots$ESS[1])
  }
  
  LWT_fleet <- cbind(LWT$Chist, LWT$C_eq, LWT$CAA, LWT$CAL, LWT$MS)
  LWT_index <- cbind(LWT$Index, LWT$IAA, LWT$IAL)
  
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
        FleetPars$Find <- apply(FleetPars$Find, 2, mean) %>% matrix(length(StockPars$R0), data$nyears, byrow = TRUE)
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
    if(is.null(dots$rescale)) {
      rescale <- 1/mean(RCMdata@Chist, na.rm = TRUE)
    } else {
      rescale <- dots$rescale
    }
    C_hist <- RCMdata@Chist
  } else {
    rescale <- 1
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
  
  if(is.null(dots$n_itF)) n_itF <- 3L else n_itF <- dots$n_itF
  if(is.null(dots$plusgroup)) plusgroup <- 1L else plusgroup <- as.integer(dots$plusgroup)
  age_only_model <- StockPars$Len_age[x, , 1:(nyears+1)] %>% apply(2, function(xx) length(xx) == n_age & max(xx) == n_age) %>% all()
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
                   sel_block = rbind(RCMdata@sel_block, RCMdata@sel_block[nyears, ]), nsel_block = RCMdata@Misc$nsel_block,
                   n_y = nyears, n_age = n_age, nfleet = nfleet, nsurvey = nsurvey,
                   M_data = if(prior$use_prior[3]) matrix(1, 1, 1) else t(StockPars$M_ageArray[x, , 1:nyears]), 
                   len_age = t(StockPars$Len_age[x, , 1:(nyears+1)]),
                   Linf = ifelse(age_only_model, n_age, StockPars$Linf[x]),
                   SD_LAA = t(StockPars$LatASD[x, , 1:nyears]), wt = t(StockPars$Wt_age[x, , 1:(nyears+1)]),
                   mat = t(StockPars$Mat_age[x, , 1:(nyears+1)]), vul_type = as.integer(selectivity),
                   ivul_type = as.integer(s_selectivity), abs_I = RCMdata@abs_I, I_units = as.integer(RCMdata@I_units), 
                   age_error = RCMdata@age_error,
                   SR_type = SR_type, LWT_fleet = LWT_fleet, LWT_index = LWT_index, comp_like = comp_like,
                   max_F = max_F, rescale = rescale, ageM = min(nyears, ceiling(StockPars$ageM[x, 1])),
                   yind_F = as.integer(rep(0.5 * nyears, nfleet)), n_itF = n_itF, plusgroup = plusgroup,
                   use_prior = prior$use_prior, prior_dist = prior$pr_matrix, nll_gr = 0L)
  
  if(SR_type == "BH") {
    transformed_h <- logit((StockPars$hs[x] - 0.2)/0.8)
  } else transformed_h <- log(StockPars$hs[x] - 0.2)
  
  LFS <- FleetPars$LFS_y[x, nyears]
  L5 <- FleetPars$L5_y[x, nyears]
  Vmaxlen <- FleetPars$Vmaxlen_y[x, nyears]
  
  if(is.null(dots$vul_par)) {
    if(any(selectivity == -2)) stop("Some (fleet) selectivity specified to be free parameters. Provide vul_par matrix to RCM.")
    vul_par <- matrix(c(LFS, L5, Vmaxlen), 3, RCMdata@Misc$nsel_block)
  } else {
    vul_par <- dots$vul_par
  }
  
  sel_check <- selectivity == -1 | selectivity == 0
  vul_par[2, sel_check] <- log(vul_par[1, sel_check] - vul_par[2, sel_check])
  vul_par[1, sel_check] <- logit(pmin(vul_par[1, sel_check]/TMB_data$Linf/0.99, 0.99))
  vul_par[3, sel_check] <- logit(pmin(vul_par[3, sel_check], 0.99))
  if(any(selectivity == -2)) vul_par[, selectivity == -2] <- logit(vul_par[, selectivity == -2], soft_bounds = TRUE)
  
  if(is.null(dots$map_vul_par)) {
    if(any(selectivity == -2)) stop("Some (fleet) selectivity specified to be free parameters. Provide map_vul_par matrix to RCM.")
    map_vul_par <- matrix(0, 3, RCMdata@Misc$nsel_block)
    map_vul_par[3, selectivity == -1] <- NA # Fix third parameter for logistic sel
    if(!is.null(dots$fix_dome) && dots$fix_dome) map_vul_par[3, selectivity == 0] <- NA # Fix dome
    
    for(ff in 1:nfleet) {
      if(all(RCMdata@CAA[,,ff] <= 0, na.rm = TRUE) && all(RCMdata@CAL[,,ff] <= 0, na.rm = TRUE)) {
        map_vul_par[, unique(RCMdata@sel_block[, ff])] <- NA # Fix sel if no comp data
      }
    }
    if(any(!is.na(map_vul_par))) map_vul_par[!is.na(map_vul_par)] <- 1:sum(!is.na(map_vul_par))
  } else {
    map_vul_par <- dots$map_vul_par
  }
  if(any(selectivity == -2)) {
    test <- dots$vul_par[, selectivity == -2] %in% c(0, 1) & is.na(map_vul_par[, selectivity == -2])
    vul_par[, selectivity == -2][test] <- logit(dots$vul_par[, selectivity == -2][test], soft_bounds = FALSE)
  }
  
  # ivul_par, and map
  if(!is.null(dots$s_vul_par)) dots$ivul_par <- dots$s_vul_par # Backwards compatibility
  if(is.null(dots$ivul_par)) {
    if(any(s_selectivity == -2)) stop("Some s_selectivity specified to be free parameters. Provide ivul_par matrix to RCM.")
    ivul_par <- matrix(c(LFS, L5, Vmaxlen), 3, nsurvey)
  } else {
    ivul_par <- dots$ivul_par
  }
  
  parametric_sel <- s_selectivity == -1 | s_selectivity == 0
  ivul_par[2, parametric_sel] <- log(ivul_par[1, parametric_sel] - ivul_par[2, parametric_sel])
  ivul_par[1, parametric_sel] <- logit(ivul_par[1, parametric_sel]/TMB_data$Linf/0.99)
  ivul_par[3, parametric_sel] <- logit(ivul_par[3, parametric_sel])
  if(any(s_selectivity == -2)) ivul_par[, s_selectivity == -2] <- logit(ivul_par[, s_selectivity == -2], soft_bounds = TRUE)
  
  if(!is.null(dots$map_s_vul_par)) dots$map_ivul_par <- dots$map_s_vul_par # Backwards compatibility
  if(is.null(dots$map_ivul_par)) {
    if(any(s_selectivity == -2)) stop("Some s_selectivity specified to be free parameters. Provide map_ivul_par matrix to RCM.")
    map_ivul_par <- matrix(0, 3, nsurvey)
    map_ivul_par[3, s_selectivity < 0] <- NA # if logistic
    for(sur in 1:nsurvey) {
      if(s_selectivity[sur] < -2 || s_selectivity[sur] > 0 ||
         (all(RCMdata@IAA[,,sur] <= 0, na.rm = TRUE) & all(RCMdata@IAL[,,sur] <= 0, na.rm = TRUE))) {
        map_ivul_par[, sur] <- NA
      }
    }
    if(any(!is.na(map_ivul_par))) map_ivul_par[!is.na(map_ivul_par)] <- 1:sum(!is.na(map_ivul_par))
  } else {
    map_ivul_par <- dots$map_ivul_par
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
  
  TMB_params <- list(R0x = ifelse(!is.na(StockPars$R0[x]), log(StockPars$R0[x] * rescale), 0),
                     transformed_h = transformed_h, log_M = log(mean(StockPars$M_ageArray[x, , nyears])),
                     vul_par = vul_par, ivul_par = ivul_par,
                     log_q_effort = rep(log(0.1), nfleet), log_F_dev = matrix(0, nyears, nfleet),
                     log_F_equilibrium = rep(log(0.05), nfleet),
                     log_CV_msize = log(RCMdata@MS_cv), log_tau = log(StockPars$procsd[x]),
                     log_early_rec_dev = log_early_rec_dev, log_rec_dev = log_rec_dev)
  if(RCMdata@Misc$condition == "catch") {
    TMB_params$log_F_dev[TMB_data$yind_F + 1, 1:nfleet] <- log(0.5 * mean(StockPars$M_ageArray[x, , nyears]))
  }
  
  map <- list()
  if(RCMdata@Misc$condition == "effort" && !sum(TMB_data$C_hist, na.rm = TRUE) && !prior$use_prior[1]) map$R0x <- factor(NA)
  if(!prior$use_prior[2]) map$transformed_h <- factor(NA)
  if(!prior$use_prior[3]) map$log_M <- factor(NA)
  map$log_tau <- factor(NA)
  map$vul_par <- factor(map_vul_par)
  map$ivul_par <- factor(map_ivul_par)
  if(RCMdata@Misc$condition != "effort") {
    map$log_q_effort <- factor(rep(NA, nfleet))
    if(any(RCMdata@C_eq == 0)) {
      map_log_F_equilibrium <- rep(NA, nfleet)
      map_log_F_equilibrium[RCMdata@C_eq > 0] <- 1:sum(RCMdata@C_eq > 0)
      map$log_F_equilibrium <- factor(map_log_F_equilibrium)
    }
  } else {
    map$log_F_equilibrium <- factor(rep(NA, nfleet))
  }
  if(RCMdata@Misc$condition != "catch") map$log_F_dev <- factor(matrix(NA, nyears, nfleet))
  map$log_CV_msize <- factor(rep(NA, nfleet))
  
  if(is.null(dots$map_log_early_rec_dev)) {
    map$log_early_rec_dev <- factor(rep(NA, n_age - 1))
  } else {
    map$log_early_rec_dev <- factor(dots$map_log_early_rec_dev)
  }
  TMB_data$est_early_rec_dev <- ifelse(is.na(map$log_early_rec_dev), 0, 1)
  
  if(is.null(dots$map_log_rec_dev)) {
    map$log_rec_dev <- factor(1:nyears)
  } else {
    map$log_rec_dev <- factor(dots$map_log_rec_dev)
  }
  TMB_data$est_rec_dev <- ifelse(is.na(map$log_rec_dev), 0, 1)
  
  if(integrate) random <- c("log_early_rec_dev", "log_rec_dev") else random <- NULL
  
  obj <- MakeADFun(data = TMB_data, parameters = TMB_params, map = map, random = random,
                   inner.control = inner.control, DLL = "SAMtool", silent = TRUE)
  
  if(RCMdata@Misc$condition == "catch2") {
    if(any(is.na(obj$report(obj$par)$F)) || any(is.infinite(obj$report(obj$par)$F))) {
      for(i in 1:10) {
        obj$par["R0x"] <- 0.5 + obj$par["R0x"]
        if(all(!is.na(obj$report(obj$par)$F)) && all(!is.infinite(obj$report(obj$par)$F))) break
      }
    }
  }
  
  mod <- optimize_TMB_model(obj, control, restart = 0)
  opt <- mod[[1]]
  SD <- mod[[2]]
  report <- obj$report(obj$env$last.par.best) %>% RCM_posthoc_adjust(obj)
  
  return(list(obj = obj, opt = opt, SD = SD, report = c(report, list(conv = !is.character(opt) && SD$pdHess))))
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
    
    new_args <- RCM_retro_subset(obj$env$data$n_y, data = obj$env$data, params = obj$env$parameters, map = obj$env$map)
    new_args$data$C_hist <- matrix(1e-8, new_args$data$n_y, new_args$data$nfleet)
    
    obj2 <- MakeADFun(data = new_args$data, parameters = new_args$params, map = new_args$map, random = obj$env$random,
                      DLL = "SAMtool", silent = TRUE)
    out <- obj2$report(par)$E
    
  } else {
    
    par[names(par) == "log_q_effort"] <- log(1e-8)
    out <- obj$report(par)$E
    
  }
  
  return(out)
}

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

  lmid <- obj$env$data$lbinmid
  nlbin <- length(lmid)
  
  report$CR <- report$Arec * report$EPR0 # Annual compensation ratio recalculated from annual EPR0
  
  age_only_model <- data$len_age %>%
    apply(1, function(x) length(x) == data$n_age && max(x) == data$n_age) %>% all()
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
  return(report)
}

get_vul_len <- function(report, selectivity, lmid, Linf) {
  vul <- matrix(NA, length(lmid), length(selectivity))
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
  ivul_len <- matrix(NA, length(lmid), length(s_selectivity)) # length-based: matrix of dimension nlbin, nsurvey
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
  SSPR_0 <- vapply(1:n_y, function(y) {
    yield_fn_SCA_int(0, M = M[y, ], mat = mat[y, ], weight = wt[y, ], vul = F_at_age[y, ]/max(F_at_age[y, ]), 
                             Arec = 1, Brec = 1, opt = FALSE)["SPR"]
  }, numeric(1))
  
  if(equilibrium) {
    SSPR_F <- vapply(1:n_y, function(y) {
      yield_fn_SCA_int(max(F_at_age[y, ]), M = M[y, ], mat = mat[y, ], weight = wt[y, ], 
                       vul = F_at_age[y, ]/max(F_at_age[y, ]), Arec = 1, Brec = 1, opt = FALSE)["SPR"]
    }, numeric(1))
  } else {
    n_age <- ncol(F_at_age)
    SSPR_F <- vapply(1:n_y, function(y) {
      if(y < n_age) {
        RR <- R_early[1:(n_age-y)] %>% rev() %>% c(R[1:y])
      } else {
        RR <- R[(y - n_age + 1):y]
      }
      sum(N_at_age[y, ] * wt[y, ] * mat[y, ]/ rev(RR))
    }, numeric(1))
  }
  SPR <- SSPR_F/SSPR_0
  return(SPR)
}
