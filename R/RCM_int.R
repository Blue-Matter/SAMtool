
RCM_int <- function(OM, RCMdata, condition = c("catch", "catch2", "effort"), selectivity = "logistic", s_selectivity = "B", LWT = list(),
                    comp_like = c("multinomial", "lognormal"), ESS = c(30, 30), prior = list(),
                    max_F = 3, cores = 1L, integrate = FALSE, mean_fit = FALSE, drop_nonconv = FALSE,
                    drop_highF = FALSE,
                    control = list(iter.max = 2e+05, eval.max = 4e+05), ...) {
  
  dots <- list(...) # can be vul_par, ivul_par, log_rec_dev, log_early_rec_dev, map_vul_par, map_ivul_par, map_log_rec_dev, map_log_early_rec_dev, rescale, plusgroup, resample, OMeff, fix_dome
  if(!is.null(dots$maxF)) max_F <- dots$maxF
  if(length(ESS) == 1) ESS <- rep(ESS, 2)
  
  comp_like <- match.arg(comp_like)
  condition <- match.arg(condition)
  
  dat_update <- check_RCMdata(RCMdata, OM, condition)
  OM <- dat_update$OM
  RCMdata <- dat_update$RCMdata
  StockPars <- dat_update$StockPars
  FleetPars <- dat_update$FleetPars
  ObsPars <- dat_update$ObsPars
  
  nsim <- OM@nsim
  proyears <- OM@proyears
  maxage <- OM@maxage
  nyears <- RCMdata@Misc$nyears
  nfleet <- RCMdata@Misc$nfleet
  nsurvey <- RCMdata@Misc$nsurvey
  
  OM@maxF <- max_F
  message("OM@maxF updated to ", max_F, ".")
  
  # No comp data
  if(!any(RCMdata@CAA > 0, na.rm = TRUE) && !any(RCMdata@CAL > 0, na.rm = TRUE)) {
    fix_sel <- TRUE
    message("No fishery length or age compositions were provided. Selectivity is fixed to values from OM.")
  } else {
    fix_sel <- FALSE
  }
  
  # Selectivity
  if(length(selectivity) == 1) selectivity <- rep(selectivity, RCMdata@Misc$nsel_block)
  if(length(selectivity) < RCMdata@Misc$nsel_block) stop("selectivity vector should be of length ", RCMdata@Misc$nsel_block, ").", call. = FALSE)
  sel <- int_sel(selectivity)
  
  # Survey selectivity
  if(nsurvey > 0) {
    if(is.null(s_selectivity)) s_selectivity <- rep("B", nsurvey)
    if(length(s_selectivity) == 1) s_selectivity <- rep(s_selectivity, nsurvey)
    s_sel <- int_s_sel(s_selectivity, nfleet)
  } else {
    s_sel <- int_s_sel("B")
  }
  
  # Likelihood weights
  RCMdata@Misc$LWT <- make_LWT(LWT, nfleet, nsurvey)
  
  # SR
  message(ifelse(OM@SRrel == 1, "Beverton-Holt", "Ricker"), " stock-recruitment relationship used.")
  
  # Generate priors
  prior <- make_prior(prior, nsurvey, OM@SRrel, dots)
  
  # Test for identical sims
  par_identical_sims <- par_identical_sims_fn(StockPars, FleetPars, ObsPars, RCMdata, dots)
  
  # Fit model
  if(!is.null(dots$resample) && dots$resample) { # Re-sample covariance matrix
    
    if(!requireNamespace("mvtnorm", quietly = TRUE)) stop("Please install the mvtnorm package.", call. = FALSE)
    
    message("\nResample = TRUE. Running mean fit model first...")
    mean_fit_output <- RCM_est(RCMdata = RCMdata, selectivity = sel, s_selectivity = s_sel,
                               LWT = RCMdata@Misc$LWT, comp_like = comp_like, prior = prior, 
                               max_F = max_F, integrate = integrate, StockPars = StockPars, ObsPars = ObsPars,
                               FleetPars = FleetPars, mean_fit = TRUE, dots = dots)
    
    if(length(mean_fit_output) > 0 && !mean_fit_output$report$conv) {
      warning("Mean fit model did not appear to converge. Will not be able to sample the covariance matrix.")
      message("Model did not converge. Returning the mean-fit model for evaluation.")
      
      drop_nonconv <- TRUE
      
      samps <- mean_fit_output$obj$env$last.par.best %>% matrix(nrow = nsim, ncol = length(mean_fit_output$obj$par), byrow = TRUE)
    } else {
      message("Model converged. Sampling covariance matrix for nsim = ", nsim, " replicates...")
      if(!all(par_identical_sims)) {
        message("Note: not all ", nsim, " replicates are identical.")
        message("These parameters should be identical among simulations: ",
                which(!par_identical_sims) %>% names() %>% paste0(collapse = " "))
      }
      samps <- mvtnorm::rmvnorm(nsim, mean_fit_output$opt$par, mean_fit_output$SD$cov.fixed)
    }
    
    report_internal_fn <- function(x, samps, obj, conv) {
      report <- obj$report(samps[x, ]) %>% RCM_posthoc_adjust(obj, par = samps[x, ])
      report$conv <- conv
      return(report)
    }
    
    res <- lapply(1:nsim, report_internal_fn, samps = samps, obj = mean_fit_output$obj, conv = mean_fit_output$report$conv)
    mod <- lapply(res, function(x) list(obj = mean_fit_output$obj, report = x))
    conv <- rep(mean_fit_output$report$conv, nsim)
  } else {
    
    if(all(par_identical_sims)) { # All identical sims detected
      message("\nAll ", nsim, " replicates are identical. Fitting once and replicating single fit...")
      
      mean_fit_output <- RCM_est(RCMdata = RCMdata, selectivity = sel, s_selectivity = s_sel,
                                 LWT = RCMdata@Misc$LWT, comp_like = comp_like, prior = prior, 
                                 max_F = max_F, integrate = integrate, StockPars = StockPars, ObsPars = ObsPars,
                                 FleetPars = FleetPars, mean_fit = TRUE, control = control, dots = dots)
      
      mod <- lapply(1:nsim, function(x) return(mean_fit_output))
    } else {
      
      message("\nFitting model (", nsim, " simulations) ...")
      if(cores > 1 && !snowfall::sfIsRunning()) MSEtool::setup(as.integer(cores))
      if(snowfall::sfIsRunning()) {
        mod <- snowfall::sfClusterApplyLB(1:nsim, RCM_est, RCMdata = RCMdata, selectivity = sel, s_selectivity = s_sel,
                                          LWT = RCMdata@Misc$LWT, comp_like = comp_like, prior = prior, 
                                          max_F = max_F, integrate = integrate, StockPars = StockPars, ObsPars = ObsPars,
                                          FleetPars = FleetPars, control = control, dots = dots)
      } else {
        mod <- lapply(1:nsim, RCM_est, RCMdata = RCMdata, selectivity = sel, s_selectivity = s_sel,
                      LWT = RCMdata@Misc$LWT, comp_like = comp_like, prior = prior, 
                      max_F = max_F, integrate = integrate, StockPars = StockPars, ObsPars = ObsPars,
                      FleetPars = FleetPars, control = control, dots = dots)
      }
      
      if(mean_fit) { ### Fit to life history means if mean_fit = TRUE
        message("Generating additional model fit from mean values of parameters in the operating model...\n")
        mean_fit_output <- RCM_est(RCMdata = RCMdata, selectivity = sel, s_selectivity = s_sel,
                                   LWT = RCMdata@Misc$LWT, comp_like = comp_like, prior = prior, 
                                   max_F = max_F, integrate = integrate, StockPars = StockPars, ObsPars = ObsPars,
                                   FleetPars = FleetPars, mean_fit = TRUE, control = control, dots = dots)
        
        if(!mean_fit_output$report$conv) warning("Mean fit model did not appear to converge.")
      } else mean_fit_output <- list()
    }
    res <- lapply(mod, getElement, "report")
    conv <- vapply(res, getElement, logical(1), name = "conv")
  }
  
  if(drop_highF) {
    highF <- vapply(res, function(x) max(getElement(x, "F"), na.rm = TRUE) >= max_F, logical(1))
    if(sum(highF)) message(sum(highF), " out of ", nsim , " model fits had F on the upper boundary (F = ", max_F, "; ", round(100*sum(highF)/nsim, 2), "% of simulations).\n")
    
    conv <- conv & !highF
  }
  if(drop_nonconv) {
    NaF <- vapply(res, function(x) getElement(x, "F") %>% is.na() %>% any(), logical(1))
    if(sum(NaF)) message(sum(NaF), " out of ", nsim , " iterations had F with NA's")
    
    conv <- conv & !NaF
  }
  if(sum(conv) < nsim) message("Non-converged iteration(s): ", paste(which(!conv), collapse = " "), "\n")
  if(sum(conv) < nsim && (drop_nonconv | drop_highF)) {
    message("Non-converged and/or highF iterations will be removed.\n")
    keep <- conv
  } else {
    keep <- !logical(OM@nsim)
  }
  
  ### R0
  OM@cpars$R0 <- vapply(res, getElement, numeric(1), "R0")
  message("Range of unfished age-0 recruitment (OM@cpars$R0): ", paste(round(range(OM@cpars$R0), 2), collapse = " - "))
  
  ### Depletion and init D - init D is only reported, OM setup for initD by adjusting rec devs
  initD <- vapply(res, function(x) x$E[1]/x$E0_SR, numeric(1))
  message("Range of initial spawning depletion: ", paste(round(range(initD), 2), collapse = " - "))
  
  OM@cpars$D <- vapply(res, function(x) x$E[length(x$E)-1]/x$E0_SR, numeric(1))
  message("Range of spawning depletion (OM@cpars$D): ", paste(round(range(OM@cpars$D), 2), collapse = " - "), "\n")
  
  ### Selectivity and F
  ### Find
  OM@isRel <- FALSE
  
  make_F <- function(x) { # Extra step to avoid apical F = 0
    apicalF <- x$F
    apicalF[apicalF < 1e-4] <- 1e-4
    F_at_age <- lapply(1:ncol(apicalF), function(xx) apicalF[, xx] * x$vul[1:nyears, , xx]) %>% 
      simplify2array() %>% apply(1:2, sum)
    return(F_at_age)
  }
  F_matrix <- lapply(res, make_F)
  apical_F <- lapply(F_matrix, function(x) apply(x, 1, max))
  
  expand_V_matrix <- function(x) {
    y <- matrix(x[nyears, ], proyears, maxage + 1, byrow = TRUE)
    rbind(x, y)
  }
  V <- Map("/", e1 = F_matrix, e2 = apical_F) %>% lapply(expand_V_matrix)
  OM@cpars$V <- simplify2array(V) %>% aperm(c(3, 2, 1))
  OM@cpars$Find <- do.call(rbind, apical_F)
  message("Historical F and selectivity trends set in OM@cpars$Find and OM@cpars$V, respectively.")
  message("Selectivity during projection period is set to that in most recent historical year.")
  
  OM@cpars$qs <- rep(1, nsim)
  Eff <- apply(OM@cpars$Find, 2, range)
  OM@EffLower <- Eff[1, ]
  OM@EffUpper <- Eff[2, ]
  if(length(OM@EffYears) != nyears) OM@EffYears <- 1:nyears
  if(length(OM@Esd) == 0 && is.null(OM@cpars$Esd)) OM@Esd <- c(0, 0)
  message("Historical effort trends set in OM@EffLower and OM@EffUpper.\n")
  
  ### Rec devs
  OM@cpars$Perr <- StockPars$procsd
  message("Recruitment standard deviation set in OM@cpars$Perr.")
  
  make_Perr <- function(x) {
    bias_corr <- ifelse(x$obj$env$data$est_rec_dev, exp(-0.5 * x$report$tau^2), 1)
    res <- exp(x$report$log_rec_dev) * bias_corr
    res[1] <- res[1] * x$report$R_eq/x$report$R0
    return(res)
  }
  Perr <- do.call(rbind, lapply(mod, make_Perr))
  
  make_early_Perr <- function(x) {
    res <- x$report$R_eq * x$report$NPR_equilibrium / x$report$R0 / x$report$NPR_unfished[1, ]
    bias_corr <- ifelse(x$obj$env$data$est_early_rec_dev, exp(-0.5 * x$report$tau^2), 1)
    early_dev <- exp(x$report$log_early_rec_dev) * bias_corr
    out <- res[-1] * early_dev
    return(rev(out))
  }
  early_Perr <- do.call(rbind, lapply(mod, make_early_Perr))
  
  OM@cpars$Perr_y <- StockPars$Perr_y
  OM@cpars$Perr_y[, 1:OM@maxage] <- early_Perr
  OM@cpars$Perr_y[, (OM@maxage+1):(OM@maxage + nyears)] <- Perr
  message("Historical recruitment set in OM@cpars$Perr_y.")
  
  log_rec_dev <- do.call(rbind, lapply(res, getElement, "log_rec_dev"))
  
  if(!all(log_rec_dev == 0)) {
    OM@cpars$AC <- apply(log_rec_dev, 1, function(x) {
      out <- acf(x, lag.max = 1, plot = FALSE)$acf[2]
      ifelse(is.na(out), 0, out)
    })
    OM@AC <- range(OM@cpars$AC)
    message("Range of recruitment autocorrelation OM@AC: ", paste(round(range(OM@AC), 2), collapse = " - "))
    
    sample_future_dev <- function() {
      if(!is.null(dots$map_log_rec_dev) && any(is.na(dots$map_log_rec_dev))) { # Sample historical rec devs for OM for most recent years
        yr_fixed_rec_dev <- which(is.na(dots$map_log_rec_dev))
        if(any(yr_fixed_rec_dev > max(dots$map_log_rec_dev, na.rm = TRUE))) {
          yr_hist_sample <- yr_fixed_rec_dev[yr_fixed_rec_dev > max(dots$map_log_rec_dev, na.rm = TRUE)] %>% min()
          samp_hist <- Map(dev_AC, AC = OM@cpars$AC, stdev = StockPars$procsd, chain_start = log_rec_dev[, yr_hist_sample - 1],
                           MoreArgs = list(n = length(yr_hist_sample:nyears), mu = 1))
          log_rec_dev[, yr_hist_sample:nyears] <<- do.call(rbind, samp_hist)
          OM@cpars$Perr_y[, (OM@maxage + yr_hist_sample):(OM@maxage + nyears)] <<- exp(log_rec_dev[, yr_hist_sample:nyears])
          message("Historical recruitment deviations sampled with autocorrelation starting in year ", yr_hist_sample, " out of OM@nyears = ", nyears)
        }
      }
      samp_proj <- Map(dev_AC, AC = OM@cpars$AC, stdev = StockPars$procsd, chain_start = log_rec_dev[, nyears],
                       MoreArgs = list(n = proyears, mu = 1))
      return(exp(do.call(rbind, samp_proj)))
    }
    OM@cpars$Perr_y[, (OM@maxage+nyears+1):ncol(OM@cpars$Perr_y)] <- sample_future_dev()
    message("Future recruitment deviations sampled with autocorrelation (in OM@cpars$Perr_y).\n")
  }
  
  ### Assign OM variables that were used in the RCM to the output
  OM@cpars$Len_age <- StockPars$Len_age
  OM@cpars$Linf <- StockPars$Linf
  OM@cpars$K <- StockPars$K
  OM@cpars$t0 <- StockPars$t0
  OM@cpars$LenCV <- StockPars$LenCV
  OM@cpars$LatASD <- StockPars$LatASD
  OM@cpars$Wt_age <- StockPars$Wt_age
  
  if(any(apply(StockPars$Mat_age, 1, function(x) all(x >= 0.5)))) { # Any simulations where all mat_age > 0.5?
    OM@cpars$L50 <- StockPars$L50
    OM@cpars$L95 <- StockPars$L95
  } else {
    OM@cpars$Mat_age <- StockPars$Mat_age
  }
  if(prior$use_prior[2]) { # May not be necessary
    OM@cpars$h <- vapply(res, getElement, numeric(1), "h")
  } else {
    OM@cpars$h <- StockPars$hs
  }
  if(prior$use_prior[3]) {
    OM@cpars$M_ageArray <- vapply(res, getElement, numeric(1), "Mest") %>% 
      array(c(nsim, maxage+1, nyears + proyears))
  } else {
    OM@cpars$M_ageArray <- StockPars$M_ageArray
  }
  
  if(any(RCMdata@CAL > 0, na.rm = TRUE) || (any(RCMdata@MS > 0, na.rm = TRUE) & RCMdata@MS_type == "length") ||
     any(RCMdata@IAL > 0, na.rm = TRUE)) {
    bw <- RCMdata@length_bin[2] - RCMdata@length_bin[1]
    OM@cpars$CAL_bins <- c(RCMdata@length_bin - 0.5 * bw, max(RCMdata@length_bin) + 0.5 * bw)
    message("RCM length bins will be added to OM.")
  }
  
  if(!is.null(dots$plusgroup) && !dots$plusgroup) OM@cpars$plusgroup <- 0L
  if(!any(RCMdata@I_sd > 0, na.rm = TRUE)) OM@cpars$Iobs <- ObsPars$Iobs
  message("Growth, maturity, natural mortality, and steepness values from RCM are set in OM@cpars.\n")
  
  ### Output list
  E <- do.call(rbind, lapply(res[keep], getElement, "E"))
  N <- array(sapply(res[keep], getElement, "N"), c(nyears+1, maxage+1, sum(keep)))
  CAA_pred <- array(sapply(res[keep], getElement, "CAApred"), c(nyears, maxage+1, nfleet, sum(keep)))
  CAL_pred <- array(sapply(res[keep], getElement, "CALpred"), c(nyears, length(RCMdata@length_bin), nfleet, sum(keep)))
  
  RCMdata@Misc$prior <- prior
  output <- new("RCModel", OM = MSEtool::SubCpars(OM, keep), SSB = E, NAA = aperm(N, c(3, 1, 2)), CAA = aperm(CAA_pred, c(4, 1:3)),
                CAL = aperm(CAL_pred, c(4, 1:3)), mean_fit = mean_fit_output, conv = conv[keep], data = RCMdata, Misc = res[keep])
  
  #if(any(!keep)) {
  #  output@data$drop_sim <- which(!keep)
  #} else {
  #  output@data$drop_sim <- numeric(0)
  #}
  
  # Data in cpars
  if(sum(RCMdata@Chist > 0, na.rm = TRUE) || nsurvey > 0) {
    real_Data <- new("Data")
    real_Data@Year <- (output@OM@CurrentYr - output@OM@nyears + 1):output@OM@CurrentYr
    real_Data@LHYear <- max(real_Data@Year)
    real_Data@MaxAge <- maxage
    if(sum(RCMdata@Chist > 0, na.rm = TRUE) && all(!is.na(RCMdata@Chist))) {
      real_Data@Cat <- matrix(rowSums(RCMdata@Chist, na.rm = TRUE), 1, nyears)
      real_Data@CV_Cat <- sqrt(exp((rowSums(RCMdata@C_sd * RCMdata@Chist)/rowSums(RCMdata@Chist))^2 - 1)) %>% 
        matrix(1, nyears)
      message("Historical catch data added to OM@cpars$Data@Cat.")
      if(nfleet > 1) message("Annual catch CV is a weighted average by fleet-specific catch.")
    }
    if(nfleet == 1) {
      if(sum(RCMdata@CAA, na.rm = TRUE)) {
        real_Data@CAA <- aperm(RCMdata@CAA, c(3, 1, 2))
        message("Age comps added to OM@cpars$Data@CAA (nfleet = 1).")
      }
      if(sum(RCMdata@CAL, na.rm = TRUE)) {
        real_Data@CAL <- aperm(RCMdata@CAL, c(3, 1, 2))
        real_Data@CAL_mids <- RCMdata@length_bin
        real_Data@CAL_bins <- c(RCMdata@length_bin - 0.5 * bw, max(RCMdata@length_bin) + 0.5 * bw)
        message("Length comps added to OM@cpars$Data@CAL (nfleet = 1).")
      }
    }
    if(.hasSlot(real_Data, "AddInd") && nsurvey > 0) {
      real_Data@AddInd <- array(RCMdata@Index, c(nyears, nsurvey, output@OM@nsim)) %>%
        aperm(perm = c(3, 2, 1))
      real_Data@CV_AddInd <- array(sqrt(exp(RCMdata@I_sd^2) - 1), c(nyears, nsurvey, output@OM@nsim)) %>%
        aperm(perm = c(3, 2, 1))
      
      # Cannot accommodate indices mirrored to fleet when nfleet > 1 and fleet has time-varying sel
      real_Data@AddIndType <- vapply(s_sel, process_AddIndType, numeric(1), nfleet = nfleet)
      real_Data@AddIndV <- lapply(1:nsurvey, process_AddIndV, Misc = output@Misc, s_sel = s_sel,
                                  n_age = maxage + 1, nfleet = nfleet, nyears = nyears) %>%
        simplify2array() %>% aperm(c(1, 3, 2))
      real_Data@AddIunits <- RCMdata@I_units
      message("Historical indices added to OM@cpars$Data@AddInd.")
    }
    output@OM@cpars$Data <- real_Data
  }
  
  # Check whether observed matches predicted
  if(RCMdata@Misc$condition == "catch2") {
    catch_check_fn <- function(x, report, data) {
      if(report[[x]]$conv) {
        catch_diff <- report[[x]]$Cpred/RCMdata@Chist - 1
        flag <- max(abs(catch_diff), na.rm = TRUE) > 0.01 | any(is.na(report[[x]]$Cpred))
      } else {
        flag <- FALSE
      }
      return(flag)
    }
    
    do_catch_check <- vapply(1:length(res), catch_check_fn, logical(1), report = res, data = data)
    if(any(do_catch_check)) {
      flag_ind <- paste(which(do_catch_check), collapse = " ")
      message("Note: there is predicted catch that deviates from observed catch by more than 1% in simulations:")
      message(flag_ind)
    }
  }
  message("Complete.")
  
  return(output)
}


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
                   IAL_hist = RCMdata@IAL, IAL_n = RCMdata@IAL_ESS, length_bin = RCMdata@length_bin, 
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
  length_bin <- obj$env$data$length_bin
  
  report$CR <- report$Arec * report$EPR0
  
  age_only_model <- data$len_age %>%
    apply(1, function(x) length(x) == data$n_age && max(x) == data$n_age) %>% all()
  if(age_only_model) {
    report$vul_len <- matrix(NA_real_, length(length_bin), data$nsel_block)
    report$ivul_len <- matrix(NA_real_, length(length_bin), dim(report$ivul)[3])
    
    report$MLpred <- matrix(NA_real_, nrow(report$F), ncol(report$F))
    report$CALpred <- array(NA_real_, dim(report$CALpred))
    report$IALpred <- array(NA_real_, dim(report$IALpred))
  } else {
    report$vul_len <- get_vul_len(report, data$vul_type, length_bin, data$Linf)
    report$ivul_len <- get_ivul_len(report, data$ivul_type, length_bin, data$Linf)
  }
  if(dynamic_SSB0) report$dynamic_SSB0 <- RCM_dynamic_SSB0(obj, par)
  return(report)
}

get_vul_len <- function(report, selectivity, length_bin, Linf) {
  vul <- matrix(NA, length(length_bin), length(selectivity))
  sel_ind  <- selectivity == 0 | selectivity == -1
  LFS <- report$LFS[sel_ind]
  L5 <- report$L5[sel_ind]
  Vmaxlen <- report$Vmaxlen[sel_ind]
  
  sls <- (LFS - L5)/sqrt(-log(0.05, 2))
  srs <- (Linf - LFS)/sqrt(-log(Vmaxlen, 2))
  asc <- Map(function(x, y) 2^-((length_bin - y)/x * (length_bin - y)/x), x = sls, y = LFS)
  dsc <- Map(function(x, y, z) ifelse(z > rep(0.99, length(length_bin)), 1, 2^-((length_bin - y)/x * (length_bin - y)/x)),
             x = srs, y = LFS, z = Vmaxlen)
  vul_out <- Map(function(x, y, z) ifelse(length_bin > x, y, z), x = LFS, y = dsc, z = asc)
  vul[, sel_ind] <- do.call(cbind, vul_out)
  return(vul)
}

get_ivul_len <- function(report, s_selectivity, length_bin, Linf) {
  ivul_len <- matrix(NA, length(length_bin), length(s_selectivity)) # length-based: matrix of dimension nlbin, nsurvey
  for(i in 1:ncol(ivul_len)) {
    if(s_selectivity[i] == -1 || s_selectivity[i] == 0) {
      sls <- (report$iLFS[i] - report$iL5[i])/sqrt(-log(0.05, 2))
      srs <- (Linf - report$iLFS[i])/sqrt(-log(report$iVmaxlen[i], 2))
      
      asc <- 2^-((length_bin - report$iLFS[i])/sls * (length_bin - report$iLFS[i])/sls)
      dsc <- ifelse(report$iVmaxlen[i] > rep(0.99, length(length_bin)), 1,
                    2^-((length_bin - report$iLFS[i])/srs * (length_bin - report$iLFS[i])/srs))
      ivul_len[, i] <- ifelse(length_bin > report$iLFS[i], dsc, asc)
    } else if(s_selectivity[i] > 0) {
      ivul_len[, i] <- report$vul_len[, s_selectivity[i]]
    }
  }
  return(ivul_len)
}

process_AddIndType <- function(s_sel, nfleet) {
  if(s_sel == -4 || s_sel == -2 || s_sel == -1 || s_sel == 0) { # -4 = B,-2 - 0 = custom sel
    return(1)
  } else if(s_sel == -3) { # SSB
    return(2)
  } else { # Fleet
    return(ifelse(nfleet > 1, 1, 3))
  }
}

process_AddIndV <- function(sur, Misc, s_sel, n_age, nfleet, nyears) { # Return a matrix of nsim x nages
  if(s_sel[sur] < -2 || (s_sel[sur] == 1 & nfleet == 1)) { # -4 = B, -3 = SSB, single-fleet VB
    out <- matrix(1, length(Misc), n_age)
  } else { # custom sel or multi-fleet
    out <- do.call(rbind, lapply(Misc, function(x) x$ivul[nyears, , sur]))
  }
  return(out)
}

RCM_retro <- function(x, nyr = 5) {
  if(length(x@mean_fit) == 0) stop("Re-run RCM() with argument `mean_fit = TRUE`", .call = FALSE)
  
  data <- x@mean_fit$obj$env$data
  params <- x@mean_fit$obj$env$parameters
  n_y <- data$n_y
  map <- x@mean_fit$obj$env$map
  
  if(data$nfleet > 1) {
    retro_ts <- array(NA, dim = c(nyr+1, n_y + 1, data$nfleet + 4))
    TS_var <- c(paste("Fleet", 1:data$nfleet, "F"), "Apical F", "SSB", "SSB_SSB0", "R")
  } else {
    retro_ts <- array(NA, dim = c(nyr+1, n_y + 1, 4))
    TS_var <- c("Apical F", "SSB", "SSB_SSB0", "R")
  }
  dimnames(retro_ts) <- list(Peel = 0:nyr, Year = (x@OM@CurrentYr - n_y):x@OM@CurrentYr + 1, Var = TS_var)
  
  new_args <- lapply(n_y - 0:nyr, RCM_retro_subset, data = data, params = params, map = map)
  lapply_fn <- function(i, new_args, x) {
    obj2 <- MakeADFun(data = new_args[[i+1]]$data, parameters = new_args[[i+1]]$params, map = new_args[[i+1]]$map,
                      random = x@mean_fit$obj$env$random, DLL = "SAMtool", silent = TRUE)
    if(new_args[[i+1]]$data$condition == "catch2") {
      if(any(is.na(obj2$report(obj2$par)$F)) || any(is.infinite(obj2$report(obj2$par)$F))) {
        for(ii in 1:10) {
          obj2$par["R0x"] <- 0.5 + obj2$par["R0x"]
          if(all(!is.na(obj2$report(obj2$par)$F)) && all(!is.infinite(obj2$report(obj2$par)$F))) break
        }
      }
    }
    mod <- optimize_TMB_model(obj2, control = list(iter.max = 2e+05, eval.max = 4e+05), restart = 0)
    opt2 <- mod[[1]]
    SD <- mod[[2]]
    
    if(!is.character(opt2)) {
      report <- obj2$report(obj2$env$last.par.best) %>% RCM_posthoc_adjust(obj2, dynamic_SSB0 = FALSE)
      
      FMort <- rbind(report$F, matrix(NA, i + 1, ncol(report$F)))
      if(data$nfleet > 1) FMort <- cbind(FMort, apply(report$F_at_age, 1, max, na.rm = TRUE) %>% c(rep(NA, i+1)))
      
      SSB <- c(report$E, rep(NA, i))
      SSB_SSB0 <- SSB/report$E0_SR
      R <- c(report$R, rep(NA, i))
      
      retro_ts[i+1, , ] <<- cbind(FMort, SSB, SSB_SSB0, R)
      
      return(SD$pdHess)
    }
    return(FALSE)
  }
  
  conv <- vapply(0:nyr, lapply_fn, logical(1), new_args = new_args, x = x)
  if(any(!conv)) warning("Peels that did not converge: ", paste0(which(!conv) - 1, collapse = " "))
  
  retro <- new("retro", Model = "RCM", Name = x@OM@Name, TS_var = TS_var, TS = retro_ts)
  if(data$nfleet > 1) {  
    attr(retro, "TS_lab") <- c(paste("Fishing mortality of Fleet", 1:data$nfleet), "Apical F",
                               "Spawning biomass", "Spawning depletion", "Recruitment")
  } else {  
    attr(retro, "TS_lab") <- c(paste("Fishing mortality of Fleet", 1:data$nfleet),
                               "Spawning biomass", "Spawning depletion", "Recruitment")
  }
  return(retro)
}


RCM_retro_subset <- function(yr, data, params, map) {
  data_out <- structure(data, check.passed = NULL)
  params_out <- lapply(params, function(x) if(!is.null(attr(x, "map"))) attr(x, "shape") else x)
  
  if(yr < data$n_y) {
    ##### Data object
    mat <- c("C_hist", "E_hist", "I_hist", "sigma_I", "sigma_C", "CAA_n", "CAL_n", "IAA_n", "IAL_n", "msize")
    if(!is.null(map$log_M)) mat <- c(mat, "M_data")
    mat_ind <- match(mat, names(data_out))
    data_out[mat_ind] <- lapply(data_out[mat_ind], function(x) x[1:yr, , drop = FALSE])
    
    mat2 <- c("len_age", "wt", "mat", "sel_block")
    mat_ind2 <- match(mat2, names(data_out))
    data_out[mat_ind2] <- lapply(data_out[mat_ind2], function(x) x[1:(yr+1), , drop = FALSE])
    
    # Update nsel_block, n_y
    data_out$nsel_block <- data_out$sel_block %>% as.vector() %>% unique() %>% length()
    data_out$n_y <- yr
    
    # Array 1:yr first index
    arr <- c("CAA_hist", "CAL_hist", "IAA_hist", "IAL_hist")
    arr_ind <- match(arr, names(data_out))
    data_out[arr_ind] <- lapply(data_out[arr_ind], function(x) x[1:yr, , , drop = FALSE])
    
    # Vector 1:yr
    data_out$est_rec_dev <- data_out$est_rec_dev[1:yr]
    
    ##### Parameters
    # Update rec devs
    params_out$log_rec_dev <- params_out$log_rec_dev[1:yr]
    if(!is.null(map$log_rec_dev)) map$log_rec_dev <- map$log_rec_dev[1:yr] %>% factor()
    
    # Update F
    if(data_out$condition == "catch") {
      params_out$log_F_dev <- params_out$log_F_dev[1:yr, , drop = FALSE]
      
      if(any(data_out$yind_F + 1 > yr)) {
        data_out$yind_F <- as.integer(0.5 * yr)
        params_out$log_F_dev[data_out$yind_F + 1, ] <- params$log_F_dev[data$yind_F + 1, ]
      }
    }
    
    ## Update nsel block if needed
    if(data$nsel_block > data_out$nsel_block) {
      sel_block_ind <- data_out$sel_block %>% as.vector() %>% unique()
      params_out$vul_par <- params_out$vul_par[, sel_block_ind, drop = FALSE]
      if(!is.null(map$vul_par)) map$vul_par <- map$vul_par[, sel_block_ind, drop = FALSE] %>% factor()
    }
  }
  
  return(list(data = data_out, params = params_out, map = map))
}


profile_likelihood_RCM <- function(x, ...) {
  dots <- list(...)
  if(!length(x@mean_fit)) stop("No model found. Re-run RCM with mean_fit = TRUE.")
  
  data <- x@mean_fit$obj$env$data
  params <- x@mean_fit$obj$env$parameters
  n_y <- data$n_y
  map <- x@mean_fit$obj$env$map
  LWT <- try(x@data@Misc$LWT, silent = TRUE)
  if(is.character(LWT)) LWT <- x@data$LWT
  
  new_args <- RCM_retro_subset(n_y, data = data, params = params, map = map)
  
  if(!is.null(dots$D)) {
    if(length(dots) > 1) message("Only doing depletion profile...")
    if(!requireNamespace("abind", quietly = TRUE)) {
      stop("Please install the abind package to run this profile.")
    }
    if(!requireNamespace("reshape2", quietly = TRUE)) {
      stop("Please install the reshape2 package to run this profile.")
    }
    profile_grid <- expand.grid(D = dots$D)
    
    new_args$data$IAA_hist <- abind::abind(new_args$data$IAA_hist, 
                                           array(NA, c(n_y, new_args$data$n_age, 1)), 
                                           along = 3)
    new_args$data$IAL_hist <- abind::abind(new_args$data$IAL_hist, 
                                           array(NA, c(n_y, length(new_args$data$length_bin), 1)), 
                                           along = 3)
    new_args$data$IAA_n <- cbind(new_args$data$IAA_n, rep(0, n_y))
    new_args$data$IAL_n <- cbind(new_args$data$IAL_n, rep(0, n_y))
    new_args$data$nsurvey <- new_args$data$nsurvey + 1
    
    new_args$data$ivul_type <- c(new_args$data$ivul_type, -3)
    new_args$data$abs_I <- c(new_args$data$abs_I, 0)
    new_args$data$I_units <- c(new_args$data$I_units, 1) 
    new_args$data$LWT_index <- rbind(new_args$data$LWT_index, rep(1, 3))
    
    new_args$data$use_prior <- c(new_args$data$use_prior, 0)
    new_args$data$prior_dist <- rbind(new_args$data$prior_dist, rep(NA_real_, 2))
    
    new_args$params$ivul_par <- cbind(new_args$params$ivul_par, rep(0, nrow(new_args$params$ivul_par)))
    if(!is.null(new_args$map$ivul_par)) {
      new_args$map$ivul_par <- factor(c(new_args$map$ivul_par, rep(NA, nrow(new_args$params$ivul_par))))
    }
    
    LWT$Index <- c(LWT$Index, 1)
    LWT$IAA <- c(LWT$IAA, 1)
    LWT$IAL <- c(LWT$IAL, 1)
    
    profile_fn <- function(i, new_args, x) {
      # Add new survey of SSB depletion
      SSB_index <- SSB_Isd <- rep(NA_real_, n_y)
      SSB_index[1] <- 1
      SSB_index[n_y] <- profile_grid$D[i]
      SSB_Isd[c(1, n_y)] <- 0.01
      
      new_args$data$I_hist <- cbind(new_args$data$I_hist, SSB_index)
      new_args$data$sigma_I <- cbind(new_args$data$sigma_I, SSB_Isd)
      
      obj2 <- MakeADFun(data = new_args$data, parameters = new_args$params, map = new_args$map,
                        random = x@mean_fit$obj$env$random, DLL = "SAMtool", silent = TRUE)
      mod <- optimize_TMB_model(obj2, control = list(iter.max = 2e+05, eval.max = 4e+05), restart = 0)
      report <- obj2$report(obj2$env$last.par.best)
      RCM_get_likelihoods(report, LWT = LWT, f_name = paste0("Fleet_", 1:new_args$data$nfleet),
                          s_name = paste0("Index_", 1:(new_args$data$nsurvey-1)) %>% c("SSB_depletion"))
    }
    
    do_profile <- lapply(1:nrow(profile_grid), profile_fn, new_args = new_args, x = x)
    profile_grid$nll <- vapply(do_profile, function(xx) xx[[1]][1, 1], numeric(1))
    
    prof <- sapply(do_profile, nll_depletion_profile) %>% t()
    prof_out <- apply(prof, 2, sum) %>% as.logical()
    
    profile_grid <- cbind(profile_grid, prof[, prof_out] %>% as.data.frame())
    
    output <- new("prof", Model = "RCM", Par = "D", 
                  MLE = x@mean_fit$report$E[n_y]/x@mean_fit$report$E0_SR, grid = profile_grid)
  } else {
    
    if(is.null(dots$R0) && is.null(dots$h)) stop("Sequence of neither D, R0, nor h was found. See help file.")
    if(!is.null(dots$R0)) {
      R0 <- dots$R0 
    } else {
      R0 <- x@mean_fit$report$R0
      profile_par <- "h"
    }
    if(!is.null(dots$h)) {
      h <- dots$h 
    } else {
      h <- x@mean_fit$report$h
      profile_par <- "R0"
    }
    
    profile_grid <- expand.grid(R0 = R0, h = h)
    joint_profile <- !exists("profile_par", inherits = FALSE)
    
    profile_fn <- function(i, new_args, x) {
      new_args$params$R0x <- log(profile_grid$R0[i]  * new_args$data$rescale)
      if(new_args$data$SR_type == "BH") {
        new_args$params$transformed_h <- logit((profile_grid$h[i] - 0.2)/0.8)
      } else {
        new_args$params$transformed_h <- log(profile_grid$h[i] - 0.2)
      }
      
      if(joint_profile) {
        new_args$map$R0x <- new_args$map$transformed_h <- factor(NA)
      } else {
        if(profile_par == "R0") new_args$map$R0x <- factor(NA) else new_args$map$transformed_h <- factor(NA)
      }
      
      obj2 <- MakeADFun(data = new_args$data, parameters = new_args$params, map = new_args$map,
                        random = x@mean_fit$obj$env$random, DLL = "SAMtool", silent = TRUE)
      
      mod <- optimize_TMB_model(obj2, control = list(iter.max = 2e+05, eval.max = 4e+05), restart = 0)
      report <- obj2$report(obj2$env$last.par.best)
      RCM_get_likelihoods(report, LWT = LWT, f_name = paste0("Fleet_", 1:new_args$data$nfleet),
                          s_name = paste0("Index_", 1:new_args$data$nsurvey))
      
      #opt2 <- optimize_TMB_model(obj2, control = list(iter.max = 2e+05, eval.max = 4e+05), restart = 0)[[1]]
      #
      #if(!is.character(opt2)) nll <- opt2$objective else nll <- NA
      #return(nll)
    }
    
    do_profile <- lapply(1:nrow(profile_grid), profile_fn, new_args = new_args, x = x)
    profile_grid$nll <- vapply(do_profile, function(xx) xx[[1]][1, 1], numeric(1))
    
    prof <- sapply(do_profile, nll_depletion_profile) %>% t()
    prof_out <- apply(prof, 2, sum) %>% as.logical()
    
    profile_grid <- cbind(profile_grid, prof[, prof_out] %>% as.data.frame())
    
    #nll <- vapply(1:nrow(profile_grid), profile_fn, numeric(1), new_args = new_args, x = x)
    #profile_grid$nll <- nll
    
    if(joint_profile) {
      pars <- c("R0", "h")
      MLE <- c(x@mean_fit$report$R0, x@mean_fit$report$h)
    } else {
      pars <- profile_par
      MLE <- getElement(x@mean_fit$report, pars)
    }
    
    output <- new("prof", Model = "RCM", Par = pars, MLE = MLE, grid = profile_grid)
  }
  return(output)
}

#' @importFrom dplyr select starts_with
nll_depletion_profile <- function(xx) {
  nll_fleet <- xx[[2]][-nrow(xx[[2]]), ] %>% select(starts_with("Fleet"))
  nll_fleet$Data <- rownames(nll_fleet)
  nll_fleet <- reshape2::melt(nll_fleet, id.var = "Data", value.var = "nll")
  
  nll_index <- xx[[4]][-nrow(xx[[4]]), ] %>% select(starts_with("Index"))
  nll_index$Data <- rownames(nll_index)
  nll_index <- reshape2::melt(nll_index, id.var = "Data", value.var = "nll")
  
  nll <- rbind(nll_fleet, nll_index)
  
  nll$Name <- paste(nll$variable, nll$Data)
  
  c(nll$value, xx[[1]][c(2, 5, 6), 1]) %>% 
    structure(names = paste(nll$variable, nll$Data) %>% c("Recruitment Deviations", "Penalty", "Prior"))
}
