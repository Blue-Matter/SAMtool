
#' @importFrom pbapply pblapply
RCM_int <- function(OM, RCMdata, condition = "catch", selectivity = "logistic", s_selectivity = "B", LWT = list(),
                    comp_like = c("multinomial", "lognormal", "mvlogistic", "dirmult1", "dirmult2"), prior = list(),
                    max_F = 3, cores = 1L, integrate = FALSE, mean_fit = FALSE, drop_nonconv = FALSE,
                    drop_highF = FALSE,
                    control = list(iter.max = 2e+05, eval.max = 4e+05), 
                    start = list(), map = list(), silent = FALSE, ...) {
  
  dots <- list(...) # can be vul_par, ivul_par, log_rec_dev, log_early_rec_dev, map_vul_par, map_ivul_par, map_log_rec_dev, map_log_early_rec_dev, rescale, plusgroup, resample, OMeff, fix_dome
  if (!is.null(dots$maxF)) max_F <- dots$maxF
  if (!is.null(dots$resample) && dots$resample && !requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Found argument: resample = TRUE. Please install the mvtnorm package.", call. = FALSE) # Re-sample covariance matrix
  }
  
  comp_like <- match.arg(comp_like)
  condition <- match.arg(condition, choices = c("catch", "catch2", "effort"), several.ok = TRUE)
  
  dat_update <- check_RCMdata(RCMdata, OM, condition, silent = silent)
  OM <- dat_update$OM
  RCMdata <- dat_update$RCMdata
  StockPars <- dat_update$StockPars
  FleetPars <- dat_update$FleetPars
  
  nsim <- OM@nsim
  proyears <- OM@proyears
  maxage <- OM@maxage
  nyears <- RCMdata@Misc$nyears
  nfleet <- RCMdata@Misc$nfleet
  nsurvey <- RCMdata@Misc$nsurvey
  
  OM@maxF <- max_F
  if (!silent) message("Maximum F in RCM will be ", max_F, ". OM@maxF is also updated.\n\n")
  
  # No comp data
  if (!any(RCMdata@CAA > 0, na.rm = TRUE) && !any(RCMdata@CAL > 0, na.rm = TRUE)) {
    fix_sel <- TRUE
    if (!silent) message_info("No fishery length or age compositions were provided. Selectivity is fixed to values from OM.\n\n")
  } else {
    fix_sel <- FALSE
  }
  
  # Selectivity
  if (length(selectivity) == 1) selectivity <- rep(selectivity, RCMdata@Misc$nsel_block)
  if (length(selectivity) < RCMdata@Misc$nsel_block) stop("selectivity vector should be of length ", RCMdata@Misc$nsel_block, ").", call. = FALSE)
  sel <- int_sel(selectivity, RCMdata, silent)
  
  # Survey selectivity
  if (nsurvey > 0) {
    if (is.null(s_selectivity)) s_selectivity <- rep("B", nsurvey)
    if (length(s_selectivity) == 1) s_selectivity <- rep(s_selectivity, nsurvey)
    if (length(s_selectivity) != nsurvey) stop("Length of s_selectivity is not equal to ", nsurvey)
    s_sel <- int_s_sel(s_selectivity, nfleet, silent)
  } else {
    s_sel <- int_s_sel(NULL, silent = silent)
  }
  
  # Likelihood weights
  RCMdata@Misc$LWT <- make_LWT(LWT, nfleet, nsurvey)
  
  # SR
  if (!silent) {
    message(switch(OM@SRrel,
                   "1" = "Beverton-Holt", 
                   "2" = "Ricker",
                   "3" = "Mesnil-Rochet"), " stock-recruitment relationship used.")
  }
  
  # Generate priors
  prior <- make_prior(prior, nsurvey, OM@SRrel, dots, msg = !silent)
  
  # Test for identical sims
  par_identical_sims <- par_identical_sims_fn(StockPars, FleetPars, RCMdata, dots)
  
  # Fit model
  if (all(par_identical_sims)) { # All identical sims detected
      
    if (!silent) message_info("All ", nsim, " replicates are identical. Fitting one model...")
    
    mean_fit_output <- RCM_est(RCMdata = RCMdata, selectivity = sel, s_selectivity = s_sel,
                               LWT = RCMdata@Misc$LWT, comp_like = comp_like, prior = prior, 
                               max_F = max_F, integrate = integrate, StockPars = StockPars,
                               FleetPars = FleetPars, mean_fit = TRUE, control = control,
                               start = start, map = map, dots = dots)
    
    if (mean_fit_output$report$conv) {
      
      if (!is.null(dots$resample) && dots$resample) { # Re-sample covariance matrix
        message_info("Sampling covariance matrix for nsim = ", nsim, " replicates...")

        samps <- mvtnorm::rmvnorm(nsim, mean_fit_output$opt$par, mean_fit_output$SD$cov.fixed)
        
        if (cores > 1 && !snowfall::sfIsRunning()) MSEtool::setup(as.integer(cores))
        res <- pblapply(1:nsim, RCM_report_samps, samps = samps, obj = mean_fit_output$obj, conv = mean_fit_output$report$conv,
                        cl = if (snowfall::sfIsRunning()) snowfall::sfGetCluster() else NULL)
      } else {
        if (!silent) message("Model converged and will be replicated for all simulations.")
        res <- list(mean_fit_output$report) 
      }
      
    } else {
      warning("Model did not converge. Returning the mean-fit model for evaluation.")
      if (!is.null(dots$resample) && dots$resample) warning("Could not sample covariance matrix.")
      res <- list(mean_fit_output$report) 
    }
    
    mod <- list(mean_fit_output) # Length 1
    conv <- rep(mean_fit_output$report$conv, nsim)
    
  } else {
    
    message_info("Fitting model (", nsim, " simulations) ...")
    if (cores > 1 && !snowfall::sfIsRunning()) MSEtool::setup(as.integer(cores))
    
    mod <- pblapply(1:nsim, RCM_est, RCMdata = RCMdata, selectivity = sel, s_selectivity = s_sel,
                    LWT = RCMdata@Misc$LWT, comp_like = comp_like, prior = prior, 
                    max_F = max_F, integrate = integrate, StockPars = StockPars,
                    FleetPars = FleetPars, control = control, start = start, map = map, dots = dots,
                    cl = if (snowfall::sfIsRunning()) snowfall::sfGetCluster() else NULL)
    
    if (mean_fit) { ### Fit to life history means if mean_fit = TRUE
      if (!silent) message_info("Generating additional model fit from mean values of parameters in the operating model...\n")
      mean_fit_output <- RCM_est(RCMdata = RCMdata, selectivity = sel, s_selectivity = s_sel,
                                 LWT = RCMdata@Misc$LWT, comp_like = comp_like, prior = prior, 
                                 max_F = max_F, integrate = integrate, StockPars = StockPars,
                                 FleetPars = FleetPars, mean_fit = TRUE, control = control, 
                                 start = start, map = map, dots = dots)
      
      if (!mean_fit_output$report$conv) warning("Mean fit model did not appear to converge.")
    } else {
      mean_fit_output <- list()
    }
    
    conv <- vapply(mod, function(x) x[["report"]][["conv"]], logical(1))
    if (!is.null(dots$resample) && dots$resample) {
      
      if (!silent) message_info("Sampling covariance matrix once for each replicate...")
      res <- pblapply(1:nsim, function(x) {
        if (conv[x]) {
          samps <- mvtnorm::rmvnorm(1, mod[[x]]$opt$par, mod[[x]]$SD$cov.fixed, checkSymmetry = FALSE)
          RCM_report_samps(1, samps = samps, obj = mod[[x]]$obj, conv = conv[x])
        } else {
          mod[[x]][["report"]]
        } 
      }, cl = NULL)
      
    } else {
      res <- lapply(mod, getElement, "report")
    }
  }
  
  highF <- vapply(res, function(x) max(getElement(x, "F"), na.rm = TRUE) >= max_F, logical(1))
  if (length(res) == 1) {
    if (highF) warning("Model had F on the upper boundary.\n")
  } else if (sum(highF)) {
    msg <- paste0(sum(highF), " out of ", nsim , " model fits had F on the upper boundary (F = ", 
                  max_F, "; ", round(100 * mean(highF), 2), "% of simulations).")
    warning(msg)
    if (drop_highF) conv <- conv & !highF
  }
  
  NaF <- vapply(res, function(x) any(is.na(x$F) | is.infinite(x$F)), logical(1))
  if (length(res) == 1) {
    if (NaF) warning("Model had F with NA's")
  } else if (sum(NaF)) {
    msg <- paste0(sum(NaF), " out of ", nsim , " iterations (", round(100 * mean(NaF), 2), "%) had F with NA's")
    warning(msg)
    if (drop_nonconv) conv <- conv & !NaF
  }
  
  keep <- !logical(OM@nsim) # Keep all simulations
  if (length(res) > 1 && sum(conv) < nsim) {
    if (!silent) message_oops("Non-converged iteration(s): ", paste(which(!conv), collapse = " "), "\n\n")
    if (!sum(conv)) {
      if (!silent) message_oops("Non-converged for all iterations. Returning all for evaluation.\n\n")
    } else if (drop_nonconv || drop_highF) {
      if (!silent) message_info("Non-converged and/or highF iterations will be removed.\n\n")
      keep <- conv
    }
  }
  
  ### Get parameters to update OM
  obj_data <- mod[[1]]$obj$env$data
  newOM <- RCM_update_OM(OM, res, StockPars, obj_data, maxage, nyears, proyears, nsim, prior, silent)
  
  ### Output S4 object
  RCMdata@Misc$prior <- prior
  output <- new("RCModel", 
                OM = MSEtool::SubCpars(newOM$OM, keep), 
                SSB = newOM$RCM_val$SSB[keep, , drop = FALSE], 
                NAA = newOM$RCM_val$NAA[keep, , , drop = FALSE], 
                CAA = newOM$RCM_val$CAA[keep, , , , drop = FALSE], 
                mean_fit = mean_fit_output, 
                conv = conv[keep], 
                data = RCMdata)
  if (!is.null(newOM$RCM_val$CAL)) output@CAL = newOM$RCM_val$CAL[keep, , , , drop = FALSE] 
  
  if (length(res) > 1) {
    output@Misc <- res[keep]
    output@config <- list(drop_sim = which(!keep))
  } else {
    output@Misc <- res
    output@config <- list(drop_sim = integer(0))
  }
  
  # Data in cpars
  if (sum(RCMdata@Chist > 0, na.rm = TRUE) || nsurvey > 0) {
    
    if (!silent) message_info("Adding some RCMdata inputs into OM@cpars$Data:\n\n")
    real_Data <- new("Data")
    real_Data@Year <- (output@OM@CurrentYr - output@OM@nyears + 1):output@OM@CurrentYr
    real_Data@LHYear <- max(real_Data@Year)
    real_Data@MaxAge <- maxage
    if (sum(RCMdata@Chist > 0, na.rm = TRUE) && all(!is.na(RCMdata@Chist))) {
      real_Data@Cat <- matrix(rowSums(RCMdata@Chist, na.rm = TRUE), 1, nyears)
      real_Data@CV_Cat <- sqrt(exp((rowSums(RCMdata@C_sd * RCMdata@Chist)/rowSums(RCMdata@Chist))^2 - 1)) %>% 
        matrix(1, nyears)
      
      if (!silent) {
        message("Historical catch data added to OM@cpars$Data@Cat.")
        if (nfleet > 1) message("Annual catch CV is a weighted average by fleet-specific catch.")
      }
      
    }
    if (nfleet == 1) {
      if (sum(RCMdata@CAA, na.rm = TRUE)) {
        real_Data@CAA <- aperm(RCMdata@CAA, c(3, 1, 2))
        if (!silent) message("Age comps added to OM@cpars$Data@CAA (nfleet = 1).")
      }
      if (sum(RCMdata@CAL, na.rm = TRUE)) {
        real_Data@CAL <- aperm(RCMdata@CAL, c(3, 1, 2))
        real_Data@CAL_mids <- RCMdata@Misc$lbinmid
        real_Data@CAL_bins <- RCMdata@Misc$lbin
        if (!silent) message("Length comps added to OM@cpars$Data@CAL (nfleet = 1).")
      }
    }
    if (.hasSlot(real_Data, "AddInd") && nsurvey > 0) {
      real_Data@AddInd <- array(RCMdata@Index, c(nyears, nsurvey, 1)) %>%
        aperm(perm = c(3, 2, 1)) %>% structure(dimnames = list(1, colnames(RCMdata@Index), real_Data@Year))
      real_Data@CV_AddInd <- array(sqrt(exp(RCMdata@I_sd^2) - 1), c(nyears, nsurvey, 1)) %>%
        aperm(perm = c(3, 2, 1))
      
      # Cannot accommodate indices mirrored to fleet when nfleet > 1 and fleet has time-varying sel
      real_Data@AddIndType <- vapply(s_sel, process_AddIndType, numeric(1), nfleet = nfleet)
      real_Data@AddIndV <- lapply(1:nsurvey, process_AddIndV, Misc = output@Misc, s_sel = s_sel,
                                  n_age = maxage + 1, nfleet = nfleet, nyears = nyears) %>%
        simplify2array() %>% aperm(c(1, 3, 2))
      real_Data@AddIunits <- RCMdata@I_units
      
      # Override to ensure index beta = 1
      output@OM@cpars$AddIbeta <- matrix(1, output@OM@nsim, nsurvey)
      if (!silent) message("Historical indices added to OM@cpars$Data@AddInd.")
    }
    output@OM@cpars$Data <- real_Data
  }
  
  # Check whether observed matches predicted
  if (any(RCMdata@Misc$condition == "catch2")) {
    catch_check_fn <- function(x, report, RCMdata) {
      if (report[[x]]$conv) {
        catch_diff <- report[[x]]$Cpred/RCMdata@Chist - 1
        catch_diff <- catch_diff[, RCMdata@Misc$condition == "catch2"]
        flag <- max(abs(catch_diff), na.rm = TRUE) > 0.01 || any(is.na(report[[x]]$Cpred))
      } else {
        flag <- FALSE
      }
      return(flag)
    }
    
    do_catch_check <- vapply(1:length(res), catch_check_fn, logical(1), report = res, RCMdata = RCMdata)
    if (any(do_catch_check)) {
      flag_ind <- paste(grep(TRUE, do_catch_check), collapse = " ")
      if (length(res) > 1) {
        warning("Note: there is predicted catch that deviates from observed catch by more than 1% in simulations: ", flag_ind)
      } else {
        warning("Note: there is predicted catch that deviates from observed catch by more than 1%")
      }
    }
  }
  if (!silent) message("Complete.")
  
  return(output)
}


RCM_report_samps <- function(x, samps, obj, conv) {
  report <- obj$report(samps[x, ]) %>% RCM_posthoc_adjust(obj, par = samps[x, ])
  report$conv <- conv
  return(report)
}

RCM_update_OM <- function(OM, report, StockPars = NULL, obj_data, maxage, nyears, proyears, nsim = length(report), prior = list(), silent = FALSE) {
  # Get parameters from RCM data and report list
  RCM_val <- .RCM_update_OM(report, obj_data, maxage, nyears, proyears, nsim)
  
  ### R0
  if (!silent) message("Updating operating model:\n\n")
  OM@cpars$R0 <- RCM_val$R0
  if (!silent) message("Range of unfished age-0 recruitment (OM@cpars$R0): ", paste(round(range(OM@cpars$R0), 2), collapse = " - "))
  
  ### Depletion and init D - init D is only reported, OM setup for initD by adjusting rec devs
  if (!silent) message("Range of initial spawning depletion: ", paste(round(range(RCM_val$initD), 2), collapse = " - "))
  
  OM@cpars$D <- RCM_val$D
  if (!silent) message("Range of spawning depletion (OM@cpars$D): ", paste(round(range(OM@cpars$D), 2), collapse = " - "), "\n")
  
  ### Selectivity and F
  ### Find
  OM@isRel <- FALSE
  OM@cpars$V <- RCM_val$V
  OM@cpars$Find <- RCM_val$Find
  
  maxFind <- max(OM@cpars$Find)
  
  if (maxFind > OM@maxF) {
    OM@maxF <- maxFind
    if (!silent) message("Updating OM@maxF to maximum F across all fleets: ", round(maxFind, 2))
  }
  
  OM@cpars$qs <- rep(1, nsim)
  if (!silent) {
    message("Historical F set with OM@cpars$Find and OM@cpars$qs.")
    message("Annual selectivity at age set in OM@cpars$V. Projection period uses selectivity of last historical year.")
  }
  
  Eff <- apply(OM@cpars$Find, 2, range)
  OM@EffLower <- Eff[1, ]
  OM@EffUpper <- Eff[2, ]
  if (length(OM@EffYears) != nyears) OM@EffYears <- 1:nyears
  if (length(OM@Esd) == 0 && is.null(OM@cpars$Esd)) OM@Esd <- c(0, 0)
  #if (!silent) message("Historical effort trends set in OM@EffLower and OM@EffUpper (N.B. not used).")
  
  #if (any(obj_data$CAL_hist > 0, na.rm = TRUE) || (any(obj_data$msize > 0, na.rm = TRUE) & obj_data$msize_type == "length") ||
  #    any(obj_data$IAL_hist > 0, na.rm = TRUE)) {
    OM@cpars$CAL_bins <- obj_data$lbin
    OM@cpars$CAL_binsmid <- obj_data$lbinmid
    if (!silent) message("RCMdata length bins will be added to OM.")
    if (!is.null(RCM_val$SLarray)) {
      OM@cpars$SLarray <- RCM_val$SLarray
      if (!silent) message("Annual selectivity at length set in OM@cpars$SLarray (used to simulate length data).")
    }
  #}
  
  ### Rec devs
  OM@cpars$Perr <- RCM_val$procsd
  if (!silent) message("Recruitment standard deviation set in OM@cpars$Perr: ", paste(round(range(OM@cpars$Perr), 2), collapse = " - "))
  
  if (!is.null(StockPars$Perr_y)) {
    Perr_y <- StockPars$Perr_y
  } else {
    Perr_y <- matrix(1, nsim, maxage + nyears + proyears)
  }
  Perr_y[, 1:maxage] <- RCM_val$early_Perr
  Perr_y[, maxage + 1:nyears] <- RCM_val$Perr
  if (!silent) message("Historical recruitment deviations set in OM@cpars$Perr_y.")
  
  # Resample future recruitment with autocorrelation
  OM@cpars$AC <- RCM_val$AC
  OM@AC <- range(RCM_val$AC)
  if (any(RCM_val$AC != 0)) {
    if (!silent) message("Range of recruitment autocorrelation OM@AC: ", paste(round(range(OM@AC), 2), collapse = " - "))
    
    OM@cpars$Perr_y <- RCM_sample_future_dev(obj_data$est_rec_dev, RCM_val$procsd, RCM_val$AC, 
                                             RCM_val$log_rec_dev, Perr_y, maxage, nyears, proyears,
                                             silent = silent)
    if (!silent) message("Future recruitment deviations in OM@cpars$Perr_y sampled with autocorrelation.")
  } else {
    OM@cpars$Perr_y <- Perr_y
  }
  
  ### Mesnil-Rochet stock recruit relationship
  if (OM@SRrel == 3 && !is.null(RCM_val$MRRmax)) {
    
    SRRfun <- function(SB, SRRpars) {
      MesnilRochet_SR(x = SB, Shinge = SRRpars$MRhinge, Rmax = SRRpars$MRRmax, gamma = SRRpars$MRgamma)
    }
    
    SRRpars <- data.frame(MRRmax = RCM_val$MRRmax,
                          MRhinge = RCM_val$MRhinge,
                          MRgamma = RCM_val$MRgamma)
    
    relRfun <- function(SSBpR, SRRpars) {
      MesnilRochet_SR(x = SSBpR, Shinge = SRRpars$MRhinge, Rmax = SRRpars$MRRmax, gamma = SRRpars$MRgamma,
                      isSB = FALSE)
    }
    
    SPRcrashfun <- function(SSBpR0, SRRpars) {
      MesnilRochet_SPRcrash(SSBpR0, Shinge = SRRpars$MRhinge, Rmax = SRRpars$MRRmax, gamma = SRRpars$MRgamma)
    }
    
    OM@cpars$SRR <- list(SRRfun = SRRfun, SRRpars = SRRpars,
                         relRfun = relRfun, SPRcrashfun = SPRcrashfun)
    if (!silent) message("Mesnil-Rochet stock recruitment relationship specified in OM@cpars$SRR.")
  }
  
  ### Assign OM variables that were used in the RCM to the output
  if (!is.null(StockPars$Len_age)) OM@cpars$Len_age <- StockPars$Len_age
  if (!is.null(StockPars$Linf)) OM@cpars$Linf <- StockPars$Linf
  if (!is.null(StockPars$K)) OM@cpars$K <- StockPars$K
  if (!is.null(StockPars$t0)) OM@cpars$t0 <- StockPars$t0
  if (!is.null(StockPars$LenCV)) OM@cpars$LenCV <- StockPars$LenCV
  if (!is.null(StockPars$LatASD)) OM@cpars$LatASD <- StockPars$LatASD
  if (!is.null(StockPars$Wt_age)) OM@cpars$Wt_age <- StockPars$Wt_age
  if (!is.null(StockPars$Mat_age)) OM@cpars$Mat_age <- StockPars$Mat_age
  
  if (!is.null(StockPars$ageMarray)) OM@cpars$ageMarray <- StockPars$ageMarray
  if (!is.null(StockPars$age95array)) OM@cpars$age95array <- StockPars$age95array
  
  if (!is.null(StockPars$L50array)) OM@cpars$L50array <- StockPars$L50array
  if (!is.null(StockPars$L95array)) OM@cpars$L95array <- StockPars$L95array
  
  if (!is.null(StockPars$Fec_Age) && !identical(StockPars$Mat_age * StockPars$Wt_age, StockPars$Fec_Age)) {
    OM@cpars$Fec_age <- StockPars$Fec_Age
  }
  
  if (!OM@SRrel == 3) {
    if (prior$use_prior[2]) {
      OM@cpars$hs <- RCM_val$h
    } else if (!is.null(StockPars$hs)) {
      OM@cpars$hs <- StockPars$hs
    }
  }
  
  if (prior$use_prior[3]) {
    OM@cpars$M_ageArray <- array(RCM_val$Mest, c(nsim, maxage+1, nyears + proyears))
  } else if (!is.null(StockPars$M_ageArray)) {
    OM@cpars$M_ageArray <- StockPars$M_ageArray
  }
  if (!silent) message("Growth, maturity, natural mortality, and stock recruit parameters from RCM are set in OM@cpars.\n\n")
  
  if (!obj_data$plusgroup) {
    OM@cpars$plusgroup <- 0L
    if (!silent) message("No plus group was used in RCM.")
  }
  return(list(OM = OM, RCM_val = RCM_val))
}

.RCM_update_OM <- function(report, obj_data, maxage, nyears, proyears, nsim = length(report)) {
  n_age <- maxage + 1
  out <- list()
  
  out$R0 <- vapply(report, getElement, numeric(1), "R0")
  out$initD <- vapply(report, function(x) x$E[1]/x$E0_SR, numeric(1))
  out$D <- vapply(report, function(x) x$E[length(x$E)-1]/x$E0_SR, numeric(1))
  
  make_F <- function(x) { # Extra step to avoid apical F = 0
    apicalF <- x$F
    apicalF[apicalF < 1e-4] <- 1e-4
    F_at_age <- lapply(1:ncol(apicalF), function(xx) apicalF[, xx] * x$vul[1:nyears, , xx]) %>% 
      simplify2array() %>% apply(1:2, sum)
    return(F_at_age)
  }
  F_matrix <- lapply(report, make_F)
  apical_F <- lapply(F_matrix, function(x) apply(x, 1, max))
  
  expand_V_matrix <- function(x) {
    y <- matrix(x[nyears, ], proyears, ncol(x), byrow = TRUE)
    rbind(x, y)
  }
  
  out$V <- Map("/", e1 = F_matrix, e2 = apical_F) %>% lapply(expand_V_matrix) %>% simplify2array() %>% aperm(3:1)
  out$Find <- do.call(rbind, apical_F)
  
  if (!is.null(report[[1]]$vul_len) && all(!is.na(report[[1]]$vul_len))) {
    make_SL <- function(x) {
      apicalF <- x$F
      apicalF[apicalF < 1e-4] <- 1e-4
      
      F_at_length <- lapply(1:ncol(apicalF), function(xx) {
        sel_block_f <- obj_data$sel_block[1:nyears, xx]
        apicalF[, xx] * t(x$vul_len[, sel_block_f])
      }) %>% 
        simplify2array() %>% 
        apply(1:2, sum)
      SL <- apply(F_at_length, 1, function(xx) xx/max(xx)) %>% t() # year x bin
      return(SL)
    }
    out$SLarray <- lapply(report, make_SL) %>% lapply(expand_V_matrix) %>% simplify2array() %>% aperm(3:1)
  }
  
  out$procsd <- vapply(report, getElement, numeric(1), "tau")
  
  make_Perr <- function(x, obj_data) {
    bias_corr <- ifelse(obj_data$est_rec_dev, exp(-0.5 * x$tau^2), 1)
    res <- exp(x$log_rec_dev) * bias_corr
    res[1] <- res[1] * x$R_eq/x$R0
    return(res)
  }
  out$Perr <- vapply(report, make_Perr, numeric(nyears), obj_data = obj_data) %>% t()
  
  make_early_Perr <- function(x, obj_data) {
    M <- if (is.null(x$Mest)) obj_data$M_data[1, ] else rep(x$Mest, n_age)
    NPR_unfished <- calc_NPR(exp(-M), length(M), obj_data$plusgroup)
    #NPR_unfished <- x$NPR_unfished[1, ]
    res <- x$R_eq * x$NPR_equilibrium / x$R0 / NPR_unfished
    bias_corr <- ifelse(obj_data$est_early_rec_dev, exp(-0.5 * x$tau^2), 1)
    early_dev <- exp(x$log_early_rec_dev) * bias_corr
    out <- res[-1] * early_dev
    return(rev(out))
  }
  out$early_Perr <- vapply(report, make_early_Perr, numeric(maxage), obj_data = obj_data) %>% t()
  
  out$log_rec_dev <- vapply(report, getElement, numeric(nyears), "log_rec_dev") %>% t()

  if (!all(out$log_rec_dev == 0)) {
    out$AC <- apply(out$log_rec_dev, 1, function(x) {
      out <- acf(x, lag.max = 1, plot = FALSE)$acf[2]
      ifelse(is.na(out), 0, out)
    })
  } else {
    out$AC <- rep(0, length(report))
  }
  
  out$h <- vapply(report, getElement, numeric(1), "h")
  out$Mest <- vapply(report, function(x) ifelse(is.null(x$Mest), NA_real_, x$Mest), numeric(1))
  
  out$SSB <- vapply(report, getElement, numeric(nyears + 1), "E") %>% t()
  out$NAA <- sapply(report, getElement, "N", simplify = "array") %>% aperm(c(3, 1, 2))
  out$CAA <- sapply(report, getElement, "CAApred", simplify = "array") %>% aperm(c(4, 1:3))
  if (!is.null(report[[1]]$CALpred)) {
    out$CAL <- sapply(report, getElement, "CALpred", simplify = "array") %>% aperm(c(4, 1:3))
  }
  
  # Mesnil Rochet parameters
  if (!is.null(report[[1]]$MR_SRR)) {
    out$MRRmax <- vapply(report, getElement, numeric(1), "MRRmax") 
    out$MRhinge <- vapply(report, getElement, numeric(1), "MRhinge") 
    out$MRgamma <- vapply(report, getElement, numeric(1), "MRgamma") 
  }
  
  if (length(report) == 1) {
    lapply(out, function(x) {
      if (is.array(x)) {
        dx <- dim(x)
        ldx <- length(dx)
        array(x, dim = c(dx[-1], nsim)) %>% aperm(c(ldx, 2:ldx - 1))
      } else {
        replicate(nsim, x)
      }
    })
  } else {
    out
  }
}


RCM_sample_future_dev <- function(est_rec_dev, procsd, AC, log_rec_dev, Perr_y, maxage, nyears, proyears, silent = FALSE) {
  if (any(!est_rec_dev)) { # Sample historical rec devs for OM for most recent years
    yr_fixed_rec_dev <- which(!est_rec_dev)
    if (any(yr_fixed_rec_dev > sum(est_rec_dev))) {
      yr_hist_sample <- yr_fixed_rec_dev[yr_fixed_rec_dev > sum(est_rec_dev)] %>% min()
      nyears <- length(est_rec_dev)
      samp_hist <- Map(dev_AC, AC = AC, stdev = procsd, 
                       chain_start = log_rec_dev[, yr_hist_sample - 1],
                       MoreArgs = list(n = length(yr_hist_sample:nyears), mu = 1))
      log_rec_dev[, yr_hist_sample:nyears] <- do.call(rbind, samp_hist)
      
      Perr_y[, (maxage + yr_hist_sample):(maxage + nyears)] <- exp(log_rec_dev[, yr_hist_sample:nyears])
      if (!silent) message("Historical recruitment deviations sampled with autocorrelation starting in year ", yr_hist_sample, " out of OM@nyears = ", nyears)
    }
  }
  samp_proj <- Map(dev_AC, AC = AC, stdev = procsd, 
                   chain_start = log_rec_dev[, nyears],
                   MoreArgs = list(n = proyears, mu = 1))
  
  Perr_y[, maxage + nyears + 1:proyears] <- exp(do.call(rbind, samp_proj))
  return(Perr_y)
}

process_AddIndType <- function(s_sel, nfleet) {
  if (s_sel == -3) { # SSB
    return(2)
  } else if (s_sel %in% c(-6, -5, -4, -2, -1, 0)) { # B (see SAMtool::int_s_sel for codes)
    return(1)
  } else if (s_sel > 0) { # Fleet
    return(ifelse(nfleet > 1, 1, 3))
  } else {
    stop("process_AddIndType() cannot process s_sel code = ", s_sel)
  }
}

process_AddIndV <- function(sur, Misc, s_sel, n_age, nfleet, nyears) { # Return a matrix of nsim x nages
  sel_B_SB <- s_sel[sur] %in% c(-3, -4) # -4 = B, -3 = SSB, coordinate with process_AddIndType() and Data@AddIndType
  sel_VB <- s_sel[sur] == 1 && nfleet == 1
  if (sel_B_SB || sel_VB) { # single-fleet VB
    out <- matrix(1, length(Misc), n_age)
  } else { # custom sel or multi-fleet
    out <- do.call(rbind, lapply(Misc, function(x) x$ivul[nyears, , sur]))
  }
  return(out)
}

#' @importFrom pbapply pbsapply
RCM_retro <- function(x, nyr = 5) {
  if (!length(x@mean_fit)) stop("Re-run RCM() with argument `mean_fit = TRUE`", .call = FALSE)
  
  data <- x@mean_fit$obj$env$data
  n_y <- data$n_y
  
  if (data$nfleet > 1) {
    retro_ts <- array(NA, dim = c(nyr+1, n_y + 1, data$nfleet + 4))
    TS_var <- c(paste("Fleet", 1:data$nfleet, "F"), "Apical F", "SSB", "SSB_SSB0", "R")
  } else {
    retro_ts <- array(NA, dim = c(nyr+1, n_y + 1, 4))
    TS_var <- c("Apical F", "SSB", "SSB_SSB0", "R")
  }
  dimnames(retro_ts) <- list(Peel = 0:nyr, Year = (x@OM@CurrentYr - n_y):x@OM@CurrentYr + 1, Var = TS_var)
  
  new_args <- lapply(n_y - 0:nyr, RCM_retro_subset, 
                     data = data, 
                     params = clean_tmb_parameters(x@mean_fit$obj), 
                     map = x@mean_fit$obj$env$map)
  lapply_fn <- function(i, new_args, x) {
    obj2 <- MakeADFun(data = new_args[[i+1]]$data, parameters = new_args[[i+1]]$params, map = new_args[[i+1]]$map,
                      random = x@mean_fit$obj$env$random, DLL = "SAMtool", silent = TRUE)
    if (any(new_args[[i+1]]$data$condition == "1L")) {
      if (any(is.na(obj2$report(obj2$par)$F)) || any(is.infinite(obj2$report(obj2$par)$F))) {
        for(ii in 1:10) {
          obj2$par["R0x"] <- 0.5 + obj2$par["R0x"]
          if (all(!is.na(obj2$report(obj2$par)$F)) && all(!is.infinite(obj2$report(obj2$par)$F))) break
        }
      }
    }
    mod <- optimize_TMB_model(obj2, control = list(iter.max = 2e+05, eval.max = 4e+05))
    opt2 <- mod[[1]]
    SD <- mod[[2]]
    
    if (!is.character(opt2)) {
      report <- obj2$report(obj2$env$last.par.best) %>% RCM_posthoc_adjust(obj2)
      
      FMort <- rbind(report$F, matrix(NA, i + 1, ncol(report$F)))
      if (data$nfleet > 1) FMort <- cbind(FMort, apply(report$F_at_age, 1, max, na.rm = TRUE) %>% c(rep(NA, i+1)))
      
      SSB <- c(report$E, rep(NA, i))
      SSB_SSB0 <- SSB/report$E0_SR
      R <- c(report$R, rep(NA, i))
      
      retro_ts[i+1, , ] <<- cbind(FMort, SSB, SSB_SSB0, R)
      
      return(SD$pdHess)
    }
    return(FALSE)
  }
  
  conv <- pbsapply(0:nyr, lapply_fn, new_args = new_args, x = x, 
                   cl = if (snowfall::sfIsRunning()) snowfall::sfGetCluster() else NULL)
  if (any(!conv)) warning("Peels that did not converge: ", paste0(which(!conv) - 1, collapse = " "))
  
  retro <- new("retro", Model = "RCM", Name = x@OM@Name, TS_var = TS_var, TS = retro_ts)
  if (data$nfleet > 1) {  
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
  params_out <- params
  
  if (yr < data$n_y) {
    ##### Data object
    mat <- c("C_hist", "E_hist", "I_hist", "sigma_I", "sigma_C", "CAA_n", "CAL_n", "IAA_n", "IAL_n", "msize")
    if (!is.null(map$log_M)) mat <- c(mat, "M_data")
    mat_ind <- match(mat, names(data_out))
    data_out[mat_ind] <- lapply(data_out[mat_ind], function(x) x[1:yr, , drop = FALSE])
    
    mat2 <- c("len_age", "wt", "mat", "sel_block", "SD_LAA", "fec")
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
    if (!is.null(map$log_rec_dev)) map$log_rec_dev <- map$log_rec_dev[1:yr] %>% factor()
    
    # Update F
    if (any(data_out$condition == 0L)) { # catch
      params_out$log_F_dev <- params_out$log_F_dev[1:yr, , drop = FALSE]
      
      if (any(data_out$yind_F + 1 > yr)) {
        data_out$yind_F <- as.integer(0.5 * yr)
        params_out$log_F_dev[data_out$yind_F + 1, ] <- params$log_F_dev[data$yind_F + 1, ]
      }
      
      if (!is.matrix(map$log_F_dev)) { # For newer versions of TMB (> 1.9.4 ?)
        map$log_F_dev <- map$log_F_dev %>%
          as.integer() %>%
          matrix(data$n_y, data$nfleet)
      } 
      map$log_F_dev <- map$log_F_dev[1:yr, , drop = FALSE] %>% factor()
    }
    
    ## Update nsel block if needed
    if (data$nsel_block > data_out$nsel_block) {
      sel_block_ind <- data_out$sel_block %>% as.vector() %>% unique()
      params_out$vul_par <- params_out$vul_par[, sel_block_ind, drop = FALSE]
      if (!is.null(map$vul_par)) map$vul_par <- map$vul_par[, sel_block_ind, drop = FALSE] %>% factor()
    }
  }
  
  return(list(data = data_out, params = params_out, map = map))
}

#' @importFrom abind abind
profile_likelihood_RCM <- function(x, ...) {
  dots <- list(...)
  if (!length(x@mean_fit)) stop("No model found. Re-run RCM with mean_fit = TRUE.")
  
  LWT <- try(x@data@Misc$LWT, silent = TRUE)
  if (is.character(LWT)) LWT <- x@data$LWT
  
  new_args <- list(data = x@mean_fit$obj$env$data,
                   params = clean_tmb_parameters(x@mean_fit$obj),
                   map = x@mean_fit$obj$env$map,
                   random = x@mean_fit$obj$env$random)
  
  if (!is.null(dots$D)) {
    if (length(dots) > 1) message("Only doing depletion profile...")
    if (!requireNamespace("reshape2", quietly = TRUE)) {
      stop("Please install the reshape2 package to run this profile.")
    }
    profile_grid <- expand.grid(D = dots$D)
    
    # Create a dummy index = depletion
    n_y <- new_args$data$n_y
    new_args$data$IAA_hist <- abind::abind(new_args$data$IAA_hist, 
                                           array(NA, c(n_y, new_args$data$n_age, 1)), 
                                           along = 3)
    new_args$data$IAL_hist <- abind::abind(new_args$data$IAL_hist, 
                                           array(NA, c(n_y, length(new_args$data$lbinmid), 1)), 
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
    if (!is.null(new_args$map$ivul_par)) {
      new_args$map$ivul_par <- factor(c(new_args$map$ivul_par, rep(NA, nrow(new_args$params$ivul_par))))
    }
    
    LWT$Index <- c(LWT$Index, 1)
    LWT$IAA <- c(LWT$IAA, 1)
    LWT$IAL <- c(LWT$IAL, 1)
    
    profile_fn <- function(i, new_args, x) { # Add new survey of SSB depletion
      SSB_index <- SSB_Isd <- rep(NA_real_, n_y)
      SSB_index[1] <- 1
      SSB_index[n_y] <- profile_grid$D[i]
      SSB_Isd[c(1, n_y)] <- 0.01
      
      new_args$data$I_hist <- cbind(new_args$data$I_hist, SSB_index)
      new_args$data$sigma_I <- cbind(new_args$data$sigma_I, SSB_Isd)
      
      obj2 <- MakeADFun(data = new_args$data, parameters = new_args$params, map = new_args$map,
                        random = x@mean_fit$obj$env$random, DLL = "SAMtool", silent = TRUE)
      mod <- optimize_TMB_model(obj2, control = list(iter.max = 2e+05, eval.max = 4e+05), do_sd = FALSE)
      report <- obj2$report(obj2$env$last.par.best)
      RCM_get_likelihoods(report, LWT = LWT, f_name = paste0("Fleet_", 1:new_args$data$nfleet),
                          s_name = paste0("Index_", 1:(new_args$data$nsurvey-1)) %>% c("SSB_depletion"))
    }
    
  } else {
    
    valid_par <- c("R0", "h", "MRRmax", "MRgamma")
    if (all(!names(dots) %in% valid_par)) stop("No valid parameters found for profiling.")
    
    profile_par <- valid_par[valid_par %in% names(dots)]
    profile_grid <- do.call(expand.grid, dots[profile_par])
    
    profile_fn <- function(i, new_args, x) {
      
      if (new_args$data$SR_type %in% c("BH", "Ricker")) {
        
        if ("R0" %in% profile_par) {
          new_args$params$R0x <- log(profile_grid$R0[i]  * new_args$data$rescale)
          new_args$map$R0x <- factor(NA)
        }
        if ("h" %in% profile_par) {
          new_args$map$transformed_h <- factor(NA)
          if (new_args$data$SR_type == "BH") {
            new_args$params$transformed_h <- logit((profile_grid$h[i] - 0.2)/0.8)
          } else if (new_args$data$SR_type == "Ricker") {
            new_args$params$transformed_h <- log(profile_grid$h[i] - 0.2)
          }
        }
      } else if (new_args$data$SR_type == "Mesnil-Rochet") {
        if ("MRRmax" %in% profile_par) {
          new_args$map$R0x <- factor(NA)
          new_args$params$R0x <- log(profile_grid$MRRmax[i]  * new_args$data$rescale)
        }
        if ("MRgamma" %in% profile_par) {
          new_args$map$MR_SRR[2] <- factor(NA)
          new_args$params$MR_SRR[2] <- log(profile_grid$MRgamma[i])
        }
      }
      
      obj2 <- MakeADFun(data = new_args$data, parameters = new_args$params, map = new_args$map,
                        random = new_args$random, DLL = "SAMtool", silent = TRUE)
      
      mod <- optimize_TMB_model(obj2, control = list(iter.max = 2e+05, eval.max = 4e+05), do_sd = FALSE)
      report <- obj2$report(obj2$env$last.par.best)
      RCM_get_likelihoods(report, LWT = LWT, f_name = paste0("Fleet_", 1:new_args$data$nfleet),
                          s_name = paste0("Index_", 1:new_args$data$nsurvey))
    }
    
  }
  
  do_profile <- pblapply(1:nrow(profile_grid), profile_fn, new_args = new_args, x = x,
                         cl = if (snowfall::sfIsRunning()) snowfall::sfGetCluster() else NULL)
  
  profile_grid$nll <- vapply(do_profile, function(xx) xx[[1]][1, 1], numeric(1))
  
  prof <- sapply(do_profile, RCM_nll_profile_summary) %>% t()
  prof_out <- apply(prof, 2, sum) %>% as.logical()
  
  profile_grid <- cbind(profile_grid, prof[, prof_out] %>% as.data.frame())
  
  if (!is.null(dots$D)) {
    MLE <- x@mean_fit$report$E[n_y]/x@mean_fit$report$E0_SR
  } else {
    MLE <- vapply(profile_par, function(xx) getElement(x@mean_fit$report, xx), numeric(1))
  }
  
  output <- new("prof", Model = "RCM", Par = profile_par, MLE = MLE, grid = profile_grid)
  return(output)
}

#' @importFrom dplyr select starts_with
RCM_nll_profile_summary <- function(x) { # Summary table of RCM likelihoods for each data-gear combination
  
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("Please install the reshape2 package.", call. = FALSE)
  }
  
  nll_fleet <- x[[2]][-nrow(x[[2]]), ] %>% select(starts_with("Fleet"))
  nll_fleet$Data <- rownames(nll_fleet)
  nll_fleet <- reshape2::melt(nll_fleet, id.var = "Data", value.var = "nll")
  
  nll_index <- x[[4]][-nrow(x[[4]]), ] %>% select(starts_with("Index"))
  nll_index$Data <- rownames(nll_index)
  nll_index <- reshape2::melt(nll_index, id.var = "Data", value.var = "nll")
  
  nll <- rbind(nll_fleet, nll_index)
  nll$Name <- paste(nll$variable, nll$Data)
  
  c(nll$value, x[[1]][c(2, 5, 6), 1]) %>% 
    structure(names = paste(nll$variable, nll$Data) %>% c("Recruitment Deviations", "Penalty", "Prior"))
}
