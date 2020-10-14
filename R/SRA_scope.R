

SRA_scope_int <- function(OM, data = list(), condition = c("catch", "catch2", "effort"), selectivity = "logistic", s_selectivity = NULL, LWT = list(),
                          comp_like = c("multinomial", "lognormal"), ESS = c(30, 30),
                          max_F = 3, cores = 1L, integrate = FALSE, mean_fit = FALSE, drop_nonconv = FALSE,
                          drop_highF = FALSE,
                          control = list(iter.max = 2e+05, eval.max = 4e+05), ...) {

  dots <- list(...) # can be vul_par, s_vul_par, map_vul_par, map_s_vul_par, map_log_rec_dev, map_early_rec_dev, rescale, plusgroup, resample, OMeff, fix_dome
  if(!is.null(dots$maxF)) max_F <- dots$maxF
  if(length(ESS) == 1) ESS <- rep(ESS, 2)

  comp_like <- match.arg(comp_like)
  condition <- match.arg(condition)

  dat_update <- update_SRA_data(data, OM, condition, dots)
  OM <- dat_update$OM
  data <- dat_update$data
  StockPars <- dat_update$StockPars
  FleetPars <- dat_update$FleetPars
  ObsPars <- dat_update$ObsPars

  nsim <- OM@nsim
  proyears <- OM@proyears
  maxage <- OM@maxage
  nyears <- data$nyears
  nfleet <- data$nfleet
  nsurvey <- data$nsurvey

  OM@maxF <- max_F
  message("OM@maxF updated to ", max_F, ".")

  # No comp data
  if(!any(data$CAA > 0) && !any(data$CAL > 0)) {
    fix_sel <- TRUE
    message("No fishery length or age compositions were provided. Selectivity is fixed to values from OM.")
  } else {
    fix_sel <- FALSE
  }

  # Selectivity
  if(length(selectivity) == 1) selectivity <- rep(selectivity, data$nsel_block)
  if(length(selectivity) < data$nsel_block) stop("selectivity vector should be of length ", data$nsel_block, ").", call. = FALSE)
  sel <- int_sel(selectivity)

  # Survey selectivity
  if(!is.null(data$I_type)) {
    message("\n\n*** WARNING: data$I_type is now obsolete. Use s_selectivity argument instead. *** \n\n")
    if(all(data$I_type != "est") && is.null(s_selectivity)) {
      message("*** For this model, data$I_type should still be compatible and moved over to s_selectivity. Attempting to run model...\n\n")
      s_selectivity <- data$I_type
      data$I_type <- NULL
    } else if(any(data$I_type == "est") && length(s_selectivity) == nsurvey) {
      s_selectivity[data$I_type != "est"] <- data$I_type[data$I_type != "est"]
      message("*** Attempting to update s_selectivity from data$I_type. Argument s_selectivity updated to:")
      message("c(", paste(s_selectivity, collapse = ", "), ")")
      message("\n\n")
    }
  }
  if(nsurvey > 0) {
    if(is.null(s_selectivity)) s_selectivity <- rep("B", nsurvey)
    if(length(s_selectivity) == 1) s_selectivity <- rep(s_selectivity, nsurvey)
    s_sel <- int_s_sel(s_selectivity, nfleet)
  } else {
    s_sel <- int_s_sel("B")
  }

  # Likelihood weights
  data$LWT <- make_LWT(LWT, nfleet, nsurvey)

  # SR
  message(ifelse(OM@SRrel == 1, "Beverton-Holt", "Ricker"), " stock-recruitment relationship used.")

  # Test for identical sims
  par_identical_sims <- par_identical_sims_fn(StockPars, FleetPars, ObsPars, data, dots)

  # Fit model
  if(!is.null(dots$resample) && dots$resample) { # Re-sample covariance matrix

    message("\nResample = TRUE. Running mean fit model first...")
    mean_fit_output <- SRA_scope_est(data = data, selectivity = sel, s_selectivity = s_sel,
                                     SR_type = ifelse(OM@SRrel == 1, "BH", "Ricker"), LWT = data$LWT, comp_like = comp_like, ESS = ESS,
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
      report <- obj$report(samps[x, ]) %>% SRA_posthoc_adjust(obj, par = samps[x, ])
      report$conv <- conv
      return(report)
    }

    res <- lapply(1:nsim, report_internal_fn, samps = samps, obj = mean_fit_output$obj, conv = mean_fit_output$report$conv)
    mod <- lapply(res, function(x) list(obj = mean_fit_output$obj, report = x))
    conv <- rep(mean_fit_output$report$conv, nsim)
  } else {

    if(all(par_identical_sims)) { # All identical sims detected
      message("\nAll ", nsim, " replicates are identical. Fitting once and replicating single fit...")

      mean_fit_output <- SRA_scope_est(data = data, selectivity = sel, s_selectivity = s_sel,
                                       SR_type = ifelse(OM@SRrel == 1, "BH", "Ricker"), LWT = data$LWT, comp_like = comp_like, ESS = ESS,
                                       max_F = max_F, integrate = integrate, StockPars = StockPars, ObsPars = ObsPars,
                                       FleetPars = FleetPars, mean_fit = TRUE, control = control, dots = dots)

      mod <- lapply(1:nsim, function(x) return(mean_fit_output))
    } else {

      message("\nFitting model (", nsim, " simulations) ...")
      if(cores > 1 && !snowfall::sfIsRunning()) DLMtool::setup(as.integer(cores))
      if(snowfall::sfIsRunning()) {
        mod <- snowfall::sfClusterApplyLB(1:nsim, SRA_scope_est, data = data, selectivity = sel, s_selectivity = s_sel,
                                          SR_type = ifelse(OM@SRrel == 1, "BH", "Ricker"), LWT = data$LWT, comp_like = comp_like, ESS = ESS,
                                          max_F = max_F, integrate = integrate, StockPars = StockPars, ObsPars = ObsPars,
                                          FleetPars = FleetPars, control = control, dots = dots)
      } else {
        mod <- lapply(1:nsim, SRA_scope_est, data = data, selectivity = sel, s_selectivity = s_sel,
                      SR_type = ifelse(OM@SRrel == 1, "BH", "Ricker"), LWT = data$LWT, comp_like = comp_like, ESS = ESS,
                      max_F = max_F, integrate = integrate, StockPars = StockPars, ObsPars = ObsPars,
                      FleetPars = FleetPars, control = control, dots = dots)
      }

      if(mean_fit) { ### Fit to life history means if mean_fit = TRUE
        message("Generating additional model fit from mean values of parameters in the operating model...\n")
        mean_fit_output <- SRA_scope_est(data = data, selectivity = sel, s_selectivity = s_sel,
                                         SR_type = ifelse(OM@SRrel == 1, "BH", "Ricker"), LWT = data$LWT, comp_like = comp_like, ESS = ESS,
                                         max_F = max_F, integrate = integrate, StockPars = StockPars, ObsPars = ObsPars,
                                         FleetPars = FleetPars, mean_fit = TRUE, control = control, dots = dots)

        if(!mean_fit_output$report$conv) warning("Mean fit model did not appear to converge.")
      } else mean_fit_output <- list()
    }
    res <- lapply(mod, getElement, "report")
    conv <- vapply(res, getElement, logical(1), name = "conv")
  }

  if(drop_highF) {
    highF <- vapply(res, function(x) getElement(x, "F") %>% max(na.rm = TRUE) >= max_F, logical(1))
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
  message("Range of unfished recruitment (OM@cpars$R0): ", paste(round(range(OM@cpars$R0), 2), collapse = " - "))

  ### Depletion and init D
  if(any(data$C_eq > 0) || any(data$E_eq > 0)) {
    initD <- vapply(res, function(x) x$E[1]/x$E0_SR, numeric(1))
    message("Estimated range in initial spawning depletion: ", paste(round(range(initD), 2), collapse = " - "))
  }

  OM@cpars$D <- vapply(res, function(x) x$E[length(x$E)-1]/x$E0_SR, numeric(1))
  message("Range of spawning depletion (OM@cpars$D): ", paste(round(range(OM@cpars$D), 2), collapse = " - "), "\n")

  ### Selectivity and F
  ### Find
  OM@isRel <- FALSE
  if(data$nsel_block > 1 || sel == -2) {
    F_matrix <- lapply(res, getElement, "F_at_age")
    apical_F <- lapply(F_matrix, function(x) apply(x, 1, max))

    expand_V_matrix <- function(x) {
      y <- matrix(x[nyears, ], proyears, maxage, byrow = TRUE)
      rbind(x, y)
    }
    V <- Map("/", e1 = F_matrix, e2 = apical_F) %>% lapply(expand_V_matrix)
    OM@cpars$V <- unlist(V) %>% array(c(nyears + proyears, maxage, nsim)) %>% aperm(c(3, 2, 1))
    OM@cpars$Find <- do.call(rbind, apical_F)
    message("Historical F and selectivity trends set in OM@cpars$Find and OM@cpars$V, respectively.")
    message("Selectivity during projection period is set to that in most recent historical year.")

  } else { # nsel_block = 1

    OM@cpars$L5 <- vapply(res, getElement, numeric(1), "L5")
    message("Range of OM@cpars$L5: ", paste(round(range(OM@cpars$L5), 2), collapse = " - "))

    OM@cpars$LFS <- vapply(res, getElement, numeric(1), "LFS")
    message("Range of OM@cpars$LFS: ", paste(round(range(OM@cpars$LFS), 2), collapse = " - "))

    if(selectivity == "logistic") {
      OM@cpars$Vmaxlen <- rep(1, nsim)
      message("With logistic selectivity, setting OM@cpars$Vmaxlen = 1")

    } else {
      OM@cpars$Vmaxlen <- vapply(res, getElement, numeric(1), "Vmaxlen")
      message("Range of OM@cpars$Vmaxlen: ", paste(round(range(OM@cpars$Vmaxlen), 2), collapse = " - "))
    }

    OM@cpars$Find <- t(do.call(cbind, lapply(res, getElement, "F")))
    message("Historical F set in OM@cpars$Find.")
  }

  if(packageVersion("DLMtool") >= "5.4.4") OM@cpars$qs <- rep(1, nsim)
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
  OM@cpars$Perr_y[, 1:(OM@maxage - 1)] <- early_Perr
  OM@cpars$Perr_y[, OM@maxage:(OM@maxage + nyears - 1)] <- Perr
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
      proc_mu <- -0.5 * StockPars$procsd^2 * (1 - OM@cpars$AC)/sqrt(1 - OM@cpars$AC^2) # http://dx.doi.org/10.1139/cjfas-2016-0167

      if(!is.null(dots$map_log_rec_dev) && any(is.na(dots$map_log_rec_dev))) { # Sample historical rec devs for OM for most recent years
        yr_fixed_rec_dev <- which(is.na(dots$map_log_rec_dev))
        if(any(yr_fixed_rec_dev > max(dots$map_log_rec_dev, na.rm = TRUE))) {
          yr_hist_sample <- yr_fixed_rec_dev[yr_fixed_rec_dev > max(dots$map_log_rec_dev, na.rm = TRUE)] %>% min()

          log_rec_dev_u <- rnorm(length(yr_hist_sample:nyears) * nsim, proc_mu, StockPars$procsd) %>%
            matrix(nsim, length(yr_hist_sample:nyears))

          for(y in 1:ncol(log_rec_dev_u)) {
            if(y == 1) {
              log_rec_dev_u[, y] <- OM@cpars$AC * log_rec_dev[, yr_hist_sample - 1] + log_rec_dev_u[, y] * sqrt(1 - OM@cpars$AC^2)
            } else {
              log_rec_dev_u[, y] <- OM@cpars$AC * log_rec_dev_u[, y-1] + log_rec_dev_u[, y] * sqrt(1 - OM@cpars$AC^2)
            }
          }
          log_rec_dev[, yr_hist_sample:nyears] <<- log_rec_dev_u
          OM@cpars$Perr_y[, (OM@maxage + yr_hist_sample - 1):(OM@maxage + nyears - 1)] <<- exp(log_rec_dev_u)

          message("Historical recruitment deviations sampled with autocorrelation starting in year ", yr_hist_sample, " out of OM@nyears = ", nyears)
        }
      }

      pro_Perr_y <- rnorm(proyears * nsim, proc_mu, StockPars$procsd) %>% matrix(nsim, proyears)
      for(y in 2:proyears) {
        if(y == 1) {
          pro_Perr_y[, y] <- OM@cpars$AC * log_rec_dev[, ncol(log_rec_dev)] + pro_Perr_y[, y] * sqrt(1 - OM@cpars$AC^2)
        } else {
          pro_Perr_y[, y] <- OM@cpars$AC * pro_Perr_y[, y-1] + pro_Perr_y[, y] * sqrt(1 - OM@cpars$AC^2)
        }
      }
      return(exp(pro_Perr_y))
    }
    OM@cpars$Perr_y[, (OM@maxage+nyears):ncol(OM@cpars$Perr_y)] <- sample_future_dev()
    message("Future recruitment deviations sampled with autocorrelation (in OM@cpars$Perr_y).\n")
  }

  ### Assign OM variables that were used in the SRA to the output
  OM@cpars$Len_age <- StockPars$Len_age
  OM@cpars$Linf <- StockPars$Linf
  OM@cpars$K <- StockPars$K
  OM@cpars$t0 <- StockPars$t0
  OM@cpars$LenCV <- StockPars$LenCV
  OM@cpars$Wt_age <- StockPars$Wt_age

  if(any(apply(StockPars$Mat_age, 1, function(x) all(x >= 0.5)))) { # Any simulations where all mat_age > 0.5?
    OM@cpars$L50 <- StockPars$L50
    OM@cpars$L95 <- StockPars$L95
  } else {
    OM@cpars$Mat_age <- StockPars$Mat_age
  }
  OM@cpars$M_ageArray <- StockPars$M_ageArray

  OM@cpars$h <- StockPars$hs

  if(any(data$CAL > 0, na.rm = TRUE) || (any(data$MS > 0, na.rm = TRUE) & data$MS_type == "length") ||
     any(data$s_CAL > 0, na.rm = TRUE)) {
    bw <- data$length_bin[2] - data$length_bin[1]
    OM@cpars$binWidth <- bw
    OM@cpars$CAL_binsmid <- data$length_bin
    OM@cpars$CAL_bins <- c(data$length_bin - 0.5 * bw, max(data$length_bin) + 0.5 * bw)
    message("SRA length bins will be added to OM.")
  }

  if(is.null(dots$plusgroup) || dots$plusgroup) OM@cpars$plusgroup <- rep(1L, nsim)
  if(is.null(data$I_sd) || !any(data$I_sd > 0, na.rm = TRUE)) OM@cpars$Iobs <- ObsPars$Iobs
  message("Growth, maturity, natural mortality, and steepness values from SRA are set in OM@cpars.\n")

  ### Output list
  E <- do.call(rbind, lapply(res[keep], getElement, "E"))
  N <- array(sapply(res[keep], getElement, "N"), c(nyears+1, maxage, sum(keep)))
  CAA_pred <- array(sapply(res[keep], getElement, "CAApred"), c(nyears, maxage, nfleet, sum(keep)))
  CAL_pred <- array(sapply(res[keep], getElement, "CALpred"), c(nyears, length(data$length_bin), nfleet, sum(keep)))

  output <- new("SRA", OM = Sub_cpars(OM, keep), SSB = E, NAA = aperm(N, c(3, 1, 2)), CAA = aperm(CAA_pred, c(4, 1:3)),
                CAL = aperm(CAL_pred, c(4, 1:3)), mean_fit = mean_fit_output, conv = conv[keep], data = data, Misc = res[keep])

  # Data in cpars
  if(sum(output@data$Chist > 0, na.rm = TRUE) || nsurvey > 0) {

    real_Data <- new("Data")
    real_Data@Year <- (output@OM@CurrentYr - output@OM@nyears + 1):output@OM@CurrentYr
    if(sum(output@data$Chist > 0, na.rm = TRUE) && all(!is.na(output@data$Chist))) {
      real_Data@Cat <- matrix(rowSums(output@data$Chist, na.rm = TRUE), 1, nyears)
      real_Data@CV_Cat <- matrix(sqrt(exp(0.01^1 - 1)), 1, nyears)
      message("Historical catch data added to OM@cpars$Data@Cat with default catch CV = 0.01.")
    }
    if(.hasSlot(real_Data, "AddInd") && nsurvey > 0) {
      real_Data@AddInd <- array(output@data$Index, c(nyears, nsurvey, output@OM@nsim)) %>%
        aperm(perm = c(3, 2, 1))
      real_Data@CV_AddInd <- array(sqrt(exp(output@data$I_sd^2) - 1), c(nyears, nsurvey, output@OM@nsim)) %>%
        aperm(perm = c(3, 2, 1))

      if(.hasSlot(real_Data, "AddIndType") && .hasSlot(real_Data, "AddIunits")) {

        # Cannot accommodate indices mirrored to fleet when nfleet > 1 and fleet has time-varying sel
        real_Data@AddIndType <- vapply(s_sel, process_AddIndType, numeric(1), nfleet = nfleet)
        real_Data@AddIndV <- lapply(1:nsurvey, process_AddIndV, Misc = output@Misc, s_sel = s_sel,
                                    maxage = maxage, nfleet = nfleet, nyears = nyears) %>%
          unlist() %>% array(c(OM@nsim, maxage, nsurvey)) %>% aperm(c(1, 3, 2))
        real_Data@AddIunits <- data$I_units

      } else {  # Backwards compatibility with DLMtool 5.4.4

        real_Data@AddIndV <- lapply(output@Misc, function(x) x$s_vul[nyears, , , drop = FALSE]) %>% unlist() %>%
          array(dim = c(maxage, nsurvey, OM@nsim)) %>% aperm(c(3, 2, 1))
        output@OM@cpars$AddIunits <- data$I_units

      }
      message("Historical indices added to OM@cpars$Data@AddInd.")
    }
    output@OM@cpars$Data <- real_Data
  }

  # Check whether observed matches predicted
  if(data$condition == "catch" || data$condition == "catch2") {
    catch_check_fn <- function(x, report, data) {
      if(report[[x]]$conv) {
        catch_diff <- report[[x]]$Cpred/data$Chist - 1
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


SRA_scope_est <- function(x = 1, data, selectivity, s_selectivity, SR_type = c("BH", "Ricker"), LWT = list(),
                          comp_like = c("multinomial", "lognormal"), ESS = c(30, 30),
                          max_F = 3, integrate = FALSE, StockPars, ObsPars, FleetPars, mean_fit = FALSE,
                          control = list(iter.max = 2e+05, eval.max = 4e+05), inner.control = list(maxit = 1e3), dots = list()) {

  SR_type <- match.arg(SR_type)
  comp_like <- match.arg(comp_like)

  nyears <- data$nyears
  nfleet <- data$nfleet
  max_age <- dim(data$CAA)[2]
  nsurvey <- ncol(data$Index)

  data$CAA <- apply(data$CAA, c(1, 3), SRA_tiny_comp) %>% aperm(c(2, 1, 3))
  data$CAL <- apply(data$CAL, c(1, 3), SRA_tiny_comp) %>% aperm(c(2, 1, 3))
  CAA_n <- apply(data$CAA, c(1, 3), sum, na.rm = TRUE)
  CAL_n <- apply(data$CAL, c(1, 3), sum, na.rm = TRUE)

  data$s_CAA <- apply(data$s_CAA, c(1, 3), SRA_tiny_comp) %>% aperm(c(2, 1, 3))
  data$s_CAL <- apply(data$s_CAL, c(1, 3), SRA_tiny_comp) %>% aperm(c(2, 1, 3))
  s_CAA_n <- apply(data$s_CAA, c(1, 3), sum, na.rm = TRUE)
  s_CAL_n <- apply(data$s_CAL, c(1, 3), sum, na.rm = TRUE)

  if(comp_like == "multinomial") {
    for(ff in 1:nfleet) { # Annual sums to effective sample size
      data$CAA[,,ff] <- data$CAA[,,ff]/CAA_n[,ff] * pmin(CAA_n[,ff], ESS[1])
      data$CAL[,,ff] <- data$CAL[,,ff]/CAL_n[,ff] * pmin(CAL_n[,ff], ESS[2])
    }
    CAA_n <- pmin(CAA_n, ESS[1])
    CAL_n <- pmin(CAL_n, ESS[2])

    for(sur in 1:nsurvey) { # Annual sums to effective sample size
      data$s_CAA[,,sur] <- data$s_CAA[,,sur]/s_CAA_n[,sur] * pmin(s_CAA_n[,sur], ESS[1])
      data$s_CAL[,,sur] <- data$s_CAL[,,sur]/s_CAL_n[,sur] * pmin(s_CAL_n[,sur], ESS[2])
    }
    s_CAA_n <- pmin(s_CAA_n, ESS[1])
    s_CAL_n <- pmin(s_CAL_n, ESS[2])
  }

  LWT_C <- matrix(c(LWT$Chist, LWT$CAA, LWT$CAL, LWT$MS, LWT$C_eq), nrow = nfleet, ncol = 5)
  LWT_Index <- cbind(LWT$Index, LWT$s_CAA, LWT$s_CAL)

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

    StockPars_ind <- match("ageM", names(StockPars))
    StockPars[StockPars_ind] <- lapply(StockPars[StockPars_ind], mean_matrix)

    if(data$condition == "effort") {
      StockPars$R0 <- mean_vector(StockPars$R0)
      if(!is.null(dots$OMeff) && dots$OMeff) {
        FleetPars$Find <- apply(FleetPars$Find, 2, mean) %>% matrix(length(StockPars$R0), data$nyears, byrow = TRUE)
      }
    }

    FleetPars_ind <- match(c("L5", "LFS", "Vmaxlen"), names(FleetPars))
    FleetPars[FleetPars_ind] <- lapply(FleetPars[FleetPars_ind], mean_matrix)

    ObsPars$Isd <- mean_vector(ObsPars$Isd)
  }

  if(is.null(data$I_sd) || !any(data$I_sd > 0, na.rm = TRUE)) data$I_sd <- matrix(sdconv(1, ObsPars$Isd[x]), nyears, nsurvey)

  if(!is.null(data$Chist) && any(data$Chist > 0, na.rm = TRUE)) {
    if(is.null(dots$rescale)) {
      rescale <- 1/mean(data$Chist, na.rm = TRUE)
    } else {
      rescale <- dots$rescale
    }
    C_hist <- data$Chist
  } else {
    rescale <- 1
    C_hist <- matrix(0, nyears, nfleet)
  }

  if(data$condition == "effort" && !is.null(dots$OMeff) && dots$OMeff) {
    data$Ehist <- matrix(FleetPars$Find[x, ], nyears, nfleet)
  }

  if(!is.null(data$Ehist) && any(data$Ehist > 0, na.rm = TRUE)) {
    rescale_effort <- 1/mean(data$Ehist, na.rm = TRUE)
    E_hist <- data$Ehist * rescale_effort
  } else {
    rescale_effort <- 1
    E_hist <- matrix(0, nyears, nfleet)
  }

  if(is.null(dots$nit_F)) nit_F <- 3L else nit_F <- dots$nit_F
  if(is.null(dots$plusgroup)) plusgroup <- 1L else plusgroup <- as.integer(dots$plusgroup)
  age_only_model <- StockPars$Len_age[x, , 1:(nyears+1)] %>% apply(2, function(xx) length(xx) == max_age & max(xx) == max_age) %>% all()
  TMB_data_all <- list(condition = data$condition,
                       nll_C = as.integer((data$condition == "effort" & nfleet > 1) || data$condition == "catch"),
                       I_hist = data$Index, sigma_I = data$I_sd, CAA_hist = data$CAA, CAA_n = pmin(CAA_n, ESS[1]),
                       CAL_hist = data$CAL, CAL_n = pmin(CAL_n, ESS[2]), s_CAA_hist = data$s_CAA, s_CAA_n = s_CAA_n,
                       s_CAL_hist = data$s_CAL, s_CAL_n = s_CAL_n, length_bin = data$length_bin, msize = data$MS, msize_type = data$MS_type,
                       sel_block = rbind(data$sel_block, data$sel_block[nyears, ]), nsel_block = data$nsel_block,
                       n_y = nyears, max_age = max_age, nfleet = nfleet, nsurvey = nsurvey,
                       M = t(StockPars$M_ageArray[x, , 1:nyears]), len_age = t(StockPars$Len_age[x, , 1:(nyears+1)]),
                       Linf = ifelse(age_only_model, max_age, StockPars$Linf[x]),
                       CV_LAA = StockPars$LenCV[x], wt = t(StockPars$Wt_age[x, , 1:(nyears+1)]),
                       mat = t(StockPars$Mat_age[x, , 1:(nyears+1)]), vul_type = as.integer(selectivity),
                       s_vul_type = as.integer(s_selectivity), abs_I = data$abs_I,
                       I_units = as.integer(data$I_units), age_error = data$age_error,
                       SR_type = SR_type, LWT_C = LWT_C, LWT_Index = LWT_Index, comp_like = comp_like,
                       max_F = max_F, rescale = rescale, ageM = min(nyears, ceiling(StockPars$ageM[x, 1])),
                       yind_F = as.integer(rep(0.5 * nyears, nfleet)), nit_F = nit_F, plusgroup = plusgroup)

  if(data$condition == "catch" || data$condition == "catch2") {
    TMB_data <- list(model = "SRA_scope", C_hist = C_hist, C_eq = data$C_eq, E_hist = E_hist, E_eq = rep(0, nfleet))
  } else {
    TMB_data <- list(model = "SRA_scope", C_hist = C_hist, C_eq = rep(0, nfleet), E_hist = E_hist, E_eq = data$E_eq * rescale_effort)
  }

  if(SR_type == "BH") {
    transformed_h <- logit((StockPars$hs[x] - 0.2)/0.8)
  } else transformed_h <- log(StockPars$hs[x] - 0.2)

  LFS <- FleetPars$LFS[nyears, x]
  L5 <- FleetPars$L5[nyears, x]
  Vmaxlen <- FleetPars$Vmaxlen[nyears, x]

  if(is.null(dots$vul_par)) {
    if(any(selectivity == -2)) stop("Some (fleet) selectivity specified to be free parameters. Provide vul_par matrix to SRA_scope.")
    vul_par <- matrix(c(LFS, L5, Vmaxlen), 3, data$nsel_block)
  } else {
    vul_par <- dots$vul_par
  }

  sel_check <- selectivity == -1 | selectivity == 0
  vul_par[2, sel_check] <- log(vul_par[1, sel_check] - vul_par[2, sel_check])
  vul_par[1, sel_check] <- logit(pmin(vul_par[1, sel_check]/TMB_data_all$Linf/0.99, 0.99))
  vul_par[3, sel_check] <- logit(pmin(vul_par[3, sel_check], 0.99))
  if(any(selectivity == -2)) vul_par[, selectivity == -2] <- logit(vul_par[, selectivity == -2], soft_bounds = TRUE)

  if(is.null(dots$map_vul_par)) {
    if(any(selectivity == -2)) stop("Some (fleet) selectivity specified to be free parameters. Provide map_vul_par matrix to SRA_scope.")
    map_vul_par <- matrix(0, 3, data$nsel_block)
    map_vul_par[3, selectivity == -1] <- NA # Fix third parameter for logistic sel
    if(!is.null(dots$fix_dome) && dots$fix_dome) map_vul_par[3, selectivity == 0] <- NA # Fix dome

    for(ff in 1:nfleet) {
      if(all(data$CAA[,,ff] <= 0, na.rm = TRUE) && all(data$CAL[,,ff] <= 0, na.rm = TRUE)) {
        map_vul_par[, unique(data$sel_block[, ff])] <- NA # Fix sel if no comp data
      }
    }
    if(any(!is.na(map_vul_par))) map_vul_par[!is.na(map_vul_par)] <- 1:sum(!is.na(map_vul_par))
  } else {
    map_vul_par <- dots$map_vul_par
  }

  # s_vul_par, and map
  if(is.null(dots$s_vul_par)) {
    if(any(s_selectivity == -2)) stop("Some s_selectivity specified to be free parameters. Provide s_vul_par matrix to SRA_scope.")
    s_vul_par <- matrix(c(LFS, L5, Vmaxlen), 3, nsurvey)
  } else {
    s_vul_par <- dots$s_vul_par
  }

  parametric_sel <- s_selectivity == -1 | s_selectivity == 0
  s_vul_par[2, parametric_sel] <- log(s_vul_par[1, parametric_sel] - s_vul_par[2, parametric_sel])
  s_vul_par[1, parametric_sel] <- logit(s_vul_par[1, parametric_sel]/TMB_data_all$Linf/0.99)
  s_vul_par[3, parametric_sel] <- logit(s_vul_par[3, parametric_sel])
  if(any(s_selectivity == -2)) s_vul_par[, s_selectivity == -2] <- logit(s_vul_par[, s_selectivity == -2], soft_bounds = TRUE)

  if(is.null(dots$map_s_vul_par)) {
    if(any(s_selectivity == -2)) stop("Some s_selectivity specified to be free parameters. Provide map_s_vul_par matrix to SRA_scope.")
    map_s_vul_par <- matrix(0, 3, nsurvey)
    map_s_vul_par[3, s_selectivity < 0] <- NA # if logistic
    for(sur in 1:nsurvey) {
      if(s_selectivity[sur] < -2 || s_selectivity[sur] > 0 ||
         (all(data$s_CAA[,,sur] <= 0, na.rm = TRUE) & all(data$s_CAL[,,sur] <= 0, na.rm = TRUE))) {
        map_s_vul_par[, sur] <- NA
      }
    }
    if(any(!is.na(map_s_vul_par))) map_s_vul_par[!is.na(map_s_vul_par)] <- 1:sum(!is.na(map_s_vul_par))
  } else {
    map_s_vul_par <- dots$map_s_vul_par
  }

  log_F_start <- matrix(0, nyears, nfleet)
  if(data$condition == "catch") log_F_start[TMB_data_all$yind_F + 1, 1:nfleet] <- log(0.5 * mean(TMB_data_all$M[nyears, ]))

  TMB_params <- list(R0x = ifelse(TMB_data_all$nll_C | data$condition == "catch2", log(StockPars$R0[x] * rescale), 0),
                     transformed_h = transformed_h, vul_par = vul_par, s_vul_par = s_vul_par,
                     log_q_effort = rep(log(0.1), nfleet), log_F = log_F_start,
                     log_F_equilibrium = rep(log(0.05), nfleet),
                     log_CV_msize = log(data$MS_cv), log_tau = log(StockPars$procsd[x]),
                     log_early_rec_dev = rep(0, max_age - 1), log_rec_dev = rep(0, nyears))

  map <- list()
  if(data$condition == "effort" && !TMB_data_all$nll_C) map$R0x <- factor(NA)
  map$transformed_h <- map$log_tau <- factor(NA)
  map$vul_par <- factor(map_vul_par)
  map$s_vul_par <- factor(map_s_vul_par)
  if(data$condition != "effort") {
    map$log_q_effort <- factor(rep(NA, nfleet))
    if(any(data$C_eq == 0)) {
      map_log_F_equilibrium <- rep(NA, nfleet)
      map_log_F_equilibrium[data$C_eq > 0] <- 1:sum(data$C_eq > 0)
      map$log_F_equilibrium <- factor(map_log_F_equilibrium)
    }
  } else {
    map$log_F_equilibrium <- factor(rep(NA, nfleet))
  }
  if(data$condition != "catch") map$log_F <- factor(matrix(NA, nyears, nfleet))
  map$log_CV_msize <- factor(rep(NA, nfleet))

  if(is.null(dots$map_log_early_rec_dev)) {
    map$log_early_rec_dev <- factor(rep(NA, max_age - 1))
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

  obj <- MakeADFun(data = c(TMB_data, TMB_data_all), parameters = TMB_params, map = map, random = random,
                   inner.control = inner.control, DLL = "MSEtool", silent = TRUE)

  if(data$condition == "catch2") {
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
  report <- obj$report(obj$env$last.par.best) %>% SRA_posthoc_adjust(obj)

  if(data$condition == "effort" && any(data$Chist > 0, na.rm = TRUE)) {
    vars_div <- c("B", "E", "C_eq_pred", "CAApred", "CALpred", "s_CAApred", "s_CALpred", "CN", "Cpred", "N", "VB",
                  "R", "R_early", "R_eq", "R0", "B0", "E0", "N0", "E0_SR")
    vars_mult <- c("Brec", "q")
    var_trans <- c("R0", "q")
    fun_trans <- c("/", "*")
    rescale <- 1/exp(mean(log(data$Chist/report$Cpred), na.rm = TRUE))
    fun_fixed <- c(NA, NA)
    rescale_report(vars_div, vars_mult, var_trans, fun_trans, fun_fixed)
  }

  return(list(obj = obj, opt = opt, SD = SD, report = c(report, list(conv = !is.character(opt) && SD$pdHess))))
}



#' @rdname SRA_scope
#' @export
Sub_cpars <- function(OM, sims = 1:OM@nsim) {

  if(is.numeric(sims)) {
    sims2 <- logical(OM@nsim)
    sims2[sims] <- TRUE
  } else if(is.logical(sims) && length(sims) == OM@nsim) {
    sims2 <- sims
  } else stop("Logical vector sims need to be of length ", OM@nsim)

  if(any(!sims2)) {
    message("Removing simulations: ", paste0(which(!sims2), collapse = " "))
    cpars <- OM@cpars

    subset_Data <- function(xx, Data, sims) {
      z <- slot(Data, xx)
      if(!all(is.na(z))) z <- z[sims, , , drop = FALSE]
      return(z)
    }
    subset_function <- function(x, sims) {
      if(is.matrix(x)) {
        return(x[sims, , drop = FALSE])
      } else if(is.array(x)) {
        if(length(dim(x)) == 3) return(x[sims, , , drop = FALSE])
        if(length(dim(x)) == 4) return(x[sims, , , , drop = FALSE])
        if(length(dim(x)) == 5) return(x[sims, , , , , drop = FALSE])
      } else if(class(x)[[1]] == "Data") {
        s_names <- c("AddIndV", "AddInd", "CV_AddInd")
        update_slots <- lapply(s_names, subset_Data, Data = x, sims = sims)
        for(i in 1:length(s_names)) slot(x, s_names[i]) <- update_slots[[i]]
        return(x)
      } else if(length(x) == OM@nsim) {
        return(x[sims])
      } else return(x)
    }
    OM@cpars <- lapply(cpars, subset_function, sims = sims2)
    OM@nsim <- sum(sims2)

    message("Set OM@nsim = ", OM@nsim)
  }

  return(OM)
}

par_identical_sims_fn <- function(StockPars, FleetPars, ObsPars, data, dots) {
  vector_fn <- function(x) sum(mean(x) - x) == 0
  array_fn <- function(x) {
    x_mean <- apply(x, 2:length(dim(x)), mean)
    all(apply(x, 1, identical, x_mean))
  }
  run_test <- function(x) if(is.null(dim(x))) vector_fn(x) else array_fn(x)

  StockPars_subset <- StockPars[c("hs", "procsd", "ageM", "M_ageArray", "Linf", "Len_age", "Wt_age", "Mat_age")]
  if(any(data$CAL > 0, na.rm = TRUE) || any(data$s_CAL > 0, na.rm = TRUE) || any(data$ML > 0, na.rm = TRUE)) {
    StockPars_subset <- c(StockPars_subset, StockPars["LenCV"])
  }
  S_test <- vapply(StockPars_subset, run_test, logical(1))

  if(data$nfleet == 1 && !any(data$CAL > 0, na.rm = TRUE) && !any(data$CAA > 0, na.rm = TRUE)) {
    FleetPars_subset <- FleetPars[c("L5", "LFS", "Vmaxlen")]
    FleetPars_subset <- lapply(FleetPars_subset, function(x) x[data$nyears, ])
    F_test_sel <- vapply(FleetPars_subset, run_test, logical(1))
  } else F_test_sel <- NULL

  if(data$condition == "effort" && !is.null(dots$OMeff) && dots$OMeff) {
    F_test_Find <- run_test(FleetPars$Find) %>% structure(names = "Find")
  } else F_test_Find <- NULL

  F_test <- c(F_test_sel, F_test_Find)

  if(data$nsurvey > 0 && !any(data$I_sd > 0, na.rm = TRUE)) {
    O_test <- run_test(ObsPars$Isd) %>% structure(names = "Isd")
  } else {
    O_test <- NULL
  }

  test_all <- c(S_test, F_test, O_test)
  return(test_all)
}


SRA_dynamic_SSB0 <- function(obj, par = obj$env$last.par.best) {

  if(obj$env$data$condition == "catch") {

    par[names(par) == "log_F" | names(par) == "log_F_equilibrium"] <- log(1e-8)
    out <- obj$report(par)$E

  } else if(obj$env$data$condition == "catch2") {

    new_args <- SRA_retro_subset(obj$env$data$n_y, data = obj$env$data, params = obj$env$parameters, map = obj$env$map)
    new_args$data$C_hist <- matrix(1e-8, new_args$data$n_y, new_args$data$nfleet)

    obj2 <- MakeADFun(data = new_args$data, parameters = new_args$params, map = new_args$map, random = obj$env$random,
                      DLL = "MSEtool", silent = TRUE)
    out <- obj2$report(par)$E

  } else {

    par[names(par) == "log_q_effort"] <- log(1e-8)
    out <- obj$report(par)$E

  }

  return(out)
}

SRA_posthoc_adjust <- function(report, obj, par = obj$env$last.par.best) {
  data <- obj$env$data
  report$F_at_age <- report$Z - data$M
  report$NPR_unfished <- do.call(rbind, report$NPR_unfished)
  length_bin <- obj$env$data$length_bin

  age_only_model <- data$len_age %>%
    apply(1, function(x) length(x) == data$max_age && max(x) == data$max_age) %>% all()
  if(age_only_model) {
    report$vul_len <- matrix(NA_real_, length(length_bin), data$nsel_block)
    report$s_vul_len <- matrix(NA_real_, length(length_bin), dim(report$s_vul)[3])

    report$MLpred <- matrix(NA_real_, nrow(report$F), ncol(report$F))
    report$CALpred <- array(NA_real_, dim(report$CALpred))
    report$s_CALpred <- array(NA_real_, dim(report$s_CALpred))
  } else {
    report$vul_len <- get_vul_len(report, data$vul_type, length_bin, data$Linf)
    report$s_vul_len <- get_s_vul_len(report, data$s_vul_type, length_bin, data$Linf)
  }
  report$dynamic_SSB0 <- SRA_dynamic_SSB0(obj, par)
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

get_s_vul_len <- function(report, s_selectivity, length_bin, Linf) {
  s_vul_len <- matrix(NA, length(length_bin), length(s_selectivity)) # length-based: matrix of dimension nlbin, nsurvey
  for(i in 1:ncol(s_vul_len)) {
    if(s_selectivity[i] == -1 || s_selectivity[i] == 0) {
      sls <- (report$s_LFS[i] - report$s_L5[i])/sqrt(-log(0.05, 2))
      srs <- (Linf - report$s_LFS[i])/sqrt(-log(report$s_Vmaxlen[i], 2))

      asc <- 2^-((length_bin - report$s_LFS[i])/sls * (length_bin - report$s_LFS[i])/sls)
      dsc <- ifelse(report$s_Vmaxlen[i] > rep(0.99, length(length_bin)), 1,
                    2^-((length_bin - report$s_LFS[i])/srs * (length_bin - report$s_LFS[i])/srs))
      s_vul_len[, i] <- ifelse(length_bin > report$s_LFS[i], dsc, asc)
    } else if(s_selectivity[i] > 0) {
      s_vul_len[, i] <- report$vul_len[, s_selectivity[i]]
    }
  }
  return(s_vul_len)
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
process_AddIndV <- function(sur, Misc, s_sel, maxage, nfleet, nyears) { # Return a matrix of nsim x nages
  if(s_sel[sur] < -2 || (s_sel[sur] == 1 & nfleet == 1)) { # -4 = B, -3 = SSB, single-fleet VB
    out <- matrix(1, length(Misc), maxage)
  } else { # custom sel or multi-fleet
    out <- do.call(rbind, lapply(Misc, function(x) x$s_vul[nyears, , sur]))
  }
  return(out)
}

SRA_retro <- function(x, nyr = 5) {
  if(length(x@mean_fit) == 0) stop("Re-run SRA_scope() with argument `mean_fit = TRUE`", .call = FALSE)

  data <- x@mean_fit$obj$env$data
  params <- x@mean_fit$obj$env$parameters
  n_y <- data$n_y
  map <- x@mean_fit$obj$env$map

  retro_ts <- array(NA, dim = c(nyr+1, n_y + 1, data$nfleet + 3))
  TS_var <- c(paste("F", 1:data$nfleet), "SSB", "SSB_SSB0", "R")
  dimnames(retro_ts) <- list(Peel = 0:nyr, Year = (x@OM@CurrentYr - n_y):x@OM@CurrentYr + 1, Var = TS_var)

  new_args <- lapply(n_y - 0:nyr, SRA_retro_subset, data = data, params = params, map = map)
  lapply_fn <- function(i, new_args, x) {
    obj2 <- MakeADFun(data = new_args[[i+1]]$data, parameters = new_args[[i+1]]$params, map = new_args[[i+1]]$map,
                      random = x@mean_fit$obj$env$random, DLL = "MSEtool", silent = TRUE)
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

    if(!is.character(opt2) && !is.character(SD)) {
      report <- obj2$report(obj2$env$last.par.best)

      if(new_args[[i+1]]$data$condition == "effort" && any(new_args[[i+1]]$data$Chist > 0, na.rm = TRUE)) {
        vars_div <- c("B", "E", "C_eq_pred", "CAApred", "CALpred", "s_CAApred", "s_CALpred", "CN", "Cpred", "N", "VB",
                      "R", "R_early", "R_eq", "R0", "B0", "E0", "N0", "E0_SR")
        vars_mult <- c("Brec", "q")
        var_trans <- c("R0", "q")
        fun_trans <- c("/", "*")
        rescale <- 1/exp(mean(log(new_args[[i+1]]$data$Chist/report$Cpred), na.rm = TRUE))
        fun_fixed <- c(NA, NA)
        rescale_report(vars_div, vars_mult, var_trans, fun_trans, fun_fixed)
      }

      FMort <- rbind(report$F, matrix(NA, i + 1, ncol(report$F)))
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

  retro <- new("retro", Model = "SRA_scope", Name = x@OM@Name, TS_var = TS_var, TS = retro_ts)
  attr(retro, "TS_lab") <- c(paste("Fishing mortality of Fleet", 1:data$nfleet),
                             "Spawning biomass", "Spawning depletion", "Recruitment")

  return(retro)
}


SRA_retro_subset <- function(yr, data, params, map) {
  ##### Data object
  data_out <- structure(data, check.passed = NULL)

  mat <- c("C_hist", "E_hist", "I_hist", "sigma_I", "CAA_n", "CAL_n", "s_CAA_n", "s_CAL_n", "msize", "M")
  mat_ind <- match(mat, names(data_out))
  data_out[mat_ind] <- lapply(data_out[mat_ind], function(x) x[1:yr, , drop = FALSE])

  mat2 <- c("len_age", "wt", "mat", "sel_block")
  mat_ind2 <- match(mat2, names(data_out))
  data_out[mat_ind2] <- lapply(data_out[mat_ind2], function(x) x[1:(yr+1), , drop = FALSE])

  # Update nsel_block, n_y
  data_out$nsel_block <- data_out$sel_block %>% as.vector() %>% unique() %>% length()
  data_out$n_y <- yr

  # Array 1:yr first index
  arr <- c("CAA_hist", "CAL_hist", "s_CAA_hist", "s_CAL_hist")
  arr_ind <- match(arr, names(data_out))
  data_out[arr_ind] <- lapply(data_out[arr_ind], function(x) x[1:yr, , , drop = FALSE])

  # Vector 1:yr
  data_out$est_rec_dev <- data_out$est_rec_dev[1:yr]

  ##### Parameters
  params_out <- structure(params, check.passed = NULL)
  for(i in 1:length(params)) {
    if(!is.null(attr(params[[i]], "map"))) params_out[[i]] <- attr(params[[i]], "shape")
  }

  # Update rec devs
  params_out$log_rec_dev <- params_out$log_rec_dev[1:yr]
  if(!is.null(map$log_rec_dev)) map$log_rec_dev <- map$log_rec_dev[1:yr] %>% factor()

  # Update F
  if(data_out$condition == "catch") {
    params_out$log_F <- params_out$log_F[1:yr, , drop = FALSE]

    if(any(data_out$yind_F + 1 > yr)) {
      data_out$yind_F <- as.integer(0.5 * yr)
      params_out$log_F[data_out$yind_F + 1, ] <- params$log_F[data$yind_F + 1, ]
    }
  }

  ## Update nsel block if needed
  if(data$nsel_block > data_out$nsel_block) {
    sel_block_ind <- data_out$sel_block %>% as.vector() %>% unique()
    params_out$vul_par <- params_out$vul_par[, sel_block_ind, drop = FALSE]
    if(!is.null(map$vul_par)) map$vul_par <- map$vul_par[, sel_block_ind, drop = FALSE] %>% factor()
  }

  return(list(data = data_out, params = params_out, map = map))
}
