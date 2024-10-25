

check_StockPars <- function(StockPars, nsim = 1, maxage, nyears) {
  
  vector_var <- c("SRrel", "R0", "Linf", "spawn_time_frac", "hs", "procsd")
  for (i in vector_var) {
    if (!length(StockPars[[i]])) stop(paste0("length(StockPars$", i, ") needs to be ", nsim))
    if (length(StockPars[[i]]) != nsim) StockPars[[i]] <- rep(StockPars[[i]], nsim)
  }
  
  vector_var <- "ageMarray"
  for (i in vector_var) {
    if (!length(StockPars[[i]])) stop(paste0("length(StockPars$", i, ") needs to be ", nsim))
    if (!is.matrix(StockPars[[i]])) StockPars[[i]] <- matrix(StockPars[[i]], nsim, 1)
  }
  
  array_var <- c("Len_age", "LatASD", "Wt_age", "Mat_age", "Fec_Age")
  for (i in array_var) {
    if (!length(StockPars[[i]])) stop(paste0("StockPars$", i, " needs to be a matrix [maxage + 1, nyears + 1]"))
    if (is.matrix(StockPars[[i]])) {
      if (!all(dim(StockPars[[i]]) == c(maxage + 1, nyears +1 ))) {
        stop(paste0("dim(StockPars$", i, ") needs to be [maxage + 1, nyears + 1]"))
      }
      StockPars[[i]] <- replicate(nsim, StockPars[[i]]) %>% aperm(c(3, 1, 2))
    }
  }
  
  array_var <- "M_ageArray"
  for (i in array_var) {
    if (!length(StockPars[[i]])) stop(paste0("StockPars$", i, " needs to be a matrix [maxage + 1, nyears]"))
    if (is.matrix(StockPars[[i]])) {
      if (!all(dim(StockPars[[i]]) == c(maxage + 1, nyears))) {
        stop(paste0("dim(StockPars$", i, ") needs to be [maxage + 1, nyears]"))
      }
      StockPars[[i]] <- replicate(nsim, StockPars[[i]]) %>% aperm(c(3, 1, 2))
    }
  }
  
  return(StockPars)
}



RCM_single_fit <- function(StockPars, RCMdata, condition = "catch", selectivity = "logistic", s_selectivity = "B", LWT = list(),
                           comp_like = c("multinomial", "lognormal", "mvlogistic", "dirmult1", "dirmult2"), prior = list(),
                           max_F = 3, integrate = FALSE,
                           control = list(iter.max = 2e+05, eval.max = 4e+05), 
                           start = list(), map = list(), silent = FALSE, ...) {
  
  dots <- list(...) # can be vul_par, ivul_par, log_rec_dev, log_early_rec_dev, map_vul_par, map_ivul_par, map_log_rec_dev, map_log_early_rec_dev, rescale, plusgroup, resample, OMeff, fix_dome
  if (!is.null(dots$maxF)) max_F <- dots$maxF
  
  comp_like <- match.arg(comp_like)
  condition <- match.arg(condition, choices = c("catch", "catch2", "effort"), several.ok = TRUE)
  
  # Check maxage
  if (length(RCMdata@CAA)) {
    if (is.matrix(RCMdata@CAA)) {
      maxage <- ncol(RCMdata@CAA) - 1
    } else if (is.array(RCM)) {
      maxage <- dim(RCMdata@CAA)[2] - 1 
    }
    RCMdata@Misc$maxage <- maxage
    if (!silent) message("Maximum age in model set to ", maxage)
  } else if (is.null(RCMdata@Misc$maxage)) {
    stop("No catch at age found to identify maximum age. Specify maxage in RCMdata@Misc$maxage.")
  }
  
  dat_update <- check_RCMdata(RCMdata, condition = condition, silent = silent)
  RCMdata <- dat_update$RCMdata
  
  maxage <- RCMdata@Misc$maxage
  nyears <- RCMdata@Misc$nyears
  nfleet <- RCMdata@Misc$nfleet
  nsurvey <- RCMdata@Misc$nsurvey
  if (is.null(RCMdata@Misc$CurrentYr)) RCMdata@Misc$CurrentYr <- nyears
  
  StockPars <- check_StockPars(StockPars, maxage = maxage, nyears = nyears)
  
  FleetPars <- list(
    L5_y = matrix(0.25 * StockPars$Linf[1], 1, RCMdata@Misc$nyears),
    LFS_y = matrix(0.5 * StockPars$Linf[1], 1, RCMdata@Misc$nyears),
    Vmaxlen_y = matrix(0.5, 1, RCMdata@Misc$nyears)
  )
  
  #nsim <- 1
  #proyears <- 0
  
  if (!silent) message("Maximum F in RCM will be ", max_F, ".\n\n")
  
  # No comp data
  if (!any(RCMdata@CAA > 0, na.rm = TRUE) && !any(RCMdata@CAL > 0, na.rm = TRUE) && is.null(start$vul_par)) {
    if (!silent) {
      message_info("No fishery length or age compositions were provided. Selectivity is fixed to default values. Use start argument to specify values.\n\n")
    }
    fix_sel <- TRUE
  } else {
    if (!silent) {
      message_info("Selectivity starting values set to generic default values. If desired, use start argument to specify values.\n\n")
    }
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
    message(switch(StockPars$SRrel[1],
                   "1" = "Beverton-Holt", 
                   "2" = "Ricker",
                   "3" = "Mesnil-Rochet"), " stock-recruitment relationship used.")
  }
  
  # Generate priors
  prior <- make_prior(prior, nsurvey, StockPars$SRrel[1], dots, msg = !silent)
  
  # Fit
  fit <- RCM_est(RCMdata = RCMdata, selectivity = sel, s_selectivity = s_sel,
                 LWT = RCMdata@Misc$LWT, comp_like = comp_like, prior = prior, 
                 max_F = max_F, integrate = integrate, StockPars = StockPars,
                 FleetPars = FleetPars, mean_fit = FALSE, control = control,
                 start = start, map = map, dots = dots)
  
  conv <- fit$report$conv
  highF <- any(fit$report$F >= max_F)
  if (highF) warning("Model had F on the upper boundary.\n")
  
  NaF <- any(is.na(fit$report$F) | is.infinite(fit$report$F))
  if (NaF) warning("Model had F with NA's")
  
  ### Output S4 object
  RCMdata@Misc$prior <- prior
  output <- new("RCModel", 
                OM = new("OM"), 
                SSB = matrix(fit$report$E, 1, nyears+1), 
                NAA = array(fit$report$N, c(1, nyears+1, maxage+1)),
                CAA = array(fit$report$CAApred, c(1, nyears+1, maxage+1, nfleet)), 
                mean_fit = fit, 
                conv = conv, 
                data = RCMdata)
  
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
    
    do_catch_check <- catch_check_fn(1, list(fit$report), RCMdata = RCMdata)
    if (do_catch_check) {
      warning("Note: there is predicted catch that deviates from observed catch by more than 1%.")
    }
  }
  if (!silent) message("Complete.")
  
  return(output)
}