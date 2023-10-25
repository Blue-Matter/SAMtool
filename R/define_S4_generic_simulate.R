

#' Class-`sim`
#'
#' An S4 class that contains output from [simulate].
#'
#' @name sim-class
#' @docType class
#'
#' @slot Model Name of the assessment model.
#' @slot data List of data from the assessment. 
#' @slot data_sim List of simulated data values. Each value returns an array.
#' @slot process_sim List of simulated process error.
#' @slot est Estimates from the original model fit.
#' @slot est_sim Estimates from the simulated data.
#' @author Q. Huynh
#' @export sim
#' @exportClass sim
sim <- setClass("sim", slots = c(Model = "character", 
                                 data = "list",
                                 data_sim = "list",
                                 process_sim = "list",
                                 est = "list",
                                 est_sim = "list"))


#' @name simulate
#' @title Generate simulated data from TMB models in SAMtool
#' @description A convenient wrapper function (`simulate`) to simulate data (and process error) from the likelihood function.
#' 
#' @param object An object of class [Assessment-class] or [RCModel-class] containing the fitted model.
#' @author Q. Huynh
#' @details Process error, e.g., recruitment deviations, will be re-sampled in the simulation.
#' @return A [sim-class] object returning the original data, simulated data, original parameters, parameters estimated
#' from simulated data, and process error used to simulate data.
#' then a nested list of model output (`opt`, `SD`, and `report`).
#' @export
setGeneric("simulate", function(object, ...) standardGeneric("simulate"))



#' @rdname simulate 
#' @aliases simulate,Assessment-method simulate.Assessment
#' @param nsim Number of simulations
#' @param seed Used for the random number generator
#' @param process_error Logical, indicates if process error is re-sampled in the simulation.
#' @param refit Logical, whether to re-fit the model for each simulated dataset.
#' @param cores The number of CPUs for parallel processing for model re-fitting if `refit = TRUE`.
#' @param ... Additional arguments
#' 
#' @importFrom stats runif
#' @exportMethod simulate
setMethod("simulate", signature(object = "Assessment"),
          function(object, nsim = 1, seed = NULL, process_error = FALSE, 
                   refit = FALSE, cores = 1, ...) {
            
            if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
            if (is.null(seed)) {
              RNGstate <- get(".Random.seed", envir = .GlobalEnv)
            } else {
              R.seed <- get(".Random.seed", envir = .GlobalEnv)
              set.seed(seed)
              RNGstate <- structure(seed, kind = as.list(RNGkind()))
              on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
            }
            
            vars <- switch(object@Model,
              "cDD" = c("C_hist", "I_hist", "MW_hist"),
              "cDD_SS" = c("C_hist", "I_hist", "MW_hist"),
              "DD_TMB" = c("C_hist", "I_hist", "MW_hist"),
              "DD_SS" = c("C_hist", "I_hist", "MW_hist"),
              "SP" = c("C_hist", "I_hist"),
              "SP_SS" = c("C_hist", "I_hist"),
              "VPA" = "I_hist",
              "SCA" = c("C_hist", "I_hist", "CAA_hist", "CAL_hist")
            )
            if (is.null(vars)) stop("Can not simulate data from model.", .call = FALSE)
            
            process_vars <- switch(object@Model,
                                   "cDD" = NULL,
                                   "cDD_SS" = "log_rec_dev_sim",
                                   "DD_TMB" = NULL,
                                   "DD_SS" = "log_rec_dev_sim",
                                   "SP" = NULL,
                                   "SP_SS" = "log_B_dev_sim",
                                   "VPA" = NULL,
                                   "SCA" = c("log_rec_dev_sim", "log_early_rec_dev_sim", "logit_M_sim", "logit_M_walk_sim")
            )
            if (process_error && is.null(process_vars)) message_info("No process error found for this model.")
            
            # Do simulation
            val <- lapply(1:nsim, Assess_sim, obj = object@obj, vars = vars, 
                          process_error = process_error, process_vars = process_vars) %>% 
              structure(names = paste0("sim_", 1:nsim))
            
            # Refit model from simulated data
            if (refit) {
              if (cores > 1 && !snowfall::sfIsRunning()) MSEtool::setup(cores)
              
              fit <- pbapply::pblapply(1:nsim, function(x) {
                newdata <- object@obj$env$data
                newdata[vars] <- val[[x]][vars]
                #newparams <- as.list(object@SD, what = "Estimate")
                newparams <- object@info$params
                obj2 <- MakeADFun(data = newdata, parameters = newparams, 
                                  random = object@obj$env$random, map = object@obj$env$map,
                                  DLL = "SAMtool", silent = TRUE)
                mod <- optimize_TMB_model(obj2, control = object@info$control, restart = 0)
                mod$report <- obj2$report(obj2$env$last.par.best)
                mod$obj <- obj2
                return(mod)
              }, cl = if (snowfall::sfIsRunning()) snowfall::sfGetCluster() else NULL) %>%
                structure(names = paste0("sim_", 1:nsim))
              
            } else {
              fit <- NULL
            }
            
            output <- create_sim_object(object, fit, val, vars, process_vars)
            attr(output, "seed") <- RNGstate
            return(output)
          })


#' @rdname simulate 
#' @aliases simulate,RCModel-method simulate.RCModel
#' @exportMethod simulate
setMethod("simulate", signature(object = "RCModel"),
          function(object, nsim = 1, seed = NULL, process_error = FALSE, refit = FALSE, cores = 1, ...) {
            
            if (!length(object@mean_fit)) stop("No TMB object was found. Re-run RCM with mean_fit = TRUE")
            
            if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
            if (is.null(seed)) {
              RNGstate <- get(".Random.seed", envir = .GlobalEnv)
            } else {
              R.seed <- get(".Random.seed", envir = .GlobalEnv)
              set.seed(seed)
              RNGstate <- structure(seed, kind = as.list(RNGkind()))
              on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
            }
            
            # Do simulation
            tmbvars <- c("C_eq", "C_hist", "I_hist", "msize")
            process_vars <- c("log_rec_dev", "log_rec_dev_sim")
            comp_vars <- c("CAA", "CAL", "IAA", "IAL")
            Rvars <- paste0(comp_vars, "_hist")
            
            val <- lapply(1:nsim, RCM_sim, obj = object@mean_fit$obj, process_error = process_error) %>% 
              structure(names = paste0("sim_", 1:nsim))
            
            # Refit model from simulated model
            if (refit) {
              
              if (cores > 1 && !snowfall::sfIsRunning()) MSEtool::setup(cores)
              
              fit <- pbapply::pblapply(1:nsim, function(x) {
                newdata <- object@mean_fit$obj$env$data
                newdata[c(tmbvars, Rvars)] <- val[[x]][c(tmbvars, Rvars)]
                newparams <- as.list(object@mean_fit$SD, what = "Estimate")
                
                obj2 <- MakeADFun(data = newdata, parameters = newparams, 
                                  random = object@mean_fit$obj$env$random, map = object@mean_fit$obj$env$map,
                                  DLL = "SAMtool", silent = TRUE)
                obj2$par[] <- object@mean_fit$obj$par
                mod <- optimize_TMB_model(obj2, control = list(iter.max = 2e+05, eval.max = 4e+05), restart = 0)
                mod$report <- obj2$report(obj2$env$last.par.best) %>% RCM_posthoc_adjust(obj2)
                mod$obj <- obj2
                return(mod)
              }, cl = if (snowfall::sfIsRunning()) snowfall::sfGetCluster() else NULL) %>% 
                structure(names = paste0("sim_", 1:nsim))
              
            } else {
              fit <- NULL
            }
            
            output <- create_sim_object(object, fit, val, 
                                        vars = c(tmbvars, Rvars), 
                                        process_vars = process_vars)
            attr(output, "seed") <- RNGstate
            return(output)
          })


# #' @param x Predicted vector
# #' @param ... multinomial: N, lognormal: tau = sqrt(0.02/p_obs), mvlogistic: tau, dirmult1: N, theta, dirmult2: beta
#' @importFrom stats rmultinom
simulate_comp <- function(x, 
                          prop = TRUE,
                          dist = c("multinomial", "lognormal", "mvlogistic", "dirmult1", "dirmult2"),
                          ...) {
  
  dist <- match.arg(dist)
  dots <- list(...)
  
  out <- numeric(length(x))
  if (dist == "multinomial") {
    
    if (is.null(dots$N)) stop("N not provided for lognormal distribution")
    out[] <- rmultinom(n = 1, size = dots$N, prob = x)
    
  } else if (dist == "lognormal") {
    
    if (is.null(dots$tau)) stop("tau not provided for lognormal distribution")
    if (!requireNamespace("mvtnorm", quietly = TRUE)) stop("Please install the mvtnorm package.", call. = FALSE)
    # tau = sqrt(0.02/p_obs)
    samp <- mvtnorm::rmvnorm(n = 1, mean = log(x), sigma = dots$tau * diag(length(x)))
    out[] <- exp(samp)
    
  } else if (dist == "mvlogistic") {
    
    if (is.null(dots$tau)) stop("tau not provided for mvlogistic distribution")
    mu <- mean(log(x))
    samp <- mvtnorm::rmvnorm(n = 1, mean = log(x) - mu, sigma = dots$tau * diag(length(x)))
    out[] <- exp(samp + mu)
    
  } else if (dist == "dirmult1") {
    
    if (is.null(dots$N)) stop("N not provided for dirmult1 distribution")
    if (is.null(dots$theta)) stop("theta not provided for dirmult1 distribution")
    if (!requireNamespace("extraDistr", quietly = TRUE)) stop("Please install the extraDistr package.", call. = FALSE)
    
    out[] <- extraDistr::rdirmnom(n = 1, size = dots$N, alpha = dots$theta * dots$N * x/sum(x))
  } else {
    
    if (is.null(dots$beta)) stop("theta not provided for dirmult2 distribution")
    if (!requireNamespace("extraDistr", quietly = TRUE)) stop("Please install the extraDistr package.", call. = FALSE)
    
    out[] <- extraDistr::rdirmnom(n = 1, size = dots$N, alpha = beta * x/sum(x))
  }
  
  if (prop) out <- out/sum(out)
  
  return(out)
  
}


RCM_sim <- function(..., obj, process_error = FALSE) {
  
  tmbvars <- c("C_eq", "C_hist", "I_hist", "msize")
  process_vars <- c("log_rec_dev", "log_rec_dev_sim")
  comp_vars <- c("CAA", "CAL", "IAA", "IAL")
  Rvars <- paste0(comp_vars, "_hist")
  
  newdata <- obj$env$data %>% structure(check.passed = NULL)
  
  if (process_error && any(names(newdata) == "sim_process_error")) {
    
    newdata$sim_process_error <- 1L
    newparams <- clean_tmb_parameters(obj)
    
    obj2 <- MakeADFun(data = newdata, parameters = newparams,
                      random = obj$env$random, map = obj$env$map,
                      DLL = "SAMtool", silent = obj$env$silent)
    report <- obj2$simulate(obj$env$last.par.best) %>% RCM_posthoc_adjust(obj)
  } else {
    report <- obj$simulate(obj$env$last.par.best) %>% RCM_posthoc_adjust(obj)
  }
  
  res <- report[c(tmbvars, process_vars)]
  res[Rvars] <- obj$env$data[Rvars]
  
  comp_like <- obj$env$data$comp_like
  
  for(y in 1:obj$env$data$n_y) {
    for(ff in 1:obj$env$data$nfleet) {
      if (obj$env$data$CAA_n[y, ff] > 0) {
        
        if (comp_like == "lognormal") {
          obs <- res$CAA_hist[y, , ff]
          dispersion_par <- sqrt(0.02/(obs/sum(obs)))
        } else {
          dispersion_par <- ifelse(is.null(report$log_compf), NA, exp(report$log_compf[ff, 1]))
        }
        res$CAA_hist[y, , ff] <- simulate_comp(x = report$CAApred[y, , ff],
                                               dist = comp_like,
                                               N = obj$env$data$CAA_n[y, ff],
                                               tau = dispersion_par,
                                               theta = dispersion_par,
                                               beta = dispersion_par)
      }
      
      if (obj$env$data$CAL_n[y, ff] > 0) {
        
        if (comp_like == "lognormal") {
          obs <- res$CAL_hist[y, , ff]
          dispersion_par <- sqrt(0.02/(obs/sum(obs)))
        } else {
          dispersion_par <- ifelse(is.null(report$log_compf), NA, exp(report$log_compf[ff, 2]))
        }
        res$CAL_hist[y, , ff] <- simulate_comp(x = report$CALpred[y, , ff],
                                               dist = comp_like,
                                               N = obj$env$data$CAL_n[y, ff],
                                               tau = dispersion_par,
                                               theta = dispersion_par,
                                               beta = dispersion_par)
      }
    }
    
    for(sur in 1:obj$env$data$nsurvey) {
      if (obj$env$data$IAA_n[y, sur] > 0) {
        
        if (comp_like == "lognormal") {
          obs <- res$IAA_hist[y, , sur]
          dispersion_par <- sqrt(0.02/(obs/sum(obs)))
        } else {
          dispersion_par <- ifelse(is.null(report$log_compi), NA, exp(report$log_compi[sur, 1]))
        }
        res$IAA_hist[y, , sur] <- simulate_comp(x = report$IAApred[y, , sur],
                                                dist = comp_like,
                                                N = obj$env$data$IAA_n[y, sur],
                                                tau = dispersion_par,
                                                theta = dispersion_par,
                                                beta = dispersion_par)
      }
      
      if (obj$env$data$IAL_n[y, ff] > 0) {
        
        if (comp_like == "lognormal") {
          obs <- res$IAL_hist[y, , sur]
          dispersion_par <- sqrt(0.02/(obs/sum(obs)))
        } else {
          dispersion_par <- ifelse(is.null(report$log_compi), NA, exp(report$log_compi[sur, 2]))
        }
        res$IAL_hist[y, , sur] <- simulate_comp(x = report$IALpred[y, , sur],
                                                dist = comp_like,
                                                N = obj$env$data$IAL_n[y, sur],
                                                tau = dispersion_par,
                                                theta = dispersion_par,
                                                beta = dispersion_par)
      }
    }
  }
  
  return(res)
}

Assess_sim <- function(..., obj, vars, process_error = FALSE, process_vars = NULL) {
  
  newdata <- obj$env$data %>% structure(check.passed = NULL)
  
  if (process_error && any(names(newdata) == "sim_process_error")) {
    
    newdata$sim_process_error <- 1L
    newparams <- clean_tmb_parameters(obj)
    
    obj2 <- MakeADFun(data = newdata, parameters = newparams,
                      random = obj$env$random, map = obj$env$map,
                      DLL = "SAMtool", silent = obj$env$silent)
    report <- obj2$simulate(obj$env$last.par.best)
    
  } else {
    report <- obj$simulate(obj$env$last.par.best) 
  }
  res <- report[c(vars, process_vars)]
  
  if (!is.null(newdata$comp_dist)) {
    comp_dist <- match.arg(newdata$comp_dist, choices = c("multinomial", "lognormal"))
  }
  if (any(vars == "CAA_hist") && is.null(res$CAA_hist)) {
    res$CAA_hist <- newdata$CAA_hist
    
    for(y in 1:newdata$n_y) {
      if (newdata$CAA_n[y] > 0) {
        
        if (comp_dist == "lognormal") {
          obs <- res$CAA_hist[y, ]
          dispersion_par <- sqrt(0.02/(obs/sum(obs)))
        } else {
          dispersion_par <- NA
        }
        res$CAA_hist[y, ] <- simulate_comp(x = report$CAApred[y, ],
                                           dist = comp_dist,
                                           N = newdata$CAA_n[y],
                                           tau = dispersion_par)
      }
    }
  }
  
  if (any(vars == "CAL_hist") && is.null(res$CAL_hist)) {
    res$CAL_hist <- newdata$CAL_hist
    
    for(y in 1:newdata$n_y) {
      if (newdata$CAL_n[y] > 0) {
        
        if (comp_dist == "lognormal") {
          obs <- res$CAL_hist[y, ]
          dispersion_par <- sqrt(0.02/(obs/sum(obs)))
        } else {
          dispersion_par <- NA
        }
        res$CAL_hist[y, ] <- simulate_comp(x = report$CALpred[y, ],
                                           dist = comp_dist,
                                           N = newdata$CAL_n[y],
                                           tau = dispersion_par)
      }
    }
  }
  
  return(res)
}

create_sim_object <- function(object, fit = NULL, val, vars, process_vars) {
  
  data_sim <- lapply(vars, function(v) {
    sapply(val, getElement, v, simplify = "array")
  }) %>% structure(names = vars)
  
  if (inherits(object, "Assessment")) {
    obj <- object@obj
  } else if (inherits(object, "RCModel")) {
    obj <- object@mean_fit$obj
  }
  sim_out <- new("sim",
                 data = obj$env$data[vars],
                 data_sim = data_sim)
  
  if (!is.null(process_vars)) {
    sim_out@process_sim <- lapply(process_vars, function(v) {
      sapply(val, getElement, v, simplify = "array")
    }) %>% structure(names = process_vars)
  }
  
  if (!is.null(fit)) {
    
    if (inherits(object, "Assessment")) {
      model <- object@Model
      TMB_report <- object@TMB_report
    } else if (inherits(object, "RCModel")) {
      model <- "RCM"
      TMB_report <- object@mean_fit$report
    }
    
    est_vars <- switch(model,
                       "cDD" = c("R0", "h", "B", "F", "R"),
                       "cDD_SS" = c("R0", "h", "B", "F", "R"),
                       "DD_TMB" = c("R0", "h", "B", "F", "R"),
                       "DD_SS" = c("R0", "h", "B", "F", "R"),
                       "SP" = c("FMSY", "MSY", "r", "K", "B", "F"),
                       "SP_SS" = c("FMSY", "MSY", "r", "K", "B", "F"),
                       "VPA" = c("B", "VB", "F"),
                       "SCA" = c("R0", "h", "B", "E", "F", "R", "VB"),
                       "RCM" = NULL
                       )
    
    sim_out@est <- TMB_report[est_vars]
    
    sim_out@est_sim <- lapply(est_vars, function(v) {
      sapply(fit, function(x) x$report[[v]], simplify = "array")
    }) %>% structure(names = est_vars)
    
  }
  
  return(sim_out)
}
