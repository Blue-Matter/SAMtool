


#' @name simulate
#' @title Generate simulated data from TMB models in SAMtool
#' @description A convenient wrapper function (\code{simulate}) to simulate data (and process error) from the likelihood function.
#' 
#' @param object An object of class \linkS4class{Assessment} or \linkS4class{RCModel} containing the fitted model.
#' @author Q. Huynh
#' @details Process error, e.g., recruitment deviations, will be re-sampled in the simulation.
#' @return A list of length `nsim` containing simulated data. If refitting the model,
#' then a nested list of model output (`opt`, `SD`, and `report`).
#' @export
setGeneric("simulate", function(object, ...) standardGeneric("simulate"))



#' @rdname simulate 
#' @aliases simulate,Assessment-method simulate.Assessment
#' @param nsim Number of simulations
#' @param seed Used for the random number generator
#' @param refit Logical, whether to re-fit the model for each simulated dataset.
#' @param cores The number of CPUs for parallel processing for model re-fitting if \code{refit = TRUE}.
#' @param ... Additional arguments
#' 
#' @importFrom stats runif
#' @exportMethod simulate
setMethod("simulate", signature(object = "Assessment"),
          function(object, nsim = 1, seed = NULL, refit = FALSE, cores = 1, ...) {
            
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
            if(is.null(vars)) stop("Can not simulate data from model.", .call = FALSE)
            
            # Do simulation
            val <- lapply(1:nsim, Assess_sim, obj = object@obj, vars = vars) %>% 
              structure(names = paste0("sim_", 1:nsim))
            
            # Refit model from simulated model
            if(refit) {
              
              if(cores > 1 && !snowfall::sfIsRunning()) MSEtool::setup(cores)
              
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
                return(mod)
              }, cl = if(snowfall::sfIsRunning()) snowfall::sfGetCluster() else NULL) %>%
                structure(names = paste0("sim_", 1:nsim))
              
              attr(fit, "seed") <- RNGstate
              return(fit)
              
            } else {
              
              attr(val, "seed") <- RNGstate
              return(val)
              
            }
            
            return(invisible())
          })


#' @rdname simulate 
#' @aliases simulate,RCModel-method simulate.RCModel
#' @exportMethod simulate
setMethod("simulate", signature(object = "RCModel"),
          function(object, nsim = 1, seed = 12, refit = FALSE, cores = 1, ...) {
            
            if(!length(object@mean_fit)) stop("No TMB object was found. Re-run RCM with mean_fit = TRUE")
            
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
            val <- lapply(1:nsim, RCM_sim, obj = object@mean_fit$obj) %>% 
              structure(names = paste0("sim_", 1:nsim))
            
            # Refit model from simulated model
            if(refit) {
              
              if(cores > 1 && !snowfall::sfIsRunning()) MSEtool::setup(cores)
              
              fit <- pbapply::pblapply(1:nsim, function(x) {
                vars <- names(val)[[1]]
                
                newdata <- object@obj$env$data
                newdata[vars] <- val[[x]]
                newparams <- as.list(object@SD, what = "Estimate")
                
                obj2 <- MakeADFun(data = newdata, parameters = newparams, 
                                  random = object@obj$env$random, map = object@obj$env$map,
                                  DLL = "SAMtool", silent = TRUE)
                mod <- optimize_TMB_model(obj2, control = list(iter.max = 2e+05, eval.max = 4e+05), restart = 0)
                mod$report <- obj2$report(obj2$env$last.par.best) %>% RCM_posthoc_adjust(obj2)
                return(mod)
              }, cl = if(snowfall::sfIsRunning()) snowfall::sfGetCluster() else NULL) %>% 
                structure(names = paste0("sim_", 1:nsim))
              
              attr(fit, "seed") <- RNGstate
              return(fit)
              
            } else {
              
              attr(val, "seed") <- RNGstate
              return(val)
              
            }
            
            
            return(invisible())
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
  if(dist == "multinomial") {
    
    if(is.null(dots$N)) stop("N not provided for lognormal distribution")
    out[] <- rmultinom(n = 1, size = dots$N, prob = x)
    
  } else if(dist == "lognormal") {
    
    if(is.null(dots$tau)) stop("tau not provided for lognormal distribution")
    if(!requireNamespace("mvtnorm", quietly = TRUE)) stop("Please install the mvtnorm package.", call. = FALSE)
    # tau = sqrt(0.02/p_obs)
    samp <- mvtnorm::rmvnorm(n = 1, mean = log(x), sigma = dots$tau * diag(length(x)))
    out[] <- exp(samp)
    
  } else if(dist == "mvlogistic") {
    
    if(is.null(dots$tau)) stop("tau not provided for mvlogistic distribution")
    mu <- mean(log(x))
    samp <- mvtnorm::rmvnorm(n = 1, mean = log(x) - mu, sigma = dots$tau * diag(length(x)))
    out[] <- exp(samp + mu)
    
  } else if(dist == "dirmult1") {
    
    if(is.null(dots$N)) stop("N not provided for dirmult1 distribution")
    if(is.null(dots$theta)) stop("theta not provided for dirmult1 distribution")
    if(!requireNamespace("extraDistr", quietly = TRUE)) stop("Please install the extraDistr package.", call. = FALSE)
    
    out[] <- extraDistr::rdirmnom(n = 1, size = dots$N, alpha = dots$theta * dots$N * x/sum(x))
  } else {
    
    if(is.null(dots$beta)) stop("theta not provided for dirmult2 distribution")
    if(!requireNamespace("extraDistr", quietly = TRUE)) stop("Please install the extraDistr package.", call. = FALSE)
    
    out[] <- extraDistr::rdirmnom(n = 1, size = dots$N, alpha = beta * x/sum(x))
  }
  
  if(prop) out <- out/sum(out)
  
  return(out)
  
}


RCM_sim <- function(..., obj) {
  
  tmbvars <- c("C_eq", "C_hist", "I_hist", "msize")
  comp_vars <- c("CAA", "CAL", "IAA", "IAL")
  Rvars <- paste0(comp_vars, "_hist")
  
  report <- obj$simulate(obj$env$last.par.best) %>% RCM_posthoc_adjust(obj)
  
  res <- report[tmbvars]
  res[Rvars] <- obj$env$data[Rvars]
  
  comp_like <- obj$env$data$comp_like
  
  for(y in 1:obj$env$data$n_y) {
    for(ff in 1:obj$env$data$nfleet) {
      if(obj$env$data$CAA_n[y, ff] > 0) {
        
        if(comp_like == "lognormal") {
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
      
      if(obj$env$data$CAL_n[y, ff] > 0) {
        
        if(comp_like == "lognormal") {
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
      if(obj$env$data$IAA_n[y, sur] > 0) {
        
        if(comp_like == "lognormal") {
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
      
      if(obj$env$data$IAL_n[y, ff] > 0) {
        
        if(comp_like == "lognormal") {
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

Assess_sim <- function(..., obj, vars) {
  report <- obj$simulate(obj$env$last.par.best)
  res <- report[vars]
  comp_dist <- match.arg(obj$env$data$comp_dist, choices = c("multinomial", "lognormal"))
  
  if(any(vars == "CAA_hist") && is.null(res$CAA_hist)) {
    res$CAA_hist <- obj$env$data$CAA_hist
    
    for(y in 1:obj$env$data$n_y) {
      if(obj$env$data$CAA_n[y] > 0) {
        
        if(comp_dist == "lognormal") {
          obs <- res$CAA_hist[y, ]
          dispersion_par <- sqrt(0.02/(obs/sum(obs)))
        } else {
          dispersion_par <- NA
        }
        res$CAA_hist[y, ] <- simulate_comp(x = report$CAApred[y, ],
                                           dist = comp_dist,
                                           N = obj$env$data$CAA_n[y],
                                           tau = dispersion_par)
      }
    }
  }
  
  if(any(vars == "CAL_hist") && is.null(res$CAL_hist)) {
    res$CAL_hist <- obj$env$data$CAL_hist
    
    for(y in 1:obj$env$data$n_y) {
      if(obj$env$data$CAL_n[y] > 0) {
        
        if(comp_dist == "lognormal") {
          obs <- res$CAL_hist[y, ]
          dispersion_par <- sqrt(0.02/(obs/sum(obs)))
        } else {
          dispersion_par <- NA
        }
        res$CAL_hist[y, ] <- simulate_comp(x = report$CALpred[y, ],
                                           dist = comp_dist,
                                           N = obj$env$data$CAL_n[y],
                                           tau = dispersion_par)
      }
    }
  }
  
  return(res)
}

