


#' @name posterior
#' @title Sample posterior of TMB models in SAMtool
#' @description A convenient wrapper function (`posterior`) to sample the posterior using MCMC in rstan 
#' and returns a `stanfit` object for diagnostics. Use `RCMstan` to update the RCM and the enclosed operating model
#' with MCMC samples..
#'
#' @param x An object of class [Assessment-class] or [RCModel-class].
#' @param ... Additional arguments to pass to `rstan::sampling` via `tmbstan::tmbstan`.
#' @author Q. Huynh
#' @return `posterior` returns an object of class `stanfit`. See `class?stanfit`.
#' 
#' `RCMstan` returns an updated `RCModel`.
#' @export
setGeneric("posterior", function(x, ...) standardGeneric("posterior"))



#' @rdname posterior 
#' @aliases posterior,RCModel-method posterior.RCModel
#' @param priors_only Logical, whether to set the likelihood to zero and sample the priors only.
#' @param laplace Logical, whether to do the Laplace approximation for random parameters.
#' @param chains The numer of MCMC chains.
#' @param iter The number of iterations for each chain, including warmup.
#' @param warmup The number of burnin iterations
#' @param thin The frequency at which iterations are kept (e.g., `5` saves every fifth iteration)
#' @param seed Seed for random number generator during the MCMC.
#' @param init The initial values of parameters for starting the MCMC chain. See `tmbstan::tmbstan`.
#' @param cores The number of cores for running in parallel, e.g., one core per MCMC chain. Used in `RCMstan`
#' for reconstructing the population.
#' 
#' @section Online Documentation:
#' A vignette on the steps to run the MCMC is available on the openMSE 
#' [website](https://openmse.com/tutorial-rcm/4-case-study-mcmc/).
#' @exportMethod posterior
setMethod("posterior", signature(x = "RCModel"),
          function(x, priors_only = FALSE, 
                   laplace = FALSE, chains = 2, iter = 2000, warmup = floor(iter/2), thin = 5,
                   seed = 34, init = "last.par.best", cores = chains, ...) {
            if (!requireNamespace("tmbstan", quietly = TRUE)) stop("Install tmbstan to use this function.")
            
            xchar <- substitute(x) %>% as.character()
            obj <- x@mean_fit$obj
            
            if (is.null(obj)) {
              stop(paste0("No TMB object was found in ", xchar, "@mean_fit$obj. Re-run RCM with mean_fit = TRUE."))
            }
            
            if (priors_only) {
              newdata <- structure(obj$env$data, check.passed = NULL)
              newdata$LWT_fleet[] <- 0
              newdata$LWT_index[] <- 0
              
              obj_new <- MakeADFun(data = newdata, 
                                   parameters = clean_tmb_parameters(obj), map = obj$env$map, 
                                   random = obj$env$random, DLL = obj$env$DLL, silent = obj$env$silent)
              
              obj_new$env$last.par <- obj$env$last.par
              obj_new$env$last.par.best <- obj$env$last.par.best
              obj_new$env$last.par.ok <- obj$env$last.par.ok
              
              obj_new$env$last.par1 <- obj$env$last.par1
              obj_new$env$last.par2 <- obj$env$last.par2
              
            } else {
              obj_new <- obj
            }
            
            lower <- rep(-Inf, length(obj$par))
            upper <- rep(Inf, length(obj$par))
            
            R0_prior <- !is.null(obj$env$data$use_prior[1]) &&
              names(obj$env$data$use_prior)[1] == "R0" && 
              obj$env$data$use_prior["R0"] > 1
            
            if (any(names(obj$par) == "R0x") && obj$env$data$use_prior["R0"] > 1) { # Uniform priors need bounds
              lower[names(obj$par) == "R0x"] <- log(obj$env$data$prior_dist[1, 1]) + log(obj$env$data$rescale)
              upper[names(obj$par) == "R0x"] <- log(obj$env$data$prior_dist[1, 2]) + log(obj$env$data$rescale)
            }
            
            tmbstan::tmbstan(obj_new, chains = chains, iter = iter, warmup = warmup, thin = thin,
                             seed = seed, init = init, cores = cores, lower = lower, upper = upper, ...)
          })


#' @rdname posterior
#' @aliases posterior,Assessment-method posterior.Assessment
#' @exportMethod posterior
setMethod("posterior", signature(x = "Assessment"),
          function(x, priors_only = FALSE, ...) {
            stop("Posterior sampling not set up for assessment models.")
            if (!requireNamespace("tmbstan", quietly = TRUE)) stop("Install tmbstan to use this function.")
            xchar <- substitute(x) %>% as.character()
            f <- get(paste0("posterior_", x@Model))
            f(x, ...)
          })


#' @rdname posterior
#' @aliases RCMstan
#' @param RCModel An object of class `RCModel`
#' @param stanfit An object of class `stanfit` returned by `posterior`.
#' @param sim A matrix of `RCModel@OM@nsim` rows and 2 columns that specifies the samples used to update the
#' operating model. The first column specifies the chain and the second columns specifies the MCMC iteration.
#' @param silent Logical to indicate if progress messages should be printed to console.
#' @export
RCMstan <- function(RCModel, stanfit, sim, cores = 1, silent = FALSE) {
  OM <- RCModel@OM
  nsim <- OM@nsim
  
  nchains <- stanfit@sim$chains
  nsim_stan <- length(stanfit@sim$samples[[1]][[1]])
  
  if (missing(sim)) {
    sim_stan <- matrix(0, nsim, 2)
    sim_stan[, 1] <- sample(1:nchains, nsim, replace = TRUE)
    sim_stan[, 2] <- sample(1:nsim_stan, nsim, replace = nsim_stan < nsim)
  } else if (!is.matrix(sim)) {
    sim_stan <- matrix(1, nsim, 2)
    sim_stan[, 2] <- sim
  } else {
    sim_stan <- sim
  }
  
  if (any(sim_stan[, 1] > nchains)) stop("There is only ", nchains, " chain(s) in the stan model.")
  if (any(sim_stan[, 2] > nsim_stan)) stop("There are only ", nsim_stan, " iterations in the stan model.")
  
  # Get samps
  if (!silent) message_info("Sampling ", nsim, " iterations for the operating model...")
  samps <- sapply(1:nsim, function(x) {
    chain <- sim_stan[x, 1]
    sim_out <- sim_stan[x, 2]
    pars <- sapply(stanfit@sim$samples[[chain]], getElement, sim_out)
  }) %>% t()
  
  # Get model output
  obj <- RCModel@mean_fit$obj
  if (is.null(obj)) stop("No TMB object found in RCModel@mean_fit.")
    
  if (!silent) message_info("Re-constructing population model from MCMC parameter samples...")
  
  if (cores > 1 && !snowfall::sfIsRunning()) {
    MSEtool::setup(as.integer(cores))
    on.exit(snowfall::sfStop())
  }   
  res <- pblapply(1:nsim, RCM_report_samps, samps = samps[, -ncol(samps)], obj = obj, conv = TRUE, 
                  cl = if (snowfall::sfIsRunning()) snowfall::sfGetCluster() else NULL)
  
  if (!silent) message("Updating operating model with MCMC samples...\n")
  newOM <- RCM_update_OM(
    OM = OM, 
    report = res, 
    StockPars = list(Perr_y = matrix(1, nsim, OM@maxage + OM@nyears + OM@proyears)), 
    obj_data = obj$env$data, 
    maxage = OM@maxage, 
    nyears = OM@nyears, 
    proyears = OM@proyears, 
    prior = list(use_prior = obj$env$data$use_prior), 
    silent = silent
  )
  
  RCModel@OM <- newOM$OM
  RCModel@SSB <- newOM$RCM_val$SSB
  RCModel@NAA <- newOM$RCM_val$NAA
  RCModel@CAA <- newOM$RCM_val$CAA
  if (!is.null(newOM$RCM_val$CAL)) RCModel@CAL <- newOM$RCM_val$CAL
  RCModel@conv <- rep(TRUE, nsim)
  RCModel@Misc <- res
  RCModel@config$drop_sim <- integer(0)
  
  if (!silent) message("Finished.")
  
  return(RCModel)
}
