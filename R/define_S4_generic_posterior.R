


#' @name posterior
#' @title Sample posterior of TMB models in SAMtool
#' @description A convenient wrapper function (\code{posterior}) to sample the posterior using MCMC in rstan 
#' and returns a \code{stanfit} object for diagnostics. Use \code{RCMstan} to update the RCM and the enclosed operating model
#' with MCMC samples..
#'
#' @param x An object of class \linkS4class{Assessment} or \linkS4class{RCModel}.
#' @param ... Additional arguments to pass to \code{rstan::sampling} via \code{tmbstan::tmbstan}.
#' @details
#'
#' @author Q. Huynh
#' @return \code{posterior} returns an object of class \code{stanfit}. See \code{class?stanfit}.
#' 
#' \code{RCMstan} returns an updated \code{RCModel}.
#' @export
setGeneric("posterior", function(x, ...) standardGeneric("posterior"))



#' @rdname posterior 
#' @aliases posterior,RCModel-method posterior.RCModel
#' @param priors_only Logical, whether to set the likelihood to zero and sample the priors only.
#' @param laplace Logical, whether to do the Laplace approximation for random parameters.
#' @param chains The numer of MCMC chains.
#' @param iter The number of iterations for each chain, including warmup.
#' @param warmup The number of burnin iterations
#' @param thin The frequency at which iterations are kept (e.g., \code{5} saves every fifth iteration)
#' @param seed Seed for random number generator during the MCMC.
#' @param init The initial values of parameters for starting the MCMC chain. See \code{tmbstan::tmbstan}.
#' @param cores The number of cores for running in parallel, e.g., one core per MCMC chain. Used in \code{RCMstan}
#' for reconstructing the population.
#' @exportMethod posterior
setMethod("posterior", signature(x = "RCModel"),
          function(x, priors_only = FALSE, 
                   laplace = FALSE, chains = 2, iter = 2000, warmup = floor(iter/2), thin = 5,
                   seed = 34, init = "last.par.best", cores = chains, ...) {
            if(!requireNamespace("tmbstan", quietly = TRUE)) stop("Install tmbstan to use this function.")
            
            xchar <- substitute(x) %>% as.character()
            
            obj <- x@mean_fit$obj
            
            if(is.null(obj)) {
              stop(paste0("No TMB object was found in ", xchar, "@mean_fit$obj. Re-run RCM with mean_fit = TRUE."))
            } else if(obj$env$data$condition != "catch2") { # With condition = "catch", infinitely diffuse prior on log_F_dev
              stop("MCMC currently only supported with condition = \"catch2\".")
            }
            
            if(priors_only) {
              
              data <- structure(obj$env$data, check.passed = NULL)
              params <- lapply(obj$env$parameters, function(x) if(!is.null(attr(x, "map"))) attr(x, "shape") else x)
              
              data$LWT_fleet[] <- 0
              data$LWT_index[] <- 0
              obj_new <- MakeADFun(data = data, parameters = params, map = obj$env$map, 
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
            
            if(!is.null(obj$env$data$use_prior[1]) && obj$env$data$use_prior[1] > 0) { # Uniform priors need bounds
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
            if(!requireNamespace("tmbstan", quietly = TRUE)) stop("Install tmbstan to use this function.")
            
            xchar <- substitute(x) %>% as.character()
            f <- get(paste0("posterior_", x@Model))
            f(x, xchar, chains, iter, warmup, seed)
          })


#' @rdname posterior
#' @aliases RCMstan
#' @param RCModel An object of class \code{RCModel}
#' @param stanfit An object of class \code{stanfit} returned by \code{posterior}.
#' @param sim A matrix of \code{RCModel@OM@nsim} rows and 2 columns that specifies the samples used to update the
#' operating model. The first column specifies the chain and the second columns specifies the MCMC iteration.
#' @export
RCMstan <- function(RCModel, stanfit, sim, cores = 1) {
  nsim <- RCModel@OM@nsim
  
  nchains <- stanfit@sim$chains
  nsim_stan <- length(stanfit@sim$samples[[1]][[1]])
  
  if(missing(sim)) {
    sim_stan <- matrix(0, nsim, 2)
    sim_stan[, 1] <- sample(1:nchains, nsim, replace = TRUE)
    sim_stan[, 2] <- sample(1:nsim_stan, nsim, replace = nsim_stan < nsim)
  } else if(!is.matrix(sim)) {
    sim_stan <- matrix(1, nsim, 2)
    sim_stan[, 2] <- sim
  }
  
  if(any(sim_stan[, 1] > nchains)) stop("There is only ", nchains, " chain(s) in the stan model.")
  if(any(sim_stan[, 2] > nsim_stan)) stop("There are only ", nsim_stan, " iterations in the stan model.")
  
  # Get samps
  message("Sampling ", nsim, " iterations for the operating model...")
  samps <- sapply(1:nsim, function(x) {
    chain <- sim_stan[x, 1]
    sim_out <- sim_stan[x, 2]
    pars <- sapply(stanfit@sim$samples[[chain]], getElement, sim_out)
  }) %>% t()
  
  # Get model output
  OM <- RCModel@OM
  maxage <- OM@maxage
  nyears <- OM@nyears
  proyears <- OM@proyears
  obj <- RCModel@mean_fit$obj
  if(is.null(obj)) stop("No TMB object found in RCModel@mean_fit.")
  
  message("Re-constructing population model...")
  
  if(cores > 1 && !snowfall::sfIsRunning()) MSEtool::setup(as.integer(cores))
  if(snowfall::sfIsRunning()) {
    res <- sfLapply(1:nsim, RCM_report_samps, samps = samps[, -ncol(samps)], obj = obj, conv = TRUE)
  } else {
    res <- lapply(1:nsim, RCM_report_samps, samps = samps[, -ncol(samps)], obj = obj, conv = TRUE)
  }
  
  OM_par <- RCM_update_OM(res, obj$env$data, maxage, nyears, proyears)
  message("Updating operating model with MCMC samples...\n")
  
  ### R0
  OM@cpars$R0 <- OM_par$R0
  message("Range of unfished age-0 recruitment (OM@cpars$R0): ", paste(round(range(OM@cpars$R0), 2), collapse = " - "))
  
  ### Depletion and init D - init D is only reported, OM setup for initD by adjusting rec devs
  message("Range of initial spawning depletion: ", paste(round(range(OM_par$initD), 2), collapse = " - "))
  
  OM@cpars$D <- OM_par$D
  message("Range of spawning depletion (OM@cpars$D): ", paste(round(range(OM@cpars$D), 2), collapse = " - "), "\n")
  
  ### Selectivity and F
  ### Find
  OM@isRel <- FALSE
  OM@cpars$V <- OM_par$V
  OM@cpars$Find <- OM_par$Find
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
  OM@cpars$Perr <- OM_par$procsd
  message("Recruitment standard deviation set in OM@cpars$Perr.")
  
  Perr_y <- OM@cpars$Perr_y
  Perr_y[, 1:maxage] <- OM_par$early_Perr
  Perr_y[, maxage + 1:nyears] <- OM_par$Perr
  message("Historical recruitment set in OM@cpars$Perr_y.")
  
  if(any(OM_par$AC != 0)) {
    OM@cpars$AC <- OM_par$AC
    OM@AC <- range(OM_par$AC)
    message("Range of recruitment autocorrelation OM@AC: ", paste(round(range(OM@AC), 2), collapse = " - "))
    
    OM@cpars$Perr_y <- RCM_sample_future_dev(obj$env$data$est_rec_dev, OM_par$procsd, OM_par$AC, 
                                             OM_par$log_rec_dev, Perr_y, maxage, nyears, proyears)
    message("Future recruitment deviations sampled with autocorrelation (in OM@cpars$Perr_y).\n")
  }
  
  prior <- list(use_prior = obj$env$data$use_prior)
  if(prior$use_prior[2]) OM@cpars$hs <- OM_par$h
  if(prior$use_prior[3]) OM@cpars$M_ageArray <- array(OM_par$Mest, c(nsim, maxage+1, nyears + proyears))
  
  RCModel@OM <- OM
  RCModel@SSB <- OM_par$SSB
  RCModel@NAA <- OM_par$NAA
  RCModel@CAA <- OM_par$CAA
  RCModel@CAL <- OM_par$CAL
  RCModel@conv <- rep(TRUE, nsim)
  RCModel@Misc <- res
  RCModel@config$drop_sim <- integer(0)
  
  return(RCModel)
}
