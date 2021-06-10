
#' Class-\code{project}
#'
#' An S4 class for the output from \link{projection}.
#'
#' @name project-class
#' @docType class
#'
#' @slot Model Name of the assessment model.
#' @slot Name Name of Data object.
#' @slot FMort A matrix of fishing mortality over \code{p_sim} rows and \code{p_years} columns.
#' @slot B An matrix of biomass with \code{p_sim} rows and \code{p_years} columns.
#' @slot SSB A matrix of spawning biomass with \code{p_sim} rows and \code{p_years} columns.
#' @slot VB A matrix of vulnerable biomass with \code{p_sim} rows and \code{p_years} columns.
#' @slot R A matrix of recruitment over \code{p_sim} rows and \code{p_years} columns.
#' @slot N A matrix of abundance over \code{p_sim} rows and \code{p_years} columns.
#' @slot Catch A matrix of simulated observed catch over \code{p_sim} rows and \code{p_years} columns.
#' @slot Index An array of simulated observed index of dimension \code{c(p_sim, p_years, nsurvey)}.
#' @slot C_at_age An array for catch-at-age with dimension \code{c(p_sim, p_years, n_age)}.
#' @seealso \link{projection}
#' @author Q. Huynh
#' @export project
#' @exportClass project
project <- setClass("project", slots = c(Model = "character", Name = "character", FMort = "matrix",
                                         B = "matrix", SSB = "matrix", VB = "matrix", R = "matrix", N = "matrix",
                                         Catch = "matrix", Index = "array", C_at_age = "array"))


#' Projections for assessment models
#'
#' This function takes an assessment model and runs a stochastic projection based on future F or catch.
#'
#' @param Assessment An object of class \linkS4class{Assessment}.
#' @param constrain Whether to project on future F or catch. By default, projects on F.
#' @param Ftarget The projection F, either of length 1 for constant F for the entirety of the projection or length p_years.
#' @param Catch The projection catch, either of length 1 for constant catch for the entirety of the projection or length p_years.
#' @param p_years Integer for the number of projection years.
#' @param p_sim Integer for the number of simulations for the projection.
#' @param obs_error A list of length two. In the first entry, a vector of length nsurvey giving the standard deviations of each future index, 
#' or alternatively an array of dimension p_sim, p_years, and nsurvey giving the deviates. The second entry is
#' the standard deviation of the projected catch. Alternatively, a matrix of 
#' simulation and year-specific error structure for the catch (p_sim rows and p_year columns; a matrix of ones indicates perfect data).
#' @param process_error Numeric, standard deviation for process error (e.g., recruitment or biomass deviates). If \code{NULL},
#' uses values from assessment model. Alternatively, a matrix of simulation and year-specific recruitment deviates (p_sim rows and p_year columns,
#' a matrix of ones indicates no recruitment deviates).
#' @param max_F The maximum allowable F if the projection is constrained on catch.
#' @param seed An integer to set the seed for the sampling observation and process error deviates.
#' @return An object of class \linkS4class{project} that contains future predicted values of F, catch, biomass, recruitment, etc.
#' @examples
#' \donttest{
#' myAssess <- SCA(Data = SimulatedData)
#' do_projection <- projection(myAssess, Ftarget = myAssess@@FMSY)
#' }
#' @export
projection <- function(Assessment, constrain = c("F", "Catch"), Ftarget, Catch, p_years = 50, p_sim = 200,
                       obs_error, process_error, max_F = 3, seed = 499) {
  constrain <- match.arg(constrain)
  if(constrain == "Catch") {
    if(missing(Catch)) stop("Need a value of argument Catch.")
    if(length(Catch) == 1) {
      Catch <- rep(Catch, p_years)
    } else {
      stop(paste0("Catch needs to be of length p_years (", p_years, ")"))
    }
    Ftarget <- NULL 
  }
  if(constrain == "F") {
    if(missing(Ftarget)) stop("Need a value of F (argument \"Ftarget\").")
    if(length(Ftarget) == 1) {
      Ftarget <- rep(Ftarget, p_years)
    } else stop(paste0("Ftarget needs to be of length p_years (", p_years, ")"))
    Catch <- NULL
  }
  if(!missing(obs_error)) {
    if(length(obs_error) < 2) {
      stop("obs_error should be a list of length 2 for the standard deviation of index and catch, respectively.")
    }
  } else {
    nsurvey <- ifelse(is.matrix(Assessment@Index), Assessment@Index %>% ncol(), 1)
    obs_error <- list(array(1, c(p_sim, p_years, nsurvey)), matrix(1, p_sim, p_years))
  }
  if(missing(process_error)) process_error <- matrix(1, p_sim, p_years)
  if(!Assessment@conv) warning("Assessment model did not appear to converge.")

  f <- get(paste0("projection_", Assessment@Model))
  out <- f(Assessment, constrain = constrain, Catch = Catch, Ftarget = Ftarget, p_years = p_years, p_sim = p_sim,
           process_error = process_error, obs_error = obs_error, max_F = max_F, seed = seed)
  return(out)
}


projection_SP <- function(Assessment, constrain = c("F", "Catch"), Ftarget, Catch, p_years = 50, p_sim = 200,
                          obs_error, process_error, max_F = 3, seed = 499) {
  constrain <- match.arg(constrain)
  TMB_report <- Assessment@TMB_report
  TMB_data <- Assessment@obj$env$data
  
  # Sample rec_devs and obs_error
  set.seed(seed)
  if(is.matrix(process_error)) {
    if(nrow(process_error) != p_sim) stop("Number of rows of process_error not equal to ", p_sim, call. = FALSE)
    if(ncol(process_error) != p_years) stop("Number of columns of process_error not equal to ", p_years, call. = FALSE)
    B_dev <- process_error
  } else {
    tau <- ifelse(!is.null(process_error), process_error[1], TMB_report$tau)
    if(!is.null(tau) && tau > 0) B_dev <- exp(matrix(rnorm(p_years * p_sim, 0, tau), p_sim, p_years) - 0.5 * tau^2)
  }
  if(!exists("B_dev", inherits = FALSE)) B_dev <- matrix(1, p_sim, p_years)
  
  nsurvey <- ifelse(is.matrix(Assessment@Index), Assessment@Index %>% ncol(), 1)
  if(is.array(obs_error[[1]])) {
    Iobs_err <- obs_error[[1]]
  } else {
    samps <- rnorm(p_years * p_sim * nsurvey, 0, obs_error[[1]]) %>% array(c(nsurvey, p_sim, p_years)) %>% 
      aperm(c(2, 3, 1))
    Iobs_err <- exp(samps)
  }
  if(is.array(obs_error[[2]])) {
    Cobs_err <- obs_error[[2]]
  } else {
    omega <- ifelse(!is.null(obs_error[[2]]), obs_error[[2]], 0)
    Cobs_err <- exp(matrix(rnorm(p_years * p_sim, 0, omega), p_sim, p_years))
  }
  
  B <- Cpred <- matrix(NA, p_sim, p_years)
  B[, 1] <- TMB_report$B[length(TMB_report$B)] * B_dev[, 1]
  Ipred <- array(NA_real_, c(p_sim, p_years, nsurvey))
  
  if(constrain == "F") {
    Fout <- matrix(pmin(Ftarget, ifelse(TMB_data$dt < 1, max_F, 1 - exp(-max_F))), p_sim, p_years, byrow = TRUE)
  } else {
    Fout <- matrix(NA, p_sim, p_years)
  }
  
  for(y in 1:p_years + 1) {
    if(constrain == "F") {
      one_ts <- lapply(B[, y-1], function(x, ...) SP_catch_solver(B = x, ...), FM = Ftarget[y-1], dt = TMB_data$dt,
                       MSY = TMB_report$MSY, K = TMB_report$K, n = TMB_report$n, n_term = TMB_report$n_term)
    } else {
      Fout[, y-1] <- vapply(B[, y-1], function(x, ...) optimize(SP_catch_solver, c(1e-8, ifelse(TMB_data$dt < 1, max_F, 1 - exp(-max_F))),
                                                                B = x, ...)$minimum, numeric(1),
                            dt = TMB_data$dt, MSY = TMB_report$MSY, K = TMB_report$K, n = TMB_report$n, n_term = TMB_report$n_term,
                            TAC = Catch[y-1])
      one_ts <- Map(SP_catch_solver, FM = Fout[, y-1], B = B[, y-1], 
                    MoreArgs = list(dt = TMB_data$dt, MSY = TMB_report$MSY, K = TMB_report$K, n = TMB_report$n, n_term = TMB_report$n_term))
    }
    Cpred[, y-1] <- vapply(one_ts, getElement, numeric(1), 1)
    Ipred[, y-1, ] <- outer(B[, y-1], TMB_report$q) * Iobs_err[, y-1, ]
    if(y <= p_years) B[, y] <- vapply(one_ts, getElement, numeric(1), 2)
  }
  
  return(new("project", Catch = Cpred * Cobs_err, Index = Ipred, B = B, VB = B, SSB = B, FMort = Fout))
}

SP_catch_solver <- function(FM, B, dt, MSY, K, n, n_term, TAC = NULL) {
  B_int <- numeric(1/dt + 1)
  B_int[1] <- B
  Cpred_int <- numeric(1/dt)

  for(j in 1:length(Cpred_int)) {
    Cpred_int[j] <- dt * FM * B_int[j]
    B_int[j+1] <- B_int[j] - Cpred_int[j] + dt *
      ifelse(n == 1, -exp(1) * MSY * B_int[j]/K * log(B_int[j]/K),
             n_term/(n-1) * MSY * (B_int[j]/K - (B_int[j]/K)^n))
  }
  if(!is.null(TAC)) {
    return((sum(Cpred_int) - TAC)^2)
  } else {
    return(c(sum(Cpred_int), B_int[length(B_int)]))
  }
}

projection_SP_SS <- projection_SP

projection_SCA <- function(Assessment, constrain = c("F", "Catch"), Ftarget, Catch,
                           p_years = 50, p_sim = 200, obs_error, process_error,
                           max_F = 3, seed = 499, ...) {
  
  constrain <- match.arg(constrain)
    
  TMB_report <- Assessment@TMB_report
  TMB_data <- Assessment@obj$env$data
  Pope <- TMB_data$catch_eq == "Pope"
  
  # Sample rec_devs and obs_error
  set.seed(seed)
  if(is.matrix(process_error)) {
    if(nrow(process_error) != p_sim) stop("Number of rows of process_error not equal to ", p_sim, call. = FALSE)
    if(ncol(process_error) != p_years) stop("Number of columns of process_error not equal to ", p_years, call. = FALSE)
    p_log_rec_dev <- process_error
  } else {
    tau <- ifelse(!is.null(process_error), process_error[1], TMB_report$tau)
    if(!is.null(tau) && tau > 0) p_log_rec_dev <- exp(matrix(rnorm(p_years * p_sim, 0, tau), p_sim, p_years) - 0.5 * tau^2)
  }
  
  if(!exists("p_log_rec_dev", inherits = FALSE)) p_log_rec_dev <- matrix(1, p_sim, p_years)
  
  nsurvey <- ifelse(is.matrix(Assessment@Index), Assessment@Index %>% ncol(), 1)
  if(is.array(obs_error[[1]])) {
    Iobs_err <- obs_error[[1]]
  } else {
    samps <- rnorm(p_years * p_sim * nsurvey, 0, obs_error[[1]]) %>% array(c(nsurvey, p_sim, p_years)) %>% 
      aperm(c(2, 3, 1))
    Iobs_err <- exp(samps)
  }
  if(is.array(obs_error[[2]])) {
    Cobs_err <- obs_error[[2]]
  } else {
    omega <- ifelse(!is.null(obs_error[[2]]), obs_error[[2]], 0)
    Cobs_err <- exp(matrix(rnorm(p_years * p_sim, 0, omega), p_sim, p_years))
  }

  p_output <- lapply(1:p_sim, function(x, ...) projection_SCA_internal(p_log_rec_dev = p_log_rec_dev[x, ], Cobs_err = Cobs_err[x, ],
                                                                       Iobs_err = Iobs_err[x, , ], ...),
                     FMort = Ftarget, Catch = Catch, constrain = constrain, TMB_report = TMB_report, TMB_data = TMB_data,
                     Pope = Pope, max_F = max_F)

  CAApred <- lapply(p_output, getElement, "CAApred") %>% simplify2array()
  Ipred <- lapply(p_output, getElement, "Ipred") %>% simplify2array()
  output <- new("project", Catch = do.call(rbind, lapply(p_output, getElement, "Cpred")),
                C_at_age = aperm(CAApred, c(3, 1, 2)), Index = aperm(Ipred, c(3, 1, 2)),
                SSB = do.call(rbind, lapply(p_output, getElement, "E")),
                R = do.call(rbind, lapply(p_output, function(x) x$N[, 1])),
                N = do.call(rbind, lapply(p_output, function(x) rowSums(x$N))), VB = do.call(rbind, lapply(p_output, getElement, "VB")),
                B = do.call(rbind, lapply(p_output, getElement, "B")), FMort = do.call(rbind, lapply(p_output, getElement, "Fout")))
  return(output)
}

projection_SCA_internal <- function(FMort, Catch, constrain, TMB_report, TMB_data, p_log_rec_dev, Cobs_err,
                                    Iobs_err, Pope = FALSE, max_F = 3) {

  weight <- TMB_data$weight
  mat <- TMB_data$mat
  vul <- TMB_report$vul
  
  tv_M <- TMB_data$tv_M
  M_bounds <- TMB_data$M_bounds

  if(constrain == "F") {
    if(Pope) {
      Fout <- pmin(FMort, 1 - exp(-max_F))
      UU <- vul %o% FMort
    } else {
      Fout <- pmin(FMort, max_F)
      FF <- vul %o% FMort
    }
  } else {
    Fout <- numeric(length(p_log_rec_dev))
    if(Pope) {
      UU <- matrix(NA, length(vul), length(p_log_rec_dev))
    } else {
      FF <- matrix(NA, length(vul), length(p_log_rec_dev))
    }
  }
  surv <- M_p <- matrix(NA, length(vul), length(p_log_rec_dev))

  N <- CAApred <- matrix(NA, length(p_log_rec_dev), ncol(TMB_report$N))
  N[1, ] <- TMB_report$N[nrow(TMB_report$N), ]
  N[1, 1] <- N[1, 1] * p_log_rec_dev[1]
  M_p[, 1] <- TMB_report$M[nrow(TMB_report$M), ]

  E <- VB <- B <- Cpred <- numeric(length(p_log_rec_dev))
  E[1] <- sum(N[1, ] * mat * weight)
  B[1] <- sum(N[1, ] * weight)
  for(y in 1:length(p_log_rec_dev)) {
    if(constrain == "Catch") {
      if(Pope) {
        VB[y] <- sum(N[y, ] * exp(-0.5 * M_p[, y]) * vul * weight)
        Fout[y] <- ifelse(Catch[y]/VB[y] > 1 - exp(-max_F), 1 - exp(-max_F), Catch[y]/VB[y])
        UU[, y] <- vul * Fout[y]
      } else {
        Fout[y] <- optimize(SCA_catch_solver, c(1e-8, max_F), N = N[y, ], weight = weight, vul = vul, M = M_p[, y], TAC = Catch[y])$minimum
        FF[, y] <- vul * Fout[y]
      }
    }
    
    if(Pope) {
      surv[, y] <- (1 - UU[, y]) * exp(-M_p[, y])
    } else {
      surv[, y] <- exp(-FF[, y] - M_p[, y])
    }

    if(y < length(p_log_rec_dev)) {
      N[y+1, 2:ncol(N)] <- N[y, 1:(ncol(N)-1)] * surv[1:(ncol(N)-1), y]
      N[y+1, ncol(N)] <- N[y+1, ncol(N)] + N[y, ncol(N)] * surv[ncol(N), y]
      
      E[y+1] <- sum(N[y+1, ] * mat * weight, na.rm = TRUE)
      B[y+1] <- sum(N[y+1, ] * weight, na.rm = TRUE)
      
      if(tv_M == "DD") {
        M_p[, y+1] <- ifelse(B[y+1] <= TMB_report$B0, M_bounds[1] + (M_bounds[2] - M_bounds[1]) * (1 - B[y+1]/TMB_report$B0), M_bounds[1])
      } else {
        M_p[, y+1] <- M_p[, 1]
      }
      N[y+1, 1] <- R_pred(E[y+1], TMB_report$h, TMB_report$R0, TMB_report$E0, TMB_data$SR_type) * p_log_rec_dev[y+1]
    }
  }
  if(Pope) {
    CAApred <- t(t(N) * exp(-0.5 * M_p) * UU)
    if(constrain == "F") VB <- colSums(t(N) * TMB_report$vul * weight * exp(-0.5 * M_p))
    Cpred <- colSums(t(CAApred) * weight)
  } else {
    out <- lapply(1:length(Fout), function(x) SCA_catch_solver(FM = Fout[x], N = N[x, ], weight = weight, vul = vul, M = M_p[, x]))
    CAApred <- do.call(rbind, lapply(out, getElement, 1))
    VB <- colSums(t(N) * TMB_report$vul * weight)
    Cpred <- vapply(out, getElement, numeric(1), 2)
  }
  
  if(!is.matrix(Iobs_err)) Iobs_err <- matrix(Iobs_err, ncol = 1)
  Ipred <- vapply(1:TMB_data$nsurvey, function(sur) {
    if(sum(TMB_data$I_vul[sur])) {
      I_vul <- TMB_data$I_vul[, sur]
    } else {
      I_vul <- vul
    }
    N_sur <- I_vul * t(N)
    if(TMB_data$I_units[sur]) {
      survey <- colSums(N_sur * weight)
    } else {
      survey <- colSums(N_sur)
    }
    Ipred <- TMB_report$q[sur] * survey * Iobs_err[, sur]
    return(Ipred)
  }, numeric(length(p_log_rec_dev)))
  
  return(list(Cpred = Cpred * Cobs_err, CAApred = CAApred, Ipred = Ipred, B = B, VB = VB, E = E, N = N, Fout = Fout))
}

Baranov <- function(sel = 1, apicalF, M, N) {
  FF <- sel * apicalF
  Z <- FF + M
  CAA <- FF / Z * (1 - exp(-Z)) * N
  return(CAA)
}

SCA_catch_solver <- function(FM, N, weight = 1, vul = 1, M, TAC = NULL) {
  CAA <- Baranov(vul, FM, M, N)
  Cw <- sum(CAA * weight)
  if(!is.null(TAC)) {
    return((Cw - TAC)^2)
  } else {
    return(list(CAA = CAA, Cpred = Cw, FM = FM))
  }
}


R_pred <- function(SSB, h, R0, SSB0, SR_type = c("BH", "Ricker", "none")) {
  SR_type <- match.arg(SR_type)
  if(SR_type == "BH") {
    den <- SSB0 * (1 - h) + (5*h - 1) * SSB
    RR <- 4 * h * R0 * SSB / den
  } else if(SR_type == "Ricker") {
    expon <- 1.25 * (1 - SSB/SSB0)
    RR <- (5*h)^expon * SSB * R0 / SSB0
  } else {
    RR <- R0
  }
  return(RR)
}

projection_cDD <- projection_cDD_SS <- function(Assessment, constrain = c("F", "Catch"), Ftarget, Catch,
                                                p_years = 50, p_sim = 200, obs_error, process_error, max_F = 3,
                                                seed = 499, ...) {
  constrain <- match.arg(constrain)
  
  TMB_report <- Assessment@TMB_report
  TMB_data <- Assessment@obj$env$data
  
  # Sample rec_devs and obs_error
  set.seed(seed)
  if(is.matrix(process_error)) {
    if(nrow(process_error) != p_sim) stop("Number of rows of process_error not equal to ", p_sim, call. = FALSE)
    if(ncol(process_error) != p_years) stop("Number of columns of process_error not equal to ", p_years, call. = FALSE)
    p_log_rec_dev <- process_error
  } else {
    tau <- ifelse(!is.null(process_error), process_error[1], TMB_report$tau)
    if(!is.null(tau) && tau > 0) p_log_rec_dev <- exp(matrix(rnorm(p_years * p_sim, 0, tau), p_sim, p_years) - 0.5 * tau^2)
  }
  
  if(!exists("p_log_rec_dev", inherits = FALSE)) p_log_rec_dev <- matrix(1, p_sim, p_years)
  
  nsurvey <- ifelse(is.matrix(Assessment@Index), Assessment@Index %>% ncol(), 1)
  if(is.array(obs_error[[1]])) {
    Iobs_err <- obs_error[[1]]
  } else {
    samps <- rnorm(p_years * p_sim * nsurvey, 0, obs_error[[1]]) %>% array(c(nsurvey, p_sim, p_years)) %>% 
      aperm(c(2, 3, 1))
    Iobs_err <- exp(samps)
  }
  if(is.array(obs_error[[2]])) {
    Cobs_err <- obs_error[[2]]
  } else {
    omega <- ifelse(!is.null(obs_error[[2]]), obs_error[[2]], 0)
    Cobs_err <- exp(matrix(rnorm(p_years * p_sim, 0, omega), p_sim, p_years))
  }
  
  B <- N <- R <- Cpred <- matrix(NA_real_, p_sim, p_years)
  Ipred <- array(NA_real_, c(p_sim, p_years, nsurvey))
  B[, 1] <- TMB_report$B[length(TMB_report$B)]
  N[, 1] <- TMB_report$N[length(TMB_report$N)]
  R[, 1:TMB_data$k] <- rep(TMB_report$R[(length(TMB_report$R) - TMB_data$k + 1):length(TMB_report$R)], each = p_sim) *
    p_log_rec_dev[, 1:TMB_data$k]

  if(constrain == "F") {
    Fout <- matrix(pmin(Ftarget, max_F), p_sim, p_years, byrow = TRUE)
  } else {
    Fout <- matrix(NA, p_sim, p_years)
  }

  for(y in 1:p_years + 1) {
    if(y+TMB_data$k-1 <= p_years) {
      R[, y+TMB_data$k-1] <- R_pred(B[, y-1], TMB_report$h, TMB_report$R0, TMB_report$B0, TMB_data$SR_type) *
        p_log_rec_dev[, y+TMB_data$k-1]
    }

    if(constrain == "F") {
      one_ts <- Map(function(x, y, z, ...) cDD_catch_solver(B = x, N = y, R = z, ...), x = B[, y-1], y = N[, y-1], z = R[, y-1],
                    MoreArgs = list(FM = Ftarget[y-1], Kappa = TMB_data$Kappa, Winf = TMB_data$Winf, wk = TMB_data$wk, M = TMB_report$M))
    } else {
      Fout[, y-1] <- mapply(function(x, y, z, ...) optimize(cDD_catch_solver, c(1e-8, max_F), B = x, N = y, R = z, ...)$minimum,
                            x = B[, y-1], y = N[, y-1], z = R[, y-1],
                            MoreArgs = list(Kappa = TMB_data$Kappa, Winf = TMB_data$Winf, wk = TMB_data$wk, M = TMB_report$M, TAC = Catch[y-1]))

      one_ts <- Map(function(x, y, z, w, ...) cDD_catch_solver(B = x, N = y, R = z, FM = w, ...), x = B[, y-1], y = N[, y-1], z = R[, y-1],
                    w = Fout[, y-1],
                    MoreArgs = list(Kappa = TMB_data$Kappa, Winf = TMB_data$Winf, wk = TMB_data$wk, M = TMB_report$M))
    }

    Cpred[, y-1] <- vapply(one_ts, getElement, numeric(1), 1)
    Ipred[, y-1, ] <- outer(B[, y-1], TMB_report$q) * Iobs_err[, y-1, ]

    if(y <= p_years) {
      N[, y] <- vapply(one_ts, getElement, numeric(1), 2)
      B[, y] <- vapply(one_ts, getElement, numeric(1), 3)
    }
  }

  new("project", Catch = Cpred * Cobs_err, Index = Ipred, B = B, VB = B, SSB = B, R = R, N = N, FMort = Fout)
}

cDD_catch_solver <- function(FM, B, N, R, M, Kappa, Winf, wk, TAC = NULL) {
  Z <- FM + M
  surv <- exp(-Z)
  ZK <- Z + Kappa
  BPR <- (Kappa * Winf/Z + wk)/ZK
  Binf <- BPR * R
  Ninf <- R/Z

  C1 <- B - Binf - (N - Ninf) * Kappa * Winf / ZK
  C2 <- 1 - exp(-ZK)
  C4 <- Binf + (N - Ninf) * Kappa * Winf / ZK
  Cpred <- FM * ((C1 * C2)/ZK + C4)

  if(!is.null(TAC)) {
    return((Cpred - TAC)^2)
  } else {

    N_next <- Ninf + (N - Ninf) * surv

    B2 <- Kappa * Winf * (N - Ninf) / ZK
    B3 <- (B - Binf - Kappa * Winf * (N - Ninf) / ZK) * exp(-ZK)
    B_next <- Binf + B2 + B3

    return(c(Cpred, N_next, B_next))
  }
}

projection_DD_TMB <- projection_DD_SS <- function(Assessment, constrain = c("F", "Catch"), Ftarget, Catch,
                                                  p_years = 50, p_sim = 200, obs_error, process_error,
                                                  max_F = 3, seed = 499, ...) {
  constrain <- match.arg(constrain)
  
  TMB_report <- Assessment@TMB_report
  TMB_data <- Assessment@obj$env$data
  
  # Sample rec_devs and obs_error
  set.seed(seed)
  if(is.matrix(process_error)) {
    if(nrow(process_error) != p_sim) stop("Number of rows of process_error not equal to ", p_sim, call. = FALSE)
    if(ncol(process_error) != p_years) stop("Number of columns of process_error not equal to ", p_years, call. = FALSE)
    p_log_rec_dev <- process_error
  } else {
    tau <- ifelse(!is.null(process_error), process_error[1], TMB_report$tau)
    if(!is.null(tau) && tau > 0) p_log_rec_dev <- exp(matrix(rnorm(p_years * p_sim, 0, tau), p_sim, p_years) - 0.5 * tau^2)
  }
  
  if(!exists("p_log_rec_dev", inherits = FALSE)) p_log_rec_dev <- matrix(1, p_sim, p_years)
  
  nsurvey <- ifelse(is.matrix(Assessment@Index), Assessment@Index %>% ncol(), 1)
  if(is.array(obs_error[[1]])) {
    Iobs_err <- obs_error[[1]]
  } else {
    samps <- rnorm(p_years * p_sim * nsurvey, 0, obs_error[[1]]) %>% array(c(nsurvey, p_sim, p_years)) %>% 
      aperm(c(2, 3, 1))
    Iobs_err <- exp(samps)
  }
  if(is.array(obs_error[[2]])) {
    Cobs_err <- obs_error[[2]]
  } else {
    omega <- ifelse(!is.null(obs_error[[2]]), obs_error[[2]], 0)
    Cobs_err <- exp(matrix(rnorm(p_years * p_sim, 0, omega), p_sim, p_years))
  }
  
  B <- N <- R <- Cpred <- matrix(NA_real_, p_sim, p_years)
  Ipred <- array(NA_real_, c(p_sim, p_years, nsurvey))
  B[, 1] <- TMB_report$B[length(TMB_report$B)]
  N[, 1] <- TMB_report$N[length(TMB_report$N)]
  R[, 1:TMB_data$k] <- rep(TMB_report$R[(length(TMB_report$R) - TMB_data$k + 1):length(TMB_report$R)], each = p_sim) *
    p_log_rec_dev[, 1:TMB_data$k]

  if(constrain == "F") {
    Fout <- matrix(pmin(Ftarget, max_F), p_sim, p_years, byrow = TRUE)
    surv <- exp(-TMB_report$M - Fout)
  } else {
    Fout <- surv <- matrix(NA_real_, p_sim, p_years)
  }

  for(y in 1:p_years + 1) {
    if(y+TMB_data$k-1 <= p_years) {
      R[, y+TMB_data$k-1] <- R_pred(B[, y-1], TMB_report$h, TMB_report$R0, TMB_report$B0, TMB_data$SR_type) *
        p_log_rec_dev[, y+TMB_data$k-1]
    }
    if(constrain == "Catch") {
      Fout[,y-1] <- vapply(1:p_sim, function(i) {
        optimize(SCA_catch_solver, c(1e-8, max_F), N = B[i, y-1], M = TMB_report$M, TAC = Catch[y-1])$minimum
      }, numeric(1))
      surv[, y-1] <- exp(-TMB_report$M - Fout[,y-1])
    }
    Cpred[, y-1] <- vapply(1:p_sim, function(i) Baranov(apicalF = Fout[i, y-1], M = TMB_report$M, N = B[i, y-1]), numeric(1))
    if(Assessment@obj$env$data$condition == "catch") {
      Ipred[, y-1, ] <- outer(B[, y-1], TMB_report$q) * Iobs_err[, y-1, ]
    } 
    if(y <= p_years) {
      B[, y] <- surv[, y-1] * (TMB_data$Alpha * N[, y-1] + TMB_data$Rho * B[, y-1]) + TMB_data$wk * R[, y]
      N[, y] <- surv[, y-1] * N[, y-1] + R[, y]
    }
  }
  if(Assessment@obj$env$data$condition == "effort") {
    Eff <- Fout/TMB_report$q/Assessment@info$E_rescale
    Ipred <- Cpred/Eff * Iobs_err
  }

  new("project", Catch = Cpred * Cobs_err, Index = Ipred, B = B, VB = B, SSB = B, R = R, N = N, FMort = Fout)
}

projection_VPA <- function(Assessment, constrain = c("F", "Catch"), Ftarget, Catch,
                           p_years = 50, p_sim = 200, obs_error, process_error,
                           max_F = 3, seed = 499, ...) {
  
  Assessment@TMB_report$meanR <- mean(Assessment@R)
  Assessment@TMB_report$vul <- Assessment@TMB_report$vul_p
  Assessment@obj$env$data$mat <- Assessment@info$LH$mat
  
  projection_SCA(Assessment, constrain, Ftarget, Catch, p_years, p_sim, obs_error, process_error,
                 max_F, seed, ...)
}


