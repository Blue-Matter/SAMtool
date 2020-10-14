
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
#' @slot Catch A matrix of observed catch over \code{p_sim} rows and \code{p_years} columns.
#' @slot Index A matrix of observed index over \code{p_sim} rows and \code{p_years} columns.
#' @slot C_at_age An array for catch-at-age with dimension \code{c(p_sim, p_years, maxage)}.
#' @seealso \link{projection}
#' @author Q. Huynh
#' @export project
#' @exportClass project
project <- setClass("project", slots = c(Model = "character", Name = "character", FMort = "matrix",
                                         B = "matrix", SSB = "matrix", VB = "matrix", R = "matrix", N = "matrix",
                                         Catch = "matrix", Index = "matrix", C_at_age = "array"))


#' Projections for assessment models
#'
#' This function takes an assessment model and runs a stochastic projection based on future F or catch.
#'
#' @param Assessment An object of class \linkS4class{Assessment}.
#' @param constrain Whether to project on future F or catch. By default, projects on F.
#' @param FMort The projection F, either of length 1 for constant F for the entirety of the projection or length p_years.
#' @param Catch The projection catch, either of length 1 for constant catch for the entirety of the projection or length p_years.
#' @param p_years Integer for the number of projection years.
#' @param p_sim Integer for the number of simulations for the projection.
#' @param obs_error Vector of length two for standard deviation of error to be added to the index and catch, respectively. If \code{NULL},
#' uses values from assessment model.
#' @param process_error Numeric, standard deviation for process error (e.g., recruitment or biomass deviates). If \code{NULL},
#' uses values from assessment model.
#' @param max_F The maximum allowable F if the projection is constrained on catch.
#' @param seed An integer to set the seed for the sampling observation and process error deviates.
#' @examples
#' \donttest{
#' myAssess <- SCA(Data = SimulatedData)
#' do_projection <- projection(myAssess, FMort = myAssess@@FMSY)
#' }
#' @export
projection <- function(Assessment, constrain = c("F", "Catch"), FMort = NULL, Catch = NULL, p_years = 50, p_sim = 200,
                       obs_error = NULL, process_error = NULL, max_F = 3, seed = 499) {
  constrain <- match.arg(constrain)
  if(constrain == "Catch") {
    if(is.null("Catch")) stop("Need a value of argument Catch.")
    if(length(Catch) == 1) {
      Catch <- rep(Catch, p_years)
    } else {
      stop(paste0("Catch needs to be of length p_years (", p_years, ")"))
    }
  }
  if(constrain == "F") {
    if(is.null(FMort)) stop("Need a value of F (argument \"FMort\").")
    if(length(FMort) == 1) {
      FMort <- rep(FMort, p_years)
    } else stop(paste0("FMort needs to be of length p_years (", p_years, ")"))
  }
  if(!is.null(obs_error) && length(obs_error) < 2) stop("obs_error should be a vector of length 2 for the standard deviation of index and catch, respectively.")

  if(!Assessment@conv) warning("Assessment model did not appear to converge.")

  f <- get(paste0("projection_", Assessment@Model))
  out <- f(Assessment, constrain = constrain, Catch = Catch, FMort = FMort, p_years = p_years, p_sim = p_sim,
           process_error = process_error, obs_error = obs_error, max_F = max_F, seed = seed)
  return(out)
}


projection_SP <- function(Assessment, constrain = c("F", "Catch"), FMort = NULL, Catch = NULL, p_years = 50, p_sim = 200,
                          obs_error = NULL, process_error = NULL, max_F = 3, seed = 499) {
  TMB_report <- Assessment@TMB_report
  TMB_data <- Assessment@obj$env$data

  # Sample rec_devs and obs_error
  set.seed(seed)
  tau <- ifelse(!is.null(process_error), process_error[1], TMB_report$tau)

  B_dev <- matrix(1, p_sim, p_years)
  if(!is.null(tau) && tau > 0) B_dev <- exp(matrix(rnorm(p_years * p_sim, 0, tau), p_sim, p_years) - 0.5 * tau^2)

  sigma <- ifelse(!is.null(obs_error), obs_error[1], TMB_report$sigma)
  omega <- ifelse(!is.null(obs_error), obs_error[2], 0)
  Iobs_err <- exp(matrix(rnorm(p_years * p_sim, 0, sigma), p_sim, p_years))
  Cobs_err <- exp(matrix(rnorm(p_years * p_sim, 0, omega), p_sim, p_years))

  B <- Cpred <- matrix(NA, p_sim, p_years)
  B[, 1] <- TMB_report$B[length(TMB_report$B)] * B_dev[, 1]

  if(constrain == "F") {
    Fout <- matrix(pmin(FMort, ifelse(TMB_data$dt < 1, max_F, 1 - exp(-max_F))), p_sim, p_years, byrow = TRUE)
  } else {
    Fout <- matrix(NA, p_sim, p_years)
  }

  for(y in 2:(p_years+1)) {
    if(constrain == "F") {
      one_ts <- lapply(B[, y-1], function(x, ...) SP_catch_solver(B = x, ...), FM = FMort[y-1], dt = TMB_data$dt,
                       MSY = TMB_report$MSY, K = TMB_report$K, n = TMB_report$n, n_term = TMB_report$n_term)
    } else {
      Fout[, y-1] <- vapply(B[, y-1], function(x, ...) optimize(SP_catch_solver, c(1e-8, ifelse(TMB_data$dt < 1, max_F, 1 - exp(-max_F))),
                                                                B = x, ...)$minimum, numeric(1),
                            dt = TMB_data$dt, MSY = TMB_report$MSY, K = TMB_report$K, n = TMB_report$n, n_term = TMB_report$n_term,
                            TAC = Catch[y-1])
      one_ts <- Map(SP_catch_solver, FM = Fout[, y-1], B = B[, y-1], MoreArgs = list(dt = TMB_data$dt,
                    MSY = TMB_report$MSY, K = TMB_report$K, n = TMB_report$n, n_term = TMB_report$n_term))
    }
    Cpred[, y-1] <- vapply(one_ts, getElement, numeric(1), 1)
    if(y <= p_years) B[, y] <- vapply(one_ts, getElement, numeric(1), 2)
  }
  Ipred <- TMB_report$q * B * Iobs_err

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

projection_SCA <- projection_SCA_Pope <- function(Assessment, constrain = c("F", "Catch"), FMort = NULL, Catch = NULL,
                                                  p_years = 50, p_sim = 200, obs_error = NULL, process_error = NULL,
                                                  max_F = 3, seed = 499, ...) {
  TMB_report <- Assessment@TMB_report
  TMB_data <- Assessment@obj$env$data
  Pope <- Assessment@Model == "SCA_Pope"

  # Sample rec_devs and obs_error
  set.seed(seed)
  tau <- ifelse(!is.null(process_error), process_error[1], TMB_report$tau)

  p_log_rec_dev <- matrix(1, p_sim, p_years)
  if(tau > 0) p_log_rec_dev <- exp(matrix(rnorm(p_years * p_sim, 0, tau), p_sim, p_years) - 0.5 * tau^2)

  sigma <- ifelse(!is.null(obs_error), obs_error[1], TMB_report$sigma)
  omega <- ifelse(!is.null(obs_error), obs_error[2], TMB_report$omega)

  Iobs_err <- exp(matrix(rnorm(p_years * p_sim, 0, sigma), p_sim, p_years))
  Cobs_err <- exp(matrix(rnorm(p_years * p_sim, 0, omega), p_sim, p_years))

  p_output <- lapply(1:p_sim, function(x, ...) projection_SCA_internal(p_log_rec_dev = p_log_rec_dev[x, ], Cobs_err = Cobs_err[x, ],
                                                                       Iobs_err = Iobs_err[x, ], ...),
                     FMort = FMort, Catch = Catch, constrain = constrain, TMB_report = TMB_report, TMB_data = TMB_data,
                     Pope = Pope, max_F = max_F)

  CAApred <- array(unlist(lapply(p_output, getElement, "CAApred")), dim = c(50, 15, p_sim))

  output <- new("project", Catch = do.call(rbind, lapply(p_output, getElement, "Cpred")),
                C_at_age = aperm(CAApred, c(3, 1, 2)),
                Index = do.call(rbind, lapply(p_output, getElement, "Ipred")),
                SSB = do.call(rbind, lapply(p_output, getElement, "E")),
                N = do.call(rbind, lapply(p_output, function(x) rowSums(x$N))), VB = do.call(rbind, lapply(p_output, getElement, "VB")),
                B = do.call(rbind, lapply(p_output, getElement, "B")), FMort = do.call(rbind, lapply(p_output, getElement, "Fout")))
  return(output)
}

projection_SCA2 <- function(...) stop("No projection function available yet for SCA2", .call = FALSE)


projection_SCA_internal <- function(FMort = NULL, Catch = NULL, constrain, TMB_report, TMB_data, p_log_rec_dev, Cobs_err,
                                    Iobs_err, Pope = FALSE, max_F = 3) {

  weight <- TMB_data$weight
  mat <- TMB_data$mat
  vul <- TMB_report$vul

  if(constrain == "F") {

    if(Pope) {
      Fout <- pmin(FMort, max_F)
      UU <- vul %o% FMort
      surv <- (1 - UU) * exp(-TMB_data$M)
    } else {
      Fout <- pmin(FMort, 1 - exp(-max_F))
      FF <- vul %o% FMort
      surv <- exp(-FF - TMB_data$M)
    }
  } else {
    Fout <- numeric(length(p_log_rec_dev))
    if(Pope) {
      surv <- UU <- matrix(NA, length(vul), length(p_log_rec_dev))
    } else {
      surv <- FF <- matrix(NA, length(vul), length(p_log_rec_dev))
    }
  }

  N <- CAApred <- matrix(NA, length(p_log_rec_dev), ncol(TMB_report$N))
  N[1, ] <- TMB_report$N[nrow(TMB_report$N), ]
  N[1, 1] <- N[1, 1] * p_log_rec_dev[1]

  E <- VB <- Cpred <- numeric(length(p_log_rec_dev))
  for(y in 1:length(p_log_rec_dev)) {
    E[y] <- sum(N[y, ] * mat * weight)

    if(constrain == "Catch") {
      if(Pope) {
        VB[y] <- sum(N[y, ] * exp(-0.5 * TMB_data$M) * vul * weight)
        Fout[y] <- ifelse(Catch[y]/VB[y] > 1 - exp(-max_F), 1 - exp(-max_F), Catch[y]/VB[y])
        UU[, y] <- vul * Fout[y]
        surv[, y] <- (1 - UU[, y]) * exp(-TMB_data$M)
      } else {
        Fout[y] <- optimize(SCA_catch_solver, c(1e-8, max_F), N = N[y, ], weight = weight, vul = vul, M = TMB_data$M, TAC = Catch[y])$minimum
        FF[, y] <- vul * Fout[y]
        surv[, y] <- exp(-FF[, y] - TMB_data$M)
      }
    }

    if(y < length(p_log_rec_dev)) {
      N[y+1, 1] <- R_pred(E[y], TMB_report$h, TMB_report$R0, TMB_report$E0, TMB_data$SR_type) * p_log_rec_dev[y+1]
      N[y+1, 2:ncol(N)] <- N[y, 1:(ncol(N)-1)] * surv[1:(ncol(N)-1), y]
      N[y+1, ncol(N)] <- N[y+1, ncol(N)] + N[y, ncol(N)] * surv[ncol(N), y]
    }
  }
  if(Pope) {
    CAApred <- t(t(N) * exp(-0.5 * TMB_data$M) * UU)
    if(constrain == "F") VB <- colSums(t(N) * TMB_report$vul * weight * exp(-0.5 * TMB_data$M))
    Cpred <- colSums(t(CAApred) * weight)
  } else {
    out <- lapply(1:length(Fout), function(x) SCA_catch_solver(FM = Fout[x], N = N[x, ], weight = weight, vul = vul, M = TMB_data$M))
    CAApred <- do.call(rbind, lapply(out, getElement, 1))
    VB <- colSums(t(N) * TMB_report$vul * weight)
    Cpred <- vapply(out, getElement, numeric(1), 2)
  }
  B <- colSums(t(N) * weight)

  if(TMB_data$I_type == "B") {
    Ipred <- TMB_report$q * B
  } else if(TMB_data$I_type == "VB") {
    Ipred <- TMB_report$q * VB
  } else {
    Ipred <- TMB_report$q * E
  }

  return(list(Cpred = Cpred * Cobs_err, CAApred = CAApred, Ipred = Ipred * Iobs_err, B = B, VB = VB, E = E, N = N, Fout = Fout))
}

SCA_catch_solver <- function(FM, N, weight, vul, M, TAC = NULL) {
  mean_N <- N * (1 - exp(-vul * FM - M)) / (vul * FM + M)
  CAA <- FM * mean_N
  Cw <- sum(CAA * weight)
  if(!is.null(TAC)) {
    return((Cw - TAC)^2)
  } else {
    return(list(CAA = CAA, Cpred = Cw, FM = FM))
  }
}


R_pred <- function(SSB, h, R0, SSB0, SR_type = c("BH", "Ricker")) {
  SR_type <- match.arg(SR_type)
  if(SR_type == "BH") {
    den <- SSB0 * (1 - h) + (5*h - 1) * SSB
    RR <- 4 * h * R0 * SSB / den
  } else {
    expon <- 1.25 * (1 - SSB/SSB0)
    RR <- (5*h)^expon * SSB * R0 / SSB0
  }
  return(RR)
}

projection_cDD <- projection_cDD_SS <- function(Assessment, constrain = c("F", "Catch"), FMort = NULL, Catch = NULL,
                                                p_years = 50, p_sim = 200, obs_error = NULL, process_error = NULL, max_F = 3,
                                                seed = 499, ...) {
  TMB_report <- Assessment@TMB_report
  TMB_data <- Assessment@obj$env$data

  # Sample rec_devs and obs_error
  set.seed(seed)
  tau <- ifelse(!is.null(process_error), process_error[1], TMB_report$tau)
  p_log_rec_dev <- matrix(1, p_sim, p_years)
  if(!is.null(tau) && tau > 0) p_log_rec_dev <- exp(matrix(rnorm(p_years * p_sim, 0, tau), p_sim, p_years) - 0.5 * tau^2)

  sigma <- ifelse(!is.null(obs_error), obs_error[1], TMB_report$sigma)
  omega <- ifelse(!is.null(obs_error), obs_error[2], 0)
  Iobs_err <- exp(matrix(rnorm(p_years * p_sim, 0, sigma), p_sim, p_years))
  Cobs_err <- exp(matrix(rnorm(p_years * p_sim, 0, omega), p_sim, p_years))

  B <- N <- R <- Cpred <- matrix(NA, p_sim, p_years)
  B[, 1] <- TMB_report$B[length(TMB_report$B)]
  N[, 1] <- TMB_report$N[length(TMB_report$N)]
  R[, 1:TMB_data$k] <- rep(TMB_report$R[(length(TMB_report$R) - TMB_data$k + 1):length(TMB_report$R)], each = p_sim) *
    p_log_rec_dev[, 1:TMB_data$k]

  if(constrain == "F") {
    Fout <- matrix(pmin(FMort, max_F), p_sim, p_years, byrow = TRUE)
  } else {
    Fout <- matrix(NA, p_sim, p_years)
  }

  for(y in 2:(p_years+1)) {
    if(y+TMB_data$k-1 <= p_years) {
      R[, y+TMB_data$k-1] <- R_pred(B[, y-1], TMB_report$h, TMB_report$R0, TMB_report$B0, TMB_data$SR_type) *
        p_log_rec_dev[, y+TMB_data$k-1]
    }

    if(constrain == "F") {
      one_ts <- Map(function(x, y, z, ...) cDD_catch_solver(B = x, N = y, R = z, ...), x = B[, y-1], y = N[, y-1], z = R[, y-1],
                    MoreArgs = list(FM = FMort[y-1], Kappa = TMB_data$Kappa, Winf = TMB_data$Winf, wk = TMB_data$wk, M = TMB_data$M))
    } else {
      Fout[, y-1] <- mapply(function(x, y, z, ...) optimize(cDD_catch_solver, c(1e-8, max_F), B = x, N = y, R = z, ...)$minimum,
                            x = B[, y-1], y = N[, y-1], z = R[, y-1],
                            MoreArgs = list(Kappa = TMB_data$Kappa, Winf = TMB_data$Winf, wk = TMB_data$wk, M = TMB_data$M, TAC = Catch[y-1]))

      one_ts <- Map(function(x, y, z, w, ...) cDD_catch_solver(B = x, N = y, R = z, FM = w, ...), x = B[, y-1], y = N[, y-1], z = R[, y-1],
                    w = Fout[, y-1],
                    MoreArgs = list(Kappa = TMB_data$Kappa, Winf = TMB_data$Winf, wk = TMB_data$wk, M = TMB_data$M))
    }

    Cpred[, y-1] <- vapply(one_ts, getElement, numeric(1), 1)

    if(y <= p_years) {
      N[, y] <- vapply(one_ts, getElement, numeric(1), 2)
      B[, y] <- vapply(one_ts, getElement, numeric(1), 3)
    }
  }
  Ipred <- TMB_report$q * B * Iobs_err

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

projection_DD_TMB <- projection_DD_SS <- function(Assessment, constrain = c("F", "Catch"), FMort = NULL, Catch = NULL,
                                                  p_years = 50, p_sim = 200, obs_error = NULL, process_error = NULL,
                                                  max_F = 3, seed = 499, ...) {
  TMB_report <- Assessment@TMB_report
  TMB_data <- Assessment@obj$env$data

  # Sample rec_devs and obs_error
  set.seed(seed)
  tau <- ifelse(!is.null(process_error), process_error[1], TMB_report$tau)
  p_log_rec_dev <- matrix(1, p_sim, p_years)
  if(!is.null(tau) && tau > 0) p_log_rec_dev <- exp(matrix(rnorm(p_years * p_sim, 0, tau), p_sim, p_years) - 0.5 * tau^2)

  sigma <- ifelse(!is.null(obs_error), obs_error[1], 0)
  omega <- ifelse(!is.null(obs_error), obs_error[2], TMB_report$omega)
  Iobs_err <- exp(matrix(rnorm(p_years * p_sim, 0, sigma), p_sim, p_years))
  Cobs_err <- exp(matrix(rnorm(p_years * p_sim, 0, omega), p_sim, p_years))

  B <- N <- R <- Cpred <- matrix(NA, p_sim, p_years)
  B[, 1] <- TMB_report$B[length(TMB_report$B)]
  N[, 1] <- TMB_report$N[length(TMB_report$N)]
  R[, 1:TMB_data$k] <- rep(TMB_report$R[(length(TMB_report$R) - TMB_data$k + 1):length(TMB_report$R)], each = p_sim) *
    p_log_rec_dev[, 1:TMB_data$k]

  if(constrain == "F") {
    Fout <- U <- matrix(pmin(FMort, 1 - exp(-max_F)), p_sim, p_years, byrow = TRUE)
    surv <- (1 - U) * TMB_data$S0
  } else {
    Fout <- U <- surv <- matrix(NA, p_sim, p_years)
  }

  for(y in 2:(p_years+1)) {
    if(y+TMB_data$k-1 <= p_years) {
      R[, y+TMB_data$k-1] <- R_pred(B[, y-1], TMB_report$h, TMB_report$R0, TMB_report$B0, TMB_data$SR_type) *
        p_log_rec_dev[, y+TMB_data$k-1]
    }
    if(constrain == "Catch") {
      Fout[,y-1] <- U[, y-1] <- ifelse(Catch[y-1]/B[, y-1] < 1 - exp(-max_F), Catch[y-1]/B[, y-1], 1 - exp(-max_F))
      surv[, y-1] <- (1 - U[, y-1]) * TMB_data$S0
    }
    Cpred[, y-1] <- U[, y-1] * B[, y-1]
    if(y <= p_years) {
      B[, y] <- surv[, y-1] * (TMB_data$Alpha * N[, y-1] + TMB_data$Rho * B[, y-1]) + TMB_data$wk * R[, y]
      N[, y] <- surv[, y-1] * N[, y-1] + R[, y]
    }
  }
  Eff <- -log(1 - U)/TMB_report$q/Assessment@info$E_rescale
  Ipred <- Cpred/Eff * Iobs_err

  new("project", Catch = Cpred * Cobs_err, Index = Ipred, B = B, VB = B, SSB = B, R = R, N = N, FMort = Fout)
}

projection_VPA <- function(...) stop("Projection function for VPA is not yet available.", call. = FALSE)
projection_SCA2 <- function(...) stop("Projection function for SCA2 is not yet available.", call. = FALSE)
