

#' Assessment emulator as a shortcut to model fitting in closed-loop simulation
#' 
#' Functions (class Assessment) that emulate a stock assessment by sampling the operating model biomass, abundance, and 
#' fishing mortality (with observation error, autocorrelation, and bias) instead of fitting a model. This output can then
#' be passed onto a harvest control rule (HCR function). \code{Shortcut} is the base function that samples the OM with an error 
#' distribution. \code{Shortcut2}, the more preferable option, fits \link{SCA} in the last historical year of the operating 
#' model, estimates the error parameters using a vector autoregressive model of the residuals, and then generates model "estimates"
#' using \link[vars]{predict.varest}. \code{Perfect} assumes no error in the assessment model and is useful for comparing the behavior of 
#' different harvest control rules. To utilize the shortcut method in closed-loop simulation, use \link{make_MP} with these functions as 
#' the Assessment model. \strong{N.B. the functions do not work with} \code{runMSE(parallel = TRUE)}.
#' 
#' @aliases Perfect
#' @param x An index for the objects in \code{Data} when running in \link[MSEtool]{runMSE}.
#' Otherwise, equals to 1 When running an assessment interactively.
#' @param Data An object of class Data.
#' @param method Indicates where the error in the OM is located. For "B", OM biomass is directly sampled with error.
#' For "N", OM abundance-at-age is sampled and biomass subsequently calculated. For "RF", recruitment and F are
#' sampled to calculate abundance and biomass. There is no error in biological parameters for "N" and "RF". By default,
#' "B" is used for \code{Shortcut} and "N" for \code{Shortcut2}.
#' @param B_err If \code{method = "B"}, a vector of length three that specifies the standard deviation (in logspace),
#' autocorrelation, and bias (1 = unbiased) for biomass.
#' @param N_err Same as B_err, but for abundance when \code{method = "N"}.
#' @param R_err Same as B_err, but for recruitment when \code{method = "RF"}.
#' @param F_err Same as B_err. Always used regardless of \code{method} to report F and selectivity for HCR.
#' @param VAR_model An object returned by \link[vars]{VAR} to generate emulated assessment error. Used by \code{Shortcut2}. 
#' @param ... Other arguments (not currently used).
#' @author Q. Huynh
#' @details Currently there is no error in FMSY (frequently the target F in the HCR).
#' 
#' See Wiedenmann et al. (2015) for guidance on the magnitude of error for the shortcut emulator.
#' @examples 
#' Shortcut_4010 <- make_MP(Shortcut, HCR40_10) 
#' Shortcut_Nerr <- make_MP(Shortcut, HCR40_10, method = "N", N_err = c(0.1, 0.1, 1)) # Highly precise!
#' 
#' # Fits SCA first and then emulate it in the projection period 
#' Shortcut2_4010 <- make_MP(Shortcut2, HCR40_10) 
#' 
#' \donttest{
#' # Compare the shortcut method vs. fitting an SCA model with a 40-10 control rule
#' MSE <- runMSE(testOM, MPs = c("Shortcut_4010", "SCA_4010"))
#' }
#' 
#' # Compare the performance of three HCRs
#' Perfect_4010 <- make_MP(Perfect, HCR40_10)
#' Perfect_6020 <- make_MP(Perfect, HCR60_20)
#' Perfect_8040MSY <- make_MP(Perfect, HCR_ramp, OCP_type = "SSB_SSBMSY", TOCP = 0.8, LOCP = 0.4)
#' 
#' \donttest{
#' MSE <- runMSE(testOM, MPs = c("Perfect_4010", "Perfect_6020", "Perfect_8040MSY"))
#' }
#' @return An object of class \linkS4class{Assessment}.
#' @references
#' Wiedenmann, J., Wilberg, M.J., Sylvia, A., and Miller, T.J. 2015. Autocorrelated error in stock assessment 
#' estimates: Implications for management strategy evaluation. Fisheries Research 172: 325-334.
#' @export
Shortcut <- function(x = 1, Data, method = c("B", "N", "RF"), B_err = c(0.3, 0.7, 1), N_err = c(0.3, 0.7, 1), 
                     R_err = c(0.3, 0.7, 1), F_err = c(0.3, 0.7, 1), VAR_model, ...) {
  OM_ind <- lapply(1:length(sys.calls()), function(xx) try(get("N_P", envir = sys.frames()[[xx]], inherits = FALSE), silent = TRUE)) %>%
    vapply(function(xx) !is.character(xx), logical(1)) %>% which()
  if(!length(OM_ind)) {
    stop("No operating model was found.")
    #return(new("Assessment", opt = "", SD = "", conv = FALSE))
  }
  
  method <- match.arg(method)
  
  n_y <- length(Data@Year) # Should be year_p + n_hist - 1
  year_p <- max(Data@Year) - Data@LHYear + 1 # You are at the beginning of this year, have the abundance but not the mortality/catches
  Hist <- Data@Misc[-c(1:nrow(Data@Cat))]
  n_hist <- dim(Hist$StockPars$SSB)[3]
  
  n_age <- Data@MaxAge + 1
  
  SSB_hist <- Hist$StockPars$SSB[x, , , ] %>% apply(2, sum)
  VB_hist <- Hist$StockPars$VBiomass[x, , , ] %>% apply(2, sum)
  B_hist <- Hist$StockPars$Biomass[x, , , ] %>% apply(2, sum)
  
  N_hist <- Hist$StockPars$N[x, , , ] %>% apply(1:2, sum)
  
  Fapical_hist <- Hist$StockPars$FM[x, , , ] %>% apply(2, max)
  
  SSB_P <- get("SSB_P", envir = sys.frames()[[OM_ind]], inherits = FALSE)[x, , , ][, 1:year_p, , drop = FALSE] %>% apply(2, sum)
  SSB <- c(SSB_hist, SSB_P)
  
  VB_P <- get("VBiomass_P", envir = sys.frames()[[OM_ind]], inherits = FALSE)[x, , , ][, 1:year_p, , drop = FALSE] %>% apply(2, sum)
  VB <- c(VB_hist, VB_P)
  
  B_P <- get("Biomass_P", envir = sys.frames()[[OM_ind]], inherits = FALSE)[x, , , ] [, 1:year_p, , drop = FALSE] %>% apply(2, sum)
  B <- c(B_hist, B_P)
  
  N_P <- get("N_P", envir = sys.frames()[[OM_ind]], inherits = FALSE)[x, , , ][, 1:year_p, , drop = FALSE] %>% apply(1:2, sum)
  N <- cbind(N_hist, N_P)
  
  if(year_p > 1) { # Need OM with equal areas
    Fapical_P <- get("FM_P", envir = sys.frames()[[OM_ind]], inherits = FALSE)[x, , , ][, 2:year_p - 1, , drop = FALSE] %>% apply(2, max)
    Fapical <- c(Fapical_hist, Fapical_P)
  } else {
    Fapical <- Fapical_hist
  }
  R <- N[1, ]
  
  if(!missing(VAR_model) && year_p > 1) {
    VAR_proj <- predict(VAR_model, n.ahead = year_p)
    F_dev <- exp(c(VAR_proj$endog[, "FM"], VAR_proj$fcst$FM[2:year_p - 1, "fcst"]))
  } else {
    stopifnot(length(F_err) >= 3)
    F_dev <- exp(dev_AC(n_y, mu = F_err[3], stdev = F_err[1], AC = F_err[2], seed = x * n_y - 5000))
  }
  F_out <- Fapical * F_dev
  
  if(method == "B") {
    
    if(!missing(VAR_model) && year_p > 1) {
      B_dev <- exp(c(VAR_proj$endog[, "B"], VAR_proj$fcst$B[, "fcst"]))
    } else {
      stopifnot(length(B_err) >= 3)
      B_dev <- exp(dev_AC(n_y+1, mu = B_err[3], stdev = B_err[1], AC = B_err[2], seed = x * n_y))
    }
    SSB_out <- SSB * B_dev
    VB_out <- VB * B_dev
    B_out <- B * B_dev
    
  } else if(method == "N") {
    
    if(!missing(VAR_model) && year_p > 1) {
      N_dev_proj <- vapply(paste0("N.", 1:n_age), function(xx) getElement(VAR_proj$fcst, xx)[, "fcst"], numeric(year_p))
      N_dev <- exp(rbind(VAR_proj$endog[, 1:n_age], N_dev_proj)) %>% t()
    } else {
      stopifnot(length(N_err) >= 3)
      
      N_dev_y <- dev_AC(n_y+1, mu = N_err[3], stdev = N_err[1], AC = N_err[2], seed = x * n_y)
      N_dev_age <- dev_AC(n_age, mu = 1, stdev = N_err[1], AC = 0, seed = x * n_y - 5000)
      N_dev <- exp(outer(N_dev_age, N_dev_y))
    }
    N_out <- N * N_dev
    
    Wt_age <- Hist$StockPars$Wt_age[x, , 0:n_y + 1]
    Mat_age <- Hist$StockPars$Mat_age[x, , 0:n_y + 1]
    V_age <- Hist$FleetPars$V_real[x, , 0:n_y + 1]
    Fec_age <- Hist$StockPars$Fec_Age[x, , 0:n_y + 1]
    
    SSB_out <- colSums(N_out * Fec_age)
    VB_out <- colSums(N_out * V_age * Wt_age)
    B_out <- colSums(N_out * Wt_age)
    
  } else if(method == "RF") {
    if(!missing(VAR_model) && year_p > 1) {
      R_dev <- exp(c(VAR_proj$endog[, "R"], VAR_proj$fcst$R[, "fcst"]))
    } else {
      stopifnot(length(R_err) >= 3)
      R_dev <- exp(dev_AC(n_y+1, mu = R_err[3], stdev = R_err[1], AC = R_err[2], seed = x * n_y))
    }
    R_out <- R * R_dev
    ASM_out <- project_ASM(x, R_out, F_out, Hist, Data)
    
    N_out <- ASM_out$N
    
    SSB_out <- ASM_out$SSB
    VB_out <- ASM_out$VB
    B_out <- ASM_out$B
  }
  
  # For assessment diagnostics
  Year <- Data@Year
  Year_plusone <- c(Year, max(Year) + 1)
  opt <- SD <- "No assessment."
  
  Assessment <- new("Assessment", 
                    Model = paste(ifelse(missing(VAR_model), "Shortcut", "Shortcut2"), method), 
                    conv = TRUE,
                    FMort = structure(F_out, names = Year),
                    SSB = structure(SSB_out, names = Year_plusone),
                    VB = structure(VB_out, names = Year_plusone),
                    B = structure(B_out, names = Year_plusone),
                    Selectivity = Hist$FleetPars$V[x, , 0:n_y + 1] %>% t(),
                    Obs_Catch = Data@Cat[x, ],
                    Obs_Index = Data@Ind[x, ],
                    Obs_C_at_age = Data@CAA[x, , ],
                    opt = "No assessment.",
                    SD = "No assessment.")
  if(exists("N_out", inherits = FALSE)) Assessment@N_at_age <- t(N_out)
  if(exists("R_out", inherits = FALSE)) {
    Assessment@R <- structure(R_out, names = Year_plusone)
  } else {
    Assessment@R <- structure(Assessment@N_at_age[, 1], names = Year_plusone)
  }
  
  if(missing(VAR_model)) {
    if(method == "B") {
      mag_bias <- B_err[3]
    } else if(method == "N") {
      mag_bias <- N_err[3]
    } else {
      mag_bias <- R_err[3]
    }
  } else {
    mag_bias <- 1
  }
  
  ref_pt <- c("N0", "B0", "SSB0", "VB0", "MSY", "SSBMSY", "BMSY", "VBMSY", "FMSY")
  lapply(ref_pt, function(xx) slot(Assessment, xx) <<- mag_bias * Hist$ReferencePoints$ByYear[[xx]][x, n_y + 1])
  Assessment@TMB_report$dynamic_SSB0 <- Hist$ReferencePoints$Dynamic_Unfished$SSB0[x, 1 + 0:n_y] %>% structure(names = Year_plusone)
  
  Assessment@R0 <- mag_bias * Hist$StockPars$R0[x]
  
  Assessment@F_FMSY <- Assessment@FMort/Assessment@FMSY
  Assessment@B_BMSY <- Assessment@B/Assessment@BMSY
  Assessment@SSB_SSBMSY <- Assessment@SSB/Assessment@SSBMSY
  Assessment@VB_VBMSY <- Assessment@VB/Assessment@VBMSY
  Assessment@B_B0 <- Assessment@B/Assessment@B0
  Assessment@SSB_SSB0 <- Assessment@SSB/Assessment@SSB0
  Assessment@VB_VB0 <- Assessment@VB/Assessment@VB0
  
  if(method == "B") {
    catch_eq <- function(Ftarget) {
      catch_equation(method = "frac", Utarget = 1 - exp(-Ftarget), B = Assessment@VB[length(Assessment@VB)])
    }
  } else {
    Assessment@info <- list(data = list(M = Hist$StockPars$M_ageArray[x, , n_y + 1], 
                                        weight = Hist$FleetPars$Wt_age_C[x, , n_y + 1]))
    catch_eq <- function(Ftarget) {
      catch_equation(method = "Baranov", Ftarget = Ftarget,
                     sel = Assessment@Selectivity[nrow(Assessment@Selectivity), ], 
                     M = Assessment@info$data$M, wt = Assessment@info$data$weight, N = Assessment@N_at_age[nrow(Assessment@N_at_age), ])
    }
  }
  
  F_SPR <- Data@Misc$ReferencePoints$ByYear$F_SPR[x, , n_y + 1] %>% rev()
  SPR <- F_SPR %>% names() %>% substr(3,4) %>% as.numeric()
  SPR <- rev(SPR/100)
  if(all(F_SPR != 0)) {
    F_SPR <- c(0, F_SPR) 
    SPR <- c(1, SPR)
  }
  Assessment@forecast <- list(per_recruit = data.frame(FM = F_SPR, SPR = SPR, 
                                                       F01 = Data@Misc$ReferencePoints$ByYear$F01_YPR[x, n_y + 1],
                                                       Fmax = Data@Misc$ReferencePoints$ByYear$Fmax_YPR[x, n_y + 1]), 
                              catch_eq = catch_eq)
  return(Assessment)
}
class(Shortcut) <- "Assess"

#' @rdname Shortcut
#' @param SCA_args Additional arguments to pass to \link{SCA}. Currently, arguments \code{SR} and \code{vulnerability}
#' are obtained from the operating model.
#' @param VAR_args Additional arguments to pass to \link[vars]{VAR}. By default, argument \code{type = "none"} 
#' (stationary time series with mean zero is assumed).
#' @importFrom vars VAR
#' @export
Shortcut2 <- function(x, Data, method = "N", SCA_args = list(), VAR_args = list(type = "none"), ...) {
  method <- match.arg(method, choices = c("N", "B", "RF"))
  
  if(!is.null(Data@Misc[[x]]$VAR_model)) {
    run_Shortcut <- Shortcut(x = x, Data = Data, method = method, VAR_model = Data@Misc[[x]]$VAR_model)
    
    run_Shortcut@info$Misc$VAR_model <- Data@Misc[[x]]$VAR_model
    return(run_Shortcut)
  } else if(max(Data@Year) == Data@LHYear) {
    SCA_args$x <- x
    SCA_args$Data <- Data
    if(is.null(SCA_args$SR) && !is.null(Data@Misc$StockPars$SRrel)) {
      SCA_args$SR <- ifelse(Data@Misc$StockPars$SRrel[x] == 1, "BH", "Ricker")
    }
    if(is.null(SCA_args$vulnerability) && !is.null(Data@OM$Vmaxlen)) {
      SCA_args$vulnerability <- ifelse(Data@OM$Vmaxlen[x] != 1, "dome", "logistic")
    }
    
    run_SCA <- do.call("SCA", SCA_args)
    
    F_est <- run_SCA@FMort
    F_OM <- Data@Misc$FleetPars$Find[x, ] * Data@Misc$FleetPars$qs[x]
    
    if(method == "B") {
      SSB_est <- run_SCA@SSB[-length(run_SCA@SSB)]
      SSB_OM <- apply(Data@Misc$StockPars$SSB[x, , , ], 2, sum)
      
      var_resid <- data.frame(B = log(SSB_est/SSB_OM), FM = log(F_est/F_OM))
      
    } else if(method == "N") {
      N_est <- run_SCA@N_at_age[-length(run_SCA@SSB), ]
      N_OM <- apply(Data@Misc$StockPars$N[x, , , ], 1:2, sum) %>% t()
      
      var_resid <- data.frame(N = log(N_est/N_OM), FM = log(F_est/F_OM))
      
    } else {
      R_est <- run_SCA@N_at_age[-length(run_SCA@SSB), 1]
      R_OM <- rowSums(Data@Misc$StockPars$N[x, 1, , ])
      
      var_resid <- data.frame(R = log(R_est/R_OM), FM = log(F_est/F_OM))
    }
    VAR_args$y <- as.matrix(var_resid)
    VAR_model <- do.call(vars::VAR, VAR_args)
    
    run_SCA@info$Misc$VAR_model <- VAR_model
    return(run_SCA)
    
  } else {
    stop("No assessment or shortcut method was run.")
  }
}
class(Shortcut2) <- "Assess"


#' @rdname Shortcut
#' @export
Perfect <- function(x, Data, ...) {
  out <- Shortcut(x, Data, method = "N", F_err = c(0, 0, 1), N_err = c(0, 0, 1))
  out@Model <- "Perfect"
  return(out)
}
class(Perfect) <- "Assess"

project_ASM <- function(x, R_out, F_out, Hist, Data) {
  n_y <- length(F_out)
  n_age <- Data@MaxAge + 1
  
  V <- Hist$FleetPars$V[x, , 1 + 0:n_y]
  M <- Hist$StockPars$M_ageArray[x, , 1:n_y]
  Wt <- Hist$StockPars$Wt_age[x, , 1 + 0:n_y]
  Mat <- Hist$StockPars$Mat_age[x, , 1 + 0:n_y]
  
  N <- matrix(NA_real_, nrow = n_age, ncol = n_y + 1)
  N[, 1] <- Hist$StockPars$N[x, , 1, ] %>% rowSums()
  for(y in 1:n_y) {
    N[1, y+1] <- R_out[y+1]
    for(a in 2:n_age - 1) N[a+1, y+1] <- N[a, y] * exp(-V[a, y] * F_out[y] - M[a, y])
    if(Hist$StockPars$plusgroup) N[n_age, y+1] <- N[n_age, y+1] + N[n_age, y] * exp(-V[n_age, y] * F_out[y] - M[n_age, y])
  }
  
  SSB <- colSums(N * Wt * Mat)
  VB <- colSums(N * Wt * V)
  B <- colSums(N * Wt)
  
  return(list(N = N, SSB = SSB, VB = VB, B = B))
}


#calc_err <- function(x_est, x_OM, AC_method) {
#  if(AC_method == "mle") {
#    opt <- nlminb(c(log(0.1), 0), err_likelihood, x_est = x_est, x_OM = x_OM, n = length(x_est))
#    c(exp(opt$par[1]), ilogit2(opt$par[2], -1, 1, 0), 1)
#  } else {
#    SD <- sd(log(x_est/x_OM))
#    AC <- acf(log(x_est/x_OM), lag.max = 1, plot = FALSE)$acf[2]
#    Bias <- mean(x_est/x_OM)
#    c(SD, AC, Bias)
#  }
#}
#
## Equation 2 of Wiedenmann et al. 2015 but drop constant terms
#err_likelihood <- function(x, x_est, x_OM, n) {
#  SD <- exp(x[1])
#  AC <- ilogit2(x[2], -1, 1, 0)
#  
#  v1 <- log(x_OM[2:n]) - AC * log(x_OM[2:n - 1]) - log(x_est[2:n]) + AC * log(x_est[2:n-1]) 
#  v2 <- log(x_OM[1]/x_est[1])
#  obj <- -n * log(SD) + 0.5 * log(1 - AC^2) - 0.5 * sum(v1^2)/SD/SD - 0.5 * (1 - AC^2) * v2^2 /SD/SD
#  return(-obj)
#}