

#' Assessment emulator as a shortcut to model fitting in closed-loop simulation
#' 
#' Functions (class Assessment) that emulate a stock assessment by sampling the operating model biomass and abundance
#' (with observation error, autocorrelation, and bias) instead of fitting a model. This output can then
#' be passed onto a harvest control rule (HCR function). To utilize the shortcut method in closed-loop simulation,
#' use \link{make_MP} with \code{Shortcut} as the Assessment function. \code{Perfect} assumes no error in the 
#' assessment model and is useful for comparing the behavior of different harvest control rules.
#' 
#' @aliases Perfect
#' @param x An index for the objects in \code{Data} when running in \link[MSEtool]{runMSE}.
#' Otherwise, equals to 1 When running an assessment interactively.
#' @param Data An object of class Data.
#' @param method Indicates where the error in the OM is located. For "B", OM biomass is directly sampled with error.
#' For "N", OM abundance-at-age is sampled and biomass subsequently calculated. For "RF", recruitment and F are
#' sampled to calculate abundance and biomass. There is no error in biological parameters for "N" and "RF".
#' @param B_err If method = "B", a vector of length three that specifies the standard deviation (in logspace),
#' autocorrelation, and bias (1 = unbiased) for biomass.
#' @param N_err Same as B_err, but for abundance when method = "N".
#' @param R_err Same as B_err, but for recruitment when method = "RF".
#' @param F_err Same as B_err, but for fishing mortality when method = "RF".
#' @param ... Other arguments (un-used).
#' 
#' @details Currently there is no error in FMSY (the target F in the HCR in SAMtool).
#' 
#' See Wiedenmann et al. (2015) for guidance on the magnitude of error for the shortcut emulator.
#' @examples 
#' Shortcut_4010 <- make_MP(Shortcut, HCR40_10) 
#' Shortcut_Nerr <- make_MP(Shortcut, HCR40_10, method = "N", N_err = c(0.1, 0.1, 1)) # Highly precise!
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
                     R_err = c(0.3, 0.7, 1), F_err = c(0.3, 0.7, 1), ...) {
  OM_ind <- lapply(1:length(sys.calls()), function(xx) try(get("N_P", envir = sys.frames()[[xx]], inherits = FALSE), silent = TRUE)) %>%
    vapply(function(xx) !is.character(xx), logical(1)) %>% which()
  if(length(OM_ind) == 0) {
    return(new("Assessment", opt = "", SD = "", conv = FALSE))
  }
  
  method <- match.arg(method)
  
  n_y <- length(Data@Year)
  year_p <- max(Data@Year) - Data@LHYear
  Hist <- Data@Misc[-c(1:nrow(Data@Cat))]
  #browser(expr = year_p == 0)
  n_age <- Data@MaxAge + 1
  
  SSB_hist <- Hist$StockPars$SSB[x, , , ] %>% apply(2, sum)
  VB_hist <- Hist$StockPars$VBiomass[x, , , ] %>% apply(2, sum)
  B_hist <- Hist$StockPars$Biomass[x, , , ] %>% apply(2, sum)
  
  N_hist <- Hist$StockPars$N[x, , , ] %>% apply(1:2, sum)
  
  Fapical_hist <- Hist$StockPars$FM[x, , , ] %>% apply(2, max)
  
  SSB_P <- get("SSB_P", envir = sys.frames()[[OM_ind]], inherits = FALSE)[x, , , ][, 1 + 0:year_p, , drop = FALSE] %>% apply(2, sum)
  SSB <- c(SSB_hist, SSB_P)
  
  VB_P <- get("VBiomass_P", envir = sys.frames()[[OM_ind]], inherits = FALSE)[x, , , ][, 1 + 0:year_p, , drop = FALSE] %>% apply(2, sum)
  VB <- c(VB_hist, VB_P)
  
  B_P <- get("Biomass_P", envir = sys.frames()[[OM_ind]], inherits = FALSE)[x, , , ][, 1 + 0:year_p, , drop = FALSE] %>% apply(2, sum)
  B <- c(B_hist, B_P)
  
  N_P <- get("N_P", envir = sys.frames()[[OM_ind]], inherits = FALSE)[x, , , ][, 1 + 0:year_p, , drop = FALSE] %>% apply(1:2, sum)
  N <- abind::abind(N_hist, N_P, along = 2)
  
  if(year_p > 0) {
    Fapical_P <- get("FM_P", envir = sys.frames()[[OM_ind]], inherits = FALSE)[x, , , ][, 1:year_p, , drop = FALSE] %>% apply(2, max)
    Fapical <- c(Fapical_hist, Fapical_P)
  } else {
    Fapical <- Fapical_hist
  }
  R <- N[1, ]
  
  if(method == "B") {
    stopifnot(length(B_err) >= 3)
    
    B_dev <- exp(dev_AC(n_y+1, mu = B_err[3], stdev = B_err[1], AC = B_err[2], seed = x * n_y))
    
    SSB_out <- SSB * B_dev
    VB_out <- VB * B_dev
    B_out <- B * B_dev
  } else if(method == "N") {
    stopifnot(length(N_err) >= 3)
    
    N_dev_y <- dev_AC(n_y+1, mu = N_err[3], stdev = N_err[1], AC = N_err[2], seed = x * n_y)
    N_dev_age <- dev_AC(n_age, mu = 1, stdev = N_err[1], AC = 0, seed = x * n_y - 5000)
    N_dev <- exp(outer(N_dev_age, N_dev_y))
    
    N_out <- N * N_dev
    
    Wt_age <- Hist$StockPars$Wt_age[x, , 0 + 1:n_y]
    Mat_age <- Hist$StockPars$Mat_age[x, , 0 + 1:n_y]
    V_age <- Hist$FleetPars$V[x, , 0 + 1:n_y]
    
    SSB_out <- colSums(N_out * Mat_age * Wt_age)
    VB_out <- colSums(N_out * V_age * Wt_age)
    B_out <- colSums(N_out * Wt_age)
    
  } else if(method == "RF") {
    stopifnot(length(R_err) >= 3, length(F_err) >= 3)
    
    R_dev <- exp(dev_AC(n_y+1, mu = R_err[3], stdev = R_err[1], AC = R_err[2], seed = x * n_y))
    F_dev <- exp(dev_AC(n_y, mu = F_err[3], stdev = F_err[1], AC = F_err[2], seed = x * n_y - 5000))
    
    R_out <- R * R_dev
    F_out <- Fapical * F_dev
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
  Assessment <- new("Assessment", Model = paste("Shortcut", method), conv = TRUE,
                    FMSY = Hist$ReferencePoints$ReferencePoints$FMSY[x],
                    SSB = structure(SSB_out, names = Year_plusone),
                    VB = structure(VB_out, names = Year_plusone),
                    B = structure(B_out, names = Year_plusone),
                    Selectivity = Hist$FleetPars$V[x, , ] %>% t(),
                    Obs_Catch = Data@Cat[x, ],
                    Obs_Index = Data@Ind[x, ],
                    Obs_C_at_age = Data@CAA[x, , ],
                    opt = "No assessment.",
                    SD = "No assessment.")
  if(exists("N_out", inherits = FALSE)) Assessment@N_at_age <- t(N_out)
  if(exists("R_out", inherits = FALSE)) Assessment@R <- structure(R_out, names = Year_plusone)
  if(exists("F_out", inherits = FALSE)) Assessment@FMort <- structure(F_out, names = Year)
  
  if(method == "B") {
    mag_bias <- B_err[3]
  } else if(method == "N") {
    mag_bias <- N_err[3]
  } else {
    mag_bias <- R_err[3]
  }
  
  ref_pt <- c("N0", "B0", "SSB0", "VB0", "MSY", "SSBMSY", "BMSY", "VBMSY")
  lapply(ref_pt, function(xx) slot(Assessment, xx) <<- mag_bias * getElement(Hist$ReferencePoints$ReferencePoints, xx)[x])
  Assessment@R0 <- mag_bias * Hist$StockPars$R0[x]
  
  Assessment@F_FMSY <- Assessment@FMort/Assessment@FMSY
  Assessment@B_BMSY <- Assessment@B/Assessment@BMSY
  Assessment@SSB_SSBMSY <- Assessment@SSB/Assessment@SSBMSY
  Assessment@VB_VBMSY <- Assessment@VB/Assessment@VBMSY
  Assessment@B_B0 <- Assessment@B/Assessment@B0
  Assessment@SSB_SSB0 <- Assessment@SSB/Assessment@SSB0
  Assessment@VB_VB0 <- Assessment@VB/Assessment@VB0
  return(Assessment)
}
class(Shortcut) <- "Assess"

#' @rdname Shortcut
#' @export
Perfect <- function(x, Data, ...) {
  out <- Shortcut(x, Data, method = "RF", R_err = c(0, 0, 1), F_err = c(0, 0, 1))
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
  
  N <- matrix(NA_real_, nrow = n_age, ncol = n_y + 1) # age x yr
  N[, 1] <- Hist$StockPars$N[x, , 1, ] %>% rowSums()
  for(y in 1:n_y) {
    N[1, y+1] <- R_out[y+1]
    for(a in 2:n_age - 1) N[a+1, y+1] <- N[a, y] * exp(-V[a, y] * F_out[y] - M[a, y])
    N[n_age, y+1] <- N[n_age, y+1] + N[n_age, y] * exp(-V[n_age, y] * F_out[y] - M[n_age, y])
  }
  
  SSB <- colSums(N * Wt * Mat)
  VB <- colSums(N * Wt * V)
  B <- colSums(N * Wt)
  
  return(list(N = N, SSB = SSB, VB = VB, B = B))
}