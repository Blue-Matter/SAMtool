#' Virtual population analysis (VPA)
#'
#' A VPA model that back-calculates abundance-at-age assuming that the catch-at-age is known without error and tuned to an index.
#' The population dynamics equations are primarily drawn from VPA-2BOX (Porch 2018). MSY reference points and per-recruit quantities
#' are then calculated from the VPA output.
#'
#' @param x A position in the Data object (by default, equal to one for assessments).
#' @param Data An object of class Data
#' @param AddInd A vector of integers or character strings indicating the indices to be used in the model. Integers assign the index to
#' the corresponding index in Data@@AddInd, "B" (or 0) represents total biomass in Data@@Ind, "VB" represents vulnerable biomass in
#' Data@@VInd, and "SSB" represents spawning stock biomass in Data@@SpInd.
#' @param expanded Whether the catch at age in \code{Data} has been expanded. If \code{FALSE}, then the catch in weight
#' should be provided in \code{Data@@Cat} so that the function can calculate annual expansion factors.
#' @param SR Stock-recruit function (either \code{"BH"} for Beverton-Holt or \code{"Ricker"}) for calculating MSY reference points.
#' @param vulnerability Whether the terminal year vulnerability is \code{"logistic"} or \code{"dome"} (double-normal). If \code{"free"},
#' independent F's are calculated in the terminal year (subject to the assumed ratio of F of the plus-group to the previous age class).
#' See details for parameterization.
#' @param start Optional list of starting values. Entries can be expressions that are evaluated in the function. See details.
#' @param fix_Fratio Logical, whether the ratio of F of the plus-group to the previous age class is fixed in the model.
#' @param fix_h Logical, whether to fix steepness to value in \code{Data@@steep}. This only affects
#' calculation of MSY and unfished reference points.
#' @param fix_Fterm Logical, whether to fix the value of the terminal F.
#' @param LWT A vector of likelihood weights for each survey.
#' @param shrinkage A named list of up to length 2 to constrain parameters:
#' \itemize{
#' \item \code{vul} - a length two vector that constrains the vulnerability-at-age in the most recent years. The first number
#' is the number of years in which vulnerability will be constrained (as a random walk in log space), the second number is the standard deviation of the random walk.
#' The default 
#' \item \code{R} - a length two vector that constrains the recruitment estimates in the most recent years. The first number
#' is the number of years in which recruitment will be constrained (as a random walk in log space), the second number is the standard deviation of the random walk.
#' }
#' @param refpt A named list of how many years to average parameters for calculating reference points, yield per recruit, and spawning potential ratio:
#' \itemize{
# #' \item \code{weight} An integer for the number of most recent years to average the weight-at-age schedule (default is 3).
#' \item \code{vul} An integer for the number of most recent years to average the vulnerability schedule (default is 3).
#' \item \code{R} A length two for the quantile used to calculate recruitment in the year following the terminal year and the number of years
#' from which that quantile is used, i.e., \code{c(0.5, 5)} is the default that calculates median recruitment from the most recent 5 years of the model.
#' }
#' @param n_itF The number of iterations for solving F in the model (via Newton's method).
#' @param min_age An integer to specify the smallest age class in the VPA. By default, the youngest age with non-zero CAA in the terminal year is used.
#' @param max_age An integer to specify the oldest age class in the VPA. By default, the oldest age with non-zero CAA for all years is used.
#' @param silent Logical, passed to \code{\link[TMB]{MakeADFun}}, whether TMB
#' will print trace information during optimization. Used for diagnostics for model convergence.
#' @param opt_hess Logical, whether the hessian function will be passed to \code{\link[stats]{nlminb}} during optimization
#' (this generally reduces the number of iterations to convergence, but is memory and time intensive and does not guarantee an increase
#' in convergence rate). Ignored if \code{integrate = TRUE}.
#' @param n_restart The number of restarts (calls to \code{\link[stats]{nlminb}}) in the optimization procedure, so long as the model
#' hasn't converged. The optimization continues from the parameters from the previous (re)start.
#' @param control A named list of arguments for optimization to be passed to
#' \code{\link[stats]{nlminb}}.
#' @param ... Other arguments to be passed.
#' @details
#' The VPA is initialized by estimating the terminal F-at-age. Parameter \code{Fterm} is the apical terminal F if
#' a functional form for vulnerability is used in the terminal year, i.e., when \code{vulnerability = "logistic"} or \code{"free"}.
#' If the terminal F-at-age are otherwise independent parameters,
#' \code{Fterm} is the F for the reference age which is half the maximum age. Once terminal-year abundance is
#' estimated, the abundance in historical years can be back-calculated. The oldest age group is a plus-group, and requires
#' an assumption regarding the ratio of F's between the plus-group and the next youngest age class. The F-ratio can
#' be fixed (default) or estimated.
#'
#' For \code{start} (optional), a named list of starting values of estimates can be provided for:
#' \itemize{
#' \item \code{Fterm} The terminal year fishing mortality. This is the apical F when \code{vulnerability = "logistic"} or \code{"free"}.
#' \item \code{Fratio} The ratio of F in the plus-group to the next youngest age. If not provided, a value of 1 is used.
#' \item \code{vul_par} Vulnerability parameters in the terminal year. This will be of length 2 vector for \code{"logistic"} or length 4 for
#' \code{"dome"}, see \link{SCA} for further documentation on parameterization. For option \code{"free"}, this will be a vector of length
#' A-2 where A is the number of age classes in the model. To estimate parameters, vulnerability is initially set to one at half the max age
#' (and subsequently re-calculated relative to the maximum F experienced in that year). Vulnerability in the plus-group is also constrained
#' by the Fratio.
#' }
#' 
#' MSY and depletion reference points are calculated by fitting the stock recruit relationship to the recruitment and SSB estimates. Per-recruit
#' quantities are also calculated, which may be used in harvest control rules.
#' 
#' @section Additional considerations:
#' The VPA tends to be finicky to implement straight out of the box. For example, zeros in plusgroup age in the catch-at-age
#' model will crash the model, as well as if the catch-at-age values are close to zero. The model sets F-at-age to 1e-4
#' if any catch-at-age value < 1e-4.
#' 
#' It is recommended to do some preliminary fits with the VPA before running simulations en masse. See example below.
#' 
#' Shrinkage, penalty functions that stabilize model estimates of recruitment and selectivity year-over-year near 
#' the end of the time series, alters the behavior of the model. This is something to tinker with in your initial
#' model fits, and worth evaluating in closed-loop simulation. 
#' 
#' @section Online Documentation:
#' Model description and equations are available on the openMSE 
#' \href{https://openmse.com/features-assessment-models/4-vpa/}{website}.
#' 
#' @return An object of class \linkS4class{Assessment}. The F vector is the apical fishing mortality experienced by any
#' age class in a given year. 
#' @examples
#' \donttest{
#' OM <- MSEtool::testOM
#' 
#' # Simulate logistic normal age comps with CV = 0.1
#' # (set CAA_ESS < 1, which is interpreted as a CV)
#' OM@CAA_ESS <- c(0.1, 0.1) 
#' Hist <- MSEtool::Simulate(OM, silent = TRUE)
#' 
#' # VPA max age is 15 (Hist@Data@MaxAge)
#' m <- VPA(x = 2, Data = Hist@Data, vulnerability = "dome")
#' 
#' # Use age-9 as the VPA max age instead
#' m9 <- VPA(x = 2, Data = Hist@Data, vulnerability = "dome", max_age = 9)
#' 
#' compare_models(m, m9)
#' }
#' @references
#' Porch, C.E. 2018. VPA-2BOX 4.01 User Guide. NOAA Tech. Memo. NMFS-SEFSC-726. 67 pp.
#' @export
VPA <- function(x = 1, Data, AddInd = "B", expanded = FALSE, SR = c("BH", "Ricker"), vulnerability = c("logistic", "dome", "free"),
                start = list(), fix_h = TRUE, fix_Fratio = TRUE, fix_Fterm = FALSE, LWT = NULL, shrinkage = list(), n_itF = 5L,
                min_age = "auto", max_age = "auto", refpt = list(),
                silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                control = list(iter.max = 2e5, eval.max = 4e5), ...) {
  dependencies <- "Data@Cat, Data@CAA, Data@Ind, Data@Mort, Data@L50, Data@L95, Data@CAA, Data@vbK, Data@vbLinf, Data@vbt0, Data@wla, Data@wlb, Data@MaxAge"
  dots <- list(...)
  if (!is.null(dots$nitF)) n_itF <- dots$nitF
  start <- lapply(start, eval, envir = environment())
  
  vulnerability <- match.arg(vulnerability)
  SR <- match.arg(SR)

  # Life history
  if (is.null(refpt$weight)) refpt$weight <- 3L
  if (is.null(refpt$vul)) refpt$vul <- 3L
  if (is.null(refpt$R)) refpt$R <- c(0.5, 5)
  n_age <- Data@MaxAge + 1
  M <- rep(Data@Mort[x], n_age)
  a <- Data@wla[x]
  b <- Data@wlb[x]
  Linf <- Data@vbLinf[x]
  K <- Data@vbK[x]
  t0 <- Data@vbt0[x]
  La <- Linf * (1 - exp(-K * (c(0:Data@MaxAge) - t0)))
  Wa <- a * La ^ b
  A50 <- min(0.5 * Data@MaxAge, iVB(t0, K, Linf, Data@L50[x]))
  A95 <- max(A50+0.5, iVB(t0, K, Linf, Data@L95[x]))
  mat_age <- 1/(1 + exp(-log(19) * (c(0:Data@MaxAge) - A50)/(A95 - A50)))
  mat_age <- mat_age/max(mat_age)

  if (any(names(dots) == "yind")) {
    yind <- eval(dots$yind)
  } else {
    yind <- 1:nrow(Data@CAA[x, , ])
  }

  Year <- Data@Year[yind]
  CAA_hist <- Data@CAA[x, yind, ]
  if (all(is.na(CAA_hist))) stop("No age composition data found in Data object.", call. = FALSE)

  if (!expanded) {
    C_hist <- Data@Cat[x, yind]
    if (any(is.na(C_hist) | C_hist < 0)) warning("Error. Catch time series is not complete.")
    expansion_factors <- C_hist/colSums(t(CAA_hist) * Wa)
    CAA_hist <- CAA_hist * expansion_factors
  }
  
  # Reduce max_age until no zero's are observed
  CAA_hist_VPA <- CAA_hist
  if (is.character(max_age) && max_age == "auto") {
    max_age <- n_age - 1
    while (any(CAA_hist_VPA[, max_age + 1] <= 0)) {
      max_age <- max_age - 1
      CAA_hist_VPA <- CAA_hist[, 0:max_age + 1]
      CAA_hist_VPA[, max_age + 1] <- rowSums(CAA_hist[, (max_age+1):ncol(CAA_hist), drop = FALSE])
    }
  } else if (is.numeric(max_age) && length(max_age) == 1) {
    CAA_hist_VPA <- CAA_hist[, 0:max_age + 1]
    CAA_hist_VPA[, max_age + 1] <- rowSums(CAA_hist[, (max_age+1):ncol(CAA_hist), drop = FALSE])
  } else {
    stop("max_age must be an integer or \"auto\".")
  }
  
  # Increase min-age until no zero's in terminal year
  if (is.character(min_age) && min_age == "auto") {
    min_age <- 0
    while (CAA_hist[nrow(CAA_hist), min_age + 1] <= 0) min_age <- min_age + 1
    #if (colSums(CAA_hist_VPA)[1] < 0.01 * sum(CAA_hist_VPA)) min_age <- min_age + 1
  } 
  if (is.numeric(min_age) && length(min_age) == 1) {
    if (min_age > 0) {
      CAA_hist_VPA <- CAA_hist_VPA[, -c(0:min_age + 1), drop = FALSE]
      CAA_hist_VPA <- CAA_hist[, 0:min_age + 1, drop = FALSE] %>% rowSums() %>% cbind(CAA_hist_VPA)
      if (ncol(CAA_hist_VPA) == 1) stop("Only one age class left after consolidating plus- and minus- groups to remove zeros.")
    }
  } else {
    stop("min_age must be an integer or \"auto\".")
  }
    
  ages <- min_age:max_age
  
  # Any missing zeros
  CAA_hist_VPA[is.na(CAA_hist_VPA) | CAA_hist_VPA < 1e-8] <- 1e-8
  
  # Avoid numerical instability with maxage and maxage-1
  maxage_ind <- CAA_hist_VPA[, ncol(CAA_hist_VPA) - 1] == 1e-8 
  CAA_hist_VPA[maxage_ind, ncol(CAA_hist_VPA) - 1] <- CAA_hist_VPA[maxage_ind, ncol(CAA_hist_VPA)]
  
  update_age_schedule <- function(x, ages) {
    xout <- x[ages + 1]
    xout[length(xout)] <- mean(x[(max(ages)+1):length(x)])
    return(xout %>% structure(names = ages))
  }

  LH <- list(LAA = update_age_schedule(La, ages), WAA = update_age_schedule(Wa, ages), 
             Linf = Linf, K = K, t0 = t0, a = a, b = b, A50 = A50, A95 = A95, 
             mat = update_age_schedule(mat_age, ages))
  
  Ind <- lapply(AddInd, Assess_I_hist, Data = Data, x = x, yind = yind)
  I_hist <- do.call(cbind, lapply(Ind, getElement, "I_hist"))
  I_sd <- do.call(cbind, lapply(Ind, getElement, "I_sd")) %>% pmax(0.05)
  I_units <- do.call(c, lapply(Ind, getElement, "I_units"))
  I_vul <- vapply(AddInd, function(xx) {
    if (xx == "B") {
      return(rep(1, length(ages)))
    } else if (xx == "SSB") {
      return(LH$mat)
    } else if (xx == "VB") {
      return(rep(0, length(ages)))
    } else {
      return(Data@AddIndV[x, suppressWarnings(as.numeric(xx)), ages + 1])
    }
  }, numeric(length(ages)))
  nsurvey <- ncol(I_hist)
  if (is.null(LWT)) LWT <- rep(1, nsurvey)
  if (length(LWT) != nsurvey) stop("LWT needs to be a vector of length ", nsurvey)
  
  # Shrinkage
  if (is.null(shrinkage$vul)) shrinkage$vul <- c(3, 0.4)
  if (is.null(shrinkage$R)) {
    sigmaR <- Data@sigmaR[x]
    if (is.null(sigmaR)) sigmaR <- 0.6
    shrinkage$R <- c(3, sigmaR)
  }
  
  data <- list(model = "VPA", I_hist = I_hist, I_sd = I_sd, I_units = I_units, I_vul = I_vul, 
               abs_I = rep(0, nsurvey), nsurvey = nsurvey, LWT = LWT,
               CAA_hist = CAA_hist_VPA, n_y = length(Year), n_age = length(ages),
               M = update_age_schedule(M, ages), weight = LH$WAA, 
               vul_type_term = vulnerability, n_itF = as.integer(n_itF),
               n_vulpen = shrinkage$vul[1], vulpen = shrinkage$vul[2], 
               n_Rpen = shrinkage$R[1], Rpen = shrinkage$R[2])

  # Starting values
  params <- list()
  if (!is.null(start)) {
    if (!is.null(start$Fterm) && is.numeric(start$Fterm)) params$log_Fterm <- log(start$Fterm)
    if (!is.null(start$Fratio) && is.numeric(start$Fratio)) params$log_Fratio <- log(start$Fratio)
    if (!is.null(start$vul_par) && is.numeric(start$vul_par)) {
      if (vulnerability == "logistic") {
        if (start$vul_par[1] > 0.75 * max_age) stop("start$vul_par[1] needs to be less than 0.75 * max_age (", max_age, ")")
        if (length(start$vul_par) < 2) stop("Two parameters needed for start$vul_par with logistic vulnerability (see help).")
        if (start$vul_par[1] <= start$vul_par[2]) stop("start$vul_par[1] needs to be greater than start$vul_par[2] (see help).")

        params$vul_par <- c(logit(start$vul_par[1]/max_age/0.75), log(start$vul_par[1] - start$vul_par[2]))
      } else if (vulnerability == "dome") {
        if (start$vul_par[1] > 0.75 * max_age) stop("start$vul_par[1] needs to be less than 0.75 * max_age (", max_age, ")")
        if (length(start$vul_par) < 4) stop("Four parameters needed for start$vul_par with dome vulnerability (see help).")
        if (start$vul_par[1] <= start$vul_par[2]) stop("start$vul_par[1] needs to be greater than start$vul_par[2] (see help).")
        if (start$vul_par[3] <= start$vul_par[1] || start$vul_par[3] >= max_age) {
          stop("start$vul_par[3] needs to be between start$vul_par[1] and max_age (", max_age, ")")
        }
        if (start$vul_par[4] <= 0 || start$vul_par[4] >= 1) stop("start$vul_par[4] needs to be between 0-1 (see help).")

        params$vul_par <- c(logit(start$vul_par[1]/max_age/0.75), log(start$vul_par[1] - start$vul_par[2]),
                            logit(1/(max_age - start$vul_par[1])), logit(start$vul_par[4]))
      } else if (vulnerability == "free") {
        if (length(start$vul_par) < length(ages)) stop(paste0("start$vul_par needs to be of length", length(ages), "."))
        params$vul_par <- log(start$vul_par[1:length(ages)])
      }
    }
  }

  if (is.null(params$log_Fterm)) params$log_Fterm <- log(0.2)
  if (is.null(params$log_Fratio)) params$log_Fratio <- log(1)
  if (is.null(params$vul_par)) {
    CAA_mode <- ages[which.max(colSums(CAA_hist, na.rm = TRUE))]
    CAA_mode <- ifelse(CAA_mode + 1 > max_age, max_age - 1, CAA_mode)
    if (vulnerability == "free") {
      params$vul_par <- ifelse(ages[-length(ages)] < CAA_mode, 0.5, 1) %>% log()
    } else {
      if ((is.na(Data@LFC[x]) && is.na(Data@LFS[x])) || (Data@LFC[x] > Linf) || (Data@LFS[x] > Linf)) {
        if (vulnerability == "logistic") params$vul_par <- c(logit(CAA_mode/(max_age - min_age)/0.75), log(1))
        if (vulnerability == "dome") {
          params$vul_par <- c(logit(CAA_mode/(max_age - min_age)/0.75), log(1), logit(1/(max_age - CAA_mode)), logit(0.5))
        }
      } else {
        A5 <- min(iVB(t0, K, Linf, Data@LFC[x]), CAA_mode-1)
        Afull <- min(iVB(t0, K, Linf, Data@LFS[x]), 0.5 * max_age)
        A5 <- min(A5, Afull - 0.5)
        A50_vul <- mean(c(A5, Afull))

        if (vulnerability == "logistic") params$vul_par <- c(logit(Afull/(max_age - min_age)/0.75), log(Afull - A50_vul))
        if (vulnerability == "dome") {
          params$vul_par <- c(logit(Afull/max_age/0.75), log(Afull - A50_vul), logit(1/(max_age - Afull)), logit(0.5))
        }
      }
    }
  }

  info <- list(Year = Year, data = data, params = params, LH = LH, SR = SR, ages = ages, control = control,
               fix_h = fix_h, refpt = refpt)

  map <- list()
  if (vulnerability == "free") {
    fixed_ind <- round(0.5 * length(ages))
    free_vul_par <- rep(log(1), length(ages) - 1)
    free_vul_par[fixed_ind] <- NA
    free_vul_par[!is.na(free_vul_par)] <- 1:(length(ages)-2)
    map$vul_par <- factor(free_vul_par)
  } else if (vulnerability == "dome") {
    map$vul_par <- c(1, 2, NA, 3) %>% factor()
  }
  if (fix_Fratio) map$log_Fratio <- factor(NA)
  if (fix_Fterm) map$log_Fterm <- factor(NA)

  obj <- MakeADFun(data = info$data, parameters = info$params, hessian = TRUE,
                   map = map, DLL = "SAMtool", silent = silent)

  mod <- optimize_TMB_model(obj, control, opt_hess, n_restart)
  opt <- mod[[1]]
  SD <- mod[[2]]
  report <- obj$report(obj$env$last.par.best) %>% VPA_posthoc(info = info)
  
  Yearplusone <- c(Year, max(Year) + 1)

  Assessment <- new("Assessment", Model = "VPA",
                    Name = Data@Name, conv = SD$pdHess,
                    FMort = structure(report$F, names = Year),
                    B = structure(report$B, names = Yearplusone), SSB = structure(report$E, names = Yearplusone),
                    VB = structure(report$VB, names = Yearplusone),
                    R = structure(report$N[, 1], names = Yearplusone), N = structure(rowSums(report$N), names = Yearplusone),
                    N_at_age = report$N, Selectivity = report$vul, h = NA_real_,
                    Obs_Catch = structure(if (expanded) colSums(t(CAA_hist_VPA) * data$weight) else C_hist, names = Year),
                    Obs_Index = structure(I_hist, dimnames = list(Year, paste0("Index_", 1:nsurvey))),
                    Obs_C_at_age = CAA_hist_VPA,
                    Catch = structure(colSums(t(report$CAApred) * data$weight), names = Year),
                    Index = structure(report$Ipred, dimnames = list(Year, paste0("Index_", 1:nsurvey))), 
                    C_at_age = report$CAApred,
                    NLL = structure(c(report$nll, report$nll_comp, report$prior, report$penalty),
                                    names = c("Total", paste0("Index_", 1:nsurvey), "Prior", "Penalty")),
                    info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                    dependencies = dependencies)

  if (Assessment@conv) {
    ref_pt <- ref_pt_VPA(E = report$E[1:(length(report$E)-min_age)], R = report$N[(min_age + 1):length(report$E), 1], 
                         weight = info$data$weight, mat = info$LH$mat, M = info$data$M, vul = report$vul_p, 
                         SR = SR, fix_h = fix_h, h = ifelse(fix_h, Data@steep[x], NA_real_)) 
    report <- c(report, ref_pt[-17])

    #Z_mat <- t(report$F) + data$M
    #VB_mid <- t(report$N[-ncol(report$N), ]) * data$weight * (1 - exp(-Z_mat))/Z_mat

    Assessment@FMSY <- report$FMSY
    #Assessment@U <- structure(Assessment@Catch/colSums(VB_mid), names = Year)
    #Assessment@UMSY <- report$MSY/report$VBMSY
    #Assessment@U_UMSY <- Assessment@U/Assessment@UMSY
    Assessment@MSY <- report$MSY
    Assessment@BMSY <- report$BMSY
    Assessment@SSBMSY <- report$EMSY
    Assessment@VBMSY <- report$VBMSY
    Assessment@B0 <- report$B0
    Assessment@R0 <- report$R0
    Assessment@N0 <- report$N0
    Assessment@SSB0 <- report$E0
    Assessment@VB0 <- report$VB0
    Assessment@h <- report$h
    Assessment@F_FMSY <- Assessment@FMort/report$FMSY
    Assessment@B_BMSY <- structure(report$B/report$BMSY, names = Yearplusone)
    Assessment@B_B0 <- structure(report$B/report$B0, names = Yearplusone)
    Assessment@SSB_SSBMSY <- structure(report$E/report$EMSY, names = Yearplusone)
    Assessment@SSB_SSB0 <- structure(report$E/report$E0, names = Yearplusone)
    Assessment@VB_VBMSY <- structure(report$VB/report$VBMSY, names = Yearplusone)
    Assessment@VB_VB0 <- structure(report$VB/report$VB0, names = Yearplusone)
    Assessment@TMB_report <- report
    
    catch_eq <- function(Ftarget) {
      catch_equation(method = "Baranov", sel = report$vul_p, M = info$data$M,
                     wt = info$data$weight, N = report@N[nrow(report$N), ])
    }
    Assessment@forecast <- list(per_recruit = ref_pt[[17]], catch_eq = catch_eq)
  }

  return(Assessment)
}

VPA_posthoc <- function(report, info) { # Calculate Terminal year + 1 abundance
  n_age <- info$data$n_age
  Y <- nrow(report$N)
  N <- numeric(n_age)
  M <- info$data$M
  FF <- report$F_at_age[Y, ]
  Z <- M + FF
  
  refpt <- info$refpt
  
  N[2:n_age] <- report$N[Y, 1:(n_age-1)] * exp(-Z[1:(n_age-1)])
  N[n_age] <- N[n_age] + report$N[Y, n_age] * exp(-Z[n_age])
  N[1] <- quantile(report$N[(Y-refpt$R[2]+1):Y, 1], refpt$R[1])
  
  report$vul_p <- report$vul[(Y-refpt$vul+1):Y, , drop = FALSE] %>% apply(2, mean)
  report$vul_p <- report$vul_p/max(report$vul_p)
  
  report$N <- rbind(report$N, N)
  report$E <- colSums(t(report$N) * info$LH$mat * info$data$weight)
  report$VB <- c(report$VB, sum(N * report$vul_p * info$data$weight))
  report$B <- c(report$B, sum(N * info$data$weight))
  
  return(report)
}


ref_pt_VPA <- function(E, R, weight, mat, M, vul, SR, fix_h, h) {
  # Per-recruit quantities
  NPR0 <- calc_NPR(exp(-M), length(M))
  EPR0 <- sum(NPR0 * weight * mat)
  
  # Fit stock-recruit curve
  Rpred <- sigmaR <- NULL
  solve_SR_par <- function(x, h = NULL) {
    R0 <- exp(x[1])
    E0 <- R0 * EPR0
    if (!fix_h) {
      if (SR == "BH") h <- 0.2 + 0.8 * ilogit(x[2])
      if (SR == "Ricker") h <- 0.2 + exp(x[2])
    }
    if (SR == "BH") Rpred <<- (0.8 * R0 * h * E)/(0.2 * EPR0 * R0 *(1-h)+(h-0.2)*E)
    if (SR == "Ricker") Rpred <<- E/EPR0 * (5*h)^(1.25 * (1 - E/E0))
    sigmaR <<- sqrt(sum((log(R/Rpred))^2)/length(Rpred))
    nLL <- -sum(dnorm(log(R/Rpred), -0.5 * sigmaR^2, sigmaR, log = TRUE))
    return(nLL)
  }
  
  if (fix_h) {
    opt <- optimize(solve_SR_par, interval = c(-10, 10), h = h)$minimum
    R0 <- exp(opt)
  } else {
    opt <- nlminb(c(10, 10), solve_SR_par)
    R0 <- exp(opt$par[1])
    if (SR == "BH") h <- 0.2 + 0.8 * ilogit(opt$par[2])
    if (SR == "Ricker") h <- 0.2 + exp(opt$par[2])
  }
  
  # Virgin reference points
  N0 <- R0 * sum(NPR0)
  E0 <- R0 * EPR0
  VB0 <- R0 * sum(NPR0 * weight * vul)
  B0 <- R0 * sum(NPR0 * weight)
  
  # Steepness
  if (SR == "BH") {
    Arec <- 4*h/(1-h)/EPR0
    Brec <- (5*h-1)/(1-h)/E0
  }
  if (SR == "Ricker") {
    Arec <- 1/EPR0 * (5*h)^1.25
    Brec <- 1.25 * log(5*h) / E0
  }
  
  opt2 <- optimize(yield_fn_SCA, interval = c(1e-4, 4), M = M, mat = mat, weight = weight, vul = vul, 
                   SR = SR, Arec = Arec, Brec = Brec)
  opt3 <- yield_fn_SCA(opt2$minimum, M = M, mat = mat, weight = weight, vul = vul, SR = SR, 
                       Arec = Arec, Brec = Brec, opt = FALSE)
  
  FMSY <- opt2$minimum
  MSY <- -1 * opt2$objective
  VBMSY <- opt3["VB"]
  RMSY <- opt3["R"]
  BMSY <- opt3["B"]
  EMSY <- opt3["E"]
  
  Fvec <- seq(0, 2.5 * FMSY, length.out = 100)
  yield <- lapply(Fvec,
                  yield_fn_SCA, M = M, mat = mat, weight = weight, vul = vul, SR = SR, 
                  Arec = Arec, Brec = Brec, opt = FALSE)
  SPR <- vapply(yield, getElement, numeric(1), "EPR")
  YPR <- vapply(yield, getElement, numeric(1), "YPR")
  
  refpt_unfished <- list(h = h, Arec = Arec, Brec = Brec, E0 = E0, R0 = R0, N0 = N0, VB0 = VB0, B0 = B0, EPR0 = EPR0, NPR0 = NPR0)
  refpt_MSY <- list(FMSY = FMSY, MSY = MSY, VBMSY = VBMSY, RMSY = RMSY, BMSY = BMSY, EMSY = EMSY,
                    per_recruit = data.frame(FM = Fvec, SPR = SPR/SPR[1], YPR = YPR))
  return(c(refpt_unfished, refpt_MSY))
}
