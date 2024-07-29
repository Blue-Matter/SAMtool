
#' Segmented harvest control rules
#' 
#' A linear segmented output control rule where the target F (used for the TAC recommendation) 
#' is a function of an operational control point (OCP) such as spawning depletion or spawning biomass. 
#' The segments of the HCR are specified by arguments `OCP` and `relF`. Beyond the range of `OCP`, the response will be flat.
#' [HCR_ramp] uses `HCR_segment` with two control points.
#' 
#' @param Assessment An object of class [Assessment-class] with estimates of
#' FMSY or UMSY, vulnerable biomass, and spawning biomass depletion in terminal year.
#' @param reps The number of stochastic samples of the TAC recommendation.
#' @param OCP_type The type of operational control points (OCPs) for the harvest control rule used to determine the reduction in F.
#' See below.
#' @param OCP Numeric vector of operational control points for the HCR (in increasing order).
#' @param Ftarget_type The type of F used for the target fishing mortality rate. See below.
#' @param relF Numeric vector of Ftarget corresponding to the values in `OCP`.
#' @param SPR_OCP The value of spawning potential ratio for the OCP if `OCP_type = "F_FSPR"`. By default, 0.4 (F40%).
#' @param SPR_targ The target value of spawning potential ratio if `Ftarget_type = "FSPR"`. By default, 0.4 (F40%). 
#' @param ... Miscellaneous arguments.
#' @return An object of class [MSEtool::Rec-class] with the TAC recommendation.
#' @author Q. Huynh
#' @details 
#' 
#' The catch advice is calculated using the catch equation of the corresponding
#' assessment. See `Assessment@forecast$catch_eq`, a function that returns the catch advice for a specified `Ftarget`.
#' 
#' *Operational control points (OCP_type)*
#' 
#' The following are the available options for harvest control rule inputs, and the source of those values
#' in the [Assessment-class] object:
#' 
#' \itemize{
#' \item **Default** `"SSB_SSB0"`: Spawning depletion. Uses the last value in  `Assessment@@SSB_SSB0` vector.
#' \item `"SSB_SSBMSY"`: Spawning biomass relative to MSY. Uses the last value in `Assessment@@SSB_SSBMSY` vector.
#' \item `"SSB_dSSB0"`: Dynamic depletion (SSB relative to the historical reconstructed biomass with F = 0). Uses the last value in
#' `Assessment@@SSB/Assessment@@TMB_report$dynamic_SSB0`.
#' \item `"F_FMSY"`: Fishing mortality relative to MSY. Uses the last value in `Assessment@@F_FMSY`.
#' \item `"F_F01"`: Fishing mortality relative to F_0.1 (yield per recruit), calculated from the data frame in
#' `Assessment@@forecast[["per_recruit"]]`.
#' \item `"F_FSPR"`: Fishing mortality relative to F_SPR% (the F that produces the spawning potential ratio specified in
#' `"SPR_OCP"`, calculated from the data frame in `Assessment@@forecast[["per_recruit"]]`.
#' }
#' 
#' *Fishing mortality target (Ftarget_type)*
#' 
#' The type of F for which the corresponding catch is calculated in the HCR is specified here. The source of those values
#' in the [Assessment-class] object is specified:
#' 
#' \itemize{
#' \item **Default** `"FMSY"`: Fishing mortality relative to MSY. Uses the value in `Assessment@@FMSY`.
#' \item `"F01"`: Fishing mortality relative to F_0.1 (yield per recruit), calculated from the data frame in
#' `Assessment@@forecast[["per_recruit"]]`.
#' \item `"Fmax"`: Fishing mortality relative to F_max (maximizing yield per recruit), calculated from the data frame in
#' `Assessment@@forecast[["per_recruit"]]`.
#' \item `"FSPR"`: Fishing mortality relative to F_SPR% (the F that produces the spawning potential ratio specified in
#' `"SPR_targ"`, calculated from data frame in `Assessment@@forecast[["per_recruit"]]`.
#' \item `"abs"`: Fishing mortality is independent of any model output and is explicitly specified in `relF`.
#' }
#' 
#' @examples 
#' # This is an MP with a 40-10 harvest control rule (using FMSY)
#' DD_40_10 <- make_MP(DD_TMB, HCR_segment, OCP_type = "SSB_SSB0", OCP = c(0.1, 0.4), relF = c(0, 1)) 
#' #' 
#' # This is an MP with a 40-10 harvest control rule with a maximum F of 0.1
#' DD_40_10 <- make_MP(DD_TMB, HCR_segment, OCP_type = "SSB_SSB0", 
#'                     Ftarget_type = "abs", OCP = c(0.1, 0.4), relF = c(0, 0.1)) 
#' @export
HCR_segment <- function(Assessment, reps = 1, OCP_type = c("SSB_SSB0", "SSB_SSBMSY", "SSB_dSSB0", "F_FMSY", "F_F01", "F_FSPR"),
                        Ftarget_type = c("FMSY", "F01", "Fmax", "FSPR", "abs"), 
                        OCP = c(0.1, 0.4), relF = c(0, 1), SPR_OCP, SPR_targ, ...) {
  dots <- list(...)
  OCP_type <- match.arg(OCP_type)
  Ftarget_type <- match.arg(Ftarget_type)
  
  n_OCP <- length(OCP)
  if (length(relF) != n_OCP) stop("Length of relF should be equal to length of OCP.")
  
  if (Assessment@conv) {
    
    if (OCP_type == "SSB_SSB0" && length(Assessment@SSB_SSB0)) {
      OCP_val <- Assessment@SSB_SSB0[length(Assessment@SSB_SSB0)]
    } else if (OCP_type == "SSB_SSBMSY" && length(Assessment@SSB_SSBMSY)) {
      OCP_val <- Assessment@SSB_SSBMSY[length(Assessment@SSB_SSBMSY)]
    } else if (OCP_type == "SSB_dSSB0" && !is.null(Assessment@TMB_report$dynamic_SSB0)) {
      OCP_val <- Assessment@SSB/Assessment@TMB_report$dynamic_SSB0
      OCP_val <- OCP_val[length(OCP_val)]
    } else if (OCP_type == "F_FMSY") {
      if (length(Assessment@U_UMSY)) {
        OCP_val <- Assessment@U_UMSY[length(Assessment@U_UMSY)]
      } else if (length(Assessment@F_FMSY)) {
        OCP_val <- Assessment@F_FMSY[length(Assessment@F_FMSY)]
      } else {
        OCP_val <- NA_real_
      }
    } else if (OCP_type == "F_F01" && length(Assessment@FMort)) {
      if (!is.null(Assessment@forecast$per_recruit$F01)) {
        F01 <- Assessment@forecast$per_recruit$F01
        OCP_val <- Assessment@FMort[length(Assessment@FMort)]/F01
      } else if (!is.null(Assessment@forecast$per_recruit$U)) {
        U01 <- get_F01(Assessment@forecast$per_recruit$U, Assessment@forecast$per_recruit$YPR)
        OCP_val <- Assessment@U[length(Assessment@U)]/U01
      } else {
        F01 <- get_F01(Assessment@forecast$per_recruit$FM, Assessment@forecast$per_recruit$YPR)
        OCP_val <- Assessment@FMort[length(Assessment@FMort)]/F01
      }
    } else if (OCP_type == "F_Fmax" && length(Assessment@FMort)) {
      if (!is.null(Assessment@forecast$per_recruit$Fmax)) {
        Fmax <- Assessment@forecast$per_recruit$Fmax
        OCP_val <- Assessment@FMort[length(Assessment@FMort)]/Fmax
      } else if (!is.null(Assessment@forecast$per_recruit$U)) {
        Umax <- get_Fmax(Assessment@forecast$per_recruit$U, Assessment@forecast$per_recruit$YPR)
        OCP_val <- Assessment@U[length(Assessment@U)]/Fmax
      } else {
        Fmax <- get_Fmax(Assessment@forecast$per_recruit$FM, Assessment@forecast$per_recruit$YPR)
        OCP_val <- Assessment@FMort[length(Assessment@FMort)]/Fmax
      }
    } else if (OCP_type == "F_FSPR" && length(Assessment@FMort)) {
      if (missing(SPR_OCP)) SPR_OCP <- 0.4
      if (!is.null(Assessment@forecast$per_recruit$U)) {
        U_SPR <- get_FSPR(Assessment@forecast$per_recruit$U, Assessment@forecast$per_recruit$SPR, target = SPR_OCP)
        OCP_val <- Assessment@U[length(Assessment@U)]/U_SPR
      } else {
        F_SPR <- get_FSPR(Assessment@forecast$per_recruit$FM, Assessment@forecast$per_recruit$SPR,
                          target = SPR_OCP)
        OCP_val <- Assessment@FMort[length(Assessment@FMort)]/F_SPR
      }
    } else {
      OCP_val <- NA_real_
    }
    
    if (!is.na(OCP_val) && OCP_val > 0) {
      alpha <- HCRlinesegment(OCP_val, OCP, relF)
      
      if (Ftarget_type == "FMSY") {
        if (length(Assessment@UMSY)) {
          Fout <- -log(1 - alpha * Assessment@UMSY)
          SE <- alpha * Assessment@SE_UMSY
        } else if (length(Assessment@FMSY)) {
          Fout <- alpha * Assessment@FMSY
          SE <- alpha * Assessment@SE_FMSY
        } 
      } else if (Ftarget_type == "F01") {
        if (!is.null(Assessment@forecast$per_recruit$F01)) {
          Fout <- alpha * Assessment@forecast$per_recruit$F01[1]
        } else if (!is.null(Assessment@forecast$per_recruit$U)) {
          U01 <- get_F01(Assessment@forecast$per_recruit$U, Assessment@forecast$per_recruit$YPR)
          Fout <- -log(1 - alpha * U01)
        } else {
          Fout <- alpha * get_F01(Assessment@forecast$per_recruit$FM, Assessment@forecast$per_recruit$YPR)
        }
      } else if (Ftarget_type == "Fmax") {
        if (!is.null(Assessment@forecast$per_recruit$Fmax)) {
          Fout <- alpha * Assessment@forecast$per_recruit$Fmax[1]
        } else if (!is.null(Assessment@forecast$per_recruit$U)) {
          Umax <- get_Fmax(Assessment@forecast$per_recruit$U, Assessment@forecast$per_recruit$YPR)
          Fout <- -log(1 - alpha * Umax)
        } else {
          Fout <- alpha * get_Fmax(Assessment@forecast$per_recruit$FM, Assessment@forecast$per_recruit$YPR)
        }
      } else if (Ftarget_type == "FSPR") {
        if (missing(SPR_targ)) {
          if (!is.null(dots$SPR)) {
            SPR_targ <- dots$SPR
          } else SPR_targ <- 0.4
        }
        if (!is.null(Assessment@forecast$per_recruit$U)) {
          U_SPR <- get_FSPR(Assessment@forecast$per_recruit$U, Assessment@forecast$per_recruit$SPR, target = SPR_targ)
          Fout <- -log(1 - alpha * U_SPR)
        } else {
          Fout <- alpha * get_FSPR(Assessment@forecast$per_recruit$FM, Assessment@forecast$per_recruit$SPR,
                                   target = SPR_targ)
        }
      } else if (Ftarget_type == "abs") {
        Fout <- alpha
      }
      
      if (exists("Fout", inherits = FALSE) && !is.na(Fout)) {
        if (Fout > 0) {
          if (!exists("SE", inherits = FALSE) || !length(SE)) SE <- 0
          FM <- trlnorm(reps, Fout, SE/Fout)
        } else {
          FM <- rep(0, reps)
        }
        TAC <- calculate_TAC(Assessment, Ftarget = FM)
      }
    }
  }
  if (!exists("TAC", inherits = FALSE)) TAC <- rep(NA_real_, reps)
  
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  return(Rec)
}
class(HCR_segment) <- "HCR"

#' Linearly ramped harvest control rules
#'
#' An output control rule with a ramp that reduces the target F (used for the TAC recommendation) linearly
#' as a function of an operational control point (OCP) such as spawning depletion or spawning biomass. The reduction in F is linear when the OCP
#' is between the target OCP (TOCP) and the limit OCP (LOCP). The target F is maximized at or above the TOCP. Below the LOCP,
#' the target F is minimized. For example, the TOCP and LOCP for 40% and 10% spawning depletion, respectively, in the 40-10 control rule.
#' Ftarget is FMSY above the TOCP and zero below the LOCP. This type of control rule can generalized with more control points (>2) in [HCR_segment].
#' Class HCR objects are typically used with function [make_MP].
#'
#' @param LOCP Numeric, the limit value for the OCP in the HCR.
#' @param TOCP Numeric, the target value for the OCP in the HCR.
#' @param relF_min The relative value of Ftarget (i.e., as a proportion) if `OCP < LOCP`.
#' @param relF_max The relative value of Ftarget if `OCP > TOCP`.
#' @inheritParams HCR_segment
#' @inherit HCR_segment details
#' @describeIn HCR_ramp Generic ramped-HCR function where user specifies OCP and corresponding limit and target
#' points, as well as minimum and maximum relative F target.
#' @return An object of class [MSEtool::Rec-class] with the TAC recommendation.
#' @author Q. Huynh & T. Carruthers
#' @references
#' Deroba, J.J. and Bence, J.R. 2008. A review of harvest policies: Understanding relative
#' performance of control rules. Fisheries Research 94:210-223.
#'
#' Edwards, C.T.T. and Dankel, D.J. (eds.). 2016. Management Science in Fisheries: an introduction
#' to simulation methods. Routledge, New York, NY. 460 pp.
#'
#' Punt, A. E, Dorn, M. W., and Haltuch, M. A. 2008. Evaluation of threshold management strategies
#' for groundfish off the U.S. West Coast. Fisheries Research 94:251-266.
#'
#' Restrepo, V.R. and Power, J.E. 1999. Precautionary control rules in US fisheries
#' management: specification and performance. ICES Journal of Marine Science 56:846-852.
#' @seealso [HCR_segment] [HCR_MSY] [HCRlin] [make_MP]
#' @examples
#' # 40-10 linear ramp
#' Brel <- seq(0, 1, length.out = 200)
#' plot(Brel, HCRlin(Brel, 0.1, 0.4), 
#'     xlab = expression("Operational control point: Estimated"~SSB/SSB[0]),
#'     ylab = expression(F[target]~~": proportion of"~~F[MSY]), 
#'     main = "40-10 harvest control rule", type = "l")
#' abline(v = c(0.1, 0.4), col = "red", lty = 2)
#'
#' # create a 40-10 MP to run in closed-loop MSE
#' DD_40_10 <- make_MP(DD_TMB, HCR40_10)
#'
#' # Alternatively,
#' DD_40_10 <- make_MP(DD_TMB, HCR_ramp, OCP_type = "SSB_SSB0", LOCP = 0.1, TOCP = 0.4)
#'
#' # An SCA with LOCP and TOCP at 0.4 and 0.8, respectively, of SSB/SSBMSY
#' SCA_80_40 <- make_MP(SCA, HCR_ramp, OCP_type = "SSB_SSBMSY", LOCP = 0.4, TOCP = 0.8)
#'
#' # A conservative HCR that fishes at 75% of FMSY at B > 80% BMSY but only reduces F
#' # to 10% of FMSY if B < 40% BMSY.
#' SCA_conservative <- make_MP(SCA, HCR_ramp, OCP_type = "SSB_SSBMSY", LOCP = 0.4, TOCP = 0.8, 
#' relF_min = 0.1, relF_max = 0.75)
#'
#' # Figure of this conservative HCR
#' Brel <- seq(0, 1, length.out = 200)
#' Frel <- HCRlin(Brel, 0.4, 0.8, 0.1, 0.75)
#' plot(Brel, Frel, 
#'     xlab = expression("Operational control point: Estimated"~SSB/SSB[MSY]),
#'     ylab = expression(F[target]~":"~~F/F[MSY]), 
#'     ylim = c(0, 1), type = "l")
#' abline(v = c(0.4, 0.8), col = "red", lty = 2)
#' 
#' # A harvest control rule as a function of BMSY, with F independent of model output, 
#' # i.e., specify F in relF argument (here maximum F of 0.1)
#' SCA_80_40 <- make_MP(SCA, HCR_ramp, OCP_type = "SSB_SSBMSY", LOCP = 0.4, TOCP = 0.8, 
#'                      relF_min = 0, relF_max = 0.1)
#'
#' @export
HCR_ramp <- function(Assessment, reps = 1, OCP_type = c("SSB_SSB0", "SSB_SSBMSY", "SSB_dSSB0", "F_FMSY", "F_F01", "F_FSPR"),
                     Ftarget_type = c("FMSY", "F01", "Fmax", "FSPR", "abs"), 
                     LOCP = 0.1, TOCP = 0.4, relF_min = 0, relF_max = 1, SPR_OCP = 0.4, SPR_targ = 0.4, ...) {
  HCR_segment(Assessment = Assessment, reps = reps, OCP_type = OCP_type, Ftarget_type = Ftarget_type,
              OCP = c(LOCP, TOCP), relF = c(relF_min, relF_max), SPR_OCP = SPR_OCP, SPR_targ = SPR_targ, ...)
}
class(HCR_ramp) <- "HCR"


#' @describeIn HCR_ramp Common U.S. west coast control rule (LOCP and TOCP of 0.1 and 0.4 spawning depletion,
#' respectively)
#' @export
HCR40_10 <- function(Assessment, reps = 1, Ftarget_type = "FMSY", SPR_targ = 0.4, ...) {
  HCR_segment(Assessment, reps, OCP = c(0.1, 0.4), Ftarget_type = Ftarget_type,
              relF = c(0, 1), SPR_targ = SPR_targ, ...)
}
class(HCR40_10) <- "HCR"



#' @describeIn HCR_ramp More conservative than `HCR40_10`, with LOCP and TOCP of 0.2 and 0.6
#' spawning depletion, respectively).
#' @export
HCR60_20 <- function(Assessment, reps = 1, Ftarget_type = "FMSY", SPR_targ = 0.4, ...) {
  HCR_segment(Assessment, reps, OCP = c(0.2, 0.6), Ftarget_type = Ftarget_type, 
              relF = c(0, 1), SPR_targ = SPR_targ, ...)
}
class(HCR60_20) <- "HCR"

#' @describeIn HCR_ramp 0.8 and 0.4 SSBMSY as the LOCP and TOCP, respectively. 
#' @export
HCR80_40MSY <- function(Assessment, reps = 1, Ftarget_type = "FMSY", SPR_targ = 0.4, ...) {
  HCR_segment(Assessment, reps, OCP_type = "SSB_SSBMSY", OCP = c(0.4, 0.8),
              Ftarget_type = Ftarget_type, relF = c(0, 1), SPR_targ = SPR_targ, ...)
}
class(HCR80_40MSY) <- "HCR"

#' Harvest control rule to fish at some fraction of maximum sustainable yield
#'
#' A simple control rule that specifies the total allowable catch (TAC) as a function of the abundance of the first 
#' projection year and some fraction of FMSY/UMSY. 
#'
#' @param Assessment An object of class [Assessment-class] with estimates of
#' FMSY or UMSY and vulnerable biomass in terminal year.
#' @param reps The number of stochastic samples of the TAC recommendation.
#' @param MSY_frac The fraction of FMSY or UMSY for calculating the TAC (e.g. MSY_frac = 0.75 fishes at 75% of FMSY).
#' @param ... Miscellaneous arguments.
#' @details 
#' The catch advice is calculated using the catch equation of the corresponding
#' assessment. See `Assessment@forecast$catch_eq`, a function that returns the catch advice for a specified `Ftarget`.
#' @return An object of class [MSEtool::Rec-class] with the TAC recommendation.
#' @author Q. Huynh
#' @references
#' Punt, A. E, Dorn, M. W., and Haltuch, M. A. 2008. Evaluation of threshold management strategies
#' for groundfish off the U.S. West Coast. Fisheries Research 94:251-266.
#' @seealso [make_MP] [HCR_ramp]
#' @examples
#' # create an MP to run in closed-loop MSE (fishes at UMSY)
#' SPMSY <- make_MP(SP, HCR_MSY)
#'
#' # The MP which fishes at 75% of FMSY
#' SP75MSY <- make_MP(SP, HCR_MSY, MSY_frac = 0.75)
#'
#' \donttest{
#' myOM <- MSEtool::runMSE(MSEtool::testOM, MPs = c("FMSYref", "SPMSY", "SP75MSY"))
#' }
#' @export
HCR_MSY <- function(Assessment, reps = 1, MSY_frac = 1, ...) {
  HCR_segment(Assessment = Assessment, reps = reps, OCP = 0, relF = MSY_frac)
}
class(HCR_MSY) <- "HCR"

#' Fixed escapement harvest control rule
#'
#' A simple control rule that allows fishing when the operational control point (OCP) is above some threshold.
#' By default, this function sets the TAC at F = 100% FMSY when spawning depletion > 0.1.
#'
#' @param Assessment An object of class [Assessment-class] with estimates of
#' FMSY or UMSY and vulnerable biomass in terminal year.
#' @param reps The number of stochastic samples of the TAC recommendation.
#' @param OCP_type The type of operational control points (OCPs) for the harvest control rule used to determine 
#' whether there is fishing. By default, use (`"SSB_SSB0"` for spawning depletion. Other biomass OCPs include `"SSB_SSBMSY"` for spawning biomass relative to MSY and
#' `"SSB_dSSB0"`, for dynamic depletion (dynamic SSB0 is the historical reconstructed biomass with F = 0).
#' For F-based OCPs, the terminal year fishing mortality relative F01 or Fmax (using yield-per-recruit) or F-SPR% (see `SPR_OCP` argument) can be used.
#' @param OCP_threshold The value of the OCP above which fishing can occur.
#' @param Ftarget_type The type of F used for the target fishing mortality rate.
#' @param relF_max The relative value of Ftarget if `OCP > OCP_treshold`.
#' @param ... Miscellaneous arguments.
#' @return An object of class [MSEtool::Rec-class] with the TAC recommendation.
#' @author Q. Huynh
#' @inherit HCR_MSY details 
#' @references
#' Deroba, J.J. and Bence, J.R. 2008. A review of harvest policies: Understanding relative
#' performance of control rules. Fisheries Research 94:210-223.
#' @seealso [make_MP] [HCR_ramp]
#' @examples
#' # create an MP to run in closed-loop MSE (fishes at FMSY when B/B0 > 0.2)
#' SP_escapement <- make_MP(SP, HCR_escapement)
#'
#' # The MP which fishes at 75% of FMSY
#' SP_escapement75 <- make_MP(SP, HCR_escapement, relF_max = 0.75)
#' 
#' # The MP which fishes at FMSY when BMSY > 0.5
#' SP_BMSY_escapement <- make_MP(SP, HCR_escapement, OCP_type = "SSB_SSBMSY", 
#'                               OCP_threshold = 0.5, relF_max = 1)
#'
#' \donttest{
#' myOM <- MSEtool::runMSE(MSEtool::testOM, MPs = c("FMSYref", "SP_escapement", "SP_BMSY_escapement"))
#' }
#' @export
HCR_escapement <- function(Assessment, reps = 1, OCP_type = "SSB_SSB0", OCP_threshold = 0.2, 
                           Ftarget_type = "FMSY", 
                           relF_max = 1, ...) {
  HCR_segment(Assessment, reps, OCP_type = OCP_type, OCP = rep(OCP_threshold, 2), 
              Ftarget_type = Ftarget_type, relF = c(0, relF_max), ...)
}
class(HCR_escapement) <- "HCR"

#' Simple fixed F harvest control rule
#'
#' A simple control rule that explicitly specifies the target apical F independent of any model.
#'
#' @param Assessment An object of class [Assessment-class] with estimates of next year's abundance or biomass.
#' @param reps The number of replicates of the TAC recommendation (not used).
#' @param Ftarget The value of F.
#' @return An object of class [MSEtool::Rec-class] with the TAC recommendation.
#' @author Q. Huynh
#' @seealso [make_MP] [HCR_ramp]#' 
#' @inherit HCR_MSY details 
#' @examples
#' # create an MP to run in closed-loop MSE (fishes at F = 0.2)
#' F0.2 <- make_MP(SP, HCR_fixedF, Ftarget = 0.2)
#' 
#' \donttest{
#' myOM <- MSEtool::runMSE(MSEtool::testOM, MPs = c("FMSYref", "F0.2"))
#' }
#' @export
HCR_fixedF <- function(Assessment, reps = 1, Ftarget = 0.1) {
  if (Assessment@conv) {
    TAC <- calculate_TAC(Assessment, Ftarget = Ftarget)
  } else {
    TAC <- NA_real_
  }
  Rec <- new("Rec")
  Rec@TAC <- rep(TAC, reps)
  return(Rec)
}
class(HCR_fixedF) <- "HCR"

#' Generic linear harvest control rule based on biomass
#'
#' A general function used by [HCR_ramp] that adjusts the output (e.g., F) by a linear ramp based on 
#' the value of the OCP relative to target and limit values.
#'
#' @param OCP_val The value of the operational control point (OCP).
#' @param LOCP Numeric, the limit value for the OCP in the HCR.
#' @param TOCP Numeric, the target value for the OCP in the HCR.
#' @param relF_min The relative maximum value (e.g. a multiple of FMSY) if `OCP < LOCP`.
#' @param relF_max The relative maximum value (e.g. a multiple of FMSY) if `OCP > TOCP`.
#' @return Numeric adjustment factor.
#' @author T. Carruthers
#' @export HCRlin
#' @examples
#' #40-10 linear ramp
#' Brel <- seq(0, 1, length.out = 200)
#' plot(Brel, HCRlin(Brel, 0.1, 0.4), xlab = "Estimated B/B0", ylab = "Relative change in F",
#'      main = "A 40-10 harvest control rule", type = 'l', col = 'blue')
#' abline(v = c(0.1,0.4), col = 'red', lty = 2)
HCRlin <- function(OCP_val, LOCP, TOCP, relF_min = 0, relF_max = 1){
  adj <- rep(relF_max, length(OCP_val))
  adj[OCP_val <= LOCP] <- relF_min
  cond <- OCP_val > LOCP & OCP_val < TOCP
  adj[cond] <- (relF_max - relF_min)/(TOCP - LOCP) * (OCP_val[cond] - LOCP) + relF_min
  return(adj)
}

HCRlinesegment <- function(OCP_val, OCP, relF) {
  OCP <- c(0, OCP, Inf)
  relF <- c(relF[1], relF, relF[length(relF)])
  OCP_ind <- findInterval(OCP_val, OCP)
  slope <- (relF[OCP_ind+1] - relF[OCP_ind])/(OCP[OCP_ind+1] - OCP[OCP_ind]) 
  x_diff <- OCP_val - OCP[OCP_ind]
  slope * x_diff + relF[OCP_ind]
}

#' A Harvest Control Rule using B/BMSY and F/FMSY to adjust TAC or TAE.
#'
#' @param Brel improper fraction: an estimate of Biomass relative to BMSY
#' @param Frel improper fraction: an estimate of Fishing mortality rate relative to FMSY
#' @param Bpow non-negative real number: controls the shape of the biomass adjustment, when zero there is no adjustment
#' @param Bgrad non-negative real number: controls the gradient of the biomass adjustment
#' @param Fpow non-negative real number: controls the adjustment speed relative to F/FMSY. When set to 1, next recommendation is FMSY. When less than 1 next recommendation is between current F and FMSY.
#' @param Fgrad improper fraction: target Fishing rate relative to FMSY
#' @return a TAC or TAE adjustment factor.
#' @author T. Carruthers
#' @references Made up for this package
#' @export
#' @examples
#' res <- 100
#' Frel <- seq(1/2, 2, length.out = res)
#' Brel <- seq(0.05, 2, length.out=res)
#' adj <- array(HCR_FB(Brel[rep(1:res, res)], Frel[rep(1:res, each = res)],
#'                     Bpow = 2, Bgrad = 1, Fpow = 1, Fgrad = 0.75), c(res, res))
#' contour(Brel, Frel, adj, nlevels = 20, xlab = "B/BMSY", ylab = "F/FMSY",
#'         main = "FBsurface TAC adjustment factor")
#' abline(h = 1, col = 'red', lty = 2)
#' abline(v = 1, col = 'red', lty = 2)
#' legend('topright', c("Bpow = 2", "Bgrad = 1", "Fpow = 1", "Fgrad = 0.75"), text.col = 'blue')
HCR_FB<-function(Brel,Frel,Bpow=2,Bgrad=1,Fpow=1,Fgrad=1){
  Fresp <- exp(log(1/Frel)*Fpow)
  Bresp <- exp(powdif(Brel-1,Bpow,Bgrad))
  Fgrad*Bresp*Fresp
}


# #' Power difference function.
# #'
# #' @param x Real number: the absolute difference between two numbers
# #' @param z Real number: the exponent the difference will be raised to
# #' @param g Real number: the gradient in the exponential difference
# #' @return a positive real number
# #' @export
# #' @examples
# #' powdif(-2,3,1)
powdif<-function(x,z,g){
  x2<-(g*abs(x))^z
  x2[x<0]<-(-x2[x<0])
  x2
}

#' Calculate MSY-based TAC from Assessment object
#'
#' A function to calculate the total allowable catch (TAC). Based on the MSY (maximum
#' sustainable yield) principle, the TAC is the product of
#' either UMSY or FMSY and the available biomass, i.e. vulnerable biomass, in terminal year.
#'
#' @param Assessment An Assessment object with estimates of UMSY or FMSY and
#' terminal year vulnerable biomass.
#' @param reps The number of stochastic draws of UMSY or FMSY.
#' @param MSY_frac The fraction of FMSY or UMSY for calculating the TAC (e.g. MSY_frac = 0.75 fishes at 75% of FMSY).
#' @note `calculate_TAC` is deprecated as of version 1.2 in favor of `TAC_MSY` because
#' the latter has a more informative name.
#' @return A vector of length `reps` of stochastic samples of TAC recommendation. Returns NA's
#' if missing either UMSY/FMSY or vulnerable biomass.
#' @seealso [HCR_MSY] [HCR40_10] [HCR60_20]
#' @aliases calculate_TAC
#' @export
TAC_MSY <- function(Assessment, reps, MSY_frac = 1) {
  if (length(Assessment@UMSY)) {
    Fout <- -log(1 - MSY_frac * Assessment@UMSY)
    SE <- MSY_frac * Assessment@SE_UMSY
  } else if (length(Assessment@FMSY)) {
    Fout <- MSY_frac * Assessment@FMSY
    SE <- MSY_frac * Assessment@SE_FMSY
  } else {
    Fout <- SE <- numeric(0)
  }
  
  if (length(Fout) && !is.na(Fout) && Fout > 0) {
    if (length(SE)) {
      FM <- trlnorm(reps, Fout, SE/Fout)
      TAC <- calculate_TAC(Assessment, Ftarget = FM)
    } else {
      TAC <- rep(NA_real_, reps)
    }
  } else {
    TAC <- rep(NA_real_, reps)
  }
  return(TAC)
}

calculate_TAC <- function(Assessment, Ftarget, Utarget) { # Vectorized for Ftarget or Utarget
  if (Assessment@conv) {
    TAC <- try(vapply(Ftarget, Assessment@forecast$catch_eq, numeric(1)), silent = TRUE)
    if (is.character(TAC)) {
      TAC <- try(mapply(catch_equation, Utarget = 1 - exp(-Ftarget), 
                        MoreArgs = list(method = "frac", B = Assessment@VB[length(Assessment@VB)])),
                 silent = TRUE)
    }
  } 
  if (!exists("TAC", inherits = FALSE) || is.character(TAC)) TAC <- rep(NA_real_, length(Ftarget))
  return(TAC)
}


catch_equation <- function(method = c("frac", "Baranov", "cDD", "SP"), ...) {
  method <- match.arg(method)
  dots <- list(...)
  
  if (method == "frac") {
    args <- dots_check(c("Utarget", "B"), dots)
    catch <- args$Utarget * args$B
    
  } else if (method == "Baranov") {
    args <- dots_check(c("sel", "Ftarget", "M", "wt", "N"), dots)
    catch <- SCA_catch_solver(FM = args$Ftarget, N = args$N, weight = args$wt, vul = args$sel, M = args$M)$Cpred
  
  } else if (method == "cDD") {
    args <- dots_check(c("Ftarget", "B", "N", "R", "M", "Kappa", "Winf", "wk"), dots)
    catch <- cDD_catch_solver(FM = args$Ftarget, B = args$B, N = args$N, R = args$R, M = args$M, Kappa = args$Kappa,
                              Winf = args$Winf, wk = args$wk)[1]
    
  } else if (method == "SP") {
    args <- dots_check(c("Ftarget", "B", "MSY", "K", "n"), dots)
    if (is.null(dots$n_seas)) {
      dt <- 1
    } else {
      dt <- 1/dots$n_seas
    }
    n_term <- ifelse(args$n == 1, exp(1), args$n^(args$n/(args$n-1)))
    
    catch <- SP_catch_solver(FM = args$Ftarget, B = args$B, dt = dt, MSY = args$MSY, 
                             K = args$K, n = args$n, n_term = n_term)[1]
  }
  if (!exists("catch", inherits = FALSE)) catch <- NA_real_
  return(catch)
}

dots_check <- function(vars, dots) {
  out <- lapply(vars, function(x) getElement(dots, x)) %>% structure(names = vars)
  check <- vapply(out, is.null, logical(1))
  if (any(check)) {
    stop(paste0(paste(vars[check], collapse = ", "), " was not found in call to catch_equation()."),
         call.= FALSE)
  }
  return(out)
}
