#' Harvest control rule to fish at some fraction of maximum sustainable yield
#'
#' A simple control rule that specifies the total allowable catch (TAC) to be the
#' product of current vulnerable biomass and UMSY.
#'
#' @param Assessment An object of class \linkS4class{Assessment} with estimates of
#' FMSY or UMSY and vulnerable biomass in terminal year.
#' @param reps The number of stochastic samples of the TAC recommendation.
#' @param MSY_frac The fraction of FMSY or UMSY for calculating the TAC (e.g. MSY_frac = 0.75 fishes at 75\% of FMSY).
#' @param ... Miscellaneous arguments.
#' @return An object of class \linkS4class{Rec} with the TAC recommendation.
#' @author Q. Huynh
#' @references
#' Punt, A. E, Dorn, M. W., and Haltuch, M. A. 2008. Evaluation of threshold management strategies
#' for groundfish off the U.S. West Coast. Fisheries Research 94:251-266.
#' @seealso \link{make_MP} \link{HCR_ramp}
#' @examples
#' # create an MP to run in closed-loop MSE (fishes at UMSY)
#' DD_MSY <- make_MP(DD_TMB, HCR_MSY)
#' class(DD_MSY)
#'
#' # The same MP which fishes at 75% of UMSY
#' DD_75MSY <- make_MP(DD_TMB, HCR_MSY, MSY_frac = 0.75)
#' class(DD_75MSY)
#'
#' \dontrun{
#' myOM <- MSEtool::runMSE(MSEtool::testOM, MPs = c("FMSYref", "DD_MSY", "DD_75MSY"))
#' }
#' @export
HCR_MSY <- function(Assessment, reps = 1, MSY_frac = 1, ...) {
  TAC <- TAC_MSY(Assessment, reps, MSY_frac)
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  return(Rec)
}
class(HCR_MSY) <- "HCR"

#' Linearly ramped harvest control rules
#'
#' An output control rule with a ramp that reduces the target F (used for the TAC recommendation) linearly
#' as a function of an operational control point (OCP) such as spawning depletion or spawning biomass relative to that at MSY. The reduction in F is linear when the OCP
#' is between the target OCP (TOCP) and the limit OCP (LOCP). Above the TOCP, the target F is maximized. Below the LOCP,
#' the target F is minimized. For example, the TOCP and LOCP for 40% and 10% spawning depletion, respectively, in the 40-10 control rule.
#' The Ftarget is FMSY above the TOCP and zero below the LOCP.
#' Class HCR objects are typically used with function \link{make_MP}.
#'
#' @param Assessment An object of class \linkS4class{Assessment} with estimates of
#' FMSY or UMSY, vulnerable biomass, and spawning biomass depletion in terminal year.
#' @param reps The number of stochastic samples of the TAC recommendation.
#' @param OCP_type The type of operational control points (OCPs) for the harvest control rule used to determine the reduction in F.
#' By default, use (\code{"SSB_SSB0"} for spawning depletion. Otherwise use \code{"SSB_SSBMSY"} for spawning biomass relative to MSY).
#' @param LOCP Numeric, the limit value for the OCP in the HCR.
#' @param TOCP Numeric, the target value for the OCP in the HCR.
#' @param Ftarget_type The type of F used for the target fishing mortality rate. Currently only FMSY is supported.
#' @param relF_min The relative value of Ftarget (i.e., as a proportion) if \code{OCP < LOCP}.
#' @param relF_max The relative value of Ftarget if \code{OCP > TOCP}.
#' @param ... Miscellaneous arguments.
#' @details \code{HCR_ramp} is the generic ramped-HCR function where user specifies OCP and corresponding limit and target
#' pointss, as well as minimum and maximum relative F target.
#'
#' \code{HCR40_10} is a common U.S. west coast control rule (LOCP and TOCP of 0.1 and 0.4 spawning depletion,
#' respectively), while \code{HCR60_20} is more conservative than 40-10, with LOCP and TOCP of 0.2 and 0.6
#' spawning depletion, respectively).
#' @return An object of class \linkS4class{Rec} with the TAC recommendation.
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
#' @seealso \link{HCR_MSY} \link{HCRlin} \link{make_MP}
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
#' \dontrun{
#' myOM <- MSEtool::runMSE(MSEtool::testOM, MPs = c("FMSYref", "DD_40_10"))
#' }
#' @export
HCR_ramp <- function(Assessment, reps = 1, OCP_type = c("SSB_SSB0", "SSB_SSBMSY"),
                     Ftarget_type = "MSY", LOCP = 0.1, TOCP = 0.4, relF_min = 0, relF_max = 1, ...) {
  OCP_type <- match.arg(OCP_type)

  if(OCP_type == "SSB_SSB0" && length(Assessment@SSB_SSB0) > 0) {
    relB <- Assessment@SSB_SSB0[length(Assessment@SSB_SSB0)]
  } else if(OCP_type == "SSB_SSBMSY" && length(Assessment@SSB_SSBMSY) > 0) {
    relB <- Assessment@SSB_SSBMSY[length(Assessment@SSB_SSBMSY)]
  } else relB <- NA_real_

  if(!is.na(relB)) {
    alpha <- HCRlin(relB, LOCP, TOCP, relF_min, relF_max)
    TAC <- TAC_MSY(Assessment, reps, MSY_frac = alpha)
  } else TAC <- rep(NA_real_, reps)

  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  return(Rec)
}
class(HCR_ramp) <- "HCR"


#' @rdname HCR_ramp
#' @export
HCR40_10 <- function(Assessment, reps = 1, ...) HCR_ramp(Assessment, reps, LOCP = 0.1, TOCP = 0.4)
class(HCR40_10) <- "HCR"


#' @rdname HCR_ramp
#' @export
HCR60_20 <- function(Assessment, reps = 1, ...) HCR_ramp(Assessment, reps, LOCP = 0.2, TOCP = 0.6)
class(HCR60_20) <- "HCR"



#' Generic linear harvest control rule based on biomass
#'
#' A general function used by \link{HCR_ramp} that adjusts the output (e.g., F) by a linear ramp based on 
#' the value of the OCP relative to target and limit values.
#'
#' @param value The value of the operational control point (OCP).
#' @param LOCP Numeric, the limit value for the OCP in the HCR.
#' @param TOCP Numeric, the target value for the OCP in the HCR.
#' @param relF_min The relative maximum value (e.g. a multiple of FMSY) if \code{OCP < LOCP}.
#' @param relF_max The relative maximum value (e.g. a multiple of FMSY) if \code{OCP > TOCP}.
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
  adj[cond] <- (relF_max - relF_min)/(TOCP - LOCP) * (adj[cond] - LOCP) + relF_min
  return(adj)
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
#' @param MSY_frac The fraction of FMSY or UMSY for calculating the TAC (e.g. MSY_frac = 0.75 fishes at 75\% of FMSY).
#' @note \code{calculate_TAC} is deprecated as of version 1.2 in favor of \code{TAC_MSY} because
#' the latter has a more informative name.
#' @return A vector of length \code{reps} of stochastic samples of TAC recommendation. Returns NA's
#' if missing either UMSY/FMSY or vulnerable biomass.
#' @seealso \link{HCR_MSY} \link{HCR40_10} \link{HCR60_20}
#' @aliases calculate_TAC
#' @export
TAC_MSY <- function(Assessment, reps, MSY_frac = 1) {
  has_UMSY <- length(Assessment@UMSY) > 0
  has_FMSY <- length(Assessment@FMSY) > 0
  has_VB <- length(Assessment@VB) > 0

  if(Assessment@conv && has_VB && (has_UMSY || has_FMSY)) {
    VB_current <- Assessment@VB[length(Assessment@VB)]
    if(has_UMSY) {
      if(length(Assessment@SE_UMSY) > 0) {
        SE_UMSY <- Assessment@SE_UMSY
      } else SE_UMSY <- 1e-8
      UMSY_vector <- trlnorm(reps, Assessment@UMSY, SE_UMSY)
      TAC <- MSY_frac * UMSY_vector * VB_current
    }
    if(has_FMSY) {
      if(length(Assessment@SE_FMSY) > 0) {
        SE_FMSY <- Assessment@SE_FMSY
      } else SE_FMSY <- 1e-8
      FMSY_vector <- trlnorm(reps, Assessment@FMSY, SE_FMSY)
      TAC <- (1 - exp(-MSY_frac * FMSY_vector)) * VB_current
    }
  } else {
    TAC <- rep(NA_real_, reps) # Missing estimates for HCR.
  }
  return(as.numeric(TAC))
}

