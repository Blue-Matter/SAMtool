% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HCRs.R
\name{HCR_segment}
\alias{HCR_segment}
\title{Segmented harvest control rules}
\usage{
HCR_segment(
  Assessment,
  reps = 1,
  OCP_type = c("SSB_SSB0", "SSB_SSBMSY", "SSB_dSSB0", "F_FMSY", "F_F01", "F_FSPR"),
  Ftarget_type = c("FMSY", "F01", "Fmax", "FSPR"),
  OCP = c(0.1, 0.4),
  relF = c(0, 1),
  SPR_OCP,
  SPR_targ,
  ...
)
}
\arguments{
\item{Assessment}{An object of class \linkS4class{Assessment} with estimates of
FMSY or UMSY, vulnerable biomass, and spawning biomass depletion in terminal year.}

\item{reps}{The number of stochastic samples of the TAC recommendation.}

\item{OCP_type}{The type of operational control points (OCPs) for the harvest control rule used to determine the reduction in F.
By default, use (\code{"SSB_SSB0"} for spawning depletion. Other biomass OCPs include \code{"SSB_SSBMSY"} for spawning biomass relative to MSY and
\code{"SSB_dSSB0"}, for dynamic depletion (dynamic SSB0 is the historical reconstructed biomass with F = 0).
For F-based OCPs, the terminal year fishing mortality relative F01 or Fmax (using yield-per-recruit) or F-SPR\% (see \code{SPR_OCP} argument) can be used.}

\item{Ftarget_type}{The type of F used for the target fishing mortality rate.}

\item{OCP}{Numeric vector of operational control points for the HCR (in increasing order).}

\item{relF}{Numeric vector of Ftarget corresponding to the values in \code{OCP}.}

\item{SPR_OCP}{The value of spawning potential ratio for the OCP if \code{OCP_type = "F_FSPR"}. By default, 0.4 (F40\%).}

\item{SPR_targ}{The target value of spawning potential ratio if \code{Ftarget_type = "FSPR"}. By default, 0.4 (F40\%).}

\item{...}{Miscellaneous arguments.}
}
\value{
An object of class \linkS4class{Rec} with the TAC recommendation.
}
\description{
A linear segmented output control rule where the target F (used for the TAC recommendation) 
is a function of an operational control point (OCP) such as spawning depletion or spawning biomass relative to that at MSY. 
The joints of the HCR are specified by arguments \code{OCP} and \code{relF}. Beyond the range of \code{OCP}, the response will be flat.
\link{HCR_ramp} uses \code{HCR_segment} with two control points.
}
\examples{
# This is an MP with a 40-10 harvest control rule
DD_40_10 <- make_MP(DD_TMB, HCR_segment, OCP_type = "SSB_SSB0", OCP = c(0.1, 0.4), relF = c(0, 1))
}
\author{
Q. Huynh
}
