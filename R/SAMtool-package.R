#' Stock Assessment Methods Toolkit 
#'
#' Simulation tools for closed-loop simulation are provided for the 'OMtool' operating model to inform data-rich fisheries. 
#' SAMtool provides an OM conditioning model, assessment models of varying complexity with standardized reporting, diagnostic tools for evaluating 
#' assessments within closed-loop simulation, and helper functions for building more complex operating models and model-based management procedures.
#'
#' @name SAMtool-package
#' @aliases SAMtool
#' @docType package
#' @author Quang Huynh \email{quang@bluematterscience.com}
#' @author Tom Carruthers \email{tom@bluematterscience.com}
#' @author Adrian Hordyk \email{adrian@bluematterscience.com}
#' @section How to use SAMtool:
#' The main features of SAMtool are the assessment models and the ability to make model-based management procedures by combining
#' assessment models with harvest control rules. Such MPs can be used and tested in management strategy evaluation
#' with OMtool operating models. An overview of these features is available in the \href{../doc/SAMtool.html}{SAMtool vignette}.
#'
#' The following assessment models are available:
#' \itemize{
#' \item \href{../doc/Surplus_production.html}{Surplus production} (\link{SP}, \link{SP_SS}, \link{SP_Fox}, and \code{spict})
#' \item \href{../doc/Delay_difference.html}{Delay difference} (\link{DD}, \link{cDD}, \link{DD_SS}, and \link{cDD_SS})
#' \item \href{../doc/SCA.html}{Statistical catch-at-age} (\link{SCA}, \link{SCA2}, and \link{SCA_Pope})
#' \item Simple Stock Synthesis (\link{SSS} which implements \link{SCA_Pope} with fixed depletion assumption)
#' \item \href{../doc/VPA.html}{Virtual population analysis} (\link{VPA})
#' }
#'
#' The \code{\link{RCM}} (Rapid Conditioning Model) can be used to condition operating models from real data. Information can be found
#' \href{../doc/RCM.html}{here}.
#'
#' All SAMtool vignettes can also be viewed by typing \code{browseVignettes("SAMtool")} into the R console or through the
#' SAMtool webpage on \href{https://cran.r-project.org/package=SAMtool}{CRAN}.
#'
#' @references
#' Carruthers, T.R., Punt, A.E., Walters, C.J., MacCall, A.,
#' McAllister, M.K., Dick, E.J., Cope, J. 2014. Evaluating methods for setting
#' catch limits in data-limited fisheries. Fisheries Research. 153: 48-68.
#'
#' Carruthers, T.R., Kell, L.T., Butterworth, D.S., Maunder, M.N., Geromont,
#' H.F., Walters, C., McAllister, M.K., Hillary, R., Levontin, P., Kitakado,
#' T., Davies, C.R. Performance review of simple management procedures. ICES
#' Journal of Marine Science. 73: 464-482.
#' @keywords  management strategy evaluation fisheries
NULL
