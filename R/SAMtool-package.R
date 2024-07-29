#' Stock Assessment Methods Toolkit 
#'
#' Simulation tools for closed-loop simulation are provided for the 'MSEtool' operating model to inform data-rich fisheries. 
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
#' with MSEtool operating models. An overview of these features is available 
#' on the [openMSE](https://openmse.com/features-assessment-models/) website.
#'
#' The [RCM()] (Rapid Conditioning Model) can be used to condition operating models from real data. 
#' 
#' The following articles are available on the openMSE website:
#' \itemize{
#' \item [Description of assessment models](https://openmse.com/features-assessment-models/)
#' \item [General overview of RCM](https://openmse.com/tutorial-rcm/)
#' }
#' 
#' The function documentation can be viewed [online](https://samtool.openmse.com/reference/index.html).
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
"_PACKAGE"
