
#' @name retrospective
#' @title Retrospective analysis of assessment models
#'
#' @description Perform a retrospective analysis, successive removals of most recent years of data to evaluate resulting
#' parameter estimates.
#'
#' @param x An S4 object of class \linkS4class{Assessment} of \linkS4class{RCModel}.
#' @param nyr The maximum number of years to remove for the retrospective analysis.
#' @param figure Indicates whether plots will be drawn.
#' @param ... More arguments.
#' @return A list with an array of model output and of model estimates from
#' the retrospective analysis.
#' @author Q. Huynh
#' @return Figures showing the time series of biomass and exploitation and parameter estimates
#' with successive number of years removed. For a variety of time series output (SSB, recruitment, etc.) and
#' estimates (R0, steepness, etc.), also returns a matrix of Mohn's rho (Mohn 1999).
#' @examples
#' \donttest{
#' output <- SP(Data = swordfish)
#' get_retro <- retrospective(output, nyr = 5, figure = FALSE)
#' }
#' @references
#' Mohn, R. 1999. The retrospective problem in sequential population analysis: an investigation using cod fishery
#' and simulated data. ICES Journal of Marine Science 56:473-488.
#' @export
setGeneric("retrospective", function(x, ...) standardGeneric("retrospective"))

#' @rdname retrospective
#' @aliases retrospective,Assessment-method
#' @exportMethod retrospective
setMethod("retrospective", signature(x = "Assessment"),
          function(x, nyr = 5, figure = TRUE) {
            if (figure) {
              old_par <- par(no.readonly = TRUE)
              on.exit(par(old_par))
            }
            
            func <- paste0('retrospective_', x@Model)
            res <- do.call2(func, list(Assessment = x, nyr = nyr))
            if (figure) plot(res)
            return(res)
          })


#' @rdname retrospective
#' @aliases retrospective,RCModel-method
#' @exportMethod retrospective
setMethod("retrospective", signature(x = "RCModel"),
          function(x, nyr = 5, figure = TRUE) {
            if (figure) {
              old_par <- par(no.readonly = TRUE)
              on.exit(par(old_par))
            }
            res <- RCM_retro(x, nyr)
            if (figure) plot(res)
            return(res)
          })

#' Class-\code{retro}
#'
#' An S4 class that contains output from \link{retrospective}.
#'
#' @name retro-class
#' @docType class
#'
#' @slot Model Name of the assessment model.
#' @slot Name Name of Data object.
#' @slot TS_var Character vector of time series variables, e.g. recruitment, biomass, from the assessment.
#' @slot TS An array of time series assessment output of dimension, indexed by: peel (the number of terminal years removed from the base assessment),
#' years, and variables (corresponding to \code{TS_var}).
#' @slot Est_var Character vector of estimated parameters, e.g. R0, steepness, in the assessment.
#' @slot Est An array for estimated parameters of dimension, indexed by: peel, variables (corresponding to \code{Est_var}), and
#' value (length 2 for estimate and standard error).
#' @seealso \link{plot.retro} \link{summary.retro} \link{plot.Assessment}
#' @author Q. Huynh
#' @export retro
#' @exportClass retro
retro <- setClass("retro", slots = c(Model = "character", Name = "character", TS_var = "character",
                                     TS = "array", Est_var = "character", Est = "array"))

#' @rdname plot.Assessment
#' @aliases plot,Assessment,retro-method
#' @exportMethod plot
setMethod("plot", signature(x = "Assessment", y = "retro"),
          function(x, y, filename = paste0("report_", x@Model), dir = tempdir(), open_file = TRUE, quiet = TRUE,
                   render_args = list(), ...) {
            report(x, y, filename = filename, dir = dir, open_file = open_file, quiet = quiet, render_args = render_args, ...)
          })

#' @name plot.retro
#' @title Methods for retro object
#' @description plot and summary functions for retro object.
#' @param object An object of class \linkS4class{retro}.
#' @param x An object of class \linkS4class{retro}.
#' @param color An optional character vector of colors for plotting.
#' @author Q. Huynh
#' @aliases plot.retro plot,retro,missing-method
#' @return A series of plots showing retrospective patterns in fishing mortality, spawning biomass, recruitment, etc.
#' @examples
#' \donttest{
#' res <- SP(Data = swordfish)
#' ret <- retrospective(res, figure = FALSE)
#'
#' summary(ret)
#' plot(ret)
#' }
#' @exportMethod plot
setMethod("plot", signature(x = "retro", y = "missing"),
          function(x, color = NULL) {
            old_par <- par(no.readonly = TRUE)
            on.exit(par(old_par))

            retro <- x
            if (is.null(color) || length(color) != dim(retro@TS)[1]) color <- rich.colors(dim(retro@TS)[1])
            Year <- as.numeric(dimnames(retro@TS)$Year)
            xlim <- range(as.numeric(dimnames(retro@TS)$Year))
            nyr_label <- dimnames(retro@TS)$Peel

            for(i in 1:length(retro@TS_var)) {
              matrix_to_plot <- t(retro@TS[, , i])
              ylim <- c(0, 1.1 * max(matrix_to_plot, na.rm = TRUE))
              ylab <- attr(retro, "TS_lab")[i]

              plot(NULL, NULL, xlim = xlim, ylim = ylim, xlab = "Year", ylab = ylab)
              abline(h = 0, col = "grey")
              if (grepl("MSY", as.character(ylab))) abline(h = 1, lty = 3)

              matlines(Year, matrix_to_plot, col = color, lty = 1)

              legend("topleft", legend = nyr_label, lwd = 1, col = color, bty = "n", title = "Years removed:")
            }
            invisible()
          })


#' @rdname plot.retro
#' @aliases summary.retro summary,retro-method
#' @exportMethod summary
setMethod("summary", signature(object = "retro"), function(object) calculate_Mohn_rho(object@TS, ts_lab = attr(object, "TS_lab")))

setMethod("show", signature(object = "retro"), function(object) print(summary(object)))
