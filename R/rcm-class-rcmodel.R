
#' Class-`rcmodel`
#'
#' Output of [RCM2()].
#'
#' @slot OM A populated [MSEtool::om-class] object conditioned on the fitted values. Each stock has the
#' fitted unfished recruitment, steepness and recruitment deviations in `@@SRR`, and each fleet has the
#' fitted apical fishing mortality (as `@@Effort` with unit catchability) and selectivity-at-age, with
#' full retention so that all catch is removed.
#' The fitted data are in `@@Data`, with index selectivity set in `@@Obs`.
#' @slot fits A list, named by stock (or by complex for hermaphroditic stocks), of [RCModel-class] objects
#' with each fit. Use [plot()] on an element for the fit diagnostics report.
#' @slot RCMdata A list, named as `fits`, of the [RCMdata-class] objects that were fitted.
#' @slot Misc A list for additional output.
#'
#' @seealso [RCM2()]
#' @name rcmodel-class
#' @rdname rcmodel
#' @exportClass rcmodel
setClass("rcmodel", slots = c(
  OM = "om", fits = "list", RCMdata = "list", Misc = "list"
))

#' @rdname rcmodel
#' @param object An `rcmodel` object.
#' @export
setMethod("show", signature(object = "rcmodel"), function(object) {
  cat("rcmodel:", length(object@fits), "stock(s),", object@OM@nSim, "simulation(s)
")
  for (st in names(object@fits)) {
    fit <- object@fits[[st]]
    R0 <- vapply(fit@report, getElement, numeric(1), "R0")
    E0 <- vapply(fit@report, getElement, numeric(1), "E0_SR")
    dep <- fit@SSB[, ncol(fit@SSB) - 1] / E0
    cat(sprintf("  %s: %d/%d converged; median R0 %s; median terminal SSB/SSB0 %s
", st, sum(fit@conv),
                length(fit@conv), signif(stats::median(R0), 3), signif(stats::median(dep), 3)))
  }
  invisible(object)
})

#' @rdname rcmodel
#' @param x An `rcmodel` object.
#' @param y Not used.
#' @param stock Name or number of the stock to report. Defaults to the first stock.
#' @param ... Arguments passed to [plot.RCModel].
#' @export
setMethod("plot", signature(x = "rcmodel", y = "missing"), function(x, y, stock = 1, ...) {
  fit <- x@fits[[stock]]
  if (is.null(fit)) cli::cli_abort("Stock {.val {stock}} is not in {.arg x}.")
  plot(fit, ...)
})
