
#' Preliminary Assessments in MSE
#'
#' Evaluates the likely performance of Assessment models in the operating model. This function will
#' apply the assessment model for Data generated during the historical period of the MSE, and report
#' the convergence rate for the model and total time elapsed in running the assessments.
#'
#' @param x Either a \code{Hist}, \code{Data} or \code{OM} object.
#' @param Assess An Assess function of class \code{Assess}.
#' @param ncpus Numeric, the number of CPUs to run the Assessment model (will run in parallel if greater than 1).
#' @param ... Arguments to be passed to \code{Assess}, e.g., model configurations.
#'
#' @return Returns invisibly a list of \linkS4class{Assessment} objects of length \code{OM@@nsim}. Messages via console.
#' @author Q. Huynh
#' @examples
#' \dontrun{
#' prelim_AM(DLMtool::testOM, DD_TMB)
#' }
#' @export
prelim_AM <- function(x, Assess, ncpus = NULL, ...) {

  # is Hist?
  if(is.list(x) && any(names(x) == "Data")) {
    Data <- x$Data
  } else
    if(inherits(x, "OM")) {
      message("Generating Hist object from OM object via runMSE:")
      runHist <- runMSE(x, Hist = TRUE)
      if(packageVersion("DLMtool") >= 5.3) Data <- runHist@Data else Data <- runHist$Data
    } else
      if(inherits(x, "Data")) {
        Data <- x
      } else {
        stop("x does not appear to be either a Hist, Data, or OM object.")
      }

  if(is.numeric(ncpus) && !snowfall::sfIsRunning()) {
    DLMtool::setup(cpus = ncpus)
    on.exit(snowfall::sfStop())
  }
  nsim <- nrow(Data@Cat)
  dots <- list(...)
  if(length(dots) > 0) message("\nAdditional arguments to be provided to ", deparse(substitute(Assess)), ":\n", paste(names(dots), collapse = "\n"))
  Assess <- match.fun(Assess)
  if(!inherits(Assess, "Assess")) stop(paste(deparse(substitute(Assess))), "does not appear to be an Assess function.")

  if(snowfall::sfIsRunning()) snowfall::sfExport(list = c("Assess", "Data"))
  timing <- proc.time()
  if(snowfall::sfIsRunning()) {
    message("Running ", deparse(substitute(Assess)), " with ", nsim, " simulations for ", deparse(substitute(x)), " on ", snowfall::sfCpus(), " CPUs.")
    res <- snowfall::sfClusterApplyLB(1:nsim, Assess, Data = Data, ...)
  } else {
    message(paste0("Running ", deparse(substitute(Assess)), " with ", nsim, " simulations for ", deparse(substitute(x)), "."))
    res <- lapply(1:nsim, Assess, Data = Data, ...)
  }
  timing2 <- (proc.time() - timing)[3]
  message("Assessments complete.")

  message(paste0("Total time to run ", nsim, " assessments: ", round(timing2, 1), " seconds"))

  nonconv <- !vapply(res, getElement, logical(1), "conv")
  message(paste0(sum(nonconv), " of ", nsim, " simulations (", round(100 *sum(nonconv)/nsim, 1), "%) failed to converge."))
  if(sum(nonconv > 0)) message(paste("See simulation number:", paste(which(nonconv), collapse = " ")))

  return(invisible(res))
}





#' diagnostic_AM (diagnostic of Assessments in MSE): Did Assess models converge during MSE?
#'
#' Diagnostic check for convergence of Assess models during MSE.
#' Assess models write output to the DLMenv environment if the MP was created with \link{make_MP}
#' with argument \code{diagnostic = TRUE}. This function summarizes and plots the diagnostic information.
#'
#' @param MSE An object of class MSE created by \code{\link[DLMtool]{runMSE}}. If no MSE object
#' is available, use argument \code{MP} instead.
#' @param MP A character vector of MPs with assessment models.
#' @param gradient_threshold The maximum magnitude (absolute value) desired for the gradient of the likelihood.
#' @param figure Logical, whether a figure will be drawn.
#' @return A matrix with diagnostic performance of assessment models in the MSE. If \code{figure = TRUE},
#' a set of figures: traffic light (red/green) plots indicating whether the model converged (defined if a positive-definite
#' Hessian matrix was obtained), the optimizer reached pre-specified iteration limits (as passed to \code{\link[stats]{nlminb}}),
#' and the maximum gradient of the likelihood in each assessment run. Also includes the number of optimization iterations
#' function evaluations reported by \code{\link[stats]{nlminb}} for each application of the assessment model.
#' @author Q. Huynh
#' @examples
#' \dontrun{
#' DD_MSY <- make_MP(DD_TMB, HCR_MSY, diagnostic = "min")
#' show(DD_MSY)
#'
#' ##### Ensure that PPD = TRUE in runMSE function
#' myMSE <- runMSE(DLMtool::testOM, MPs = "DD_MSY", PPD = TRUE)
#' diagnostic_AM(myMSE)
#' }
#' @importFrom graphics layout
#' @seealso \link{retrospective_AM}
#' @export
diagnostic_AM <- function(MSE, MP = NULL, gradient_threshold = 0.1, figure = TRUE) {
  if(!inherits(MSE, "MSE")) stop("No object of class MSE was provided.")
  if(packageVersion("DLMtool") >= 5.3) {
    if(length(MSE@Misc$Data) == 0) stop("Nothing found in MSE@Misc$Data. Use an MP created by 'make_MP(diagnostic = 'min')' and set 'runMSE(PPD = TRUE)'.")
  } else {
    if(length(MSE@Misc) == 0) stop("Nothing found in MSE@Misc. Use an MP created by 'make_MP(diagnostic = 'min')' and set 'runMSE(PPD = TRUE)'.")
  }

  if(figure) {
    old_par <- par(no.readonly = TRUE)
    on.exit(layout(matrix(1)))
    on.exit(par(old_par), add = TRUE)

    par(mar = c(5, 4, 1, 1), oma = c(0, 0, 8, 0))
  }

  MPs <- MSE@MPs
  has_diagnostic_fn <- function(Data) {
    Misc <- Data@Misc
    all(vapply(Misc, function(y) any(names(y) == "diagnostic"), logical(1)))
  }

  if(packageVersion("DLMtool") >= 5.3) {
    has_diagnostic <- vapply(MSE@Misc$Data, has_diagnostic_fn, logical(1))
    if(all(!has_diagnostic)) stop("No diagnostic information found in MSE@Misc$Data for any MP. Use an MP created by 'make_MP(diagnostic = 'min')' and set 'runMSE(PPD = TRUE)'.")

  } else {
    has_diagnostic <- vapply(MSE@Misc, has_diagnostic_fn, logical(1))
    if(all(!has_diagnostic)) stop("No diagnostic information found in MSE@Misc for any MP. Use an MP created by 'make_MP(diagnostic = 'min')' and set 'runMSE(PPD = TRUE)'.")
  }

  MPs <- MPs[has_diagnostic]

  if(!is.null(MP)) {
    match_ind <- pmatch(MP, MPs)
    if(is.na(match_ind)) stop(paste(MP, "MP was not found in the MSE object. Available options are:", paste(MPs, collapse = " ")))
    MPs <- MPs[match_ind]
  }

  message(paste0("Creating plots for MP:\n", paste(MPs, collapse = "\n")))
  res_mat <- matrix(NA, ncol = length(MPs), nrow = 5)
  get_code <- function(x, y) vapply(x, getElement, numeric(1), y)
  get_code_char <- function(x, y) vapply(x, getElement, character(1), y)
  for(i in 1:length(MPs)) {
    if(packageVersion("DLMtool") >= 5.3) {
      objects <- MSE@Misc$Data[[which(MPs[i] == MSE@MPs)]]
    } else objects <- MSE@Misc[[which(MPs[i] == MSE@MPs)]]
    diagnostic <- lapply(objects@Misc, getElement, "diagnostic")

    if(!is.null(diagnostic)) {
      hessian_code <- lapply(diagnostic, get_code, y = "hess")
      msg <- lapply(diagnostic, get_code_char, y = "msg")
      max_gr <- lapply(diagnostic, get_code, y = "maxgrad")
      iter <- lapply(diagnostic, get_code, y = "iter")
      fn_eval <- lapply(diagnostic, get_code, y = "fn_eval")
      Year <- lapply(diagnostic, get_code, y = "Year")

      if(figure) {
        layout(matrix(c(1, 2, 3, 4, 4, 5), ncol = 3, byrow = TRUE))

        plot_hessian_convergence(hessian_code, Year)
        plot_msg(msg, Year)
        plot_max_gr(max_gr, Year, gradient_threshold)

        plot_iter(iter, Year, "line", "Optimization iterations")
        plot_iter(fn_eval, Year, "hist", "Function evaluations")

        title(paste(MPs[i], "management procedure"), outer = TRUE)
      }

      hessian_code <- do.call(c, hessian_code)
      msg <- do.call(c, msg)
      msg2 <- msg == "function evaluation limit reached without convergence (9)" |
        msg == "iteration limit reached without convergence (10)"
      max_gr <- do.call(c, max_gr)
      iter <- do.call(c, iter)
      fn_eval <- do.call(c, fn_eval)

      res_mat[, i] <- c(100 * sum(hessian_code)/length(hessian_code),
                        100 * sum(msg2)/length(msg2),
                        100 * sum(abs(max_gr) <= gradient_threshold)/length(max_gr),
                        median(iter, na.rm = TRUE), median(fn_eval, na.rm = TRUE))
    }
  }

  res_mat <- round(res_mat, 2)
  colnames(res_mat) <- MPs
  rownames(res_mat) <- c("Percent positive-definite Hessian", "Percent iteration limit reached",
                         paste("Percent max. gradient <", gradient_threshold),
                         "Median iterations", "Median function evaluations")
  return(res_mat)
}

plot_hessian_convergence <- function(convergence_code, Year) {
  nsim <- length(convergence_code)
  yr_range <- range(do.call(c, Year))

  plot(x = NULL, y = NULL, xlim = c(-0.5, 0.5) + yr_range, ylim = c(0.5, nsim + 0.5),
       xlab = "MSE year", ylab = "Simulation #", xaxs = "i", yaxs = "i")
  mtext('Is Hessian matrix \n positive definite? \n (green: yes, red: no)')
  for(i in 1:nsim) {
    for(j in 1:length(convergence_code[[i]])) {
      color <- ifelse(convergence_code[[i]][j], 'green', 'red')
      xcoord <- rep(c(Year[[i]][j] - 0.5, Year[[i]][j] + 0.5), each = 2)
      ycoord <- c(i - 0.5, i + 0.5, i + 0.5, i - 0.5)
      polygon(x = xcoord, y = ycoord, border = 'black', col = color)
    }
  }
  return(invisible())
}

plot_msg <- function(convergence_code, Year) {
  nsim <- length(convergence_code)
  yr_range <- range(do.call(c, Year))

  plot(x = NULL, y = NULL, xlim = c(-0.5, 0.5) + yr_range, ylim = c(0.5, nsim + 0.5),
       xlab = "MSE year", ylab = "Simulation #", xaxs = "i", yaxs = "i")
  mtext('Did nlminb reach iteration limit? \n (red: yes)')
  for(i in 1:nsim) {
    for(j in 1:length(convergence_code[[i]])) {
      if(convergence_code[[i]][j] == "function evaluation limit reached without convergence (9)" ||
         convergence_code[[i]][j] == "iteration limit reached without convergence (10)") {
        color <- "red"
      } else color <- NA
      xcoord <- rep(c(Year[[i]][j] - 0.5, Year[[i]][j] + 0.5), each = 2)
      ycoord <- c(i - 0.5, i + 0.5, i + 0.5, i - 0.5)
      polygon(x = xcoord, y = ycoord, border = 'black', col = color)
    }
  }
  return(invisible())
}

#' @importFrom gplots rich.colors
plot_iter <- function(x, Year, plot_type = c('line', 'hist'), lab = c("Optimization iterations", "Function evaluations")) {
  plot_type <- match.arg(plot_type)
  lab <- match.arg(lab)

  if(plot_type == "hist") {
    hist(do.call(c, x), main = "", xlab = lab)
    mtext(paste("Median:", round(median(do.call(c, x), na.rm = TRUE), 2)))
  }

  if(plot_type == "line") {
    nsim <- length(x)
    yr_range <- range(do.call(c, Year))

    color.vec <- rich.colors(nsim)
    plot(x = NULL, y = NULL, xlim = yr_range,
         ylim = c(0.9, 1.1) * range(do.call(c, x), na.rm = TRUE),
         xlab = "MSE Year", ylab = lab)
    mtext(lab)

    for(i in 1:nsim) {
      points(Year[[i]], x[[i]], col = color.vec[i], typ = 'l')
      text(Year[[i]], x[[i]], labels = i, col = color.vec[i])
    }
  }

  return(invisible())
}

plot_max_gr <- function(max_gr, Year, threshold = 1) {
  nsim <- length(max_gr)
  yr_range <- range(do.call(c, Year))

  plot(x = NULL, y = NULL, xlim = c(-0.5, 0.5) + yr_range, ylim = c(0.5, nsim + 0.5),
       xlab = "MSE Year", ylab = "Simulation #", xaxs = "i", yaxs = "i")
  mtext(paste0('Maximum gradient of likelihood \n (green: < ', threshold, ', red: > ', threshold, ')'))

  for(i in 1:nsim) {
    for(j in 1:length(max_gr[[i]])) {
      color <- ifelse(max_gr[[i]][j] < threshold, 'green', 'red')
      xcoord <- rep(c(Year[[i]][j] - 0.5, Year[[i]][j] + 0.5), each = 2)
      ycoord <- c(i - 0.5, i + 0.5, i + 0.5, i - 0.5)
      polygon(x = xcoord, y = ycoord, border = 'black', col = color)
    }
  }
  return(invisible())
}



# Call in MP created by make_MP when diagnostic = "min" or "full"
Assess_diagnostic <- function(x, Data, Assessment, include_assessment = TRUE) {

  # Update reporting objects
  Year <- Data@Year[length(Data@Year)]
  if(inherits(Assessment, "Assessment")) {
    msg <- ifelse(is.character(Assessment@opt), Assessment@opt, Assessment@opt$message)
    hess <- ifelse(is.character(Assessment@SD), FALSE, Assessment@SD$pdHess)
    maxgrad <- ifelse(is.character(Assessment@SD), 1e10, max(abs(Assessment@SD$gradient.fixed)))
    iter <- ifelse(is.character(Assessment@opt), NA, Assessment@opt$iterations)
    fn_eval <- ifelse(is.character(Assessment@opt), NA, Assessment@opt$evaluations[1])

    dg <- list(hess = hess, msg = msg, maxgrad = maxgrad, iter = iter, fn_eval = fn_eval, Year = Year)

  } else {
    dg <- list(hess = FALSE, msg = msg, maxgrad = 1e10, iter = NA, fn_eval = NA, Year = Year)
  }

  # Assign report objects to output
  if(length(Data@Misc) == 0) Data@Misc <- vector("list", length(Data@Mort))
  diagnostic <- Data@Misc[[x]]$diagnostic
  len_diag <- length(diagnostic)
  diagnostic[[len_diag + 1]] <- dg

  output <- list(diagnostic = diagnostic)

  if(include_assessment) {
    # Remove some objects to save memory/disk space
    if(inherits(Assessment, "Assessment")) {
      if(dg$hess) Assessment@obj <- list()
      Assessment@info <- Assessment@TMB_report <- list()
      Assessment@N_at_age <- Assessment@C_at_age <- Assessment@Obs_C_at_age <- Assessment@Selectivity <- array(dim = c(0, 0))
    }

    Assessment_report <- Data@Misc[[x]]$Assessment_report
    len_Assess <- length(Assessment_report)
    if(len_Assess > 0) {
      Assessment_report[[len_Assess + 1]] <- Assessment
    } else Assessment_report <- list(Assessment)
    output$Assessment_report <- Assessment_report
  }

  return(output)
}

