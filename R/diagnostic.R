
#' Preliminary Assessments in MSE
#'
#' Evaluates the likely performance of Assessment models in the operating model. This function will
#' apply the assessment model for Data generated during the historical period of the MSE, and report
#' the convergence rate for the model and total time elapsed in running the assessments.
#'
#' @param x Either a `Hist`, `Data` or `OM` object.
#' @param Assess An Assess function of class `Assess`.
#' @param ncpus Numeric, the number of CPUs to run the Assessment model (will run in parallel if greater than 1).
#' @param ... Arguments to be passed to `Assess`, e.g., model configurations.
#'
#' @return Returns invisibly a list of [Assessment-class] objects of length `OM@@nsim`. Messages via console.
#' @author Q. Huynh
#' @examples
#' \donttest{
#' prelim_AM(MSEtool::SimulatedData, SP)
#' }
#' @export
prelim_AM <- function(x, Assess, ncpus = NULL, ...) {
  
  if (is.numeric(ncpus) && !snowfall::sfIsRunning()) {
    MSEtool::setup(cpus = ncpus)
    on.exit(snowfall::sfStop())
  }
  
  # is Hist?
  if (inherits(x, "Hist")) {
    Data <- x@Data
  } else if (inherits(x, "OM")) {
    message("Generating Hist object from OM via runMSE...")
    runHist <- runMSE(x, Hist = TRUE, silent = TRUE, parallel = snowfall::sfIsRunning())
    Data <- runHist@Data
  } else if (inherits(x, "Data")) {
    Data <- x
  } else {
    stop("x does not appear to be either a Hist, Data, or OM object.")
  }
  
  Assess_char <- as.character(substitute(Assess))
  
  nsim <- nrow(Data@Cat)
  dots <- list(...)
  if (length(dots) > 0) message("\nAdditional arguments to be provided to ", deparse(substitute(Assess)), ":\n", paste(names(dots), collapse = "\n"))
  Assess <- match.fun(Assess)
  if (!inherits(Assess, "Assess")) stop(paste(deparse(substitute(Assess))), "does not appear to be an Assess function.")

  if (snowfall::sfIsRunning()) snowfall::sfExport(list = c("Assess", "Data"))
  timing <- proc.time()
  message(paste0("Running ", Assess_char, " with ", nsim, " simulations for ", deparse(substitute(x)), "."))
  res <- pblapply(1:nsim, Assess, Data = Data, ..., cl = if (snowfall::sfIsRunning()) snowfall::sfGetCluster() else NULL)
  timing2 <- (proc.time() - timing)[3]
  message("Assessments complete.")

  message(paste0("Total time to run ", nsim, " assessments: ", round(timing2, 1), " seconds"))

  nonconv <- !vapply(res, getElement, logical(1), "conv")
  message(paste0(sum(nonconv), " of ", nsim, " simulations (", round(100 *sum(nonconv)/nsim, 1), "%) failed to converge."))
  if (sum(nonconv > 0)) message(paste("See simulation number:", paste(which(nonconv), collapse = " ")))

  return(invisible(res))
}





#' Diagnostic of assessments in MSE: did Assess models converge during MSE?
#'
#' Diagnostic check for convergence of Assess models during closed-loop simulation. Use when the MP was 
#' created with [make_MP] with argument `diagnostic = "min"` or `"full"`. 
#' This function summarizes and plots the diagnostic information.
#'
#' @param MSE An object of class MSE created by [MSEtool::runMSE()].
#' @param MP Optional, a character vector of MPs that use assessment models.
#' @param gradient_threshold The maximum magnitude (absolute value) desired for the gradient of the likelihood.
#' @param figure Logical, whether a figure will be drawn.
#' @param ... Arguments to pass to `diagnostic`.
#' @return A matrix with diagnostic performance of assessment models in the MSE. If `figure = TRUE`,
#' a set of figures: traffic light (red/green) plots indicating whether the model converged (defined if a positive-definite
#' Hessian matrix was obtained), the optimizer reached pre-specified iteration limits (as passed to [stats::nlminb()]),
#' and the maximum gradient of the likelihood in each assessment run. Also includes the number of optimization iterations
#' function evaluations reported by [stats::nlminb()] for each application of the assessment model.
#' @author Q. Huynh
#' @aliases diagnostic_AM
#' @examples 
#' \donttest{
#' OM <- MSEtool::testOM; OM@@proyears <- 20
#' myMSE <- runMSE(OM, MPs = "SCA_4010")
#' diagnostic(myMSE)
#' 
#' # How to get all the reporting
#' library(dplyr)
#' conv_statistics <- lapply(1:myMSE@@nMPs, function(m) {
#'   lapply(1:myMSE@@nsim, function(x) {
#'     myMSE@@PPD[[m]]@@Misc[[x]]$diagnostic %>%
#'       mutate(MP = myMSE@@MPs[m], Simulation = x)
#'  }) %>% bind_rows()
#' }) %>% bind_rows()
#' }
#' @seealso [retrospective_AM]
#' @export
diagnostic <- function(MSE, MP, gradient_threshold = 0.1, figure = TRUE) {
  if (!inherits(MSE, "MSE")) stop("No object of class MSE was provided.")

  if (figure) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
  }

  MPs <- MSE@MPs
  has_diagnostic_fn <- function(Data) {
    Misc <- Data@Misc[1:MSE@nsim]
    all(vapply(Misc, function(y) any(names(y) == "diagnostic"), logical(1)))
  }

  has_diagnostic <- vapply(MSE@PPD, has_diagnostic_fn, logical(1))
  if (all(!has_diagnostic)) stop("No diagnostic information found in MSE@PPD for any MP. Use an MP created by: make_MP(diagnostic = \"min\").")
  
  MPs <- MPs[has_diagnostic]

  if (!missing(MP)) {
    match_ind <- pmatch(MP, MPs)
    if (is.na(match_ind)) stop(paste(MP, "MP was not found in the MSE object. Available options are:", paste(MPs, collapse = " ")))
    MPs <- MPs[match_ind]
  }

  message(paste0("Creating plots for MP:\n", paste(MPs, collapse = "\n")))
  res_mat <- matrix(NA, ncol = length(MPs), nrow = 5)
  for(i in 1:length(MPs)) {
    objects <- MSE@PPD[[grep(MPs[i], MSE@MPs)]]
    diagnostic <- lapply(objects@Misc[1:MSE@nsim], getElement, "diagnostic")

    if (!is.null(diagnostic)) {
      if (is.data.frame(diagnostic[[1]])) {
        
        hessian_code <- lapply(diagnostic, getElement, "hess")
        msg <- lapply(diagnostic, getElement, "msg")
        max_gr <- lapply(diagnostic, getElement, "maxgrad")
        iter <- lapply(diagnostic, getElement, "iter")
        fn_eval <- lapply(diagnostic, getElement, "fn_eval")
        Year <- lapply(diagnostic, getElement, "Year")
        
      } else { # Backwards compatibility to version < 1.5.0
        
        get_code <- function(x, y) sapply(x, getElement, y)
        hessian_code <- lapply(diagnostic, get_code, y = "hess")
        msg <- lapply(diagnostic, get_code, y = "msg")
        max_gr <- lapply(diagnostic, get_code, y = "maxgrad")
        iter <- lapply(diagnostic, get_code, y = "iter")
        fn_eval <- lapply(diagnostic, get_code, y = "fn_eval")
        Year <- lapply(diagnostic, get_code, y = "Year")
        
      }

      if (figure) {
        layout(matrix(c(1, 2, 3, 4, 4, 5), ncol = 3, byrow = TRUE))
        par(mar = c(5, 4, 1, 1), oma = c(0, 0, 8, 0))

        plot_hessian_convergence(hessian_code, Year)
        plot_msg(msg, Year)
        plot_max_gr(max_gr, Year, gradient_threshold)

        plot_iter(iter, Year, "line", "nlminb iterations")
        plot_iter(fn_eval, Year, "hist", "nlminb function evaluations")

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

#' @rdname diagnostic
#' @export
diagnostic_AM <- function(...) .Deprecated("diagnostic")

plot_hessian_convergence <- function(convergence_code, Year) {
  nsim <- length(convergence_code)
  yr_range <- range(do.call(c, Year))

  plot(x = NULL, y = NULL, xlim = c(-0.5, 0.5) + yr_range, ylim = c(0.5, nsim + 0.5),
       xlab = "MSE year", ylab = "Simulation #", xaxs = "i", yaxs = "i", axes = FALSE)
  axis(1, at = pretty(yr_range[1]:yr_range[2]) %>% round() %>% unique())
  axis(2, at = pretty(1:nsim) %>% round() %>% unique())
  box()
  mtext('Is Hessian matrix \n positive definite? \n (green: yes, red: no)')
  
  lapply(1:nsim, function(i) {
    lapply(1:length(convergence_code[[i]]), function(j) {
      xcoord <- rep(c(Year[[i]][j] - 0.5, Year[[i]][j] + 0.5), each = 2)
      ycoord <- c(i - 0.5, i + 0.5, i + 0.5, i - 0.5)
      polygon(x = xcoord, y = ycoord, border = 'black', col = ifelse(convergence_code[[i]][j], 'green', 'red'))
    })
  })
  
  return(invisible())
}

plot_msg <- function(convergence_code, Year) {
  nsim <- length(convergence_code)
  yr_range <- range(do.call(c, Year))

  plot(x = NULL, y = NULL, xlim = c(-0.5, 0.5) + yr_range, ylim = c(0.5, nsim + 0.5),
       xlab = "MSE year", ylab = "Simulation #", xaxs = "i", yaxs = "i", axes = FALSE)
  axis(1, at = pretty(yr_range[1]:yr_range[2]) %>% round() %>% unique())
  axis(2, at = pretty(1:nsim) %>% round() %>% unique())
  box()
  mtext('Did nlminb reach iteration limit? \n (red: yes)')
  
  lapply(1:nsim, function(i) {
    lapply(1:length(convergence_code[[i]]), function(j) {
      if (convergence_code[[i]][j] == "function evaluation limit reached without convergence (9)" ||
         convergence_code[[i]][j] == "iteration limit reached without convergence (10)") {
        color <- "red"
      } else color <- NA
      xcoord <- rep(c(Year[[i]][j] - 0.5, Year[[i]][j] + 0.5), each = 2)
      ycoord <- c(i - 0.5, i + 0.5, i + 0.5, i - 0.5)
      polygon(x = xcoord, y = ycoord, border = 'black', col = color)
    })
  })
  
  return(invisible())
}

#' @importFrom gplots rich.colors
plot_iter <- function(x, Year, plot_type = c('line', 'hist'), 
                      lab = c("nlminb iterations", "nlminb function evaluations")) {
  plot_type <- match.arg(plot_type)
  lab <- match.arg(lab)
  x_plot <- do.call(c, x)
  if (plot_type == "hist") {
    if (all(is.na(x_plot))) { # For Perfect()
      hist(0, ylim = c(0, 1), main = "", xlab = lab)
    } else {
      hist(x_plot[!is.na(x_plot)], main = "", xlab = lab)
    }
    mtext(paste("Median:", round(median(do.call(c, x), na.rm = TRUE), 2)))
  }

  if (plot_type == "line") {
    nsim <- length(x)
    yr_range <- range(do.call(c, Year))

    color.vec <- rich.colors(nsim)
    if (all(is.na(x_plot))) {
      ylim <- c(0, 1)
    } else {
      ylim <- c(0.9, 1.1) * range(x_plot, na.rm = TRUE)
    }
    plot(x = NULL, y = NULL, xlim = yr_range, ylim = ylim, xlab = "MSE Year", ylab = lab)
    mtext(lab)
    
    lapply(1:nsim, function(i) {
      points(Year[[i]], x[[i]], col = color.vec[i], typ = 'l')
      text(Year[[i]], x[[i]], labels = i, col = color.vec[i])
    })
  }

  return(invisible())
}

plot_max_gr <- function(max_gr, Year, threshold = 1) {
  nsim <- length(max_gr)
  yr_range <- range(do.call(c, Year))

  plot(x = NULL, y = NULL, xlim = c(-0.5, 0.5) + yr_range, ylim = c(0.5, nsim + 0.5),
       xlab = "MSE Year", ylab = "Simulation #", xaxs = "i", yaxs = "i", axes = FALSE)
  axis(1, at = pretty(yr_range[1]:yr_range[2]) %>% round() %>% unique())
  axis(2, at = pretty(1:nsim) %>% round() %>% unique())
  box()
  mtext(paste0('Maximum gradient of likelihood \n (green: < ', threshold, ', red: > ', threshold, ')'))
  
  lapply(1:nsim, function(i) {
    lapply(1:length(max_gr[[i]]), function(j) {
      xcoord <- rep(c(Year[[i]][j] - 0.5, Year[[i]][j] + 0.5), each = 2)
      ycoord <- c(i - 0.5, i + 0.5, i + 0.5, i - 0.5)
      polygon(x = xcoord, y = ycoord, border = 'black', col = ifelse(max_gr[[i]][j] < threshold, 'green', 'red'))
    })
  })
  
  return(invisible())
}



# Call in MP created by make_MP when diagnostic = "min" or "full"
#' @importFrom dplyr bind_rows
Assess_diagnostic <- function(x, Data, Assessment, include_assessment = TRUE) {

  # Update reporting objects
  Year <- Data@Year[length(Data@Year)]
  if (inherits(Assessment, "Assessment")) {
    msg <- ifelse(is.character(Assessment@opt), Assessment@opt, Assessment@opt$message)
    hess <- ifelse(is.character(Assessment@SD), FALSE, Assessment@SD$pdHess)
    maxgrad <- ifelse(is.character(Assessment@SD), 1e10, max(abs(Assessment@SD$gradient.fixed)))
    iter <- ifelse(is.character(Assessment@opt), NA, Assessment@opt$iterations)
    fn_eval <- ifelse(is.character(Assessment@opt), NA, Assessment@opt$evaluations[1])
    
    dg <- data.frame(Year = Year, hess = hess, msg = msg, maxgrad = maxgrad, iter = iter, fn_eval = fn_eval)
  } else {
    dg <- data.frame(Year = Year, hess = FALSE, msg = "", maxgrad = 1e10, iter = NA, fn_eval = NA)
  }

  # Assign report objects to output
  if (!length(Data@Misc)) Data@Misc <- vector("list", length(Data@Mort))
  output <- list()
  if (is.null(Data@Misc[[x]]$diagnostic) || !length(Data@Misc[[x]]$diagnostic)) {
    output$diagnostic <- dg
  } else {
    output$diagnostic <- rbind(Data@Misc[[x]]$diagnostic, dg)
  }

  if (include_assessment) {
    if (inherits(Assessment, "Assessment")) { # Remove some objects to save memory/disk space
      vars <- c("FMort", "SSB", "R", "U", "UMSY", "FMSY", "SSBMSY", "R0", "h", "SSB0", "MSY", "VB", "conv", 
                "dynamic_SSB0", "Index", "q", "M")
      vars_TMBreport <- c("dynamic_SSB0", "q", "M")
      Assessment_report <- lapply(vars, function(x) {
        
        if (x %in% vars_TMBreport) {
          v <- Assessment@TMB_report[[x]]
          if (x == "M" && is.matrix(v)) {
            v <- v[1:length(Data@Year), 1] %>% structure(names = Data@Year)
          }
        } else {
          v <- slot(Assessment, x)
        }
        if (length(v)) {
          if (is.matrix(v)) {
            if (!requireNamespace("reshape2", quietly = TRUE)) {
              stop("Please install the reshape2 package.", call. = FALSE)
            }
            vmelt <- reshape2::melt(v)
            data.frame(Year_assess = Year,
                       Year_est = vmelt$Var1,
                       variable = vmelt$Var2,
                       value = as.numeric(vmelt$value))
            
          } else {
            Year_est <- suppressWarnings(as.numeric(names(v)))
            xout <- x
            if (is.null(Year_est) || !length(Year_est)) {
              Year_est <- NA_real_
              if (length(v) > 1) xout <- paste0(x, "", 1:length(v))
            } 
            data.frame(Year_assess = Year,
                       Year_est = Year_est,
                       variable = xout,
                       value = as.numeric(v)) # for conv
          }
         
        } else {
          data.frame()
        }
      }) %>% bind_rows()
    } else {
      Assessment_report <- data.frame()
    }
    if(is.null(Data@Misc[[x]]$Assessment_report) || !length(Data@Misc[[x]]$Assessment_report)) {
      output$Assessment_report <- Assessment_report
    } else {
      output$Assessment_report <- rbind(Data@Misc[[x]]$Assessment_report, Assessment_report)
    }
  }
  
  return(output)
}


