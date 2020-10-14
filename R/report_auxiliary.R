

#' Compare output from several assessment models
#'
#' Plot biomass, recruitment, and fishing mortality time series from several . This function can be used to compare outputs among
#' different assessment models from the same Data object.
#'
#' @param ... Objects of class \linkS4class{Assessment}.
#' @param label A character vector of the models for the legend.
#' @param color A vector of colors for each assessment model.
#' @author Q. Huynh
#' @return A set of figures of biomass, recruitment, and fishing mortality estimates among the models.
#' @examples
#' res <- cDD_SS(Data = DLMtool::SimulatedData)
#' res2 <- SCA(Data = DLMtool::SimulatedData)
#' res3 <- SCA2(Data = DLMtool::SimulatedData)
#' res4 <- VPA(Data = DLMtool::SimulatedData)
#'
#' compare_models(res, res2, res3)
#' @importFrom gplots rich.colors
#' @export
compare_models <- function(..., label = NULL, color = NULL) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  dots <- list(...)
  class_dots <- vapply(dots, inherits, logical(1), what = "Assessment")
  if(!all(class_dots)) stop("Some objects provided were not of class Assessment", call. = FALSE)

  n_assess <- length(dots)
  if(n_assess <= 1) stop("Need more than one assessment model for this function.", call. = FALSE)

  if(is.null(label) || length(label) != n_assess) {
    label <- vapply(dots, slot, character(1), name = "Model")
    if(length(unique(label)) != length(dots)) label <- as.character(substitute(list(...)))[-1]
  }
  if(is.null(color) || length(color) != n_assess) {
    color <- rich.colors(n_assess)
  }

  par(mfrow = c(3, 2), mar = c(5, 4, 1, 1), oma = c(2, 0, 0, 0))

  # F
  slot_F <- function(x) {
    if(length(x@FMort) == 0) {
      x@U
    } else {
      x@FMort
    }
  }
  FM <- do.call(rbind, lapply(dots, slot_F))
  ts_matplot(FM, "Fishing Mortality", color = color)

  # F/FMSY
  slot_FMSY <- function(x) {
    if(length(x@F_FMSY) == 0) {
      x@U_UMSY
    } else {
      x@F_FMSY
    }
  }
  F_FMSY <- do.call(rbind, lapply(dots, slot_FMSY))
  ts_matplot(F_FMSY, expression(F/F[MSY]), color = color, dotted_one = TRUE)

  # B/BMSY
  B_BMSY <- do.call(rbind, lapply(dots, slot, name = "SSB_SSBMSY"))
  ts_matplot(B_BMSY, expression(SSB/SSB[MSY]), color = color, dotted_one = TRUE)

  # B/B0
  B_B0 <- do.call(rbind, lapply(dots, slot, name = "SSB_SSB0"))
  ts_matplot(B_B0, expression(SSB/SSB[0]), color = color)

  # VB
  VB <- do.call(rbind, lapply(dots, slot, name = "VB"))
  ts_matplot(VB, "Vulnerable Biomass", color = color)

  # R
  RR <- lapply(dots, slot, name = "R")
  R <- match_R_years(RR)
  if(!all(is.na(R))) ts_matplot(R, "Recruitment", color = color)

  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", label, col = color, xpd = TRUE, horiz = TRUE, bty = "n", lwd = 2)

  invisible()
}

#' @importFrom graphics matlines
ts_matplot <- function(m, ylab, color, dotted_one = FALSE) {
  m <- t(m)
  x <- matrix(as.numeric(rownames(m)), ncol = ncol(m), nrow = nrow(m))
  plot(NULL, NULL, xlim = range(as.numeric(rownames(m))), ylim = c(0, 1.1 * max(m, na.rm = TRUE)), xlab = "Year", ylab = ylab)
  abline(h = 0, col = "grey")
  matlines(x, m, type = "l", col = color, lty = 1, lwd = 2)
  if(dotted_one) abline(h = 1, lty = 3)
}

match_R_years <- function(RR) {
  yrs <- do.call(c, lapply(RR, function(x) if(length(x) > 0) as.numeric(names(x)) else NA))
  if(all(is.na(yrs))) return(NA)

  yrs <- range(yrs, na.rm = TRUE)

  R <- matrix(NA, nrow = length(RR), ncol = diff(yrs) + 1)
  R_yrs <- seq(min(yrs), max(yrs))
  for(i in 1:length(RR)) {
    if(!all(is.na(RR[[i]]))) {
      ind <- match(as.numeric(names(RR[[i]])), R_yrs)
      R[i, ind] <- RR[[i]]
    }
  }
  colnames(R) <- R_yrs
  return(R)
}


#' @importFrom rmarkdown render
report <- function(Assessment, retro = NULL, filename = paste0("report_", Assessment@Model), dir = tempdir(), open_file = TRUE, quiet = TRUE,
                   render_args = list(), ...) {
  name <- ifelse(nchar(Assessment@Name) > 0, Assessment@Name, substitute(Assessment))

  # Generate markdown report
  filename_html <- paste0(filename, ".html")
  filename_rmd <- paste0(filename, ".Rmd")

  if(!dir.exists(dir)) {
    message("Creating directory: \n", dir)
    dir.create(dir)
  }
  message("Writing markdown file: ", file.path(dir, filename_rmd))

  if(Assessment@Model == "SCA2") Assessment@info$data$SR_type <- Assessment@info$SR
  f <- get(paste0("rmd_", Assessment@Model))
  rmd_model <- f(Assessment, ...)

  if(!is.null(retro)) {
    rmd_ret <- rmd_retrospective()
  } else rmd_ret <- ""

  rmd <- c(rmd_head(name), rmd_model, rmd_ret, rmd_footer())
  write(rmd, file = file.path(dir, filename_rmd))

  # Render markdown file
  assign_Assessment_slots(Assessment)

  if(is.null(render_args$input)) render_args$input <- file.path(dir, filename_rmd)
  if(is.null(render_args$output_format)) render_args$output_format <- "html_document"
  if(is.null(render_args$output_options) && render_args$output_format == "html_document") {
    render_args$output_options <- list(df_print = "paged")
  }
  render_args$quiet <- quiet

  message("Rendering markdown file...")
  output_filename <- do.call(rmarkdown::render, render_args)
  message("Rendered file: ", output_filename)

  if(open_file) browseURL(output_filename)
  invisible(output_filename)
}






#' Plots a lognormal variable
#'
#' Plots the probability distribution function of a lognormal variable from the
#' mean and standard deviation in either transformed (normal) or untransformed space.
#'
#' @param m A vector of means of the distribution.
#' @param sd A vector of standard deviations of the distribution.
#' @param label Name of the variable to be used as x-axis label.
#' @param logtransform Indicates whether the mean and standard deviation are in
#' transformed (normal) or untransformed space.
#' @param color A vector of colors.
#' @return A plot of the probability distribution function. Vertical dotted line
#' indicates mean of distribution. This function can plot multiple curves when multiple means
#' and standard deviations are provided.
#' @author Q. Huynh
#' @export plot_lognormalvar
#' @seealso \code{\link{plot_betavar}} \code{\link{plot_steepness}}
#' @examples
#' mu <- 0.5
#' stddev <- 0.1
#' plot_lognormalvar(mu, stddev) # mean of plot should be 0.5
#'
#' #logtransformed parameters
#' mu <- 0
#' stddev <- 0.1
#' plot_lognormalvar(mu, stddev, logtransform = TRUE) # mean of plot should be 1
plot_lognormalvar <- function(m, sd, label = NULL, logtransform = FALSE, color = "black") {
  # plots life history parameters: Linf, K, t0, M, FMSY_M
  ncurve <- length(m)
  if(!logtransform) {
    true.m <- m
    if(all(m < 0)) m <- -1 * m # special case needed when t0 < 0
    mu <- mconv(m, sd)
    sdlog <- sdconv(m, sd)
    support <- seq(0.001, max(m + 5*sdlog), length.out = 1e3)

    dist <- matrix(NA, nrow = length(support), ncol = ncurve)
    for(i in 1:ncurve) dist[, i] <- dlnorm(support, mu[i], sdlog[i])
    dist[is.infinite(dist)] <- NA

    dist.max <- max(dist, na.rm = TRUE)
    tails <- apply(dist, 1, function(x) all(x < 0.001 * dist.max))
    tails <- which(!tails)
    ind.tails <- c(tails[1]:tails[length(tails)])

    support <- support[ind.tails]
    dist <- as.matrix(dist[ind.tails, ])

    if(all(true.m < 0)) {
      support <- -1 * support
      xlim_truncated <- range(pretty(support))
      plot(support, dist[, 1], typ = 'l', xlab = label,
           ylab = 'Probability density function', xlim = xlim_truncated,
           ylim = c(0, 1.1 * max(dist, na.rm = TRUE)), col = color[1])
      if(ncurve > 1) {
        for(i in 2:ncurve) lines(support, dist[, i], col = color[i])
      }
    }
    if(all(true.m > 0)) {
      xlim_truncated <- range(pretty(support))
      plot(support, dist[, 1], typ = 'l', xlab = label,
           ylab = 'Probability density function', xlim = xlim_truncated,
           ylim = c(0, 1.1 * max(dist, na.rm = TRUE)), col = color[1])
      if(ncurve > 1) {
        for(i in 2:ncurve) lines(support, dist[, i], col = color[i])
      }
    }
    abline(h = 0, col = 'grey')
    abline(v = true.m, lty = 2, col = color)
  }

  if(logtransform) {
    #f_y(y) = f_x(g-1(y)) * abs(d/dy[g-1(y)])
    #where f is the pdf of distribution, g(y) = exp(X) is the transformation
    #y is the lognormal variable, x is a normal variable
    support.norm <- seq(min(m - 5*sd, na.rm = TRUE), max(m+5*sd, na.rm = TRUE),
                        length.out = 1e3)
    support <- exp(support.norm)

    dist <- matrix(NA, nrow = length(support), ncol = ncurve)
    for(i in 1:ncurve) dist[, i] <- dnorm(support.norm, m[i], sd[i])/abs(support)
    dist[is.infinite(dist)] <- NA

    dist.max <- max(dist, na.rm = TRUE)
    tails <- apply(dist, 1, function(x) all(x < 0.001 * dist.max))
    tails <- which(!tails)
    ind.tails <- c(tails[1]:tails[length(tails)])

    support <- support[ind.tails]
    dist <- as.matrix(dist[ind.tails, ])

    xlim_truncated <- range(pretty(support), finite = TRUE, na.rm = TRUE)

    plot(support, dist[, 1], typ = 'l', xlab = label, xlim = xlim_truncated,
         ylab = 'Probability density function',
         ylim = c(0, 1.1 * max(dist, na.rm = TRUE)), col = color[1])
    if(ncurve > 1) {
      for(i in 2:ncurve) lines(support, dist[, i], col = color[i])
    }
    abline(h = 0, col = 'grey')
    abline(v = exp(m), lty = 2, col = color)
  }

  invisible()
}


plot_normalvar <- function(m, sd, label = NULL, color = "black") {
  ncurve <- length(m)
  support <- seq(max(0, min(m - 5 * sd)), max(m + 5 * sd), length.out = 1e3)
  dist <- matrix(NA, nrow = length(support), ncol = ncurve)
  for(i in 1:ncurve) dist[, i] <- dnorm(support, m[i], sd[i])
  dist[is.infinite(dist)] <- NA

  dist.max <- max(dist, na.rm = TRUE)
  tails <- apply(dist, 1, function(x) all(x < 0.001 * dist.max))
  tails <- which(!tails)
  ind.tails <- c(tails[1]:tails[length(tails)])

  support <- support[ind.tails]
  dist <- as.matrix(dist[ind.tails, ])

  xlim_truncated <- range(pretty(support))
  plot(support, dist[, 1], typ = 'l', xlab = label, ylab = 'Probability density function',
       xlim = xlim_truncated, ylim = c(0, 1.1 * max(dist, na.rm = TRUE)), col = color[1])
  if(ncurve > 1) {
    for(i in 2:ncurve) lines(support, dist[, i], col = color[i])
  }
  abline(h = 0, col = 'grey')
  abline(v = m, lty = 2, col = color)

  invisible()
}

#' Plots a beta variable
#'
#' Plots the probability distribution function of a beta variable from the
#' mean and standard deviation in either transformed (logit) or untransformed space.
#'
#' @param m A vector of means of the distribution.
#' @param sd A vector of standard deviations of the distribution.
#' @param label Name of the variable to be used as x-axis label.
#' @param is_logit Logical that indicates whether the means and standard deviations are in
#' transformed (logit) or untransformed space.
#' @param color A vector of colors.
#' @return A plot of the probability distribution function. Vertical dotted line
#' indicates mean of distribution. This function can plot multiple curves when multiple means
#' and standard deviations are provided.
#' @author Q. Huynh
#' @export plot_betavar
#' @seealso \code{\link{plot_lognormalvar}} \code{\link{plot_steepness}}
#' @examples
#' mu <- 0.5
#' stddev <- 0.1
#' plot_betavar(mu, stddev) # mean of plot should be 0.5
#'
#' #logit parameters
#' mu <- 0
#' stddev <- 0.1
#' plot_betavar(mu, stddev, is_logit = TRUE) # mean of plot should be 0.5
plot_betavar <- function(m, sd, label = NULL, is_logit = FALSE, color = "black") {
  support <- seq(0.01, 0.99, length.out = 1e3)
  ncurve <- length(m)
  dist <- matrix(NA, nrow = length(support), ncol = ncurve)
  if(!is_logit) {
    a <- alphaconv(m, sd)
    b <- betaconv(m, sd)
    for(i in 1:ncurve) dist[, i] <- dbeta(support, a[i], b[i])
  }
  if(is_logit) {
    for(i in 1:ncurve) {
      #f_y(y) = f_x(g-1(y)) * abs(d/dy[g-1(y)])
      #where f is the pdf of distribution, g(y) = 1/(1 + exp(-X)) is the transformation
      #y is the beta variable, x is a normal variable
      dist[, i] <- dnorm(log(support/(1-support)), m[i], sd[i])/abs(support * (1-support))
      m[i] <- ilogit(m[i])
    }
  }
  dist[is.infinite(dist)] <- NA

  dist.max <- max(dist, na.rm = TRUE)
  tails <- apply(dist, 1, function(x) all(x < 0.001 * dist.max))
  tails <- which(!tails)
  ind.tails <- c(tails[1]:tails[length(tails)])

  support <- support[ind.tails]
  dist <- as.matrix(dist[ind.tails, ])
  xlim_truncated <- range(pretty(support), finite = TRUE)

  plot(support, dist[, 1], typ = 'l', xlab = label, ylab = 'Probability density function',
       xlim = xlim_truncated, ylim = c(0, 1.1 * max(dist, na.rm = TRUE)), col = color[1])
  if(ncurve > 1) {
    for(i in 2:ncurve) lines(support, dist[, i], col = color[i])
  }
  abline(h = 0, col = 'grey')
  abline(v = m, lty = 2, col = color)

  invisible()
}

#' Plots probability distribution function of stock-recruit steepness
#'
#' Plots the probability distribution function of steepness from the
#' mean and standard deviation.
#'
#' @param m The mean of the distribution (vectorized).
#' @param sd The standard deviation of the distribution (vectorized).
#' @param is_transform Logical, whether the mean and standard deviation are in normal space
#' (FALSE) or transformed space.
#' @param SR The stock recruitment relationship (determines the range and, if relevant, transformation of
#' steepness).
#' @param color A vector of colors.
#' @return A plot of the probability distribution function. Vertical dotted line
#' indicates mean of distribution.
#' @note The function samples from a beta distribution with parameters alpha and beta
#' that are converted from the mean and standard deviation. Then, the distribution is
#' transformed from 0 - 1 to 0.2 - 1.
#' @author Q. Huynh
#' @export
#' @seealso \code{\link{plot_lognormalvar}} \code{\link{plot_betavar}}
#' @examples
#' mu <- DLMtool::Simulation_1@@steep
#' stddev <- DLMtool::Simulation_1@@steep * DLMtool::Simulation_1@@CV_steep
#' plot_steepness(mu, stddev)
plot_steepness <- function(m, sd, is_transform = FALSE, SR = c("BH", "Ricker"), color = "black") {
  SR <- match.arg(SR)
  ncurve <- length(m)

  if(SR == "BH") {
    support <- seq(0.201, 0.999, 0.001)
    dist <- matrix(NA, nrow = length(support), ncol = ncurve)
    if(is_transform) {
      #f_y(y) = f_x(g-1(y)) * abs(d/dy[g-1(y)])
      #where f is the pdf of distribution, g(y) = 0.2 + 0.8/(1 + exp(-X)) is the transformation
      #y is steepness, x is a normal variable
      z <- (support - 0.2)/0.8
      for(i in 1:ncurve) {
        dist[, i] <- dnorm(logit(z), m[i], sd[i]) * (1/z + 1/(1-z)) * 1.25 * support
      }
      m <- ilogit(m) * 0.8 + 0.2
    } else {
      #f_y(y) = f_x(g-1(y)) * abs(d/dy[g-1(y)])
      #where f is the pdf of distribution, g(y) = 0.8X + 0.2 is the transformation
      #y is steepness, x is a beta variable
      m.transformed <- (m - 0.2)/0.8
      a <- alphaconv(m = m.transformed, sd = sd/0.8)
      b <- betaconv(m = m.transformed, sd = sd/0.8)

      for(i in 1:ncurve) {
        if(a[i] > 0 && b[i] > 0) dist[, i] <- dbeta((support - 0.2)/0.8, a[i], b[i]) / 0.8
        else {
          dist[, i] <- dnorm(support, m[i], sd[i])
          #warning("Transformation not possible with value of m and sd for steepness. Typically, sd is too high given the value of m, resulting in negative beta distribution parameters).")
        }
      }
    }
  }

  if(SR == "Ricker") {
    if(is_transform) {
      #f_y(y) = f_x(g-1(y)) * abs(d/dy[g-1(y)])
      #where f is the pdf of distribution, g(y) = 0.2 + 0.8/(1 + exp(-X)) is the transformation
      #y is steepness, x is a normal variable
      support.norm <- seq(min(m - 5*sd, na.rm = TRUE), max(m+5*sd, na.rm = TRUE),
                          length.out = 1e3)
      support <- exp(support.norm) + 0.2
      dist <- matrix(NA, nrow = length(support), ncol = ncurve)
      for(i in 1:ncurve) dist[, i] <- dnorm(support.norm, m[i], sd[i]) / (support - 0.2)

      dist[is.infinite(dist)] <- NA

      dist.max <- max(dist, na.rm = TRUE)
      tails <- apply(dist, 1, function(x) all(x < 0.001 * dist.max))
      tails <- which(!tails)
      ind.tails <- c(tails[1]:tails[length(tails)])

      support <- support[ind.tails]
      dist <- as.matrix(dist[ind.tails, ])
      m <- exp(m) + 0.2

    } else {
      #f_y(y) = f_x(g-1(y)) * abs(d/dy[g-1(y)])
      #where f is the pdf of distribution, g(y) = exp(X) + 0.2 is the transformation
      #y is steepness, x is a lognormal variable
      mulog <- mconv(m = m - 0.2, sd = sd)
      sdlog <- sdconv(m = m - 0.2, sd = sd)
      support <- seq(0.001, max(m + 5*sdlog), length.out = 1e3)

      dist <- matrix(NA, nrow = length(support), ncol = ncurve)
      for(i in 1:ncurve) dist[, i] <- dlnorm(support, mulog[i], sdlog[i])

      dist[is.infinite(dist)] <- NA

      dist.max <- max(dist, na.rm = TRUE)
      tails <- apply(dist, 1, function(x) all(x < 0.001 * dist.max))
      tails <- which(!tails)
      ind.tails <- c(tails[1]:tails[length(tails)])

      support <- support[ind.tails] + 0.2
      dist <- as.matrix(dist[ind.tails, ])
    }
  }

  plot(support, dist[, 1], typ = 'l', xlab = 'Steepness (h)',
       ylab = 'Probability density function', xlim = range(pretty(support)),
       ylim = c(0, 1.1 * max(dist, na.rm = TRUE)), col = color[1])
  if(ncurve > 1) {
    for(i in 2:ncurve) lines(support, dist[, i], col = color[i])
  }
  abline(h = 0, col = 'grey')
  abline(v = m, lty = 2, col = color)

  invisible()
}





#' Plot time series of data
#'
#' Plot time series of observed (with lognormally-distributed error bars) vs.
#' predicted data.
#'
#' @param Year A vector of years for the data.
#' @param obs A vector of observed data.
#' @param fit A vector of predicted data (e.g., from an assessment model).
#' @param obs_CV A vector of year-specific coefficient of variation in the observed data.
#' @param obs_CV_CI The confidence interval for the error bars based for \code{obs_CV}.
#' @param obs_upper A vector of year-specific upper bounds for the error bars of the observed data (in lieu of argument \code{obs_CV}).
#' @param obs_lower A vector of year-specific lower bounds for the error bars of the observed data (in lieu of argument \code{obs_CV}).
#' @param obs_ind_blue Indices of \code{obs} for which the plotted points and error bars will be blue.
#' @param fit_linewidth Argument \code{lwd} for fitted line.
#' @param fit_color Color of fitted line.
#' @param label Character string that describes the data to label the y-axis.
#' @author Q. Huynh
#' @seealso \code{\link{plot_residuals}}
#' @examples
#' data(Red_snapper)
#' plot_timeseries(Red_snapper@@Year, Red_snapper@@Cat[1, ],
#' obs_CV = Red_snapper@@CV_Cat, label = "Catch")
#' @export plot_timeseries
plot_timeseries <- function(Year, obs, fit = NULL, obs_CV = NULL, obs_CV_CI = 0.95,
                            obs_upper = NULL, obs_lower = NULL, obs_ind_blue = NULL, fit_linewidth = 3,
                            fit_color = "red", label = "Observed data") {
  old.warning <- options()$warn
  options(warn = -1)
  on.exit(options(warn = old.warning))

  # Without CV interval
  if(is.null(obs_CV)) {
    y.max <- max(c(obs, fit), na.rm = TRUE)
    if(is.null(obs_ind_blue)) {
      plot(Year, obs, typ = 'o', ylab = label, ylim = c(0, 1.1 * y.max))
    } else {
      plot(Year, obs, typ = 'n', ylab = label, ylim = c(0, 1.1 * y.max))
      lines(Year[-obs_ind_blue], obs[-obs_ind_blue], typ = 'o')
      lines(Year[obs_ind_blue], obs[obs_ind_blue], typ = 'o', col = "blue")
    }
  }

  # With CV interval
  if(!is.null(obs_CV) || (!is.null(obs_upper) & !is.null(obs_lower))) {
    sigma <- sdconv(1, obs_CV)
    if(is.null(obs_upper))
      obs_upper <- exp(log(obs) + qnorm(1-0.5*(1-obs_CV_CI)) * sigma)
    if(is.null(obs_lower))
      obs_lower <- exp(log(obs) + qnorm(0.5*(1-obs_CV_CI)) * sigma)
    y.max <- max(c(obs_lower, obs_upper, obs, fit), na.rm = TRUE)

    if(is.null(obs_ind_blue)) {
      plot(Year, obs, typ = 'o', ylab = label, ylim = c(0, 1.1 * y.max))
      arrows(Year, obs_lower, Year, obs_upper, length = 0.025, angle = 90, code = 3)
    } else {
      plot(Year, obs, typ = 'n', ylab = label, ylim = c(0, 1.1 * y.max))
      lines(Year[-obs_ind_blue], obs[-obs_ind_blue], typ = 'o')
      arrows(Year[-obs_ind_blue], obs_lower[-obs_ind_blue], Year[-obs_ind_blue],
             obs_upper[-obs_ind_blue], length = 0.025, angle = 90, code = 3, col = "grey30")

      lines(Year[obs_ind_blue], obs[obs_ind_blue], typ = 'o', col = "blue")
      arrows(Year[obs_ind_blue], obs_lower[obs_ind_blue], Year[obs_ind_blue],
             obs_upper[obs_ind_blue], length = 0.025, angle = 90, code = 3, col = "blue")
    }
  }
  if(!is.null(fit)) lines(Year, fit, lwd = fit_linewidth, col = fit_color)
  abline(h = 0, col = 'grey')

  invisible()
}

#' Plot residuals
#'
#' Plots figure of residuals (or any time series with predicted mean of zero).
#'
#' @param Year A vector of years for the data.
#' @param res A vector of residuals.
#' @param res_sd A vector of year specific standard deviation for \code{res}.
#' @param res_sd_CI The confidence interval for the error bars based for \code{res_sd}.
#' @param res_upper A vector of year-specific upper bounds for the error bars of the residual (in lieu of argument \code{res_CV}).
#' @param res_lower A vector of year-specific lower bounds for the error bars of the residual (in lieu of argument \code{res_CV}).
#' @param res_ind_blue Indices of \code{obs} for which the plotted residuals and error bars will be blue.
#' @param draw_zero Indicates whether a horizontal line should be drawn at zero.
#' @param zero_linetype Passes argument \code{lty} (e.g. solid line = 1, dotted = 2) to \code{draw_zero}.
#' @param label Character string that describes the data to label the y-axis.
#' @author Q. Huynh
#' @seealso \code{\link{plot_timeseries}}
#' @export plot_residuals
plot_residuals <- function(Year, res, res_sd = NULL, res_sd_CI = 0.95,
                           res_upper = NULL, res_lower = NULL, res_ind_blue = NULL, draw_zero = TRUE,
                           zero_linetype = 2, label = "Residual") {
  old.warning <- options()$warn
  options(warn = -1)
  on.exit(options(warn = old.warning))

  # Without sd interval
  if(is.null(res_sd)) {
    res.lim <- max(abs(res), na.rm = TRUE)
    if(is.null(res_ind_blue) || all(!res_ind_blue)) {
      plot(Year, res, typ = 'o', ylab = label, ylim = c(-res.lim, res.lim))
    } else {
      plot(Year, res, typ = 'n', ylab = label, ylim = c(-res.lim, res.lim))
      lines(Year[-res_ind_blue], res[-res_ind_blue], typ = 'o')
      lines(Year[res_ind_blue], res[res_ind_blue], typ = 'o', col = "blue")
    }
  }

  # With CV interval
  if(!is.null(res_sd) || (!is.null(res_upper) & !is.null(res_lower))) {
    if(is.null(res_upper)) res_upper <- res + qnorm(1-0.5*(1-res_sd_CI)) * res_sd
    if(is.null(res_lower)) res_lower <- res + qnorm(0.5*(1-res_sd_CI)) * res_sd
    res.lim <- max(abs(c(res_lower, res_upper, res)), na.rm = TRUE)

    if(is.null(res_ind_blue) || all(!res_ind_blue)) {
      plot(Year, res, typ = 'o', ylab = label, ylim = c(-res.lim, res.lim))
      arrows(Year, res_lower, Year, res_upper, length = 0.025, angle = 90,
             code = 3, col = 'grey30')
    } else {
      plot(Year, res, typ = 'n', ylab = label, ylim = c(-res.lim, res.lim))
      lines(Year[-res_ind_blue], res[-res_ind_blue], typ = 'o')
      arrows(Year[-res_ind_blue], res_lower[-res_ind_blue], Year[-res_ind_blue],
             res_upper[-res_ind_blue], length = 0.025, angle = 90, code = 3, col = 'grey30')

      lines(Year[res_ind_blue], res[res_ind_blue], typ = 'o', col = "blue")
      arrows(Year[res_ind_blue], res_lower[res_ind_blue], Year[res_ind_blue],
             res_upper[res_ind_blue], length = 0.025, angle = 90, code = 3, col = 'blue')
    }
  }
  if(draw_zero) abline(h = 0, lty = zero_linetype)
  invisible()
}




#' Plot composition data
#'
#' Plots annual length or age composition data.
#'
#' @param Year A vector of years.
#' @param obs A matrix of either length or age composition data. For lengths, rows and columns
#' should index years and length bin, respectively. For ages, rows and columns should index
#' years and age, respectively.
#' @param fit A matrix of predicted length or age composition from an assessment model.
#' Same dimensions as obs.
#' @param plot_type Indicates which plots to create. Options include annual distributions,
#' bubble plot of the data, and bubble plot of the residuals, and annual means.
#' @param N Annual sample sizes. Vector of length \code{nrow(obs)}.
#' @param CAL_bins A vector of lengths corresponding to the columns in \code{obs}.
#' and \code{fit}. Ignored for age data.
#' @param ages An optional vector of ages corresponding to the columns in \code{obs}.
#' @param ind A numeric vector for plotting a subset of rows (which indexes year) of \code{obs} and \code{fit}.
#' @param annual_ylab Character string for y-axis label when \code{plot_type = "annual"}.
#' @param annual_yscale For annual composition plots (\code{plot_type = "annual"}), whether the raw values
#' ("raw") or frequencies ("proportions") are plotted.
#' @param bubble_adj Numeric, for adjusting the relative size of bubbles in bubble plots
#' (larger number = larger bubbles).
#' @param fit_linewidth Argument \code{lwd} for fitted line.
#' @param fit_color Color of fitted line.
#' @param bubble_color Colors for negative and positive residuals, respectively, for bubble plots.
#' @return Plots depending on \code{plot_type}.
#' @author Q. Huynh
#' @export plot_composition
#' @examples
#' \donttest{
#' data(Red_snapper)
#' plot_composition(obs = Red_snapper@@CAA[1, , ], plot_type = "annual")
#' plot_composition(obs = Red_snapper@@CAA[1, , ], plot_type = "bubble_data")
#'
#' plot_composition(obs = Red_snapper@@CAL[1, , ], plot_type = "annual", Red_snapper@@CAL_bins[1:43])
#' plot_composition(obs = Red_snapper@@CAL[1, , ], plot_type = "bubble_data",
#' CAL_bins = Red_snapper@@CAL_bins[1:43])
#' }
plot_composition <- function(Year = 1:nrow(obs), obs, fit = NULL, plot_type = c('annual', 'bubble_data', 'bubble_residuals', 'mean'),
                             N = rowSums(obs), CAL_bins = NULL, ages = NULL, ind = 1:nrow(obs),
                             annual_ylab = "Frequency", annual_yscale = c("proportions", "raw"),
                             bubble_adj = 5, bubble_color = c("black", "white"), fit_linewidth = 3, fit_color = "red") {
  plot_type <- match.arg(plot_type)
  annual_yscale <- match.arg(annual_yscale)
  if(is.null(CAL_bins)) data_type <- "age" else data_type <- "length"

  if(!is.null(fit) && !all(dim(fit) == dim(obs))) stop("Dimensions of 'obs' and 'fit' do not match.")

  if(data_type == 'length') {
    data_val <- CAL_bins
    data_lab <- "Length"
  }
  if(data_type == 'age') {
    data_val <- if(is.null(ages)) 1:ncol(obs) else ages
    data_lab <- "Age"
  }
  if(!is.null(N)) N <- round(N, 1)

  if(annual_yscale == "proportions") {
    obs_prob_all <- obs/rowSums(obs, na.rm = TRUE)
    if(!is.null(fit)) fit_prob_all <- fit/rowSums(fit, na.rm = TRUE)
    else fit_prob_all <- NULL
  }
  if(annual_yscale == "raw") {
    obs_prob_all <- obs
    if(!is.null(fit)) fit_prob_all <- fit
    else fit_prob_all <- NULL
  }

  # subset
  #ind <- rowSums(obs, na.rm = TRUE) > 0
  Year <- Year[ind]
  obs <- obs[ind, , drop = FALSE]
  if(!is.null(fit)) fit <- fit[ind, , drop = FALSE]
  if(!is.null(N)) N <- N[ind]

  # Bubble plot (obs)
  if('bubble_data' %in% plot_type) {
    range_obs <- pretty(obs, n = 6)
    n1 <- range_obs[2]
    n2 <- pretty(quantile(obs[obs > 0], na.rm = TRUE, probs = 0.9))[2]
    if(n2 < n1) n1 <- 0.5 * n2
    diameter_max <- bubble_adj / n2
    plot(NULL, NULL, typ = 'n', xlim = range(Year), xlab = "Year",
         ylim = range(data_val), ylab = data_lab)
    for(i in 1:length(Year)) {
      for(j in 1:length(data_val)) {
        points(Year[i], data_val[j], cex = 0.5 * diameter_max * pmin(obs[i, j], n2), pch = 21, bg = "white")
      }
    }
    legend("topleft", legend = c(n1, paste0(">", n2)), pt.cex = 0.5 * diameter_max * c(n1, n2),
           pt.bg = "white", pch = 21, horiz = TRUE)
    return(invisible())
  }
  # Bubble plot (residuals if applicable)
  if('bubble_residuals' %in% plot_type) {
    if(is.null(fit)) stop("No fitted data available.")

    obs_prob <- obs/rowSums(obs, na.rm = TRUE)
    fit_prob <- fit/rowSums(fit, na.rm = TRUE)

    resid <- N * (obs_prob - fit_prob) / sqrt(N * fit_prob * (1 - fit_prob))
    diameter_max <- bubble_adj / pmin(10, max(abs(resid), na.rm = TRUE))
    plot(NULL, NULL, typ = 'n', xlim = range(Year), xlab = "Year", ylim = range(data_val), ylab = data_lab)

    Year_mat <- matrix(Year, ncol = ncol(resid), nrow = nrow(resid))
    data_mat <- matrix(data_val, ncol = ncol(resid), nrow = nrow(resid), byrow = TRUE)
    isPositive <- resid > 0
    points(Year_mat[!isPositive], data_mat[!isPositive], cex = pmin(0.5 * diameter_max * abs(resid[!isPositive]), diameter_max), pch = 21, bg = bubble_color[1])
    points(Year_mat[isPositive], data_mat[isPositive], cex = pmin(0.5 * diameter_max * resid[isPositive], diameter_max), pch = 21, bg = bubble_color[2])
    legend("topleft", legend = c("<-10", "-1", "1", ">10"),
           pt.cex = c(diameter_max, 0.5 * diameter_max, 0.5 * diameter_max, diameter_max),
           pt.bg = rep(bubble_color, each = 2), pch = 21, horiz = TRUE)

    return(invisible())
  }
  # Mean length or age over time
  if('mean' %in% plot_type) {
    mu <- mupred <- numeric(length = length(Year))
    for(i in 1:length(mu)) {
      mu[i] <- weighted.mean(data_val, obs[i, ], na.rm = TRUE)
      if(!is.null(fit)) mupred[i] <- weighted.mean(data_val, fit[i, ], na.rm = TRUE)
    }
    ind2 <- (which(!is.na(mu) & mu > 0)[1]):length(mu)
    plot(Year[ind2], mu[ind2], xlab = "Year", ylab = paste0("Mean ", data_type), typ = "o")
    if(!is.null(fit)) lines(Year[ind2], mupred[ind2], lwd = fit_linewidth, col = fit_color)

    return(invisible())
  }

  # Annual comps (obs vs. fitted if available)
  if("annual" %in% plot_type) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mfcol = c(4, 4), mar = rep(0, 4), oma = c(5.1, 5.1, 2.1, 2.1))
    ylim <- c(0, 1.1 * max(obs_prob_all, fit_prob_all, na.rm = TRUE))

    if(annual_yscale == "proportions") {
      obs_prob <- obs/rowSums(obs, na.rm = TRUE)
      if(!is.null(fit)) fit_prob <- fit/rowSums(fit, na.rm = TRUE)
      else fit_prob <- NULL
    }
    if(annual_yscale == "raw") {
      obs_prob <- obs
      fit_prob <- fit
    }

    yaxp <- c(0, max(pretty(ylim, n = 4)), 4)
    if(max(obs_prob_all, fit_prob_all, na.rm = TRUE) == 1) yaxp <- c(0, 1, 4)

    las <- 1
    for(i in 1:length(Year)) {
      yaxt <- ifelse(i %% 16 %in% c(1:4), "s", "n") # TRUE = first column
      xaxt <- ifelse(i < length(Year) & i %% 4 %in% c(1:3), "n", "s") # TRUE = first three rows

      plot(data_val, obs_prob[i, ], typ = "n", ylim = ylim, yaxp = yaxp, xaxt = xaxt, yaxt = yaxt, las = las)
      abline(h = 0, col = "grey")
      lines(data_val, obs_prob[i, ], typ = "o")
      if(!is.null(fit)) lines(data_val, fit_prob[i, ], lwd = fit_linewidth, col = fit_color)
      legend("topright", legend = c(Year[i], ifelse(is.null(N), "", paste0("N = ", N[i]))), bty = "n", xjust = 1)

      if(i %% 16 == 1) {
        mtext(data_lab, side = 1, line = 3, outer = TRUE)
        mtext(annual_ylab, side = 2, line = 3.5, outer = TRUE)
      }
    }
    return(invisible())

  }

  invisible()
}







#' @importFrom graphics arrows
plot_surplus_production <- function(B, B0 = NULL, C, yield_fn = NULL, arrow_size = 0.07, xlab = NULL) {
  old.warning <- options()$warn
  options(warn = -1)
  on.exit(options(warn = old.warning))

  if(!is.null(B0)) {
    B <- B/B0
    if(is.null(xlab)) xlab <- expression(B/B[0])
  } else {
    if(is.null(xlab)) xlab <- "Biomass"
  }
  B_now <- B[1:(length(B)-1)]
  B_next <- B[2:length(B)]
  SP_now <- B_next - B_now + C

  xlim <- c(0, max(B))
  ylim <- c(min(0, min(SP_now)), max(SP_now))
  plot(B_now, SP_now, typ = "n", xlab = xlab, xlim = xlim, ylim = ylim,
       ylab = "Surplus production")
  if(!is.null(yield_fn)) {
    if(!is.null(B0)) lines(yield_fn$B_B0, yield_fn$Yield) else {
      lines(yield_fn$B, yield_fn$Yield)
    }
  }
  arrows(x0 = B_now, y0 = SP_now[1:(length(B)-1)], x1 = B_next, y1 = SP_now[2:length(B)],
         length = arrow_size)
  abline(h = 0, col = "grey")

  invisible()
}


plot_Kobe <- function(biomass, exploit, arrow_size = 0.07, color = TRUE, xlab = expression(B/B[MSY]), ylab = expression(F/F[MSY])) {
  old.warning <- options()$warn
  on.exit(options(warn = old.warning))
  options(warn = -1)

  n.arrows <- length(exploit)
  if(length(biomass) > n.arrows) biomass <- biomass[1:n.arrows]

  x.max <- max(biomass, 1)
  y.max <- max(exploit, 1)
  plot(NULL, NULL, typ = "n", xlab = xlab, ylab = ylab, xlim = c(0, max(1.1, 1.1 * x.max)), ylim = c(0, max(1.1, 1.1 * y.max)))
  if(color) {
    # Colors from https://www.rapidtables.com/web/color/html-color-codes.html
    green <- "#228B22"    #forestgreen
    yellow <- "#F0E68C"   #khaki
    red <- "#CD5C5C"      #indianred
    polygon(x = c(1, 1, 10*x.max, 10*x.max), y = c(1, -5, -5, 1), col = green, border = NA)
    polygon(x = c(-5, -5, 1, 1), y = c(1, -5, -5, 1), col = yellow, border = NA)
    polygon(x = c(1, 1, 10*x.max, 10*x.max), y = c(10*y.max, 1, 1, 10*y.max),
            col = yellow, border = NA)
    polygon(x = c(-5, -5, 1, 1), y = c(10*y.max, 1, 1, 10*y.max),
            col = red, border = NA)
    box()
  }
  arrows(x0 = biomass[1:(n.arrows-1)], y0 = exploit[1:(n.arrows-1)],
         x1 = biomass[2:n.arrows], y1 = exploit[2:n.arrows], length = arrow_size)
  abline(h = 0, col = "grey")
  abline(v = 0, col = "grey")
  if(!color) {
    abline(h = 1, lty = 2)
    abline(v = 1, lty = 2)
  }
  invisible()
}


#' Plot stock-recruitment function
#'
#' Plot stock-recruitment (with recruitment deviations if estimated).
#'
#' @param Spawners A vector of the number of the spawners (x-axis).
#' @param expectedR A vector of the expected recruitment (from the
#' stock-recruit function) corresponding to values of \code{Spawners}.
#' @param R0 Virgin recruitment.
#' @param S0 Virgin spawners.
#' @param rec_dev If recruitment deviations are estimated, a vector of estimated recruitment
#' (in normal space) corresponding to values of \code{Spawners}.
#' @param trajectory Indicates whether arrows will be drawn showing the trajectory of
#' spawners and recruitment deviations over time.
#' @param y_zoom If recruitment deviations are plotted, the y-axis limit relative to
#' maximum expected recruitment \code{expectedR}. If \code{NULL}, all recruitments are plotted.
#' @param ylab Character string for label on y-axis.
#' @author Q. Huynh
#' @return A stock-recruit plot
#' @export plot_SR
plot_SR <- function(Spawners, expectedR, R0 = NULL, S0 = NULL, rec_dev = NULL, trajectory = FALSE,
                    y_zoom = NULL, ylab = "Recruitment") {
  if(is.null(rec_dev)) {
    R.max <- 1.1 * max(expectedR)
  } else {
    if(is.null(y_zoom)) R.max <- 1.1 * max(rec_dev)
    else R.max <- y_zoom * max(expectedR)
  }
  S.max <- 1.1 * max(c(Spawners, S0))
  plot(Spawners[order(Spawners)], expectedR[order(Spawners)], typ = "l", xlim = c(0, 1.05 * S.max), ylim = c(0, 1.1 * R.max),
       xlab = "Spawning Stock Biomass (SSB)", ylab = ylab)
  if(!trajectory) {
    if(is.null(rec_dev)) points(Spawners, expectedR)
    if(!is.null(rec_dev)) points(Spawners, rec_dev)
  }
  if(trajectory) {
    old.warning <- options()$warn
    on.exit(options(warn = old.warning), add = TRUE)
    options(warn = -1)
    n.arrows <- length(Spawners)

    arrows(x0 = Spawners[1:(n.arrows-1)], y0 = rec_dev[1:(n.arrows-1)],
           x1 = Spawners[2:n.arrows], y1 = rec_dev[2:n.arrows], length = 0.07)
  }
  if(!is.null(R0) && !is.null(S0)) points(S0, R0, col = "red", pch = 16)
  abline(h = 0, col = "grey")
  abline(v = 0, col = "grey")
}


plot_generic_at_age <- function(Age, quantity, label, ymax = 1.1 * max(quantity)) {
  plot(Age, quantity, ylab = label, typ = "n", ylim = c(0, ymax))
  abline(h = 0, col = "grey")
  lines(Age, quantity, typ = "o")

  invisible()

}

plot_ogive <- function(Age, ogive, label = "Selectivity") {
  plot_generic_at_age(Age = Age, quantity = ogive, label = label, ymax = 1.1)

  invisible()
}


calculate_Mohn_rho <- function(ts, est = NULL, ts_lab, est_lab = NULL) {
  rho_ts <- apply(ts, 3, Mohn_rho)
  rho_est <- if(!is.null(est)) apply(est, 2, Mohn_rho, type = "est") else NULL

  ans <- matrix(c(rho_ts, rho_est), ncol = 1)
  dimnames(ans) <- list(c(ts_lab, est_lab), "Mohn's rho")
  return(round(ans, 3))
}


Mohn_rho <- function(x, type = c("ts", "est")) {
  type <- match.arg(type)
  if(type == "ts") { # let x be a matrix of nrow = n_peel + 1 and ncol = nyear + 1
    terminal_ind <- apply(x[-1, , drop = FALSE], 1, function(y) sum(!is.na(y)))
    n_peel <- length(terminal_ind)
    rho <- diag(x[1:n_peel + 1, terminal_ind, drop = FALSE])/x[1, terminal_ind] - 1
  } else { # let x be a vector of length n_peel + 1
    n_peel <- length(x) - 1
    rho <- x[1:n_peel + 1]/x[1] - 1
  }
  return(mean(rho))
}
