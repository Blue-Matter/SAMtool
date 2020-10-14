#' retrospective_AM (retrospective of Assessment model in MSE)
#'
#' Plots the true retrospective of an assessment model during the MSE. A series of time series estimates of SSB, F, and VB
#' are plotted over the course of the MSE are plotted against the operating model (true) values (in black).
#'
#' @param MSE An object of class MSE created by \code{\link[DLMtool]{runMSE}} with \code{PPD = TRUE}.
#' @param sim Integer between 1 and MSE@@nsim. The simulation number for which the retrospectives will be plotted.
#' @param MP Character. The name of the management procedure created by \code{\link{make_MP}} containing the asssessment model.
#' @param MSE_Hist Optional. The list containing historical data for the MSE, created by \code{\link[DLMtool]{runMSE}} with argument \code{Hist = TRUE}.
#' Currently only used to plot operating model vulnerable biomass in historical period.
#' @param plot_legend Logical. Whether to plot legend to reference year of assessment in the MSE.
#' @author Q. Huynh
#' @details For assessment models that utilize annual harvest rates (u), the instantaneous fishing mortality rates
#' are obtained as F = -log(1 - u).
#' @note This function only plots retrospectives from a single simulation in the MSE. Results from one figure
#' may not be indicative of general assessment behavior and performance overall.
#'
#' For \link{SP} and \link{SP_SS} assessment models don't model SSB. Instead, the estimated vulnerable biomass is plotted.
#' @return A series of figures for spawning stock biomass
#' (SSB, including absolute magnitude and relative to MSY and virgin), fishing mortality (F, including absolute
#' magnitude and relative to MSY), and vulnerable biomass (VB) estimates over the course of the MSE are plotted
#' against the operating model (true) values (in black).
#' @examples
#' \dontrun{
#' DD_MSY <- makeMP(DD_TMB, HCR_MSY, diagnostic = "full")
#' myMSE_hist <- DLMtool::runMSE(DLMtool::testOM, Hist = TRUE)
#' myMSE <- DLMtool::runMSE(DLMtool::testOM, MPs = "DD_MSY", PPD = TRUE)
#' retrospective_AM(myMSE, sim = 1, MP = "DD_MSY")
#' retrospective_AM(myMSE, sim = 1, MP = "DD_MSY", Hist = myMSE_hist)
#' }
#' @seealso \link{diagnostic_AM}
#' @importFrom gplots rich.colors
#' @export
retrospective_AM <- function(MSE, sim = 1, MP, MSE_Hist = NULL, plot_legend = FALSE) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  if(!inherits(MSE, "MSE")) stop("No object of class MSE was provided.")
  if(packageVersion("DLMtool") >= 5.3) {
    if(length(MSE@Misc$Data) == 0) stop("Nothing found in MSE@Misc$Data. Use an MP created by 'make_MP(diagnostic = 'full')' and set 'runMSE(PPD = TRUE)'.")
  } else {
    if(length(MSE@Misc) == 0) stop("Nothing found in MSE@Misc. Use an MP created by 'make_MP(diagnostic = 'full')' and set 'runMSE(PPD = TRUE)'.")
  }
  if(length(sim) > 1 || sim > MSE@nsim) stop(paste0(sim, " should be a number between 1 and ", MSE@nsim, "."))

  MPs <- MSE@MPs
  match_ind <- match(MP, MPs)
  if(is.na(match_ind)) stop(paste(MP, "MP was not found in the MSE object. Available options are:", paste(MPs, collapse = " ")))

  has_Assess_fn <- function(Data) {
    Misc <- Data@Misc
    all(vapply(Misc, function(y) any(names(y) == "Assessment_report"), logical(1)))
  }

  if(packageVersion("DLMtool") >= 5.3) {
    has_Assess <- has_Assess_fn(MSE@Misc$Data[[match_ind]])
    if(!has_Assess) stop("No Assessment objects were found in MSE@Misc$Data for any MP. Use an MP created by 'make_MP(diagnostic = 'full')' and set 'runMSE(PPD = TRUE)'.")
    Assessment_report <- lapply(MSE@Misc$Data[[match_ind]]@Misc, getElement, "Assessment_report")[[sim]]
  } else {
    has_Assess <- has_Assess_fn(MSE@Misc[[match_ind]])
    if(!has_Assess) stop("No Assessment objects were found in MSE@Misc for any MP. Use an MP created by 'make_MP(diagnostic = 'full')' and set 'runMSE(PPD = TRUE)'.")
    Assessment_report <- lapply(MSE@Misc[[match_ind]]@Misc, getElement, "Assessment_report")[[sim]]
  }

  isSP <- grepl("SP", Assessment_report[[1]]@Model)

  color.vec <- rich.colors(length(Assessment_report))
  Yr_MSE <- 1:(MSE@nyears + MSE@proyears)
  Assess <- lapply(Assessment_report, slot, "VB")
  End_Assess_Yr <- vapply(Assess, function(x) as.numeric(names(x))[length(x)-1], numeric(1))

  plot_type = c("SSB", "F", "SSB_SSBMSY", "F_FMSY", "SSB_SSB0", "VB")

  par(mfcol = c(2, 3), mar = c(5, 4, 1, 1), oma = c(0, 0, 2.5, 0))
  for(i in 1:length(plot_type)) {
    if(plot_type[i] == "SSB_SSBMSY") {
      ylab <- expression(SSB/SSB[MSY])

      Hist <- apply(MSE@SSB_hist, c(1, 3), sum)[sim, ]/MSE@OM$SSBMSY[sim]
      Proj <- MSE@B_BMSY[sim, match_ind, ]
      if(!isSP) {
        Assess <- lapply(Assessment_report, slot, "SSB_SSBMSY")
      } else {
        Assess <- lapply(Assessment_report, slot, "VB_VBMSY")
      }
    }
    if(plot_type[i] == "F_FMSY") {
      ylab <- expression(F/F[MSY])
      Hist <- apply(MSE@FM_hist, c(1, 3), max)[sim, ]/MSE@OM$FMSY[sim]
      Proj <- MSE@F_FMSY[sim, match_ind, ]
      Assess <- lapply(Assessment_report, slot, "F_FMSY")
      if(length(do.call(c, Assess)) == 0) {
        AssessU <- lapply(Assessment_report, slot, "U")
        AssessUMSY <- lapply(Assessment_report, slot, "UMSY")
        Assess <- Map(function(x, y) log(1-x)/log(1-y), x = AssessU, y = AssessUMSY)
      }
    }
    if(plot_type[i] == "SSB") {
      ylab <- "SSB"
      Hist <- apply(MSE@SSB_hist, c(1, 3), sum)[sim, ]
      Proj <- MSE@SSB[sim, match_ind, ]
      if(!isSP) {
        Assess <- lapply(Assessment_report, slot, "SSB")
      } else {
        Assess <- lapply(Assessment_report, slot, "VB")
      }
    }
    if(plot_type[i] == "F") {
      ylab <- "F"
      Hist <- apply(MSE@FM_hist, c(1, 3), max)[sim, ]
      Proj <- MSE@FM[sim, match_ind, ]
      Assess <- lapply(Assessment_report, slot, "FMort")
      if(length(do.call(c, Assess)) == 0) {
        AssessU <- lapply(Assessment_report, slot, "U")
        Assess <- lapply(AssessU, function(x) -log(1 - x))
      }
    }
    if(plot_type[i] == "SSB_SSB0") {
      ylab <- expression(SSB/SSB[0])
      Hist <- apply(MSE@SSB_hist, c(1, 3), sum)[sim, ]/MSE@OM$SSB0[sim]
      Proj <- MSE@SSB[sim, match_ind, ]/MSE@OM$SSB0[sim]
      if(isSP) {
        Assess <- lapply(Assessment_report, slot, "VB_VB0")
      } else {
        Assess <- lapply(Assessment_report, slot, "SSB_SSB0")
      }
    }
    if(plot_type[i] == "VB") {
      ylab <- "Vulnerable biomass"
      if(!is.null(MSE_Hist)) {
        Hist <- MSE_Hist$TSdata$VB[, sim]
      } else {
        Hist <- rep(NA, MSE@nyears)
        message("Provide Hist object in order to plot the historical vulnerable biomass in the operating model.")
      }
      Proj <- MSE@VB[sim, match_ind, ]
      Assess <- lapply(Assessment_report, slot, "VB")
    }

    if(plot_type[i] == "Recruit") {
      ylab <- "Recruitment"
      if(!is.null(MSE_Hist)) {
        Hist <- MSE_Hist$AtAge$Nage[sim, 1, ]
      } else {
        Hist <- rep(NA, MSE@nyears)
        message("Provide Hist object in order to plot simulated recruitment in the operating model.")
      }
      Proj <- rep(NA, MSE@proyears)
      Assess <- lapply(Assessment_report, slot, "R")
    }

    if(plot_type[i] != "MSY") {
      Assess_Yr <- lapply(Assess, function(x) as.numeric(names(x)))
      xlimits <- c(1, MSE@nyears + MSE@proyears)

      converged_assessments <- vapply(Assessment_report, slot, logical(1), "conv")
      ylimits <- c(0, 1.1 * max(c(Hist, Proj, do.call(c, Assess[converged_assessments])), na.rm = TRUE))

      if(all(!is.na(ylimits))) {
        plot(Yr_MSE, c(Hist, Proj), xlab = "MSE year", ylab = ylab, xlim = xlimits, ylim = ylimits, lwd = 2, typ = 'l')
        for(j in length(Assess):1) {
          if(Assessment_report[[j]]@conv && length(Assess[[j]]) > 0) lines(Assess_Yr[[j]], Assess[[j]], col = color.vec[j])
        }
        if(plot_type[i] == "SSB_SSBMSY" || plot_type[i] == "F_FMSY") abline(h = 1)
        abline(h = 0, col = "grey")
        abline(v = MSE@nyears, lty = 2)
        if(plot_legend && i == 1) legend("topleft", c("OM", End_Assess_Yr), col = c("black", color.vec), lwd = c(2, rep(1, length(Assessment_report))))
      } else {
        message(paste0("Skipped plot for ", plot_type[i], "."))
      }
    }

    if(plot_type[i] == "MSY") {
      Hist <- c(MSE@OM$FMSY[sim], MSE@OM$MSY[sim])
      Assess <- vapply(Assessment_report, function(x) c(-log(1 - slot(x, "UMSY")), slot(x, "MSY")), numeric(2))

      xc <- c(Hist[1], Assess[1, ])
      xc <- xc[xc < 2]
      xlimits <- c(max(0, mean(xc) - 2*sd(xc)), mean(xc + 2*sd(xc))) # FMSY
      yc <- c(Hist[2], Assess[2, ])
      ylimits <- c(max(0, mean(yc) - 2*sd(yc)), mean(yc + 2*sd(yc))) # MSY
      plot(x = Assess[1, ], y = Assess[2, ], xlim = xlimits, ylim = ylimits,
           xlab = "FMSY estimates", ylab = "MSY estimates", col = color.vec, pch = 16,
           cex = 2)
      j <- ncol(Assess) - 1
      arrows(x0 = Assess[1, 1:j], y0 = Assess[2, 1:j], x1 = Assess[1, 2:(j+1)], y1 = Assess[2, 2:(j+1)],
             length = 0.1)
      points(Hist[1], Hist[2], col = "grey", cex = 2, pch = 0)
    }

  }
  title(paste0(MP, " management procedure \n Simulation #", sim), outer = TRUE)
  invisible()
}
