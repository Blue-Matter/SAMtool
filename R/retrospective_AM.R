#' retrospective_AM (retrospective of Assessment model in MSE)
#'
#' Plots the true retrospective of an assessment model during the closed-loop simulation. A series of time series estimates of SSB, F, and VB
#' are plotted over the course of the MSE are plotted against the operating model (true) values (in black).
#'
#' @param MSE An object of class \linkS4class{MSE}.
#' @param MP Character. The name of the management procedure created by \code{\link{make_MP}} containing the asssessment model. 
#' @param sim Integer between 1 and MSE@@nsim. The simulation number for which the retrospectives will be plotted.
#' @param plot_legend Logical. Whether to plot legend to reference year of assessment in the MSE.
#' @author Q. Huynh
#' @details For assessment models that utilize annual harvest rates (u), the instantaneous fishing mortality rates
#' are obtained as F = -log(1 - u).
#' @note This function only plots retrospectives from a single simulation in the MSE. Results from one figure
#' may not be indicative of general assessment behavior and performance overall.
#' @return A series of figures for spawning stock biomass
#' (SSB, including absolute magnitude and relative to MSY and unfished), fishing mortality (F, including absolute
#' magnitude and relative to MSY), and vulnerable biomass (VB) estimates over the course of the MSE are plotted
#' against the operating model (true) values (in black).
#' @examples
#' \dontrun{
#' DD_MSY <- makeMP(DD_TMB, HCR_MSY, diagnostic = "full")
#' myMSE <- MSEtool::runMSE(OM = MSEtool::testOM, MPs = "DD_MSY")
#' retrospective_AM(myMSE, MP = "DD_MSY", sim = 1)
#' }
#' @seealso \link{diagnostic}
#' @importFrom gplots rich.colors
#' @export
retrospective_AM <- function(MSE, MP, sim = 1, plot_legend = FALSE) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  if(!inherits(MSE, "MSE")) stop("No object of class MSE was provided.")
  if(length(sim) > 1 || sim > MSE@nsim) stop(paste0(sim, " should be a number between 1 and ", MSE@nsim, "."))

  MPs <- MSE@MPs
  match_ind <- match(MP, MPs)
  if(is.na(match_ind)) stop(paste(MP, "MP was not found in the MSE object. Available options are:", paste(MPs, collapse = " ")))

  has_Assess_fn <- function(Data) {
    Misc <- Data@Misc[1:MSE@nsim]
    all(vapply(Misc, function(y) any(names(y) == "Assessment_report"), logical(1)))
  }

  has_Assess <- has_Assess_fn(MSE@PPD[[match_ind]])
  if(!has_Assess) stop("No Assessment objects were found in MSE@PPD for ", MP, ". Use an MP created by: make_MP(diagnostic = \"full\").")
  Assessment_report <- lapply(MSE@PPD[[match_ind]]@Misc[1:MSE@nsim], getElement, "Assessment_report")[[sim]]
  
  color.vec <- rich.colors(length(Assessment_report))
  Yr_MSE <- 1:(MSE@nyears + MSE@proyears)
  Assess <- lapply(Assessment_report, slot, "VB")
  End_Assess_Yr <- vapply(Assess, function(x) as.numeric(names(x))[length(x)-1], numeric(1))

  plot_type = c("SSB", "F", "SSB_SSBMSY", "F_FMSY", "SSB_SSB0", "VB")

  par(mfcol = c(2, 3), mar = c(5, 4, 1, 1), oma = c(0, 0, 2.5, 0))
  for(i in 1:length(plot_type)) {
    if(plot_type[i] == "SSB_SSBMSY") {
      ylab <- expression(SSB/SSB[MSY])

      Hist <- apply(MSE@Hist@TSdata$SBiomass[sim, , ], 1, sum)/MSE@OM$SSBMSY[sim]
      Proj <- MSE@SB_SBMSY[sim, match_ind, ]
      Assess <- lapply(Assessment_report, slot, "SSB_SSBMSY")
    }
    if(plot_type[i] == "F_FMSY") {
      ylab <- expression(F/F[MSY])
      Hist <- apply(MSE@Hist@AtAge$F.Mortality[sim, , , ], 2, max)/MSE@OM$FMSY[sim]
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
      Hist <- apply(MSE@Hist@TSdata$SBiomass[sim, , ], 1, sum)
      Proj <- MSE@SSB[sim, match_ind, ]
      Assess <- lapply(Assessment_report, slot, "SSB")
    }
    if(plot_type[i] == "F") {
      ylab <- "F"
      Hist <- apply(MSE@Hist@AtAge$F.Mortality[sim, , , ], 2, max)
      Proj <- MSE@FM[sim, match_ind, ]
      Assess <- lapply(Assessment_report, slot, "FMort")
      if(length(do.call(c, Assess)) == 0) {
        AssessU <- lapply(Assessment_report, slot, "U")
        Assess <- lapply(AssessU, function(x) -log(1 - x))
      }
    }
    if(plot_type[i] == "SSB_SSB0") {
      ylab <- expression(SSB/SSB[0])
      Hist <- apply(MSE@Hist@TSdata$SBiomass[sim, , ], 1, sum)/MSE@OM$SSB0[sim]
      Proj <- MSE@SSB[sim, match_ind, ]/MSE@OM$SSB0[sim]
      Assess <- lapply(Assessment_report, slot, "SSB_SSB0")
    }
    if(plot_type[i] == "VB") {
      ylab <- "Vulnerable biomass"
      Hist <- apply(MSE@Hist@TSdata$VBiomass[sim, , ], 1, sum)
      Proj <- MSE@VB[sim, match_ind, ]
      Assess <- lapply(Assessment_report, slot, "VB")
    }
    #if(plot_type[i] == "Recruit") {
    #  ylab <- "Recruitment"
    #  if(!is.null(Hist)) {
    #    Hist_ts <- Hist@AtAge$Number[sim, 1, , ] %>% rowSums()
    #  } else {
    #    Hist_ts <- rep(NA, MSE@nyears)
    #    message("Provide Hist object in order to plot simulated recruitment in the operating model.")
    #  }
    #  Proj <- rep(NA, MSE@proyears)
    #  Assess <- lapply(Assessment_report, slot, "R")
    #}
    
    Assess_Yr <- lapply(Assess, function(x) as.numeric(names(x)))
    xlimits <- c(1, MSE@nyears + MSE@proyears)
    
    converged_assessments <- vapply(Assessment_report, slot, logical(1), "conv")
    ylimits <- c(0, 1.1 * max(c(Hist, Proj, do.call(c, Assess[converged_assessments])), na.rm = TRUE))
    
    if(all(!is.na(ylimits))) {
      plot(NULL, NULL, xlab = "MSE year", ylab = ylab, xlim = xlimits, ylim = ylimits)
      for(j in length(Assess):1) {
        if(Assessment_report[[j]]@conv && length(Assess[[j]]) > 0) {
          lines(Assess_Yr[[j]], Assess[[j]], col = color.vec[j])
          points(max(Assess_Yr[[j]]), Assess[[j]][length(Assess[[j]])], col = color.vec[j], pch = 16, cex = 1.5)
        }
      }
      if(plot_type[i] == "SSB_SSBMSY" || plot_type[i] == "F_FMSY") abline(h = 1)
      lines(Yr_MSE, c(Hist, Proj), xlab = "MSE year", ylab = ylab, xlim = xlimits, ylim = ylimits, lwd = 2)
      abline(h = 0, col = "grey")
      abline(v = MSE@nyears, lty = 2)
      if(plot_legend && i == 1) legend("topleft", c("OM", End_Assess_Yr), col = c("black", color.vec), lwd = c(2, rep(1, length(Assessment_report))))
    } else {
      message(paste0("Skipped plot for ", plot_type[i], "."))
    }
  }
  title(paste0("Simulation #", sim, " of ", MP, "\n(OM in black, MP in colors)"), outer = TRUE)
  invisible()
}
