
#' Class-\code{RCModel}
#'
#' An S4 class for the output from \link{RCM}.
#'
#' @name RCModel-class
#' @docType class
#'
#' @slot OM An updated operating model, class \linkS4class{OM}.
#' @slot SSB A matrix of estimated spawning biomass with \code{OM@@nsim} rows and \code{OM@@nyears+1} columns.
#' @slot NAA An array for the predicted numbers at age with dimension \code{OM@@nsim}, \code{OM@@nyears+1}, and \code{OM@@maxage+1}.
#' @slot CAA An array for the predicted catch at age with dimension \code{OM@@nsim}, \code{OM@@nyears}, \code{OM@@maxage}, and nfleet.
#' @slot CAL An array for the predicted catch at length with dimension \code{OM@@nsim}, \code{OM@@nyears}, length bins, and nfleet.
#' @slot conv A logical vector of length \code{OM@@nsim} indicating convergence of the RCM in the i-th simulation.
#' @slot Misc A list of length \code{OM@@nsim} with more output from the fitted RCM. Within each simulation, a named list containing items of interest include:
#'
#' \itemize{
#' \item B - total biomass - vector of length nyears+1
#' \item EPR0 - annual unfished spawners per recruit - vector of length nyears
#' \item ageM - age of 50\% maturity - integer
#' \item EPR0_SR - unfished spawners per recruit for the stock-recruit relationship (mean EPR0 over the first \code{ageM} years) - numeric
#' \item R0 - unfished recruitment for the stock-recruit relationship - numeric
#' \item h - steepness for the stock-recruit relationship - numeric
#' \item Arec - stock-recruit alpha parameter - numeric
#' \item Brec - stock-recruit beta parameter - numeric
#' \item E0_SR - unfished spawning biomass for the stock-recruit relationship (product of EPR0_SR and R0) - numeric
#' \item CR_SR - compensation ratio, the product of Arec and EPR0_SR - numeric
#' \item E0 - annual unfished spawning biomass (intersection of stock-recruit relationship and unfished spawners per recruit) - vector of length nyears
#' \item R0_annual - annual unfished recruitment (annual ratio of E0 and EPR0) - vector of length nyears
#' \item h_annual - annual steepness (calculated from EPR0 and Arec) - vector of length nyears
#' \item CR - annual compensation ratio, the product of alpha and annual unfished spawners per recruit (EPR0) - vector of length nyears
#' \item R - recruitment - vector of length nyears+1
#' \item R_early - recruitment for the cohorts in first year of the model - vector n_age-1 (where n_age = maxage + 1)
#' \item VB - vulnerable biomass - matrix of nyears x nfleet
#' \item N - abundance at age - matrix of nyears+1 x n_age
#' \item F - apical fishing mortality - matrix of nyears x nfleet
#' \item F_at_age - fishing mortality at age - matrix of nyears x n_age
#' \item F_equilibrium - equilibrium fishing mortality prior to first year - vector of length nfleet
#' \item M - natural mortality - matrix of nyears x n_age
#' \item Z - total mortality - matrix of nyears x n_age
#' \item q - index catchability - vector of length nsurvey
#' \item ivul - index selectivity at age - array of dim nyears+1, n_age, nsurvey
#' \item ivul_len - corresponding index selectivity at length - matrix of nbins x nsurvey
#' \item Ipred - predicted index values - matrix of nyears x nsurvey
#' \item IAApred - predicted index catch at age - array of dim nyears, n_age, nsurvey
#' \item vul - fleet selectivity at age - array of dim nyears+1, n_age, nfleet (or nsel_block)
#' \item vul_len - corresponding fleet selectivity at length - matrix of nbins x nfleet (or nsel_block)
#' \item IALpred - predicted index catch at length - array of dim nyears, nbins, nsurvey
#' \item MLpred - predicted mean length - matrix of nyears x nfleet
#' \item MWpred - predicted mean weight - matrix of nyears x nfleet
#' \item CAApred - predicted catch at age - array of nyears, n_age, nfleet
#' \item CALpred - predicted catch at length - array of nyears, nbins, nfleet
#' \item Cpred - predicted catch in weight - matrix of nyears x nfleet
#' \item CN - predicted catch in numbers - matrix of nyears x nfleet
#' \item dynamic_SSB0 - the dynamic unfished spawning biomass calcaluated by projecting the historical model with zero catches - vector of length nyears+1
#' \item SPR_eq - equilibrium spawning potential ratio calculated from annual F-at-age - vector of length nyears
#' \item SPR_dyn - dynamic (transitional) spawning potential ratio calculated from cumulative survival of cohorts - vector of length nyears
#' \item nll - total objective function of the model - numeric
#' \item nll_fleet - objective function values for each annual data point(s) from fleets - array of nyears x nfleet x 5 (for Catch, equilibrium catch, CAA, CAL, and mean size)
#' \item nll_index - objective function values for each annual data point(s) in the index - array of nyears x nsurvey x 3 (for Index, IAA, and IAL)
#' \item prior - penalty value added to the objective function from priors - numeric
#' \item penalty - additional penalty values added to the objective function due to high F - numeric
#' \item conv - whether the model converged (whether a positive-definite Hessian was obtained) - logical
#' }
#'
#' @slot mean_fit A list of output from fit to mean values of life history parameters in the operating model. The named list consists of:
#'
#' \itemize{
#' \item obj - a list with components returned from \code{\link[TMB]{MakeADFun}}.
#' \item opt - a list with components from calling \code{\link[stats]{nlminb}} to \code{obj}.
#' \item SD - a list (class sdreport) with parameter estimates and their standard errors, obtained from
#' \code{\link[TMB]{sdreport}}.
#' \item report - a list of model output reported from the TMB executable, i.e. \code{obj$report()}. See Misc.
#' }
#' @slot data A \linkS4class{RCMdata} object containing data inputs for the RCM.
#' @slot config A list describing configuration of the RCM:
#' \itemize{
#' \item drop_sim - a vector of simulations that were dropped for the output
#' }
#'
#' @seealso \link{plot.RCModel} \link{RCM}
#' @author Q. Huynh
#' @export RCModel
#' @exportClass RCModel
RCModel <- setClass("RCModel", slots = c(OM = "ANY", SSB = "matrix", NAA = "array",
                                         CAA = "array", CAL = "array", conv = "logical", Misc = "list", mean_fit = "list",
                                         data = "ANY", config = "ANY"))


setMethod("initialize", "RCModel", function(.Object, ...) {
  dots <- list(...)
  if (length(dots)) {
    for(i in names(dots)) slot(.Object, i) <- dots[[i]]
  }
  attr(.Object, "version") <- paste("SAMtool", packageVersion("SAMtool"), "with MSEtool", packageVersion("MSEtool"))
  attr(.Object, "date") <- date()
  attr(.Object, "R.version") <- getRversion()
  
  return(.Object)
})

#' @name plot.RCModel
#' @aliases plot,RCModel,missing-method
#' @title Plot RCM scope output
#' @description Produces HTML file (via markdown) figures of parameter estimates and output from an \linkS4class{Assessment} object.
#' Plots histograms of operating model parameters that are updated by the RCM scoping function, as well as diagnostic plots
#' for the fits to the RCM for each simulation. \code{compare_RCM} plots a short report that compares output from multiple RCM objects,
#' assuming the same model structure, i.e., identical matrix and array dimensions among models, but different data weightings, data omissions, etc.
#'
#' @param x An object of class \linkS4class{RCModel} (output from \link{RCM}).
#' @param compare Logical, if TRUE, the function will run \code{runMSE} to compare the historical period of the operating model
#' and the RCM output.
#' @param filename Character string for the name of the markdown and HTML files.
#' @param dir The directory in which the markdown and HTML files will be saved.
#' @param sims A logical vector of length \code{x@@OM@@nsim} or a numeric vector indicating which simulations to keep.
#' @param Year Optional, a vector of years for the historical period for plotting.
#' @param f_name Character vector for fleet names.
#' @param s_name Character vector for survey names.
#' @param MSY_ref A numeric vector for reference horizontal lines for B/BMSY plots.
#' @param bubble_adj A number to adjust the size of bubble plots (for residuals of age and length comps).
#' @param scenario Optional, a named list to label each simulation in the RCM for plotting, e.g.:
#' \code{list(names = c("low M", "high M"), col = c("blue", "red"))}.
#' @param title Optional character string for an alternative title for the markdown report.
#' @param open_file Logical, whether the HTML document is opened after it is rendered.
#' @param quiet Logical, whether to silence the markdown rendering function.
#' @param render_args A list of other arguments to pass to \link[rmarkdown]{render}.
#' @param ... For \code{compare_RCM}, multiple RCM objects for comparison.
#' @return Returns invisibly the output from \link[rmarkdown]{render}.
#' @importFrom rmarkdown render
#' @seealso \linkS4class{RCModel} \link{RCM}
#' @exportMethod plot
setMethod("plot", signature(x = "RCModel", y = "missing"),
          function(x, compare = FALSE, filename = "RCM", dir = tempdir(), sims = 1:x@OM@nsim, Year = NULL,
                   f_name = NULL, s_name = NULL, MSY_ref = c(0.5, 1), bubble_adj = 10, scenario = list(), title = NULL,
                   open_file = TRUE, quiet = TRUE, render_args, ...) {

            # Run retrospective (dots$retrospective = TRUE with dots$nyr)
            # Or directly provide a retrospective object (dots$retro)
            dots <- list(...)
            if (!is.null(dots$retrospective) && dots$retrospective) {
              message("Running retrospective on mean fit object...")
              if (is.null(dots$nyr)) {
                retro <- retrospective(x)
              } else retro <- retrospective(x, dots$nyr)
            } else if (!is.null(dots$retro)) {
              retro <- dots$retro
            }

            # Update scenario
            if (is.null(scenario$col)) {
              scenario$col <- "red"
              scenario$col2 <- "black"
            } else {
              scenario$col2 <- scenario$col
            }

            if (is.null(scenario$lwd)) scenario$lwd <- 1
            if (is.null(scenario$lty)) scenario$lty <- 1

            ####### Function arguments for rmarkdown::render
            filename_rmd <- paste0(filename, ".Rmd")

            if (missing(render_args)) render_args <- list()
            render_args$input <- file.path(dir, filename_rmd)
            if (is.null(render_args$output_format)) {
              render_args$output_format <- "html_document"
            }
            if (is.null(render_args$output_options)) {
              if (render_args$output_format == "html_document") {
                render_args$output_options <- list(df_print = "paged")
              } else {
                render_args$output_options <- list(toc = TRUE, df_print = "kable")
              }
            }
            if (is.null(render_args$quiet)) render_args$quiet <- quiet

            ####### Assign variables
            OM <- MSEtool::SubCpars(x@OM, sims)
            if (length(x@Misc) == 1) {
              report_list <- x@Misc
            } else {
              report_list <- x@Misc[sims]
            }

            nsim <- OM@nsim
            RCMdata <- x@data
            stopifnot(inherits(RCMdata, "RCMdata"))

            max_age <- OM@maxage
            age <- 0:max_age
            nyears <- OM@nyears
            proyears <- OM@proyears
            if (is.null(Year)) Year <- (OM@CurrentYr - nyears + 1):OM@CurrentYr
            Yearplusone <- c(Year, max(Year) + 1)

            nfleet <- RCMdata@Misc$nfleet
            nsel_block <- RCMdata@Misc$nsel_block
            nsurvey <- RCMdata@Misc$nsurvey
            length_bin <- RCMdata@Misc$lbinmid

            if (is.null(f_name)) f_name <- paste("Fleet", 1:nfleet)
            if (is.null(s_name)) {
              if (nsurvey > 0) {
                s_name <- paste("Index", 1:nsurvey)
              } else {
                s_name <- "Index"
              }
            }

            ####### Document header
            if (is.null(title)) title <- "Operating model (OM) conditioning report for `r ifelse(nchar(OM@Name) > 0, OM@Name, substitute(x))`"
            header <- c("---",
                        paste0("title: \"", title, "\""),
                        "subtitle: Output from SAMtool Rapid Conditioning Model (RCM)",
                        "date: \"`r Sys.Date()`\"",
                        "---",
                        "<style type=\"text/css\">",
                        "h1 { /* Header 1 */",
                        "  font-size: 24px;",
                        "}",
                        "</style>",
                        "",
                        "```{r setup, include = FALSE, echo = FALSE}",
                        "  knitr::opts_chunk$set(collapse = TRUE, echo = FALSE, message = FALSE,",
                        paste0("  fig.width = 6, fig.height = 4.5, ", ifelse(render_args$output_format == "html_document", "", "dpi = 400, "),
                               "out.width = \"650px\", comment = \"#>\")"),
                        "```\n")

            ####### Updated OM
            OM_update <- c("# Summary {.tabset}\n",
                           "## Operating model {.tabset}\n",
                           "### Updated parameters\n", 
                           rmd_RCM_R0(), rmd_RCM_D(), rmd_RCM_Perr(), rmd_RCM_Find(), rmd_RCM_sel())
            
            vary_R0 <- !is.null(OM@cpars$R0) && length(unique(OM@cpars$R0)) > 1 
            vary_D <- !is.null(OM@cpars$D) && length(unique(OM@cpars$R0)) > 1
            vary_hs <- !is.null(OM@cpars$hs) && length(unique(OM@cpars$hs)) > 1
            vary_Mest <- RCMdata@Misc$prior$use_prior[3] && length(report_list) > 1 # Only plot there is a prior
            
            if (sum(vary_R0, vary_D, vary_hs, vary_Mest) > 1) {
              vars <- c("R0", "D", "hs", "Mest")
              var_labs <- c(R0 = "expression(R[0])", D = "\"Depletion\"", hs = "\"Steepness\"", Mest = "\"Natural mortality\"")
              var_names <- c(R0 = "unfished recruitment", D = "depletion", hs = "steepness", Mest = "natural mortality")
              
              if (vary_Mest) OM@cpars[["Mest"]] <- vapply(report_list, getElement, numeric(1), "Mest") # For plotting only
              
              corr_series <- lapply(1:3, function(i) data.frame(x = vars[i], y = vars[(i+1):4])) %>%
                bind_rows()
              corr_rmd <- lapply(1:nrow(corr_series), function(i) {
                x <- corr_series$x[i]
                y <- corr_series$y[i]
                if (eval(as.symbol(paste0("vary_", x))) && eval(as.symbol(paste0("vary_", y)))) { 
                  rmd_corr(paste0("OM@cpars[[\"", x, "\"]]"), paste0("OM@cpars[[\"", y, "\"]]"), 
                           var_labs[x], var_labs[y], 
                           paste0("Correlation between ", var_names[x], " and ", var_names[y], " in operating model.")
                  )
                } else {
                  ""
                }
              })
              
              OM_update <- c(OM_update, "### Parameter correlations\n", do.call(c, corr_rmd))
            }
            
            if (compare) {
              message("Getting Hist object from runMSE...")
              Hist <- runMSE(OM, Hist = TRUE, silent = TRUE, parallel = snowfall::sfIsRunning())
              compare_rmd <- rmd_RCM_Hist_compare()
            } else {
              compare_rmd <- c("### RCM-OM comparison\n",
                               "Re-run `plot()` function with argument `compare = TRUE`.\n\n")
            }
            OM_update <- c(OM_update, compare_rmd)

            ####### Output from all simulations {.tabset}
            fleet_output <- lapply(1:nfleet, rmd_RCM_fleet_output, f_name = f_name)

            if (any(RCMdata@Index > 0, na.rm = TRUE)) {
              index_output <- lapply(1:nsurvey, rmd_RCM_index_output, s_name = s_name)
            } else index_output <- NULL

            all_sims_output <- c(fleet_output, index_output, "### Model predictions\n",
                                 rmd_RCM_initD(), rmd_RCM_R_output(), rmd_RCM_SSB_output(), rmd_log_rec_dev(), 
                                 rmd_RCM_SPR())

            ####### Fit to mean inputs from operating model
            # Generate summary table (parameter estimates)

            if (length(x@mean_fit)) {
              report <- x@mean_fit$report
              conv <- report$conv
              data_mean_fit <- x@mean_fit$obj$env$data
              SD <- x@mean_fit$SD

              if (render_args$output_format == "html_document") {
                sumry <- c("## Fit to mean parameters of the OM {.tabset}\n",
                           "### RCM Estimates\n",
                           "`r sdreport_int(SD) %>% signif(3) %>% as.data.frame()`\n\n")
              } else {
                sumry <- c("## Fit to mean parameters of the OM {.tabset}\n",
                           "### RCM Estimates\n",
                           "`r sdreport_int(SD) %>% signif(3) %>% as.data.frame() %>% knitr::kable(format = \"markdown\")`\n\n")
              }

              # Life History section
              LH_varies_fn <- function(x) {
                n_unique <- apply(x, 2, function(y) length(unique(y)))
                any(n_unique > 1)
              }
              if (any(data_mean_fit$MLpred > 0) && any(data_mean_fit$CALpred > 0)) {
                SD_LAA <- "data_mean_fit$SD_LAA[nyears, ]"
              } else {
                SD_LAA <- ""
              }
              LAA <- rmd_LAA(age = "age", LAA = "data_mean_fit$len_age[nyears, ]", header = "### Life History\n", 
                             SD_LAA = SD_LAA, fig.cap = "Length-at-age in last historical year.")
              if (LH_varies_fn(data_mean_fit$len_age)) {
                LAA_persp <- rmd_persp_plot(x = "Year", y = "age", z = "data_mean_fit$len_age[1:nyears, ]", xlab = "Year", ylab = "Age",
                                            zlab = "Length-at-age", phi = 35, theta = 45, expand = 0.55, fig.cap = "Annual length-at-age.")
              } else LAA_persp <- NULL

              mat <- rmd_mat(age = "age", "data_mean_fit$mat[nyears, ]", fig.cap = "Maturity-at-age in last historical year.")
              if (LH_varies_fn(data_mean_fit$mat)) {
                mat_persp <- rmd_persp_plot(x = "Year", y = "age", z = "data_mean_fit$mat[1:nyears, ]", xlab = "Year", ylab = "Age",
                                            zlab = "Maturity-at-age", phi = 35, theta = 45, expand = 0.55, fig.cap = "Annual maturity-at-age.")
              } else mat_persp <- NULL

              if (data_mean_fit$use_prior[3]) {
                NatM <- rmd_at_age(age = "age", "rep(report$Mest, length(age))", fig.cap = "Natural mortality (time constant).", label = "Natural mortality")
              } else {
                NatM <- rmd_at_age(age = "age", "data_mean_fit$M[nyears, ]", fig.cap = "Natural mortality in last historical year.", label = "Natural mortality")
              }
              if (!data_mean_fit$use_prior[3] && LH_varies_fn(data_mean_fit$M)) {
                NatM_persp <- rmd_persp_plot(x = "Year", y = "age", z = "data_mean_fit$M[1:nyears, ]", xlab = "Year", ylab = "Age",
                                             zlab = "Natural mortality", phi = 35, theta = 45, expand = 0.55, fig.cap = "Annual M-at-age.")
              } else NatM_persp <- NULL

              LH_section <- c(LAA, LAA_persp, mat, mat_persp, NatM, NatM_persp)

              # Data and fit section
              individual_matrix_fn <- function(i, obs, pred, fig.cap, label, resids = FALSE, match = FALSE) {
                if (resids) {
                  rmd_assess_resid2("Year", paste0(obs, "[, ", i, "]"), paste0(pred, "[, ", i, "]"),
                                    fig.cap = paste(fig.cap, i), label = label[i])
                } else {
                  rmd_assess_fit2("Year", paste0(obs, "[, ", i, "]"), paste0(pred, "[, ", i, "]"),
                                  fig.cap = paste(fig.cap, i), label = label[i], match = match)
                }
              }
              individual_array_fn <- function(i, obs, pred, comps = c("age", "length"), label, plot_mean = TRUE) {
                comps <- match.arg(comps)
                obs2 <- paste0(obs, "[, , ", i, "]")
                pred2 <- paste0(pred, "[, , ", i, "]")
                fig.cap <- paste0("Observed (black) and predicted (red) ", comps, " composition from ", label[i], ".")
                fig.cap2 <- paste0("Residuals for ", comps, " composition from ", label[i], ".")
                fig.cap3 <- paste0("Observed (black) and predicted (red) mean ", comps, " from the composition data for ", 
                                   label[i], ".")
                if (comps == "age") {
                  #rr <- rmd_fit_comps("Year", obs2, pred2, type = c("annual", "bubble_residuals", "mean"), 
                  #                    ages = "age", fig.cap = fig.cap)
                  rr <- rmd_fit_comps("Year", obs2, pred2, type = "annual", ages = "age", fig.cap = fig.cap)
                  rr2 <- rmd_fit_comps("Year", obs2, pred2, type = "bubble_residuals", ages = "age", fig.cap = fig.cap2)
                  rr3 <- rmd_fit_comps("Year", obs2, pred2, type = "mean", ages = "age", fig.cap = fig.cap3)
                } else {
                  rr <- rmd_fit_comps("Year", obs2, pred2, type = "annual", CAL_bins = "RCMdata@Misc$lbinmid", fig.cap = fig.cap)
                  rr2 <- rmd_fit_comps("Year", obs2, pred2, type = "bubble_residuals", CAL_bins = "RCMdata@Misc$lbinmid", fig.cap = fig.cap2)
                  rr3 <- rmd_fit_comps("Year", obs2, pred2, type = "mean", CAL_bins = "RCMdata@Misc$lbinmid", fig.cap = fig.cap3)
                }
                c(rr, rr2, rr3)
              }

              if (any(RCMdata@Chist > 0, na.rm = TRUE)) {
                C_matplot <- rmd_matplot(x = "Year", y = "RCMdata@Chist", col = "rich.colors(nfleet)",
                                         xlab = "Year", ylab = "Catch", legend.lab = "f_name",
                                         fig.cap = "Catch time series.", header = "### Data and Fit {.tabset}\n#### Catch \n")

                if (data_mean_fit$condition == "effort" || ncol(RCMdata@Chist) > 1) {
                  C_plots <- lapply(1:nfleet, individual_matrix_fn, obs = "RCMdata@Chist", pred = "report$Cpred",
                                    fig.cap = "catch from fleet", label = f_name, match = data_mean_fit$condition == "catch2")
                } else C_plots <- NULL
              } else C_matplot <- C_plots <- NULL

              if (any(RCMdata@Ehist > 0, na.rm = TRUE)) {
                if (is.null(C_matplot)) {
                  E_header <- "### Data and Fit {.tabset}\n#### Effort \n"
                } else {
                  E_header <- "#### Effort \n"
                }
                E_matplot <- rmd_matplot(x = "Year", y = "data_mean_fit$E_hist", col = "rich.colors(nfleet)",
                                         xlab = "Year", ylab = "Effort", legend.lab = "f_name", fig.cap = "Effort time series.", header = E_header)
              } else E_matplot <- NULL

              if (any(RCMdata@Index > 0, na.rm = TRUE)) {
                I_plots <- c("#### Index \n",
                             lapply(1:nsurvey, individual_matrix_fn, obs = "RCMdata@Index", pred = "report$Ipred",
                                    fig.cap = "index from survey", label = s_name),
                             lapply(1:nsurvey, individual_matrix_fn, obs = "RCMdata@Index", pred = "report$Ipred",
                                    fig.cap = "index from survey", label = paste(s_name, "Residuals"), resids = TRUE))
              } else I_plots <- NULL

              if (any(RCMdata@CAA > 0, na.rm = TRUE)) {
                CAA_plots <- c("#### Age comps \n",
                               lapply(1:nfleet, individual_array_fn, obs = "RCMdata@CAA", pred = "report$CAApred", comps = "age", label = f_name))
              } else CAA_plots <- NULL

              if (any(RCMdata@CAL > 0, na.rm = TRUE)) {
                CAL_plots <- c("#### Length comps \n",
                               lapply(1:nfleet, individual_array_fn, obs = "RCMdata@CAL", pred = "report$CALpred", comps = "length", label = f_name))
              } else CAL_plots <- NULL

              if (any(RCMdata@MS > 0, na.rm = TRUE)) {
                if (RCMdata@MS_type == "length") {
                  MS_label <- paste("Mean Length from", f_name)
                } else {
                  MS_label <- paste("Mean Weight from", f_name)
                }
                MS_plots <- c(paste0("#### Mean ", RCMdata@MS_type, "\n"),
                              lapply(1:nfleet, individual_matrix_fn, obs = "RCMdata@MS",
                                     pred = ifelse(RCMdata@MS_type == "length", "report$MLpred", "report$MWpred"),
                                     fig.cap = paste("mean", RCMdata@MS_type, "from fleet"), label = MS_label))
              } else MS_plots <- NULL

              if (any(RCMdata@IAA > 0, na.rm = TRUE)) {
                IAA_plots <- c("#### Index age comps \n",
                                 lapply(1:nsurvey, individual_array_fn, obs = "RCMdata@IAA", pred = "report$IAApred", comps = "age", label = s_name))
              } else IAA_plots <- NULL

              if (any(RCMdata@IAL > 0, na.rm = TRUE)) {
                IAL_plots <- c("#### Index length comps \n",
                                 lapply(1:nsurvey, individual_array_fn, obs = "RCMdata@IAL", pred = "report$IALpred", comps = "length", label = s_name))
              } else IAL_plots <- NULL

              data_section <- c(C_matplot, C_plots, E_matplot, I_plots, CAA_plots, CAL_plots, MS_plots, IAA_plots, IAL_plots)

              # Model output
              sel_matplot <- rmd_matplot(x = "age", y = "matrix(report$vul[nyears, , ], max_age + 1, nfleet)", col = "rich.colors(nfleet)",
                                         xlab = "Age", ylab = "Selectivity", legend.lab = "f_name",
                                         fig.cap = "Terminal year selectivity by fleet.", header = "### Output \n")

              F_matplot <- rmd_matplot(x = "Year", y = "report$F", col = "rich.colors(nfleet)",
                                       xlab = "Year", ylab = "Fishing Mortality (F)", legend.lab = "f_name",
                                       fig.cap = "Time series of fishing mortality by fleet.")

              if (length(unique(report$E0)) > 1) {
                SSB0_eq_plot <- rmd_assess_timeseries("structure(report$E0, names = Year)", 
                                                      "equilibrium unfished spawning biomass (growth and/or M are time-varying)",
                                                      "expression(Equilibrium~SSB[0])")
                SSB_SSB0_plot <- rmd_SSB_SSB0(FALSE, "structure(report$E/report$E0_SR, names = Yearplusone)",
                                              fig.cap = "spawning depletion, using the unfished biomass from the replacement line at the beginning of time series")
              } else {
                SSB0_eq_plot <- NULL
                SSB_SSB0_plot <- rmd_SSB_SSB0(FALSE, "structure(report$E/report$E0_SR, names = Yearplusone)")
              }

              N_bubble <- rmd_bubble("Yearplusone", "report$N", ages = "age", fig.cap = "Predicted abundance-at-age.")
              CAA_bubble <- rmd_bubble("Year", "apply(report$CAApred, 1:2, sum)", ages = "age",
                                       fig.cap = "Predicted catch-at-age (summed over all fleets).", 
                                       bubble_adj = as.character(bubble_adj))

              if (sum(report$CALpred, na.rm = TRUE) > 0) {
                CAL_bubble <- rmd_bubble("Year", "apply(report$CALpred, 1:2, sum)", CAL_bins = "RCMdata@Misc$lbinmid",
                                         fig.cap = "Predicted catch-at-length (summed over all fleets).", bubble_adj = as.character(bubble_adj))
              } else CAL_bubble <- NULL

              ts_output <- c(sel_matplot, F_matplot, rmd_SSB("structure(report$E, names = Yearplusone)"), 
                             SSB0_eq_plot, SSB_SSB0_plot, 
                             rmd_dynamic_SSB0("structure(report$dynamic_SSB0, names = Yearplusone)"), 
                             rmd_R("structure(report$R, names = Yearplusone)"), 
                             rmd_RCM_SR(), rmd_RCM_SPR2(),
                             rmd_residual("structure(report$log_rec_dev, names = Year)", 
                                          fig.cap = "Time series of recruitment deviations.", 
                                          label = "log-Recruitment deviations"),
                             rmd_residual("structure(report$log_rec_dev, names = Year)", 
                                          "ifelse(data_mean_fit$est_rec_dev == 1, as.list(SD, \"Std. Error\")$log_rec_dev, 0)", 
                                          fig.cap = "Time series of recruitment deviations with 95% confidence intervals.",
                                          label = "log-Recruitment deviations", 
                                          conv_check = TRUE),
                             rmd_N("structure(rowSums(report$N), names = Yearplusone)"), 
                             N_bubble, 
                             CAA_bubble, 
                             CAL_bubble)

              nll <- RCM_get_likelihoods(report, RCMdata@Misc$LWT, f_name, s_name)
              
              nll_table <- c("### Likelihood components\n",
                             "#### Summary\n",
                             "`r nll[[1]] ", ifelse(render_args$output_format == "html_document", "", "%>% knitr::kable(format = \"markdown\")"), "`\n\n",
                             "#### Fleet likelihoods\n",
                             "`r nll[[2]] ", ifelse(render_args$output_format == "html_document", "", "%>% knitr::kable(format = \"markdown\")"), "`\n\n",
                             "#### Fleet likelihood weights\n",
                             "`r nll[[3]] ", ifelse(render_args$output_format == "html_document", "", "%>% knitr::kable(format = \"markdown\")"), "`\n\n",
                             "#### Survey likelihoods\n",
                             "`r nll[[4]] ", ifelse(render_args$output_format == "html_document", "", "%>% knitr::kable(format = \"markdown\")"), "`\n\n",
                             "#### Survey likelihood weights\n",
                             "`r nll[[5]] ", ifelse(render_args$output_format == "html_document", "", "%>% knitr::kable(format = \"markdown\")"), "`\n\n")
              
              corr_matrix <- c("### Correlation matrix\n",
                               "`r SD$env$corr.fixed %>% structure(dimnames = list(make_unique_names(rownames(.)), make_unique_names(colnames(.)))) %>% as.data.frame()`\n\n")
              
              if (exists("retro", inherits = FALSE)) {
                ret <- rmd_RCM_retrospective(render_args)
              } else ret <- NULL
              

              mean_fit_rmd <- c(sumry, LH_section, data_section, ts_output, nll_table, corr_matrix, ret)
            } else {
              mean_fit_rmd <- c("## Fit to mean parameters of OM {.tabset}\n",
                                "No model found. Re-run `RCM()` with `mean_fit = TRUE`.\n\n")
            }

            rmd <- c(header, OM_update, all_sims_output, mean_fit_rmd, rmd_footer())
            if (is.list(rmd)) rmd <- do.call(c, rmd)

            # Generate markdown report
            if (!dir.exists(dir)) {
              message("Creating directory: \n", dir)
              dir.create(dir)
            }
            write(rmd, file = file.path(dir, filename_rmd))
            message("Generated markdown file: ", file.path(dir, filename_rmd))

            # Rendering markdown file
            message("Rendering markdown file...")
            output_filename <- do.call(rmarkdown::render, render_args)
            message("Rendered file: ", output_filename)

            if (open_file) browseURL(output_filename)
            invisible(output_filename)
          })





#' @rdname plot.RCModel
#' @export
compare_RCM <- function(..., compare = FALSE, filename = "compare_RCM", dir = tempdir(), Year = NULL,
                        f_name = NULL, s_name = NULL, MSY_ref = c(0.5, 1), bubble_adj = 10, scenario = list(), title = NULL,
                        open_file = TRUE, quiet = TRUE, render_args) {

  dots <- list(...)
  test <- vapply(dots, inherits, logical(1), what = "RCModel") %>% all()
  if (!test) stop("Objects provided are not of class RCModel.", call. = FALSE)
  
  # Update scenario
  if (is.null(scenario$names)) scenario$names <- paste("Scenario", 1:length(dots))
  if (is.null(scenario$col)) {
    scenario$col <- gplots::rich.colors(length(dots))
  }
  scenario$col2 <- scenario$col

  if (is.null(scenario$lwd)) scenario$lwd <- 1
  if (is.null(scenario$lty)) scenario$lty <- 1

  ####### Function arguments for rmarkdown::render
  filename_rmd <- paste0(filename, ".Rmd")

  if (missing(render_args)) render_args <- list()
  render_args$input <- file.path(dir, filename_rmd)
  if (is.null(render_args$output_format)) {
    render_args$output_format <- "html_document"
  }
  if (is.null(render_args$output_options)) {
    if (render_args$output_format == "html_document") {
      render_args$output_options <- list(df_print = "paged")
    } else {
      render_args$output_options <- list(toc = TRUE, df_print = "kable")
    }
  }
  if (is.null(render_args$quiet)) render_args$quiet <- quiet

  ####### Assign variables
  x <- dots[[1]] # Dummy variable
  report_list <- lapply(dots, function(xx) if (length(xx@mean_fit)) return(xx@mean_fit$report) else stop("Error in RCM objects."))

  nsim <- length(report_list)
  RCMdata <- dots[[1]]@data

  max_age <- dots[[1]]@OM@maxage
  age <- 0:max_age
  nyears <- dots[[1]]@OM@nyears
  if (is.null(Year)) Year <- (dots[[1]]@OM@CurrentYr - nyears + 1):dots[[1]]@OM@CurrentYr
  Yearplusone <- c(Year, max(Year) + 1)

  nfleet <- vapply(dots, function(xx) xx@data@Misc$nfleet, numeric(1)) %>% unique()
  nsurvey <- vapply(dots, function(xx) xx@data@Misc$nsurvey, numeric(1)) %>% unique()
  length_bin <- RCMdata@Misc$lbinmid

  if (is.null(f_name)) f_name <- paste("Fleet", 1:nfleet)
  if (is.null(s_name)) s_name <- paste("Index", 1:nsurvey)

  ####### Document header
  if (is.null(title)) title <- "Comparisons of OM conditioning"
  header <- c("---",
              paste0("title: \"", title, "\""),
              "subtitle: Output from Rapid Conditioning Model (RCM)",
              "date: \"`r Sys.Date()`\"",
              "---",
              "<style type=\"text/css\">",
              "h1 { /* Header 1 */",
              "  font-size: 24px;",
              "}",
              "</style>",
              "",
              "```{r setup, include = FALSE, echo = FALSE}",
              "  knitr::opts_chunk$set(collapse = TRUE, echo = FALSE, message = FALSE,",
              paste0("  fig.width = 6, fig.height = 4.5, ", ifelse(render_args$output_format == "html_document", "", "dpi = 400, "),
                     "out.width = \"650px\", comment = \"#>\")"),
              "```\n")


  ####### Output from all simulations {.tabset}
  fleet_output <- lapply(1:nfleet, rmd_RCM_fleet_output, f_name = f_name)

  if (any(RCMdata@Index > 0, na.rm = TRUE)) {
    survey_output <- lapply(1:nsurvey, rmd_RCM_index_output, s_name = s_name)
  } else survey_output <- NULL

  #### MSY comparisons
  if (compare) {
    message("Running runMSE() to get MSY reference points...")
    if (snowfall::sfIsRunning()) {
      Hist <- snowfall::sfClusterApplyLB(dots, function(xx) runMSE(xx@OM, Hist = TRUE, silent = TRUE))
    } else {
      Hist <- lapply(dots, function(xx) runMSE(xx@OM, Hist = TRUE, silent = TRUE))
    }

    SSB_MSY <- c("```{r, fig.cap = \"Spawning biomass (SSB) relative to MSY from the operating model.\"}",
                 "SSB <- do.call(rbind, lapply(report_list, function(x) x$E[1:nyears]))",
                 "SSBMSY <- vapply(Hist, function(x) mean(x@Ref$ReferencePoints$SSBMSY), numeric(1))",
                 "SSB_SSBMSY <- t(SSB/SSBMSY)",
                 "matplot(Year, SSB_SSBMSY, typ = \"n\", xlab = \"Year\", ylab = expression(SSB/SSB[MSY]), ylim = c(0, 1.1 * max(SSB_SSBMSY)))",
                 "matlines(Year, SSB_SSBMSY, col = scenario$col2)",
                 "abline(h = c(0, MSY_ref), col = \"grey\")",
                 "if (!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col2, lty = scenario$lty)",
                 "```\n")

    ref_pt_fn <- function(xx) c(mean(xx@Ref$ReferencePoints$FMSY), mean(xx@Ref$ReferencePoints$MSY), mean(xx@Ref$ReferencePoints$SSBMSY_SSB0))
    ref_pt <- do.call(cbind, lapply(Hist, ref_pt_fn)) %>%
      structure(dimnames = list(c("FMSY", "MSY", "Spawning depletion at MSY"), scenario$names)) %>% as.data.frame()

    if (render_args$output_format == "html_document") {
      rmd_ref_pt <- paste0("## Reference points \n",
                           "`r ref_pt`\n\n")
    } else {
      rmd_ref_pt <- paste0("## Reference points \n",
                           "`r ref_pt %>% knitr::kable(format = \"markdown\")`\n\n")
    }

  } else {
    SSB_MSY <- rmd_ref_pt <- ""
  }

  all_sims_output <- c("# Summary {.tabset}\n\n", fleet_output, survey_output, "### Model predictions\n",
                       rmd_RCM_initD(), rmd_RCM_R_output(), rmd_RCM_SSB_output(), SSB_MSY, rmd_log_rec_dev())

  #### Likelihoods
  nll <- Map(RCM_get_likelihoods, x = report_list, LWT = lapply(dots, function(xx) xx@data@Misc$LWT),
             MoreArgs = list(f_name = f_name, s_name = s_name))

  summary_nll <- vapply(nll, function(xx) xx[[1]] %>% as.matrix(), numeric(6)) %>%
    structure(dimnames = list(rownames(nll[[1]][[1]]), scenario$names)) %>% as.data.frame()

  if (render_args$output_format == "html_document") {

    nll_table_fn <- function(i, ii) {
      c(paste0("#### ", scenario$names[i]),
        paste0("`r nll[[", i, "]][[", ii, "]]`\n\n"))
    }

    nll_table <- c("## Likelihood components\n",
                   "### Summary\n",
                   "`r summary_nll`\n\n",
                   "### Fleet likelihoods {.tabset}\n",
                   do.call(c, lapply(1:length(dots), nll_table_fn, ii = 2)),
                   "### Fleet weights {.tabset}\n",
                   do.call(c, lapply(1:length(dots), nll_table_fn, ii = 3)),
                   "### Survey likelihoods {.tabset}\n",
                   do.call(c, lapply(1:length(dots), nll_table_fn, ii = 4)),
                   "### Survey weights {.tabset}\n",
                   do.call(c, lapply(1:length(dots), nll_table_fn, ii = 5)))
  } else {

    nll_table_fn <- function(i, ii) {
      c(paste0("#### ", scenario$names[i]),
        paste0("`r nll[[", i, "]][[", ii, "]] %>% knitr::kable(format = \"markdown\")`\n\n"))
    }
    nll_table <- c("## Likelihood components\n",
                   "### Summary\n",
                   "`r summary_nll %>% knitr::kable(format = \"markdown\")`\n\n",
                   "### Fleet likelihoods\n",
                   do.call(c, lapply(1:length(dots), nll_table_fn, ii = 2)),
                   "### Fleet weights\n",
                   do.call(c, lapply(1:length(dots), nll_table_fn, ii = 3)),
                   "### Survey likelihoods\n",
                   do.call(c, lapply(1:length(dots), nll_table_fn, ii = 4)),
                   "### Survey weights\n",
                   do.call(c, lapply(1:length(dots), nll_table_fn, ii = 5)))
  }

  rmd <- c(header, all_sims_output, nll_table, rmd_ref_pt, rmd_footer())
  if (is.list(rmd)) rmd <- do.call(c, rmd)

  # Generate markdown report
  if (!dir.exists(dir)) {
    message("Creating directory: \n", dir)
    dir.create(dir)
  }
  write(rmd, file = file.path(dir, filename_rmd))
  message("Generated markdown file: ", file.path(dir, filename_rmd))

  # Rendering markdown file
  message("Rendering markdown file...")
  output_filename <- do.call(rmarkdown::render, render_args)
  message("Rendered file: ", output_filename)

  if (open_file) browseURL(output_filename)
  invisible(output_filename)
}



