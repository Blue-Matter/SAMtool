
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
#' @slot Misc A list of length \code{OM@@nsim} with more output from the fitted RCM. Within each simulation, items of interest include:
#'
#' \itemize{
#' \item B - total biomass - vector of length nyears+1
#' \item E0 - annual unfished spawning biomass - vector of length nyears
#' \item E0_SR - unfished spawning biomass for the stock-recruit relationship - numeric
#' \item CR - annual compensation ratio - vector of length nyears
#' \item Arec - alpha parameter of the stock-recruit relationship - numeric
#' \item Brec - beta parameter of the stock-recruit relationship - numeric
#' \item R - recruitment - vector of length nyears+1
#' \item R_early - recruitment for the cohorts in first year of the model - vector maxage-1
#' \item VB - vulnerable biomass - matrix of nyears x nfleet
#' \item N - abundance at age - matrix of nyears+1 x maxage
#' \item F - apical fishing mortality - matrix of nyears x nfleet
#' \item F_at_age - fishing mortality at age - array of nyears x maxage x nfleet
#' \item F_equilibrium - equilibrium fishing mortality prior to first year - vector of length nfleet
#' \item M - natural mortality - matrix of nyears x maxage
#' \item Z - total mortality - matrix of nyears x maxage
#' \item q - survey catchability - vector of length nsurvey
#' \item s_vul - survey selectivity at age - array of dim nyears+1, maxage, nsurvey
#' \item s_vul_len - corresponding survey selectivity at length - matrix of nbins x nsurvey
#' \item Ipred - predicted index values - matrix of nyears x nsurvey
#' \item s_CAApred - predicted survey catch at age - array of dim nyears, maxage, nsurvey
#' \item vul - fleet selectivity at age - array of dim nyears+1, maxage, nfleet (or nsel_block)
#' \item vul_len - corresponding fleet selectivity at length - matrix of nbins x nfleet (or nsel_block)
#' \item s_CALpred - predicted survey catch at length - array of dim nyears, nbins, nsurvey
#' \item MLpred - predicted mean length - matrix of nyears x nfleet
#' \item MWpred - predicted mean weight - matrix of nyears x nfleet
#' \item CAApred - predicted catch at age - array of nyears, maxage, nfleet
#' \item CALpred - predicted catch at length - array of nyears, nbins, nfleet
#' \item Cpred - predicted catch in weight - matrix of nyears x nfleet
#' \item CN - predicted catch in numbers - matrix of nyears x nfleet
#' \item nll - Total objective function of the model - numeric
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
#' @slot data A list of the data inputs for the RCM.
#' @slot config A data frame describing configuration of the RCM (not currently used).
#'
#' @seealso \link{plot.RCModel} \link{RCM}
#' @author Q. Huynh
#' @export RCModel
#' @exportClass RCModel
RCModel <- setClass("RCModel", slots = c(OM = "ANY", SSB = "matrix", NAA = "array",
                                         CAA = "array", CAL = "array", conv = "logical", Misc = "list", mean_fit = "list",
                                         data = "list", config = "data.frame"))


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
          function(x, compare = TRUE, filename = "RCM", dir = tempdir(), sims = 1:x@OM@nsim, Year = NULL,
                   f_name = NULL, s_name = NULL, MSY_ref = c(0.5, 1), bubble_adj = 10, scenario = list(), title = NULL,
                   open_file = TRUE, quiet = TRUE, render_args, ...) {

            # Run retrospective (dots$retrospective = TRUE with dots$nyr)
            # Or directly provide a retrospective object (dots$retro)
            dots <- list(...)
            if(!is.null(dots$retrospective) && dots$retrospective) {
              message("Running retrospective on mean fit object...")
              if(is.null(dots$nyr)) {
                retro <- retrospective(x)
              } else retro <- retrospective(x, dots$nyr)
            } else if(!is.null(dots$retro)) {
              retro <- dots$retro
            }

            # Update scenario
            if(is.null(scenario$col)) {
              scenario$col <- "red"
              scenario$col2 <- "black"
            } else {
              scenario$col2 <- scenario$col
            }

            if(is.null(scenario$lwd)) scenario$lwd <- 1
            if(is.null(scenario$lty)) scenario$lty <- 1

            ####### Function arguments for rmarkdown::render
            filename_rmd <- paste0(filename, ".Rmd")

            if(missing(render_args)) render_args <- list()
            render_args$input <- file.path(dir, filename_rmd)
            if(is.null(render_args$output_format)) {
              render_args$output_format <- "html_document"
            }
            if(is.null(render_args$output_options)) {
              if(render_args$output_format == "html_document") {
                render_args$output_options <- list(df_print = "paged")
              } else {
                render_args$output_options <- list(toc = TRUE, df_print = "kable")
              }
            }
            if(is.null(render_args$quiet)) render_args$quiet <- quiet

            ####### Assign variables
            OM <- MSEtool::SubCpars(x@OM, sims)
            mean_fit <- x@mean_fit
            report_list <- x@Misc[sims]

            nsim <- OM@nsim
            data <- x@data

            max_age <- OM@maxage
            age <- 0:max_age
            nyears <- OM@nyears
            proyears <- OM@proyears
            if(is.null(Year)) Year <- (OM@CurrentYr - nyears + 1):OM@CurrentYr
            Yearplusone <- c(Year, max(Year) + 1)

            # Backwards compatibility
            if(is.null(data$nsel_block)) {
              data$nsel_block <- data$nfleet
              data$sel_block <- matrix(1:data$nfleet, nyears, data$nfleet, byrow = TRUE)
            }

            nfleet <- data$nfleet
            nsel_block <- data$nsel_block
            nsurvey <- data$nsurvey
            length_bin <- data$length_bin

            if(is.null(f_name)) f_name <- paste("Fleet", 1:nfleet)
            if(is.null(s_name)) {
              if(nsurvey > 0) {
                s_name <- paste("Survey", 1:nsurvey)
              } else {
                s_name <- "Survey"
              }
            }

            ####### Document header
            if(is.null(title)) title <- "Operating model (OM) conditioning report for `r ifelse(nchar(OM@Name) > 0, OM@Name, substitute(OM))`"
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

            ####### Updated historical OM parameters
            #Year_matrix <- matrix(Year, ncol = nsim, nrow = nyears)
            #Yearplusone_matrix <- matrix(Yearplusone, ncol = nsim, nrow = nyears+1)

            OM_update <- c("# Summary {.tabset}\n",
                           "## Updated historical OM parameters\n", rmd_RCM_R0(),
                           rmd_RCM_D(), rmd_RCM_Perr(), rmd_RCM_Find(), rmd_RCM_sel())

            ####### Output from all simulations {.tabset}
            fleet_output <- lapply(1:nfleet, rmd_RCM_fleet_output, f_name = f_name)

            if(any(data$Index > 0, na.rm = TRUE)) {
              survey_output <- lapply(1:nsurvey, rmd_RCM_survey_output, s_name = s_name)
            } else survey_output <- NULL

            all_sims_output <- c(fleet_output, survey_output, "### Model predictions\n",
                                 rmd_RCM_initD(), rmd_RCM_R_output(), rmd_RCM_SSB_output(), rmd_log_rec_dev())

            ####### Fit to mean inputs from operating model
            # Generate summary table (parameter estimates)

            if(length(x@mean_fit) > 0) {
              SD <- x@mean_fit$SD
              report <- x@mean_fit$report
              data_mean_fit <- x@mean_fit$obj$env$data
              length_bin <- data_mean_fit$length_bin

              # Backwards compatibility
              if(is.null(data_mean_fit$nsel_block)) {
                data_mean_fit$nsel_block <- data$nfleet
                data_mean_fit$sel_block <- matrix(1:data$nfleet, nyears, data$nfleet, byrow = TRUE)
              }

              conv <- report$conv

              SD2 <- sdreport_int(SD) %>% signif(4) %>% format() %>% as.data.frame()
              if(render_args$output_format == "html_document") {
                sumry <- c("## Fit to mean parameters of the OM {.tabset}\n",
                           "### RCM Estimates\n",
                           "`r SD2`\n\n")
              } else {
                sumry <- c("## Fit to mean parameters of the OM {.tabset}\n",
                           "### RCM Estimates\n",
                           "`r SD2 %>% knitr::kable(format = \"markdown\")`\n\n")
              }

              # Life History section
              LH_varies_fn <- function(x) {
                n_unique <- apply(x, 2, function(y) length(unique(y)))
                any(n_unique > 1)
              }
              LAA <- rmd_at_age(age, data_mean_fit$len_age[nyears, ], header = "### Life History\n", fig.cap = "Length-at-age in last historical year.",
                                label = "Mean Length-at-age")
              if(LH_varies_fn(data_mean_fit$len_age)) {
                LAA_persp <- rmd_persp_plot(x = "Year", y = "age", z = "data_mean_fit$len_age[1:nyears, ]", xlab = "Year", ylab = "Age",
                                            zlab = "Length-at-age", phi = 35, theta = 45, expand = 0.55, fig.cap = "Annual length-at-age.")
              } else LAA_persp <- NULL

              mat <- rmd_mat(age, data_mean_fit$mat[nyears, ], fig.cap = "Maturity-at-age in last historical year.")
              if(LH_varies_fn(data_mean_fit$mat)) {
                mat_persp <- rmd_persp_plot(x = "Year", y = "age", z = "data_mean_fit$mat[1:nyears, ]", xlab = "Year", ylab = "Age",
                                            zlab = "Maturity-at-age", phi = 35, theta = 45, expand = 0.55, fig.cap = "Annual maturity-at-age.")
              } else mat_persp <- NULL

              if(data_mean_fit$use_prior[3]) {
                NatM <- rmd_at_age(age, rep(report$Mest, length(age)), fig.cap = "Natural mortality (time constant).", label = "Natural mortality")
              } else {
                NatM <- rmd_at_age(age, data_mean_fit$M[nyears, ], fig.cap = "Natural mortality in last historical year.", label = "Natural mortality")
              }
              if(!data_mean_fit$use_prior[3] && LH_varies_fn(data_mean_fit$M)) {
                NatM_persp <- rmd_persp_plot(x = "Year", y = "age", z = "data_mean_fit$M[1:nyears, ]", xlab = "Year", ylab = "Age",
                                             zlab = "Natural mortality", phi = 35, theta = 45, expand = 0.55, fig.cap = "Annual M-at-age.")
              } else NatM_persp <- NULL

              LH_section <- c(LAA, LAA_persp, mat, mat_persp, NatM, NatM_persp)

              # Data and fit section
              individual_matrix_fn <- function(i, obs, pred, fig.cap, label, resids = FALSE, match = FALSE) {
                if(resids) {
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
                if(comps == "age") {
                  #rr <- rmd_fit_comps("Year", obs2, pred2, type = c("annual", "bubble_residuals", "mean"), 
                  #                    ages = "age", fig.cap = fig.cap)
                  rr <- rmd_fit_comps("Year", obs2, pred2, type = "annual", ages = "age", fig.cap = fig.cap)
                  rr2 <- rmd_fit_comps("Year", obs2, pred2, type = "bubble_residuals", ages = "age", fig.cap = fig.cap2)
                  rr3 <- rmd_fit_comps("Year", obs2, pred2, type = "mean", ages = "age", fig.cap = fig.cap3)
                } else {
                  rr <- rmd_fit_comps("Year", obs2, pred2, type = "annual", CAL_bins = "data$length_bin", fig.cap = fig.cap)
                  rr2 <- rmd_fit_comps("Year", obs2, pred2, type = "bubble_residuals", CAL_bins = "data$length_bin", fig.cap = fig.cap2)
                  rr3 <- rmd_fit_comps("Year", obs2, pred2, type = "mean", CAL_bins = "data$length_bin", fig.cap = fig.cap3)
                }
                c(rr, rr2, rr3)
              }

              if(any(data$Chist > 0, na.rm = TRUE)) {
                C_matplot <- rmd_matplot(x = "matrix(Year, nyears, nfleet)", y = "data$Chist", col = "rich.colors(nfleet)",
                                         xlab = "Year", ylab = "Catch", legend.lab = "f_name",
                                         fig.cap = "Catch time series.", header = "### Data and Fit {.tabset}\n#### Catch \n")

                if(data_mean_fit$condition == "effort" || ncol(data$Chist) > 1) {
                  C_plots <- lapply(1:nfleet, individual_matrix_fn, obs = "data$Chist", pred = "report$Cpred",
                                    fig.cap = "catch from fleet", label = f_name, match = data_mean_fit$condition == "catch2")
                } else C_plots <- NULL
              } else C_matplot <- C_plots <- NULL

              if(any(data$Ehist > 0, na.rm = TRUE)) {
                if(is.null(C_matplot)) {
                  E_header <- "### Data and Fit {.tabset}\n#### Effort \n"
                } else {
                  E_header <- "#### Effort \n"
                }
                E_matplot <- rmd_matplot(x = "matrix(Year, nyears, nfleet)", y = "data_mean_fit$E_hist", col = "rich.colors(nfleet)",
                                         xlab = "Year", ylab = "Effort", legend.lab = "f_name", fig.cap = "Effort time series.", header = E_header)
              } else E_matplot <- NULL

              if(any(data$Index > 0, na.rm = TRUE)) {
                I_plots <- c("#### Surveys \n",
                             lapply(1:nsurvey, individual_matrix_fn, obs = "data$Index", pred = "report$Ipred",
                                    fig.cap = "index from survey", label = s_name),
                             lapply(1:nsurvey, individual_matrix_fn, obs = "data$Index", pred = "report$Ipred",
                                    fig.cap = "index from survey", label = paste(s_name, "Residuals"), resids = TRUE))
              } else I_plots <- NULL

              if(any(data$CAA > 0, na.rm = TRUE)) {
                CAA_plots <- c("#### Age comps \n",
                               lapply(1:nfleet, individual_array_fn, obs = "data$CAA", pred = "report$CAApred", comps = "age", label = f_name))
              } else CAA_plots <- NULL

              if(any(data$CAL > 0, na.rm = TRUE)) {
                CAL_plots <- c("#### Length comps \n",
                               lapply(1:nfleet, individual_array_fn, obs = "data$CAL", pred = "report$CALpred", comps = "length", label = f_name))
              } else CAL_plots <- NULL

              if(any(data$MS > 0, na.rm = TRUE)) {
                if(data$MS_type == "length") {
                  MS_label <- paste("Mean Length from", f_name)
                } else {
                  MS_label <- paste("Mean Weight from", f_name)
                }
                MS_plots <- c(paste0("#### Mean ", data$MS_type, "\n"),
                              lapply(1:nfleet, individual_matrix_fn, obs = "data$MS",
                                     pred = ifelse(data$MS_type == "length", "report$MLpred", "report$MWpred"),
                                     fig.cap = paste("mean", data$MS_type, "from fleet"), label = MS_label))
              } else MS_plots <- NULL

              if(any(data$s_CAA > 0, na.rm = TRUE)) {
                s_CAA_plots <- c("#### Survey age comps \n",
                                 lapply(1:nsurvey, individual_array_fn, obs = "data$s_CAA", pred = "report$s_CAApred", comps = "age", label = s_name))
              } else s_CAA_plots <- NULL

              if(any(data$s_CAL > 0, na.rm = TRUE)) {
                s_CAL_plots <- c("#### Survey length comps \n",
                                 lapply(1:nsurvey, individual_array_fn, obs = "data$s_CAL", pred = "report$s_CALpred", comps = "length", label = s_name))
              } else s_CAL_plots <- NULL

              data_section <- c(C_matplot, E_matplot, C_plots, I_plots, CAA_plots, CAL_plots, MS_plots, s_CAA_plots, s_CAL_plots)

              # Model output
              sel_matplot <- rmd_matplot(x = "matrix(age, max_age + 1, nfleet)", y = "matrix(report$vul[nyears, , ], max_age + 1, nfleet)", col = "rich.colors(nfleet)",
                                         xlab = "Age", ylab = "Selectivity", legend.lab = "f_name",
                                         fig.cap = "Terminal year selectivity by fleet.", header = "### Output \n")

              F_matplot <- rmd_matplot(x = "matrix(Year, nyears, nfleet)", y = "report$F", col = "rich.colors(nfleet)",
                                       xlab = "Year", ylab = "Fishing Mortality (F)", legend.lab = "f_name",
                                       fig.cap = "Time series of fishing mortality by fleet.")

              SSB <- structure(report$E, names = c(Year, max(Year) + 1))
              dynamic_SSB0 <- structure(report$dynamic_SSB0, names = names(SSB))
              SSB0 <- structure(report$E0, names = Year)
              if(length(unique(SSB0)) > 1) {
                SSB_plot <- rmd_assess_timeseries("SSB0", "unfished spawning depletion (growth and/or M are time-varying)",
                                                  "expression(SSB[0])")
              } else SSB_plot <- NULL

              SSB_SSB0 <- structure(report$E/report$E0_SR, names = c(Year, max(Year) + 1))

              R <- structure(report$R, names = c(Year, max(Year) + 1))
              N_at_age <- report$N
              N <- structure(rowSums(N_at_age), names = c(Year, max(Year) + 1))

              log_rec_dev <- structure(report$log_rec_dev, names = Year)
              log_rec_dev_SE <- summary(SD)[, 2]
              log_rec_dev_SE <- log_rec_dev_SE[names(log_rec_dev_SE) == "log_rec_dev"]
              log_rec_dev_SE <- ifelse(data_mean_fit$est_rec_dev, log_rec_dev_SE, 0)

              N_bubble <- rmd_bubble("c(Year, max(Year)+1)", "N_at_age", fig.cap = "Predicted abundance-at-age.")

              CAA_all <- apply(report$CAApred, c(1, 2), sum)
              CAA_bubble <- rmd_bubble("Year", "CAA_all", fig.cap = "Predicted catch-at-age (summed over all fleets).", bubble_adj = as.character(bubble_adj))

              CAL_all <- apply(report$CALpred, c(1, 2), sum)
              if(sum(CAL_all, na.rm = TRUE) > 0) {
                CAL_bubble <- rmd_bubble("Year", "CAL_all", CAL_bins = "data_mean_fit$length_bin",
                                         fig.cap = "Predicted catch-at-length (summed over all fleets).", bubble_adj = as.character(bubble_adj))
              } else CAL_bubble <- ""

              ts_output <- c(sel_matplot, F_matplot, rmd_SSB(), SSB_plot, rmd_SSB_SSB0(FALSE), rmd_dynamic_SSB0(), rmd_R(), rmd_RCM_SR(),
                             rmd_residual("log_rec_dev", fig.cap = "Time series of recruitment deviations.", label = "log-Recruitment deviations"),
                             rmd_residual("log_rec_dev", "log_rec_dev_SE", fig.cap = "Time series of recruitment deviations with 95% confidence intervals.",
                                          label = "log-Recruitment deviations", conv_check = TRUE),
                             rmd_N(), N_bubble, CAA_bubble, CAL_bubble)

              nll <- RCM_get_likelihoods(report, data$LWT, f_name, s_name)
              
              nll_table <- c("### Likelihood components\n",
                             "#### Summary\n",
                             "`r nll[[1]] ", ifelse(render_args$output_format == "html_document", "", "%>% knitr::kable(format = \"markdown\")"), "`\n\n",
                             "#### Fleet likelihoods\n",
                             "`r nll[[2]] ", ifelse(render_args$output_format == "html_document", "", "%>% knitr::kable(format = \"markdown\")"), "`\n\n",
                             "#### Fleet weights\n",
                             "`r nll[[3]] ", ifelse(render_args$output_format == "html_document", "", "%>% knitr::kable(format = \"markdown\")"), "`\n\n",
                             "#### Survey likelihoods\n",
                             "`r nll[[4]] ", ifelse(render_args$output_format == "html_document", "", "%>% knitr::kable(format = \"markdown\")"), "`\n\n",
                             "#### Survey weights\n",
                             "`r nll[[5]] ", ifelse(render_args$output_format == "html_document", "", "%>% knitr::kable(format = \"markdown\")"), "`\n\n")
              
              corr_matrix <- c("### Correlation matrix\n",
                               "`r SD$env$corr.fixed %>% as.data.frame()`\n\n")

              if(exists("retro", inherits = FALSE)) {
                ret <- rmd_RCM_retrospective(render_args)
              } else ret <- NULL

              mean_fit_rmd <- c(sumry, LH_section, data_section, ts_output, nll_table, corr_matrix, ret)
            } else mean_fit_rmd <- c("## Fit to mean parameters of OM {.tabset}\n",
                                     "No model found. Re-run `RCM()` with `mean_fit = TRUE`.\n\n")

            if(compare) {
              message("Getting Hist object from runMSE...")
              Hist <- runMSE(OM, Hist = TRUE, silent = TRUE, parallel = OM@nsim >= 48 & snowfall::sfIsRunning())

              compare_rmd <- c("## Updated OM {.tabset}\n",
                               "### OM historical period\n\n",
                               "```{r, fig.cap = \"Apical F from the operating model.\"}",
                               "Hist_F <- apply(Hist@AtAge$F.Mortality, c(1, 3), max, na.rm = TRUE)",
                               "matplot(Year, t(Hist_F), typ = \"l\", col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd, xlab = \"Year\", ylab = \"OM Apical F\", ylim = c(0, 1.1 * max(Hist_F)))",
                               "abline(h = 0, col = \"grey\")",
                               "if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd)",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Spawning biomass (SSB) from the operating model.\"}",
                               "Hist_SSB <- apply(Hist@TSdata$SBiomass, 1:2, sum)",
                               "matplot(Year, t(Hist_SSB), typ = \"l\", col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd, xlab = \"Year\", ylab = \"OM SSB\", ylim = c(0, 1.1 * max(Hist_SSB)))",
                               "abline(h = 0, col = \"grey\")",
                               "if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd)",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Spawning biomass (SSB) relative to MSY from the operating model.\"}",
                               "matplot(Year, t(Hist_SSB/Hist@Ref$ReferencePoints$SSBMSY), col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd, typ = \"l\", xlab = \"Year\", ylab = expression(OM~~SSB/SSB[MSY]), ylim = c(0, 1.1 * max(Hist_SSB/Hist@Ref$ReferencePoints$SSBMSY)))",
                               "abline(h = c(0, MSY_ref), col = \"grey\")",
                               "if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd)",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Spawning biomass (SSB) relative to MSY in the most recent decade.\"}",
                               "if(length(Year) > 10) {",
                               "  Yr_ind <- Year > max(Year) - 10",
                               "  matplot(Year[Yr_ind], t(Hist_SSB[, Yr_ind, drop = FALSE]/Hist@Ref$ReferencePoints$SSBMSY), typ = \"l\", col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd, xlab = \"Year\", ylab = expression(OM~~SSB/SSB[MSY]), ylim = c(0, 1.1 * max(Hist_SSB[, Yr_ind, drop = FALSE]/Hist@Ref$ReferencePoints$SSBMSY)))",
                               "  abline(h = c(0, MSY_ref), col = \"grey\")",
                               "  if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd)",
                               "}",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Spawning depletion from the operating model.\"}",
                               "matplot(Year, t(Hist_SSB/Hist@Ref$ReferencePoints$SSB0), typ = \"l\", col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd, xlab = \"Year\", ylab = expression(OM~~SSB/SSB[0]), ylim = c(0, 1.1 * max(Hist_SSB/Hist@Ref$ReferencePoints$SSB0)))",
                               "abline(h = 0, col = \"grey\")",
                               "if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd)",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Recruitment (age-0) from the operating model.\"}",
                               "Hist_R <- Hist@AtAge$Number[, 1, , ] %>% apply(c(1, 2), sum)",
                               "matplot(Year, t(Hist_R), typ = \"l\", col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd, xlab = \"Year\", ylab = \"OM Recruitment\", ylim = c(0, 1.1 * max(Hist_R)))",
                               "abline(h = 0, col = \"grey\")",
                               "if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd)",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Catch (total removals, including discards) from the operating model.\"}",
                               "Hist_C <- apply(Hist@TSdata$Removals, 1:2, sum)",
                               "matplot(Year, t(Hist_C), typ = \"l\", col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd, xlab = \"Year\", ylab = \"OM Catch\", ylim = c(0, 1.1 * max(Hist_C)))",
                               "abline(h = 0, col = \"grey\")",
                               "if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd)",
                               "```\n",
                               "",
                               "### OM/RCM Comparison\n\n",
                               "```{r, fig.cap = \"Apical F comparison between the OM and RCM.\"}",
                               "matplot(Year, t(Hist_F), typ = \"o\", pch = 16, col = \"red\", xlab = \"Year\", ylab = \"Apical F\", ylim = c(0, 1.1 * max(c(Hist_F, OM@cpars$Find))))",
                               "matlines(Year, t(OM@cpars$Find), col = \"black\")",
                               "abline(h = 0, col = \"grey\")",
                               "legend(\"topleft\", c(\"RCM\", \"OM\"), col = c(\"black\", \"red\"), pch = c(NA_integer_, 16), lwd = c(1, 1), bty = \"n\")",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Difference in apical F between the OM and RCM. Positive values indicate higher F in the OM.\"}",
                               "matplot(Year, t(Hist_F - OM@cpars$Find), typ = \"n\", xlab = \"Year\", ylab = \"Difference in apical F\")",
                               "abline(h = 0, col = \"grey\")",
                               "matlines(Year, t(Hist_F - OM@cpars$Find), col = \"black\")",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"SSB comparison between the OM and RCM.\"}",
                               "matplot(Year, t(Hist_SSB), typ = \"o\", col = \"red\", pch = 16, xlab = \"Year\", ylab = \"SSB\",",
                               "        ylim = c(0, 1.1 * max(c(Hist_SSB, x@SSB[sims, 1:OM@nyears]))))",
                               "matlines(Year, t(x@SSB[sims, 1:OM@nyears, drop = FALSE]), col = \"black\")",
                               "abline(h = 0, col = \"grey\")",
                               "legend(\"topleft\", c(\"RCM\", \"OM\"), col = c(\"black\", \"red\"), pch = c(NA_integer_, 16), lwd = c(1, 1), bty = \"n\")",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Difference in spawning biomass (SSB), relative to SSB0, between the OM and RCM, calculated as $(SSB^{OM}_y - SSB^{RCM}_y)/SSB^{OM}_0$. Positive values indicate higher SSB in the OM.\"}",
                               "matplot(Year, t((Hist_SSB - x@SSB[sims, 1:OM@nyears, drop = FALSE])/Hist@Ref$ReferencePoints$SSB0), typ = \"n\", xlab = \"Year\", ylab = \"Difference in relative SSB\")",
                               "abline(h = 0, col = \"grey\")",
                               "matlines(Year, t((Hist_SSB - x@SSB[sims, 1:OM@nyears, drop = FALSE])/Hist@Ref$ReferencePoints$SSB0), col = \"black\")",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Recruitment comparison between the OM and RCM.\"}",
                               "matplot(Year, t(Hist_R), typ = \"o\", col = \"red\", pch = 16, xlab = \"Year\", ylab = \"Recruitment\",",
                               "        ylim = c(0, 1.1 * max(c(Hist_R, x@NAA[sims, 1:OM@nyears, 1]))))",
                               "matlines(Year, t(x@NAA[, 1:OM@nyears, 1][sims, , drop = FALSE]), col = \"black\")",
                               "abline(h = 0, col = \"grey\")",
                               "legend(\"topleft\", c(\"RCM\", \"OM\"), col = c(\"black\", \"red\"), pch = c(NA_integer_, 16), lwd = c(1, 1), bty = \"n\")",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Difference in recruitment (relative to R0) between the OM and RCM, calculated as $(R^{OM}_y - R^{RCM}_y)/R^{OM}_0$. Positive values indicate higher recruitment in the OM.\"}",
                               "matplot(Year, t(Hist_R/OM@cpars$R0 - x@NAA[, 1:OM@nyears, 1][sims, , drop = FALSE]/OM@cpars$R0), typ = \"n\", xlab = \"Year\", ylab = \"Difference in relative recruitment\")",
                               "abline(h = 0, col = \"grey\")",
                               "matlines(Year, t(Hist_R/OM@cpars$R0 - x@NAA[, 1:OM@nyears, 1][sims, , drop = FALSE]/OM@cpars$R0),",
                               "         col = \"black\")",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Comparison of total removals between the OM and RCM.\"}",
                               "matplot(Year, t(Hist_C), typ = \"o\", col = \"red\", pch = 16, xlab = \"Year\", ylab = \"Total removals\",",
                               "        ylim = c(0, 1.1 * max(c(Hist_C, data$Chist, na.rm = TRUE))))",
                               "lines(Year, rowSums(data$Chist, na.rm = TRUE), col = \"black\")",
                               "abline(h = 0, col = \"grey\")",
                               "legend(\"topleft\", c(\"RCM\", \"OM\"), col = c(\"black\", \"red\"), pch = c(NA_integer_, 16), lwd = c(1, 1), bty = \"n\")",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Difference in annual catch (relative to observed), calculated as $C^{OM}_y/C^{obs}_y - 1$. Positive values indicate higher catch in the OM. Catch in the OM is the total removals (both landings and discards).\"}",
                               "if(any(data$Chist > 0, na.rm = TRUE)) {",
                               "Catch_difference <- t(Hist_C)/rowSums(data$Chist, na.rm = TRUE) - 1",
                               "Catch_difference[is.infinite(Catch_difference)] <- 0",
                               "matplot(Year, Catch_difference, typ = \"n\", xlab = \"Year\", ylab = \"Difference in relative catch\")",
                               "abline(h = 0, col = \"grey\")",
                               "matlines(Year, Catch_difference, col = \"black\")",
                               "}",
                               "```\n")

            } else compare_rmd <- c("## Updated OM\n",
                                    "Re-run `plot()` function with argument `compare = TRUE`.\n\n")

            rmd <- c(header, OM_update, all_sims_output, mean_fit_rmd, compare_rmd, rmd_footer())
            if(is.list(rmd)) rmd <- do.call(c, rmd)

            # Generate markdown report
            if(!dir.exists(dir)) {
              message("Creating directory: \n", dir)
              dir.create(dir)
            }
            write(rmd, file = file.path(dir, filename_rmd))
            message("Generated markdown file: ", file.path(dir, filename_rmd))

            # Rendering markdown file
            message("Rendering markdown file...")
            output_filename <- do.call(rmarkdown::render, render_args)
            message("Rendered file: ", output_filename)

            if(open_file) browseURL(output_filename)
            invisible(output_filename)
          })



rmd_persp_plot <- function(x, y, z, xlab, ylab, zlab, phi, theta, expand, fig.cap, header = NULL) {
  ans <- c(paste0("```{r, fig.cap = \"", fig.cap, "\"}"),
           paste0("persp(x = ", x, ", y = ", y, ", z = ", z, ", theta = ", theta, ", phi = ", phi, ", expand = ", expand, ", xlab = \"", xlab, "\",
                   ylab = \"", ylab, "\", zlab = \"", zlab, "\", ticktype = \"detailed\")"),
           " ```\n")
  if(!is.null(header)) ans <- c(header, ans)
  return(ans)
}

rmd_matplot <- function(x, y, col, xlab, ylab, legend.lab, type = "l", lty = 1, fig.cap, header = NULL) {
  ans <- c(paste0("```{r, fig.cap = \"", fig.cap, "\"}"),
           paste0("xx <- ", x, "; yy <- ", y),
           paste0("matplot(xx, yy, type = \"", type, "\", lty = ", lty, ", col = ", col,
                  ", ylim = c(0, 1.1 * max(yy, na.rm = TRUE)), xlab = \"", xlab, "\", ylab = \"", ylab, "\")"),
           "abline(h = 0, col = \"grey\")",
           paste0("if(ncol(xx) > 1) legend(\"topleft\", ", legend.lab, ", text.col = ", col, ")"),
           " ```\n")
  if(!is.null(header)) ans <- c(header, ans)
  return(ans)
}

# For RCM function
rmd_assess_fit2 <- function(year, obs, fit, fig.cap, label = fig.cap, match = FALSE) {
  fig.cap2 <- paste0("Observed (black) and predicted (red) ", fig.cap, ".")
  if(match) fig.cap2 <- paste(fig.cap2, "Predicted", fig.cap, "should match observed in this model.")

  c(paste0("```{r, fig.cap = \"", fig.cap2, "\"}"),
    paste0("plot_timeseries(", year, ", ", obs, ", ", fit, ", label = \"", label, "\")"),
    "```\n")
}

# For RCM function
rmd_assess_resid2 <- function(year, obs, fit, fig.cap, label = fig.cap) {
  fig.cap2 <- paste0("Index residuals (in log space) for ", fig.cap, ".")
  
  c(paste0("```{r, fig.cap = \"", fig.cap2, "\"}"),
    paste0("istart <- which(!is.na(", obs, "))[1]"),
    paste0("istop <- which(!is.na(", obs, ")) %>% max()"),
    paste0("plot_residuals(", year, "[istart:istop], log(", obs, "[istart:istop]/", fit, "[istart:istop]), label = \"", label, "\")"),
    "```\n")
}



rmd_fit_comps <- function(year, obs, fit, type = c("bubble_data", "annual", "bubble_residuals", "mean"), ages = "NULL", CAL_bins = "NULL", fig.cap,
                          bubble_adj = "10") {
  type <- match.arg(type)
  arg <- paste0("\"", type, "\", CAL_bins = ", CAL_bins, ", ages = ", ages)
  c(paste0("```{r, fig.cap = \"", fig.cap, "\"}"),
    paste0("ind_valid <- rowSums(", obs, ", na.rm = TRUE) > 0"),
    paste0("if(any(ind_valid)) plot_composition(", year, "[ind_valid], ", obs, "[ind_valid, , drop = FALSE], ", fit, "[ind_valid, , drop = FALSE], plot_type = ", arg, ", bubble_adj = ", bubble_adj, ")"),
    "```\n")
}

rmd_RCM_R0 <- function(fig.cap = "Histogram of R0 (unfished recruitment).") {
  c(paste0("```{r, fig.cap = \"", fig.cap, "\"}"),
    "if(!is.null(OM@cpars$R0)) hist(OM@cpars$R0, main = \"\", xlab = expression(R[0]))",
    "```\n")
}

rmd_RCM_D <- function(fig.cap = "Histogram of historical depletion.") {
  c(paste0("```{r, fig.cap = \"", fig.cap, "\"}"),
    "if(!is.null(OM@cpars$D)) hist(OM@cpars$D, main = \"\", xlab = \"Depletion\")",
    "```\n")
}

rmd_RCM_Perr <- function(fig.cap = "Historical recruitment deviations among simulations.") {
  c(paste0("```{r, fig.cap = \"", fig.cap, "\"}"),
    "Perr <- OM@cpars$Perr_y[, (max_age+1):(max_age+nyears), drop = FALSE]",
    "matplot(Year, t(Perr), type = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Recruitment deviations\",",
    "        ylim = c(0, 1.1 * max(Perr)))",
    "abline(h = 0, col = \"grey\")",
    "```\n",
    "",
    "```{r, fig.cap = \"Future recruitment deviations (up to 5 simulations).\"}",
    "Perr_future <- OM@cpars$Perr_y[, (max_age+nyears+1):(max_age+nyears+proyears)]",
    "matplot(Year, t(Perr), type = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Recruitment deviations\",",
    "        xlim = c(min(Year), max(Year) + proyears), ylim = c(0, 1.1 * max(c(Perr, Perr_future))))",
    "matlines(max(Year) + 1:proyears, t(Perr_future[1:min(5, nrow(OM@cpars$Perr_y)), ]), type = \"l\")",
    "abline(h = 0, col = \"grey\")",
    "abline(v = max(Year), lty = 3)",
    "```\n",
    "",
    "```{r, fig.cap = \"Annual mean and median of future recruitment deviations.\"}",
    "matplot(Year, t(Perr), type = \"n\", xlab = \"Year\", ylab = \"Recruitment deviations\",",
    "        xlim = c(min(Year), max(Year) + proyears), ylim = c(0, 1.1 * max(c(Perr, Perr_future))))",
    "abline(h = c(0, 1), col = \"grey\")",
    "abline(v = max(Year), lty = 3)",
    "matlines(Year, t(Perr), type = \"l\", col = \"black\")",
    "lines(max(Year) + 1:proyears, apply(Perr_future, 2, mean), col = \"red\")",
    "lines(max(Year) + 1:proyears, apply(Perr_future, 2, median), lty = 2)",
    "legend(\"topleft\", c(\"Mean\", \"Median\"), col = c(\"red\", \"black\"), lty = c(1, 2))",
    "```\n",
    "```{r, fig.cap = \"Histogram of recruitment autocorrelation.\"}",
    "if(!is.null(OM@cpars$AC)) hist(OM@cpars$AC, main = \"\", xlab = \"Recruitment Autocorrelation\")",
    "```\n")
}

rmd_RCM_Find <- function(fig.cap = "Apical F from RCM model. These values may be subsequently re-scaled in the operating model in order to match the specified depletion") {
  c(paste0("```{r, fig.cap = \"", fig.cap, "\"}"),
    "matplot(Year, t(OM@cpars$Find), type = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Apical F\")",
    "abline(h = 0, col = \"grey\")",
    "```\n")
}

rmd_RCM_sel <- function(fig.cap = "Operating model selectivity among simulations.") {
  c(paste0("```{r, fig.cap = \"", fig.cap, "\"}"),
    "if(nsel_block == 1) {",
    "  vul <- do.call(cbind, lapply(report_list, getElement, \"vul_len\"))",
    "  matplot(length_bin, vul, type = \"l\", col = \"black\",",
    "          xlab = \"Length\", ylab = \"Selectivity (last historical year)\", ylim = c(0, 1.1))",
    "} else {",
    "  if(nsim == 1) V_plot <- matrix(OM@cpars$V[, , nyears], 1, byrow = TRUE) else V_plot <- OM@cpars$V[, , nyears]",
    "  matplot(age, t(V_plot), type = \"l\", col = \"black\",",
    "          xlab = \"Age\", ylab = \"Selectivity (last historical year)\", ylim = c(0, 1.1))",
    "}",
    "abline(h = 0, col = \"grey\")",
    "```\n")
}

rmd_RCM_fleet_output <- function(ff, f_name) {
  if(ff == 1) header <- "## RCM output {.tabset}\n" else header <- NULL
  ans <- c(paste("### ", f_name[ff], "\n"),
           paste0("```{r, fig.cap = \"Selectivity of ", f_name[ff], ".\"}"),
           paste0("bl <- unique(data$sel_block[, ", ff, "])"),
           "vul_bb <- list()",
           "bl_col <- gplots::rich.colors(length(bl))",
           "Year_legend <- character(length(bl))",
           "for(bb in 1:length(bl)) {",
           "  vul_bb[[bb]] <- do.call(cbind, lapply(report_list, function(x) x$vul_len[, bl[bb]]))",
           paste0("  Year_legend[bb] <- Year[data$sel_block[, ", ff, "] == bl[bb]] %>% range() %>% paste(collapse = \"-\")"),
           "}",
           "test <- vapply(vul_bb, function(x) all(!is.na(x)), logical(1))",
           "if(all(test)) {",
           paste0("  matplot(length_bin, length_bin, type = \"n\", xlab = \"Length\", ylim = c(0, 1), ylab = \"Selectivity of Fleet ", ff, "\")"),
           "  abline(h = 0, col = \"grey\")",
           "  for(bb in 1:length(bl)) {",
           "    matlines(length_bin, vul_bb[[bb]], type = \"l\", col = bl_col[bb], lty = scenario$lty, lwd = scenario$lwd)",
           "  }",
           "  if(length(bl) > 1) legend(\"topright\", Year_legend, col = bl_col, lwd = 1)",
           #"if(!is.null(scenario$names)) legend("topleft", scenario$names, col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd)",
           "}",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Corresponding age-based selectivity of ", f_name[ff], ".\"}"),
           paste0("matplot(age, age, type = \"n\", xlab = \"Age\", ylim = c(0, 1), ylab = \"Selectivity of Fleet ", ff, "\")"),
           "abline(h = 0, col = \"grey\")",
           "",
           "for(bb in 1:length(bl)) {",
           paste0("  vul_bb_age <- do.call(rbind, lapply(report_list, function(x) x$vul[data$sel_block[, ", ff, "] == bl[bb], , ", ff, "])) %>% t()"),
           "  matlines(age, vul_bb_age, type = \"l\", col = bl_col[bb], lty = scenario$lty, lwd = scenario$lwd)",
           "}",
           "if(length(bl) > 1) legend(\"topleft\", Year_legend, col = bl_col, lwd = 1)",
           #"if(!is.null(scenario$names)) legend("topleft", scenario$names, col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd)",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Fishing Mortality of ", f_name[ff], ".\"}"),
           paste0("FM <- do.call(cbind, lapply(report_list, function(x) x$F[, ", ff, "]))"),
           paste0("matplot(Year, FM, type = \"l\", col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd, ylim = c(0, 1.1 * max(FM, na.rm = TRUE)), xlab = \"Year\", ylab = \"Fishing Mortality of ", f_name[ff], "\")"),
           "abline(h = 0, col = \"grey\")",
           "if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd)",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (black) and predicted (red) catch from ", f_name[ff], ".\"}"),
           paste0("if(any(data$Chist[, ", ff, "] > 0)) {"),
           paste0("  Cpred <- do.call(cbind, lapply(report_list, function(x) x$Cpred[, ", ff, "]))"),
           paste0("  Chist <- data$Chist[, ", ff, "]"),
           "  ylim <- c(0.9, 1.1) * range(c(Cpred, Chist), na.rm = TRUE)",
           paste0("  plot(Year, Chist, type = \"o\", xlab = \"Year\", ylab = \"Catch of ", f_name[ff], "\", ylim = ylim)"),
           paste0("  matlines(Year, Cpred, col = scenario$col, lty = scenario$lty, lwd = scenario$lwd)"),
           "} else {",
           paste0("  Cpred <- do.call(cbind, lapply(report_list, function(x) x$Cpred[, ", ff, "]/mean(x$Cpred[, ", ff, "])))"),
           paste0("  matplot(Year, Cpred, col = scenario$col, type = \"l\", lty = scenario$lty, lwd = scenario$lwd, xlab = \"Year\", ylab = \"Predicted relative catch of ", f_name[ff], "\")"),
           "}",
           "abline(h = 0, col = \"grey\")",
           "if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col, lty = scenario$lty, lwd = scenario$lwd)",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (black) and predicted (red) mean ages from ", f_name[ff], ".\"}"),
           paste0("MApred <- do.call(cbind, lapply(report_list, function(x) x$CAApred[, , ", ff, "] %*% age/x$CN[, ", ff, "]))"),
           paste0("MAobs <- (data$CAA[, , ", ff, "] %*% age)/rowSums(data$CAA[, , ", ff, "], na.rm = TRUE)"),
           "ylim <- c(0.9, 1.1) * range(c(MApred, MAobs), na.rm = TRUE)",
           "matplot(Year, MApred, type = \"l\", col = scenario$col, lty = scenario$lty, lwd = scenario$lwd, xlab = \"Year\", ylab = \"Mean age\", ylim = ylim)",
           paste0("if(any(data$CAA[, , ", ff, "] > 0, na.rm = TRUE)) {"),
           paste0("  lines(Year, MAobs, col = \"black\", typ = \"o\")"),
           "}",
           "if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col, lty = scenario$lty, lwd = scenario$lwd)",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (black) and predicted (red) mean lengths from ", f_name[ff], ".\"}"),
           paste0("MLpred <- do.call(cbind, lapply(report_list, function(x) x$MLpred[, ", ff, "]))"),
           paste0("if(any(data$CAL[, , ", ff, "] > 0, na.rm = TRUE)) {"),
           paste0("  MLobs <- (data$CAL[, , ", ff, "] %*% length_bin)/rowSums(data$CAL[, , ", ff, "], na.rm = TRUE)"),
           paste0("} else if(data$MS_type == \"length\" && any(data$MS[, ", ff, "] > 0, na.rm = TRUE)) MLobs <- data$MS[, ", ff, "] else MLobs <- NA"),
           "if(!all(is.na(MLpred))) {",
           "  ylim <- c(0.9, 1.1) * range(c(MLpred, MLobs), na.rm = TRUE)",
           "  matplot(Year, MLpred, type = \"l\", col = scenario$col, lty = scenario$lty, lwd = scenario$lwd, xlab = \"Year\", ylab = \"Mean length\", ylim = ylim)",
           "  if(!all(is.na(MLobs))) lines(Year, MLobs, col = \"black\", typ = \"o\")",
           "  if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col, lty = scenario$lty, lwd = scenario$lwd)",
           "}",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (black) and predicted (red) mean weights from ", f_name[ff], ".\"}"),
           paste0("if(data$MS_type == \"weight\" && any(data$MS[, ", ff, "] > 0, na.rm = TRUE)) {"),
           paste0("  MWobs <- data$MS[, ", ff, "]"),
           paste0("} else MWobs <- NA"),
           "if(!all(is.na(MWobs))) {",
           paste0("  MWpred <- do.call(cbind, lapply(report_list, function(x) x$MWpred[, ", ff, "]))"),
           "  ylim <- c(0.9, 1.1) * range(c(MWpred, MWobs), na.rm = TRUE)",
           "  matplot(Year, MWpred, type = \"l\", col = scenario$col, lty = scenario$lty, lwd = scenario$lwd, xlab = \"Year\", ylab = \"Mean weight\", ylim = ylim)",
           "  lines(Year, MWobs, col = \"black\", typ = \"o\")",
           "  if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col, lty = scenario$lty, lwd = scenario$lwd)",
           "}",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (black) and predicted (red) age composition from ", f_name[ff], ".\"}"),
           paste0("if(any(data$CAA[, , ", ff, "] > 0, na.rm = TRUE)) {"),
           paste0("if(nsim == 1) CAA_plot <- array(x@CAA[, , , ", ff, "], c(1, nyears, max_age + 1)) else CAA_plot <- x@CAA[, , , ", ff, "]"),
           paste0("plot_composition_RCM(Year, CAA_plot, data$CAA[, , ", ff, "], ages = age, dat_col = scenario$col)"),
           "}",
           "```\n",
           paste0("```{r, fig.cap = \"Predicted age composition from fleet ", ff, ".\"}"),
           paste0("if(any(data$CAA[, , ", ff, "] > 0, na.rm = TRUE)) {"),
           paste0("plot_composition_RCM(Year, CAA_plot, ages = age, dat_col = scenario$col)"),
           "}",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (black) and predicted (red) length composition from ", f_name[ff], ".\"}"),
           paste0("if(any(data$CAL[, , ", ff, "] > 0, na.rm = TRUE)) {"),
           paste0("if(nsim == 1) CAL_plot <- array(x@CAL[, , , ", ff, "], c(1, nyears, length(data$length_bin))) else CAL_plot <- x@CAL[, , , ", ff, "]"),
           paste0("plot_composition_RCM(Year, CAL_plot, data$CAL[, , ", ff, "], CAL_bins = data$length_bin, dat_col = scenario$col)"),
           "}",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Predicted length composition from ", f_name[ff], ".\"}"),
           paste0("if(any(data$CAL[, , ", ff, "] > 0, na.rm = TRUE)) {"),
           paste0("plot_composition_RCM(Year, CAL_plot, CAL_bins = data$length_bin, dat_col = scenario$col)"),
           "}",
           "```\n")

  c(header, ans)
}

rmd_RCM_survey_output <- function(sur, s_name) {
  ans <- c(paste0("### ", s_name[sur], " \n"),
           "",
           paste0("```{r, fig.cap = \"Selectivity of ", s_name[sur], " in last historical year.\"}"),
           "if(!is.null(report_list[[1]]$s_vul)) {",
           paste0("s_vul_ff_age <- do.call(cbind, lapply(report_list, function(x) x$s_vul[nyears, , ", sur, "]))"),
           paste0("matplot(age, s_vul_ff_age, type = \"l\", col = scenario$col2, xlab = \"Age\", ylim = c(0, 1), ylab = \"Selectivity of ", s_name[sur], "\")"),
           "abline(h = 0, col = \"grey\")",
           "if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd)",
           "}",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (black) and predicted (red) index values for ", s_name[sur], ".\"}"),
           paste0("Ipred <- do.call(cbind, lapply(report_list, function(x) x$Ipred[, ", sur, "]))"),
           paste0("matplot(Year, Ipred, type = \"l\", col = scenario$col, lty = scenario$lty, lwd = scenario$lwd, ylim = c(0, 1.1 * max(c(Ipred, data$Index[, ", sur, "]), na.rm = TRUE)), xlab = \"Year\", ylab = \"", s_name[sur], "\")"),
           paste0("lines(Year, data$Index[, ", sur, "], col = \"black\", typ = \"o\")"),
           "abline(h = 0, col = \"grey\")",
           "if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col, lty = scenario$lty, lwd = scenario$lwd)",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (black) and predicted (red) index values for ", s_name[sur], ". Error bars indicate 95% confidence intervals for observed values.\"}"),
           "if(!is.null(data$I_sd)) {",
           paste0("  II <- data$Index[, ", sur, "]"),
           "  ind <- seq(min(which(!is.na(II))), max(which(!is.na(II))), 1)",
           paste0("  err <- exp(log(II) + outer(data$I_sd[, ", sur, "], c(-1.96, 1.96)))"),
           paste0("  matplot(Year[ind], Ipred[ind, ], type = \"l\", col = scenario$col, lty = scenario$lty, lwd = scenario$lwd, ylim = c(0, 1.1 * max(c(Ipred[ind, ], II[ind], err[ind, ]), na.rm = TRUE)), xlab = \"Year\", ylab = \"", s_name[sur], "\")"),
           "  points(Year[ind], II[ind], lwd = 3, pch = 16)",
           "  arrows(Year[ind], y0 = err[ind, 1], y1 = err[ind, 2], length = 0, lwd = 1.5)",
           "  abline(h = 0, col = \"grey\")",
           "  if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col, lty = scenario$lty, lwd = scenario$lwd)",
           "}",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (black) and predicted (red) mean ages from ", s_name[sur], ".\"}"),
           paste0("if(!is.null(data$s_CAA) && any(data$s_CAA[, , ", sur, "] > 0, na.rm = TRUE)) {"),
           paste0("MApred <- do.call(cbind, lapply(report_list, function(x) x$s_CAApred[, , ", sur, "] %*% age/rowSums(x$N[1:nyears, ] * x$s_vul[1:nyears, , ", sur, "])))"),
           paste0("MAobs <- (data$s_CAA[, , ", sur, "] %*% age)/rowSums(data$s_CAA[, , ", sur, "], na.rm = TRUE)"),
           "ylim <- c(0.9, 1.1) * range(c(MApred, MAobs), na.rm = TRUE)",
           "matplot(Year, MApred, type = \"l\", col = scenario$col, lty = scenario$lty, lwd = scenario$lwd, xlab = \"Year\", ylab = \"Mean age\", ylim = ylim)",
           paste0("if(any(data$s_CAA[, , ", sur, "] > 0, na.rm = TRUE)) {"),
           paste0("  lines(Year, MAobs, col = \"black\", typ = \"o\")"),
           "}",
           "if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col, lty = scenario$lty, lwd = scenario$lwd)",
           "}",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (black) and predicted (red) mean lengths from ", s_name[sur], ".\"}"),
           paste0("if(!is.null(data$s_CAL) && any(data$s_CAL[, , ", sur, "] > 0, na.rm = TRUE)) {"),
           paste0("MLpred <- do.call(cbind, lapply(report_list, function(x) x$s_CALpred[, , ", sur, "] %*% length_bin/rowSums(x$N[1:nyears, ] * x$s_vul[1:nyears, , ", sur, "])))"),
           paste0("MLobs <- (data$s_CAL[, , ", sur, "] %*% length_bin)/rowSums(data$s_CAL[, , ", sur, "], na.rm = TRUE)"),
           "ylim <- c(0.9, 1.1) * range(c(MLpred, MLobs), na.rm = TRUE)",
           "matplot(Year, MLpred, type = \"l\", col = scenario$col, lty = scenario$lty, lwd = scenario$lwd, xlab = \"Year\", ylab = \"Mean length\", ylim = ylim)",
           paste0("if(any(data$s_CAL[, , ", sur, "] > 0, na.rm = TRUE)) {"),
           paste0("  lines(Year, MLobs, col = \"black\", typ = \"o\")"),
           "}",
           "if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col, lty = scenario$lty, lwd = scenario$lwd)",
           "}",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (black) and predicted (red) age composition from ", s_name[sur], ".\"}"),
           paste0("if(!is.null(data$s_CAA) && any(data$s_CAA[, , ", sur, "] > 0, na.rm = TRUE)) {"),
           paste0("pred_sCAA <- lapply(report_list, function(x) x$s_CAA[,, ", sur, "]) %>% simplify2array() %>% aperm(perm = c(3, 1, 2))"),
           paste0("plot_composition_RCM(Year, pred_sCAA, data$s_CAA[, , ", sur, "], ages = age, dat_col = scenario$col)"),
           "}",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (black) and predicted (red) length composition from ", s_name[sur], ".\"}"),
           paste0("if(!is.null(data$s_CAL) && any(data$s_CAL[, , ", sur, "] > 0, na.rm = TRUE)) {"),
           paste0("pred_sCAL <- lapply(report_list, function(x) x$s_CAL[,, ", sur, "]) %>% simplify2array() %>% aperm(perm = c(3, 1, 2))"),
           paste0("plot_composition_RCM(Year, pred_sCAL, data$s_CAL[, , ", sur, "], CAL_bins = data$length_bin, dat_col = scenario$col)"),
           "}",
           "```\n"
           )
  ans
}

rmd_RCM_initD <- function(fig.cap = "Histogram of initial depletion among all simulations.") {
  c(paste0("```{r, fig.cap = \"", fig.cap, "\"}"),
    "initD <- vapply(report_list, function(x) x$E[1]/x$E0[1], numeric(1))",
    "hist(initD, main = \"\", xlab = \"Initial depletion\")",
    "```\n")
}

rmd_RCM_R_output <- function() {
  c("```{r, fig.cap = \"Estimated recruitment among all simulations.\"}",
    "R_out <- do.call(cbind, lapply(report_list, getElement, \"R\"))",
    "matplot(Yearplusone, R_out, ylim = c(0, 1.1 * max(R_out, na.rm = TRUE)), type = \"l\", col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd, xlab = \"Year\", ylab = \"Recruitment\")",
    "abline(h = 0, col = \"grey\")",
	"if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd)",
    "```\n")
}

rmd_RCM_SSB_output <- function() {
  c("```{r, fig.cap = \"Estimated spawning biomass among all simulations.\"}",
    "E <- do.call(cbind, lapply(report_list, getElement, \"E\"))",
    "matplot(Yearplusone, E, ylim = c(0, 1.1 * max(E, na.rm = TRUE)), type = \"l\", col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd, xlab = \"Year\", ylab = \"Spawning biomass\")",
    "abline(h = 0, col = \"grey\")",
    "if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd)",
    "```\n",
    "",
    "```{r, fig.cap = \"Estimated spawning depletion among all simulations. Unfished spawning biomass is the value calculated from first year life history parameters.\"}",
    "E_E0 <- do.call(cbind, lapply(report_list, function(x) x$E/x$E0_SR))",
    "matplot(Yearplusone, E_E0, ylim = c(0, 1.1 * max(E_E0, na.rm = TRUE)), type = \"l\", col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd, xlab = \"Year\", ylab = \"Spawning depletion\")",
    "abline(h = 0, col = \"grey\")",
    "if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd)",
    "```\n",
    "",
    "```{r, fig.cap = \"Dynamic SSB0 among all simulations. Model is re-run assuming no historical catches.\"}",
    "dyn_SSB0 <- do.call(cbind, lapply(report_list, function(x) x$dynamic_SSB0))",
    "matplot(Yearplusone, dyn_SSB0, ylim = c(0, 1.1 * max(dyn_SSB0, na.rm = TRUE)), type = \"l\", col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd, xlab = \"Year\", ylab = expression(Dynamic~~SSB[0]))",
    "abline(h = 0, col = \"grey\")",
    "if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd)",
    "```\n")
}

rmd_log_rec_dev <- function() {
  c("```{r, fig.cap = \"Estimated recruitment deviations among all simulations.\"}",
    "log_rec_dev2 <- do.call(cbind, lapply(report_list, getElement, \"log_rec_dev\"))",
    "matplot(Year, log_rec_dev2, type = \"n\", xlab = \"Year\", ylab = \"log-recruitment deviations\")",
    "abline(h = 0, col = \"grey\")",
    "matlines(Year, log_rec_dev2, col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd)",
	"if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd)",
    "```\n")
}

rmd_RCM_SR <- function() {
  c("```{r, fig.cap = \"Stock-recruit relationship and estimated recruitment.\"}",
    "if(OM@SRrel == 1) {",
    "  expectedR <- report$Arec * report$E[1:nyears] / (1 + report$Brec * report$E[1:nyears])",
    "} else {",
    "  expectedR <- report$Arec * report$E[1:nyears] * exp(-report$Brec * report$E[1:nyears])",
    "}",
    "plot_SR(report$E[1:nyears], expectedR, report$R0, report$E0_SR, report$R[2:(nyears+1)])",
    "```\n")
}

rmd_RCM_retrospective <- function(render_args) {
  if(render_args$output_format == "html_document") {
    x <- "summary(retro) %>% as.data.frame()"
  } else {
    x <- "summary(retro) %>% as.data.frame() %>% knitr::kable(format = \"markdown\")"
  }

  c("### Retrospective\n",
    "```{r}",
    x,
    "plot(retro)",
    "```\n")
}

plot_composition_RCM <- function(Year, fit, dat = NULL, CAL_bins = NULL, ages = NULL, annual_ylab = "Frequency",
                                 annual_yscale = c("proportions", "raw"), N = if(is.null(dat)) NULL else round(rowSums(dat)), dat_linewidth = 2, dat_color = "black") {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfcol = c(4, 4), mar = rep(0, 4), oma = c(5.1, 5.1, 2.1, 2.1))

  annual_yscale <- match.arg(annual_yscale)
  if(is.null(CAL_bins)) data_type <- "age" else data_type <- "length"

  if(data_type == 'length') {
    data_val <- CAL_bins
    data_lab <- "Length"
  }
  if(data_type == 'age') {
    data_val <- if(is.null(ages)) 1:dim(fit)[3] else ages
    data_lab <- "Age"
  }

  # Annual comps (fit vs. dat if available)
  # Dim of
  fit_plot <- fit
  dat_plot <- dat
  if(annual_yscale == "proportions") {
    for(i in 1:length(Year)) {
      fit_plot[, i, ] <- fit[, i, ]/rowSums(fit[, i, ])
      if(!is.null(dat)) dat_plot[i, ] <- dat[i, ]/sum(dat[i, ])
    }
  }
  ylim <- c(0, 1.1 * max(fit_plot, dat_plot, na.rm = TRUE))
  yaxp <- c(0, max(pretty(ylim, n = 4)), 4)
  if(max(c(fit_plot, dat_plot), na.rm = TRUE) == 1) yaxp <- c(0, 1, 4)

  las <- 1

  for(i in 1:length(Year)) {
    yaxt <- ifelse(i %% 16 %in% c(1:4), "s", "n") # TRUE = first column
    xaxt <- ifelse(i < length(Year) & i %% 4 %in% c(1:3), "n", "s") # TRUE = first three rows

    if(dim(fit_plot)[1] == 1) {
      plot(data_val, fit_plot[, i, ], type = "n", ylim = ylim, yaxp = yaxp, xaxt = xaxt, yaxt = yaxt, las = las)
      abline(h = 0, col = 'grey')
      lines(data_val, fit_plot[, i, ], col = dat_color)
    } else {
      matplot(data_val, t(fit_plot[, i, ]), type = "n", ylim = ylim, yaxp = yaxp, xaxt = xaxt, yaxt = yaxt, las = las)
      abline(h = 0, col = 'grey')
      matlines(data_val, t(fit_plot[, i, ]), col = dat_color)
    }
    abline(h = 0, col = 'grey')
    if(!is.null(dat)) lines(data_val, dat_plot[i, ], lwd = 1.5)
    legend("topright", legend = c(Year[i], ifelse(is.null(N) | is.na(N[i]), "", paste("N =", N[i]))), bty = "n", xjust = 1)

    if(i %% 16 == 1) {
      mtext(data_lab, side = 1, line = 3, outer = TRUE)
      mtext(annual_ylab, side = 2, line = 3.5, outer = TRUE)
    }
  }

  invisible()
}

RCM_get_likelihoods <- function(x, LWT, f_name, s_name) {
  nll_fleet <- x$nll_fleet %>% t()
  nll_fleet[is.na(nll_fleet)] <- 0
  nll_fleet <- cbind(nll_fleet, rowSums(nll_fleet))
  nll_fleet <- rbind(nll_fleet, colSums(nll_fleet))
  colnames(nll_fleet) <- c(f_name, "Sum")
  rownames(nll_fleet) <- c("Catch", "Equilibrium Catch", "CAA", "CAL", "Mean Size", "Sum")

  wt_fleet <- rbind(LWT$Chist, LWT$C_eq, LWT$CAA, LWT$CAL, LWT$MS) %>% structure(dimnames = list(rownames(nll_fleet)[1:5], f_name))
  
  nll_survey <- x$nll_survey %>% t()
  nll_survey[is.na(nll_survey)] <- 0
  nll_survey <- cbind(nll_survey, rowSums(nll_survey))
  nll_survey <- rbind(nll_survey, colSums(nll_survey))
  colnames(nll_survey) <- c(s_name, "Sum")
  rownames(nll_survey) <- c("Index", "CAA", "CAL", "Sum")

  wt_survey <- rbind(LWT$Index, LWT$s_CAA, LWT$s_CAL) %>% structure(dimnames = list(rownames(nll_survey)[1:3], s_name))

  tot <- c(x$nll, x$nll_log_rec_dev, nll_fleet[6, length(f_name) + 1], nll_survey[4, length(s_name) + 1], x$penalty, x$prior) %>% matrix(ncol = 1)
  dimnames(tot) <- list(c("Total", "Recruitment Deviations", "Fleets", "Surveys", "Penalty (High F)", "Priors"), 
                        "Negative log-likelihood")

  res <- list(tot, nll_fleet, wt_fleet, nll_survey, wt_survey) %>% lapply(FUN = function(xx) xx %>% round(2) %>% as.data.frame())
  return(res)
}

#' @rdname plot.RCModel
#' @export
compare_RCM <- function(..., compare = TRUE, filename = "compare_RCM", dir = tempdir(), Year = NULL,
                        f_name = NULL, s_name = NULL, MSY_ref = c(0.5, 1), bubble_adj = 10, scenario = list(), title = NULL,
                        open_file = TRUE, quiet = TRUE, render_args) {

  dots <- list(...)
  test <- vapply(dots, inherits, character(1), what = "RCModel") %>% all()
  if(!test) stop("Objects provided are not of class RCModel.", call. = FALSE)
  
  # Update scenario
  if(is.null(scenario$names)) scenario$names <- paste("Scenario", 1:length(dots))
  if(is.null(scenario$col)) {
    scenario$col <- gplots::rich.colors(length(dots))
  }
  scenario$col2 <- scenario$col

  if(is.null(scenario$lwd)) scenario$lwd <- 1
  if(is.null(scenario$lty)) scenario$lty <- 1

  ####### Function arguments for rmarkdown::render
  filename_rmd <- paste0(filename, ".Rmd")

  if(missing(render_args)) render_args <- list()
  render_args$input <- file.path(dir, filename_rmd)
  if(is.null(render_args$output_format)) {
    render_args$output_format <- "html_document"
  }
  if(is.null(render_args$output_options)) {
    if(render_args$output_format == "html_document") {
      render_args$output_options <- list(df_print = "paged")
    } else {
      render_args$output_options <- list(toc = TRUE, df_print = "kable")
    }
  }
  if(is.null(render_args$quiet)) render_args$quiet <- quiet

  ####### Assign variables
  x <- dots[[1]] # Dummy variable
  report_list <- lapply(dots, function(xx) if(length(xx@mean_fit) > 0) return(xx@mean_fit$report) else stop("Error in RCM objects."))

  nsim <- length(report_list)
  data <- dots[[1]]@data

  max_age <- dots[[1]]@OM@maxage
  age <- 0:max_age
  nyears <- dots[[1]]@OM@nyears
  if(is.null(Year)) Year <- (dots[[1]]@OM@CurrentYr - nyears + 1):dots[[1]]@OM@CurrentYr
  Yearplusone <- c(Year, max(Year) + 1)

  nfleet <- vapply(dots, function(xx) xx@data$nfleet, numeric(1)) %>% unique()
  nsurvey <- vapply(dots, function(xx) xx@data$nsurvey, numeric(1)) %>% unique()
  length_bin <- dots[[1]]@data$length_bin

  # Backwards compatibility
  if(is.null(data$nsel_block)) {
    data$nsel_block <- data$nfleet
    data$sel_block <- matrix(1:data$nfleet, nyears, data$nfleet, byrow = TRUE)
  }

  if(is.null(f_name)) f_name <- paste("Fleet", 1:nfleet)
  if(is.null(s_name)) s_name <- paste("Survey", 1:nsurvey)

  ####### Document header
  if(is.null(title)) title <- "Comparisons of OM conditioning"
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

  if(any(data$Index > 0, na.rm = TRUE)) {
    survey_output <- lapply(1:nsurvey, rmd_RCM_survey_output, s_name = s_name)
  } else survey_output <- NULL

  #### MSY comparisons
  if(compare) {
    message("Running runMSE() to get MSY reference points...")
    if(snowfall::sfIsRunning()) {
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
                 "if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col2, lty = scenario$lty)",
                 "```\n")

    ref_pt_fn <- function(xx) c(mean(xx@Ref$ReferencePoints$FMSY), mean(xx@Ref$ReferencePoints$MSY), mean(xx@Ref$ReferencePoints$SSBMSY_SSB0))
    ref_pt <- do.call(cbind, lapply(Hist, ref_pt_fn)) %>%
      structure(dimnames = list(c("FMSY", "MSY", "Spawning depletion at MSY"), scenario$names)) %>% as.data.frame()

    if(render_args$output_format == "html_document") {
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
  nll <- Map(RCM_get_likelihoods, x = report_list, LWT = lapply(dots, function(xx) xx@data$LWT),
             MoreArgs = list(f_name = f_name, s_name = s_name))

  summary_nll <- vapply(nll, function(xx) xx[[1]] %>% as.matrix(), numeric(4)) %>%
    structure(dimnames = list(rownames(nll[[1]][[1]]), scenario$names)) %>% as.data.frame()

  if(render_args$output_format == "html_document") {

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
  if(is.list(rmd)) rmd <- do.call(c, rmd)

  # Generate markdown report
  if(!dir.exists(dir)) {
    message("Creating directory: \n", dir)
    dir.create(dir)
  }
  write(rmd, file = file.path(dir, filename_rmd))
  message("Generated markdown file: ", file.path(dir, filename_rmd))

  # Rendering markdown file
  message("Rendering markdown file...")
  output_filename <- do.call(rmarkdown::render, render_args)
  message("Rendered file: ", output_filename)

  if(open_file) browseURL(output_filename)
  invisible(output_filename)
}

