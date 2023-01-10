
summary_SSS <- function(Assessment) {
  output <- summary_SCA(Assessment)
  output$model <- "Simple Stock Synthesis (SSS)"
  return(output)
}

rmd_SSS <- function(Assessment, ...) {
  ss <- rmd_summary("Simple Stock Synthesis (SSS)")

  # Life History
  LH_section <- c(rmd_LAA(header = "## Life History\n"), rmd_WAA(), rmd_LW(),
                  rmd_mat(fig.cap = "Maturity at age. Length-based maturity parameters were converted to the corresponding ages."),
                  rmd_sel(fig.cap = "Selectivity at age."))

  # Data section
  data_section <- rmd_data_timeseries("Catch", header = "## Data\n")

  # Assessment
  #### Pars and Fit
  assess_fit <- c(rmd_R0(header = "## Assessment {.tabset}\n### Estimates and Model Fit\n"),
                  rmd_assess_fit_series(nsets = 1))

  #### Time Series
  ts_output <- c(rmd_U(header = "### Time Series Output\n"), rmd_U_UMSY(), rmd_SSB(), 
                 rmd_dynamic_SSB0("TMB_report$dynamic_SSB0"), rmd_SSB_SSBMSY(), rmd_SSB_SSB0(), 
                 rmd_Kobe("SSB_SSBMSY", "U_UMSY", xlab = "expression(SSB/SSB[MSY])", ylab = "expression(U/U[MSY])"), rmd_R(),
                 rmd_N(), rmd_N_at_age())

  # Productivity
  SR_calc <- c("SSB_SR <- SSB",
               "Rest <- R_SR <- R[as.numeric(names(R)) >= info$Year[1]]")
  productivity <- c(rmd_SR(header = "### Productivity\n\n\n", SR_calc = SR_calc),
                    rmd_SR(fig.cap = "Stock-recruit relationship (trajectory plot).", trajectory = TRUE),
                    rmd_yield_U("SCA_Pope"), rmd_yield_depletion("SCA_Pope"), rmd_sp(), rmd_SPR(), rmd_YPR())

  return(c(ss, LH_section, data_section, assess_fit, ts_output, productivity))
}


profile_likelihood_SSS <- function(Assessment, ...) {
  dots <- list(...)
  if (!"R0" %in% names(dots)) stop("Sequence of R0 is needed for profile.")

  nll <- vapply(log(dots$R0 * Assessment@obj$env$data$rescale), function(xx) Assessment@obj$fn(xx), numeric(1))
  output <- new("prof", Model = Assessment@Model, Name = Assessment@Name, Par = "R0", MLE = Assessment@R0,
                grid = data.frame(R0 = dots$R0, nll = nll))
  return(output)
}


retrospective_SSS <- function(Assessment, nyr) retrospective_SCA(Assessment, nyr)
