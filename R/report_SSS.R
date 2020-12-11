
summary_SSS <- function(Assessment) {
  output <- summary_SCA_Pope(Assessment)
  output$model <- "Simple Stock Synthesis (SSS)"
  return(output)
}

rmd_SSS <- function(Assessment, ...) {
  ss <- rmd_summary("Simple Stock Synthesis (SSS)")

  # Life History
  age <- 0:(Assessment@info$data$n_age - 1)
  LH_section <- c(rmd_LAA(age, Assessment@info$LH$LAA, header = "## Life History\n"), rmd_WAA(age, Assessment@info$LH$WAA),
                  rmd_LW(Assessment@info$LH$LAA, Assessment@info$LH$WAA),
                  rmd_mat(age, Assessment@info$data$mat,
                          fig.cap = "Maturity at age. Length-based maturity parameters were converted to the corresponding ages."),
                  rmd_sel(age, Assessment@Selectivity[nrow(Assessment@Selectivity), ], fig.cap = "Selectivity at age."))

  # Data section
  data_section <- rmd_data_timeseries("Catch", header = "## Data\n")

  # Assessment
  #### Pars and Fit
  assess_fit <- c(rmd_R0(header = "## Assessment {.tabset}\n### Estimates and Model Fit\n"),
                  rmd_assess_fit("Index", "index"))

  #### Time Series
  ts_output <- c(rmd_U(header = "### Time Series Output\n"), rmd_U_UMSY(), rmd_SSB(), rmd_SSB_SSBMSY(),
                 rmd_SSB_SSB0(), rmd_Kobe("SSB_SSBMSY", "U_UMSY", xlab = "expression(SSB/SSB[MSY])", ylab = "expression(U/U[MSY])"), rmd_R(),
                 rmd_N(), rmd_N_at_age())

  # Productivity
  SSB <- Assessment@SSB
  expectedR <- estR <- Assessment@R[as.numeric(names(Assessment@R)) >= Assessment@info$Year[1]]

  productivity <- c(rmd_SR(SSB, expectedR, estR, header = "### Productivity\n\n\n"),
                    rmd_SR(SSB, expectedR, estR, fig.cap = "Stock-recruit relationship (trajectory plot).", trajectory = TRUE),
                    rmd_yield_U("SCA_Pope"), rmd_yield_depletion("SCA_Pope"), rmd_sp(), rmd_SPR(), rmd_YPR())

  return(c(ss, LH_section, data_section, assess_fit, ts_output, productivity))
}


profile_likelihood_SSS <- function(Assessment, ...) {
  dots <- list(...)
  if(!"R0" %in% names(dots)) stop("Sequence of R0 is needed for profile.")

  nll <- vapply(log(dots$R0 * Assessment@obj$env$data$rescale), function(xx) Assessment@obj$fn(xx), numeric(1))
  output <- new("prof", Model = Assessment@Model, Name = Assessment@Name, Par = "R0", MLE = Assessment@R0,
                grid = data.frame(R0 = dots$R0, nll = nll))
  return(output)
}


retrospective_SSS <- function(Assessment, nyr) retrospective_SCA_Pope(Assessment, nyr)
