
summary_SCA_RWM <- function(Assessment) {
  assign_Assessment_slots(Assessment)

  if(conv) current_status <- c(F_FMSY[length(F_FMSY)], SSB_SSBMSY[length(SSB_SSBMSY)], SSB_SSB0[length(SSB_SSB0)])
  else current_status <- c(NA, NA, SSB_SSB0[length(SSB_SSB0)])
  current_status <- data.frame(Value = current_status)
  rownames(current_status) <- c("F/FMSY", "SSB/SSBMSY", "SSB/SSB0")

  Value <- c(h, info$data$n_age - 1, info$LH$Linf, info$LH$K, info$LH$t0,
             info$LH$a * info$LH$Linf ^ info$LH$b, info$LH$A50, info$LH$A95)
  Description = c("Stock-recruit steepness", "Maximum age (plus-group)", "Asymptotic length", "Growth coefficient",
                  "Age at length-zero", "Asymptotic weight", "Age of 50% maturity", "Age of 95% maturity")
  rownam <- c("h", "maxage", "Linf", "K", "t0", "Winf", "A50", "A95")
  input_parameters <- data.frame(Value = Value, Description = Description, stringsAsFactors = FALSE)
  rownames(input_parameters) <- rownam
  if(!"transformed_h" %in% names(obj$env$map)) input_parameters <- input_parameters[-1, ]

  if(conv) Value <- c(VB0, SSB0[length(SSB0)], MSY, FMSY, VBMSY, SSBMSY, SSBMSY/SSB0[length(SSB0)])
  else Value <- rep(NA, 7)

  Description <- c("Unfished vulnerable biomass",
                   "Unfished spawning stock biomass (SSB)", "Maximum sustainable yield (MSY)", "Fishing mortality at MSY",
                   "Vulnerable biomass at MSY", "SSB at MSY", "Spawning depletion at MSY")
  derived <- data.frame(Value = Value, Description = Description, stringsAsFactors = FALSE)
  rownames(derived) <- c("VB0", "SSB0", "MSY", "FMSY", "VBMSY", "SSBMSY", "SSBMSY/SSB0")

  model_estimates <- sdreport_int(SD)
  if(!is.character(model_estimates)) {
    rownames(model_estimates)[rownames(model_estimates) == "log_F_"] <- paste0("log_F_dev_", names(FMort))
    rownames(model_estimates)[rownames(model_estimates) == "log_rec_dev"] <- paste0("log_rec_dev_", names(FMort)[as.logical(obj$env$data$est_rec_dev)])
    rownames(model_estimates)[rownames(model_estimates) == "logit_M"] <- paste0("logit_M_", names(FMort))
  }

  output <- list(model = "Statistical Catch-at-Age with Random Walk in M",
                 current_status = current_status, input_parameters = input_parameters,
                 derived_quantities = derived, model_estimates = model_estimates,
                 log_likelihood = matrix(NLL, ncol = 1, dimnames = list(names(NLL), "Neg.LL")))
  return(output)
}

rmd_SCA_RWM <- function(Assessment, ...) {
  ss <- rmd_summary("Statistical Catch-at-Age with Random Walk in M")

  # Life History
  age <- 0:(Assessment@info$data$n_age - 1)
  LH_section <- c(rmd_LAA(age, Assessment@info$LH$LAA, header = "## Life History\n"), rmd_WAA(age, Assessment@info$LH$WAA),
                  rmd_LW(Assessment@info$LH$LAA, Assessment@info$LH$WAA),
                  rmd_mat(age, Assessment@info$data$mat,
                          fig.cap = "Maturity at age. Length-based maturity parameters were converted to the corresponding ages."))
  # Data section
  data_section <- c(rmd_data_timeseries("Catch", header = "## Data\n"),
                    rmd_data_timeseries("Index", is_matrix = is.matrix(Assessment@Obs_Index), nsets = ncol(Assessment@Obs_Index)),
                    rmd_data_age_comps("bubble"), rmd_data_age_comps("annual"))

  # Assessment
  #### Pars and Fit
  lead_par <- c(rmd_R0(header = "## Assessment {.tabset}\n### Estimates and Model Fit\n"), rmd_h())

  assess_fit <- c(lead_par,
                  rmd_sel(age, Assessment@Selectivity[nrow(Assessment@Selectivity), ], fig.cap = "Estimated selectivity at age."),
                  rmd_assess_fit("Catch", "catch"), rmd_assess_resid("Catch"), rmd_assess_qq("Catch", "catch"),
                  rmd_assess_fit_series(nsets = ncol(Assessment@Index)),
                  rmd_fit_age_comps("bubble"), rmd_fit_age_comps("annual"),
                  rmd_residual("Dev", fig.cap = "Time series of recruitment deviations.", label = Assessment@Dev_type,
                               blue = any(as.numeric(names(Assessment@Dev)) < Assessment@info$Year[1])),
                  rmd_residual("Dev", "SE_Dev", fig.cap = "Time series of recruitment deviations with 95% confidence intervals.",
                               label = Assessment@Dev_type, conv_check = TRUE, blue = any(as.numeric(names(Assessment@Dev)) < Assessment@info$Year[1])),
                  rmd_M_rw())
  
  #### Time Series
  ts_output <- c(rmd_F(header = "### Time Series Output\n"), rmd_F_FMSY(),  rmd_M(), rmd_SSB(),
                 rmd_dynamic_SSB0("TMB_report$dynamic_SSB0"), rmd_SSB_SSBMSY(),
                 rmd_SSB_SSB0(), rmd_Kobe("SSB_SSBMSY", xlab = "expression(SSB/SSB[MSY])"), rmd_R(),
                 rmd_N(), rmd_N_at_age(), rmd_C_at_age(), rmd_C_mean_age())

  # Productivity
  Arec <- Assessment@TMB_report$Arec
  Brec <- Assessment@TMB_report$Brec
  SSB <- Assessment@SSB

  SR <- Assessment@info$data$SR_type
  if(SR == "BH") {
    expectedR <- Arec * SSB / (1 + Brec * SSB)
  } else {
    expectedR <- Arec * SSB * exp(-Brec * SSB)
  }
  estR <- Assessment@R[as.numeric(names(Assessment@R)) >= Assessment@info$Year[1]]

  productivity <- c(rmd_SR(SSB, expectedR, estR, header = "### Productivity\n\n\n", unfished = FALSE),
                    rmd_SR(SSB, expectedR, estR, fig.cap = "Stock-recruit relationship (trajectory plot).", 
                           trajectory = TRUE, unfished = FALSE),
                    rmd_SPR(), rmd_YPR())

  return(c(ss, LH_section, data_section, assess_fit, ts_output, productivity))
}

retrospective_SCA_RWM <- function(Assessment, nyr, ...) retrospective_SCA(Assessment, nyr)

profile_likelihood_SCA_RWM <- function(Assessment, ...) profile_likelihood_SCA(Assessment, ...)

