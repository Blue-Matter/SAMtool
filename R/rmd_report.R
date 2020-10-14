rmd_head <- function(name, assessment = TRUE) {
   if(assessment) name <- paste("Assessment summary for", name)
   nametag <- c("---", paste0("title: \"", name, "\""))
   ans <- c(nametag,
            "subtitle: Tables and Figures",
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
            "  fig.width = 6, fig.height = 4.5, out.width = \"650px\", comment = \"#>\")",
            "```\n")
   return(ans)
}

rmd_summary <- function(modname) {
  nametag <- paste0("# ", modname, "{.tabset}\n")
  ans <- c(nametag,
           "## Summary Tables {.tabset}\n",
           "```{r}",
           "  sx <- summary(Assessment)[-1]",
           "  for(i in 1:length(sx)) {",
           "    dat <- as.data.frame(sx[[i]])",
           "    for(j in 1:ncol(dat)) if(nrow(dat) > 0 && is.numeric(dat[, j])) dat[, j] <- ifelse(dat[, j] > 1e3, round(dat[, j], 0), signif(dat[, j], 3))",
           "    sx[[i]] <- dat",
           "  }",
           "```\n\n",
           "### Current Status",
           "`r sx[[1]]`\n",
           "### Input Parameters",
           "`r sx[[2]]`\n",
           "### Derived Quantities",
           "`r sx[[3]]`\n",
           "### Model Estimates",
           "`r sx[[4]]`\n",
           "### Likelihoods",
           "`r sx[[5]]`\n")
  return(ans)
}



### Life history
vector2char <- function(x) paste0("c(", paste0(x, collapse = ", "), ")")

rmd_at_age <- function(age, y_var, fig.cap, label, header = NULL) {
  age <- vector2char(age)
  y_var <- vector2char(y_var)

  ans <- c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
           paste0("plot_generic_at_age(", age, ", ", y_var, ", label = \"", label, "\")"),
           " ```\n")
  if(!is.null(header)) ans <- c(header, ans)
  return(ans)
}

rmd_LAA <- function(age, LAA, header = NULL) {
  rmd_at_age(age, LAA, fig.cap = "Mean length-at-age from Data object.", label = "Mean Length-at-age", header = header)
}

rmd_WAA <- function(age, WAA, header = NULL) {
  rmd_at_age(age, WAA, fig.cap = "Mean weight-at-age from Data object.", label = "Mean Weight-at-age", header = header)
}

rmd_LW <- function(LAA, WAA) {
  LAA_char <- paste0("c(", paste0(LAA, collapse = ", "), ")")
  WAA_char <- paste0("c(", paste0(WAA, collapse = ", "), ")")

  return(c("```{r, fig.cap=\"Length-weight relationship.\"}",
           paste0("plot(", LAA_char, ", ", WAA_char, ", typ = \"o\", xlab = \"Length\", ylab = \"Weight\")"),
           "abline(h = 0, col = \"grey\")",
           "```\n"))
}

rmd_mat <- function(age, mat, fig.cap) rmd_at_age(age, mat, fig.cap, "Maturity")




#### Data
rmd_data_timeseries <- function(type, header = NULL, is_matrix = FALSE, nsets = 1) {
  if(is_matrix) {
    lapply_fn <- function(x) {
      c(paste0("```{r, fig.cap=\"", type, " ", x, " time series.\"}"),
        paste0("plot_timeseries(as.numeric(rownames(Obs_", type, ")), Obs_", type, "[, ", x, "], label = paste(\"", type, "\", ", x, "))"),
        "```\n")
    }
    ans <- do.call(c, lapply(1:nsets, lapply_fn))
  } else {
    ans <- c(paste0("```{r, fig.cap=\"", type, " time series.\"}"),
             paste0("plot_timeseries(as.numeric(names(Obs_", type, ")), Obs_", type, ", label = \"", type, "\")"),
             "```\n")
  }
  if(!is.null(header)) ans <- c(header, ans)
  return(ans)
}


rmd_data_age_comps <- function(type = c("bubble", "annual"), ages = "NULL", annual_yscale = "\"proportions\"",
                               annual_ylab = "\"Frequency\"")  {
  type <- match.arg(type)
  if(type == "bubble") {
    arg <- "\"bubble_data\""
    fig.cap = "Age composition bubble plot."
  } else {
    arg <- "\"annual\""
    fig.cap <- "Annual age compositions."
  }
  c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
    "ind_valid <- rowSums(Obs_C_at_age, na.rm = TRUE) > 0",
    paste0("plot_composition(info$Year[ind_valid], Obs_C_at_age[ind_valid, ], ages = ", ages, ", plot_type = ", arg, ","),
    paste0("                 annual_yscale = ", annual_yscale, ", annual_ylab = ", annual_ylab, ")"),
    "```\n")
}

#### Assessment parameters
rmd_R0 <- function(header = NULL) {
  fig.cap <- "Estimate of R0, distribution based on normal approximation of estimated covariance matrix."
  ans <- c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
           "if(conv) {",
           "  ind <- names(SD$par.fixed) == \"R0x\"",
           "  mu <- SD$par.fixed[ind] - log(obj$env$data$rescale)",
           "  sig <- sqrt(diag(SD$cov.fixed)[ind])",
           "  plot_lognormalvar(mu, sig, label = expression(Unfished~~recruitment~~(R[0])), logtransform = TRUE)",
           "}",
           "```\n")
  if(!is.null(header)) ans <- c(header, ans)
  return(ans)
}

rmd_meanR <- function(header = NULL) {
  fig.cap <- "Estimate of mean recruitment, distribution based on normal approximation of estimated covariance matrix."
  ans <- c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
           "if(conv) {",
           "  ind <- names(SD$par.fixed) == \"meanRx\"",
           "  mu <- SD$par.fixed[ind] - log(obj$env$data$rescale)",
           "  sig <- sqrt(diag(SD$cov.fixed)[ind])",
           "  plot_lognormalvar(mu, sig, label = \"Mean recruitment\", logtransform = TRUE)",
           "}",
           "```\n")
  if(!is.null(header)) ans <- c(header, ans)
  return(ans)
}


rmd_h <- function() {
  fig.cap <- "Estimate of steepness, distribution based on normal approximation of estimated covariance matrix."
  c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
    "if(conv && !\"transformed_h\" %in% names(obj$env$map)) {",
    "  ind <- names(SD$par.fixed) == \"transformed_h\"",
    "  plot_steepness(SD$par.fixed[ind], sqrt(diag(SD$cov.fixed)[ind]), is_transform = TRUE, SR = info$data$SR_type)",
    "}",
    "```\n")
}

rmd_FMSY <- function(header = NULL) {
  fig.cap <- "Estimate of FMSY, distribution based on normal approximation of estimated covariance matrix."
  ans <- c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
           "if(conv) {",
           "  Fmsy.ind <- names(SD$par.fixed) == \"log_FMSY\"",
           "  log.Fmsy <- SD$par.fixed[Fmsy.ind]",
           "  log.Fmsy.sd <- sqrt(diag(SD$cov.fixed)[Fmsy.ind])",
           "  plot_lognormalvar(log.Fmsy, log.Fmsy.sd, label = expression(hat(F)[MSY]), logtransform = TRUE)",
           "}",
           "```\n")
  if(!is.null(header)) ans <- c(header, ans)
  return(ans)
}

rmd_MSY <- function(par = "MSYx") {
  fig.cap <- "Estimate of MSY, distribution based on normal approximation of estimated covariance matrix."
  c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
    "if(conv) {",
    paste0("  msy.ind <- names(SD$par.fixed) == \"", par, "\""),
    "  log_rescale <- ifelse(is.null(obj$env$data$rescale), 0, log(obj$env$data$rescale))",
    "  log.msy <- SD$par.fixed[msy.ind] - log_rescale",
    "  log.msy.sd <- sqrt(diag(SD$cov.fixed)[msy.ind])",
    "  plot_lognormalvar(log.msy, log.msy.sd, logtransform = TRUE, label = expression(widehat(MSY)))",
    "}",
    "```\n")
}

rmd_F_FMSY_terminal <- function() {
  fig.cap <- paste0("Estimate of F/FMSY in terminal year, distribution based on normal approximation of estimated covariance matrix.")
  c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
    "if(conv) {",
    "  Fy <- names(F_FMSY)[length(F_FMSY)]",
    "  plot_normalvar(F_FMSY[length(F_FMSY)], SE_F_FMSY_final, label = bquote(F[.(Fy)]/F[MSY]))",
    "}",
    "```\n")
}

rmd_B_BMSY_terminal <- function() {
  fig.cap <- paste0("Estimate of B/BMSY in terminal year, distribution based on normal approximation of estimated covariance matrix.")
  c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
    "if(conv) {",
    "  By <- names(B_BMSY)[length(B_BMSY)]",
    "  plot_normalvar(B_BMSY[length(B_BMSY)], SE_B_BMSY_final, label = bquote(B[.(By)]/B[MSY]))",
    "}",
    "```\n")
}

rmd_B_B0_terminal <- function() {
  fig.cap <- paste0("Estimate of biomass depletion in terminal year, distribution based on normal approximation of estimated covariance matrix.")
  c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
    "if(conv) {",
    "  By <- names(B_B0)[length(B_B0)]",
    "  plot_normalvar(B_B0[length(B_B0)], SE_B_B0_final, label = bquote(B[.(By)]/B[0]))",
    "}",
    "```\n")
}


rmd_sel <- function(age, sel, fig.cap) {
  age <- vector2char(age)
  sel <- vector2char(sel)
  c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
    paste0("plot_ogive(", age, ", ", sel, ")"),
    "```\n")
}

#' @importFrom graphics persp
rmd_sel_persp <- function(age, sel = "Selectivity", fig.cap = "Perspective plot of selectivity.") {
  age <- vector2char(age)
  c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
    paste0("persp(info$Year, ", age, ", ", sel, ", expand = 0.35, ticktype = \"detailed\", phi = 25,"),
    paste0("      theta = 45, xlab = \"Year\", ylab = \"Age\", zlab = \"Selectivity\")"),
    "```\n")
}

rmd_sel_annual <- function(age, sel = "Selectivity", fig.cap = "Annual selectivity.") {
  age <- vector2char(age)
  c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
    paste0("plot_composition(info$Year, Selectivity, plot_type = \"annual\", ages = ", age, ", annual_yscale = \"raw\", annual_ylab = \"Selectivity\")"),
    "```\n")
}


#### Assessment time series
rmd_assess_fit <- function(par, fig.cap, label = par, match = FALSE) {
  fig.cap2 <- paste0("Observed (black) and predicted (red) ", fig.cap, ".")
  if(match) fig.cap2 <- paste(fig.cap2, "Predicted", fig.cap, "should match observed in this model.")

  c(paste0("```{r, fig.cap=\"", fig.cap2, "\"}"),
    paste0("plot_timeseries(as.numeric(names(", par, ")), Obs_", par, ", ", par, ", label = \"", label, "\")"),
    "```\n")

}

rmd_assess_resid <- function(par, fig.cap = paste(par, "residuals in log-space."), label = paste0("log(", par, ") Residual")) {
  c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
    paste0("plot_residuals(as.numeric(names(", par, ")), log(Obs_", par, "/", par, "), label = \"", label, "\")"),
    "```\n")
}

#' @importFrom stats qqline qqnorm
rmd_assess_qq <- function(par, fig.cap) {
  fig.cap2 <- paste("QQ-plot of", fig.cap, "residuals in log-space.")
  c(paste0("```{r, fig.cap=\"", fig.cap2, "\"}"),
    paste0("qqnorm(log(Obs_", par, "/", par, "), main = \"\")"),
    paste0("qqline(log(Obs_", par, "/", par, "))"),
    "```\n")
}

rmd_assess_fit_series <- function(par = "Index", fig.cap = "index", label = par, nsets = 1) {

  lapply_fn <- function(x) {
    c(paste0("```{r, fig.cap=\"Observed (black) and predicted (red) values of ", fig.cap, " ", x, ".\"}"),
      paste0("plot_timeseries(as.numeric(rownames(", par, ")), Obs_", par, "[, ", x, "], ", par, "[, ", x, "], label = paste(\"", label, "\", ", x, "))"),
      "```\n",
      paste0("```{r, fig.cap=\"", par, " ", x, " residuals in log-space.\"}"),
      paste0("plot_residuals(as.numeric(rownames(", par, ")), log(Obs_", par, "[, ", x, "]/", par, "[, ", x, "]), "),
      paste0("               label = paste(\"", label, "\", ", x, ", \"log-residuals\"))"),
      "```\n",
      paste0("```{r, fig.cap=\"QQ-plot of log-residuals from ", label, " ", x, ".\"}"),
      paste0("qqnorm(log(Obs_", par, "[, ", x, "]/", par, "[, ", x, "]), main = \"\")"),
      paste0("qqline(log(Obs_", par, "[, ", x, "]/", par, "[, ", x, "]))"),
      "```\n")
  }
  do.call(c, lapply(1:nsets, lapply_fn))

}


rmd_assess_timeseries <- function(par, fig.cap, label, header = NULL, conv_check = FALSE, one_line = FALSE) {
  if(conv_check) {
    conv <- "if(conv) {\n  "
    conv2 <- "}"
  } else conv <- conv2 <- ""
  if(one_line) ab_one <- "abline(h = 1, lty = 3)" else ab_one <- ""
  ans <- c(paste0("```{r, fig.cap=\"Time series of ", fig.cap, ".\"}"),
           paste0(conv, "plot_timeseries(as.numeric(names(", par, ")),", par, ", label = ", label, ")"),
           ab_one,
           conv2,
           "```\n")
  if(!is.null(header)) ans <- c(header, ans)
  return(ans)
}

rmd_residual <- function(par, se = "NULL", fig.cap, label, conv_check = FALSE, blue = FALSE) {
  if(conv_check) conv <- "if(conv) " else conv <- ""
  if(blue) {
    blue_arg <- ", res_ind_blue = as.numeric(names(Dev)) < info$Year[1]"
    fig.cap <- paste(fig.cap, "Deviations prior to the first year of the model are in blue.")
  } else {
    blue_arg <- ""
  }

  c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
    paste0(conv, "plot_residuals(as.numeric(names(", par, ")), ", par, " , res_sd = ", se, blue_arg, ", label = \"", label, "\")"),
    "```\n")
}

rmd_fit_age_comps <- function(type = c("bubble", "annual"), ages = "NULL", match = FALSE)  {
  type <- match.arg(type)
  if(type == "bubble") {
    arg <- "\"bubble_residuals\", bubble_adj = 20"
    fig.cap = "Pearson residual bubble plot of age compositions (grey bubbles are negative, white are positive)."
  } else {
    arg <- paste("\"annual\", ages =", ages)

    if(!match) {
      arg <- paste0(arg, ", N = info$data$CAA_n[ind_valid]")
      fig.cap <- "Annual observed (black) and predicted (red) age compositions."
    } else {
      fig.cap <- "Annual observed (black) and predicted (red) age compositions. Predicted should match observed in this model."
    }
  }
  c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
    "ind_valid <- rowSums(Obs_C_at_age, na.rm = TRUE) > 0",
    paste0("plot_composition(info$Year[ind_valid], Obs_C_at_age[ind_valid, ], C_at_age[ind_valid, ], plot_type = ", arg, ")"),
    "```\n")
}

rmd_bubble <- function(year, par, CAL_bins = "NULL", ages = "NULL", fig.cap, bubble_adj = "5") {
  c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
    paste0("plot_composition(", year, ", ", par, ", CAL_bins = ", CAL_bins, ", ages = ", ages, ", plot_type = \"bubble_data\", bubble_adj = ", bubble_adj, ")"),
    "```\n")
}


rmd_F <- function(header = NULL, fig.cap = "fishing mortality") {
  rmd_assess_timeseries("FMort", "fishing mortality", "\"Fishing Mortality (F)\"", header = header)
}

rmd_F_FMSY <- function(conv_check = TRUE) {
  rmd_assess_timeseries("F_FMSY", "F/FMSY", "expression(F/F[MSY])", conv_check = conv_check, one_line = TRUE)
}

rmd_U <- function(header = NULL, fig.cap = "harvest rate") {
  rmd_assess_timeseries("U", fig.cap, "\"Harvest rate (U)\"", header = header)
}

rmd_U_UMSY <- function(conv_check = TRUE, fig.cap = "U/UMSY") {
  rmd_assess_timeseries("U_UMSY", fig.cap, "expression(U/U[MSY])", conv_check = conv_check, one_line = TRUE)
}

rmd_SSB <- function() rmd_assess_timeseries("SSB", "spawning biomass", "\"Spawning biomass\"")

rmd_dynamic_SSB0 <- function() rmd_assess_timeseries("dynamic_SSB0", "dynamic SSB0", "\"Dynamic SSB0\"")

rmd_SSB_SSBMSY <- function(conv_check = TRUE) {
  rmd_assess_timeseries("SSB_SSBMSY", "SSB/SSBMSY", "expression(SSB/SSB[MSY])", conv_check = conv_check, one_line = TRUE)
}

rmd_SSB_SSB0 <- function(conv_check = TRUE) {
  rmd_assess_timeseries("SSB_SSB0", "spawning depletion", "expression(SSB/SSB[0])", conv_check = conv_check)
}

rmd_B <- function() rmd_assess_timeseries("B", "biomass", "\"Biomass\"")

rmd_B_BMSY <- function(conv_check = TRUE) {
  rmd_assess_timeseries("B_BMSY", "B/BMSY", "expression(B/B[MSY])", conv_check = conv_check, one_line = TRUE)
}

rmd_B_B0 <- function(conv_check = TRUE) {
  rmd_assess_timeseries("B_B0", "depletion", "expression(B/B[0])", conv_check = conv_check)
}

rmd_R <- function() rmd_assess_timeseries("R", "recruitment", "\"Recruitment (R)\"")

rmd_N <- function() rmd_assess_timeseries("N", "abundance", "\"Abundance (N)\"")

rmd_N_at_age <- function() rmd_bubble("c(info$Year, max(info$Year)+1)", "N_at_age", fig.cap = "Abundance-at-age bubble plot.")

rmd_C_at_age <- function() rmd_bubble("info$Year", "C_at_age", fig.cap = "Predicted catch-at-age bubble plot.")

rmd_C_mean_age <- function() {
  c(paste0("```{r, fig.cap=\"Observed (black) and predicted (red) mean age of the composition data.\"}"),
    "plot_composition(info$Year, Obs_C_at_age, C_at_age, plot_type = \"mean\")",
    "```\n")
}

rmd_Kobe <- function(Bvar = "B_BMSY", Fvar = "F_FMSY", xlab = "expression(B/B[MSY])", ylab = "expression(F/F[MSY])",
                     conv_check = TRUE) {
  if(conv_check) conv <- "if(conv) " else conv <- ""
  c("```{r, fig.cap=\"Kobe plot trajectory.\"}",
    paste0(conv, "plot_Kobe(", Bvar, ", ", Fvar, ", xlab = ", xlab, ", ylab = ", ylab, ")"),
    "```\n")
}


#### Productivity
rmd_SR <- function(Bvec, expectedR, Rpoints, fig.cap = "Stock-recruit relationship.", trajectory = FALSE, ylab = "Recruitment",
                   conv_check = FALSE, header = NULL) {
  Bvec <- vector2char(Bvec)
  expectedR <- vector2char(expectedR)
  Rpoints <- vector2char(Rpoints)

  if(conv_check) conv <- "if(conv) " else conv <- ""
  ans <- c(paste0("```{r fig.cap=\"", fig.cap, "\"}"),
           paste0(conv, "plot_SR(", Bvec, ", ", expectedR, ", rec_dev = ", Rpoints, ","),
           paste0("R0 = R0, S0 = SSB0, ylab = \"", ylab, "\", trajectory = ", as.character(trajectory), ")"),
           "```\n")
  if(!is.null(header)) ans <- c(header, ans)
  return(ans)
}

rmd_yield_F <- function(model, conv_check = TRUE, header = NULL) {
  if(conv_check) conv <- "if(conv) " else conv <- ""
  if(model == "VPA") extra_code <- "info$data$SR_type <- info$SR; TMB_report$vul <- TMB_report$vul_refpt" else extra_code <- ""

  ans <- c("```{r, fig.cap=\"Yield plot relative to fishing mortality.\"}",
           extra_code,
           paste0(conv, "plot_yield_", model, "(info$data, TMB_report, FMSY, MSY, xaxis = \"F\")"),
           "```\n")
  if(!is.null(header)) ans <- c(header, ans)
  return(ans)
}

rmd_yield_U <- function(model, conv_check = TRUE, header = NULL) {
  if(conv_check) conv <- "if(conv) " else conv <- ""

  ans <- c("```{r, fig.cap=\"Yield plot relative to harvest rate.\"}",
           paste0(conv, "plot_yield_", model, "(info$data, TMB_report, UMSY, MSY, xaxis = \"U\")"),
           "```\n")
  if(!is.null(header)) ans <- c(header, ans)
  return(ans)
}

rmd_yield_depletion <- function(model, conv_check = TRUE) {
  if(conv_check) conv <- "if(conv) " else conv <- ""
  if(model == "SCA_Pope") rate <- "UMSY" else rate <- "FMSY"

  c("```{r, fig.cap=\"Yield plot relative to depletion.\"}",
    paste0(conv, "yield_fn <- plot_yield_", model, "(info$data, TMB_report, ", rate, ", MSY, xaxis = \"Depletion\")"),
    "```\n")
}

rmd_sp <- function(conv_check = TRUE, depletion = TRUE) {
  if(conv_check) conv <- "if(conv) " else conv <- ""
  if(depletion) B0 <- "B0" else B0 <- "NULL"
  c("```{r, fig.cap=\"Comparison of historical surplus production and estimated yield curve.\"}",
    paste0(conv, "plot_surplus_production(B, B0 = ", B0, ", Catch, yield_fn = yield_fn)"),
    "```\n")
}

rmd_retrospective <- function() {
  c("## Retrospective\n",
    "```{r}",
    "as.data.frame(summary(retro))",
    "plot(retro)",
    "```\n")
}

#### Footer
rmd_footer <- function() {
  c("## About\n",
    "This report was generated on: `r Sys.time()`<br />",
    "MSEtool R package version `r packageVersion(\"MSEtool\")`<br />",
    "`r R.version.string`<br />\n")
}


