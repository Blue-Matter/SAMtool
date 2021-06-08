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
           paste0("if(ncol(yy) > 1) legend(\"topleft\", ", legend.lab, ", text.col = ", col, ")"),
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
           paste0("bl <- unique(RCMdata@sel_block[, ", ff, "])"),
           "vul_bb <- list()",
           "bl_col <- gplots::rich.colors(length(bl))",
           "Year_legend <- character(length(bl))",
           "for(bb in 1:length(bl)) {",
           "  vul_bb[[bb]] <- do.call(cbind, lapply(report_list, function(x) x$vul_len[, bl[bb]]))",
           paste0("  Year_legend[bb] <- Year[RCMdata@sel_block[, ", ff, "] == bl[bb]] %>% range() %>% paste(collapse = \"-\")"),
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
           paste0("  vul_bb_age <- do.call(rbind, lapply(report_list, function(x) x$vul[RCMdata@sel_block[, ", ff, "] == bl[bb], , ", ff, "])) %>% t()"),
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
           paste0("if(any(RCMdata@Chist[, ", ff, "] > 0)) {"),
           paste0("  Cpred <- do.call(cbind, lapply(report_list, function(x) x$Cpred[, ", ff, "]))"),
           paste0("  Chist <- RCMdata@Chist[, ", ff, "]"),
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
           paste0("MAobs <- (RCMdata@CAA[, , ", ff, "] %*% age)/rowSums(RCMdata@CAA[, , ", ff, "], na.rm = TRUE)"),
           "ylim <- c(0.9, 1.1) * range(c(MApred, MAobs), na.rm = TRUE)",
           "matplot(Year, MApred, type = \"l\", col = scenario$col, lty = scenario$lty, lwd = scenario$lwd, xlab = \"Year\", ylab = \"Mean age\", ylim = ylim)",
           paste0("if(any(RCMdata@CAA[, , ", ff, "] > 0, na.rm = TRUE)) {"),
           paste0("  lines(Year, MAobs, col = \"black\", typ = \"o\")"),
           "}",
           "if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col, lty = scenario$lty, lwd = scenario$lwd)",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (black) and predicted (red) mean lengths from ", f_name[ff], ".\"}"),
           paste0("MLpred <- do.call(cbind, lapply(report_list, function(x) x$MLpred[, ", ff, "]))"),
           paste0("if(any(RCMdata@CAL[, , ", ff, "] > 0, na.rm = TRUE)) {"),
           paste0("  MLobs <- (RCMdata@CAL[, , ", ff, "] %*% length_bin)/rowSums(RCMdata@CAL[, , ", ff, "], na.rm = TRUE)"),
           paste0("} else if(RCMdata@MS_type == \"length\" && any(RCMdata@MS[, ", ff, "] > 0, na.rm = TRUE)) MLobs <- RCMdata@MS[, ", ff, "] else MLobs <- NA"),
           "if(!all(is.na(MLpred))) {",
           "  ylim <- c(0.9, 1.1) * range(c(MLpred, MLobs), na.rm = TRUE)",
           "  matplot(Year, MLpred, type = \"l\", col = scenario$col, lty = scenario$lty, lwd = scenario$lwd, xlab = \"Year\", ylab = \"Mean length\", ylim = ylim)",
           "  if(!all(is.na(MLobs))) lines(Year, MLobs, col = \"black\", typ = \"o\")",
           "  if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col, lty = scenario$lty, lwd = scenario$lwd)",
           "}",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (black) and predicted (red) mean weights from ", f_name[ff], ".\"}"),
           paste0("if(RCMdata@MS_type == \"weight\" && any(RCMdata@MS[, ", ff, "] > 0, na.rm = TRUE)) {"),
           paste0("  MWobs <- RCMdata@MS[, ", ff, "]"),
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
           paste0("if(any(RCMdata@CAA[, , ", ff, "] > 0, na.rm = TRUE)) {"),
           paste0("if(nsim == 1) CAA_plot <- array(x@CAA[, , , ", ff, "], c(1, nyears, max_age + 1)) else CAA_plot <- x@CAA[, , , ", ff, "]"),
           paste0("plot_composition_RCM(Year, CAA_plot, RCMdata@CAA[, , ", ff, "], ages = age, dat_col = scenario$col)"),
           "}",
           "```\n",
           paste0("```{r, fig.cap = \"Predicted age composition from fleet ", ff, ".\"}"),
           paste0("if(any(RCMdata@CAA[, , ", ff, "] > 0, na.rm = TRUE)) {"),
           paste0("plot_composition_RCM(Year, CAA_plot, ages = age, dat_col = scenario$col)"),
           "}",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (black) and predicted (red) length composition from ", f_name[ff], ".\"}"),
           paste0("if(any(RCMdata@CAL[, , ", ff, "] > 0, na.rm = TRUE)) {"),
           paste0("if(nsim == 1) CAL_plot <- array(x@CAL[, , , ", ff, "], c(1, nyears, length(RCMdata@length_bin))) else CAL_plot <- x@CAL[, , , ", ff, "]"),
           paste0("plot_composition_RCM(Year, CAL_plot, RCMdata@CAL[, , ", ff, "], CAL_bins = RCMdata@length_bin, dat_col = scenario$col)"),
           "}",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Predicted length composition from ", f_name[ff], ".\"}"),
           paste0("if(any(RCMdata@CAL[, , ", ff, "] > 0, na.rm = TRUE)) {"),
           paste0("plot_composition_RCM(Year, CAL_plot, CAL_bins = RCMdata@length_bin, dat_col = scenario$col)"),
           "}",
           "```\n")
  
  c(header, ans)
}

rmd_RCM_index_output <- function(sur, s_name) {
  ans <- c(paste0("### ", s_name[sur], " \n"),
           "",
           paste0("```{r, fig.cap = \"Selectivity of ", s_name[sur], " in last historical year.\"}"),
           "if(!is.null(report_list[[1]]$ivul)) {",
           paste0("ivul_ff_age <- do.call(cbind, lapply(report_list, function(x) x$ivul[nyears, , ", sur, "]))"),
           paste0("matplot(age, ivul_ff_age, type = \"l\", col = scenario$col2, xlab = \"Age\", ylim = c(0, 1), ylab = \"Selectivity of ", s_name[sur], "\")"),
           "abline(h = 0, col = \"grey\")",
           "if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col2, lty = scenario$lty, lwd = scenario$lwd)",
           "}",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (black) and predicted (red) index values for ", s_name[sur], ".\"}"),
           paste0("Ipred <- do.call(cbind, lapply(report_list, function(x) x$Ipred[, ", sur, "]))"),
           paste0("matplot(Year, Ipred, type = \"l\", col = scenario$col, lty = scenario$lty, lwd = scenario$lwd, ylim = c(0, 1.1 * max(c(Ipred, RCMdata@Index[, ", sur, "]), na.rm = TRUE)), xlab = \"Year\", ylab = \"", s_name[sur], "\")"),
           paste0("lines(Year, RCMdata@Index[, ", sur, "], col = \"black\", typ = \"o\")"),
           "abline(h = 0, col = \"grey\")",
           "if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col, lty = scenario$lty, lwd = scenario$lwd)",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (black) and predicted (red) index values for ", s_name[sur], ". Error bars indicate 95% confidence intervals for observed values.\"}"),
           "if(length(RCMdata@I_sd) && any(RCMdata@I_sd > 0, na.rm = TRUE)) {",
           paste0("  II <- RCMdata@Index[, ", sur, "]"),
           "  ind <- seq(min(which(!is.na(II))), max(which(!is.na(II))), 1)",
           paste0("  err <- exp(log(II) + outer(RCMdata@I_sd[, ", sur, "], c(-1.96, 1.96)))"),
           paste0("  matplot(Year[ind], Ipred[ind, ], type = \"l\", col = scenario$col, lty = scenario$lty, lwd = scenario$lwd, ylim = c(0, 1.1 * max(c(Ipred[ind, ], II[ind], err[ind, ]), na.rm = TRUE)), xlab = \"Year\", ylab = \"", s_name[sur], "\")"),
           "  points(Year[ind], II[ind], lwd = 3, pch = 16)",
           "  arrows(Year[ind], y0 = err[ind, 1], y1 = err[ind, 2], length = 0, lwd = 1.5)",
           "  abline(h = 0, col = \"grey\")",
           "  if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col, lty = scenario$lty, lwd = scenario$lwd)",
           "}",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (black) and predicted (red) mean ages from ", s_name[sur], ".\"}"),
           paste0("if(length(RCMdata@IAA) && any(RCMdata@IAA[, , ", sur, "] > 0, na.rm = TRUE)) {"),
           paste0("MApred <- do.call(cbind, lapply(report_list, function(x) x$IAApred[, , ", sur, "] %*% age/rowSums(x$N[1:nyears, ] * x$ivul[1:nyears, , ", sur, "])))"),
           paste0("MAobs <- (RCMdata@IAA[, , ", sur, "] %*% age)/rowSums(RCMdata@IAA[, , ", sur, "], na.rm = TRUE)"),
           "ylim <- c(0.9, 1.1) * range(c(MApred, MAobs), na.rm = TRUE)",
           "matplot(Year, MApred, type = \"l\", col = scenario$col, lty = scenario$lty, lwd = scenario$lwd, xlab = \"Year\", ylab = \"Mean age\", ylim = ylim)",
           paste0("if(any(RCMdata@IAA[, , ", sur, "] > 0, na.rm = TRUE)) {"),
           paste0("  lines(Year, MAobs, col = \"black\", typ = \"o\")"),
           "}",
           "if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col, lty = scenario$lty, lwd = scenario$lwd)",
           "}",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (black) and predicted (red) mean lengths from ", s_name[sur], ".\"}"),
           paste0("if(length(RCMdata@IAL) && any(RCMdata@IAL[, , ", sur, "] > 0, na.rm = TRUE)) {"),
           paste0("MLpred <- do.call(cbind, lapply(report_list, function(x) x$IALpred[, , ", sur, "] %*% length_bin/rowSums(x$N[1:nyears, ] * x$ivul[1:nyears, , ", sur, "])))"),
           paste0("MLobs <- (RCMdata@IAL[, , ", sur, "] %*% length_bin)/rowSums(RCMdata@IAL[, , ", sur, "], na.rm = TRUE)"),
           "ylim <- c(0.9, 1.1) * range(c(MLpred, MLobs), na.rm = TRUE)",
           "matplot(Year, MLpred, type = \"l\", col = scenario$col, lty = scenario$lty, lwd = scenario$lwd, xlab = \"Year\", ylab = \"Mean length\", ylim = ylim)",
           paste0("if(any(RCMdata@IAL[, , ", sur, "] > 0, na.rm = TRUE)) {"),
           paste0("  lines(Year, MLobs, col = \"black\", typ = \"o\")"),
           "}",
           "if(!is.null(scenario$names)) legend(\"topleft\", scenario$names, col = scenario$col, lty = scenario$lty, lwd = scenario$lwd)",
           "}",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (black) and predicted (red) age composition from ", s_name[sur], ".\"}"),
           paste0("if(length(RCMdata@IAA) && any(RCMdata@IAA[, , ", sur, "] > 0, na.rm = TRUE)) {"),
           paste0("pred_sCAA <- lapply(report_list, function(x) x$IAA[,, ", sur, "]) %>% simplify2array() %>% aperm(perm = c(3, 1, 2))"),
           paste0("plot_composition_RCM(Year, pred_sCAA, RCMdata@IAA[, , ", sur, "], ages = age, dat_col = scenario$col)"),
           "}",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (black) and predicted (red) length composition from ", s_name[sur], ".\"}"),
           paste0("if(length(RCMdata@IAL) && any(RCMdata@IAL[, , ", sur, "] > 0, na.rm = TRUE)) {"),
           paste0("pred_sCAL <- lapply(report_list, function(x) x$IAL[,, ", sur, "]) %>% simplify2array() %>% aperm(perm = c(3, 1, 2))"),
           paste0("plot_composition_RCM(Year, pred_sCAL, RCMdata@IAL[, , ", sur, "], CAL_bins = RCMdata@length_bin, dat_col = scenario$col)"),
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
    "  expectedR <- report$Arec * report$E / (1 + report$Brec * report$E)",
    "} else {",
    "  expectedR <- report$Arec * report$E * exp(-report$Brec * report$E)",
    "}",
    "plot_SR(report$E, expectedR, report$R0, report$E0_SR, report$R)",
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

rmd_RCM_Hist_compare <- function() {
  c("## Updated OM {.tabset}\n",
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
    "        ylim = c(0, 1.1 * max(c(Hist_C, RCMdata@Chist, na.rm = TRUE))))",
    "lines(Year, rowSums(RCMdata@Chist, na.rm = TRUE), col = \"black\")",
    "abline(h = 0, col = \"grey\")",
    "legend(\"topleft\", c(\"RCM\", \"OM\"), col = c(\"black\", \"red\"), pch = c(NA_integer_, 16), lwd = c(1, 1), bty = \"n\")",
    "```\n",
    "",
    "```{r, fig.cap = \"Difference in annual catch (relative to observed), calculated as $C^{OM}_y/C^{obs}_y - 1$. Positive values indicate higher catch in the OM. Catch in the OM is the total removals (both landings and discards).\"}",
    "if(any(RCMdata@Chist > 0, na.rm = TRUE)) {",
    "Catch_difference <- t(Hist_C)/rowSums(RCMdata@Chist, na.rm = TRUE) - 1",
    "Catch_difference[is.infinite(Catch_difference)] <- 0",
    "matplot(Year, Catch_difference, typ = \"n\", xlab = \"Year\", ylab = \"Difference in relative catch\")",
    "abline(h = 0, col = \"grey\")",
    "matlines(Year, Catch_difference, col = \"black\")",
    "}",
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
  if(inherits(x$nll_fleet, "array")) {
    nll_fleet <- apply(x$nll_fleet, 2:3, sum) %>% t()
  } else {
    nll_fleet <- x$nll_fleet %>% t()
  }
  nll_fleet[is.na(nll_fleet)] <- 0
  nll_fleet <- cbind(nll_fleet, rowSums(nll_fleet))
  nll_fleet <- rbind(nll_fleet, colSums(nll_fleet))
  colnames(nll_fleet) <- c(f_name, "Sum")
  rownames(nll_fleet) <- c("Catch", "Equilibrium Catch", "CAA", "CAL", "Mean Size", "Sum")
  
  wt_fleet <- rbind(LWT$Chist, LWT$C_eq, LWT$CAA, LWT$CAL, LWT$MS) %>% structure(dimnames = list(rownames(nll_fleet)[1:5], f_name))
  
  if(inherits(x$nll_index, "array")) {
    nll_index <- apply(x$nll_index, 2:3, sum) %>% t()
  } else {
    nll_index <- x$nll_index %>% t()
  }
  nll_index[is.na(nll_index)] <- 0
  nll_index <- cbind(nll_index, rowSums(nll_index))
  nll_index <- rbind(nll_index, colSums(nll_index))
  colnames(nll_index) <- c(s_name, "Sum")
  rownames(nll_index) <- c("Index", "CAA", "CAL", "Sum")
  
  wt_index <- rbind(LWT$Index, LWT$IAA, LWT$IAL) %>% structure(dimnames = list(rownames(nll_index)[1:3], s_name))
  
  tot <- c(x$nll, x$nll_log_rec_dev, nll_fleet[6, length(f_name) + 1], nll_index[4, length(s_name) + 1], x$penalty, x$prior) %>% matrix(ncol = 1)
  dimnames(tot) <- list(c("Total", "Recruitment Deviations", "Fleets", "Indices", "Penalty (High F)", "Priors"), 
                        "Negative log-likelihood")
  
  res <- list(tot, nll_fleet, wt_fleet, nll_index, wt_index) %>% lapply(FUN = function(xx) xx %>% round(2) %>% as.data.frame())
  return(res)
}


rmd_RCM_likelihood_gradients <- function(f_name, s_name, do_index) {
  header <- c("```{r}",
              "obj <- x@mean_fit$obj",
              "new_dat <- structure(obj$env$data, check.passed = NULL)",
              "new_dat$nll_gr <- 1L",
              "new_par <- as.list(SD, \"Estimate\") %>% structure(what = NULL)",
              "",
              "obj2 <- MakeADFun(data = new_dat, parameters = new_par, map = obj$env$map, random = obj$env$random, ADreport = TRUE,",
              "                  DLL = obj$env$DLL, silent = obj$env$silent)",
              "gr <- obj2$gr() %>% structure(dimnames = list(rownames = names(obj2$fn()), colnames = names(SD$par.fixed)))",
              "par_names <- data.frame(par_type = colnames(gr), par = make_unique_names(names(SD$par.fixed)))",
              "```\n\n")
  
  if(!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("reshape2", quietly = TRUE)) {
    body <- c("#### NULL \n\nInstall the ggplot2 and reshape2 packages to plot gradients.\n\n")
  } else {
    
    fleet_lapply_fn <- function(ff) {
      c(paste("####", f_name[ff], "\n"),
        "```{r, fig.cap = \"Likelihood gradients (annual values by data type in columns) with respect to model parameters (rows).\"}",
        paste0("gr_plot <- dplyr::filter(gr_fleet, Fleet == \"", f_name[ff], "\")"),
        "if(nrow(gr_plot)) {",
        "  ggplot2::ggplot(gr_plot, ggplot2::aes(Year, Gradient, group = par, colour = par)) + ggplot2::facet_grid(par_type ~ data_type, scales = \"free_y\") +",
        paste0("  ggplot2::geom_hline(yintercept = 0, linetype = 3) + ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::theme(legend.position = \"none\") + ggplot2::ggtitle(\"", f_name[ff], "\")"),
        "}",
        "```\n\n")
    }
    index_lapply_fn <- function(sur) {
      c(paste("####", s_name[sur], "\n"),
        "```{r, fig.cap = \"Likelihood gradients (annual values by data type in columns) with respect to model parameters (rows).\"}",
        paste0("gr_plot <- dplyr::filter(gr_index, Index == \"", s_name[sur], "\")"),
        "if(nrow(gr_plot)) {",
        "  ggplot2::ggplot(gr_plot, ggplot2::aes(Year, Gradient, group = par, colour = par)) + ggplot2::facet_grid(par_type ~ data_type, scales = \"free_y\") +",
        paste0("  ggplot2::geom_hline(yintercept = 0, linetype = 3) + ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::theme(legend.position = \"none\") + ggplot2::ggtitle(\"", s_name[sur], "\")"),
        "}",
        "```\n\n")
    }
    
    f_plots <- lapply(1:length(f_name), fleet_lapply_fn) %>% unlist()
    if(do_index) {
      s_plots <- lapply(1:length(s_name), index_lapply_fn) %>% unlist()
    } else {
      s_plots <- NULL
    }
    body <- c("```{r}",
              "gr_fleet <- gr[rownames(gr) == \"nll_fleet\", ] %>% array(dim = dim(report$nll_fleet) %>% c(length(SD$par.fixed))) %>%",
              "  structure(dimnames = list(Year = Year, Fleet = f_name, data_type = c(\"Catch\", \"Equilibrium Catch\", \"CAA\", \"CAL\", \"Mean size\"),",
              "                            par = par_names$par)) %>%", 
              "  reshape2::melt(value.name = \"Gradient\") %>% dplyr::left_join(par_names, by = \"par\") %>%",
              "  dplyr::group_by(data_type) %>% dplyr::filter(any(Gradient != 0))",
              "",
              "gr_index <- gr[rownames(gr) == \"nll_index\", ] %>% array(dim = dim(report$nll_index) %>% c(length(SD$par.fixed))) %>%",
              "  structure(dimnames = list(Year = Year, Index = s_name, data_type = c(\"Index\", \"CAA\", \"CAL\"),",
              "                            par = par_names$par)) %>%", 
              "  reshape2::melt(value.name = \"Gradient\") %>% dplyr::left_join(par_names, by = \"par\") %>%",
              "  dplyr::group_by(data_type) %>% dplyr::filter(any(Gradient != 0))",
              "```\n\n", f_plots, s_plots)
  }
  
  if(requireNamespace("caret", quietly = TRUE)) {
    jac <- c("#### Linear combos\n",
             "",
             "```{r}",
             "gr_combo <- gr[rownames(gr) %in% c(\"nll_fleet\", \"nll_index\"), ] %>% caret::findLinearCombos()",
             "if(is.null(gr_combo$remove)) {",
             "  print(\"Jacobian matrix is of full rank, according to caret::findLinearCombos().\")",
             "} else {",
             "  combo_report <- data.frame(`Parameter Number` = gr_combo$remove, `Parameter Name` = par_names$par[gr_combo$remove])",
             "  print(\"Jacobian matrix is not of full rank, according to caret::findLinearCombos(). Reduce number of model parameters?\")",
             "  print(\"See table below:\")",
             "}",
             "```\n\n",
             "```{r} \nif(!is.null(gr_combo$remove)) combo_report\n```\n\n",
             "Output of caret::findLinearCombos():\n",
             "```{r}\n gr_combo\n```\n\n")
  } else {
    jac <-  c("#### Linear combos\n",
              "",
              "Install the caret package to evaluate if Jacobian matrix is of full rank.\n\n")
  }
  c(header, body, jac)
}

make_unique_names <- function(par_names) {
  unique_names <- par_names %>% unique()
  par_new <- lapply(unique_names, function(x) {
    ind <- par_names == x
    if(sum(ind) > 1) {
      paste0(x, "_", 1:sum(ind))
    } else x
  })
  do.call(c, par_new)
}
