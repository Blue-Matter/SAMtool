
# Write RCM_int() fits back into an om (fitted values, then PopulateOM())

.SRRHSlot <- function(Model) {
  nm <- .SRRModelName(Model)
  if (!is.na(nm) && identical(nm, "Ricker")) "hR" else "h"
}

# Sim x Year apical F for one fleet
.FitApicalF <- function(report, fl, nyears) {
  out <- vapply(report, function(x) x$F[1:nyears, fl], numeric(nyears))
  matrix(t(out), nrow = length(report))
}

# Sim x Age x Year selectivity for one fleet
.FitVulArray <- function(report, fl, nyears, nAge) {
  vul_list <- lapply(report, function(x) t(matrix(x$vul[1:nyears, , fl], nyears, nAge)))
  arr <- array(unlist(vul_list), c(nAge, nyears, length(report)))
  aperm(arr, c(3, 1, 2))
}

# Recruitment deviations from the legacy Perr_y (bias-corrected, initial equilibrium included)
.FitRecDevs <- function(fit, nAge, Years) {
  perr_y <- fit@OM@cpars$Perr_y
  if (is.null(perr_y)) return(NULL)
  maxage <- nAge - 1
  nyears <- length(Years$Hist)
  nproj <- min(length(Years$Proj), ncol(perr_y) - maxage - nyears)
  sims <- seq_len(nrow(perr_y))
  out <- list(
    Init = matrix(perr_y[, maxage:1], nrow(perr_y), dimnames = list(Sim = sims, Age = seq_len(maxage))),
    Hist = matrix(perr_y[, maxage + seq_len(nyears)], nrow(perr_y), dimnames = list(Sim = sims, Year = Years$Hist))
  )
  if (nproj > 0 && nproj == length(Years$Proj)) {
    out$Proj <- matrix(perr_y[, maxage + nyears + seq_len(nproj)], nrow(perr_y),
                       dimnames = list(Sim = sims, Year = Years$Proj))
  }
  out
}

.fit2om <- function(OM, fits, units, Data = NULL, silent = FALSE) {
  if (!inherits(OM, "om")) cli::cli_abort("{.arg OM} must be an {.cls om} object.")

  Years <- list(Hist = MSEtool::Years(OM, "Historical"), Proj = MSEtool::Years(OM, "Projection"))
  stock_index <- function(s) match(s, MSEtool::StockNames(OM))

  for (u in names(units)) {
    fit <- fits[[u]]
    if (!inherits(fit, "RCModel")) cli::cli_abort("{.arg fits}[[{.val {u}}]] must be an {.cls RCModel} object.")

    report <- fit@report
    if (length(report) == 1 && nrow(fit@SSB) > 1) report <- rep(report, nrow(fit@SSB))
    nyears <- fit@data@Misc$nyears
    nfleet <- fit@data@Misc$nfleet
    nAge <- ncol(report[[1]]$N)
    members <- units[[u]]$stocks
    devs <- .FitRecDevs(fit, nAge, Years)

    for (s in members) {
      si <- stock_index(s)
      stock <- OM@Stock[[si]]
      if (s == members[1]) {
        stock@SRR@R0 <- vapply(report, getElement, numeric(1), "R0")
        stock@SRR@Pars[[.SRRHSlot(stock@SRR@Model)]] <- vapply(report, getElement, numeric(1), "h")
        stock@SRR@SD <- vapply(report, getElement, numeric(1), "tau")
        if (!is.null(fit@OM@cpars$AC)) stock@SRR@AC <- fit@OM@cpars$AC
      }
      stock@Depletion@Final <- numeric()
      # Members share the cohort deviations so initial numbers follow the recruiting stock
      if (!is.null(devs)) {
        dimnames(devs$Init)$Age <- as.character(stock@Ages@Classes[-1])
        stock@SRR@RecDevInit <- devs$Init
        stock@SRR@RecDevHist <- devs$Hist
        if (!is.null(devs$Proj)) stock@SRR@RecDevProj <- devs$Proj
      }
      OM@Stock[[si]] <- stock

      fleets <- OM@Fleet[[si]]
      if (inherits(fleets, "fleet")) fleets <- stats::setNames(list(fleets), fleets@Name)
      fleet_names <- names(fleets) %||% paste0("Fleet", seq_len(nfleet))
      for (fi in seq_along(fleet_names)) {
        fl <- fleet_names[fi]
        blank <- MSEtool::Fleet(Name = fl)
        fleet <- fleets[[fl]] %||% blank
        fleet@Effort <- blank@Effort
        fleet@Effort@Effort <- .FitApicalF(report, fi, nyears)
        fleet@Catchability <- blank@Catchability
        fleet@Catchability@Efficiency <- array(1, c(1, 1))
        fleet@Selectivity <- blank@Selectivity
        fleet@Selectivity@MeanAtAge <- .FitVulArray(report, fi, nyears, nAge)
        fleet@Selectivity@isAtLength <- FALSE
        # Fitted F is for total removals
        fleet@Retention <- blank@Retention
        fleet@DiscardMortality <- blank@DiscardMortality
        fleets[[fl]] <- fleet
      }
      OM@Fleet[[si]] <- fleets
    }

    if (!is.null(Data[[u]])) {
      k <- units[[u]]$complex
      if (!length(OM@Data)) OM@Data <- list()
      OM@Data[[k]] <- Data[[u]]
      if (is.null(names(OM@Data)) || !nzchar(names(OM@Data)[k])) names(OM@Data)[k] <- u
    }
  }

  OM <- MSEtool::PopulateOM(OM, silent = silent)
  for (u in names(units)) {
    if (!is.null(Data[[u]])) OM <- .SetIndexSel(OM, fits[[u]], Data[[u]], units[[u]])
  }
  OM
}

# Index selectivity on the om from the RCM index selectivity type and fit
.SetIndexSel <- function(OM, fit, data, unit) {
  idx <- c("CPUE", "Survey")
  nms <- unlist(lapply(idx, function(sl) methods::slot(data, sl)@Name))
  type <- unlist(lapply(idx, function(sl) rep(sl, length(methods::slot(data, sl)@Name))))
  if (!length(nms)) return(OM)

  s_sel <- fit@data@Misc$s_selectivity %||% "B"
  if (length(s_sel) == 1) s_sel <- rep(s_sel, length(nms))
  Complex <- unit$complex
  si <- match(unit$stocks[1], MSEtool::StockNames(OM))
  Ages <- as.character(OM@Stock[[si]]@Ages@Classes)
  nyears <- fit@data@Misc$nyears
  fleet_names <- names(OM@Fleet[[si]])

  for (i in seq_along(nms)) {
    # Combined spawning output is not an om index type, so hermaphroditic SSB indices use the fitted array
    sel <- switch(as.character(s_sel[i]), B = "Biomass", SSB = if (!unit$herm) "SBiomass", NULL)
    fleet_match <- suppressWarnings(as.integer(s_sel[i]))
    same_fleet <- !is.na(fleet_match) && identical(fleet_names[fleet_match], nms[i])
    if (is.null(sel) && !same_fleet) {
      sel <- do.call(rbind, lapply(rep(fit@report, length.out = nrow(fit@SSB)), function(x) x$ivul[nyears, , i]))
      dimnames(sel) <- list(Sim = seq_len(nrow(sel)), Age = Ages)
    }
    if (!is.null(sel)) MSEtool::IndexSelectivity(OM, nms[i], Type = type[i], Complex = Complex) <- sel
  }
  OM
}
