
# Engine for RCM2(): fit every stock in an om as one population with RCM_int(), return an rcmodel

.ListOf <- function(x) if (isS4(x)) list(x) else as.list(x)

.AsStockList <- function(Stock) {
  if (inherits(Stock, "stock")) {
    nm <- Stock@Name
    if (is.null(nm) || !nzchar(nm)) nm <- "Stock"
    return(stats::setNames(list(Stock), nm))
  }
  if (is.list(Stock)) {
    if (is.null(names(Stock)) || any(!nzchar(names(Stock)))) {
      cli::cli_abort("A list of {.cls stock} objects must be named.")
    }
    return(Stock)
  }
  cli::cli_abort("{.arg Stock} must be a {.cls stock} object or a named list of {.cls stock} objects.")
}

# The one data object, from a data object or a list of length one
.OneData <- function(data, arg = "data") {
  if (inherits(data, "data")) return(data)
  if (is.list(data) && length(data) == 1L && inherits(data[[1]], "data")) return(data[[1]])
  if (is.list(data) && length(data) > 1L) {
    cli::cli_abort(c("{.fn RCM2} fits one data set.",
                     "x" = "{.arg {arg}} has {length(data)} {.cls data} objects.",
                     "i" = "Sum the catches and indices over the stocks into one {.cls data} object."))
  }
  cli::cli_abort("{.arg {arg}} must be a {.cls data} object.")
}

# Fleet names from the data, then the Fleet list, then Fleet1, Fleet2, ...
.ResolveFleetNames <- function(data, Fleet) {
  nm <- data@Landings@Name
  if (!length(nm)) nm <- data@Discards@Name
  if (!length(nm) && is.list(Fleet) && !inherits(Fleet, "fleet")) nm <- names(Fleet)
  if (!length(nm)) {
    nfleet <- 1L
    if (length(data@Landings@Value)) {
      nfleet <- ncol(as.matrix(data@Landings@Value))
    } else if (length(data@Discards@Value)) {
      nfleet <- ncol(as.matrix(data@Discards@Value))
    } else if (length(data@Effort@Value)) {
      nfleet <- ncol(as.matrix(data@Effort@Value))
    } else if (is.list(Fleet) && !inherits(Fleet, "fleet")) {
      nfleet <- length(Fleet)
    }
    nm <- paste0("Fleet", seq_len(nfleet))
  }
  nm
}

.AlignFleetList <- function(Fleet, fleet_names) {
  out <- stats::setNames(vector("list", length(fleet_names)), fleet_names)
  if (is.null(Fleet)) return(out)
  if (inherits(Fleet, "fleet")) {
    out[[1]] <- Fleet
    return(out)
  }
  for (nm in intersect(fleet_names, names(Fleet))) out[[nm]] <- Fleet[[nm]]
  out
}

.IsPopulated <- function(Stock) length(Stock@nYear) > 0 && length(Stock@Length@MeanAtAge) > 0

.DataYears <- function(data) {
  yrs <- data@Years
  if (!length(yrs)) {
    for (sl in c("Landings", "Discards", "Effort")) {
      v <- methods::slot(data, sl)@Value
      if (length(v) && length(rownames(v))) return(as.numeric(rownames(v)))
    }
  }
  yrs
}

# Populate an unpopulated stock over the years in its data
.PrepareStock <- function(Stock, data, om_meta, silent) {
  if (.IsPopulated(Stock)) return(Stock)
  yrs <- .DataYears(data)
  nYear <- om_meta$nYear %||% length(yrs)
  CurrentYear <- om_meta$CurrentYear %||% if (length(yrs)) max(yrs)
  if (!length(nYear) || !nYear) {
    cli::cli_abort("Supply {.arg nYear} or {.code data@Years} to populate stock {.val {Stock@Name}}.")
  }
  MSEtool::PopulateStock(Stock, nYear = nYear, pYear = om_meta$pYear %||% 0, CurrentYear = CurrentYear,
                         nSim = om_meta$nSim %||% 48, seed = om_meta$Seed %||% 101, silent = silent)
}

.CheckDataYears <- function(data, Stock, st) {
  for (sl in c("Landings", "Discards", "Effort")) {
    v <- methods::slot(data, sl)@Value
    if (length(v) && NROW(v) != Stock@nYear) {
      cli::cli_abort("{.code data@{sl}@Value} for {.val {st}} has {NROW(v)} row{?s} but the stock has {Stock@nYear} historical year{?s}.")
    }
  }
  invisible(TRUE)
}


# Fit every stock in an om as one population to OM@Data with RCM_int(), then write the fit back into the om
.RCMCore <- function(OM, condition = "catch", selectivity = "logistic", s_selectivity = NULL,
                     LWT = list(), comp_like = c("multinomial", "lognormal", "mvlogistic", "dirmult1", "dirmult2"),
                     prior = list(), max_F = 3, integrate = FALSE, mean_fit = FALSE,
                     drop_nonconv = FALSE, drop_highF = FALSE,
                     control = list(iter.max = 2e5, eval.max = 4e5), start = list(), map = list(),
                     cores = 1L, silent = FALSE, ...) {
  if (!is.list(comp_like)) comp_like <- match.arg(comp_like)
  if (OM@Seasons > 1) cli::cli_abort("RCM fits annual models only; {.code OM@Seasons} must be 1.")
  if (!length(OM@Data)) cli::cli_abort("{.code OM@Data} is empty. Add fishery data to {.arg OM} first.")
  data <- .OneData(OM@Data, "OM@Data")

  OM <- .PrepareComplex(OM)
  all_stocks <- MSEtool::StockNames(OM)
  Stock <- stats::setNames(.ListOf(OM@Stock), all_stocks)
  om_meta <- list(nSim = OM@nSim, nYear = OM@nYear, pYear = OM@pYear, CurrentYear = OM@CurrentYear, Seed = OM@Seed)
  for (s in all_stocks) Stock[[s]] <- .PrepareStock(Stock[[s]], data, om_meta, silent)

  units <- .FitUnit(OM, Stock)
  name <- names(units)
  unit <- units[[1]]
  stk <- if (unit$combined) .CombineStocks(OM, Stock, unit, name) else Stock[[unit$stocks[1]]]
  Fleet <- OM@Fleet[[match(unit$stocks[1], all_stocks)]]
  nSim <- OM@nSim

  fleet_names <- .ResolveFleetNames(data, Fleet)
  dat <- .FillDataNames(data, fleet_names)
  .CheckDataYears(dat, stk, name)
  fl_aligned <- .AlignFleetList(Fleet, fleet_names)
  fl_for_convert <- if (length(fl_aligned) == 1L) fl_aligned[[1]] else fl_aligned

  RCMdata <- .data2RCMdata(dat, condition = condition)
  OM_fit <- .stock2OM(stk, fl_for_convert, nSim = nSim, pYear = OM@pYear, CurrentYear = OM@CurrentYear)

  if (is.null(start$vul_par) && length(fleet_names) > 1L) {
    start$vul_par <- .FleetSelStart(fl_for_convert, nSim)$vul_par
  }

  fit <- RCM_int(OM = OM_fit, RCMdata = RCMdata, condition = condition,
                 selectivity = selectivity, s_selectivity = s_selectivity,
                 LWT = LWT, comp_like = comp_like, prior = prior,
                 max_F = max_F, cores = cores, integrate = integrate,
                 mean_fit = mean_fit, drop_nonconv = drop_nonconv,
                 drop_highF = drop_highF, control = control,
                 start = start, map = map, silent = silent, ...)

  if (nrow(fit@SSB) != nSim) {
    cli::cli_abort(c("{.val {name}} kept {nrow(fit@SSB)} of {nSim} simulations.",
                     "i" = "Set {.arg drop_nonconv} and {.arg drop_highF} to {.code FALSE} to keep every simulation in the operating model."))
  }

  # Keep the biology that was fitted
  for (s in all_stocks) OM@Stock[[match(s, all_stocks)]] <- Stock[[s]]
  OM@Data <- list()

  fits <- stats::setNames(list(fit), name)
  OM_out <- .fit2om(OM, fits, units, Data = stats::setNames(list(dat), name), silent = silent)
  methods::new("rcmodel", OM = OM_out, fits = fits, RCMdata = lapply(fits, methods::slot, "data"))
}
