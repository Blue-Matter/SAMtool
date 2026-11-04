
# Engine for RCM2(): fit each stock or hermaphroditic complex with RCM_int(), return an rcmodel

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
  cli::cli_abort("{.arg OM} must be a {.cls stock} object or a named list of {.cls stock} objects.")
}

# Key data/fleet input by stock name; a list named by stocks is kept as is
.AsStockKeyedList <- function(x, stock_names, cls) {
  if (is.null(x)) return(stats::setNames(vector("list", length(stock_names)), stock_names))
  if (inherits(x, cls)) return(stats::setNames(list(x), stock_names[1]))
  if (is.list(x)) {
    if (length(names(x)) && any(names(x) %in% stock_names)) return(x)
    if (is.null(names(x)) && length(x) == length(stock_names) && all(vapply(x, inherits, logical(1), cls))) {
      return(stats::setNames(x, stock_names))
    }
    return(stats::setNames(list(x), stock_names[1]))
  }
  cli::cli_abort("Expected a {.cls {cls}} object or a list of them named by stock.")
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

# Per-stock value from a list named by stock, otherwise the same value for every stock
.PerStockArg <- function(arg, stock_names, default) {
  is_per_stock <- is.list(arg) && length(names(arg)) > 0 && any(names(arg) %in% stock_names)
  out <- if (is_per_stock) {
    lapply(stock_names, function(st) if (st %in% names(arg)) arg[[st]] else default)
  } else {
    lapply(stock_names, function(st) arg)
  }
  stats::setNames(out, stock_names)
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
      cli::cli_abort("{.code data@{sl}@Value} for stock {.val {st}} has {NROW(v)} row{?s} but the stock has {Stock@nYear} historical year{?s}.")
    }
  }
  invisible(TRUE)
}


# Fit each stock or hermaphroditic complex in an om with RCM_int(), then write the fits back into the om
.RCMCore <- function(OM, condition = "catch", selectivity = "logistic", s_selectivity = NULL,
                     LWT = list(), comp_like = c("multinomial", "lognormal", "mvlogistic", "dirmult1", "dirmult2"),
                     prior = list(), max_F = 3, integrate = FALSE, mean_fit = FALSE,
                     drop_nonconv = FALSE, drop_highF = FALSE,
                     control = list(iter.max = 2e5, eval.max = 4e5), start = list(), map = list(),
                     cores = 1L, stocks = NULL, silent = FALSE, ...) {
  if (!is.list(comp_like)) comp_like <- match.arg(comp_like)
  if (OM@Seasons > 1) cli::cli_abort("RCM fits annual models only; {.code OM@Seasons} must be 1.")
  if (!length(OM@Data)) cli::cli_abort("{.code OM@Data} is empty. Add fishery data to {.arg OM} first.")

  units <- .FitUnits(OM)
  Stock <- stats::setNames(.ListOf(OM@Stock), MSEtool::StockNames(OM))
  all_stocks <- names(Stock)
  Fleet <- .AsStockKeyedList(OM@Fleet, all_stocks, "fleet")
  om_data <- .ListOf(OM@Data)
  data <- lapply(seq_along(units), function(k) {
    u <- names(units)[k]
    if (length(names(om_data)) && u %in% names(om_data)) om_data[[u]] else if (k <= length(om_data)) om_data[[k]]
  })
  names(data) <- names(units)
  om_meta <- list(nSim = OM@nSim, nYear = OM@nYear, pYear = OM@pYear, CurrentYear = OM@CurrentYear, Seed = OM@Seed)

  keep <- vapply(names(units), function(u) is.null(stocks) || any(c(u, units[[u]]$stocks) %in% stocks), logical(1))
  stock_names <- names(units)[keep]
  if (!length(stock_names)) cli::cli_abort("None of {.arg stocks} match the stock names {.val {all_stocks}}.")

  missing_data <- stock_names[vapply(stock_names, function(st) is.null(data[[st]]), logical(1))]
  if (length(missing_data)) cli::cli_abort("No {.cls data} object for {.val {missing_data}}.")

  for (st in stock_names) {
    for (s in units[[st]]$stocks) Stock[[s]] <- .PrepareStock(Stock[[s]], data[[st]], om_meta, silent)
  }
  FitStock <- lapply(stats::setNames(stock_names, stock_names), function(st) {
    u <- units[[st]]
    if (u$herm) .CombineHermStock(OM, Stock, u, st) else Stock[[u$stocks[1]]]
  })
  Fleet <- lapply(stats::setNames(stock_names, stock_names), function(st) Fleet[[units[[st]]$stocks[1]]])
  nSim <- OM@nSim

  condition_l <- .PerStockArg(condition, stock_names, "catch")
  selectivity_l <- .PerStockArg(selectivity, stock_names, "logistic")
  s_selectivity_l <- .PerStockArg(s_selectivity, stock_names, NULL)
  LWT_l <- .PerStockArg(LWT, stock_names, list())
  comp_like_l <- .PerStockArg(comp_like, stock_names, "multinomial")
  prior_l <- .PerStockArg(prior, stock_names, list())
  max_F_l <- .PerStockArg(max_F, stock_names, 3)
  integrate_l <- .PerStockArg(integrate, stock_names, FALSE)
  mean_fit_l <- .PerStockArg(mean_fit, stock_names, FALSE)
  drop_nonconv_l <- .PerStockArg(drop_nonconv, stock_names, FALSE)
  drop_highF_l <- .PerStockArg(drop_highF, stock_names, FALSE)
  control_l <- .PerStockArg(control, stock_names, list(iter.max = 2e5, eval.max = 4e5))
  start_l <- .PerStockArg(start, stock_names, list())
  map_l <- .PerStockArg(map, stock_names, list())

  # Parallel over stocks (forked, so not on Windows), otherwise over simulations in RCM_int()
  cores_stock <- if (.Platform$OS.type == "windows") 1L else max(1L, min(length(stock_names), cores))
  cores_sim <- if (cores_stock > 1L) 1L else cores

  fit_one <- function(st) {
    stk <- FitStock[[st]]
    fleet_names_st <- .ResolveFleetNames(data[[st]], Fleet[[st]])
    dat <- .FillDataNames(data[[st]], fleet_names_st)
    .CheckDataYears(dat, stk, st)
    fl_aligned <- .AlignFleetList(Fleet[[st]], fleet_names_st)
    fl_for_convert <- if (length(fl_aligned) == 1L) fl_aligned[[1]] else fl_aligned

    RCMdata_st <- .data2RCMdata(dat, condition = condition_l[[st]])
    OM_st <- .stock2OM(stk, fl_for_convert, nSim = nSim, pYear = OM@pYear, CurrentYear = OM@CurrentYear)

    start_st <- start_l[[st]]
    if (is.null(start_st$vul_par) && length(fleet_names_st) > 1L) {
      start_st$vul_par <- .FleetSelStart(fl_for_convert, nSim)$vul_par
    }

    fit <- RCM_int(OM = OM_st, RCMdata = RCMdata_st, condition = condition_l[[st]],
                   selectivity = selectivity_l[[st]], s_selectivity = s_selectivity_l[[st]],
                   LWT = LWT_l[[st]], comp_like = comp_like_l[[st]], prior = prior_l[[st]],
                   max_F = max_F_l[[st]], cores = cores_sim, integrate = integrate_l[[st]],
                   mean_fit = mean_fit_l[[st]], drop_nonconv = drop_nonconv_l[[st]],
                   drop_highF = drop_highF_l[[st]], control = control_l[[st]],
                   start = start_st, map = map_l[[st]], silent = silent, ...)
    list(fit = fit, data = dat)
  }

  results <- if (cores_stock > 1L) {
    parallel::mclapply(stock_names, fit_one, mc.cores = cores_stock)
  } else {
    lapply(stock_names, fit_one)
  }
  names(results) <- stock_names
  fits <- lapply(results, getElement, "fit")

  nsim_fit <- vapply(fits, function(x) nrow(x@SSB), integer(1))
  if (any(nsim_fit != nSim)) {
    cli::cli_abort(c("{.val {stock_names[nsim_fit != nSim][1]}} kept {nsim_fit[nsim_fit != nSim][1]} of {nSim} simulations.",
                     "i" = "Set {.arg drop_nonconv} and {.arg drop_highF} to {.code FALSE} to keep every simulation in the operating model."))
  }

  # Keep the biology that was fitted
  for (s in unlist(lapply(units[stock_names], getElement, "stocks"))) OM@Stock[[match(s, all_stocks)]] <- Stock[[s]]

  OM_out <- .fit2om(OM, fits, units[stock_names], Data = lapply(results, getElement, "data"), silent = silent)
  methods::new("rcmodel", OM = OM_out, fits = fits, RCMdata = lapply(fits, methods::slot, "data"))
}
