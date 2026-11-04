
#' Rapid Conditioning Model for MSEtool operating models
#'
#' Conditions an operating model on fishery data with the Rapid Conditioning Model (see [RCM()]) and returns
#' a complete [MSEtool::om-class] object. The life history comes from [MSEtool::stock-class] objects and the
#' fishery data from [MSEtool::data-class] objects. `RCM2()` builds the operating model, fits each stock (or
#' hermaphroditic complex) and writes the fitted values into the operating model.
#'
#' @param Stock A [MSEtool::stock-class] object, or a named list of them. Alternatively, an
#' [MSEtool::om-class] object with fishery data in `OM@@Data`.
#' @param data A [MSEtool::data-class] object, or a list of them named by stock (or by complex, see
#' "Hermaphroditic stocks"). Not used when `Stock` is an [MSEtool::om-class] object with data.
#' @param Fleet Optional [MSEtool::fleet-class] object, or a list of them named by fleet, used for every stock.
#' Can also be a list of these named by stock. Fleet names are matched to `data@@Landings@@Name`. Fleet
#' selectivity parameters are used as starting values, and other fleet settings are kept in the returned
#' operating model.
#' @param Herm Optional [MSEtool::stocktransition-class] object (see [MSEtool::Herm()]), or a list of them,
#' for hermaphroditic stocks.
#' @param nSim Number of simulations. Defaults to the value in populated stock objects, otherwise 48.
#' @param pYear Number of projection years. Defaults to the value in populated stock objects, otherwise 50.
#' @param CurrentYear Last historical year. Defaults to the last year in `data@@Years`.
#' @param Obs Optional observation model for the returned operating model. See [MSEtool::OM()].
#' @param Imp Optional implementation model for the returned operating model. See [MSEtool::OM()].
#' @param Name Name of the returned operating model.
#' @param Seed Random seed for populating the stock objects.
#' @param stocks Optional character vector of stock or complex names to fit. Other stocks are returned unchanged.
#' @inheritParams RCM
#' @param mean_fit Logical, whether to run an additional fit with mean values of the life history parameters.
#' @param drop_nonconv,drop_highF Must be `FALSE`; every simulation is kept in the operating model.
#' @param ... Other arguments passed to [RCM()] for starting values and fixed parameters.
#'
#' @section Data:
#' Each stock (or hermaphroditic complex) is fitted to its own [MSEtool::data-class] object. Fleets correspond to
#' the columns of `data@@Landings@@Value`, plus `data@@Discards@@Value` if provided, and each fleet has its own
#' fitted selectivity and fishing mortality. Catch CVs in `@@CV` are optional.
#' Indices in `data@@CPUE` and `data@@Survey` are fitted in that order, so `s_selectivity` refers to the CPUE
#' indices first. Index CVs are required. Age and length compositions are taken from `data@@LandingsAtAge` and
#' `data@@LandingsAtSize`. The number of historical years and `CurrentYear` are taken from `data@@Years`, and the
#' first year of data is the first historical year.
#'
#' Arguments such as `condition`, `selectivity`, `s_selectivity`, `LWT`, `prior`, `start` and `map` apply to
#' every stock, or can be a list named by stock (or complex) to set values for each.
#'
#' @section Hermaphroditic stocks:
#' A hermaphroditic stock is two or more [MSEtool::stock-class] objects, one for each sex, linked by `Herm`.
#' The stocks share their age classes. Recruitment is to the stock that individuals start in (the `From` stock
#' that is not the `To` stock of any transition). `RCM2()` sets `SRR@@R0 = 0` for the other stocks and, if not
#' set, their `SRR@@SPFrom` to the recruiting stock. The linked stocks form one complex, named by joining the
#' stock names with `"-"`, and are fitted to one [MSEtool::data-class] object with catches and indices summed
#' over the stocks. In a list of data objects, this element is named by the complex or by any of its stocks.
#'
#' The complex is fitted as one population. Natural mortality, weight and length-at-age are averaged over the
#' stocks at age, weighted by the unfished numbers-at-age of each stock. Spawning output is that of the stocks
#' in `SRR@@SPFrom` of the recruiting stock. The fitted unfished recruitment, steepness and recruitment
#' deviations are assigned to the recruiting stock. Every stock in the complex gets the fitted fishing mortality
#' and selectivity-at-age of each fleet, so the proportion of each stock at age is the same as in the unfished
#' population.
#'
#' @section Operating model:
#' The returned operating model has the fitted unfished recruitment, steepness and recruitment deviations,
#' the fitted fishing mortality and selectivity-at-age of each fleet, and full retention of the catch. Index
#' selectivity is set in `@@Obs` from `s_selectivity`. The operating model can be used directly with
#' [MSEtool::Simulate()] and [MSEtool::runMSE()]. RCM is a single-area model, so operating models with more than
#' one area simulate a history close to, but not identical to, the fitted one.
#'
#' Assessment-based management procedures need the life-history data generated by
#' [MSEtool::LifeHistoryObs()] in the observation model.
#'
#' @return An [rcmodel-class] object with the conditioned operating model in `@@OM` and the [RCModel-class] fit for
#' each stock or complex in `@@fits`.
#'
#' @examples
#' \dontrun{
#' # Female and male stocks of a protogynous species, with data summed over both sexes
#' fit <- RCM2(list(Female = female, Male = male), data = dat,
#'             Herm = MSEtool::Herm(From = "Female", To = "Male", Frac = frac),
#'             nSim = 48, pYear = 50, selectivity = "logistic_length", s_selectivity = "B")
#' fit
#' plot(fit)
#' Hist <- MSEtool::Simulate(fit@OM)
#' }
#' @seealso [RCM()] [rcmodel-class]
#' @export
RCM2 <- function(Stock, data = NULL, Fleet = NULL, Herm = NULL, nSim = NULL, pYear = NULL, CurrentYear = NULL,
                 Obs = NULL, Imp = NULL, Name = "RCM", Seed = 101, stocks = NULL,
                 condition = "catch", selectivity = "logistic", s_selectivity = NULL, LWT = list(),
                 comp_like = c("multinomial", "lognormal", "mvlogistic", "dirmult1", "dirmult2"), prior = list(),
                 max_F = 3, integrate = FALSE, mean_fit = FALSE, drop_nonconv = FALSE, drop_highF = FALSE,
                 control = list(iter.max = 2e+05, eval.max = 4e+05), start = list(), map = list(),
                 cores = 1L, silent = FALSE, ...) {
  if (inherits(Stock, "om")) {
    OM <- Stock
    if (!is.null(data)) OM@Data <- if (inherits(data, "data")) list(data) else data
  } else {
    OM <- .BuildRCMOM(Stock, data, Fleet = Fleet, Herm = Herm, nSim = nSim, pYear = pYear,
                      CurrentYear = CurrentYear, Obs = Obs, Imp = Imp, Name = Name, Seed = Seed)
  }
  .RCMCore(OM, condition = condition, selectivity = selectivity, s_selectivity = s_selectivity, LWT = LWT,
           comp_like = comp_like, prior = prior, max_F = max_F, integrate = integrate, mean_fit = mean_fit,
           drop_nonconv = drop_nonconv, drop_highF = drop_highF, control = control, start = start, map = map,
           cores = cores, stocks = stocks, silent = silent, ...)
}

# Groups of stocks linked by Herm, in stock order
.HermGroups <- function(herm, stock_names) {
  grp <- stats::setNames(seq_along(stock_names), stock_names)
  for (h in herm) grp[grp == grp[[h@To]]] <- grp[[h@From]]
  unname(split(stock_names, factor(grp, levels = unique(grp))))
}

# Data object for each complex, from one data object or a list named by complex, stock or position
.ComplexData <- function(data, comps, cx_names) {
  if (inherits(data, "data")) {
    if (length(comps) > 1) cli::cli_abort("Supply {.arg data} as a list named by stock for {length(comps)} stocks.")
    return(stats::setNames(list(data), cx_names))
  }
  if (!is.list(data) || !length(data)) cli::cli_abort("{.arg data} must be a {.cls data} object or a list of them.")
  out <- lapply(seq_along(comps), function(k) {
    nm <- intersect(c(cx_names[k], comps[[k]]), names(data))
    if (length(nm)) return(data[[nm[1]]])
    if (is.null(names(data)) && length(data) == length(comps)) return(data[[k]])
    cli::cli_abort("No {.cls data} object for {.val {cx_names[k]}}.")
  })
  stats::setNames(out, cx_names)
}

# A complete om from stock, fleet, data and Herm objects
.BuildRCMOM <- function(Stock, data, Fleet = NULL, Herm = NULL, nSim = NULL, pYear = NULL, CurrentYear = NULL,
                        Obs = NULL, Imp = NULL, Name = "RCM", Seed = 101) {
  Stock <- .AsStockList(Stock)
  stock_names <- names(Stock)
  for (s in stock_names) Stock[[s]]@Name <- s

  herm <- if (inherits(Herm, "stocktransition")) list(Herm) else as.list(Herm)
  herm <- lapply(herm, function(h) {
    if (!is.character(h@From)) h@From <- stock_names[h@From]
    if (!is.character(h@To)) h@To <- stock_names[h@To]
    bad <- setdiff(c(h@From, h@To), stock_names)
    if (length(bad)) cli::cli_abort("{.arg Herm} stock{?s} {.val {bad}} not in {.arg Stock} ({.val {stock_names}}).")
    h
  })
  comps <- .HermGroups(herm, stock_names)
  cx_names <- vapply(comps, paste, character(1), collapse = "-")

  to <- vapply(herm, function(h) h@To, character(1))
  for (members in comps[lengths(comps) > 1]) {
    rec <- setdiff(members, to)
    if (length(rec) != 1) {
      cli::cli_abort("The {.arg Herm} transitions among {.val {members}} must start from one stock.")
    }
    for (s in setdiff(members, rec)) {
      Stock[[s]]@SRR@R0 <- 0
      if (!length(Stock[[s]]@SRR@SPFrom)) Stock[[s]]@SRR@SPFrom <- rec
    }
  }

  data <- .ComplexData(data, comps, cx_names)
  yrs <- .DataYears(data[[1]])
  pop <- Filter(.IsPopulated, Stock)
  if (!length(yrs)) {
    if (!length(pop)) cli::cli_abort("Set {.code data@@Years} to define the historical years.")
    nYear <- pop[[1]]@nYear
    CurrentYear <- CurrentYear %||% pop[[1]]@CurrentYear
  } else {
    nYear <- length(yrs)
    CurrentYear <- CurrentYear %||% max(yrs)
  }
  if (is.null(nSim)) nSim <- if (length(pop)) pop[[1]]@nSim else 48
  if (is.null(pYear)) pYear <- if (length(pop)) pop[[1]]@pYear else 50

  per_stock <- is.list(Fleet) && !inherits(Fleet, "fleet") && length(names(Fleet)) && all(names(Fleet) %in% stock_names)
  FleetList <- lapply(stock_names, function(s) {
    k <- which(vapply(comps, function(m) s %in% m, logical(1)))
    fl_in <- if (per_stock) Fleet[[s]] %||% Fleet[[comps[[k]][1]]] else Fleet
    fleet_names <- .ResolveFleetNames(data[[k]], fl_in)
    fl <- .AlignFleetList(fl_in, fleet_names)
    stats::setNames(lapply(fleet_names, function(nm) fl[[nm]] %||% MSEtool::Fleet(Name = nm)), fleet_names)
  })
  names(FleetList) <- stock_names

  args <- list(Name = Name, nSim = nSim, nYear = nYear, pYear = pYear, CurrentYear = CurrentYear,
               Stock = Stock, Fleet = FleetList, Data = data,
               Complexes = stats::setNames(lapply(comps, match, stock_names), cx_names), Seed = Seed)
  if (length(herm)) args$Herm <- herm
  if (!is.null(Obs)) args$Obs <- Obs
  if (!is.null(Imp)) args$Imp <- Imp
  do.call(MSEtool::OM, args)
}
