
#' Rapid Conditioning Model for MSEtool operating models
#'
#' Conditions an operating model on fishery data with the Rapid Conditioning Model (see [RCM()]) and returns
#' a complete [MSEtool::om-class] object. The life history comes from one or more [MSEtool::stock-class] objects
#' and the fishery data from one [MSEtool::data-class] object. `RCM2()` builds the operating model, fits the
#' stocks as one population to the data, and writes the fitted values into the operating model.
#'
#' @param Stock A [MSEtool::stock-class] object, or a named list of them. The names become the stock names.
#' Alternatively, an [MSEtool::om-class] object with fishery data in `OM@@Data`.
#' @param data A [MSEtool::data-class] object with catches and indices summed over all stocks.
#' Not used when `Stock` is an [MSEtool::om-class] object with data.
#' @param Fleet Optional [MSEtool::fleet-class] object, or a list of them named by fleet, used for every stock.
#' Fleet names are matched to `data@@Landings@@Name`. Fleet selectivity parameters are used as starting values,
#' and other fleet settings are kept in the returned operating model.
#' @param Herm Optional [MSEtool::stocktransition-class] object (see [MSEtool::Herm()]), or a list of them,
#' for hermaphroditic stocks. `From` and `To` are names (or positions) in `Stock`.
#' @param nSim Number of simulations. Defaults to the value in populated stock objects, otherwise 48.
#' @param pYear Number of projection years. Defaults to the value in populated stock objects, otherwise 50.
#' @param CurrentYear Last historical year. Defaults to the last year in `data@@Years`.
#' @param Obs Optional observation model for the returned operating model. See [MSEtool::OM()].
#' @param Imp Optional implementation model for the returned operating model. See [MSEtool::OM()].
#' @param Name Name of the returned operating model.
#' @param Seed Random seed for populating the stock objects.
#' @inheritParams RCM
#' @param mean_fit Logical, whether to run an additional fit with mean values of the life history parameters.
#' @param drop_nonconv,drop_highF Must be `FALSE`; every simulation is kept in the operating model.
#' @param ... Other arguments passed to [RCM()] for starting values and fixed parameters.
#'
#' @section Data:
#' `RCM2()` fits one [MSEtool::data-class] object. Fleets correspond to the columns of `data@@Landings@@Value`,
#' plus `data@@Discards@@Value` if provided, and each fleet has its own fitted selectivity and fishing mortality.
#' Catch CVs in `@@CV` are optional. Indices in `data@@CPUE` and `data@@Survey` are fitted in that order, so
#' `s_selectivity` refers to the CPUE indices first. Index CVs are required. Age and length compositions are
#' taken from `data@@LandingsAtAge` and `data@@LandingsAtSize`. The number of historical years and `CurrentYear`
#' are taken from `data@@Years`, and the first year of data is the first historical year.
#'
#' @section Several stocks:
#' Several [MSEtool::stock-class] objects, e.g., the sexes of a stock or its components, are fitted together as
#' one population (a complex) to the one data object, with catches and indices summed over the stocks. The
#' stocks share their age classes. The complex is named by joining the stock names with `"-"`, and the returned
#' operating model has every stock in that complex in `@@Complexes`.
#'
#' Natural mortality, weight and length-at-age are averaged over the stocks at age, weighted by the unfished
#' numbers-at-age of each stock. Spawning output is set by `SRR@@SPFrom` of the recruiting stocks, as in MSEtool: a
#' stock name, a position in `Stock`, or a weight vector named by stock (name or position). An empty `SRR@@SPFrom` is
#' the stock itself. Recruiting stocks with the same `SRR@@SPFrom` recruit from the same spawning output. Recruiting
#' stocks with an empty `SRR@@SPFrom` each recruit from their own spawning output, and the fit uses their sum.
#' Every stock gets the fitted fishing mortality and selectivity-at-age of each fleet, so the proportion of each
#' stock at age is the same as in the unfished population.
#'
#' Recruitment is to the stocks that are not the `To` stock of any `Herm` transition. `RCM2()` sets
#' `SRR@@R0 = 0` for the other stocks. Each recruiting stock gets the fitted unfished recruitment times its share of
#' the input `SRR@@R0` (equal shares if no stock has `SRR@@R0`), and the fitted steepness and recruitment
#' deviations. An empty `SRR@@SPFrom` is set to the recruiting stocks. The recruiting stocks must have the same
#' `SRR@@SPFrom` and `SRR@@Model`.
#'
#' For a hermaphroditic stock, link the stocks (one for each sex) with `Herm`. Recruitment is then to the stock
#' that individuals start in, from the spawning output of the stocks in its `SRR@@SPFrom` (itself if empty).
#' `RCM2()` does not change a `SRR@@SPFrom` that is set. For a protogynous stock (`Herm(From = "Female", To = "Male")`),
#' the female stock recruits from its own spawning output by default. For a protandrous stock
#' (`Herm(From = "Male", To = "Female")`), set `SRR@@SPFrom = "Female"` in the male stock so that recruitment comes
#' from female spawning output.
#'
#' When `Stock` is a named list, the names replace `Stock@@Name`. Stock names in `SRR@@SPFrom`, `Herm@@From` and
#' `Herm@@To` are updated to match.
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
#' @return An [rcmodel-class] object with the conditioned operating model in `@@OM` and the [RCModel-class] fit in
#' `@@fits`, named by the stock or complex.
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
                 Obs = NULL, Imp = NULL, Name = "RCM", Seed = 101,
                 condition = "catch", selectivity = "logistic", s_selectivity = NULL, LWT = list(),
                 comp_like = c("multinomial", "lognormal", "mvlogistic", "dirmult1", "dirmult2"), prior = list(),
                 max_F = 3, integrate = FALSE, mean_fit = FALSE, drop_nonconv = FALSE, drop_highF = FALSE,
                 control = list(iter.max = 2e+05, eval.max = 4e+05), start = list(), map = list(),
                 cores = 1L, silent = FALSE, ...) {
  if (inherits(Stock, "om")) {
    OM <- Stock
    if (!is.null(data)) OM@Data <- list(.OneData(data))
  } else {
    OM <- .BuildRCMOM(Stock, data, Fleet = Fleet, Herm = Herm, nSim = nSim, pYear = pYear,
                      CurrentYear = CurrentYear, Obs = Obs, Imp = Imp, Name = Name, Seed = Seed)
  }
  .RCMCore(OM, condition = condition, selectivity = selectivity, s_selectivity = s_selectivity, LWT = LWT,
           comp_like = comp_like, prior = prior, max_F = max_F, integrate = integrate, mean_fit = mean_fit,
           drop_nonconv = drop_nonconv, drop_highF = drop_highF, control = control, start = start, map = map,
           cores = cores, silent = silent, ...)
}

# Replace old stock names with new ones in SRR@SPFrom (values, or names of a weight vector)
.RenameSPFrom <- function(SPFrom, renames) {
  if (!length(SPFrom) || !length(renames)) return(SPFrom)
  swap <- function(x) ifelse(x %in% names(renames), renames[x], x)
  if (length(names(SPFrom)) && any(nzchar(names(SPFrom)))) {
    names(SPFrom) <- swap(names(SPFrom))
  } else if (is.character(SPFrom)) {
    SPFrom <- unname(swap(SPFrom))
  }
  SPFrom
}

# A complete om from stock, fleet, data and Herm objects
.BuildRCMOM <- function(Stock, data, Fleet = NULL, Herm = NULL, nSim = NULL, pYear = NULL, CurrentYear = NULL,
                        Obs = NULL, Imp = NULL, Name = "RCM", Seed = 101) {
  if (is.null(data)) cli::cli_abort("Supply {.arg data}.")
  data <- .OneData(data)
  Stock <- .AsStockList(Stock)
  stock_names <- names(Stock)

  # Old names that the list names replace
  old <- vapply(Stock, function(x) if (length(x@Name)) x@Name[1] else NA_character_, character(1))
  renames <- stats::setNames(stock_names, old)[!is.na(old) & nzchar(old) & old != stock_names]
  renames <- renames[!names(renames) %in% stock_names]
  for (s in stock_names) {
    Stock[[s]]@Name <- s
    Stock[[s]]@SRR@SPFrom <- .RenameSPFrom(Stock[[s]]@SRR@SPFrom, renames)
  }

  herm <- if (inherits(Herm, "stocktransition")) list(Herm) else as.list(Herm)
  herm <- lapply(herm, function(h) {
    for (sl in c("From", "To")) {
      x <- methods::slot(h, sl)
      x <- if (is.character(x)) unname(ifelse(x %in% names(renames), renames[x], x)) else stock_names[x]
      methods::slot(h, sl) <- x
    }
    bad <- setdiff(c(h@From, h@To), stock_names)
    if (length(bad)) cli::cli_abort("{.arg Herm} stock{?s} {.val {bad}} not in {.arg Stock} ({.val {stock_names}}).")
    h
  })

  yrs <- .DataYears(data)
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

  fleet_names <- .ResolveFleetNames(data, Fleet)
  fl <- .AlignFleetList(Fleet, fleet_names)
  fl <- stats::setNames(lapply(fleet_names, function(nm) fl[[nm]] %||% MSEtool::Fleet(Name = nm)), fleet_names)
  FleetList <- stats::setNames(rep(list(fl), length(stock_names)), stock_names)

  cx_name <- paste(stock_names, collapse = "-")
  args <- list(Name = Name, nSim = nSim, nYear = nYear, pYear = pYear, CurrentYear = CurrentYear,
               Stock = Stock, Fleet = FleetList, Data = stats::setNames(list(data), cx_name),
               Complexes = stats::setNames(list(seq_along(stock_names)), cx_name), Seed = Seed)
  if (length(herm)) args$Herm <- herm
  if (!is.null(Obs)) args$Obs <- Obs
  if (!is.null(Imp)) args$Imp <- Imp
  .PrepareComplex(do.call(MSEtool::OM, args))
}
