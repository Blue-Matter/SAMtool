
# The fitting unit for RCM2(): every stock in the om, fitted as one population to one data set

# Herm pairs as stock names
.HermPairNames <- function(OM, stock_names) {
  herm <- OM@Herm
  if (inherits(herm, "stocktransition")) herm <- list(herm)
  nm <- function(x) if (is.character(x)) x else stock_names[x]
  lapply(herm, function(h) {
    h@From <- nm(h@From)
    h@To <- nm(h@To)
    h
  })
}

# Stocks that recruit: those that are not the `To` stock of any Herm transition
.RecruitingStocks <- function(OM, stock_names = MSEtool::StockNames(OM)) {
  to <- vapply(.HermPairNames(OM, stock_names), function(h) h@To, character(1))
  setdiff(stock_names, to)
}

# The recruiting stock that each stock comes from, following the Herm transitions back
.RecruitingSource <- function(OM, stock_names, rec) {
  herm <- .HermPairNames(OM, stock_names)
  from <- stats::setNames(vapply(herm, function(h) h@From, character(1)), vapply(herm, function(h) h@To, character(1)))
  vapply(stock_names, function(s) {
    seen <- s
    while (!s %in% rec && s %in% names(from)) {
      s <- from[[s]]
      if (s %in% seen) break
      seen <- c(seen, s)
    }
    s
  }, character(1))
}

# SRR@SPFrom with stock names in place of positions, as MSEtool's populate step does: a numeric position, or
# names of a weight vector that are numeric strings
.NameSPFrom <- function(SPFrom, stock_names, stock) {
  if (!length(SPFrom)) return(SPFrom)
  pos <- function(x) vapply(as.character(x), function(v) {
    i <- suppressWarnings(as.numeric(v))
    if (!is.na(i) && i >= 1 && i <= length(stock_names) && i == round(i)) stock_names[i] else v
  }, character(1), USE.NAMES = FALSE)
  if (length(names(SPFrom)) && any(nzchar(names(SPFrom)))) {
    names(SPFrom) <- pos(names(SPFrom))
    src <- names(SPFrom)
  } else {
    SPFrom <- pos(SPFrom)
    src <- SPFrom
  }
  bad <- setdiff(src, stock_names)
  if (length(bad)) {
    cli::cli_abort("{.code SRR@SPFrom} of stock {.val {stock}} names {.val {bad}}, which {?is/are} not in the stocks ({.val {stock_names}}).")
  }
  SPFrom
}

# Spawning sources of the fitted population as named weights, from SRR@SPFrom of the recruiting stocks.
# Recruitment depends on spawning output relative to unfished, so the weights set the relative contribution of
# each source
.SpawningSources <- function(OM, rec, stock_names) {
  raw <- lapply(rec, function(s) .SPFromWeights(OM@Stock[[match(s, stock_names)]]@SRR@SPFrom, s, stock_names))
  w <- lapply(raw, function(x) {
    x <- x[order(names(x))]
    x / sum(x)
  })
  same <- vapply(w, function(x) identical(names(x), names(w[[1]])) && isTRUE(all.equal(unname(x), unname(w[[1]]))), logical(1))
  if (all(same)) return(raw[[1]])
  # Each recruits from its own spawning output: the population spawns from all of them
  self <- vapply(seq_along(rec), function(k) identical(names(w[[k]]), rec[k]), logical(1))
  if (all(self)) return(stats::setNames(rep(1, length(rec)), rec))
  cli::cli_abort(c("The recruiting stocks {.val {rec}} are fitted with one stock-recruit relationship.",
                   "x" = "Their {.code SRR@SPFrom} name different spawning stocks.",
                   "i" = "Give them the same {.code SRR@SPFrom}, or none so that each recruits from its own spawning output."))
}

# Every stock in one complex. Stocks that do not recruit get SRR@R0 = 0, and an empty SRR@SPFrom is set to
# the recruiting stock they come from
.PrepareComplex <- function(OM) {
  stock_names <- MSEtool::StockNames(OM)
  cx <- OM@Complexes
  if (length(cx) > 1L) {
    cli::cli_abort(c("{.fn RCM2} fits one data set, so every stock must be in one complex.",
                     "x" = "{.code OM@Complexes} has {length(cx)} complexes."))
  }
  name <- if (length(cx)) names(cx)[1] else NULL
  if (length(cx)) {
    members <- cx[[1]]
    if (!is.character(members)) members <- stock_names[members]
    if (!setequal(members, stock_names)) {
      cli::cli_abort("{.fn RCM2} fits one data set, so {.code OM@Complexes} must include every stock.")
    }
  }
  if (is.null(name) || !nzchar(name)) name <- paste(stock_names, collapse = "-")
  OM@Complexes <- stats::setNames(list(seq_along(stock_names)), name)
  for (i in seq_along(stock_names)) {
    OM@Stock[[i]]@SRR@SPFrom <- .NameSPFrom(OM@Stock[[i]]@SRR@SPFrom, stock_names, stock_names[i])
  }

  rec <- .RecruitingStocks(OM, stock_names)
  if (!length(rec)) cli::cli_abort("Every stock is the {.code To} stock of a {.arg Herm} transition; one stock must recruit.")
  source <- .RecruitingSource(OM, stock_names, rec)
  # Stocks whose spawning output drives recruitment keep their own spawning output
  spawners <- unique(unlist(lapply(rec, function(s) {
    names(.SPFromWeights(OM@Stock[[match(s, stock_names)]]@SRR@SPFrom, s, stock_names))
  })))
  for (s in setdiff(stock_names, rec)) {
    i <- match(s, stock_names)
    OM@Stock[[i]]@SRR@R0 <- 0
    if (!length(OM@Stock[[i]]@SRR@SPFrom) && !s %in% spawners) OM@Stock[[i]]@SRR@SPFrom <- source[[s]]
  }

  # The recruiting stocks are fitted with one stock-recruit relationship
  .SpawningSources(OM, rec, stock_names)
  if (length(rec) > 1L) {
    models <- vapply(rec, function(s) .SRRModelName(OM@Stock[[match(s, stock_names)]]@SRR@Model), character(1))
    if (length(unique(models)) > 1L) {
      cli::cli_abort("The recruiting stocks {.val {rec}} must have the same {.code SRR@Model}.")
    }
  }
  OM
}

# Share of the fitted unfished recruitment for each recruiting stock, by simulation, from the input SRR@R0
.R0Shares <- function(Stocks, rec, nSim) {
  R0 <- lapply(Stocks[rec], function(x) x@SRR@R0)
  given <- lengths(R0) > 0
  if (!any(given)) {
    R0 <- lapply(rec, function(s) rep(1, nSim))
  } else if (!all(given)) {
    cli::cli_abort("Set {.code SRR@R0} for every recruiting stock ({.val {rec[!given]}} {?has/have} none) to split recruitment among them.")
  } else {
    R0 <- lapply(R0, function(x) rep(as.numeric(x), length.out = nSim))
  }
  tot <- Reduce(`+`, R0)
  if (any(tot <= 0)) cli::cli_abort("{.code SRR@R0} of the recruiting stocks {.val {rec}} must sum to more than zero.")
  stats::setNames(lapply(R0, function(x) x / tot), rec)
}

# The fitting unit, named by the complex: member stocks with the recruiting stocks first
.FitUnit <- function(OM, Stocks = NULL) {
  stock_names <- MSEtool::StockNames(OM)
  if (is.null(Stocks)) Stocks <- stats::setNames(.ListOf(OM@Stock), stock_names)
  rec <- .RecruitingStocks(OM, stock_names)
  unit <- list(stocks = c(rec, setdiff(stock_names, rec)), rec = rec,
               share = .R0Shares(Stocks, rec, OM@nSim), spawn = .SpawningSources(OM, rec, stock_names),
               complex = 1L,
               combined = length(stock_names) > 1L)
  stats::setNames(list(unit), names(OM@Complexes)[1])
}

# Array of dims Sim x Age x Year from a scalar, Sim vector or collapsed array
.AgeArray <- function(x, dims) {
  if (is.null(x) || !length(x)) return(NULL)
  if (is.null(dim(x))) {
    if (length(x) %in% c(1L, dims[1])) return(array(x, dims))
    cli::cli_abort("Cannot expand a vector of length {length(x)} to Sim x Age x Year.")
  }
  out <- .BroadcastArray(x, dims)
  if (!identical(dim(out), as.integer(dims))) cli::cli_abort("Cannot expand an array of dimension {dim(x)} to {dims}.")
  out
}

# Spawning sources of a stock as named weights; an empty SPFrom is the stock itself
.SPFromWeights <- function(SPFrom, stock, stock_names) {
  if (is.null(SPFrom) || !length(SPFrom)) return(stats::setNames(1, stock))
  SPFrom <- .NameSPFrom(SPFrom, stock_names, stock)
  if (length(names(SPFrom)) && any(nzchar(names(SPFrom)))) return(stats::setNames(as.numeric(SPFrom), names(SPFrom)))
  stats::setNames(rep(1, length(SPFrom)), SPFrom)
}

.PlaceholderFleet <- function(OM, nAge) {
  fl <- MSEtool::Fleet(Name = "Placeholder")
  fl@Effort@Effort <- array(0, c(1, OM@nYear))
  fl@Catchability@Efficiency <- array(1, c(1, 1))
  fl@Selectivity@MeanAtAge <- array(1, c(1, nAge))
  fl@Selectivity@isAtLength <- FALSE
  fl
}

# One stock with biology averaged over the members at the unfished proportion of each stock by age
.CombineStocks <- function(OM, Stocks, unit, name) {
  members <- unit$stocks
  rec <- members[1]
  stock_names <- MSEtool::StockNames(OM)

  tmp <- OM
  tmp@Stock <- Stocks[members]
  for (s in members) tmp@Stock[[s]]@SRR@R0 <- if (s %in% unit$rec) unit$share[[s]] else 0
  tmp@Herm <- .HermPairNames(OM, stock_names)
  # Unfished numbers need stocks only; placeholder fleets let the om populate
  tmp@Fleet <- lapply(stats::setNames(members, members), function(s) list(Placeholder = .PlaceholderFleet(OM, length(Stocks[[s]]@Ages@Classes))))
  tmp@Complexes <- stats::setNames(list(seq_along(members)), name)
  tmp@Data <- list()
  tmp@Obs <- list()
  tmp@Imp <- list()
  N <- MSEtool::CalcUnfishedNumber(MSEtool::PopulateOM(tmp, silent = TRUE))
  if (is.null(names(N))) names(N) <- members

  stk <- Stocks[[rec]]
  nSim <- OM@nSim
  nYear <- OM@nYear
  nTot <- nYear + OM@pYear
  Ages <- stk@Ages@Classes
  dims <- c(nSim, length(Ages), nTot)
  dn <- list(Sim = seq_len(nSim), Age = Ages, Year = c(MSEtool::Years(OM, "Historical"), MSEtool::Years(OM, "Projection")))

  yr <- c(seq_len(nYear), rep(nYear, nTot - nYear))
  Ntot <- Reduce(`+`, N[members])
  p <- lapply(N[members], function(x) {
    out <- x / Ntot
    out[!is.finite(out)] <- 0
    .AgeArray(out, c(nSim, length(Ages), nYear))[, , yr, drop = FALSE]
  })
  empty <- .AgeArray(Ntot, c(nSim, length(Ages), nYear))[, , yr, drop = FALSE] <= 0
  p[[rec]][empty] <- 1

  get <- function(s, f) .AgeArray(f(Stocks[[s]]), dims)
  wsum <- function(f) Reduce(`+`, lapply(members, function(s) p[[s]] * f(s)))

  M <- -log(wsum(function(s) exp(-get(s, function(x) x@NaturalMortality@MeanAtAge))))
  W <- wsum(function(s) get(s, function(x) x@Weight@MeanAtAge))
  L <- wsum(function(s) get(s, function(x) x@Length@MeanAtAge))
  L2 <- wsum(function(s) {
    len <- get(s, function(x) x@Length@MeanAtAge)
    cv <- get(s, function(x) x@Length@CVatAge) %||% array(0.1, dims)
    (cv * len)^2 + len^2
  })
  CV <- sqrt(pmax(L2 - L^2, 0)) / L
  CV[!is.finite(CV)] <- 0

  # Spawning output of the population: the stocks in SRR@SPFrom of the recruiting stocks
  src <- unit$spawn %||% .SPFromWeights(stk@SRR@SPFrom, rec, stock_names)
  fec <- function(s) {
    out <- get(s, function(x) x@Fecundity@MeanAtAge)
    out %||% (get(s, function(x) x@Weight@MeanAtAge) * get(s, function(x) x@Maturity@MeanAtAge))
  }
  Fec <- Reduce(`+`, lapply(names(src), function(s) src[[s]] * p[[s]] * fec(s)))
  # Mature spawner biomass as a fraction of total biomass at age
  Mat <- Reduce(`+`, lapply(names(src), function(s) {
    src[[s]] * p[[s]] * get(s, function(x) x@Weight@MeanAtAge) * get(s, function(x) x@Maturity@MeanAtAge)
  })) / W
  Mat[!is.finite(Mat)] <- 0

  setdn <- function(x) {
    dimnames(x) <- dn
    x
  }
  stk@Name <- name
  stk@NaturalMortality@MeanAtAge <- setdn(M)
  stk@Weight@MeanAtAge <- setdn(W)
  stk@Length@MeanAtAge <- setdn(L)
  stk@Length@CVatAge <- setdn(CV)
  stk@Fecundity@MeanAtAge <- setdn(Fec)
  stk@Maturity@MeanAtAge <- setdn(pmin(Mat, 1))
  stk
}
