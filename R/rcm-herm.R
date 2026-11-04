
# Fitting units for RCM(om): single stocks, or hermaphroditic complexes fitted as one population

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

# Units named by complex (or stock): member stocks with the recruiting stock first, and complex index
.FitUnits <- function(OM) {
  stock_names <- MSEtool::StockNames(OM)
  herm <- .HermPairNames(OM, stock_names)
  cx <- OM@Complexes
  if (!length(cx)) {
    if (length(herm)) cli::cli_abort("Set {.code OM@Complexes} so that each pair of stocks in {.code OM@Herm} is one complex.")
    cx <- stats::setNames(as.list(stock_names), stock_names)
  }
  cx <- lapply(cx, function(x) if (is.character(x)) x else stock_names[x])
  if (is.null(names(cx))) names(cx) <- vapply(cx, paste, character(1), collapse = "-")

  units <- lapply(seq_along(cx), function(k) {
    members <- cx[[k]]
    if (length(members) == 1L) return(list(stocks = members, complex = k, herm = FALSE))
    pairs <- Filter(function(h) h@From %in% members || h@To %in% members, herm)
    linked <- unique(unlist(lapply(pairs, function(h) c(h@From, h@To))))
    if (!length(pairs) || !setequal(linked, members)) {
      cli::cli_abort(c("RCM fits single stocks or hermaphroditic complexes.",
                       "x" = "Complex {.val {names(cx)[k]}} has stocks that are not linked by {.code OM@Herm}."))
    }
    R0 <- vapply(members, function(s) {
      R0 <- OM@Stock[[match(s, stock_names)]]@SRR@R0
      if (length(R0)) sum(R0) else 1
    }, numeric(1))
    rec <- members[R0 > 0]
    if (length(rec) != 1L) {
      cli::cli_abort("Complex {.val {names(cx)[k]}} must have one stock with recruitment ({.code SRR@R0 > 0}); the other stocks need {.code SRR@R0 = 0}.")
    }
    list(stocks = c(rec, setdiff(members, rec)), complex = k, herm = TRUE)
  })
  stats::setNames(units, names(cx))
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

# Spawning sources of the recruiting stock as named weights
.SPFromWeights <- function(SPFrom, rec, stock_names) {
  if (is.null(SPFrom) || !length(SPFrom)) return(stats::setNames(1, rec))
  if (length(names(SPFrom)) && any(nzchar(names(SPFrom)))) return(SPFrom)
  src <- if (is.character(SPFrom)) SPFrom else stock_names[SPFrom]
  stats::setNames(rep(1, length(src)), src)
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
.CombineHermStock <- function(OM, Stocks, unit, name) {
  members <- unit$stocks
  rec <- members[1]
  stock_names <- MSEtool::StockNames(OM)

  tmp <- OM
  tmp@Stock <- Stocks[members]
  for (s in members) tmp@Stock[[s]]@SRR@R0 <- as.numeric(s == rec)
  tmp@Herm <- Filter(function(h) h@From %in% members, .HermPairNames(OM, stock_names))
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

  src <- .SPFromWeights(stk@SRR@SPFrom, rec, stock_names)
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
