
# Convert stock/fleet objects to a legacy OM (values in cpars) for RCM_int()

# Expand dimensions that PopulateStock() collapsed to length 1
.BroadcastArray <- function(arr, dims) {
  if (is.null(arr)) return(NULL)
  d <- dim(arr)
  if (is.null(d) || length(d) != length(dims)) return(arr)
  if (all(d == dims)) return(arr)
  idx <- lapply(seq_along(d), function(i) {
    if (d[i] == 1 && dims[i] > 1) rep(1L, dims[i]) else seq_len(d[i])
  })
  out <- do.call(`[`, c(list(arr), idx, list(drop = FALSE)))
  dim(out) <- dims
  out
}

.ToSimVec <- function(x, nSim) {
  if (is.null(x) || !length(x)) return(NULL)
  if (is.array(x) || is.matrix(x)) {
    d <- dim(x)
    if (length(d) && d[1] == nSim) {
      out <- apply(x, 1, function(v) v[1])
      return(as.numeric(out))
    }
    return(rep(mean(x, na.rm = TRUE), nSim))
  }
  if (length(x) == nSim) return(as.numeric(x))
  if (length(x) == 1) return(rep(as.numeric(x), nSim))
  rep(as.numeric(x[1]), nSim)
}

.ModelName <- function(Model) {
  if (is.character(Model) && length(Model)) return(Model[1])
  if (is.function(Model)) {
    if (identical(Model, MSEtool::Ricker)) return("Ricker")
    if (identical(Model, MSEtool::BevertonHolt)) return("BevertonHolt")
    if (identical(Model, MSEtool::HockeyStick)) return("HockeyStick")
    if (identical(Model, MSEtool::DoubleNormal)) return("DoubleNormal")
  }
  NA_character_
}
.SRRModelName <- .ModelName

# Legacy SRrel code: 1 Beverton-Holt, 2 Ricker, 3 hockey stick
.SRRelCode <- function(Model) {
  nm <- .SRRModelName(Model)
  if (is.na(nm)) return(1L)
  switch(nm, Ricker = 2L, HockeyStick = 3L, 1L)
}

# Steepness is Pars$hR for Ricker, Pars$h otherwise
.SRRSteepness <- function(SRR) {
  pars <- SRR@Pars
  if (!length(pars)) return(NULL)
  if (identical(.SRRModelName(SRR@Model), "Ricker")) return(pars$hR)
  pars$h
}

# Starting L5/LFS/Vmaxlen for one fleet's selectivity
.SingleFleetSelStart <- function(Fleet, nSim) {
  default <- list(L5 = rep(1, nSim), LFS = rep(2, nSim), Vmaxlen = rep(1, nSim))
  if (is.null(Fleet)) return(default)

  Sel <- Fleet@Selectivity
  pars <- Sel@Pars
  if (length(pars) && all(c("L5", "LFS", "Vmaxlen") %in% names(pars))) {
    return(list(
      L5 = .ToSimVec(pars$L5, nSim),
      LFS = .ToSimVec(pars$LFS, nSim),
      Vmaxlen = .ToSimVec(pars$Vmaxlen, nSim)
    ))
  }

  MAL <- Sel@MeanAtLength
  Classes <- Sel@Classes
  if (!is.null(MAL) && length(Classes)) {
    d <- dim(MAL)
    sel1 <- if (length(d) >= 3) MAL[1, , 1] else MAL[1, ]
    mids <- Classes
    vmax <- max(sel1, na.rm = TRUE)
    if (vmax <= 0) return(default)
    rel <- sel1 / vmax
    L5 <- tryCatch(stats::approx(rel, mids, xout = 0.05, ties = "ordered")$y, error = function(e) NA)
    LFS_idx <- which.max(rel)
    LFS <- mids[LFS_idx]
    Vmaxlen <- rel[length(rel)]
    if (is.na(L5)) L5 <- min(mids, na.rm = TRUE)
    return(list(L5 = rep(L5, nSim), LFS = rep(LFS, nSim), Vmaxlen = rep(Vmaxlen, nSim)))
  }

  default
}

# Starting selectivity for the first fleet, and vul_par (LFS, L5, Vmaxlen by fleet) for all fleets
.FleetSelStart <- function(Fleet, nSim) {
  fleet_list <- if (is.null(Fleet) || inherits(Fleet, "fleet")) list(Fleet) else Fleet
  per_fleet <- lapply(fleet_list, .SingleFleetSelStart, nSim = nSim)

  vul_par <- vapply(per_fleet, function(s) {
    c(LFS = s$LFS[1], L5 = s$L5[1], Vmaxlen = s$Vmaxlen[1])
  }, numeric(3))
  dimnames(vul_par) <- list(c("LFS", "L5", "Vmaxlen"), names(fleet_list))

  first <- per_fleet[[1]]
  list(L5 = first$L5, LFS = first$LFS, Vmaxlen = first$Vmaxlen, vul_par = vul_par)
}

.stock2OM <- function(Stock, Fleet = NULL, nSim = NULL, pYear = NULL, CurrentYear = NULL, seed = 1) {
  if (!inherits(Stock, "stock")) cli::cli_abort("{.arg Stock} must be a {.cls stock} object.")

  nSim <- nSim %||% Stock@nSim
  pYear <- pYear %||% Stock@pYear %||% 0
  CurrentYear <- CurrentYear %||% Stock@CurrentYear
  nyears <- Stock@nYear
  maxage <- Stock@Ages@MaxAge

  legacyStock <- methods::new("Stock")
  legacyStock@Name <- Stock@Name %||% "Stock"
  legacyStock@maxage <- maxage
  legacyStock@SRrel <- .SRRelCode(Stock@SRR@Model)
  legacyStock@Msd <- c(0, 0)
  legacyStock@Linfsd <- c(0, 0)
  legacyStock@Ksd <- c(0, 0)
  legacyStock@AC <- c(0, 0)
  legacyStock@D <- c(0, 0)
  legacyStock@LenCV <- c(0.1, 0.1)
  legacyStock@L50 <- c(0, 0)
  legacyStock@L50_95 <- c(0, 0)
  legacyStock@Fdisc <- c(0, 0)
  legacyStock@Size_area_1 <- legacyStock@Frac_area_1 <- legacyStock@Prob_staying <- c(0.5, 0.5)
  # Placeholders; values by simulation are in OM@cpars
  legacyStock@M <- c(0.2, 0.2)
  legacyStock@h <- c(0.999, 0.999)
  legacyStock@Linf <- c(1, 1)
  legacyStock@K <- c(1, 1)
  legacyStock@t0 <- c(0, 0)
  legacyStock@a <- 1
  legacyStock@b <- 1
  legacyStock@Perr <- c(0.4, 0.4)

  legacyFleet <- methods::new("Fleet")
  legacyFleet@CurrentYr <- CurrentYear
  legacyFleet@nyears <- nyears
  legacyFleet@EffYears <- c(1, nyears)
  legacyFleet@EffLower <- legacyFleet@EffUpper <- c(0, 1)
  legacyFleet@Esd <- legacyFleet@qinc <- legacyFleet@qcv <- c(0, 0)
  legacyFleet@DR <- c(0, 0)
  legacyFleet@Spat_targ <- c(1, 1)
  legacyFleet@isRel <- FALSE
  legacyFleet@MPA <- FALSE

  sel_start <- .FleetSelStart(Fleet, nSim)
  legacyFleet@L5 <- range(sel_start$L5)
  legacyFleet@LFS <- range(sel_start$LFS)
  legacyFleet@Vmaxlen <- range(sel_start$Vmaxlen)

  OM <- methods::new("OM", Stock = legacyStock, Fleet = legacyFleet,
                      Obs = MSEtool::Generic_Obs, Imp = MSEtool::Perfect_Imp,
                      nsim = nSim, proyears = pYear)
  OM@seed <- seed
  OM@CurrentYr <- CurrentYear
  OM@nyears <- nyears

  nAge <- length(Stock@Ages@Classes)
  nYearTot <- if (length(Stock@Years)) length(Stock@Years) else nyears + pYear
  target_dim <- c(nSim, nAge, nYearTot)

  OM@cpars$Len_age <- .BroadcastArray(Stock@Length@MeanAtAge, target_dim)
  OM@cpars$Wt_age <- .BroadcastArray(Stock@Weight@MeanAtAge, target_dim)
  OM@cpars$Mat_age <- .BroadcastArray(Stock@Maturity@MeanAtAge, target_dim)
  OM@cpars$M_ageArray <- .BroadcastArray(Stock@NaturalMortality@MeanAtAge, target_dim)
  if (length(Stock@Fecundity@MeanAtAge)) OM@cpars$Fec_age <- .BroadcastArray(Stock@Fecundity@MeanAtAge, target_dim)

  CV <- Stock@Length@CVatAge
  if (length(dim(CV)) == 3) {
    OM@cpars$LatASD <- .BroadcastArray(CV, target_dim) * OM@cpars$Len_age
  } else if (length(CV)) {
    OM@cpars$LenCV <- .ToSimVec(CV, nSim)
  }

  h_val <- .SRRSteepness(Stock@SRR)
  if (!is.null(h_val)) OM@cpars$h <- .ToSimVec(h_val, nSim)

  if (length(Stock@SRR@SD)) OM@cpars$Perr <- .ToSimVec(Stock@SRR@SD, nSim)
  if (length(Stock@SRR@AC)) OM@cpars$AC <- .ToSimVec(Stock@SRR@AC, nSim)

  R0val <- Stock@SRR@R0
  if (is.null(R0val) || !length(R0val)) R0val <- 1000
  OM@cpars$R0 <- .ToSimVec(R0val, nSim)

  if (length(Stock@Depletion@Final)) OM@cpars$D <- .ToSimVec(Stock@Depletion@Final, nSim)

  OM@cpars$L5 <- sel_start$L5
  OM@cpars$LFS <- sel_start$LFS
  OM@cpars$Vmaxlen <- sel_start$Vmaxlen

  OM
}
