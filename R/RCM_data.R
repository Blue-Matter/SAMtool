

vec_slot_fn <- function(x, Data, err = FALSE) {
  res <- slot(Data, x)
  if(length(res) && !all(is.na(res))) {
    return(res[1, ])
  } else {
    if(err) stop(paste0("Nothing found in Data@", x), call. = FALSE)
    return(NULL)
  }
}

matrix_slot_fn <- function(x, Data) {
  res <- slot(Data, x)
  if(length(res) && !all(is.na(res))) return(res[1, , ]) else return(NULL)
}

pull_Ind <- function(Data, maxage) {
  Ind_name <- c("Ind", "SpInd", "VInd")
  s_sel_Ind <- c("B", "SSB", 1)
  s_sel_AddIndType <- 1:3
  lapply_fn <- function(x) {
    Index <- vec_slot_fn(x = x, Data = Data)
    if(!is.null(Index)) {
      ICV <- vec_slot_fn(paste0("CV_", x), Data)
      if(is.null(ICV)) {
        I_sd <- rep(NA_real_, length(Index))
      } else {
        if(sum(!is.na(ICV)) == 1) I_sd <- rep(ICV[1], length(Index))
        I_sd <- sdconv(1, ICV)
      }
      s_sel <- s_sel_Ind[match(x, Ind_name)]
      slotname <- x
    } else {
      I_sd <- s_sel <- slotname <- NULL
    }
    return(list(Index = Index, I_sd = I_sd, s_sel = s_sel, slotname = slotname))
  }

  get_Ind <- lapply(Ind_name, lapply_fn)
  out <- list(Index = do.call(cbind, lapply(get_Ind, getElement, "Index")),
              I_sd = do.call(cbind, lapply(get_Ind, getElement, "I_sd")),
              s_sel = do.call(c, lapply(get_Ind, getElement, "s_sel")),
              slotname = do.call(c, lapply(get_Ind, getElement, "slotname")))
  if(!is.null(out$Index)) {
    out$V <- matrix(NA_real_, maxage + 1, ncol(out$Index))
    out$I_units <- rep(1, ncol(out$Index))
  } else {
    out$V <- out$I_units <- NULL
  }
  return(out)
}

pull_AddInd <- function(Data, maxage) {
  Ind_name <- c("Ind", "SpInd", "VInd")
  s_sel_Ind <- c("B", "SSB", 1)
  s_sel_AddIndType <- 1:3
  if(!all(is.na(Data@AddInd))) {
    nindex <- dim(Data@AddInd)[2]
    nyears <- dim(Data@AddInd)[3]
    Index <- Data@AddInd[1, , ] %>% matrix(nyears, nindex, byrow = TRUE)
    if(!all(is.na(Data@CV_AddInd[1, , ]))) {
      I_sd <- sdconv(1, Data@CV_AddInd[1, , ]) %>% matrix(nyears, nindex, byrow = TRUE)
    } else {
      I_sd <- array(NA_real_, dim(Index))
    }

    V <- matrix(NA_real_, maxage + 1, nindex)
    s_sel <- rep(NA_character_, nindex)
    for(i in 1:nindex) {
      if(!all(is.na(Data@AddIndV[1, i, ]))) {
        V[, i] <- Data@AddIndV[1, i, ]
        s_sel[i] <- "free"
      } else if(!is.na(Data@AddIndType[i])) {
        sel_arg <- match(Data@AddIndType[i], s_sel_AddIndType)
        if(is.na(sel_arg)) {
          message("Data@AddIndType[", i, "] is undefined.")
          s_sel[i] <- 1
          message("Selectivity for survey number ", i, " assumed to be a vulnerable biomass survey.")
        } else {
          s_sel[i] <- s_sel_Ind[sel_arg]
        }
      } else {
        s_sel[i] <- 1
        message("Selectivity for survey number ", i, " in Data@AddInd is undefined. Assuming to be a vulnerable biomass survey.")
      }
    }
    if(all(!is.na(Data@AddIunits)) && length(Data@AddIunits) == nindex) {
      I_units <- Data@AddIunits
    } else {
      I_units <- rep(1, nindex)
    }
    return(list(Index = Index, I_sd = I_sd, s_sel = s_sel, slotname = rep("AddInd", ncol(Index)), V = V,
                I_units = I_units))
  } else {
    return(list(Index = NULL, I_sd = NULL, s_sel = NULL, slotname = NULL, V = NULL, I_units = NULL))
  }
}


pull_Index <- function(Data, maxage) {
  Ind <- pull_Ind(Data, maxage)
  AddInd <- pull_AddInd(Data, maxage)

  if(all(is.na(Ind$I_sd)) && all(is.na(AddInd$I_sd))) {
    I_sd <- NULL
  } else {
    I_sd <- cbind(Ind$I_sd, AddInd$I_sd)
  }
  return(list(Index = cbind(Ind$Index, AddInd$Index), I_sd = I_sd, s_selectivity = c(Ind$s_sel, AddInd$s_sel),
              slotname = c(Ind$slotname, AddInd$slotname), V = cbind(Ind$V, AddInd$V),
              I_units = c(Ind$I_units, AddInd$I_units)))
}



int_s_sel <- function(s_selectivity, RCMdata) {
  if(is.null(s_selectivity)) return(-4)
  
  nfleet <- RCMdata@Misc$nfleet
  nsurvey <- length(s_selectivity)
  stopifnot(length(s_selectivity) == RCMdata@Misc$nsurvey)
  
  s_sel <- suppressWarnings(as.numeric(s_selectivity)) # Numbers match fleets, otherwise see next lines
  s_sel[s_selectivity == "B"] <- -4
  s_sel[s_selectivity == "SSB"] <- -3
  s_sel[s_selectivity == "free"] <- -2
  s_sel[s_selectivity == "logistic"] <- -1
  s_sel[s_selectivity == "dome"] <- 0

  if(any(s_sel > nfleet, na.rm = TRUE)) {
    stop(paste("There are undefined fishing fleets in s_selectivity (for indices). There are only", nfleet, "fleets."),
         call. = FALSE)
  }

  if(any(is.na(s_sel))) {
    stop("Character entries for s_selectivity (for indices) must be either: \"B\", \"SSB\", \"logistic\", \"dome\", or \"free\"", call. = FALSE)
  }
  
  message("\nIndex selectivity setup:")
  for(sur in 1:nsurvey) {
    if(s_sel[sur] > 0) {
      sout <- paste("fishery fleet", s_sel[sur])
    } else {
      sout <- switch(s_sel[sur] %>% as.character(),
                     "-4" = "total biomass",
                     "-3" = "spawning biomass",
                     "-2" = "individual parameters at age (free)",
                     "-1" = "logistic function",
                     "0" = "dome function")
    }
    message("Index ", sur, ": ", sout)
  }
  return(s_sel)
}


int_sel <- function(selectivity, RCMdata) {
  sel <- suppressWarnings(as.numeric(selectivity))
  sel[selectivity == "free"] <- -2
  sel[selectivity == "logistic"] <- -1
  sel[selectivity == "dome"] <- 0
  
  if(any(is.na(sel))) {
    stop("Character entries for selectivity (for fleets) must be either: \"logistic\", \"dome\", or \"free\"", call. = FALSE)
  }
  
  message("\nFishery selectivity setup:")
  Yr <- RCMdata@Misc$CurrentYr - RCMdata@Misc$nyears:1 + 1
  no_blocks <- apply(RCMdata@sel_block, 2, function(x) length(unique(x)) == 1) %>% all()
  for(bb in 1:length(sel)) {
    fout <- switch(sel[bb] %>% as.character(),
                   "-2" = "individual parameters at age (free)",
                   "-1" = "logistic function",
                   "0" = "dome function")
    if(no_blocks) {
      message("Fleet ", bb, ": ", fout)
    } else {
      fleet <- lapply(1:ncol(RCMdata@sel_block), function(ff) {
        y <- Yr[RCMdata@sel_block[, ff] == bb]
        if(length(y)) {
          if(all(diff(y) == 1)) {
            paste0(ff, " (", range(y) %>% paste(collapse = "-"), ")")
          } else {
            paste0(ff, " (", range(y) %>% paste(collapse = "-"), ", with gaps)")
          }
        } else {
          NULL
        }
      })
      message("Block ", bb, " (", fout, ") assigned to fishery:\n", do.call(c, fleet) %>% paste(collapse = "\n"))
    }
  }
  
  return(sel)
}


make_LWT <- function(LWT, nfleet, nsurvey) {

  if(is.null(LWT$Chist)) {
    LWT$Chist <- rep(1, nfleet)
  } else if(length(LWT$Chist) == 1 && nfleet > 1) {
    LWT$Chist <- rep(LWT$Chist, nfleet)
  }
  if(length(LWT$Chist) != nfleet) stop("LWT$Chist should be a vector of length ", nfleet, ".")

  if(is.null(LWT$Index)) {
    LWT$Index <- rep(1, max(1, nsurvey))
  } else if(length(LWT$Index) == 1 && nsurvey > 1) {
    LWT$Index <- rep(LWT$Index, nsurvey)
  }
  if(length(LWT$Index) != max(1, nsurvey)) stop("LWT$Index should be a vector of length ", nsurvey, ".")

  if(is.null(LWT$CAA)) {
    LWT$CAA <- rep(1, nfleet)
  } else if(length(LWT$CAA) == 1 && nfleet > 1) {
    LWT$CAA <- rep(LWT$CAA, nfleet)
  }
  if(length(LWT$CAA) != nfleet) stop("LWT$CAA should be a vector of length ", nfleet, ".")

  if(is.null(LWT$CAL)) {
    LWT$CAL <- rep(1, nfleet)
  } else if(length(LWT$CAL) == 1 && nfleet > 1) {
    LWT$CAL <- rep(LWT$CAL, nfleet)
  }
  if(length(LWT$CAL) != nfleet) stop("LWT$CAL should be a vector of length ", nfleet, ".")

  if(is.null(LWT$MS)) {
    LWT$MS <- rep(1, nfleet)
  } else if(length(LWT$MS) == 1 && nfleet > 1) {
    LWT$MS <- rep(LWT$MS, nfleet)
  }
  if(length(LWT$MS) != nfleet) stop("LWT$MS should be a vector of length ", nfleet, ".")

  if(is.null(LWT$C_eq)) {
    LWT$C_eq <- rep(1, max(1, nfleet))
  } else if(length(LWT$C_eq) == 1 && nfleet > 1) {
    LWT$C_eq <- rep(LWT$C_eq, nfleet)
  }
  if(length(LWT$C_eq) != nfleet) stop("LWT$C_eq should be a vector of length ", nfleet, ".")

  if(is.null(LWT$IAA)) {
    LWT$IAA <- rep(1, max(1, nsurvey))
  } else if(length(LWT$IAA) == 1 && nsurvey > 1) {
    LWT$IAA <- rep(LWT$IAA, nsurvey)
  }
  if(length(LWT$IAA) != max(1, nsurvey)) stop("LWT$IAA should be a vector of length ", nsurvey, ".")

  if(is.null(LWT$IAL)) {
    LWT$IAL <- rep(1, max(1, nsurvey))
  } else if(length(LWT$IAL) == 1 && nsurvey > 1) {
    LWT$IAL <- rep(LWT$IAL, nsurvey)
  }
  if(length(LWT$IAL) != max(1, nsurvey)) stop("LWT$IAL should be a vector of length ", nsurvey, ".")

  return(LWT)
}



#' @rdname RCM
#' @param RCMdata An \linkS4class{RCMdata} object.
#' @export
check_RCMdata <- function(RCMdata, OM, condition = c("catch", "catch2", "effort")) {

  message("\nChecking OM and data...\n")
  condition <- match.arg(condition)

  # Preliminary OM check for basics
  if(length(OM@nyears) == 0) stop("OM@nyears is needed.", call. = FALSE)
  if(length(OM@maxage) == 0) stop("OM@maxage is needed.", call. = FALSE)

  if(!length(RCMdata@Chist) && length(RCMdata@Ehist) && condition != "effort") {
    message("No catch found. Only effort found. Switching condition = \"effort\".")
    RCMdata@Misc$condition <- "effort"
  }

  if(condition == "catch" || condition == "catch2") {
    if(!length(RCMdata@Chist)) {
      stop("Full time series of catch is needed.", call. = FALSE)
    } else {
      if(any(is.na(RCMdata@Chist))) {
        stop("One or more of the historical annual catch observations is missing. Suggestion: use linear interpolation to fill these data.", call. = FALSE)
      }
      if(any(RCMdata@Chist < 0)) stop("All catch values should be zero or greater.", call. = FALSE)

      # Convert single fleet inputs to multiple fleet, e.g. matrices to arrays
      if(!is.matrix(RCMdata@Chist)) RCMdata@Chist <- matrix(RCMdata@Chist, ncol = 1)

      RCMdata@Misc$nyears <- nrow(RCMdata@Chist)
      RCMdata@Misc$nfleet <- ncol(RCMdata@Chist)
      RCMdata@Ehist <- matrix(0, RCMdata@Misc$nyears, RCMdata@Misc$nfleet)

      if(!length(RCMdata@Index) && !length(RCMdata@CAA) && !length(RCMdata@CAL) && 
         !length(RCMdata@MS) && !length(RCMdata@Ehist)) {
        message("No data other than Chist is provided. Model will switch to conditioning on equilibrium effort.")
        RCMdata@Misc$condition <- "effort"
        RCMdata@Ehist <- matrix(1, RCMdata@Misc$nyears, RCMdata@Misc$nfleet)
        RCMdata@E_eq <- rep(1, RCMdata@Misc$nfleet)
      } else {
        RCMdata@Misc$condition <- condition
      }
    }
  }
  
  if(condition == "effort") {
    RCMdata@Misc$condition <- "effort"
    if(!length(RCMdata@Ehist)) {
      stop("Full time series of effort is needed.")
    } else {
      if(any(is.na(RCMdata@Ehist))) stop("Effort time series is not complete (contains NA's)")
      if(any(RCMdata@Ehist < 0)) stop("All effort values should be positive.")

      if(!is.matrix(RCMdata@Ehist)) RCMdata@Ehist <- matrix(RCMdata@Ehist, ncol = 1)

      RCMdata@Misc$nyears <- nrow(RCMdata@Ehist)
      RCMdata@Misc$nfleet <- ncol(RCMdata@Ehist)

      if(length(RCMdata@Chist) && !is.matrix(RCMdata@Chist)) RCMdata@Chist <- matrix(RCMdata@Chist, ncol = 1)
      if(!length(RCMdata@Chist)) RCMdata@Chist <- matrix(0, RCMdata@Misc$nyears, RCMdata@Misc$nfleet)
    }
  }

  message("RCM model is conditioned on: ", RCMdata@Misc$condition)
  message(RCMdata@Misc$nfleet, " fleet(s) detected.")
  message(RCMdata@Misc$nyears, " years of data detected.")

  # Match number of historical years of catch/effort to OM
  if(OM@nyears != RCMdata@Misc$nyears) {
    cpars_cond <- length(OM@cpars) && any(vapply(OM@cpars, function(x) inherits(x, "matrix") || inherits(x, "array"), logical(1)))
    if(cpars_cond) {
      stmt <- paste0("OM@nyears != length(", ifelse(grepl("catch", RCMdata@Misc$condition), "Chist", "Ehist"), "). ",
                     "There will be indexing errors in your custom parameters (OM@cpars).")
      stop(stmt, call. = FALSE)
    } else {
      message("OM@nyears was updated to length(", ifelse(grepl("catch", RCMdata@Misc$condition), "Chist", "Ehist"), "): ", RCMdata@Misc$nyears)
      OM@nyears <- RCMdata@Misc$nyears
    }
  }
  if(!length(OM@CurrentYr)) OM@CurrentYr <- RCMdata@Misc$nyears
  
  # C_sd
  if(length(RCMdata@C_sd)) {
    if(is.vector(RCMdata@C_sd)) {
      if(length(RCMdata@C_sd) != RCMdata@Misc$nyears) stop("Length of C_sd vector does not equal nyears (", RCMdata@Misc$nyears, ").", call. = FALSE)
      RCMdata@C_sd <- matrix(RCMdata@C_sd, ncol = 1)
    } else if(is.matrix(RCMdata@C_sd)) {
      if(nrow(RCMdata@C_sd) != RCMdata@Misc$nyears) stop("Number of rows of C_sd matrix does not equal nyears (", RCMdata@Misc$nyears, "). NAs are acceptable.", call. = FALSE)
      if(ncol(RCMdata@C_sd) != RCMdata@Misc$nfleet) stop("Number of columns of C_sd matrix does not equal nfleet (", RCMdata@Misc$nfleet, ").", call. = FALSE)
    }
  } else {
    RCMdata@C_sd <- matrix(0.01, RCMdata@Misc$nyears, RCMdata@Misc$nfleet)
  }

  # Indices
  if(length(RCMdata@Index)) {
    if(is.vector(RCMdata@Index)) {
      if(length(RCMdata@Index) != RCMdata@Misc$nyears) stop("Length of Index vector does not equal nyears (", RCMdata@Misc$nyears, "). NAs are acceptable.", call. = FALSE)
      RCMdata@Index <- matrix(RCMdata@Index, ncol = 1)
    } else if(is.matrix(RCMdata@Index)) {
      if(nrow(RCMdata@Index) != RCMdata@Misc$nyears) stop("Number of rows of Index matrix does not equal nyears (", RCMdata@Misc$nyears, "). NAs are acceptable.", call. = FALSE)
    } else stop("Index is neither a vector nor a matrix.", call. = FALSE)

    RCMdata@Misc$nsurvey <- ncol(RCMdata@Index)
    
    if(length(RCMdata@I_sd)) {
      if(is.vector(RCMdata@I_sd)) {
        if(length(RCMdata@I_sd) != RCMdata@Misc$nyears) stop("Length of I_sd vector does not equal nyears (", RCMdata@Misc$nyears, ").", call. = FALSE)
        RCMdata@I_sd <- matrix(RCMdata@I_sd, ncol = 1)
      } else if(is.matrix(RCMdata@I_sd)) {
        if(nrow(RCMdata@I_sd) != RCMdata@Misc$nyears) stop("Number of rows of I_sd matrix does not equal nyears (", RCMdata@Misc$nyears, "). NAs are acceptable.", call. = FALSE)
        if(ncol(RCMdata@I_sd) != RCMdata@Misc$nsurvey) stop("Number of columns of I_sd matrix does not equal nsurvey (", RCMdata@Misc$nsurvey, ").", call. = FALSE)
      }
      
      SD_NA <- is.na(RCMdata@I_sd)
      if(sum(SD_NA)) {
        SD_out <- !is.na(RCMdata@Index[SD_NA])
        if(any(SD_out)) stop("There are NA's in data@I_sd for years associated with survey values in data@Index.", call. = FALSE)
      }
    }
  } else {
    RCMdata@Misc$nsurvey <- 0
    RCMdata@Index <- RCMdata@I_sd <- matrix(NA, ncol = 1, nrow = RCMdata@Misc$nyears)
  }
  message(RCMdata@Misc$nsurvey, " survey(s) detected.")

  # Process age comps
  if(length(RCMdata@CAA)) {

    if(is.matrix(RCMdata@CAA)) RCMdata@CAA <- array(RCMdata@CAA, c(dim(RCMdata@CAA), 1))

    if(dim(RCMdata@CAA)[1] != RCMdata@Misc$nyears) {
      stop("Number of CAA rows (", dim(RCMdata@CAA)[1], ") does not equal nyears (", RCMdata@Misc$nyears, "). NAs are acceptable.", call. = FALSE)
    }
    if(dim(RCMdata@CAA)[2] < OM@maxage + 1) {
      message("Number of CAA columns (", dim(RCMdata@CAA)[2], ") does not equal OM@maxage + 1 (", OM@maxage + 1, ").")
      message("Assuming no observations for ages greater than 0 - ", dim(RCMdata@CAA)[2] - 1, " and filling with zeros.")
      add_ages <- OM@maxage + 1 - dim(RCMdata@CAA)[2]
      CAA_new <- array(0, c(RCMdata@Misc$nyears, OM@maxage + 1, RCMdata@Misc$nfleet))
      CAA_new[, 1:dim(RCMdata@CAA)[2], ] <- RCMdata@CAA
      RCMdata@CAA <- CAA_new
    }
    if(dim(RCMdata@CAA)[2] > OM@maxage + 1) {
      OM@maxage <- dim(RCMdata@CAA)[2] - 1
      message("Increasing OM@maxage to ", OM@maxage, ".")
    }
    if(dim(RCMdata@CAA)[3] != RCMdata@Misc$nfleet) {
      stop("Number of CAA slices (", dim(RCMdata@CAA)[3], ") does not equal nfleet (", RCMdata@Misc$nfleet, "). NAs are acceptable.", call. = FALSE)
    }
    message("Fleet age comps (CAA) processed, assuming ages 0 - ", OM@maxage, " in array.")
    
    if(!length(RCMdata@CAA_ESS)) {
      RCMdata@CAA_ESS <- apply(RCMdata@CAA, c(1, 3), sum, na.rm = TRUE)
    }
    if(is.vector(RCMdata@CAA_ESS)) {
      if(length(RCMdata@CAA_ESS) != RCMdata@Misc$nyears) stop("Length of CAA_ESS vector does not equal nyears (", RCMdata@Misc$nyears, "). NAs are acceptable.", call. = FALSE)
      RCMdata@CAA_ESS <- matrix(RCMdata@CAA_ESS, ncol = 1)
    } else if(is.matrix(RCMdata@CAA_ESS)) {
      if(nrow(RCMdata@CAA_ESS) != RCMdata@Misc$nyears) stop("Number of rows of CAA_ESS matrix does not equal nyears (", RCMdata@Misc$nyears, "). NAs are acceptable.", call. = FALSE)
      if(ncol(RCMdata@CAA_ESS) != RCMdata@Misc$nfleet) stop("Number of columns of CAA_ESS matrix does not equal nfleet (", RCMdata@Misc$nfleet, "). NAs are acceptable.", call. = FALSE)
    } else stop("CAA_ESS is neither a vector nor a matrix.", call. = FALSE)
    
    # Check if CAA_ESS > 0 if there are no data
    RCMdata@CAA_ESS[apply(RCMdata@CAA, c(1, 3), sum, na.rm = TRUE) == 0] <- 0
    
  } else {
    RCMdata@CAA <- array(0, c(RCMdata@Misc$nyears, OM@maxage + 1, RCMdata@Misc$nfleet))
    RCMdata@CAA_ESS <- matrix(0, RCMdata@Misc$nyears, RCMdata@Misc$nfleet)
  }
  RCMdata@CAA <- apply(RCMdata@CAA, c(1, 3), find_na) %>% aperm(c(2, 1, 3))

  OM_samp <- check_OM_for_sampling(OM, RCMdata) # Sample life history, selectivity, and obs parameters

  set.seed(OM@seed)
  message("Getting biological parameters from OM...")
  suppressMessages({
    StockPars <- MSEtool::SampleStockPars(OM_samp, msg = FALSE)
    ObsPars <- MSEtool::SampleObsPars(OM_samp)
    FleetPars <- MSEtool::SampleFleetPars(OM_samp, msg = FALSE)
  })

  # Process length comps
  if(length(RCMdata@CAL)) {
    if(is.matrix(RCMdata@CAL)) RCMdata@CAL <- array(RCMdata@CAL, c(dim(RCMdata@CAL), 1))
    
    if(dim(RCMdata@CAL)[1] != RCMdata@Misc$nyears) {
      stop("Number of CAL rows (", dim(RCMdata@CAL)[1], ") does not equal nyears (", RCMdata@Misc$nyears, "). NAs are acceptable.", call. = FALSE)
    }
    message(dim(RCMdata@CAL)[2], " length bins detected in CAL.")
    if(dim(RCMdata@CAL)[3] != RCMdata@Misc$nfleet) {
      stop("Number of CAL slices (", dim(RCMdata@CAA)[3], ") does not equal nfleet (", RCMdata@Misc$nfleet, "). NAs are acceptable.", call. = FALSE)
    }
    if(!length(RCMdata@CAL_ESS)) {
      RCMdata@CAL_ESS <- apply(RCMdata@CAL, c(1, 3), sum, na.rm = TRUE)
    }
    if(is.vector(RCMdata@CAL_ESS)) {
      if(length(RCMdata@CAL_ESS) != RCMdata@Misc$nyears) stop("Length of CAL_ESS vector does not equal nyears (", RCMdata@Misc$nyears, "). NAs are acceptable.", call. = FALSE)
      RCMdata@CAL_ESS <- matrix(RCMdata@CAL_ESS, ncol = 1)
    } else if(is.matrix(RCMdata@CAL_ESS)) {
      if(nrow(RCMdata@CAL_ESS) != RCMdata@Misc$nyears) stop("Number of rows of CAL_ESS matrix does not equal nyears (", RCMdata@Misc$nyears, "). NAs are acceptable.", call. = FALSE)
      if(ncol(RCMdata@CAL_ESS) != RCMdata@Misc$nfleet) stop("Number of columns of CAL_ESS matrix does not equal nfleet (", RCMdata@Misc$nfleet, "). NAs are acceptable.", call. = FALSE)
    } else stop("CAL_ESS is neither a vector nor a matrix.", call. = FALSE)
    
    # Check if CAL_ESS > 0 if there are no data
    RCMdata@CAL_ESS[apply(RCMdata@CAL, c(1, 3), sum, na.rm = TRUE) == 0] <- 0
  } else {
    RCMdata@CAL_ESS <- matrix(0, RCMdata@Misc$nyears, RCMdata@Misc$nfleet)
  }

  # Process mean size
  if(length(RCMdata@MS)) {
    if(!length(RCMdata@MS_type) || !nchar(RCMdata@MS_type)) {
      message("Mean size (RCMdata@MS) found, but not type (RCMdata@MS_type). Assuming it's mean length.")
      RCMdata@MS_type <- "length"
    } else {
      RCMdata@MS_type <- match.arg(RCMdata@MS_type, choices = c("length", "weight"))
      message("Mean ", RCMdata@MS_type, " data found.")
    }
    if(is.vector(RCMdata@MS)) {
      if(length(RCMdata@MS) != RCMdata@Misc$nyears) stop("Mean size vector (MS) must be of length ", RCMdata@Misc$nyears, ".", call. = FALSE)
      RCMdata@MS <- matrix(RCMdata@MS, ncol = 1)
    }
    if(nrow(RCMdata@MS) != RCMdata@Misc$nyears) stop("Number of MS rows (", nrow(RCMdata@MS), ") does not equal nyears (", RCMdata@Misc$nyears, "). NAs are acceptable.", call. = FALSE)
    if(ncol(RCMdata@MS) != RCMdata@Misc$nfleet) stop("Number of MS columns (", ncol(RCMdata@MS), ") does not equal nfleet (", RCMdata@Misc$nfleet, "). NAs are acceptable.", call. = FALSE)

    if(!length(RCMdata@MS_cv)) {
      RCMdata@MS_cv <- rep(0.2, RCMdata@Misc$nfleet)
    } else if(length(RCMdata@MS_cv) == 1) {
      RCMdata@MS_cv <- rep(RCMdata@MS_cv, RCMdata@Misc$nfleet)
    }
    if(length(RCMdata@MS_cv) != RCMdata@Misc$nfleet) stop("Mean size CV vector (MS_cv) must be of length ", RCMdata@Misc$nfleet, ".", call. = FALSE)
  } else {
    RCMdata@MS <- matrix(NA, nrow = RCMdata@Misc$nyears, ncol = RCMdata@Misc$nfleet)
    RCMdata@MS_cv <- rep(0.2, RCMdata@Misc$nfleet)
    RCMdata@MS_type <- "length"
  }

  # Process equilibrium catch/effort - C_eq
  if(!length(RCMdata@C_eq)) RCMdata@C_eq <- rep(0, RCMdata@Misc$nfleet)
  if(grepl("catch", RCMdata@Misc$condition)) {
    if(length(RCMdata@C_eq) == 1) RCMdata@C_eq <- rep(RCMdata@C_eq, RCMdata@Misc$nfleet)
    if(length(RCMdata@C_eq) < RCMdata@Misc$nfleet) stop("C_eq needs to be of length nfleet (", RCMdata@Misc$nfleet, ").", call. = FALSE)
  }
  
  if(!length(RCMdata@C_eq_sd)) {
    RCMdata@C_eq_sd <- rep(0.01, RCMdata@Misc$nfleet)
  } else if(length(RCMdata@C_eq_sd) == 1) {
    RCMdata@C_eq_sd <- rep(RCMdata@C_eq_sd, RCMdata@Misc$nfleet)
  }
  if(length(RCMdata@C_eq_sd) != RCMdata@Misc$nfleet) stop("C_eq_sd needs to be of length nfleet (", RCMdata@Misc$nfleet, ").", call. = FALSE)
  
  if(RCMdata@Misc$condition != "effort" && any(RCMdata@C_eq > 0)) {
    message("Equilibrium catch was detected. The corresponding equilibrium F will be estimated.")
  }

  if(!length(RCMdata@E_eq)) RCMdata@E_eq <- rep(0, RCMdata@Misc$nfleet)
  if(RCMdata@Misc$condition == "effort") {
    if(length(RCMdata@E_eq) == 1) RCMdata@E_eq <- rep(RCMdata@E_eq, RCMdata@Misc$nfleet)
    if(length(RCMdata@E_eq) < RCMdata@Misc$nfleet) stop("E_eq needs to be of length nfleet (", RCMdata@Misc$nfleet, ").", call. = FALSE)
    if(any(RCMdata@E_eq > 0)) {
      message("Equilibrium effort was detected. The corresponding equilibrium F will be estimated.")
    }
  }

  # Process survey age comps
  if(length(RCMdata@IAA)) {
    if(is.matrix(RCMdata@IAA)) RCMdata@IAA <- array(RCMdata@IAA, c(dim(RCMdata@IAA), 1))
    if(dim(RCMdata@IAA)[1] != RCMdata@Misc$nyears) {
      stop("Number of IAA rows (", dim(RCMdata@IAA)[1], ") does not equal nyears (", RCMdata@Misc$nyears, "). NAs are acceptable.", call. = FALSE)
    }
    if(dim(RCMdata@IAA)[2] < OM@maxage + 1) {
      message("Number of IAA columns (", dim(RCMdata@IAA)[2], ") does not equal OM@maxage + 1 (", OM@maxage + 1, ").")
      message("Assuming no observations for ages greater than 0 - ", dim(RCMdata@IAA)[2] - 1, " and filling with zeros.")
      add_ages <- OM@maxage + 1 - dim(RCMdata@IAA)[2]
      IAA_new <- array(0, c(RCMdata@Misc$nyears, OM@maxage, RCMdata@Misc$nsurvey))
      IAA_new[, 1:dim(RCMdata@IAA)[2], ] <- RCMdata@IAA
      RCMdata@IAA <- IAA_new
    }
    if(dim(RCMdata@IAA)[2] > OM@maxage + 1) {
      stop("Error in age dimension of IAA.", call. = FALSE)
    }
    if(dim(RCMdata@IAA)[3] != RCMdata@Misc$nsurvey) {
      stop("Number of CAA slices (", dim(RCMdata@IAA)[3], ") does not equal nsurvey (", RCMdata@Misc$nsurvey, "). NAs are acceptable.", call. = FALSE)
    }
    message("Index age comps (IAA) processed, assuming ages 0 - ", OM@maxage, " in array.")
    
    if(!length(RCMdata@IAA_ESS)) {
      RCMdata@IAA_ESS <- apply(RCMdata@IAA, c(1, 3), sum, na.rm = TRUE)
    }
    if(is.vector(RCMdata@IAA_ESS)) {
      if(length(RCMdata@IAA_ESS) != RCMdata@Misc$nyears) stop("Length of IAA_ESS vector does not equal nyears (", RCMdata@Misc$nyears, "). NAs are acceptable.", call. = FALSE)
      RCMdata@IAA_ESS <- matrix(RCMdata@IAA_ESS, ncol = 1)
    } else if(is.matrix(RCMdata@IAA_ESS)) {
      if(nrow(RCMdata@IAA_ESS) != RCMdata@Misc$nyears) stop("Number of rows of IAA_ESS matrix does not equal nyears (", RCMdata@Misc$nyears, "). NAs are acceptable.", call. = FALSE)
      if(ncol(RCMdata@IAA_ESS) != RCMdata@Misc$nsurvey) stop("Number of columns of IAA_ESS matrix does not equal nsurvey (", RCMdata@Misc$nsurvey, "). NAs are acceptable.", call. = FALSE)
    } else stop("IAA_ESS is neither a vector nor a matrix.", call. = FALSE)
    
    # Check if IAA_ESS > 0 if there are no data
    RCMdata@IAA_ESS[apply(RCMdata@IAA, c(1, 3), sum, na.rm = TRUE) == 0] <- 0
  } else {
    RCMdata@IAA <- array(0, c(RCMdata@Misc$nyears, OM@maxage + 1, ncol(RCMdata@Index)))
    RCMdata@IAA_ESS <- array(0, dim(RCMdata@Index))
  }
  RCMdata@IAA <- apply(RCMdata@IAA, c(1, 3), find_na) %>% aperm(c(2, 1, 3))

  # Process survey length comps
  if(length(RCMdata@IAL)) {
    if(is.matrix(RCMdata@IAL)) RCMdata@IAL <- array(RCMdata@IAL, c(dim(RCMdata@IAL), 1))
    if(dim(RCMdata@IAL)[1] != RCMdata@Misc$nyears) {
      stop("Number of IAL rows (", dim(RCMdata@IAL)[1], ") does not equal nyears (", RCMdata@Misc$nyears, "). NAs are acceptable.", call. = FALSE)
    }
    message(dim(RCMdata@IAL)[2], " length bins detected in CAL.")
    if(dim(RCMdata@IAL)[3] != RCMdata@Misc$nsurvey) {
      stop("Number of IAL slices (", dim(RCMdata@IAL)[3], ") does not equal nsurvey (", RCMdata@Misc$nsurvey, "). NAs are acceptable.", call. = FALSE)
    }
    if(!length(RCMdata@IAL_ESS)) {
      RCMdata@IAL_ESS <- apply(RCMdata@IAL, c(1, 3), sum, na.rm = TRUE)
    }
    if(is.vector(RCMdata@IAL_ESS)) {
      if(length(RCMdata@IAL_ESS) != RCMdata@Misc$nyears) stop("Length of IAL_ESS vector does not equal nyears (", RCMdata@Misc$nyears, "). NAs are acceptable.", call. = FALSE)
      RCMdata@IAL_ESS <- matrix(RCMdata@IAL_ESS, ncol = 1)
    } else if(is.matrix(RCMdata@IAL_ESS)) {
      if(nrow(RCMdata@IAL_ESS) != RCMdata@Misc$nyears) stop("Number of rows of IAL_ESS matrix does not equal nyears (", RCMdata@Misc$nyears, "). NAs are acceptable.", call. = FALSE)
      if(ncol(RCMdata@IAL_ESS) != RCMdata@Misc$nsurvey) stop("Number of columns of IAL_ESS matrix does not equal nfleet (", RCMdata@Misc$nsurvey, "). NAs are acceptable.", call. = FALSE)
    } else stop("IAL_ESS is neither a vector nor a matrix.", call. = FALSE)
    
    # Check if IAL_ESS > 0 if there are no data
    RCMdata@IAL_ESS[apply(RCMdata@IAL, c(1, 3), sum, na.rm = TRUE) == 0] <- 0
  } else {
    RCMdata@IAL_ESS <- array(0, dim(RCMdata@Index))
  }
  
  # Length bin
  if(!sum(RCMdata@CAL_ESS) && !sum(RCMdata@IAL_ESS)) { # No length/index data
    RCMdata@Misc$lbin <- StockPars$CAL_bins
    RCMdata@Misc$lbinmid <- StockPars$CAL_binsmid
    RCMdata@Misc$nlbin <- length(RCMdata@Misc$lbinmid)
    RCMdata@CAL <- array(0, c(RCMdata@Misc$nyears, RCMdata@Misc$nlbin, RCMdata@Misc$nfleet))
    RCMdata@IAL <- array(0, c(RCMdata@Misc$nyears, RCMdata@Misc$nlbin, ncol(RCMdata@Index)))
  } else {
    if(sum(RCMdata@CAL_ESS) && !sum(RCMdata@IAL_ESS)) { # CAL only
      RCMdata@Misc$nlbin <- dim(RCMdata@CAL)[2]
      RCMdata@IAL <- array(0, c(RCMdata@Misc$nyears, RCMdata@Misc$nlbin, ncol(RCMdata@Index)))
    } else if(!sum(RCMdata@CAL_ESS) && sum(RCMdata@IAL_ESS)) { # IAL only
      RCMdata@Misc$nlbin <- dim(RCMdata@IAL)[2]
      RCMdata@CAL <- array(0, c(RCMdata@Misc$nyears, RCMdata@Misc$nlbin, RCMdata@Misc$nfleet))
    } else if(dim(RCMdata@CAL)[2] == dim(RCMdata@IAL)[2]) { # Both CAL/IAL. Check nlbin for both
      RCMdata@Misc$nlbin <- dim(RCMdata@CAL)[2]
    } else {
      stop("Number of length bins in CAL is not equal to those in IAL.", call. = FALSE)
    }
    
    if(!length(RCMdata@length_bin)) { # No length bins
      stop("You must specify length_bin for your length composition.", call. = FALSE)
      
    } else if(length(RCMdata@length_bin) == RCMdata@Misc$nlbin) { # Even length bins
      binWidth <- unique(diff(RCMdata@length_bin))
      if(length(binWidth) == 1) {
        RCMdata@Misc$lbinmid <- RCMdata@length_bin
        RCMdata@Misc$lbin <- c(RCMdata@Misc$lbinmid - 0.5 * binWidth, max(RCMdata@Misc$lbinmid) + 0.5 * binWidth)
      } else {
        stop("Uneven length bins detected. Provide a vector of length n_bin + 1 of the boundaries of all length bins", call. = FALSE)
      }
      
    } else if(length(RCMdata@length_bin) == RCMdata@Misc$nlbin + 1)  { # Uneven length bins
      RCMdata@Misc$lbin <- RCMdata@length_bin 
      RCMdata@Misc$lbinmid <- 0.5 * (RCMdata@Misc$lbin[1:RCMdata@Misc$nlbin + 1] - RCMdata@Misc$lbin[1:RCMdata@Misc$nlbin])
    } else {
      stop("Check vector of length_bin vs. the dimensions of the length compositions.", call. = FALSE)
    }
    RCMdata@CAL <- apply(RCMdata@CAL, c(1, 3), find_na) %>% aperm(c(2, 1, 3))
    RCMdata@IAL <- apply(RCMdata@IAL, c(1, 3), find_na) %>% aperm(c(2, 1, 3))
  }

  # Absolute survey
  if(RCMdata@Misc$nsurvey > 0) {
    if(!length(RCMdata@abs_I)) RCMdata@abs_I <- rep(0L, RCMdata@Misc$nsurvey)
    if(length(RCMdata@abs_I) < RCMdata@Misc$nsurvey) stop("abs_I should be of length", RCMdata@Misc$nsurvey, call. = FALSE)
  } else {
    RCMdata@abs_I <- 0L
  }

  # Index units - biomass/abundance
  if(RCMdata@Misc$nsurvey > 0) {
    if(!length(RCMdata@I_units)) RCMdata@I_units <- rep(1L, RCMdata@Misc$nsurvey)
    if(length(RCMdata@I_units) < RCMdata@Misc$nsurvey) stop("I_basis should be of length", RCMdata@Misc$nsurvey, call. = FALSE)
  } else {
    RCMdata@I_units <- 1L
  }

  # Ageing error
  if(!length(RCMdata@age_error)) RCMdata@age_error <- diag(OM@maxage + 1)
  if(any(dim(RCMdata@age_error) != OM@maxage + 1)) stop("age_error should be a square matrix of OM@maxage + 1 rows and columns", call. = FALSE)

  # Sel_block dummy fleets
  if(!length(RCMdata@sel_block)) {
    RCMdata@sel_block <- matrix(1:RCMdata@Misc$nfleet, nrow = RCMdata@Misc$nyears, ncol = RCMdata@Misc$nfleet, byrow = TRUE)
  } else {
    if(nrow(RCMdata@sel_block) != RCMdata@Misc$nyears) {
      stop(paste("sel_block should be a matrix of", RCMdata@Misc$nyears, "rows."), call. = FALSE)
    }
    if(ncol(RCMdata@sel_block) != RCMdata@Misc$nfleet) {
      stop(paste("sel_block should be a matrix of", RCMdata@Misc$nfleet, "columns."), call. = FALSE)
    }
  }
  RCMdata@Misc$nsel_block <- as.numeric(RCMdata@sel_block) %>% unique() %>% length()
  RCMdata@Misc$CurrentYr <- OM@CurrentYr

  return(list(RCMdata = RCMdata, OM = OM, StockPars = StockPars, ObsPars = ObsPars, FleetPars = FleetPars))
}

check_OM_for_sampling <- function(OM, RCMdata) {
  if(length(OM@nsim) == 0) stop("OM@nsim is needed.", call. = FALSE)
  if(length(OM@proyears) == 0) stop("OM@proyears is needed.", call. = FALSE)
  if(length(OM@seed) == 0) stop("OM@seed is needed.", call. = FALSE)

  cpars <- OM@cpars

  ###### Stock parameters
  # Len_at_age
  len_check <- !is.null(cpars$Len_age)
  if(!len_check) {
    Linf_check <- length(OM@Linf) == 2 || !is.null(cpars$Linf) || !is.null(cpars$Linfarray)
    K_check <- length(OM@K) == 2 || !is.null(cpars$K) || !is.null(cpars$Karray)
    t0_check <- length(OM@t0) == 2 || !is.null(cpars$t0)

    if(!Linf_check && !K_check && !t0_check) stop("Length-at-age not found in OM.", call. = FALSE)
  }
  # LenCV
  LenCV_check <- length(OM@LenCV) == 2 || !is.null(cpars$LenCV) || !is.null(cpars$LatASD)
  if(!LenCV_check) {
    any_CAL <- !is.null(RCMdata@CAL) && any(RCMdata@CAL > 0, na.rm = TRUE)
    any_ML <- RCMdata@MS_type == "length" && any(RCMdata@MS > 0, na.rm = TRUE)
    any_IAL <- !is.null(RCMdata@IAL) && any(RCMdata@IAL > 0, na.rm = TRUE)
    if(any_CAL || any_ML || any_IAL) {
      stop("OM@LenCV not found in OM.", call. = FALSE)
    } else {
      OM@LenCV <- c(0.1, 0.1)
    }
  }

  # Weight_at_age
  wt_check <- !is.null(cpars$Wt_age)
  if(!wt_check) {
    wt_check2 <- length(OM@a) > 0 && length(OM@b) > 0
    if(!wt_check2) stop("Weight-at-age not found in OM.", call. = FALSE)
  }

  # Maturity
  mat_check <- !is.null(cpars$Mat_age) ||
    (!is.null(cpars$ageM) & !is.null(cpars$age95)) ||
    (!is.null(cpars$L50) & (!is.null(cpars$L95) | !is.null(cpars$L50_95)))
  if(!mat_check) {
    mat_check2 <- length(OM@L50) == 2 && length(OM@L50_95) == 2
    if(!mat_check2) stop("Maturity-at-age not found in OM.", call. = FALSE)
  }

  # Natural mortality
  M_check <- !is.null(cpars$M_at_Length) || !is.null(cpars$Mage) || !is.null(cpars$M) || !is.null(cpars$Marray) || !is.null(cpars$M_ageArray)
  if(!M_check) {
    M_check2 <- length(OM@M) >= 2
    if(!M_check2) stop("Natural mortality not found in OM.", call. = FALSE)
  }

  # Msd - placeholder
  Msd_check <- length(OM@Msd) == 2 || !is.null(cpars$Msd)
  if(!Msd_check) OM@Msd <- c(0, 0)

  # Steepness
  h_check <- length(OM@h) == 2 || !is.null(cpars$h)
  if(!h_check) stop("Steepness (OM@h) not found.", call. = FALSE)

  # procsd
  procsd_check <- !is.null(cpars$Perr) || length(OM@Perr) == 2
  if(!procsd_check) stop("OM@Perr not found.", call. = FALSE)

  # autocorrelation - placeholder
  AC_check <- length(OM@AC) == 2 || !is.null(cpars$AC)
  if(!AC_check) OM@AC <- c(0, 0)

  # Depletion - placeholder
  D_check <- length(OM@D) == 2 || !is.null(cpars$D)
  if(!D_check) OM@D <- c(0, 0)

  # Stock recruit relationship
  SR_check <- length(OM@SRrel) == 1
  if(!SR_check) stop("Stock-recruit relationship (OM@SRrel) not found.", call. = FALSE)

  # R0 check
  R0_check <- length(OM@R0) > 0 || !is.null(cpars$R0)
  if(!R0_check) {
    if(RCMdata@Misc$condition == "effort") {
      OM@R0 <- 1
    } else {
      OM@R0 <- 1e3
      message("OM@R0 is used as the starting value for R0 in RCM, but was not found. By default, using 1000.")
    }
  }

  # Mvt parameters
  OM@Size_area_1 <- OM@Frac_area_1 <- OM@Prob_staying <- c(0.5, 0.5)

  # Fdisc placeholder
  Fdisc_check <- length(OM@Fdisc) == 2 || !is.null(cpars$Fdisc)
  if(!Fdisc_check) OM@Fdisc <- c(0, 0)

  ###### Fleet Parameters
  sel_check <- (length(OM@L5) == 2 | !is.null(cpars$L5)) &&
    (length(OM@LFS) == 2 | !is.null(cpars$LFS)) &&
    (length(OM@Vmaxlen) == 2 | !is.null(cpars$Vmaxlen))
  if(!sel_check) {
    stop("Selectivity parameters (OM@L5, OM@LFS, OM@Vmaxlen) not found. These are starting values for selectivity in the model", call. = FALSE)
  }

  # More placeholders
  OM@Esd <- OM@qinc <- OM@qcv <- OM@DR <- c(0, 0)
  OM@Spat_targ <- c(1, 1)

  OM@EffYears <- c(1, OM@nyears)
  OM@EffLower <- OM@EffUpper <- c(0, 1)

  ###### Observation Parameters - Iobs
  if(any(RCMdata@Index > 0, na.rm = TRUE) && !any(RCMdata@I_sd > 0, na.rm = TRUE)) {
    Isd_check <- length(OM@Iobs) == 2 || !is.null(cpars$Iobs)
    if(!Isd_check) stop("OM@Iobs is needed.", call. = FALSE)
  }
  Iobs <- OM@Iobs
  OM <- MSEtool::Replace(OM, MSEtool::Generic_Obs, silent = TRUE)
  if(length(Iobs) == 2) OM@Iobs <- Iobs

  ###### Imp
  OM <- MSEtool::Replace(OM, MSEtool::Perfect_Imp, silent = TRUE)

  return(OM)
}
