

vec_slot_fn <- function(x, Data, err = FALSE) {
  res <- slot(Data, x)
  if(!all(is.na(res))) {
    return(res[1, ])
  } else {
    if(err) stop(paste0("Nothing found in Data@", x), call. = FALSE)
    return(NULL)
  }
}

matrix_slot_fn <- function(x, Data) {
  res <- slot(Data, x)
  if(!all(is.na(res))) return(res[1, , ]) else return(NULL)
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
    out$V <- matrix(NA_real_, maxage, ncol(out$Index))
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

    V <- matrix(NA_real_, maxage, nindex)
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


SRA_tiny_comp <- function(x) {
  all_zero <- all(is.na(x)) | sum(x, na.rm = TRUE) == 0
  if(!all_zero) {
    ind <- is.na(x) | x == 0
    if(any(ind)) x[ind] <- 1e-8
  }
  return(x)
}

int_s_sel <- function(s_selectivity, nfleet = 1) {
  s_sel <- suppressWarnings(as.numeric(s_selectivity)) # Numbers match fleets, otherwise see next lines
  s_sel[s_selectivity == "B"] <- -4
  s_sel[s_selectivity == "SSB"] <- -3
  s_sel[s_selectivity == "free"] <- -2
  s_sel[s_selectivity == "logistic"] <- -1
  s_sel[s_selectivity == "dome"] <- 0

  if(any(s_sel > nfleet, na.rm = TRUE)) {
    stop(paste("There are undefined fishing fleets in s_selectivity (for surveys). There are only", nfleet, "fleets."),
         call. = FALSE)
  }

  if(any(is.na(s_sel))) {
    stop("Character entries for s_selectivity (for surveys) must be either: \"B\", \"SSB\", \"logistic\", \"dome\", or \"free\"", call. = FALSE)
  }

  return(s_sel)
}

int_sel <- function(selectivity) {
  sel <- suppressWarnings(as.numeric(selectivity))
  sel[selectivity == "free"] <- -2
  sel[selectivity == "logistic"] <- -1
  sel[selectivity == "dome"] <- 0

  if(any(is.na(sel))) {
    stop("Character entries for s_selectivity (for fleets) must be either: \"logistic\", \"dome\", or \"free\"", call. = FALSE)
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
  if(length(LWT$Index) != max(1, nsurvey)) stop("LWT$Index should be a vector of length ", data$nsurvey, ".")

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

  if(is.null(LWT$s_CAA)) {
    LWT$s_CAA <- rep(1, max(1, nsurvey))
  } else if(length(LWT$s_CAA) == 1 && nsurvey > 1) {
    LWT$s_CAA <- rep(LWT$s_CAA, nsurvey)
  }
  if(length(LWT$s_CAA) != max(1, nsurvey)) stop("LWT$s_CAA should be a vector of length ", nsurvey, ".")

  if(is.null(LWT$s_CAL)) {
    LWT$s_CAL <- rep(1, max(1, nsurvey))
  } else if(length(LWT$s_CAL) == 1 && nsurvey > 1) {
    LWT$s_CAL <- rep(LWT$s_CAL, nsurvey)
  }
  if(length(LWT$s_CAL) != max(1, nsurvey)) stop("LWT$s_CAL should be a vector of length ", nsurvey, ".")

  return(LWT)
}


update_SRA_data <- function(data, OM, condition, dots) {

  message("\nChecking OM and data...\n")

  # Preliminary OM check for basics
  if(length(OM@nyears) == 0) stop("OM@nyears is needed.", call. = FALSE)
  if(length(OM@maxage) == 0) stop("OM@maxage is needed.", call. = FALSE)

  assign_for_compatibility <- function(x) {
    if(is.null(data[[x]]) && !is.null(dots[[x]])) {
      data[[x]] <<- getElement(dots, x)
      dots[[x]] <<- NULL
    }
    invisible()
  }

  dat_names <- c("Chist", "Ehist", "Index", "I_sd", "I_type", "CAA", "CAL", "ML", "ML_sd", "C_eq", "E_eq", "s_CAA", "s_CAL", "length_bin")

  lapply(dat_names, assign_for_compatibility)

  if(is.null(data$Chist) && !is.null(data$Ehist)) condition <- "effort"

  if(condition == "catch" || condition == "catch2") {
    if(is.null(data$Chist)) {
      stop("Full time series of catch is needed.", call. = FALSE)
    } else {
      if(any(is.na(data$Chist))) {
        stop("One or more of the historical annual catch observations is missing. Suggestion: use linear interpolation to fill these data.", call. = FALSE)
      }
      if(any(data$Chist < 0)) stop("All catch values should be zero or greater.", call. = FALSE)

      # Convert single fleet inputs to multiple fleet, e.g. matrices to arrays
      if(!is.matrix(data$Chist)) data$Chist <- matrix(data$Chist, ncol = 1)

      data$nyears <- nrow(data$Chist)
      data$nfleet <- ncol(data$Chist)
      data$Ehist <- matrix(0, data$nyears, data$nfleet)

      if(is.null(data$Index) && is.null(data$CAA) && is.null(data$CAL) && is.null(data$ML) && is.null(data$Ehist)) {
        message("No data other than Chist is provided. Model will switch to conditioning on equilibrium effort.")
        data$condition <- "effort"
        data$Ehist <- matrix(1, data$nyears, data$nfleet)
        data$E_eq <- rep(1, data$nfleet)
      } else {
        data$condition <- condition
      }
    }
  }

  if(condition == "effort") {
    data$condition <- "effort"
    if(is.null(data$Ehist)) {
      stop("Full time series of effort is needed.")
    } else {
      if(any(is.na(data$Ehist))) stop("Effort time series is not complete (contains NA's")
      if(any(data$Ehist < 0)) stop("All effort values should be positive.")

      if(!is.matrix(data$Ehist)) data$Ehist <- matrix(data$Ehist, ncol = 1)

      data$nyears <- nrow(data$Ehist)
      data$nfleet <- ncol(data$Ehist)

      if(!is.null(data$Chist) && !is.matrix(data$Chist)) data$Chist <- matrix(data$Chist, ncol = 1)
      if(is.null(data$Chist)) data$Chist <- matrix(0, data$nyears, data$nfleet)
    }
  }

  message("SRA model is conditioned on: ", data$condition)
  message(data$nfleet, " fleet(s) detected.")
  message(data$nyears, " years of data detected.")

  # Match number of historical years of catch/effort to OM
  if(OM@nyears != data$nyears) {
    cpars_cond <- length(OM@cpars) > 0 && any(vapply(OM@cpars, function(x) class(x) == "matrix" || class(x) == "array", logical(1)))
    if(cpars_cond) {
      stmt <- paste0("OM@nyears != length(", ifelse(data$condition == "catch" || data$condition == "catch2", "Chist", "Ehist"), "). ",
                     "There will be indexing errors in your custom parameters (OM@cpars).")
      stop(stmt, call. = FALSE)
    } else {
      message("OM@nyears was updated to length(", ifelse(data$condition == "catch" || data$condition == "catch2", "Chist", "Ehist"), "): ", data$nyears)
      OM@nyears <- data$nyears
    }
  }
  if(length(OM@CurrentYr) == 0) OM@CurrentYr <- data$nyears

  # Indices
  if(!is.null(data$Index)) {
    if(is.vector(data$Index)) {
      if(length(data$Index) != data$nyears) stop("Length of Index vector does not equal nyears (", data$nyears, "). NAs are acceptable.", call. = FALSE)
      data$Index <- matrix(data$Index, ncol = 1)
    } else if(is.matrix(data$Index)) {
      if(nrow(data$Index) != data$nyears) stop("Number of rows of Index matrix does not equal nyears (", data$nyears, "). NAs are acceptable.", call. = FALSE)
    } else stop("Index is neither a vector nor a matrix.", call. = FALSE)

    data$nsurvey <- ncol(data$Index)
  } else {
    data$nsurvey <- 0
    data$Index <- matrix(NA, ncol = 1, nrow = data$nyears)
  }

  if(!is.null(data$I_sd)) {
    if(is.vector(data$I_sd)) {
      if(length(data$I_sd) != data$nyears) stop("Length of I_sd vector does not equal nyears (", data$nyears, ").", call. = FALSE)
      data$I_sd <- matrix(data$I_sd, ncol = 1)
    } else if(is.matrix(data$I_sd)) {
      if(nrow(data$I_sd) != data$nyears) stop("Number of rows of I_sd matrix does not equal nyears (", data$nyears, "). NAs are acceptable.", call. = FALSE)
      if(ncol(data$I_sd) != data$nsurvey) stop("Number of columns of I_sd matrix does not equal nsurvey (", data$nsurvey, ").", call. = FALSE)
    }

    SD_NA <- is.na(data$I_sd)
    if(sum(SD_NA)) {
      SD_out <- !is.na(data$Index[SD_NA])
      if(any(SD_out)) stop("There are NA's in data$I_sd for years associated with survey values in data$Index.", call. = FALSE)
    }
  }
  message(data$nsurvey, " survey(s) detected.")

  # Process age comps
  if(!is.null(data$CAA)) {

    if(is.matrix(data$CAA)) data$CAA <- array(data$CAA, c(dim(data$CAA), 1))

    if(dim(data$CAA)[1] != data$nyears) {
      stop("Number of CAA rows (", dim(data$CAA)[1], ") does not equal nyears (", data$nyears, "). NAs are acceptable.", call. = FALSE)
    }
    if(dim(data$CAA)[2] < OM@maxage) {
      message("Number of CAA columns (", dim(data$CAA)[2], ") does not equal OM@maxage (", OM@maxage, ").")
      message("Assuming no observations for ages greater than ", dim(data$CAA)[2], " and filling with zeros.")
      add_ages <- OM@maxage - dim(data$CAA)[2]
      CAA_new <- array(0, c(data$nyears, OM@maxage, data$nfleet))
      CAA_new[, 1:dim(data$CAA)[2], ] <- data$CAA
      data$CAA <- CAA_new
    }
    if(dim(data$CAA)[2] > OM@maxage) {
      OM@maxage <- dim(data$CAA)[2]
      message("Increasing OM@maxage to ", OM@maxage, ".")
    }
    if(dim(data$CAA)[3] != data$nfleet) {
      stop("Number of CAA slices (", dim(data$CAA)[3], ") does not equal nfleet (", data$nfleet, "). NAs are acceptable.", call. = FALSE)
    }

  } else {
    data$CAA <- array(0, c(data$nyears, OM@maxage, data$nfleet))
  }

  # Sample life history, selectivity, and obs parameters
  OM_samp <- check_OM_for_sampling(OM, data)

  old_warning <- options()$warn
  options(warn = -1)
  on.exit(options(warn = old_warning))

  set.seed(OM@seed)
  message("Getting biological parameters from OM...")
  suppressMessages({
    StockPars <- SampleStockPars(OM_samp, msg = FALSE)
    ObsPars <- SampleObsPars(OM_samp)
    FleetPars <- SampleFleetPars(OM_samp, msg = FALSE)
  })

  # Process length comps
  if(!is.null(data$CAL)) {
    if(is.matrix(data$CAL)) data$CAL <- array(data$CAL, c(dim(data$CAL), 1))

    if(is.null(data$length_bin)) {
      stop("You must specify length_bin, which is the mean length of each length bin (columns) of the CAL data.", call. = FALSE)
    }
    if(dim(data$CAL)[1] != data$nyears) {
      stop("Number of CAL rows (", dim(data$CAL)[1], ") does not equal nyears (", data$nyears, "). NAs are acceptable.", call. = FALSE)
    }
    if(dim(data$CAL)[2] != length(data$length_bin)) {
      stop("Number of CAL columns (", dim(data$CAL)[2], ") does not equal length(length_bin) (", length(data$length_bin), ").", call. = FALSE)
    }
    if(dim(data$CAL)[3] != data$nfleet) {
      stop("Number of CAL slices (", dim(data$CAA)[3], ") does not equal nfleet (", data$nfleet, "). NAs are acceptable.", call. = FALSE)
    }
  } else {
    data$CAL <- array(0, c(data$nyears, length(StockPars$CAL_binsmid), data$nfleet))
    data$length_bin <- StockPars$CAL_binsmid
  }

  # Process mean size (either lengths or weights)
  if(!is.null(data$ML) && is.null(data$MS)) { # Backwards compatibility
    message("\n\n ** data$ML is no longer used. Use data$MS (mean size) instead. ** \n\n")
    data$MS <- data$ML
  }
  if(!is.null(data$MS)) {
    if(is.null(data$MS_type)) {
      message("Mean size (data$MS) found, but not type (data$MS_type). Assuming it's mean length.")
      data$MS_type <- "length"
    } else {
      data$MS_type <- match.arg(data$MS_type, choices = c("length", "weight"))
      message("Mean ", data$MS_type, " data found.")
    }
    if(is.vector(data$MS)) {
      if(length(data$MS) != data$nyears) stop("Mean size vector (MS) must be of length ", data$nyears, ".", call. = FALSE)
      data$MS <- matrix(data$MS, ncol = 1)
    }
    if(nrow(data$MS) != data$nyears) stop("Number of MS rows (", nrow(data$ML), ") does not equal nyears (", data$nyears, "). NAs are acceptable.", call. = FALSE)
    if(ncol(data$MS) != data$nfleet) stop("Number of MS columns (", ncol(data$ML), ") does not equal nfleet (", data$nfleet, "). NAs are acceptable.", call. = FALSE)

    if(!is.null(data$ML_sd) && is.null(data$MS_cv)) { # Backwards compatibility
      message("\n\n ** data$ML_sd is no longer used. Use data$MS_cv instead. ** \n\n")
      data$MS_cv <- data$ML_sd/apply(data$MS, 2, mean, na.rm = TRUE)
      message("data$MS_cv = ", paste(data$MS_cv, collapse = ", "))
    }
    if(is.null(data$MS_cv)) {
      data$MS_cv <- rep(0.2, data$nfleet)
    } else if(length(data$MS_cv) == 1) data$MS_cv <- rep(data$MS_cv, data$nfleet)
    if(length(data$MS_cv) != data$nfleet) stop("Mean size CV vector (MS_cv) must be of length ", data$nfleet, ".", call. = FALSE)
  } else {
    data$MS <- matrix(NA, nrow = data$nyears, ncol = data$nfleet)
    data$MS_cv <- rep(0.2, data$nfleet)
    data$MS_type <- "length"
  }

  # Process equilibrium catch/effort - Ceq
  if(is.null(data$C_eq)) data$C_eq <- rep(0, data$nfleet)
  if(data$condition == "catch" || data$condition == "catch2") {
    if(length(data$C_eq) == 1) data$C_eq <- rep(data$C_eq, data$nfleet)
    if(length(data$C_eq) < data$nfleet) stop("C_eq needs to be of length nfleet (", data$nfleet, ").", call. = FALSE)
  }

  if(data$condition == "catch2" && any(data$C_eq > 0)) {
    message("Equilibrium catch was detected. The corresponding equilibrium F will be estimated.")
  }

  if(is.null(data$E_eq)) data$E_eq <- rep(0, data$nfleet)
  if(data$condition == "effort") {
    if(length(data$E_eq) == 1) data$E_eq <- rep(data$E_eq, data$nfleet)
    if(length(data$E_eq) < data$nfleet) stop("E_eq needs to be of length nfleet (", data$nfleet, ").", call. = FALSE)
  }

  # Process survey age comps
  if(!is.null(data$s_CAA)) {

    if(is.matrix(data$s_CAA)) data$s_CAA <- array(data$s_CAA, c(dim(data$s_CAA), 1))

    if(dim(data$s_CAA)[1] != data$nyears) {
      stop("Number of s_CAA rows (", dim(data$s_CAA)[1], ") does not equal nyears (", data$nyears, "). NAs are acceptable.", call. = FALSE)
    }
    if(dim(data$s_CAA)[2] < OM@maxage) {
      message("Number of s_CAA columns (", dim(data$s_CAA)[2], ") does not equal OM@maxage (", OM@maxage, ").")
      message("Assuming no observations for ages greater than ", dim(data$s_CAA)[2], " and filling with zeros.")
      add_ages <- OM@maxage - dim(data$s_CAA)[2]
      CAA_new <- array(0, c(data$nyears, OM@maxage, data$survey))
      CAA_new[, 1:dim(data$s_CAA)[2], ] <- data$s_CAA
      data$s_CAA <- CAA_new
    }
    if(dim(data$s_CAA)[2] > OM@maxage) {
      stop("Error in age dimension of s_CAA.", call. = FALSE)
    }
    if(dim(data$s_CAA)[3] != data$nsurvey) {
      stop("Number of CAA slices (", dim(data$s_CAA)[3], ") does not equal nsurvey (", data$nsurvey, "). NAs are acceptable.", call. = FALSE)
    }

  } else {
    data$s_CAA <- array(0, c(data$nyears, OM@maxage, ncol(data$Index)))
  }

  # Process survey length comps
  if(!is.null(data$s_CAL)) {
    if(is.matrix(data$s_CAL)) data$s_CAL <- array(data$s_CAL, c(dim(data$s_CAL), 1))

    if(is.null(data$length_bin)) {
      stop("You must specify length_bin, which is the mean length of each length bin (columns) of the CAL data.", call. = FALSE)
    }
    if(dim(data$s_CAL)[1] != data$nyears) {
      stop("Number of s_CAL rows (", dim(data$s_CAL)[1], ") does not equal nyears (", data$nyears, "). NAs are acceptable.", call. = FALSE)
    }
    if(dim(data$s_CAL)[2] != length(data$length_bin)) {
      stop("Number of s_CAL columns (", dim(data$s_CAL)[2], ") does not equal length(length_bin) (", length(data$length_bin), ").", call. = FALSE)
    }
    if(dim(data$s_CAL)[3] != data$nsurvey) {
      stop("Number of s_CAL slices (", dim(data$s_CAL)[3], ") does not equal nsurvey (", data$nsurvey, "). NAs are acceptable.", call. = FALSE)
    }
  } else {
    data$s_CAL <- array(0, c(data$nyears, length(data$length_bin), ncol(data$Index)))
  }

  # Absolute survey
  if(data$nsurvey > 0) {
    if(is.null(data$abs_I)) data$abs_I <- rep(0L, data$nsurvey)
    if(length(data$abs_I) < data$nsurvey) stop("abs_I should be of length", data$nsurvey, call. = FALSE)
  } else {
    data$abs_I <- 0L
  }

  # Index units - biomass/abundance
  if(data$nsurvey > 0) {
    if(!is.null(data$I_basis)) data$I_units <- data$I_basis # Backwards compatability with 1.4.4
    if(is.null(data$I_units)) data$I_units <- rep(1L, data$nsurvey)
    if(length(data$I_units) < data$nsurvey) stop("I_basis should be of length", data$nsurvey, call. = FALSE)
  } else {
    data$I_units <- 1L
  }

  # Ageing error
  if(is.null(data$age_error)) data$age_error <- diag(OM@maxage)
  if(any(dim(data$age_error) != OM@maxage)) stop("data$age_error should be a square matrix of OM@maxage rows and columns", call. = FALSE)

  # Sel_block dummy fleets
  if(is.null(data$sel_block)) {
    data$sel_block <- matrix(1:data$nfleet, nrow = data$nyears, ncol = data$nfleet, byrow = TRUE)
  } else {
    if(nrow(data$sel_block) != data$nyears) {
      stop(paste("data$sel_block should be a matrix of", data$nyears, "rows."), call. = FALSE)
    }
    if(ncol(data$sel_block) != data$nfleet) {
      stop(paste("data$sel_block should be a matrix of", data$nfleet, "columns."), call. = FALSE)
    }
  }
  data$nsel_block <- as.numeric(data$sel_block) %>% unique() %>% length()

  return(list(data = data, OM = OM, StockPars = StockPars, ObsPars = ObsPars, FleetPars = FleetPars))
}

check_OM_for_sampling <- function(OM, data) {
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
    any_CAL <- !is.null(data$CAL) && any(data$CAL > 0, na.rm = TRUE)
    any_ML <- !is.null(data$ML) && any(data$ML > 0, na.rm = TRUE)
    any_s_CAL <- !is.null(data$s_CAL) && any(data$s_CAL > 0, na.rm = TRUE)
    if(any_CAL || any_ML || any_s_CAL) {
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
    if(data$condition == "effort") {
      OM@R0 <- 1
    } else {
      OM@R0 <- 1e3
      message("OM@R0 is used as the starting value for R0 in SRA_scope, but was not found. By default, using 1000.")
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
  if(any(data$Index > 0, na.rm = TRUE)) {
    Isd_check <- !is.null(data$I_sd) && any(data$I_sd > 0, na.rm = TRUE)
    if(!Isd_check) {
      Isd_check2 <- length(OM@Iobs) == 2 || !is.null(cpars$Iobs)
      if(!Isd_check2) stop("OM@Iobs is needed.", call. = FALSE)
    }
  }
  Iobs <- OM@Iobs
  OM <- Replace(OM, DLMtool::Generic_Obs, silent = TRUE)
  OM@Iobs <- Iobs

  ###### Imp
  OM <- Replace(OM, DLMtool::Perfect_Imp, silent = TRUE)

  return(OM)
}
