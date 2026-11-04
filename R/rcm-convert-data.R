
# Convert a data object to RCMdata, keeping each fleet separate

# data arrays are Year x Fleet x Class; RCMdata uses Year x Class x Fleet
.CompData2Array <- function(compdata) {
  if (is.null(compdata) || !length(compdata@Value)) return(NULL)
  aperm(compdata@Value, c(1, 3, 2))
}

.SDFromCV <- function(CV, dims) {
  if (is.null(CV) || !length(CV)) return(array(NA_real_, dims))
  array(MSEtool::sdconv(1, CV), dims)
}

.data2RCMdata <- function(data, condition = "catch") {
  if (!inherits(data, "data")) cli::cli_abort("{.arg data} must be a {.cls data} object.")

  RCMdata <- methods::new("RCMdata")
  Land <- data@Landings
  Disc <- data@Discards

  if (any(condition == "effort") && length(data@Effort@Value)) {
    RCMdata@Ehist <- as.matrix(data@Effort@Value)
  }

  if (length(Land@Value) || length(Disc@Value)) {
    Chist <- if (length(Land@Value)) as.matrix(Land@Value) else NULL
    if (length(Disc@Value)) Chist <- if (is.null(Chist)) as.matrix(Disc@Value) else Chist + as.matrix(Disc@Value)
    RCMdata@Chist <- Chist
    if (length(Land@CV)) RCMdata@C_sd <- .SDFromCV(Land@CV, dim(Chist))
  }

  idx_objs <- Filter(function(x) length(x@Value) > 0, list(CPUE = data@CPUE, Survey = data@Survey))
  no_cv <- names(idx_objs)[vapply(idx_objs, function(x) !length(x@CV) || all(is.na(x@CV)), logical(1))]
  if (length(no_cv)) cli::cli_abort("Index CVs are required in {.code data@{no_cv}@CV}.")
  if (length(idx_objs)) {
    RCMdata@Index <- do.call(cbind, lapply(idx_objs, function(x) as.matrix(x@Value)))
    RCMdata@I_sd <- do.call(cbind, lapply(idx_objs, function(x) .SDFromCV(x@CV, dim(as.matrix(x@Value)))))
    RCMdata@I_units <- unlist(lapply(idx_objs, function(x) {
      units <- x@Units
      if (!length(units)) units <- rep("Biomass", ncol(as.matrix(x@Value)))
      ifelse(units == "Number", 0L, 1L)
    }))
  }

  CAA <- .CompData2Array(data@LandingsAtAge)
  if (!is.null(CAA)) {
    RCMdata@CAA <- CAA
    RCMdata@CAA_ESS <- apply(CAA, c(1, 3), sum, na.rm = TRUE)
  }

  CAL <- .CompData2Array(data@LandingsAtSize)
  if (!is.null(CAL)) {
    RCMdata@CAL <- CAL
    RCMdata@CAL_ESS <- apply(CAL, c(1, 3), sum, na.rm = TRUE)
    cls <- data@LandingsAtSize@Classes
    if (is.list(cls)) cls <- cls[[1]]
    RCMdata@length_bin <- cls
  }

  RCMdata
}

# Default names for unnamed catch, effort and index columns
.FillDataNames <- function(data, fleet_names) {
  for (sl in c("Landings", "Discards", "Effort")) {
    obj <- methods::slot(data, sl)
    if (length(obj@Value) && !length(obj@Name)) {
      obj@Name <- fleet_names[seq_len(ncol(as.matrix(obj@Value)))]
      methods::slot(data, sl) <- obj
    }
  }
  for (sl in c("CPUE", "Survey")) {
    obj <- methods::slot(data, sl)
    if (length(obj@Value) && !length(obj@Name)) {
      obj@Name <- paste0(sl, seq_len(ncol(as.matrix(obj@Value))))
      methods::slot(data, sl) <- obj
    }
  }
  data
}
