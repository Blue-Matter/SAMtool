#' Convert RCM to a multi-fleet operating model (MOM)
#' 
#' The RCM (Rapid Conditioning Model) returns a single-fleet operating model, implying constant effort among fleets for projections. 
#' Here, we convert the single-fleet OM to a multi-fleet OM, preserving the multiple fleet structure used in the conditioning model
#' for projections. This allows for testing management procedures that explicitly specify fleet allocation in the management advice.
#' 
#' @param RCModel Output from \link{RCM}, a class \linkS4class{RCModel} object.
#' @return A class \linkS4class{MOM} object.
#' @author Q. Huynh
#' @export
RCM2MOM <- function(RCModel) {
  if(!requireNamespace("abind", quietly = TRUE)) stop("Install the abind package to use this function.")
  
  MOM <- suppressMessages(new("MOM"))
  
  nf <- ncol(RCModel@data@Chist)
  
  slot_intersect <- intersect(slotNames("MOM"), slotNames("OM"))
  for(i in slot_intersect) slot(MOM, i) <- slot(RCModel@OM, i)
  
  cpars <- lapply(1:nf, function(f) {
    cp <- RCModel@OM@cpars
    cp$Find <- sapply(1:MOM@nsim, function(x) RCModel@Misc[[x]]$F[, f]) %>% t()
    Vhist <- sapply(1:MOM@nsim, function(x) RCModel@Misc[[x]]$vul[1:RCModel@OM@nyears, , f], simplify = "array") %>%
      aperm(3:1) 
    Vpro <- array(Vhist[, , dim(Vhist)[3]], c(dim(Vhist)[1:2], RCModel@OM@proyears))
    cp$V <- abind::abind(Vhist, Vpro, along = 3)
    
    if(!is.null(cp$Data)) {
      if (sum(RCModel@data@Chist[, f] > 0, na.rm = TRUE)) {
        cp$Data@Cat <- matrix(RCModel@data@Chist[, f], 1, RCModel@OM@nyears)
        cp$Data@CV_Cat <- sqrt(exp(RCModel@data@C_sd[, f]^2 - 1)) %>% matrix(1, RCModel@OM@nyears)
      } else if(!all(is.na(cp$Data@Cat))) {
        cp$Data@Cat <- new("Data")@Cat
        cp$Data@CV_Cat <- new("Data")@CV_Cat
      }
      if (nf > 1 && sum(RCModel@data@CAA[, , f] > 0, na.rm = TRUE)) {
        cp$Data@CAA <- aperm(RCModel@data@CAA[, , f, drop = FALSE], c(3, 1, 2))
      }
      if (nf > 1 && sum(RCModel@data@CAL[, , f] > 0, na.rm = TRUE)) {
        cp$Data@CAL <- aperm(RCModel@data@CAL[, , f, drop = FALSE], c(3, 1, 2))
        cp$Data@CAL_mids <- RCModel@data@Misc$lbinmid
        cp$Data@CAL_bins <- RCModel@data@Misc$lbin
      }
    }
    return(cp)
  })
  MOM@cpars <- list(cpars)
    
  MOM@Stocks <- list(SubOM(RCModel@OM, "Stock"))
  Fleets <- lapply(1:nf, function(f) {
    ff <- SubOM(RCModel@OM, "Fleet")
    ff@Name <- paste(ff@Name, "Fleet", f)
    ff
  })
  MOM@Fleets <- list(Fleets)
  
  MOM@Obs <- list(lapply(1:nf, function(f) SubOM(RCModel@OM, "Obs")))
  MOM@Imps <- list(lapply(1:nf, function(f) SubOM(RCModel@OM, "Imp")))
  
  CatchFrac <- sapply(1:MOM@nsim, function(x) {
    cvec <- RCModel@Misc[[x]]$Cpred
    cvec[nrow(cvec), ]/sum(cvec[nrow(cvec), ])
  }) %>% t()
  
  MOM@CatchFrac <- list(CatchFrac)
  
  return(MOM)
}

