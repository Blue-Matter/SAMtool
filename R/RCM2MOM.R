#' Convert RCM to a multi-fleet operating model (MOM)
#' 
#' The RCM (Rapid Conditioning Model) returns a single-fleet operating model, implying constant effort among fleets for projections. 
#' Here, we convert the single-fleet OM to a multi-fleet OM, preserving the multiple fleet structure used in the conditioning model
#' for projections. This allows for testing management procedures that explicitly specify fleet allocation in the management advice.
#' 
#' @param RCModel Output from [RCM], a class [RCModel-class] object.
#' @return A class [MSEtool::MOM-class] object.
#' @author Q. Huynh
#' @examples
#' \donttest{
#' data(pcod) 
#' mat_ogive <- pcod$OM@cpars$Mat_age[1, , 1]
#' OM <- MSEtool::SubCpars(pcod$OM, 1:3)
#' out <- RCM(OM = pcod$OM, data = pcod$data, 
#'            condition = "catch", mean_fit = TRUE,
#'            selectivity = "free", s_selectivity = rep("SSB", ncol(pcod$data@Index)),
#'            start = list(vul_par = matrix(mat_ogive, length(mat_ogive), 1)),
#'            map = list(vul_par = matrix(NA, length(mat_ogive), 1),
#'                       log_early_rec_dev = rep(1, pcod$OM@maxage)),
#'            prior = pcod$prior)
#' MOM <- RCM2MOM(out)
#' }
#' 
#' @importFrom abind abind
#' @export
RCM2MOM <- function(RCModel) {
  MOM <- suppressMessages(new("MOM"))
  
  nf <- ncol(RCModel@data@Chist)
  
  slot_intersect <- intersect(slotNames("MOM"), slotNames("OM"))
  for(i in slot_intersect) slot(MOM, i) <- slot(RCModel@OM, i)
  
  report <- RCModel@Misc
  if (length(report) == 1) report <- lapply(1:MOM@nsim, function(...) RCModel@Misc[[1]])
  cpars <- lapply(1:nf, function(f) {
    cp <- RCModel@OM@cpars
    
    F_age <- sapply(
      1:MOM@nsim, 
      function(x) report[[x]]$F[, f] * report[[x]]$vul[1:RCModel@OM@nyears, , f],
      simplify = "array"
    ) %>% 
      aperm(3:1)
    
    cp$Find <- apply(F_age, c(1, 3), max)
    Vhist <- apply(F_age, c(1, 3), function(x) x/max(x))
    Vpro <- array(Vhist[, , dim(Vhist)[3]], c(dim(Vhist)[1:2], RCModel@OM@proyears))
    cp$V <- abind::abind(Vhist, Vpro, along = 3) %>%
      aperm(c(2, 1, 3))
    
    if (!is.null(cp$SLarray)) {
      cp$SLarray <- lapply(report, make_SL, sel_block = RCModel@data@sel_block) %>% 
        lapply(expand_V_matrix, nyears = RCModel@OM@nyears, proyears = RCModel@OM@proyears) %>% 
        simplify2array() %>% aperm(c(3, 1, 2))
    }
    
    if (!is.null(cp$Data)) {
      if (sum(RCModel@data@Chist[, f] > 0, na.rm = TRUE)) {
        cp$Data@Cat <- matrix(RCModel@data@Chist[, f], 1, RCModel@OM@nyears)
        cp$Data@CV_Cat <- sqrt(exp(RCModel@data@C_sd[, f]^2 - 1)) %>% matrix(1, RCModel@OM@nyears)
      } else if (!all(is.na(cp$Data@Cat))) {
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
  
  if (nf > 1) {
    CatchFrac <- sapply(1:MOM@nsim, function(x) {
      cvec <- report[[x]]$Cpred
      cvec[nrow(cvec), ]/sum(cvec[nrow(cvec), ])
    }) %>% t()
  } else {
    CatchFrac <- matrix(1, MOM@nsim, 1)
  }
  
  MOM@CatchFrac <- list(CatchFrac)
  
  return(MOM)
}

