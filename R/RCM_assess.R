
#' The rapid conditioning model as an assessment function
#' 
#' In beta testing. A function that uses RCM as an assessment function for use in MPs. More function arguments will be added
#' to tinker with model settings and data inputs. 
#' 
#' @param x A position in the Data object (by default, equal to one for assessments).
#' @param Data An object of class Data
#' @param AddInd A vector of integers or character strings indicating the indices to be used in the model. Integers assign the index to
#' the corresponding index in Data@@AddInd, "B" (or 0) represents total biomass in Data@@Ind, "VB" represents vulnerable biomass in
#' Data@@VInd, and "SSB" represents spawning stock biomass in Data@@SpInd. Vulnerability to the survey is fixed in the model.
#' @param SR Stock-recruit function (either \code{"BH"} for Beverton-Holt or \code{"Ricker"}).
#' @param selectivity Whether to model "logistic" or "dome" selectivity for the fishery.
#' @param CAA_multiplier Numeric for data weighting of catch-at-age matrix. If greater than 1, then this is the maximum multinomial sample size
#' in any year. If less than one, then the multinomial sample size is this fraction of the sample size. 
#' @param prior A named list for the parameters of any priors to be added to the model. See documentation in \link{SCA}.
#' @param LWT A named list (Index, CAA, Catch) of likelihood weights for the data components. For the index, a vector of length survey. For
#' CAL and Catch, a single value.
#' @param StockPars Either a string ("Data" or "OM") to indicate whether to grab biological parameters from the Data object,
#' or operating model. Alternatively, a named list to provide custom parameters for the assessment.
#' @param ... Additional arguments (to be added).
#' 
#' @section Data:
#' Currently uses catch, CAA, and indices of abundance in the corresponding slots in the Data object.
#' 
#' @section StockPars:
#' Biological parameters can be used from the (1) Data object, (2) operating model, or (3) provided directly in the
#' \code{StockPars} argument.
#' 
#' Options 2 and 3 allow for time-varying growth, maturity, and natural mortality. Natural mortality can also be age-varying.
#' 
#' \code{StockPars} can be a named list of parameters used to provide inputs to the assessment model:
#' 
#' \itemize{
#' \item \code{Wt_age} - annual weight at age, array [sim, ages, year]
#' \item \code{Mat_age} - annual maturity at age, array [sim, ages, year]
#' \item \code{hs} - Stock-recruit steepness, vector of length [sim]
#' \item \code{M_ageArray} - annual natural mortality, array [sim, ages, year]  
#' }
#' 
#' @examples  
#' r <- RCM_assess(Data = SimulatedData)
#' myMP <- make_MP(RCM_assess, HCR_MSY)
#' myMP(x = 1, Data = SimulatedData)
#' @export
RCM_assess <- function(x = 1, Data, AddInd = "B", SR = c("BH", "Ricker"), 
                       selectivity = c("logistic", "dome"), CAA_multiplier = 50, 
                       prior = list(), LWT = list(), StockPars = "Data", ...) {
  dots <- list(...)
  SR <- match.arg(SR)
  selectivity <- match.arg(selectivity)
  
  if (any(names(dots) == "yind")) {
    yind <- eval(dots$yind)
  } else {
    yind <- which(!is.na(Data@Cat[x, ]))[1]
    yind <- yind:length(Data@Cat[x, ])
  }
  
  RCMdata <- new("RCMdata")
  RCMdata@Misc$nyears <- nyears <- length(yind)
  RCMdata@Misc$nfleet <- nfleet <- 1
  RCMdata@Misc$condition <- "catch"
  
  n_age <- Data@MaxAge + 1
  
  # Age based model only
  RCMdata@Misc$lbinmid <- lbin <- c(5, 10)
  RCMdata@Misc$nlbin <- nlbin <- length(lbin)
  
  RCMdata@Chist <- matrix(Data@Cat[x, yind], nyears, nfleet)
  RCMdata@C_sd <- sdconv(1, Data@CV_Cat[x, yind]) %>% matrix(nyears, nfleet)
  
  RCMdata@Ehist <- array(NA_real_, dim(RCMdata@Chist))
  
  RCMdata@CAA <- Data@CAA[x, yind, , drop = FALSE] %>% aperm(c(2, 3, 1))
  RCMdata@CAA_ESS <- local({
    N <- apply(RCMdata@CAA, c(1, 3), sum, na.rm = TRUE)
    if (CAA_multiplier > 1) {
      pmin(N, CAA_multiplier)
    } else {
      N * CAA_multiplier
    }
  })
  
  sel <- int_sel(selectivity, silent = TRUE)
  
  # No length comps or mean size
  RCMdata@CAL <- array(NA_real_, c(nyears, max(c(nlbin, 2)), nfleet))
  RCMdata@CAL_ESS <- array(0, c(nyears, nfleet))
  RCMdata@MS <- array(NA_real_, c(nyears, nfleet))
  RCMdata@MS_type <- "length"
  RCMdata@MS_cv <- 0.2
  
  # Indices
  RCMdata@Misc$nsurvey <- nsurvey <- length(AddInd)
  Index <- lapply(AddInd, Assess_I_hist, Data, x, yind)
  RCMdata@Index <- vapply(Index, getElement, numeric(nyears), "I_hist")
  RCMdata@I_sd <- vapply(Index, getElement, numeric(nyears), "I_sd")
  RCMdata@I_units <- vapply(Index, getElement, numeric(1), "I_units")
  
  s_sel <- local({
    vout <- ifelse(AddInd == "VB", 1, AddInd) # AddInd - fix selectivity to values in Data@AddIndV
    int_s_sel(vout, nfleet, silent = TRUE)
  })
  
  # No index age/length comps
  RCMdata@IAA <- array(NA_real_, c(nyears, n_age, nsurvey))
  RCMdata@IAL <- array(NA_real_, c(nyears, max(c(nlbin,2)), nsurvey))
  RCMdata@IAA_ESS <- RCMdata@IAL_ESS <- matrix(0, nyears, nsurvey)
  
  # Misc arguments that need to be filled in - no equilibrium catch right now
  RCMdata@C_eq <- RCMdata@E_eq <- rep(0, nfleet)
  RCMdata@C_eq_sd <- rep(0.01, nfleet)
  
  RCMdata@sel_block <- matrix(1, nyears, nfleet)
  RCMdata@Misc$nsel_block <- 1
  
  RCMdata@abs_I <- rep(0, nsurvey)
  RCMdata@age_error <- diag(n_age)
  
  LWT <- make_LWT(list(), nfleet, nsurvey)
  prior <- make_prior(prior, nsurvey, msg = FALSE)
  
  Ages <- 0:Data@MaxAge
  n_age <- length(Ages)
  nsim <- 1
  
  # Create StockPars_out from Data, OM, or custom list
  if (is.character(StockPars)) {
    if (StockPars == "Data") {
      StockPars_out <- RCM_assess_StockPars(x, Data, list(), n_age, nyears, nsim, SR)
    } else if (StockPars == "OM") {
      StockPars_out <- RCM_assess_StockPars(x, Data, Data@Misc$StockPars, n_age, nyears, nsim, SR)
    } else {
      stop("StockPars should be either \"Data\", \"OM\", or a list of parameters")
    }
  } else {
    StockPars_out <- RCM_assess_StockPars(x, Data, StockPars, n_age, nyears, nsim, SR)
  }
  
  # Create FleetPars from Data
  FleetPars <- list()
  #FleetPars$LFS_y <- array(Data@LFS[x], dim=c(1, nyears+1))
  #FleetPars$L5_y <- array(Data@LFC[x], dim=c(1, nyears+1))
  FleetPars$LFS_y <- apply(RCMdata@CAA, 2, sum, na.rm = TRUE) %>% which.max() %>% pmax(1) %>% pmin(0.75 * n_age) %>%
    array(c(1, nyears + 1))
  FleetPars$L5_y <- FleetPars$LFS_y - 1
  FleetPars$Vmaxlen_y <- ifelse(is.null(Data@Vmaxlen[x]), 0.5, Data@Vmaxlen[x]) %>% array(c(1, nyears+1))
  
  RCM_out <- RCM_est(x = 1, RCMdata = RCMdata, selectivity = sel, s_selectivity = s_sel, LWT = LWT,
                     comp_like = "multinomial", prior = prior, StockPars = StockPars_out,
                     FleetPars = FleetPars, mean_fit = FALSE)
  obj <- RCM_out$obj
  opt <- RCM_out$opt
  SD <- RCM_out$SD
  report <- RCM_out$report
  conv <- SD$pdHess
  
  Year <- Data@Year[yind]
  Yearplusone <- YearR <- c(Year, max(Year) + 1)
  Assessment <- new("Assessment", Model = "RCM", Name = Data@Name, conv = conv,
                    h = report$h, 
                    FMort = structure(apply(report$F_at_age, 1, max), names = Year), # Annual apical F
                    B = structure(report$B, names = Yearplusone),
                    SSB = structure(report$E, names = Yearplusone),
                    R = structure(report$R, names = YearR),
                    N = structure(rowSums(report$N), names = Yearplusone),
                    N_at_age = report$N,
                    Selectivity = report$F_at_age/apply(report$F_at_age, 1, max), # Aggregate selectivity
                    Dev = structure(report$log_rec_dev, names = Year), Dev_type = "log-Recruitment deviations",
                    NLL = ifelse(is.character(opt), NA_real_, opt$objective),
                    obj = obj, opt = opt, SD = SD, TMB_report = report)
  
  # Calculate annual reference points
  if (conv) {
    ref_pt <- RCM_assess_ref(obj, report, yref = 1:nyears)
    
    report$FMSY <- sapply(ref_pt, getElement, "FMSY")
    tv_ref_pt <- length(unique(report$FMSY)) > 1
    
    report$MSY <- sapply(ref_pt, getElement, "MSY")
    report$VBMSY <- sapply(ref_pt, getElement, "VBMSY")
    report$RMSY <- sapply(ref_pt, getElement, "RMSY")
    report$BMSY <- sapply(ref_pt, getElement, "BMSY")
    report$EMSY <- sapply(ref_pt, getElement, "EMSY")
    
    refyear <- nyears # Year of reference points for reporting 
    
    report$new_B0 <- Assessment@B0 <- ref_pt[[refyear]]$new_B0
    report$new_E0 <- Assessment@SSB0 <- ref_pt[[refyear]]$new_E0
    report$new_VB0 <- Assessment@VB0 <- ref_pt[[refyear]]$new_VB0
    report$new_R0 <- Assessment@R0 <- ref_pt[[refyear]]$new_R0
    report$new_h <- Assessment@h <- ref_pt[[refyear]]$new_h
    
    Assessment@B_B0 <- Assessment@B/Assessment@B0
    Assessment@SSB_SSB0 <- Assessment@SSB/Assessment@SSB0
    Assessment@VB_VB0 <- Assessment@VB/Assessment@VB0
    
    Assessment@FMSY <- report$FMSY[refyear]
    Assessment@F_FMSY <- structure(Assessment@FMort/Assessment@FMSY, names = Year)
    
    Assessment@MSY <- report$MSY[refyear]
    Assessment@BMSY <- report$BMSY[refyear]
    Assessment@SSBMSY <- report$EMSY[refyear]
    Assessment@VBMSY <- report$VBMSY[refyear]
    Assessment@B_BMSY <- structure(Assessment@B/Assessment@BMSY, names = Yearplusone)
    Assessment@SSB_SSBMSY <- structure(Assessment@SSB/Assessment@SSBMSY, names = Yearplusone)
    #Assessment@VB_VBMSY <- structure(Assessment@VB/Assessment@VBMSY, names = Yearplusone)
    Assessment@TMB_report <- report
    
    catch_eq_fn <- function(Ftarget) {
      catch_equation(method = "Baranov", Ftarget = Ftarget, 
                     M = obj$env$data$M_data[nyears, ],
                     wt = obj$env$data$wt[nyears, ],
                     N = report$N[nyears + 1, ],
                     sel = report$F_at_age[nyears, ]/max(report$F_at_age[nyears, ]))
    }
    Assessment@forecast <- list(per_recruit = ref_pt[[refyear]][["per_recruit"]], catch_eq = catch_eq_fn)
  }
  return(Assessment)
}
class(RCM_assess) <- "Assess"


RCM_assess_StockPars <- function(x, Data, StockPars = list(), n_age, nyears, nsim, SR) {
  out <- list()
  Ages <- 1:n_age - 1
  
  # Age-only model for now
  if (is.null(StockPars$Len_age)) {
    Len_age <- Data@vbLinf[x]*(1-exp(-Data@vbK[x]*(Ages-Data@vbt0[x])))
    out$Len_age <- local({
      Len_age <- Data@vbLinf[x]*(1-exp(-Data@vbK[x]*(Ages-Data@vbt0[x])))
      array(Len_age, dim = c(1, n_age, nyears+1))
    })
  } else {
    out$Len_age <- StockPars$Len_age[x, 1:n_age, 1:(nyears+1), drop = FALSE]
  }
  
  if (is.null(StockPars$Linf)) {
    out$Linf <- Data@vbLinf[x]
  } else {
    out$Linf <- StockPars$Linf[x]
  }
  
  if (is.null(StockPars$LatASD)) {
    Len_age <- Data@vbLinf[x]*(1-exp(-Data@vbK[x]*(Ages-Data@vbt0[x])))
    out$LatASD <- out$Len_age * Data@LenCV[x]
  } else {
    out$LatASD <- StockPars$LatASD[x, 1:n_age, 1:(nyears+1), drop = FALSE]
  }
  
  if (is.null(StockPars$Wt_age)) {
    out$Wt_age <- local({
      Len_age <- Data@vbLinf[x]*(1-exp(-Data@vbK[x]*(Ages-Data@vbt0[x])))
      array(Data@wla[x] * Len_age ^ Data@wlb[x], dim = c(1, n_age, nyears+1))
    })
  } else {
    out$Wt_age <- StockPars$Wt_age[x, 1:n_age, 1:(nyears+1), drop = FALSE]
  }
  
  if (is.null(StockPars$ageM)) {
    out$ageM <- min(0.5 * Data@MaxAge, iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x], Data@L50[x])) %>% matrix(1, 1)
  } else {
    out$ageM <- StockPars$ageM[x]
  }
  if (!is.matrix(out$ageM)) out$ageM <- matrix(out$ageM, 1, 1)
  
  if (is.null(StockPars$Mat_age)) {
    out$Mat_age <- local({
      A50 <- min(0.5 * Data@MaxAge, iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x], Data@L50[x]))
      A95 <- max(A50+0.5, iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x], Data@L95[x]))
      m <- c(0, 1/(1 + exp(-log(19) * (c(1:Data@MaxAge) - A50)/(A95 - A50)))) # Age-0 is immature
      array(m, dim = c(1, n_age, nyears+1))
    })
  } else {
    out$Mat_age <- StockPars$Mat_age[x, 1:n_age, 1:(nyears+1), drop = FALSE]
  }
  
  if (is.null(StockPars$SRrel)) {
    out$SRrel <- ifelse(SR == "BH", 1, 2)
  } else {
    out$SRrel <- StockPars$SRrel[x]
  }
  
  if (is.null(StockPars$procsd)) {
    out$procsd <- Data@sigmaR[x]
  } else {
    out$procsd <- StockPars$procsd
  }
  
  if (is.null(StockPars$R0)) {
    out$R0 <- 2 * mean(Data@Cat[x, ])
  } else {
    out$R0 <- StockPars$R0[x]
  }
  
  if (is.null(StockPars$hs)) {
    out$hs <- Data@steep[x]
  } else {
    out$hs <- StockPars$hs[x]
  }
  
  if (is.null(StockPars$M_ageArray)) {
    out$M_ageArray <- array(Data@Mort[x], dim = c(1, n_age, nyears+1))
  } else {
    out$M_ageArray <- StockPars$M_ageArray[x, 1:n_age, 1:(nyears+1), drop = FALSE]
  }
  
  # Fec_age
  
  check <- vapply(out, function(y) any(is.na(y)), logical(1))
  if (any(check)) stop("Input parameters not found for RCM_assess: ", paste(names(check)[check], collapse = ", "))
  return(out)
}


RCM_assess_ref <- function(obj, yref = 1:obj$env$data$n_y) {
  
  ref_pt <- lapply(yref, function(y) {
    M <- obj$env$data$M_data[y, ]
    mat <- obj$env$data$mat[y, ]
    weight <- obj$env$data$wt[y, ]
    fec <- obj$env$data$fec[y, ]
    spawn_time_frac <- obj$env$data$spawn_time_frac
    
    vul <- report$F_at_age[y, ]/max(report$F_at_age[y, ])
    SR <- obj$env$data$SR_type
    catch_eq <- "Baranov"
    tv_M <- "none"
    Arec <- report$Arec
    Brec <- report$Brec
    
    # Optimize for MSY
    opt2 <- optimize(yield_fn_SCA, interval = c(1e-4, 4), M = M, mat = mat, weight = weight, fec = fec, 
                     vul = vul, SR = SR, Arec = Arec, Brec = Brec, catch_eq = catch_eq, tv_M = tv_M, 
                     spawn_time_frac = spawn_time_frac)
    opt3 <- yield_fn_SCA(opt2$minimum, M = M, mat = mat, weight = weight, fec = fec, vul = vul, SR = SR, 
                         Arec = Arec, Brec = Brec, opt = FALSE, catch_eq = catch_eq, tv_M = tv_M,
                         spawn_time_frac = spawn_time_frac)
    
    FMSY <- opt2$minimum
    MSY <- -1 * opt2$objective
    VBMSY <- opt3["VB"]
    RMSY <- opt3["R"]
    BMSY <- opt3["B"]
    EMSY <- opt3["E"]
    
    Fvec <- seq(0, 2.5 * FMSY, length.out = 100)
    
    # Yield curve
    yield <- lapply(Fvec, yield_fn_SCA, M = M, mat = mat, weight = weight, fec = fec, vul = vul, SR = SR, 
                    Arec = Arec, Brec = Brec, opt = FALSE, catch_eq = catch_eq, tv_M = tv_M, 
                    spawn_time_frac = spawn_time_frac)
    EPR <- vapply(yield, getElement, numeric(1), "EPR")
    YPR <- vapply(yield, getElement, numeric(1), "YPR")
    
    new_B0 <- yield[[1]]["B"] # New due to change in M
    new_E0 <- yield[[1]]["E"]
    new_VB0 <- yield[[1]]["VB"]
    new_R0 <- yield[[1]]["R"]
    if (SR == "BH") {
      new_h <- Arec * EPR[1]/ (4 + Arec * EPR[1])
    } else {
      new_h <- 0.2 * (Arec * EPR[1])^0.8
    }
    
    return(list(FMSY = FMSY, MSY = MSY, VBMSY = VBMSY, RMSY = RMSY, BMSY = BMSY, EMSY = EMSY,
                per_recruit = data.frame(FM = Fvec, SPR = EPR/EPR[1], YPR = YPR), SR_par = SR,
                new_B0 = new_B0, new_E0 = new_B0, new_VB0 = new_VB0, new_R0 = new_R0, new_h = new_h))
  })
  
  return(ref_pt)
}
  
