
#' The rapid conditioning model as an assessment function
#' 
#' In beta testing. A function that uses RCM as an assessment function for use in MPs. Many more function arguments will be added
#' to tinker with model settings and data inputs. Currently uses catch, CAA, and biomass index from Data object, as well as biological
#' parameters in \code{Data@@Misc$StockPars}.
#' 
#' @param x A position in the Data object (by default, equal to one for assessments).
#' @param Data An object of class Data
#' @param selectivity Whether to model "logistic" or "dome" selectivity for the fishery.
#' @param CAA_ESS The annual multinomial sample size for the fishery age composition.
#' @param ... Additional arguments (to be added).
#' 
#' @examples  
#' r <- RCM_assess(Data = SimulatedData)
#' myMP <- make_MP(RCM_assess, HCR_MSY)
#' myMP(x = 1, Data = SimulatedData)
#' @export
RCM_assess <- function(x = 1, Data, 
                       selectivity = c("logistic", "dome"), CAA_ESS = 50, ...) {
  RCMdata <- new("RCMdata")
  RCMdata@Misc$nyears <- nyears <- ncol(Data@Cat)
  RCMdata@Misc$nfleet <- nfleet <- 1
  RCMdata@Misc$condition <- "catch"
  
  n_age <- Data@MaxAge + 1
  RCMdata@Misc$lbinmid <- lbin <- Data@CAL_mids
  RCMdata@Misc$nlbin <- nlbin <- length(lbin)
  
  RCMdata@Chist <- matrix(Data@Cat[x, ], nyears, nfleet)
  RCMdata@C_sd <- array(0.1, dim(RCMdata@Chist))
  
  RCMdata@Ehist <- array(NA_real_, dim(RCMdata@Chist))
  
  RCMdata@CAA <- Data@CAA[x, , , drop = FALSE] %>% aperm(c(2, 3, 1))
  RCMdata@CAA_ESS <- array(CAA_ESS, dim(RCMdata@CAA)[c(1, 3)])
  
  selectivity <- match.arg(selectivity)
  sel <- int_sel(selectivity)
  
  # No length comps or mean size
  RCMdata@CAL <- array(0, c(nyears, nlbin, nfleet))
  RCMdata@CAL_ESS <- RCMdata@MS <- array(NA_real_, c(nyears, nfleet))
  RCMdata@MS_type <- "length"
  RCMdata@MS_cv <- 0.2
  
  # Only one survey for now
  AddInd <- "B"
  s_sel <- int_s_sel(AddInd)
  RCMdata@Misc$nsurvey <- nsurvey <- length(s_sel)
  Index <- Assess_I_hist(AddInd, Data, x, 1:nyears)
  RCMdata@Index <- matrix(Index$I_hist, nyears, nsurvey)
  RCMdata@I_sd <- matrix(Index$I_sd, nyears, nsurvey)
  RCMdata@I_units <- Index$I_units
  
  # No index age/length comps
  RCMdata@IAA <- array(NA_real_, c(nyears, n_age, nsurvey))
  RCMdata@IAL <- array(NA_real_, c(nyears, nlbin, nsurvey))
  RCMdata@IAA_ESS <- RCMdata@IAL_ESS <- matrix(0, nyears, nsurvey)
  
  # Misc arguments that need to be filled in
  RCMdata@C_eq <- RCMdata@E_eq <- rep(0, nfleet)
  RCMdata@C_eq_sd <- rep(0.01, nfleet)
  RCMdata@sel_block <- matrix(1, nyears, nfleet)
  RCMdata@Misc$nsel_block <- 1
  RCMdata@abs_I <- rep(0, nsurvey)
  RCMdata@age_error <- diag(n_age)
  
  LWT <- make_LWT(list(), nfleet, nsurvey)
  prior <- make_prior(list(), nsurvey, msg = FALSE)
  
  RCM_out <- RCM_est(x = x, RCMdata = RCMdata, selectivity = sel, s_selectivity = s_sel, LWT = LWT,
                     comp_like = "multinomial", prior = prior, StockPars = Data@Misc$StockPars, ObsPars = list(Isd = Data@Obs$Isd),
                     FleetPars = Data@Misc$FleetPars, mean_fit = FALSE)
  obj <- RCM_out$obj
  opt <- RCM_out$opt
  SD <- RCM_out$SD
  report <- RCM_out$report
  conv <- SD$pdHess
  
  Year <- Data@Year
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
  if(conv) {
    ref_pt <- lapply(1:nyears, function(y) {
      M <- obj$env$data$M_data[y, ]
      mat <- obj$env$data$mat[y, ]
      weight <- obj$env$data$wt[y, ]
      vul <- report$F_at_age[y, ]/max(report$F_at_age[y, ])
      SR <- obj$env$data$SR_type
      catch_eq <- "Baranov"
      tv_M <- "none"
      Arec <- report$Arec
      Brec <- report$Brec
      
      # Optimize for MSY
      opt2 <- optimize(yield_fn_SCA, interval = c(1e-4, 4), M = M, mat = mat, weight = weight, vul = vul, 
                       SR = SR, Arec = Arec, Brec = Brec, catch_eq = catch_eq, tv_M = tv_M)
      opt3 <- yield_fn_SCA(opt2$minimum, M = M, mat = mat, weight = weight, vul = vul, SR = SR, 
                           Arec = Arec, Brec = Brec, opt = FALSE, catch_eq = catch_eq, tv_M = tv_M)
      FMSY <- opt2$minimum
      MSY <- -1 * opt2$objective
      VBMSY <- opt3["VB"]
      RMSY <- opt3["R"]
      BMSY <- opt3["B"]
      EMSY <- opt3["E"]
      
      Fvec <- seq(0, 2.5 * FMSY, length.out = 100)
      
      # Yield curve
      yield <- lapply(Fvec, yield_fn_SCA, M = M, mat = mat, weight = weight, vul = vul, SR = SR, 
                      Arec = Arec, Brec = Brec, opt = FALSE, catch_eq = catch_eq, tv_M = tv_M)
      EPR <- vapply(yield, getElement, numeric(1), "EPR")
      YPR <- vapply(yield, getElement, numeric(1), "YPR")
      
      new_B0 <- yield[[1]]["B"] # New due to change in M
      new_E0 <- yield[[1]]["E"]
      new_VB0 <- yield[[1]]["VB"]
      new_R0 <- yield[[1]]["R"]
      if(SR == "BH") {
        new_h <- Arec * EPR[1]/ (4 + Arec * EPR[1])
      } else {
        new_h <- 0.2 * (Arec * EPR[1])^0.8
      }
      
      return(list(FMSY = FMSY, MSY = MSY, VBMSY = VBMSY, RMSY = RMSY, BMSY = BMSY, EMSY = EMSY,
                  per_recruit = data.frame(FM = Fvec, SPR = EPR/EPR[1], YPR = YPR), SR_par = SR,
                  new_B0 = new_B0, new_E0 = new_B0, new_VB0 = new_VB0, new_R0 = new_R0, new_h = new_h))
    })
    
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
