
#' @rdname make_MP
#' @aliases make_interim_MP
#' @param AddInd A vector of integers or character strings indicating the indices to be used in the assessment model. 
#' Integers assign the index to the corresponding index in Data@@AddInd, "B" (or 0) represents total biomass in Data@@Ind, 
#' "VB" represents vulnerable biomass in Data@@VInd, and "SSB" represents spawning stock biomass in Data@@SpInd. For the interim
#' procedure, the function will use the first index in \code{AddInd}.
#' @param assessment_interval The time interval for when the assesment model is applied (number of years). In all other years, the 
#' interim procedure is applied.
#' @param type How the index is used to calculate the TAC in the interim procedure. See details.
#' @param type_par A control parameter for the interim procedure. See details.
#' 
#' @details 
#' \code{make_interim_MP} creates an MP that runs the interim procedure (updating the TAC according to index observations in between periodic 
#' assessment intervals. \strong{Always ensure to set:} \code{OM@@interval <- 1}. The assessment frequency is specified in argument 
#' \code{assessment_interval}.
#' 
#' In the year when the assessment is applied, the TAC is set by fitting the model and then running the harvest control rule. Between assessments,
#' the TAC is updated as
#' \deqn{
#' TAC_{y+1} = Cref (I_y + b \times s)/(Iref + b \times s)
#' }
#' where \code{Cref} is the TAC calculated from the most recent assessment, \code{Iref} is the value of the index when \code{Cref} was calculated
#' (see Equations 6 and 7 of Huynh et al. 2020). The value of \code{I_y} depends on \code{type}, with \code{b} and \code{s} equal zero unless
#' \code{type = "buffer"}:
#' 
#' \itemize{
#' \item \code{"buffer"} - \code{I_y} is the most recent index with \code{b} is specifed by \code{type_par} (default = 1), and \code{s} is 
#' the standard deviation of index residuals from the most recent assessment. 
#' \item \code{"mean"} - \code{I_y} is the mean value of the index over the most recent \code{type_par} years (default = 3).
#' \item \code{"loess"} - \code{I_y} is the most recent index predicted by a \link[stats]{loess} smoother applied over the entire time series of the index.
#' Use \code{type_par} to adjust the \code{span} parameter (default = 0.75).
#' \item \code{"none"} - \code{I_y} is the most recent index. Index values are not adjusted in the interim procedure.
#' }
#' 
#' @references 
#' Huynh et al. 2020. The interim management procedure approach for assessed stocks: Responsive management advice and lower assessment
#' frequency. Fish Fish. 21:663–679. \doi{10.1111/faf.12453}
#' @examples 
#' \donttest{
#' # Interim MPs
#' MP_buffer_5 <- make_interim_MP(assessment_interval = 5)
#' MP_buffer_10 <- make_interim_MP(assessment_interval = 10)
#' OM <- MSEtool::testOM
#' OM@@interval <- 1
#' 
#' MSE <- MSEtool::runMSE(OM, MPs = c("MP_buffer_5", "MP_buffer_10")) 
#' }
#' @export
make_interim_MP <- function(.Assess = "SCA", .HCR = "HCR_MSY", AddInd = "VB", assessment_interval = 5, type = c("buffer", "mean", "loess", "none"),
                            type_par = NULL, diagnostic = c("min", "full", "none"), ...) {
  type <- match.arg(type)
  diagnostic <- match.arg(diagnostic)
  
  if(is.character(.Assess)) {
    Assess_char <- .Assess
    .Assess <- as.symbol(.Assess)
  } else {
    .Assess <- substitute(.Assess)
    Assess_char <- as.character(.Assess)
  }
  if(is.character(.HCR)) {
    .HCR <- as.symbol(.HCR)
  } else {
    .HCR <- substitute(.HCR)
  }
  if(!inherits(eval(.Assess), "Assess")) {
    stop(paste(.Assess, "does not belong to class 'Assess'. Use: avail('Assess') to find eligible objects."))
  }
  if(!inherits(eval(.HCR), "HCR")) {
    stop(paste(.HCR, "does not belong to class 'HCR.' Use: avail('HCR') to find eligible objects."))
  }
  
  dots <- list(...)
  dots$AddInd <- AddInd
  dots_in_Assess <- dots[match(names(formals(eval(.Assess))), names(dots), nomatch = 0)]
  dots_in_HCR <- dots[match(names(formals(eval(.HCR))), names(dots), nomatch = 0)]
  
  Assess_call <- as.call(c(.Assess, x = quote(x), Data = quote(Data), dots_in_Assess))
  HCR_call <- as.call(c(.HCR, Assessment = quote(do_Assessment), reps = quote(reps), dots_in_HCR))
  
  MP_body <- bquote({
    dependencies <- .(get_dependencies(Assess_char, dots))
    
    ny <- length(Data@Year)
    Current_Yr <- Data@Year[ny]
    first_assess_yr <- Current_Yr == Data@LHYear
    run_assessment <- first_assess_yr || Current_Yr == Data@Misc[[x]]$interim$next_assess_yr
    run_interim <- !run_assessment
    
    if(run_assessment) {
      do_Assessment <- .(Assess_call)
      
      if(diagnostic != "none") {
        Assess_diag_output <- Assess_diagnostic(x, Data, do_Assessment, include_assessment = .(diagnostic == "full"))
      }
      
      if(do_Assessment@conv) { # Assessment converged. Run the HCR and report parameters for future interim procedure
        Rec <- .(HCR_call)
        if(diagnostic != "none") Rec@Misc$diagnostic <- Assess_diag_output$diagnostic
        if(!is.null(do_Assessment@info$Misc)) Rec@Misc <- c(Rec@Misc, do_Assessment@info$Misc)
        
        Rec@Misc$interim <- list(Cref = Rec@TAC, # A vector of nsim reps
                                 Iref = do_Assessment@Index[ny, 1],
                                 next_assess_yr = Current_Yr + assessment_interval)
        if(type == "buffer") {
          Rec@Misc$interim$s <- log(do_Assessment@Obs_Index[, 1]/do_Assessment@Index[, 1]) %>% sd(na.rm = TRUE)
        } else {
          Rec@Misc$interim$s <- 0
        }
        return(Rec)
        
      } else { # Assessment did not converge. Try the assessment again next year
        
        next_assess_yr <- Current_Yr + 1
        
        if(first_assess_yr || is.null(Data@Misc[[x]]$interim)) { # No assessment has been run in the past, can't do interim procedure, return TAC = NA
          
          Rec <- new("Rec")
          Rec@TAC <- rep(NA_real_, reps)
          
          if(diagnostic != "none") Rec@Misc$diagnostic <- Assess_diag_output$diagnostic
          if(!is.null(do_Assessment@info$Misc)) Rec@Misc <- c(Rec@Misc, do_Assessment@info$Misc)
          
          Rec@Misc$interim <- list(Cref = NA_real_, Iref = NA_real_, next_assess_yr = next_assess_yr)
          return(Rec)
          
        } else { # There should be a previous assessment to continue the interim procedure for the current year
          run_interim <- TRUE
        }
      }
    }
    
    if(!run_assessment || run_interim) {
      Cref <- Data@Misc[[x]]$interim$Cref
      Iref <- Data@Misc[[x]]$interim$Iref
      I_y <- switch(type,
                    "buffer" = interim_get_index(AddInd[1], x, Data, ny),
                    "mean" =  local({
                      if(is.null(type_par)) type_par <- 3
                      interim_get_index(AddInd[1], x, Data, seq(ny - type_par + 1, ny)) %>% mean(na.rm = TRUE)
                    }),
                    "loess" = local({
                      I_df <- data.frame(Year = Data@Year, Ind = interim_get_index(AddInd[1], x, Data, 1:ny))
                      if(is.null(type_par)) type_par <- formals(loess)$span
                      fit <- loess(Ind ~ Year, data = I_df, span = type_par)
                      fit$fitted[length(fit$fitted)]
                    }),
                    "none" = interim_get_index(AddInd[1], x, Data, ny)
      )
      
      if(type == "buffer") {
        if(is.null(type_par)) {
          b <- 1
        } else {
          b <- type_par
        }
        s <- Data@Misc[[x]]$interim$s
      } else {
        b <- s <- 0
      }
      
      Rec <- new("Rec")
      Rec@TAC <- Cref * (I_y + b * s)/(Iref + b * s)
      
      Rec@Misc <- Data@Misc[[x]]
      if(exists("next_assess_yr", inherits = FALSE)) Rec@Misc$interim$next_assess_yr <- next_assess_yr
      if(exists("Assess_diag_output", inherits = FALSE)) Rec@Misc$diagnostic <- Assess_diag_output$diagnostic
    }
    return(Rec)
  })
  
  custom_MP <- eval(call("function", as.pairlist(alist(x = 1, Data = , reps = 1)), MP_body))
  formals(custom_MP)$assessment_interval <- assessment_interval
  formals(custom_MP)$AddInd <- AddInd
  formals(custom_MP)$type <- type
  formals(custom_MP)$type_par <- type_par
  formals(custom_MP)$diagnostic <- diagnostic
  
  return(structure(custom_MP, class = "MP"))
}


interim_get_index <- function(xx, x, Data, y) {
  if(xx == "B") {
    Data@Ind[x, y]
  } else if(xx == "SSB") {
    Data@SpInd[x, y]
  } else if(xx == "VB") {
    Data@VInd[x, y]
  } else {
    Data@AddInd[x, suppressWarnings(as.numeric(xx)), y]
  }
}