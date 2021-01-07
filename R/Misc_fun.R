

TACfilter <- function(TAC) {
  TAC[TAC < 0] <- NA_real_ 
  TAC[TAC > (mean(TAC, na.rm = TRUE) + 5 * sd(TAC, na.rm = TRUE))] <- NA_real_  # remove very large TAC samples
  return(TAC)
}

#' Get the SAMtool vignettes
#'
#' A convenient function to open a web browser with the SAMtool package vignettes
#' @examples
#' \dontrun{
#' userguide()
#' }
#' @export
userguide <- function() browseVignettes("SAMtool")


squeeze <- function(x) (1 - .Machine$double.eps) * (x - 0.5) + 0.5
iVB <- function(t0, K, Linf, L) max(1, ((-log(1 - L/Linf))/K + t0))  # Inverse Von-B

logit <- function(p, soft_bounds = TRUE, minp = 0.01, maxp = 0.99) { #log(p/(1 - p))
  p <- squeeze(p)
  if(soft_bounds) {
    p <- pmax(p, minp)
    p <- pmin(p, maxp)
  }
  qlogis(p)
}

ilogit <- function(x) plogis(x) #1/(1 + exp(-x))

ilogitm <- function(x) {
  if(inherits(x, "matrix")) {
    return(exp(x)/apply(exp(x), 1, sum))
  } else {
    return(exp(x)/sum(exp(x)))
  }
}

ilogit2 <- function(x, ymin = 0, ymax = 1, y0 = 0.5, scale = 1) {
  location <- scale * log((ymax - ymin)/(y0 - ymin) - 1)
  return((ymax - ymin) * plogis(x, location, scale) + ymin)
}

logit2 <- function(v, ymin = 0, ymax = 1, y0 = 0.5, scale = 1) {
  location <- scale * log((ymax - ymin)/(y0 - ymin) - 1)
  p <- (v - ymin)/(ymax - ymin)
  return(qlogis(p, location, scale))
}

get_F01 <- function(FM, YPR) {
  if(is.null(FM) && is.null(YPR)) stop("F01 can not be used.")
  stopifnot(length(FM) == length(YPR))
  dY_dF <- (YPR[2:length(YPR)] - YPR[2:length(YPR) - 1])/(FM[2:length(FM)] - FM[2:length(FM) - 1])
  LinInterp(dY_dF, FM[-length(FM)], xlev = 0.1 * dY_dF[1])
}

get_Fmax <- function(FM, YPR) {
  if(is.null(FM) && is.null(YPR)) stop("Fmax can not be used.")
  FM[which.max(YPR)]
}

get_FSPR <- function(FM, SPR, target = 0.4) {
  if(is.null(FM) && is.null(SPR)) stop("SPR can not be used.")
  stopifnot(length(FM) == length(SPR))
  LinInterp(SPR, FM, xlev = target)
}

LinInterp <- function(x,y,xlev,ascending=F,zeroint=F){
  if(zeroint){
    x<-c(0,x)
    y<-c(0,y)
  }
  if(ascending){
    cond<-(1:length(x))<which.max(x)
  }else{
    cond<-rep(TRUE,length(x))
  }
  
  close<-which.min((x[cond]-xlev)^2)
  ind<-c(close,close+(x[close]<xlev)*2-1)
  ind <- ind[ind <= length(x)]
  if (length(ind)==1) ind <- c(ind, ind-1)
  ind<-ind[order(ind)]
  pos<-(xlev-x[ind[1]])/(x[ind[2]]-x[ind[1]])
  y[ind[1]]+pos*(y[ind[2]]-y[ind[1]])
  
}


optimize_TMB_model <- function(obj, control = list(), use_hessian = FALSE, restart = 1) {
  
  old_warn <- options()$warn
  options(warn = -1)
  on.exit(options(warn = old_warn))
  
  restart <- as.integer(restart)
  if(is.null(obj$env$random) && use_hessian) h <- obj$he else h <- NULL
  low <- rep(-Inf, length(obj$par))
  if(any(c("U_equilibrium", "F_equilibrium") %in% names(obj$par))) {
    low[match(c("U_equilibrium", "F_equilibrium"), names(obj$par))] <- 0
  }
  opt <- try(nlminb(obj$par, obj$fn, obj$gr, h, control = control, lower = low), silent = TRUE)
  if(is.character(opt) && all(is.na(obj$gr()))) {
    opt <- try(nlminb(obj$par, obj$fn, hessian = h, control = control, lower = low), silent = TRUE)
  }
  SD <- get_sdreport(obj, opt)

  if(restart > 0 && !is.character(SD) && !SD$pdHess) {
    if(!is.character(opt)) obj$par <- opt$par
    obj$par <- obj$par * exp(rnorm(length(obj$par), 0, 1e-3))
    Recall(obj, control, use_hessian, restart - 1)
  } else {
    return(list(opt = opt, SD = SD))
  }
}


get_sdreport <- function(obj, opt) {
  
  old_warn <- options()$warn
  options(warn = -1)
  on.exit(options(warn = old_warn))
  
  if(is.character(opt)) par.fixed <- NULL else par.fixed <- opt$par
  if(is.null(obj$env$random) && !is.character(opt)) h <- obj$he(opt$par) else h <- NULL

  res <- try(sdreport(obj, par.fixed = par.fixed, hessian.fixed = h, getReportCovariance = FALSE), silent = TRUE)

  if(!is.character(res) && !res$pdHess && !is.character(opt) && !is.null(h)) {
    h <- optimHess(opt$par, obj$fn, obj$gr)
    res <- try(sdreport(obj, par.fixed = par.fixed, hessian.fixed = h, getReportCovariance = FALSE), silent = TRUE)
  }

  if(!is.character(res) && res$pdHess && all(is.nan(res$cov.fixed))) {
    if(is.null(h)) h <- optimHess(opt$par, obj$fn, obj$gr)
    if(!is.character(try(chol(h), silent = TRUE))) res$cov.fixed <- chol2inv(chol(h))
  }

  if(!is.character(res) && !is.null(par.fixed)) {
    obj2 <- MakeADFun(obj$env$data, obj$env$parameters, type = "ADFun",
                      ADreport = TRUE, DLL = obj$env$DLL, silent = obj$env$silent)
    gr <- obj2$gr(obj$env$last.par.best)
    if(any(is.na(gr))) {
      res$env$gradient.AD <- rep(NA_real_, nrow(gr))
    } else {
      inv_gr <- try(gr %>% pseudoinverse(tol = 1e-4), silent = TRUE)
      if(is.character(inv_gr)) {
        res$env$gradient.AD <- rep(NA_real_, nrow(gr))
      } else {
        if(!is.null(obj$env$random)) inv_gr <- inv_gr[-obj$env$random, , drop = FALSE]
        res$env$gradient.AD <- colSums(inv_gr * as.vector(res$gradient.fixed))
      }
    }
  }
  if(!is.character(res)) {
    res$env$corr.fixed <- res$cov.fixed %>% cov2cor() %>% round(3) %>% 
      structure(dimnames = list(names(res$par.fixed), names(res$par.fixed)))
  }
  return(res)
}

sdreport_int <- function(object, select = c("all", "fixed", "random", "report"), p.value = FALSE, ...) {
  if(is.character(object)) return(object)
  select <- match.arg(select, several.ok = TRUE)
  if("report" %in% select) {
    gradient.AD <- object$env$gradient.AD %>% as.vector()
    if(is.null(gradient.AD)) {
      gradient.AD <- rep(NA_real_, length(object$value))
    } else {
      gradient.AD <- ifelse(object$sd, gradient.AD, 0)
    }
    AD <- TMB::summary.sdreport(object, "report", p.value = p.value) %>% cbind("Gradient" = gradient.AD)
  } else AD <- NULL

  if("fixed" %in% select) {
    fix <- TMB::summary.sdreport(object, "fixed", p.value = p.value) %>% cbind("Gradient" = as.vector(object$gradient.fixed))
  } else fix <- NULL

  if(!is.null(object$par.random) && "random" %in% select) {
    random <- TMB::summary.sdreport(object, "random", p.value = p.value) %>% cbind("Gradient" = rep(NA_real_, length(object$par.random)))
  } else {
    random <- NULL
  }

  out <- rbind(AD, fix, random)
  out <- cbind(out, "CV" = ifelse(abs(out[, "Estimate"]) > 0, out[, "Std. Error"]/abs(out[, "Estimate"]), NA_real_))
  return(out)
}

# Call from inside generate_plots() and summary.Assessment
assign_Assessment_slots <- function(Assessment = NULL) {
  if(is.null(Assessment)) Assessment <- get("Assessment", envir = parent.frame(), inherits = FALSE)
  Nslots <- length(slotNames(Assessment))
  for(i in 1:Nslots) {
    assign(slotNames(Assessment)[i], slot(Assessment, slotNames(Assessment)[i]), envir = parent.frame())
  }
  invisible()
}

# For SCA, if there are fewer years of CAA/CAL than Year, add NAs to matrix
expand_comp_matrix <- function(Data, comp_type = c("CAA", "CAL")) {
  comp_type <- match.arg(comp_type)
  ny <- length(Data@Year)

  comp <- slot(Data, comp_type)
  dim_comp <- dim(comp)
  ny_comp <- dim_comp[2]
  if(ny_comp < ny) {
    newcomp <- array(NA, dim = c(1, ny, dim_comp[3]))
    ind_new <- ny - ny_comp + 1
    newcomp[ , ind_new:ny, ] <- comp
    slot(Data, comp_type) <- newcomp
  }
  if(ny_comp > ny) {
    ind_new <- ny_comp - ny + 1
    newcomp <- comp[, ind_new:ny, ]
    slot(Data, comp_type) <- newcomp
  }
  return(Data)
}

# var_div - report variables which are divided by the catch rescale
# var_mult - report variables which are multiplied by the catch rescale
# var_trans - transformed variables which need to be rescaled
# fun_trans - the function for rescaling the transformed variables (usually either "*" or "/")
# fun_fixed - the transformation from the output variable to the estimated variable indicated in var_trans (e.g. log, logit, NULL)
rescale_report <- function(var_div, var_mult, var_trans = NULL, fun_trans = NULL, fun_fixed = NULL) {
  output <- mget(c("report", "rescale", "SD"), envir = parent.frame(), ifnotfound = list(NULL))
  report <- output$report

  if(!is.null(var_div)) report[var_div] <- lapply(report[var_div], "/", output$rescale)
  if(!is.null(var_mult)) report[var_mult] <- lapply(report[var_mult], "*", output$rescale)
  assign("report", report, envir = parent.frame())

  if(!is.null(output$SD) && !is.character(output$SD)) {
    SD <- output$SD
    if(!is.null(var_trans)) {
      for(i in 1:length(var_trans)) {
        var_trans2 <- var_trans[i]
        fun_trans2 <- fun_trans[i]
        fun_fixed2 <- fun_fixed[i]

        ind <- var_trans2 == names(SD$value)
        if(any(ind)) {
          SD$value[ind] <- do.call(match.fun(fun_trans2), list(SD$value[ind], output$rescale))
          SD$sd[ind] <- do.call(match.fun(fun_trans2), list(SD$sd[ind], output$rescale))
        }
        if(!is.na(fun_fixed2)) {
          fixed_name <- paste0(fun_fixed2, "_", var_trans2)
          ind_fixed <- fixed_name == names(SD$par.fixed)
          if(any(ind_fixed)) SD$par.fixed[ind_fixed] <- do.call(match.fun(fun_fixed2), list(SD$value[var_trans2]))
        }
      }
    }
    assign("SD", SD, envir = parent.frame())
  }
  invisible()
}


sample_steepness3 <- function(n, mu, cv, SR_type = c("BH", "Ricker")) {
  if(n == 1) {
    return(mu)
  } else if(SR_type == "BH") {
    sigma <- mu * cv
    mu.beta.dist <- (mu - 0.2)/0.8
    sigma.beta.dist <- sigma/0.8
    beta.par <- MSEtool::derive_beta_par(mu.beta.dist, sigma.beta.dist)
    h.transformed <- rbeta(n, beta.par[1], beta.par[2])
    h <- 0.8 * h.transformed + 0.2
    h[h > 0.99] <- 0.99
    h[h < 0.2] <- 0.2
    return(h)
  } else {
    sigma <- mu * cv
    mu.lognorm.dist <- mconv(mu - 0.2)
    sigma.lognorm.dist <- sigma

    h.transformed <- trlnorm(n, mconv(mu.lognorm.dist, sigma.lognorm.dist), sdconv(mu.lognorm.dist, sigma.lognorm.dist))
    h <- h.transformed + 0.2
    h[h < 0.2] <- 0.2
    return(h)
  }
}

Assess_I_hist <- function(xx, Data, x, yind) {
  if(xx == 0 || xx == "B") {

    I_hist <- Data@Ind[x, yind]
    I_sd <- sdconv(1, Data@CV_Ind[x, yind])
    I_units <- 1L

  } else if(xx == "SSB" && .hasSlot(Data, "SpInd")) {

    I_hist <- Data@SpInd[x, yind]
    I_sd <- sdconv(1, Data@CV_SpInd[x, yind])
    I_units <- 1L

  } else if(xx == "VB" && .hasSlot(Data, "VInd")) {

    I_hist <- Data@VInd[x, yind]
    I_sd <- sdconv(1, Data@CV_VInd[x, yind])
    I_units <- 1L

  } else if(is.numeric(xx) && xx > 0 && .hasSlot(Data, "AddInd") && xx <= dim(Data@AddInd)[2]) {

    I_hist <- Data@AddInd[x, xx, yind]
    I_sd <- sdconv(1, Data@CV_AddInd[x, xx, yind])
    if(.hasSlot(Data, "AddIunits") && !is.na(Data@AddIunits[xx])) {
      I_units <- Data@AddIunits[xx]
    } else I_units <- 1L

  }

  if(exists("I_hist", inherits = FALSE)) {
    I_hist[I_hist <= 0] <- NA
  } else {
    I_hist <- I_sd <- I_units <- NULL
  }
  return(list(I_hist = I_hist, I_sd = I_sd, I_units = I_units))
}


dev_AC <- function(n, mu = 1, stdev, AC, seed, chain_start) {
  if(!missing(seed)) set.seed(seed)
  
  log_mean <- log(mu) - 0.5 * stdev^2 * (1 - AC/sqrt(1 - AC^2)) #http://dx.doi.org/10.1139/cjfas-2016-0167
  samp <- rnorm(n, log_mean, stdev)
  out <- numeric(n)
  if(missing(chain_start)) {
    out[1] <- samp[1]
  } else {
    out[1] <- chain_start * AC + samp[1] * sqrt(1 - AC^2)
  }
  for(i in 2:n) out[i] <- out[i-1] * AC + samp[i] * sqrt(1 - AC^2)
  return(out)
}
