
max <- function(..., na.rm = TRUE) suppressWarnings(base::max(..., na.rm = na.rm))

arrows <- function(...) suppressWarnings(graphics::arrows(...))

hist.numeric <- function(x, ...) {
  if(all(!diff(x))) {
    breaks <- c(0.99, 1.01) * x[1]
    xlim <- c(0.8, 1.2) * x[1]
    graphics::hist.default(x, breaks = breaks, xlim = xlim, ...)
  } else {
    graphics::hist.default(x, ...)
  }
}

#' Get the SAMtool vignettes
#'
#' A convenient function to open a web browser with the SAMtool package vignettes
#' @examples
#' userguide()
#' 
#' @return Displays a browser webpage to the package vignette.
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

tiny_comp <- function(x) {
  all_zero <- all(is.na(x)) | sum(x, na.rm = TRUE) == 0
  if(!all_zero) {
    x_out <- x/sum(x, na.rm = TRUE)
    ind <- is.na(x) | x == 0
    if(any(ind)) x_out[ind] <- 1e-8
  } else {
    x_out <- x
  }
  return(x_out)
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
  restart <- as.integer(restart)
  if(is.null(obj$env$random) && use_hessian) h <- obj$he else h <- NULL
  low <- rep(-Inf, length(obj$par))
  if(any(c("U_equilibrium", "F_equilibrium") %in% names(obj$par))) {
    low[match(c("U_equilibrium", "F_equilibrium"), names(obj$par))] <- 0
  }
  opt <- tryCatch(suppressWarnings(nlminb(obj$par, obj$fn, obj$gr, h, control = control, lower = low)),
                  error = function(e) as.character(e))
  SD <- get_sdreport(obj)

  if(!SD$pdHess && restart > 0) {
    if(!is.character(opt)) obj$par <- opt$par * exp(rnorm(length(opt$par), 0, 1e-3))
    Recall(obj, control, use_hessian, restart - 1)
  } else {
    return(list(opt = opt, SD = SD))
  }
}


check_det <- function(h, abs_val = 0.1, is_null = TRUE) {
  if(is.null(h)) return(is_null)
  det_h <- det(h) %>% abs()
  !is.na(det_h) && det_h < abs_val
}

#' @importFrom corpcor pseudoinverse
get_sdreport <- function(obj, getReportCovariance = FALSE) {
  par.fixed <- obj$env$last.par.best
  if(is.null(obj$env$random)) {
    h <- obj$he(par.fixed)
    if(any(is.na(h)) || any(is.infinite(h)) || det(h) <= 0) {
      h <- NULL
    } else {
      res <- suppressWarnings(sdreport(obj, par.fixed = par.fixed, hessian.fixed = h, getReportCovariance = getReportCovariance))
      if(!res$pdHess) h <- NULL
    }
  } else {
    par.fixed <- par.fixed[-obj$env$random] 
    h <- NULL
  }
  
  if(check_det(h)) {  # If hessian doesn't exist or marginal positive-definite cases, with -0.1 < det(h) <= 0
    h <- suppressWarnings(optimHess(par.fixed, obj$fn, obj$gr))
    res <- suppressWarnings(sdreport(obj, par.fixed = par.fixed, hessian.fixed = h, getReportCovariance = getReportCovariance))
  }
  
  if(check_det(h) && !res$pdHess && requireNamespace("numDeriv", quietly = TRUE)) {
    h <- numDeriv::jacobian(obj$gr, par.fixed)
    res <- suppressWarnings(sdreport(obj, par.fixed = par.fixed, hessian.fixed = h, getReportCovariance = getReportCovariance))
  }
  
  if(all(is.na(res$cov.fixed)) && res$pdHess) {
    if(!is.character(try(chol(h), silent = TRUE))) res$cov.fixed <- chol2inv(chol(h))
  }

  obj2 <- MakeADFun(obj$env$data, obj$env$parameters, type = "ADFun", ADreport = TRUE, 
                    DLL = obj$env$DLL, silent = obj$env$silent)
  gr <- obj2$gr(obj$env$last.par.best)
  if(any(is.na(gr))) {
    res$env$gradient.AD <- rep(NA_real_, nrow(gr))
  } else {
    inv_gr <- try(pseudoinverse(gr, tol = 1e-4), silent = TRUE)
    if(is.character(inv_gr)) {
      res$env$gradient.AD <- rep(NA_real_, nrow(gr))
    } else {
      if(!is.null(obj$env$random)) inv_gr <- inv_gr[-obj$env$random, , drop = FALSE]
      res$env$gradient.AD <- colSums(inv_gr * as.vector(res$gradient.fixed))
    }
  }
  
  res$env$corr.fixed <- suppressWarnings(cov2cor(res$cov.fixed) %>% round(3)) %>% 
    structure(dimnames = list(names(res$par.fixed), names(res$par.fixed)))
  
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
  rownames(out) <- make_unique_names(rownames(out))
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
    I_hist[!is.na(I_hist) & I_hist <= 0] <- NA_real_
    if(all(is.na(I_sd)) || all(!I_sd, na.rm = TRUE)) {
      if(!is.null(Data@Obs$Isd[x])) {
        I_sd <- rep(Data@Obs$Isd[x], length(yind))
      } else {
        I_sd <- rep(0.2, length(yind))
      }
    }
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


make_prior <- function(prior, nsurvey, SR_rel, dots = list(), msg = TRUE) { # log_R0, log_M, h, q
  if(length(prior) == 0 && !is.null(dots$priors)) prior <- dots$priors
  
  no_index <- nsurvey == 0
  if(no_index) nsurvey <- 1 # Use only on next two lines
  use_prior <- rep(0L, nsurvey + 3)
  pr_matrix <- matrix(NA_real_, nsurvey + 3, 2) %>% 
    structure(dimnames = list(c("log_R0", "h", "log_M", paste0("q_", 1:nsurvey)), c("par1", "par2")))
  
  if(!is.null(prior$R0)) {
    if(msg) message("Prior for log_R0 found.")
    use_prior[1] <- 1L
    pr_matrix[1, ] <- c(log(prior$R0[1]), prior$R0[2])
  }
  if(!is.null(prior$h)) {
    if(msg) message("Prior for steepness (h) found.")
    use_prior[2] <- 1L
    if(SR_rel == 1) { #BH
      a <- MSEtool::alphaconv(1.25 * prior$h[1] - 0.25, 1.25 * prior$h[2])
      b <- MSEtool::betaconv(1.25 * prior$h[1] - 0.25, 1.25 * prior$h[2])
      
      if(a <= 0) stop("The alpha parameter < 0 (beta distribution) for the steepness prior. Try reducing the prior SD.", call. = FALSE)
      if(b <= 0) stop("The beta parameter < 0 (beta distribution) for the steepness prior. Try reducing the prior SD.", call. = FALSE)
      pr_matrix[2, ] <- c(a, b)
    } else { #Ricker
      pr_matrix[2, ] <- prior$h
    }
  }
  if(!is.null(prior$M)) {
    if(msg) message("Prior for log_M found.")
    use_prior[3] <- 1L
    pr_matrix[3, ] <- c(log(prior$M[1]), prior$M[2])
  }
  if(!no_index && !is.null(prior$q)) {
    if(msg) message("Prior for q found.")
    if(nsurvey == 1) {
      if(is.matrix(prior$q)) {
        pr_matrix[4, ] <- prior$q[1, ]
      } else {
        pr_matrix[4, ] <- prior$q
      }
      use_prior[4] <- !is.na(pr_matrix[4, 1])
    } else if(is.matrix(prior$q) && nrow(prior$q) == nsurvey) {
      pr_matrix[4:(4+nsurvey-1), ] <- prior$q[1:nsurvey, ]
      use_prior[4:(4+nsurvey-1)] <- !is.na(pr_matrix[4:(4+nsurvey-1), 1])
    } else {
      stop("prior$q should be a matrix of ", nsurvey, "rows and 2 columns.")
    }
    if(any(pr_matrix[4:(4+nsurvey-1), 1] <= 0, na.rm = TRUE)) stop("Ensure q prior mean is greater than > 0.")
    pr_matrix[4:(4+nsurvey-1), 1] <- log(pr_matrix[4:(4+nsurvey-1), 1])
  }
  return(list(use_prior = use_prior, pr_matrix = pr_matrix))
}
