#' What objects of this class are available
#'
#' Generic class finder
#'
#' Finds objects of the specified class in the global environment or in the
#' MSEtool and DLMtool packages. This function is an addendum to the \code{\link[DLMtool]{avail}}
#' function in DLMtool.
#'
#' @param classy A class of object (character string, e.g. 'Fleet')
#' @param all_avail Logical. If TRUE, function will return all objects of class \code{classy} available to user.
#' If FALSE, returns only those objects included in MSEtool.
#' @author Q. Huynh
#' @examples
#' avail("Assess")
#' avail("HCR")
#' avail("Stock")
#' avail("MP")
#' avail("MP", all_avail = FALSE)
#' @export
avail <- function(classy, all_avail = TRUE) {
  temp <- try(class(classy), silent = TRUE)
  if(class(temp) == "try-error") classy <- deparse(substitute(classy))
  if(temp == "function") classy <- deparse(substitute(classy))

  temp <- ls("package:MSEtool")[vapply(ls("package:MSEtool"), getclass, logical(1), classy = classy)]

  if(all_avail) {
    temp_globalenv <- ls(envir = .GlobalEnv)[vapply(ls(envir = .GlobalEnv), getclass, logical(1), classy = classy)]
    temp <- c(temp, temp_globalenv)

    temp_DLMtool <- try(DLMtool::avail(classy), silent = TRUE)
    if(!inherits(temp_DLMtool, "try-error")) temp <- unique(c(temp, temp_DLMtool))
  }

  if(length(temp) < 1) stop("No objects of class '", classy, "' found", call. = FALSE)
  return(temp)
}

getclass <- function(x, classy) any(inherits(get(x), classy))

#' Get the MSEtool vignettes
#'
#' A convenient function to open a web browser with the MSEtool package vignettes
#' @examples
#' \dontrun{
#' MSEtool::userguide()
#' DLMtool::userguide()
#' }
#' @seealso \link[DLMtool]{userguide}
#' @export
userguide <- function() {
  message("For the DLMtool user guide, type in \"DLMtool::userguide()\" to the console.")
  browseVignettes("MSEtool")
}


# Internal DLMtool functions that are also needed for MSEtool
iVB <- function(t0, K, Linf, L) max(1, ((-log(1 - L/Linf))/K + t0))
mconv <- function (m, sd) log(m) - 0.5 * log(1 + ((sd^2)/(m^2)))
squeeze <- function(x) (1 - .Machine$double.eps) * (x - 0.5) + 0.5

logit <- function(p, soft_bounds = TRUE, minp = 0.01, maxp = 0.99) {
  p <- squeeze(p)
  if(soft_bounds) {
    p <- pmax(minp, p)
    p <- pmin(maxp, p)
  }
  log(p/(1 - p))
}
ilogit <- function(x) 1/(1 + exp(-x))
ilogitm <- function(x) exp(x)/apply(exp(x), 1, sum)

#' @importFrom stats nlminb
optimize_TMB_model <- function(obj, control = list(), use_hessian = FALSE, restart = 1) {
  restart <- as.integer(restart)
  if(is.null(obj$env$random) && use_hessian) h <- obj$he else h <- NULL
  low <- rep(-Inf, length(obj$par))
  if(any(c("U_equilibrium", "F_equilibrium") %in% names(obj$par))) {
    low[match(c("U_equilibrium", "F_equilibrium"), names(obj$par))] <- 0
  }
  opt <- suppressWarnings(tryCatch(nlminb(obj$par, obj$fn, obj$gr, h, control = control, lower = low), error = as.character))
  if(is.character(opt) && all(is.na(obj$gr()))) {
    opt <- suppressWarnings(tryCatch(nlminb(obj$par, obj$fn, hessian = h, control = control, lower = low), error = as.character))
  }
  SD <- get_sdreport(obj, opt)

  if((is.character(SD) || !SD$pdHess) && !is.character(opt) && restart > 0) {
    obj$par <- opt$par
    Recall(obj, control, use_hessian, restart - 1)
  } else {
    res <- list(opt = opt, SD = SD)
    return(res)
  }
}

#' @importFrom stats optimHess
get_sdreport <- function(obj, opt) {
  if(is.character(opt)) par.fixed <- NULL else par.fixed <- opt$par
  if(is.null(obj$env$random) && !is.character(opt)) h <- obj$he(opt$par) else h <- NULL

  res <- tryCatch(sdreport(obj, par.fixed = par.fixed, hessian.fixed = h, getReportCovariance = FALSE), error = as.character)

  if(!is.character(res) && !is.character(opt) && !is.null(h) && !res$pdHess) {
    h <- optimHess(opt$par, obj$fn, obj$gr)
    res <- tryCatch(sdreport(obj, par.fixed = par.fixed, hessian.fixed = h, getReportCovariance = FALSE), error = as.character)
  }

  if(inherits(res, "sdreport") && res$pdHess && all(is.nan(res$cov.fixed))) {
    if(is.null(h)) h <- optimHess(opt$par, obj$fn, obj$gr)
    if(!is.character(try(chol(h), silent = TRUE))) res$cov.fixed <- chol2inv(chol(h))
  }

  if(inherits(res, "sdreport") && !is.null(par.fixed)) {
    obj2 <- MakeADFun(obj$env$data, obj$env$parameters, type = "ADFun",
                      ADreport = TRUE, DLL = obj$env$DLL, silent = obj$env$silent)
    gr <- obj2$gr(obj$env$last.par.best)
    if(any(is.na(gr))) {
      res$env$gradient.AD <- rep(NA_real_, nrow(gr))
    } else {
      inv_gr <- gr %>% pseudoinverse(tol = 1e-4)
      if(!is.null(obj$env$random)) inv_gr <- inv_gr[-obj$env$random, , drop = FALSE]
      res$env$gradient.AD <- colSums(inv_gr * as.vector(res$gradient.fixed))
    }
  }
  return(res)
}

sdreport_int <- function(object, select = c("all", "fixed", "random", "report"), p.value = FALSE, ...) {
  if(is.character(object)) return(object)
  select <- match.arg(select, several.ok = TRUE)
  if("report" %in% select) {
    gradient.AD <- object$env$gradient.AD %>% as.vector()
    if(is.null(gradient.AD)) gradient.AD <- rep(NA_real_, length(object$value))

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
    beta.par <- derive_beta_par(mu.beta.dist, sigma.beta.dist)
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

