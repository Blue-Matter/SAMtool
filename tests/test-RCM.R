
#### For data-limited situations - see notes below

library(SAMtool)
#library(testthat)

#### Functions

make_OM <- function(x, Hist) {
  require(dplyr)
  
  Stock <- new("Stock")
  Stock@maxage <- Hist@SampPars$Stock$maxage
  Stock@R0 <- Hist@SampPars$Stock$R0[x]
  Stock@M <- Hist@SampPars$Stock$M[x] %>% rep(2)
  Stock@h <- Hist@SampPars$Stock$h[x] %>% rep(2)
  Stock@SRrel <- 1L
  Stock@Perr <- Hist@SampPars$Stock$procsd[x] %>% rep(2)
  Stock@AC <- c(0, 0)
  Stock@Linf <- Hist@SampPars$Stock$Linf[x] %>% rep(2)
  Stock@K <- Hist@SampPars$Stock$K[x] %>% rep(2)
  Stock@t0 <- Hist@SampPars$Stock$t0[x] %>% rep(2)
  Stock@LenCV <- Hist@SampPars$Stock$LenCV[x] %>% rep(2)
  Stock@a <- Hist@SampPars$Stock$a
  Stock@b <- Hist@SampPars$Stock$b
  Stock@Msd <- Stock@Linfsd <- Stock@Ksd <- c(0, 0)
  Stock@Size_area_1 <- Stock@Frac_area_1 <- Stock@Prob_staying <- c(0.5, 0.5)
  Stock@Fdisc <- c(0, 0)
  Stock@L50 <- Hist@SampPars$Stock$L50[x] %>% rep(2)
  Stock@L50_95 <- Hist@SampPars$Stock$L50_95[x] %>% rep(2)
  
  # Placeholders
  Stock@D <- c(0.9, 0.9)
  
  Fleet <- new("Fleet")
  Fleet@CurrentYr <- Fleet@nyears <- Hist@OM@nyears
  Fleet@EffYears <- 1:Fleet@nyears
  Fleet@EffLower <- rep(1e-4, Fleet@nyears)
  Fleet@EffUpper <- rep(1e-3, Fleet@nyears)
  Fleet@Esd <- Fleet@qinc <- Fleet@qcv <- c(0, 0)
  
  Fleet@L5 <- Hist@SampPars$Fleet$L5_y[x, 1] %>% rep(2)
  Fleet@LFS <- Hist@SampPars$Fleet$LFS_y[x, 1] %>% rep(2)
  Fleet@Vmaxlen <- Hist@SampPars$Fleet$Vmaxlen_y[x, 1] %>% rep(2)
  Fleet@isRel <- FALSE
  Fleet@DR <- c(0, 0)
  Fleet@Spat_targ <- c(1, 1)
  Fleet@MPA <- FALSE
  
  OM <- new("OM", Stock = Stock, Fleet = Fleet, Obs = MSEtool::Generic_Obs, Imp = MSEtool::Perfect_Imp)
  OM@nsim <- 2
  OM@proyears <- 5
  OM@interval <- 2
  OM@maxF <- 3
  
  return(OM)
}

lapply2_fn <- function(x, i = 1, data_master, Hist) {
  
  OM <- make_OM(x, Hist)
  RCMdata <- new("RCMdata")
  
  args <- list()
  
  if(data_master$cond[i] == "catch") { # Condition on catch
    args$condition <- "catch2"
    RCMdata@Chist <- Hist@Data@Cat[x, ]
  } else {
    if(data_master$Catch_with_effort[i]) { # Conditioned on effort, but may include some catches
      RCMdata@Chist <- Hist@Data@Cat[x, ]
      #RCMdata@Chist[1:40, ] <- NA
    }
    args$condition <- "effort"
    RCMdata@Ehist <- Hist@TSdata$Find[x, ]
  }
  
  if(data_master$Index[i]) { # Only recent ten years index
    RCMdata@Index <- Hist@Data@Ind[1, ]
    #RCMdata@Index[c(1:20)] <- NA
    RCMdata@I_sd <- rep(Hist@SampPars$Obs$Isd[x], length(RCMdata@Index))
    args$s_selectivity <- "B"
  }
  
  if(data_master$ML[i]) { # Only recent mean length
    RCMdata@MS <- Hist@Data@ML[x, ]
    #args$data$MS[c(1:40)] <- NA
    RCMdata@MS_type <- "length"
  } else if(data_master$MW[i]) {
    WAL <- Hist@Data@wla[x] * Hist@Data@CAL_mids^Hist@Data@wlb[x]
    
    MW <- colSums(t(Hist@Data@CAL[x, , ]) * WAL) / rowSums(Hist@Data@CAL[x, , ])
    RCMdata@MS <- MW
    RCMdata@MS_type <- "weight"
  } else if(data_master$CAA[i]) {
    RCMdata@CAA <- Hist@Data@CAA[x, , ]
    #args$data$CAA[1:35, ] <- NA
  } else if(data_master$CAL[i]) {
    RCMdata@CAL <- Hist@Data@CAL[x, , ]
    #args$data$CAL[1:35, ] <- NA
    RCMdata@length_bin <- Hist@Data@CAL_mids
  }
  args$mean_fit <- TRUE
  
  args$OM <- OM
  args$data <- RCMdata
  args$selectivity <- ifelse(mean(Hist@SampPars$Fleet$Vmaxlen_y[x, ]) < 1, "dome", "logistic")
  
  SRA <- do.call(RCM, args)
  
  return(c(R0 = SRA@mean_fit$report$R0, dep = SRA@OM@cpars$D %>% mean()))
}

lapply_fn <- function(x, data_master, Hist) {
  rr <- try(lapply(1:nrow(data_master), function(xx) lapply2_fn(x = x, i = xx, data_master = data_master, Hist = Hist)),
            silent = TRUE)
  return(rr)
}





# Run simulation
OM <- testOM
OM@nsim <- 100
Hist <- runMSE(OM, Hist = TRUE)

bool <- c(TRUE, FALSE)

#### Function to run RCM in parallel (with list)
data_master <- expand.grid(cond = c("catch2", "effort"), Catch_with_effort = bool,
                           Index = bool, ML = bool, MW = bool, CAL = bool, CAA = bool)
data_master$MW[data_master$ML] <- FALSE
data_master$CAA[data_master$ML | data_master$MW] <- data_master$CAL[data_master$ML | data_master$MW] <- FALSE

#data_master <- data_master[1:3, ]

setup(8)
sfExportAll()
sfLibrary(dplyr)

#### Run RCM
res <- sfClusterApplyLB(1:OM@nsim, lapply_fn, data_master = data_master, Hist = Hist)
#res <- lapply(1:2, lapply_fn, data_master = data_master, Hist = Hist)

saveRDS(res, file = "sim_RCM.rds")

#res <- lapply(5, lapply_fn, data_master = data_master, Hist = Hist, OM = OM)
#
## Test MW
#res <- lapply_fn(9, data_master = data_master, Hist = Hist, OM = OM)
#
#
##### Test whether plot function works
#lapply_fn_plot <- function(x, res) {
#  try(R.utils::withTimeout(plot(res[[x]], compare = FALSE, open_file = FALSE, filename = as.character(x)),
#                           timeout = 5),
#      silent = TRUE)
#}
#
##sfLibrary(R.utils)
##sfExport(list = c("lapply_fn_plot", "res"))
#test_plot <- lapply(1:nrow(data_master), lapply_fn_plot, res = res)
#saveRDS(test_plot, file = "tests/test_plot.rds")
#
#error_fn <- function(x) {
#  length(grep("elapsed time limit", x[1])) > 0 | length(grep("C:/", x[1])) > 0
#}
#
#
#
#data_master$conv <- lapply(res, function(x) sum(x@conv)/length(x@conv))
#data_master$figure <- vapply(test_plot, error_fn, numeric(1))
#View(data_master)
#
#