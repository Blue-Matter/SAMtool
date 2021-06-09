
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
  
  OM <- new("OM", Stock = Stock, Fleet = Fleet, Obs = MSEtool::Precise_Unbiased, Imp = MSEtool::Perfect_Imp)
  OM@nsim <- 2
  OM@proyears <- 5
  OM@interval <- 2
  OM@maxF <- 3
  
  return(OM)
}

lapply2_fn <- function(x, i = 1, data_master, Hist) {
  
  OM <- make_OM(x, Hist)
  RCMdata <- new("RCMdata")
  
  if(data_master$cond[i] == "catch2") { # Condition on catch
    condition <- "catch2"
    RCMdata@Chist <- Hist@Data@Cat[x, ]
    OM@R0 <- 2 * OM@R0
  } else {
    condition <- "effort"
    RCMdata@Ehist <- Hist@Data@Cat[x, ]/Hist@Data@Ind[x, ] #Hist@TSdata$Find[x, ]
    if(data_master$cond[i] == "effort_with_catch") { # Conditioned on effort, but may include some catches
      RCMdata@Chist <- Hist@Data@Cat[x, ]
      #RCMdata@Chist[1:40, ] <- NA
    }
  }
  
  if(data_master$Index[i]) { # Only recent ten years index
    RCMdata@Index <- Hist@Data@Ind[1, ]
    #RCMdata@Index[c(1:20)] <- NA
    RCMdata@I_sd <- rep(Hist@SampPars$Obs$Isd[x], length(RCMdata@Index))
  }
  
  if(data_master$size[i] == "ML") { # Only recent mean length
    RCMdata@MS <- Hist@Data@ML[x, ]
    #args$data$MS[c(1:40)] <- NA
    RCMdata@MS_type <- "length"
  } else if(data_master$size[i] == "MW") {
    RCMdata@MS <- apply(Hist@Data@CAL[x, , ], 1, function(xx) {
      weighted.mean(x = Hist@Data@wla[x] * Hist@Data@CAL_mids ^ Hist@Data@wlb[x], w = xx, na.rm = TRUE)
    })
    #args$data$MS[c(1:40)] <- NA
    RCMdata@MS_type <- "weight"
  } else if(data_master$size[i] == "CAA") {
    RCMdata@CAA <- Hist@Data@CAA[x, , ]
    RCMdata@CAA_ESS <- rep(Hist@Data@Obs$CAA_ESS[x], length(Hist@Data@Year)) %>% pmin(50)
    #args$data$CAA[1:35, ] <- NA
  } else if(data_master$size[i] == "CAL") {
    RCMdata@CAL <- Hist@Data@CAL[x, , ]
    #args$data$CAL[1:35, ] <- NA
    RCMdata@length_bin <- Hist@Data@CAL_mids
    RCMdata@CAL_ESS <- rep(Hist@Data@Obs$CAL_ESS[x], length(Hist@Data@Year)) %>% pmin(50)
  }
  
  selectivity <- ifelse(mean(Hist@SampPars$Fleet$Vmaxlen_y[x, ]) < 1, 0, -1)
  s_selectivity <- -4
  
  do_check <- check_RCMdata(RCMdata, OM, condition)
  
  out <- SAMtool:::RCM_est(x = 1, RCMdata = do_check$RCMdata, selectivity = selectivity,
                           s_selectivity = s_selectivity, LWT = SAMtool:::make_LWT(list(), 1, 1),
                           comp_like = "multinomial", prior = SAMtool:::make_prior(list(), 1), 
                           StockPars = do_check$StockPars, FleetPars = do_check$FleetPars,
                           ObsPars = do_check$ObsPars)
  return(c(R0 = out$report$R0, dep = out$report$E[length(out$report$E)-1]/out$report$E0_SR, conv = as.numeric(out$SD$pdHess)))
}

lapply_fn <- function(x, data_master, Hist) {
  rr <- try(sapply(1:nrow(data_master), function(xx) lapply2_fn(x = x, i = xx, data_master = data_master, Hist = Hist)),
            silent = TRUE)
  return(rr)
}


# Run simulation
OM <- testOM %>% Replace(from = Precise_Unbiased)
OM@Prob_staying <- OM@Frac_area_1 <- OM@Size_area_1 <- rep(0.5, 2)

OM@nsim <- 100
Hist <- runMSE(OM, Hist = TRUE)

#### Function to run RCM in parallel (with list)
data_master <- expand.grid(cond = c("catch2", "effort", "effort_with_catch"),
                           Index = c(TRUE, FALSE), size = c("ML", "CAA", "CAL", "MW", "none")) %>%
  dplyr::filter(!(!Index & size == "none"))

setup(12)
sfExportAll()
sfLibrary(dplyr)

#### Run RCM
res <- sfClusterApplyLB(1:OM@nsim, lapply_fn, data_master = data_master, Hist = Hist)
sfStop()
#res <- suppressMessages(lapply(1, lapply_fn, data_master = data_master, Hist = Hist))
saveRDS(res, file = "tests/sim_RCM.rds")

res <- readRDS("tests/sim_RCM.rds")
#

out <- lapply(1:nrow(data_master), function(x) {
  res_x <- simplify2array(res)[, x, ]
  conv <- res_x[3, ] == 1
  
  R0 <- Hist@SampPars$Stock$R0[conv]
  D <- Hist@SampPars$Stock$D[conv]
  
  dev_R0 <- res_x[1, conv]/Hist@SampPars$Stock$R0[conv] - 1
  dev_D <- res_x[2, conv]/Hist@SampPars$Stock$D[conv] - 1
  
  bias_R0 <- mean(dev_R0)
  bias_D <- mean(dev_D)
  
  rmse_R0 <- sqrt(mean(dev_R0^2))
  rmse_D <- sqrt(mean(dev_D^2))
  
  data.frame(Value = c(bias_R0, bias_D, rmse_R0, rmse_D),
             Var = c("R0", "D", "R0", "D"),
             Type = c("bias", "bias", "rmse", "rmse"),
             combo = x, 
             no_catch = ifelse(rmse_R0 <= 1e-8, TRUE, FALSE))
}) %>% bind_rows()

out %>% #dplyr::filter(abs(Value) <= 1e5)
  ggplot(aes(combo, Value, shape = no_catch)) + geom_point() + 
  facet_grid(Type ~ Var, scales = "free_y") + theme_bw() +
  scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 16)) +
  coord_cartesian(ylim = c(-2, 2))
#



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