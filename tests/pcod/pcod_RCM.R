


library(SAMtool)

make_pcod_OM <- function(nsim = 48, proyears = 50, FYr = 1956, maxage = 10, stochastic = FALSE) {
  
  #rho <- 0.059
  #theta2 <- 1.471
  #sigmaI <- sqrt(rho/theta2) # 0.2
  #
  #sigmaR <- sqrt((1 - rho)/theta2) # 0.8
  
  Stock <- new("Stock")
  Stock@maxage <- maxage
  Stock@R0 <- 12000
  Stock@M <- c(0.31, 0.31) # Remember to add a prior to RCM
  Stock@h <- c(0.73, 0.73) # Remember to add a prior to RCM
  Stock@SRrel <- 1L
  Stock@Perr <- c(0.8, 0.8)
  Stock@AC <- c(0, 0)
  Stock@Linf <- c(95.51, 95.51)
  Stock@K <- c(0.19, 0.19)
  Stock@t0 <- c(-0.81, -0.81)
  Stock@LenCV <- c(0.1, 0.1) # Not available
  Stock@a <- 6.72e-6
  Stock@b <- 3.11
  Stock@Msd <- Stock@Linfsd <- Stock@Ksd <- c(0, 0)
  Stock@Size_area_1 <- Stock@Frac_area_1 <- Stock@Prob_staying <- c(0.5, 0.5)
  Stock@Fdisc <- c(0, 0)
  
  # Placeholders
  Stock@D <- c(0.9, 0.9)
  Stock@L50 <- Stock@L50_95 <- c(0, 0)
  
  Fleet <- new("Fleet")
  Fleet@CurrentYr <- 2020
  
  Fleet@nyears <- length(FYr:Fleet@CurrentYr)
  Fleet@EffYears <- FYr:Fleet@CurrentYr
  Fleet@EffLower <- rep(1e-4, Fleet@nyears)
  Fleet@EffUpper <- rep(1e-3, Fleet@nyears)
  Fleet@Esd <- Fleet@qinc <- Fleet@qcv <- c(0, 0)
  
  Fleet@L5 <- c(20, 20)
  Fleet@LFS <- c(30, 30)
  Fleet@Vmaxlen <- c(1, 1)
  Fleet@isRel <- FALSE
  Fleet@DR <- c(0, 0)
  Fleet@Spat_targ <- c(1, 1)
  Fleet@MPA <- FALSE
  
  OM <- new("OM", Stock = Stock, Fleet = Fleet, Obs = Generic_Obs, Imp = Perfect_Imp)
  OM@nsim <- nsim
  OM@proyears <- proyears
  OM@interval <- 2
  OM@maxF <- 3
  
  if(stochastic) {
    set.seed(24)
    OM@cpars$M <- rlnorm(nsim, log(0.31), sdconv(1, 0.1))
    
    h_a <- alphaconv((0.73 - 0.2)/0.8, 0.12/0.8)
    h_b <- betaconv((0.73 - 0.2)/0.8, 0.12/0.8)
    #median(rbeta(1e5, h_a, h_b) * 0.8 + 0.2)
    OM@cpars$h <- rbeta(nsim, h_a, h_b) * 0.8 + 0.2
  }
  
  # Maturity = selectivity
  Mat_age <- ifelse(0:OM@maxage >= 2, 1, 0)
  
  OM@cpars$Mat_age <- OM@cpars$V <- 
    array(Mat_age, c(OM@maxage + 1, OM@nyears + OM@proyears, OM@nsim)) %>% aperm(c(3, 1, 2))
  
  return(OM)
}


OM <- make_pcod_OM(nsim = 2)
OM <- make_pcod_OM(nsim = 48, stochastic = TRUE)
mat_ogive <- OM@cpars$Mat_age[1,,1]

pcod_data <- readRDS("data-raw/pcod-data-5ABCD.rds")
RCM_data <- pcod_data[c(1, 3, 4, 6, 7, 8, 18, 23)]
RCM_data$C_sd <- matrix(0.05, OM@nyears, 1)
RCM_data$MS_cv <- 0.2

survey_name <- c("HSMSAS", "SYN QCS", "SYN HS", "CPUE Hist", "CPUE Modern")

             # mu q      # sd q
#prior_q <- c(0.99908020, 0.25323739,
#             0.40820819, 0.12518667,
#             0.06548742, 0.02011711,
#             1.00009286, 0.25413549,
#             0.99990984, 0.25414730) %>% matrix(5, 2, byrow = TRUE) %>% 
#  structure(dimnames = list(survey_name, c("mu", "sd")))

prior_q <- c(NA, NA,
             0.40820819, 0.3, # in 2020
             0.06548742, 0.3, # in 2020
             NA, NA,
             NA, NA) %>% 
  matrix(5, 2, byrow = TRUE) %>% 
  structure(dimnames = list(survey_name, c("mu", "sd")))
prior <- list( #M = c(0.5, 0.1), h = c(0.7, 0.15), 
              q = prior_q)

out <- RCM(OM, RCM_data, 
           condition = "catch", mean_fit = TRUE, #resample = TRUE,
           # Next five lines to fix selectivity = maturity,
           selectivity = "free", s_selectivity = rep("SSB", ncol(pcod_data$Index)),
           vul_par = matrix(mat_ogive, length(mat_ogive), 1),
           map_vul_par = matrix(NA, length(mat_ogive), 1),
           #map_log_early_rec_dev = 1:OM@maxage,
           map_log_early_rec_dev = rep(1, OM@maxage), # Estimate R-init
           prior = prior)

saveRDS(out, file = "tests/pcod/pcod-RCM2.rds")

out <- readRDS("tests/pcod/pcod-RCM2.rds")
summary(out@mean_fit$SD) # Parameter estimates
plot(1956:2021, out@mean_fit$report$E, typ = 'o', xlab = "Year", ylab = "SSB")
abline(h = out@mean_fit$report$E0_SR) # SSB0
out@mean_fit$report$CR


out <- RCM(OM, RCM_data, 
           condition = "catch", mean_fit = TRUE, #resample = TRUE,
           # Next five lines to fix selectivity = maturity,
           selectivity = "free", s_selectivity = rep("SSB", ncol(pcod_data$Index)),
           vul_par = matrix(mat_ogive, length(mat_ogive), 1),
           map_vul_par = matrix(NA, length(mat_ogive), 1),
           #map_log_early_rec_dev = 1:OM@maxage,
           map_log_early_rec_dev = rep(1, OM@maxage), # Estimate R-init
           prior = prior)

matplot(1956:2021, t(out@SSB/sapply(out@Misc, function(x) x$E[1956:2021 == 2000])), typ = 'l', col = 1, lty = 1,
        xlab = "Year", ylab = "SSB/LRP", ylim = c(0, 7))
abline(h = 1, lty = 3)
abline(h = 0, col = "grey")
matplot(1956:2021, t(out@SSB)/1e3, typ = 'l', col = 1, lty = 1,
        xlab = "Year", ylab = "SSB (mt)", ylim = c(0, 1e2))
abline(h = 0, col = "grey")
LPR_ind <- 1956:2021 == 2020

saveRDS(out, file = "tests/pcod/pcod-RCM2.rds")

out <- readRDS("tests/pcod/pcod-RCM2.rds")
summary(out@mean_fit$SD) # Parameter estimates
plot(1956:2021, out@mean_fit$report$E, typ = 'o', xlab = "Year", ylab = "SSB")
abline(h = out@mean_fit$report$E0_SR) # SSB0
out@mean_fit$report$CR

plot(out, compare = FALSE, s_name = survey_name)

rr <- profile(out, D = seq(0.1, 0.6, 0.025))
rr <- profile(out, R0 = seq(18000, 22000, 1000))
rr <- profile(out, h = seq(0.7, 0.9, 0.05))
rr <- profile(out, R0 = seq(18000, 22000, 1000), h = seq(0.7, 0.9, 0.05))
plot(rr)


