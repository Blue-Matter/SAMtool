


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
  Stock@D <- c(0.5, 0.5)
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
  
  if (stochastic) {
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


#OM <- make_pcod_OM(nsim = 2)
OM <- make_pcod_OM(nsim = 48, stochastic = TRUE)
mat_ogive <- OM@cpars$Mat_age[1,,1]

make_RCMdata <- function(OM) {
  
  pcod_data <- readRDS("data-raw/pcod-data-5ABCD.rds")
  RCMdata <- new("RCMdata")
  
  RCMdata@Chist <- pcod_data$Chist
  RCMdata@C_sd <- matrix(0.05, OM@nyears, 1)
  RCMdata@Index <- pcod_data$Index %>% 
    structure(dimnames = list(Year = 1956:2020, 
                              Survey = c("HSMSAS", "SYN QCS", "SYN HS", "CPUE Hist", "CPUE Modern")))
  RCMdata@I_sd <- pcod_data$I_sd
  RCMdata@MS <- pcod_data$MS
  RCMdata@MS_type <- "weight"
  RCMdata@MS_cv <- 0.2
  RCMdata@I_units <- pcod_data$I_units
  
  return(RCMdata)
}
RCM_data <- make_RCMdata(OM)

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
  structure(dimnames = list(colnames(RCM_data@Index), c("Mean", "SD")))
prior <- list(#M = c(0.5, 0.1), h = c(0.7, 0.15), 
              q = prior_q)

pcod <- list(OM = OM, data = RCM_data, prior = prior)
usethis::use_data(pcod, overwrite = TRUE)
