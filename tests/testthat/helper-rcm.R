
# Shared single-stock fixtures for the RCM tests

n_sim <- 3
nyears <- 20
p_year <- 2
maxage <- 8
current_year <- 2020

Chist_true <- round(100 + 20 * sin(seq(0, 3 * pi, length.out = nyears)) + seq(0, 30, length.out = nyears), 1)
Ind_true <- round(pmax(50 - 1.2 * seq_len(nyears) + 5 * cos(seq(0, 2 * pi, length.out = nyears)), 1), 2)

make_legacy_fixture <- function() {
  legacyStock <- methods::new("Stock")
  legacyStock@maxage <- maxage
  legacyStock@R0 <- 1000
  legacyStock@M <- c(0.2, 0.2)
  legacyStock@h <- c(0.7, 0.7)
  legacyStock@SRrel <- 1L
  legacyStock@Perr <- c(0.4, 0.4)
  legacyStock@AC <- c(0, 0)
  legacyStock@Linf <- c(100, 100)
  legacyStock@K <- c(0.2, 0.2)
  legacyStock@t0 <- c(-0.5, -0.5)
  legacyStock@LenCV <- c(0.1, 0.1)
  legacyStock@a <- 1e-5
  legacyStock@b <- 3
  legacyStock@Msd <- legacyStock@Linfsd <- legacyStock@Ksd <- c(0, 0)
  legacyStock@Size_area_1 <- legacyStock@Frac_area_1 <- legacyStock@Prob_staying <- c(0.5, 0.5)
  legacyStock@Fdisc <- c(0, 0)
  legacyStock@D <- c(0, 0)
  legacyStock@L50 <- c(40, 40)
  legacyStock@L50_95 <- c(5, 5)

  legacyFleet <- methods::new("Fleet")
  legacyFleet@CurrentYr <- current_year
  legacyFleet@nyears <- nyears
  legacyFleet@EffYears <- c(1, nyears)
  legacyFleet@EffLower <- legacyFleet@EffUpper <- c(0, 1)
  legacyFleet@Esd <- legacyFleet@qinc <- legacyFleet@qcv <- c(0, 0)
  legacyFleet@L5 <- c(30, 30)
  legacyFleet@LFS <- c(50, 50)
  legacyFleet@Vmaxlen <- c(1, 1)
  legacyFleet@isRel <- FALSE
  legacyFleet@DR <- c(0, 0)
  legacyFleet@Spat_targ <- c(1, 1)
  legacyFleet@MPA <- FALSE

  legacyOM <- methods::new("OM", Stock = legacyStock, Fleet = legacyFleet,
                            Obs = MSEtool::Generic_Obs, Imp = MSEtool::Perfect_Imp,
                            nsim = n_sim, proyears = p_year)
  legacyOM@seed <- 1
  legacyOM@maxF <- 3

  legacyRCMdata <- methods::new("RCMdata")
  legacyRCMdata@Chist <- Chist_true
  legacyRCMdata@C_sd <- rep(0.01, nyears)
  legacyRCMdata@Index <- Ind_true
  legacyRCMdata@I_sd <- rep(0.2, nyears)

  list(OM = legacyOM, RCMdata = legacyRCMdata)
}

make_new_fixture <- function() {
  ages <- MSEtool::Ages(MaxAge = maxage, MinAge = 0, Units = "year", PlusGroup = TRUE)
  length_obj <- MSEtool::Length(Pars = list(Linf = 100, K = 0.2, t0 = -0.5), Model = "vonBert")
  weight_obj <- MSEtool::Weight(Pars = list(Alpha = 1e-5, Beta = 3), Model = "WeightatLength")
  natm <- MSEtool::NaturalMortality(Pars = list(M = 0.2), Model = "MortalityAtAge")
  mat <- MSEtool::Maturity(Pars = list(L50 = 40, L50_95 = 5), Model = "MaturityAtLength")
  srr <- MSEtool::SRR(Pars = list(h = 0.7), Model = "BevertonHolt", R0 = 1000, SD = 0.4, AC = 0)

  stock <- MSEtool::Stock(Name = "TestStock", Ages = ages, Length = length_obj, Weight = weight_obj,
                           NaturalMortality = natm, Maturity = mat, SRR = srr)
  stock <- MSEtool::PopulateStock(stock, nYear = nyears, pYear = p_year, nSim = n_sim,
                                   CurrentYear = current_year, seed = 1)

  sel <- MSEtool::Selectivity(Pars = list(L5 = 30, LFS = 50, Vmaxlen = 1), Model = "DoubleNormal")
  fleet <- MSEtool::Fleet(Name = "Fleet1", Selectivity = sel)

  dat <- MSEtool::Data(Years = (current_year - nyears + 1):current_year)
  dat@Landings <- MSEtool::CatchData(Value = matrix(Chist_true, ncol = 1), CV = matrix(0.01, nyears, 1))
  dat@CPUE <- MSEtool::IndicesData(Value = matrix(Ind_true, ncol = 1), CV = matrix(0.2, nyears, 1))

  list(Stock = stock, Fleet = fleet, data = dat)
}
