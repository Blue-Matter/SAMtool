
fit_fixture <- function(fx = make_new_fixture(), ...) {
  RCM2(fx$Stock, data = fx$data, Fleet = fx$Fleet, nSim = n_sim, pYear = p_year,
      CurrentYear = current_year, s_selectivity = "B", silent = TRUE, ...)
}

test_that("the conditioned om simulates the fitted SSB, recruitment and catch", {
  rc <- fit_fixture()
  hist <- MSEtool::Simulate(rc@OM, silent = TRUE)
  fit <- rc@fits[[1]]

  sb <- MSEtool::SBiomass(hist)
  sb <- matrix(sb$Value[sb$Period == "Historical"], n_sim, nyears, byrow = TRUE)
  expect_equal(sb, unname(fit@SSB[, 1:nyears]), tolerance = 1e-4)

  rec <- MSEtool::Recruits(hist)
  rec <- matrix(rec$Value[rec$Period == "Historical"], n_sim, nyears, byrow = TRUE)
  rec_fit <- t(vapply(rep(fit@report, length.out = n_sim), function(x) x$R[1:nyears], numeric(nyears)))
  expect_equal(rec, rec_fit, tolerance = 1e-4)

  removals <- MSEtool::Removals(hist)
  removals <- removals$Value[removals$Sim == 1 & removals$Period == "Historical"]
  expect_equal(removals, Chist_true, tolerance = 1e-3)
})

test_that("fitted recruitment deviations include the bias correction", {
  rc <- fit_fixture()
  r <- rc@fits[[1]]@report[[1]]
  dev <- as.numeric(rc@OM@Stock[[1]]@SRR@RecDevHist[1, ])
  expect_equal(dev[-1], exp(r$log_rec_dev[-1] - 0.5 * r$tau^2))
})

test_that("unnamed indices get default names and the RCM index selectivity", {
  rc <- fit_fixture()
  expect_equal(rc@OM@Data[[1]]@CPUE@Name, "CPUE1")
  expect_equal(rc@OM@Data[[1]]@Landings@Name, "Fleet1")
  expect_equal(MSEtool::IndexSelectivity(rc@OM, "CPUE1", Type = "CPUE"), "Biomass")
})

test_that("RCM2(om) conditions an om with data in OM@Data and keeps its settings", {
  fx <- make_new_fixture()
  om <- MSEtool::OM(Name = "TestOM", nSim = n_sim, nYear = nyears, pYear = p_year, CurrentYear = current_year,
                    Stock = list(TestStock = fx$Stock), Fleet = list(TestStock = list(Fleet1 = fx$Fleet)),
                    Data = list(TestStock = fx$data), maxF = 2.5)
  rc <- RCM2(om, s_selectivity = "B", silent = TRUE)

  expect_s4_class(rc, "rcmodel")
  expect_equal(rc@OM@Name, "TestOM")
  expect_equal(rc@OM@maxF, 2.5)
  R0 <- vapply(rc@fits[[1]]@report, getElement, numeric(1), "R0")
  expect_equal(as.numeric(rc@OM@Stock[[1]]@SRR@R0), rep(R0, length.out = n_sim))
  expect_s4_class(MSEtool::Simulate(rc@OM, silent = TRUE), "hist")
})

test_that("RCM2() populates an unpopulated stock from the data years", {
  fx <- make_new_fixture()
  stk <- MSEtool::Stock(Name = "TestStock", Ages = fx$Stock@Ages, Length = MSEtool::Length(Pars = list(Linf = 100, K = 0.2, t0 = -0.5), Model = "vonBert"),
                        Weight = MSEtool::Weight(Pars = list(Alpha = 1e-5, Beta = 3), Model = "WeightatLength"),
                        NaturalMortality = MSEtool::NaturalMortality(Pars = list(M = 0.2), Model = "MortalityAtAge"),
                        Maturity = MSEtool::Maturity(Pars = list(L50 = 40, L50_95 = 5), Model = "MaturityAtLength"),
                        SRR = MSEtool::SRR(Pars = list(h = 0.7), Model = "BevertonHolt", R0 = 1000, SD = 0.4, AC = 0))
  rc <- RCM2(stk, data = fx$data, Fleet = fx$Fleet, nSim = n_sim, pYear = p_year,
            s_selectivity = "B", silent = TRUE)
  expect_equal(rc@OM@nYear, nyears)
  expect_equal(rc@OM@CurrentYear, current_year)
  expect_true(all(rc@fits[[1]]@conv))
})

test_that("RCM2() fits catch without CVs and requires index CVs", {
  fx <- make_new_fixture()
  fx$data@Landings@CV <- NULL
  rc <- fit_fixture(fx)
  expect_true(all(rc@fits[[1]]@conv))

  fx$data@CPUE@CV <- NULL
  expect_error(fit_fixture(fx), "Index CVs")
})

test_that("RCM2() rejects seasonal models and mismatched data years", {
  fx <- make_new_fixture()
  om <- .BuildRCMOM(fx$Stock, fx$data, Fleet = fx$Fleet, nSim = n_sim, pYear = p_year)
  om@Seasons <- 4
  expect_error(RCM2(om, silent = TRUE), "annual")

  fx$data@Landings@Value <- fx$data@Landings@Value[-1, , drop = FALSE]
  fx$data@Landings@CV <- fx$data@Landings@CV[-1, , drop = FALSE]
  expect_error(fit_fixture(fx), "historical year")
})

test_that("rcmodel show and plot methods are registered", {
  expect_true(methods::existsMethod("show", "rcmodel"))
  expect_true(methods::existsMethod("plot", methods::signature("rcmodel", "missing")))
  rc <- fit_fixture()
  expect_output(methods::show(rc), "TestStock: 3/3 converged")
})

test_that("simulations with different natural mortality are fitted and simulated separately", {
  fx <- make_new_fixture()
  stk <- fx$Stock
  stk@NaturalMortality <- MSEtool::NaturalMortality(Pars = list(M = c(0.18, 0.2, 0.22)), Model = "MortalityAtAge")
  stk <- MSEtool::PopulateStock(stk, nYear = nyears, pYear = p_year, nSim = n_sim, CurrentYear = current_year,
                                seed = 1, force = TRUE)
  rc <- RCM2(stk, data = fx$data, Fleet = fx$Fleet, s_selectivity = "B", silent = TRUE)
  fit <- rc@fits[[1]]
  expect_length(fit@report, n_sim)
  R0 <- vapply(fit@report, getElement, numeric(1), "R0")
  expect_gt(stats::sd(R0), 0)

  hist <- MSEtool::Simulate(rc@OM, silent = TRUE)
  sb <- MSEtool::SBiomass(hist)
  sb <- matrix(sb$Value[sb$Period == "Historical"], n_sim, nyears, byrow = TRUE)
  expect_equal(sb, unname(fit@SSB[, 1:nyears]), tolerance = 1e-4)
})

test_that("the legacy RCM(OM, RCMdata) method still returns an RCModel", {
  fx <- make_legacy_fixture()
  out <- RCM(fx$OM, fx$RCMdata, condition = "catch", s_selectivity = "B", silent = TRUE)
  expect_s4_class(out, "RCModel")
  expect_s4_class(out@OM, "OM")
  expect_true(all(out@conv))
})

test_that("fleets in the conditioned om retain all catch", {
  fx <- make_new_fixture()
  fx$Fleet@Retention <- MSEtool::Retention(MeanAtAge = c(0, 0, 0, 0.5, rep(1, maxage - 3)))
  fx$Fleet@DiscardMortality <- MSEtool::DiscardMortality(MeanAtAge = 0.5)
  rc <- fit_fixture(fx)
  hist <- MSEtool::Simulate(rc@OM, silent = TRUE)

  removals <- MSEtool::Removals(hist)
  removals <- removals$Value[removals$Sim == 1 & removals$Period == "Historical"]
  expect_equal(removals, Chist_true, tolerance = 1e-3)
})
