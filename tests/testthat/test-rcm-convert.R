
# New-style stock/fleet/data fixtures matched to an equivalent legacy OM/RCMdata

test_that(".stock2OM() produces an OM that passes check_OM_for_sampling()", {
  new_fixture <- make_new_fixture()
  convOM <- .stock2OM(new_fixture$Stock, new_fixture$Fleet, seed = 1)

  expect_s4_class(convOM, "OM")
  expect_equal(convOM@nsim, nSim)
  expect_equal(convOM@maxage, maxage)
  expect_equal(dim(convOM@cpars$Len_age), c(nSim, maxage + 1, nyears + pYear))

  dummyRCMdata <- methods::new("RCMdata")
  dummyRCMdata@Chist <- Chist_true
  dummyRCMdata@Index <- Ind_true
  dummyRCMdata@I_sd <- rep(0.2, nyears)

  expect_error(SAMtool:::check_OM_for_sampling(convOM, dummyRCMdata), NA)
})

test_that(".data2RCMdata() produces an RCMdata object that passes check_RCMdata()", {
  new_fixture <- make_new_fixture()
  convOM <- .stock2OM(new_fixture$Stock, new_fixture$Fleet, seed = 1)
  convRCMdata <- .data2RCMdata(new_fixture$data, condition = "catch")

  expect_s4_class(convRCMdata, "RCMdata")
  expect_equal(as.numeric(convRCMdata@Chist), Chist_true)
  expect_equal(as.numeric(convRCMdata@Index), Ind_true)

  expect_error(SAMtool:::check_RCMdata(convRCMdata, convOM, condition = "catch", silent = TRUE), NA)
})

test_that("new stock/fleet/data pathway reproduces legacy OM/RCMdata RCM_int() fit", {
  legacy_fixture <- make_legacy_fixture()
  new_fixture <- make_new_fixture()

  out_legacy <- SAMtool:::RCM_int(OM = legacy_fixture$OM, RCMdata = legacy_fixture$RCMdata,
                                  condition = "catch", selectivity = "logistic", s_selectivity = "B",
                                  silent = TRUE)

  convOM <- .stock2OM(new_fixture$Stock, new_fixture$Fleet, seed = 1)
  convRCMdata <- .data2RCMdata(new_fixture$data, condition = "catch")

  out_new <- SAMtool:::RCM_int(OM = convOM, RCMdata = convRCMdata,
                               condition = "catch", selectivity = "logistic", s_selectivity = "B",
                               silent = TRUE)

  expect_true(all(out_legacy@conv))
  expect_true(all(out_new@conv))

  R0_legacy <- vapply(out_legacy@report, function(x) x$R0, numeric(1))
  R0_new <- vapply(out_new@report, function(x) x$R0, numeric(1))
  expect_equal(R0_new, R0_legacy, tolerance = 5e-3)

  expect_equal(out_new@SSB, out_legacy@SSB, tolerance = 5e-3)

  F_legacy <- vapply(out_legacy@report, function(x) x$F, matrix(0, nyears, 1))
  F_new <- vapply(out_new@report, function(x) x$F, matrix(0, nyears, 1))
  expect_equal(F_new, F_legacy, tolerance = 5e-3)

  vul_legacy <- vapply(out_legacy@report, function(x) x$vul[nyears + 1, , 1], numeric(maxage + 1))
  vul_new <- vapply(out_new@report, function(x) x$vul[nyears + 1, , 1], numeric(maxage + 1))
  expect_equal(vul_new, vul_legacy, tolerance = 5e-3)
})

test_that(".fit2om() writes the fit back into a valid, non-default om/rcmodel", {
  new_fixture <- make_new_fixture()

  convOM <- .stock2OM(new_fixture$Stock, new_fixture$Fleet, seed = 1)
  convRCMdata <- .data2RCMdata(new_fixture$data, condition = "catch")

  fit <- SAMtool:::RCM_int(OM = convOM, RCMdata = convRCMdata,
                           condition = "catch", selectivity = "logistic", s_selectivity = "B",
                           silent = TRUE)
  expect_true(all(fit@conv))

  om_skeleton <- MSEtool::OM(Name = "TestOM", Stock = list(TestStock = new_fixture$Stock),
                             Fleet = list(TestStock = list(Fleet1 = new_fixture$Fleet)),
                             nSim = nSim, nYear = nyears, pYear = pYear, CurrentYear = CurrentYear)

  fits <- list(TestStock = fit)
  units <- list(TestStock = list(stocks = "TestStock", complex = 1, herm = FALSE))
  om_fitted <- .fit2om(om_skeleton, fits, units,
                          Data = list(TestStock = new_fixture$data), silent = TRUE)

  expect_error(methods::validObject(om_fitted), NA)

  rcmodel_obj <- methods::new("rcmodel", OM = om_fitted, fits = fits,
                              RCMdata = list(TestStock = convRCMdata), Misc = list())
  expect_error(methods::validObject(rcmodel_obj), NA)

  R0_report <- vapply(fit@report, getElement, numeric(1), "R0")
  expect_equal(as.numeric(om_fitted@Stock$TestStock@SRR@R0), rep(R0_report, length.out = nSim))
  expect_true(all(om_fitted@Stock$TestStock@SRR@R0 != 1000)) # not the pre-fit starting value

  expect_length(om_fitted@Stock$TestStock@Depletion@Final, 0)

  fitted_fleet <- om_fitted@Fleet$TestStock$Fleet1
  F_report <- t(vapply(rep(fit@report, length.out = nSim), function(x) x$F[1:nyears, 1], numeric(nyears)))
  expect_equal(matrix(fitted_fleet@Effort@Effort, nrow(F_report), ncol(F_report)), F_report)

  sel_arr <- fitted_fleet@Selectivity@MeanAtAge
  sel_terminal <- sel_arr[1, , dim(sel_arr)[3], 1]
  expect_equal(unname(sel_terminal), unname(fit@report[[1]]$vul[nyears, , 1]))
  expect_equal(unname(sel_terminal[length(sel_terminal)]), 1) # fully selected plus-group age
})

test_that("RCM2(stock, data = data) reproduces the RCM_int() fit", {
  new_fixture <- make_new_fixture()

  rc <- RCM2(new_fixture$Stock, data = new_fixture$data, Fleet = new_fixture$Fleet,
           nSim = nSim, pYear = pYear, CurrentYear = CurrentYear,
           condition = "catch", selectivity = "logistic", s_selectivity = "B", silent = TRUE)

  expect_s4_class(rc, "rcmodel")
  expect_error(methods::validObject(rc), NA)
  expect_error(methods::validObject(rc@OM), NA)
  expect_equal(names(rc@fits), "TestStock")
  expect_equal(names(rc@RCMdata), "TestStock")

  convOM <- .stock2OM(new_fixture$Stock, new_fixture$Fleet, seed = 1)
  convRCMdata <- .data2RCMdata(new_fixture$data, condition = "catch")
  ref_fit <- SAMtool:::RCM_int(OM = convOM, RCMdata = convRCMdata,
                               condition = "catch", selectivity = "logistic", s_selectivity = "B",
                               silent = TRUE)
  R0_ref <- vapply(ref_fit@report, getElement, numeric(1), "R0")

  expect_equal(as.numeric(rc@OM@Stock$TestStock@SRR@R0), rep(R0_ref, length.out = nSim))
})

test_that("RCM() keeps only the legacy methods", {
  expect_false(methods::existsMethod("RCM", methods::signature("stock", "data")))
  expect_false(methods::existsMethod("RCM", methods::signature("om", "missing")))
  expect_true(methods::existsMethod("RCM", methods::signature("OM", "RCMdata")))
  expect_true(methods::existsMethod("RCM", methods::signature("OM", "list")))
  expect_true(methods::existsMethod("RCM", methods::signature("OM", "Data")))
  expect_true(methods::existsMethod("RCM", methods::signature("list", "RCMdata")))
})

test_that("RCM2(list, data = list) fits two stocks", {
  new_fixture <- make_new_fixture()
  stock2 <- new_fixture$Stock
  stock2@Name <- "TestStock2"

  rc2 <- RCM2(list(TestStock = new_fixture$Stock, TestStock2 = stock2),
            data = list(TestStock = new_fixture$data, TestStock2 = new_fixture$data),
            Fleet = list(TestStock = new_fixture$Fleet, TestStock2 = new_fixture$Fleet),
            nSim = nSim, pYear = pYear, CurrentYear = CurrentYear,
            condition = "catch", selectivity = "logistic", s_selectivity = "B", silent = TRUE)

  expect_s4_class(rc2, "rcmodel")
  expect_error(methods::validObject(rc2), NA)
  expect_error(methods::validObject(rc2@OM), NA)
  expect_equal(sort(names(rc2@fits)), c("TestStock", "TestStock2"))
  expect_true(all(vapply(rc2@fits, function(x) all(x@conv), logical(1))))
})

test_that("RCM2() fits two stocks with two fleets each independently", {
  mk_stock <- function(name, M, Linf, h) {
    ages <- MSEtool::Ages(MaxAge = maxage, MinAge = 0, Units = "year", PlusGroup = TRUE)
    length_obj <- MSEtool::Length(Pars = list(Linf = Linf, K = 0.2, t0 = -0.5), Model = "vonBert")
    weight_obj <- MSEtool::Weight(Pars = list(Alpha = 1e-5, Beta = 3), Model = "WeightatLength")
    natm <- MSEtool::NaturalMortality(Pars = list(M = M), Model = "MortalityAtAge")
    mat <- MSEtool::Maturity(Pars = list(L50 = 40, L50_95 = 5), Model = "MaturityAtLength")
    srr <- MSEtool::SRR(Pars = list(h = h), Model = "BevertonHolt", R0 = 1000, SD = 0.4, AC = 0)
    stock <- MSEtool::Stock(Name = name, Ages = ages, Length = length_obj, Weight = weight_obj,
                             NaturalMortality = natm, Maturity = mat, SRR = srr)
    MSEtool::PopulateStock(stock, nYear = nyears, pYear = pYear, nSim = nSim,
                            CurrentYear = CurrentYear, seed = 1)
  }
  mk_fleet <- function(name, L5, LFS, Vmaxlen = 1) {
    sel <- MSEtool::Selectivity(Pars = list(L5 = L5, LFS = LFS, Vmaxlen = Vmaxlen), Model = "DoubleNormal")
    MSEtool::Fleet(Name = name, Selectivity = sel)
  }
  mk_data <- function() {
    Chist1 <- round(100 + 20 * sin(seq(0, 3 * pi, length.out = nyears)) + seq(0, 30, length.out = nyears), 1)
    Chist2 <- round(60 + 10 * cos(seq(0, 2 * pi, length.out = nyears)) + seq(0, 15, length.out = nyears), 1)
    Ind_true2 <- round(pmax(50 - 1.2 * seq_len(nyears) + 5 * cos(seq(0, 2 * pi, length.out = nyears)), 1), 2)
    dat <- MSEtool::Data(Years = (CurrentYear - nyears + 1):CurrentYear)
    dat@Landings <- MSEtool::CatchData(Name = c("Small", "Large"), Value = unname(cbind(Chist1, Chist2)),
                                        CV = matrix(0.01, nyears, 2))
    dat@CPUE <- MSEtool::IndicesData(Value = matrix(Ind_true2, ncol = 1), CV = matrix(0.2, nyears, 1))
    dat
  }

  # Different biology per stock, different selectivity per fleet within a stock.
  stockA <- mk_stock("StockA", M = 0.2, Linf = 100, h = 0.7)
  stockB <- mk_stock("StockB", M = 0.35, Linf = 80, h = 0.85)
  fleetA <- list(Small = mk_fleet("Small", L5 = 20, LFS = 35), Large = mk_fleet("Large", L5 = 50, LFS = 70))
  fleetB <- list(Small = mk_fleet("Small", L5 = 15, LFS = 30), Large = mk_fleet("Large", L5 = 40, LFS = 60))
  dataA <- mk_data()
  dataB <- mk_data()

  rc <- RCM2(list(StockA = stockA, StockB = stockB),
            data = list(StockA = dataA, StockB = dataB),
            Fleet = list(StockA = fleetA, StockB = fleetB),
            nSim = nSim, pYear = pYear, CurrentYear = CurrentYear,
            condition = "catch", selectivity = "logistic", s_selectivity = "B", silent = TRUE)

  expect_s4_class(rc, "rcmodel")
  expect_error(methods::validObject(rc), NA)
  expect_error(methods::validObject(rc@OM), NA)
  expect_equal(sort(names(rc@fits)), c("StockA", "StockB"))
  expect_true(all(rc@fits$StockA@conv))
  expect_true(all(rc@fits$StockB@conv))

  # No cross-contamination between stocks: fitted R0 differs substantially.
  R0_A <- mean(rc@OM@Stock$StockA@SRR@R0)
  R0_B <- mean(rc@OM@Stock$StockB@SRR@R0)
  expect_false(isTRUE(all.equal(R0_A, R0_B)))
  expect_gt(abs(R0_A - R0_B) / R0_A, 0.1)

  # Independent per-fleet selectivity/F within a stock (no collapse to one aggregate fleet).
  selA_small <- rc@OM@Fleet$StockA$Small@Selectivity@MeanAtAge
  selA_large <- rc@OM@Fleet$StockA$Large@Selectivity@MeanAtAge
  nAge <- dim(selA_small)[2]
  expect_false(isTRUE(all.equal(
    selA_small[1, , dim(selA_small)[3], 1],
    selA_large[1, , dim(selA_large)[3], 1]
  )))

  FA_small <- rc@OM@Fleet$StockA$Small@Effort@Effort
  FA_large <- rc@OM@Fleet$StockA$Large@Effort@Effort
  expect_false(isTRUE(all.equal(as.numeric(FA_small[1, ]), as.numeric(FA_large[1, ]))))
})
