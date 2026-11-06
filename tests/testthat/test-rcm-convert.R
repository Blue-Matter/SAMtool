
# New-style stock/fleet/data fixtures matched to an equivalent legacy OM/RCMdata

test_that(".stock2OM() produces an OM that passes check_OM_for_sampling()", {
  new_fixture <- make_new_fixture()
  convOM <- .stock2OM(new_fixture$Stock, new_fixture$Fleet, seed = 1)

  expect_s4_class(convOM, "OM")
  expect_equal(convOM@nsim, n_sim)
  expect_equal(convOM@maxage, maxage)
  expect_equal(dim(convOM@cpars$Len_age), c(n_sim, maxage + 1, nyears + p_year))

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
                             nSim = n_sim, nYear = nyears, pYear = p_year, CurrentYear = current_year)

  fits <- list(TestStock = fit)
  units <- list(TestStock = list(stocks = "TestStock", rec = "TestStock", share = list(TestStock = 1),
                                complex = 1, combined = FALSE))
  om_fitted <- .fit2om(om_skeleton, fits, units,
                          Data = list(TestStock = new_fixture$data), silent = TRUE)

  expect_error(methods::validObject(om_fitted), NA)

  rcmodel_obj <- methods::new("rcmodel", OM = om_fitted, fits = fits,
                              RCMdata = list(TestStock = convRCMdata), Misc = list())
  expect_error(methods::validObject(rcmodel_obj), NA)

  R0_report <- vapply(fit@report, getElement, numeric(1), "R0")
  expect_equal(as.numeric(om_fitted@Stock$TestStock@SRR@R0), rep(R0_report, length.out = n_sim))
  expect_true(all(om_fitted@Stock$TestStock@SRR@R0 != 1000)) # not the pre-fit starting value

  expect_length(om_fitted@Stock$TestStock@Depletion@Final, 0)

  fitted_fleet <- om_fitted@Fleet$TestStock$Fleet1
  F_report <- t(vapply(rep(fit@report, length.out = n_sim), function(x) x$F[1:nyears, 1], numeric(nyears)))
  expect_equal(matrix(fitted_fleet@Effort@Effort, nrow(F_report), ncol(F_report)), F_report)

  sel_arr <- fitted_fleet@Selectivity@MeanAtAge
  sel_terminal <- sel_arr[1, , dim(sel_arr)[3], 1]
  expect_equal(unname(sel_terminal), unname(fit@report[[1]]$vul[nyears, , 1]))
  expect_equal(unname(sel_terminal[length(sel_terminal)]), 1) # fully selected plus-group age
})

test_that("RCM2(stock, data = data) reproduces the RCM_int() fit", {
  new_fixture <- make_new_fixture()

  rc <- RCM2(new_fixture$Stock, data = new_fixture$data, Fleet = new_fixture$Fleet,
           nSim = n_sim, pYear = p_year, CurrentYear = current_year,
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

  expect_equal(as.numeric(rc@OM@Stock$TestStock@SRR@R0), rep(R0_ref, length.out = n_sim))
})

test_that("RCM() keeps only the legacy methods", {
  expect_false(methods::existsMethod("RCM", methods::signature("stock", "data")))
  expect_false(methods::existsMethod("RCM", methods::signature("om", "missing")))
  expect_true(methods::existsMethod("RCM", methods::signature("OM", "RCMdata")))
  expect_true(methods::existsMethod("RCM", methods::signature("OM", "list")))
  expect_true(methods::existsMethod("RCM", methods::signature("OM", "Data")))
  expect_true(methods::existsMethod("RCM", methods::signature("list", "RCMdata")))
})

test_that("RCM2() fits several stocks as one population to one data set", {
  fx <- make_new_fixture()
  stockA <- stockB <- fx$Stock
  stockA@SRR@R0 <- 600
  stockB@SRR@R0 <- 400

  rc <- RCM2(list(StockA = stockA, StockB = stockB), data = fx$data, Fleet = fx$Fleet,
             nSim = n_sim, pYear = p_year, s_selectivity = "B", silent = TRUE)

  expect_named(rc@fits, "StockA-StockB")
  expect_true(all(rc@fits[[1]]@conv))
  expect_equal(names(rc@OM@Complexes), "StockA-StockB")

  # The fitted R0 is split by the input R0
  R0 <- vapply(rc@fits[[1]]@report, getElement, numeric(1), "R0")
  expect_equal(as.numeric(rc@OM@Stock$StockA@SRR@R0), 0.6 * rep(R0, length.out = n_sim))
  expect_equal(as.numeric(rc@OM@Stock$StockB@SRR@R0), 0.4 * rep(R0, length.out = n_sim))

  # Same biology, so the summed history matches the fit
  hist <- MSEtool::Simulate(rc@OM, silent = TRUE)
  sb <- MSEtool::SBiomass(hist)
  sb <- sb[sb$Period == "Historical", ]
  sb <- tapply(sb$Value, list(sb$Sim, sb$Year), sum)
  # Identical simulations are simulated once
  ssb_fit <- rc@fits[[1]]@SSB[seq_len(nrow(sb)), 1:nyears, drop = FALSE]
  expect_equal(unname(sb), unname(ssb_fit), tolerance = 1e-4)

  r <- MSEtool::Removals(hist)
  r <- r[r$Sim == 1 & r$Period == "Historical", ]
  expect_equal(as.numeric(tapply(r$Value, r$Year, sum)), Chist_true, tolerance = 1e-3)
})

test_that("RCM2() splits recruitment equally without input R0", {
  fx <- make_new_fixture()
  stockA <- stockB <- fx$Stock
  stockA@SRR@R0 <- stockB@SRR@R0 <- numeric()
  rc <- RCM2(list(StockA = stockA, StockB = stockB), data = fx$data, Fleet = fx$Fleet,
             nSim = n_sim, pYear = p_year, s_selectivity = "B", silent = TRUE)
  expect_equal(rc@OM@Stock$StockA@SRR@R0, rc@OM@Stock$StockB@SRR@R0)
})

test_that("RCM2() fits one data set", {
  fx <- make_new_fixture()
  stockB <- fx$Stock
  expect_error(RCM2(list(StockA = fx$Stock, StockB = stockB), data = list(fx$data, fx$data),
                    Fleet = fx$Fleet, nSim = n_sim, pYear = p_year, silent = TRUE), "one data set")
  expect_error(RCM2(fx$Stock, Fleet = fx$Fleet, nSim = n_sim, pYear = p_year, silent = TRUE), "data")

  # A list of one data object is accepted
  rc <- RCM2(fx$Stock, data = list(fx$data), Fleet = fx$Fleet, nSim = n_sim, pYear = p_year,
             s_selectivity = "B", silent = TRUE)
  expect_named(rc@fits, "TestStock")
})
