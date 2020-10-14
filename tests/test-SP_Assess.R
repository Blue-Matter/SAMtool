
library(testthat)
library(MSEtool)


context("Surplus production in assessment mode")

tdir <- tempdir()

test_that("SP nonconvergence", {

  res <- expect_s4_class(SP(Data = swordfish, start = list(MSY = 1e9)), "Assessment")
  expect_equivalent(plot(res, save_dir = tdir), invisible())

})

test_that("SP_SS nonconvergence", {

  res <- expect_s4_class(SP_SS(Data = swordfish, fix_sigma = FALSE), "Assessment")
  expect_equivalent(plot(res, save_dir = tdir), invisible())

})


test_that("SP assess model", {
  expect_s4_class(SP(Data = swordfish), "Assessment")
  res <- expect_s4_class(SP(Data = swordfish, start = list(dep = 0.95)), "Assessment")

  expect_equivalent(plot(res, save_figure = FALSE), invisible())

  pro <- profile(res, UMSY = seq(0.05, 0.5, 0.005),
                            MSY = seq(0.1, 2, 0.01) * 1e4, figure = FALSE)
  expect_type(pro, "list")
  expect_true(is.data.frame(pro))

  expect_type(retrospective(res, 5, save_figure = FALSE), "list")

})


test_that("SP_SS assess model", {
  res <- expect_s4_class(SP_SS(Data = swordfish, start = list(dep = 0.95, tau = 0.3), fix_sigma = TRUE), "Assessment")

  expect_equivalent(plot(res, save_dir = tdir), invisible())

  pro <- profile(res, UMSY = seq(0.15, 0.2, 0.005),
                            MSY = seq(1.4, 1.5, 0.01) * 1e4, save_dir = tdir)
  expect_type(pro, "list")
  expect_true(is.data.frame(pro))

  expect_type(retrospective(res, 5, save_dir = tdir), "list")

})

