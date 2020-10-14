
library(testthat)
library(MSEtool)

context("SCA in assessment mode")

load('tests_results/Data_from_Hist.RData')

tdir <- tempdir()

test_that("SCA nonconvergence", {

  res <- expect_s4_class(SCA(Data = Data_from_Hist[[1]]), "Assessment")
  expect_warning(plot(res, save_dir = tdir))

})


test_that("SCA2 nonconvergence", {

  res <- expect_s4_class(SCA2(x = 12, Data = Data_from_Hist[[1]]), "Assessment")
  expect_warning(plot(res, save_dir = tdir))

})


test_that("SCA assess model", {
  res <- expect_s4_class(SCA(Data = Simulation_1), "Assessment")
  expect_s4_class(SCA(Data = Simulation_1, start = list(R0 = 1)), "Assessment")
  expect_s4_class(SCA(Data = Simulation_1, integrate = TRUE, fix_tau = F), "Assessment")
  expect_s4_class(SCA(Data = Simulation_1, SR = "Ricker"), "Assessment")
  expect_s4_class(SCA(Data = Red_snapper, SR = "Ricker"), "Assessment")
  expect_s4_class(SCA(Data = Red_snapper, start = list(R0 = 0.25)), "Assessment")
  expect_s4_class(SCA(Data = Simulation_1, fix_h = TRUE), "Assessment")

  expect_s4_class(SCA(Data = Simulation_1, early_dev = "comp", integrate = TRUE), "Assessment")
  expect_s4_class(SCA(Data = Simulation_1, early_dev = "comp"), "Assessment")

  expect_s4_class(SCA(Data = Simulation_1, late_dev = 5, integrate = TRUE), "Assessment")
  expect_s4_class(SCA(Data = Simulation_1, late_dev = 5), "Assessment")

  expect_s4_class(SCA(Data = Simulation_1, CAA_multiplier = 0.4), "Assessment")

  expect_equivalent(plot(res, save_dir = tdir), invisible())

  pro <- profile(res, R0 = seq(600, 700, 10), h = seq(0.7, 0.99, 0.025), save_dir = tdir)
  expect_type(pro, "list")
  expect_true(is.data.frame(pro))

  expect_type(retrospective(res, 5, save_dir = tdir), "list")

})


test_that("SCA assess model fix_h", {
  res <- expect_s4_class(SCA(Data = Simulation_1, fix_h = TRUE), "Assessment")

  expect_equivalent(plot(res, save_dir = tdir), invisible())

  pro <- profile(res, R0 = seq(600, 700, 10), figure = FALSE)
  expect_type(pro, "list")
  expect_true(is.data.frame(pro))

  expect_type(retrospective(res, 5, save_dir = tdir), "list")

})


test_that("SCA2 assess model", {
  expect_s4_class(SCA2(Data = Simulation_1), "Assessment")
  expect_s4_class(SCA2(Data = Simulation_1, start = list(R0 = 1)), "Assessment")
  expect_s4_class(SCA2(Data = Simulation_1, integrate = TRUE, fix_tau = F), "Assessment")
  expect_s4_class(SCA2(Data = Simulation_1, SR = "Ricker"), "Assessment")

  expect_s4_class(SCA2(Data = Red_snapper, SR = "Ricker"), "Assessment")
  expect_s4_class(SCA2(Data = Simulation_1, fix_h = TRUE), "Assessment")
  expect_s4_class(SCA2(Data = Simulation_1, early_dev = "comp"), "Assessment")
  expect_s4_class(SCA2(Data = Simulation_1, common_dev = 5), "Assessment")
  expect_s4_class(SCA2(Data = Simulation_1, CAA_multiplier = 0.4), "Assessment")

  res <- expect_s4_class(SCA2(Data = Red_snapper), "Assessment")
  expect_s4_class(SCA2(Data = Red_snapper, SR = "Ricker"), "Assessment")

  expect_equivalent(plot(res, save_dir = tdir), invisible())

  pro <- profile(res, meanR = seq(1.5, 2.5, 0.1), figure = FALSE)
  expect_type(pro, "list")
  expect_true(is.data.frame(pro))

  expect_type(retrospective(res, 5, save_dir = tdir), "list")

})

test_that("SCA2 assess model fix_h", {

  res <- expect_s4_class(SCA2(Data = Simulation_1, fix_h = TRUE), "Assessment")

  expect_equivalent(plot(res, save_dir = tdir), invisible())

  pro <- profile(res, meanR = seq(1.5, 2.5, 0.1), figure = FALSE)
  expect_type(pro, "list")
  expect_true(is.data.frame(pro))

  expect_type(retrospective(res, 5, save_dir = tdir), "list")

})

