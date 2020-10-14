
library(testthat)
library(MSEtool)

context("Delay-Difference in assessment mode")
Red_snapper@sigmaR <- 0.3

tdir <- tempdir()

if(exists("res")) rm(res)
if(exists("pro")) rm(pro)

load('tests_results/Data_from_Hist.RData')

test_that("DD_TMB nonconvergence", {

  res <- expect_s4_class(DD_TMB(x = 21, Data = Data_from_Hist[[1]]), "Assessment")
  expect_warning(plot(res, save_dir = tdir))

})

if(exists("res")) rm(res)
if(exists("pro")) rm(pro)

test_that("DD_TMB nonconvergence", {

  res <- expect_s4_class(DD_SS(x = 21, Data = Data_from_Hist[[1]]), "Assessment")
  expect_warning(plot(res, save_dir = tdir))

})

if(exists("res")) rm(res)
if(exists("pro")) rm(pro)

test_that("DD_TMB assess model", {
  expect_s4_class(DD_TMB(Data = SimulatedData, start = list(R0 = 1)), "Assessment")
  expect_s4_class(DD_TMB(Data = SimulatedData), "Assessment")
  expect_s4_class(DD_TMB(Data = SimulatedData, SR = "Ricker"), "Assessment")

  res <- expect_s4_class(DD_TMB(Data = Red_snapper), "Assessment")

  expect_s4_class(DD_TMB(Data = Red_snapper, start = list(R0 = 0.25)), "Assessment")
  expect_s4_class(DD_TMB(Data = Red_snapper, SR = "Ricker"), "Assessment")

  expect_equivalent(plot(res, save_dir = tdir), invisible())

  pro <- profile(res, R0 = seq(0.75, 1.25, 0.025), h = seq(0.95, 1, 0.0025), save_dir = tdir)
  expect_type(pro, "list")
  expect_true(is.data.frame(pro))

  expect_type(retrospective(res, 5, save_dir = tdir), "list")

})


if(exists("res")) rm(res)
if(exists("pro")) rm(pro)

test_that("DD_SS assess model", {
  res <- expect_s4_class(DD_SS(Data = Red_snapper), "Assessment")

  expect_s4_class(DD_SS(Data = Red_snapper, start = list(R0 = 0.25)), "Assessment")

  expect_s4_class(DD_SS(Data = SimulatedData, integrate = TRUE), "Assessment")
  expect_s4_class(DD_SS(Data = SimulatedData, SR = "Ricker"), "Assessment")
  expect_s4_class(DD_SS(Data = SimulatedData, start = list(R0 = 1)), "Assessment")

  expect_s4_class(DD_SS(Data = Red_snapper), "Assessment")
  expect_s4_class(DD_SS(Data = Red_snapper, integrate = TRUE), "Assessment")
  expect_s4_class(DD_SS(Data = Red_snapper, SR = "Ricker"), "Assessment")

  expect_equivalent(plot(res, save_dir = tdir), invisible())

  pro <- profile(res, R0 = seq(0.75, 1.25, 0.025), h = seq(0.9, 0.99, 0.01), save_dir = tdir)
  expect_type(pro, "list")
  expect_true(is.data.frame(pro))

  expect_type(retrospective(res, 5, save_dir = tdir), "list")

})


