
library(testthat)
library(MSEtool)

# Tests for certain features:
#
# Default (one season):                                  GOM Red snapper
# Multiple seasons within year:                          ICCAT bigeye tuna
# Season-as-years (i.e. seasonal time step as year):     Indian Ocean yellowfin tuna
# Ricker S-R function:
# Discard mortality:                                     Something GOM, maybe greater amberjack, Spanish mackerel
# No rec devs:                                           GOM vermilion snapper
# Multiple areas:                                        ICCAT bigeye tuna, IO yellowfin tuna
# Age-varying growth (deviations in K)                   ICCAT yellowfin tuna
# Time-varying growth
# Maturity specified in length vs. age:
# Then, do both for both SS2OM and SS2Data

context("Test SS2OM")

get_nsim <- function(x) {
  if(is.vector(x)) nsim <- length(x)
  if(is.matrix(x) || is.array(x)) nsim <- dim(x)[1]
  return(nsim)
}


DLM_dev_home <- "C:/~/DLMDev"

test_that("GOM red snapper", {
  dir <- file.path(DLM_dev_home, "DLMextra_Data/buildOMfromSS/Red_snapper_GOM")
  res <- expect_s4_class(SS2OM(dir), "OM")

  out <- plot(res)
  expect_type(out, "list")

  out <- runMSE(res, Hist = TRUE)
  expect_type(out, "list")

  # cpars objects are all of the correct dimension
  expect_true(all(vapply(res@cpars, get_nsim, numeric(1)) == res@nsim))

  expect_true(all(dim(res@cpars$Len_age) == c(res@nsim, res@maxage, res@proyears + res@nyears)))
  expect_true(all(res@cpars$Len_age > 0))

  expect_true(all(dim(res@cpars$Wt_age) == c(res@nsim, res@maxage, res@proyears + res@nyears)))
  expect_true(all(res@cpars$Wt_age > 0))

  expect_true(length(res@cpars$hs) == res@nsim)
  hmax <- ifelse(res@SRrel == 1L, 1, Inf)
  expect_true(all(res@cpars$hs >= 0.2 & res@cpars$h <= hmax))

  expect_true(all(dim(res@cpars$Perr_y) == c(res@nsim, res@maxage - 1 + res@proyears + res@nyears)))
  expect_true(all(res@cpars$Perr_y > 0))

  expect_true(all(dim(res@cpars$V) == c(res@nsim, res@maxage, res@proyears + res@nyears)))
  expect_true(all(res@cpars$V >= 0 & res@cpars$V <= 1))

  expect_true(all(dim(res@cpars$Find) == c(res@nsim, res@nyears)))
  expect_true(all(res@cpars$Find >= 0))

  expect_equal(DLMtool::cparscheck(res@cpars), res@nsim)

  expect_s4_class(runMSE(res, Hist = TRUE, silent = TRUE), "Hist")
})


test_that("Indian Ocean yellowfin tuna", {
  dir <- file.path(DLM_dev_home, "DLMextra_Data/buildOMfromSS/Yellowfin_tuna_IO")
  res <- expect_s4_class(SS2OM(dir), "OM")

  out <- plot(res)
  expect_type(out, "list")

  # cpars objects are all of the correct dimension
  expect_true(all(vapply(res@cpars, get_nsim, numeric(1)) == res@nsim))

  expect_true(all(dim(res@cpars$Len_age) == c(res@nsim, res@maxage, res@proyears + res@nyears)))
  expect_true(all(res@cpars$Len_age > 0))

  expect_true(all(dim(res@cpars$Wt_age) == c(res@nsim, res@maxage, res@proyears + res@nyears)))
  expect_true(all(res@cpars$Wt_age > 0))

  expect_true(length(res@cpars$hs) == res@nsim)
  hmax <- ifelse(res@SRrel == 1L, 1, Inf)
  expect_true(all(res@cpars$hs >= 0.2 & res@cpars$h <= hmax))

  expect_true(all(dim(res@cpars$Perr_y) == c(res@nsim, res@maxage - 1 + res@proyears + res@nyears)))
  expect_true(all(res@cpars$Perr_y > 0))

  expect_true(all(dim(res@cpars$V) == c(res@nsim, res@maxage, res@proyears + res@nyears)))
  expect_true(all(res@cpars$V >= 0 & res@cpars$V <= 1))

  expect_true(all(dim(res@cpars$Find) == c(res@nsim, res@nyears)))
  expect_true(all(res@cpars$Find >= 0))

  expect_equal(DLMtool::cparscheck(res@cpars), res@nsim)

  expect_s4_class(runMSE(res, Hist = TRUE, silent = TRUE), "Hist")
})


test_that("GOM vermilion snapper", {
  dir <- file.path(DLM_dev_home, "DLMextra_Data/buildOMfromSS/Vermilion_Snapper_GOM")
  res <- expect_s4_class(SS2OM(dir), "OM")

  out <- plot(res)
  expect_type(out, "list")

  # cpars objects are all of the correct dimension
  expect_true(all(vapply(res@cpars, get_nsim, numeric(1)) == res@nsim))

  expect_true(all(dim(res@cpars$Len_age) == c(res@nsim, res@maxage, res@proyears + res@nyears)))
  expect_true(all(res@cpars$Len_age > 0))

  expect_true(all(dim(res@cpars$Wt_age) == c(res@nsim, res@maxage, res@proyears + res@nyears)))
  expect_true(all(res@cpars$Wt_age > 0))

  expect_true(length(res@cpars$hs) == res@nsim)
  hmax <- ifelse(res@SRrel == 1L, 1, Inf)
  expect_true(all(res@cpars$hs >= 0.2 & res@cpars$h <= hmax))

  expect_true(all(dim(res@cpars$Perr_y) == c(res@nsim, res@maxage - 1 + res@proyears + res@nyears)))
  expect_true(all(res@cpars$Perr_y > 0))

  expect_true(all(dim(res@cpars$V) == c(res@nsim, res@maxage, res@proyears + res@nyears)))
  expect_true(all(res@cpars$V >= 0 & res@cpars$V <= 1))

  expect_true(all(dim(res@cpars$Find) == c(res@nsim, res@nyears)))
  expect_true(all(res@cpars$Find >= 0))

  expect_equal(DLMtool::cparscheck(res@cpars), res@nsim)

  expect_s4_class(runMSE(res, Hist = TRUE, silent = TRUE), "Hist")
})



test_that("ICCAT yellowfin tuna", {
  dir <- file.path(DLM_dev_home, "DLMextra_Data/buildOMfromSS/Yellowfin_tuna_ICCAT")
  res <- expect_s4_class(SS2OM(dir), "OM")

  out <- plot(res)
  expect_type(out, "list")

  # cpars objects are all of the correct dimension
  expect_true(all(vapply(res@cpars, get_nsim, numeric(1)) == res@nsim))

  expect_true(all(dim(res@cpars$Len_age) == c(res@nsim, res@maxage, res@proyears + res@nyears)))
  expect_true(all(res@cpars$Len_age > 0))

  expect_true(all(dim(res@cpars$Wt_age) == c(res@nsim, res@maxage, res@proyears + res@nyears)))
  expect_true(all(res@cpars$Wt_age > 0))

  expect_true(length(res@cpars$hs) == res@nsim)
  hmax <- ifelse(res@SRrel == 1L, 1, Inf)
  expect_true(all(res@cpars$hs >= 0.2 & res@cpars$h <= hmax))

  expect_true(all(dim(res@cpars$Perr_y) == c(res@nsim, res@maxage - 1 + res@proyears + res@nyears)))
  expect_true(all(res@cpars$Perr_y > 0))

  expect_true(all(dim(res@cpars$V) == c(res@nsim, res@maxage, res@proyears + res@nyears)))
  expect_true(all(res@cpars$V >= 0 & res@cpars$V <= 1))

  expect_true(all(dim(res@cpars$Find) == c(res@nsim, res@nyears)))
  expect_true(all(res@cpars$Find >= 0))

  expect_equal(DLMtool::cparscheck(res@cpars), res@nsim)

  expect_s4_class(runMSE(res, Hist = TRUE, silent = TRUE), "Hist")
})
