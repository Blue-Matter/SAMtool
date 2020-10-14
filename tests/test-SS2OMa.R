
library(testthat)
library(MSEtool)


context("SS2OM Case studies")


DLM_dev_home <- "C:/~/DLMDev/Case_Studies"

dir_vec <- c("Blue_Shark_IO_IOTC/data/SS", "Bluefin_Tuna_EAtl_ICCAT/data/SS/EAST/Run82",
             "Bluefin_Tuna_WAtl_ICCAT/data/SS/WEST/WBFT12", "Bluefin_Tuna_WAtl_ICCAT/data/SS/WEST/WBFT13",
             "Red_Snapper_GOM_NOAA/data/SS", "Swordfish_NAtl_DFO/data/SS/SWO_AMO", "Yellowfin_Tuna_IO_IOTC/data/SS")

for(i in 1:length(dir_vec)) {
  test_that(dir_vec[i], {
    dir <- file.path(DLM_dev_home, dir_vec[i])
    res <- expect_s4_class(SS2OM(dir), "OM")

    #out <- plot(res)
    #expect_type(out, "list")

    # cpars objects are all of the correct dimension

    expect_true(all(dim(res@cpars$Len_age) == c(res@nsim, res@maxage, res@proyears + res@nyears)))
    expect_true(all(res@cpars$Len_age > 0))

    expect_true(all(dim(res@cpars$Wt_age) == c(res@nsim, res@maxage, res@proyears + res@nyears)))
    expect_true(all(res@cpars$Wt_age > 0))

    #expect_true(length(res@cpars$hs) == res@nsim)
    #hmax <- ifelse(res@SRrel == 1L, 1, Inf)
    #expect_true(all(res@cpars$hs >= 0.2 & res@cpars$h <= hmax))

    expect_true(all(dim(res@cpars$Perr_y) == c(res@nsim, res@maxage - 1 + res@proyears + res@nyears)))
    expect_true(all(res@cpars$Perr_y > 0))

    expect_true(all(dim(res@cpars$V) == c(res@nsim, res@maxage, res@proyears + res@nyears)))
    expect_true(all(res@cpars$V >= 0 & res@cpars$V <= 1))

    expect_true(all(dim(res@cpars$Find) == c(res@nsim, res@nyears)))
    expect_true(all(res@cpars$Find >= 0))

    #expect_true(length(res@EffYears) == length(res@EffUpper))
    #expect_true(length(res@EffYears) == length(res@EffLower))

    expect_true(length(res@M) == length(res@M2))

    expect_equal(DLMtool::cparscheck(res@cpars), res@nsim)


  })
}
