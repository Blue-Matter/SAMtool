
library(testthat)
library(MSEtool)


context("SS2Data Case studies")

DLM_dev_home <- "C:/~/DLMDev/Case_Studies"

dir_vec <- c("Blue_Shark_IO_IOTC/data/SS", "Bluefin_Tuna_EAtl_ICCAT/data/SS/EAST/Run82",
             "Bluefin_Tuna_WAtl_ICCAT/data/SS/WEST/WBFT12", "Bluefin_Tuna_WAtl_ICCAT/data/SS/WEST/WBFT13",
             "Red_Snapper_GOM_NOAA/data/SS", "Swordfish_NAtl_DFO/data/SS/SWO_AMO", "Yellowfin_Tuna_IO_IOTC/data/SS")

for(i in 1:length(dir_vec)) {
  test_that(dir_vec[i], {
    dir <- file.path(DLM_dev_home, dir_vec[i])
    res <- expect_s4_class(SS2Data(dir), "Data")

    nyears <- length(res@Year)

    expect_message(plot(res))

    expect_true(dim(res@Cat)[2] == nyears)

    nindex <- length(res@AddIndType)
    expect_true(all(dim(res@AddInd[1, , ]) == c(nindex, nyears)))
    expect_true(all(dim(res@CV_AddInd[1, , ]) == c(nindex, nyears)))
    expect_true(all(dim(res@AddIndV[1, , ]) == c(nindex, res@MaxAge)))
    expect_true(length(res@AddIunits) == nindex)

    if(!all(is.na(res@CAL))) {
      nbins <- length(res@CAL_mids)
      expect_true(!all(is.na(res@CAL_mids)))
      expect_true(!all(is.na(res@CAL_bins)))

      expect_true(length(res@CAL_mids) == nbins)
      expect_true(length(res@CAL_bins) == nbins + 1)
      expect_true(all(dim(res@CAL[1, , ]) == c(nyears, nbins)))
    }

    if(!all(is.na(res@CAA))) expect_true(all(dim(res@CAA[1, , ]) == c(nyears, res@MaxAge)))

    expect_true(dim(res@Rec)[2] == nyears)
    expect_true(res@t == nyears)

    expect_true(length(res@Mort) == 1 & res@Mort > 0)

    expect_true(length(res@FMSY_M) == 1 & !is.na(res@FMSY_M) & res@FMSY_M > 0)
    expect_true(length(res@BMSY_B0) == 1 & !is.na(res@BMSY_B0) & res@BMSY_B0 > 0)

    expect_true(length(res@L50) == 1 & !is.na(res@L50) & res@L50 < res@vbLinf)
    expect_true(length(res@L95) == 1 & !is.na(res@L95) & res@L95 < res@vbLinf)

    expect_true(res@L50 < res@L95)

    expect_true(all(is.na(res@ML)) | dim(res@ML)[2] == nyears)
    expect_true(all(is.na(res@Lbar)) | dim(res@Lbar)[2] == nyears)
    expect_true(all(is.na(res@Lc)) | dim(res@Lc)[2] == nyears)

    expect_true(length(res@LFC) == 1 & !is.na(res@LFC) & res@LFC < res@vbLinf)
    expect_true(length(res@LFS) == 1 & !is.na(res@LFS) & res@LFS < res@vbLinf)

    expect_true(res@LFC < res@LFS)
    expect_true(length(res@steep) == 1 & !is.na(res@steep) & res@steep >= 0.2)

  })
}
