

# Test that default configurations are very, very robust when running MSEs.
library(testthat)
library(MSEtool)


DD_MP <- make_MP(DD_TMB, HCR_MSY, diagnostic = 'min')
DDSS_MP <- make_MP(DD_SS, HCR_MSY, diagnostic = 'min')
SP_MP <- make_MP(SP, HCR_MSY, diagnostic = 'min')
SPSS_MP <- make_MP(SP_SS, HCR_MSY, diagnostic = 'min')
SCA_MP_dome <- make_MP(SCA, HCR_MSY, diagnostic = 'min', vulnerability = 'dome')
SCA_MP <- make_MP(SCA, HCR_MSY, diagnostic = 'min')

MP_vec <- c("DD_MP", "DDSS_MP", "SP_MP", "SPSS_MP", "SCA_MP_dome", "SCA_MP")

testOM@nsim <- 200


for(i in MP_vec) {
  setup(12)
  tim <- proc.time()
  res <- runMSE(testOM, MPs = i, PPD = TRUE, parallel = TRUE)
  message(paste0("Run time: ", (proc.time() - tim)[3], " seconds"))
  sfStop()
  cat(diagnostic_AM(res, figure = FALSE))
  save(res, file = paste0("tests_results/", i, ".RData"))
}
