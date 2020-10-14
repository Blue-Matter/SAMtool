

# Test that default configurations are very, very robust when running Assessments.
library(testthat)
library(MSEtool)
library(DLMextra)

OMs <- avail("OM")

# Make hist objects once
#setup(12)
#sfLibrary(DLMextra)
#Hist <- sfClusterApplyLB(OMs, function(x) runMSE(get(x), Hist = TRUE))
#sfStop()
#Data_from_Hist <- lapply(Hist, getElement, "Data")
#save(Data_from_Hist, file = "tests_results/Data_from_Hist.RData")

# Run test
context("Run Assess models from Hist object")

load("tests_results/Data_from_Hist.RData")

# Other OMs need to be evaluated separately
for(i in 1:length(Data_from_Hist)) {

  #test_that(OMs[i], {
    message(paste("OM:", OMs[i], "\n"))

    Data <- Data_from_Hist[[i]]

    # Main call: change desired Assess model as necessary here.
    res <- prelim_AM(Data, SCA, 12)

  #})

}

sfStop()

