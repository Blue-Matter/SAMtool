
#### For data-limited situations - see notes below

library(MSEtool)
library(testthat)

Data <- new("Data", "tests/Data_files/Data_ORH7A.csv")
Data <- new("Data", "tests/Data_files/GulfMenhaden_Data.csv")
Data <- new("Data", "tests/Data_files/MERA_TokyoBayData.csv")
Data <- new("Data", "tests/Data_files/Tiger_flathead_data.csv")
Data <- new("Data", "tests/Data_files/Yellow_striped_MPerf.csv")

testOM@maxage <- Data@MaxAge
SRA <- SRA_scope(testOM, Data)

generate_OM <- function(OM, Data) {
  Stock <- new("Stock")
  Stock@maxage <- Data@MaxAge
  Stock@R0 <- 1e4
  Stock@M <- rep(Data@Mort, 2)
  Msd <- rep(0, 2)
  h <- rep(Data@steep, 2)
  SRrel <- 1L
  Perr <- rep(Data@sigmaR, 2)
  Linf <- rep(Data@vbLinf, 2)
  K <- rep(Data@vbK, 2)
  t0 <- rep(Data@vbt0, 2)
  LenCV <- rep(Data@LenCV, 2)
  Ksd <- rep(0, 2)
  Linfsd <- rep(0, 2)




  return(OM)

}
