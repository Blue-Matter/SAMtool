
library(SAMtool)

#assess <- c("DD_TMB", "DD_SS", "cDD", "cDD_SS", "SCA", "SCA2", "SCA_Pope", "SCA_RWM", "SP", "SP_SS")
#
#assess_sim <- lapply(assess, do.call, what = MSEtool::SimulatedData)
#names(assess_sim) <- "assess"
#saveRDS(assess_sim, file = "assess_sim.rds")
#
#assess_rs <- lapply(assess[!grepl("SCA", assess)], do.call, what = MSEtool::Red_snapper)
#names(assess_rs) <- assess[!grepl("SCA", assess)]
#saveRDS(assess_rs, file = "assess_rs.rds")

assess_sim <- readRDS("assess_sim.rds")
for(i in 1:length(assess_sim)) {
  print(names(assess_sim)[i])
  assess_sim@obj$retape()
  assess_sim@obj$fn()
}

assess_rs <- readRDS("assess_rs.rds")
for(i in 1:length(assess_rs)) {
  print(names(assess_rs)[i])
  assess_rs@obj$retape()
  assess_sim@obj$fn()
}