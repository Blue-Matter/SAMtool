
library(SAMtool)

assess <- c("DD_TMB", "DD_SS", "cDD", "cDD_SS", "SCA", "SCA2", "SCA_Pope", "SCA_RWM", "SP", "SP_SS")

assess_sim <- lapply(assess, do.call, what = MSEtool::SimulatedData)
names(assess_sim) <- "assess"
saveRDS(assess_sim, file = "assess_sim.rds")

assess_rs <- lapply(assess[!grepl("SCA", assess)], do.call, what = MSEtool::Red_snapper)
names(assess_rs) <- assess[!grepl("SCA", assess)]
saveRDS(assess_rs, file = "assess_rs.rds")
