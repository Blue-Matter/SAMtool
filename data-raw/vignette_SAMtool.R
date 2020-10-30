
library(SAMtool)

# Figure 1. MSEtool_design.png is a Powerpoint slide.
# Figure 2. SAMtool_design.png is a Powerpoint slide.

SCA_assessment <- SCA(x = 3, Data = MSEtool::SimulatedData, fix_h = FALSE)
summary(SCA_assessment)

png("vignettes/fig_assess/retrospective.png", units = "in", res = 400, width = 6, height = 4)
par(mar = c(5, 4, 1, 1))
retrospective(SCA_assessment)
dev.off()

png("vignettes/fig_assess/profile.png", units = "in", res = 400, width = 6, height = 4)
profile(SCA_assessment, R0 = seq(600, 1000, 50), h = seq(0.75, 0.95, 0.02))
dev.off()


SCA_MSY_ <- make_MP(SCA, HCR_MSY, diagnostic = "full")
SCA_4010_ <- make_MP(SCA, HCR40_10, diagnostic = "full")
#myMSE <- MSEtool::runMSE(OM = MSEtool::testOM, MPs = c("FMSYref", "AvC", "SCA_MSY_", "SCA_4010_"))

saveRDS(myMSE, file = "data-raw/myMSE.rds")

rr <- Tplot(myMSE)
rr$Plots[[3]]
ggplot2::ggsave("vignettes/fig_assess/Tplot.png", width = 5, height = 3.5)

myMSE <- readRDS("data-raw/myMSE.rds")
png("vignettes/fig_assess/diagnostic_AM.png", units = "in", res = 400, width = 9, height = 6)
diagnostic_AM(myMSE, "SCA_MSY_")
dev.off()

png("vignettes/fig_assess/retrospective_AM.png", units = "in", res = 400, width = 8.5, height = 5.5)
retrospective_AM(myMSE, "SCA_MSY_", sim = 3)
dev.off()

rr <- TAC(SimulatedData, c("SCA_MSY_", "SCA_4010_"))
