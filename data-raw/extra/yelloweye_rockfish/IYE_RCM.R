library(SAMtool)

############ Condition operating models with SRA_scope and data
#SRA_data <- readRDS("tests/Data_files/IYE_data.rds")
#SRA_data$s_CAA[!is.na(SRA_data$s_CAA[, 80, 1]), 81, 1] <- 0
#saveRDS(SRA_data, file = "tests/Data_files/IYE_data.rds")

#data_names <- c("Chist", "Index", "I_sd", "I_type", "length_bin", "s_CAA", "CAA", "CAL", "I_units")
#data_ind <- match(data_names, names(SRA_data))
#OM_condition <- readRDS("tests/Data_files/IYE_OM_2sim.rds")

# Base
SRA_data <- readRDS("tests/yelloweye_rockfish/IYE_data.rds")
OM_condition <- readRDS("tests/yelloweye_rockfish/IYE_OM_2sim.rds")

# Update OM to include Eobs slots
#OM <- new("OM")
#for(i in 1:length(slotNames(OM_condition))) {
#  if(!slotNames(OM_condition)[i] %in% c("Eobs", "Ebiascv")) slot(OM, slotNames(OM_condition)[i]) <- slot(OM_condition, slotNames(OM_condition)[i])
#}
#OM@Eobs <- c(0.1, 0.3)
#OM@Ebiascv <- 0.1

# IYE nsim = 250
OM_250 <- OM_condition
OM_250@nsim <- 250
set.seed(91283)
M_samps <- rlnorm(250, log(0.045) - 0.5 * 0.2^2, 0.2)
OM_250@cpars$M <- M_samps

# Sample steepness h ~ transformed beta with mean = 0.71, sd = 0.1.
# x = (h - 0.2)/0.8 where x ~ Beta(mean = 0.6375, sd = 0.12)
h_alpha <- alphaconv(0.6375, 0.12)
h_beta <- betaconv(0.6375, 0.12)
set.seed(65423)

h_samps <- rbeta(250, h_alpha, h_beta)
h_samps <- 0.8 * h_samps + 0.2
OM_250@cpars$h <- h_samps
OM_250@cpars$Mat_age <- OM_condition@cpars$Mat_age[1, , 1] %>% array(c(81, 250, 202)) %>% aperm(c(2, 1, 3))
saveRDS(OM_250, file = "tests/Data_files/IYE_OM_250sim.rds")


data_names <- c("Chist", "Index", "I_sd", "length_bin", "s_CAA", "CAA", "CAL", "I_units")

SRA <- RCM(OM = OM_condition, data = SRA_data[match(data_names, names(SRA_data))], condition = "catch2", selectivity = rep("free", 2),
           s_selectivity = rep("logistic", 5), #resample = TRUE, #cores = 1,
           vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 81, 2),
           map_ivul_par = SRA_data$map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
           #prior = list(M = c(0.02, 0.05)),
           LWT = list(CAL = 0, CAA = 0))
saveRDS(SRA, "tests/Data_files/IYE_SRA.rds")

SRA <- readRDS("tests/Data_files/IYE_SRA.rds")
ret <- retrospective(SRA, 3)

saveRDS(ret, "tests/Data_files/IYE_ret.rds")
ret <- readRDS("tests/Data_files/IYE_ret.rds")

plot(SRA)
plot(SRA, retro = ret)

Hist <- runMSE(SRA@OM, Hist = TRUE)

N <- SRA@mean_fit$report$N[1:102, ]
NOM <- Hist@AtAge$Number[1, , , ] %>% apply(c(1, 2), sum) %>% t()

# nsim = 250
OM_condition <- readRDS("tests/Data_files/IYE_OM_250sim.rds")

data_names <- c("Chist", "Index", "I_sd", "I_type", "length_bin", "s_CAA", "CAA", "CAL", "I_units")
data_ind <- match(data_names, names(SRA_data))

SRA <- RCM(OM = OM_condition, data = SRA_data[data_ind], condition = "catch2", selectivity = rep("free", 2),
           s_selectivity = rep("logistic", 5), cores = 4,
           vul_par = SRA_data$vul_par, map_vul_par = matrix(NA, 81, 2),
           map_s_vul_par = SRA_data$map_s_vul_par, map_log_rec_dev = SRA_data$map_log_rec_dev,
           LWT = list(CAL = 0, CAA = 0))
SRA@OM@MPA <- FALSE

set.seed(24)
SRA <- SRA_for_selectivity
nsim <- 250
AddIerr <- rnorm(nsim * 5 * (SRA@OM@nyears + SRA@OM@proyears), -0.5 * 0.25^2, 0.25) %>% exp() %>%
  array(dim = c(nsim, 5, SRA@OM@nyears + SRA@OM@proyears))

add_Ierr <- function(SRA) {
  AddIbeta <- matrix(1, nsim, 5)
  SRA@OM@cpars$AddIbeta <- AddIbeta
  SRA@OM@cpars$AddIerr <- AddIerr
  return(SRA)
}
SRA2 <- add_Ierr(SRA)
saveRDS(SRA2, file = "tests/Data_files/IYE_RCM_250sim.rds")

saveRDS(SRA2@OM, file = "yelloweye_OM.rds")