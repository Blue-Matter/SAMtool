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
SRA@OM@MPA <- FALSE
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
# Let's play with some priors






########
# Cod
#args <- readRDS("C:/~/NE_retro/GoM_cod/args_cod_M02.rds")
#newOM <- readRDS("tests/Data_files/cod_OM.rds")
#for(i in 1:length(slotNames(newOM))) if(slotNames(newOM)[i] != "MPA") slot(newOM, slotNames(newOM)[i]) <- slot(args$OM, slotNames(newOM)[i])
#
#args$OM <- newOM
#args$OM@DR <- rep(0, 2)
#args$OM@MPA <- FALSE
#
#args$data$CAA <- array(0, c(37, 1, 1)) %>% abind::abind(args$data$CAA, along = 2)
#args$data$s_CAA <- array(0, c(37, 1, 3)) %>% abind::abind(args$data$s_CAA, along = 2)
#
#args$OM@cpars$M_ageArray <- array(0.2, c(100, 1, 87)) %>% abind::abind(args$OM@cpars$M_ageArray, along = 2)
#args$OM@cpars$Mat_age <- array(0, c(100, 1, 87)) %>% abind::abind(args$OM@cpars$Mat_age, along = 2)
#
#args$OM@cpars$Len_age <- abind::abind(args$OM@cpars$Len_age, array(10, c(100, 1, 87)), along = 2)
#args$OM@cpars$Wt_age <- array(0.1, c(100, 1, 87)) %>% abind::abind(args$OM@cpars$Wt_age, along = 2)
#args$map_log_early_rec_dev <- 1:9
#args$s_vul_par <- matrix(0, 1, 3) %>% rbind(args$s_vul_par)
#args$map_s_vul_par <- matrix(NA, 1, 3) %>% rbind(args$map_s_vul_par)

#saveRDS(args, file = "tests/Data_files/cod_args.rds")

args <- readRDS(file = "tests/Data_files/cod_args.rds")
args <- readRDS(file = "tests/Data_files/cod_args_MRAMP.rds")
#args$resample <- TRUE

SRA <- do.call(RCM, args)
saveRDS(SRA, "tests/Data_files/cod_SRA.rds")
SRA <- readRDS("tests/Data_files/cod_SRA.rds")


ret <- retrospective(SRA, 7, figure = FALSE)
plot(ret)

plot(SRA)

Hist <- runMSE(SRA@OM, Hist = TRUE)


N <- SRA@mean_fit$report$N[1:37, ]
NOM <- Hist@AtAge$Number[1, , , ] %>% apply(c(1, 2), sum) %>% t()

View(N - NOM)

plot(N[, 1], typ = 'o')
lines(NOM[, 1], col = 'red')

plot(N[, 1]/NOM[, 1] - 1)

Hist@Data@Misc$StockPars$R0a

Hist@Data@Misc$StockPars$SSBpR





###### Condition on effort with patchy catch
datS<-readRDS("tests/snapper/datS.rda")
OM<-readRDS("tests/snapper/OM.rda")

input<-list()
input$ESS<-100
input$wt_comp<-1
input$max_F<-5
input$C_eq_val=0
selectivity="logistic"
ncpus=1
condition="effort"

out<-RCM(OM,datS,
         ESS = rep(input$ESS,2),
         LWT = list(CAA=input$Wt_comp,CAL=input$Wt_comp),
         max_F = input$max_F,
         C_eq=input$C_eq_val,
         selectivity=selectivity,
         condition=condition,
         cores = ncpus,
         mean_fit=TRUE,
         #drop_nonconv=TRUE,
         control=list(eval.max=1E4, iter.max=1E4, abs.tol=1e-6))

plot(out)

### Compare RCM to OM
#OM <- out@OM
#Hist <- runMSE(out@OM, Hist = TRUE)
#
#N_RCM <- lapply(out@Misc, getElement, "N") %>% simplify2array()
#N_Hist <- Hist@AtAge$Number %>% apply(c(1, 2, 3), sum)
#
#range(OM@cpars$Find - Hist@SampPars$Fleet$Find * Hist@SampPars$Fleet$qs)
#range(OM@cpars$V - Hist@AtAge$Select)
#OM@cpars$R0 - Hist@SampPars$Stock$R0
#
#SSB_RCM <- lapply(out@Misc, getElement, "E") %>% simplify2array()
#SSB_Hist <- Hist@AtAge$SBiomass %>% apply(c(1, 3), sum)
#
#OM2 <- out@OM
#OM2@cpars$initD <- NULL
#OM2@cpars$qs <- NULL
#OM2@cpars$Fdisc <- rep(0, 24)
#OM2@cpars$Asize <- OM2@cpars$control <- NULL
#OM2@Prob_staying <- OM2@Size_area_1 <- OM2@Frac_area_1 <- rep(0.5, 2)
#OM2@cpars$mov <- OM2@cpars$MPA <- NULL
#
#out2 <- out
#out2@OM <- OM2
#
#plot(out2)
#Hist2 <- runMSE(out@OM, Hist = TRUE)



###### Issue with catch2?
dat <- readRDS("tests/tiger_flathead/dat.rda")
OM <- readRDS("tests/tiger_flathead/OM.rda")
out <- RCM(OM, dat, mean_fit = TRUE, cores = 4, condition = 'catch2')
