



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

