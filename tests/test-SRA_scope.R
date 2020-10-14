
#### For data-limited situations - see notes below

library(MSEtool)
library(testthat)

OM <- testOM
OM@nsim = 2
OM@Msd <- OM@Ksd <- OM@Linfsd <- c(0, 0)

Hist = runMSE(OM, Hist = TRUE)

bool <- c(TRUE, FALSE)

#### Function to run SRA_scope in parallel (with list)
data_master <- expand.grid(cond = c("catch", "effort"), Catch_with_effort = bool,
                           Index = bool, ML = bool, MW = bool, CAL = bool, CAA = bool)

data_master$MW[data_master$ML] <- NA
data_master$CAA[data_master$ML | data_master$MW] <- data_master$CAL[data_master$ML | data_master$MW] <- NA

lapply_fn <- function(i, data_master, OM, Hist) {
  args <- list(OM = OM)
  args$data <- list()

  if(data_master$cond[i] == "catch") { # Condition on catch
    args$condition <- "catch2"
    args$data$Chist <- Hist@Data@Cat[1, ]
  } else {
    if(data_master$Catch_with_effort[i]) { # Conditioned on effort, but may include some catches
      args$data$Chist <- Hist@Data@Cat[1, ]
      args$data$Chist[c(1:40)] <- NA
    }
    args$condition <- "effort"
    args$data$Ehist <- Hist@TSdata$Find[1, ]
  }

  if(data_master$Index[i]) { # Only recent ten years index
    args$data$Index <- Hist@Data@Ind[1, ]
    args$data$Index[c(1:20)] <- NA
    args$data$I_sd <- c(rep(NA, 20), rep(0.3, 30))
    args$s_selectivity <- "B"
  }

  if(data_master$ML[i]) { # Only recent mean length
    args$data$MS <- Hist@Data@ML[1, ]
    args$data$MS[c(1:40)] <- NA
    args$data$MS_type <- "length"
  } else if(data_master$MW[i]) {
    Data <- Hist@Data
    a <- Data@wla[1]
    b <- Data@wlb[1]
    WAL <- a * Data@CAL_mids^b

    MW <- colSums(t(Data@CAL[1, , ]) * WAL) / rowSums(Data@CAL[1, , ])
    args$data$MS <- MW
    args$data$MS_cv <- 0.2
    args$data$MS_type <- "weight"
  } else if(data_master$CAA[i]) {
    args$data$CAA <- Hist@Data@CAA[1, , ]
    args$data$CAA[1:35, ] <- NA
  } else if(data_master$CAL[i]) {

    args$data$CAL <- Hist@Data@CAL[1, , ]
    args$data$CAL[1:35, ] <- NA

    args$data$length_bin <- Hist@Data@CAL_mids
  }
  args$mean_fit <- TRUE
  SRA <- do.call(SRA_scope, args)

  return(SRA)
}


#### Function to run SRA_scope in parallel (with Data object)
data_master <- expand.grid(cond = c("catch", "effort"), Catch_with_effort = bool,
                           Index = bool, ML = bool, CAL = bool, CAA = bool)

data_master$CAA[data_master$ML] <- data_master$CAL[data_master$ML] <- NA

lapply_fn <- function(i, data_master, OM, Hist) {
  args <- list(OM = OM)
  args$data <- new("Data")

  if(data_master$cond[i] == "catch") { # Condition on catch
    args$condition <- "catch2"
    args$data@Cat <- Hist@Data@Cat[1, , drop = FALSE]
  } else {
    if(data_master$Catch_with_effort[i]) { # Conditioned on effort, but may include some catches
      args$data@Cat <- Hist@Data@Cat[1, , drop = FALSE]
      args$data@Cat[1, 1:40] <- NA
    }
    args$condition <- "effort"
    args$data@Effort <- Hist@TSdata$Find[1, , drop = FALSE]
  }

  if(data_master$Index[i]) { # Only recent ten years index
    args$data@Ind <- Hist@Data@Ind[1, , drop = FALSE]
    args$data@Ind[1, 1:20] <- NA
    args$data@CV_Ind <- c(rep(NA, 20), rep(0.3, 30)) %>% matrix(nrow = 1)
    args$s_selectivity <- "B"
  }

  if(data_master$ML[i]) { # Only recent mean length
    args$data@ML <- Hist@Data@ML[1, , drop = FALSE]
    args$data@ML[1, 1:40] <- NA
  } else if(data_master$CAA[i]) {
    args$data@CAA <- Hist@Data@CAA[1, , , drop = FALSE]
    args$data@CAA[, 1:35, ] <- NA
  } else if(data_master$CAL[i]) {

    args$data@CAL <- Hist@Data@CAL[1, , , drop = FALSE]
    args$data@CAL[, 1:35, ] <- NA

    args$data@CAL_mids <- Hist@Data@CAL_mids
  }
  args$mean_fit <- TRUE
  SRA <- do.call(SRA_scope, args)

  return(SRA)
}

DLMtool::setup(4)
sfExportAll()
sfLibrary(dplyr)



#### Run SRA
res <- sfClusterApplyLB(1:nrow(data_master), lapply_fn, data_master = data_master, Hist = Hist, OM = OM)
res <- lapply(5, lapply_fn, data_master = data_master, Hist = Hist, OM = OM)

# Test MW
res <- lapply_fn(9, data_master = data_master, Hist = Hist, OM = OM)


#### Test whether plot function works
lapply_fn_plot <- function(x, res) {
  try(R.utils::withTimeout(plot(res[[x]], compare = FALSE, open_file = FALSE, filename = as.character(x)),
                           timeout = 5),
      silent = TRUE)
}

#sfLibrary(R.utils)
#sfExport(list = c("lapply_fn_plot", "res"))
test_plot <- lapply(1:nrow(data_master), lapply_fn_plot, res = res)
saveRDS(test_plot, file = "tests/test_plot.rds")

error_fn <- function(x) {
  length(grep("elapsed time limit", x[1])) > 0 | length(grep("C:/", x[1])) > 0
}



data_master$conv <- lapply(res, function(x) sum(x@conv)/length(x@conv))
data_master$figure <- vapply(test_plot, error_fn, numeric(1))
View(data_master)

