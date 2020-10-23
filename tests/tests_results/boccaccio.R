library(SAMtool)
load('boccaccio_OM.RData')

sfInit(parallel = TRUE, cpus = 8)
sfLibrary(SAMtool)

start_time <- proc.time()
res <- sfLapply(1:48, SCA, Data = Data, start = list(R0 = 10 * OM@R0), vulnerability = "dome", early_dev = "comp")
message(paste("Timing:", round((proc.time() - start_time)[3], 2), "seconds"))

conv <- vapply(res, function(x) if(is.character(x@SD)) return(FALSE) else return(x@SD$pdHess), logical(1))

is_na <- sum(conv)
message(paste0(is_na, " didn't converge out of ", OM@nsim, " (", round(100 * is_na/OM@nsim, 1), " %)"))

sfStop()

save(res, file = "boccaccio5.RData")

# 1 = Use hessian for converence
# 2 = No hessian, 3 restarts
# 3 = Dome
# 4 = Dome, no estimate early_dev
# 5 = Dome, re-do log_rec_dev
