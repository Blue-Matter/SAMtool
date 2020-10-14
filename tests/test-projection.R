
gen_plots <- function() {
  matplot(t(pro$Cpred), typ = 'l', ylim = c(0, 1.1 * max(pro$Cpred)))
  matplot(t(pro$F), typ = 'l', ylim = c(0, 1.1 * max(pro$F)))
  matplot(t(pro$B), typ = 'l', ylim = c(0, 1.1 * max(pro$B)))
  title(mods[i])
}

mods <- avail("Assess")

for(i in 1:length(mods)) {
  res <- do.call(mods[i], list(x = 1, Data = SimulatedData))

  par(mfrow = c(3,3))

  pro <- MSEtool:::projection(res, FMort = ifelse(length(res@FMSY) > 0, res@FMSY, res@UMSY))
  gen_plots()

  pro <- MSEtool:::projection(res, FMort = 0.25 * ifelse(length(res@FMSY) > 0, res@FMSY, res@UMSY))
  gen_plots()

  pro <- MSEtool:::projection(res, FMort = 2.5 * ifelse(length(res@FMSY) > 0, res@FMSY, res@UMSY))
  gen_plots()


  pro <- MSEtool:::projection(res, Catch = 10, constrain = "Catch")
  gen_plots()

  pro <- MSEtool:::projection(res, Catch = 100, constrain = "Catch")
  gen_plots()

  pro <- MSEtool:::projection(res, Catch = 100000, constrain = "Catch")
  gen_plots()
}


