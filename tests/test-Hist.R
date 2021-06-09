

# Test that default configurations are very, very robust when running Assessments.
library(SAMtool)
OM <- testOM
OM@nsim <- 48
Hist <- runMSE(OM, Hist = TRUE)

assess <- avail("Assess")

MSEtool::setup(10)
sfExportAll()

# Generic run all assessments in default
message("Default assessment settings")
for(i in 1:length(assess)) {
  message(assess[i])
  run_mod <- sfClusterApplyLB(1:OM@nsim, assess[i], Data = Hist@Data)
  conv <- sapply(run_mod, getElement, 'conv')/length(run_mod)
  message("Percent converged... ", 100 * sum(conv)/length(run_mod), "\n")
}

# Run with depletion
message("\n\nAssessment dep 0.5")
do_dep <- sapply(assess, function(x) any(grepl("dep", names(formals(x)))))
for(i in 1:length(assess)) {
  if(do_dep[i]) {
    message(assess[i])
    run_mod <- sfClusterApplyLB(1:OM@nsim, assess[i], Data = Hist@Data, dep = 0.5, start = list(dep = 0.5))
    conv <- sapply(run_mod, getElement, 'conv')/length(run_mod)
    message("Percent converged... ", 100 * sum(conv)/length(run_mod), "\n")
  }
}

# Run with Ricker
message("\n\nDo Ricker")
do_Ricker <- sapply(assess, function(x) any(grepl("SR", names(formals(x)))))
for(i in 1:length(assess)) {
  if(do_Ricker[i]) {
    message(assess[i])
    run_mod <- sfClusterApplyLB(1:OM@nsim, assess[i], Data = Hist@Data, SR = "Ricker")
    conv <- sapply(run_mod, getElement, 'conv')/length(run_mod)
    message("Percent converged... ", 100 * sum(conv)/length(run_mod), "\n")
  }
}

# Shortened Index
message("\n\nShorten index")
Hist@Data@Ind[, 1:30] <- NA
for(i in 1:length(assess)) {
  message(assess[i])
  run_mod <- sfClusterApplyLB(1:OM@nsim, assess[i], Data = Hist@Data)
  conv <- sapply(run_mod, getElement, 'conv')/length(run_mod)
  message("Percent converged... ", 100 * sum(conv)/length(run_mod), "\n")
}

# Zero catch scenario

# SCA data weights

#  fix/est pars

# retro

# projections

# profile


