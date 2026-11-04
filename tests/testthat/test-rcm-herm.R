# Hermaphroditic complexes: fitted as one population, written back to each stock

# Single-area HermOM with data from one simulated history
make_herm_fixture <- function(nArea = 1) {
  om <- MSEtool::HermOM
  om@nSim <- 3
  if (nArea == 1) for (s in 1:2) om@Stock[[s]]@Spatial <- MSEtool::Spatial()
  hist <- MSEtool::Simulate(om, silent = TRUE)
  yrs <- MSEtool::Years(om, "Historical")

  r <- MSEtool::Removals(hist)
  r <- r[r$Sim == 1 & r$Period == "Historical", ]
  C <- tapply(r$Value, list(r$Year, r$Fleet), sum)
  b <- MSEtool::Biomass(hist)
  b <- b[b$Sim == 1 & b$Period == "Historical", ]
  B <- tapply(b$Value, b$Year, sum)

  dat <- MSEtool::Data(Years = yrs)
  dat@Landings <- MSEtool::CatchData(Name = colnames(C), Value = unname(C), CV = matrix(0.01, length(yrs), ncol(C)))
  dat@Survey <- MSEtool::IndicesData(Name = "Survey1", Value = matrix(B, ncol = 1), CV = matrix(0.1, length(yrs), 1))
  om@Data <- list(Grouper = dat)
  om
}

herm_om <- make_herm_fixture()
herm_rc <- suppressWarnings(RCM2(herm_om, s_selectivity = "B", silent = TRUE))

test_that("a hermaphroditic complex is fitted once and named by complex", {
  expect_s4_class(herm_rc, "rcmodel")
  expect_named(herm_rc@fits, "Grouper")
  expect_true(all(herm_rc@fits$Grouper@conv))
})

test_that("recruitment goes to the recruiting stock and fleets apply to every stock", {
  fit <- herm_rc@fits$Grouper
  R0 <- vapply(fit@report, getElement, numeric(1), "R0")
  expect_equal(as.numeric(herm_rc@OM@Stock[[1]]@SRR@R0), rep(R0, length.out = 3))
  expect_true(all(herm_rc@OM@Stock[[2]]@SRR@R0 == 0))

  sel1 <- herm_rc@OM@Fleet[[1]][[1]]@Selectivity@MeanAtAge
  sel2 <- herm_rc@OM@Fleet[[2]][[1]]@Selectivity@MeanAtAge
  expect_equal(sel1, sel2)
})

test_that("the conditioned hermaphroditic om simulates the fitted spawning output, numbers and catch", {
  fit <- herm_rc@fits$Grouper
  hist <- MSEtool::Simulate(herm_rc@OM, silent = TRUE)
  nyears <- herm_om@nYear

  sb <- MSEtool::SBiomass(hist)
  sb <- sb[sb$Period == "Historical" & sb$Stock == levels(sb$Stock)[1], ]
  sb <- tapply(sb$Value, list(sb$Sim, sb$Year), sum)
  expect_equal(unname(sb), unname(fit@SSB[, 1:nyears]), tolerance = 1e-4)

  N <- hist@Number
  Ntot <- apply(N[[1]][1, , , , drop = FALSE], 2:3, sum) + apply(N[[2]][1, , , , drop = FALSE], 2:3, sum)
  expect_equal(unname(t(Ntot)), unname(fit@report[[1]]$N[1:nyears, ]), tolerance = 1e-4)

  r <- MSEtool::Removals(hist)
  r <- r[r$Sim == 1 & r$Period == "Historical", ]
  C <- tapply(r$Value, list(r$Year, r$Fleet), sum)
  Cpred <- fit@report[[1]]$Cpred[1:nyears, ]
  expect_equal(unname(C[-1, ]), unname(Cpred[-1, ]), tolerance = 1e-3)
})

test_that("a hermaphroditic SSB index uses the fitted selectivity-at-age", {
  rc <- suppressWarnings(RCM2(herm_om, s_selectivity = "SSB", silent = TRUE))
  sel <- MSEtool::IndexSelectivity(rc@OM, "Survey1", Type = "Survey")
  expect_true(is.array(sel))
  expect_equal(unname(sel[1, ]), unname(rc@fits[[1]]@report[[1]]$ivul[herm_om@nYear, , 1]))
})

test_that("a multi-area hermaphroditic om is conditioned and simulates", {
  om <- make_herm_fixture(nArea = 2)
  rc <- suppressWarnings(RCM2(om, s_selectivity = "B", silent = TRUE))
  expect_s4_class(MSEtool::Simulate(rc@OM, silent = TRUE), "hist")
})

test_that("RCM2() builds the hermaphroditic om from stocks, data and Herm", {
  female <- herm_om@Stock[[1]]
  male <- herm_om@Stock[[2]]
  male@SRR@R0 <- 5000
  male@SRR@SPFrom <- NULL
  h <- herm_om@Herm[[1]]
  h@From <- "Female"
  h@To <- "Male"

  rc <- suppressWarnings(RCM2(list(Female = female, Male = male), data = herm_om@Data[[1]],
                              Fleet = herm_om@Fleet[[1]], Herm = h, nSim = 3, pYear = 5,
                              s_selectivity = "B", silent = TRUE))
  expect_named(rc@fits, "Female-Male")
  expect_true(all(rc@fits[[1]]@conv))
  expect_true(all(rc@OM@Stock$Male@SRR@R0 == 0))
  expect_equal(rc@OM@Stock$Male@SRR@SPFrom, "Female")
  expect_equal(rc@OM@nYear, herm_om@nYear)
  expect_equal(rc@OM@CurrentYear, herm_om@CurrentYear)

  hist <- MSEtool::Simulate(rc@OM, silent = TRUE)
  sb <- MSEtool::SBiomass(hist)
  sb <- sb[sb$Period == "Historical" & sb$Stock == "Female", ]
  sb <- tapply(sb$Value, list(sb$Sim, sb$Year), sum)
  expect_equal(unname(sb), unname(rc@fits[[1]]@SSB[, 1:herm_om@nYear]), tolerance = 1e-4)
})

test_that("complexes need Herm links and one recruiting stock", {
  om <- herm_om
  om@Herm <- list()
  expect_error(RCM2(om, silent = TRUE), "not linked")

  om <- herm_om
  om@Stock[[2]]@SRR@R0 <- 100
  expect_error(RCM2(om, silent = TRUE), "one stock with recruitment")

  om <- herm_om
  om@Complexes <- list()
  expect_error(RCM2(om, silent = TRUE), "OM@Complexes")
})
