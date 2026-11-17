# btime <- Sys.time(); devtools::test(filter = "mohns_rho"); etime <- Sys.time(); runtime = etime - btime; runtime;

context("Mohn's rho")
test_that("Mohn's rho uses peel terminal years and regional fishing mortality", {
  report <- list(
    SSB = matrix(seq_len(10) * 100, nrow = 5),
    Fbar = cbind(rep(100, 5), rep(200, 5), 1:5, 6:10),
    NAA = array(seq_len(60), dim = c(2, 2, 5, 3))
  )
  report$NAA[2, 1, , 2] <- 0
  # Non-terminal years and fleet Fbar columns must not enter the averages.
  peel1 <- list(rep = lapply(report, function(x) x * 10))
  peel2 <- list(rep = lapply(report, function(x) x * 10))
  peel1$rep$SSB[4, ] <- report$SSB[4, ] * c(1.2, 0.8)
  peel2$rep$SSB[3, ] <- report$SSB[3, ] * c(0.9, 1.4)
  peel1$rep$Fbar[4, 3:4] <- report$Fbar[4, 3:4] * c(1.4, 0.7)
  peel2$rep$Fbar[3, 3:4] <- report$Fbar[3, 3:4] * c(0.8, 1.1)
  peel1$rep$NAA[, , 4, ] <- report$NAA[, , 4, ] * 1.2
  peel2$rep$NAA[, , 3, ] <- report$NAA[, , 3, ] * 0.9
  naa_where <- array(1, dim = c(2, 2, 3))
  naa_where[2, 1, 2] <- 0
  model <- list(
    env = list(data = list(
      n_years_model = 5,
      n_regions = 2,
      n_stocks = 2,
      n_fleets = 2,
      n_ages = 3,
      NAA_where = naa_where
    )),
    rep = report,
    peels = list(peel1, peel2),
    ages.lab = c("1", "2", "3+")
  )
  rho <- mohns_rho(model)
  expect_named(rho, c("SSB", "Fbar", "naa"), ignore.order = TRUE)
  expect_equal(rho$SSB, c(0.05, 0.1))
  expect_equal(rho$Fbar, c(0.1, -0.1))
  expected_naa <- array(0.05, dim = c(2, 2, 3),
                        dimnames = list(NULL, NULL, c("R", "N2", "N3+")))
  expected_naa[2, 1, 2] <- NA_real_
  expect_equal(rho$naa, expected_naa)
  expect_identical(unname(rho$naa[2, 1, 2]), NA_real_)
})
test_that("Mohn's rho rejects models without retrospective peels", {
  expect_error(mohns_rho(list(peels = list())),
               "There are no peels in this model", fixed = TRUE)
  expect_error(mohns_rho(list()),
               "There are no peels in this model", fixed = TRUE)
})
