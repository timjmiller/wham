test_that("aic reports the marginal value and its metadata", {
  model <- list(
    TMB_version = "test",
    opt = list(objective = 12),
    par = c(a = 1, b = 2),
    env = list(data = list(obs = list(val = c(1, NA, 2, NA, 3))))
  )

  result <- aic(model)

  expect_equal(unclass(result), 28)
  expect_identical(attr(result, "df"), 2L)
  expect_identical(attr(result, "n"), 3L)
  expect_identical(attr(result, "type"), "marginal")
})

test_that("aic rejects objects without TMB model metadata", {
  expect_error(aic(list()), "TMB_version does not exist")
})

test_that("mohns_rho summarizes peels and masks unavailable abundance states", {
  base_naa <- array(10, c(1L, 1L, 5L, 3L))
  first_naa <- base_naa
  first_naa[1L, 1L, 4L, 1L] <- 11
  first_naa[1L, 1L, 4L, 2L] <- 9
  second_naa <- base_naa
  second_naa[1L, 1L, 3L, 1L] <- 8
  second_naa[1L, 1L, 3L, 2L] <- 12

  base_fbar <- matrix(2, 5L, 2L)
  first_fbar <- base_fbar
  first_fbar[4L, 2L] <- 2.2
  second_fbar <- base_fbar
  second_fbar[3L, 2L] <- 1.8
  model <- list(
    env = list(data = list(
      n_years_model = 5L,
      n_regions = 1L,
      n_stocks = 1L,
      n_fleets = 1L,
      n_ages = 3L,
      NAA_where = array(c(1L, 1L, 0L), c(1L, 1L, 3L))
    )),
    rep = list(SSB = matrix(100, 5L, 1L), Fbar = base_fbar, NAA = base_naa),
    peels = list(
      list(rep = list(
        SSB = matrix(100, 5L, 1L),
        Fbar = first_fbar,
        NAA = first_naa
      )),
      list(rep = list(
        SSB = matrix(100, 5L, 1L),
        Fbar = second_fbar,
        NAA = second_naa
      ))
    ),
    ages.lab = c("1", "2", "3")
  )
  model$peels[[1L]]$rep$SSB[4L, 1L] <- 110
  model$peels[[2L]]$rep$SSB[3L, 1L] <- 80

  result <- mohns_rho(model)

  expect_equal(result$SSB, -0.05)
  expect_equal(result$Fbar, 0)
  expect_equal(result$naa[1L, 1L, 1:2], c(-0.05, 0.1))
  expect_true(is.na(result$naa[1L, 1L, 3L]))
  expect_identical(dimnames(result$naa)[[3L]], c("R", "N2", "N3"))
})

test_that("mohns_rho requires at least one peel", {
  model <- list(
    env = list(data = list(
      n_years_model = 2L,
      n_regions = 1L,
      n_stocks = 1L,
      n_fleets = 1L,
      n_ages = 1L
    )),
    peels = list()
  )

  expect_error(mohns_rho(model), "There are no peels")
})
