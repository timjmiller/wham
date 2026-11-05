test_that("legacy model codes are translated for the current WHAM objective", {
  input <- list(data = list(recruit_model = 3L))
  data <- list(
    M_model = 3L,
    recruit_model = 3L
  )

  migrated <- migrate_legacy_refit_model_codes(input, data)

  expect_identical(migrated$data$M_model, 2L)
  expect_identical(migrated$data$recruit_model, 3L)
})

test_that("weight-based mortality retains shared annual deviations", {
  pars <- list(
    M_re = matrix(c(0.1, 0.2, 0.3, 0.1, 0.2, 0.3), nrow = 3L)
  )
  old <- list(map = list(M_re = factor(c(1:3, 1:3))))
  input <- list(
    data = list(n_M_re = matrix(0L, 1L, 1L), log_b_model = 0L),
    par = list(M_re = array(0, c(1L, 1L, 3L, 2L))),
    map = list()
  )

  migrated <- migrate_legacy_weight_mortality(
    input, old, pars, list(n_ages = 2L)
  )

  expect_equal(migrated$par$M_re[1L, 1L, , 1L], pars$M_re[, 1L])
  expect_identical(migrated$data$n_M_re, matrix(1L, 1L, 1L))
  expect_identical(migrated$data$log_b_model, 1L)
  expect_identical(
    as.character(migrated$map$M_re)[1:3],
    as.character(old$map$M_re[1:3])
  )
  expect_true(all(is.na(migrated$map$M_re[4:6])))
})

test_that("legacy mortality comparisons use the original baseline", {
  pars <- list(Mpars = c(0.5, 0.6), M_a = c(0.1, 0.2))
  fit <- list(
    input = list(data = list(M_model = 3L), par = list(M0 = 0.3)),
    parList = list()
  )

  compared <- refit_legacy_mortality_parameters(pars, fit)

  expect_identical(compared$M0, 0.5)
  expect_equal(compared$M_a, c(-0.4, -0.3))
})
