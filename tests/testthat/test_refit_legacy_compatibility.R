test_that("legacy model codes are translated for the current WHAM objective", {
  input <- list(data = list(recruit_model = 3L))
  data <- list(
    M_model = 3L,
    use_steepness = 1L,
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

test_that("legacy steepness transforms to and from current recruitment models", {
  for (model in 3:4) {
    spr <- 0.2 * exp(-0.5 * 0.2) +
      exp(-0.2) / (1 - exp(-0.1)) * 0.8 * exp(-0.5 * 0.1)
    steepness <- 0.65
    r0 <- 100
    legacy_pars <- if (model == 3L) {
      c(qlogis((steepness - 0.2) / 0.8), log(r0))
    } else {
      c(log(steepness - 0.2), log(r0))
    }
    data <- list(
      n_stocks = 1L, n_regions = 1L, n_seasons = 1L,
      n_years_model = 3L, n_ages = 2L, n_M_re = 1L,
      M_re_model = 1L, M_model = 1L, use_steepness = 1L,
      recruit_model = model, do_proj = 0L,
      mature = array(rep(c(0.2, 0.8), 3L), c(1L, 3L, 2L)),
      fracyr_SSB = matrix(0.5, 3L, 1L)
    )
    waa <- array(1, c(1L, 3L, 2L))
    mortality <- array(
      rep(log(c(0.2, 0.1)), 3L),
      c(1L, 1L, 3L, 2L)
    )
    input <- list(
      data = list(
        n_stocks = 1L,
        n_regions = 1L,
        n_seasons = 1L,
        fracyr_SSB = matrix(0.5, 3L, 1L),
        waa = waa,
        waa_pointer_ssb = 1L,
        waa_pointer_M = 1L,
        Ecov_how_M = array(0L, c(1L, 1L, 2L, 1L))
      ),
      par = list(mean_rec_pars = matrix(legacy_pars, 1L))
    )
    fit <- list(rep = list(log_M = mortality))

    migrated <- migrate_legacy_steepness(input, fit, data, list(
      mean_rec_pars = matrix(legacy_pars, 1L)
    ))
    restored <- refit_legacy_steepness_parameters(
      migrated$par, migrated$legacy_steepness
    )

    expect_identical(migrated$data$recruit_model, model)
    expect_equal(migrated$legacy_steepness$spr, spr, tolerance = 1e-12)
    expect_equal(as.vector(restored$mean_rec_pars), legacy_pars, tolerance = 1e-12)
  }

  data$n_M_re <- 1L
  data$M_re_model <- 2L
  expect_error(
    migrate_legacy_steepness(input, fit, data, list(
      mean_rec_pars = matrix(c(0, log(100)), 1L)
    )),
    "time-varying mortality random effects"
  )
})
