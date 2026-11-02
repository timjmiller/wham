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

test_that("legacy steepness retains native parameters without cached reports", {
  for (model in 3:4) {
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
      mature = array(rep(c(0.2, 0.8), each = 3L), c(1L, 3L, 2L)),
      fracyr_SSB = matrix(0.5, 3L, 1L)
    )
    waa <- array(1, c(1L, 3L, 2L))
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
    input$map <- list(mean_rec_pars = factor(c(1L, NA_integer_)))
    fit <- list()

    migrated <- migrate_legacy_steepness(input, fit, data, list(
      mean_rec_pars = matrix(legacy_pars, 1L)
    ))
    restored <- refit_legacy_steepness_parameters(
      migrated$par, migrated$legacy_steepness
    )

    expect_identical(migrated$data$recruit_model, model + 2L)
    expect_identical(migrated$map, input$map)
    expect_equal(as.vector(restored$mean_rec_pars), legacy_pars, tolerance = 1e-12)
    data$M_re_model <- 2L
    expect_identical(migrate_legacy_steepness(input, fit, data, list(
      mean_rec_pars = matrix(legacy_pars, 1L)
    ))$par, input$par)
  }

  data$do_proj <- 1L
  expect_error(
    migrate_legacy_steepness(input, fit, data, list(
      mean_rec_pars = matrix(c(0, log(100)), 1L)
    )),
    "without projections"
  )
})

test_that("native steepness preserves historical likelihoods and annual SPR", {
  for (case in c("steepness_bh", "steepness_ricker")) {
    path <- test_path("..", "..", "oldfitobjects", "release_refits",
      "v1.0.6", "features", paste0(case, "_original.rds"))
    skip_if_not(file.exists(path), "Historical steepness fixtures are unavailable")
    fit <- readRDS(path)
    input <- migrate_legacy_refit_input(fit)
    objective <- TMB::MakeADFun(input$data, input$par, random = input$random,
      map = input$map, DLL = "wham", silent = TRUE)
    nll <- objective$fn(objective$par)
    TMB::FreeADFun(objective)
    expect_equal(as.numeric(nll), as.numeric(fit$opt$objective),
      tolerance = 1e-6, scale = 1)

    for (shift in c(0, 0.15)) {
      pars <- input$par
      pars$Mpars[] <- pars$Mpars + shift
      data <- input$data
      if (shift > 0) {
        data$mature[1L, 2L, ] <- 0.9 * data$mature[1L, 2L, ]
        data$fracyr_SSB[2L, 1L] <- data$fracyr_SSB[2L, 1L] + 0.02
      }
      objective <- TMB::MakeADFun(data, pars,
        map = input$map, DLL = "wham", silent = TRUE)
      objective$fn(objective$par)
      report <- objective$report()
      TMB::FreeADFun(objective)
      spr <- vapply(seq_len(data$n_years_model), function(year) {
        mortality <- exp(report$log_M[1L, 1L, year, ])
        ages <- length(mortality)
        survival <- cumprod(c(1, exp(-mortality[-ages])))
        survival[ages] <- survival[ages] / (1 - exp(-mortality[ages]))
        sum(survival * data$mature[1L, year, ] *
          data$waa[data$waa_pointer_ssb[1L], year, ] *
          exp(-data$fracyr_SSB[year, 1L] * mortality))
      }, numeric(1))
      coefficients <- pars$mean_rec_pars[1L, ]
      steepness <- if (case == "steepness_bh") {
        0.2 + 0.8 * plogis(coefficients[1L])
      } else 0.2 + exp(coefficients[1L])
      if (case == "steepness_bh") {
        expected_a <- log(4 * steepness / ((1 - steepness) * spr))
        expected_b <- log((5 * steepness - 1) / ((1 - steepness) * spr)) - coefficients[2L]
      } else {
        expected_a <- 1.25 * log(5 * steepness) - log(spr)
        expected_b <- log(1.25 * log(5 * steepness) / spr) - coefficients[2L]
      }
      expect_equal(as.vector(report$log_SR_a), expected_a, tolerance = 1e-10)
      expect_equal(as.vector(report$log_SR_b), expected_b, tolerance = 1e-10)
      if (shift > 0) expect_gt(abs(expected_a[2L] - expected_a[1L]), 1e-6)
    }
  }
})
