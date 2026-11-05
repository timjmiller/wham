test_that("legacy year-age output is normalized for comparison", {
  values <- matrix(seq_len(6), nrow = 2)
  mod <- list(
    is_sdrep = FALSE,
    env = list(data = list(n_ages = 3, percentSPR = 40)),
    years = 2001:2002,
    years_full = 2001:2002,
    rep = list(
      MAA = values,
      NAA = values,
      log_NAA_rep = log(values),
      log_SSB = log(c(10, 20)),
      log_F = log(matrix(c(0.1, 0.2)))
    )
  )
  result <- read_wham_fit(mod)
  expect_identical(dim(result$MAA), c(1L, 1L, 2L, 3L))
  expect_equal(result$MAA[1, 1, , ], values)
  for (name in c("NAA", "log_NAA_rep")) {
    expect_identical(dim(result[[name]]$est), c(1L, 1L, 2L, 3L))
    expect_equal(result[[name]]$est[1, 1, , ], mod$rep[[name]])
    expect_identical(dim(result[[name]]$se), dim(result[[name]]$est))
    expect_identical(dim(result[[name]]$ci$lo), dim(result[[name]]$est))
    expect_identical(dim(result[[name]]$ci$hi), dim(result[[name]]$est))
  }
  expect_equal(result$log_SSB_all$est, log(c(10, 20)))
  expect_equal(result$log_F_tot$est, log(c(0.1, 0.2)))
  expect_identical(dim(mod$rep$MAA), c(2L, 3L))

  mod$rep$MAA <- array(values, c(1, 1, 2, 3))
  expect_identical(read_wham_fit(mod)$MAA, mod$rep$MAA)
})

test_that("natural-scale reports supply missing log estimates without sdreport", {
  naa <- array(seq_len(12), c(1L, 1L, 2L, 6L))
  mod <- list(
    is_sdrep = FALSE,
    env = list(data = list(n_ages = 6, percentSPR = 40)),
    years = 2001:2002,
    years_full = 2001:2002,
    rep = list(
      MAA = naa * 0.01,
      NAA = naa,
      SSB = cbind(c(10, 20), c(30, 40)),
      log_F_tot = log(c(0.1, 0.2))
    )
  )
  result <- read_wham_fit(mod)
  expect_equal(result$log_SSB_all$est, log(c(40, 60)))
  expect_equal(result$log_NAA_rep$est, log(naa))
  for (name in c("log_SSB_all", "log_NAA_rep")) {
    zero <- result[[name]]$est
    zero[] <- 0
    expect_equal(result[[name]]$se, zero)
    expect_equal(result[[name]]$ci, list(lo = zero, hi = zero))
  }
  expect_identical(dim(result$log_NAA_rep$est), dim(naa))
  expect_null(mod$rep$log_NAA_rep)

  mod$rep$log_SSB_all <- log(c(100, 200))
  mod$rep$log_NAA_rep <- log(naa + 1)
  result <- read_wham_fit(mod)
  expect_equal(result$log_SSB_all$est, mod$rep$log_SSB_all)
  expect_equal(result$log_NAA_rep$est, mod$rep$log_NAA_rep)

  mod$rep$log_SSB_all <- mod$rep$log_NAA_rep <- NULL
  mod$rep$MAA <- matrix(0.2, 2, 6)
  mod$rep$NAA <- matrix(seq_len(12), 2, 6)
  result <- read_wham_fit(mod)
  expect_identical(dim(result$log_NAA_rep$est), c(1L, 1L, 2L, 6L))
  expect_equal(result$log_NAA_rep$est[1, 1, , ], log(mod$rep$NAA))
})

test_that("saved legacy fits retain status estimates and covariance", {
  for (model in c("Butterfishfit", "GBK_Cod_2024MTfit")) {
    path <- test_path("..", "..", "oldfitobjects", paste0(model, ".rds"))
    skip_if_not(file.exists(path), "Local legacy fit is unavailable")
    mod <- readRDS(path)
    result <- read_wham_fit(mod)
    summary <- summary(mod$sdrep, "report")
    indices <- c(
      which(rownames(summary) == "log_SSB")[1],
      which(rownames(summary) == "log_SSB_FXSPR")[1],
      which(rownames(summary) == "log_F")[1],
      which(rownames(summary) == "log_FXSPR")[1]
    )
    contrast <- cbind(c(1, -1, 0, 0), c(0, 0, 1, -1))
    expect_equal(
      result$log_rel_ssb_F_cov[[1]][[1]],
      t(contrast) %*% summary[indices, 1]
    )
    expect_equal(
      result$log_rel_ssb_F_cov[[1]][[2]],
      t(contrast) %*% mod$sdrep$cov[indices, indices] %*% contrast
    )
    expect_equal(result$log_SSB_all$est, as.vector(result$log_SSB$est))
    expect_equal(result$log_SSB_all$se, as.vector(result$log_SSB$se))
    expect_equal(result$log_SSB_all$ci$lo, as.vector(result$log_SSB$ci$lo))
    expect_equal(result$log_SSB_all$ci$hi, as.vector(result$log_SSB$ci$hi))
    expect_equal(result$log_F_tot$est, as.vector(result$log_F$est))
    expect_equal(result$log_F_tot$se, as.vector(result$log_F$se))
    expect_equal(
      result$log_NAA_rep$est[1, 1, , ],
      as.list(mod$sdrep, report = TRUE, what = "Est")$log_NAA_rep
    )
  }
})
