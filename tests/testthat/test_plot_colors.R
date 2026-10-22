context("plot color palettes")

test_that("WHAM plot color modes return the expected palettes", {
  old.options <- options(wham.colors = "default")
  on.exit(options(old.options), add = TRUE)

  expect_equal(wham:::wham_palette(4), viridisLite::viridis(4))

  options(wham.colors = "bw")
  bw <- grDevices::col2rgb(wham:::wham_palette(4))
  expect_equal(bw[1, ], bw[2, ])
  expect_equal(bw[2, ], bw[3, ])

  options(wham.colors = "gray")
  gray <- grDevices::col2rgb(wham:::wham_palette(4))
  expect_equal(gray[1, ], gray[2, ])
  expect_equal(gray[2, ], gray[3, ])

  options(wham.colors = "cividis")
  expect_equal(wham:::wham_palette(4), viridisLite::cividis(4))

  options(wham.colors = "highcontrast")
  expect_length(wham:::wham_palette(4), 4)

  options(wham.colors = "pastel")
  expect_length(wham:::wham_palette(4), 4)

  options(wham.colors = "turbo")
  expect_equal(wham:::wham_palette(4), viridisLite::turbo(4, begin = 0.2, end = 0.8))
})