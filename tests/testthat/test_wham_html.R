# btime <- Sys.time(); devtools::test(filter = "wham_html"); etime <- Sys.time(); runtime = etime - btime; runtime;

context("HTML output")

test_that("HTML output links plot tabs and respects the title and width", {
  output_dir <- tempfile("wham-html-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
  tabs <- c("diagnostics", "input_data", "misc", "ref_points", "results", "retro")
  labels <- c("Diagnostics", "Input Data", "Misc", "Reference Points",
              "Results", "Retrospective")
  for (tab in tabs) {
    dir.create(file.path(output_dir, "plots_png", tab), recursive = TRUE)
    file.create(file.path(output_dir, "plots_png", tab, "plot.png"))
  }
  file.create(file.path(output_dir, "plots_png", "diagnostics", "summary_text.png"))
  opened <- character()
  old_options <- options(browser = function(url, ...) opened <<- c(opened, url))
  on.exit(options(old_options), add = TRUE)
  wham:::wham_html(output_dir, title = "Test assessment", width = 400,
                  openfile = FALSE)
  home_path <- file.path(output_dir, "wham_output.html")
  home <- paste(readLines(home_path, warn = FALSE), collapse = "\n")
  expect_match(home, "<title>Test assessment</title>", fixed = TRUE)
  expect_match(home, "summary_text.png", fixed = TRUE)
  expect_match(home, "width=680", fixed = TRUE)
  expect_length(opened, 0)
  for (i in seq_along(tabs)) {
    tab_path <- file.path(output_dir, "html", paste0(tabs[i], ".html"))
    expect_true(file.exists(tab_path))
    expect_match(home, paste0(">", labels[i], "</a>"), fixed = TRUE)
    page <- paste(readLines(tab_path, warn = FALSE), collapse = "\n")
    expect_match(page, home_path, fixed = TRUE)
    expect_match(page, file.path(output_dir, "plots_png", tabs[i], "plot.png"),
                 fixed = TRUE)
    expect_match(page, "width=400", fixed = TRUE)
    expect_match(page, "</html>", fixed = TRUE)
  }
  expect_message(wham:::wham_html(output_dir), "Opening HTML file")
  expect_identical(opened, home_path)
})
