# Plot WHAM output

Generates many output plots and tables for a fit WHAM model.

## Usage

``` r
plot_wham_output(
  mod,
  dir.main = getwd(),
  out.type = "html",
  res = 72,
  plot.opts = NULL
)
```

## Arguments

- mod:

  output from
  [`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md)

- dir.main:

  character, directory to save plots to (default =
  [`getwd()`](https://rdrr.io/r/base/getwd.html))

- out.type:

  character, either `'html'`, `'pdf'`, or `'png'` (default = `'html'`)

- res:

  resolution to save .png files (dpi)

- plot.opts:

  (optional) list of plot modifications

## Details

`out.type = 'html'` (default) makes a html file for viewing plot .png
files and html tables of parameter estimates in a browser.
`out.type = 'pdf'` makes one pdf file of all plots and tables.
`out.type = 'png'` creates a subdirectory \`plots_png“ in `dir.main` and
saves .png files within. `out.type = 'pdf' or 'png'` makes LaTeX and pdf
files of tables of parameter estimates. (tabs: 'input data',
'diagnostics', 'results', 'ref_points', 'retro', and 'misc').

`plot.opts` holds optional arguments to modify plots:

- `$ages.lab`:

  Character vector, will change age labels in plots (default is
  `1:n.ages`).

- `$font.family`:

  Font family, e.g. `"Times"`.

- `$browse`:

  T/F whether to open the html file in a browser. Default = T.

Plot functions are located in `wham_plots_tables.R` Table function is
located in `par_tables_fn.R`

## See also

[`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md),
[`wham_html`](https://timjmiller.github.io/wham/reference/wham_html.md),
`wham_plots_tables`

## Examples

``` r
if (FALSE) { # \dontrun{
data("input4_SNEMAYT") # load fit wham model
mod <- fit_wham(input4_SNEMAYT)
plot_wham_output(mod)
} # }
```
