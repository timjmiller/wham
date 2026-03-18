# Compare multiple WHAM (or ASAP) models

After fitting multiple WHAM (or ASAP) models, `compare_wham_models`
produces plots and a table of AIC and Mohn's rho to aid model
comparison.

## Usage

``` r
compare_wham_models(
  mods,
  do.table = TRUE,
  do.plot = TRUE,
  fdir = getwd(),
  compare.opts = NULL,
  table.opts = NULL,
  plot.opts = NULL,
  fname = NULL,
  sort = NULL,
  calc.rho = NULL,
  calc.aic = NULL,
  do.print = NULL
)
```

## Arguments

- mods:

  (named) list of fit WHAM/ASAP models. To read in ASAP model output,
  use
  [`read_asap3_fit`](https://timjmiller.github.io/wham/reference/read_asap3_fit.md).
  If no names are given, 'm1', 'm2', ... will be used.

- do.table:

  T/F, produce table of AIC and/or Mohn's rho? Default = TRUE.

- do.plot:

  T/F, produce plots? Default = TRUE.

- fdir:

  character, path to directory to save table and/or plots. Default =
  [`getwd()`](https://rdrr.io/r/base/getwd.html).

- compare.opts:

  list of options to generate comparison results:

  `$stock`

  :   integer, which stock to include in results. Default = 1.

  `$region`

  :   integer, which region to include in results. Default = 1.

- table.opts:

  list of options for AIC/rho table:

  `$fname`

  :   character, filename to save CSV results table (.csv will be
      appended). Default = `'model_comparison'`.

  `$sort`

  :   T/F, sort by AIC? Default = TRUE.

  `$calc.rho`

  :   T/F, calculate Mohn's rho? Retrospective analysis must have been
      run for all modes. Default = TRUE.

  `$calc.aic`

  :   T/F, calculate AIC? Default = TRUE.

  `$print`

  :   T/F, print table to console? Default = TRUE.

  `$save.csv`

  :   T/F, save table as a CSV file? Default = FALSE.

- plot.opts:

  list of options for plots:

  `$out.type`

  :   character, either `'pdf'` or `'png'` (default = `'png'` because I
      am not sure `system('pdftk')` will work across platforms.)

  `$ci`

  :   vector of T/F, length = 1 (applied to all models) or number of
      models

  `$years`

  :   vector, which years to plot? Default = all (model and projection
      years).

  `$which`

  :   vector, which plots to make? Default = all. See details.

  `$relative.to`

  :   character, name of "base" model to plot differences relative to.

  `$alpha`

  :   scalar, (1-alpha)% confidence intervals will be plotted. Default =
      0.05 for 95% CI.

  `$ages.lab`

  :   vector, overwrite model age labels.

  `$kobe.yr`

  :   integer vector (length = 1 (applied to all models or number of
      models, which year to use in Kobe plot (relative status). Default
      = terminal model year(s).

  `$M.age`

  :   integer, which age to use in M time-series plot. Default =
      `max(data$which_F_age)` (max age of F to use for full total F).

  `$return.ggplot`

  :   T/F, return a list of ggplot2 objects for later modification?
      Default = TRUE.

  `$kobe.prob`

  :   T/F, print probabilities for each model in each quadrant of Kobe
      plot? Default = TRUE.

  `$refpt`

  :   "XSPR" or "MSY", which reference point to use. Default = "XSPR".

  `$browse`

  :   Open html document in web browser (if \$out.type = `'png'`)
      (default = TRUE).

## Value

a list with the following components:

- `daic`:

  Vector of delta-AIC by model (if `do.table=T` and
  `table.opts$calc.aic=T`)

- `aic`:

  Vector of AIC by model (if `do.table=T` and `table.opts$calc.aic=T`)

- `rho`:

  Matrix of Mohn's rho by model (if `do.table=T` and
  `table.opts$calc.rho=T`)

- `best`:

  Name of best model (lowest AIC) (if `do.table=T` and
  `table.opts$calc.aic=T`)

- `tab`:

  Results table of AIC and Mohn's rho (if `do.table=T`)

- `g`:

  List of ggplot2 objects for later modification (if `do.plot=T` and
  `plot.opts$return.ggplot=T`)

## Details

`plot.opts$which` specifies which plots to make:

- 1:

  3-panel of SSB (spawning stock biomass), F (fully-selected fishing
  mortality), and Recruitment

- 2:

  CV (coefficient of variation) for SSB, F, and Recruitment

- 3:

  Fleet selectivity (by block, averaged across years)

- 4:

  Index selectivity (by block, averaged across years)

- 5:

  Selectivity tile (fleets + indices, useful for time-varying random
  effects)

- 6:

  M time series (natural mortality, can specify which age with
  plot.opts\$M.age)

- 7:

  M tile (useful for time-varying random effects)

- 8:

  3-panel of F X% SPR, SSB at F_X%SPR, and yield at F_X%SPR

- 9:

  2-panel of relative status (SSB / SSB at F_X%SPR and F / F_X%SPR)

- 10:

  Kobe status (relative SSB vs. relative F)

If `plot.opts$return.ggplot = TRUE`, a list `g` is returned holding the
above ggplot2 objects for later modification. `g[[i]]` holds the plot
corresponding to `i` above, e.g. `g[[2]]` is the CV plot.

## See also

[`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md),
[`read_asap3_fit`](https://timjmiller.github.io/wham/reference/read_asap3_fit.md)`, `[`read_wham_fit`](https://timjmiller.github.io/wham/reference/read_wham_fit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
base <- read_asap3_fit()
m1 <- fit_wham(input1)
m2 <- fit_wham(input2)
mods <- list(base=base, m1=m1, m2=m2)
res <- compare_wham_models(mods)
} # }
```
