# Read WHAM fit

Gets output from a fit WHAM model for plotting with other models.
Internal function, called within
[`compare_wham_models`](https://timjmiller.github.io/wham/reference/compare_wham_models.md).

## Usage

``` r
read_wham_fit(mod, alphaCI = 0.05)
```

## Arguments

- mod:

  output from
  [`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md)

- alphaCI:

  (1-alpha)% confidence intervals will be calculated. Default = 0.05 for
  95% CI.

## Value

a named list with the following elements:

- `$years`:

  numeric vector, model years only, e.g. `1972:2020`

- `$years_full`:

  numeric vector, model + proj years, e.g. `1972:2022`

- `$selAA`:

  list of length(n_selblocks), first the fleet blocks then indices, i.e.
  if 4 fleet blocks and 3 indices, `selAA[[5]]` is for index 1. Each
  element is a matrix, years (rows) x ages (cols), selectivity at age

- `$selblock_pointer_fleets`:

  matrix, years x fleets, indices of selAA used by each fleet in each
  year

- `$selblock_pointer_indices`:

  matrix, years x indices, indices of selAA used by each index in each
  year

- `$MAA`:

  array, stocks x regions x years x ages, natural mortality

- `$log_SSB`:

  matrix, years x 2, log-scale spawning stock biomass. 1st col = MLE,
  2nd col = SE.

- `$log_F`:

  matrix, years x 2, log-scale fully-selected F. 1st col = MLE, 2nd col
  = SE.

- `$log_NAA_rep`:

  array, stocks x regions x years x ages, numbers at age

- `$NAA_CV`:

  array, stocks x regions x years x ages, CV of numbers at age

- `$percentSPR`:

  scalar, X% SPR used to calculate reference points, default = 40

- `$log_Y_FXSPR`:

  matrix, years x 2, log-scale yield at FXSPR. 1st col = MLE, 2nd col =
  SE.

- `$log_FXSPR`:

  matrix, years x 2, log-scale FXSPR. 1st col = MLE, 2nd col = SE.

- `$log_SSB_FXSPR`:

  matrix, years x 2, log-scale SSB at FXSPR. 1st col = MLE, 2nd col =
  SE.

- `$log_rel_ssb_F_cov`:

  list, length n_years, each element is a 2x2 covariance matrix with
  SSB/SSB_FXSPR first and F/F_FXSPR second

## See also

[`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md),
[`read_asap3_fit`](https://timjmiller.github.io/wham/reference/read_asap3_fit.md),
[`compare_wham_models`](https://timjmiller.github.io/wham/reference/compare_wham_models.md)
