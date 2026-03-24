# Read ASAP3 fit

Gets output from a fit ASAP3 model for plotting with WHAM models.

## Usage

``` r
read_asap3_fit(wd, asap.name, pSPR = 40)
```

## Arguments

- wd:

  character, directory where ASAP3 output files are located (ex:
  'C:/MY/file/directories/model/'). 5 files are needed: `.rdat`, `.dat`,
  `.std`, `.cor`, and `.par`.

- asap.name:

  character, base name of original .dat file (i.e. without the .dat
  extension)

- pSPR:

  scalar, user-specified percent SPR to use for reference points,
  expressed as 100\*SSBPR(Fspr)/SSBPR(F=0). Default = 40.

## Value

a named list with the following elements:

- `$years`:

  numeric vector, model years only, e.g. `1972:2020`

- `$years_full`:

  numeric vector, model + proj years, e.g. `1972:2022`. For ASAP this
  will be the same as `$years`.

- `$selAA`:

  list of length(n_selblocks), first the fleet blocks then indices, i.e.
  if 4 fleet blocks and 3 indices, `selAA[[5]]` is for index 1. Each
  element is a matrix, years (rows) x ages (cols), selectivity at age

- `$selblock_pointer_fleets`:

  matrix, n_years x n_fleets, indices of selAA used by each fleet in
  each year

- `$selblock_pointer_indices`:

  matrix, n_years x n_indices, indices of selAA used by each index in
  each year

- `$MAA`:

  matrix, n_years x n_ages, natural mortality

- `$log_SSB`:

  matrix, n_years x 2, log-scale spawning stock biomass. 1st col = MLE,
  2nd col = SE (from .std file in ADMB).

- `$log_F`:

  matrix, n_years x 2, log-scale fully-selected F. 1st col = MLE, 2nd
  col = SE (from .std file in ADMB).

- `$log_NAA`:

  matrix, n_years x n_ages, numbers at age

- `$NAA_CV`:

  matrix, n_years x n_ages, CV of numbers at age

- `$percentSPR`:

  scalar, X% SPR used to calculate reference points, default = 40

- `$log_Y_FXSPR`:

  matrix, n_years x 2, log-scale yield at FXSPR. 1st col = MLE, 2nd col
  = SE.

- `$log_FXSPR`:

  matrix, n_years x 2, log-scale FXSPR. 1st col = MLE, 2nd col = SE.

- `$log_SSB_FXSPR`:

  matrix, n_years x 2, log-scale SSB at FXSPR, i.e. annual numerator of
  SPR(FXSPR) \* Recruits. 1st col = MLE, 2nd col = SE.

- `$log_rel_ssb_F_cov`:

  list, length n_years, each element is a 2x2 covariance matrix with SSB
  / SSB_FXSPR first and F / F_FXSPR second

## See also

[`compare_wham_models`](https://timjmiller.github.io/wham/reference/compare_wham_models.md),
[`read_wham_fit`](https://timjmiller.github.io/wham/reference/read_wham_fit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
base <- read_asap3_fit(wd=file.path(getwd(),'asap_results'), asap.name='BASE_5C.DAT', pSPR=40)
m1 <- fit_wham(input1)
m2 <- fit_wham(input2)
mods <- list(base=base, m1=m1, m2=m2)
res <- compare_wham_models(mods)
} # }
```
