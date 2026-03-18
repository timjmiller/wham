# Fit retrospective peels and add them to the fitted model object

Function to add fitted retrospective peels to fitted model returned by
[`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md).
Calls [`retro`](https://timjmiller.github.io/wham/reference/retro.md).
to fit the model peeling off `1, 2, ..., n.peels` years of data. This is
just a wrapper for retro that instead of returning just the list of
peels, returns the fitted model with the peels.

## Usage

``` r
do_retro_peels(
  model,
  n.peels = 7,
  ran = NULL,
  use.mle = TRUE,
  do.sdrep = FALSE,
  n.newton = 0,
  MakeADFun.silent = FALSE,
  retro.silent = FALSE,
  save.input = FALSE,
  do.brps = FALSE,
  check.version = TRUE,
  save.sdrep = FALSE
)
```

## Arguments

- model:

  Optimized TMB model, output from
  [`fit_tmb`](https://timjmiller.github.io/wham/reference/fit_tmb.md).

- n.peels:

  Integer, number of peels to use in retrospective analysis. Default =
  `7`.

- ran:

  Character, specifies which parameters to treat as random effects.
  Default = `"model$input$random"`.

- use.mle:

  T/F, use MLEs from full model fit as initial values for each peel? If
  not, the initial values from full model input are used. Default =
  `TRUE`.

- do.sdrep:

  T/F, calculate standard deviations of model parameters for each peel?
  Default = `FALSE`.

- n.newton:

  integer, number of additional Newton steps after optimization for each
  peel. Default = `0`.

- MakeADFun.silent:

  T/F, Passed to silent argument of
  [`TMB::MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html).
  Default = `FALSE`.

- retro.silent:

  T/F, Passed to argument of internal fit_peel function. Determines
  whether peel number is printed to screen. Default = `FALSE`.

- save.input:

  T/F, should modified input list be saved for every peel? Necessary to
  project from a peel but increases model object size. Default =
  `FALSE`.

- do.brps:

  T/F, calculate and report biological reference points

- check.version:

  T/F, whether to verify the wham package commit and version for the
  fitted model are the same as the currently used package.

- save.sdrep:

  T/F, save the full
  [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html) object?
  If `FALSE`, only save
  [`summary.sdreport`](https://rdrr.io/pkg/TMB/man/summary.sdreport.html)
  to reduce model object file size. Default = `FALSE`.

## Value

`peels`, a list of length `n.peels`, where entry *i* is a model fit by
peeling off *i* years of data.

## See also

[`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md),
[`retro`](https://timjmiller.github.io/wham/reference/retro.md),
[`fit_peel`](https://timjmiller.github.io/wham/reference/fit_peel.md)
