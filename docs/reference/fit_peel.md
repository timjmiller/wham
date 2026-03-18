# Fit model peeling off *i* years of data

Internal function called by
[`retro`](https://timjmiller.github.io/wham/reference/retro.md) for *i*
in 1–`n.peels`. Fits the model peeling off *i* years of data (calls
[`fit_tmb`](https://timjmiller.github.io/wham/reference/fit_tmb.md)).

## Usage

``` r
fit_peel(
  peel,
  input,
  do.sdrep = FALSE,
  n.newton = 3,
  MakeADFun.silent = FALSE,
  retro.silent = FALSE,
  save.input = FALSE
)
```

## Arguments

- peel:

  Integer, number of years of data to remove before model fitting.

- input:

  input with same structure as that provided by
  [`prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md).
  May want to use input\$par = model\$parList to start at MLEs.

- do.sdrep:

  T/F, calculate standard deviations of model parameters? Default =
  `FALSE`.

- n.newton:

  integer, number of additional Newton steps after optimizafit_tmbtion
  for each peel. Default = `3`.

- MakeADFun.silent:

  T/F, Passed to silent argument of
  [`TMB::MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html).
  Default = `FALSE`.

- retro.silent:

  T/F, Passed to argument of internal fit_peel function. Determines
  whether peel number is printed to screen. Default = `FALSE`.

- save.input:

  T/F, should modified input list be saved? Necessary to project from a
  peel but increases model object size. Default = `FALSE`.

## Value

`out`, output of
[`fit_tmb`](https://timjmiller.github.io/wham/reference/fit_tmb.md) for
peel *i*

## See also

[`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md),
[`retro`](https://timjmiller.github.io/wham/reference/retro.md),
[`fit_tmb`](https://timjmiller.github.io/wham/reference/fit_tmb.md)
