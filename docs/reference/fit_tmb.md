# Fit TMB model using nlminb

Runs optimization on the TMB model using
[`stats::nlminb`](https://rdrr.io/r/stats/nlminb.html). If specified,
takes additional Newton steps and calculates standard deviations.
Internal function called by
[`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md).

## Usage

``` r
fit_tmb(
  model,
  n.newton = 3,
  do.sdrep = TRUE,
  do.check = FALSE,
  save.sdrep = FALSE,
  use.optim = FALSE,
  opt.control = NULL
)
```

## Arguments

- model:

  Output from
  [`TMB::MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html).

- n.newton:

  Integer, number of additional Newton steps after optimization. Default
  = `3`.

- do.sdrep:

  T/F, calculate standard deviations of model parameters? See
  [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html). Default
  = `TRUE`.

- do.check:

  T/F, check if model parameters are identifiable? Runs internal
  [`check_estimability`](https://timjmiller.github.io/wham/reference/check_estimability.md),
  originally provided by the
  [TMBhelper](https://github.com/kaskr/TMB_contrib_R/tree/master/TMBhelper)
  package. Default = `TRUE`.

- save.sdrep:

  T/F, save the full
  [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html) object?
  If `FALSE`, only save
  [`summary.sdreport)`](https://rdrr.io/pkg/TMB/man/summary.sdreport.html)
  to reduce model object file size. Default = `FALSE`.

- use.optim:

  T/F, use [`stats::optim`](https://rdrr.io/r/stats/optim.html) instead
  of [`stats::nlminb`](https://rdrr.io/r/stats/nlminb.html)? Default =
  `FALSE`.

- opt.control:

  list of control parameters to pass to optimizer. For nlminb default =
  list(iter.max = 1000, eval.max = 1000). For optim default =
  list(maxit=1000).

## Value

`model`, appends the following:

- `model$opt`:

  Output from [`stats::nlminb`](https://rdrr.io/r/stats/nlminb.html)

- `model$date`:

  System date

- `model$dir`:

  Current working directory

- `model$rep`:

  model\$report(model\$env\$last.par.best)

- `model$TMB_version`:

  Version of TMB installed

- `model$parList`:

  List of parameters,
  `model$env$parList(x = model$opt$par, par = model$env$last.par.best)`

- `model$final_gradient`:

  Final gradient, `model$gr(model$opt$par)`

- `model$sdrep`:

  Estimated standard deviations for model parameters,
  [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html) or
  [`summary.sdreport)`](https://rdrr.io/pkg/TMB/man/summary.sdreport.html)

## See also

[`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md),
[`retro`](https://timjmiller.github.io/wham/reference/retro.md),
`check_estimability`
