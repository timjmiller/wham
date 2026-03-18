# Add TMB sdreport object to WHAM model

Runs [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html) and
adds the object to the fitted (and projected) model list. E.g.,
fit\$sdrep.

## Usage

``` r
do_sdreport(
  model,
  save.sdrep = TRUE,
  TMB.bias.correct = FALSE,
  TMB.jointPrecision = FALSE
)
```

## Arguments

- model:

  a fitted WHAM model object returned by fit_wham or project_wham.

- save.sdrep:

  T/F, save the full
  [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html) object?
  If `FALSE`, only save
  [`summary.sdreport`](https://rdrr.io/pkg/TMB/man/summary.sdreport.html)
  to reduce model object file size. Default = `TRUE`.

- TMB.bias.correct:

  T/F whether to use the bias.correct feature of TMB::sdreport. Default
  = `FALSE`.

- TMB.jointPrecision:

  T/F whether to return the joint precision matrix for the fixed and
  random effects from TMB::sdreport. Default = `FALSE`.

## See also

[`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md),
[`project_wham`](https://timjmiller.github.io/wham/reference/project_wham.md)
