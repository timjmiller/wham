# Add reporting of biological reference points to WHAM model

Changes internal flags to do the extra calculations and reporting for
reference points.

## Usage

``` r
do_reference_points(model, do.sdrep = FALSE, save.sdrep = TRUE)
```

## Arguments

- model:

  a fitted WHAM model object returned by fit_wham or project_wham.

- do.sdrep:

  T/F, calculate standard deviations of model parameters? See
  [`sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html). Default =
  `FALSE`.

- save.sdrep:

  T/F, save the full
  [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html) object?
  If `FALSE`, only save
  [`summary.sdreport`](https://rdrr.io/pkg/TMB/man/summary.sdreport.html)
  to reduce model object file size. Default = `TRUE`.

## See also

[`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md),
[`project_wham`](https://timjmiller.github.io/wham/reference/project_wham.md)
