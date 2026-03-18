# Calculate one-step-ahead residuals

Standard residuals are not appropriate for models with random effects.
Instead, one-step-ahead (OSA) residuals can be used for evaluating model
goodness-of-fit ([Thygeson et al.
(2017)](https://doi.org/10.1007/s10651-017-0372-4), implemented in
[`TMB::oneStepPredict`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)).
OSA residual options are passed to
[`TMB::oneStepPredict`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
in a list `osa.opts`. Current options are method: oneStepGaussianOffMode
(default), oneStepGaussian, or oneStepGeneric, and parallel: TRUE/FALSE.
See
[`TMB::oneStepPredict`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
for further details. It is not recommended to run this function (or
[`TMB::oneStepPredict`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html))
with any random effects and mvtweedie age composition likelihoods due to
extensive computational demand. An error will be thrown in such cases.
See [Trijoulet et al.
(2023)](https://doi.org/10.1016/j.fishres.2022.106487) for OSA methods
for age composition OSA residuals.

## Usage

``` r
make_osa_residuals(
  model,
  osa.opts = list(method = "oneStepGaussianOffMode", parallel = TRUE),
  sdrep_required = TRUE
)
```

## Arguments

- model:

  A fit WHAM model, output from
  [`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md).

## Value

the same fit TMB model with additional elements for osa residuals:

- `$OSA.Ecov`:

  data.frame returned by
  [`TMB::oneStepPredict`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  for environmental observations, if applicable.

- `$OSA.agregate`:

  data.frame returned by
  [`TMB::oneStepPredict`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  for aggregate catch and index observations conditional on any
  environmental observations, if applicable.

- `$OSA.agecomp`:

  data.frame returned by
  [`TMB::oneStepPredict`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  for age composition observations conditional on any aggregate catch or
  index, or environmental observations, if applicable.

- `$osa`:

  One-step-ahead residuals (if `do.osa=TRUE`)

## See also

[`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data("input4_SNEMAYT") # load SNEMA yellowtail flounder data and parameter settings
mod <- fit_wham(input4_SNEMAYT, do.osa =FALSE, do.retro =FALSE)
mod <- make_osa_residuals(mod) # calculate Mohn's rho
plot_wham_output(mod)
} # }
```
