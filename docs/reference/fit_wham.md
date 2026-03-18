# Fit WHAM model

Fits the compiled WHAM model using
[`TMB::MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html) and
[`stats::nlminb`](https://rdrr.io/r/stats/nlminb.html). Runs
retrospective analysis if specified.

## Usage

``` r
fit_wham(
  input,
  n.newton = 3,
  do.sdrep = TRUE,
  do.retro = TRUE,
  n.peels = 7,
  do.osa = TRUE,
  osa.opts = list(method = "oneStepGaussianOffMode", parallel = TRUE),
  do.post.samp = TRUE,
  model = NULL,
  do.check = FALSE,
  MakeADFun.silent = FALSE,
  retro.silent = FALSE,
  do.proj = FALSE,
  proj.opts = list(n.yrs = 3, use.last.F = TRUE, use.avg.F = FALSE, use.FXSPR = FALSE,
    proj.F = NULL, proj.catch = NULL, avg.yrs = NULL, cont.ecov = TRUE, use.last.ecov =
    FALSE, avg.ecov.yrs = NULL, proj.ecov = NULL, cont.Mre = NULL, avg.rec.yrs = NULL,
    percentFXSPR = 100),
  do.fit = TRUE,
  save.sdrep = TRUE,
  do.brps = TRUE,
  fit.tmb.control = NULL,
  TMB.bias.correct = FALSE,
  TMB.jointPrecision = FALSE
)
```

## Arguments

- input:

  Named list with components:

  `$data`

  :   Data to fit the assessment model to.

  `$par`

  :   Parameters, a list of all parameter objects required by the user
      template (both random and fixed effects). See
      [`MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html).

  `$map`

  :   Map, a mechanism for collecting and fixing parameters. See
      [`MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html).

  `$random`

  :   Character vector defining the parameters to treat as random
      effect. See
      [`MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html).

  `$years`

  :   Numeric vector of the years which the model spans. Not important
      for model fitting, but useful for plotting.

  `$model_name`

  :   Character, name of the model, e.g. `"Yellowtail flounder"`

  `$ages.lab`

  :   Character vector of the age labels, e.g. `c("1","2","3","4+").`

- n.newton:

  integer, number of additional Newton steps after optimization. Passed
  to
  [`fit_tmb`](https://timjmiller.github.io/wham/reference/fit_tmb.md).
  Default = `3`.

- do.sdrep:

  T/F, calculate standard deviations of model parameters? See
  [`sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html). Default =
  `TRUE`.

- do.retro:

  T/F, do retrospective analysis? Default = `TRUE`.

- n.peels:

  integer, number of peels to use in retrospective analysis. Default =
  `7`.

- do.osa:

  T/F, calculate one-step-ahead (OSA) residuals? Default = `TRUE`. See
  details. Returned as `mod$osa$residual`.

- osa.opts:

  list of 2 options (method, parallel) for calculating OSA residuals,
  passed to
  [`TMB::oneStepPredict`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html).
  Default:
  `osa.opts = list(method="oneStepGaussianOffMode", parallel=TRUE)`. See
  [`make_osa_residuals`](https://timjmiller.github.io/wham/reference/make_osa_residuals.md).

- do.post.samp:

  T/F, obtain sample from posterior of random effects? Default = `TRUE`.
  NOT YET IMPLEMENTED.

- model:

  (optional), a previously fit wham model.

- do.check:

  T/F, check if model parameters are identifiable? Passed to
  [`fit_tmb`](https://timjmiller.github.io/wham/reference/fit_tmb.md).
  Runs internal function `check_estimability`, originally provided by
  https://github.com/kaskr/TMB_contrib_R/TMBhelper. Default = `TRUE`.

- MakeADFun.silent:

  T/F, Passed to silent argument of
  [`TMB::MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html).
  Default = `FALSE`.

- retro.silent:

  T/F, Passed to argument of internal retro function. Determines whether
  peel number is printed to screen. Default = `FALSE`.

- do.proj:

  T/F, do projections? Default = `FALSE`. If true, runs
  [`project_wham`](https://timjmiller.github.io/wham/reference/project_wham.md).

- proj.opts:

  a named list with the following components:

  - `$n.yrs` (integer), number of years to project/forecast. Default =
    `3`.

  - `$use.last.F` (T/F), use terminal year F for projections. Default =
    `TRUE`.

  - `$use.FXSPR` (T/F), calculate and use F at X

  - `$use.FMSY` (T/F), calculate and use FMSY for projections.

  - `$proj.F` (vector), user-specified fishing mortality for
    projections. Length must equal `n.yrs`.

  - `$proj.catch` (vector), user-specified aggregate catch for
    projections. Length must equal `n.yrs`.

  - `$avg.yrs` (vector), specify which years to average over for
    calculating reference points. Default = last 5 model years,
    `tail(model$years, 5)`.

  - `$cont.ecov` (T/F), continue ecov process (e.g. random walk or AR1)
    for projections. Default = `TRUE`.

  - `$use.last.ecov` (T/F), use terminal year ecov for projections.

  - `$avg.ecov.yrs` (vector), specify which years to average over the
    environmental covariate(s) for projections.

  - `$proj.ecov` (matrix), user-specified environmental covariate(s) for
    projections. `n.yrs x n_Ecov`.

  - `$cont.Mre` (T/F), continue M random effects (i.e. AR1_y or 2D AR1)
    for projections. Default = `TRUE`. If `FALSE`, M will be averaged
    over `$avg.yrs` (which defaults to last 5 model years).

  - `$avg.rec.yrs` (vector), specify which years to calculate the CDF of
    recruitment for use in projections. Default = all model years.

  - `$percentFXSPR` (scalar), percent of F_XSPR to use for calculating
    catch in projections, only used if proj.opts\$use.FXSPR = TRUE. For
    example, GOM cod uses F = 75

  - `$percentFMSY` (scalar), percent of F_MSY to use for calculating
    catch in projections, only used if \$use.FMSY = TRUE.

- do.fit:

  T/F, fit the model using `fit_tmb`. Default = `TRUE`.

- save.sdrep:

  T/F, save the full
  [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html) object?
  If `FALSE`, only save
  [`summary.sdreport`](https://rdrr.io/pkg/TMB/man/summary.sdreport.html)
  to reduce model object file size. Default = `TRUE`.

- do.brps:

  T/F, calculate and report biological reference points. Default =
  `TRUE`.

- fit.tmb.control:

  list of optimizer controlling attributes passed to
  [`fit_tmb`](https://timjmiller.github.io/wham/reference/fit_tmb.md).
  Default is
  `list(use.optim = FALSE, opt.control = list(iter.max = 1000, eval.max = 1000))`,
  so stats::nlminb is used to opitmize.

- TMB.bias.correct:

  T/F whether to use the bias.correct feature of TMB::sdreport. Default
  = `FALSE`.

- TMB.jointPrecision:

  T/F whether TMB::sdreport should return the joint precision matrix for
  the fixed and random effects. Default = `FALSE`.

## Value

a fit TMB model with additional output if specified:

- `$rep`:

  List of derived quantity estimates (see examples)

- `$sdrep`:

  Parameter estimates (and standard errors if `do.sdrep=TRUE`)

- `$peels`:

  Retrospective analysis (if `do.retro=TRUE`)

- `$osa`:

  One-step-ahead residuals (if `do.osa=TRUE`)

## Details

Standard residuals are not appropriate for models with random effects.
Instead, one-step-ahead (OSA) residuals can be used for evaluating model
goodness-of-fit ([Thygeson et al.
(2017)](https://link.springer.com/article/10.1007/s10651-017-0372-4),
implemented in
[`TMB::oneStepPredict`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)).
Additional OSA residual options are passed to
[`TMB::oneStepPredict`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
in a list `osa.opts`. For example, to use the (much faster, ~1 sec
instead of 2 min) full Gaussian approximation instead of the (default)
generic method, you can use `osa.opts=list(method="fullGaussian")`.

## See also

[`fit_tmb`](https://timjmiller.github.io/wham/reference/fit_tmb.md),
[`retro`](https://timjmiller.github.io/wham/reference/retro.md),
[`TMB::oneStepPredict`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html),
[`project_wham`](https://timjmiller.github.io/wham/reference/project_wham.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data("input4_SNEMAYT") # load SNEMA yellowtail flounder data and parameter settings
mod = fit_wham(input4_SNEMAYT) # using default values
mod = fit_wham(input4_SNEMAYT, do.retro=FALSE, osa.opts=list(method="oneStepGeneric")) # slower OSA method. 

names(mod$rep) # list of derived quantities
mod$rep$SSB # get SSB estimates (weight, not numbers)
m1$rep$NAA[,1] # get recruitment estimates (numbers, first column of numbers-at-age matrix)
m1$rep$F[,1] # get F estimates for fleet 1
} # }
```
