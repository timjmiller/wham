# Prepare input data and parameters to project an already fit WHAM model

`prepare_projection` is an internal function called by
[`project_wham`](https://timjmiller.github.io/wham/reference/project_wham.md),
which in turn is called by
[`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md) if
`do.proj = TRUE`.

## Usage

``` r
prepare_projection(model, proj.opts, check.version = FALSE)
```

## Arguments

- model:

  a previously fit wham model

- proj.opts:

  a named list with the following components:

  - `$n.yrs` (integer), number of years to project/forecast. Default =
    `3`.

  - `$use.last.F` (T/F), use terminal year F for projections. Default =
    `TRUE`.

  - `$use.avg.F` (T/F), use average of F over certain years for
    projections. Default = `FALSE`. Years to average over determined by
    \$avg.yrs defined below.

  - `$use.FXSPR` (T/F), calculate and use F at X% SPR for projections.
    Default = `FALSE`.

  - `$use.FMSY` (T/F), calculate and use FMSY for projections. Default =
    `FALSE`.

  - `$proj.F` (vector), user-specified fishing mortality for
    projections. Length must equal `n.yrs`.

  - `$proj.catch` (vector), user-specified aggregate catch for
    projections. Length must equal `n.yrs`.

  - `$avg.yrs` (vector), specify which years to use to average
    population attributes (MAA,FAA,WAA,maturity,movement) in projection
    years. Any BRPs calculated in projection years will also use these.
    Default = last 5 years, `tail(model$years, 5)`.

  - `$cont.ecov` (T/F), continue ecov process (e.g. random walk or AR1)
    for projections. Default = `TRUE`.

  - `$use.last.ecov` (T/F), use terminal year ecov for projections.

  - `$avg.ecov.yrs` (vector), specify which years to average the
    environmental covariate(s) over for projections.

  - `$proj.ecov` (matrix), user-specified environmental covariate(s) for
    projections. `n.yrs x n.ecov`.

  - `$cont.M.re` (T/F), continue M random effects (i.e. AR1_y or 2D AR1)
    for projections. Default = `FALSE`. If `FALSE`, M will be averaged
    over `$avg.yrs.M` (which defaults to last 5 model years).

  - `$cont.move.re` (T/F), continue any movement random effects for
    projections. Default = `FALSE`. If `FALSE`, movement parameters will
    be averaged over `$avg.yrs.move` (which defaults to last 5 model
    years).

  - `$cont.L.re` (T/F), continue any L ("extra mortality rate") random
    effects for projections. Default = `FALSE`. If `FALSE`, L parameters
    will be averaged over `$avg.yrs.L` (which defaults to last 5 model
    years).

  - `$avg.rec.yrs` (vector), specify which years to calculate the CDF of
    recruitment for use in projections. Default = all model years. Only
    used when recruitment is estimated as fixed effects (SCAA).

  - `$percentFXSPR` (scalar), percent of F_XSPR to use for projections,
    only used if \$use.FXSPR = TRUE. For example, to project with F =
    75% F_40%SPR, `proj.opts$percentFXSPR = 75`. Default = 100.

  - `$percentFMSY` (scalar), percent of F_MSY to use for projections,
    only used if \$use.FMSY = TRUE and a stock-recruit relationship is
    assumed. Default = 100.

  - `$proj_F_opt` (vector), integers specifying how to configure each
    year of the projection: 1: use terminal F, 2: use average F, 3: use
    F at X% SPR, 4: use specified F, 5: use specified catch, 6: use
    Fmsy. Overrides any of the above specifications.

  - `$proj_Fcatch` (vector or matrix), catch or F values to use each
    projection year: values are not used when using Fmsy, FXSPR,
    terminal F or average F. Overrides any of the above specifications
    of proj.F or proj.catch. if vector, total catch or F is supplied
    else matrix columns should be fleets for fleet-specific F to be
    found/used (`n.yrs` x 1 or n_fleets).

  - `$proj_mature` (array), user-supplied maturity values for the
    projection years with dimensions (n_stocks x `n.yrs` x n_ages).

  - `$proj_waa` (3-d array), user-supplied waa values for the projection
    years with first and third dimensions equal to that of
    `model$input$data$waa` (waa source x `n.yrs` x n_ages).

  - `$proj_R_opt` (integer), 1: continue any RE processes for
    recruitment, 2: make projected recruitment consistent with average
    recruitment in SPR reference points and cancel any bias correction
    for NAA in projection years. 3: average recruitment deviations over
    \$avg.yrs.R (if \$sigma = "rec") 4: no recruitment deviations (if
    \$sigma = "rec").

  - `$proj_NAA_opt` (integer), 1: continue any RE processes for NAA, 2:
    average NAA deviations over \$avg.yrs.NAA. 3: no NAA deviations.

  - `$proj_NAA_init` (scalar), the default starting value for all NAA
    random effects in projection years is exp(10), which may not be
    large enough for some catch specification. Use this to change the
    default if a call to project_wham suggests it.

  - `$proj_F_init` which F to initialize internal newton search for
    annual projected F for a given user-specifed catch. Default is 0.1

  - `$avg.yrs.sel` list (length = n_fleets), years to average
    selectivity or FAA for each fleet for projection years. Any BRPs
    calculated in projection years will also use this. Default = last 5
    years, `tail(model$years, 5)`.

  - `$avg.yrs.waacatch` list (length = n_fleets), years to average
    weight at age for each fleet for projection years (if \$proj_waa is
    NULL). Any BRPs calculated in projection years will also use this.
    Default = last 5 years, `tail(model$years, 5)`.

  - `$avg.yrs.waassb` list (length = n_stocks), years to average weight
    at age for each stock SSB for projection years (if \$proj_waa is
    NULL). Any BRPs calculated in projection years will also use this.
    Default = last 5 years, `tail(model$years, 5)`.

  - `$avg.yrs.mature` list (length = n_stocks), years to average
    maturity at age for each stock for projection years (if
    \$proj_mature is NULL). Any BRPs calculated in projection years will
    also use this. Default = last 5 years, `tail(model$years, 5)`.

  - `$avg.yrs.L` list (length = n_regions), years to average extra
    mortality at age for each region for projection years (if
    \$cond.L.re = FALSE). Any BRPs calculated in projection years will
    also use this. Default = last 5 years, `tail(model$years, 5)`.

  - `$avg.yrs.M` list (length = n_stocks, each is a list with length =
    n_regions), years to average natural mortality at age for each stock
    and region for projection years (if \$cont.M.re = FALSE). Any BRPs
    calculated in projection years will also use this. Default = last 5
    years, `tail(model$years, 5)`.

  - `$avg.yrs.move` list (length = n_stocks, each is a list with length
    = n_regions), years to average movement rates at age and season for
    each stock and region (at beginning of interval) for projection
    years (if \$cont.move.re = FALSE). Any BRPs calculated in projection
    years will also use this. Default = last 5 years,
    `tail(model$years, 5)`.

  - `$avg.yrs.R` list (length = n_stocks), years to average recruitment
    deviations for each stock and region for projection years (if
    \$proj_R_opt = 3). Any BRPs calculated in projection years will also
    use this. Default = last 5 years, `tail(model$years, 5)`.

  - `$avg.yrs.NAA` list (length = n_stocks, each is a list with length =
    n_regions), years to average NAA deviations for each stock and
    region for projection years (if \$proj_NAA_opt = 2). Any BRPs
    calculated in projection years will also use this. Default = last 5
    years, `tail(model$years, 5)`.

- check.version:

  T/F check whether version WHAM and TMB for fitted model match that of
  the version of WHAM using for projections. Default = `TRUE`.

## Value

same as
[`prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md),
a list ready for
[`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md):

- `data`:

  Named list of data, passed to
  [`TMB::MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html)

- `par`:

  Named list of parameters, passed to
  [`TMB::MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html)

- `map`:

  Named list of factors that determine which parameters are estimated,
  passed to
  [`TMB::MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html)

- `random`:

  Character vector of parameters to treat as random effects, passed to
  [`TMB::MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html)

- `years`:

  Numeric vector of representing (non-projection) model years of WHAM
  model

- `years_full`:

  Numeric vector of representing all model and projection years of WHAM
  model

- `ages.lab`:

  Character vector of age labels, ending with plus-group

- `model_name`:

  Character, name of stock/model (specified in call to
  `prepare_wham_input`)

## See also

[`prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md),
[`project_wham`](https://timjmiller.github.io/wham/reference/project_wham.md)
