# Specify model and parameter configuration for selectivity

Specify model and parameter configuration for selectivity

## Usage

``` r
set_selectivity(input, selectivity)
```

## Arguments

- input:

  list containing data, parameters, map, and random elements (output
  from
  [`prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md))

- selectivity:

  (optional) list specifying options for selectivity blocks, models,
  initial parameter values, parameter fixing/mapping, and random effects
  (see details)

  `set_selectivity` specifies options for selectivity and allows you to
  overwrite existing options in the `input` list or as specified in the
  ASAP data file. If `selectivity = NULL`, selectivity options from
  `input` are used.

  [`prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md)`(..., selectivity=selectivity)`
  calls `set_selectivity(..., selectivity=selectivity)`. If you already
  have created `input` with `prepare_wham_input`, you can also use
  `set_selectivity(input, selectivity=selectivity)` to modify the
  selectivity specification.

  `selectivity` is a list with the following entries:

  \$model

  :   Selectivity model for each block. Vector with length = number of
      selectivity blocks. Each entry must be one of: "age-specific",
      "logistic", "double-logistic", or "decreasing-logistic".

  \$re

  :   Time-varying (random effects) for each block. Vector with length =
      number of selectivity blocks. If `NULL`, selectivity parameters in
      all blocks are constant over time and uncorrelated. Each entry of
      `selectivity$re` must be one of the following options, where the
      selectivity parameters are:

      "none"

      :   (default) are constant and uncorrelated

      "iid"

      :   vary by year and age/par, but uncorrelated

      "ar1"

      :   correlated by age/par (AR1), but not year

      "ar1_y"

      :   correlated by year (AR1), but not age/par

      "2dar1"

      :   correlated by year and age/par (2D AR1)

  \$initial_pars

  :   Initial parameter values for each block. List of length = number
      of selectivity blocks. Each entry must be a vector of length \#
      parameters in the block, i.e. `c(2,0.2)` for logistic (a50 and
      1/slope) or `c(0.5,0.5,0.5,1,1,0.5)` for age-specific parameters
      when there are 6 ages. Default is to set at middle of parameter
      range. This is 0.5 for age-specific and n.ages/2 or logistic,
      double-logistic, and decreasing-logistic.

  \$fix_pars

  :   Alternative to `$map_pars` for specifying which selectivity
      parameters (only fixed effects) to fix at initial values. List of
      length = number of selectivity blocks. E.g. model with 3
      age-specific blocks and 6 ages, `list(4:5, 4, 2:4))` will fix ages
      4 and 5 in block 1, age 4 in block 2, and ages 2, 3, and 4 in
      block 3. Use NULL to not fix any parameters for a block, e.g.
      list(NULL, 4, 2) does not fix any pars in block 1.

  \$par_min

  :   The lower bound for selectivity parameters and is used to populate
      `data$selpars_lower`. List of length = number of selectivity
      blocks, where each item is a vector of length = number of
      selectivity parameters (age-specific: n.ages, logistic: 2,
      double-logistic: 4).

  \$par_max

  :   The upper bound for selectivity parameters and is used to populate
      `data$selpars_upper`. List of length = number of selectivity
      blocks, where each item is a vector of length = number of
      selectivity parameters (age-specific: n.ages, logistic: 2,
      double-logistic: 4).

  \$map_pars

  :   Alternative to `$fix_pars` for specifying how to fix selectivity
      parameters (only fixed effects), corresponds to
      `map$logit_selpars`. List of length = number of selectivity
      blocks, where each item is a vector of length = number of
      selectivity parameters (age-specific: n.ages, logistic: 2,
      double-logistic: 4). Use `NA` to fix a parameter and integers to
      estimate. Use the same integer for multiple ages or fleets/indices
      to estimate a shared parameter. E.g. for a model with 3
      age-specific blocks (1 fleet, 2 indices) and 6 ages,
      `$map_pars = list(c(1,2,3,NA,NA,4), c(5,6,7,NA,8,8), c(1,2,3,NA,NA,4))`
      will estimate ages 1-3 and 6 in block 1 (fleet), ages 1-3 and 4-5
      (shared) in block 2 (index 1), and then set the index 2 (block 3)
      selectivity equal to the fleet.

  \$sigma_vals

  :   Initial standard deviation values to use for the random effect
      deviations. Must be a vector with length = number of blocks. Use
      natural (not log) scale, must be positive. `par$sel_repars[,1]`
      will be estimated on log-scale. Not used if `re = 'none'` for all
      blocks.

  \$map_sigma

  :   Specify which SD parameters to fix for the random effect
      deviations. Must be a vector with length = number of blocks. Use
      `NA` to fix a parameter and integers to estimate. Use the same
      integer for multiple blocks to estimate a shared SD parameter. Not
      used if `re = 'none'` for all blocks.

  \$cor_vals

  :   Initial correlation values to use for the random effect
      deviations. Must be a n_selblocks x 2 integer matrix. Columns
      correspond to correlation by age and year, respectively. If
      `re = 'ar1'` or `re = 'ar1_y'` only the corresponding values are
      used. Values must be between -1 and 1, but parameters are
      estimated on a logit transformed scale internally. Not used if
      `re = 'none'` or `re = 'iid'` for all blocks.

  \$map_cor

  :   Specify which correlation parameters to fix for the random effect
      deviations. Must be a n_selblocks x 2 matrix. Columns correspond
      to correlation by age and year,respectively. Parameters can be
      shared by setting corresponding values of `map_cor` to the same
      integer. Use `NA` to fix a parameter. If `re = 'ar1'` or
      `re = 'ar1_y'`, only the column for the corresponding correlation
      are used. Not used if `re = 'none'` or `re = 'iid'` for all
      blocks.

  \$n_selblocks

  :   How many selectivity blocks. Optional. If unspecified and no asap3
      object, then this is set to the number of fleets + indices. If
      specified, ensure other components of `selectivity` are
      consistent.

## Value

a named list with same elements as the input provided with selectivity
options modified.

## See also

[`prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md)

## Examples

``` r
if (FALSE) { # \dontrun{
wham.dir <- find.package("wham")
path_to_examples <- system.file("extdata", package="wham")
asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
input <- prepare_wham_input(asap3, NAA_re = list(sigma = "rec"))
sel <- list(model=rep("logistic",input$data$n_selblocks),
   initial_pars=rep(list(c(3,3)),input$data$n_selblocks),
   fix_pars=rep(list(NULL),input$data$n_selblocks))
input <- set_selectivity(input, selectivity = sel) #logistic selectivity for all selectivity blocks
} # }
```
