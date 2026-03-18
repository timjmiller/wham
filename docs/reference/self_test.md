# Perform self-test simulation and estimation of a fitted WHAM model

Generate simulated data from fitted model and refit the model to the
simulated data

## Usage

``` r
self_test(
  fit_RDS = NULL,
  n = 10,
  seeds = NULL,
  which_seeds = NULL,
  conditional = TRUE,
  map_change = NULL,
  do_parallel = TRUE,
  n_cores = NULL,
  res_dir = NULL,
  wham_location = NULL,
  test_dir = NULL,
  save_inputs = FALSE
)
```

## Arguments

- fit_RDS:

  (required) location of RDS file with fitted WHAM model.

- n:

  the number of simulated data sets to create and fit. Default = 10.

- seeds:

  (optional) vector of seeds to use to generate each simulated data set.

- which_seeds:

  (optional) which `seeds` to use for simualted data set. Useful if
  doing self test in stages.

- conditional:

  T/F whether to fix rather than simulate estimated random effects.
  Deafult = TRUE.

- map_change:

  (optional) list of input\$map elements for altering mapping
  assumptions of fitted model.

- do_parallel:

  T/F whether to do self-test fits in parallel. Requires snowfall and
  parallel packages to be installed. Default = TRUE.

- n_cores:

  (optional) the number of cores to use for parallel fitting.

- res_dir:

  directory where to save individual files for each self test fit.
  Useful if doing self test in stages. If not provided, no files will be
  saved.

- wham_location:

  (optional) location of WHAM package. Useful if not using the WHAM
  installation in the standard library location.

- test_dir:

  (optional) directory for package repository. To be used when the
  function is being called during package testing rather than an
  installed version of WHAM.

- save_inputs:

  T/F whether to save the simulated inputs in res_dir. Default = FALSE.

## Value

a list of two elements. First is the results which is a list (length =
n) of lists with 6 elements: minimized negative log-likelihood, MLEs,
gradient, SSB, F, abundance at age. Second element is the vector of
seeds used for self-test simulations.

## See also

[`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md)
