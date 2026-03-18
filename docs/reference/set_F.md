# Specify configuration for fully-selected fishing mortality

Specify configuration for fully-selected fishing mortality

## Usage

``` r
set_F(input, F_opts = NULL)
```

## Arguments

- input:

  list containing data, parameters, map, and random elements (output
  from
  [`wham::prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md))

- F_opts:

  (optional) named list of initial values for annual fully-selected
  fishing mortality and configuration method for estimation.

  `F_opts` specifies a few as well as the effect on the population.
  Environmental covariate data need not span the same years as the
  fisheries data. It can be `NULL` if no environmental data are to be
  fit. Otherwise, it must be a named list with the following components:

  \$F

  :   matrix (n_years x n_fleets) of (initial) values for fully-selected
      fishing morality.

  \$F_config

  :   integer 1: (default) configure F parameters (on log scale) as an F
      in the initial year and then deviations from one year to the next,
      or 2: configure F parameters as (log) annual values.

  \$map_F

  :   Specify whether to fix any fully-selected F parameters,
      corresponds to `map$F_pars`. integer matrix (n_years x n_fleets).
      Use NA to fix parameters and common integers will make those
      parameters equal. E.g., for a given column (fleet) values
      NA,1,2,... will fix the first year and estimate other years. Use
      unique values for all distinct parameters for each year and fleet.
