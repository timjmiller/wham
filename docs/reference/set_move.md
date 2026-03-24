# Specify model and parameter configuration for movement when input\$data\$n_regions \> 1

Specify model and parameter configuration for movement when
input\$data\$n_regions \> 1

## Usage

``` r
set_move(input, move)
```

## Arguments

- input:

  list containing data, parameters, map, and random elements (output
  from
  [`wham::prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md))

- move:

  (optional) list specifying movement options: model, random effects,
  initial values, and parameters to fix (see details)

  `move` specifies estimation options for movement. If `NULL`, no
  movement will occur. If there are multiple regions, each stock will be
  modeled separately in different regions without movement. `move` is a
  list with the following entries:

  \$stock_move

  :   length = n_stocks, T/F whether each stock can move. If not
      provided then movement will be defined below for all stocks.

  \$separable

  :   length = n_stocks, T/F whether movement should be modeled
      separably from mortality or both occuring simultaneously.

  \$mean_model

  :   matrix (n_regions x (n_regions-1)): model options for fixed
      effects (mean and possibly variance) for each movement parameter
      are:

      "none"

      :   (default) no movement between regions.

      "constant"

      :   estimate a single movement rate to each region shared across
          all stocks, seasons, ages, years

      "season"

      :   estimate movement rates to each region for each season shared
          across all stocks, ages, years

      "stock_constant"

      :   estimate a movement rate for each stock to each region shared
          across all seasons, ages, years

      "stock_season"

      :   estimate a movement rate for each stock each season to each
          region shared across all ages, years

  \$age_re

  :   matrix (n_regions x (n_regions-1)): options for age random effects
      (for each mean parameter defined in `move$mean_model`):

      "none"

      :   (default) no movement rate random effects by age.

      "iid"

      :   independent movement rate random effects by age.

      "ar1"

      :   allow first order autoregressive correlation of movement rate
          random effects by age.

  \$year_re

  :   matrix (n_regions x (n_regions-1)): options for yearly random
      effects (for each mean parameter defined in `move$mean_model`):

      "none"

      :   (default) no movement rate random effects by year.

      "iid"

      :   independent movement rate random effects by year.

      "ar1"

      :   allow first order autoregressive correlation of movement rate
          random effects by year.

  \$prior_sigma

  :   array (n_stocks x n_seasons x n_regions x n_regions - 1) of sd
      parameters for normal priors on mean movement parameters on
      transformed scale (-Inf,Inf)

  \$use_prior

  :   array (n_stocks x n_seasons x n_regions x n_regions - 1) 0/1
      indicator whether to include prior for mean movement parameters in
      joint log-likelihood.

  \$can_move

  :   array (n_stocks x n_seasons x n_regions x n_regions) 0/1 indicator
      whether movement can occur from one region to another.

  \$must_move

  :   array (n_stocks x n_seasons x n_regions) 0/1 indicator whether
      movement from region must occur.

  \$mean_vals

  :   array (n_stocks x n_seasons x n_regions x n_regions-1) of initial
      movement rate parameters \*from\* each region. Usage depends on
      `move$mean_model`.

  \$sigma_vals

  :   array (n_stocks x n_seasons x n_regions x n_regions -1) of initial
      standard deviations to use for random effects. Usage depends on
      `move$age_re` and `move$year_re`.

  \$cor_vals

  :   array (n_stocks x n_seasons x n_regions x n_regions - 1x 2) of
      initial correlation values to use for random effects. Usage
      depends on `move$age_re` and `move$year_re`. cor_vals\[,,,,1\] is
      for correlation with age, and cor_vals\[,,,,2\] is for correlation
      with year.
