# Set up observation vector that is used by the model for likelihood calculations and one-step-ahead residuals.

Set up observation vector that is used by the model for likelihood
calculations and one-step-ahead residuals.

## Usage

``` r
set_osa_obs(input)
```

## Arguments

- input:

  list containing data, parameters, map, and random elements (output
  from
  [`wham::prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md))

## Value

the same input list as provided, but with \$obs and \$obsvec configured.
This is run after any changes have been made to the data
