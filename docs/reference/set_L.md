# Specify model and parameter configuration for "extra" mortality not directly attributed to natural mortality

Specify model and parameter configuration for "extra" mortality not
directly attributed to natural mortality

## Usage

``` r
set_L(input, L)
```

## Arguments

- input:

  list containing data, parameters, map, and random elements (output
  from
  [`wham::prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md))

- L:

  (optional) list specifying "extra" mortality options: model, random
  effects, initial values, and parameters to fix (see details)

  `L` specifies estimation options for "extra" mortality. If `NULL`,
  This mortality source is not used. `L` is a list with the following
  entries:

  \$model

  :   length = n_regions. "extra" mortality model options are:

      "none"

      :   (default) no extra mortality for this region.

      "constant"

      :   estimate a single mean mortality for the region shared across
          all ages

      "iid_re"

      :   estimate independent random effects over years, for the region

      "ar1_re"

      :   estimate random effect correlated over years, for the region

  \$initial_means

  :   Initial/mean L-at-region

  \$sigma_vals

  :   Initial standard deviation by region value to use for the L random
      effects. Values are not used if `L$model` = "none".

  \$cor_vals

  :   Initial correlation values to use for the L random effects. If
      unspecified all initial values are 0
