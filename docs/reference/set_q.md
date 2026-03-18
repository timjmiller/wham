# Specify model and parameter configuration for catchability

Specify model and parameter configuration for catchability

## Usage

``` r
set_q(input, catchability = NULL)
```

## Arguments

- input:

  list containing data, parameters, map, and random elements (output
  from
  [`prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md))

- catchability:

  (optional) list specifying options for numbers-at-age random effects,
  initial parameter values, and recruitment model (see details)

  `catchability` specifies options for catchability. If `NULL` and
  `asap3` is not NULL, a single catchability parameter for each index is
  used with initial values specified in ASAP file. If both are NULL,
  initial catchabilities for all indices = 0.3. Otherwise, it is a list
  with the following optional entries:

  \$re

  :   Time-varying (random effects) for each index. Vector with length =
      number of indices. Each entry of `catchability$re` must be one of
      the following options:

      "none"

      :   (default) are constant

      "iid"

      :   vary by year and age/par, but uncorrelated

      "ar1"

      :   correlated by year (AR1)

  \$initial_q

  :   Initial catchabilities for each index. vector length = number of
      indices. Will override values provided in `basic_info$q`. If
      `basic_info$q` and `asap3` are not provided, default q values are
      0.3.

  \$q_lower

  :   Lower bound for catchabilities for each index. vector length =
      number of indices. For indices with NULL components default lower
      values are 0.

  \$q_upper

  :   Upper bound for catchabilities for each index. vector length =
      number of indices. For indices with NULL components default lower
      values are 1000.

  \$prior_sd

  :   vector of NA and standard deviations to use for gaussian prior on
      logit transform of catchability parameter. Length = number of
      indices. Indices with non-NA values will have mean logit q as a
      random effect with mean determined by logit transform of
      `catchability$initial_q` and sigma as standard error.

  \$sigma_val

  :   Vector of initial standard deviation values to use for annual
      random effects for each index. Values are not used if `q$re` =
      "none". Otherwise, a single value for all indices.

  \$sigma_map

  :   Specify which sigma parameters to fix for the random effect
      deviations. Must be a vector with length = number of indices. Use
      `NA` to fix a parameter and integers to estimate. Use the same
      integer for multiple indices to share the same sigma parameter.
      Not used if `re = 'none'` for all indices.

  \$cor_vals

  :   Vector of initial correlation values to use for annual random
      effects for each index. If unspecified all initial values are 0.
      Only used if `q$re` = "ar1"

  \$cor_map

  :   Specify which ar1 correlation parameters to fix for the random
      effect deviations. Must be a vector with length = number of
      indices. If `re = 'ar1'`, each element (index) must be a single
      value. Use `NA` to fix a parameter and integers to estimate. Use
      the same integer for multiple indices to share the same
      correlation parameter. Not used if `re = 'none'` or `re = 'iid'`
      for all indices.

## Value

a named list with same elements as the input provided with catchability
options modified.

## See also

[`prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md)

## Examples

``` r
if (FALSE) { # \dontrun{
wham.dir <- find.package("wham")
path_to_examples <- system.file("extdata", package="wham")
asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
input <- prepare_wham_input(asap3)
catchability <- list(re = c("iid", "none"))
input <- set_q(input, catchability = catchability) #independent time-varying random effects on q for first survey.
} # }
```
