# Specify model and parameter configuration for natural mortality

Specify model and parameter configuration for natural mortality

## Usage

``` r
set_M(input, M)
```

## Arguments

- input:

  list containing data, parameters, map, and random elements (output
  from
  [`wham::prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md))

- M:

  (optional) list specifying natural mortality options: model, random
  effects, initial values, and parameters to fix (see details)

  `M` specifies estimation options for natural mortality and can
  overwrite M-at-age values specified in the ASAP data file. If `NULL`,
  the M-at-age matrix from the ASAP data file is used (M fixed, not
  estimated). To estimate M-at-age shared/mirrored among some but not
  all ages, modify `M$means_map` (see vignette for more details). `M` is
  a list with the following entries:

  \$mean_model

  :   Character describing the type of model for M stock and regional
      models for natural mortality. Options are:

      "fixed-M"

      :   Use initial values from ASAP3 dat files or `$initial_means`
          for (mean) M as fixed values. If no ASAP3 files and
          `$initial_means` is not provided, default is M = 0.2 for all
          stocks, regions and ages

      "estimate-M"

      :   estimate one or more (mean) M parameters. Default is to
          estimate a single M shared across all stocks and ages, but use
          `$means_map` to fix or estimate parameters for specific
          stocks, regions, ages.

      "weight-at-age"

      :   specifies M as a function of weight-at-age, \\M_y,a = exp(b0 +
          b1\*log(W_y,a))\\, as in [Lorenzen
          (1996)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1095-8649.1996.tb00060.x)
          and [Miller & Hyun
          (2018)](https://www.nrcresearchpress.com/doi/10.1139/cjfas-2017-0035).
          Default is to estimate a single model shared across all stocks
          and regions, but use `$means_map[s,r,1]` to fix or estimate
          the intercept for specific stocks, regions. See also
          `$logb_prior` and `$initial_b` configuring the slope on log
          scale.

  \$initial_means

  :   array (n_stocks x n_regions x n_ages) of initial/mean M by stock,
      region and age. If `NULL`, initial mean M-at-age values for a
      given stock and region are taken from the first row of the MAA
      matrix in the ASAP data file. If no ASAP data file, M = 0.2 is the
      default. If `$mean_model` is "weight-at-age" only elements for the
      first age (`$initial_means[,,1]`) are used (for the intercept of
      log(M)).

  \$means_map

  :   array (n_stocks x n_regions x n_ages) of NA or integers ( 0 \<=
      max \<= n_stocks \* n_regions \* n_ages) indicating which ages to
      estimate (mean) M and whether to set any ages to be identical.
      E.g. in a model with 2 stock, 2 regions and 6 ages with constant M
      estimated for each stock across regions and ages
      `$M_ages_map[1,,] = 1` and `$M_ages_map[2,,] = 2`.
      `$M_ages_map[1,1,] = c(NA,1,1,2,2,3)` will fix M for age 1 at the
      initial value, and estimates for ages 2 and 3 are identical as are
      those for ages 4 and 5 and different from age 6+ for stock 1 and
      region 1. If `NULL`, specifies all ages fixed at
      `M$initial_means`. If `$mean_model` is "weight-at-age" these are
      used for all stocks and regions and only the elements for the
      first age (`$M_ages_map[,,1]`) are used (for the intercept of
      log(M)).

  \$intial_MAA

  :   array (n_stocks x n_regions x n_years x n_ages) of initial values
      for M at age. Intended to be uses when nothing pertaining to M
      estimated.

  \$b_model

  :   "constant","stock","region", "stock-region" defining whether
      parameter is constant, stock-specific, region-specific, stock- and
      region-specific. Only used if `M$mean_model` = "weight-at-age".

  \$b_prior

  :   T/F, should a N(mu, 0.08) prior (where mu = log(0.305) by default)
      be used on log_b? Based on Fig. 1 and Table 1 (marine fish) in
      [Lorenzen
      (1996)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1095-8649.1996.tb00060.x).
      (Only used if `$mean_model` is "weight-at-age").

  \$intial_b

  :   if any elements of `$mean_model` is "weight-at-age", initial value
      for mean b for weight-at-age model.

  \$re_model

  :   Character matrix (n_stocks x n_regions) of options for time- and
      age-varying (random effects) on M by stock and region. Possible
      values are:

      "none"

      :   (default) No random effects by age or year.

      "iid_a"

      :   uncorrelated M by age, constant in time.

      "iid_y"

      :   uncorrelated M by year, constant all ages.

      "ar1_a"

      :   M correlated by age (AR1), constant in time.

      "ar1_y"

      :   M correlated by year (AR1), constant all ages.

      "iid_ay"

      :   M uncorrelated by year and age (2D).

      "ar1_ay"

      :   M correlated by year and age (2D AR1), as in [Cadigan
          (2016)](https://www.nrcresearchpress.com/doi/10.1139/cjfas-2015-0047).

  \$re_map

  :   array (n_stocks x n_regions x n_ages) of NA and integers (1 \<=
      max \<= n_ages) indicating which ages, for a given stock and
      region, have random effects (not NA) and whether to set RE for any
      ages to be identical. E.g. in a model with 2 stock, 2 regions and
      6 ages, `$re_map[2,1,] = c(NA,1,1,2,2,3)` will not estimate RE for
      age 1, and those for ages 2 and 3 are identical as are those for
      ages 4 and 5 and different from age 6+ for stock 2 and region 1.
      If `NULL`, and `$re_model` specifies M random effects at age, at
      least two ages must be specified for correlation among ages to be
      estimated.

  \$sigma_vals

  :   n_stocks x n_regions matrix Initial standard deviation value to
      use for the M random effects. Values are not used if `M$re_model`
      = "none". Otherwise, a single value. If unspecified all values are
      0.1.

  \$cor_vals

  :   n_stocks x n_regions x 2 array of initial correlation values to
      use for the M deviations. If unspecified all initial values are 0.
      When `M$re_model` =

      "iid_a", "iid_y", "iid_ay" or "none"

      :   values are not used.

      "ar1_a"

      :   first value cor_vals\[s,r,1\] is used.

      "ar1_y"

      :   second value cor_vals\[s,r,2\] is used.

      "ar1_ay"

      :   First is for "age", second is for "year".

  \$sigma_map

  :   n_stocks x n_region matrix of NA or integers indicating which
      random effects sd is estimated and whether to set any to be
      identical. If not supplied a single sd will be estimated for any
      stock and region where \$re_model is other than "none".

  \$cor_map

  :   n_stocks x n_region matrix x 2 array of NA or integers indicating
      which random effects correlation parameters are estimated and
      whether to set any to be identical. If not supplied a single value
      for age and/or year will be estimated for any stock and region
      where \$re_model is other than "none", "iid_a", "iid_y".

## Value

a named list with same elements as the input provided with natural
mortality options modified.

## See also

[`prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md)

## Examples

``` r
if (FALSE) { # \dontrun{
wham.dir <- find.package("wham")
path_to_examples <- system.file("extdata", package="wham")
asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
input <- prepare_wham_input(asap3)
M = list(mean_model = "estimate-M")
input <- set_q(input, M = M) #estimate a constant M parameters
} # }
```
