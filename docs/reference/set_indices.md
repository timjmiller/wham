# Specify index selectivity blocks and aggregate and age composition observations for indices

Specify index selectivity blocks and aggregate and age composition
observations for indices

## Usage

``` r
set_indices(input, index_info = NULL)
```

## Arguments

- input:

  list containing data, parameters, map, and random elements (output
  from
  [`prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md))

- index_info:

  (optional) list specifying various aspects about catch by indices (see
  details)

  `index_info` specifies observations, and various configuration options
  for index-specific catch observations and will overwrite attributes
  specified in the ASAP data file. If `NULL`, all settings from the ASAP
  data file or basic_info are used. `index_info` is a list with any of
  the following entries:

  \$n_indices

  :   number of indices

  \$index_regions

  :   vector (n_indices) of regions where each fleet operates.

  \$index_seasons

  :   vector (n_indices) of 0/1 values flagging which seasons each index
      occurs.

  \$index_names

  :   character vector (n_indices) of names for indices Used for naming
      results in plots and tables.

  \$agg_indices

  :   matrix (n_years_model x n_indices) of annual aggregate index
      catches.

  \$agg_index_cv

  :   matrix (n_years_model x n_indices) of CVs for annual aggregate
      index catches.

  \$fracyr_indices

  :   matrix (n_years_model x n_indices) of fractions of year at which
      index occurs within the season (difference between time of survey
      and time at start of season).

  \$use_indices

  :   matrix (n_years_model x n_indices) of 0/1 values flagging whether
      to use aggregate observations.

  \$units_indices

  :   matrix (n_years_model x n_indices) of 1/2 values flagging whether
      aggregate observations are biomass (1) or numbers (2).

  \$index_paa

  :   array (n_indices x n_years_model x n_ages) of annual catch
      proportions at age by index.

  \$use_index_paa

  :   matrix (n_years_model x n_indices) of 0/1 values flagging whether
      to use proportions at age observations.

  \$units_index_paa

  :   matrix (n_years_model x n_indices) of 1/2 values flagging whether
      composition observations are biomass (1) or numbers (2).

  \$index_Neff

  :   matrix (n_years_model x n_indices) of effective sample sizes for
      proportions at age observations.

  \$waa_pointer_indices

  :   vector (n_indices) of itegers indicated waa to use for each index.

  \$selblock_pointer_indices

  :   matrix (n_years_model x n_indices) of itegers indicated selblocks
      to use.

  \$initial_index_sd_scale

  :   vector (n_indices) of scalar multipliers of annual log-observation
      standard deviation. Default = 1.

  \$map_index_sd_scale

  :   integer vector (n_indices) specifying which sd scalar parameters
      to fix. Use `NA` to fix a parameter and integers to estimate. Use
      the same integer for multiple indices to estimate a shared scalar
      parameter.
