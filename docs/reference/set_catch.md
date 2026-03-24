# Specify catch selectivity blocks and aggregate and age composition observations for catch

Specify catch selectivity blocks and aggregate and age composition
observations for catch

## Usage

``` r
set_catch(input, catch_info = NULL)
```

## Arguments

- input:

  list containing data, parameters, map, and random elements (output
  from
  [`prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md))

- catch_info:

  (optional) list specifying various aspects about catch by fleet (see
  details)

  `catch_info` specifies observations, and various configuration options
  for fleet-specific catch observations and will overwrite attributes
  specified in the ASAP data file. If `NULL`, all settings from the ASAP
  data file or basic_info are used. `catch_info` is a list with any of
  the following entries:

  \$n_fleets

  :   number of fleets

  \$fleet_regions

  :   vector (n_fleets) of regions where each fleet operates.

  \$fleet_seasons

  :   matrix (n_fleets x n_seasons) of 0/1 values flagging which seasons
      each fleet operates.

  \$fleet_names

  :   character vector (n_fleets) of names for fleets. Used for naming
      results in plots and tables.

  \$agg_catch

  :   matrix (n_years_model x n_fleets) of annual aggregate catches by
      fleet.

  \$agg_catch_cv

  :   matrix (n_years_model x n_fleets) of CVs for annual aggregate
      catches by fleet.

  \$catch_paa

  :   array (n_fleets x n_years_model x n_ages) of annual catch
      proportions at age by fleet.

  \$use_catch_paa

  :   matrix (n_years_model x n_fleets) of 0/1 values flagging whether
      to use proportions at age observations.

  \$catch_Neff

  :   matrix (n_years_model x n_fleets) of effective sample sizes for
      proportions at age observations.

  \$waa_pointer_fleets

  :   vector (n_fleets) of itegers indicated waa to use for each fleet.

  \$selblock_pointer_fleets

  :   matrix (n_years_model x n_fleets) of itegers indicated selblocks
      to use.

  \$initial_catch_sd_scale

  :   vector (n_fleets) of scalar multipliers of annual log-observation
      standard deviation. Default = 1.

  \$map_catch_sd_scale

  :   integer vector (n_fleets) specifying which sd scalar parameters to
      fix. Use `NA` to fix a parameter and integers to estimate. Use the
      same integer for multiple fleets to estimate a shared scalar
      parameter.

## Value

a named list with same elements as the input provided with catch
observations and fleet options modified.

## See also

[`prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md)

## Examples

``` r
if (FALSE) { # \dontrun{
wham.dir <- find.package("wham")
path_to_examples <- system.file("extdata", package="wham")
asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
input <- prepare_wham_input(asap3)
input <- set_catch(input, catch_info = list(agg_catch = newcatch)) #constant catch of 500 mt
} # }
```
