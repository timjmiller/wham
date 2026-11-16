# Prepare input data and parameters for WHAM model

Prepares data and parameter settings for
[`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md),
optionally using an ASAP3 data file read into R by
[`read_asap3_dat`](https://timjmiller.github.io/wham/reference/read_asap3_dat.md).
By default, this will set up a SCAA version like
[ASAP](https://www.nefsc.noaa.gov/nft/ASAP.html) without any penalties
or random effects. If any asap3, basic_info, catch_info, and index_info
arguments are not provided, some respective arbitrary assumptions are
made. The various options described below, such as `NAA_re` and
`selectivity`, can still be used without the asap3 object.

## Usage

``` r
prepare_wham_input(
  asap3 = NULL,
  model_name = "WHAM for unnamed stock",
  recruit_model = 2,
  ecov = NULL,
  selectivity = NULL,
  M = NULL,
  NAA_re = NULL,
  catchability = NULL,
  age_comp = NULL,
  move = NULL,
  L = NULL,
  F = NULL,
  catch_info = NULL,
  index_info = NULL,
  basic_info = NULL
)
```

## Arguments

- asap3:

  (optional) list containing data and parameters (output from
  [`read_asap3_dat`](https://timjmiller.github.io/wham/reference/read_asap3_dat.md))

- model_name:

  character, name of stock/model

- recruit_model:

  numeric, option to specify stock-recruit model (see details)

- ecov:

  (optional) named list of environmental covariate data and parameters
  (see
  [`set_ecov`](https://timjmiller.github.io/wham/reference/set_ecov.md)
  for full details)

- selectivity:

  (optional) list specifying selectivity options by block: models,
  initial values, parameters to fix, and random effects (see
  [`set_selectivity`](https://timjmiller.github.io/wham/reference/set_selectivity.md)
  for full details)

- M:

  (optional) list specifying natural mortality options: model, random
  effects, initial values, and parameters to fix (see
  [`set_M`](https://timjmiller.github.io/wham/reference/set_M.md) for
  full details)

- NAA_re:

  (optional) list specifying options for random effect on
  numbers-at-age, initial numbers at age, and recruitment model (see
  [`set_NAA`](https://timjmiller.github.io/wham/reference/set_NAA.md)
  for full details)

- catchability:

  (optional) list specifying options for priors and random effects on
  catchability (see
  [`set_q`](https://timjmiller.github.io/wham/reference/set_q.md) for
  full details)

- age_comp:

  (optional) character or named list, specifies age composition model
  for fleet(s) and indices (see
  [`set_age_comp`](https://timjmiller.github.io/wham/reference/set_age_comp.md)
  for full details)

- move:

  (optional) list specifying movement/migration options for models with
  more than 1 region (see
  [`set_move`](https://timjmiller.github.io/wham/reference/set_move.md)
  for full details)

- L:

  (optional) list specifying "extra" mortality options (see
  [`set_L`](https://timjmiller.github.io/wham/reference/set_L.md) for
  full details)

- F:

  (optional) list specifying fishing mortality options (see
  [`set_F`](https://timjmiller.github.io/wham/reference/set_F.md) for
  full details)

- catch_info:

  (optional) list specifying catch information (see
  [`set_catch`](https://timjmiller.github.io/wham/reference/set_catch.md)
  for full details)

- index_info:

  (optional) list specifying index informaiton (see
  [`set_indices`](https://timjmiller.github.io/wham/reference/set_indices.md)
  for full details)

- basic_info:

  (optional) list of basic population information for use when asap3 is
  not provided (see details)

## Value

a named list with the following components:

- `data`:

  Named list of data, passed to
  [`TMB::MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html)

- `par`:

  Named list of parameters, passed to
  [`TMB::MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html)

- `map`:

  Named list defining how to optionally collect and fix parameters,
  passed to
  [`TMB::MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html)

- `random`:

  Character vector of parameters to treat as random effects, passed to
  [`TMB::MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html)

- `years`:

  Numeric vector of years to fit WHAM model (specified in ASAP3 .dat
  file)

- `ages.lab`:

  Character vector of age labels, ending with plus-group (specified in
  ASAP3 .dat file)

- `model_name`:

  Character, name of stock/model (specified in call to
  `prepare_wham_input`)

- `log`:

  list of character strings attempting to describe the input and what
  the model assumes.

## Details

`recruit_model` specifies the stock-recruit model. See `wham.cpp` for
implementation.

- = 1:

  SCAA (without NAA_re option specified) or Random walk (if
  NAA_re\$sigma specified), i.e. predicted recruitment in year i =
  recruitment in year i-1

- = 2:

  (default) Random effects about constant mean

- = 3:

  Beverton-Holt

- = 4:

  Ricker

Note:
[`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md)
allows fitting a SCAA model (`NAA_re = NULL`), which estimates
recruitment in every year as separate fixed effect parameters, but in
that case no stock-recruit function is estimated. A warning message is
printed if `recruit_model > 2` and `NAA_re = NULL`. If you wish to use a
stock-recruit function for expected recruitment, choose recruitment
deviations as random effects, either only age-1
(`NAA_re = list(sigma="rec")`) or all ages
(`NAA_re = list(sigma="rec+1")`, "full state-space" model). See below
for details on `NAA_re` specification.

`ecov` specifies any environmental covariate data and models. See
[`set_ecov`](https://timjmiller.github.io/wham/reference/set_ecov.md)
for full details.

`selectivity` specifies options for selectivity, to overwrite existing
options specified in the ASAP data file. See
[`set_selectivity`](https://timjmiller.github.io/wham/reference/set_selectivity.md)
for full details.

`M` specifies estimation options for natural mortality and can overwrite
M-at-age values specified in the ASAP data file. See
[`set_M`](https://timjmiller.github.io/wham/reference/set_M.md) for full
details.

`NAA_re` specifies options for random effects on numbers-at-age (NAA,
i.e. state-space model or not). See
[`set_NAA`](https://timjmiller.github.io/wham/reference/set_NAA.md) for
full details. If `NULL`, a traditional statistical catch-at-age model is
fit.

`catchability` specifies options for catchability. See
[`set_q`](https://timjmiller.github.io/wham/reference/set_q.md) for full
details. If `NULL` and `asap3` is not NULL, a single catchability
parameter for each index is used with initial values specified in ASAP
file. If both are NULL, initial catchabilities for all indices = 0.3.

`move` specifies options for movement if there are multiple regions. See
[`set_move`](https://timjmiller.github.io/wham/reference/set_move.md)
for full details.

`age_comp` specifies the age composition models for fleet(s) and
indices. See
[`set_age_comp`](https://timjmiller.github.io/wham/reference/set_age_comp.md)
for full details.

`catch_info` is an optional list of fishery catch information that can
be used to set up these types of observations when there is no asap3
file given. See
[`set_catch`](https://timjmiller.github.io/wham/reference/set_catch.md)
for full details. Useful for setting up an operating model to simulate
population processes and observations. Also can be useful for setting up
the structure of assessment model when asap3 has not been used.

`index_info` is an optional list of survey/index information that can be
used to set up these types of observations when there is no asap3 file
given. See
[`set_indices`](https://timjmiller.github.io/wham/reference/set_indices.md)
for full details. Useful for setting up an operating model to simulate
population processes and observations. Also can be useful for setting up
the structure of assessment model when asap3 has not been used.

`basic_info` is an optional list of information that can be used to set
up the population and types of observations when there is no asap3 file
given. Particularly useful for setting up an operating model to simulate
population processes and observations. Also can be useful for setting up
the structure of assessment model when asap3 has not been used. The
current options are:

- \$ages:

  integer vector of ages (years) with the last being a plus group

- \$years:

  integer vector of years that the population model spans.

- \$n_seasons:

  number of seasons within year.

- \$n_fleets:

  number of fishing fleets.

- \$fracyr_seasons:

  proportions of year for each season within year (sums to 1).

- \$F:

  matrix (length(years) x n_fleets) of annual fishing mortality rates
  for each fleet to initialize the model.

- \$waa:

  array ((n_fleets + n_indices + n_stocks) x length(years) x
  length(ages)) of annual weight at at age for each fleet, each index,
  and spawning biomass for each stock.

- \$maturity:

  array (n_stocks x length(years) x length(ages)) of annual maturity at
  age for estimating spawning biomass for each stock.

- \$fracyr_SSB:

  matrix (n_years x n_stocks) (1 or length(years)) of yearly proportions
  (0-1) of the year at which to calculate spawning biomass.

- \$spawn_seasons:

  vector (n_stocks) of seasons in which each stock spawns.

- \$spawn_regions:

  vector (n_stocks) of regions in which each stock spawns.

- \$NAA_where:

  array (n_stocks x n_regions x n_ages) of 0/1 indicating where
  individuals of each stock may exist on January 1 of each year.

- \$Fbar_ages:

  integer vector of ages to use to average F at age for reported "Fbar"
  across all fleets, by regions and by fleet.

- \$q:

  vector (length(n_indices)) of catchabilities for each of the indices
  to initialize the model.

- \$percentSPR:

  (0-100) percentage of unfished spawning biomass per recruit for
  determining equilibrium fishing mortality reference point

- \$XSPR_input_average_years:

  which years to average inputs to per recruit calculation (selectivity,
  M, WAA, maturity) for static (or prevailing) SPR-based reference
  points. Default is last 5 years (utils::tail(1:length(years),5)).
  Below specific values will override this option.

- \$average_years_L:

  list (length = n_regions), which years to average L (extra mortality)
  for per recruit calculation for static (or prevailing) reference
  points. Default is last 5 years (utils::tail(1:length(years),5))

- \$average_years_M:

  list (length = n_stocks, each is a list with length = n_regions),
  which years to average M for per recruit calculation for static (or
  prevailing) reference points. Default is last 5 years
  (utils::tail(1:length(years),5))

- \$average_years_move:

  list (length = n_stocks, each is a list with length = n_regions),
  which years to average movement parameters for per recruit calculation
  for static (or prevailing) reference points. Default is last 5 years
  (utils::tail(1:length(years),5))

- \$average_years_mat:

  list (length = n_stocks), which years to average maturity for per
  recruit calculation for static (or prevailing) reference points.
  Default is last 5 years (utils::tail(1:length(years),5))

- \$average_years_waassb:

  list (length = n_stocks), which years to average ssb waa for per
  recruit calculation for static (or prevailing) reference points.
  Default is last 5 years (utils::tail(1:length(years),5))

- \$average_years_sel:

  list (length = n_fleets), which years to average selectivity for per
  recruit calculation for static (or prevailing) reference points.
  Default is last 5 years (utils::tail(1:length(years),5))

- \$average_years_waacatch:

  list (length = n_fleets), which years to average catch waa for per
  recruit calculation for static (or prevailing) reference points.
  Default is last 5 years (utils::tail(1:length(years),5))

- \$average_years_SRR:

  list (length = n_stocks), which years to average SRR parameters for
  static (or prevailing) MSY-based reference points. Default is last 5
  years (utils::tail(1:length(years),5))

- \$XSPR_R_avg_yrs:

  which years to average recruitments for calculating SPR-based SSB
  reference points. Default is 1:length(years)

- \$XSPR_R_opt:

  1(3): use annual R estimates(predictions) for annual SSB_XSPR, 2(4):
  use average R estimates(predictions). 5: use bias-corrected expected
  recruitment. For long-term projections, may be important to use
  certain years for XSPR_R_avg_yrs

- \$FXSPR_init:

  which F to initialize internal newton search for annual and static F
  at X percent SPR. Default is 0.5

- \$FMSY_init:

  which F to initialize internal newton search for annual and static
  Fmsy (if a stock-recruit model is assumed). Default is 0.5

- \$simulate_process_error:

  T/F vector (length = 9). When simulating from the model, whether to
  simulate any process errors for (NAA, M, selectivity, q, movement,
  unidentified mortality, q priors, movement priors, Ecov). Only used
  for applicable random effects.

- \$simulate_observation_error:

  T/F vector (length = 3). When simulating from the model, whether to
  simulate catch, index, and ecov observations.

- \$simulate_period:

  T/F vector (length = 2). When simulating from the model, whether to
  simulate base period (model years) and projection period.

- \$bias_correct_process:

  T/F. Perform bias correction of log-normal random effects for NAA.

- \$bias_correct_BRPs:

  T/F. Perform bias correction of analytic SSB/R and Y/R when there is
  bias correction of log-normal NAA. May want to use XSPR_R_opt = 5 for
  long-term projections.

If other arguments to `prepare_wham_input` are provided such as
`selectivity`, `M`, and `age_comp`, the information provided there must
be consistent with `basic_info`. For example the dimensions for number
of years, ages, fleets, and indices.

## See also

[`read_asap3_dat`](https://timjmiller.github.io/wham/reference/read_asap3_dat.md),
[`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md),
[ASAP](https://www.nefsc.noaa.gov/nft/ASAP.html), [Iles & Beverton
(1998)](https://www.sciencedirect.com/science/article/pii/S1385110197000221)

## Examples

``` r
# \donttest{
asap3 = read_asap3_dat(file.path(system.file("extdata", package="wham"), "ex1_SNEMAYT.dat"))
input = prepare_wham_input(asap3)
#> 
#> --Creating input---------------------------------------------------------------------------------------------------------------------
#> basic_info done
#> catch done
#> indices done
#> WAA done
#> NAA done
#> q done
#> selectivity done
#> age_comp done
#> F done
#> M done
#> move done
#> osa_obs done
#> --Input Complete--------------------------------------------------------------------------------------------------------------------
#> -------------------------------------------------------------------------------------------------------------------------------------
#> WHAM version 2.0.0 and forward: 
#> 1) makes major changes to the structure of some data, parameters, and reported objects.
#> 2) decouples random effects for recruitment and random effects for older ages by default. To obtain results from previous versions set NAA_re$decouple_recruitment = FALSE.
#> 3) does not bias correct any log-normal process or observation errors. To configure these, set basic_info$bias_correct_process = TRUE and/or basic_info$bias_correct_observation = TRUE.
#> -------------------------------------------------------------------------------------------------------------------------------------
#> 
#> --USING ASAP3 Inputs-----------------------------------------------------------------------------------------------------------------
#> 1 asap3 dat file was processed. A single stock and region will be assumed. 
#> -------------------------------------------------------------------------------------------------------------------------------------
#> 
#> --Catch------------------------------------------------------------------------------------------------------------------------------ 
#> waa_pointer_fleets determined from ASAP file(s). 
#> Selectivity blocks for each fleet:
#> fleet_1: 1
#> 
#> -------------------------------------------------------------------------------------------------------------------------------------
#> 
#> --Indices----------------------------------------------------------------------------------------------------------------------------
#> waa_pointer_indices determined from ASAP file(s). 
#> Selectivity blocks for each index:
#> index_1: 2
#> index_2: 3
#> 
#> -------------------------------------------------------------------------------------------------------------------------------------
#> 
#> --WAA--------------------------------------------------------------------------------------------------------------------------------
#> basic_info$waa_pointer_M was not provided, so waa_pointer_ssb will be used. 
#> -------------------------------------------------------------------------------------------------------------------------------------
#> 
#> --Selectivity------------------------------------------------------------------------------------------------------------------------
#> (Mean) selectivity block models are:
#> Block 1: logistic
#> Block 2: logistic
#> Block 3: logistic
#> 
#> Random effects options for each selectivity block are:
#> Block 1: none
#> Block 2: none
#> Block 3: none
#> 
#> -------------------------------------------------------------------------------------------------------------------------------------
#> 
#> --Reference points-------------------------------------------------------------------------------------------------------------------
#> For annual SPR-based reference points, corresponding annual conditionally expected recruitments are used. 
#> -------------------------------------------------------------------------------------------------------------------------------------
#> 
mod <- fit_wham(input, do.osa = FALSE, do.retro = FALSE)
#> Order of parameters:
#>  [1] "log_catch_sig_scale"  "log_index_sig_scale"  "log_N1"              
#>  [4] "N1_repars"            "log_NAA_sigma"        "trans_NAA_rho"       
#>  [7] "log_NAA"              "mean_rec_pars"        "logit_q"             
#> [10] "q_re"                 "q_repars"             "q_prior_re"          
#> [13] "logit_selpars"        "selpars_re"           "sel_repars"          
#> [16] "catch_paa_pars"       "index_paa_pars"       "F_pars"              
#> [19] "Mpars"                "M_re"                 "M_repars"            
#> [22] "log_b"                "logR_proj"            "mu_prior_re"         
#> [25] "trans_mu"             "mu_re"                "mu_repars"           
#> [28] "L_re"                 "L_repars"             "Ecov_obs_logsigma"   
#> [31] "Ecov_obs_logsigma_re" "Ecov_obs_sigma_par"   "Ecov_process_pars"   
#> [34] "Ecov_re"              "Ecov_beta_R"          "Ecov_beta_q"         
#> [37] "Ecov_beta_M"          "Ecov_beta_mu"        
#> Not matching template order:
#>  [1] "mean_rec_pars"        "logit_q"              "q_prior_re"          
#>  [4] "q_re"                 "q_repars"             "F_pars"              
#>  [7] "mu_prior_re"          "trans_mu"             "mu_re"               
#> [10] "mu_repars"            "N1_repars"            "log_N1"              
#> [13] "log_NAA_sigma"        "trans_NAA_rho"        "log_NAA"             
#> [16] "logR_proj"            "logit_selpars"        "selpars_re"          
#> [19] "sel_repars"           "catch_paa_pars"       "index_paa_pars"      
#> [22] "log_catch_sig_scale"  "log_index_sig_scale"  "Mpars"               
#> [25] "M_re"                 "M_repars"             "log_b"               
#> [28] "L_re"                 "L_repars"             "Ecov_re"             
#> [31] "Ecov_beta_R"          "Ecov_beta_M"          "Ecov_beta_mu"        
#> [34] "Ecov_beta_q"          "Ecov_process_pars"    "Ecov_obs_logsigma"   
#> [37] "Ecov_obs_logsigma_re" "Ecov_obs_sigma_par"  
#> Your parameter list has been re-ordered.
#> (Disable this warning with checkParameterOrder=FALSE)
#> Using stats::nlminb for optimization with these control parameters:
#> iter.max: 1000
#> eval.max: 1000
#> outer mgc:  1471.042 
#> outer mgc:  1385.041 
#> outer mgc:  1051.061 
#> outer mgc:  1924.258 
#> outer mgc:  863.8018 
#> outer mgc:  750.4628 
#> outer mgc:  604.2889 
#> outer mgc:  829.776 
#> outer mgc:  939.6738 
#> outer mgc:  1460.33 
#> outer mgc:  2076.431 
#> outer mgc:  767.5407 
#> outer mgc:  758.4262 
#> outer mgc:  540.7402 
#> outer mgc:  468.1114 
#> outer mgc:  637.78 
#> outer mgc:  402.3814 
#> outer mgc:  355.8152 
#> outer mgc:  642.6439 
#> outer mgc:  279.326 
#> outer mgc:  308.7317 
#> outer mgc:  346.621 
#> outer mgc:  390.9605 
#> outer mgc:  529.728 
#> outer mgc:  196.518 
#> outer mgc:  322.462 
#> outer mgc:  295.1015 
#> outer mgc:  243.8621 
#> outer mgc:  198.209 
#> outer mgc:  673.7442 
#> outer mgc:  661.9501 
#> outer mgc:  128.5663 
#> outer mgc:  331.3173 
#> outer mgc:  206.7063 
#> outer mgc:  247.115 
#> outer mgc:  221.07 
#> outer mgc:  197.5152 
#> outer mgc:  312.1423 
#> outer mgc:  589.655 
#> outer mgc:  594.1248 
#> outer mgc:  411.8637 
#> outer mgc:  165.7469 
#> outer mgc:  438.8875 
#> outer mgc:  181.2036 
#> outer mgc:  105.4714 
#> outer mgc:  170.7329 
#> outer mgc:  140.4304 
#> outer mgc:  95.78204 
#> outer mgc:  144.2011 
#> outer mgc:  146.4235 
#> outer mgc:  78.11377 
#> outer mgc:  155.4508 
#> outer mgc:  69.23239 
#> outer mgc:  238.2647 
#> outer mgc:  204.1208 
#> outer mgc:  416.8242 
#> outer mgc:  233.6249 
#> outer mgc:  170.9533 
#> outer mgc:  195.3753 
#> outer mgc:  64.74556 
#> outer mgc:  45.96703 
#> outer mgc:  63.55537 
#> outer mgc:  59.57929 
#> outer mgc:  76.98854 
#> outer mgc:  42.94833 
#> outer mgc:  80.84296 
#> outer mgc:  30.07479 
#> outer mgc:  35.75543 
#> outer mgc:  45.79472 
#> outer mgc:  27.00725 
#> outer mgc:  63.60967 
#> outer mgc:  80.48628 
#> outer mgc:  87.66053 
#> outer mgc:  89.87816 
#> outer mgc:  109.8868 
#> outer mgc:  135.1671 
#> outer mgc:  67.46743 
#> outer mgc:  58.35881 
#> outer mgc:  26.45785 
#> outer mgc:  114.2341 
#> outer mgc:  53.61375 
#> outer mgc:  32.38649 
#> outer mgc:  36.94307 
#> outer mgc:  22.49244 
#> outer mgc:  31.94953 
#> outer mgc:  20.70378 
#> outer mgc:  19.66949 
#> outer mgc:  20.27573 
#> outer mgc:  49.65985 
#> outer mgc:  41.68128 
#> outer mgc:  36.22574 
#> outer mgc:  50.33551 
#> outer mgc:  19.55701 
#> outer mgc:  19.54934 
#> outer mgc:  19.00862 
#> outer mgc:  19.1817 
#> outer mgc:  18.31261 
#> outer mgc:  19.07406 
#> outer mgc:  49.35664 
#> outer mgc:  20.18388 
#> outer mgc:  28.72863 
#> outer mgc:  13.61682 
#> outer mgc:  27.1636 
#> outer mgc:  25.1287 
#> outer mgc:  33.52772 
#> outer mgc:  23.47264 
#> outer mgc:  11.39304 
#> outer mgc:  24.67729 
#> outer mgc:  32.74759 
#> outer mgc:  10.18052 
#> outer mgc:  11.68369 
#> outer mgc:  11.24285 
#> outer mgc:  9.944645 
#> outer mgc:  10.35671 
#> outer mgc:  8.116827 
#> outer mgc:  9.595103 
#> outer mgc:  6.186934 
#> outer mgc:  12.04335 
#> outer mgc:  5.589237 
#> outer mgc:  11.20889 
#> outer mgc:  5.340422 
#> outer mgc:  7.953319 
#> outer mgc:  8.25513 
#> outer mgc:  5.994672 
#> outer mgc:  5.037356 
#> outer mgc:  9.712536 
#> outer mgc:  15.88604 
#> outer mgc:  12.33504 
#> outer mgc:  7.695554 
#> outer mgc:  12.56381 
#> outer mgc:  12.06773 
#> outer mgc:  9.213528 
#> outer mgc:  6.201596 
#> outer mgc:  10.37895 
#> outer mgc:  17.19026 
#> outer mgc:  8.346575 
#> outer mgc:  3.865654 
#> outer mgc:  6.126388 
#> outer mgc:  8.160044 
#> outer mgc:  14.81724 
#> outer mgc:  12.78784 
#> outer mgc:  2.317833 
#> outer mgc:  4.373007 
#> outer mgc:  5.114099 
#> outer mgc:  3.147958 
#> outer mgc:  2.909466 
#> outer mgc:  3.07079 
#> outer mgc:  2.402695 
#> outer mgc:  4.027797 
#> outer mgc:  3.162727 
#> outer mgc:  8.711973 
#> outer mgc:  12.64145 
#> outer mgc:  10.04501 
#> outer mgc:  2.225802 
#> outer mgc:  1.571785 
#> outer mgc:  4.656341 
#> outer mgc:  1.71837 
#> outer mgc:  3.98229 
#> outer mgc:  1.443449 
#> outer mgc:  2.89838 
#> outer mgc:  1.894677 
#> outer mgc:  2.917418 
#> outer mgc:  4.287873 
#> outer mgc:  4.681034 
#> outer mgc:  3.900503 
#> outer mgc:  3.075215 
#> outer mgc:  2.414769 
#> outer mgc:  1.875776 
#> outer mgc:  3.980331 
#> outer mgc:  1.728553 
#> outer mgc:  0.5477011 
#> outer mgc:  1.600941 
#> outer mgc:  0.5939936 
#> outer mgc:  0.8511864 
#> outer mgc:  0.4902308 
#> outer mgc:  0.9273049 
#> outer mgc:  1.771194 
#> outer mgc:  1.318009 
#> outer mgc:  0.5481822 
#> outer mgc:  0.613916 
#> outer mgc:  1.577432 
#> outer mgc:  1.33428 
#> outer mgc:  0.6577855 
#> outer mgc:  0.577304 
#> outer mgc:  1.279212 
#> outer mgc:  1.314878 
#> outer mgc:  0.4671219 
#> outer mgc:  0.2004487 
#> outer mgc:  0.7427812 
#> outer mgc:  0.812929 
#> outer mgc:  0.2301668 
#> outer mgc:  0.5115683 
#> outer mgc:  0.5412341 
#> outer mgc:  0.3167261 
#> outer mgc:  0.264794 
#> outer mgc:  0.2659858 
#> outer mgc:  0.2253876 
#> outer mgc:  0.3522257 
#> outer mgc:  0.1350663 
#> outer mgc:  0.6270517 
#> outer mgc:  0.5535786 
#> outer mgc:  0.1350758 
#> outer mgc:  0.4422347 
#> outer mgc:  0.1316514 
#> outer mgc:  0.4088483 
#> outer mgc:  0.1142766 
#> outer mgc:  0.2550084 
#> outer mgc:  0.09382741 
#> outer mgc:  0.1426867 
#> outer mgc:  0.1411335 
#> outer mgc:  0.1100854 
#> outer mgc:  0.1661936 
#> outer mgc:  0.06371693 
#> outer mgc:  0.2694342 
#> outer mgc:  0.07112256 
#> outer mgc:  0.07067425 
#> outer mgc:  0.07000295 
#> outer mgc:  0.09297061 
#> outer mgc:  0.11482 
#> outer mgc:  0.0882267 
#> outer mgc:  0.1185607 
#> outer mgc:  0.1799139 
#> outer mgc:  0.1360937 
#> outer mgc:  0.05101703 
#> outer mgc:  0.08865691 
#> outer mgc:  0.1463756 
#> outer mgc:  0.09396539 
#> outer mgc:  0.03961496 
#> outer mgc:  0.07369442 
#> outer mgc:  0.02277367 
#> outer mgc:  0.05871531 
#> outer mgc:  0.02204421 
#> outer mgc:  0.05910879 
#> outer mgc:  0.05779574 
#> outer mgc:  0.05959536 
#> outer mgc:  0.1038497 
#> outer mgc:  0.06119139 
#> outer mgc:  0.02899737 
#> outer mgc:  0.02179275 
#> outer mgc:  0.03434646 
#> outer mgc:  0.006702836 
#> outer mgc:  0.02355213 
#> outer mgc:  0.007522799 
#> outer mgc:  0.01888813 
#> outer mgc:  0.006849522 
#> outer mgc:  0.03485359 
#> outer mgc:  0.01254535 
#> outer mgc:  0.01063946 
#> outer mgc:  0.003240917 
#> outer mgc:  0.008887901 
#> outer mgc:  0.003054112 
#> outer mgc:  0.007051713 
#> outer mgc:  0.002210686 
#> outer mgc:  0.01120369 
#> outer mgc:  0.001781682 
#> outer mgc:  0.001781682 
#> outer mgc:  0.001781682 
#> outer mgc:  1.197527 
#> outer mgc:  1.197523 
#> outer mgc:  1.119649 
#> outer mgc:  1.122586 
#> outer mgc:  7.410608 
#> outer mgc:  7.401089 
#> outer mgc:  7.211279 
#> outer mgc:  7.20193 
#> outer mgc:  6.843592 
#> outer mgc:  6.834585 
#> outer mgc:  6.690354 
#> outer mgc:  6.683041 
#> outer mgc:  6.599985 
#> outer mgc:  6.59297 
#> outer mgc:  6.514313 
#> outer mgc:  6.50752 
#> outer mgc:  6.399902 
#> outer mgc:  6.392474 
#> outer mgc:  6.279856 
#> outer mgc:  6.273186 
#> outer mgc:  6.204693 
#> outer mgc:  6.198558 
#> outer mgc:  6.114578 
#> outer mgc:  6.108244 
#> outer mgc:  5.96009 
#> outer mgc:  5.952826 
#> outer mgc:  5.648514 
#> outer mgc:  5.641436 
#> outer mgc:  5.295934 
#> outer mgc:  5.289254 
#> outer mgc:  5.120881 
#> outer mgc:  5.114493 
#> outer mgc:  4.943325 
#> outer mgc:  4.937952 
#> outer mgc:  4.869697 
#> outer mgc:  4.863477 
#> outer mgc:  4.769168 
#> outer mgc:  4.762919 
#> outer mgc:  4.537699 
#> outer mgc:  4.53153 
#> outer mgc:  3.809091 
#> outer mgc:  3.802717 
#> outer mgc:  3.017181 
#> outer mgc:  3.010983 
#> outer mgc:  2.521934 
#> outer mgc:  2.518352 
#> outer mgc:  2.425796 
#> outer mgc:  2.422768 
#> outer mgc:  2.350047 
#> outer mgc:  2.347102 
#> outer mgc:  2.268373 
#> outer mgc:  2.264794 
#> outer mgc:  2.176436 
#> outer mgc:  2.172897 
#> outer mgc:  2.069524 
#> outer mgc:  2.066605 
#> outer mgc:  1.96522 
#> outer mgc:  1.96186 
#> outer mgc:  1.843508 
#> outer mgc:  1.839035 
#> outer mgc:  1.687804 
#> outer mgc:  1.683443 
#> outer mgc:  1.445602 
#> outer mgc:  1.441502 
#> outer mgc:  1.27097 
#> outer mgc:  1.268244 
#> outer mgc:  1.146029 
#> outer mgc:  1.1444 
#> outer mgc:  1.075999 
#> outer mgc:  1.074276 
#> outer mgc:  1.011703 
#> outer mgc:  1.01092 
#> outer mgc:  0.9416134 
#> outer mgc:  0.9403846 
#> outer mgc:  0.868526 
#> outer mgc:  0.8674787 
#> outer mgc:  0.8049732 
#> outer mgc:  0.8045814 
#> outer mgc:  0.7387655 
#> outer mgc:  0.7376168 
#> outer mgc:  0.651532 
#> outer mgc:  0.6503983 
#> outer mgc:  0.5801737 
#> outer mgc:  0.5771658 
#> outer mgc:  0.4625025 
#> outer mgc:  0.4597228 
#> outer mgc:  0.2702918 
#> outer mgc:  0.2681184 
#> outer mgc:  0.09519635 
#> outer mgc:  0.09534295 
#> outer mgc:  0.07127049 
#> outer mgc:  0.07111639 
#> outer mgc:  0.1042693 
#> outer mgc:  0.1035835 
#> outer mgc:  0.07449826 
#> outer mgc:  0.07367042 
#> outer mgc:  0.1334469 
#> outer mgc:  0.1331797 
#> outer mgc:  0.06756141 
#> outer mgc:  0.06919739 
#> outer mgc:  0.03405096 
#> outer mgc:  0.0342973 
#> outer mgc:  0.03935506 
#> outer mgc:  0.03919657 
#> outer mgc:  0.08030692 
#> outer mgc:  0.08036545 
#> outer mgc:  0.2373811 
#> outer mgc:  0.2382832 
#> outer mgc:  0.09577961 
#> outer mgc:  0.0949735 
#> outer mgc:  0.2416543 
#> outer mgc:  0.241139 
#> outer mgc:  0.2150332 
#> outer mgc:  0.2155873 
#> outer mgc:  0.16172 
#> outer mgc:  0.1614195 
#> outer mgc:  0.1598343 
#> outer mgc:  0.1610084 
#> outer mgc:  0.2946625 
#> outer mgc:  0.2951873 
#> outer mgc:  0.169666 
#> outer mgc:  0.170069 
#> outer mgc:  0.07787973 
#> outer mgc:  0.07966162 
#> outer mgc:  0.1617097 
#> outer mgc:  0.1621756 
#> outer mgc:  0.2327487 
#> outer mgc:  0.2317454 
#> outer mgc:  0.07811188 
#> outer mgc:  0.07821535 
#> outer mgc:  0.1081498 
#> outer mgc:  0.1092107 
#> outer mgc:  0.4658821 
#> outer mgc:  0.4659442 
#> outer mgc:  0.2498095 
#> outer mgc:  0.2529273 
#> outer mgc:  0.1808717 
#> outer mgc:  0.1840382 
#> outer mgc:  0.1037288 
#> outer mgc:  0.103585 
#> outer mgc:  0.1240407 
#> outer mgc:  0.1246314 
#> outer mgc:  0.1654354 
#> outer mgc:  0.1670617 
#> outer mgc:  0.2529841 
#> outer mgc:  0.2504515 
#> outer mgc:  0.1635201 
#> outer mgc:  0.1609387 
#> outer mgc:  0.1224878 
#> outer mgc:  0.1198909 
#> outer mgc:  0.2352855 
#> outer mgc:  0.235608 
#> outer mgc:  0.1509494 
#> outer mgc:  0.1510012 
#> outer mgc:  0.247227 
#> outer mgc:  0.2465334 
#> outer mgc:  0.1897589 
#> outer mgc:  0.1894565 
#> outer mgc:  0.222673 
#> outer mgc:  0.2230381 
#> outer mgc:  0.1104867 
#> outer mgc:  0.1086092 
#> outer mgc:  0.1301088 
#> outer mgc:  0.1303686 
#> outer mgc:  0.1381907 
#> outer mgc:  0.1387184 
#> outer mgc:  0.3427163 
#> outer mgc:  0.3423138 
#> outer mgc:  0.2628599 
#> outer mgc:  0.2641642 
#> outer mgc:  0.1197283 
#> outer mgc:  0.119822 
#> outer mgc:  0.1955166 
#> outer mgc:  0.1945391 
#> outer mgc:  0.1717368 
#> outer mgc:  0.1730437 
#> outer mgc:  0.1305349 
#> outer mgc:  0.1314873 
#> outer mgc:  0.3040046 
#> outer mgc:  0.3045023 
#> outer mgc:  0.1332542 
#> outer mgc:  0.132391 
#> outer mgc:  0.1445515 
#> outer mgc:  0.1448826 
#> outer mgc:  0.04225243 
#> outer mgc:  0.04237443 
#> outer mgc:  0.02390502 
#> outer mgc:  0.02338735 
#> outer mgc:  0.02708104 
#> outer mgc:  0.02375456 
#> outer mgc:  12.07012 
#> outer mgc:  12.06639 
#> outer mgc:  4.760808 
#> outer mgc:  4.743826 
#> outer mgc:  1.164905 
#> outer mgc:  1.160987 
#> outer mgc:  4.328867 
#> outer mgc:  4.332539 
#> outer mgc:  1.34368 
#> outer mgc:  1.344355 
#> outer mgc:  0.3112966 
#> outer mgc:  0.3126904 
#> outer mgc:  1.45783e-08 
#> outer mgc:  1.197524 
#> outer mgc:  1.197525 
#> outer mgc:  1.121117 
#> outer mgc:  1.121118 
#> outer mgc:  7.408827 
#> outer mgc:  7.402871 
#> outer mgc:  7.209499 
#> outer mgc:  7.203713 
#> outer mgc:  6.841811 
#> outer mgc:  6.836367 
#> outer mgc:  6.689297 
#> outer mgc:  6.684099 
#> outer mgc:  6.599046 
#> outer mgc:  6.593913 
#> outer mgc:  6.51344 
#> outer mgc:  6.508396 
#> outer mgc:  6.398347 
#> outer mgc:  6.393416 
#> outer mgc:  6.278886 
#> outer mgc:  6.274158 
#> outer mgc:  6.203962 
#> outer mgc:  6.199291 
#> outer mgc:  6.11365 
#> outer mgc:  6.108977 
#> outer mgc:  5.95831 
#> outer mgc:  5.953703 
#> outer mgc:  5.646733 
#> outer mgc:  5.642266 
#> outer mgc:  5.294153 
#> outer mgc:  5.289815 
#> outer mgc:  5.119571 
#> outer mgc:  5.115395 
#> outer mgc:  4.942569 
#> outer mgc:  4.938708 
#> outer mgc:  4.868461 
#> outer mgc:  4.864712 
#> outer mgc:  4.767933 
#> outer mgc:  4.764154 
#> outer mgc:  4.535917 
#> outer mgc:  4.532054 
#> outer mgc:  3.807311 
#> outer mgc:  3.803243 
#> outer mgc:  3.0154 
#> outer mgc:  3.011508 
#> outer mgc:  2.521802 
#> outer mgc:  2.518485 
#> outer mgc:  2.425874 
#> outer mgc:  2.422691 
#> outer mgc:  2.350108 
#> outer mgc:  2.347042 
#> outer mgc:  2.268106 
#> outer mgc:  2.265062 
#> outer mgc:  2.176169 
#> outer mgc:  2.173165 
#> outer mgc:  2.06952 
#> outer mgc:  2.066609 
#> outer mgc:  1.964956 
#> outer mgc:  1.962124 
#> outer mgc:  1.841726 
#> outer mgc:  1.83897 
#> outer mgc:  1.686023 
#> outer mgc:  1.683379 
#> outer mgc:  1.443821 
#> outer mgc:  1.441427 
#> outer mgc:  1.270621 
#> outer mgc:  1.268466 
#> outer mgc:  1.146145 
#> outer mgc:  1.144285 
#> outer mgc:  1.075986 
#> outer mgc:  1.074291 
#> outer mgc:  1.012134 
#> outer mgc:  1.01049 
#> outer mgc:  0.9417863 
#> outer mgc:  0.940213 
#> outer mgc:  0.8687442 
#> outer mgc:  0.8672618 
#> outer mgc:  0.8054812 
#> outer mgc:  0.8040746 
#> outer mgc:  0.7388802 
#> outer mgc:  0.7375033 
#> outer mgc:  0.6516466 
#> outer mgc:  0.6502848 
#> outer mgc:  0.5783927 
#> outer mgc:  0.5770684 
#> outer mgc:  0.4607211 
#> outer mgc:  0.4596283 
#> outer mgc:  0.2685103 
#> outer mgc:  0.2678799 
#> outer mgc:  0.09526644 
#> outer mgc:  0.09527266 
#> outer mgc:  0.07117687 
#> outer mgc:  0.07120998 
#> outer mgc:  0.1039717 
#> outer mgc:  0.1038806 
#> outer mgc:  0.0741186 
#> outer mgc:  0.07404936 
#> outer mgc:  0.1333716 
#> outer mgc:  0.1332548 
#> outer mgc:  0.0684155 
#> outer mgc:  0.06834509 
#> outer mgc:  0.03419476 
#> outer mgc:  0.03415384 
#> outer mgc:  0.03929953 
#> outer mgc:  0.03925182 
#> outer mgc:  0.08037272 
#> outer mgc:  0.08029977 
#> outer mgc:  0.237906 
#> outer mgc:  0.237759 
#> outer mgc:  0.09541799 
#> outer mgc:  0.09533448 
#> outer mgc:  0.2414781 
#> outer mgc:  0.2413149 
#> outer mgc:  0.2153896 
#> outer mgc:  0.2152315 
#> outer mgc:  0.1616341 
#> outer mgc:  0.1615052 
#> outer mgc:  0.1604875 
#> outer mgc:  0.1603562 
#> outer mgc:  0.2950222 
#> outer mgc:  0.2948274 
#> outer mgc:  0.1699339 
#> outer mgc:  0.1698014 
#> outer mgc:  0.0788078 
#> outer mgc:  0.07873531 
#> outer mgc:  0.1620047 
#> outer mgc:  0.1618812 
#> outer mgc:  0.2323192 
#> outer mgc:  0.2321741 
#> outer mgc:  0.07819723 
#> outer mgc:  0.07813013 
#> outer mgc:  0.108723 
#> outer mgc:  0.1086384 
#> outer mgc:  0.4660404 
#> outer mgc:  0.4657863 
#> outer mgc:  0.2514608 
#> outer mgc:  0.2512743 
#> outer mgc:  0.1825242 
#> outer mgc:  0.1823864 
#> outer mgc:  0.1037028 
#> outer mgc:  0.103611 
#> outer mgc:  0.124389 
#> outer mgc:  0.1242837 
#> outer mgc:  0.1663181 
#> outer mgc:  0.1661806 
#> outer mgc:  0.2518084 
#> outer mgc:  0.2516256 
#> outer mgc:  0.1622967 
#> outer mgc:  0.16216 
#> outer mgc:  0.1212414 
#> outer mgc:  0.1211351 
#> outer mgc:  0.2355342 
#> outer mgc:  0.2353601 
#> outer mgc:  0.1510407 
#> outer mgc:  0.1509102 
#> outer mgc:  0.2469748 
#> outer mgc:  0.2467853 
#> outer mgc:  0.1896901 
#> outer mgc:  0.1895251 
#> outer mgc:  0.222946 
#> outer mgc:  0.2227656 
#> outer mgc:  0.1095966 
#> outer mgc:  0.1094977 
#> outer mgc:  0.1302944 
#> outer mgc:  0.1301833 
#> outer mgc:  0.1385157 
#> outer mgc:  0.138394 
#> outer mgc:  0.3426467 
#> outer mgc:  0.3423831 
#> outer mgc:  0.2636258 
#> outer mgc:  0.2633995 
#> outer mgc:  0.1198319 
#> outer mgc:  0.1197185 
#> outer mgc:  0.1951158 
#> outer mgc:  0.1949391 
#> outer mgc:  0.1724703 
#> outer mgc:  0.1723115 
#> outer mgc:  0.1310729 
#> outer mgc:  0.1309502 
#> outer mgc:  0.3043745 
#> outer mgc:  0.3041331 
#> outer mgc:  0.1328805 
#> outer mgc:  0.1327641 
#> outer mgc:  0.1447738 
#> outer mgc:  0.1446607 
#> outer mgc:  0.04233245 
#> outer mgc:  0.04229455 
#> outer mgc:  0.02365667 
#> outer mgc:  0.02363525 
#> outer mgc:  0.02542887 
#> outer mgc:  0.02540673 
#> outer mgc:  12.07176 
#> outer mgc:  12.06473 
#> outer mgc:  4.760808 
#> outer mgc:  4.743822 
#> outer mgc:  1.16482 
#> outer mgc:  1.161071 
#> outer mgc:  4.327215 
#> outer mgc:  4.334191 
#> outer mgc:  1.343677 
#> outer mgc:  1.344357 
#> outer mgc:  0.3113822 
#> outer mgc:  0.312607 
#> outer mgc:  3.371151e-12 
#> outer mgc:  1.197524 
#> outer mgc:  1.197525 
#> outer mgc:  1.121117 
#> outer mgc:  1.121118 
#> outer mgc:  7.408827 
#> outer mgc:  7.402871 
#> outer mgc:  7.209499 
#> outer mgc:  7.203713 
#> outer mgc:  6.841811 
#> outer mgc:  6.836367 
#> outer mgc:  6.689297 
#> outer mgc:  6.684099 
#> outer mgc:  6.599045 
#> outer mgc:  6.593913 
#> outer mgc:  6.51344 
#> outer mgc:  6.508396 
#> outer mgc:  6.398347 
#> outer mgc:  6.393416 
#> outer mgc:  6.278886 
#> outer mgc:  6.274158 
#> outer mgc:  6.203962 
#> outer mgc:  6.199291 
#> outer mgc:  6.11365 
#> outer mgc:  6.108977 
#> outer mgc:  5.95831 
#> outer mgc:  5.953703 
#> outer mgc:  5.646733 
#> outer mgc:  5.642266 
#> outer mgc:  5.294153 
#> outer mgc:  5.289815 
#> outer mgc:  5.119571 
#> outer mgc:  5.115395 
#> outer mgc:  4.942569 
#> outer mgc:  4.938708 
#> outer mgc:  4.868461 
#> outer mgc:  4.864712 
#> outer mgc:  4.767932 
#> outer mgc:  4.764154 
#> outer mgc:  4.535917 
#> outer mgc:  4.532054 
#> outer mgc:  3.807311 
#> outer mgc:  3.803243 
#> outer mgc:  3.0154 
#> outer mgc:  3.011508 
#> outer mgc:  2.521802 
#> outer mgc:  2.518485 
#> outer mgc:  2.425874 
#> outer mgc:  2.422691 
#> outer mgc:  2.350108 
#> outer mgc:  2.347042 
#> outer mgc:  2.268106 
#> outer mgc:  2.265062 
#> outer mgc:  2.176169 
#> outer mgc:  2.173165 
#> outer mgc:  2.06952 
#> outer mgc:  2.066609 
#> outer mgc:  1.964956 
#> outer mgc:  1.962124 
#> outer mgc:  1.841726 
#> outer mgc:  1.83897 
#> outer mgc:  1.686023 
#> outer mgc:  1.683379 
#> outer mgc:  1.443821 
#> outer mgc:  1.441427 
#> outer mgc:  1.270621 
#> outer mgc:  1.268466 
#> outer mgc:  1.146145 
#> outer mgc:  1.144285 
#> outer mgc:  1.075986 
#> outer mgc:  1.074291 
#> outer mgc:  1.012134 
#> outer mgc:  1.01049 
#> outer mgc:  0.9417863 
#> outer mgc:  0.940213 
#> outer mgc:  0.8687442 
#> outer mgc:  0.8672618 
#> outer mgc:  0.8054812 
#> outer mgc:  0.8040746 
#> outer mgc:  0.7388802 
#> outer mgc:  0.7375033 
#> outer mgc:  0.6516466 
#> outer mgc:  0.6502849 
#> outer mgc:  0.5783927 
#> outer mgc:  0.5770684 
#> outer mgc:  0.4607211 
#> outer mgc:  0.4596283 
#> outer mgc:  0.2685103 
#> outer mgc:  0.2678799 
#> outer mgc:  0.09526644 
#> outer mgc:  0.09527266 
#> outer mgc:  0.07117687 
#> outer mgc:  0.07120998 
#> outer mgc:  0.1039717 
#> outer mgc:  0.1038806 
#> outer mgc:  0.0741186 
#> outer mgc:  0.07404936 
#> outer mgc:  0.1333716 
#> outer mgc:  0.1332548 
#> outer mgc:  0.06841549 
#> outer mgc:  0.0683451 
#> outer mgc:  0.03419476 
#> outer mgc:  0.03415384 
#> outer mgc:  0.03929953 
#> outer mgc:  0.03925182 
#> outer mgc:  0.08037272 
#> outer mgc:  0.08029977 
#> outer mgc:  0.237906 
#> outer mgc:  0.237759 
#> outer mgc:  0.09541799 
#> outer mgc:  0.09533448 
#> outer mgc:  0.2414781 
#> outer mgc:  0.2413149 
#> outer mgc:  0.2153896 
#> outer mgc:  0.2152315 
#> outer mgc:  0.1616341 
#> outer mgc:  0.1615052 
#> outer mgc:  0.1604875 
#> outer mgc:  0.1603562 
#> outer mgc:  0.2950222 
#> outer mgc:  0.2948274 
#> outer mgc:  0.1699339 
#> outer mgc:  0.1698014 
#> outer mgc:  0.07880779 
#> outer mgc:  0.07873531 
#> outer mgc:  0.1620047 
#> outer mgc:  0.1618812 
#> outer mgc:  0.2323192 
#> outer mgc:  0.2321741 
#> outer mgc:  0.07819723 
#> outer mgc:  0.07813013 
#> outer mgc:  0.108723 
#> outer mgc:  0.1086384 
#> outer mgc:  0.4660404 
#> outer mgc:  0.4657863 
#> outer mgc:  0.2514608 
#> outer mgc:  0.2512743 
#> outer mgc:  0.1825242 
#> outer mgc:  0.1823864 
#> outer mgc:  0.1037028 
#> outer mgc:  0.103611 
#> outer mgc:  0.124389 
#> outer mgc:  0.1242837 
#> outer mgc:  0.1663181 
#> outer mgc:  0.1661806 
#> outer mgc:  0.2518084 
#> outer mgc:  0.2516256 
#> outer mgc:  0.1622967 
#> outer mgc:  0.16216 
#> outer mgc:  0.1212414 
#> outer mgc:  0.1211351 
#> outer mgc:  0.2355342 
#> outer mgc:  0.2353601 
#> outer mgc:  0.1510407 
#> outer mgc:  0.1509102 
#> outer mgc:  0.2469748 
#> outer mgc:  0.2467853 
#> outer mgc:  0.1896901 
#> outer mgc:  0.1895251 
#> outer mgc:  0.222946 
#> outer mgc:  0.2227656 
#> outer mgc:  0.1095966 
#> outer mgc:  0.1094977 
#> outer mgc:  0.1302944 
#> outer mgc:  0.1301833 
#> outer mgc:  0.1385157 
#> outer mgc:  0.138394 
#> outer mgc:  0.3426467 
#> outer mgc:  0.3423831 
#> outer mgc:  0.2636258 
#> outer mgc:  0.2633995 
#> outer mgc:  0.1198319 
#> outer mgc:  0.1197185 
#> outer mgc:  0.1951158 
#> outer mgc:  0.1949391 
#> outer mgc:  0.1724703 
#> outer mgc:  0.1723115 
#> outer mgc:  0.1310729 
#> outer mgc:  0.1309502 
#> outer mgc:  0.3043745 
#> outer mgc:  0.3041331 
#> outer mgc:  0.1328805 
#> outer mgc:  0.1327641 
#> outer mgc:  0.1447738 
#> outer mgc:  0.1446607 
#> outer mgc:  0.04233245 
#> outer mgc:  0.04229455 
#> outer mgc:  0.02365667 
#> outer mgc:  0.02363525 
#> outer mgc:  0.02542886 
#> outer mgc:  0.02540674 
#> outer mgc:  12.07176 
#> outer mgc:  12.06473 
#> outer mgc:  4.760808 
#> outer mgc:  4.743822 
#> outer mgc:  1.16482 
#> outer mgc:  1.161071 
#> outer mgc:  4.327215 
#> outer mgc:  4.334191 
#> outer mgc:  1.343677 
#> outer mgc:  1.344357 
#> outer mgc:  0.3113822 
#> outer mgc:  0.312607 
#> outer mgc:  1.40215e-12 
#> outer mgc:  1.40215e-12 
#> outer mgc:  1.40215e-12 
#> outer mgc:  1.197524 
#> outer mgc:  1.197525 
#> outer mgc:  1.121117 
#> outer mgc:  1.121118 
#> outer mgc:  7.408827 
#> outer mgc:  7.402871 
#> outer mgc:  7.209499 
#> outer mgc:  7.203713 
#> outer mgc:  6.841811 
#> outer mgc:  6.836367 
#> outer mgc:  6.689297 
#> outer mgc:  6.684099 
#> outer mgc:  6.599045 
#> outer mgc:  6.593913 
#> outer mgc:  6.51344 
#> outer mgc:  6.508396 
#> outer mgc:  6.398347 
#> outer mgc:  6.393416 
#> outer mgc:  6.278886 
#> outer mgc:  6.274158 
#> outer mgc:  6.203962 
#> outer mgc:  6.199291 
#> outer mgc:  6.11365 
#> outer mgc:  6.108977 
#> outer mgc:  5.95831 
#> outer mgc:  5.953703 
#> outer mgc:  5.646733 
#> outer mgc:  5.642266 
#> outer mgc:  5.294153 
#> outer mgc:  5.289815 
#> outer mgc:  5.119571 
#> outer mgc:  5.115395 
#> outer mgc:  4.942569 
#> outer mgc:  4.938708 
#> outer mgc:  4.868461 
#> outer mgc:  4.864712 
#> outer mgc:  4.767932 
#> outer mgc:  4.764154 
#> outer mgc:  4.535917 
#> outer mgc:  4.532054 
#> outer mgc:  3.807311 
#> outer mgc:  3.803243 
#> outer mgc:  3.0154 
#> outer mgc:  3.011508 
#> outer mgc:  2.521802 
#> outer mgc:  2.518485 
#> outer mgc:  2.425874 
#> outer mgc:  2.422691 
#> outer mgc:  2.350108 
#> outer mgc:  2.347042 
#> outer mgc:  2.268106 
#> outer mgc:  2.265062 
#> outer mgc:  2.176169 
#> outer mgc:  2.173165 
#> outer mgc:  2.06952 
#> outer mgc:  2.066609 
#> outer mgc:  1.964956 
#> outer mgc:  1.962124 
#> outer mgc:  1.841726 
#> outer mgc:  1.83897 
#> outer mgc:  1.686023 
#> outer mgc:  1.683379 
#> outer mgc:  1.443821 
#> outer mgc:  1.441427 
#> outer mgc:  1.270621 
#> outer mgc:  1.268466 
#> outer mgc:  1.146145 
#> outer mgc:  1.144285 
#> outer mgc:  1.075986 
#> outer mgc:  1.074291 
#> outer mgc:  1.012134 
#> outer mgc:  1.01049 
#> outer mgc:  0.9417863 
#> outer mgc:  0.940213 
#> outer mgc:  0.8687442 
#> outer mgc:  0.8672618 
#> outer mgc:  0.8054812 
#> outer mgc:  0.8040746 
#> outer mgc:  0.7388802 
#> outer mgc:  0.7375033 
#> outer mgc:  0.6516466 
#> outer mgc:  0.6502849 
#> outer mgc:  0.5783927 
#> outer mgc:  0.5770684 
#> outer mgc:  0.4607211 
#> outer mgc:  0.4596283 
#> outer mgc:  0.2685103 
#> outer mgc:  0.2678799 
#> outer mgc:  0.09526644 
#> outer mgc:  0.09527266 
#> outer mgc:  0.07117687 
#> outer mgc:  0.07120998 
#> outer mgc:  0.1039717 
#> outer mgc:  0.1038806 
#> outer mgc:  0.0741186 
#> outer mgc:  0.07404936 
#> outer mgc:  0.1333716 
#> outer mgc:  0.1332548 
#> outer mgc:  0.06841549 
#> outer mgc:  0.0683451 
#> outer mgc:  0.03419476 
#> outer mgc:  0.03415384 
#> outer mgc:  0.03929953 
#> outer mgc:  0.03925182 
#> outer mgc:  0.08037272 
#> outer mgc:  0.08029977 
#> outer mgc:  0.237906 
#> outer mgc:  0.237759 
#> outer mgc:  0.09541799 
#> outer mgc:  0.09533448 
#> outer mgc:  0.2414781 
#> outer mgc:  0.2413149 
#> outer mgc:  0.2153896 
#> outer mgc:  0.2152315 
#> outer mgc:  0.1616341 
#> outer mgc:  0.1615052 
#> outer mgc:  0.1604875 
#> outer mgc:  0.1603562 
#> outer mgc:  0.2950222 
#> outer mgc:  0.2948274 
#> outer mgc:  0.1699339 
#> outer mgc:  0.1698014 
#> outer mgc:  0.07880779 
#> outer mgc:  0.07873531 
#> outer mgc:  0.1620047 
#> outer mgc:  0.1618812 
#> outer mgc:  0.2323192 
#> outer mgc:  0.2321741 
#> outer mgc:  0.07819723 
#> outer mgc:  0.07813013 
#> outer mgc:  0.108723 
#> outer mgc:  0.1086384 
#> outer mgc:  0.4660404 
#> outer mgc:  0.4657863 
#> outer mgc:  0.2514608 
#> outer mgc:  0.2512743 
#> outer mgc:  0.1825242 
#> outer mgc:  0.1823864 
#> outer mgc:  0.1037028 
#> outer mgc:  0.103611 
#> outer mgc:  0.124389 
#> outer mgc:  0.1242837 
#> outer mgc:  0.1663181 
#> outer mgc:  0.1661806 
#> outer mgc:  0.2518084 
#> outer mgc:  0.2516256 
#> outer mgc:  0.1622967 
#> outer mgc:  0.16216 
#> outer mgc:  0.1212414 
#> outer mgc:  0.1211351 
#> outer mgc:  0.2355342 
#> outer mgc:  0.2353601 
#> outer mgc:  0.1510407 
#> outer mgc:  0.1509102 
#> outer mgc:  0.2469748 
#> outer mgc:  0.2467853 
#> outer mgc:  0.1896901 
#> outer mgc:  0.1895251 
#> outer mgc:  0.222946 
#> outer mgc:  0.2227656 
#> outer mgc:  0.1095966 
#> outer mgc:  0.1094977 
#> outer mgc:  0.1302944 
#> outer mgc:  0.1301833 
#> outer mgc:  0.1385157 
#> outer mgc:  0.138394 
#> outer mgc:  0.3426467 
#> outer mgc:  0.3423831 
#> outer mgc:  0.2636258 
#> outer mgc:  0.2633995 
#> outer mgc:  0.1198319 
#> outer mgc:  0.1197185 
#> outer mgc:  0.1951158 
#> outer mgc:  0.1949391 
#> outer mgc:  0.1724703 
#> outer mgc:  0.1723115 
#> outer mgc:  0.1310729 
#> outer mgc:  0.1309502 
#> outer mgc:  0.3043745 
#> outer mgc:  0.3041331 
#> outer mgc:  0.1328805 
#> outer mgc:  0.1327641 
#> outer mgc:  0.1447738 
#> outer mgc:  0.1446607 
#> outer mgc:  0.04233245 
#> outer mgc:  0.04229455 
#> outer mgc:  0.02365667 
#> outer mgc:  0.02363525 
#> outer mgc:  0.02542886 
#> outer mgc:  0.02540674 
#> outer mgc:  12.07176 
#> outer mgc:  12.06473 
#> outer mgc:  4.760808 
#> outer mgc:  4.743822 
#> outer mgc:  1.16482 
#> outer mgc:  1.161071 
#> outer mgc:  4.327215 
#> outer mgc:  4.334191 
#> outer mgc:  1.343677 
#> outer mgc:  1.344357 
#> outer mgc:  0.3113822 
#> outer mgc:  0.312607 
#> outer mgc:  6.886272 
# }
```
