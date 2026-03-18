# Specify configuration for environmental covariates, effects on the population, and parameter values

Specify configuration for environmental covariates, effects on the
population, and parameter values

## Usage

``` r
set_ecov(input, ecov)
```

## Arguments

- input:

  list containing data, parameters, map, and random elements (output
  from
  [`wham::prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md))

- ecov:

  (optional) named list of environmental covariate data and parameters
  (see details)

  `ecov` specifies any environmental covariate data and model as well as
  the effect on the population. Environmental covariate data need not
  span the same years as the fisheries data. It can be `NULL` if no
  environmental data are to be fit. Otherwise, it must be a named list
  with the following components:

  \$label

  :   Name(s) of the environmental covariate(s). Used in printing.

  \$mean

  :   Mean observations (matrix). number of years x number of
      covariates. Missing values = NA.

  \$logsigma

  :   Configure observation standard errors. Options:

      Matrix of \\log\\ standard errors with same dimensions as `$mean`

      :   Specified values for each time step

      log standard errors for each covariate, numeric vector or matrix w/ dim 1 x n.ecov

      :   Specified value the same for all time steps

      estimation option (for all covariates). character string:

      :   `"est_1"`: Estimated, one value shared among time steps.
          `"est_re"`: Estimated value for each time step as random
          effects with two parameters (mean, var)

      list of two elements.

      :   First is the matrix of log standard errors or the vector of
          single values for each covariate as above. Second is a
          character vector of estimation options (`NA`,
          `"est_1"`,`"est_re"`) for each covariate. For covariates with
          non-NA values, values in the first element are ignored.

  \$year

  :   Years corresponding to observations (vector of same length as
      `$mean` and `$logsigma`)

  \$use_obs

  :   T/F (or 1/0) vector/matrix of the same dimension as `$mean` and
      `$logsigma`. Use the observation? Can be used to ignore subsets of
      the ecov without changing data files.

  \$process_model

  :   Process model for the ecov time-series. `"rw"` = random walk,
      `"ar1"` = 1st order autoregressive, `NA` = do not fit

  \$process_mean_vals

  :   vector of (initial) mean values for the ecov time-series.

  \$process_sig_vals

  :   vector of (initial) standard deviation values for the ecov
      time-series.

  \$process_cor_vals

  :   vector of (initial) correlation values for the ecov time-series.

  \$recruitment_how

  :   character matrix (n_Ecov x n_stocks) indicating how each ecov
      affects recruitment for each stock. Options are based on (see
      [Iles & Beverton
      (1998)](https://www.sciencedirect.com/science/article/pii/S1385110197000221))
      combined with the order of orthogonal polynomial of the covariate
      and has the form "type-lag-order". "type" can be:

      = "none"

      :   no effect.

      = "controlling"

      :   pre-recruit density-independent mortality.

      = "limiting"

      :   maximum recruitment, e.g. ecov determines amount of suitable
          habitat)

      = "lethal"

      :   threshold, i.e. R –\> 0 at some ecov value.

      = "masking"

      :   metabolic/growth, decreases dR/dS

      = "directive"

      :   e.g. behavioral

      for type other than "none", "lag" can be:

      = "lag-n"

      :   lag = n which can be 0,1,2,.... lag-1 implies the covariate in
          year y affects recruitment in year y+1.

      for "type" being other than "none", "order" can be:

      = "linear"

      :   the covariate effect is linear on the transformed recruitment
          parameter (e.g., log).

      = "poly-n"

      :   orthogonal polynomial where n = 1 (same as "linear"),2,...

      so "limiting-lag-1-poly-2" would model the covariate affecting
      recruitment the next year (lag = 1) as a second order orthogonal
      polynomial (\\b_0 + b_1\*ecov + b_2\*ecov^2 + ...\\) limiting
      effect.

  \$M_how

  :   character array (n_Ecov x n_stocks x n_ages x n_regions)
      indicating how each ecov affects M by age,stock,region and has the
      form "lag-order". "lag" can be:

      = "none"

      :   no effect.

      = "lag-n"

      :   lag = n which can be 0,1,2,.... lag-1 implies the covariate in
          year y affects M in year y+1.

      for "lag" being other than "none", "order" can be:

      = "linear"

      :   the covariate effect is linear on the transformed M parameter
          (e.g., log).

      = "poly-n"

      :   orthogonal polynomial where n = 1 (same as "linear"),2,...

  \$M_effect_map

  :   integer array (n_stocks x n_ages x n_regions x n_Ecov) indicating
      which estimated effects are common by age,stock,region. If not
      specified there the same effect is estimated for all M where
      \$M_how is other than "none" for each covariate.

  \$q_how

  :   character matrix (n_Ecov x n_indices) indicating whether each ecov
      affects catchability for each index. and has the form "lag-order".
      "lag" can be:

      = "none"

      :   no effect.

      = "lag-n"

      :   lag = n which can be 0,1,2,.... lag-1 implies the covariate in
          year y affects catchability in year y+1.

      for "lag" being other than "none", "order" can be:

      = "linear"

      :   the covariate effect is linear on the transformed catchability
          parameter (e.g., log).

      = "poly-n"

      :   orthogonal polynomial where n = 1 (same as "linear"),2,...

  \$move_how

  :   character array (n_Ecov x n_stocks x n_ages x n_seasons x
      n_regions x n_regions - 1) indicating whether each ecov affects
      movement from one region to the others by stock,age,season. and
      has the form "lag-order". "lag" can be:

      = "none"

      :   no effect.

      = "lag-n"

      :   lag = n which can be 0,1,2,.... lag-1 implies the covariate in
          year y affects a movement parameter in year y+1.

      for "lag" being other than "none", "order" can be:

      = "linear"

      :   the covariate effect is linear on the transformed movement
          parameter (e.g., log).

      = "poly-n"

      :   orthogonal polynomial where n = 1 (same as "linear"),2,...

  \$move_effect_map

  :   integer array (n_stocks x n_ages x n_seasons x n_regions x
      n_regions-1 x n_Ecov) indicating which estimated effects are
      common by age,stock,region, season etc. If not specified the same
      effect is estimated for all movement parameters where \$move_how
      is other than "none" for each covariate.

  \$beta_R_vals

  :   n_stocks x n_ecov x max(n_poly_R) array of initial values for
      effects on recruitment.

  \$beta_M_vals

  :   n_stocks x n_ages x n_regions x n_ecov x max(n_poly_M) array of
      initial values for effects on natural mortality.

  \$beta_q_vals

  :   n_indices x n_ecov x max(n_poly_q) array of initial values for
      effects on catchability.

  \$beta_mu_vals

  :   n_stocks x n_ages x n_seasons x n_regions x n_regions - 1 x n_ecov
      x max(n_poly_move) array of initial values for effects on movement
      parameters.

## Value

a named list with same elements as the input provided with environmental
covariate observations, effects, and model options modified.

## See also

[`prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md)

## Examples

``` r
if (FALSE) { # \dontrun{
wham.dir <- find.package("wham")
path_to_examples <- system.file("extdata", package="wham")
asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
env.dat <- read.csv(file.path(path_to_examples,"GSI.csv"), header=T)
input <- prepare_wham_input(asap3, NAA_re = list(sigma = "rec"))
ecov <- list(
 label = "GSI",
 mean = as.matrix(env.dat$GSI),
 logsigma = 'est_1', # estimate obs sigma, 1 value shared across years
 year = env.dat$year,
 use_obs = matrix(1, ncol=1, nrow=dim(env.dat)[1]), # use all obs (=1)
 process_model = 'ar1') # "rw" or "ar1"
input <- set_ecov(input, ecov = ecov) #GSI in the model without any effects
} # }
```
