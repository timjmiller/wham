# Specify the age composition models for fleet(s) and indices.

Specify the age composition models for fleet(s) and indices.

## Usage

``` r
set_age_comp(input, age_comp)
```

## Arguments

- input:

  list containing data, parameters, map, and random elements (output
  from
  [`wham::prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md))

- age_comp:

  specifies the age composition models for fleet(s) and indices. If
  `NULL`, the multinomial is used because this was the only option in
  ASAP.

  The age composition models available are:

  `"multinomial"`

  :   Multinomial. This is the default because it was the only option in
      ASAP. 0 parameters.

  `"dir-mult"`

  :   Saturating Dirichlet-multinomial, parameterized such that
      effective-sample-size is a nonlinear and saturating function with
      respect to input-sample-size. 1 parameter. Effective sample size
      is estimated by the model ([Candy
      2008](https://www.ccamlr.org/es/system/files/science_journal_papers/07candy.pdf))

  `"dirichlet-pool0"`

  :   Dirichlet, pooling zero observations with adjacent age classes. 1.
      parameter. See [Francis
      2014](https://www.sciencedirect.com/science/article/abs/pii/S0165783613003093)
      and [Albertsen et al.
      2016](https://cdnsciencepub.com/doi/abs/10.1139/cjfas-2015-0532)

  `"dirichlet-miss0"`

  :   

  `"logistic-normal-miss0"`

  :   Logistic normal, treating zero observations as missing. 1
      parameter.

  `"logistic-normal-ar1-miss0"`

  :   Logistic normal, treating zero observations as missing. 1
      parameter.

  `"logistic-normal-pool0"`

  :   Logistic normal, pooling zero observations with adjacent age
      classes. 1 parameter. See [Schnute and Haigh
      (2007)](https://doi.org/10.1093/icesjms/fsl024) and [Francis
      (2014)](https://doi.org/10.1016/j.fishres.2013.12.015)

  `"logistic-normal-01-infl"`

  :   Zero-or-one inflated logistic normal. Inspired by zero-one
      inflated beta in [Ospina and Ferrari
      (2012)](https://www.sciencedirect.com/science/article/abs/pii/S0167947311003628).
      3 parameters. . No OSA residuals.

  `"logistic-normal-01-infl-2par"`

  :   Zero-one inflated logistic normal where p0 is a function of
      binomial sample size. 2 parameters. No OSA residuals.

  `"mvtweedie"`

  :   Multivariate-tweedie, where the product of composition proportions
      and input sample sizes follows a distribution with mean equal to
      the product of predicted proportions and input sample size, and
      other parameters define the ratio of effective to input sample
      size (with is bounded 0 to Inf) and the probability of zeros. 2
      parameters. No OSA residuals.

  `"dir-mult-linear"`

  :   Linear Dirichlet-multinomial, parameterized such that
      effective-sample-size is a linear function with respect to
      input-sample-size, estimating 1 parameter, \\log(\theta)\\, where
      the ratio of effective and input sample size is approximately
      \\\theta / (1+\theta)\\, i.e., the logistic transformation of the
      estimated parameter \\log(\theta)\\. ([Thorson et al.
      2017](https://doi.org/10.1016/j.fishres.2016.06.005))

  The two Dirichlet-multinomial options will only differ when
  input-sample-size differs among years. In these cases, the
  linear-Dirichlet multinomial is designed to decrease the effective
  sample size in each year by approximately the same proportion, while
  the saturating-Dirichlet multinomial will decrease the years with
  highest input-sample-size much more than those with lower
  input-sample-size. One-step-ahead residuals will be calculated for all
  but options 8-10 when `do.osa=TRUE` (Nielsen et al. in prep.). An age
  composition model needs to be specified for each fleet and index. If
  you would like all fleets and indices to use the same age composition
  likelihood, you can simply specify one of the strings above, i.e.
  `age_comp = "logistic-normal-miss0"`. If you do not want the same age
  composition model to be used for all fleets and indices, you must
  specify a named list with the following entries:

  \$fleets

  :   A vector of the above strings with length = the number of fleets.

  \$indices

  :   A vector of the above strings with length = the number of indices.

## Value

a named list with same elements as the input provided with age
composition likelihood options modified.

## See also

[`prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md)

## Examples

``` r
if (FALSE) { # \dontrun{
wham.dir <- find.package("wham")
path_to_examples <- system.file("extdata", package="wham")
asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
input <- prepare_wham_input(asap3)
input <- set_age_comp(input, age_comp = "logistic-normal-miss0") #no longer multinomial
} # }
```
