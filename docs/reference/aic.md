# Calculate Akaike's Information Criterion for a WHAM model

Calculates AIC from a fitted WHAM model's negative log-likelihood and
number of estimated parameters.

## Usage

``` r
aic(mod, conditional = FALSE)
```

## Arguments

- mod:

  A fitted WHAM model object returned by
  [`fit_wham()`](https://timjmiller.github.io/wham/reference/fit_wham.md).

- conditional:

  (TRUE/FALSE) When the model includes random effects, the default
  (`conditional = FALSE`) calculation uses the marginal likelihood and
  number of fixed effects parameters. If `conditional = TRUE`, the joint
  likelihood of the data conditional on the estimated random effects is
  used with an estimated effective degress of freedom that is calculated
  using the approach descsribed by [Zhang et al.
  2024](https://doi.org/10.48550/arXiv.2411.14185) and code provided by
  Noel Cadigan.

## Value

A numeric AIC value with attributes denoting the degrees of freedom,
number of observations, and the type (marginal or conditional).

## Examples

``` r
if (FALSE) { # \dontrun{
mod <- fit_wham(input)
aic(mod)
aic(mod, conditional = TRUE)
} # }
```
