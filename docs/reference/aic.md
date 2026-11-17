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
# \donttest{
input <- readRDS(file.path(system.file("extdata", package="wham"), "ytf_rec1.RDS"))
#> Warning: cannot open compressed file 'C:/Users/timothy.j.miller/AppData/Local/Programs/R/R-4.4.3/library/wham/extdata/ytf_rec1.RDS', probable reason 'No such file or directory'
#> Error in gzfile(file, "rb"): cannot open the connection
mod <- fit_wham(input, do.osa = FALSE, do.retro = FALSE)
#> Error: object 'input' not found
aic(mod)
#> Error: object 'mod' not found
aic(mod, conditional = TRUE)
#> Error: object 'mod' not found
# }
```
