# Calculate Mohn's rho for a WHAM model with peels

Calculate Mohn's rho for a WHAM model with peels

## Usage

``` r
mohns_rho(model)
```

## Arguments

- model:

  A fit WHAM model, output from
  [`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md)
  with `do.retro = TRUE`.

## Value

`rho`, a vector of Mohn's rho

## See also

[`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md),
[`retro`](https://timjmiller.github.io/wham/reference/retro.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data("input4_SNEMAYT") # load SNEMA yellowtail flounder data and parameter settings
mod <- fit_wham(input4_SNEMAYT) # using default values: do.retro = T, n.peels = 7
mohns_rho(mod) # calculate Mohn's rho
} # }
```
