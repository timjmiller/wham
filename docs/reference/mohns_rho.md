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
# \donttest{
input <- readRDS(file.path(system.file("extdata", package="wham"), "ytf_rec1.RDS"))
#> Warning: cannot open compressed file 'C:/Users/timothy.j.miller/AppData/Local/Programs/R/R-4.4.3/library/wham/extdata/ytf_rec1.RDS', probable reason 'No such file or directory'
#> Error in gzfile(file, "rb"): cannot open the connection
mod <- fit_wham(input, do.osa = FALSE, do.retro = TRUE)
#> Error: object 'input' not found
mohns_rho(mod) # calculate Mohn's rho
#> Error: object 'mod' not found
# }
```
