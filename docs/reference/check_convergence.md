# Check convergence of WHAM model

Access quick convergence checks from \`TMB\` and \`nlminb\`.

## Usage

``` r
check_convergence(mod, ret = FALSE)
```

## Arguments

- mod:

  output from
  [`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md)

- ret:

  T/F, return list? Default = FALSE, just prints to console

## Value

a list with at least the first three of these components:

- `$convergence`:

  From [`stats::nlminb`](https://rdrr.io/r/stats/nlminb.html), "0
  indicates successful convergence for nlminb"

- `$maxgr`:

  Max absolute gradient value, from \`max(abs(mod\$gr(mod\$opt\$par)))\`

- `$maxgr_par`:

  Name of parameter with max gradient

- `$is_sdrep`:

  If [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html) was
  performed for this model, this indicates whether it performed without
  error

- `$na_sdrep`:

  If [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html) was
  performed without error for this model, this indicates which (if any)
  components of the diagonal of the inverted hessian were returned as NA

## See also

[`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md),
[`fit_tmb`](https://timjmiller.github.io/wham/reference/fit_tmb.md),
[`stats::nlminb`](https://rdrr.io/r/stats/nlminb.html)

## Examples

``` r
# \donttest{
input <- readRDS(file.path(system.file("extdata", package="wham"), "ytf_rec1.RDS"))
#> Warning: cannot open compressed file 'C:/Users/timothy.j.miller/AppData/Local/Programs/R/R-4.4.3/library/wham/extdata/ytf_rec1.RDS', probable reason 'No such file or directory'
#> Error in gzfile(file, "rb"): cannot open the connection
mod <- fit_wham(input, do.osa = FALSE, do.retro = FALSE)
#> Error: object 'input' not found
check_convergence(mod)
#> Error: object 'mod' not found
# }
```
