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
if (FALSE) { # \dontrun{
data("input4_SNEMAYT") # load SNEMA yellowtail flounder data and parameter settings
mod <- fit_wham(input4_SNEMAYT) # using default values
check_convergence(mod)
} # }
```
