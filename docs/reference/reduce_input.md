# Reduce the years of the model

Internal function called by
[`fit_peel`](https://timjmiller.github.io/wham/reference/fit_peel.md)
for *i* in 1–`n.peels`. Creates the input for the model peeling off *i*
years of data (calls
[`fit_tmb`](https://timjmiller.github.io/wham/reference/fit_tmb.md)).
Reduces the year dimension of all data, parameters, maps, so that one
can project correctly from the peeled model.

## Usage

``` r
reduce_input(input, years_peeled, retro = TRUE)
```

## Arguments

- input:

  list containing data, parameters, map, and random elements (output
  from
  [`prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md)).
  NOT from prepare_projection or project_wham

- years_peeled:

  which of input\$years to peel from the model input

- retro:

  (T/F) whether this is for a retro peel (Default = TRUE)
