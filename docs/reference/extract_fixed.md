# Extract fixed effects Originally provided by the [TMBhelper](https://github.com/kaskr/TMB_contrib_R/tree/master/TMBhelper) package. Internal function called by [`check_estimability`](https://timjmiller.github.io/wham/reference/check_estimability.md).

`extract_fixed` extracts the best previous value of fixed effects, in a
way that works for both mixed and fixed effect models

## Usage

``` r
extract_fixed(obj)
```

## Arguments

- obj, :

  The compiled object

## Value

A vector of fixed-effect estimates
