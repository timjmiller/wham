# Check for identifiability of fixed effects Originally provided by the [TMBhelper](https://github.com/kaskr/TMB_contrib_R/tree/master/TMBhelper) package. Internal function called by [`fit_tmb`](https://timjmiller.github.io/wham/reference/fit_tmb.md).

`check_estimability` calculates the matrix of second-derivatives of the
marginal likelihood w.r.t. fixed effects, to see if any linear
combinations are not estimable (i.e. cannot be uniquely estimated
conditional upon model structure and available data, e.g., resulting in
a likelihood ridge and singular, non-invertable Hessian matrix)

## Usage

``` r
check_estimability(obj, h)
```

## Arguments

- obj:

  The compiled object

- h:

  optional argument containing pre-computed Hessian matrix

## Value

A tagged list of the hessian and the message
