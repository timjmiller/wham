# Read an ASAP3 .dat file into R

WHAM is built on ASAP (\[Legault and Restrepo
(1999)\](http://sedarweb.org/docs/wsupp/S12RD06%20ASAPdoc.pdf)) and this
function provides functionality to use a preexisting ASAP3 input data
file. The output of `read_asap3_dat` should then be passed to
[`prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md).
If you are not familiar with ASAP3 input files, see the ASAP
[documentation](https://github.com/cmlegault/ASAPplots/tree/master/pdf)
and [code](https://nmfs-fish-tools.github.io/ASAP/).

## Usage

``` r
read_asap3_dat(filename)
```

## Arguments

- filename:

  character vector, names of ASAP3 .dat files. The file either needs to
  be in the current working directory, or `filename` can include the
  path. If multipile files, a multi-stock model will be assumed.

## Value

a named list with the following components:

- `dat`:

  Named list of input data and parameters

- `comments`:

  Comments at top of ASAP3 .dat file (indicated by "#")

## See also

[`prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md),
[`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md),
[ASAP
documentation](https://github.com/cmlegault/ASAPplots/tree/master/pdf)

## Examples

``` r
if (FALSE) { # \dontrun{
asap3 = read_asap3_dat("ASAP_SNEMAYT.dat")
input = prepare_wham_input(asap3)
mod = fit_wham(input)
} # }
```
