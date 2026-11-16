# Read an ASAP3 .dat file into R

WHAM was originally built to fit models very much like
[ASAP](https://nmfs-ost.github.io/noaa-fit/ASAP) ([Legault and Restrepo
(1999)](https://www.iccat.int/Documents/CVSP/CV049_1999/n_2/CV049020246.pdf))
and extend from there to include other features. This function provides
functionality to use a preexisting ASAP3 input data file. The output of
`read_asap3_dat` should then be passed to
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
# \donttest{
asap3 = read_asap3_dat(file.path(system.file("extdata", package="wham"), "ex1_SNEMAYT.dat"))
sinput = prepare_wham_input(asap3)
#> 
#> --Creating input---------------------------------------------------------------------------------------------------------------------
#> basic_info done
#> catch done
#> indices done
#> WAA done
#> NAA done
#> q done
#> selectivity done
#> age_comp done
#> F done
#> M done
#> move done
#> osa_obs done
#> --Input Complete--------------------------------------------------------------------------------------------------------------------
#> -------------------------------------------------------------------------------------------------------------------------------------
#> WHAM version 2.0.0 and forward: 
#> 1) makes major changes to the structure of some data, parameters, and reported objects.
#> 2) decouples random effects for recruitment and random effects for older ages by default. To obtain results from previous versions set NAA_re$decouple_recruitment = FALSE.
#> 3) does not bias correct any log-normal process or observation errors. To configure these, set basic_info$bias_correct_process = TRUE and/or basic_info$bias_correct_observation = TRUE.
#> -------------------------------------------------------------------------------------------------------------------------------------
#> 
#> --USING ASAP3 Inputs-----------------------------------------------------------------------------------------------------------------
#> 1 asap3 dat file was processed. A single stock and region will be assumed. 
#> -------------------------------------------------------------------------------------------------------------------------------------
#> 
#> --Catch------------------------------------------------------------------------------------------------------------------------------ 
#> waa_pointer_fleets determined from ASAP file(s). 
#> Selectivity blocks for each fleet:
#> fleet_1: 1
#> 
#> -------------------------------------------------------------------------------------------------------------------------------------
#> 
#> --Indices----------------------------------------------------------------------------------------------------------------------------
#> waa_pointer_indices determined from ASAP file(s). 
#> Selectivity blocks for each index:
#> index_1: 2
#> index_2: 3
#> 
#> -------------------------------------------------------------------------------------------------------------------------------------
#> 
#> --WAA--------------------------------------------------------------------------------------------------------------------------------
#> basic_info$waa_pointer_M was not provided, so waa_pointer_ssb will be used. 
#> -------------------------------------------------------------------------------------------------------------------------------------
#> 
#> --Selectivity------------------------------------------------------------------------------------------------------------------------
#> (Mean) selectivity block models are:
#> Block 1: logistic
#> Block 2: logistic
#> Block 3: logistic
#> 
#> Random effects options for each selectivity block are:
#> Block 1: none
#> Block 2: none
#> Block 3: none
#> 
#> -------------------------------------------------------------------------------------------------------------------------------------
#> 
#> --Reference points-------------------------------------------------------------------------------------------------------------------
#> For annual SPR-based reference points, corresponding annual conditionally expected recruitments are used. 
#> -------------------------------------------------------------------------------------------------------------------------------------
#> 
mod <- fit_wham(input, do.osa = FALSE, do.retro = FALSE)
#> Error: object 'input' not found
# }
```
