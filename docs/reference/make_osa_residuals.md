# Calculate one-step-ahead residuals

Standard residuals are not appropriate for models with random effects.
Instead, one-step-ahead (OSA) residuals can be used for evaluating model
goodness-of-fit ([Thygeson et al.
(2017)](https://doi.org/10.1007/s10651-017-0372-4), implemented in
[`TMB::oneStepPredict`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)).
OSA residual options are passed to
[`TMB::oneStepPredict`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
in a list `osa.opts`. Current options are method: oneStepGaussianOffMode
(default), oneStepGaussian, or oneStepGeneric, and parallel: TRUE/FALSE.
See
[`TMB::oneStepPredict`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
for further details. It is not recommended to run this function (or
[`TMB::oneStepPredict`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html))
with any random effects and mvtweedie age composition likelihoods due to
extensive computational demand. An error will be thrown in such cases.
See [Trijoulet et al.
(2023)](https://doi.org/10.1016/j.fishres.2022.106487) for OSA methods
for age composition OSA residuals.

## Usage

``` r
make_osa_residuals(
  model,
  osa.opts = list(method = "oneStepGaussianOffMode", parallel = TRUE),
  sdrep_required = TRUE
)
```

## Arguments

- model:

  A fit WHAM model, output from
  [`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md).

- osa.opts:

  a list of options: `method` = ("oneStepGaussianOffMode" (default),
  "oneStepGaussian", or "oneStepGeneric") and `parallel` (T (default)/F)
  to pass to
  [`TMB::oneStepPredict`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html).

- sdrep_required:

  (T/F) whether the osa residuals should be calculated without the
  [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html) being
  performed on the model.

## Value

the same fit TMB model with additional elements for osa residuals:

- `$OSA.Ecov`:

  data.frame returned by
  [`TMB::oneStepPredict`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  for environmental observations, if applicable.

- `$OSA.agregate`:

  data.frame returned by
  [`TMB::oneStepPredict`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  for aggregate catch and index observations conditional on any
  environmental observations, if applicable.

- `$OSA.agecomp`:

  data.frame returned by
  [`TMB::oneStepPredict`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  for age composition observations conditional on any aggregate catch or
  index, or environmental observations, if applicable.

- `$osa`:

  One-step-ahead residuals (if `do.osa=TRUE`)

## See also

[`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.md)

## Examples

``` r
# \donttest{
input <- readRDS(file.path(system.file("extdata", package="wham"), "ytf_rec1.RDS"))
mod <- fit_wham(input, do.osa = FALSE, do.retro = FALSE)
#> Order of parameters:
#>  [1] "log_catch_sig_scale"  "log_index_sig_scale"  "log_N1"              
#>  [4] "N1_repars"            "log_NAA_sigma"        "trans_NAA_rho"       
#>  [7] "log_NAA"              "mean_rec_pars"        "logit_q"             
#> [10] "q_re"                 "q_repars"             "q_prior_re"          
#> [13] "logit_selpars"        "selpars_re"           "sel_repars"          
#> [16] "catch_paa_pars"       "index_paa_pars"       "F_pars"              
#> [19] "Mpars"                "M_re"                 "M_repars"            
#> [22] "log_b"                "logR_proj"            "mu_prior_re"         
#> [25] "trans_mu"             "mu_re"                "mu_repars"           
#> [28] "L_re"                 "L_repars"             "Ecov_obs_logsigma"   
#> [31] "Ecov_obs_logsigma_re" "Ecov_obs_sigma_par"   "Ecov_process_pars"   
#> [34] "Ecov_re"              "Ecov_beta_R"          "Ecov_beta_q"         
#> [37] "Ecov_beta_M"          "Ecov_beta_mu"        
#> Not matching template order:
#>  [1] "mean_rec_pars"        "logit_q"              "q_prior_re"          
#>  [4] "q_re"                 "q_repars"             "F_pars"              
#>  [7] "mu_prior_re"          "trans_mu"             "mu_re"               
#> [10] "mu_repars"            "N1_repars"            "log_N1"              
#> [13] "log_NAA_sigma"        "trans_NAA_rho"        "log_NAA"             
#> [16] "logR_proj"            "logit_selpars"        "selpars_re"          
#> [19] "sel_repars"           "catch_paa_pars"       "index_paa_pars"      
#> [22] "log_catch_sig_scale"  "log_index_sig_scale"  "Mpars"               
#> [25] "M_re"                 "M_repars"             "log_b"               
#> [28] "L_re"                 "L_repars"             "Ecov_re"             
#> [31] "Ecov_beta_R"          "Ecov_beta_M"          "Ecov_beta_mu"        
#> [34] "Ecov_beta_q"          "Ecov_process_pars"    "Ecov_obs_logsigma"   
#> [37] "Ecov_obs_logsigma_re" "Ecov_obs_sigma_par"  
#> Your parameter list has been re-ordered.
#> (Disable this warning with checkParameterOrder=FALSE)
#> iter: 1  value: 7882.877 mgc: 159.82 ustep: 0.01720205 
#> iter: 2  value: 4915.778 mgc: 95.24501 ustep: 0.06792212 
#> iter: 3  value: 4319.932 mgc: 86.22232 ustep: 0.2606927 
#> iter: 4  value: 4231.836 mgc: 10.41171 ustep: 0.5106296 
#> iter: 5  value: 4218.468 mgc: 2.827659 ustep: 0.7146121 
#> iter: 6  value: 4217.41 mgc: 0.7701239 ustep: 0.8453628 
#> iter: 7  value: 4217.387 mgc: 0.1044678 ustep: 0.9194442 
#> iter: 8  value: 4217.387 mgc: 0.008092813 ustep: 0.9588806 
#> iter: 9  value: 4217.387 mgc: 0.0003825941 ustep: 0.9792266 
#> iter: 10  value: 4217.387 mgc: 8.535063e-06 ustep: 0.9895598 
#> iter: 11  value: 4217.387 mgc: 9.355553e-08 ustep: 0.9947667 
#> iter: 12  mgc: 5.097112e-10 
#> iter: 1  mgc: 5.097112e-10 
#> Using stats::nlminb for optimization with these control parameters:
#> iter.max: 1000
#> eval.max: 1000
#> iter: 1  mgc: 5.097112e-10 
#> iter: 1  mgc: 5.097112e-10 
#> Matching hessian patterns... Done
#> outer mgc:  838.4045 
#> iter: 1  value: 23056.44 mgc: 263.4276 ustep: 0.006061026 
#> iter: 2  value: 18235.24 mgc: 221.5431 ustep: 0.00873895 
#> iter: 3  value: 15787.24 mgc: 130.5354 ustep: 0.01630886 
#> iter: 4  value: 13747.11 mgc: 65.39022 ustep: 0.03359784 
#> iter: 5  value: 12881.8 mgc: 42.39166 ustep: 0.09631942 
#> iter: 6  value: 11861.02 mgc: 354.3893 ustep: 0.08871461 
#> iter: 7  value: 11586.84 mgc: 94.41374 ustep: 0.1620818 
#> iter: 8  value: 11472.83 mgc: 21.00957 ustep: 0.2270954 
#> iter: 9  value: 11419.62 mgc: 43.31709 ustep: 0.149401 
#> iter: 10  value: 11396.08 mgc: 11.64454 ustep: 0.114986 
#> iter: 11  value: 11381.67 mgc: 7.618771 ustep: 0.09835932 
#> iter: 12  value: 11365.44 mgc: 8.602592 ustep: 0.3136915 
#> iter: 13  value: 11362.96 mgc: 2.142472 ustep: 0.3367361 
#> iter: 14  value: 11350.85 mgc: 21.84388 ustep: 0.5803317 
#> iter: 15  value: 11349.35 mgc: 10.00032 ustep: 0.7618189 
#> iter: 16  value: 11349.2 mgc: 1.374198 ustep: 0.8728351 
#> iter: 17  value: 11349.18 mgc: 0.4708787 ustep: 0.934263 
#> iter: 18  value: 11349.18 mgc: 0.02871575 ustep: 0.9665762 
#> iter: 19  value: 11349.18 mgc: 0.001317889 ustep: 0.9831477 
#> iter: 20  value: 11349.18 mgc: 3.459377e-05 ustep: 0.9915389 
#> iter: 21  value: 11349.18 mgc: 4.982199e-07 ustep: 0.9957609 
#> iter: 22  mgc: 3.663683e-09 
#> iter: 1  value: 3978.725 mgc: 43.15786 ustep: 1 
#> iter: 2  value: 3923.777 mgc: 53.51634 ustep: 1 
#> iter: 3  value: 3914.555 mgc: 16.92558 ustep: 1 
#> iter: 4  value: 3912.124 mgc: 6.327685 ustep: 1 
#> iter: 5  value: 3911.738 mgc: 1.869159 ustep: 1 
#> iter: 6  value: 3911.724 mgc: 0.307486 ustep: 1 
#> iter: 7  value: 3911.724 mgc: 0.01167152 ustep: 1 
#> iter: 8  value: 3911.724 mgc: 1.81751e-05 ustep: 1 
#> iter: 9  mgc: 4.422622e-11 
#> iter: 1  mgc: 4.422622e-11 
#> outer mgc:  749.5992 
#> iter: 1  value: 3737.958 mgc: 35.62372 ustep: 1 
#> iter: 2  value: 3737.721 mgc: 0.8122296 ustep: 1 
#> iter: 3  value: 3737.72 mgc: 0.01599528 ustep: 1 
#> iter: 4  value: 3737.72 mgc: 2.332566e-05 ustep: 1 
#> iter: 5  mgc: 4.983219e-11 
#> iter: 1  mgc: 4.983219e-11 
#> outer mgc:  584.0328 
#> iter: 1  value: 4266.603 mgc: 58.36646 ustep: 0.4367355 
#> iter: 2  value: 3691.738 mgc: 302.4568 ustep: 0.6608936 
#> iter: 3  value: 3630.969 mgc: 60.79021 ustep: 0.8129723 
#> iter: 4  value: 3624.955 mgc: 14.30922 ustep: 0.9016596 
#> iter: 5  value: 3624.711 mgc: 2.744502 ustep: 0.9495626 
#> iter: 6  value: 3624.71 mgc: 0.1620449 ustep: 0.9744576 
#> iter: 7  value: 3624.71 mgc: 0.0008345152 ustep: 0.9871475 
#> iter: 8  value: 3624.71 mgc: 4.905308e-07 ustep: 0.9935536 
#> iter: 9  mgc: 2.359464e-10 
#> iter: 1  mgc: 2.359464e-10 
#> outer mgc:  860.2559 
#> iter: 1  value: 3286.005 mgc: 56.07818 ustep: 1 
#> iter: 2  value: 3283.968 mgc: 2.934293 ustep: 1 
#> iter: 3  value: 3283.965 mgc: 0.08807975 ustep: 1 
#> iter: 4  value: 3283.965 mgc: 0.0005970503 ustep: 1 
#> iter: 5  value: 3283.965 mgc: 4.78873e-08 ustep: 1 
#> iter: 6  mgc: 7.597994e-14 
#> iter: 1  mgc: 7.597994e-14 
#> outer mgc:  302.9141 
#> iter: 1  value: 3128.213 mgc: 16.21896 ustep: 1 
#> iter: 2  value: 3128.009 mgc: 0.8443066 ustep: 1 
#> iter: 3  value: 3128.009 mgc: 0.0290793 ustep: 1 
#> iter: 4  value: 3128.009 mgc: 0.0001255467 ustep: 1 
#> iter: 5  mgc: 2.253581e-09 
#> iter: 1  mgc: 2.253581e-09 
#> outer mgc:  338.9048 
#> iter: 1  value: 2927.787 mgc: 54.81738 ustep: 1 
#> iter: 2  value: 2927.521 mgc: 1.391304 ustep: 1 
#> iter: 3  value: 2927.52 mgc: 0.02409551 ustep: 1 
#> iter: 4  value: 2927.52 mgc: 2.35672e-05 ustep: 1 
#> iter: 5  mgc: 2.947366e-11 
#> iter: 1  mgc: 2.947366e-11 
#> outer mgc:  344.0037 
#> iter: 1  value: 2721.669 mgc: 24.57158 ustep: 1 
#> iter: 2  value: 2721.22 mgc: 1.38444 ustep: 1 
#> iter: 3  value: 2721.219 mgc: 0.03909399 ustep: 1 
#> iter: 4  value: 2721.219 mgc: 0.0002287833 ustep: 1 
#> iter: 5  mgc: 7.557662e-09 
#> iter: 1  mgc: 7.557662e-09 
#> outer mgc:  437.3271 
#> iter: 1  value: 2717.633 mgc: 64.2714 ustep: 1 
#> iter: 2  value: 2714.479 mgc: 3.580022 ustep: 1 
#> iter: 3  value: 2714.475 mgc: 0.1198489 ustep: 1 
#> iter: 4  value: 2714.475 mgc: 0.0005246523 ustep: 1 
#> iter: 5  value: 2714.475 mgc: 1.544549e-08 ustep: 1 
#> mgc: 1.074597e-13 
#> iter: 1  mgc: 1.074597e-13 
#> outer mgc:  572.6003 
#> iter: 1  value: 2461.609 mgc: 31.1633 ustep: 1 
#> iter: 2  value: 2456.623 mgc: 3.531759 ustep: 1 
#> iter: 3  value: 2456.525 mgc: 0.7551975 ustep: 1 
#> iter: 4  value: 2456.525 mgc: 0.06118095 ustep: 1 
#> iter: 5  value: 2456.525 mgc: 0.0004533908 ustep: 1 
#> iter: 6  value: 2456.525 mgc: 2.49853e-08 ustep: 1 
#> mgc: 1.296254e-13 
#> iter: 1  mgc: 1.296254e-13 
#> outer mgc:  226.0325 
#> iter: 1  value: 2379.718 mgc: 29.1523 ustep: 1 
#> iter: 2  value: 2379.646 mgc: 0.7022036 ustep: 1 
#> iter: 3  value: 2379.646 mgc: 0.006580331 ustep: 1 
#> iter: 4  value: 2379.646 mgc: 1.573864e-06 ustep: 1 
#> iter: 5  mgc: 1.616882e-13 
#> iter: 1  mgc: 1.616882e-13 
#> outer mgc:  192.4525 
#> iter: 1  value: 2254.159 mgc: 16.60862 ustep: 1 
#> iter: 2  value: 2254.127 mgc: 0.9175179 ustep: 1 
#> iter: 3  value: 2254.127 mgc: 0.002398864 ustep: 1 
#> iter: 4  value: 2254.127 mgc: 6.430957e-07 ustep: 1 
#> iter: 5  mgc: 1.056046e-13 
#> iter: 1  mgc: 1.056046e-13 
#> outer mgc:  228.5998 
#> iter: 1  value: 2162.375 mgc: 36.33658 ustep: 1 
#> iter: 2  value: 2162.159 mgc: 0.9327455 ustep: 1 
#> iter: 3  value: 2162.159 mgc: 0.01199258 ustep: 1 
#> iter: 4  value: 2162.159 mgc: 1.669359e-05 ustep: 1 
#> iter: 5  mgc: 2.99184e-11 
#> iter: 1  mgc: 2.99184e-11 
#> outer mgc:  152.9878 
#> iter: 1  value: 2762.865 mgc: 61.9729 ustep: 1 
#> iter: 2  value: 2298.677 mgc: 380.5837 ustep: 1 
#> iter: 3  value: 2257.922 mgc: 72.55175 ustep: 1 
#> iter: 4  value: 2247.118 mgc: 23.77544 ustep: 1 
#> iter: 5  value: 2243.864 mgc: 8.798231 ustep: 1 
#> iter: 6  value: 2243.243 mgc: 2.818696 ustep: 1 
#> iter: 7  value: 2243.211 mgc: 0.5499186 ustep: 1 
#> iter: 8  value: 2243.211 mgc: 0.03011872 ustep: 1 
#> iter: 9  value: 2243.211 mgc: 9.980867e-05 ustep: 1 
#> iter: 10  mgc: 1.102562e-09 
#> iter: 1  value: 2121.279 mgc: 24.59156 ustep: 1 
#> iter: 2  value: 2120.075 mgc: 2.585009 ustep: 1 
#> iter: 3  value: 2120.04 mgc: 0.5259633 ustep: 1 
#> iter: 4  value: 2120.039 mgc: 0.02857164 ustep: 1 
#> iter: 5  value: 2120.039 mgc: 8.959776e-05 ustep: 1 
#> iter: 6  mgc: 8.753081e-10 
#> iter: 1  mgc: 8.753081e-10 
#> outer mgc:  229.6082 
#> iter: 1  value: 2143.844 mgc: 48.71411 ustep: 1 
#> iter: 2  value: 2142.826 mgc: 2.409764 ustep: 1 
#> iter: 3  value: 2142.826 mgc: 0.03355966 ustep: 1 
#> iter: 4  value: 2142.826 mgc: 2.57527e-05 ustep: 1 
#> iter: 5  mgc: 5.132127e-11 
#> iter: 1  value: 2129.353 mgc: 32.21959 ustep: 1 
#> iter: 2  value: 2128.924 mgc: 1.294651 ustep: 1 
#> iter: 3  value: 2128.924 mgc: 0.01367994 ustep: 1 
#> iter: 4  value: 2128.924 mgc: 2.127558e-05 ustep: 1 
#> iter: 5  mgc: 4.752264e-11 
#> iter: 1  value: 2131.669 mgc: 25.88899 ustep: 1 
#> iter: 2  value: 2131.402 mgc: 0.8097926 ustep: 1 
#> iter: 3  value: 2131.402 mgc: 0.01443238 ustep: 1 
#> iter: 4  value: 2131.402 mgc: 2.353806e-05 ustep: 1 
#> iter: 5  mgc: 5.807028e-11 
#> iter: 1  value: 2098.876 mgc: 12.45824 ustep: 1 
#> iter: 2  value: 2098.856 mgc: 0.2154628 ustep: 1 
#> iter: 3  value: 2098.856 mgc: 0.001670923 ustep: 1 
#> iter: 4  value: 2098.856 mgc: 3.225386e-07 ustep: 1 
#> iter: 5  mgc: 1.500958e-13 
#> iter: 1  mgc: 1.500958e-13 
#> outer mgc:  129.9969 
#> iter: 1  value: 2090.604 mgc: 6.42427 ustep: 1 
#> iter: 2  value: 2090.602 mgc: 0.1000729 ustep: 1 
#> iter: 3  value: 2090.602 mgc: 0.000144679 ustep: 1 
#> iter: 4  mgc: 1.186227e-09 
#> iter: 1  mgc: 1.186227e-09 
#> outer mgc:  110.1973 
#> iter: 1  value: 2073.914 mgc: 11.16234 ustep: 1 
#> iter: 2  value: 2073.91 mgc: 0.1087578 ustep: 1 
#> iter: 3  value: 2073.91 mgc: 0.0003175944 ustep: 1 
#> iter: 4  value: 2073.91 mgc: 1.170354e-08 ustep: 1 
#> mgc: 1.047497e-13 
#> iter: 1  mgc: 1.047497e-13 
#> outer mgc:  96.10531 
#> iter: 1  value: 2053.169 mgc: 2.719103 ustep: 1 
#> iter: 2  value: 2053.169 mgc: 0.03703428 ustep: 1 
#> iter: 3  value: 2053.169 mgc: 8.407329e-05 ustep: 1 
#> iter: 4  mgc: 8.131854e-10 
#> iter: 1  mgc: 8.131854e-10 
#> outer mgc:  106.9058 
#> iter: 1  value: 2036.572 mgc: 14.24249 ustep: 1 
#> iter: 2  value: 2036.563 mgc: 0.1840021 ustep: 1 
#> iter: 3  value: 2036.563 mgc: 0.0005497474 ustep: 1 
#> iter: 4  value: 2036.563 mgc: 1.636641e-08 ustep: 1 
#> mgc: 1.574266e-13 
#> iter: 1  mgc: 1.574266e-13 
#> outer mgc:  84.62055 
#> iter: 1  value: 2019.159 mgc: 8.954107 ustep: 1 
#> iter: 2  value: 2019.159 mgc: 0.02503877 ustep: 1 
#> iter: 3  value: 2019.159 mgc: 1.11622e-05 ustep: 1 
#> iter: 4  mgc: 1.974742e-12 
#> iter: 1  mgc: 1.974742e-12 
#> outer mgc:  62.11866 
#> iter: 1  value: 1994.336 mgc: 17.2447 ustep: 1 
#> iter: 2  value: 1994.32 mgc: 0.2536224 ustep: 1 
#> iter: 3  value: 1994.32 mgc: 0.0004679887 ustep: 1 
#> iter: 4  value: 1994.32 mgc: 1.305022e-08 ustep: 1 
#> mgc: 1.644872e-13 
#> iter: 1  value: 1976.665 mgc: 11.97185 ustep: 1 
#> iter: 2  value: 1976.662 mgc: 0.1033863 ustep: 1 
#> iter: 3  value: 1976.662 mgc: 0.0001516678 ustep: 1 
#> iter: 4  mgc: 3.629762e-09 
#> iter: 1  mgc: 3.629762e-09 
#> outer mgc:  88.37673 
#> iter: 1  value: 1944.302 mgc: 58.14252 ustep: 1 
#> iter: 2  value: 1943.757 mgc: 1.527867 ustep: 1 
#> iter: 3  value: 1943.757 mgc: 0.03781426 ustep: 1 
#> iter: 4  value: 1943.757 mgc: 5.057419e-05 ustep: 1 
#> iter: 5  mgc: 2.480793e-10 
#> iter: 1  mgc: 2.480793e-10 
#> outer mgc:  114.4297 
#> iter: 1  value: 3637.214 mgc: 89.87506 ustep: 1 
#> iter: 2  value: 3600.204 mgc: 6.549042 ustep: 1 
#> iter: 3  value: 3546.542 mgc: 53.18353 ustep: 1 
#> iter: 4  value: 3537.283 mgc: 14.29943 ustep: 0.8994163 
#> iter: 5  value: 3529.784 mgc: 6.183712 ustep: 0.9483808 
#> iter: 6  value: 3529.747 mgc: 0.3292379 ustep: 0.973851 
#> iter: 7  value: 3529.747 mgc: 0.00720933 ustep: 0.9868402 
#> iter: 8  value: 3529.747 mgc: 0.0001751717 ustep: 0.993399 
#> iter: 9  value: 3529.747 mgc: 2.934617e-06 ustep: 0.9966944 
#> iter: 10  value: 3529.747 mgc: 2.940453e-08 ustep: 0.998346 
#> iter: 11  mgc: 7.614122e-10 
#> iter: 1  value: 2051.219 mgc: 31.79301 ustep: 1 
#> iter: 2  value: 2050.718 mgc: 1.227465 ustep: 1 
#> iter: 3  value: 2050.717 mgc: 0.01776177 ustep: 1 
#> iter: 4  value: 2050.717 mgc: 2.244918e-05 ustep: 1 
#> iter: 5  mgc: 5.520548e-11 
#> iter: 1  value: 1938.062 mgc: 5.948224 ustep: 1 
#> iter: 2  value: 1938.061 mgc: 0.05167908 ustep: 1 
#> iter: 3  value: 1938.061 mgc: 4.536121e-05 ustep: 1 
#> iter: 4  mgc: 1.132661e-10 
#> iter: 1  mgc: 1.132661e-10 
#> outer mgc:  112.1901 
#> iter: 1  value: 1932.372 mgc: 5.132122 ustep: 1 
#> iter: 2  value: 1932.371 mgc: 0.04316491 ustep: 1 
#> iter: 3  value: 1932.371 mgc: 8.169268e-05 ustep: 1 
#> iter: 4  mgc: 3.597879e-10 
#> iter: 1  mgc: 3.597879e-10 
#> outer mgc:  104.0288 
#> iter: 1  value: 1924.536 mgc: 7.046888 ustep: 1 
#> iter: 2  value: 1924.534 mgc: 0.05514962 ustep: 1 
#> iter: 3  value: 1924.534 mgc: 0.0001011492 ustep: 1 
#> iter: 4  mgc: 7.261118e-10 
#> iter: 1  mgc: 7.261118e-10 
#> outer mgc:  110.5486 
#> iter: 1  value: 1913.302 mgc: 3.734862 ustep: 1 
#> iter: 2  value: 1913.302 mgc: 0.03964341 ustep: 1 
#> iter: 3  value: 1913.302 mgc: 4.993686e-05 ustep: 1 
#> iter: 4  mgc: 1.338027e-10 
#> iter: 1  mgc: 1.338027e-10 
#> outer mgc:  93.86252 
#> iter: 1  value: 1898.091 mgc: 12.59814 ustep: 1 
#> iter: 2  value: 1898.084 mgc: 0.1290158 ustep: 1 
#> iter: 3  value: 1898.084 mgc: 0.0006211282 ustep: 1 
#> iter: 4  value: 1898.084 mgc: 2.304735e-08 ustep: 1 
#> iter: 5  mgc: 1.38459e-13 
#> iter: 1  mgc: 1.38459e-13 
#> outer mgc:  104.3391 
#> iter: 1  value: 1878.909 mgc: 12.76383 ustep: 1 
#> iter: 2  value: 1878.892 mgc: 0.2070918 ustep: 1 
#> iter: 3  value: 1878.892 mgc: 0.001002897 ustep: 1 
#> iter: 4  value: 1878.892 mgc: 7.108416e-08 ustep: 1 
#> iter: 5  mgc: 9.977728e-14 
#> iter: 1  mgc: 9.977728e-14 
#> outer mgc:  76.84346 
#> iter: 1  value: 1861.935 mgc: 20.89779 ustep: 1 
#> iter: 2  value: 1861.919 mgc: 0.1238725 ustep: 1 
#> iter: 3  value: 1861.919 mgc: 0.001943234 ustep: 1 
#> iter: 4  value: 1861.919 mgc: 3.186256e-07 ustep: 1 
#> iter: 5  mgc: 1.071791e-13 
#> iter: 1  mgc: 1.071791e-13 
#> outer mgc:  85.44911 
#> iter: 1  value: 1893.261 mgc: 26.24283 ustep: 1 
#> iter: 2  value: 1893.19 mgc: 0.3755634 ustep: 1 
#> iter: 3  value: 1893.19 mgc: 0.009246496 ustep: 1 
#> iter: 4  value: 1893.19 mgc: 9.280566e-06 ustep: 1 
#> iter: 5  mgc: 8.241155e-12 
#> iter: 1  value: 1866.944 mgc: 9.563519 ustep: 1 
#> iter: 2  value: 1866.93 mgc: 0.2523388 ustep: 1 
#> iter: 3  value: 1866.93 mgc: 0.003429696 ustep: 1 
#> iter: 4  value: 1866.93 mgc: 5.904504e-07 ustep: 1 
#> iter: 5  mgc: 1.455492e-13 
#> iter: 1  value: 1869.308 mgc: 10.25118 ustep: 1 
#> iter: 2  value: 1869.292 mgc: 0.221746 ustep: 1 
#> iter: 3  value: 1869.292 mgc: 0.002713186 ustep: 1 
#> iter: 4  value: 1869.292 mgc: 3.76011e-07 ustep: 1 
#> iter: 5  mgc: 1.955997e-13 
#> iter: 1  value: 1871.712 mgc: 10.25268 ustep: 1 
#> iter: 2  value: 1871.692 mgc: 0.1828084 ustep: 1 
#> iter: 3  value: 1871.692 mgc: 0.001886039 ustep: 1 
#> iter: 4  value: 1871.692 mgc: 1.854816e-07 ustep: 1 
#> iter: 5  mgc: 1.730684e-13 
#> iter: 1  value: 1859.861 mgc: 2.852587 ustep: 1 
#> iter: 2  value: 1859.861 mgc: 0.01259044 ustep: 1 
#> iter: 3  value: 1859.861 mgc: 7.859349e-06 ustep: 1 
#> iter: 4  mgc: 3.234206e-12 
#> iter: 1  mgc: 3.234206e-12 
#> outer mgc:  42.27835 
#> iter: 1  value: 1858.528 mgc: 2.462126 ustep: 1 
#> iter: 2  value: 1858.528 mgc: 0.0102943 ustep: 1 
#> iter: 3  value: 1858.528 mgc: 1.184747e-06 ustep: 1 
#> iter: 4  mgc: 1.561881e-13 
#> iter: 1  mgc: 1.561881e-13 
#> outer mgc:  55.88646 
#> iter: 1  value: 1856.674 mgc: 2.71141 ustep: 1 
#> iter: 2  value: 1856.674 mgc: 0.01205462 ustep: 1 
#> iter: 3  value: 1856.674 mgc: 8.254081e-06 ustep: 1 
#> iter: 4  mgc: 3.485335e-12 
#> iter: 1  mgc: 3.485335e-12 
#> outer mgc:  55.58047 
#> iter: 1  value: 1854.416 mgc: 4.145365 ustep: 1 
#> iter: 2  value: 1854.416 mgc: 0.01859906 ustep: 1 
#> iter: 3  value: 1854.416 mgc: 2.431592e-05 ustep: 1 
#> iter: 4  mgc: 5.275121e-11 
#> iter: 1  mgc: 5.275121e-11 
#> outer mgc:  46.49094 
#> iter: 1  value: 1851.636 mgc: 2.798182 ustep: 1 
#> iter: 2  value: 1851.636 mgc: 0.01697202 ustep: 1 
#> iter: 3  value: 1851.636 mgc: 1.551575e-05 ustep: 1 
#> iter: 4  mgc: 1.235264e-11 
#> iter: 1  mgc: 1.235264e-11 
#> outer mgc:  52.85448 
#> iter: 1  value: 1846.848 mgc: 5.317538 ustep: 1 
#> iter: 2  value: 1846.847 mgc: 0.06561445 ustep: 1 
#> iter: 3  value: 1846.847 mgc: 5.145749e-05 ustep: 1 
#> iter: 4  mgc: 3.112587e-11 
#> iter: 1  mgc: 3.112587e-11 
#> outer mgc:  33.28613 
#> iter: 1  value: 1847.662 mgc: 5.069742 ustep: 1 
#> iter: 2  value: 1847.661 mgc: 0.0515119 ustep: 1 
#> iter: 3  value: 1847.661 mgc: 0.0001350596 ustep: 1 
#> iter: 4  mgc: 8.352017e-10 
#> iter: 1  mgc: 8.352017e-10 
#> outer mgc:  29.35665 
#> iter: 1  value: 1868.228 mgc: 20.52224 ustep: 1 
#> iter: 2  value: 1868.095 mgc: 0.3327715 ustep: 1 
#> iter: 3  value: 1868.095 mgc: 0.01253892 ustep: 1 
#> iter: 4  value: 1868.095 mgc: 1.763502e-05 ustep: 1 
#> iter: 5  mgc: 4.845466e-11 
#> iter: 1  value: 1844.257 mgc: 5.895744 ustep: 1 
#> iter: 2  value: 1844.256 mgc: 0.0202269 ustep: 1 
#> iter: 3  value: 1844.256 mgc: 6.430089e-05 ustep: 1 
#> iter: 4  mgc: 5.86481e-10 
#> iter: 1  mgc: 5.86481e-10 
#> outer mgc:  61.09649 
#> iter: 1  value: 1848.957 mgc: 6.050627 ustep: 1 
#> iter: 2  value: 1848.953 mgc: 0.08268862 ustep: 1 
#> iter: 3  value: 1848.953 mgc: 0.0003803282 ustep: 1 
#> iter: 4  mgc: 7.460035e-09 
#> iter: 1  value: 1844.117 mgc: 3.478968 ustep: 1 
#> iter: 2  value: 1844.117 mgc: 0.02026256 ustep: 1 
#> iter: 3  value: 1844.117 mgc: 2.225271e-05 ustep: 1 
#> iter: 4  mgc: 3.112106e-11 
#> iter: 1  value: 1843.375 mgc: 1.648868 ustep: 1 
#> iter: 2  value: 1843.375 mgc: 0.004416191 ustep: 1 
#> iter: 3  value: 1843.375 mgc: 1.025642e-06 ustep: 1 
#> iter: 4  mgc: 1.274835e-13 
#> iter: 1  mgc: 1.274835e-13 
#> outer mgc:  20.6319 
#> iter: 1  value: 1842.9 mgc: 1.491342 ustep: 1 
#> iter: 2  value: 1842.9 mgc: 0.004144313 ustep: 1 
#> iter: 3  value: 1842.9 mgc: 2.140215e-07 ustep: 1 
#> iter: 4  mgc: 1.729785e-13 
#> iter: 1  mgc: 1.729785e-13 
#> outer mgc:  33.42338 
#> iter: 1  value: 1842.251 mgc: 1.285974 ustep: 1 
#> iter: 2  value: 1842.251 mgc: 0.003592713 ustep: 1 
#> iter: 3  value: 1842.251 mgc: 7.299581e-07 ustep: 1 
#> iter: 4  mgc: 1.384392e-13 
#> iter: 1  mgc: 1.384392e-13 
#> outer mgc:  34.38701 
#> iter: 1  value: 1841.352 mgc: 1.796275 ustep: 1 
#> iter: 2  value: 1841.352 mgc: 0.004743027 ustep: 1 
#> iter: 3  value: 1841.352 mgc: 4.571805e-07 ustep: 1 
#> iter: 4  mgc: 1.239958e-13 
#> iter: 1  mgc: 1.239958e-13 
#> outer mgc:  25.36312 
#> iter: 1  value: 1839.709 mgc: 1.532969 ustep: 1 
#> iter: 2  value: 1839.709 mgc: 0.008305357 ustep: 1 
#> iter: 3  value: 1839.709 mgc: 8.212436e-07 ustep: 1 
#> iter: 4  mgc: 1.09668e-13 
#> iter: 1  mgc: 1.09668e-13 
#> outer mgc:  27.75771 
#> iter: 1  value: 1838.432 mgc: 3.330485 ustep: 1 
#> iter: 2  value: 1838.432 mgc: 0.03000592 ustep: 1 
#> iter: 3  value: 1838.432 mgc: 2.599976e-06 ustep: 1 
#> iter: 4  mgc: 1.40779e-13 
#> iter: 1  mgc: 1.40779e-13 
#> outer mgc:  27.70165 
#> iter: 1  value: 1835.827 mgc: 1.9724 ustep: 1 
#> iter: 2  value: 1835.827 mgc: 0.01355414 ustep: 1 
#> iter: 3  value: 1835.827 mgc: 5.042175e-06 ustep: 1 
#> iter: 4  mgc: 1.315375e-12 
#> iter: 1  mgc: 1.315375e-12 
#> outer mgc:  15.90539 
#> iter: 1  value: 1838.352 mgc: 12.44564 ustep: 1 
#> iter: 2  value: 1838.348 mgc: 0.1406089 ustep: 1 
#> iter: 3  value: 1838.348 mgc: 0.0001371599 ustep: 1 
#> iter: 4  mgc: 1.013433e-09 
#> iter: 1  value: 1837.556 mgc: 9.174303 ustep: 1 
#> iter: 2  value: 1837.555 mgc: 0.1143099 ustep: 1 
#> iter: 3  value: 1837.555 mgc: 0.0001043729 ustep: 1 
#> iter: 4  mgc: 2.953143e-10 
#> iter: 1  value: 1837.135 mgc: 6.158556 ustep: 1 
#> iter: 2  value: 1837.134 mgc: 0.06781265 ustep: 1 
#> iter: 3  value: 1837.134 mgc: 4.67191e-05 ustep: 1 
#> iter: 4  mgc: 5.228672e-11 
#> iter: 1  value: 1835.561 mgc: 1.8753 ustep: 1 
#> iter: 2  value: 1835.561 mgc: 0.00656595 ustep: 1 
#> iter: 3  value: 1835.561 mgc: 4.588384e-07 ustep: 1 
#> iter: 4  mgc: 1.096555e-13 
#> iter: 1  mgc: 1.096555e-13 
#> outer mgc:  26.66539 
#> iter: 1  value: 1835.222 mgc: 1.056829 ustep: 1 
#> iter: 2  value: 1835.222 mgc: 0.002179779 ustep: 1 
#> iter: 3  value: 1835.222 mgc: 2.783386e-07 ustep: 1 
#> iter: 4  mgc: 1.091283e-13 
#> iter: 1  mgc: 1.091283e-13 
#> outer mgc:  15.99609 
#> iter: 1  value: 1835.073 mgc: 1.054527 ustep: 1 
#> iter: 2  value: 1835.073 mgc: 0.001669695 ustep: 1 
#> iter: 3  value: 1835.073 mgc: 3.652194e-08 ustep: 1 
#> iter: 4  mgc: 1.167081e-13 
#> iter: 1  mgc: 1.167081e-13 
#> outer mgc:  15.94446 
#> iter: 1  value: 1834.733 mgc: 0.4857804 ustep: 1 
#> iter: 2  value: 1834.733 mgc: 0.0007098264 ustep: 1 
#> iter: 3  value: 1834.733 mgc: 2.880114e-08 ustep: 1 
#> mgc: 2.003813e-13 
#> iter: 1  mgc: 2.003813e-13 
#> outer mgc:  16.25764 
#> iter: 1  value: 1834.166 mgc: 1.626817 ustep: 1 
#> iter: 2  value: 1834.166 mgc: 0.003577436 ustep: 1 
#> iter: 3  value: 1834.166 mgc: 1.538106e-07 ustep: 1 
#> iter: 4  mgc: 1.5648e-13 
#> iter: 1  mgc: 1.5648e-13 
#> outer mgc:  15.93048 
#> iter: 1  value: 1835.186 mgc: 2.030324 ustep: 1 
#> iter: 2  value: 1835.186 mgc: 0.01109607 ustep: 1 
#> iter: 3  value: 1835.186 mgc: 7.829109e-07 ustep: 1 
#> iter: 4  mgc: 1.745239e-13 
#> iter: 1  mgc: 1.745239e-13 
#> outer mgc:  16.46464 
#> iter: 1  value: 1833.78 mgc: 2.938971 ustep: 1 
#> iter: 2  value: 1833.78 mgc: 0.01553498 ustep: 1 
#> iter: 3  value: 1833.78 mgc: 6.03116e-06 ustep: 1 
#> iter: 4  mgc: 2.07819e-12 
#> iter: 1  mgc: 2.07819e-12 
#> outer mgc:  67.29718 
#> iter: 1  value: 1830.717 mgc: 5.046202 ustep: 1 
#> iter: 2  value: 1830.716 mgc: 0.07654968 ustep: 1 
#> iter: 3  value: 1830.716 mgc: 6.048742e-05 ustep: 1 
#> iter: 4  mgc: 2.038375e-10 
#> iter: 1  mgc: 2.038375e-10 
#> outer mgc:  31.48924 
#> iter: 1  value: 1829.834 mgc: 1.649223 ustep: 1 
#> iter: 2  value: 1829.834 mgc: 0.004311186 ustep: 1 
#> iter: 3  value: 1829.834 mgc: 1.976916e-07 ustep: 1 
#> iter: 4  mgc: 1.061229e-13 
#> iter: 1  mgc: 1.061229e-13 
#> outer mgc:  10.69557 
#> iter: 1  value: 1828.104 mgc: 1.300638 ustep: 1 
#> iter: 2  value: 1828.104 mgc: 0.003074822 ustep: 1 
#> iter: 3  value: 1828.104 mgc: 2.897004e-07 ustep: 1 
#> iter: 4  mgc: 1.427847e-13 
#> iter: 1  mgc: 1.427847e-13 
#> outer mgc:  9.609329 
#> iter: 1  value: 1830.135 mgc: 4.164722 ustep: 1 
#> iter: 2  value: 1830.134 mgc: 0.01767023 ustep: 1 
#> iter: 3  value: 1830.134 mgc: 7.105836e-06 ustep: 1 
#> iter: 4  mgc: 5.601172e-12 
#> iter: 1  value: 1829.578 mgc: 3.449795 ustep: 1 
#> iter: 2  value: 1829.578 mgc: 0.01547805 ustep: 1 
#> iter: 3  value: 1829.578 mgc: 5.548626e-06 ustep: 1 
#> iter: 4  mgc: 2.502131e-12 
#> iter: 1  value: 1828.473 mgc: 1.818426 ustep: 1 
#> iter: 2  value: 1828.473 mgc: 0.006482108 ustep: 1 
#> iter: 3  value: 1828.473 mgc: 7.290146e-07 ustep: 1 
#> iter: 4  mgc: 9.691131e-14 
#> iter: 1  value: 1828.139 mgc: 0.6892311 ustep: 1 
#> iter: 2  value: 1828.139 mgc: 0.0009423018 ustep: 1 
#> iter: 3  value: 1828.139 mgc: 1.584e-08 ustep: 1 
#> mgc: 2.03129e-13 
#> iter: 1  mgc: 2.03129e-13 
#> outer mgc:  11.94862 
#> iter: 1  value: 1828.118 mgc: 0.542276 ustep: 1 
#> iter: 2  value: 1828.118 mgc: 0.0007299609 ustep: 1 
#> iter: 3  value: 1828.118 mgc: 1.246083e-08 ustep: 1 
#> mgc: 2.429514e-13 
#> iter: 1  mgc: 2.429514e-13 
#> outer mgc:  9.421351 
#> iter: 1  value: 1828.247 mgc: 0.6218052 ustep: 1 
#> iter: 2  value: 1828.247 mgc: 0.0004552423 ustep: 1 
#> iter: 3  mgc: 4.589138e-09 
#> iter: 1  mgc: 4.589138e-09 
#> outer mgc:  9.222598 
#> iter: 1  value: 1828.345 mgc: 0.6221731 ustep: 1 
#> iter: 2  value: 1828.345 mgc: 0.001475258 ustep: 1 
#> iter: 3  value: 1828.345 mgc: 3.374602e-08 ustep: 1 
#> mgc: 2.032237e-13 
#> iter: 1  mgc: 2.032237e-13 
#> outer mgc:  9.112927 
#> iter: 1  value: 1827.872 mgc: 0.6226218 ustep: 1 
#> iter: 2  value: 1827.872 mgc: 0.001141691 ustep: 1 
#> iter: 3  value: 1827.872 mgc: 3.653024e-08 ustep: 1 
#> mgc: 1.358107e-13 
#> iter: 1  value: 1827.303 mgc: 1.092165 ustep: 1 
#> iter: 2  value: 1827.303 mgc: 0.003276104 ustep: 1 
#> iter: 3  value: 1827.303 mgc: 2.420489e-07 ustep: 1 
#> iter: 4  mgc: 9.115292e-14 
#> iter: 1  mgc: 9.115292e-14 
#> outer mgc:  8.187901 
#> iter: 1  value: 1827.468 mgc: 2.66328 ustep: 1 
#> iter: 2  value: 1827.468 mgc: 0.03075308 ustep: 1 
#> iter: 3  value: 1827.468 mgc: 8.359944e-06 ustep: 1 
#> iter: 4  mgc: 4.415704e-13 
#> iter: 1  mgc: 4.415704e-13 
#> outer mgc:  19.00968 
#> iter: 1  value: 1828.127 mgc: 2.933561 ustep: 1 
#> iter: 2  value: 1828.127 mgc: 0.009092746 ustep: 1 
#> iter: 3  value: 1828.127 mgc: 1.342777e-06 ustep: 1 
#> iter: 4  mgc: 1.990448e-13 
#> iter: 1  mgc: 1.990448e-13 
#> outer mgc:  30.24515 
#> iter: 1  value: 1827.723 mgc: 2.790517 ustep: 1 
#> iter: 2  value: 1827.723 mgc: 0.03569617 ustep: 1 
#> iter: 3  value: 1827.723 mgc: 1.121054e-05 ustep: 1 
#> iter: 4  mgc: 1.52965e-12 
#> iter: 1  mgc: 1.52965e-12 
#> outer mgc:  10.17029 
#> iter: 1  value: 1827.374 mgc: 2.109408 ustep: 1 
#> iter: 2  value: 1827.374 mgc: 0.02083945 ustep: 1 
#> iter: 3  value: 1827.374 mgc: 3.839992e-06 ustep: 1 
#> iter: 4  mgc: 2.134243e-13 
#> iter: 1  mgc: 2.134243e-13 
#> outer mgc:  8.460503 
#> iter: 1  value: 1827.248 mgc: 2.629251 ustep: 1 
#> iter: 2  value: 1827.247 mgc: 0.0284613 ustep: 1 
#> iter: 3  value: 1827.247 mgc: 1.913786e-05 ustep: 1 
#> iter: 4  mgc: 1.903345e-11 
#> iter: 1  mgc: 1.903345e-11 
#> outer mgc:  35.58816 
#> iter: 1  value: 1826.751 mgc: 2.313529 ustep: 1 
#> iter: 2  value: 1826.751 mgc: 0.008842839 ustep: 1 
#> iter: 3  value: 1826.751 mgc: 1.285305e-06 ustep: 1 
#> iter: 4  mgc: 1.367969e-13 
#> iter: 1  mgc: 1.367969e-13 
#> outer mgc:  30.22536 
#> iter: 1  value: 1826.065 mgc: 1.800481 ustep: 1 
#> iter: 2  value: 1826.065 mgc: 0.02050504 ustep: 1 
#> iter: 3  value: 1826.065 mgc: 4.937588e-06 ustep: 1 
#> iter: 4  mgc: 6.121957e-13 
#> iter: 1  mgc: 6.121957e-13 
#> outer mgc:  5.399189 
#> iter: 1  value: 1826.386 mgc: 1.871116 ustep: 1 
#> iter: 2  value: 1826.386 mgc: 0.01528655 ustep: 1 
#> iter: 3  value: 1826.386 mgc: 1.323202e-05 ustep: 1 
#> iter: 4  mgc: 1.135865e-11 
#> iter: 1  value: 1826.372 mgc: 2.175956 ustep: 1 
#> iter: 2  value: 1826.372 mgc: 0.005173616 ustep: 1 
#> iter: 3  value: 1826.372 mgc: 3.497958e-06 ustep: 1 
#> iter: 4  mgc: 1.978876e-12 
#> iter: 1  value: 1825.928 mgc: 0.9755895 ustep: 1 
#> iter: 2  value: 1825.928 mgc: 0.001486281 ustep: 1 
#> iter: 3  value: 1825.928 mgc: 6.105186e-08 ustep: 1 
#> iter: 4  mgc: 1.47272e-13 
#> iter: 1  mgc: 1.47272e-13 
#> outer mgc:  6.61011 
#> iter: 1  value: 1825.939 mgc: 0.7054563 ustep: 1 
#> iter: 2  value: 1825.939 mgc: 0.0002093981 ustep: 1 
#> iter: 3  mgc: 4.521354e-09 
#> iter: 1  mgc: 4.521354e-09 
#> outer mgc:  7.361866 
#> iter: 1  value: 1825.898 mgc: 0.8302428 ustep: 1 
#> iter: 2  value: 1825.898 mgc: 0.003174869 ustep: 1 
#> iter: 3  value: 1825.898 mgc: 8.857776e-08 ustep: 1 
#> iter: 4  mgc: 1.494103e-13 
#> iter: 1  mgc: 1.494103e-13 
#> outer mgc:  4.663853 
#> iter: 1  value: 1825.875 mgc: 0.2744759 ustep: 1 
#> iter: 2  value: 1825.875 mgc: 0.0002750709 ustep: 1 
#> iter: 3  mgc: 2.059932e-09 
#> iter: 1  value: 1825.836 mgc: 0.5523232 ustep: 1 
#> iter: 2  value: 1825.836 mgc: 0.002462281 ustep: 1 
#> iter: 3  value: 1825.836 mgc: 6.696386e-08 ustep: 1 
#> mgc: 1.32559e-13 
#> iter: 1  value: 1825.848 mgc: 1.734083 ustep: 1 
#> iter: 2  value: 1825.848 mgc: 0.02502708 ustep: 1 
#> iter: 3  value: 1825.848 mgc: 6.843314e-06 ustep: 1 
#> iter: 4  mgc: 1.060804e-12 
#> iter: 1  mgc: 1.060804e-12 
#> outer mgc:  7.001753 
#> iter: 1  value: 1825.976 mgc: 3.421943 ustep: 1 
#> iter: 2  value: 1825.976 mgc: 0.09431103 ustep: 1 
#> iter: 3  value: 1825.976 mgc: 7.177912e-05 ustep: 1 
#> iter: 4  mgc: 3.973408e-11 
#> iter: 1  mgc: 3.973408e-11 
#> outer mgc:  10.73881 
#> iter: 1  value: 1825.273 mgc: 4.122913 ustep: 1 
#> iter: 2  value: 1825.272 mgc: 0.1281336 ustep: 1 
#> iter: 3  value: 1825.272 mgc: 0.0001485703 ustep: 1 
#> iter: 4  mgc: 5.73337e-10 
#> iter: 1  mgc: 5.73337e-10 
#> outer mgc:  30.10412 
#> iter: 1  value: 1824.687 mgc: 4.53048 ustep: 1 
#> iter: 2  value: 1824.686 mgc: 0.1300817 ustep: 1 
#> iter: 3  value: 1824.686 mgc: 0.0001479498 ustep: 1 
#> iter: 4  mgc: 1.590132e-10 
#> iter: 1  mgc: 1.590132e-10 
#> outer mgc:  10.09854 
#> iter: 1  value: 1824.262 mgc: 1.909847 ustep: 1 
#> iter: 2  value: 1824.262 mgc: 0.01029808 ustep: 1 
#> iter: 3  value: 1824.262 mgc: 1.719248e-06 ustep: 1 
#> iter: 4  mgc: 1.400375e-13 
#> iter: 1  mgc: 1.400375e-13 
#> outer mgc:  5.280114 
#> iter: 1  value: 1824.931 mgc: 1.706874 ustep: 1 
#> iter: 2  value: 1824.931 mgc: 0.009552343 ustep: 1 
#> iter: 3  value: 1824.931 mgc: 9.831731e-07 ustep: 1 
#> iter: 4  mgc: 1.242698e-13 
#> iter: 1  value: 1825.039 mgc: 2.334561 ustep: 1 
#> iter: 2  value: 1825.039 mgc: 0.01030263 ustep: 1 
#> iter: 3  value: 1825.039 mgc: 1.516602e-06 ustep: 1 
#> iter: 4  mgc: 1.669137e-13 
#> iter: 1  value: 1824.42 mgc: 1.238205 ustep: 1 
#> iter: 2  value: 1824.42 mgc: 0.003445966 ustep: 1 
#> iter: 3  value: 1824.42 mgc: 1.518063e-07 ustep: 1 
#> iter: 4  mgc: 1.422285e-13 
#> iter: 1  value: 1824.26 mgc: 0.3311691 ustep: 1 
#> iter: 2  value: 1824.26 mgc: 0.0002491291 ustep: 1 
#> iter: 3  mgc: 8.040823e-10 
#> iter: 1  mgc: 8.040823e-10 
#> outer mgc:  5.712031 
#> iter: 1  value: 1824.24 mgc: 0.2465093 ustep: 1 
#> iter: 2  value: 1824.24 mgc: 9.76585e-05 ustep: 1 
#> iter: 3  mgc: 5.948406e-10 
#> iter: 1  mgc: 5.948406e-10 
#> outer mgc:  4.982434 
#> iter: 1  value: 1824.236 mgc: 0.2485692 ustep: 1 
#> iter: 2  value: 1824.236 mgc: 0.0001223164 ustep: 1 
#> iter: 3  mgc: 2.299439e-10 
#> iter: 1  mgc: 2.299439e-10 
#> outer mgc:  5.353427 
#> iter: 1  value: 1824.206 mgc: 0.3066426 ustep: 1 
#> iter: 2  value: 1824.206 mgc: 0.0001618631 ustep: 1 
#> iter: 3  mgc: 1.602892e-09 
#> iter: 1  mgc: 1.602892e-09 
#> outer mgc:  4.635515 
#> iter: 1  value: 1824.19 mgc: 0.216181 ustep: 1 
#> iter: 2  value: 1824.19 mgc: 0.0001065236 ustep: 1 
#> iter: 3  mgc: 1.229275e-10 
#> iter: 1  mgc: 1.229275e-10 
#> outer mgc:  5.23998 
#> iter: 1  value: 1824.172 mgc: 0.3095059 ustep: 1 
#> iter: 2  value: 1824.172 mgc: 0.0005411666 ustep: 1 
#> iter: 3  mgc: 8.007697e-09 
#> iter: 1  mgc: 8.007697e-09 
#> outer mgc:  4.476694 
#> iter: 1  value: 1824.146 mgc: 0.2347218 ustep: 1 
#> iter: 2  value: 1824.146 mgc: 0.000130208 ustep: 1 
#> iter: 3  mgc: 3.336395e-10 
#> iter: 1  mgc: 3.336395e-10 
#> outer mgc:  4.02874 
#> iter: 1  value: 1823.959 mgc: 0.7543415 ustep: 1 
#> iter: 2  value: 1823.959 mgc: 0.003277073 ustep: 1 
#> iter: 3  value: 1823.959 mgc: 1.164473e-07 ustep: 1 
#> iter: 4  mgc: 1.46394e-13 
#> iter: 1  mgc: 1.46394e-13 
#> outer mgc:  2.83525 
#> iter: 1  value: 1824.063 mgc: 0.8584383 ustep: 1 
#> iter: 2  value: 1824.062 mgc: 0.007151257 ustep: 1 
#> iter: 3  value: 1824.062 mgc: 1.384888e-06 ustep: 1 
#> iter: 4  mgc: 1.331818e-13 
#> iter: 1  mgc: 1.331818e-13 
#> outer mgc:  2.264246 
#> iter: 1  value: 1824.285 mgc: 2.873974 ustep: 1 
#> iter: 2  value: 1824.285 mgc: 0.03604026 ustep: 1 
#> iter: 3  value: 1824.285 mgc: 1.163727e-05 ustep: 1 
#> iter: 4  mgc: 9.913225e-13 
#> iter: 1  value: 1824.23 mgc: 2.002931 ustep: 1 
#> iter: 2  value: 1824.23 mgc: 0.01444298 ustep: 1 
#> iter: 3  value: 1824.23 mgc: 1.875276e-06 ustep: 1 
#> iter: 4  mgc: 1.129319e-13 
#> iter: 1  value: 1824.217 mgc: 1.459604 ustep: 1 
#> iter: 2  value: 1824.217 mgc: 0.005727753 ustep: 1 
#> iter: 3  value: 1824.217 mgc: 3.224819e-07 ustep: 1 
#> iter: 4  mgc: 1.029058e-13 
#> iter: 1  value: 1824.078 mgc: 0.492325 ustep: 1 
#> iter: 2  value: 1824.078 mgc: 0.0005195313 ustep: 1 
#> iter: 3  mgc: 4.561323e-09 
#> iter: 1  value: 1824.063 mgc: 0.2286455 ustep: 1 
#> iter: 2  value: 1824.063 mgc: 0.000111706 ustep: 1 
#> iter: 3  mgc: 2.093717e-10 
#> iter: 1  mgc: 2.093717e-10 
#> outer mgc:  4.272081 
#> iter: 1  value: 1824.068 mgc: 0.1530367 ustep: 1 
#> iter: 2  value: 1824.068 mgc: 4.778486e-05 ustep: 1 
#> iter: 3  mgc: 5.580381e-11 
#> iter: 1  mgc: 5.580381e-11 
#> outer mgc:  2.265795 
#> iter: 1  value: 1824.067 mgc: 0.1729325 ustep: 1 
#> iter: 2  value: 1824.067 mgc: 0.0001018334 ustep: 1 
#> iter: 3  mgc: 9.294731e-11 
#> iter: 1  value: 1824.069 mgc: 0.06500758 ustep: 1 
#> iter: 2  value: 1824.069 mgc: 2.287013e-05 ustep: 1 
#> iter: 3  mgc: 4.613868e-12 
#> iter: 1  value: 1824.078 mgc: 0.2310041 ustep: 1 
#> iter: 2  value: 1824.078 mgc: 0.0002899211 ustep: 1 
#> iter: 3  mgc: 7.364161e-10 
#> iter: 1  value: 1824.13 mgc: 0.7688327 ustep: 1 
#> iter: 2  value: 1824.13 mgc: 0.003252788 ustep: 1 
#> iter: 3  value: 1824.13 mgc: 9.271975e-08 ustep: 1 
#> mgc: 1.308734e-13 
#> iter: 1  mgc: 1.308734e-13 
#> outer mgc:  3.660808 
#> iter: 1  value: 1824.154 mgc: 0.495964 ustep: 1 
#> iter: 2  value: 1824.154 mgc: 0.00193172 ustep: 1 
#> iter: 3  value: 1824.154 mgc: 1.032694e-07 ustep: 1 
#> iter: 4  mgc: 1.208653e-13 
#> iter: 1  mgc: 1.208653e-13 
#> outer mgc:  1.787063 
#> iter: 1  value: 1824.028 mgc: 1.155827 ustep: 1 
#> iter: 2  value: 1824.028 mgc: 0.007003203 ustep: 1 
#> iter: 3  value: 1824.028 mgc: 4.38671e-07 ustep: 1 
#> iter: 4  mgc: 1.139473e-13 
#> iter: 1  mgc: 1.139473e-13 
#> outer mgc:  3.130018 
#> iter: 1  value: 1824.108 mgc: 0.9680725 ustep: 1 
#> iter: 2  value: 1824.108 mgc: 0.00107494 ustep: 1 
#> iter: 3  value: 1824.108 mgc: 2.356953e-08 ustep: 1 
#> iter: 4  mgc: 1.143239e-13 
#> iter: 1  mgc: 1.143239e-13 
#> outer mgc:  7.599048 
#> iter: 1  value: 1824.013 mgc: 1.415297 ustep: 1 
#> iter: 2  value: 1824.013 mgc: 0.01046189 ustep: 1 
#> iter: 3  value: 1824.013 mgc: 9.756765e-07 ustep: 1 
#> iter: 4  mgc: 1.620627e-13 
#> iter: 1  mgc: 1.620627e-13 
#> outer mgc:  3.802895 
#> iter: 1  value: 1823.947 mgc: 0.3322677 ustep: 1 
#> iter: 2  value: 1823.947 mgc: 0.001075045 ustep: 1 
#> iter: 3  value: 1823.947 mgc: 2.815618e-08 ustep: 1 
#> mgc: 2.260007e-13 
#> iter: 1  mgc: 2.260007e-13 
#> outer mgc:  3.23066 
#> iter: 1  value: 1823.95 mgc: 1.015586 ustep: 1 
#> iter: 2  value: 1823.95 mgc: 0.005203593 ustep: 1 
#> iter: 3  value: 1823.95 mgc: 2.427958e-07 ustep: 1 
#> iter: 4  mgc: 1.115336e-13 
#> iter: 1  mgc: 1.115336e-13 
#> outer mgc:  4.354267 
#> iter: 1  value: 1823.908 mgc: 0.4703953 ustep: 1 
#> iter: 2  value: 1823.908 mgc: 0.001374505 ustep: 1 
#> iter: 3  value: 1823.908 mgc: 4.540691e-08 ustep: 1 
#> mgc: 1.347449e-13 
#> iter: 1  mgc: 1.347449e-13 
#> outer mgc:  2.750474 
#> iter: 1  value: 1824.053 mgc: 1.241126 ustep: 1 
#> iter: 2  value: 1824.053 mgc: 0.001013025 ustep: 1 
#> iter: 3  mgc: 4.703549e-09 
#> iter: 1  mgc: 4.703549e-09 
#> outer mgc:  1.825814 
#> iter: 1  value: 1824.122 mgc: 0.4601509 ustep: 1 
#> iter: 2  value: 1824.122 mgc: 0.001319868 ustep: 1 
#> iter: 3  value: 1824.122 mgc: 1.50286e-08 ustep: 1 
#> mgc: 1.824102e-13 
#> iter: 1  mgc: 1.824102e-13 
#> outer mgc:  4.812017 
#> iter: 1  value: 1823.979 mgc: 1.084329 ustep: 1 
#> iter: 2  value: 1823.979 mgc: 0.005230436 ustep: 1 
#> iter: 3  value: 1823.979 mgc: 2.456124e-07 ustep: 1 
#> iter: 4  mgc: 1.902868e-13 
#> iter: 1  mgc: 1.902868e-13 
#> outer mgc:  3.896554 
#> iter: 1  value: 1824.062 mgc: 0.2994343 ustep: 1 
#> iter: 2  value: 1824.062 mgc: 0.0004425831 ustep: 1 
#> iter: 3  mgc: 5.006242e-09 
#> iter: 1  mgc: 5.006242e-09 
#> outer mgc:  2.063378 
#> iter: 1  value: 1824.08 mgc: 0.8079386 ustep: 1 
#> iter: 2  value: 1824.08 mgc: 0.004137006 ustep: 1 
#> iter: 3  value: 1824.08 mgc: 4.998607e-07 ustep: 1 
#> iter: 4  mgc: 1.896218e-13 
#> iter: 1  value: 1824.067 mgc: 0.5391846 ustep: 1 
#> iter: 2  value: 1824.067 mgc: 0.0006020665 ustep: 1 
#> iter: 3  mgc: 5.547206e-09 
#> iter: 1  value: 1824.08 mgc: 0.4888943 ustep: 1 
#> iter: 2  value: 1824.08 mgc: 0.0003738202 ustep: 1 
#> iter: 3  mgc: 1.520584e-09 
#> iter: 1  value: 1824.087 mgc: 0.4665453 ustep: 1 
#> iter: 2  value: 1824.087 mgc: 0.0003781937 ustep: 1 
#> iter: 3  mgc: 3.410246e-09 
#> iter: 1  value: 1824.069 mgc: 0.2667973 ustep: 1 
#> iter: 2  value: 1824.069 mgc: 0.0001606028 ustep: 1 
#> iter: 3  mgc: 4.791398e-10 
#> iter: 1  value: 1824.062 mgc: 0.1051795 ustep: 1 
#> iter: 2  value: 1824.062 mgc: 2.50091e-05 ustep: 1 
#> iter: 3  mgc: 1.167874e-11 
#> iter: 1  mgc: 1.167874e-11 
#> outer mgc:  1.971203 
#> iter: 1  value: 1824.055 mgc: 0.06863991 ustep: 1 
#> iter: 2  value: 1824.055 mgc: 4.350973e-06 ustep: 1 
#> iter: 3  mgc: 1.209361e-12 
#> iter: 1  mgc: 1.209361e-12 
#> outer mgc:  1.240925 
#> iter: 1  value: 1824.055 mgc: 0.07051201 ustep: 1 
#> iter: 2  value: 1824.055 mgc: 7.975547e-06 ustep: 1 
#> iter: 3  mgc: 1.616753e-12 
#> iter: 1  mgc: 1.616753e-12 
#> outer mgc:  1.736895 
#> iter: 1  value: 1824.049 mgc: 0.05996945 ustep: 1 
#> iter: 2  value: 1824.049 mgc: 9.736303e-06 ustep: 1 
#> iter: 3  mgc: 2.2672e-12 
#> iter: 1  mgc: 2.2672e-12 
#> outer mgc:  0.8471655 
#> iter: 1  value: 1824.046 mgc: 0.04906514 ustep: 1 
#> iter: 2  value: 1824.046 mgc: 4.94029e-06 ustep: 1 
#> iter: 3  mgc: 3.826327e-13 
#> iter: 1  value: 1824.043 mgc: 0.03895909 ustep: 1 
#> iter: 2  value: 1824.043 mgc: 1.18523e-05 ustep: 1 
#> iter: 3  mgc: 3.848372e-12 
#> iter: 1  value: 1824.033 mgc: 0.1211028 ustep: 1 
#> iter: 2  value: 1824.033 mgc: 0.000114526 ustep: 1 
#> iter: 3  mgc: 3.57758e-10 
#> iter: 1  mgc: 3.57758e-10 
#> outer mgc:  1.691278 
#> iter: 1  value: 1824.029 mgc: 0.1407216 ustep: 1 
#> iter: 2  value: 1824.029 mgc: 4.689854e-05 ustep: 1 
#> iter: 3  mgc: 5.458795e-11 
#> iter: 1  mgc: 5.458795e-11 
#> outer mgc:  0.5872276 
#> iter: 1  value: 1823.993 mgc: 0.2206256 ustep: 1 
#> iter: 2  value: 1823.993 mgc: 0.0003744227 ustep: 1 
#> iter: 3  mgc: 1.886336e-09 
#> iter: 1  mgc: 1.886336e-09 
#> outer mgc:  0.6532466 
#> iter: 1  value: 1823.955 mgc: 0.2009748 ustep: 1 
#> iter: 2  value: 1823.955 mgc: 7.224823e-05 ustep: 1 
#> iter: 3  mgc: 5.979143e-11 
#> iter: 1  mgc: 5.979143e-11 
#> outer mgc:  1.952536 
#> iter: 1  value: 1823.908 mgc: 0.2459614 ustep: 1 
#> iter: 2  value: 1823.908 mgc: 9.121918e-05 ustep: 1 
#> iter: 3  mgc: 8.797386e-11 
#> iter: 1  mgc: 8.797386e-11 
#> outer mgc:  1.551127 
#> iter: 1  value: 1823.976 mgc: 0.2762101 ustep: 1 
#> iter: 2  value: 1823.976 mgc: 9.54326e-05 ustep: 1 
#> iter: 3  mgc: 1.931405e-10 
#> iter: 1  mgc: 1.931405e-10 
#> outer mgc:  0.6309545 
#> iter: 1  value: 1823.965 mgc: 0.2200542 ustep: 1 
#> iter: 2  value: 1823.965 mgc: 4.27141e-05 ustep: 1 
#> iter: 3  mgc: 4.251205e-11 
#> iter: 1  mgc: 4.251205e-11 
#> outer mgc:  0.8189808 
#> iter: 1  value: 1824.006 mgc: 0.1643511 ustep: 1 
#> iter: 2  value: 1824.006 mgc: 0.0001091399 ustep: 1 
#> iter: 3  mgc: 1.177274e-10 
#> iter: 1  mgc: 1.177274e-10 
#> outer mgc:  1.319348 
#> iter: 1  value: 1824.009 mgc: 0.3711326 ustep: 1 
#> iter: 2  value: 1824.009 mgc: 0.0003254007 ustep: 1 
#> iter: 3  mgc: 9.253024e-10 
#> iter: 1  mgc: 9.253024e-10 
#> outer mgc:  1.999095 
#> iter: 1  value: 1824.001 mgc: 0.2443711 ustep: 1 
#> iter: 2  value: 1824.001 mgc: 9.427507e-05 ustep: 1 
#> iter: 3  mgc: 6.109856e-11 
#> iter: 1  mgc: 6.109856e-11 
#> outer mgc:  0.7188707 
#> iter: 1  value: 1824.038 mgc: 0.2032706 ustep: 1 
#> iter: 2  value: 1824.038 mgc: 0.0001546213 ustep: 1 
#> iter: 3  mgc: 6.335301e-10 
#> iter: 1  mgc: 6.335301e-10 
#> outer mgc:  0.5654962 
#> iter: 1  value: 1824.008 mgc: 0.2546046 ustep: 1 
#> iter: 2  value: 1824.008 mgc: 5.899977e-05 ustep: 1 
#> iter: 3  mgc: 5.013881e-11 
#> iter: 1  mgc: 5.013881e-11 
#> outer mgc:  1.074703 
#> iter: 1  value: 1823.975 mgc: 0.1594867 ustep: 1 
#> iter: 2  value: 1823.975 mgc: 0.0001492134 ustep: 1 
#> iter: 3  mgc: 6.025408e-10 
#> iter: 1  mgc: 6.025408e-10 
#> outer mgc:  1.38664 
#> iter: 1  value: 1823.945 mgc: 0.1449545 ustep: 1 
#> iter: 2  value: 1823.945 mgc: 6.799124e-05 ustep: 1 
#> iter: 3  mgc: 4.769614e-11 
#> iter: 1  mgc: 4.769614e-11 
#> outer mgc:  1.218415 
#> iter: 1  value: 1823.922 mgc: 0.1163433 ustep: 1 
#> iter: 2  value: 1823.922 mgc: 8.614095e-05 ustep: 1 
#> iter: 3  mgc: 1.458135e-10 
#> iter: 1  mgc: 1.458135e-10 
#> outer mgc:  0.3674795 
#> iter: 1  value: 1823.951 mgc: 0.2283347 ustep: 1 
#> iter: 2  value: 1823.951 mgc: 9.131184e-05 ustep: 1 
#> iter: 3  mgc: 2.145214e-10 
#> iter: 1  mgc: 2.145214e-10 
#> outer mgc:  1.144934 
#> iter: 1  value: 1823.986 mgc: 0.09489399 ustep: 1 
#> iter: 2  value: 1823.986 mgc: 3.427643e-05 ustep: 1 
#> iter: 3  mgc: 2.074076e-11 
#> iter: 1  mgc: 2.074076e-11 
#> outer mgc:  0.4240457 
#> iter: 1  value: 1824.017 mgc: 0.08433126 ustep: 1 
#> iter: 2  value: 1824.017 mgc: 3.270083e-05 ustep: 1 
#> iter: 3  mgc: 9.899874e-12 
#> iter: 1  mgc: 9.899874e-12 
#> outer mgc:  1.069755 
#> iter: 1  value: 1824.013 mgc: 0.2256685 ustep: 1 
#> iter: 2  value: 1824.013 mgc: 0.0001749566 ustep: 1 
#> iter: 3  mgc: 2.691457e-10 
#> iter: 1  mgc: 2.691457e-10 
#> outer mgc:  0.5165124 
#> iter: 1  value: 1823.993 mgc: 0.05543808 ustep: 1 
#> iter: 2  value: 1823.993 mgc: 1.862418e-05 ustep: 1 
#> iter: 3  mgc: 7.676922e-12 
#> iter: 1  mgc: 7.676922e-12 
#> outer mgc:  0.2053157 
#> iter: 1  value: 1823.979 mgc: 0.07684606 ustep: 1 
#> iter: 2  value: 1823.979 mgc: 1.359223e-05 ustep: 1 
#> iter: 3  mgc: 1.999144e-12 
#> iter: 1  mgc: 1.999144e-12 
#> outer mgc:  0.5977271 
#> iter: 1  value: 1823.978 mgc: 0.08496553 ustep: 1 
#> iter: 2  value: 1823.978 mgc: 1.731873e-05 ustep: 1 
#> iter: 3  mgc: 5.120062e-12 
#> iter: 1  mgc: 5.120062e-12 
#> outer mgc:  0.4309044 
#> iter: 1  value: 1823.97 mgc: 0.07374451 ustep: 1 
#> iter: 2  value: 1823.97 mgc: 2.058229e-05 ustep: 1 
#> iter: 3  mgc: 4.687613e-12 
#> iter: 1  mgc: 4.687613e-12 
#> outer mgc:  0.1918016 
#> iter: 1  value: 1823.979 mgc: 0.07530742 ustep: 1 
#> iter: 2  value: 1823.979 mgc: 1.86208e-05 ustep: 1 
#> iter: 3  mgc: 4.950928e-12 
#> iter: 1  mgc: 4.950928e-12 
#> outer mgc:  0.2855536 
#> iter: 1  value: 1823.99 mgc: 0.0535193 ustep: 1 
#> iter: 2  value: 1823.99 mgc: 4.991182e-06 ustep: 1 
#> iter: 3  mgc: 2.066622e-13 
#> iter: 1  mgc: 2.066622e-13 
#> outer mgc:  0.3807738 
#> iter: 1  value: 1823.998 mgc: 0.06641819 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 1.242823e-05 ustep: 1 
#> iter: 3  mgc: 7.233067e-13 
#> iter: 1  mgc: 7.233067e-13 
#> outer mgc:  0.2868459 
#> iter: 1  value: 1824.002 mgc: 0.1102196 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 1.949244e-05 ustep: 1 
#> iter: 3  mgc: 4.383996e-12 
#> iter: 1  mgc: 4.383996e-12 
#> outer mgc:  0.3071323 
#> iter: 1  value: 1824.01 mgc: 0.05341373 ustep: 1 
#> iter: 2  value: 1824.01 mgc: 7.842989e-06 ustep: 1 
#> iter: 3  mgc: 2.325446e-12 
#> iter: 1  mgc: 2.325446e-12 
#> outer mgc:  0.3815191 
#> iter: 1  value: 1824 mgc: 0.08089564 ustep: 1 
#> iter: 2  value: 1824 mgc: 1.452026e-05 ustep: 1 
#> iter: 3  mgc: 1.042973e-12 
#> iter: 1  mgc: 1.042973e-12 
#> outer mgc:  0.3803674 
#> iter: 1  value: 1823.986 mgc: 0.1844588 ustep: 1 
#> iter: 2  value: 1823.986 mgc: 2.762839e-05 ustep: 1 
#> iter: 3  mgc: 1.297607e-11 
#> iter: 1  value: 1823.995 mgc: 0.05832573 ustep: 1 
#> iter: 2  value: 1823.995 mgc: 2.588334e-06 ustep: 1 
#> iter: 3  mgc: 1.565805e-13 
#> iter: 1  mgc: 1.565805e-13 
#> outer mgc:  0.2007175 
#> iter: 1  value: 1823.991 mgc: 0.07168108 ustep: 1 
#> iter: 2  value: 1823.991 mgc: 1.455729e-05 ustep: 1 
#> iter: 3  mgc: 1.025624e-12 
#> iter: 1  mgc: 1.025624e-12 
#> outer mgc:  0.2501556 
#> iter: 1  value: 1823.987 mgc: 0.07548632 ustep: 1 
#> iter: 2  value: 1823.987 mgc: 3.47845e-06 ustep: 1 
#> iter: 3  mgc: 1.976236e-13 
#> iter: 1  mgc: 1.976236e-13 
#> outer mgc:  0.1138994 
#> iter: 1  value: 1823.997 mgc: 0.0422356 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.420186e-06 ustep: 1 
#> iter: 3  mgc: 2.850502e-13 
#> iter: 1  mgc: 2.850502e-13 
#> outer mgc:  0.1729511 
#> iter: 1  value: 1824.004 mgc: 0.02738386 ustep: 1 
#> iter: 2  value: 1824.004 mgc: 2.145204e-06 ustep: 1 
#> iter: 3  mgc: 1.107321e-13 
#> iter: 1  mgc: 1.107321e-13 
#> outer mgc:  0.116876 
#> iter: 1  value: 1824.002 mgc: 0.02910713 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 1.995929e-06 ustep: 1 
#> iter: 3  mgc: 1.513417e-13 
#> iter: 1  mgc: 1.513417e-13 
#> outer mgc:  0.1396068 
#> iter: 1  value: 1824 mgc: 0.03914133 ustep: 1 
#> iter: 2  value: 1824 mgc: 1.575583e-06 ustep: 1 
#> iter: 3  mgc: 2.270539e-13 
#> iter: 1  mgc: 2.270539e-13 
#> outer mgc:  0.1242514 
#> iter: 1  value: 1823.993 mgc: 0.02802661 ustep: 1 
#> iter: 2  value: 1823.993 mgc: 3.628643e-06 ustep: 1 
#> iter: 3  mgc: 1.271175e-13 
#> iter: 1  value: 1823.998 mgc: 0.01303762 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 5.482558e-07 ustep: 1 
#> iter: 3  mgc: 1.234563e-13 
#> iter: 1  mgc: 1.234563e-13 
#> outer mgc:  0.07665125 
#> iter: 1  value: 1823.997 mgc: 0.05239491 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 1.482708e-06 ustep: 1 
#> iter: 3  mgc: 1.846945e-13 
#> iter: 1  value: 1823.997 mgc: 0.01886457 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 1.420144e-07 ustep: 1 
#> iter: 3  mgc: 2.193166e-13 
#> iter: 1  mgc: 2.193166e-13 
#> outer mgc:  0.08158812 
#> iter: 1  value: 1824 mgc: 0.01339481 ustep: 1 
#> iter: 2  value: 1824 mgc: 5.102133e-07 ustep: 1 
#> iter: 3  mgc: 1.056448e-13 
#> iter: 1  mgc: 1.056448e-13 
#> outer mgc:  0.02006747 
#> iter: 1  value: 1824 mgc: 0.005626234 ustep: 1 
#> iter: 2  value: 1824 mgc: 7.943239e-08 ustep: 1 
#> mgc: 1.634923e-13 
#> iter: 1  mgc: 1.634923e-13 
#> outer mgc:  0.02466537 
#> iter: 1  value: 1824 mgc: 0.007413305 ustep: 1 
#> iter: 2  value: 1824 mgc: 6.030482e-08 ustep: 1 
#> iter: 3  mgc: 1.135483e-13 
#> iter: 1  mgc: 1.135483e-13 
#> outer mgc:  0.0229146 
#> iter: 1  value: 1823.998 mgc: 0.007382067 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 2.252063e-08 ustep: 1 
#> mgc: 1.338481e-13 
#> iter: 1  mgc: 1.338481e-13 
#> outer mgc:  0.007657609 
#> iter: 1  value: 1823.998 mgc: 0.003027093 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 2.895632e-08 ustep: 1 
#> mgc: 1.591077e-13 
#> iter: 1  mgc: 1.591077e-13 
#> outer mgc:  0.01446768 
#> iter: 1  value: 1823.999 mgc: 0.003144647 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 3.19317e-08 ustep: 1 
#> mgc: 1.678375e-13 
#> iter: 1  mgc: 1.678375e-13 
#> outer mgc:  0.01478877 
#> iter: 1  value: 1823.999 mgc: 0.00376451 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.91846e-08 ustep: 1 
#> mgc: 1.413485e-13 
#> iter: 1  mgc: 1.413485e-13 
#> outer mgc:  0.006681925 
#> iter: 1  value: 1823.999 mgc: 0.007683415 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.053505e-07 ustep: 1 
#> mgc: 1.273945e-13 
#> iter: 1  value: 1823.999 mgc: 0.001970661 ustep: 1 
#> iter: 2  mgc: 7.0311e-09 
#> iter: 1  mgc: 7.0311e-09 
#> outer mgc:  0.008194309 
#> iter: 1  value: 1823.999 mgc: 0.002561813 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 2.911621e-08 ustep: 1 
#> mgc: 1.256115e-13 
#> iter: 1  value: 1823.999 mgc: 0.001221812 ustep: 1 
#> iter: 2  mgc: 6.580191e-09 
#> iter: 1  mgc: 6.580191e-09 
#> outer mgc:  0.005153914 
#> iter: 1  value: 1823.999 mgc: 0.002608169 ustep: 1 
#> iter: 2  mgc: 5.509315e-09 
#> iter: 1  mgc: 5.509315e-09 
#> outer mgc:  0.004660017 
#> iter: 1  value: 1823.999 mgc: 0.0007406081 ustep: 1 
#> iter: 2  mgc: 2.093475e-09 
#> iter: 1  mgc: 2.093475e-09 
#> outer mgc:  0.00161958 
#> iter: 1  value: 1823.999 mgc: 0.0005019671 ustep: 1 
#> iter: 2  mgc: 7.919174e-10 
#> iter: 1  mgc: 7.919174e-10 
#> outer mgc:  0.001849121 
#> iter: 1  value: 1823.999 mgc: 0.0005250924 ustep: 1 
#> iter: 2  mgc: 6.837104e-10 
#> iter: 1  mgc: 6.837104e-10 
#> outer mgc:  0.001513977 
#> iter: 1  mgc: 6.837104e-10 
#> iter: 1  mgc: 6.837104e-10 
#> outer mgc:  0.001513977 
#> iter: 1  mgc: 6.837104e-10 
#> outer mgc:  0.001513977 
#> iter: 1  value: 1823.998 mgc: 0.0005963965 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 1.751766e-08 ustep: 1 
#> mgc: 1.474277e-13 
#> outer mgc:  0.02452112 
#> iter: 1  value: 1823.999 mgc: 0.0005963965 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.751628e-08 ustep: 1 
#> mgc: 1.29097e-13 
#> outer mgc:  0.02448971 
#> iter: 1  value: 1823.999 mgc: 0.02240033 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.617244e-07 ustep: 1 
#> iter: 3  mgc: 1.526483e-13 
#> outer mgc:  0.9293774 
#> iter: 1  value: 1823.999 mgc: 0.02240034 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.618802e-07 ustep: 1 
#> iter: 3  mgc: 1.246429e-13 
#> outer mgc:  0.9300413 
#> iter: 1  value: 1823.999 mgc: 0.02307942 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.599763e-06 ustep: 1 
#> iter: 3  mgc: 1.811233e-13 
#> outer mgc:  0.8468434 
#> iter: 1  value: 1823.999 mgc: 0.02307943 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.600994e-06 ustep: 1 
#> iter: 3  mgc: 1.14482e-13 
#> outer mgc:  0.8479677 
#> iter: 1  value: 1824.004 mgc: 0.0549207 ustep: 1 
#> iter: 2  value: 1824.004 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 6.876829e-13 
#> outer mgc:  1.674441 
#> iter: 1  value: 1823.996 mgc: 0.05493975 ustep: 1 
#> iter: 2  value: 1823.996 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 6.977626e-13 
#> outer mgc:  1.671133 
#> iter: 1  value: 1824.004 mgc: 0.0549207 ustep: 1 
#> iter: 2  value: 1824.004 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 6.876829e-13 
#> outer mgc:  1.600416 
#> iter: 1  value: 1823.996 mgc: 0.05493975 ustep: 1 
#> iter: 2  value: 1823.996 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 6.977626e-13 
#> outer mgc:  1.597182 
#> iter: 1  value: 1824.004 mgc: 0.0549207 ustep: 1 
#> iter: 2  value: 1824.004 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 6.776913e-13 
#> outer mgc:  1.549866 
#> iter: 1  value: 1823.996 mgc: 0.05493975 ustep: 1 
#> iter: 2  value: 1823.996 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 7.128525e-13 
#> outer mgc:  1.54649 
#> iter: 1  value: 1824.003 mgc: 0.0549207 ustep: 1 
#> iter: 2  value: 1824.003 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 7.272149e-13 
#> outer mgc:  1.507062 
#> iter: 1  value: 1823.996 mgc: 0.05493975 ustep: 1 
#> iter: 2  value: 1823.996 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 6.413206e-13 
#> outer mgc:  1.503727 
#> iter: 1  value: 1824.003 mgc: 0.0549207 ustep: 1 
#> iter: 2  value: 1824.003 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 6.930407e-13 
#> outer mgc:  1.470266 
#> iter: 1  value: 1823.996 mgc: 0.05493975 ustep: 1 
#> iter: 2  value: 1823.996 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 7.099506e-13 
#> outer mgc:  1.466814 
#> iter: 1  value: 1824.003 mgc: 0.0549207 ustep: 1 
#> iter: 2  value: 1824.003 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 7.020187e-13 
#> outer mgc:  1.43296 
#> iter: 1  value: 1823.996 mgc: 0.05493975 ustep: 1 
#> iter: 2  value: 1823.996 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 7.177857e-13 
#> outer mgc:  1.429296 
#> iter: 1  value: 1824.003 mgc: 0.0549207 ustep: 1 
#> iter: 2  value: 1824.003 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 7.160648e-13 
#> outer mgc:  1.392624 
#> iter: 1  value: 1823.997 mgc: 0.05493975 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 6.784649e-13 
#> outer mgc:  1.389152 
#> iter: 1  value: 1824.002 mgc: 0.0549207 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 6.885517e-13 
#> outer mgc:  1.344637 
#> iter: 1  value: 1823.997 mgc: 0.05493975 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 6.789002e-13 
#> outer mgc:  1.341398 
#> iter: 1  value: 1824.002 mgc: 0.0549207 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 6.464132e-13 
#> outer mgc:  1.309439 
#> iter: 1  value: 1823.997 mgc: 0.05493975 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 7.022605e-13 
#> outer mgc:  1.306184 
#> iter: 1  value: 1824.002 mgc: 0.0549207 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 6.484405e-13 
#> outer mgc:  1.277735 
#> iter: 1  value: 1823.997 mgc: 0.05493975 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 6.844138e-13 
#> outer mgc:  1.274352 
#> iter: 1  value: 1824.002 mgc: 0.0549207 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 6.926063e-13 
#> outer mgc:  1.241826 
#> iter: 1  value: 1823.997 mgc: 0.05493975 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 6.681632e-13 
#> outer mgc:  1.239 
#> iter: 1  value: 1824.002 mgc: 0.0549207 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 6.856556e-13 
#> outer mgc:  1.200885 
#> iter: 1  value: 1823.997 mgc: 0.05493975 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 7.164799e-13 
#> outer mgc:  1.197448 
#> iter: 1  value: 1824.002 mgc: 0.0549207 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 6.784153e-13 
#> outer mgc:  1.145698 
#> iter: 1  value: 1823.997 mgc: 0.05493975 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 6.54234e-13 
#> outer mgc:  1.142732 
#> iter: 1  value: 1824.001 mgc: 0.0549207 ustep: 1 
#> iter: 2  value: 1824.001 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 6.997018e-13 
#> outer mgc:  1.091362 
#> iter: 1  value: 1823.998 mgc: 0.05493975 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 7.115466e-13 
#> outer mgc:  1.088671 
#> iter: 1  value: 1824.001 mgc: 0.0549207 ustep: 1 
#> iter: 2  value: 1824.001 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 6.901446e-13 
#> outer mgc:  1.042834 
#> iter: 1  value: 1823.998 mgc: 0.05493975 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 6.864451e-13 
#> outer mgc:  1.040233 
#> iter: 1  value: 1824.001 mgc: 0.0549207 ustep: 1 
#> iter: 2  value: 1824.001 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 7.098382e-13 
#> outer mgc:  0.9930856 
#> iter: 1  value: 1823.998 mgc: 0.05493975 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 6.89347e-13 
#> outer mgc:  0.990343 
#> iter: 1  value: 1824 mgc: 0.05480459 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 7.001364e-13 
#> outer mgc:  0.9543303 
#> iter: 1  value: 1823.998 mgc: 0.05482375 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 6.648258e-13 
#> outer mgc:  0.9520287 
#> iter: 1  value: 1824 mgc: 0.04565736 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 6.895657e-13 
#> outer mgc:  0.9084202 
#> iter: 1  value: 1823.998 mgc: 0.0456655 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 6.886212e-13 
#> outer mgc:  0.9057299 
#> iter: 1  value: 1824 mgc: 0.04565736 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 7.149072e-13 
#> outer mgc:  0.8527281 
#> iter: 1  value: 1823.998 mgc: 0.0456655 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 6.934089e-13 
#> outer mgc:  0.8502493 
#> iter: 1  value: 1824 mgc: 0.04565736 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 6.688605e-13 
#> outer mgc:  0.8046 
#> iter: 1  value: 1823.999 mgc: 0.0456655 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 6.635178e-13 
#> outer mgc:  0.8022738 
#> iter: 1  value: 1824 mgc: 0.04565736 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 6.690108e-13 
#> outer mgc:  0.7580404 
#> iter: 1  value: 1823.999 mgc: 0.0456655 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 6.330427e-13 
#> outer mgc:  0.7554032 
#> iter: 1  value: 1824 mgc: 0.04565736 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 5.208209e-13 
#> outer mgc:  0.7239925 
#> iter: 1  value: 1823.999 mgc: 0.0456655 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 5.332254e-13 
#> outer mgc:  0.7221721 
#> iter: 1  value: 1824 mgc: 0.04565736 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 4.861912e-13 
#> outer mgc:  0.689842 
#> iter: 1  value: 1823.999 mgc: 0.0456655 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 5.462523e-13 
#> outer mgc:  0.6879226 
#> iter: 1  value: 1824 mgc: 0.04565736 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 5.468245e-13 
#> outer mgc:  0.6572433 
#> iter: 1  value: 1823.999 mgc: 0.0456655 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 5.146872e-13 
#> outer mgc:  0.6555897 
#> iter: 1  value: 1823.999 mgc: 0.04565736 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 5.675773e-13 
#> outer mgc:  0.6270639 
#> iter: 1  value: 1823.999 mgc: 0.0456655 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 5.478807e-13 
#> outer mgc:  0.6251849 
#> iter: 1  value: 1823.999 mgc: 0.04565736 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 5.688274e-13 
#> outer mgc:  0.5953489 
#> iter: 1  value: 1823.999 mgc: 0.0456655 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 5.429956e-13 
#> outer mgc:  0.5937548 
#> iter: 1  value: 1823.999 mgc: 0.04565736 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 5.295721e-13 
#> outer mgc:  0.5649231 
#> iter: 1  value: 1823.999 mgc: 0.0456655 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 5.159397e-13 
#> outer mgc:  0.5631785 
#> iter: 1  value: 1823.999 mgc: 0.04565736 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 5.429489e-13 
#> outer mgc:  0.5368733 
#> iter: 1  value: 1823.999 mgc: 0.0456655 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 5.16566e-13 
#> outer mgc:  0.5353996 
#> iter: 1  value: 1823.999 mgc: 0.04565736 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 5.365731e-13 
#> outer mgc:  0.5046088 
#> iter: 1  value: 1823.999 mgc: 0.0456655 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 5.16566e-13 
#> outer mgc:  0.5030921 
#> iter: 1  value: 1823.999 mgc: 0.04565736 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 5.296972e-13 
#> outer mgc:  0.4757983 
#> iter: 1  value: 1823.999 mgc: 0.0456655 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 5.190711e-13 
#> outer mgc:  0.4744888 
#> iter: 1  value: 1823.999 mgc: 0.04565736 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 5.106948e-13 
#> outer mgc:  0.442752 
#> iter: 1  value: 1823.999 mgc: 0.0456655 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 5.196972e-13 
#> outer mgc:  0.4413106 
#> iter: 1  value: 1823.999 mgc: 0.04565736 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 5.025695e-13 
#> outer mgc:  0.4114196 
#> iter: 1  value: 1823.999 mgc: 0.0456655 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 5.06419e-13 
#> outer mgc:  0.4098316 
#> iter: 1  value: 1823.999 mgc: 0.04565736 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 4.545731e-13 
#> outer mgc:  0.3801085 
#> iter: 1  value: 1823.999 mgc: 0.0456655 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 4.378929e-13 
#> outer mgc:  0.3789883 
#> iter: 1  value: 1823.999 mgc: 0.04565736 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 2.865944e-13 
#> outer mgc:  0.3516918 
#> iter: 1  value: 1823.999 mgc: 0.0456655 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 2.775254e-13 
#> outer mgc:  0.3504236 
#> iter: 1  value: 1823.999 mgc: 0.0455201 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 2.607383e-13 
#> outer mgc:  0.3210384 
#> iter: 1  value: 1823.999 mgc: 0.04552838 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 2.568064e-13 
#> outer mgc:  0.3197553 
#> iter: 1  value: 1823.999 mgc: 0.04315113 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 2.640044e-13 
#> outer mgc:  0.2889471 
#> iter: 1  value: 1823.999 mgc: 0.04314931 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 2.412672e-13 
#> outer mgc:  0.2881966 
#> iter: 1  value: 1823.999 mgc: 0.04315113 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 2.773407e-13 
#> outer mgc:  0.2569161 
#> iter: 1  value: 1823.999 mgc: 0.04314931 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 2.889753e-13 
#> outer mgc:  0.2559145 
#> iter: 1  value: 1823.999 mgc: 0.04315113 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200859e-06 ustep: 1 
#> iter: 3  mgc: 2.101149e-13 
#> outer mgc:  0.2262013 
#> iter: 1  value: 1823.999 mgc: 0.04314931 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204881e-06 ustep: 1 
#> iter: 3  mgc: 2.472646e-13 
#> outer mgc:  0.2251206 
#> iter: 1  value: 1823.999 mgc: 0.04315113 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200858e-06 ustep: 1 
#> iter: 3  mgc: 2.6455e-13 
#> outer mgc:  0.1938352 
#> iter: 1  value: 1823.999 mgc: 0.04314931 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.20488e-06 ustep: 1 
#> iter: 3  mgc: 2.824311e-13 
#> outer mgc:  0.1931947 
#> iter: 1  value: 1823.999 mgc: 0.04315113 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200841e-06 ustep: 1 
#> iter: 3  mgc: 1.896405e-13 
#> outer mgc:  0.1626156 
#> iter: 1  value: 1823.999 mgc: 0.04314931 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204863e-06 ustep: 1 
#> iter: 3  mgc: 2.042591e-13 
#> outer mgc:  0.1616746 
#> iter: 1  value: 1823.999 mgc: 0.04315113 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200621e-06 ustep: 1 
#> iter: 3  mgc: 1.796031e-13 
#> outer mgc:  0.130544 
#> iter: 1  value: 1823.999 mgc: 0.04314931 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204642e-06 ustep: 1 
#> iter: 3  mgc: 1.972826e-13 
#> outer mgc:  0.1300223 
#> iter: 1  value: 1823.999 mgc: 0.04315113 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.198923e-06 ustep: 1 
#> iter: 3  mgc: 1.796031e-13 
#> outer mgc:  0.09802413 
#> iter: 1  value: 1823.999 mgc: 0.04314931 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.202936e-06 ustep: 1 
#> iter: 3  mgc: 1.796031e-13 
#> outer mgc:  0.09771398 
#> iter: 1  value: 1823.999 mgc: 0.04315113 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.209272e-06 ustep: 1 
#> iter: 3  mgc: 1.796031e-13 
#> outer mgc:  0.06679092 
#> iter: 1  value: 1823.999 mgc: 0.04314931 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.213279e-06 ustep: 1 
#> iter: 3  mgc: 1.864297e-13 
#> outer mgc:  0.06647934 
#> iter: 1  value: 1823.999 mgc: 0.0425037 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.508467e-06 ustep: 1 
#> iter: 3  mgc: 1.796031e-13 
#> outer mgc:  0.03415435 
#> iter: 1  value: 1823.999 mgc: 0.04250246 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.512771e-06 ustep: 1 
#> iter: 3  mgc: 1.796031e-13 
#> outer mgc:  0.03415554 
#> iter: 1  value: 1823.999 mgc: 0.002082818 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 2.20398e-08 ustep: 1 
#> mgc: 1.28849e-13 
#> outer mgc:  0.009844456 
#> iter: 1  value: 1823.999 mgc: 0.002082818 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 2.203963e-08 ustep: 1 
#> mgc: 1.496118e-13 
#> outer mgc:  0.009775134 
#> iter: 1  value: 1823.999 mgc: 0.002082818 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 6.106181e-08 ustep: 1 
#> mgc: 1.28849e-13 
#> outer mgc:  0.02690621 
#> iter: 1  value: 1823.999 mgc: 0.002082818 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 6.10593e-08 ustep: 1 
#> mgc: 1.591433e-13 
#> outer mgc:  0.02689242 
#> iter: 1  value: 1823.999 mgc: 0.002082818 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 3.462897e-08 ustep: 1 
#> mgc: 1.803873e-13 
#> outer mgc:  0.06113827 
#> iter: 1  value: 1823.999 mgc: 0.002082818 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 3.462833e-08 ustep: 1 
#> mgc: 1.28849e-13 
#> outer mgc:  0.06117702 
#> iter: 1  value: 1823.999 mgc: 0.002082818 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.334945e-08 ustep: 1 
#> mgc: 1.28849e-13 
#> outer mgc:  0.0403763 
#> iter: 1  value: 1823.999 mgc: 0.002082818 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.334609e-08 ustep: 1 
#> mgc: 1.577817e-13 
#> outer mgc:  0.04011876 
#> iter: 1  value: 1823.999 mgc: 0.0003812371 ustep: 1 
#> iter: 2  mgc: 2.530021e-09 
#> outer mgc:  0.02753921 
#> iter: 1  value: 1823.999 mgc: 0.0003809271 ustep: 1 
#> iter: 2  mgc: 2.525876e-09 
#> outer mgc:  0.0275529 
#> iter: 1  value: 1823.999 mgc: 0.001701891 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.041954e-08 ustep: 1 
#> mgc: 1.591945e-13 
#> outer mgc:  0.03126143 
#> iter: 1  value: 1823.999 mgc: 0.001701581 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.039918e-08 ustep: 1 
#> mgc: 1.28849e-13 
#> outer mgc:  0.03098093 
#> iter: 1  value: 1824.001 mgc: 0.004491842 ustep: 1 
#> iter: 2  value: 1824.001 mgc: 8.35069e-07 ustep: 1 
#> iter: 3  mgc: 1.283574e-13 
#> outer mgc:  0.07725241 
#> iter: 1  value: 1823.997 mgc: 0.004500835 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 8.381584e-07 ustep: 1 
#> iter: 3  mgc: 1.007952e-13 
#> outer mgc:  0.07733868 
#> iter: 1  value: 1824.05 mgc: 0.01041406 ustep: 1 
#> iter: 2  value: 1824.05 mgc: 1.94114e-06 ustep: 1 
#> iter: 3  mgc: 1.971431e-13 
#> outer mgc:  0.261927 
#> iter: 1  value: 1823.948 mgc: 0.01043491 ustep: 1 
#> iter: 2  value: 1823.948 mgc: 1.9466e-06 ustep: 1 
#> iter: 3  mgc: 1.973671e-13 
#> outer mgc:  0.2626136 
#> iter: 1  value: 1823.998 mgc: 0.03732136 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 3.831758e-06 ustep: 1 
#> iter: 3  mgc: 1.252838e-13 
#> outer mgc:  0.1673553 
#> iter: 1  value: 1824 mgc: 0.03729913 ustep: 1 
#> iter: 2  value: 1824 mgc: 3.830337e-06 ustep: 1 
#> iter: 3  mgc: 1.490151e-13 
#> outer mgc:  0.16719 
#> iter: 1  value: 1823.998 mgc: 0.01626404 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 7.795472e-07 ustep: 1 
#> iter: 3  mgc: 1.274286e-13 
#> outer mgc:  0.08018727 
#> iter: 1  value: 1823.999 mgc: 0.01625515 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.790815e-07 ustep: 1 
#> iter: 3  mgc: 1.654977e-13 
#> outer mgc:  0.08009587 
#> iter: 1  value: 1823.999 mgc: 0.01549459 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.3409e-06 ustep: 1 
#> iter: 3  mgc: 1.147474e-13 
#> outer mgc:  0.08593942 
#> iter: 1  value: 1823.999 mgc: 0.01549902 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.342397e-06 ustep: 1 
#> iter: 3  mgc: 1.258316e-13 
#> outer mgc:  0.0871117 
#> iter: 1  value: 1823.999 mgc: 0.03268821 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.802601e-06 ustep: 1 
#> iter: 3  mgc: 1.513596e-13 
#> outer mgc:  0.2980045 
#> iter: 1  value: 1823.999 mgc: 0.03270147 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.803399e-06 ustep: 1 
#> iter: 3  mgc: 1.615154e-13 
#> outer mgc:  0.2954736 
#> iter: 1  value: 1823.999 mgc: 0.01221777 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.321694e-07 ustep: 1 
#> iter: 3  mgc: 1.991232e-13 
#> outer mgc:  0.1484056 
#> iter: 1  value: 1823.999 mgc: 0.01222341 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.327079e-07 ustep: 1 
#> iter: 3  mgc: 1.152843e-13 
#> outer mgc:  0.1491272 
#> iter: 1  value: 1823.999 mgc: 0.007496946 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.737444e-07 ustep: 1 
#> iter: 3  mgc: 1.62647e-13 
#> outer mgc:  0.09372831 
#> iter: 1  value: 1823.999 mgc: 0.007502524 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.740449e-07 ustep: 1 
#> iter: 3  mgc: 7.825048e-14 
#> outer mgc:  0.09122588 
#> iter: 1  value: 1823.999 mgc: 0.004455654 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.281382e-08 ustep: 1 
#> mgc: 1.251872e-13 
#> outer mgc:  0.05325928 
#> iter: 1  value: 1823.999 mgc: 0.004459095 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.292884e-08 ustep: 1 
#> mgc: 1.451657e-13 
#> outer mgc:  0.05396189 
#> iter: 1  value: 1823.999 mgc: 0.0001850322 ustep: 1 
#> iter: 2  mgc: 2.96752e-10 
#> outer mgc:  0.0007248966 
#> iter: 1  value: 1823.999 mgc: 0.000185213 ustep: 1 
#> iter: 2  mgc: 2.973281e-10 
#> outer mgc:  0.002478066 
#> iter: 1  value: 1823.999 mgc: 0.005196846 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.849452e-07 ustep: 1 
#> mgc: 1.631416e-13 
#> outer mgc:  0.01502116 
#> iter: 1  value: 1823.999 mgc: 0.005197538 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.850349e-07 ustep: 1 
#> mgc: 1.596924e-13 
#> outer mgc:  0.01759443 
#> iter: 1  value: 1823.999 mgc: 0.02860371 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.503674e-06 ustep: 1 
#> iter: 3  mgc: 2.605348e-13 
#> outer mgc:  0.05548238 
#> iter: 1  value: 1823.999 mgc: 0.02859277 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.508132e-06 ustep: 1 
#> iter: 3  mgc: 2.534154e-13 
#> outer mgc:  0.05292655 
#> iter: 1  value: 1823.999 mgc: 0.009316193 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.207936e-07 ustep: 1 
#> iter: 3  mgc: 1.286577e-13 
#> outer mgc:  0.0337857 
#> iter: 1  value: 1823.999 mgc: 0.009313524 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.206916e-07 ustep: 1 
#> iter: 3  mgc: 1.167274e-13 
#> outer mgc:  0.0338489 
#> iter: 1  value: 1823.999 mgc: 0.005284572 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.20923e-07 ustep: 1 
#> iter: 3  mgc: 1.613561e-13 
#> outer mgc:  0.02418199 
#> iter: 1  value: 1823.999 mgc: 0.005282106 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.207502e-07 ustep: 1 
#> iter: 3  mgc: 1.057528e-13 
#> outer mgc:  0.02430963 
#> iter: 1  value: 1823.999 mgc: 0.0005188545 ustep: 1 
#> iter: 2  mgc: 7.048389e-10 
#> iter: 1  mgc: 7.048389e-10 
#> outer mgc:  1.853157e-08 
#> iter: 1  value: 1823.998 mgc: 0.0005963968 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 1.751769e-08 ustep: 1 
#> mgc: 9.996312e-14 
#> outer mgc:  0.02450545 
#> iter: 1  value: 1824 mgc: 0.0005963975 ustep: 1 
#> iter: 2  value: 1824 mgc: 1.751632e-08 ustep: 1 
#> mgc: 1.40884e-13 
#> outer mgc:  0.0245054 
#> iter: 1  value: 1823.999 mgc: 0.02240033 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.617242e-07 ustep: 1 
#> iter: 3  mgc: 1.344393e-13 
#> outer mgc:  0.9297086 
#> iter: 1  value: 1823.999 mgc: 0.02240034 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.618802e-07 ustep: 1 
#> iter: 3  mgc: 9.095005e-14 
#> outer mgc:  0.9297101 
#> iter: 1  value: 1823.999 mgc: 0.02307942 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.599763e-06 ustep: 1 
#> iter: 3  mgc: 1.282595e-13 
#> outer mgc:  0.8474113 
#> iter: 1  value: 1823.999 mgc: 0.02307943 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.600995e-06 ustep: 1 
#> iter: 3  mgc: 2.034351e-13 
#> outer mgc:  0.8473999 
#> iter: 1  value: 1824.004 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.004 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 7.052072e-13 
#> outer mgc:  1.673155 
#> iter: 1  value: 1823.996 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.996 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.488681e-13 
#> outer mgc:  1.672417 
#> iter: 1  value: 1824.004 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.004 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 7.052072e-13 
#> outer mgc:  1.599179 
#> iter: 1  value: 1823.996 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.996 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.488681e-13 
#> outer mgc:  1.598417 
#> iter: 1  value: 1824.004 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.004 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 7.020215e-13 
#> outer mgc:  1.548532 
#> iter: 1  value: 1823.996 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.996 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.646836e-13 
#> outer mgc:  1.547822 
#> iter: 1  value: 1824.003 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.003 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.94926e-13 
#> outer mgc:  1.505729 
#> iter: 1  value: 1823.996 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.996 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.700521e-13 
#> outer mgc:  1.505059 
#> iter: 1  value: 1824.003 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.003 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 7.133164e-13 
#> outer mgc:  1.468869 
#> iter: 1  value: 1823.996 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.996 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.83546e-13 
#> outer mgc:  1.468209 
#> iter: 1  value: 1824.003 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.003 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.674127e-13 
#> outer mgc:  1.431445 
#> iter: 1  value: 1823.996 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.996 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.741148e-13 
#> outer mgc:  1.430809 
#> iter: 1  value: 1824.003 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.003 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.839207e-13 
#> outer mgc:  1.391194 
#> iter: 1  value: 1823.997 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 7.151769e-13 
#> outer mgc:  1.390581 
#> iter: 1  value: 1824.002 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 7.302587e-13 
#> outer mgc:  1.343305 
#> iter: 1  value: 1823.997 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.97185e-13 
#> outer mgc:  1.342729 
#> iter: 1  value: 1824.002 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 7.157781e-13 
#> outer mgc:  1.308091 
#> iter: 1  value: 1823.997 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.842715e-13 
#> outer mgc:  1.30753 
#> iter: 1  value: 1824.002 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.840655e-13 
#> outer mgc:  1.276321 
#> iter: 1  value: 1823.997 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 7.005222e-13 
#> outer mgc:  1.275764 
#> iter: 1  value: 1824.002 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.653855e-13 
#> outer mgc:  1.240686 
#> iter: 1  value: 1823.997 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.915263e-13 
#> outer mgc:  1.240138 
#> iter: 1  value: 1824.002 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 7.285211e-13 
#> outer mgc:  1.199428 
#> iter: 1  value: 1823.997 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 7.470979e-13 
#> outer mgc:  1.198904 
#> iter: 1  value: 1824.002 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 7.123027e-13 
#> outer mgc:  1.144443 
#> iter: 1  value: 1823.997 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.626522e-13 
#> outer mgc:  1.143986 
#> iter: 1  value: 1824.001 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.001 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.898577e-13 
#> outer mgc:  1.090205 
#> iter: 1  value: 1823.998 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 7.193846e-13 
#> outer mgc:  1.089827 
#> iter: 1  value: 1824.001 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.001 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 7.256249e-13 
#> outer mgc:  1.041701 
#> iter: 1  value: 1823.998 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 7.127102e-13 
#> outer mgc:  1.041365 
#> iter: 1  value: 1824.001 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.001 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 7.273627e-13 
#> outer mgc:  0.9918555 
#> iter: 1  value: 1823.998 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 7.064711e-13 
#> outer mgc:  0.9915724 
#> iter: 1  value: 1824 mgc: 0.0548045 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.742189e-13 
#> outer mgc:  0.9533095 
#> iter: 1  value: 1823.998 mgc: 0.05482366 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.712126e-13 
#> outer mgc:  0.9530489 
#> iter: 1  value: 1824 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.918854e-13 
#> outer mgc:  0.9071931 
#> iter: 1  value: 1823.998 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.982003e-13 
#> outer mgc:  0.9069563 
#> iter: 1  value: 1824 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.94782e-13 
#> outer mgc:  0.8515715 
#> iter: 1  value: 1823.998 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.831099e-13 
#> outer mgc:  0.8514055 
#> iter: 1  value: 1824 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.807374e-13 
#> outer mgc:  0.8034923 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.849945e-13 
#> outer mgc:  0.8033813 
#> iter: 1  value: 1824 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.358524e-13 
#> outer mgc:  0.7567484 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.433469e-13 
#> outer mgc:  0.7566949 
#> iter: 1  value: 1824 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.712833e-13 
#> outer mgc:  0.7231028 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.201482e-13 
#> outer mgc:  0.7230615 
#> iter: 1  value: 1824 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.869563e-13 
#> outer mgc:  0.6888978 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.493851e-13 
#> outer mgc:  0.6888667 
#> iter: 1  value: 1824 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.474509e-13 
#> outer mgc:  0.6564308 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.467547e-13 
#> outer mgc:  0.656402 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.444505e-13 
#> outer mgc:  0.6261387 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.771925e-13 
#> outer mgc:  0.6261099 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.126961e-13 
#> outer mgc:  0.5945653 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.434979e-13 
#> outer mgc:  0.5945382 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.768299e-13 
#> outer mgc:  0.5640639 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.323499e-13 
#> outer mgc:  0.5640375 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 4.791914e-13 
#> outer mgc:  0.5361501 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 4.926428e-13 
#> outer mgc:  0.5361225 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.224474e-13 
#> outer mgc:  0.5038627 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.51765e-13 
#> outer mgc:  0.5038381 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.573273e-13 
#> outer mgc:  0.4751546 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.411179e-13 
#> outer mgc:  0.4751323 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.254481e-13 
#> outer mgc:  0.442041 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 4.855029e-13 
#> outer mgc:  0.4420214 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.343252e-13 
#> outer mgc:  0.4106333 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.429957e-13 
#> outer mgc:  0.4106178 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 4.29195e-13 
#> outer mgc:  0.3795521 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 4.323827e-13 
#> outer mgc:  0.3795446 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 2.376052e-13 
#> outer mgc:  0.3510614 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 2.712567e-13 
#> outer mgc:  0.351054 
#> iter: 1  value: 1823.999 mgc: 0.04552017 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 2.585624e-13 
#> outer mgc:  0.3204003 
#> iter: 1  value: 1823.999 mgc: 0.04552845 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 2.908853e-13 
#> outer mgc:  0.3203934 
#> iter: 1  value: 1823.999 mgc: 0.04315111 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 2.373331e-13 
#> outer mgc:  0.2885755 
#> iter: 1  value: 1823.999 mgc: 0.04314929 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 2.478114e-13 
#> outer mgc:  0.2885682 
#> iter: 1  value: 1823.999 mgc: 0.04315111 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 2.621006e-13 
#> outer mgc:  0.2564191 
#> iter: 1  value: 1823.999 mgc: 0.04314929 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 2.570804e-13 
#> outer mgc:  0.2564115 
#> iter: 1  value: 1823.999 mgc: 0.04315111 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 2.479479e-13 
#> outer mgc:  0.225665 
#> iter: 1  value: 1823.999 mgc: 0.04314929 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 2.824339e-13 
#> outer mgc:  0.2256568 
#> iter: 1  value: 1823.999 mgc: 0.04315111 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200853e-06 ustep: 1 
#> iter: 3  mgc: 2.237256e-13 
#> outer mgc:  0.1935189 
#> iter: 1  value: 1823.999 mgc: 0.04314929 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204875e-06 ustep: 1 
#> iter: 3  mgc: 2.532625e-13 
#> outer mgc:  0.1935109 
#> iter: 1  value: 1823.999 mgc: 0.04315111 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200836e-06 ustep: 1 
#> iter: 3  mgc: 1.849966e-13 
#> outer mgc:  0.1621492 
#> iter: 1  value: 1823.999 mgc: 0.04314929 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204858e-06 ustep: 1 
#> iter: 3  mgc: 1.611744e-13 
#> outer mgc:  0.162141 
#> iter: 1  value: 1823.999 mgc: 0.04315111 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200616e-06 ustep: 1 
#> iter: 3  mgc: 1.551959e-13 
#> outer mgc:  0.1302867 
#> iter: 1  value: 1823.999 mgc: 0.04314929 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204637e-06 ustep: 1 
#> iter: 3  mgc: 1.620123e-13 
#> outer mgc:  0.1302795 
#> iter: 1  value: 1823.999 mgc: 0.04315111 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.198918e-06 ustep: 1 
#> iter: 3  mgc: 1.416501e-13 
#> outer mgc:  0.09787191 
#> iter: 1  value: 1823.999 mgc: 0.04314929 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.202931e-06 ustep: 1 
#> iter: 3  mgc: 1.321679e-13 
#> outer mgc:  0.09786615 
#> iter: 1  value: 1823.999 mgc: 0.04315111 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.209267e-06 ustep: 1 
#> iter: 3  mgc: 1.989298e-13 
#> outer mgc:  0.06663676 
#> iter: 1  value: 1823.999 mgc: 0.04314929 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.213274e-06 ustep: 1 
#> iter: 3  mgc: 1.875928e-13 
#> outer mgc:  0.06663347 
#> iter: 1  value: 1823.999 mgc: 0.04250369 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.508461e-06 ustep: 1 
#> iter: 3  mgc: 1.776487e-13 
#> outer mgc:  0.03415457 
#> iter: 1  value: 1823.999 mgc: 0.04250244 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.512765e-06 ustep: 1 
#> iter: 3  mgc: 1.775104e-13 
#> outer mgc:  0.03415532 
#> iter: 1  value: 1823.999 mgc: 0.002082815 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 2.203977e-08 ustep: 1 
#> mgc: 1.223321e-13 
#> outer mgc:  0.009813275 
#> iter: 1  value: 1823.999 mgc: 0.002082815 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 2.203954e-08 ustep: 1 
#> mgc: 9.750411e-14 
#> outer mgc:  0.009806254 
#> iter: 1  value: 1823.999 mgc: 0.002082815 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 6.106167e-08 ustep: 1 
#> mgc: 8.761404e-14 
#> outer mgc:  0.02691008 
#> iter: 1  value: 1823.999 mgc: 0.002082815 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 6.105912e-08 ustep: 1 
#> mgc: 1.002545e-13 
#> outer mgc:  0.02688856 
#> iter: 1  value: 1823.999 mgc: 0.002082815 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 3.462886e-08 ustep: 1 
#> mgc: 1.581194e-13 
#> outer mgc:  0.06118104 
#> iter: 1  value: 1823.999 mgc: 0.002082815 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 3.462823e-08 ustep: 1 
#> mgc: 1.781536e-13 
#> outer mgc:  0.06113435 
#> iter: 1  value: 1823.999 mgc: 0.002082815 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.334922e-08 ustep: 1 
#> mgc: 7.996366e-14 
#> outer mgc:  0.04026473 
#> iter: 1  value: 1823.999 mgc: 0.002082815 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.334586e-08 ustep: 1 
#> mgc: 1.350163e-13 
#> outer mgc:  0.04023014 
#> iter: 1  value: 1823.999 mgc: 0.0003812381 ustep: 1 
#> iter: 2  mgc: 2.530023e-09 
#> outer mgc:  0.02755915 
#> iter: 1  value: 1823.999 mgc: 0.0003809268 ustep: 1 
#> iter: 2  mgc: 2.525873e-09 
#> outer mgc:  0.02753301 
#> iter: 1  value: 1823.999 mgc: 0.001701888 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.041944e-08 ustep: 1 
#> mgc: 1.258522e-13 
#> outer mgc:  0.03113579 
#> iter: 1  value: 1823.999 mgc: 0.001701577 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.03989e-08 ustep: 1 
#> mgc: 1.119228e-13 
#> outer mgc:  0.03110633 
#> iter: 1  value: 1824.001 mgc: 0.004491843 ustep: 1 
#> iter: 2  value: 1824.001 mgc: 8.350706e-07 ustep: 1 
#> iter: 3  mgc: 1.619249e-13 
#> outer mgc:  0.07723292 
#> iter: 1  value: 1823.997 mgc: 0.004500836 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 8.3816e-07 ustep: 1 
#> iter: 3  mgc: 1.070924e-13 
#> outer mgc:  0.07735822 
#> iter: 1  value: 1824.05 mgc: 0.01041407 ustep: 1 
#> iter: 2  value: 1824.05 mgc: 1.941144e-06 ustep: 1 
#> iter: 3  mgc: 1.904939e-13 
#> outer mgc:  0.2621763 
#> iter: 1  value: 1823.948 mgc: 0.01043492 ustep: 1 
#> iter: 2  value: 1823.948 mgc: 1.946604e-06 ustep: 1 
#> iter: 3  mgc: 1.956775e-13 
#> outer mgc:  0.2623638 
#> iter: 1  value: 1823.998 mgc: 0.03732134 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 3.831755e-06 ustep: 1 
#> iter: 3  mgc: 1.756824e-13 
#> outer mgc:  0.1672986 
#> iter: 1  value: 1824 mgc: 0.03729911 ustep: 1 
#> iter: 2  value: 1824 mgc: 3.830334e-06 ustep: 1 
#> iter: 3  mgc: 1.179474e-13 
#> outer mgc:  0.1672466 
#> iter: 1  value: 1823.998 mgc: 0.01626405 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 7.795488e-07 ustep: 1 
#> iter: 3  mgc: 1.615644e-13 
#> outer mgc:  0.08016426 
#> iter: 1  value: 1823.999 mgc: 0.01625516 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.790832e-07 ustep: 1 
#> iter: 3  mgc: 9.610775e-14 
#> outer mgc:  0.08011886 
#> iter: 1  value: 1823.999 mgc: 0.01549462 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.340908e-06 ustep: 1 
#> iter: 3  mgc: 1.560419e-13 
#> outer mgc:  0.08650743 
#> iter: 1  value: 1823.999 mgc: 0.01549906 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.342404e-06 ustep: 1 
#> iter: 3  mgc: 1.411376e-13 
#> outer mgc:  0.08654404 
#> iter: 1  value: 1823.999 mgc: 0.03268818 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.802603e-06 ustep: 1 
#> iter: 3  mgc: 1.292489e-13 
#> outer mgc:  0.2967189 
#> iter: 1  value: 1823.999 mgc: 0.03270143 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.8034e-06 ustep: 1 
#> iter: 3  mgc: 1.829283e-13 
#> outer mgc:  0.2967591 
#> iter: 1  value: 1823.999 mgc: 0.01221777 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.321695e-07 ustep: 1 
#> iter: 3  mgc: 1.538724e-13 
#> outer mgc:  0.1487368 
#> iter: 1  value: 1823.999 mgc: 0.0122234 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.32708e-07 ustep: 1 
#> iter: 3  mgc: 1.254132e-13 
#> outer mgc:  0.148796 
#> iter: 1  value: 1823.999 mgc: 0.00749696 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.737449e-07 ustep: 1 
#> iter: 3  mgc: 1.219585e-13 
#> outer mgc:  0.09244283 
#> iter: 1  value: 1823.999 mgc: 0.007502538 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.740454e-07 ustep: 1 
#> iter: 3  mgc: 9.455804e-14 
#> outer mgc:  0.09251155 
#> iter: 1  value: 1823.999 mgc: 0.004455653 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.281374e-08 ustep: 1 
#> mgc: 1.309568e-13 
#> outer mgc:  0.05359049 
#> iter: 1  value: 1823.999 mgc: 0.004459095 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.292879e-08 ustep: 1 
#> mgc: 1.701062e-13 
#> outer mgc:  0.05363067 
#> iter: 1  value: 1823.999 mgc: 0.0001850018 ustep: 1 
#> iter: 2  mgc: 2.966691e-10 
#> outer mgc:  0.001092835 
#> iter: 1  value: 1823.999 mgc: 0.0001851827 ustep: 1 
#> iter: 2  mgc: 2.972631e-10 
#> outer mgc:  0.001093943 
#> iter: 1  value: 1823.999 mgc: 0.005196852 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.849457e-07 ustep: 1 
#> mgc: 2.00017e-13 
#> outer mgc:  0.01630672 
#> iter: 1  value: 1823.999 mgc: 0.005197544 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.850353e-07 ustep: 1 
#> mgc: 1.387958e-13 
#> outer mgc:  0.01630885 
#> iter: 1  value: 1823.999 mgc: 0.02860371 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.503677e-06 ustep: 1 
#> iter: 3  mgc: 2.489215e-13 
#> outer mgc:  0.05419683 
#> iter: 1  value: 1823.999 mgc: 0.02859278 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.508135e-06 ustep: 1 
#> iter: 3  mgc: 2.383039e-13 
#> outer mgc:  0.05421214 
#> iter: 1  value: 1823.999 mgc: 0.009316197 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.207944e-07 ustep: 1 
#> iter: 3  mgc: 1.265395e-13 
#> outer mgc:  0.0338114 
#> iter: 1  value: 1823.999 mgc: 0.009313529 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.206923e-07 ustep: 1 
#> iter: 3  mgc: 1.193826e-13 
#> outer mgc:  0.03382319 
#> iter: 1  value: 1823.999 mgc: 0.00528458 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.20925e-07 ustep: 1 
#> iter: 3  mgc: 1.69848e-13 
#> outer mgc:  0.0242449 
#> iter: 1  value: 1823.999 mgc: 0.005282114 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.207523e-07 ustep: 1 
#> iter: 3  mgc: 1.65611e-13 
#> outer mgc:  0.02424674 
#> iter: 1  value: 1823.999 mgc: 1.657263e-08 ustep: 1 
#> mgc: 1.868681e-13 
#> iter: 1  mgc: 1.868681e-13 
#> outer mgc:  1.829692e-11 
#> iter: 1  value: 1823.998 mgc: 0.0005963968 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 1.751769e-08 ustep: 1 
#> mgc: 1.426255e-13 
#> outer mgc:  0.02450545 
#> iter: 1  value: 1824 mgc: 0.0005963968 ustep: 1 
#> iter: 2  value: 1824 mgc: 1.751631e-08 ustep: 1 
#> mgc: 1.551606e-13 
#> outer mgc:  0.0245054 
#> iter: 1  value: 1823.999 mgc: 0.02240033 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.617242e-07 ustep: 1 
#> iter: 3  mgc: 1.317665e-13 
#> outer mgc:  0.9297086 
#> iter: 1  value: 1823.999 mgc: 0.02240034 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.618801e-07 ustep: 1 
#> iter: 3  mgc: 1.301361e-13 
#> outer mgc:  0.9297101 
#> iter: 1  value: 1823.999 mgc: 0.02307942 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.599763e-06 ustep: 1 
#> iter: 3  mgc: 1.321198e-13 
#> outer mgc:  0.8474113 
#> iter: 1  value: 1823.999 mgc: 0.02307943 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.600995e-06 ustep: 1 
#> iter: 3  mgc: 1.316435e-13 
#> outer mgc:  0.8473999 
#> iter: 1  value: 1824.004 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.004 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.81459e-13 
#> outer mgc:  1.673155 
#> iter: 1  value: 1823.996 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.996 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.755657e-13 
#> outer mgc:  1.672417 
#> iter: 1  value: 1824.004 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.004 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.81459e-13 
#> outer mgc:  1.599179 
#> iter: 1  value: 1823.996 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.996 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.755657e-13 
#> outer mgc:  1.598417 
#> iter: 1  value: 1824.004 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.004 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.960844e-13 
#> outer mgc:  1.548532 
#> iter: 1  value: 1823.996 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.996 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.598954e-13 
#> outer mgc:  1.547822 
#> iter: 1  value: 1824.003 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.003 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.978221e-13 
#> outer mgc:  1.505729 
#> iter: 1  value: 1823.996 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.996 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.939929e-13 
#> outer mgc:  1.505059 
#> iter: 1  value: 1824.003 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.003 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.842103e-13 
#> outer mgc:  1.468869 
#> iter: 1  value: 1823.996 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.996 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 7.060358e-13 
#> outer mgc:  1.468209 
#> iter: 1  value: 1824.003 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.003 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 7.354718e-13 
#> outer mgc:  1.431445 
#> iter: 1  value: 1823.996 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.996 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 7.286708e-13 
#> outer mgc:  1.430809 
#> iter: 1  value: 1824.003 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.003 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.742186e-13 
#> outer mgc:  1.391194 
#> iter: 1  value: 1823.997 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 7.100985e-13 
#> outer mgc:  1.390581 
#> iter: 1  value: 1824.002 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 7.073793e-13 
#> outer mgc:  1.343305 
#> iter: 1  value: 1823.997 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.736795e-13 
#> outer mgc:  1.342729 
#> iter: 1  value: 1824.002 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 7.153437e-13 
#> outer mgc:  1.308091 
#> iter: 1  value: 1823.997 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.797735e-13 
#> outer mgc:  1.30753 
#> iter: 1  value: 1824.002 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.491671e-13 
#> outer mgc:  1.276321 
#> iter: 1  value: 1823.997 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.567033e-13 
#> outer mgc:  1.275764 
#> iter: 1  value: 1824.002 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.824726e-13 
#> outer mgc:  1.240686 
#> iter: 1  value: 1823.997 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 7.00232e-13 
#> outer mgc:  1.240138 
#> iter: 1  value: 1824.002 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.520633e-13 
#> outer mgc:  1.199428 
#> iter: 1  value: 1823.997 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 7.115495e-13 
#> outer mgc:  1.198904 
#> iter: 1  value: 1824.002 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 7.267834e-13 
#> outer mgc:  1.144443 
#> iter: 1  value: 1823.997 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 7.00232e-13 
#> outer mgc:  1.143986 
#> iter: 1  value: 1824.001 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.001 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.513392e-13 
#> outer mgc:  1.090205 
#> iter: 1  value: 1823.998 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 7.33604e-13 
#> outer mgc:  1.089827 
#> iter: 1  value: 1824.001 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.001 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.94926e-13 
#> outer mgc:  1.041701 
#> iter: 1  value: 1823.998 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.780323e-13 
#> outer mgc:  1.041365 
#> iter: 1  value: 1824.001 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.001 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 7.065105e-13 
#> outer mgc:  0.9918555 
#> iter: 1  value: 1823.998 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.844165e-13 
#> outer mgc:  0.9915724 
#> iter: 1  value: 1824 mgc: 0.0548045 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.999944e-13 
#> outer mgc:  0.9533095 
#> iter: 1  value: 1823.998 mgc: 0.05482366 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.890594e-13 
#> outer mgc:  0.953049 
#> iter: 1  value: 1824 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.978225e-13 
#> outer mgc:  0.9071931 
#> iter: 1  value: 1823.998 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.725183e-13 
#> outer mgc:  0.9069563 
#> iter: 1  value: 1824 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.830527e-13 
#> outer mgc:  0.8515715 
#> iter: 1  value: 1823.998 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 7.031331e-13 
#> outer mgc:  0.8514055 
#> iter: 1  value: 1824 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.666911e-13 
#> outer mgc:  0.8034923 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 7.083549e-13 
#> outer mgc:  0.8033813 
#> iter: 1  value: 1824 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.945821e-13 
#> outer mgc:  0.7567484 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.284022e-13 
#> outer mgc:  0.7566949 
#> iter: 1  value: 1824 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.631738e-13 
#> outer mgc:  0.7231028 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.310299e-13 
#> outer mgc:  0.7230616 
#> iter: 1  value: 1824 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.780801e-13 
#> outer mgc:  0.6888978 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.126843e-13 
#> outer mgc:  0.6888667 
#> iter: 1  value: 1824 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.543268e-13 
#> outer mgc:  0.6564308 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.27089e-13 
#> outer mgc:  0.656402 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 4.916931e-13 
#> outer mgc:  0.6261387 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.456273e-13 
#> outer mgc:  0.6261099 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.044449e-13 
#> outer mgc:  0.5945653 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.184462e-13 
#> outer mgc:  0.5945382 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.151964e-13 
#> outer mgc:  0.5640639 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.379865e-13 
#> outer mgc:  0.5640375 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.219473e-13 
#> outer mgc:  0.5361501 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.507629e-13 
#> outer mgc:  0.5361225 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.211972e-13 
#> outer mgc:  0.5038627 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.356066e-13 
#> outer mgc:  0.5038381 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.439505e-13 
#> outer mgc:  0.4751546 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.342287e-13 
#> outer mgc:  0.4751323 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.193222e-13 
#> outer mgc:  0.442041 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.38362e-13 
#> outer mgc:  0.4420214 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.435765e-13 
#> outer mgc:  0.4106333 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 4.722246e-13 
#> outer mgc:  0.4106178 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 4.65951e-13 
#> outer mgc:  0.3795521 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 4.244916e-13 
#> outer mgc:  0.3795446 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 2.253575e-13 
#> outer mgc:  0.3510614 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 2.45903e-13 
#> outer mgc:  0.351054 
#> iter: 1  value: 1823.999 mgc: 0.04552017 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 2.294401e-13 
#> outer mgc:  0.3204003 
#> iter: 1  value: 1823.999 mgc: 0.04552845 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 2.688031e-13 
#> outer mgc:  0.3203934 
#> iter: 1  value: 1823.999 mgc: 0.04315111 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 2.816969e-13 
#> outer mgc:  0.2885755 
#> iter: 1  value: 1823.999 mgc: 0.04314929 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 2.42359e-13 
#> outer mgc:  0.2885682 
#> iter: 1  value: 1823.999 mgc: 0.04315111 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 2.34067e-13 
#> outer mgc:  0.2564191 
#> iter: 1  value: 1823.999 mgc: 0.04314929 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 2.42359e-13 
#> outer mgc:  0.2564115 
#> iter: 1  value: 1823.999 mgc: 0.04315111 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 2.313454e-13 
#> outer mgc:  0.225665 
#> iter: 1  value: 1823.999 mgc: 0.04314929 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 2.344529e-13 
#> outer mgc:  0.2256568 
#> iter: 1  value: 1823.999 mgc: 0.04315111 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200853e-06 ustep: 1 
#> iter: 3  mgc: 2.057622e-13 
#> outer mgc:  0.1935189 
#> iter: 1  value: 1823.999 mgc: 0.04314929 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204875e-06 ustep: 1 
#> iter: 3  mgc: 2.096436e-13 
#> outer mgc:  0.1935109 
#> iter: 1  value: 1823.999 mgc: 0.04315111 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200836e-06 ustep: 1 
#> iter: 3  mgc: 1.967922e-13 
#> outer mgc:  0.1621492 
#> iter: 1  value: 1823.999 mgc: 0.04314929 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204858e-06 ustep: 1 
#> iter: 3  mgc: 1.863973e-13 
#> outer mgc:  0.162141 
#> iter: 1  value: 1823.999 mgc: 0.04315111 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200616e-06 ustep: 1 
#> iter: 3  mgc: 1.702898e-13 
#> outer mgc:  0.1302867 
#> iter: 1  value: 1823.999 mgc: 0.04314929 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204637e-06 ustep: 1 
#> iter: 3  mgc: 1.478752e-13 
#> outer mgc:  0.1302795 
#> iter: 1  value: 1823.999 mgc: 0.04315111 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.198918e-06 ustep: 1 
#> iter: 3  mgc: 2.124752e-13 
#> outer mgc:  0.09787191 
#> iter: 1  value: 1823.999 mgc: 0.04314929 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.202931e-06 ustep: 1 
#> iter: 3  mgc: 1.478752e-13 
#> outer mgc:  0.09786615 
#> iter: 1  value: 1823.999 mgc: 0.04315111 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.209267e-06 ustep: 1 
#> iter: 3  mgc: 1.478752e-13 
#> outer mgc:  0.06663676 
#> iter: 1  value: 1823.999 mgc: 0.04314929 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.213274e-06 ustep: 1 
#> iter: 3  mgc: 1.693761e-13 
#> outer mgc:  0.06663347 
#> iter: 1  value: 1823.999 mgc: 0.04250369 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.508461e-06 ustep: 1 
#> iter: 3  mgc: 1.621673e-13 
#> outer mgc:  0.03415457 
#> iter: 1  value: 1823.999 mgc: 0.04250244 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.512765e-06 ustep: 1 
#> iter: 3  mgc: 1.972768e-13 
#> outer mgc:  0.03415532 
#> iter: 1  value: 1823.999 mgc: 0.002082815 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 2.203967e-08 ustep: 1 
#> mgc: 1.321414e-13 
#> outer mgc:  0.009813275 
#> iter: 1  value: 1823.999 mgc: 0.002082815 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 2.203958e-08 ustep: 1 
#> mgc: 1.478752e-13 
#> outer mgc:  0.009806254 
#> iter: 1  value: 1823.999 mgc: 0.002082815 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 6.106163e-08 ustep: 1 
#> mgc: 9.623622e-14 
#> outer mgc:  0.02691008 
#> iter: 1  value: 1823.999 mgc: 0.002082815 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 6.105912e-08 ustep: 1 
#> mgc: 1.478752e-13 
#> outer mgc:  0.02688856 
#> iter: 1  value: 1823.999 mgc: 0.002082815 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 3.462877e-08 ustep: 1 
#> mgc: 1.478752e-13 
#> outer mgc:  0.06118104 
#> iter: 1  value: 1823.999 mgc: 0.002082815 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 3.462823e-08 ustep: 1 
#> mgc: 1.825388e-13 
#> outer mgc:  0.06113435 
#> iter: 1  value: 1823.999 mgc: 0.002082815 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.334923e-08 ustep: 1 
#> mgc: 1.234596e-13 
#> outer mgc:  0.04026473 
#> iter: 1  value: 1823.999 mgc: 0.002082815 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.334585e-08 ustep: 1 
#> mgc: 1.478752e-13 
#> outer mgc:  0.04023014 
#> iter: 1  value: 1823.999 mgc: 0.0003812381 ustep: 1 
#> iter: 2  mgc: 2.530001e-09 
#> outer mgc:  0.02755915 
#> iter: 1  value: 1823.999 mgc: 0.0003809268 ustep: 1 
#> iter: 2  mgc: 2.52584e-09 
#> outer mgc:  0.02753301 
#> iter: 1  value: 1823.999 mgc: 0.001701888 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.041941e-08 ustep: 1 
#> mgc: 1.075649e-13 
#> outer mgc:  0.03113579 
#> iter: 1  value: 1823.999 mgc: 0.001701577 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.039894e-08 ustep: 1 
#> mgc: 1.267543e-13 
#> outer mgc:  0.03110633 
#> iter: 1  value: 1824.001 mgc: 0.004491843 ustep: 1 
#> iter: 2  value: 1824.001 mgc: 8.350706e-07 ustep: 1 
#> iter: 3  mgc: 1.462851e-13 
#> outer mgc:  0.07723292 
#> iter: 1  value: 1823.997 mgc: 0.004500836 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 8.3816e-07 ustep: 1 
#> iter: 3  mgc: 1.110568e-13 
#> outer mgc:  0.07735822 
#> iter: 1  value: 1824.05 mgc: 0.01041407 ustep: 1 
#> iter: 2  value: 1824.05 mgc: 1.941144e-06 ustep: 1 
#> iter: 3  mgc: 2.014104e-13 
#> outer mgc:  0.2621763 
#> iter: 1  value: 1823.948 mgc: 0.01043492 ustep: 1 
#> iter: 2  value: 1823.948 mgc: 1.946604e-06 ustep: 1 
#> iter: 3  mgc: 1.937397e-13 
#> outer mgc:  0.2623638 
#> iter: 1  value: 1823.998 mgc: 0.03732134 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 3.831755e-06 ustep: 1 
#> iter: 3  mgc: 1.285053e-13 
#> outer mgc:  0.1672986 
#> iter: 1  value: 1824 mgc: 0.03729911 ustep: 1 
#> iter: 2  value: 1824 mgc: 3.830335e-06 ustep: 1 
#> iter: 3  mgc: 2.511539e-13 
#> outer mgc:  0.1672466 
#> iter: 1  value: 1823.998 mgc: 0.01626405 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 7.795488e-07 ustep: 1 
#> iter: 3  mgc: 1.007051e-13 
#> outer mgc:  0.08016426 
#> iter: 1  value: 1823.999 mgc: 0.01625516 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.790832e-07 ustep: 1 
#> iter: 3  mgc: 1.160943e-13 
#> outer mgc:  0.08011886 
#> iter: 1  value: 1823.999 mgc: 0.01549462 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.340907e-06 ustep: 1 
#> iter: 3  mgc: 2.452374e-13 
#> outer mgc:  0.08650743 
#> iter: 1  value: 1823.999 mgc: 0.01549906 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.342404e-06 ustep: 1 
#> iter: 3  mgc: 1.549316e-13 
#> outer mgc:  0.08654403 
#> iter: 1  value: 1823.999 mgc: 0.03268818 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.802603e-06 ustep: 1 
#> iter: 3  mgc: 1.852903e-13 
#> outer mgc:  0.2967189 
#> iter: 1  value: 1823.999 mgc: 0.03270143 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.8034e-06 ustep: 1 
#> iter: 3  mgc: 1.5744e-13 
#> outer mgc:  0.2967592 
#> iter: 1  value: 1823.999 mgc: 0.01221777 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.321695e-07 ustep: 1 
#> iter: 3  mgc: 2.357535e-13 
#> outer mgc:  0.1487368 
#> iter: 1  value: 1823.999 mgc: 0.0122234 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.327079e-07 ustep: 1 
#> iter: 3  mgc: 1.176059e-13 
#> outer mgc:  0.148796 
#> iter: 1  value: 1823.999 mgc: 0.00749696 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.73745e-07 ustep: 1 
#> iter: 3  mgc: 1.43204e-13 
#> outer mgc:  0.09244281 
#> iter: 1  value: 1823.999 mgc: 0.007502538 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.740454e-07 ustep: 1 
#> iter: 3  mgc: 1.384875e-13 
#> outer mgc:  0.09251157 
#> iter: 1  value: 1823.999 mgc: 0.004455653 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.281372e-08 ustep: 1 
#> mgc: 1.218342e-13 
#> outer mgc:  0.0535905 
#> iter: 1  value: 1823.999 mgc: 0.004459095 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.292878e-08 ustep: 1 
#> mgc: 2.189405e-13 
#> outer mgc:  0.05363066 
#> iter: 1  value: 1823.999 mgc: 0.0001850018 ustep: 1 
#> iter: 2  mgc: 2.96685e-10 
#> outer mgc:  0.001092853 
#> iter: 1  value: 1823.999 mgc: 0.0001851826 ustep: 1 
#> iter: 2  mgc: 2.972531e-10 
#> outer mgc:  0.001093924 
#> iter: 1  value: 1823.999 mgc: 0.005196852 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.849457e-07 ustep: 1 
#> mgc: 1.249092e-13 
#> outer mgc:  0.01630674 
#> iter: 1  value: 1823.999 mgc: 0.005197544 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.850353e-07 ustep: 1 
#> mgc: 1.797382e-13 
#> outer mgc:  0.01630883 
#> iter: 1  value: 1823.999 mgc: 0.02860371 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.503677e-06 ustep: 1 
#> iter: 3  mgc: 2.845371e-13 
#> outer mgc:  0.05419681 
#> iter: 1  value: 1823.999 mgc: 0.02859278 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.508135e-06 ustep: 1 
#> iter: 3  mgc: 2.507035e-13 
#> outer mgc:  0.05421216 
#> iter: 1  value: 1823.999 mgc: 0.009316197 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.207943e-07 ustep: 1 
#> iter: 3  mgc: 1.982589e-13 
#> outer mgc:  0.0338114 
#> iter: 1  value: 1823.999 mgc: 0.009313529 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.206923e-07 ustep: 1 
#> iter: 3  mgc: 1.822912e-13 
#> outer mgc:  0.03382319 
#> iter: 1  value: 1823.999 mgc: 0.00528458 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.20925e-07 ustep: 1 
#> iter: 3  mgc: 1.238725e-13 
#> outer mgc:  0.0242449 
#> iter: 1  value: 1823.999 mgc: 0.005282114 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.207523e-07 ustep: 1 
#> iter: 3  mgc: 1.691226e-13 
#> outer mgc:  0.02424674 
#> iter: 1  mgc: 1.435393e-11 
#> iter: 1  mgc: 1.435393e-11 
#> outer mgc:  2.092992e-12 
#> iter: 1  mgc: 1.435393e-11 
#> outer mgc:  2.092992e-12 
#> iter: 1  mgc: 1.435393e-11 
#> iter: 1  mgc: 1.435393e-11 
#> iter: 1  mgc: 1.435393e-11 
#> iter: 1  mgc: 1.435393e-11 
#> Matching hessian patterns... Done
#> outer mgc:  2.092992e-12 
#> iter: 1  value: 1823.998 mgc: 0.0005963968 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 1.75177e-08 ustep: 1 
#> mgc: 1.330076e-13 
#> outer mgc:  0.02450545 
#> iter: 1  value: 1824 mgc: 0.0005963968 ustep: 1 
#> iter: 2  value: 1824 mgc: 1.751633e-08 ustep: 1 
#> mgc: 1.183669e-13 
#> outer mgc:  0.0245054 
#> iter: 1  value: 1823.999 mgc: 0.02240033 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.617242e-07 ustep: 1 
#> iter: 3  mgc: 1.524042e-13 
#> outer mgc:  0.9297086 
#> iter: 1  value: 1823.999 mgc: 0.02240034 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.618802e-07 ustep: 1 
#> iter: 3  mgc: 1.044475e-13 
#> outer mgc:  0.9297101 
#> iter: 1  value: 1823.999 mgc: 0.02307942 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.599763e-06 ustep: 1 
#> iter: 3  mgc: 1.305613e-13 
#> outer mgc:  0.8474113 
#> iter: 1  value: 1823.999 mgc: 0.02307943 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.600995e-06 ustep: 1 
#> iter: 3  mgc: 1.617181e-13 
#> outer mgc:  0.8473999 
#> iter: 1  value: 1824.004 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.004 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.955052e-13 
#> outer mgc:  1.673155 
#> iter: 1  value: 1823.996 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.996 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.861577e-13 
#> outer mgc:  1.672417 
#> iter: 1  value: 1824.004 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.004 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.955052e-13 
#> outer mgc:  1.599179 
#> iter: 1  value: 1823.996 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.996 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.861577e-13 
#> outer mgc:  1.598417 
#> iter: 1  value: 1824.004 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.004 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.936227e-13 
#> outer mgc:  1.548532 
#> iter: 1  value: 1823.996 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.996 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.793382e-13 
#> outer mgc:  1.547822 
#> iter: 1  value: 1824.003 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.003 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.99415e-13 
#> outer mgc:  1.505729 
#> iter: 1  value: 1823.996 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.996 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 7.096632e-13 
#> outer mgc:  1.505059 
#> iter: 1  value: 1824.003 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.003 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.788525e-13 
#> outer mgc:  1.468869 
#> iter: 1  value: 1823.996 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.996 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 7.019732e-13 
#> outer mgc:  1.468209 
#> iter: 1  value: 1824.003 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.003 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 7.267834e-13 
#> outer mgc:  1.431445 
#> iter: 1  value: 1823.996 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.996 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 7.026986e-13 
#> outer mgc:  1.430809 
#> iter: 1  value: 1824.003 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.003 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.60462e-13 
#> outer mgc:  1.391194 
#> iter: 1  value: 1823.997 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.680208e-13 
#> outer mgc:  1.390581 
#> iter: 1  value: 1824.002 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.510496e-13 
#> outer mgc:  1.343305 
#> iter: 1  value: 1823.997 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.556876e-13 
#> outer mgc:  1.342729 
#> iter: 1  value: 1824.002 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.64227e-13 
#> outer mgc:  1.308091 
#> iter: 1  value: 1823.997 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.778873e-13 
#> outer mgc:  1.30753 
#> iter: 1  value: 1824.002 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.9565e-13 
#> outer mgc:  1.276321 
#> iter: 1  value: 1823.997 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 7.006673e-13 
#> outer mgc:  1.275764 
#> iter: 1  value: 1824.002 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.598828e-13 
#> outer mgc:  1.240686 
#> iter: 1  value: 1823.997 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 7.231571e-13 
#> outer mgc:  1.240138 
#> iter: 1  value: 1824.002 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 7.105651e-13 
#> outer mgc:  1.199428 
#> iter: 1  value: 1823.997 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 7.025535e-13 
#> outer mgc:  1.198904 
#> iter: 1  value: 1824.002 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.002 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 7.037592e-13 
#> outer mgc:  1.144443 
#> iter: 1  value: 1823.997 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 7.257689e-13 
#> outer mgc:  1.143986 
#> iter: 1  value: 1824.001 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.001 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 7.09841e-13 
#> outer mgc:  1.090205 
#> iter: 1  value: 1823.998 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.983458e-13 
#> outer mgc:  1.089827 
#> iter: 1  value: 1824.001 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.001 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.885545e-13 
#> outer mgc:  1.041701 
#> iter: 1  value: 1823.998 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.947184e-13 
#> outer mgc:  1.041365 
#> iter: 1  value: 1824.001 mgc: 0.05492061 ustep: 1 
#> iter: 2  value: 1824.001 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.94926e-13 
#> outer mgc:  0.9918555 
#> iter: 1  value: 1823.998 mgc: 0.05493966 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.730991e-13 
#> outer mgc:  0.9915724 
#> iter: 1  value: 1824 mgc: 0.0548045 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 7.26494e-13 
#> outer mgc:  0.9533095 
#> iter: 1  value: 1823.998 mgc: 0.05482366 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.968946e-13 
#> outer mgc:  0.953049 
#> iter: 1  value: 1824 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.953608e-13 
#> outer mgc:  0.9071931 
#> iter: 1  value: 1823.998 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 7.169176e-13 
#> outer mgc:  0.9069563 
#> iter: 1  value: 1824 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.869624e-13 
#> outer mgc:  0.8515715 
#> iter: 1  value: 1823.998 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.967489e-13 
#> outer mgc:  0.8514055 
#> iter: 1  value: 1824 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.732074e-13 
#> outer mgc:  0.8034923 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.784652e-13 
#> outer mgc:  0.8033813 
#> iter: 1  value: 1824 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 6.571392e-13 
#> outer mgc:  0.7567484 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 6.456684e-13 
#> outer mgc:  0.7566949 
#> iter: 1  value: 1824 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.213232e-13 
#> outer mgc:  0.7231028 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.523582e-13 
#> outer mgc:  0.7230616 
#> iter: 1  value: 1824 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.098207e-13 
#> outer mgc:  0.6888978 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.386128e-13 
#> outer mgc:  0.6888667 
#> iter: 1  value: 1824 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1824 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.544518e-13 
#> outer mgc:  0.6564308 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.272143e-13 
#> outer mgc:  0.656402 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.323238e-13 
#> outer mgc:  0.6261387 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.614099e-13 
#> outer mgc:  0.6261099 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.537017e-13 
#> outer mgc:  0.5945653 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.262122e-13 
#> outer mgc:  0.5945382 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.7733e-13 
#> outer mgc:  0.5640639 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 4.910145e-13 
#> outer mgc:  0.5640375 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.208222e-13 
#> outer mgc:  0.5361501 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.563996e-13 
#> outer mgc:  0.5361225 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.514514e-13 
#> outer mgc:  0.5038627 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.194482e-13 
#> outer mgc:  0.5038381 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.000694e-13 
#> outer mgc:  0.4751546 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.744367e-13 
#> outer mgc:  0.4751323 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.324491e-13 
#> outer mgc:  0.442041 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.434977e-13 
#> outer mgc:  0.4420214 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 5.166977e-13 
#> outer mgc:  0.4106333 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 5.371085e-13 
#> outer mgc:  0.4106178 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 4.586999e-13 
#> outer mgc:  0.3795521 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 3.993152e-13 
#> outer mgc:  0.3795446 
#> iter: 1  value: 1823.999 mgc: 0.04565743 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 2.45226e-13 
#> outer mgc:  0.3510614 
#> iter: 1  value: 1823.999 mgc: 0.04566557 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 2.671674e-13 
#> outer mgc:  0.351054 
#> iter: 1  value: 1823.999 mgc: 0.04552017 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 2.422321e-13 
#> outer mgc:  0.3204003 
#> iter: 1  value: 1823.999 mgc: 0.04552845 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 2.606245e-13 
#> outer mgc:  0.3203934 
#> iter: 1  value: 1823.999 mgc: 0.04315111 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 2.743483e-13 
#> outer mgc:  0.2885755 
#> iter: 1  value: 1823.999 mgc: 0.04314929 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 2.655316e-13 
#> outer mgc:  0.2885682 
#> iter: 1  value: 1823.999 mgc: 0.04315111 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 2.305288e-13 
#> outer mgc:  0.2564191 
#> iter: 1  value: 1823.999 mgc: 0.04314929 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 2.548995e-13 
#> outer mgc:  0.2564115 
#> iter: 1  value: 1823.999 mgc: 0.04315111 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200854e-06 ustep: 1 
#> iter: 3  mgc: 2.601956e-13 
#> outer mgc:  0.225665 
#> iter: 1  value: 1823.999 mgc: 0.04314929 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204876e-06 ustep: 1 
#> iter: 3  mgc: 2.679851e-13 
#> outer mgc:  0.2256568 
#> iter: 1  value: 1823.999 mgc: 0.04315111 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200853e-06 ustep: 1 
#> iter: 3  mgc: 2.536646e-13 
#> outer mgc:  0.1935189 
#> iter: 1  value: 1823.999 mgc: 0.04314929 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204875e-06 ustep: 1 
#> iter: 3  mgc: 2.137328e-13 
#> outer mgc:  0.1935109 
#> iter: 1  value: 1823.999 mgc: 0.04315111 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200836e-06 ustep: 1 
#> iter: 3  mgc: 1.950592e-13 
#> outer mgc:  0.1621492 
#> iter: 1  value: 1823.999 mgc: 0.04314929 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204858e-06 ustep: 1 
#> iter: 3  mgc: 1.913668e-13 
#> outer mgc:  0.162141 
#> iter: 1  value: 1823.999 mgc: 0.04315111 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.200616e-06 ustep: 1 
#> iter: 3  mgc: 1.625493e-13 
#> outer mgc:  0.1302867 
#> iter: 1  value: 1823.999 mgc: 0.04314929 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.204637e-06 ustep: 1 
#> iter: 3  mgc: 2.038719e-13 
#> outer mgc:  0.1302795 
#> iter: 1  value: 1823.999 mgc: 0.04315111 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.198918e-06 ustep: 1 
#> iter: 3  mgc: 1.5597e-13 
#> outer mgc:  0.09787191 
#> iter: 1  value: 1823.999 mgc: 0.04314929 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.202931e-06 ustep: 1 
#> iter: 3  mgc: 1.49997e-13 
#> outer mgc:  0.09786615 
#> iter: 1  value: 1823.999 mgc: 0.04315111 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.209267e-06 ustep: 1 
#> iter: 3  mgc: 1.571314e-13 
#> outer mgc:  0.06663676 
#> iter: 1  value: 1823.999 mgc: 0.04314929 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.213274e-06 ustep: 1 
#> iter: 3  mgc: 1.550354e-13 
#> outer mgc:  0.06663347 
#> iter: 1  value: 1823.999 mgc: 0.04250369 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.508461e-06 ustep: 1 
#> iter: 3  mgc: 1.900338e-13 
#> outer mgc:  0.03415457 
#> iter: 1  value: 1823.999 mgc: 0.04250244 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.512765e-06 ustep: 1 
#> iter: 3  mgc: 1.565812e-13 
#> outer mgc:  0.03415532 
#> iter: 1  value: 1823.999 mgc: 0.002082815 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 2.203978e-08 ustep: 1 
#> mgc: 9.952664e-14 
#> outer mgc:  0.009813275 
#> iter: 1  value: 1823.999 mgc: 0.002082815 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 2.203947e-08 ustep: 1 
#> mgc: 1.434724e-13 
#> outer mgc:  0.009806254 
#> iter: 1  value: 1823.999 mgc: 0.002082815 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 6.106164e-08 ustep: 1 
#> mgc: 1.866952e-13 
#> outer mgc:  0.02691008 
#> iter: 1  value: 1823.999 mgc: 0.002082815 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 6.105912e-08 ustep: 1 
#> mgc: 1.330432e-13 
#> outer mgc:  0.02688856 
#> iter: 1  value: 1823.999 mgc: 0.002082815 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 3.462882e-08 ustep: 1 
#> mgc: 1.037414e-13 
#> outer mgc:  0.06118104 
#> iter: 1  value: 1823.999 mgc: 0.002082815 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 3.462819e-08 ustep: 1 
#> mgc: 2.277364e-13 
#> outer mgc:  0.06113435 
#> iter: 1  value: 1823.999 mgc: 0.002082815 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.334924e-08 ustep: 1 
#> mgc: 1.161309e-13 
#> outer mgc:  0.04026473 
#> iter: 1  value: 1823.999 mgc: 0.002082815 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.334582e-08 ustep: 1 
#> mgc: 9.066273e-14 
#> outer mgc:  0.04023014 
#> iter: 1  value: 1823.999 mgc: 0.0003812382 ustep: 1 
#> iter: 2  mgc: 2.52998e-09 
#> outer mgc:  0.02755915 
#> iter: 1  value: 1823.999 mgc: 0.0003809268 ustep: 1 
#> iter: 2  mgc: 2.525815e-09 
#> outer mgc:  0.02753301 
#> iter: 1  value: 1823.999 mgc: 0.001701888 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.041935e-08 ustep: 1 
#> mgc: 9.228726e-14 
#> outer mgc:  0.03113579 
#> iter: 1  value: 1823.999 mgc: 0.001701577 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.039888e-08 ustep: 1 
#> mgc: 1.739145e-13 
#> outer mgc:  0.03110633 
#> iter: 1  value: 1824.001 mgc: 0.004491843 ustep: 1 
#> iter: 2  value: 1824.001 mgc: 8.350706e-07 ustep: 1 
#> iter: 3  mgc: 1.136406e-13 
#> outer mgc:  0.07723292 
#> iter: 1  value: 1823.997 mgc: 0.004500836 ustep: 1 
#> iter: 2  value: 1823.997 mgc: 8.3816e-07 ustep: 1 
#> iter: 3  mgc: 1.71975e-13 
#> outer mgc:  0.07735822 
#> iter: 1  value: 1824.05 mgc: 0.01041407 ustep: 1 
#> iter: 2  value: 1824.05 mgc: 1.941144e-06 ustep: 1 
#> iter: 3  mgc: 1.863257e-13 
#> outer mgc:  0.2621763 
#> iter: 1  value: 1823.948 mgc: 0.01043492 ustep: 1 
#> iter: 2  value: 1823.948 mgc: 1.946604e-06 ustep: 1 
#> iter: 3  mgc: 1.982614e-13 
#> outer mgc:  0.2623638 
#> iter: 1  value: 1823.998 mgc: 0.03732134 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 3.831755e-06 ustep: 1 
#> iter: 3  mgc: 1.413181e-13 
#> outer mgc:  0.1672986 
#> iter: 1  value: 1824 mgc: 0.03729911 ustep: 1 
#> iter: 2  value: 1824 mgc: 3.830334e-06 ustep: 1 
#> iter: 3  mgc: 2.047165e-13 
#> outer mgc:  0.1672466 
#> iter: 1  value: 1823.998 mgc: 0.01626405 ustep: 1 
#> iter: 2  value: 1823.998 mgc: 7.795487e-07 ustep: 1 
#> iter: 3  mgc: 1.218607e-13 
#> outer mgc:  0.08016426 
#> iter: 1  value: 1823.999 mgc: 0.01625516 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.790833e-07 ustep: 1 
#> iter: 3  mgc: 1.282274e-13 
#> outer mgc:  0.08011886 
#> iter: 1  value: 1823.999 mgc: 0.01549462 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.340907e-06 ustep: 1 
#> iter: 3  mgc: 1.045547e-13 
#> outer mgc:  0.08650743 
#> iter: 1  value: 1823.999 mgc: 0.01549906 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.342404e-06 ustep: 1 
#> iter: 3  mgc: 1.266333e-13 
#> outer mgc:  0.08654403 
#> iter: 1  value: 1823.999 mgc: 0.03268818 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.802603e-06 ustep: 1 
#> iter: 3  mgc: 1.662701e-13 
#> outer mgc:  0.2967189 
#> iter: 1  value: 1823.999 mgc: 0.03270143 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.8034e-06 ustep: 1 
#> iter: 3  mgc: 1.557287e-13 
#> outer mgc:  0.2967592 
#> iter: 1  value: 1823.999 mgc: 0.01221777 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.321694e-07 ustep: 1 
#> iter: 3  mgc: 1.278293e-13 
#> outer mgc:  0.1487368 
#> iter: 1  value: 1823.999 mgc: 0.0122234 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.32708e-07 ustep: 1 
#> iter: 3  mgc: 1.096548e-13 
#> outer mgc:  0.148796 
#> iter: 1  value: 1823.999 mgc: 0.00749696 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.73745e-07 ustep: 1 
#> iter: 3  mgc: 1.325905e-13 
#> outer mgc:  0.09244281 
#> iter: 1  value: 1823.999 mgc: 0.007502538 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.740454e-07 ustep: 1 
#> iter: 3  mgc: 1.065266e-13 
#> outer mgc:  0.09251157 
#> iter: 1  value: 1823.999 mgc: 0.004455653 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.281368e-08 ustep: 1 
#> mgc: 1.356779e-13 
#> outer mgc:  0.0535905 
#> iter: 1  value: 1823.999 mgc: 0.004459095 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.292878e-08 ustep: 1 
#> mgc: 1.11083e-13 
#> outer mgc:  0.05363066 
#> iter: 1  value: 1823.999 mgc: 0.0001850018 ustep: 1 
#> iter: 2  mgc: 2.96667e-10 
#> outer mgc:  0.001092853 
#> iter: 1  value: 1823.999 mgc: 0.0001851826 ustep: 1 
#> iter: 2  mgc: 2.972671e-10 
#> outer mgc:  0.001093924 
#> iter: 1  value: 1823.999 mgc: 0.005196852 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.849457e-07 ustep: 1 
#> mgc: 1.430126e-13 
#> outer mgc:  0.01630674 
#> iter: 1  value: 1823.999 mgc: 0.005197544 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 1.850353e-07 ustep: 1 
#> mgc: 1.253846e-13 
#> outer mgc:  0.01630883 
#> iter: 1  value: 1823.999 mgc: 0.02860371 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.503677e-06 ustep: 1 
#> iter: 3  mgc: 2.915053e-13 
#> outer mgc:  0.05419681 
#> iter: 1  value: 1823.999 mgc: 0.02859278 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 5.508135e-06 ustep: 1 
#> iter: 3  mgc: 2.344291e-13 
#> outer mgc:  0.05421216 
#> iter: 1  value: 1823.999 mgc: 0.009316197 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.207944e-07 ustep: 1 
#> iter: 3  mgc: 1.110513e-13 
#> outer mgc:  0.0338114 
#> iter: 1  value: 1823.999 mgc: 0.009313529 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 7.206923e-07 ustep: 1 
#> iter: 3  mgc: 1.084194e-13 
#> outer mgc:  0.03382319 
#> iter: 1  value: 1823.999 mgc: 0.00528458 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.20925e-07 ustep: 1 
#> iter: 3  mgc: 1.30977e-13 
#> outer mgc:  0.0242449 
#> iter: 1  value: 1823.999 mgc: 0.005282114 ustep: 1 
#> iter: 2  value: 1823.999 mgc: 4.207523e-07 ustep: 1 
#> iter: 3  mgc: 1.101087e-13 
#> outer mgc:  0.02424674 
#> outer mgc:  1 
#> iter: 1  mgc: 1.435393e-11 
mod <- make_osa_residuals(mod)
#> Doing OSA residuals...
#> Doing OSA for aggregate catch and index observations...
#> Loading required namespace: parallel
#> [1] 132
#> [1] 131
#> [1] 130
#> [1] 129
#> [1] 128
#> [1] 127
#> [1] 126
#> [1] 125
#> [1] 124
#> [1] 123
#> [1] 122
#> [1] 121
#> [1] 120
#> [1] 119
#> [1] 118
#> [1] 117
#> [1] 116
#> [1] 115
#> [1] 114
#> [1] 113
#> [1] 112
#> [1] 111
#> [1] 110
#> [1] 109
#> [1] 108
#> [1] 107
#> [1] 106
#> [1] 105
#> [1] 104
#> [1] 103
#> [1] 102
#> [1] 101
#> [1] 100
#> [1] 99
#> [1] 98
#> [1] 97
#> [1] 96
#> [1] 95
#> [1] 94
#> [1] 93
#> [1] 92
#> [1] 91
#> [1] 90
#> [1] 89
#> [1] 88
#> [1] 87
#> [1] 86
#> [1] 85
#> [1] 84
#> [1] 83
#> [1] 82
#> [1] 81
#> [1] 80
#> [1] 79
#> [1] 78
#> [1] 77
#> [1] 76
#> [1] 75
#> [1] 74
#> [1] 73
#> [1] 72
#> [1] 71
#> [1] 70
#> [1] 69
#> [1] 68
#> [1] 67
#> [1] 66
#> [1] 65
#> [1] 64
#> [1] 63
#> [1] 62
#> [1] 61
#> [1] 60
#> [1] 59
#> [1] 58
#> [1] 57
#> [1] 56
#> [1] 55
#> [1] 54
#> [1] 53
#> [1] 52
#> [1] 51
#> [1] 50
#> [1] 49
#> [1] 48
#> [1] 47
#> [1] 46
#> [1] 45
#> [1] 44
#> [1] 43
#> [1] 42
#> [1] 41
#> [1] 40
#> [1] 39
#> [1] 38
#> [1] 37
#> [1] 36
#> [1] 35
#> [1] 34
#> [1] 33
#> [1] 32
#> [1] 31
#> [1] 30
#> [1] 29
#> [1] 28
#> [1] 27
#> [1] 26
#> [1] 25
#> [1] 24
#> [1] 23
#> [1] 22
#> [1] 21
#> [1] 20
#> [1] 19
#> [1] 18
#> [1] 17
#> [1] 16
#> [1] 15
#> [1] 14
#> [1] 13
#> [1] 12
#> [1] 11
#> [1] 10
#> [1] 9
#> [1] 8
#> [1] 7
#> [1] 6
#> [1] 5
#> [1] 4
#> [1] 3
#> [1] 2
#> [1] 1
#> Doing OSA for catch and index age comp observations...
#> Doing OSA for age comp observations...
#> [1] 792
#> [1] 791
#> [1] 790
#> [1] 789
#> [1] 788
#> [1] 787
#> [1] 786
#> [1] 785
#> [1] 784
#> [1] 783
#> [1] 782
#> [1] 781
#> [1] 780
#> [1] 779
#> [1] 778
#> [1] 777
#> [1] 776
#> [1] 775
#> [1] 774
#> [1] 773
#> [1] 772
#> [1] 771
#> [1] 770
#> [1] 769
#> [1] 768
#> [1] 767
#> [1] 766
#> [1] 765
#> [1] 764
#> [1] 763
#> [1] 762
#> [1] 761
#> [1] 760
#> [1] 759
#> [1] 758
#> [1] 757
#> [1] 756
#> [1] 755
#> [1] 754
#> [1] 753
#> [1] 752
#> [1] 751
#> [1] 750
#> [1] 749
#> [1] 748
#> [1] 747
#> [1] 746
#> [1] 745
#> [1] 744
#> [1] 743
#> [1] 742
#> [1] 741
#> [1] 740
#> [1] 739
#> [1] 738
#> [1] 737
#> [1] 736
#> [1] 735
#> [1] 734
#> [1] 733
#> [1] 732
#> [1] 731
#> [1] 730
#> [1] 729
#> [1] 728
#> [1] 727
#> [1] 726
#> [1] 725
#> [1] 724
#> [1] 723
#> [1] 722
#> [1] 721
#> [1] 720
#> [1] 719
#> [1] 718
#> [1] 717
#> [1] 716
#> [1] 715
#> [1] 714
#> [1] 713
#> [1] 712
#> [1] 711
#> [1] 710
#> [1] 709
#> [1] 708
#> [1] 707
#> [1] 706
#> [1] 705
#> [1] 704
#> [1] 703
#> [1] 702
#> [1] 701
#> [1] 700
#> [1] 699
#> [1] 698
#> [1] 697
#> [1] 696
#> [1] 695
#> [1] 694
#> [1] 693
#> [1] 692
#> [1] 691
#> [1] 690
#> [1] 689
#> [1] 688
#> [1] 687
#> [1] 686
#> [1] 685
#> [1] 684
#> [1] 683
#> [1] 682
#> [1] 681
#> [1] 680
#> [1] 679
#> [1] 678
#> [1] 677
#> [1] 676
#> [1] 675
#> [1] 674
#> [1] 673
#> [1] 672
#> [1] 671
#> [1] 670
#> [1] 669
#> [1] 668
#> [1] 667
#> [1] 666
#> [1] 665
#> [1] 664
#> [1] 663
#> [1] 662
#> [1] 661
#> [1] 660
#> [1] 659
#> [1] 658
#> [1] 657
#> [1] 656
#> [1] 655
#> [1] 654
#> [1] 653
#> [1] 652
#> [1] 651
#> [1] 650
#> [1] 649
#> [1] 648
#> [1] 647
#> [1] 646
#> [1] 645
#> [1] 644
#> [1] 643
#> [1] 642
#> [1] 641
#> [1] 640
#> [1] 639
#> [1] 638
#> [1] 637
#> [1] 636
#> [1] 635
#> [1] 634
#> [1] 633
#> [1] 632
#> [1] 631
#> [1] 630
#> [1] 629
#> [1] 628
#> [1] 627
#> [1] 626
#> [1] 625
#> [1] 624
#> [1] 623
#> [1] 622
#> [1] 621
#> [1] 620
#> [1] 619
#> [1] 618
#> [1] 617
#> [1] 616
#> [1] 615
#> [1] 614
#> [1] 613
#> [1] 612
#> [1] 611
#> [1] 610
#> [1] 609
#> [1] 608
#> [1] 607
#> [1] 606
#> [1] 605
#> [1] 604
#> [1] 603
#> [1] 602
#> [1] 601
#> [1] 600
#> [1] 599
#> [1] 598
#> [1] 597
#> [1] 596
#> [1] 595
#> [1] 594
#> [1] 593
#> [1] 592
#> [1] 591
#> [1] 590
#> [1] 589
#> [1] 588
#> [1] 587
#> [1] 586
#> [1] 585
#> [1] 584
#> [1] 583
#> [1] 582
#> [1] 581
#> [1] 580
#> [1] 579
#> [1] 578
#> [1] 577
#> [1] 576
#> [1] 575
#> [1] 574
#> [1] 573
#> [1] 572
#> [1] 571
#> [1] 570
#> [1] 569
#> [1] 568
#> [1] 567
#> [1] 566
#> [1] 565
#> [1] 564
#> [1] 563
#> [1] 562
#> [1] 561
#> [1] 560
#> [1] 559
#> [1] 558
#> [1] 557
#> [1] 556
#> [1] 555
#> [1] 554
#> [1] 553
#> [1] 552
#> [1] 551
#> [1] 550
#> [1] 549
#> [1] 548
#> [1] 547
#> [1] 546
#> [1] 545
#> [1] 544
#> [1] 543
#> [1] 542
#> [1] 541
#> [1] 540
#> [1] 539
#> [1] 538
#> [1] 537
#> [1] 536
#> [1] 535
#> [1] 534
#> [1] 533
#> [1] 532
#> [1] 531
#> [1] 530
#> [1] 529
#> [1] 528
#> [1] 527
#> [1] 526
#> [1] 525
#> [1] 524
#> [1] 523
#> [1] 522
#> [1] 521
#> [1] 520
#> [1] 519
#> [1] 518
#> [1] 517
#> [1] 516
#> [1] 515
#> [1] 514
#> [1] 513
#> [1] 512
#> [1] 511
#> [1] 510
#> [1] 509
#> [1] 508
#> [1] 507
#> [1] 506
#> [1] 505
#> [1] 504
#> [1] 503
#> [1] 502
#> [1] 501
#> [1] 500
#> [1] 499
#> [1] 498
#> [1] 497
#> [1] 496
#> [1] 495
#> [1] 494
#> [1] 493
#> [1] 492
#> [1] 491
#> [1] 490
#> [1] 489
#> [1] 488
#> [1] 487
#> [1] 486
#> [1] 485
#> [1] 484
#> [1] 483
#> [1] 482
#> [1] 481
#> [1] 480
#> [1] 479
#> [1] 478
#> [1] 477
#> [1] 476
#> [1] 475
#> [1] 474
#> [1] 473
#> [1] 472
#> [1] 471
#> [1] 470
#> [1] 469
#> [1] 468
#> [1] 467
#> [1] 466
#> [1] 465
#> [1] 464
#> [1] 463
#> [1] 462
#> [1] 461
#> [1] 460
#> [1] 459
#> [1] 458
#> [1] 457
#> [1] 456
#> [1] 455
#> [1] 454
#> [1] 453
#> [1] 452
#> [1] 451
#> [1] 450
#> [1] 449
#> [1] 448
#> [1] 447
#> [1] 446
#> [1] 445
#> [1] 444
#> [1] 443
#> [1] 442
#> [1] 441
#> [1] 440
#> [1] 439
#> [1] 438
#> [1] 437
#> [1] 436
#> [1] 435
#> [1] 434
#> [1] 433
#> [1] 432
#> [1] 431
#> [1] 430
#> [1] 429
#> [1] 428
#> [1] 427
#> [1] 426
#> [1] 425
#> [1] 424
#> [1] 423
#> [1] 422
#> [1] 421
#> [1] 420
#> [1] 419
#> [1] 418
#> [1] 417
#> [1] 416
#> [1] 415
#> [1] 414
#> [1] 413
#> [1] 412
#> [1] 411
#> [1] 410
#> [1] 409
#> [1] 408
#> [1] 407
#> [1] 406
#> [1] 405
#> [1] 404
#> [1] 403
#> [1] 402
#> [1] 401
#> [1] 400
#> [1] 399
#> [1] 398
#> [1] 397
#> [1] 396
#> [1] 395
#> [1] 394
#> [1] 393
#> [1] 392
#> [1] 391
#> [1] 390
#> [1] 389
#> [1] 388
#> [1] 387
#> [1] 386
#> [1] 385
#> [1] 384
#> [1] 383
#> [1] 382
#> [1] 381
#> [1] 380
#> [1] 379
#> [1] 378
#> [1] 377
#> [1] 376
#> [1] 375
#> [1] 374
#> [1] 373
#> [1] 372
#> [1] 371
#> [1] 370
#> [1] 369
#> [1] 368
#> [1] 367
#> [1] 366
#> [1] 365
#> [1] 364
#> [1] 363
#> [1] 362
#> [1] 361
#> [1] 360
#> [1] 359
#> [1] 358
#> [1] 357
#> [1] 356
#> [1] 355
#> [1] 354
#> [1] 353
#> [1] 352
#> [1] 351
#> [1] 350
#> [1] 349
#> [1] 348
#> [1] 347
#> [1] 346
#> [1] 345
#> [1] 344
#> [1] 343
#> [1] 342
#> [1] 341
#> [1] 340
#> [1] 339
#> [1] 338
#> [1] 337
#> [1] 336
#> [1] 335
#> [1] 334
#> [1] 333
#> [1] 332
#> [1] 331
#> [1] 330
#> [1] 329
#> [1] 328
#> [1] 327
#> [1] 326
#> [1] 325
#> [1] 324
#> [1] 323
#> [1] 322
#> [1] 321
#> [1] 320
#> [1] 319
#> [1] 318
#> [1] 317
#> [1] 316
#> [1] 315
#> [1] 314
#> [1] 313
#> [1] 312
#> [1] 311
#> [1] 310
#> [1] 309
#> [1] 308
#> [1] 307
#> [1] 306
#> [1] 305
#> [1] 304
#> [1] 303
#> [1] 302
#> [1] 301
#> [1] 300
#> [1] 299
#> [1] 298
#> [1] 297
#> [1] 296
#> [1] 295
#> [1] 294
#> [1] 293
#> [1] 292
#> [1] 291
#> [1] 290
#> [1] 289
#> [1] 288
#> [1] 287
#> [1] 286
#> [1] 285
#> [1] 284
#> [1] 283
#> [1] 282
#> [1] 281
#> [1] 280
#> [1] 279
#> [1] 278
#> [1] 277
#> [1] 276
#> [1] 275
#> [1] 274
#> [1] 273
#> [1] 272
#> [1] 271
#> [1] 270
#> [1] 269
#> [1] 268
#> [1] 267
#> [1] 266
#> [1] 265
#> [1] 264
#> [1] 263
#> [1] 262
#> [1] 261
#> [1] 260
#> [1] 259
#> [1] 258
#> [1] 257
#> [1] 256
#> [1] 255
#> [1] 254
#> [1] 253
#> [1] 252
#> [1] 251
#> [1] 250
#> [1] 249
#> [1] 248
#> [1] 247
#> [1] 246
#> [1] 245
#> [1] 244
#> [1] 243
#> [1] 242
#> [1] 241
#> [1] 240
#> [1] 239
#> [1] 238
#> [1] 237
#> [1] 236
#> [1] 235
#> [1] 234
#> [1] 233
#> [1] 232
#> [1] 231
#> [1] 230
#> [1] 229
#> [1] 228
#> [1] 227
#> [1] 226
#> [1] 225
#> [1] 224
#> [1] 223
#> [1] 222
#> [1] 221
#> [1] 220
#> [1] 219
#> [1] 218
#> [1] 217
#> [1] 216
#> [1] 215
#> [1] 214
#> [1] 213
#> [1] 212
#> [1] 211
#> [1] 210
#> [1] 209
#> [1] 208
#> [1] 207
#> [1] 206
#> [1] 205
#> [1] 204
#> [1] 203
#> [1] 202
#> [1] 201
#> [1] 200
#> [1] 199
#> [1] 198
#> [1] 197
#> [1] 196
#> [1] 195
#> [1] 194
#> [1] 193
#> [1] 192
#> [1] 191
#> [1] 190
#> [1] 189
#> [1] 188
#> [1] 187
#> [1] 186
#> [1] 185
#> [1] 184
#> [1] 183
#> [1] 182
#> [1] 181
#> [1] 180
#> [1] 179
#> [1] 178
#> [1] 177
#> [1] 176
#> [1] 175
#> [1] 174
#> [1] 173
#> [1] 172
#> [1] 171
#> [1] 170
#> [1] 169
#> [1] 168
#> [1] 167
#> [1] 166
#> [1] 165
#> [1] 164
#> [1] 163
#> [1] 162
#> [1] 161
#> [1] 160
#> [1] 159
#> [1] 158
#> [1] 157
#> [1] 156
#> [1] 155
#> [1] 154
#> [1] 153
#> [1] 152
#> [1] 151
#> [1] 150
#> [1] 149
#> [1] 148
#> [1] 147
#> [1] 146
#> [1] 145
#> [1] 144
#> [1] 143
#> [1] 142
#> [1] 141
#> [1] 140
#> [1] 139
#> [1] 138
#> [1] 137
#> [1] 136
#> [1] 135
#> [1] 134
#> [1] 133
#> [1] 132
#> [1] 131
#> [1] 130
#> [1] 129
#> [1] 128
#> [1] 127
#> [1] 126
#> [1] 125
#> [1] 124
#> [1] 123
#> [1] 122
#> [1] 121
#> [1] 120
#> [1] 119
#> [1] 118
#> [1] 117
#> [1] 116
#> [1] 115
#> [1] 114
#> [1] 113
#> [1] 112
#> [1] 111
#> [1] 110
#> [1] 109
#> [1] 108
#> [1] 107
#> [1] 106
#> [1] 105
#> [1] 104
#> [1] 103
#> [1] 102
#> [1] 101
#> [1] 100
#> [1] 99
#> [1] 98
#> [1] 97
#> [1] 96
#> [1] 95
#> [1] 94
#> [1] 93
#> [1] 92
#> [1] 91
#> [1] 90
#> [1] 89
#> [1] 88
#> [1] 87
#> [1] 86
#> [1] 85
#> [1] 84
#> [1] 83
#> [1] 82
#> [1] 81
#> [1] 80
#> [1] 79
#> [1] 78
#> [1] 77
#> [1] 76
#> [1] 75
#> [1] 74
#> [1] 73
#> [1] 72
#> [1] 71
#> [1] 70
#> [1] 69
#> [1] 68
#> [1] 67
#> [1] 66
#> [1] 65
#> [1] 64
#> [1] 63
#> [1] 62
#> [1] 61
#> [1] 60
#> [1] 59
#> [1] 58
#> [1] 57
#> [1] 56
#> [1] 55
#> [1] 54
#> [1] 53
#> [1] 52
#> [1] 51
#> [1] 50
#> [1] 49
#> [1] 48
#> [1] 47
#> [1] 46
#> [1] 45
#> [1] 44
#> [1] 43
#> [1] 42
#> [1] 41
#> [1] 40
#> [1] 39
#> [1] 38
#> [1] 37
#> [1] 36
#> [1] 35
#> [1] 34
#> [1] 33
#> [1] 32
#> [1] 31
#> [1] 30
#> [1] 29
#> [1] 28
#> [1] 27
#> [1] 26
#> [1] 25
#> [1] 24
#> [1] 23
#> [1] 22
#> [1] 21
#> [1] 20
#> [1] 19
#> [1] 18
#> [1] 17
#> [1] 16
#> [1] 15
#> [1] 14
#> [1] 13
#> [1] 12
#> [1] 11
#> [1] 10
#> [1] 9
#> [1] 8
#> [1] 7
#> [1] 6
#> [1] 5
#> [1] 4
#> [1] 3
#> [1] 2
#> [1] 1
plot_wham_output(mod, dir.main = tempdir(check = TRUE))
#> Generating plot files... Please wait ~30 seconds...
#> html output works best with Google Chrome browser
#> Making HTML document.
#> Opening HTML document in your default web-browser.
# }
```
