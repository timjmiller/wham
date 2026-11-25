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
mod <- fit_wham(input, do.osa = FALSE, do.retro = TRUE)
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
#> Retro Peel: 1
#> Using stats::nlminb for optimization with these control parameters:
#> iter.max: 1000
#> eval.max: 1000
#> iter: 1  value: 1766.986 mgc: 2.184558 ustep: 1 
#> iter: 2  value: 1766.979 mgc: 0.2280424 ustep: 1 
#> iter: 3  value: 1766.979 mgc: 0.004926529 ustep: 1 
#> iter: 4  value: 1766.979 mgc: 2.438648e-06 ustep: 1 
#> iter: 5  mgc: 6.143575e-13 
#> iter: 1  mgc: 6.143575e-13 
#> Matching hessian patterns... Done
#> outer mgc:  9.528688 
#> iter: 1  value: 111388.9 mgc: 414.3604 ustep: 1 
#> iter: 2  value: 59609.96 mgc: 6.516144e+14 ustep: 1 
#> iter: 3  value: 59236.59 mgc: 6.516144e+14 ustep: 0.008923143 
#> iter: 4  value: 58408.16 mgc: 6.516144e+14 ustep: 0.09455295 
#> iter: 5  value: 57055.94 mgc: 6.516144e+14 ustep: 0.3075639 
#> iter: 6  value: 57055.94 mgc: 6.516144e+14 ustep: 2.931278e-12 
#> Error in newton(par = c(log_NAA = 9.83305460908798, log_NAA = 11.17239355608,  : 
#>   Newton failed to find minimum.
#> Warning: NA/NaN function evaluation
#> iter: 1  value: 1856.071 mgc: 20.98766 ustep: 1 
#> iter: 2  value: 1855.899 mgc: 0.7322968 ustep: 1 
#> iter: 3  value: 1855.899 mgc: 0.0133447 ustep: 1 
#> iter: 4  value: 1855.899 mgc: 5.814163e-06 ustep: 1 
#> iter: 5  mgc: 1.990791e-12 
#> iter: 1  value: 1767.595 mgc: 2.204095 ustep: 1 
#> iter: 2  value: 1767.595 mgc: 0.00911527 ustep: 1 
#> iter: 3  value: 1767.595 mgc: 2.286575e-06 ustep: 1 
#> iter: 4  mgc: 3.138341e-13 
#> iter: 1  value: 1766.958 mgc: 0.4190731 ustep: 1 
#> iter: 2  value: 1766.958 mgc: 0.0003340103 ustep: 1 
#> iter: 3  mgc: 3.235548e-09 
#> iter: 1  mgc: 3.235548e-09 
#> outer mgc:  9.424451 
#> iter: 1  value: 1766.913 mgc: 0.2527061 ustep: 1 
#> iter: 2  value: 1766.913 mgc: 0.0001146605 ustep: 1 
#> iter: 3  mgc: 8.016803e-10 
#> iter: 1  mgc: 8.016803e-10 
#> outer mgc:  9.415443 
#> iter: 1  value: 1766.868 mgc: 0.4872252 ustep: 1 
#> iter: 2  value: 1766.868 mgc: 0.0005225452 ustep: 1 
#> iter: 3  mgc: 7.428955e-09 
#> iter: 1  mgc: 7.428955e-09 
#> outer mgc:  9.198791 
#> iter: 1  value: 1766.717 mgc: 0.4755746 ustep: 1 
#> iter: 2  value: 1766.717 mgc: 0.0004659157 ustep: 1 
#> iter: 3  value: 1766.717 mgc: 1.068987e-08 ustep: 1 
#> mgc: 1.488155e-13 
#> iter: 1  mgc: 1.488155e-13 
#> outer mgc:  8.858308 
#> iter: 1  value: 1766.558 mgc: 0.7427654 ustep: 1 
#> iter: 2  value: 1766.558 mgc: 0.0003810082 ustep: 1 
#> iter: 3  mgc: 5.981427e-09 
#> iter: 1  mgc: 5.981427e-09 
#> outer mgc:  8.075582 
#> iter: 1  value: 1766.332 mgc: 0.734432 ustep: 1 
#> iter: 2  value: 1766.332 mgc: 0.0006245926 ustep: 1 
#> iter: 3  value: 1766.332 mgc: 3.606838e-08 ustep: 1 
#> iter: 4  mgc: 1.187889e-13 
#> iter: 1  mgc: 1.187889e-13 
#> outer mgc:  6.956993 
#> iter: 1  value: 1766.332 mgc: 0.8195664 ustep: 1 
#> iter: 2  value: 1766.332 mgc: 0.000873629 ustep: 1 
#> iter: 3  value: 1766.332 mgc: 2.470658e-08 ustep: 1 
#> mgc: 1.72928e-13 
#> iter: 1  mgc: 1.72928e-13 
#> outer mgc:  14.5483 
#> iter: 1  value: 1766.226 mgc: 1.263281 ustep: 1 
#> iter: 2  value: 1766.226 mgc: 0.004724337 ustep: 1 
#> iter: 3  value: 1766.226 mgc: 6.101892e-07 ustep: 1 
#> iter: 4  mgc: 1.400155e-13 
#> iter: 1  mgc: 1.400155e-13 
#> outer mgc:  9.866229 
#> iter: 1  value: 1766.118 mgc: 0.8568123 ustep: 1 
#> iter: 2  value: 1766.118 mgc: 0.001160133 ustep: 1 
#> iter: 3  value: 1766.118 mgc: 5.37793e-08 ustep: 1 
#> iter: 4  mgc: 1.216902e-13 
#> iter: 1  mgc: 1.216902e-13 
#> outer mgc:  5.037795 
#> iter: 1  value: 1765.955 mgc: 1.035154 ustep: 1 
#> iter: 2  value: 1765.955 mgc: 0.000935144 ustep: 1 
#> iter: 3  value: 1765.955 mgc: 5.801286e-08 ustep: 1 
#> iter: 4  mgc: 1.197012e-13 
#> iter: 1  mgc: 1.197012e-13 
#> outer mgc:  11.1547 
#> iter: 1  value: 1765.963 mgc: 0.9886176 ustep: 1 
#> iter: 2  value: 1765.963 mgc: 0.001600483 ustep: 1 
#> iter: 3  value: 1765.963 mgc: 9.576748e-08 ustep: 1 
#> iter: 4  mgc: 1.101101e-13 
#> iter: 1  mgc: 1.101101e-13 
#> outer mgc:  15.52213 
#> iter: 1  value: 1765.795 mgc: 0.8473629 ustep: 1 
#> iter: 2  value: 1765.795 mgc: 0.001974146 ustep: 1 
#> iter: 3  value: 1765.795 mgc: 1.069479e-07 ustep: 1 
#> iter: 4  mgc: 1.51312e-13 
#> iter: 1  mgc: 1.51312e-13 
#> outer mgc:  4.218777 
#> iter: 1  value: 1765.776 mgc: 0.4269459 ustep: 1 
#> iter: 2  value: 1765.776 mgc: 0.0002450382 ustep: 1 
#> iter: 3  mgc: 1.637616e-09 
#> iter: 1  mgc: 1.637616e-09 
#> outer mgc:  3.74796 
#> iter: 1  value: 1765.744 mgc: 0.2804007 ustep: 1 
#> iter: 2  value: 1765.744 mgc: 0.0003045478 ustep: 1 
#> iter: 3  mgc: 2.537718e-09 
#> iter: 1  mgc: 2.537718e-09 
#> outer mgc:  3.782393 
#> iter: 1  value: 1765.763 mgc: 0.5196797 ustep: 1 
#> iter: 2  value: 1765.763 mgc: 0.0003700313 ustep: 1 
#> iter: 3  mgc: 6.765612e-09 
#> iter: 1  mgc: 6.765612e-09 
#> outer mgc:  7.663253 
#> iter: 1  value: 1765.716 mgc: 0.4991434 ustep: 1 
#> iter: 2  value: 1765.716 mgc: 0.0004012688 ustep: 1 
#> iter: 3  mgc: 5.13732e-09 
#> iter: 1  mgc: 5.13732e-09 
#> outer mgc:  3.176274 
#> iter: 1  value: 1765.675 mgc: 0.3739928 ustep: 1 
#> iter: 2  value: 1765.675 mgc: 0.0003023854 ustep: 1 
#> iter: 3  mgc: 2.356349e-09 
#> iter: 1  mgc: 2.356349e-09 
#> outer mgc:  3.77326 
#> iter: 1  value: 1765.574 mgc: 0.4032009 ustep: 1 
#> iter: 2  value: 1765.574 mgc: 0.0006625033 ustep: 1 
#> iter: 3  value: 1765.574 mgc: 1.220841e-08 ustep: 1 
#> mgc: 1.627367e-13 
#> iter: 1  mgc: 1.627367e-13 
#> outer mgc:  4.356212 
#> iter: 1  value: 1765.609 mgc: 0.4867177 ustep: 1 
#> iter: 2  value: 1765.609 mgc: 0.0003903856 ustep: 1 
#> iter: 3  mgc: 4.918723e-09 
#> iter: 1  mgc: 4.918723e-09 
#> outer mgc:  6.462967 
#> iter: 1  value: 1765.698 mgc: 0.7083293 ustep: 1 
#> iter: 2  value: 1765.698 mgc: 0.001445863 ustep: 1 
#> iter: 3  value: 1765.698 mgc: 8.742725e-08 ustep: 1 
#> iter: 4  mgc: 1.228125e-13 
#> iter: 1  mgc: 1.228125e-13 
#> outer mgc:  5.645805 
#> iter: 1  value: 1765.658 mgc: 0.5277898 ustep: 1 
#> iter: 2  value: 1765.658 mgc: 0.0003983283 ustep: 1 
#> iter: 3  mgc: 3.547584e-09 
#> iter: 1  mgc: 3.547584e-09 
#> outer mgc:  4.46758 
#> iter: 1  value: 1765.553 mgc: 0.3844711 ustep: 1 
#> iter: 2  value: 1765.553 mgc: 0.0002245129 ustep: 1 
#> iter: 3  mgc: 6.261734e-09 
#> iter: 1  mgc: 6.261734e-09 
#> outer mgc:  4.184503 
#> iter: 1  value: 1765.495 mgc: 0.3374441 ustep: 1 
#> iter: 2  value: 1765.495 mgc: 0.0002134897 ustep: 1 
#> iter: 3  mgc: 3.970326e-10 
#> iter: 1  mgc: 3.970326e-10 
#> outer mgc:  3.519053 
#> iter: 1  value: 1765.549 mgc: 0.2499822 ustep: 1 
#> iter: 2  value: 1765.549 mgc: 0.0001430123 ustep: 1 
#> iter: 3  mgc: 7.496264e-10 
#> iter: 1  mgc: 7.496264e-10 
#> outer mgc:  1.645751 
#> iter: 1  value: 1765.516 mgc: 0.1632596 ustep: 1 
#> iter: 2  value: 1765.516 mgc: 4.026521e-05 ustep: 1 
#> iter: 3  mgc: 6.363648e-11 
#> iter: 1  mgc: 6.363648e-11 
#> outer mgc:  2.25326 
#> iter: 1  value: 1765.416 mgc: 0.6269917 ustep: 1 
#> iter: 2  value: 1765.416 mgc: 0.0006022789 ustep: 1 
#> iter: 3  value: 1765.416 mgc: 2.392096e-08 ustep: 1 
#> mgc: 1.053212e-13 
#> iter: 1  mgc: 1.053212e-13 
#> outer mgc:  3.178952 
#> iter: 1  value: 1765.54 mgc: 0.3675012 ustep: 1 
#> iter: 2  value: 1765.54 mgc: 0.00020896 ustep: 1 
#> iter: 3  mgc: 6.926366e-10 
#> iter: 1  mgc: 6.926366e-10 
#> outer mgc:  3.383963 
#> iter: 1  value: 1765.353 mgc: 0.1785147 ustep: 1 
#> iter: 2  value: 1765.353 mgc: 0.0003250767 ustep: 1 
#> iter: 3  mgc: 2.427416e-09 
#> iter: 1  mgc: 2.427416e-09 
#> outer mgc:  1.569301 
#> iter: 1  value: 1765.48 mgc: 0.3418484 ustep: 1 
#> iter: 2  value: 1765.48 mgc: 0.0002865162 ustep: 1 
#> iter: 3  mgc: 2.393264e-09 
#> iter: 1  mgc: 2.393264e-09 
#> outer mgc:  3.243142 
#> iter: 1  value: 1765.534 mgc: 0.8203853 ustep: 1 
#> iter: 2  value: 1765.534 mgc: 0.00148365 ustep: 1 
#> iter: 3  value: 1765.534 mgc: 1.391784e-07 ustep: 1 
#> iter: 4  mgc: 2.375949e-13 
#> iter: 1  value: 1765.487 mgc: 0.4861612 ustep: 1 
#> iter: 2  value: 1765.487 mgc: 0.0004361664 ustep: 1 
#> iter: 3  value: 1765.487 mgc: 1.138362e-08 ustep: 1 
#> mgc: 1.681643e-13 
#> iter: 1  value: 1765.468 mgc: 0.2513764 ustep: 1 
#> iter: 2  value: 1765.468 mgc: 0.0001045781 ustep: 1 
#> iter: 3  mgc: 4.719875e-10 
#> iter: 1  value: 1765.471 mgc: 0.1144857 ustep: 1 
#> iter: 2  value: 1765.471 mgc: 2.16501e-05 ustep: 1 
#> iter: 3  mgc: 2.017153e-11 
#> iter: 1  mgc: 2.017153e-11 
#> outer mgc:  2.331709 
#> iter: 1  value: 1765.468 mgc: 0.06926102 ustep: 1 
#> iter: 2  value: 1765.468 mgc: 1.075016e-05 ustep: 1 
#> iter: 3  mgc: 2.720258e-12 
#> iter: 1  mgc: 2.720258e-12 
#> outer mgc:  2.164195 
#> iter: 1  value: 1765.451 mgc: 0.1167452 ustep: 1 
#> iter: 2  value: 1765.451 mgc: 2.490192e-05 ustep: 1 
#> iter: 3  mgc: 3.643474e-11 
#> iter: 1  mgc: 3.643474e-11 
#> outer mgc:  2.332353 
#> iter: 1  value: 1765.433 mgc: 0.1157988 ustep: 1 
#> iter: 2  value: 1765.433 mgc: 1.892502e-05 ustep: 1 
#> iter: 3  mgc: 6.140134e-12 
#> iter: 1  mgc: 6.140134e-12 
#> outer mgc:  1.7519 
#> iter: 1  value: 1765.5 mgc: 0.07797266 ustep: 1 
#> iter: 2  value: 1765.5 mgc: 5.007906e-05 ustep: 1 
#> iter: 3  mgc: 6.899689e-11 
#> iter: 1  mgc: 6.899689e-11 
#> outer mgc:  1.398223 
#> iter: 1  value: 1765.287 mgc: 0.06208525 ustep: 1 
#> iter: 2  value: 1765.287 mgc: 7.465155e-05 ustep: 1 
#> iter: 3  mgc: 1.931737e-10 
#> iter: 1  mgc: 1.931737e-10 
#> outer mgc:  0.8698238 
#> iter: 1  value: 1765.556 mgc: 1.006424 ustep: 1 
#> iter: 2  value: 1765.556 mgc: 0.002110858 ustep: 1 
#> iter: 3  value: 1765.556 mgc: 8.726807e-08 ustep: 1 
#> iter: 4  mgc: 1.163e-13 
#> iter: 1  value: 1765.3 mgc: 0.1319088 ustep: 1 
#> iter: 2  value: 1765.3 mgc: 3.912118e-05 ustep: 1 
#> iter: 3  mgc: 3.657821e-11 
#> iter: 1  value: 1765.291 mgc: 0.05229008 ustep: 1 
#> iter: 2  value: 1765.291 mgc: 6.151298e-06 ustep: 1 
#> iter: 3  mgc: 9.15599e-13 
#> iter: 1  mgc: 9.15599e-13 
#> outer mgc:  1.075854 
#> iter: 1  value: 1765.291 mgc: 0.05026656 ustep: 1 
#> iter: 2  value: 1765.291 mgc: 4.495113e-06 ustep: 1 
#> iter: 3  mgc: 1.089943e-12 
#> iter: 1  mgc: 1.089943e-12 
#> outer mgc:  0.7504038 
#> iter: 1  value: 1765.297 mgc: 0.03896931 ustep: 1 
#> iter: 2  value: 1765.297 mgc: 2.585859e-06 ustep: 1 
#> iter: 3  mgc: 1.561382e-13 
#> iter: 1  mgc: 1.561382e-13 
#> outer mgc:  0.8560606 
#> iter: 1  value: 1765.3 mgc: 0.05088056 ustep: 1 
#> iter: 2  value: 1765.3 mgc: 9.537238e-06 ustep: 1 
#> iter: 3  mgc: 3.533026e-12 
#> iter: 1  mgc: 3.533026e-12 
#> outer mgc:  0.6602805 
#> iter: 1  value: 1765.303 mgc: 0.03534248 ustep: 1 
#> iter: 2  value: 1765.303 mgc: 5.274714e-06 ustep: 1 
#> iter: 3  mgc: 6.457352e-13 
#> iter: 1  value: 1765.303 mgc: 0.03520403 ustep: 1 
#> iter: 2  value: 1765.303 mgc: 2.63186e-05 ustep: 1 
#> iter: 3  mgc: 1.941113e-11 
#> iter: 1  value: 1765.303 mgc: 0.07766006 ustep: 1 
#> iter: 2  value: 1765.303 mgc: 0.0001284821 ustep: 1 
#> iter: 3  mgc: 4.632348e-10 
#> iter: 1  mgc: 4.632348e-10 
#> outer mgc:  0.8170882 
#> iter: 1  value: 1765.375 mgc: 0.1343158 ustep: 1 
#> iter: 2  value: 1765.375 mgc: 0.0003993496 ustep: 1 
#> iter: 3  mgc: 3.720684e-09 
#> iter: 1  mgc: 3.720684e-09 
#> outer mgc:  1.139541 
#> iter: 1  value: 1765.448 mgc: 0.3057516 ustep: 1 
#> iter: 2  value: 1765.448 mgc: 0.0007345471 ustep: 1 
#> iter: 3  value: 1765.448 mgc: 1.401094e-08 ustep: 1 
#> mgc: 1.692246e-13 
#> iter: 1  value: 1765.413 mgc: 0.2077391 ustep: 1 
#> iter: 2  value: 1765.413 mgc: 0.0001580496 ustep: 1 
#> iter: 3  mgc: 5.42146e-10 
#> iter: 1  value: 1765.398 mgc: 0.1765244 ustep: 1 
#> iter: 2  value: 1765.398 mgc: 6.219135e-05 ustep: 1 
#> iter: 3  mgc: 9.394932e-11 
#> iter: 1  value: 1765.389 mgc: 0.1610626 ustep: 1 
#> iter: 2  value: 1765.389 mgc: 4.764453e-05 ustep: 1 
#> iter: 3  mgc: 3.313516e-11 
#> iter: 1  value: 1765.385 mgc: 0.1675172 ustep: 1 
#> iter: 2  value: 1765.385 mgc: 5.118527e-05 ustep: 1 
#> iter: 3  mgc: 4.657351e-11 
#> iter: 1  value: 1765.378 mgc: 0.07760322 ustep: 1 
#> iter: 2  value: 1765.378 mgc: 1.126071e-05 ustep: 1 
#> iter: 3  mgc: 3.414649e-12 
#> iter: 1  mgc: 3.414649e-12 
#> outer mgc:  1.241751 
#> iter: 1  value: 1765.375 mgc: 0.04924157 ustep: 1 
#> iter: 2  value: 1765.375 mgc: 4.447551e-06 ustep: 1 
#> iter: 3  mgc: 1.100974e-12 
#> iter: 1  mgc: 1.100974e-12 
#> outer mgc:  0.5266185 
#> iter: 1  value: 1765.377 mgc: 0.05085016 ustep: 1 
#> iter: 2  value: 1765.377 mgc: 4.601111e-06 ustep: 1 
#> iter: 3  mgc: 2.481176e-13 
#> iter: 1  mgc: 2.481176e-13 
#> outer mgc:  1.047724 
#> iter: 1  value: 1765.377 mgc: 0.04878314 ustep: 1 
#> iter: 2  value: 1765.377 mgc: 1.29154e-05 ustep: 1 
#> iter: 3  mgc: 4.654844e-12 
#> iter: 1  mgc: 4.654844e-12 
#> outer mgc:  0.5315693 
#> iter: 1  value: 1765.382 mgc: 0.03243681 ustep: 1 
#> iter: 2  value: 1765.382 mgc: 8.223577e-06 ustep: 1 
#> iter: 3  mgc: 1.937135e-12 
#> iter: 1  value: 1765.388 mgc: 0.03458879 ustep: 1 
#> iter: 2  value: 1765.388 mgc: 2.591981e-05 ustep: 1 
#> iter: 3  mgc: 1.580596e-11 
#> iter: 1  value: 1765.416 mgc: 0.158649 ustep: 1 
#> iter: 2  value: 1765.416 mgc: 0.0005497434 ustep: 1 
#> iter: 3  mgc: 7.10663e-09 
#> iter: 1  value: 1765.503 mgc: 0.4588634 ustep: 1 
#> iter: 2  value: 1765.503 mgc: 0.00469284 ustep: 1 
#> iter: 3  value: 1765.503 mgc: 5.19299e-07 ustep: 1 
#> iter: 4  mgc: 1.400829e-13 
#> iter: 1  mgc: 1.400829e-13 
#> outer mgc:  1.157058 
#> iter: 1  value: 1765.373 mgc: 0.6699892 ustep: 1 
#> iter: 2  value: 1765.373 mgc: 0.009678618 ustep: 1 
#> iter: 3  value: 1765.373 mgc: 2.252336e-06 ustep: 1 
#> iter: 4  mgc: 1.339052e-13 
#> iter: 1  value: 1765.43 mgc: 0.3053284 ustep: 1 
#> iter: 2  value: 1765.43 mgc: 0.001085217 ustep: 1 
#> iter: 3  value: 1765.43 mgc: 2.840396e-08 ustep: 1 
#> mgc: 1.361472e-13 
#> iter: 1  mgc: 1.361472e-13 
#> outer mgc:  2.525733 
#> iter: 1  value: 1765.391 mgc: 0.347906 ustep: 1 
#> iter: 2  value: 1765.391 mgc: 0.0002763969 ustep: 1 
#> iter: 3  mgc: 2.234315e-09 
#> iter: 1  value: 1765.417 mgc: 0.1648232 ustep: 1 
#> iter: 2  value: 1765.417 mgc: 4.492798e-05 ustep: 1 
#> iter: 3  mgc: 2.856253e-11 
#> iter: 1  mgc: 2.856253e-11 
#> outer mgc:  1.215701 
#> iter: 1  value: 1765.401 mgc: 0.2125464 ustep: 1 
#> iter: 2  value: 1765.401 mgc: 0.000199875 ustep: 1 
#> iter: 3  mgc: 1.286264e-09 
#> iter: 1  mgc: 1.286264e-09 
#> outer mgc:  0.4627719 
#> iter: 1  value: 1765.384 mgc: 0.1517581 ustep: 1 
#> iter: 2  value: 1765.384 mgc: 0.0001215166 ustep: 1 
#> iter: 3  mgc: 3.968132e-10 
#> iter: 1  mgc: 3.968132e-10 
#> outer mgc:  1.272037 
#> iter: 1  value: 1765.365 mgc: 0.07583412 ustep: 1 
#> iter: 2  value: 1765.365 mgc: 0.0001060447 ustep: 1 
#> iter: 3  mgc: 2.51384e-10 
#> iter: 1  mgc: 2.51384e-10 
#> outer mgc:  0.7290172 
#> iter: 1  value: 1765.366 mgc: 0.1510197 ustep: 1 
#> iter: 2  value: 1765.366 mgc: 0.0001239162 ustep: 1 
#> iter: 3  mgc: 3.559399e-10 
#> iter: 1  mgc: 3.559399e-10 
#> outer mgc:  0.9772766 
#> iter: 1  value: 1765.368 mgc: 0.08640158 ustep: 1 
#> iter: 2  value: 1765.368 mgc: 1.642214e-05 ustep: 1 
#> iter: 3  mgc: 5.218376e-12 
#> iter: 1  mgc: 5.218376e-12 
#> outer mgc:  1.32024 
#> iter: 1  value: 1765.356 mgc: 0.06894081 ustep: 1 
#> iter: 2  value: 1765.356 mgc: 4.098269e-05 ustep: 1 
#> iter: 3  mgc: 3.889259e-11 
#> iter: 1  mgc: 3.889259e-11 
#> outer mgc:  0.2849028 
#> iter: 1  value: 1765.353 mgc: 0.03340888 ustep: 1 
#> iter: 2  value: 1765.353 mgc: 1.695824e-05 ustep: 1 
#> iter: 3  mgc: 8.531075e-12 
#> iter: 1  mgc: 8.531075e-12 
#> outer mgc:  0.4639085 
#> iter: 1  value: 1765.338 mgc: 0.07346127 ustep: 1 
#> iter: 2  value: 1765.338 mgc: 0.0001061226 ustep: 1 
#> iter: 3  mgc: 2.693655e-10 
#> iter: 1  mgc: 2.693655e-10 
#> outer mgc:  0.3321509 
#> iter: 1  value: 1765.324 mgc: 0.1932226 ustep: 1 
#> iter: 2  value: 1765.324 mgc: 0.0003228569 ustep: 1 
#> iter: 3  mgc: 2.933532e-09 
#> iter: 1  value: 1765.333 mgc: 0.1621402 ustep: 1 
#> iter: 2  value: 1765.333 mgc: 4.541596e-05 ustep: 1 
#> iter: 3  mgc: 7.05263e-11 
#> iter: 1  value: 1765.337 mgc: 0.01620925 ustep: 1 
#> iter: 2  value: 1765.337 mgc: 4.534244e-07 ustep: 1 
#> iter: 3  mgc: 1.281735e-13 
#> iter: 1  mgc: 1.281735e-13 
#> outer mgc:  0.3039353 
#> iter: 1  value: 1765.337 mgc: 0.01035995 ustep: 1 
#> iter: 2  value: 1765.337 mgc: 1.791824e-07 ustep: 1 
#> iter: 3  mgc: 1.250063e-13 
#> iter: 1  mgc: 1.250063e-13 
#> outer mgc:  0.2755773 
#> iter: 1  value: 1765.336 mgc: 0.01105899 ustep: 1 
#> iter: 2  value: 1765.336 mgc: 4.067341e-07 ustep: 1 
#> iter: 3  mgc: 1.221945e-13 
#> iter: 1  mgc: 1.221945e-13 
#> outer mgc:  0.2555108 
#> iter: 1  value: 1765.335 mgc: 0.009876275 ustep: 1 
#> iter: 2  value: 1765.335 mgc: 1.466727e-07 ustep: 1 
#> iter: 3  mgc: 1.385917e-13 
#> iter: 1  mgc: 1.385917e-13 
#> outer mgc:  0.2435127 
#> iter: 1  value: 1765.334 mgc: 0.008268048 ustep: 1 
#> iter: 2  value: 1765.334 mgc: 1.13493e-06 ustep: 1 
#> iter: 3  mgc: 1.504422e-13 
#> iter: 1  mgc: 1.504422e-13 
#> outer mgc:  0.2394328 
#> iter: 1  value: 1765.334 mgc: 0.01067326 ustep: 1 
#> iter: 2  value: 1765.334 mgc: 1.735254e-06 ustep: 1 
#> iter: 3  mgc: 1.685418e-13 
#> iter: 1  value: 1765.333 mgc: 0.01647968 ustep: 1 
#> iter: 2  value: 1765.333 mgc: 5.912265e-06 ustep: 1 
#> iter: 3  mgc: 9.555008e-13 
#> iter: 1  value: 1765.331 mgc: 0.05889145 ustep: 1 
#> iter: 2  value: 1765.331 mgc: 7.569638e-05 ustep: 1 
#> iter: 3  mgc: 1.583031e-10 
#> iter: 1  mgc: 1.583031e-10 
#> outer mgc:  0.3734858 
#> iter: 1  value: 1765.325 mgc: 0.1160424 ustep: 1 
#> iter: 2  value: 1765.325 mgc: 0.0002745802 ustep: 1 
#> iter: 3  mgc: 1.817322e-09 
#> iter: 1  mgc: 1.817322e-09 
#> outer mgc:  0.3870352 
#> iter: 1  value: 1765.309 mgc: 0.1755718 ustep: 1 
#> iter: 2  value: 1765.309 mgc: 0.0006911503 ustep: 1 
#> iter: 3  value: 1765.309 mgc: 1.050644e-08 ustep: 1 
#> mgc: 1.064804e-13 
#> iter: 1  value: 1765.317 mgc: 0.07355184 ustep: 1 
#> iter: 2  value: 1765.317 mgc: 0.0001140887 ustep: 1 
#> iter: 3  mgc: 2.811283e-10 
#> iter: 1  mgc: 2.811283e-10 
#> outer mgc:  0.4668187 
#> iter: 1  value: 1765.314 mgc: 0.1166348 ustep: 1 
#> iter: 2  value: 1765.314 mgc: 0.00022384 ustep: 1 
#> iter: 3  mgc: 1.302316e-09 
#> iter: 1  value: 1765.314 mgc: 0.0733237 ustep: 1 
#> iter: 2  value: 1765.314 mgc: 6.129513e-05 ustep: 1 
#> iter: 3  mgc: 1.049345e-10 
#> iter: 1  value: 1765.315 mgc: 0.06333323 ustep: 1 
#> iter: 2  value: 1765.315 mgc: 2.382458e-05 ustep: 1 
#> iter: 3  mgc: 1.588821e-11 
#> iter: 1  value: 1765.315 mgc: 0.06528417 ustep: 1 
#> iter: 2  value: 1765.315 mgc: 1.123031e-05 ustep: 1 
#> iter: 3  mgc: 3.51871e-12 
#> iter: 1  value: 1765.316 mgc: 0.04088738 ustep: 1 
#> iter: 2  value: 1765.316 mgc: 2.857533e-06 ustep: 1 
#> iter: 3  mgc: 2.53651e-13 
#> iter: 1  value: 1765.317 mgc: 0.01854782 ustep: 1 
#> iter: 2  value: 1765.317 mgc: 5.879383e-07 ustep: 1 
#> iter: 3  mgc: 1.40819e-13 
#> iter: 1  mgc: 1.40819e-13 
#> outer mgc:  0.2490362 
#> iter: 1  value: 1765.316 mgc: 0.007818602 ustep: 1 
#> iter: 2  value: 1765.316 mgc: 8.607404e-08 ustep: 1 
#> iter: 3  mgc: 1.076017e-13 
#> iter: 1  value: 1765.316 mgc: 0.009463376 ustep: 1 
#> iter: 2  value: 1765.316 mgc: 1.261003e-07 ustep: 1 
#> iter: 3  mgc: 1.125653e-13 
#> iter: 1  mgc: 1.125653e-13 
#> outer mgc:  0.2905964 
#> iter: 1  value: 1765.315 mgc: 0.01777026 ustep: 1 
#> iter: 2  value: 1765.315 mgc: 1.242688e-06 ustep: 1 
#> iter: 3  mgc: 1.033164e-13 
#> iter: 1  mgc: 1.033164e-13 
#> outer mgc:  0.1817502 
#> iter: 1  value: 1765.315 mgc: 0.005623075 ustep: 1 
#> iter: 2  value: 1765.315 mgc: 9.594674e-08 ustep: 1 
#> iter: 3  mgc: 1.207165e-13 
#> iter: 1  mgc: 1.207165e-13 
#> outer mgc:  0.1705634 
#> iter: 1  value: 1765.315 mgc: 0.009372844 ustep: 1 
#> iter: 2  value: 1765.315 mgc: 1.427829e-06 ustep: 1 
#> iter: 3  mgc: 1.452728e-13 
#> iter: 1  value: 1765.315 mgc: 0.0061513 ustep: 1 
#> iter: 2  value: 1765.315 mgc: 8.066584e-07 ustep: 1 
#> iter: 3  mgc: 9.894233e-14 
#> iter: 1  value: 1765.315 mgc: 0.01770912 ustep: 1 
#> iter: 2  value: 1765.315 mgc: 6.691189e-06 ustep: 1 
#> iter: 3  mgc: 1.244245e-12 
#> iter: 1  mgc: 1.244245e-12 
#> outer mgc:  0.1573843 
#> iter: 1  value: 1765.31 mgc: 0.01957497 ustep: 1 
#> iter: 2  value: 1765.31 mgc: 1.113873e-05 ustep: 1 
#> iter: 3  mgc: 3.82163e-12 
#> iter: 1  mgc: 3.82163e-12 
#> outer mgc:  0.07185227 
#> iter: 1  value: 1765.312 mgc: 0.0256689 ustep: 1 
#> iter: 2  value: 1765.312 mgc: 1.254649e-05 ustep: 1 
#> iter: 3  mgc: 4.403019e-12 
#> iter: 1  mgc: 4.403019e-12 
#> outer mgc:  0.3636361 
#> iter: 1  value: 1765.308 mgc: 0.04217257 ustep: 1 
#> iter: 2  value: 1765.308 mgc: 2.208678e-05 ustep: 1 
#> iter: 3  mgc: 1.517585e-11 
#> iter: 1  mgc: 1.517585e-11 
#> outer mgc:  0.2236307 
#> iter: 1  value: 1765.304 mgc: 0.02606666 ustep: 1 
#> iter: 2  value: 1765.304 mgc: 1.15797e-05 ustep: 1 
#> iter: 3  mgc: 2.997181e-12 
#> iter: 1  mgc: 2.997181e-12 
#> outer mgc:  0.09967383 
#> iter: 1  value: 1765.303 mgc: 0.02795538 ustep: 1 
#> iter: 2  value: 1765.303 mgc: 4.031116e-06 ustep: 1 
#> iter: 3  mgc: 3.006015e-13 
#> iter: 1  mgc: 3.006015e-13 
#> outer mgc:  0.3535299 
#> iter: 1  value: 1765.302 mgc: 0.02303054 ustep: 1 
#> iter: 2  value: 1765.302 mgc: 9.493817e-06 ustep: 1 
#> iter: 3  mgc: 2.503539e-12 
#> iter: 1  mgc: 2.503539e-12 
#> outer mgc:  0.1488761 
#> iter: 1  value: 1765.297 mgc: 0.0222649 ustep: 1 
#> iter: 2  value: 1765.297 mgc: 1.025659e-05 ustep: 1 
#> iter: 3  mgc: 2.832394e-12 
#> iter: 1  mgc: 2.832394e-12 
#> outer mgc:  0.08519223 
#> iter: 1  value: 1765.3 mgc: 0.03203506 ustep: 1 
#> iter: 2  value: 1765.3 mgc: 5.438893e-06 ustep: 1 
#> iter: 3  mgc: 7.120351e-13 
#> iter: 1  mgc: 7.120351e-13 
#> outer mgc:  0.3172291 
#> iter: 1  value: 1765.298 mgc: 0.04496141 ustep: 1 
#> iter: 2  value: 1765.298 mgc: 1.185747e-05 ustep: 1 
#> iter: 3  mgc: 6.060474e-12 
#> iter: 1  mgc: 6.060474e-12 
#> outer mgc:  0.2648864 
#> iter: 1  value: 1765.303 mgc: 0.01447752 ustep: 1 
#> iter: 2  value: 1765.303 mgc: 2.632428e-06 ustep: 1 
#> iter: 3  mgc: 3.055419e-13 
#> iter: 1  mgc: 3.055419e-13 
#> outer mgc:  0.04816669 
#> iter: 1  value: 1765.299 mgc: 0.02311489 ustep: 1 
#> iter: 2  value: 1765.299 mgc: 2.966615e-06 ustep: 1 
#> iter: 3  mgc: 2.095623e-13 
#> iter: 1  mgc: 2.095623e-13 
#> outer mgc:  0.05701545 
#> iter: 1  value: 1765.301 mgc: 0.01743997 ustep: 1 
#> iter: 2  value: 1765.301 mgc: 3.591471e-06 ustep: 1 
#> iter: 3  mgc: 4.887807e-13 
#> iter: 1  mgc: 4.887807e-13 
#> outer mgc:  0.1632678 
#> iter: 1  value: 1765.298 mgc: 0.0207219 ustep: 1 
#> iter: 2  value: 1765.298 mgc: 1.339074e-06 ustep: 1 
#> iter: 3  mgc: 1.062347e-13 
#> iter: 1  mgc: 1.062347e-13 
#> outer mgc:  0.1856145 
#> iter: 1  value: 1765.296 mgc: 0.02050006 ustep: 1 
#> iter: 2  value: 1765.296 mgc: 2.436957e-06 ustep: 1 
#> iter: 3  mgc: 1.938538e-13 
#> iter: 1  mgc: 1.938538e-13 
#> outer mgc:  0.06427167 
#> iter: 1  value: 1765.296 mgc: 0.02694347 ustep: 1 
#> iter: 2  value: 1765.296 mgc: 8.722626e-07 ustep: 1 
#> iter: 3  mgc: 1.400138e-13 
#> iter: 1  mgc: 1.400138e-13 
#> outer mgc:  0.03725954 
#> iter: 1  value: 1765.296 mgc: 0.01572778 ustep: 1 
#> iter: 2  value: 1765.296 mgc: 6.921888e-07 ustep: 1 
#> iter: 3  mgc: 1.159506e-13 
#> iter: 1  mgc: 1.159506e-13 
#> outer mgc:  0.1013685 
#> iter: 1  value: 1765.295 mgc: 0.01064041 ustep: 1 
#> iter: 2  value: 1765.295 mgc: 6.901539e-07 ustep: 1 
#> iter: 3  mgc: 1.405118e-13 
#> iter: 1  mgc: 1.405118e-13 
#> outer mgc:  0.08450795 
#> iter: 1  value: 1765.295 mgc: 0.01156514 ustep: 1 
#> iter: 2  value: 1765.295 mgc: 3.553966e-07 ustep: 1 
#> iter: 3  mgc: 1.480498e-13 
#> iter: 1  value: 1765.295 mgc: 0.007713633 ustep: 1 
#> iter: 2  value: 1765.295 mgc: 4.935134e-07 ustep: 1 
#> iter: 3  mgc: 1.186354e-13 
#> iter: 1  mgc: 1.186354e-13 
#> outer mgc:  0.07288165 
#> iter: 1  value: 1765.294 mgc: 0.0171952 ustep: 1 
#> iter: 2  value: 1765.294 mgc: 2.361183e-06 ustep: 1 
#> iter: 3  mgc: 1.418645e-13 
#> iter: 1  value: 1765.292 mgc: 0.03625975 ustep: 1 
#> iter: 2  value: 1765.292 mgc: 1.959461e-05 ustep: 1 
#> iter: 3  mgc: 7.449925e-12 
#> iter: 1  mgc: 7.449925e-12 
#> outer mgc:  0.3469453 
#> iter: 1  value: 1765.286 mgc: 0.05267145 ustep: 1 
#> iter: 2  value: 1765.286 mgc: 4.837228e-05 ustep: 1 
#> iter: 3  mgc: 5.387163e-11 
#> iter: 1  mgc: 5.387163e-11 
#> outer mgc:  0.2884319 
#> iter: 1  value: 1765.28 mgc: 0.03297092 ustep: 1 
#> iter: 2  value: 1765.28 mgc: 1.090374e-05 ustep: 1 
#> iter: 3  mgc: 3.055175e-12 
#> iter: 1  mgc: 3.055175e-12 
#> outer mgc:  0.2278035 
#> iter: 1  value: 1765.276 mgc: 0.04357029 ustep: 1 
#> iter: 2  value: 1765.276 mgc: 1.446779e-05 ustep: 1 
#> iter: 3  mgc: 4.578192e-12 
#> iter: 1  mgc: 4.578192e-12 
#> outer mgc:  0.2943947 
#> iter: 1  value: 1765.279 mgc: 0.04713349 ustep: 1 
#> iter: 2  value: 1765.279 mgc: 3.447196e-06 ustep: 1 
#> iter: 3  mgc: 1.793848e-13 
#> iter: 1  mgc: 1.793848e-13 
#> outer mgc:  0.3424616 
#> iter: 1  value: 1765.281 mgc: 0.03290462 ustep: 1 
#> iter: 2  value: 1765.281 mgc: 2.996876e-06 ustep: 1 
#> iter: 3  mgc: 3.608271e-13 
#> iter: 1  value: 1765.286 mgc: 0.02590797 ustep: 1 
#> iter: 2  value: 1765.286 mgc: 2.515685e-06 ustep: 1 
#> iter: 3  mgc: 1.135819e-13 
#> iter: 1  mgc: 1.135819e-13 
#> outer mgc:  0.0862411 
#> iter: 1  value: 1765.273 mgc: 0.08587882 ustep: 1 
#> iter: 2  value: 1765.273 mgc: 1.222744e-05 ustep: 1 
#> iter: 3  mgc: 3.497115e-12 
#> iter: 1  mgc: 3.497115e-12 
#> outer mgc:  0.1740468 
#> iter: 1  value: 1765.274 mgc: 0.05350919 ustep: 1 
#> iter: 2  value: 1765.274 mgc: 1.142335e-05 ustep: 1 
#> iter: 3  mgc: 3.279912e-12 
#> iter: 1  mgc: 3.279912e-12 
#> outer mgc:  0.3473863 
#> iter: 1  value: 1765.277 mgc: 0.05621909 ustep: 1 
#> iter: 2  value: 1765.277 mgc: 1.782216e-05 ustep: 1 
#> iter: 3  mgc: 7.277884e-12 
#> iter: 1  mgc: 7.277884e-12 
#> outer mgc:  0.3144528 
#> iter: 1  value: 1765.272 mgc: 0.06237658 ustep: 1 
#> iter: 2  value: 1765.272 mgc: 3.454397e-06 ustep: 1 
#> iter: 3  mgc: 4.696932e-13 
#> iter: 1  mgc: 4.696932e-13 
#> outer mgc:  0.462884 
#> iter: 1  value: 1765.262 mgc: 0.1430556 ustep: 1 
#> iter: 2  value: 1765.262 mgc: 2.797213e-05 ustep: 1 
#> iter: 3  mgc: 6.952118e-12 
#> iter: 1  mgc: 6.952118e-12 
#> outer mgc:  0.2675018 
#> iter: 1  value: 1765.264 mgc: 0.1204855 ustep: 1 
#> iter: 2  value: 1765.264 mgc: 6.145607e-05 ustep: 1 
#> iter: 3  mgc: 7.957871e-11 
#> iter: 1  value: 1765.263 mgc: 0.04240716 ustep: 1 
#> iter: 2  value: 1765.263 mgc: 1.322042e-05 ustep: 1 
#> iter: 3  mgc: 3.783929e-12 
#> iter: 1  mgc: 3.783929e-12 
#> outer mgc:  0.2852488 
#> iter: 1  value: 1765.268 mgc: 0.05885619 ustep: 1 
#> iter: 2  value: 1765.268 mgc: 3.803551e-05 ustep: 1 
#> iter: 3  mgc: 3.323488e-11 
#> iter: 1  mgc: 3.323488e-11 
#> outer mgc:  0.1706452 
#> iter: 1  value: 1765.272 mgc: 0.06362509 ustep: 1 
#> iter: 2  value: 1765.272 mgc: 8.885756e-05 ustep: 1 
#> iter: 3  mgc: 2.092255e-10 
#> iter: 1  mgc: 2.092255e-10 
#> outer mgc:  0.4766987 
#> iter: 1  value: 1765.266 mgc: 0.09390373 ustep: 1 
#> iter: 2  value: 1765.266 mgc: 8.743365e-05 ustep: 1 
#> iter: 3  mgc: 2.438032e-10 
#> iter: 1  value: 1765.271 mgc: 0.05285753 ustep: 1 
#> iter: 2  value: 1765.271 mgc: 1.44397e-05 ustep: 1 
#> iter: 3  mgc: 5.702957e-12 
#> iter: 1  mgc: 5.702957e-12 
#> outer mgc:  0.3890927 
#> iter: 1  value: 1765.271 mgc: 0.02581307 ustep: 1 
#> iter: 2  value: 1765.271 mgc: 2.944955e-06 ustep: 1 
#> iter: 3  mgc: 2.038652e-13 
#> iter: 1  mgc: 2.038652e-13 
#> outer mgc:  0.05681168 
#> iter: 1  value: 1765.273 mgc: 0.01710828 ustep: 1 
#> iter: 2  value: 1765.273 mgc: 4.414586e-07 ustep: 1 
#> iter: 3  mgc: 1.299704e-13 
#> iter: 1  mgc: 1.299704e-13 
#> outer mgc:  0.1728743 
#> iter: 1  value: 1765.276 mgc: 0.0168127 ustep: 1 
#> iter: 2  value: 1765.276 mgc: 6.792492e-07 ustep: 1 
#> iter: 3  mgc: 2.263801e-13 
#> iter: 1  mgc: 2.263801e-13 
#> outer mgc:  0.06255569 
#> iter: 1  value: 1765.275 mgc: 0.01773209 ustep: 1 
#> iter: 2  value: 1765.275 mgc: 2.477632e-07 ustep: 1 
#> iter: 3  mgc: 1.126213e-13 
#> iter: 1  mgc: 1.126213e-13 
#> outer mgc:  0.04208672 
#> iter: 1  value: 1765.271 mgc: 0.02746398 ustep: 1 
#> iter: 2  value: 1765.271 mgc: 1.433615e-06 ustep: 1 
#> iter: 3  mgc: 1.112212e-13 
#> iter: 1  mgc: 1.112212e-13 
#> outer mgc:  0.09994126 
#> iter: 1  value: 1765.27 mgc: 0.009028785 ustep: 1 
#> iter: 2  value: 1765.27 mgc: 1.201163e-07 ustep: 1 
#> iter: 3  mgc: 1.43181e-13 
#> iter: 1  mgc: 1.43181e-13 
#> outer mgc:  0.0828494 
#> iter: 1  value: 1765.274 mgc: 0.02435679 ustep: 1 
#> iter: 2  value: 1765.274 mgc: 1.789035e-06 ustep: 1 
#> iter: 3  mgc: 1.337103e-13 
#> iter: 1  value: 1765.272 mgc: 0.01146854 ustep: 1 
#> iter: 2  value: 1765.272 mgc: 2.863238e-07 ustep: 1 
#> iter: 3  mgc: 2.283595e-13 
#> iter: 1  mgc: 2.283595e-13 
#> outer mgc:  0.03426501 
#> iter: 1  value: 1765.278 mgc: 0.02772881 ustep: 1 
#> iter: 2  value: 1765.278 mgc: 2.178577e-06 ustep: 1 
#> iter: 3  mgc: 1.469102e-13 
#> iter: 1  value: 1765.273 mgc: 0.009317575 ustep: 1 
#> iter: 2  value: 1765.273 mgc: 2.228941e-07 ustep: 1 
#> iter: 3  mgc: 8.783145e-14 
#> iter: 1  mgc: 8.783145e-14 
#> outer mgc:  0.08741192 
#> iter: 1  value: 1765.274 mgc: 0.009971616 ustep: 1 
#> iter: 2  value: 1765.274 mgc: 6.14683e-08 ustep: 1 
#> iter: 3  mgc: 1.693484e-13 
#> iter: 1  mgc: 1.693484e-13 
#> outer mgc:  0.008390059 
#> iter: 1  value: 1765.275 mgc: 0.004953808 ustep: 1 
#> iter: 2  value: 1765.275 mgc: 5.463539e-08 ustep: 1 
#> mgc: 1.773528e-13 
#> iter: 1  mgc: 1.773528e-13 
#> outer mgc:  0.02374389 
#> iter: 1  value: 1765.274 mgc: 0.003798467 ustep: 1 
#> iter: 2  value: 1765.274 mgc: 3.730581e-08 ustep: 1 
#> mgc: 1.782558e-13 
#> iter: 1  mgc: 1.782558e-13 
#> outer mgc:  0.00579776 
#> iter: 1  value: 1765.274 mgc: 0.00488562 ustep: 1 
#> iter: 2  value: 1765.274 mgc: 3.233728e-08 ustep: 1 
#> mgc: 1.743618e-13 
#> iter: 1  mgc: 1.743618e-13 
#> outer mgc:  0.009287108 
#> iter: 1  value: 1765.275 mgc: 0.007266457 ustep: 1 
#> iter: 2  value: 1765.275 mgc: 9.296279e-08 ustep: 1 
#> iter: 3  mgc: 1.409807e-13 
#> iter: 1  value: 1765.274 mgc: 0.002457041 ustep: 1 
#> iter: 2  mgc: 8.036517e-09 
#> iter: 1  mgc: 8.036517e-09 
#> outer mgc:  0.008850122 
#> iter: 1  value: 1765.275 mgc: 0.001469307 ustep: 1 
#> iter: 2  mgc: 2.554801e-09 
#> iter: 1  mgc: 2.554801e-09 
#> outer mgc:  0.006308953 
#> iter: 1  value: 1765.275 mgc: 0.001277407 ustep: 1 
#> iter: 2  mgc: 4.770493e-09 
#> iter: 1  mgc: 4.770493e-09 
#> outer mgc:  0.01223007 
#> iter: 1  value: 1765.275 mgc: 0.001026876 ustep: 1 
#> iter: 2  mgc: 2.767143e-09 
#> iter: 1  mgc: 2.767143e-09 
#> outer mgc:  0.004139453 
#> iter: 1  value: 1765.275 mgc: 0.001923468 ustep: 1 
#> iter: 2  mgc: 3.763437e-09 
#> iter: 1  value: 1765.275 mgc: 0.0009570344 ustep: 1 
#> iter: 2  mgc: 8.989763e-10 
#> iter: 1  mgc: 8.989763e-10 
#> outer mgc:  0.001981581 
#> iter: 1  value: 1765.275 mgc: 0.0006948528 ustep: 1 
#> iter: 2  mgc: 3.045922e-09 
#> iter: 1  mgc: 3.045922e-09 
#> outer mgc:  0.001180739 
#> iter: 1  value: 1765.275 mgc: 0.0008149843 ustep: 1 
#> iter: 2  mgc: 5.813061e-10 
#> iter: 1  mgc: 3.045922e-09 
#> iter: 1  mgc: 3.045922e-09 
#> outer mgc:  0.001180739 
#> iter: 1  mgc: 3.045922e-09 
#> outer mgc:  0.001180739 
#> Retro Peel: 2
#> Using stats::nlminb for optimization with these control parameters:
#> iter.max: 1000
#> eval.max: 1000
#> iter: 1  value: 1701.289 mgc: 2.24493 ustep: 1 
#> iter: 2  value: 1701.277 mgc: 0.4243491 ustep: 1 
#> iter: 3  value: 1701.277 mgc: 0.009577269 ustep: 1 
#> iter: 4  value: 1701.277 mgc: 5.187526e-06 ustep: 1 
#> iter: 5  mgc: 1.521464e-12 
#> iter: 1  mgc: 1.521464e-12 
#> Matching hessian patterns... Done
#> outer mgc:  10.42982 
#> iter: 1  value: 11413.54 mgc: 199.2832 ustep: 0.01133084 
#> iter: 2  value: 9959.712 mgc: 286.7779 ustep: 0.02073357 
#> iter: 3  value: 7850.763 mgc: 517.0925 ustep: 0.07483153 
#> iter: 4  value: 7533.769 mgc: 170.4746 ustep: 0.2736262 
#> iter: 5  value: 6974.225 mgc: 161.1333 ustep: 0.5231406 
#> iter: 6  value: 6936.443 mgc: 32.18893 ustep: 0.7233123 
#> iter: 7  value: 6933.95 mgc: 7.704587 ustep: 0.8504926 
#> iter: 8  value: 6933.912 mgc: 0.8732686 ustep: 0.9222293 
#> iter: 9  value: 6933.912 mgc: 0.04382198 ustep: 0.9603317 
#> iter: 10  value: 6933.912 mgc: 0.0003588063 ustep: 0.9799672 
#> iter: 11  value: 6933.912 mgc: 1.686221e-06 ustep: 0.9899339 
#> iter: 12  value: 6933.912 mgc: 1.497872e-08 ustep: 0.9949547 
#> iter: 13  mgc: 6.629555e-11 
#> iter: 1  value: 1751.603 mgc: 17.96582 ustep: 1 
#> iter: 2  value: 1751.341 mgc: 0.9902706 ustep: 1 
#> iter: 3  value: 1751.34 mgc: 0.03551239 ustep: 1 
#> iter: 4  value: 1751.34 mgc: 7.522222e-05 ustep: 1 
#> iter: 5  mgc: 3.276637e-10 
#> iter: 1  value: 1701.416 mgc: 1.748164 ustep: 1 
#> iter: 2  value: 1701.416 mgc: 0.006489923 ustep: 1 
#> iter: 3  value: 1701.416 mgc: 1.244565e-06 ustep: 1 
#> iter: 4  mgc: 1.352935e-13 
#> iter: 1  value: 1701.205 mgc: 0.4990443 ustep: 1 
#> iter: 2  value: 1701.205 mgc: 0.0005181898 ustep: 1 
#> iter: 3  mgc: 7.817852e-09 
#> iter: 1  mgc: 7.817852e-09 
#> outer mgc:  10.34374 
#> iter: 1  value: 1701.168 mgc: 0.4647433 ustep: 1 
#> iter: 2  value: 1701.168 mgc: 0.0003611706 ustep: 1 
#> iter: 3  mgc: 5.682437e-09 
#> iter: 1  mgc: 5.682437e-09 
#> outer mgc:  10.10271 
#> iter: 1  value: 1701.015 mgc: 0.5242656 ustep: 1 
#> iter: 2  value: 1701.015 mgc: 0.0007197376 ustep: 1 
#> iter: 3  value: 1701.015 mgc: 1.432968e-08 ustep: 1 
#> mgc: 1.342912e-13 
#> iter: 1  mgc: 1.342912e-13 
#> outer mgc:  9.788548 
#> iter: 1  value: 1700.772 mgc: 0.6264267 ustep: 1 
#> iter: 2  value: 1700.772 mgc: 0.000647725 ustep: 1 
#> iter: 3  value: 1700.772 mgc: 1.343476e-08 ustep: 1 
#> mgc: 1.798688e-13 
#> iter: 1  mgc: 1.798688e-13 
#> outer mgc:  8.926683 
#> iter: 1  value: 1700.274 mgc: 0.6987772 ustep: 1 
#> iter: 2  value: 1700.274 mgc: 0.0008092353 ustep: 1 
#> iter: 3  value: 1700.274 mgc: 4.277535e-08 ustep: 1 
#> iter: 4  mgc: 1.246299e-13 
#> iter: 1  mgc: 1.246299e-13 
#> outer mgc:  9.061904 
#> iter: 1  value: 1699.935 mgc: 1.15555 ustep: 1 
#> iter: 2  value: 1699.935 mgc: 0.001072188 ustep: 1 
#> iter: 3  value: 1699.935 mgc: 7.702413e-08 ustep: 1 
#> iter: 4  mgc: 1.488982e-13 
#> iter: 1  mgc: 1.488982e-13 
#> outer mgc:  8.345005 
#> iter: 1  value: 1699.66 mgc: 0.9279407 ustep: 1 
#> iter: 2  value: 1699.66 mgc: 0.001011822 ustep: 1 
#> iter: 3  value: 1699.66 mgc: 8.697905e-08 ustep: 1 
#> iter: 4  mgc: 1.73456e-13 
#> iter: 1  mgc: 1.73456e-13 
#> outer mgc:  17.6443 
#> iter: 1  value: 1699.659 mgc: 1.982858 ustep: 1 
#> iter: 2  value: 1699.659 mgc: 0.008290213 ustep: 1 
#> iter: 3  value: 1699.659 mgc: 1.446452e-06 ustep: 1 
#> iter: 4  mgc: 1.29258e-13 
#> iter: 1  mgc: 1.29258e-13 
#> outer mgc:  19.34016 
#> iter: 1  value: 1699.124 mgc: 1.206889 ustep: 1 
#> iter: 2  value: 1699.124 mgc: 0.00323432 ustep: 1 
#> iter: 3  value: 1699.124 mgc: 3.206204e-07 ustep: 1 
#> iter: 4  mgc: 1.524944e-13 
#> iter: 1  mgc: 1.524944e-13 
#> outer mgc:  6.693382 
#> iter: 1  value: 1698.906 mgc: 2.215957 ustep: 1 
#> iter: 2  value: 1698.906 mgc: 0.006426323 ustep: 1 
#> iter: 3  value: 1698.906 mgc: 1.208156e-06 ustep: 1 
#> iter: 4  mgc: 2.501766e-13 
#> iter: 1  mgc: 2.501766e-13 
#> outer mgc:  17.91834 
#> iter: 1  value: 1698.294 mgc: 2.177819 ustep: 1 
#> iter: 2  value: 1698.294 mgc: 0.00623088 ustep: 1 
#> iter: 3  value: 1698.294 mgc: 2.541163e-06 ustep: 1 
#> iter: 4  mgc: 4.374311e-13 
#> iter: 1  mgc: 4.374311e-13 
#> outer mgc:  29.18718 
#> iter: 1  value: 1697.684 mgc: 2.086194 ustep: 1 
#> iter: 2  value: 1697.684 mgc: 0.007326227 ustep: 1 
#> iter: 3  value: 1697.684 mgc: 1.728837e-06 ustep: 1 
#> iter: 4  mgc: 1.895862e-13 
#> iter: 1  mgc: 1.895862e-13 
#> outer mgc:  16.48461 
#> iter: 1  value: 1697.408 mgc: 1.088284 ustep: 1 
#> iter: 2  value: 1697.408 mgc: 0.001497805 ustep: 1 
#> iter: 3  value: 1697.408 mgc: 1.032283e-07 ustep: 1 
#> iter: 4  mgc: 1.263294e-13 
#> iter: 1  mgc: 1.263294e-13 
#> outer mgc:  3.393495 
#> iter: 1  value: 1697.238 mgc: 1.154505 ustep: 1 
#> iter: 2  value: 1697.238 mgc: 0.002632011 ustep: 1 
#> iter: 3  value: 1697.238 mgc: 2.739166e-07 ustep: 1 
#> iter: 4  mgc: 1.460244e-13 
#> iter: 1  mgc: 1.460244e-13 
#> outer mgc:  7.86768 
#> iter: 1  value: 1697.1 mgc: 1.263226 ustep: 1 
#> iter: 2  value: 1697.1 mgc: 0.002001521 ustep: 1 
#> iter: 3  value: 1697.1 mgc: 2.793175e-07 ustep: 1 
#> iter: 4  mgc: 1.643259e-13 
#> iter: 1  mgc: 1.643259e-13 
#> outer mgc:  6.714232 
#> iter: 1  value: 1697.603 mgc: 2.86577 ustep: 1 
#> iter: 2  value: 1697.603 mgc: 0.01347118 ustep: 1 
#> iter: 3  value: 1697.603 mgc: 5.710969e-06 ustep: 1 
#> iter: 4  mgc: 1.894805e-12 
#> iter: 1  value: 1697.124 mgc: 1.029027 ustep: 1 
#> iter: 2  value: 1697.124 mgc: 0.001731381 ustep: 1 
#> iter: 3  value: 1697.124 mgc: 9.797577e-08 ustep: 1 
#> iter: 4  mgc: 2.215128e-13 
#> iter: 1  value: 1697.077 mgc: 0.4701656 ustep: 1 
#> iter: 2  value: 1697.077 mgc: 0.0003628041 ustep: 1 
#> iter: 3  mgc: 4.370548e-09 
#> iter: 1  mgc: 4.370548e-09 
#> outer mgc:  7.161414 
#> iter: 1  value: 1697.038 mgc: 0.4269573 ustep: 1 
#> iter: 2  value: 1697.038 mgc: 0.0002059502 ustep: 1 
#> iter: 3  mgc: 2.278644e-09 
#> iter: 1  mgc: 2.278644e-09 
#> outer mgc:  5.690939 
#> iter: 1  value: 1696.993 mgc: 0.4002689 ustep: 1 
#> iter: 2  value: 1696.993 mgc: 0.0002946723 ustep: 1 
#> iter: 3  mgc: 2.77741e-09 
#> iter: 1  mgc: 2.77741e-09 
#> outer mgc:  4.2529 
#> iter: 1  value: 1696.9 mgc: 0.3292621 ustep: 1 
#> iter: 2  value: 1696.9 mgc: 0.0001037443 ustep: 1 
#> iter: 3  mgc: 7.212658e-10 
#> iter: 1  mgc: 7.212658e-10 
#> outer mgc:  3.565075 
#> iter: 1  value: 1696.822 mgc: 0.2450684 ustep: 1 
#> iter: 2  value: 1696.822 mgc: 0.0002151001 ustep: 1 
#> iter: 3  mgc: 7.226579e-10 
#> iter: 1  mgc: 7.226579e-10 
#> outer mgc:  5.191759 
#> iter: 1  value: 1696.774 mgc: 0.2898761 ustep: 1 
#> iter: 2  value: 1696.774 mgc: 0.0001096982 ustep: 1 
#> iter: 3  mgc: 6.255191e-10 
#> iter: 1  value: 1696.744 mgc: 0.1932944 ustep: 1 
#> iter: 2  value: 1696.744 mgc: 0.0001978824 ustep: 1 
#> iter: 3  mgc: 8.042676e-10 
#> iter: 1  value: 1696.747 mgc: 0.3442022 ustep: 1 
#> iter: 2  value: 1696.747 mgc: 0.0006206656 ustep: 1 
#> iter: 3  mgc: 7.913165e-09 
#> iter: 1  mgc: 7.913165e-09 
#> outer mgc:  3.382649 
#> iter: 1  value: 1697.086 mgc: 2.565992 ustep: 1 
#> iter: 2  value: 1697.086 mgc: 0.01139972 ustep: 1 
#> iter: 3  value: 1697.086 mgc: 3.591457e-06 ustep: 1 
#> iter: 4  mgc: 4.446014e-13 
#> iter: 1  value: 1696.929 mgc: 1.921621 ustep: 1 
#> iter: 2  value: 1696.929 mgc: 0.005277704 ustep: 1 
#> iter: 3  value: 1696.929 mgc: 7.663052e-07 ustep: 1 
#> iter: 4  mgc: 1.198816e-13 
#> iter: 1  value: 1696.723 mgc: 0.9825782 ustep: 1 
#> iter: 2  value: 1696.723 mgc: 0.001286937 ustep: 1 
#> iter: 3  value: 1696.723 mgc: 5.620622e-08 ustep: 1 
#> iter: 4  mgc: 1.14038e-13 
#> iter: 1  value: 1696.697 mgc: 0.3959742 ustep: 1 
#> iter: 2  value: 1696.697 mgc: 0.0002086006 ustep: 1 
#> iter: 3  mgc: 1.486493e-09 
#> iter: 1  mgc: 1.486493e-09 
#> outer mgc:  6.887789 
#> iter: 1  value: 1696.688 mgc: 0.3229579 ustep: 1 
#> iter: 2  value: 1696.688 mgc: 0.0001891618 ustep: 1 
#> iter: 3  mgc: 1.344068e-09 
#> iter: 1  mgc: 1.344068e-09 
#> outer mgc:  3.225123 
#> iter: 1  value: 1696.665 mgc: 0.2654109 ustep: 1 
#> iter: 2  value: 1696.665 mgc: 0.0001455687 ustep: 1 
#> iter: 3  mgc: 5.841813e-10 
#> iter: 1  mgc: 5.841813e-10 
#> outer mgc:  4.029473 
#> iter: 1  value: 1696.659 mgc: 0.2490691 ustep: 1 
#> iter: 2  value: 1696.659 mgc: 0.0001431153 ustep: 1 
#> iter: 3  mgc: 7.10517e-10 
#> iter: 1  mgc: 7.10517e-10 
#> outer mgc:  4.394951 
#> iter: 1  value: 1696.597 mgc: 0.5064491 ustep: 1 
#> iter: 2  value: 1696.597 mgc: 0.0003601707 ustep: 1 
#> iter: 3  mgc: 6.1488e-09 
#> iter: 1  mgc: 6.1488e-09 
#> outer mgc:  4.794136 
#> iter: 1  value: 1696.525 mgc: 0.4937333 ustep: 1 
#> iter: 2  value: 1696.525 mgc: 0.0006727195 ustep: 1 
#> iter: 3  mgc: 7.600346e-09 
#> iter: 1  mgc: 7.600346e-09 
#> outer mgc:  3.441706 
#> iter: 1  value: 1696.085 mgc: 0.9678049 ustep: 1 
#> iter: 2  value: 1696.085 mgc: 0.00291554 ustep: 1 
#> iter: 3  value: 1696.085 mgc: 1.172282e-07 ustep: 1 
#> iter: 4  mgc: 1.715032e-13 
#> iter: 1  value: 1696.282 mgc: 0.4741865 ustep: 1 
#> iter: 2  value: 1696.282 mgc: 0.0002789752 ustep: 1 
#> iter: 3  mgc: 5.810617e-09 
#> iter: 1  mgc: 5.810617e-09 
#> outer mgc:  4.457041 
#> iter: 1  value: 1696.346 mgc: 0.3832927 ustep: 1 
#> iter: 2  value: 1696.346 mgc: 0.0009759744 ustep: 1 
#> iter: 3  value: 1696.346 mgc: 1.950276e-08 ustep: 1 
#> mgc: 1.890324e-13 
#> iter: 1  value: 1696.344 mgc: 0.5290949 ustep: 1 
#> iter: 2  value: 1696.344 mgc: 0.000772743 ustep: 1 
#> iter: 3  value: 1696.344 mgc: 1.737516e-08 ustep: 1 
#> mgc: 1.171872e-13 
#> iter: 1  value: 1696.349 mgc: 0.6387235 ustep: 1 
#> iter: 2  value: 1696.349 mgc: 0.0007550151 ustep: 1 
#> iter: 3  value: 1696.349 mgc: 1.958828e-08 ustep: 1 
#> mgc: 1.060114e-13 
#> iter: 1  value: 1696.316 mgc: 0.5312538 ustep: 1 
#> iter: 2  value: 1696.316 mgc: 0.0004569483 ustep: 1 
#> iter: 3  mgc: 8.18035e-09 
#> iter: 1  value: 1696.282 mgc: 0.2135618 ustep: 1 
#> iter: 2  value: 1696.282 mgc: 7.416756e-05 ustep: 1 
#> iter: 3  mgc: 2.165784e-10 
#> iter: 1  mgc: 2.165784e-10 
#> outer mgc:  1.811551 
#> iter: 1  value: 1696.274 mgc: 0.09675683 ustep: 1 
#> iter: 2  value: 1696.274 mgc: 2.132903e-05 ustep: 1 
#> iter: 3  mgc: 1.076807e-11 
#> iter: 1  value: 1696.27 mgc: 0.1092528 ustep: 1 
#> iter: 2  value: 1696.27 mgc: 2.718564e-05 ustep: 1 
#> iter: 3  mgc: 1.740127e-11 
#> iter: 1  mgc: 1.740127e-11 
#> outer mgc:  3.244007 
#> iter: 1  value: 1696.282 mgc: 0.1483032 ustep: 1 
#> iter: 2  value: 1696.282 mgc: 5.010571e-05 ustep: 1 
#> iter: 3  mgc: 9.032144e-11 
#> iter: 1  mgc: 9.032144e-11 
#> outer mgc:  1.739481 
#> iter: 1  value: 1696.289 mgc: 0.2149008 ustep: 1 
#> iter: 2  value: 1696.289 mgc: 7.506885e-05 ustep: 1 
#> iter: 3  mgc: 1.227377e-10 
#> iter: 1  value: 1696.317 mgc: 0.2605345 ustep: 1 
#> iter: 2  value: 1696.317 mgc: 0.0001513058 ustep: 1 
#> iter: 3  mgc: 6.970737e-11 
#> iter: 1  value: 1696.41 mgc: 0.6761913 ustep: 1 
#> iter: 2  value: 1696.41 mgc: 0.001024445 ustep: 1 
#> iter: 3  mgc: 3.219298e-09 
#> iter: 1  mgc: 3.219298e-09 
#> outer mgc:  2.689213 
#> iter: 1  value: 1696.111 mgc: 0.4444629 ustep: 1 
#> iter: 2  value: 1696.111 mgc: 0.0006640875 ustep: 1 
#> iter: 3  mgc: 7.343489e-09 
#> iter: 1  mgc: 7.343489e-09 
#> outer mgc:  6.799823 
#> iter: 1  value: 1696.629 mgc: 2.314341 ustep: 1 
#> iter: 2  value: 1696.629 mgc: 0.008823474 ustep: 1 
#> iter: 3  value: 1696.629 mgc: 8.812472e-07 ustep: 1 
#> iter: 4  mgc: 1.256301e-13 
#> iter: 1  value: 1696.379 mgc: 1.65668 ustep: 1 
#> iter: 2  value: 1696.379 mgc: 0.003720527 ustep: 1 
#> iter: 3  value: 1696.379 mgc: 4.69283e-07 ustep: 1 
#> iter: 4  mgc: 1.297101e-13 
#> iter: 1  value: 1696.107 mgc: 0.3384739 ustep: 1 
#> iter: 2  value: 1696.107 mgc: 0.0001566837 ustep: 1 
#> iter: 3  mgc: 8.631024e-10 
#> iter: 1  mgc: 8.631024e-10 
#> outer mgc:  6.748154 
#> iter: 1  value: 1696.096 mgc: 0.2582655 ustep: 1 
#> iter: 2  value: 1696.096 mgc: 0.0001193228 ustep: 1 
#> iter: 3  mgc: 6.342328e-10 
#> iter: 1  mgc: 6.342328e-10 
#> outer mgc:  3.701549 
#> iter: 1  value: 1696.107 mgc: 0.1611495 ustep: 1 
#> iter: 2  value: 1696.107 mgc: 4.059142e-05 ustep: 1 
#> iter: 3  mgc: 3.509796e-11 
#> iter: 1  mgc: 3.509796e-11 
#> outer mgc:  4.681224 
#> iter: 1  value: 1696.108 mgc: 0.0765647 ustep: 1 
#> iter: 2  value: 1696.108 mgc: 2.137721e-05 ustep: 1 
#> iter: 3  mgc: 1.979387e-11 
#> iter: 1  mgc: 1.979387e-11 
#> outer mgc:  2.90821 
#> iter: 1  value: 1696.08 mgc: 0.2430232 ustep: 1 
#> iter: 2  value: 1696.08 mgc: 0.0001148101 ustep: 1 
#> iter: 3  mgc: 1.531365e-10 
#> iter: 1  value: 1696.05 mgc: 0.1560305 ustep: 1 
#> iter: 2  value: 1696.05 mgc: 4.104691e-05 ustep: 1 
#> iter: 3  mgc: 3.069488e-11 
#> iter: 1  mgc: 3.069488e-11 
#> outer mgc:  4.403251 
#> iter: 1  value: 1696.337 mgc: 0.4301988 ustep: 1 
#> iter: 2  value: 1696.337 mgc: 0.0004769217 ustep: 1 
#> iter: 3  mgc: 9.413765e-09 
#> iter: 1  mgc: 9.413765e-09 
#> outer mgc:  2.788843 
#> iter: 1  value: 1695.826 mgc: 0.6693616 ustep: 1 
#> iter: 2  value: 1695.826 mgc: 0.0011576 ustep: 1 
#> iter: 3  value: 1695.826 mgc: 2.333234e-08 ustep: 1 
#> mgc: 1.868188e-13 
#> iter: 1  mgc: 1.868188e-13 
#> outer mgc:  6.696837 
#> iter: 1  value: 1695.764 mgc: 0.7954295 ustep: 1 
#> iter: 2  value: 1695.764 mgc: 0.00107599 ustep: 1 
#> iter: 3  value: 1695.764 mgc: 4.950378e-08 ustep: 1 
#> iter: 4  mgc: 1.317616e-13 
#> iter: 1  value: 1695.804 mgc: 0.7928675 ustep: 1 
#> iter: 2  value: 1695.804 mgc: 0.0009849975 ustep: 1 
#> iter: 3  value: 1695.804 mgc: 4.338913e-08 ustep: 1 
#> iter: 4  mgc: 1.54986e-13 
#> iter: 1  value: 1695.802 mgc: 0.5750542 ustep: 1 
#> iter: 2  value: 1695.802 mgc: 0.0005296485 ustep: 1 
#> iter: 3  value: 1695.802 mgc: 1.108314e-08 ustep: 1 
#> mgc: 1.551775e-13 
#> iter: 1  mgc: 1.551775e-13 
#> outer mgc:  6.348813 
#> iter: 1  value: 1695.828 mgc: 0.2695554 ustep: 1 
#> iter: 2  value: 1695.828 mgc: 0.0001493381 ustep: 1 
#> iter: 3  mgc: 3.772885e-10 
#> iter: 1  mgc: 3.772885e-10 
#> outer mgc:  1.15881 
#> iter: 1  value: 1695.896 mgc: 0.3349468 ustep: 1 
#> iter: 2  value: 1695.896 mgc: 0.0002249345 ustep: 1 
#> iter: 3  mgc: 2.078811e-09 
#> iter: 1  mgc: 2.078811e-09 
#> outer mgc:  3.15298 
#> iter: 1  value: 1696.01 mgc: 0.3615988 ustep: 1 
#> iter: 2  value: 1696.01 mgc: 0.0003615575 ustep: 1 
#> iter: 3  mgc: 1.019563e-09 
#> iter: 1  mgc: 1.019563e-09 
#> outer mgc:  1.979917 
#> iter: 1  value: 1695.997 mgc: 0.1425814 ustep: 1 
#> iter: 2  value: 1695.997 mgc: 2.52912e-05 ustep: 1 
#> iter: 3  mgc: 2.61184e-11 
#> iter: 1  mgc: 2.61184e-11 
#> outer mgc:  1.412818 
#> iter: 1  value: 1695.887 mgc: 0.08579106 ustep: 1 
#> iter: 2  value: 1695.887 mgc: 3.289438e-05 ustep: 1 
#> iter: 3  mgc: 1.898305e-11 
#> iter: 1  value: 1695.759 mgc: 0.0883479 ustep: 1 
#> iter: 2  value: 1695.759 mgc: 1.993522e-05 ustep: 1 
#> iter: 3  mgc: 1.040526e-11 
#> iter: 1  mgc: 1.040526e-11 
#> outer mgc:  1.326337 
#> iter: 1  value: 1695.844 mgc: 0.3734058 ustep: 1 
#> iter: 2  value: 1695.844 mgc: 0.0001675438 ustep: 1 
#> iter: 3  mgc: 2.369664e-10 
#> iter: 1  mgc: 2.369664e-10 
#> outer mgc:  0.9565213 
#> iter: 1  value: 1695.715 mgc: 1.102704 ustep: 1 
#> iter: 2  value: 1695.715 mgc: 0.001255269 ustep: 1 
#> iter: 3  value: 1695.715 mgc: 2.935187e-08 ustep: 1 
#> iter: 4  mgc: 1.232279e-13 
#> iter: 1  value: 1695.807 mgc: 0.5819494 ustep: 1 
#> iter: 2  value: 1695.807 mgc: 0.0004716344 ustep: 1 
#> iter: 3  mgc: 4.123347e-09 
#> iter: 1  value: 1695.847 mgc: 0.4058934 ustep: 1 
#> iter: 2  value: 1695.847 mgc: 0.000291903 ustep: 1 
#> iter: 3  mgc: 1.527188e-09 
#> iter: 1  value: 1695.875 mgc: 0.2935855 ustep: 1 
#> iter: 2  value: 1695.875 mgc: 0.0001988445 ustep: 1 
#> iter: 3  mgc: 6.693649e-10 
#> iter: 1  value: 1695.854 mgc: 0.1045127 ustep: 1 
#> iter: 2  value: 1695.854 mgc: 2.168513e-05 ustep: 1 
#> iter: 3  mgc: 1.034532e-11 
#> iter: 1  value: 1695.848 mgc: 0.04590152 ustep: 1 
#> iter: 2  value: 1695.848 mgc: 4.184792e-06 ustep: 1 
#> iter: 3  mgc: 4.189068e-13 
#> iter: 1  mgc: 4.189068e-13 
#> outer mgc:  0.946744 
#> iter: 1  value: 1695.849 mgc: 0.03765838 ustep: 1 
#> iter: 2  value: 1695.849 mgc: 2.213675e-06 ustep: 1 
#> iter: 3  mgc: 2.408439e-13 
#> iter: 1  mgc: 2.408439e-13 
#> outer mgc:  0.9444858 
#> iter: 1  value: 1695.853 mgc: 0.04710124 ustep: 1 
#> iter: 2  value: 1695.853 mgc: 5.173485e-06 ustep: 1 
#> iter: 3  mgc: 4.460852e-13 
#> iter: 1  mgc: 4.460852e-13 
#> outer mgc:  0.9335024 
#> iter: 1  value: 1695.848 mgc: 0.05075827 ustep: 1 
#> iter: 2  value: 1695.848 mgc: 4.118374e-06 ustep: 1 
#> iter: 3  mgc: 1.064452e-12 
#> iter: 1  mgc: 1.064452e-12 
#> outer mgc:  0.9379458 
#> iter: 1  value: 1695.848 mgc: 0.07275599 ustep: 1 
#> iter: 2  value: 1695.848 mgc: 1.023387e-05 ustep: 1 
#> iter: 3  mgc: 1.505351e-12 
#> iter: 1  value: 1695.841 mgc: 0.04545835 ustep: 1 
#> iter: 2  value: 1695.841 mgc: 2.956342e-06 ustep: 1 
#> iter: 3  mgc: 1.874779e-13 
#> iter: 1  value: 1695.824 mgc: 0.1229001 ustep: 1 
#> iter: 2  value: 1695.824 mgc: 2.157751e-05 ustep: 1 
#> iter: 3  mgc: 9.984197e-12 
#> iter: 1  mgc: 9.984197e-12 
#> outer mgc:  1.047556 
#> iter: 1  value: 1695.932 mgc: 0.2327646 ustep: 1 
#> iter: 2  value: 1695.932 mgc: 5.46594e-05 ustep: 1 
#> iter: 3  mgc: 1.366687e-10 
#> iter: 1  mgc: 1.366687e-10 
#> outer mgc:  0.9834283 
#> iter: 1  value: 1695.955 mgc: 0.1395329 ustep: 1 
#> iter: 2  value: 1695.955 mgc: 6.66365e-05 ustep: 1 
#> iter: 3  mgc: 9.267578e-11 
#> iter: 1  mgc: 9.267578e-11 
#> outer mgc:  2.331877 
#> iter: 1  value: 1695.99 mgc: 0.2384952 ustep: 1 
#> iter: 2  value: 1695.99 mgc: 7.125333e-05 ustep: 1 
#> iter: 3  mgc: 5.827624e-11 
#> iter: 1  mgc: 5.827624e-11 
#> outer mgc:  0.8449921 
#> iter: 1  value: 1695.989 mgc: 0.1026626 ustep: 1 
#> iter: 2  value: 1695.989 mgc: 4.103115e-05 ustep: 1 
#> iter: 3  mgc: 6.427323e-11 
#> iter: 1  mgc: 6.427323e-11 
#> outer mgc:  0.6824409 
#> iter: 1  value: 1695.927 mgc: 0.3492303 ustep: 1 
#> iter: 2  value: 1695.927 mgc: 0.0002732261 ustep: 1 
#> iter: 3  mgc: 2.088769e-09 
#> iter: 1  mgc: 2.088769e-09 
#> outer mgc:  4.243316 
#> iter: 1  value: 1695.931 mgc: 0.220395 ustep: 1 
#> iter: 2  value: 1695.931 mgc: 0.0001077068 ustep: 1 
#> iter: 3  mgc: 1.576021e-10 
#> iter: 1  mgc: 1.576021e-10 
#> outer mgc:  2.570943 
#> iter: 1  value: 1695.9 mgc: 0.08711177 ustep: 1 
#> iter: 2  value: 1695.9 mgc: 1.569297e-05 ustep: 1 
#> iter: 3  mgc: 4.186564e-12 
#> iter: 1  mgc: 4.186564e-12 
#> outer mgc:  0.6174168 
#> iter: 1  value: 1695.896 mgc: 0.09022735 ustep: 1 
#> iter: 2  value: 1695.896 mgc: 1.519956e-05 ustep: 1 
#> iter: 3  mgc: 9.556675e-12 
#> iter: 1  value: 1695.891 mgc: 0.07956603 ustep: 1 
#> iter: 2  value: 1695.891 mgc: 1.081435e-05 ustep: 1 
#> iter: 3  mgc: 4.665188e-12 
#> iter: 1  mgc: 4.665188e-12 
#> outer mgc:  1.310563 
#> iter: 1  value: 1695.878 mgc: 0.1443992 ustep: 1 
#> iter: 2  value: 1695.878 mgc: 6.490784e-05 ustep: 1 
#> iter: 3  mgc: 8.787726e-11 
#> iter: 1  mgc: 8.787726e-11 
#> outer mgc:  0.9678046 
#> iter: 1  value: 1695.86 mgc: 0.1875428 ustep: 1 
#> iter: 2  value: 1695.86 mgc: 4.015437e-05 ustep: 1 
#> iter: 3  mgc: 3.689833e-11 
#> iter: 1  mgc: 3.689833e-11 
#> outer mgc:  1.492537 
#> iter: 1  value: 1695.809 mgc: 0.1903263 ustep: 1 
#> iter: 2  value: 1695.809 mgc: 8.769977e-05 ustep: 1 
#> iter: 3  mgc: 1.226694e-10 
#> iter: 1  mgc: 1.226694e-10 
#> outer mgc:  0.5284049 
#> iter: 1  value: 1695.821 mgc: 0.1244973 ustep: 1 
#> iter: 2  value: 1695.821 mgc: 4.807329e-05 ustep: 1 
#> iter: 3  mgc: 1.359076e-10 
#> iter: 1  mgc: 1.359076e-10 
#> outer mgc:  0.8363515 
#> iter: 1  value: 1695.907 mgc: 0.4120925 ustep: 1 
#> iter: 2  value: 1695.907 mgc: 0.0003176143 ustep: 1 
#> iter: 3  mgc: 2.13291e-09 
#> iter: 1  mgc: 2.13291e-09 
#> outer mgc:  3.40151 
#> iter: 1  value: 1695.993 mgc: 0.3262709 ustep: 1 
#> iter: 2  value: 1695.993 mgc: 0.000189821 ustep: 1 
#> iter: 3  mgc: 4.669816e-10 
#> iter: 1  mgc: 4.669816e-10 
#> outer mgc:  5.217506 
#> iter: 1  value: 1696.01 mgc: 0.2416832 ustep: 1 
#> iter: 2  value: 1696.01 mgc: 0.0001226025 ustep: 1 
#> iter: 3  mgc: 7.116882e-10 
#> iter: 1  mgc: 7.116882e-10 
#> outer mgc:  0.9031873 
#> iter: 1  value: 1695.938 mgc: 0.1819976 ustep: 1 
#> iter: 2  value: 1695.938 mgc: 7.741535e-05 ustep: 1 
#> iter: 3  mgc: 2.59925e-10 
#> iter: 1  mgc: 2.59925e-10 
#> outer mgc:  0.9202696 
#> iter: 1  value: 1695.916 mgc: 0.1863377 ustep: 1 
#> iter: 2  value: 1695.916 mgc: 8.05955e-05 ustep: 1 
#> iter: 3  mgc: 1.440547e-10 
#> iter: 1  mgc: 1.440547e-10 
#> outer mgc:  0.8769405 
#> iter: 1  value: 1695.842 mgc: 0.1629374 ustep: 1 
#> iter: 2  value: 1695.842 mgc: 4.647424e-05 ustep: 1 
#> iter: 3  mgc: 8.063218e-11 
#> iter: 1  mgc: 8.063218e-11 
#> outer mgc:  0.9102165 
#> iter: 1  value: 1695.851 mgc: 0.1483937 ustep: 1 
#> iter: 2  value: 1695.851 mgc: 7.937678e-05 ustep: 1 
#> iter: 3  mgc: 1.012012e-10 
#> iter: 1  mgc: 1.012012e-10 
#> outer mgc:  2.010992 
#> iter: 1  value: 1695.885 mgc: 0.1808217 ustep: 1 
#> iter: 2  value: 1695.885 mgc: 3.519845e-05 ustep: 1 
#> iter: 3  mgc: 7.084103e-11 
#> iter: 1  mgc: 7.084103e-11 
#> outer mgc:  1.637215 
#> iter: 1  value: 1695.901 mgc: 0.1713348 ustep: 1 
#> iter: 2  value: 1695.901 mgc: 8.359214e-05 ustep: 1 
#> iter: 3  mgc: 1.647175e-10 
#> iter: 1  mgc: 1.647175e-10 
#> outer mgc:  0.8961098 
#> iter: 1  value: 1695.893 mgc: 0.183847 ustep: 1 
#> iter: 2  value: 1695.893 mgc: 4.285095e-05 ustep: 1 
#> iter: 3  mgc: 7.297855e-11 
#> iter: 1  mgc: 7.297855e-11 
#> outer mgc:  1.097992 
#> iter: 1  value: 1695.869 mgc: 0.1031574 ustep: 1 
#> iter: 2  value: 1695.869 mgc: 4.041863e-05 ustep: 1 
#> iter: 3  mgc: 2.774214e-11 
#> iter: 1  mgc: 2.774214e-11 
#> outer mgc:  0.7980964 
#> iter: 1  value: 1695.861 mgc: 0.09906088 ustep: 1 
#> iter: 2  value: 1695.861 mgc: 2.385571e-05 ustep: 1 
#> iter: 3  mgc: 1.668712e-11 
#> iter: 1  value: 1695.863 mgc: 0.079692 ustep: 1 
#> iter: 2  value: 1695.863 mgc: 5.769127e-06 ustep: 1 
#> iter: 3  mgc: 1.137893e-12 
#> iter: 1  mgc: 1.137893e-12 
#> outer mgc:  1.144463 
#> iter: 1  value: 1695.877 mgc: 0.08938931 ustep: 1 
#> iter: 2  value: 1695.877 mgc: 1.24444e-05 ustep: 1 
#> iter: 3  mgc: 2.619683e-12 
#> iter: 1  mgc: 2.619683e-12 
#> outer mgc:  0.271055 
#> iter: 1  value: 1695.874 mgc: 0.06107327 ustep: 1 
#> iter: 2  value: 1695.874 mgc: 7.301751e-06 ustep: 1 
#> iter: 3  mgc: 1.865275e-12 
#> iter: 1  mgc: 1.865275e-12 
#> outer mgc:  0.1664725 
#> iter: 1  value: 1695.877 mgc: 0.09203479 ustep: 1 
#> iter: 2  value: 1695.877 mgc: 1.06793e-05 ustep: 1 
#> iter: 3  mgc: 8.122066e-13 
#> iter: 1  mgc: 8.122066e-13 
#> outer mgc:  0.3178567 
#> iter: 1  value: 1695.877 mgc: 0.03695496 ustep: 1 
#> iter: 2  value: 1695.877 mgc: 9.394619e-06 ustep: 1 
#> iter: 3  mgc: 5.94353e-13 
#> iter: 1  mgc: 5.94353e-13 
#> outer mgc:  0.3995466 
#> iter: 1  value: 1695.873 mgc: 0.09080528 ustep: 1 
#> iter: 2  value: 1695.873 mgc: 3.116462e-05 ustep: 1 
#> iter: 3  mgc: 1.091019e-11 
#> iter: 1  mgc: 1.091019e-11 
#> outer mgc:  0.553819 
#> iter: 1  value: 1695.865 mgc: 0.07137024 ustep: 1 
#> iter: 2  value: 1695.865 mgc: 7.386196e-06 ustep: 1 
#> iter: 3  mgc: 2.240467e-12 
#> iter: 1  mgc: 2.240467e-12 
#> outer mgc:  0.7064942 
#> iter: 1  value: 1695.863 mgc: 0.07188282 ustep: 1 
#> iter: 2  value: 1695.863 mgc: 6.094381e-06 ustep: 1 
#> iter: 3  mgc: 8.928504e-13 
#> iter: 1  mgc: 8.928504e-13 
#> outer mgc:  0.1227844 
#> iter: 1  value: 1695.866 mgc: 0.05491323 ustep: 1 
#> iter: 2  value: 1695.866 mgc: 4.231252e-06 ustep: 1 
#> iter: 3  mgc: 1.094523e-13 
#> iter: 1  mgc: 1.094523e-13 
#> outer mgc:  0.1175261 
#> iter: 1  value: 1695.869 mgc: 0.03277608 ustep: 1 
#> iter: 2  value: 1695.869 mgc: 3.136704e-06 ustep: 1 
#> iter: 3  mgc: 1.253814e-13 
#> iter: 1  mgc: 1.253814e-13 
#> outer mgc:  0.1378319 
#> iter: 1  value: 1695.874 mgc: 0.01928638 ustep: 1 
#> iter: 2  value: 1695.874 mgc: 9.441424e-07 ustep: 1 
#> iter: 3  mgc: 1.594448e-13 
#> iter: 1  mgc: 1.594448e-13 
#> outer mgc:  0.1239274 
#> iter: 1  value: 1695.873 mgc: 0.03087581 ustep: 1 
#> iter: 2  value: 1695.873 mgc: 1.52432e-06 ustep: 1 
#> iter: 3  mgc: 1.211903e-13 
#> iter: 1  mgc: 1.211903e-13 
#> outer mgc:  0.06010895 
#> iter: 1  value: 1695.869 mgc: 0.02484718 ustep: 1 
#> iter: 2  value: 1695.869 mgc: 1.642547e-06 ustep: 1 
#> iter: 3  mgc: 1.11279e-13 
#> iter: 1  mgc: 1.11279e-13 
#> outer mgc:  0.1063055 
#> iter: 1  value: 1695.864 mgc: 0.043083 ustep: 1 
#> iter: 2  value: 1695.864 mgc: 2.207224e-06 ustep: 1 
#> iter: 3  mgc: 1.756761e-13 
#> iter: 1  value: 1695.867 mgc: 0.02818134 ustep: 1 
#> iter: 2  value: 1695.867 mgc: 7.318877e-07 ustep: 1 
#> iter: 3  mgc: 1.688434e-13 
#> iter: 1  value: 1695.868 mgc: 0.02540987 ustep: 1 
#> iter: 2  value: 1695.868 mgc: 6.695718e-07 ustep: 1 
#> iter: 3  mgc: 1.5862e-13 
#> iter: 1  value: 1695.868 mgc: 0.02441649 ustep: 1 
#> iter: 2  value: 1695.868 mgc: 8.374311e-07 ustep: 1 
#> iter: 3  mgc: 2.411869e-13 
#> iter: 1  value: 1695.869 mgc: 0.01430421 ustep: 1 
#> iter: 2  value: 1695.869 mgc: 3.610448e-07 ustep: 1 
#> iter: 3  mgc: 1.488625e-13 
#> iter: 1  value: 1695.869 mgc: 0.006081541 ustep: 1 
#> iter: 2  value: 1695.869 mgc: 6.525807e-08 ustep: 1 
#> iter: 3  mgc: 1.380323e-13 
#> iter: 1  mgc: 1.380323e-13 
#> outer mgc:  0.08905375 
#> iter: 1  value: 1695.869 mgc: 0.005325096 ustep: 1 
#> iter: 2  value: 1695.869 mgc: 3.098394e-08 ustep: 1 
#> iter: 3  mgc: 1.492854e-13 
#> iter: 1  mgc: 1.492854e-13 
#> outer mgc:  0.09549226 
#> iter: 1  value: 1695.869 mgc: 0.004975144 ustep: 1 
#> iter: 2  value: 1695.869 mgc: 2.9506e-08 ustep: 1 
#> iter: 3  mgc: 1.705447e-13 
#> iter: 1  mgc: 1.705447e-13 
#> outer mgc:  0.05401707 
#> iter: 1  value: 1695.869 mgc: 0.002962756 ustep: 1 
#> iter: 2  value: 1695.869 mgc: 1.897197e-08 ustep: 1 
#> mgc: 1.682714e-13 
#> iter: 1  mgc: 1.682714e-13 
#> outer mgc:  0.06664724 
#> iter: 1  value: 1695.868 mgc: 0.005027396 ustep: 1 
#> iter: 2  value: 1695.868 mgc: 1.867223e-08 ustep: 1 
#> mgc: 1.164645e-13 
#> iter: 1  mgc: 1.164645e-13 
#> outer mgc:  0.04932767 
#> iter: 1  value: 1695.868 mgc: 0.003286759 ustep: 1 
#> iter: 2  value: 1695.868 mgc: 3.732002e-08 ustep: 1 
#> mgc: 1.256966e-13 
#> iter: 1  mgc: 1.256966e-13 
#> outer mgc:  0.04425 
#> iter: 1  value: 1695.867 mgc: 0.00543139 ustep: 1 
#> iter: 2  value: 1695.867 mgc: 3.530448e-08 ustep: 1 
#> mgc: 1.291498e-13 
#> iter: 1  mgc: 1.291498e-13 
#> outer mgc:  0.03096208 
#> iter: 1  value: 1695.867 mgc: 0.008236922 ustep: 1 
#> iter: 2  value: 1695.867 mgc: 2.108305e-07 ustep: 1 
#> iter: 3  mgc: 1.199765e-13 
#> iter: 1  value: 1695.867 mgc: 0.009884794 ustep: 1 
#> iter: 2  value: 1695.867 mgc: 2.360177e-07 ustep: 1 
#> iter: 3  mgc: 1.444761e-13 
#> iter: 1  value: 1695.867 mgc: 0.004167996 ustep: 1 
#> iter: 2  value: 1695.867 mgc: 3.504718e-08 ustep: 1 
#> mgc: 1.08946e-13 
#> iter: 1  value: 1695.867 mgc: 0.001776975 ustep: 1 
#> iter: 2  mgc: 6.370167e-09 
#> iter: 1  mgc: 6.370167e-09 
#> outer mgc:  0.03259532 
#> iter: 1  value: 1695.867 mgc: 0.001413225 ustep: 1 
#> iter: 2  mgc: 3.046377e-09 
#> iter: 1  mgc: 3.046377e-09 
#> outer mgc:  0.02862367 
#> iter: 1  value: 1695.867 mgc: 0.0009821306 ustep: 1 
#> iter: 2  mgc: 2.363802e-09 
#> iter: 1  mgc: 2.363802e-09 
#> outer mgc:  0.02312058 
#> iter: 1  value: 1695.867 mgc: 0.001436316 ustep: 1 
#> iter: 2  mgc: 2.560227e-09 
#> iter: 1  mgc: 2.560227e-09 
#> outer mgc:  0.02699949 
#> iter: 1  value: 1695.867 mgc: 0.0009721213 ustep: 1 
#> iter: 2  mgc: 1.991617e-09 
#> iter: 1  value: 1695.867 mgc: 0.0009437432 ustep: 1 
#> iter: 2  mgc: 2.077106e-09 
#> iter: 1  value: 1695.867 mgc: 0.002293316 ustep: 1 
#> iter: 2  value: 1695.867 mgc: 1.226548e-08 ustep: 1 
#> mgc: 1.450886e-13 
#> iter: 1  mgc: 1.450886e-13 
#> outer mgc:  0.02665073 
#> iter: 1  value: 1695.868 mgc: 0.003222897 ustep: 1 
#> iter: 2  value: 1695.868 mgc: 1.927361e-08 ustep: 1 
#> mgc: 1.099871e-13 
#> iter: 1  value: 1695.868 mgc: 0.00764209 ustep: 1 
#> iter: 2  value: 1695.868 mgc: 9.5078e-08 ustep: 1 
#> iter: 3  mgc: 1.323951e-13 
#> iter: 1  mgc: 1.323951e-13 
#> outer mgc:  0.0315523 
#> iter: 1  value: 1695.867 mgc: 0.007773133 ustep: 1 
#> iter: 2  value: 1695.867 mgc: 2.1467e-07 ustep: 1 
#> iter: 3  mgc: 9.948893e-14 
#> iter: 1  mgc: 9.948893e-14 
#> outer mgc:  0.04009792 
#> iter: 1  value: 1695.867 mgc: 0.006160984 ustep: 1 
#> iter: 2  value: 1695.867 mgc: 1.019523e-07 ustep: 1 
#> mgc: 2.021887e-13 
#> iter: 1  mgc: 2.021887e-13 
#> outer mgc:  0.01656306 
#> iter: 1  value: 1695.866 mgc: 0.004172335 ustep: 1 
#> iter: 2  value: 1695.866 mgc: 2.781009e-08 ustep: 1 
#> mgc: 1.649287e-13 
#> iter: 1  mgc: 1.649287e-13 
#> outer mgc:  0.02285762 
#> iter: 1  value: 1695.866 mgc: 0.003566006 ustep: 1 
#> iter: 2  value: 1695.866 mgc: 3.991536e-08 ustep: 1 
#> mgc: 1.286193e-13 
#> iter: 1  mgc: 1.286193e-13 
#> outer mgc:  0.02875903 
#> iter: 1  value: 1695.866 mgc: 0.002551471 ustep: 1 
#> iter: 2  value: 1695.866 mgc: 2.236137e-08 ustep: 1 
#> mgc: 1.098757e-13 
#> iter: 1  mgc: 1.098757e-13 
#> outer mgc:  0.01249567 
#> iter: 1  value: 1695.867 mgc: 0.002968237 ustep: 1 
#> iter: 2  value: 1695.867 mgc: 1.366321e-08 ustep: 1 
#> mgc: 9.265548e-14 
#> iter: 1  value: 1695.867 mgc: 0.003949987 ustep: 1 
#> iter: 2  value: 1695.867 mgc: 2.839723e-08 ustep: 1 
#> mgc: 9.970738e-14 
#> iter: 1  mgc: 9.265548e-14 
#> outer mgc:  0.01857601 
#> iter: 1  value: 1695.867 mgc: 0.00394889 ustep: 1 
#> iter: 2  value: 1695.867 mgc: 1.506997e-08 ustep: 1 
#> mgc: 1.470128e-13 
#> iter: 1  mgc: 1.470128e-13 
#> outer mgc:  0.02416407 
#> iter: 1  value: 1695.867 mgc: 0.003167257 ustep: 1 
#> iter: 2  value: 1695.867 mgc: 2.670868e-08 ustep: 1 
#> mgc: 1.181081e-13 
#> iter: 1  value: 1695.867 mgc: 0.003521271 ustep: 1 
#> iter: 2  value: 1695.867 mgc: 3.643047e-08 ustep: 1 
#> mgc: 1.440255e-13 
#> iter: 1  mgc: 1.440255e-13 
#> outer mgc:  0.02215541 
#> iter: 1  value: 1695.867 mgc: 0.006848257 ustep: 1 
#> iter: 2  value: 1695.867 mgc: 1.308287e-07 ustep: 1 
#> iter: 3  mgc: 1.22913e-13 
#> iter: 1  mgc: 1.22913e-13 
#> outer mgc:  0.0344119 
#> iter: 1  value: 1695.866 mgc: 0.004905411 ustep: 1 
#> iter: 2  value: 1695.866 mgc: 7.635996e-08 ustep: 1 
#> iter: 3  mgc: 1.477108e-13 
#> iter: 1  value: 1695.865 mgc: 0.009795729 ustep: 1 
#> iter: 2  value: 1695.865 mgc: 2.202768e-07 ustep: 1 
#> iter: 3  mgc: 1.435683e-13 
#> iter: 1  mgc: 1.435683e-13 
#> outer mgc:  0.08378627 
#> iter: 1  value: 1695.865 mgc: 0.008645525 ustep: 1 
#> iter: 2  value: 1695.865 mgc: 2.770362e-07 ustep: 1 
#> iter: 3  mgc: 1.251185e-13 
#> iter: 1  mgc: 1.251185e-13 
#> outer mgc:  0.04249942 
#> iter: 1  value: 1695.868 mgc: 0.006470028 ustep: 1 
#> iter: 2  value: 1695.868 mgc: 3.022633e-07 ustep: 1 
#> iter: 3  mgc: 1.177234e-13 
#> iter: 1  mgc: 1.177234e-13 
#> outer mgc:  0.03373898 
#> iter: 1  value: 1695.868 mgc: 0.006903179 ustep: 1 
#> iter: 2  value: 1695.868 mgc: 3.21453e-07 ustep: 1 
#> iter: 3  mgc: 1.828656e-13 
#> iter: 1  mgc: 1.828656e-13 
#> outer mgc:  0.036156 
#> iter: 1  value: 1695.867 mgc: 0.01062004 ustep: 1 
#> iter: 2  value: 1695.867 mgc: 3.087449e-07 ustep: 1 
#> iter: 3  mgc: 1.400357e-13 
#> iter: 1  mgc: 1.400357e-13 
#> outer mgc:  0.03451156 
#> iter: 1  value: 1695.866 mgc: 0.005199234 ustep: 1 
#> iter: 2  value: 1695.866 mgc: 1.463053e-07 ustep: 1 
#> iter: 3  mgc: 3.124534e-13 
#> iter: 1  value: 1695.864 mgc: 0.006436426 ustep: 1 
#> iter: 2  value: 1695.864 mgc: 2.038882e-07 ustep: 1 
#> iter: 3  mgc: 1.786626e-13 
#> iter: 1  mgc: 1.786626e-13 
#> outer mgc:  0.09317833 
#> iter: 1  value: 1695.864 mgc: 0.01460317 ustep: 1 
#> iter: 2  value: 1695.864 mgc: 6.167653e-07 ustep: 1 
#> iter: 3  mgc: 1.139058e-13 
#> iter: 1  mgc: 1.139058e-13 
#> outer mgc:  0.03732152 
#> iter: 1  value: 1695.865 mgc: 0.02665187 ustep: 1 
#> iter: 2  value: 1695.865 mgc: 4.44951e-06 ustep: 1 
#> iter: 3  mgc: 2.176942e-13 
#> iter: 1  value: 1695.865 mgc: 0.009621612 ustep: 1 
#> iter: 2  value: 1695.865 mgc: 2.123059e-07 ustep: 1 
#> iter: 3  mgc: 1.30427e-13 
#> iter: 1  mgc: 1.30427e-13 
#> outer mgc:  0.02029616 
#> iter: 1  value: 1695.865 mgc: 0.008482834 ustep: 1 
#> iter: 2  value: 1695.865 mgc: 2.272297e-07 ustep: 1 
#> iter: 3  mgc: 1.227189e-13 
#> iter: 1  mgc: 1.227189e-13 
#> outer mgc:  0.04165243 
#> iter: 1  value: 1695.865 mgc: 0.01652861 ustep: 1 
#> iter: 2  value: 1695.865 mgc: 3.596513e-07 ustep: 1 
#> iter: 3  mgc: 1.13435e-13 
#> iter: 1  value: 1695.865 mgc: 0.006992384 ustep: 1 
#> iter: 2  value: 1695.865 mgc: 7.21805e-08 ustep: 1 
#> iter: 3  mgc: 1.363659e-13 
#> iter: 1  mgc: 1.363659e-13 
#> outer mgc:  0.05244265 
#> iter: 1  value: 1695.865 mgc: 0.01200238 ustep: 1 
#> iter: 2  value: 1695.865 mgc: 1.967934e-07 ustep: 1 
#> iter: 3  mgc: 1.463101e-13 
#> iter: 1  value: 1695.865 mgc: 0.005011623 ustep: 1 
#> iter: 2  value: 1695.865 mgc: 3.255457e-08 ustep: 1 
#> mgc: 1.323247e-13 
#> iter: 1  mgc: 1.323247e-13 
#> outer mgc:  0.01211795 
#> iter: 1  value: 1695.864 mgc: 0.004601805 ustep: 1 
#> iter: 2  value: 1695.864 mgc: 6.777845e-08 ustep: 1 
#> mgc: 1.295935e-13 
#> iter: 1  mgc: 1.295935e-13 
#> outer mgc:  0.02245857 
#> iter: 1  value: 1695.864 mgc: 0.002286314 ustep: 1 
#> iter: 2  mgc: 8.577714e-09 
#> iter: 1  mgc: 8.577714e-09 
#> outer mgc:  0.006653785 
#> iter: 1  value: 1695.864 mgc: 0.003052986 ustep: 1 
#> iter: 2  value: 1695.864 mgc: 2.767648e-08 ustep: 1 
#> mgc: 1.240032e-13 
#> iter: 1  mgc: 1.240032e-13 
#> outer mgc:  0.006577472 
#> iter: 1  value: 1695.864 mgc: 0.001978969 ustep: 1 
#> iter: 2  value: 1695.864 mgc: 1.967498e-08 ustep: 1 
#> mgc: 1.410793e-13 
#> iter: 1  mgc: 1.410793e-13 
#> outer mgc:  0.01049392 
#> iter: 1  value: 1695.864 mgc: 0.003238365 ustep: 1 
#> iter: 2  value: 1695.864 mgc: 1.951677e-08 ustep: 1 
#> mgc: 1.283481e-13 
#> iter: 1  value: 1695.864 mgc: 0.001540859 ustep: 1 
#> iter: 2  mgc: 4.874413e-09 
#> iter: 1  mgc: 4.874413e-09 
#> outer mgc:  0.008904417 
#> iter: 1  value: 1695.864 mgc: 0.001223607 ustep: 1 
#> iter: 2  mgc: 2.037173e-09 
#> iter: 1  mgc: 2.037173e-09 
#> outer mgc:  0.00499252 
#> iter: 1  value: 1695.864 mgc: 0.002603208 ustep: 1 
#> iter: 2  mgc: 9.72495e-09 
#> iter: 1  value: 1695.864 mgc: 0.0008686503 ustep: 1 
#> iter: 2  mgc: 1.138724e-09 
#> iter: 1  mgc: 1.138724e-09 
#> outer mgc:  0.002614829 
#> iter: 1  value: 1695.864 mgc: 0.00047369 ustep: 1 
#> iter: 2  mgc: 3.877679e-10 
#> iter: 1  mgc: 3.877679e-10 
#> outer mgc:  0.001477719 
#> iter: 1  mgc: 3.877679e-10 
#> iter: 1  mgc: 3.877679e-10 
#> outer mgc:  0.001477719 
#> iter: 1  mgc: 3.877679e-10 
#> outer mgc:  0.001477719 
#> Retro Peel: 3
#> Using stats::nlminb for optimization with these control parameters:
#> iter.max: 1000
#> eval.max: 1000
#> iter: 1  value: 1660.162 mgc: 1.739049 ustep: 1 
#> iter: 2  value: 1660.161 mgc: 0.06256312 ustep: 1 
#> iter: 3  value: 1660.161 mgc: 0.0001051489 ustep: 1 
#> iter: 4  mgc: 2.995186e-10 
#> iter: 1  mgc: 2.995186e-10 
#> Matching hessian patterns... Done
#> outer mgc:  24.38513 
#> iter: 1  value: 9914.48 mgc: 181.7162 ustep: 0.006770861 
#> iter: 2  value: 8515.004 mgc: 130.3481 ustep: 0.01905181 
#> iter: 3  value: 6408.633 mgc: 545.1853 ustep: 0.07161494 
#> iter: 4  value: 5837.708 mgc: 186.3223 ustep: 0.2676829 
#> iter: 5  value: 5753.828 mgc: 40.02684 ustep: 0.5174291 
#> iter: 6  value: 5747.26 mgc: 5.862228 ustep: 0.7193535 
#> iter: 7  value: 5746.979 mgc: 2.082706 ustep: 0.8481623 
#> iter: 8  value: 5746.977 mgc: 0.1070772 ustep: 0.9209652 
#> iter: 9  value: 5746.977 mgc: 0.002726731 ustep: 0.9596733 
#> iter: 10  value: 5746.977 mgc: 9.148099e-05 ustep: 0.9796312 
#> iter: 11  value: 5746.977 mgc: 1.640038e-06 ustep: 0.9897642 
#> iter: 12  value: 5746.977 mgc: 1.457471e-08 ustep: 0.9948695 
#> iter: 13  mgc: 6.462196e-11 
#> iter: 1  value: 1697.613 mgc: 16.47185 ustep: 1 
#> iter: 2  value: 1697.46 mgc: 0.7683627 ustep: 1 
#> iter: 3  value: 1697.46 mgc: 0.02115546 ustep: 1 
#> iter: 4  value: 1697.46 mgc: 2.681852e-05 ustep: 1 
#> iter: 5  mgc: 4.187131e-11 
#> iter: 1  value: 1660.019 mgc: 1.607038 ustep: 1 
#> iter: 2  value: 1660.019 mgc: 0.005259189 ustep: 1 
#> iter: 3  value: 1660.019 mgc: 8.200153e-07 ustep: 1 
#> iter: 4  mgc: 1.180096e-13 
#> iter: 1  mgc: 1.180096e-13 
#> outer mgc:  28.09783 
#> iter: 1  value: 1659.566 mgc: 1.317059 ustep: 1 
#> iter: 2  value: 1659.566 mgc: 0.003004383 ustep: 1 
#> iter: 3  value: 1659.566 mgc: 3.924482e-07 ustep: 1 
#> iter: 4  mgc: 1.407232e-13 
#> iter: 1  mgc: 1.407232e-13 
#> outer mgc:  16.85105 
#> iter: 1  value: 1659.261 mgc: 0.7818364 ustep: 1 
#> iter: 2  value: 1659.261 mgc: 0.001502417 ustep: 1 
#> iter: 3  value: 1659.261 mgc: 6.091987e-08 ustep: 1 
#> iter: 4  mgc: 1.422694e-13 
#> iter: 1  mgc: 1.422694e-13 
#> outer mgc:  11.2041 
#> iter: 1  value: 1658.793 mgc: 1.227136 ustep: 1 
#> iter: 2  value: 1658.793 mgc: 0.002501899 ustep: 1 
#> iter: 3  value: 1658.793 mgc: 2.367107e-07 ustep: 1 
#> iter: 4  mgc: 1.716087e-13 
#> iter: 1  mgc: 1.716087e-13 
#> outer mgc:  11.62133 
#> iter: 1  value: 1658.253 mgc: 0.9730207 ustep: 1 
#> iter: 2  value: 1658.253 mgc: 0.002304107 ustep: 1 
#> iter: 3  value: 1658.253 mgc: 1.540546e-07 ustep: 1 
#> iter: 4  mgc: 1.863664e-13 
#> iter: 1  mgc: 1.863664e-13 
#> outer mgc:  10.18111 
#> iter: 1  value: 1657.898 mgc: 0.9681651 ustep: 1 
#> iter: 2  value: 1657.898 mgc: 0.001420973 ustep: 1 
#> iter: 3  value: 1657.898 mgc: 4.600111e-08 ustep: 1 
#> iter: 4  mgc: 1.46762e-13 
#> iter: 1  value: 1657.615 mgc: 0.5483193 ustep: 1 
#> iter: 2  value: 1657.615 mgc: 0.0005030225 ustep: 1 
#> iter: 3  mgc: 3.104234e-09 
#> iter: 1  mgc: 3.104234e-09 
#> outer mgc:  12.44792 
#> iter: 1  value: 1657.071 mgc: 1.696902 ustep: 1 
#> iter: 2  value: 1657.071 mgc: 0.00652284 ustep: 1 
#> iter: 3  value: 1657.071 mgc: 1.214098e-06 ustep: 1 
#> iter: 4  mgc: 1.625557e-13 
#> iter: 1  mgc: 1.625557e-13 
#> outer mgc:  12.57314 
#> iter: 1  value: 1656.697 mgc: 1.45376 ustep: 1 
#> iter: 2  value: 1656.697 mgc: 0.004248581 ustep: 1 
#> iter: 3  value: 1656.697 mgc: 4.885412e-07 ustep: 1 
#> iter: 4  mgc: 1.714849e-13 
#> iter: 1  mgc: 1.714849e-13 
#> outer mgc:  20.23588 
#> iter: 1  value: 1658.15 mgc: 3.982197 ustep: 1 
#> iter: 2  value: 1658.149 mgc: 0.01929418 ustep: 1 
#> iter: 3  value: 1658.149 mgc: 1.771065e-05 ustep: 1 
#> iter: 4  mgc: 1.888089e-11 
#> iter: 1  value: 1656.739 mgc: 2.21396 ustep: 1 
#> iter: 2  value: 1656.739 mgc: 0.008188205 ustep: 1 
#> iter: 3  value: 1656.739 mgc: 3.078447e-06 ustep: 1 
#> iter: 4  mgc: 5.836924e-13 
#> iter: 1  value: 1656.468 mgc: 1.087394 ustep: 1 
#> iter: 2  value: 1656.468 mgc: 0.001938834 ustep: 1 
#> iter: 3  value: 1656.468 mgc: 1.708627e-07 ustep: 1 
#> iter: 4  mgc: 1.793267e-13 
#> iter: 1  mgc: 1.793267e-13 
#> outer mgc:  6.688342 
#> iter: 1  value: 1656.414 mgc: 0.2338211 ustep: 1 
#> iter: 2  value: 1656.414 mgc: 0.0001120189 ustep: 1 
#> iter: 3  mgc: 5.19666e-10 
#> iter: 1  value: 1656.351 mgc: 0.08355742 ustep: 1 
#> iter: 2  value: 1656.351 mgc: 9.704277e-06 ustep: 1 
#> iter: 3  mgc: 2.831383e-12 
#> iter: 1  value: 1656.116 mgc: 0.3432962 ustep: 1 
#> iter: 2  value: 1656.116 mgc: 0.0001634379 ustep: 1 
#> iter: 3  mgc: 8.046814e-10 
#> iter: 1  value: 1655.657 mgc: 0.992183 ustep: 1 
#> iter: 2  value: 1655.657 mgc: 0.001352052 ustep: 1 
#> iter: 3  value: 1655.657 mgc: 5.546517e-08 ustep: 1 
#> iter: 4  mgc: 1.114297e-13 
#> iter: 1  mgc: 1.114297e-13 
#> outer mgc:  10.01052 
#> iter: 1  value: 1654.766 mgc: 1.876758 ustep: 1 
#> iter: 2  value: 1654.766 mgc: 0.004514297 ustep: 1 
#> iter: 3  value: 1654.766 mgc: 1.121159e-06 ustep: 1 
#> iter: 4  mgc: 2.621609e-13 
#> iter: 1  mgc: 2.621609e-13 
#> outer mgc:  19.17236 
#> iter: 1  value: 1654.65 mgc: 2.926877 ustep: 1 
#> iter: 2  value: 1654.65 mgc: 0.02221587 ustep: 1 
#> iter: 3  value: 1654.65 mgc: 9.412008e-06 ustep: 1 
#> iter: 4  mgc: 2.127044e-12 
#> iter: 1  mgc: 2.127044e-12 
#> outer mgc:  36.56584 
#> iter: 1  value: 1654.262 mgc: 2.897268 ustep: 1 
#> iter: 2  value: 1654.262 mgc: 0.01513136 ustep: 1 
#> iter: 3  value: 1654.262 mgc: 9.945592e-06 ustep: 1 
#> iter: 4  mgc: 3.982552e-12 
#> iter: 1  mgc: 3.982552e-12 
#> outer mgc:  18.9275 
#> iter: 1  value: 1654.132 mgc: 2.188687 ustep: 1 
#> iter: 2  value: 1654.132 mgc: 0.004804131 ustep: 1 
#> iter: 3  value: 1654.132 mgc: 2.727771e-06 ustep: 1 
#> iter: 4  mgc: 7.622305e-13 
#> iter: 1  mgc: 7.622305e-13 
#> outer mgc:  5.544 
#> iter: 1  value: 1654.475 mgc: 3.009488 ustep: 1 
#> iter: 2  value: 1654.475 mgc: 0.01627534 ustep: 1 
#> iter: 3  value: 1654.475 mgc: 1.563649e-05 ustep: 1 
#> iter: 4  mgc: 2.289556e-11 
#> iter: 1  value: 1653.896 mgc: 0.8071452 ustep: 1 
#> iter: 2  value: 1653.896 mgc: 0.001150155 ustep: 1 
#> iter: 3  value: 1653.896 mgc: 7.477588e-08 ustep: 1 
#> iter: 4  mgc: 1.750517e-13 
#> iter: 1  mgc: 1.750517e-13 
#> outer mgc:  6.691793 
#> iter: 1  value: 1654.288 mgc: 1.863785 ustep: 1 
#> iter: 2  value: 1654.288 mgc: 0.004024192 ustep: 1 
#> iter: 3  value: 1654.288 mgc: 6.100577e-07 ustep: 1 
#> iter: 4  mgc: 1.052181e-13 
#> iter: 1  value: 1653.9 mgc: 0.4887092 ustep: 1 
#> iter: 2  value: 1653.9 mgc: 0.0002893209 ustep: 1 
#> iter: 3  mgc: 3.339936e-09 
#> iter: 1  mgc: 3.339936e-09 
#> outer mgc:  8.2992 
#> iter: 1  value: 1653.833 mgc: 0.3924005 ustep: 1 
#> iter: 2  value: 1653.833 mgc: 0.0003037501 ustep: 1 
#> iter: 3  mgc: 2.393367e-09 
#> iter: 1  mgc: 2.393367e-09 
#> outer mgc:  3.811512 
#> iter: 1  value: 1653.775 mgc: 0.3748823 ustep: 1 
#> iter: 2  value: 1653.775 mgc: 0.0001564886 ustep: 1 
#> iter: 3  mgc: 1.634251e-09 
#> iter: 1  mgc: 1.634251e-09 
#> outer mgc:  5.34693 
#> iter: 1  value: 1653.657 mgc: 0.3647941 ustep: 1 
#> iter: 2  value: 1653.657 mgc: 0.0002515794 ustep: 1 
#> iter: 3  mgc: 2.3762e-09 
#> iter: 1  mgc: 2.3762e-09 
#> outer mgc:  5.043068 
#> iter: 1  value: 1653.602 mgc: 0.19408 ustep: 1 
#> iter: 2  value: 1653.602 mgc: 0.0001316347 ustep: 1 
#> iter: 3  mgc: 3.004188e-10 
#> iter: 1  mgc: 3.004188e-10 
#> outer mgc:  3.264645 
#> iter: 1  value: 1653.473 mgc: 0.3086963 ustep: 1 
#> iter: 2  value: 1653.473 mgc: 0.0001553384 ustep: 1 
#> iter: 3  mgc: 1.091291e-09 
#> iter: 1  value: 1653.265 mgc: 0.3538891 ustep: 1 
#> iter: 2  value: 1653.265 mgc: 0.0002874586 ustep: 1 
#> iter: 3  mgc: 1.233766e-09 
#> iter: 1  mgc: 1.233766e-09 
#> outer mgc:  4.945986 
#> iter: 1  value: 1653.37 mgc: 0.5725921 ustep: 1 
#> iter: 2  value: 1653.37 mgc: 0.001067374 ustep: 1 
#> iter: 3  value: 1653.37 mgc: 1.487372e-08 ustep: 1 
#> mgc: 1.852862e-13 
#> iter: 1  mgc: 1.852862e-13 
#> outer mgc:  4.271681 
#> iter: 1  value: 1652.636 mgc: 1.275969 ustep: 1 
#> iter: 2  value: 1652.636 mgc: 0.003269631 ustep: 1 
#> iter: 3  value: 1652.636 mgc: 3.414535e-07 ustep: 1 
#> iter: 4  mgc: 1.559747e-13 
#> iter: 1  mgc: 1.559747e-13 
#> outer mgc:  14.37489 
#> iter: 1  value: 1652.818 mgc: 1.987861 ustep: 1 
#> iter: 2  value: 1652.818 mgc: 0.01086056 ustep: 1 
#> iter: 3  value: 1652.818 mgc: 3.104379e-06 ustep: 1 
#> iter: 4  mgc: 4.104632e-13 
#> iter: 1  value: 1652.693 mgc: 1.463631 ustep: 1 
#> iter: 2  value: 1652.693 mgc: 0.004216778 ustep: 1 
#> iter: 3  value: 1652.693 mgc: 5.411768e-07 ustep: 1 
#> iter: 4  mgc: 2.260689e-13 
#> iter: 1  value: 1652.718 mgc: 1.335029 ustep: 1 
#> iter: 2  value: 1652.718 mgc: 0.002993286 ustep: 1 
#> iter: 3  value: 1652.718 mgc: 3.225041e-07 ustep: 1 
#> iter: 4  mgc: 1.856752e-13 
#> iter: 1  value: 1652.621 mgc: 0.7811304 ustep: 1 
#> iter: 2  value: 1652.621 mgc: 0.000971665 ustep: 1 
#> iter: 3  value: 1652.621 mgc: 3.588911e-08 ustep: 1 
#> iter: 4  mgc: 1.257135e-13 
#> iter: 1  mgc: 1.257135e-13 
#> outer mgc:  6.724772 
#> iter: 1  value: 1652.629 mgc: 0.2809591 ustep: 1 
#> iter: 2  value: 1652.629 mgc: 0.0001338346 ustep: 1 
#> iter: 3  mgc: 7.9891e-10 
#> iter: 1  mgc: 7.9891e-10 
#> outer mgc:  3.12776 
#> iter: 1  value: 1652.705 mgc: 0.3350782 ustep: 1 
#> iter: 2  value: 1652.705 mgc: 0.0002660571 ustep: 1 
#> iter: 3  mgc: 1.428098e-09 
#> iter: 1  mgc: 1.428098e-09 
#> outer mgc:  5.992131 
#> iter: 1  value: 1652.783 mgc: 0.3923535 ustep: 1 
#> iter: 2  value: 1652.783 mgc: 0.0003242147 ustep: 1 
#> iter: 3  mgc: 3.091384e-09 
#> iter: 1  mgc: 3.091384e-09 
#> outer mgc:  3.050125 
#> iter: 1  value: 1652.651 mgc: 0.2088892 ustep: 1 
#> iter: 2  value: 1652.651 mgc: 0.0001270426 ustep: 1 
#> iter: 3  mgc: 4.082906e-10 
#> iter: 1  value: 1652.511 mgc: 0.1314302 ustep: 1 
#> iter: 2  value: 1652.511 mgc: 3.670278e-05 ustep: 1 
#> iter: 3  mgc: 3.858801e-11 
#> iter: 1  mgc: 3.858801e-11 
#> outer mgc:  1.828394 
#> iter: 1  value: 1652.681 mgc: 0.546229 ustep: 1 
#> iter: 2  value: 1652.681 mgc: 0.0005652798 ustep: 1 
#> iter: 3  value: 1652.681 mgc: 1.848865e-08 ustep: 1 
#> mgc: 1.857617e-13 
#> iter: 1  mgc: 1.857617e-13 
#> outer mgc:  5.313852 
#> iter: 1  value: 1653.025 mgc: 1.08047 ustep: 1 
#> iter: 2  value: 1653.025 mgc: 0.002292144 ustep: 1 
#> iter: 3  value: 1653.025 mgc: 2.410886e-07 ustep: 1 
#> iter: 4  mgc: 1.557304e-13 
#> iter: 1  value: 1652.857 mgc: 0.8057333 ustep: 1 
#> iter: 2  value: 1652.857 mgc: 0.00114898 ustep: 1 
#> iter: 3  value: 1652.857 mgc: 5.910634e-08 ustep: 1 
#> iter: 4  mgc: 9.915563e-14 
#> iter: 1  value: 1652.79 mgc: 0.7025991 ustep: 1 
#> iter: 2  value: 1652.79 mgc: 0.000858111 ustep: 1 
#> iter: 3  value: 1652.79 mgc: 2.618768e-08 ustep: 1 
#> iter: 4  mgc: 1.303693e-13 
#> iter: 1  value: 1652.747 mgc: 0.6649784 ustep: 1 
#> iter: 2  value: 1652.747 mgc: 0.0007038633 ustep: 1 
#> iter: 3  value: 1652.747 mgc: 1.987463e-08 ustep: 1 
#> mgc: 1.108709e-13 
#> iter: 1  value: 1652.681 mgc: 0.2154127 ustep: 1 
#> iter: 2  value: 1652.681 mgc: 7.431528e-05 ustep: 1 
#> iter: 3  mgc: 2.23098e-10 
#> iter: 1  mgc: 2.23098e-10 
#> outer mgc:  1.730168 
#> iter: 1  value: 1652.676 mgc: 0.111575 ustep: 1 
#> iter: 2  value: 1652.676 mgc: 2.308636e-05 ustep: 1 
#> iter: 3  mgc: 1.415781e-11 
#> iter: 1  mgc: 1.415781e-11 
#> outer mgc:  1.717088 
#> iter: 1  value: 1652.684 mgc: 0.1438614 ustep: 1 
#> iter: 2  value: 1652.684 mgc: 2.925683e-05 ustep: 1 
#> iter: 3  mgc: 4.176897e-11 
#> iter: 1  mgc: 4.176897e-11 
#> outer mgc:  2.480301 
#> iter: 1  value: 1652.678 mgc: 0.1322145 ustep: 1 
#> iter: 2  value: 1652.678 mgc: 2.564255e-05 ustep: 1 
#> iter: 3  mgc: 2.70718e-11 
#> iter: 1  mgc: 2.70718e-11 
#> outer mgc:  1.708255 
#> iter: 1  value: 1652.683 mgc: 0.104819 ustep: 1 
#> iter: 2  value: 1652.683 mgc: 1.453595e-05 ustep: 1 
#> iter: 3  mgc: 1.552664e-11 
#> iter: 1  value: 1652.691 mgc: 0.02359542 ustep: 1 
#> iter: 2  value: 1652.691 mgc: 2.327997e-06 ustep: 1 
#> iter: 3  mgc: 2.574619e-13 
#> iter: 1  value: 1652.727 mgc: 0.1035312 ustep: 1 
#> iter: 2  value: 1652.727 mgc: 4.482924e-05 ustep: 1 
#> iter: 3  mgc: 8.156104e-11 
#> iter: 1  value: 1652.867 mgc: 0.2905005 ustep: 1 
#> iter: 2  value: 1652.867 mgc: 0.0003532522 ustep: 1 
#> iter: 3  mgc: 5.083919e-09 
#> iter: 1  mgc: 5.083919e-09 
#> outer mgc:  4.000389 
#> iter: 1  value: 1653.632 mgc: 0.9635297 ustep: 1 
#> iter: 2  value: 1653.632 mgc: 0.0008505488 ustep: 1 
#> iter: 3  value: 1653.632 mgc: 2.354733e-08 ustep: 1 
#> mgc: 1.421835e-13 
#> iter: 1  mgc: 1.421835e-13 
#> outer mgc:  5.563771 
#> iter: 1  value: 1653.003 mgc: 1.653677 ustep: 1 
#> iter: 2  value: 1653.003 mgc: 0.003579809 ustep: 1 
#> iter: 3  value: 1653.003 mgc: 1.87433e-07 ustep: 1 
#> iter: 4  mgc: 1.632062e-13 
#> iter: 1  mgc: 1.632062e-13 
#> outer mgc:  12.00874 
#> iter: 1  value: 1652.503 mgc: 1.376193 ustep: 1 
#> iter: 2  value: 1652.503 mgc: 0.003287836 ustep: 1 
#> iter: 3  value: 1652.503 mgc: 1.017647e-07 ustep: 1 
#> iter: 4  mgc: 1.743036e-13 
#> iter: 1  value: 1652.744 mgc: 1.322865 ustep: 1 
#> iter: 2  value: 1652.744 mgc: 0.001669545 ustep: 1 
#> iter: 3  value: 1652.744 mgc: 1.073554e-07 ustep: 1 
#> iter: 4  mgc: 1.368203e-13 
#> iter: 1  value: 1652.924 mgc: 1.472701 ustep: 1 
#> iter: 2  value: 1652.924 mgc: 0.002328418 ustep: 1 
#> iter: 3  value: 1652.924 mgc: 3.222527e-07 ustep: 1 
#> iter: 4  mgc: 1.160869e-13 
#> iter: 1  value: 1653.054 mgc: 1.574607 ustep: 1 
#> iter: 2  value: 1653.054 mgc: 0.003657503 ustep: 1 
#> iter: 3  value: 1653.054 mgc: 5.895973e-07 ustep: 1 
#> iter: 4  mgc: 1.500388e-13 
#> iter: 1  value: 1652.919 mgc: 0.6522838 ustep: 1 
#> iter: 2  value: 1652.919 mgc: 0.0006613937 ustep: 1 
#> iter: 3  value: 1652.919 mgc: 1.80324e-08 ustep: 1 
#> mgc: 1.374883e-13 
#> iter: 1  mgc: 1.374883e-13 
#> outer mgc:  5.90405 
#> iter: 1  value: 1652.896 mgc: 0.1950827 ustep: 1 
#> iter: 2  value: 1652.896 mgc: 8.168462e-05 ustep: 1 
#> iter: 3  mgc: 2.007958e-10 
#> iter: 1  mgc: 2.007958e-10 
#> outer mgc:  3.396458 
#> iter: 1  value: 1652.832 mgc: 0.2428223 ustep: 1 
#> iter: 2  value: 1652.832 mgc: 5.631016e-05 ustep: 1 
#> iter: 3  mgc: 1.41347e-10 
#> iter: 1  mgc: 1.41347e-10 
#> outer mgc:  5.420406 
#> iter: 1  value: 1652.762 mgc: 0.2367603 ustep: 1 
#> iter: 2  value: 1652.762 mgc: 0.0001356777 ustep: 1 
#> iter: 3  mgc: 6.801721e-10 
#> iter: 1  mgc: 6.801721e-10 
#> outer mgc:  3.445669 
#> iter: 1  value: 1652.615 mgc: 0.1720933 ustep: 1 
#> iter: 2  value: 1652.615 mgc: 5.893036e-05 ustep: 1 
#> iter: 3  mgc: 5.685957e-11 
#> iter: 1  mgc: 5.685957e-11 
#> outer mgc:  4.297491 
#> iter: 1  value: 1652.656 mgc: 0.246465 ustep: 1 
#> iter: 2  value: 1652.656 mgc: 0.000186586 ustep: 1 
#> iter: 3  mgc: 1.557722e-09 
#> iter: 1  mgc: 1.557722e-09 
#> outer mgc:  2.452996 
#> iter: 1  value: 1652.293 mgc: 0.9194044 ustep: 1 
#> iter: 2  value: 1652.293 mgc: 0.001526045 ustep: 1 
#> iter: 3  mgc: 7.427058e-09 
#> iter: 1  mgc: 7.427058e-09 
#> outer mgc:  2.226423 
#> iter: 1  value: 1652.901 mgc: 1.793044 ustep: 1 
#> iter: 2  value: 1652.901 mgc: 0.006134042 ustep: 1 
#> iter: 3  value: 1652.901 mgc: 7.113759e-07 ustep: 1 
#> iter: 4  mgc: 1.232455e-13 
#> iter: 1  value: 1652.292 mgc: 0.3644724 ustep: 1 
#> iter: 2  value: 1652.292 mgc: 0.0002258483 ustep: 1 
#> iter: 3  mgc: 1.741309e-09 
#> iter: 1  value: 1652.287 mgc: 0.0815314 ustep: 1 
#> iter: 2  value: 1652.287 mgc: 1.127903e-05 ustep: 1 
#> iter: 3  mgc: 4.367961e-12 
#> iter: 1  mgc: 4.367961e-12 
#> outer mgc:  2.474477 
#> iter: 1  value: 1652.289 mgc: 0.08541084 ustep: 1 
#> iter: 2  value: 1652.289 mgc: 1.205316e-05 ustep: 1 
#> iter: 3  mgc: 5.409703e-12 
#> iter: 1  mgc: 5.409703e-12 
#> outer mgc:  1.983266 
#> iter: 1  value: 1652.286 mgc: 0.08030372 ustep: 1 
#> iter: 2  value: 1652.286 mgc: 1.248667e-05 ustep: 1 
#> iter: 3  mgc: 3.0281e-12 
#> iter: 1  mgc: 3.0281e-12 
#> outer mgc:  2.502289 
#> iter: 1  value: 1652.291 mgc: 0.1036738 ustep: 1 
#> iter: 2  value: 1652.291 mgc: 1.713074e-05 ustep: 1 
#> iter: 3  mgc: 1.503813e-11 
#> iter: 1  mgc: 1.503813e-11 
#> outer mgc:  1.887525 
#> iter: 1  value: 1652.287 mgc: 0.09131208 ustep: 1 
#> iter: 2  value: 1652.287 mgc: 1.092975e-05 ustep: 1 
#> iter: 3  mgc: 1.153982e-12 
#> iter: 1  mgc: 1.153982e-12 
#> outer mgc:  2.48907 
#> iter: 1  value: 1652.263 mgc: 0.2156709 ustep: 1 
#> iter: 2  value: 1652.263 mgc: 8.34942e-05 ustep: 1 
#> iter: 3  mgc: 2.512247e-10 
#> iter: 1  mgc: 2.512247e-10 
#> outer mgc:  1.80821 
#> iter: 1  value: 1652.273 mgc: 0.3839344 ustep: 1 
#> iter: 2  value: 1652.273 mgc: 8.274326e-05 ustep: 1 
#> iter: 3  mgc: 9.034577e-11 
#> iter: 1  mgc: 9.034577e-11 
#> outer mgc:  1.71563 
#> iter: 1  value: 1652.305 mgc: 0.6779814 ustep: 1 
#> iter: 2  value: 1652.305 mgc: 0.0009065435 ustep: 1 
#> iter: 3  value: 1652.305 mgc: 2.255441e-08 ustep: 1 
#> mgc: 1.235686e-13 
#> iter: 1  value: 1652.257 mgc: 0.1589682 ustep: 1 
#> iter: 2  value: 1652.257 mgc: 4.936603e-05 ustep: 1 
#> iter: 3  mgc: 6.685252e-11 
#> iter: 1  mgc: 6.685252e-11 
#> outer mgc:  3.645786 
#> iter: 1  value: 1652.239 mgc: 0.2375568 ustep: 1 
#> iter: 2  value: 1652.239 mgc: 0.000100368 ustep: 1 
#> iter: 3  mgc: 3.791695e-10 
#> iter: 1  mgc: 3.791695e-10 
#> outer mgc:  2.117883 
#> iter: 1  value: 1652.212 mgc: 0.05898058 ustep: 1 
#> iter: 2  value: 1652.212 mgc: 2.669921e-06 ustep: 1 
#> iter: 3  mgc: 9.830287e-13 
#> iter: 1  mgc: 9.830287e-13 
#> outer mgc:  1.326765 
#> iter: 1  value: 1652.18 mgc: 0.1129541 ustep: 1 
#> iter: 2  value: 1652.18 mgc: 5.457035e-05 ustep: 1 
#> iter: 3  mgc: 8.084323e-11 
#> iter: 1  mgc: 8.084323e-11 
#> outer mgc:  1.800182 
#> iter: 1  value: 1652.137 mgc: 0.1572963 ustep: 1 
#> iter: 2  value: 1652.137 mgc: 8.13345e-05 ustep: 1 
#> iter: 3  mgc: 3.081614e-10 
#> iter: 1  mgc: 3.081614e-10 
#> outer mgc:  1.201041 
#> iter: 1  value: 1652.222 mgc: 0.7667194 ustep: 1 
#> iter: 2  value: 1652.222 mgc: 0.0008119677 ustep: 1 
#> iter: 3  value: 1652.222 mgc: 1.42913e-08 ustep: 1 
#> mgc: 1.288085e-13 
#> iter: 1  value: 1652.15 mgc: 0.278987 ustep: 1 
#> iter: 2  value: 1652.15 mgc: 0.0001223572 ustep: 1 
#> iter: 3  mgc: 5.667221e-10 
#> iter: 1  value: 1652.139 mgc: 0.08231648 ustep: 1 
#> iter: 2  value: 1652.139 mgc: 1.063785e-05 ustep: 1 
#> iter: 3  mgc: 4.321987e-12 
#> iter: 1  mgc: 4.321987e-12 
#> outer mgc:  2.183032 
#> iter: 1  value: 1652.142 mgc: 0.07288558 ustep: 1 
#> iter: 2  value: 1652.142 mgc: 9.301721e-06 ustep: 1 
#> iter: 3  mgc: 3.852374e-12 
#> iter: 1  mgc: 3.852374e-12 
#> outer mgc:  1.098786 
#> iter: 1  value: 1652.142 mgc: 0.06272221 ustep: 1 
#> iter: 2  value: 1652.142 mgc: 5.396433e-06 ustep: 1 
#> iter: 3  mgc: 5.851056e-13 
#> iter: 1  mgc: 5.851056e-13 
#> outer mgc:  1.747328 
#> iter: 1  value: 1652.144 mgc: 0.05842115 ustep: 1 
#> iter: 2  value: 1652.144 mgc: 1.441796e-05 ustep: 1 
#> iter: 3  mgc: 5.640327e-12 
#> iter: 1  value: 1652.144 mgc: 0.01559234 ustep: 1 
#> iter: 2  value: 1652.144 mgc: 2.557527e-06 ustep: 1 
#> iter: 3  mgc: 1.648922e-13 
#> iter: 1  value: 1652.142 mgc: 0.0572453 ustep: 1 
#> iter: 2  value: 1652.142 mgc: 3.443407e-05 ustep: 1 
#> iter: 3  mgc: 3.068997e-11 
#> iter: 1  value: 1652.144 mgc: 0.1674263 ustep: 1 
#> iter: 2  value: 1652.144 mgc: 0.000293589 ustep: 1 
#> iter: 3  mgc: 2.22861e-09 
#> iter: 1  mgc: 2.22861e-09 
#> outer mgc:  1.846488 
#> iter: 1  value: 1652.235 mgc: 0.7923631 ustep: 1 
#> iter: 2  value: 1652.235 mgc: 0.001376173 ustep: 1 
#> iter: 3  value: 1652.235 mgc: 5.301756e-08 ustep: 1 
#> mgc: 1.299384e-13 
#> iter: 1  mgc: 1.299384e-13 
#> outer mgc:  2.813165 
#> iter: 1  value: 1652.278 mgc: 0.2143628 ustep: 1 
#> iter: 2  value: 1652.278 mgc: 0.0002218389 ustep: 1 
#> iter: 3  mgc: 1.142298e-09 
#> iter: 1  mgc: 1.142298e-09 
#> outer mgc:  3.826962 
#> iter: 1  value: 1652.261 mgc: 0.3969603 ustep: 1 
#> iter: 2  value: 1652.261 mgc: 9.713983e-05 ustep: 1 
#> iter: 3  mgc: 2.234214e-10 
#> iter: 1  mgc: 2.234214e-10 
#> outer mgc:  4.468326 
#> iter: 1  value: 1652.187 mgc: 0.4408782 ustep: 1 
#> iter: 2  value: 1652.187 mgc: 0.0001103767 ustep: 1 
#> iter: 3  mgc: 3.33011e-10 
#> iter: 1  mgc: 3.33011e-10 
#> outer mgc:  1.173397 
#> iter: 1  value: 1652.212 mgc: 0.6505903 ustep: 1 
#> iter: 2  value: 1652.212 mgc: 0.001185388 ustep: 1 
#> iter: 3  value: 1652.212 mgc: 4.538318e-08 ustep: 1 
#> mgc: 1.343576e-13 
#> iter: 1  value: 1652.207 mgc: 0.4999843 ustep: 1 
#> iter: 2  value: 1652.207 mgc: 0.0005279975 ustep: 1 
#> iter: 3  mgc: 8.515819e-09 
#> iter: 1  value: 1652.202 mgc: 0.4151393 ustep: 1 
#> iter: 2  value: 1652.202 mgc: 0.0002629252 ustep: 1 
#> iter: 3  mgc: 2.294855e-09 
#> iter: 1  value: 1652.188 mgc: 0.1234052 ustep: 1 
#> iter: 2  value: 1652.188 mgc: 2.327715e-05 ustep: 1 
#> iter: 3  mgc: 1.80811e-11 
#> iter: 1  mgc: 1.80811e-11 
#> outer mgc:  1.757067 
#> iter: 1  value: 1652.184 mgc: 0.06735299 ustep: 1 
#> iter: 2  value: 1652.184 mgc: 6.862922e-06 ustep: 1 
#> iter: 3  mgc: 1.638849e-12 
#> iter: 1  mgc: 1.638849e-12 
#> outer mgc:  0.9636 
#> iter: 1  value: 1652.185 mgc: 0.06093884 ustep: 1 
#> iter: 2  value: 1652.185 mgc: 8.453835e-06 ustep: 1 
#> iter: 3  mgc: 3.0202e-12 
#> iter: 1  mgc: 3.0202e-12 
#> outer mgc:  1.015449 
#> iter: 1  value: 1652.185 mgc: 0.04465896 ustep: 1 
#> iter: 2  value: 1652.185 mgc: 2.098148e-06 ustep: 1 
#> iter: 3  mgc: 1.327376e-13 
#> iter: 1  value: 1652.188 mgc: 0.02603233 ustep: 1 
#> iter: 2  value: 1652.188 mgc: 4.953271e-06 ustep: 1 
#> iter: 3  mgc: 8.822106e-13 
#> iter: 1  value: 1652.199 mgc: 0.08416724 ustep: 1 
#> iter: 2  value: 1652.199 mgc: 5.172701e-05 ustep: 1 
#> iter: 3  mgc: 9.56538e-11 
#> iter: 1  mgc: 9.56538e-11 
#> outer mgc:  1.069055 
#> iter: 1  value: 1652.235 mgc: 0.2525451 ustep: 1 
#> iter: 2  value: 1652.235 mgc: 0.000197924 ustep: 1 
#> iter: 3  mgc: 4.588663e-10 
#> iter: 1  mgc: 4.588663e-10 
#> outer mgc:  1.084323 
#> iter: 1  value: 1652.286 mgc: 0.2853324 ustep: 1 
#> iter: 2  value: 1652.286 mgc: 0.0002283913 ustep: 1 
#> iter: 3  mgc: 1.206575e-09 
#> iter: 1  value: 1652.243 mgc: 0.1589106 ustep: 1 
#> iter: 2  value: 1652.243 mgc: 5.357856e-05 ustep: 1 
#> iter: 3  mgc: 1.039492e-11 
#> iter: 1  mgc: 1.039492e-11 
#> outer mgc:  1.57241 
#> iter: 1  value: 1652.238 mgc: 0.1638282 ustep: 1 
#> iter: 2  value: 1652.238 mgc: 0.000127579 ustep: 1 
#> iter: 3  mgc: 7.4691e-10 
#> iter: 1  mgc: 7.4691e-10 
#> outer mgc:  1.343223 
#> iter: 1  value: 1652.155 mgc: 0.1553342 ustep: 1 
#> iter: 2  value: 1652.155 mgc: 4.478357e-05 ustep: 1 
#> iter: 3  mgc: 1.063226e-10 
#> iter: 1  mgc: 1.063226e-10 
#> outer mgc:  0.6995568 
#> iter: 1  value: 1652.175 mgc: 0.3347661 ustep: 1 
#> iter: 2  value: 1652.175 mgc: 0.0002127408 ustep: 1 
#> iter: 3  mgc: 1.633525e-09 
#> iter: 1  value: 1652.162 mgc: 0.2887328 ustep: 1 
#> iter: 2  value: 1652.162 mgc: 0.0001472538 ustep: 1 
#> iter: 3  mgc: 8.573059e-10 
#> iter: 1  value: 1652.154 mgc: 0.151824 ustep: 1 
#> iter: 2  value: 1652.154 mgc: 3.884484e-05 ustep: 1 
#> iter: 3  mgc: 4.961815e-11 
#> iter: 1  value: 1652.154 mgc: 0.03795164 ustep: 1 
#> iter: 2  value: 1652.154 mgc: 2.42533e-06 ustep: 1 
#> iter: 3  mgc: 1.930998e-13 
#> iter: 1  mgc: 1.930998e-13 
#> outer mgc:  0.7012656 
#> iter: 1  value: 1652.156 mgc: 0.0241489 ustep: 1 
#> iter: 2  value: 1652.156 mgc: 9.755904e-07 ustep: 1 
#> iter: 3  mgc: 1.782512e-13 
#> iter: 1  mgc: 1.782512e-13 
#> outer mgc:  0.503434 
#> iter: 1  value: 1652.156 mgc: 0.02966718 ustep: 1 
#> iter: 2  value: 1652.156 mgc: 1.56747e-06 ustep: 1 
#> iter: 3  mgc: 1.415593e-13 
#> iter: 1  mgc: 1.415593e-13 
#> outer mgc:  0.6130194 
#> iter: 1  value: 1652.158 mgc: 0.01996093 ustep: 1 
#> iter: 2  value: 1652.158 mgc: 5.667311e-07 ustep: 1 
#> iter: 3  mgc: 1.231973e-13 
#> iter: 1  mgc: 1.231973e-13 
#> outer mgc:  0.4172714 
#> iter: 1  value: 1652.159 mgc: 0.03035046 ustep: 1 
#> iter: 2  value: 1652.159 mgc: 1.423427e-06 ustep: 1 
#> iter: 3  mgc: 1.447631e-13 
#> iter: 1  value: 1652.161 mgc: 0.02064476 ustep: 1 
#> iter: 2  value: 1652.161 mgc: 5.386402e-07 ustep: 1 
#> iter: 3  mgc: 1.767812e-13 
#> iter: 1  value: 1652.166 mgc: 0.05114471 ustep: 1 
#> iter: 2  value: 1652.166 mgc: 3.306085e-06 ustep: 1 
#> iter: 3  mgc: 2.925172e-13 
#> iter: 1  mgc: 2.925172e-13 
#> outer mgc:  0.78118 
#> iter: 1  value: 1652.176 mgc: 0.1270073 ustep: 1 
#> iter: 2  value: 1652.176 mgc: 1.521169e-05 ustep: 1 
#> iter: 3  mgc: 1.442996e-11 
#> iter: 1  mgc: 1.442996e-11 
#> outer mgc:  0.3123026 
#> iter: 1  value: 1652.203 mgc: 0.1658159 ustep: 1 
#> iter: 2  value: 1652.203 mgc: 5.672341e-05 ustep: 1 
#> iter: 3  mgc: 9.615865e-11 
#> iter: 1  mgc: 9.615865e-11 
#> outer mgc:  0.2229302 
#> iter: 1  value: 1652.208 mgc: 0.1760762 ustep: 1 
#> iter: 2  value: 1652.208 mgc: 2.882263e-05 ustep: 1 
#> iter: 3  mgc: 5.882996e-12 
#> iter: 1  mgc: 5.882996e-12 
#> outer mgc:  1.218025 
#> iter: 1  value: 1652.192 mgc: 0.1038036 ustep: 1 
#> iter: 2  value: 1652.192 mgc: 5.640718e-05 ustep: 1 
#> iter: 3  mgc: 7.633978e-11 
#> iter: 1  mgc: 7.633978e-11 
#> outer mgc:  1.006259 
#> iter: 1  value: 1652.193 mgc: 0.1326361 ustep: 1 
#> iter: 2  value: 1652.193 mgc: 5.280318e-05 ustep: 1 
#> iter: 3  mgc: 1.709665e-11 
#> iter: 1  mgc: 1.709665e-11 
#> outer mgc:  0.8270409 
#> iter: 1  value: 1652.18 mgc: 0.2504228 ustep: 1 
#> iter: 2  value: 1652.18 mgc: 0.0001594643 ustep: 1 
#> iter: 3  mgc: 7.689719e-11 
#> iter: 1  value: 1652.186 mgc: 0.1508837 ustep: 1 
#> iter: 2  value: 1652.186 mgc: 4.12191e-05 ustep: 1 
#> iter: 3  mgc: 4.709313e-11 
#> iter: 1  value: 1652.188 mgc: 0.1582533 ustep: 1 
#> iter: 2  value: 1652.188 mgc: 3.280544e-05 ustep: 1 
#> iter: 3  mgc: 5.774775e-11 
#> iter: 1  value: 1652.191 mgc: 0.0333297 ustep: 1 
#> iter: 2  value: 1652.191 mgc: 1.457864e-06 ustep: 1 
#> iter: 3  mgc: 1.43568e-13 
#> iter: 1  mgc: 1.43568e-13 
#> outer mgc:  0.3707 
#> iter: 1  value: 1652.192 mgc: 0.03131916 ustep: 1 
#> iter: 2  value: 1652.192 mgc: 1.767039e-06 ustep: 1 
#> iter: 3  mgc: 1.946355e-13 
#> iter: 1  mgc: 1.946355e-13 
#> outer mgc:  0.4569612 
#> iter: 1  value: 1652.191 mgc: 0.01916226 ustep: 1 
#> iter: 2  value: 1652.191 mgc: 5.321929e-07 ustep: 1 
#> iter: 3  mgc: 8.439706e-14 
#> iter: 1  mgc: 8.439706e-14 
#> outer mgc:  0.243618 
#> iter: 1  value: 1652.191 mgc: 0.02206295 ustep: 1 
#> iter: 2  value: 1652.191 mgc: 1.210738e-06 ustep: 1 
#> iter: 3  mgc: 1.755624e-13 
#> iter: 1  mgc: 1.755624e-13 
#> outer mgc:  0.2722199 
#> iter: 1  value: 1652.19 mgc: 0.02480247 ustep: 1 
#> iter: 2  value: 1652.19 mgc: 1.53563e-06 ustep: 1 
#> iter: 3  mgc: 1.39054e-13 
#> iter: 1  mgc: 1.39054e-13 
#> outer mgc:  0.2010623 
#> iter: 1  value: 1652.191 mgc: 0.01971292 ustep: 1 
#> iter: 2  value: 1652.191 mgc: 1.705005e-06 ustep: 1 
#> iter: 3  mgc: 1.31356e-13 
#> iter: 1  value: 1652.192 mgc: 0.03130598 ustep: 1 
#> iter: 2  value: 1652.192 mgc: 2.919067e-06 ustep: 1 
#> iter: 3  mgc: 2.053079e-13 
#> iter: 1  mgc: 2.053079e-13 
#> outer mgc:  0.1608554 
#> iter: 1  value: 1652.186 mgc: 0.03319492 ustep: 1 
#> iter: 2  value: 1652.186 mgc: 3.879925e-06 ustep: 1 
#> iter: 3  mgc: 2.501301e-13 
#> iter: 1  value: 1652.179 mgc: 0.05960217 ustep: 1 
#> iter: 2  value: 1652.179 mgc: 1.013507e-05 ustep: 1 
#> iter: 3  mgc: 1.733876e-12 
#> iter: 1  mgc: 1.733876e-12 
#> outer mgc:  0.5617927 
#> iter: 1  value: 1652.184 mgc: 0.05695843 ustep: 1 
#> iter: 2  value: 1652.184 mgc: 6.080493e-06 ustep: 1 
#> iter: 3  mgc: 1.258124e-12 
#> iter: 1  mgc: 1.258124e-12 
#> outer mgc:  0.129946 
#> iter: 1  value: 1652.186 mgc: 0.06922619 ustep: 1 
#> iter: 2  value: 1652.186 mgc: 9.713426e-06 ustep: 1 
#> iter: 3  mgc: 1.295118e-12 
#> iter: 1  mgc: 1.295118e-12 
#> outer mgc:  0.8515329 
#> iter: 1  value: 1652.198 mgc: 0.1086317 ustep: 1 
#> iter: 2  value: 1652.198 mgc: 1.932167e-05 ustep: 1 
#> iter: 3  mgc: 1.661972e-12 
#> iter: 1  value: 1652.193 mgc: 0.05846386 ustep: 1 
#> iter: 2  value: 1652.193 mgc: 9.578518e-06 ustep: 1 
#> iter: 3  mgc: 4.451049e-13 
#> iter: 1  mgc: 4.451049e-13 
#> outer mgc:  0.2112131 
#> iter: 1  value: 1652.193 mgc: 0.06432241 ustep: 1 
#> iter: 2  value: 1652.193 mgc: 5.038647e-06 ustep: 1 
#> iter: 3  mgc: 1.882006e-13 
#> iter: 1  mgc: 1.882006e-13 
#> outer mgc:  0.2769597 
#> iter: 1  value: 1652.189 mgc: 0.08643684 ustep: 1 
#> iter: 2  value: 1652.189 mgc: 9.653457e-06 ustep: 1 
#> iter: 3  mgc: 4.075472e-13 
#> iter: 1  mgc: 4.075472e-13 
#> outer mgc:  0.3603231 
#> iter: 1  value: 1652.185 mgc: 0.07390164 ustep: 1 
#> iter: 2  value: 1652.185 mgc: 7.850679e-06 ustep: 1 
#> iter: 3  mgc: 2.171764e-13 
#> iter: 1  mgc: 2.171764e-13 
#> outer mgc:  0.1657288 
#> iter: 1  value: 1652.187 mgc: 0.1042255 ustep: 1 
#> iter: 2  value: 1652.187 mgc: 1.563459e-05 ustep: 1 
#> iter: 3  mgc: 1.549735e-12 
#> iter: 1  mgc: 1.549735e-12 
#> outer mgc:  0.6554731 
#> iter: 1  value: 1652.186 mgc: 0.03394853 ustep: 1 
#> iter: 2  value: 1652.186 mgc: 5.729469e-06 ustep: 1 
#> iter: 3  mgc: 1.602509e-12 
#> iter: 1  mgc: 1.602509e-12 
#> outer mgc:  0.4005079 
#> iter: 1  value: 1652.187 mgc: 0.03633926 ustep: 1 
#> iter: 2  value: 1652.187 mgc: 4.279425e-06 ustep: 1 
#> iter: 3  mgc: 2.602754e-13 
#> iter: 1  mgc: 2.602754e-13 
#> outer mgc:  0.2140237 
#> iter: 1  value: 1652.183 mgc: 0.02927618 ustep: 1 
#> iter: 2  value: 1652.183 mgc: 2.132542e-06 ustep: 1 
#> iter: 3  mgc: 2.055986e-13 
#> iter: 1  mgc: 2.055986e-13 
#> outer mgc:  0.09527991 
#> iter: 1  value: 1652.186 mgc: 0.0367471 ustep: 1 
#> iter: 2  value: 1652.186 mgc: 2.236627e-06 ustep: 1 
#> iter: 3  mgc: 1.33748e-13 
#> iter: 1  mgc: 1.33748e-13 
#> outer mgc:  0.1941599 
#> iter: 1  value: 1652.19 mgc: 0.01826039 ustep: 1 
#> iter: 2  value: 1652.19 mgc: 4.532311e-07 ustep: 1 
#> iter: 3  mgc: 1.250873e-13 
#> iter: 1  value: 1652.198 mgc: 0.03022021 ustep: 1 
#> iter: 2  value: 1652.198 mgc: 2.539112e-06 ustep: 1 
#> iter: 3  mgc: 1.574493e-13 
#> iter: 1  mgc: 1.574493e-13 
#> outer mgc:  0.4905098 
#> iter: 1  value: 1652.202 mgc: 0.07140331 ustep: 1 
#> iter: 2  value: 1652.202 mgc: 6.384532e-06 ustep: 1 
#> iter: 3  mgc: 2.699847e-13 
#> iter: 1  mgc: 2.699847e-13 
#> outer mgc:  0.1434465 
#> iter: 1  value: 1652.196 mgc: 0.04153952 ustep: 1 
#> iter: 2  value: 1652.196 mgc: 2.185733e-06 ustep: 1 
#> iter: 3  mgc: 1.228518e-13 
#> iter: 1  mgc: 1.228518e-13 
#> outer mgc:  0.1869753 
#> iter: 1  value: 1652.193 mgc: 0.04172082 ustep: 1 
#> iter: 2  value: 1652.193 mgc: 2.818107e-06 ustep: 1 
#> iter: 3  mgc: 2.363116e-13 
#> iter: 1  mgc: 2.363116e-13 
#> outer mgc:  0.5464544 
#> iter: 1  value: 1652.193 mgc: 0.07935972 ustep: 1 
#> iter: 2  value: 1652.193 mgc: 2.282359e-06 ustep: 1 
#> iter: 3  mgc: 1.867733e-13 
#> iter: 1  mgc: 1.867733e-13 
#> outer mgc:  0.06635358 
#> iter: 1  value: 1652.189 mgc: 0.03206009 ustep: 1 
#> iter: 2  value: 1652.189 mgc: 6.127467e-07 ustep: 1 
#> iter: 3  mgc: 1.513203e-13 
#> iter: 1  value: 1652.185 mgc: 0.04860156 ustep: 1 
#> iter: 2  value: 1652.185 mgc: 2.410471e-06 ustep: 1 
#> iter: 3  mgc: 2.196759e-13 
#> iter: 1  mgc: 2.196759e-13 
#> outer mgc:  0.1849262 
#> iter: 1  value: 1652.183 mgc: 0.0424482 ustep: 1 
#> iter: 2  value: 1652.183 mgc: 4.025106e-06 ustep: 1 
#> iter: 3  mgc: 2.829043e-13 
#> iter: 1  value: 1652.182 mgc: 0.07289067 ustep: 1 
#> iter: 2  value: 1652.182 mgc: 1.8462e-05 ustep: 1 
#> iter: 3  mgc: 5.850619e-12 
#> iter: 1  mgc: 5.850619e-12 
#> outer mgc:  0.5107202 
#> iter: 1  value: 1652.193 mgc: 0.04935571 ustep: 1 
#> iter: 2  value: 1652.193 mgc: 2.504222e-05 ustep: 1 
#> iter: 3  mgc: 1.96182e-11 
#> iter: 1  mgc: 1.96182e-11 
#> outer mgc:  0.1987704 
#> iter: 1  value: 1652.171 mgc: 0.2839283 ustep: 1 
#> iter: 2  value: 1652.171 mgc: 6.941948e-05 ustep: 1 
#> iter: 3  mgc: 6.465401e-11 
#> iter: 1  mgc: 6.465401e-11 
#> outer mgc:  0.7395688 
#> iter: 1  value: 1652.215 mgc: 0.2280627 ustep: 1 
#> iter: 2  value: 1652.215 mgc: 9.034461e-05 ustep: 1 
#> iter: 3  mgc: 6.308521e-11 
#> iter: 1  value: 1652.183 mgc: 0.09882282 ustep: 1 
#> iter: 2  value: 1652.183 mgc: 1.070867e-05 ustep: 1 
#> iter: 3  mgc: 6.865039e-13 
#> iter: 1  mgc: 6.865039e-13 
#> outer mgc:  0.1553915 
#> iter: 1  value: 1652.203 mgc: 0.1077465 ustep: 1 
#> iter: 2  value: 1652.203 mgc: 5.136313e-05 ustep: 1 
#> iter: 3  mgc: 1.832173e-11 
#> iter: 1  mgc: 1.832173e-11 
#> outer mgc:  0.3734948 
#> iter: 1  value: 1652.18 mgc: 0.06918065 ustep: 1 
#> iter: 2  value: 1652.18 mgc: 1.260897e-05 ustep: 1 
#> iter: 3  mgc: 5.059508e-12 
#> iter: 1  mgc: 5.059508e-12 
#> outer mgc:  0.177382 
#> iter: 1  value: 1652.175 mgc: 0.1183499 ustep: 1 
#> iter: 2  value: 1652.175 mgc: 1.958742e-05 ustep: 1 
#> iter: 3  mgc: 1.931966e-12 
#> iter: 1  mgc: 1.931966e-12 
#> outer mgc:  0.492481 
#> iter: 1  value: 1652.189 mgc: 0.06614923 ustep: 1 
#> iter: 2  value: 1652.189 mgc: 1.898584e-05 ustep: 1 
#> iter: 3  mgc: 6.706303e-12 
#> iter: 1  mgc: 6.706303e-12 
#> outer mgc:  0.2745685 
#> iter: 1  value: 1652.197 mgc: 0.1129288 ustep: 1 
#> iter: 2  value: 1652.197 mgc: 5.559535e-05 ustep: 1 
#> iter: 3  mgc: 2.892387e-11 
#> iter: 1  mgc: 2.892387e-11 
#> outer mgc:  0.3849058 
#> iter: 1  value: 1652.217 mgc: 0.2258848 ustep: 1 
#> iter: 2  value: 1652.217 mgc: 2.242091e-05 ustep: 1 
#> iter: 3  mgc: 7.908697e-12 
#> iter: 1  value: 1652.202 mgc: 0.06377008 ustep: 1 
#> iter: 2  value: 1652.202 mgc: 5.538496e-06 ustep: 1 
#> iter: 3  mgc: 1.37278e-13 
#> iter: 1  mgc: 1.37278e-13 
#> outer mgc:  0.2537189 
#> iter: 1  value: 1652.199 mgc: 0.04252176 ustep: 1 
#> iter: 2  value: 1652.199 mgc: 3.611732e-06 ustep: 1 
#> iter: 3  mgc: 2.895324e-13 
#> iter: 1  mgc: 2.895324e-13 
#> outer mgc:  0.2848004 
#> iter: 1  value: 1652.193 mgc: 0.04368312 ustep: 1 
#> iter: 2  value: 1652.193 mgc: 5.759186e-06 ustep: 1 
#> iter: 3  mgc: 2.620518e-13 
#> iter: 1  mgc: 2.620518e-13 
#> outer mgc:  0.2073611 
#> iter: 1  value: 1652.198 mgc: 0.04433854 ustep: 1 
#> iter: 2  value: 1652.198 mgc: 3.066503e-06 ustep: 1 
#> iter: 3  mgc: 2.153304e-13 
#> iter: 1  mgc: 2.153304e-13 
#> outer mgc:  0.192915 
#> iter: 1  value: 1652.195 mgc: 0.02353755 ustep: 1 
#> iter: 2  value: 1652.195 mgc: 1.46814e-06 ustep: 1 
#> iter: 3  mgc: 1.377547e-13 
#> iter: 1  mgc: 1.377547e-13 
#> outer mgc:  0.1285838 
#> iter: 1  value: 1652.191 mgc: 0.02241641 ustep: 1 
#> iter: 2  value: 1652.191 mgc: 1.424537e-06 ustep: 1 
#> iter: 3  mgc: 1.431063e-13 
#> iter: 1  mgc: 1.431063e-13 
#> outer mgc:  0.06937493 
#> iter: 1  value: 1652.192 mgc: 0.02086474 ustep: 1 
#> iter: 2  value: 1652.192 mgc: 9.136497e-07 ustep: 1 
#> iter: 3  mgc: 1.597758e-13 
#> iter: 1  mgc: 1.597758e-13 
#> outer mgc:  0.0312859 
#> iter: 1  value: 1652.189 mgc: 0.01960782 ustep: 1 
#> iter: 2  value: 1652.189 mgc: 8.817464e-07 ustep: 1 
#> iter: 3  mgc: 1.241482e-13 
#> iter: 1  value: 1652.184 mgc: 0.0349373 ustep: 1 
#> iter: 2  value: 1652.184 mgc: 2.670181e-06 ustep: 1 
#> iter: 3  mgc: 1.40412e-13 
#> iter: 1  mgc: 1.40412e-13 
#> outer mgc:  0.1031543 
#> iter: 1  value: 1652.184 mgc: 0.02972131 ustep: 1 
#> iter: 2  value: 1652.184 mgc: 1.089602e-06 ustep: 1 
#> iter: 3  mgc: 9.6816e-14 
#> iter: 1  mgc: 9.6816e-14 
#> outer mgc:  0.05715969 
#> iter: 1  value: 1652.187 mgc: 0.03888713 ustep: 1 
#> iter: 2  value: 1652.187 mgc: 3.711405e-06 ustep: 1 
#> iter: 3  mgc: 1.809154e-13 
#> iter: 1  mgc: 1.809154e-13 
#> outer mgc:  0.0710849 
#> iter: 1  value: 1652.194 mgc: 0.03261919 ustep: 1 
#> iter: 2  value: 1652.194 mgc: 4.467765e-06 ustep: 1 
#> iter: 3  mgc: 1.260909e-13 
#> iter: 1  mgc: 1.260909e-13 
#> outer mgc:  0.08456901 
#> iter: 1  value: 1652.191 mgc: 0.02228853 ustep: 1 
#> iter: 2  value: 1652.191 mgc: 1.729281e-06 ustep: 1 
#> iter: 3  mgc: 2.54745e-13 
#> iter: 1  mgc: 2.54745e-13 
#> outer mgc:  0.09019485 
#> iter: 1  value: 1652.191 mgc: 0.01998514 ustep: 1 
#> iter: 2  value: 1652.191 mgc: 5.958531e-07 ustep: 1 
#> iter: 3  mgc: 1.314006e-13 
#> iter: 1  mgc: 1.314006e-13 
#> outer mgc:  0.03631214 
#> iter: 1  value: 1652.191 mgc: 0.03453007 ustep: 1 
#> iter: 2  value: 1652.191 mgc: 3.680646e-06 ustep: 1 
#> iter: 3  mgc: 1.364667e-13 
#> iter: 1  mgc: 1.364667e-13 
#> outer mgc:  0.05853598 
#> iter: 1  value: 1652.189 mgc: 0.0128476 ustep: 1 
#> iter: 2  value: 1652.189 mgc: 6.759782e-07 ustep: 1 
#> iter: 3  mgc: 1.998477e-13 
#> iter: 1  mgc: 1.998477e-13 
#> outer mgc:  0.09498922 
#> iter: 1  value: 1652.188 mgc: 0.01429986 ustep: 1 
#> iter: 2  value: 1652.188 mgc: 5.85792e-07 ustep: 1 
#> iter: 3  mgc: 1.120637e-13 
#> iter: 1  mgc: 1.120637e-13 
#> outer mgc:  0.02118464 
#> iter: 1  value: 1652.189 mgc: 0.02918216 ustep: 1 
#> iter: 2  value: 1652.189 mgc: 6.661432e-07 ustep: 1 
#> iter: 3  mgc: 1.642417e-13 
#> iter: 1  mgc: 1.642417e-13 
#> outer mgc:  0.03806839 
#> iter: 1  value: 1652.19 mgc: 0.01785493 ustep: 1 
#> iter: 2  value: 1652.19 mgc: 6.921572e-07 ustep: 1 
#> iter: 3  mgc: 1.411761e-13 
#> iter: 1  mgc: 1.411761e-13 
#> outer mgc:  0.09966607 
#> iter: 1  value: 1652.189 mgc: 0.01425133 ustep: 1 
#> iter: 2  value: 1652.189 mgc: 2.419194e-07 ustep: 1 
#> iter: 3  mgc: 1.23748e-13 
#> iter: 1  mgc: 1.23748e-13 
#> outer mgc:  0.06793899 
#> iter: 1  value: 1652.188 mgc: 0.0153502 ustep: 1 
#> iter: 2  value: 1652.188 mgc: 3.312795e-07 ustep: 1 
#> iter: 3  mgc: 9.856538e-14 
#> iter: 1  mgc: 9.856538e-14 
#> outer mgc:  0.08748618 
#> iter: 1  value: 1652.19 mgc: 0.0191307 ustep: 1 
#> iter: 2  value: 1652.19 mgc: 2.350275e-07 ustep: 1 
#> iter: 3  mgc: 1.166486e-13 
#> iter: 1  mgc: 1.166486e-13 
#> outer mgc:  0.03204475 
#> iter: 1  value: 1652.181 mgc: 0.01304472 ustep: 1 
#> iter: 2  value: 1652.181 mgc: 2.4577e-07 ustep: 1 
#> iter: 3  mgc: 8.863661e-14 
#> iter: 1  mgc: 8.863661e-14 
#> outer mgc:  0.156434 
#> iter: 1  value: 1652.186 mgc: 0.02246862 ustep: 1 
#> iter: 2  value: 1652.186 mgc: 7.612242e-07 ustep: 1 
#> iter: 3  mgc: 1.688006e-13 
#> iter: 1  mgc: 1.688006e-13 
#> outer mgc:  0.1230355 
#> iter: 1  value: 1652.191 mgc: 0.03506742 ustep: 1 
#> iter: 2  value: 1652.191 mgc: 1.326418e-06 ustep: 1 
#> iter: 3  mgc: 1.099434e-13 
#> iter: 1  value: 1652.188 mgc: 0.01436017 ustep: 1 
#> iter: 2  value: 1652.188 mgc: 2.236788e-07 ustep: 1 
#> iter: 3  mgc: 1.958441e-13 
#> iter: 1  mgc: 1.958441e-13 
#> outer mgc:  0.01346274 
#> iter: 1  value: 1652.189 mgc: 0.03062508 ustep: 1 
#> iter: 2  value: 1652.189 mgc: 1.171972e-06 ustep: 1 
#> iter: 3  mgc: 1.458172e-13 
#> iter: 1  value: 1652.188 mgc: 0.005830903 ustep: 1 
#> iter: 2  value: 1652.188 mgc: 3.829039e-08 ustep: 1 
#> mgc: 1.136399e-13 
#> iter: 1  mgc: 1.136399e-13 
#> outer mgc:  0.03140772 
#> iter: 1  value: 1652.189 mgc: 0.003163808 ustep: 1 
#> iter: 2  value: 1652.189 mgc: 1.1449e-08 ustep: 1 
#> mgc: 1.470844e-13 
#> iter: 1  value: 1652.19 mgc: 0.003875426 ustep: 1 
#> iter: 2  value: 1652.19 mgc: 2.364882e-08 ustep: 1 
#> mgc: 1.397242e-13 
#> iter: 1  mgc: 1.397242e-13 
#> outer mgc:  0.03750259 
#> iter: 1  value: 1652.188 mgc: 0.008547299 ustep: 1 
#> iter: 2  value: 1652.188 mgc: 1.405653e-07 ustep: 1 
#> iter: 3  mgc: 1.699426e-13 
#> iter: 1  value: 1652.189 mgc: 0.003105972 ustep: 1 
#> iter: 2  value: 1652.189 mgc: 2.022984e-08 ustep: 1 
#> mgc: 1.226553e-13 
#> iter: 1  mgc: 1.226553e-13 
#> outer mgc:  0.01670406 
#> iter: 1  value: 1652.188 mgc: 0.002535161 ustep: 1 
#> iter: 2  value: 1652.188 mgc: 1.901971e-08 ustep: 1 
#> mgc: 1.335898e-13 
#> iter: 1  mgc: 1.335898e-13 
#> outer mgc:  0.01401928 
#> iter: 1  value: 1652.188 mgc: 0.003321928 ustep: 1 
#> iter: 2  value: 1652.188 mgc: 2.272657e-08 ustep: 1 
#> mgc: 1.41406e-13 
#> iter: 1  mgc: 1.41406e-13 
#> outer mgc:  0.006524178 
#> iter: 1  value: 1652.188 mgc: 0.0008807942 ustep: 1 
#> iter: 2  mgc: 2.535052e-09 
#> iter: 1  value: 1652.187 mgc: 0.002571138 ustep: 1 
#> iter: 2  value: 1652.187 mgc: 1.46172e-08 ustep: 1 
#> mgc: 1.604418e-13 
#> iter: 1  mgc: 1.604418e-13 
#> outer mgc:  0.01438176 
#> iter: 1  value: 1652.19 mgc: 0.006376819 ustep: 1 
#> iter: 2  value: 1652.19 mgc: 6.394009e-08 ustep: 1 
#> iter: 3  mgc: 1.125369e-13 
#> iter: 1  value: 1652.188 mgc: 0.002447156 ustep: 1 
#> iter: 2  mgc: 8.704249e-09 
#> iter: 1  mgc: 8.704249e-09 
#> outer mgc:  0.0112644 
#> iter: 1  value: 1652.188 mgc: 0.002309059 ustep: 1 
#> iter: 2  value: 1652.188 mgc: 1.189401e-08 ustep: 1 
#> mgc: 1.193524e-13 
#> iter: 1  mgc: 1.193524e-13 
#> outer mgc:  0.01144651 
#> iter: 1  value: 1652.189 mgc: 0.004800762 ustep: 1 
#> iter: 2  value: 1652.189 mgc: 5.849632e-08 ustep: 1 
#> mgc: 1.533671e-13 
#> iter: 1  mgc: 1.533671e-13 
#> outer mgc:  0.0161467 
#> iter: 1  value: 1652.188 mgc: 0.001487384 ustep: 1 
#> iter: 2  mgc: 5.066629e-09 
#> iter: 1  mgc: 5.066629e-09 
#> outer mgc:  0.009179983 
#> iter: 1  value: 1652.187 mgc: 0.004774273 ustep: 1 
#> iter: 2  value: 1652.187 mgc: 6.04354e-08 ustep: 1 
#> mgc: 1.175916e-13 
#> iter: 1  value: 1652.188 mgc: 0.002198816 ustep: 1 
#> iter: 2  value: 1652.188 mgc: 1.23913e-08 ustep: 1 
#> mgc: 1.748822e-13 
#> iter: 1  mgc: 1.748822e-13 
#> outer mgc:  0.01028277 
#> iter: 1  value: 1652.188 mgc: 0.001454386 ustep: 1 
#> iter: 2  mgc: 2.407349e-09 
#> iter: 1  mgc: 2.407349e-09 
#> outer mgc:  0.007262259 
#> iter: 1  value: 1652.189 mgc: 0.006650331 ustep: 1 
#> iter: 2  value: 1652.189 mgc: 4.01941e-08 ustep: 1 
#> mgc: 1.766042e-13 
#> iter: 1  value: 1652.188 mgc: 0.003029319 ustep: 1 
#> iter: 2  mgc: 7.826102e-09 
#> iter: 1  mgc: 7.826102e-09 
#> outer mgc:  0.004266858 
#> iter: 1  value: 1652.188 mgc: 0.0005607363 ustep: 1 
#> iter: 2  mgc: 2.744491e-10 
#> iter: 1  value: 1652.188 mgc: 0.0006457043 ustep: 1 
#> iter: 2  mgc: 3.635055e-10 
#> iter: 1  mgc: 3.635055e-10 
#> outer mgc:  0.001527805 
#> iter: 1  value: 1652.188 mgc: 0.00193015 ustep: 1 
#> iter: 2  mgc: 7.209364e-09 
#> iter: 1  value: 1652.188 mgc: 0.0006805881 ustep: 1 
#> iter: 2  mgc: 9.567342e-10 
#> iter: 1  mgc: 9.567342e-10 
#> outer mgc:  0.001076584 
#> iter: 1  value: 1652.188 mgc: 0.0004016501 ustep: 1 
#> iter: 2  mgc: 2.670182e-10 
#> iter: 1  value: 1652.188 mgc: 0.0007773238 ustep: 1 
#> iter: 2  mgc: 1.001261e-09 
#> iter: 1  mgc: 1.001261e-09 
#> outer mgc:  0.006174534 
#> iter: 1  value: 1652.188 mgc: 0.001884348 ustep: 1 
#> iter: 2  mgc: 6.318216e-09 
#> iter: 1  value: 1652.188 mgc: 0.0007727261 ustep: 1 
#> iter: 2  mgc: 9.956843e-10 
#> iter: 1  mgc: 9.956843e-10 
#> outer mgc:  0.001864663 
#> iter: 1  value: 1652.188 mgc: 0.0003758684 ustep: 1 
#> iter: 2  mgc: 2.100717e-10 
#> iter: 1  mgc: 2.100717e-10 
#> outer mgc:  0.003302829 
#> iter: 1  mgc: 2.100717e-10 
#> iter: 1  mgc: 2.100717e-10 
#> outer mgc:  0.003302829 
#> iter: 1  mgc: 2.100717e-10 
#> outer mgc:  0.003302829 
#> Retro Peel: 4
#> Using stats::nlminb for optimization with these control parameters:
#> iter.max: 1000
#> eval.max: 1000
#> iter: 1  value: 1622.143 mgc: 1.240068 ustep: 1 
#> iter: 2  value: 1622.142 mgc: 0.1109194 ustep: 1 
#> iter: 3  value: 1622.142 mgc: 0.0008764108 ustep: 1 
#> iter: 4  value: 1622.142 mgc: 5.583464e-08 ustep: 1 
#> mgc: 1.056351e-13 
#> iter: 1  mgc: 1.056351e-13 
#> Matching hessian patterns... Done
#> outer mgc:  26.66262 
#> iter: 1  value: 4640.057 mgc: 128.183 ustep: 0.0247161 
#> iter: 2  value: 3965.363 mgc: 48.48902 ustep: 0.08200369 
#> iter: 3  value: 3753.493 mgc: 79.62734 ustep: 0.2864342 
#> iter: 4  value: 3728.18 mgc: 21.32512 ustep: 0.535242 
#> iter: 5  value: 3726.415 mgc: 4.126907 ustep: 0.7316292 
#> iter: 6  value: 3726.367 mgc: 0.3942301 ustep: 0.8553677 
#> iter: 7  value: 3726.367 mgc: 0.02100893 ustep: 0.9248684 
#> iter: 8  value: 3726.367 mgc: 0.001385888 ustep: 0.9617046 
#> iter: 9  value: 3726.367 mgc: 4.635616e-05 ustep: 0.9806673 
#> iter: 10  value: 3726.367 mgc: 7.738801e-07 ustep: 0.9902875 
#> iter: 11  mgc: 6.447928e-09 
#> iter: 1  value: 1639.315 mgc: 11.68621 ustep: 1 
#> iter: 2  value: 1639.286 mgc: 0.348682 ustep: 1 
#> iter: 3  value: 1639.286 mgc: 0.004007792 ustep: 1 
#> iter: 4  value: 1639.286 mgc: 9.611023e-07 ustep: 1 
#> iter: 5  mgc: 1.490332e-13 
#> iter: 1  value: 1621.832 mgc: 1.21626 ustep: 1 
#> iter: 2  value: 1621.832 mgc: 0.00308581 ustep: 1 
#> iter: 3  value: 1621.832 mgc: 2.734484e-07 ustep: 1 
#> iter: 4  mgc: 1.141009e-13 
#> iter: 1  mgc: 1.141009e-13 
#> outer mgc:  27.08885 
#> iter: 1  value: 1621.495 mgc: 1.293726 ustep: 1 
#> iter: 2  value: 1621.495 mgc: 0.002839951 ustep: 1 
#> iter: 3  value: 1621.495 mgc: 3.464622e-07 ustep: 1 
#> iter: 4  mgc: 1.104812e-13 
#> iter: 1  mgc: 1.104812e-13 
#> outer mgc:  20.70388 
#> iter: 1  value: 1621.136 mgc: 0.9834165 ustep: 1 
#> iter: 2  value: 1621.136 mgc: 0.002196314 ustep: 1 
#> iter: 3  value: 1621.136 mgc: 1.244508e-07 ustep: 1 
#> iter: 4  mgc: 1.276862e-13 
#> iter: 1  mgc: 1.276862e-13 
#> outer mgc:  12.49662 
#> iter: 1  value: 1620.726 mgc: 1.042695 ustep: 1 
#> iter: 2  value: 1620.726 mgc: 0.001841844 ustep: 1 
#> iter: 3  value: 1620.726 mgc: 1.205642e-07 ustep: 1 
#> iter: 4  mgc: 9.407996e-14 
#> iter: 1  mgc: 9.407996e-14 
#> outer mgc:  13.44897 
#> iter: 1  value: 1620.159 mgc: 0.7731599 ustep: 1 
#> iter: 2  value: 1620.159 mgc: 0.00226626 ustep: 1 
#> iter: 3  value: 1620.159 mgc: 1.307487e-07 ustep: 1 
#> iter: 4  mgc: 1.16743e-13 
#> iter: 1  mgc: 1.16743e-13 
#> outer mgc:  9.962375 
#> iter: 1  value: 1619.795 mgc: 0.9426379 ustep: 1 
#> iter: 2  value: 1619.795 mgc: 0.0005083878 ustep: 1 
#> iter: 3  value: 1619.795 mgc: 1.740623e-08 ustep: 1 
#> iter: 4  mgc: 1.714052e-13 
#> iter: 1  mgc: 1.714052e-13 
#> outer mgc:  9.497003 
#> iter: 1  value: 1619.467 mgc: 0.9590822 ustep: 1 
#> iter: 2  value: 1619.467 mgc: 0.001018496 ustep: 1 
#> iter: 3  value: 1619.467 mgc: 8.829534e-08 ustep: 1 
#> iter: 4  mgc: 1.64354e-13 
#> iter: 1  mgc: 1.64354e-13 
#> outer mgc:  8.605111 
#> iter: 1  value: 1619.168 mgc: 1.349124 ustep: 1 
#> iter: 2  value: 1619.168 mgc: 0.001355162 ustep: 1 
#> iter: 3  value: 1619.168 mgc: 1.065115e-07 ustep: 1 
#> iter: 4  mgc: 1.729225e-13 
#> iter: 1  mgc: 1.729225e-13 
#> outer mgc:  10.10989 
#> iter: 1  value: 1620.38 mgc: 3.945866 ustep: 1 
#> iter: 2  value: 1620.38 mgc: 0.03012134 ustep: 1 
#> iter: 3  value: 1620.38 mgc: 4.002935e-05 ustep: 1 
#> iter: 4  mgc: 9.343552e-11 
#> iter: 1  value: 1619.461 mgc: 1.767046 ustep: 1 
#> iter: 2  value: 1619.461 mgc: 0.005391741 ustep: 1 
#> iter: 3  value: 1619.461 mgc: 1.378297e-06 ustep: 1 
#> iter: 4  mgc: 1.196314e-13 
#> iter: 1  value: 1619.098 mgc: 0.4201845 ustep: 1 
#> iter: 2  value: 1619.098 mgc: 0.0002991872 ustep: 1 
#> iter: 3  mgc: 4.14279e-09 
#> iter: 1  mgc: 4.14279e-09 
#> outer mgc:  8.290095 
#> iter: 1  value: 1619.07 mgc: 0.3519557 ustep: 1 
#> iter: 2  value: 1619.07 mgc: 0.0001878201 ustep: 1 
#> iter: 3  mgc: 1.516091e-09 
#> iter: 1  mgc: 1.516091e-09 
#> outer mgc:  8.226486 
#> iter: 1  value: 1618.988 mgc: 0.4017443 ustep: 1 
#> iter: 2  value: 1618.988 mgc: 0.0003120362 ustep: 1 
#> iter: 3  mgc: 3.887748e-09 
#> iter: 1  mgc: 3.887748e-09 
#> outer mgc:  8.123276 
#> iter: 1  value: 1618.872 mgc: 0.4116906 ustep: 1 
#> iter: 2  value: 1618.872 mgc: 0.0001635126 ustep: 1 
#> iter: 3  mgc: 1.814289e-09 
#> iter: 1  value: 1618.721 mgc: 0.2052786 ustep: 1 
#> iter: 2  value: 1618.721 mgc: 0.000131157 ustep: 1 
#> iter: 3  mgc: 4.290066e-10 
#> iter: 1  value: 1618.299 mgc: 0.7849065 ustep: 1 
#> iter: 2  value: 1618.299 mgc: 0.001899449 ustep: 1 
#> iter: 3  value: 1618.299 mgc: 8.997602e-08 ustep: 1 
#> iter: 4  mgc: 1.370319e-13 
#> iter: 1  mgc: 1.370319e-13 
#> outer mgc:  9.077422 
#> iter: 1  value: 1617.154 mgc: 2.587962 ustep: 1 
#> iter: 2  value: 1617.154 mgc: 0.005230356 ustep: 1 
#> iter: 3  value: 1617.154 mgc: 2.526968e-06 ustep: 1 
#> iter: 4  mgc: 8.147844e-13 
#> iter: 1  mgc: 8.147844e-13 
#> outer mgc:  13.09304 
#> iter: 1  value: 1618.145 mgc: 4.199589 ustep: 1 
#> iter: 2  value: 1618.144 mgc: 0.03289188 ustep: 1 
#> iter: 3  value: 1618.144 mgc: 4.001819e-05 ustep: 1 
#> iter: 4  mgc: 9.850792e-11 
#> iter: 1  value: 1617.923 mgc: 2.674038 ustep: 1 
#> iter: 2  value: 1617.923 mgc: 0.01260021 ustep: 1 
#> iter: 3  value: 1617.923 mgc: 7.41293e-06 ustep: 1 
#> iter: 4  mgc: 3.394008e-12 
#> iter: 1  value: 1617.188 mgc: 1.055792 ustep: 1 
#> iter: 2  value: 1617.188 mgc: 0.001680542 ustep: 1 
#> iter: 3  value: 1617.188 mgc: 1.494543e-07 ustep: 1 
#> iter: 4  mgc: 1.844951e-13 
#> iter: 1  mgc: 1.844951e-13 
#> outer mgc:  12.79813 
#> iter: 1  value: 1617.078 mgc: 0.6421369 ustep: 1 
#> iter: 2  value: 1617.078 mgc: 0.0006752367 ustep: 1 
#> iter: 3  value: 1617.078 mgc: 1.892587e-08 ustep: 1 
#> mgc: 3.32639e-13 
#> iter: 1  mgc: 3.32639e-13 
#> outer mgc:  5.005499 
#> iter: 1  value: 1617.014 mgc: 0.452719 ustep: 1 
#> iter: 2  value: 1617.014 mgc: 0.0005090622 ustep: 1 
#> iter: 3  mgc: 6.355419e-09 
#> iter: 1  mgc: 6.355419e-09 
#> outer mgc:  4.86102 
#> iter: 1  value: 1616.948 mgc: 0.6471741 ustep: 1 
#> iter: 2  value: 1616.948 mgc: 0.0007583549 ustep: 1 
#> iter: 3  value: 1616.948 mgc: 2.524393e-08 ustep: 1 
#> mgc: 1.237803e-13 
#> iter: 1  mgc: 1.237803e-13 
#> outer mgc:  6.24787 
#> iter: 1  value: 1616.838 mgc: 0.5592112 ustep: 1 
#> iter: 2  value: 1616.838 mgc: 0.0005290289 ustep: 1 
#> iter: 3  value: 1616.838 mgc: 1.147151e-08 ustep: 1 
#> mgc: 1.461378e-13 
#> iter: 1  mgc: 1.461378e-13 
#> outer mgc:  5.803196 
#> iter: 1  value: 1616.777 mgc: 0.283093 ustep: 1 
#> iter: 2  value: 1616.777 mgc: 0.0001732431 ustep: 1 
#> iter: 3  mgc: 7.118903e-10 
#> iter: 1  value: 1616.728 mgc: 0.2072493 ustep: 1 
#> iter: 2  value: 1616.728 mgc: 5.773689e-05 ustep: 1 
#> iter: 3  mgc: 2.247094e-10 
#> iter: 1  mgc: 2.247094e-10 
#> outer mgc:  4.575768 
#> iter: 1  value: 1616.558 mgc: 0.4865542 ustep: 1 
#> iter: 2  value: 1616.558 mgc: 0.0008008138 ustep: 1 
#> iter: 3  value: 1616.558 mgc: 1.456382e-08 ustep: 1 
#> mgc: 1.197983e-13 
#> iter: 1  value: 1616.411 mgc: 0.4898763 ustep: 1 
#> iter: 2  value: 1616.411 mgc: 0.0002796812 ustep: 1 
#> iter: 3  mgc: 5.957279e-09 
#> iter: 1  value: 1616.319 mgc: 1.146441 ustep: 1 
#> iter: 2  value: 1616.319 mgc: 0.001521317 ustep: 1 
#> iter: 3  value: 1616.319 mgc: 1.725361e-07 ustep: 1 
#> iter: 4  mgc: 1.065931e-13 
#> iter: 1  mgc: 1.065931e-13 
#> outer mgc:  11.14487 
#> iter: 1  value: 1616.641 mgc: 3.733982 ustep: 1 
#> iter: 2  value: 1616.64 mgc: 0.02040742 ustep: 1 
#> iter: 3  value: 1616.64 mgc: 2.199509e-05 ustep: 1 
#> iter: 4  mgc: 4.682386e-11 
#> iter: 1  value: 1617.171 mgc: 3.323334 ustep: 1 
#> iter: 2  value: 1617.171 mgc: 0.01856999 ustep: 1 
#> iter: 3  value: 1617.171 mgc: 1.661525e-05 ustep: 1 
#> iter: 4  mgc: 1.638197e-11 
#> iter: 1  value: 1616.544 mgc: 2.069915 ustep: 1 
#> iter: 2  value: 1616.544 mgc: 0.007206665 ustep: 1 
#> iter: 3  value: 1616.544 mgc: 2.232178e-06 ustep: 1 
#> iter: 4  mgc: 2.880336e-13 
#> iter: 1  value: 1616.231 mgc: 0.6018736 ustep: 1 
#> iter: 2  value: 1616.231 mgc: 0.0005968943 ustep: 1 
#> iter: 3  value: 1616.231 mgc: 1.493776e-08 ustep: 1 
#> mgc: 1.084987e-13 
#> iter: 1  mgc: 1.084987e-13 
#> outer mgc:  5.671484 
#> iter: 1  value: 1616.201 mgc: 0.2987458 ustep: 1 
#> iter: 2  value: 1616.201 mgc: 8.856553e-05 ustep: 1 
#> iter: 3  mgc: 2.540897e-10 
#> iter: 1  mgc: 2.540897e-10 
#> outer mgc:  5.793292 
#> iter: 1  value: 1616.119 mgc: 0.4096592 ustep: 1 
#> iter: 2  value: 1616.119 mgc: 0.0002716785 ustep: 1 
#> iter: 3  mgc: 3.358346e-09 
#> iter: 1  mgc: 3.358346e-09 
#> outer mgc:  7.513054 
#> iter: 1  value: 1616.019 mgc: 0.3908799 ustep: 1 
#> iter: 2  value: 1616.019 mgc: 0.0004283896 ustep: 1 
#> iter: 3  mgc: 3.757082e-09 
#> iter: 1  mgc: 3.757082e-09 
#> outer mgc:  4.947338 
#> iter: 1  value: 1615.821 mgc: 0.5025112 ustep: 1 
#> iter: 2  value: 1615.821 mgc: 0.0004452 ustep: 1 
#> iter: 3  mgc: 9.048895e-09 
#> iter: 1  mgc: 9.048895e-09 
#> outer mgc:  9.628935 
#> iter: 1  value: 1615.617 mgc: 0.4580412 ustep: 1 
#> iter: 2  value: 1615.617 mgc: 0.0007557801 ustep: 1 
#> iter: 3  mgc: 7.513577e-09 
#> iter: 1  mgc: 7.513577e-09 
#> outer mgc:  6.572817 
#> iter: 1  value: 1615.437 mgc: 0.5144713 ustep: 1 
#> iter: 2  value: 1615.437 mgc: 0.0005768245 ustep: 1 
#> iter: 3  value: 1615.437 mgc: 1.850644e-08 ustep: 1 
#> mgc: 1.158765e-13 
#> iter: 1  mgc: 1.158765e-13 
#> outer mgc:  5.042435 
#> iter: 1  value: 1615.361 mgc: 1.092009 ustep: 1 
#> iter: 2  value: 1615.361 mgc: 0.001316983 ustep: 1 
#> iter: 3  value: 1615.361 mgc: 6.631952e-08 ustep: 1 
#> iter: 4  mgc: 9.828403e-14 
#> iter: 1  mgc: 9.828403e-14 
#> outer mgc:  12.80468 
#> iter: 1  value: 1615.536 mgc: 1.920979 ustep: 1 
#> iter: 2  value: 1615.536 mgc: 0.005634246 ustep: 1 
#> iter: 3  value: 1615.536 mgc: 1.347963e-06 ustep: 1 
#> iter: 4  mgc: 1.667366e-13 
#> iter: 1  value: 1615.448 mgc: 1.505926 ustep: 1 
#> iter: 2  value: 1615.448 mgc: 0.003657931 ustep: 1 
#> iter: 3  value: 1615.448 mgc: 5.676151e-07 ustep: 1 
#> iter: 4  mgc: 1.351828e-13 
#> iter: 1  value: 1615.327 mgc: 0.9854289 ustep: 1 
#> iter: 2  value: 1615.327 mgc: 0.001614323 ustep: 1 
#> iter: 3  value: 1615.327 mgc: 1.116159e-07 ustep: 1 
#> iter: 4  mgc: 1.3424e-13 
#> iter: 1  value: 1615.284 mgc: 0.447674 ustep: 1 
#> iter: 2  value: 1615.284 mgc: 0.0003306892 ustep: 1 
#> iter: 3  mgc: 4.661588e-09 
#> iter: 1  mgc: 4.661588e-09 
#> outer mgc:  6.625753 
#> iter: 1  value: 1615.257 mgc: 0.3743288 ustep: 1 
#> iter: 2  value: 1615.257 mgc: 0.0002325613 ustep: 1 
#> iter: 3  mgc: 1.515766e-09 
#> iter: 1  mgc: 1.515766e-09 
#> outer mgc:  7.739144 
#> iter: 1  value: 1615.209 mgc: 0.2877302 ustep: 1 
#> iter: 2  value: 1615.209 mgc: 0.0001230005 ustep: 1 
#> iter: 3  mgc: 9.788335e-10 
#> iter: 1  mgc: 9.788335e-10 
#> outer mgc:  4.627784 
#> iter: 1  value: 1615.177 mgc: 0.3549278 ustep: 1 
#> iter: 2  value: 1615.177 mgc: 0.0002927294 ustep: 1 
#> iter: 3  mgc: 2.010945e-09 
#> iter: 1  mgc: 2.010945e-09 
#> outer mgc:  4.995482 
#> iter: 1  value: 1615.09 mgc: 0.3560827 ustep: 1 
#> iter: 2  value: 1615.09 mgc: 0.0002728925 ustep: 1 
#> iter: 3  mgc: 3.188879e-10 
#> iter: 1  value: 1614.986 mgc: 0.3976703 ustep: 1 
#> iter: 2  value: 1614.986 mgc: 0.0004195994 ustep: 1 
#> iter: 3  mgc: 6.504731e-10 
#> iter: 1  value: 1614.873 mgc: 0.8341979 ustep: 1 
#> iter: 2  value: 1614.873 mgc: 0.001856529 ustep: 1 
#> iter: 3  value: 1614.873 mgc: 1.209156e-08 ustep: 1 
#> mgc: 1.503904e-13 
#> iter: 1  mgc: 1.503904e-13 
#> outer mgc:  5.508638 
#> iter: 1  value: 1614.507 mgc: 1.422921 ustep: 1 
#> iter: 2  value: 1614.507 mgc: 0.001959661 ustep: 1 
#> iter: 3  value: 1614.507 mgc: 5.797692e-08 ustep: 1 
#> mgc: 1.937585e-13 
#> iter: 1  mgc: 1.937585e-13 
#> outer mgc:  14.62512 
#> iter: 1  value: 1614.456 mgc: 2.548334 ustep: 1 
#> iter: 2  value: 1614.456 mgc: 0.005975194 ustep: 1 
#> iter: 3  value: 1614.456 mgc: 1.148728e-06 ustep: 1 
#> iter: 4  mgc: 9.279108e-14 
#> iter: 1  value: 1614.355 mgc: 1.348751 ustep: 1 
#> iter: 2  value: 1614.355 mgc: 0.002608658 ustep: 1 
#> iter: 3  value: 1614.355 mgc: 3.365886e-07 ustep: 1 
#> iter: 4  mgc: 1.382295e-13 
#> iter: 1  mgc: 1.382295e-13 
#> outer mgc:  8.673705 
#> iter: 1  value: 1615.1 mgc: 1.000164 ustep: 1 
#> iter: 2  value: 1615.1 mgc: 0.001905277 ustep: 1 
#> iter: 3  value: 1615.1 mgc: 9.962182e-08 ustep: 1 
#> iter: 4  mgc: 1.956382e-13 
#> iter: 1  value: 1614.638 mgc: 0.5477015 ustep: 1 
#> iter: 2  value: 1614.638 mgc: 0.0004015864 ustep: 1 
#> iter: 3  mgc: 4.843848e-09 
#> iter: 1  mgc: 4.843848e-09 
#> outer mgc:  2.405777 
#> iter: 1  value: 1614.066 mgc: 0.742041 ustep: 1 
#> iter: 2  value: 1614.066 mgc: 0.0004587516 ustep: 1 
#> iter: 3  mgc: 9.052872e-09 
#> iter: 1  mgc: 9.052872e-09 
#> outer mgc:  4.091925 
#> iter: 1  value: 1614.324 mgc: 1.331996 ustep: 1 
#> iter: 2  value: 1614.324 mgc: 0.003275272 ustep: 1 
#> iter: 3  value: 1614.324 mgc: 3.160259e-07 ustep: 1 
#> iter: 4  mgc: 9.784443e-14 
#> iter: 1  value: 1614.169 mgc: 1.046122 ustep: 1 
#> iter: 2  value: 1614.169 mgc: 0.001844012 ustep: 1 
#> iter: 3  value: 1614.169 mgc: 1.20924e-07 ustep: 1 
#> iter: 4  mgc: 1.193126e-13 
#> iter: 1  value: 1614.073 mgc: 0.6712784 ustep: 1 
#> iter: 2  value: 1614.073 mgc: 0.0007162332 ustep: 1 
#> iter: 3  value: 1614.073 mgc: 2.042988e-08 ustep: 1 
#> mgc: 1.148988e-13 
#> iter: 1  value: 1614.052 mgc: 0.240573 ustep: 1 
#> iter: 2  value: 1614.052 mgc: 9.145614e-05 ustep: 1 
#> iter: 3  mgc: 3.309281e-10 
#> iter: 1  mgc: 3.309281e-10 
#> outer mgc:  3.570221 
#> iter: 1  value: 1614.06 mgc: 0.1220372 ustep: 1 
#> iter: 2  value: 1614.06 mgc: 1.202516e-05 ustep: 1 
#> iter: 3  mgc: 1.111098e-11 
#> iter: 1  mgc: 1.111098e-11 
#> outer mgc:  2.496218 
#> iter: 1  value: 1614.086 mgc: 0.215128 ustep: 1 
#> iter: 2  value: 1614.086 mgc: 8.135331e-05 ustep: 1 
#> iter: 3  mgc: 2.635887e-10 
#> iter: 1  mgc: 2.635887e-10 
#> outer mgc:  4.394545 
#> iter: 1  value: 1614.158 mgc: 0.1591216 ustep: 1 
#> iter: 2  value: 1614.158 mgc: 8.651644e-05 ustep: 1 
#> iter: 3  mgc: 1.62704e-10 
#> iter: 1  mgc: 1.62704e-10 
#> outer mgc:  2.029511 
#> iter: 1  value: 1614.188 mgc: 0.1412701 ustep: 1 
#> iter: 2  value: 1614.188 mgc: 3.00549e-05 ustep: 1 
#> iter: 3  mgc: 2.928572e-11 
#> iter: 1  mgc: 2.928572e-11 
#> outer mgc:  3.406757 
#> iter: 1  value: 1613.999 mgc: 0.1456153 ustep: 1 
#> iter: 2  value: 1613.999 mgc: 0.0001176331 ustep: 1 
#> iter: 3  mgc: 5.730755e-10 
#> iter: 1  mgc: 5.730755e-10 
#> outer mgc:  1.365693 
#> iter: 1  value: 1614.211 mgc: 0.3327001 ustep: 1 
#> iter: 2  value: 1614.211 mgc: 0.0001359334 ustep: 1 
#> iter: 3  mgc: 3.150091e-10 
#> iter: 1  mgc: 3.150091e-10 
#> outer mgc:  2.673506 
#> iter: 1  value: 1614.384 mgc: 0.4985833 ustep: 1 
#> iter: 2  value: 1614.384 mgc: 0.001277996 ustep: 1 
#> iter: 3  value: 1614.384 mgc: 4.295315e-08 ustep: 1 
#> iter: 4  mgc: 1.025876e-13 
#> iter: 1  mgc: 1.025876e-13 
#> outer mgc:  3.534957 
#> iter: 1  value: 1614.325 mgc: 0.4356296 ustep: 1 
#> iter: 2  value: 1614.325 mgc: 0.0003746008 ustep: 1 
#> iter: 3  mgc: 3.531967e-09 
#> iter: 1  value: 1614.349 mgc: 0.2897222 ustep: 1 
#> iter: 2  value: 1614.349 mgc: 0.0002131594 ustep: 1 
#> iter: 3  mgc: 1.222364e-09 
#> iter: 1  value: 1614.362 mgc: 0.399994 ustep: 1 
#> iter: 2  value: 1614.362 mgc: 0.0002915966 ustep: 1 
#> iter: 3  mgc: 2.706244e-09 
#> iter: 1  value: 1614.369 mgc: 0.4478594 ustep: 1 
#> iter: 2  value: 1614.369 mgc: 0.0003187671 ustep: 1 
#> iter: 3  mgc: 4.005393e-09 
#> iter: 1  value: 1614.369 mgc: 0.1774214 ustep: 1 
#> iter: 2  value: 1614.369 mgc: 4.984546e-05 ustep: 1 
#> iter: 3  mgc: 9.754354e-11 
#> iter: 1  mgc: 9.754354e-11 
#> outer mgc:  1.115958 
#> iter: 1  value: 1614.359 mgc: 0.0427746 ustep: 1 
#> iter: 2  value: 1614.359 mgc: 1.310693e-06 ustep: 1 
#> iter: 3  mgc: 1.580518e-13 
#> iter: 1  value: 1614.341 mgc: 0.08560522 ustep: 1 
#> iter: 2  value: 1614.341 mgc: 5.258771e-06 ustep: 1 
#> iter: 3  mgc: 1.8619e-12 
#> iter: 1  mgc: 1.8619e-12 
#> outer mgc:  2.289005 
#> iter: 1  value: 1614.348 mgc: 0.06048832 ustep: 1 
#> iter: 2  value: 1614.348 mgc: 1.156491e-05 ustep: 1 
#> iter: 3  mgc: 3.306847e-12 
#> iter: 1  mgc: 3.306847e-12 
#> outer mgc:  0.8783779 
#> iter: 1  value: 1614.341 mgc: 0.09957895 ustep: 1 
#> iter: 2  value: 1614.341 mgc: 1.808043e-05 ustep: 1 
#> iter: 3  mgc: 6.674108e-12 
#> iter: 1  value: 1614.329 mgc: 0.1204114 ustep: 1 
#> iter: 2  value: 1614.329 mgc: 2.315128e-05 ustep: 1 
#> iter: 3  mgc: 6.653254e-12 
#> iter: 1  value: 1614.303 mgc: 0.3360782 ustep: 1 
#> iter: 2  value: 1614.303 mgc: 0.0001795284 ustep: 1 
#> iter: 3  mgc: 4.031525e-10 
#> iter: 1  mgc: 4.031525e-10 
#> outer mgc:  1.79433 
#> iter: 1  value: 1614.259 mgc: 0.57806 ustep: 1 
#> iter: 2  value: 1614.259 mgc: 0.0003810695 ustep: 1 
#> iter: 3  mgc: 1.358068e-09 
#> iter: 1  mgc: 1.358068e-09 
#> outer mgc:  1.087256 
#> iter: 1  value: 1613.921 mgc: 0.5302605 ustep: 1 
#> iter: 2  value: 1613.921 mgc: 0.0004455584 ustep: 1 
#> iter: 3  value: 1613.921 mgc: 1.102674e-08 ustep: 1 
#> mgc: 1.60217e-13 
#> iter: 1  mgc: 1.60217e-13 
#> outer mgc:  1.702549 
#> iter: 1  value: 1614.102 mgc: 0.7676994 ustep: 1 
#> iter: 2  value: 1614.102 mgc: 0.001801709 ustep: 1 
#> iter: 3  value: 1614.102 mgc: 2.07809e-07 ustep: 1 
#> iter: 4  mgc: 1.970951e-13 
#> iter: 1  mgc: 1.970951e-13 
#> outer mgc:  6.917593 
#> iter: 1  value: 1613.94 mgc: 1.164654 ustep: 1 
#> iter: 2  value: 1613.94 mgc: 0.0007664447 ustep: 1 
#> iter: 3  value: 1613.94 mgc: 2.029376e-08 ustep: 1 
#> iter: 4  mgc: 1.2316e-13 
#> iter: 1  mgc: 1.2316e-13 
#> outer mgc:  6.629717 
#> iter: 1  value: 1614.014 mgc: 0.6471991 ustep: 1 
#> iter: 2  value: 1614.014 mgc: 0.0008879938 ustep: 1 
#> iter: 3  value: 1614.014 mgc: 1.488217e-08 ustep: 1 
#> mgc: 1.077642e-13 
#> iter: 1  mgc: 1.077642e-13 
#> outer mgc:  3.784432 
#> iter: 1  value: 1614.253 mgc: 1.27213 ustep: 1 
#> iter: 2  value: 1614.253 mgc: 0.001006752 ustep: 1 
#> iter: 3  value: 1614.253 mgc: 2.63473e-08 ustep: 1 
#> iter: 4  mgc: 1.498298e-13 
#> iter: 1  mgc: 1.498298e-13 
#> outer mgc:  2.154671 
#> iter: 1  value: 1614.771 mgc: 1.655733 ustep: 1 
#> iter: 2  value: 1614.771 mgc: 0.00768938 ustep: 1 
#> iter: 3  value: 1614.771 mgc: 3.805211e-07 ustep: 1 
#> iter: 4  mgc: 1.511422e-13 
#> iter: 1  value: 1614.253 mgc: 0.7759821 ustep: 1 
#> iter: 2  value: 1614.253 mgc: 0.001128339 ustep: 1 
#> iter: 3  value: 1614.253 mgc: 2.625795e-08 ustep: 1 
#> mgc: 1.386754e-13 
#> iter: 1  value: 1614.24 mgc: 0.2528296 ustep: 1 
#> iter: 2  value: 1614.24 mgc: 0.0001193306 ustep: 1 
#> iter: 3  mgc: 2.896954e-10 
#> iter: 1  mgc: 2.896954e-10 
#> outer mgc:  4.191198 
#> iter: 1  value: 1614.245 mgc: 0.2264067 ustep: 1 
#> iter: 2  value: 1614.245 mgc: 8.701854e-05 ustep: 1 
#> iter: 3  mgc: 3.404281e-10 
#> iter: 1  mgc: 3.404281e-10 
#> outer mgc:  2.687523 
#> iter: 1  value: 1614.233 mgc: 0.1136186 ustep: 1 
#> iter: 2  value: 1614.233 mgc: 2.040434e-05 ustep: 1 
#> iter: 3  mgc: 6.891497e-12 
#> iter: 1  mgc: 6.891497e-12 
#> outer mgc:  1.768083 
#> iter: 1  value: 1614.226 mgc: 0.1038225 ustep: 1 
#> iter: 2  value: 1614.226 mgc: 2.837425e-05 ustep: 1 
#> iter: 3  mgc: 4.91339e-11 
#> iter: 1  mgc: 4.91339e-11 
#> outer mgc:  2.038557 
#> iter: 1  value: 1614.211 mgc: 0.07004892 ustep: 1 
#> iter: 2  value: 1614.211 mgc: 2.294697e-05 ustep: 1 
#> iter: 3  mgc: 1.290219e-11 
#> iter: 1  mgc: 1.290219e-11 
#> outer mgc:  1.488061 
#> iter: 1  value: 1614.209 mgc: 0.06237598 ustep: 1 
#> iter: 2  value: 1614.209 mgc: 6.404507e-05 ustep: 1 
#> iter: 3  mgc: 1.084616e-10 
#> iter: 1  value: 1614.207 mgc: 0.04867672 ustep: 1 
#> iter: 2  value: 1614.207 mgc: 5.004231e-05 ustep: 1 
#> iter: 3  mgc: 6.560196e-11 
#> iter: 1  value: 1614.205 mgc: 0.1125428 ustep: 1 
#> iter: 2  value: 1614.205 mgc: 0.0002663551 ustep: 1 
#> iter: 3  mgc: 1.858233e-09 
#> iter: 1  mgc: 1.858233e-09 
#> outer mgc:  1.549049 
#> iter: 1  value: 1614.254 mgc: 0.2659016 ustep: 1 
#> iter: 2  value: 1614.254 mgc: 0.0004935649 ustep: 1 
#> iter: 3  mgc: 7.347228e-09 
#> iter: 1  mgc: 7.347228e-09 
#> outer mgc:  1.312038 
#> iter: 1  value: 1614.129 mgc: 0.3596223 ustep: 1 
#> iter: 2  value: 1614.129 mgc: 0.0005338198 ustep: 1 
#> iter: 3  mgc: 8.280878e-09 
#> iter: 1  value: 1614.175 mgc: 0.35508 ustep: 1 
#> iter: 2  value: 1614.175 mgc: 0.0001791325 ustep: 1 
#> iter: 3  mgc: 2.993681e-10 
#> iter: 1  value: 1614.2 mgc: 0.3529457 ustep: 1 
#> iter: 2  value: 1614.2 mgc: 0.0002653368 ustep: 1 
#> iter: 3  mgc: 1.100403e-09 
#> iter: 1  value: 1614.235 mgc: 0.09701069 ustep: 1 
#> iter: 2  value: 1614.235 mgc: 2.055532e-05 ustep: 1 
#> iter: 3  mgc: 6.754879e-12 
#> iter: 1  mgc: 6.754879e-12 
#> outer mgc:  1.877336 
#> iter: 1  value: 1614.234 mgc: 0.0918631 ustep: 1 
#> iter: 2  value: 1614.234 mgc: 1.66418e-05 ustep: 1 
#> iter: 3  mgc: 1.203839e-11 
#> iter: 1  mgc: 1.203839e-11 
#> outer mgc:  1.040446 
#> iter: 1  value: 1614.224 mgc: 0.04001964 ustep: 1 
#> iter: 2  value: 1614.224 mgc: 1.911031e-06 ustep: 1 
#> iter: 3  mgc: 1.326952e-13 
#> iter: 1  mgc: 1.326952e-13 
#> outer mgc:  1.173544 
#> iter: 1  value: 1614.217 mgc: 0.07723138 ustep: 1 
#> iter: 2  value: 1614.217 mgc: 2.181261e-05 ustep: 1 
#> iter: 3  mgc: 2.597778e-11 
#> iter: 1  mgc: 2.597778e-11 
#> outer mgc:  0.924759 
#> iter: 1  value: 1614.205 mgc: 0.06761091 ustep: 1 
#> iter: 2  value: 1614.205 mgc: 2.379642e-05 ustep: 1 
#> iter: 3  mgc: 1.530205e-11 
#> iter: 1  mgc: 1.530205e-11 
#> outer mgc:  0.8688727 
#> iter: 1  value: 1614.174 mgc: 0.1817813 ustep: 1 
#> iter: 2  value: 1614.174 mgc: 0.0002153833 ustep: 1 
#> iter: 3  mgc: 1.218621e-09 
#> iter: 1  mgc: 1.218621e-09 
#> outer mgc:  0.9247954 
#> iter: 1  value: 1614.132 mgc: 0.2481048 ustep: 1 
#> iter: 2  value: 1614.132 mgc: 0.0002482946 ustep: 1 
#> iter: 3  mgc: 2.917804e-09 
#> iter: 1  value: 1614.167 mgc: 0.2472567 ustep: 1 
#> iter: 2  value: 1614.167 mgc: 0.0001342212 ustep: 1 
#> iter: 3  mgc: 8.719386e-10 
#> iter: 1  value: 1614.184 mgc: 0.2468795 ustep: 1 
#> iter: 2  value: 1614.184 mgc: 9.609066e-05 ustep: 1 
#> iter: 3  mgc: 3.945317e-10 
#> iter: 1  value: 1614.174 mgc: 0.04187558 ustep: 1 
#> iter: 2  value: 1614.174 mgc: 2.677084e-06 ustep: 1 
#> iter: 3  mgc: 3.209954e-13 
#> iter: 1  mgc: 3.209954e-13 
#> outer mgc:  0.566333 
#> iter: 1  value: 1614.17 mgc: 0.03149628 ustep: 1 
#> iter: 2  value: 1614.17 mgc: 1.638205e-06 ustep: 1 
#> iter: 3  mgc: 1.511664e-13 
#> iter: 1  mgc: 1.511664e-13 
#> outer mgc:  0.6311447 
#> iter: 1  value: 1614.169 mgc: 0.02869432 ustep: 1 
#> iter: 2  value: 1614.169 mgc: 1.466435e-06 ustep: 1 
#> iter: 3  mgc: 1.269146e-13 
#> iter: 1  mgc: 1.269146e-13 
#> outer mgc:  0.540462 
#> iter: 1  value: 1614.163 mgc: 0.02836622 ustep: 1 
#> iter: 2  value: 1614.163 mgc: 1.022782e-06 ustep: 1 
#> iter: 3  mgc: 1.305392e-13 
#> iter: 1  mgc: 1.305392e-13 
#> outer mgc:  0.6140303 
#> iter: 1  value: 1614.158 mgc: 0.03310483 ustep: 1 
#> iter: 2  value: 1614.158 mgc: 3.485115e-06 ustep: 1 
#> iter: 3  mgc: 5.342128e-13 
#> iter: 1  mgc: 5.342128e-13 
#> outer mgc:  0.5028427 
#> iter: 1  value: 1614.143 mgc: 0.05039197 ustep: 1 
#> iter: 2  value: 1614.143 mgc: 5.662388e-06 ustep: 1 
#> iter: 3  mgc: 2.634044e-13 
#> iter: 1  value: 1614.125 mgc: 0.05382127 ustep: 1 
#> iter: 2  value: 1614.125 mgc: 1.011534e-05 ustep: 1 
#> iter: 3  mgc: 2.723448e-12 
#> iter: 1  value: 1614.077 mgc: 0.1365855 ustep: 1 
#> iter: 2  value: 1614.077 mgc: 6.509057e-05 ustep: 1 
#> iter: 3  mgc: 1.119508e-10 
#> iter: 1  mgc: 1.119508e-10 
#> outer mgc:  0.6323809 
#> iter: 1  value: 1614.061 mgc: 0.1180383 ustep: 1 
#> iter: 2  value: 1614.061 mgc: 0.0001208917 ustep: 1 
#> iter: 3  mgc: 3.359823e-10 
#> iter: 1  mgc: 3.359823e-10 
#> outer mgc:  0.5687066 
#> iter: 1  value: 1614.083 mgc: 0.1184759 ustep: 1 
#> iter: 2  value: 1614.083 mgc: 0.0002658636 ustep: 1 
#> iter: 3  mgc: 1.861641e-09 
#> iter: 1  mgc: 1.861641e-09 
#> outer mgc:  1.443301 
#> iter: 1  value: 1614.08 mgc: 0.1601523 ustep: 1 
#> iter: 2  value: 1614.08 mgc: 5.7237e-05 ustep: 1 
#> iter: 3  mgc: 4.148243e-11 
#> iter: 1  mgc: 4.148243e-11 
#> outer mgc:  0.5146168 
#> iter: 1  value: 1614.11 mgc: 0.1328292 ustep: 1 
#> iter: 2  value: 1614.11 mgc: 0.0001370566 ustep: 1 
#> iter: 3  mgc: 5.279723e-10 
#> iter: 1  value: 1614.087 mgc: 0.1045586 ustep: 1 
#> iter: 2  value: 1614.087 mgc: 1.389853e-05 ustep: 1 
#> iter: 3  mgc: 4.666399e-12 
#> iter: 1  value: 1614.078 mgc: 0.09860992 ustep: 1 
#> iter: 2  value: 1614.078 mgc: 1.516138e-05 ustep: 1 
#> iter: 3  mgc: 5.762253e-12 
#> iter: 1  value: 1614.074 mgc: 0.1060057 ustep: 1 
#> iter: 2  value: 1614.074 mgc: 1.793234e-05 ustep: 1 
#> iter: 3  mgc: 1.379419e-11 
#> iter: 1  value: 1614.078 mgc: 0.01932861 ustep: 1 
#> iter: 2  value: 1614.078 mgc: 5.97368e-07 ustep: 1 
#> iter: 3  mgc: 2.07165e-13 
#> iter: 1  mgc: 2.07165e-13 
#> outer mgc:  0.3979248 
#> iter: 1  value: 1614.079 mgc: 0.02368113 ustep: 1 
#> iter: 2  value: 1614.079 mgc: 8.571763e-07 ustep: 1 
#> iter: 3  mgc: 1.521263e-13 
#> iter: 1  mgc: 1.521263e-13 
#> outer mgc:  0.375674 
#> iter: 1  value: 1614.078 mgc: 0.01517029 ustep: 1 
#> iter: 2  value: 1614.078 mgc: 3.847037e-07 ustep: 1 
#> iter: 3  mgc: 1.407633e-13 
#> iter: 1  mgc: 1.407633e-13 
#> outer mgc:  0.3646386 
#> iter: 1  value: 1614.079 mgc: 0.02313274 ustep: 1 
#> iter: 2  value: 1614.079 mgc: 9.460726e-07 ustep: 1 
#> iter: 3  mgc: 1.206525e-13 
#> iter: 1  mgc: 1.206525e-13 
#> outer mgc:  0.3280124 
#> iter: 1  value: 1614.079 mgc: 0.01211217 ustep: 1 
#> iter: 2  value: 1614.079 mgc: 2.065279e-07 ustep: 1 
#> iter: 3  mgc: 1.645028e-13 
#> iter: 1  value: 1614.08 mgc: 0.007829627 ustep: 1 
#> iter: 2  value: 1614.08 mgc: 4.896886e-07 ustep: 1 
#> iter: 3  mgc: 1.22801e-13 
#> iter: 1  value: 1614.085 mgc: 0.02591391 ustep: 1 
#> iter: 2  value: 1614.085 mgc: 5.359622e-06 ustep: 1 
#> iter: 3  mgc: 8.663344e-13 
#> iter: 1  value: 1614.094 mgc: 0.04429383 ustep: 1 
#> iter: 2  value: 1614.094 mgc: 1.563666e-05 ustep: 1 
#> iter: 3  mgc: 7.451093e-12 
#> iter: 1  mgc: 7.451093e-12 
#> outer mgc:  0.4460225 
#> iter: 1  value: 1614.072 mgc: 0.1460398 ustep: 1 
#> iter: 2  value: 1614.072 mgc: 2.658683e-05 ustep: 1 
#> iter: 3  mgc: 1.912723e-11 
#> iter: 1  mgc: 1.912723e-11 
#> outer mgc:  0.3754897 
#> iter: 1  value: 1614.065 mgc: 0.07385256 ustep: 1 
#> iter: 2  value: 1614.065 mgc: 4.185132e-05 ustep: 1 
#> iter: 3  mgc: 4.513081e-11 
#> iter: 1  mgc: 4.513081e-11 
#> outer mgc:  0.6767156 
#> iter: 1  value: 1614.054 mgc: 0.04890178 ustep: 1 
#> iter: 2  value: 1614.054 mgc: 5.218084e-06 ustep: 1 
#> iter: 3  mgc: 3.218922e-13 
#> iter: 1  mgc: 3.218922e-13 
#> outer mgc:  0.2545858 
#> iter: 1  value: 1614.044 mgc: 0.1635897 ustep: 1 
#> iter: 2  value: 1614.044 mgc: 2.908891e-05 ustep: 1 
#> iter: 3  mgc: 5.166145e-12 
#> iter: 1  mgc: 5.166145e-12 
#> outer mgc:  0.8007753 
#> iter: 1  value: 1614.046 mgc: 0.07456245 ustep: 1 
#> iter: 2  value: 1614.046 mgc: 9.951172e-06 ustep: 1 
#> iter: 3  mgc: 2.42364e-12 
#> iter: 1  mgc: 2.42364e-12 
#> outer mgc:  0.6559651 
#> iter: 1  value: 1614.054 mgc: 0.02228304 ustep: 1 
#> iter: 2  value: 1614.054 mgc: 8.887523e-07 ustep: 1 
#> iter: 3  mgc: 1.032527e-13 
#> iter: 1  mgc: 1.032527e-13 
#> outer mgc:  0.1716163 
#> iter: 1  value: 1614.049 mgc: 0.06512011 ustep: 1 
#> iter: 2  value: 1614.049 mgc: 1.321107e-05 ustep: 1 
#> iter: 3  mgc: 4.674634e-12 
#> iter: 1  value: 1614.052 mgc: 0.05510411 ustep: 1 
#> iter: 2  value: 1614.052 mgc: 1.028325e-05 ustep: 1 
#> iter: 3  mgc: 2.792769e-12 
#> iter: 1  value: 1614.054 mgc: 0.05804022 ustep: 1 
#> iter: 2  value: 1614.054 mgc: 8.771562e-06 ustep: 1 
#> iter: 3  mgc: 1.989515e-12 
#> iter: 1  value: 1614.055 mgc: 0.04858152 ustep: 1 
#> iter: 2  value: 1614.055 mgc: 4.854663e-06 ustep: 1 
#> iter: 3  mgc: 4.076947e-13 
#> iter: 1  mgc: 4.076947e-13 
#> outer mgc:  0.241449 
#> iter: 1  value: 1614.055 mgc: 0.03010831 ustep: 1 
#> iter: 2  value: 1614.055 mgc: 1.325427e-06 ustep: 1 
#> iter: 3  mgc: 1.243589e-13 
#> iter: 1  mgc: 1.243589e-13 
#> outer mgc:  0.395639 
#> iter: 1  value: 1614.057 mgc: 0.01310436 ustep: 1 
#> iter: 2  value: 1614.057 mgc: 3.120086e-07 ustep: 1 
#> iter: 3  mgc: 1.585373e-13 
#> iter: 1  mgc: 1.585373e-13 
#> outer mgc:  0.1143429 
#> iter: 1  value: 1614.058 mgc: 0.01542697 ustep: 1 
#> iter: 2  value: 1614.058 mgc: 2.911443e-07 ustep: 1 
#> iter: 3  mgc: 1.547745e-13 
#> iter: 1  value: 1614.058 mgc: 0.006147306 ustep: 1 
#> iter: 2  value: 1614.058 mgc: 7.948492e-08 ustep: 1 
#> iter: 3  mgc: 1.177462e-13 
#> iter: 1  value: 1614.06 mgc: 0.0307316 ustep: 1 
#> iter: 2  value: 1614.06 mgc: 1.986672e-06 ustep: 1 
#> iter: 3  mgc: 1.111748e-13 
#> iter: 1  mgc: 1.111748e-13 
#> outer mgc:  0.1212159 
#> iter: 1  value: 1614.071 mgc: 0.03661984 ustep: 1 
#> iter: 2  value: 1614.071 mgc: 1.832339e-06 ustep: 1 
#> iter: 3  mgc: 1.22138e-13 
#> iter: 1  mgc: 1.22138e-13 
#> outer mgc:  0.08884075 
#> iter: 1  value: 1614.066 mgc: 0.0279514 ustep: 1 
#> iter: 2  value: 1614.066 mgc: 1.154431e-06 ustep: 1 
#> iter: 3  mgc: 1.312649e-13 
#> iter: 1  mgc: 1.312649e-13 
#> outer mgc:  0.1283212 
#> iter: 1  value: 1614.055 mgc: 0.02821218 ustep: 1 
#> iter: 2  value: 1614.055 mgc: 5.02551e-06 ustep: 1 
#> iter: 3  mgc: 6.529542e-13 
#> iter: 1  mgc: 6.529542e-13 
#> outer mgc:  0.2514076 
#> iter: 1  value: 1614.064 mgc: 0.04145233 ustep: 1 
#> iter: 2  value: 1614.064 mgc: 1.770003e-06 ustep: 1 
#> iter: 3  mgc: 1.227866e-13 
#> iter: 1  mgc: 1.227866e-13 
#> outer mgc:  0.262389 
#> iter: 1  value: 1614.064 mgc: 0.03575335 ustep: 1 
#> iter: 2  value: 1614.064 mgc: 3.612264e-06 ustep: 1 
#> iter: 3  mgc: 4.70624e-13 
#> iter: 1  mgc: 4.70624e-13 
#> outer mgc:  0.5442447 
#> iter: 1  value: 1614.065 mgc: 0.02385087 ustep: 1 
#> iter: 2  value: 1614.065 mgc: 1.467137e-06 ustep: 1 
#> iter: 3  mgc: 1.853962e-13 
#> iter: 1  value: 1614.067 mgc: 0.01535381 ustep: 1 
#> iter: 2  value: 1614.067 mgc: 8.625905e-07 ustep: 1 
#> iter: 3  mgc: 1.423937e-13 
#> iter: 1  value: 1614.076 mgc: 0.06058918 ustep: 1 
#> iter: 2  value: 1614.076 mgc: 1.077889e-05 ustep: 1 
#> iter: 3  mgc: 1.838305e-12 
#> iter: 1  mgc: 1.838305e-12 
#> outer mgc:  0.3563372 
#> iter: 1  value: 1614.08 mgc: 0.1734499 ustep: 1 
#> iter: 2  value: 1614.08 mgc: 6.75956e-05 ustep: 1 
#> iter: 3  mgc: 1.8716e-11 
#> iter: 1  mgc: 1.8716e-11 
#> outer mgc:  0.5047855 
#> iter: 1  value: 1614.096 mgc: 0.187817 ustep: 1 
#> iter: 2  value: 1614.096 mgc: 6.263481e-05 ustep: 1 
#> iter: 3  mgc: 3.288869e-11 
#> iter: 1  value: 1614.086 mgc: 0.06409769 ustep: 1 
#> iter: 2  value: 1614.086 mgc: 1.205407e-05 ustep: 1 
#> iter: 3  mgc: 1.198761e-12 
#> iter: 1  mgc: 1.198761e-12 
#> outer mgc:  0.3806277 
#> iter: 1  value: 1614.097 mgc: 0.1364811 ustep: 1 
#> iter: 2  value: 1614.097 mgc: 5.768815e-05 ustep: 1 
#> iter: 3  mgc: 1.308614e-11 
#> iter: 1  mgc: 1.308614e-11 
#> outer mgc:  0.7369581 
#> iter: 1  value: 1614.099 mgc: 0.1881159 ustep: 1 
#> iter: 2  value: 1614.099 mgc: 5.966882e-05 ustep: 1 
#> iter: 3  mgc: 1.272793e-11 
#> iter: 1  mgc: 1.272793e-11 
#> outer mgc:  0.5913147 
#> iter: 1  value: 1614.089 mgc: 0.1178623 ustep: 1 
#> iter: 2  value: 1614.089 mgc: 2.999463e-05 ustep: 1 
#> iter: 3  mgc: 9.072779e-12 
#> iter: 1  mgc: 9.072779e-12 
#> outer mgc:  1.511201 
#> iter: 1  value: 1614.11 mgc: 0.2210938 ustep: 1 
#> iter: 2  value: 1614.11 mgc: 0.0002053176 ustep: 1 
#> iter: 3  mgc: 1.71165e-10 
#> iter: 1  mgc: 1.71165e-10 
#> outer mgc:  0.9783013 
#> iter: 1  value: 1614.022 mgc: 0.4936438 ustep: 1 
#> iter: 2  value: 1614.022 mgc: 0.000287019 ustep: 1 
#> iter: 3  mgc: 4.630255e-10 
#> iter: 1  value: 1614.08 mgc: 0.2280535 ustep: 1 
#> iter: 2  value: 1614.08 mgc: 3.970482e-05 ustep: 1 
#> iter: 3  mgc: 2.257388e-11 
#> iter: 1  value: 1614.098 mgc: 0.1728734 ustep: 1 
#> iter: 2  value: 1614.098 mgc: 2.050868e-05 ustep: 1 
#> iter: 3  mgc: 8.802289e-12 
#> iter: 1  value: 1614.108 mgc: 0.1426416 ustep: 1 
#> iter: 2  value: 1614.108 mgc: 1.423436e-05 ustep: 1 
#> iter: 3  mgc: 9.147153e-12 
#> iter: 1  value: 1614.113 mgc: 0.1273814 ustep: 1 
#> iter: 2  value: 1614.113 mgc: 1.52009e-05 ustep: 1 
#> iter: 3  mgc: 1.20797e-11 
#> iter: 1  value: 1614.116 mgc: 0.1209428 ustep: 1 
#> iter: 2  value: 1614.116 mgc: 1.611352e-05 ustep: 1 
#> iter: 3  mgc: 1.345488e-11 
#> iter: 1  value: 1614.113 mgc: 0.06818923 ustep: 1 
#> iter: 2  value: 1614.113 mgc: 5.453331e-06 ustep: 1 
#> iter: 3  mgc: 1.552989e-12 
#> iter: 1  mgc: 1.552989e-12 
#> outer mgc:  0.8967773 
#> iter: 1  value: 1614.111 mgc: 0.04306502 ustep: 1 
#> iter: 2  value: 1614.111 mgc: 2.86647e-06 ustep: 1 
#> iter: 3  mgc: 3.175405e-13 
#> iter: 1  mgc: 3.175405e-13 
#> outer mgc:  0.2630561 
#> iter: 1  value: 1614.109 mgc: 0.04367731 ustep: 1 
#> iter: 2  value: 1614.109 mgc: 3.309279e-06 ustep: 1 
#> iter: 3  mgc: 2.694166e-13 
#> iter: 1  mgc: 2.694166e-13 
#> outer mgc:  0.5229502 
#> iter: 1  value: 1614.109 mgc: 0.0244025 ustep: 1 
#> iter: 2  value: 1614.109 mgc: 8.495555e-07 ustep: 1 
#> iter: 3  mgc: 1.8802e-13 
#> iter: 1  mgc: 1.8802e-13 
#> outer mgc:  0.204734 
#> iter: 1  value: 1614.107 mgc: 0.009456172 ustep: 1 
#> iter: 2  value: 1614.107 mgc: 2.508943e-07 ustep: 1 
#> iter: 3  mgc: 2.001904e-13 
#> iter: 1  value: 1614.105 mgc: 0.006837211 ustep: 1 
#> iter: 2  value: 1614.105 mgc: 7.324922e-08 ustep: 1 
#> mgc: 1.709898e-13 
#> iter: 1  value: 1614.096 mgc: 0.03262768 ustep: 1 
#> iter: 2  value: 1614.096 mgc: 1.667715e-06 ustep: 1 
#> iter: 3  mgc: 9.762995e-14 
#> iter: 1  value: 1614.073 mgc: 0.08839202 ustep: 1 
#> iter: 2  value: 1614.073 mgc: 1.223473e-05 ustep: 1 
#> iter: 3  mgc: 7.618088e-13 
#> iter: 1  mgc: 7.618088e-13 
#> outer mgc:  0.3131471 
#> iter: 1  value: 1614.052 mgc: 0.08820583 ustep: 1 
#> iter: 2  value: 1614.052 mgc: 1.908976e-05 ustep: 1 
#> iter: 3  mgc: 1.140167e-12 
#> iter: 1  mgc: 1.140167e-12 
#> outer mgc:  0.2124163 
#> iter: 1  value: 1614.062 mgc: 0.1224332 ustep: 1 
#> iter: 2  value: 1614.062 mgc: 1.297222e-05 ustep: 1 
#> iter: 3  mgc: 6.604859e-13 
#> iter: 1  mgc: 6.604859e-13 
#> outer mgc:  0.316721 
#> iter: 1  value: 1614.07 mgc: 0.04094853 ustep: 1 
#> iter: 2  value: 1614.07 mgc: 3.713248e-06 ustep: 1 
#> iter: 3  mgc: 2.751336e-13 
#> iter: 1  value: 1614.087 mgc: 0.07136541 ustep: 1 
#> iter: 2  value: 1614.087 mgc: 1.513224e-05 ustep: 1 
#> iter: 3  mgc: 1.144918e-12 
#> iter: 1  mgc: 1.144918e-12 
#> outer mgc:  0.9412445 
#> iter: 1  value: 1614.09 mgc: 0.1458562 ustep: 1 
#> iter: 2  value: 1614.09 mgc: 1.167439e-05 ustep: 1 
#> iter: 3  mgc: 3.490197e-12 
#> iter: 1  mgc: 3.490197e-12 
#> outer mgc:  0.4634587 
#> iter: 1  value: 1614.078 mgc: 0.2049858 ustep: 1 
#> iter: 2  value: 1614.078 mgc: 3.320564e-05 ustep: 1 
#> iter: 3  mgc: 1.173793e-11 
#> iter: 1  mgc: 1.173793e-11 
#> outer mgc:  1.101056 
#> iter: 1  value: 1614.001 mgc: 0.1847122 ustep: 1 
#> iter: 2  value: 1614.001 mgc: 8.376627e-05 ustep: 1 
#> iter: 3  mgc: 2.745134e-10 
#> iter: 1  value: 1614.047 mgc: 0.1087677 ustep: 1 
#> iter: 2  value: 1614.047 mgc: 2.96317e-05 ustep: 1 
#> iter: 3  mgc: 2.324414e-11 
#> iter: 1  mgc: 2.324414e-11 
#> outer mgc:  0.4119623 
#> iter: 1  value: 1614.09 mgc: 0.2003477 ustep: 1 
#> iter: 2  value: 1614.09 mgc: 5.815092e-05 ustep: 1 
#> iter: 3  mgc: 6.628989e-11 
#> iter: 1  value: 1614.06 mgc: 0.1016337 ustep: 1 
#> iter: 2  value: 1614.06 mgc: 2.022691e-05 ustep: 1 
#> iter: 3  mgc: 1.003742e-11 
#> iter: 1  mgc: 1.003742e-11 
#> outer mgc:  1.329494 
#> iter: 1  value: 1614.049 mgc: 0.06683877 ustep: 1 
#> iter: 2  value: 1614.049 mgc: 6.47854e-06 ustep: 1 
#> iter: 3  mgc: 2.468868e-12 
#> iter: 1  mgc: 2.468868e-12 
#> outer mgc:  0.3602504 
#> iter: 1  value: 1614.045 mgc: 0.06811031 ustep: 1 
#> iter: 2  value: 1614.045 mgc: 5.085904e-06 ustep: 1 
#> iter: 3  mgc: 3.065979e-13 
#> iter: 1  mgc: 3.065979e-13 
#> outer mgc:  0.1256921 
#> iter: 1  value: 1614.051 mgc: 0.05142512 ustep: 1 
#> iter: 2  value: 1614.051 mgc: 1.225595e-06 ustep: 1 
#> iter: 3  mgc: 1.287353e-13 
#> iter: 1  value: 1614.056 mgc: 0.04017754 ustep: 1 
#> iter: 2  value: 1614.056 mgc: 1.447275e-06 ustep: 1 
#> iter: 3  mgc: 1.091383e-13 
#> iter: 1  mgc: 1.091383e-13 
#> outer mgc:  0.1601482 
#> iter: 1  value: 1614.061 mgc: 0.04558458 ustep: 1 
#> iter: 2  value: 1614.061 mgc: 1.457994e-06 ustep: 1 
#> iter: 3  mgc: 1.261775e-13 
#> iter: 1  value: 1614.067 mgc: 0.06458891 ustep: 1 
#> iter: 2  value: 1614.067 mgc: 2.977796e-06 ustep: 1 
#> iter: 3  mgc: 1.895549e-13 
#> iter: 1  mgc: 1.895549e-13 
#> outer mgc:  0.8056979 
#> iter: 1  value: 1614.069 mgc: 0.1123197 ustep: 1 
#> iter: 2  value: 1614.069 mgc: 1.603164e-05 ustep: 1 
#> iter: 3  mgc: 1.038184e-11 
#> iter: 1  mgc: 1.038184e-11 
#> outer mgc:  0.5628634 
#> iter: 1  value: 1614.055 mgc: 0.162702 ustep: 1 
#> iter: 2  value: 1614.055 mgc: 7.392771e-05 ustep: 1 
#> iter: 3  mgc: 3.114461e-11 
#> iter: 1  value: 1614.061 mgc: 0.08974145 ustep: 1 
#> iter: 2  value: 1614.061 mgc: 1.778146e-05 ustep: 1 
#> iter: 3  mgc: 1.024459e-12 
#> iter: 1  mgc: 1.024459e-12 
#> outer mgc:  0.5432756 
#> iter: 1  value: 1614.045 mgc: 0.08187977 ustep: 1 
#> iter: 2  value: 1614.045 mgc: 1.158038e-05 ustep: 1 
#> iter: 3  mgc: 3.532753e-12 
#> iter: 1  mgc: 3.532753e-12 
#> outer mgc:  0.5661408 
#> iter: 1  value: 1614.033 mgc: 0.1058041 ustep: 1 
#> iter: 2  value: 1614.033 mgc: 2.542567e-05 ustep: 1 
#> iter: 3  mgc: 1.469595e-11 
#> iter: 1  mgc: 1.469595e-11 
#> outer mgc:  0.8697601 
#> iter: 1  value: 1614.032 mgc: 0.1402606 ustep: 1 
#> iter: 2  value: 1614.032 mgc: 3.285767e-05 ustep: 1 
#> iter: 3  mgc: 1.05139e-11 
#> iter: 1  mgc: 1.05139e-11 
#> outer mgc:  1.161477 
#> iter: 1  value: 1614.053 mgc: 0.07570412 ustep: 1 
#> iter: 2  value: 1614.053 mgc: 1.293503e-05 ustep: 1 
#> iter: 3  mgc: 1.383932e-12 
#> iter: 1  mgc: 1.383932e-12 
#> outer mgc:  0.5363327 
#> iter: 1  value: 1614.117 mgc: 0.286133 ustep: 1 
#> iter: 2  value: 1614.117 mgc: 0.0001296671 ustep: 1 
#> iter: 3  mgc: 9.5635e-10 
#> iter: 1  value: 1614.067 mgc: 0.05215914 ustep: 1 
#> iter: 2  value: 1614.067 mgc: 5.043812e-06 ustep: 1 
#> iter: 3  mgc: 1.411885e-12 
#> iter: 1  mgc: 1.411885e-12 
#> outer mgc:  0.3314373 
#> iter: 1  value: 1614.073 mgc: 0.0617403 ustep: 1 
#> iter: 2  value: 1614.073 mgc: 2.949005e-06 ustep: 1 
#> iter: 3  mgc: 2.070924e-13 
#> iter: 1  mgc: 2.070924e-13 
#> outer mgc:  0.6857884 
#> iter: 1  value: 1614.068 mgc: 0.06068772 ustep: 1 
#> iter: 2  value: 1614.068 mgc: 5.084279e-06 ustep: 1 
#> iter: 3  mgc: 6.215078e-13 
#> iter: 1  mgc: 6.215078e-13 
#> outer mgc:  0.1879357 
#> iter: 1  value: 1614.066 mgc: 0.02102097 ustep: 1 
#> iter: 2  value: 1614.066 mgc: 7.217481e-07 ustep: 1 
#> iter: 3  mgc: 1.295147e-13 
#> iter: 1  value: 1614.065 mgc: 0.02201382 ustep: 1 
#> iter: 2  value: 1614.065 mgc: 1.363459e-06 ustep: 1 
#> iter: 3  mgc: 1.045773e-13 
#> iter: 1  mgc: 1.045773e-13 
#> outer mgc:  0.2268318 
#> iter: 1  value: 1614.06 mgc: 0.03668543 ustep: 1 
#> iter: 2  value: 1614.06 mgc: 3.699832e-06 ustep: 1 
#> iter: 3  mgc: 2.01773e-13 
#> iter: 1  mgc: 2.01773e-13 
#> outer mgc:  0.1096741 
#> iter: 1  value: 1614.052 mgc: 0.04176593 ustep: 1 
#> iter: 2  value: 1614.052 mgc: 3.210391e-06 ustep: 1 
#> iter: 3  mgc: 1.854589e-13 
#> iter: 1  mgc: 1.854589e-13 
#> outer mgc:  0.2354827 
#> iter: 1  value: 1614.062 mgc: 0.03075105 ustep: 1 
#> iter: 2  value: 1614.062 mgc: 1.575443e-06 ustep: 1 
#> iter: 3  mgc: 2.009324e-13 
#> iter: 1  mgc: 2.009324e-13 
#> outer mgc:  0.3513357 
#> iter: 1  value: 1614.067 mgc: 0.03866224 ustep: 1 
#> iter: 2  value: 1614.067 mgc: 3.722307e-06 ustep: 1 
#> iter: 3  mgc: 3.45864e-13 
#> iter: 1  mgc: 3.45864e-13 
#> outer mgc:  0.07821328 
#> iter: 1  value: 1614.074 mgc: 0.04738708 ustep: 1 
#> iter: 2  value: 1614.074 mgc: 6.165297e-06 ustep: 1 
#> iter: 3  mgc: 1.797303e-13 
#> iter: 1  mgc: 1.797303e-13 
#> outer mgc:  0.3159 
#> iter: 1  value: 1614.067 mgc: 0.02896228 ustep: 1 
#> iter: 2  value: 1614.067 mgc: 1.632462e-06 ustep: 1 
#> iter: 3  mgc: 1.645155e-13 
#> iter: 1  mgc: 1.645155e-13 
#> outer mgc:  0.1913792 
#> iter: 1  value: 1614.056 mgc: 0.02543623 ustep: 1 
#> iter: 2  value: 1614.056 mgc: 1.595976e-06 ustep: 1 
#> iter: 3  mgc: 1.229525e-13 
#> iter: 1  mgc: 1.229525e-13 
#> outer mgc:  0.2173089 
#> iter: 1  value: 1614.056 mgc: 0.06376203 ustep: 1 
#> iter: 2  value: 1614.056 mgc: 5.895317e-06 ustep: 1 
#> iter: 3  mgc: 2.073476e-13 
#> iter: 1  value: 1614.056 mgc: 0.0232221 ustep: 1 
#> iter: 2  value: 1614.056 mgc: 9.347786e-07 ustep: 1 
#> iter: 3  mgc: 1.272436e-13 
#> iter: 1  mgc: 1.272436e-13 
#> outer mgc:  0.2051575 
#> iter: 1  value: 1614.055 mgc: 0.01329611 ustep: 1 
#> iter: 2  value: 1614.055 mgc: 3.939803e-07 ustep: 1 
#> iter: 3  mgc: 1.593524e-13 
#> iter: 1  mgc: 1.593524e-13 
#> outer mgc:  0.09081538 
#> iter: 1  value: 1614.059 mgc: 0.01091293 ustep: 1 
#> iter: 2  value: 1614.059 mgc: 2.854787e-07 ustep: 1 
#> iter: 3  mgc: 1.179776e-13 
#> iter: 1  mgc: 1.179776e-13 
#> outer mgc:  0.05770042 
#> iter: 1  value: 1614.063 mgc: 0.0188166 ustep: 1 
#> iter: 2  value: 1614.063 mgc: 7.740451e-07 ustep: 1 
#> iter: 3  mgc: 1.335699e-13 
#> iter: 1  mgc: 1.335699e-13 
#> outer mgc:  0.1090521 
#> iter: 1  value: 1614.064 mgc: 0.0203131 ustep: 1 
#> iter: 2  value: 1614.064 mgc: 2.872574e-07 ustep: 1 
#> iter: 3  mgc: 1.199084e-13 
#> iter: 1  mgc: 1.199084e-13 
#> outer mgc:  0.1383342 
#> iter: 1  value: 1614.058 mgc: 0.02074397 ustep: 1 
#> iter: 2  value: 1614.058 mgc: 6.179703e-07 ustep: 1 
#> iter: 3  mgc: 1.027086e-13 
#> iter: 1  mgc: 1.027086e-13 
#> outer mgc:  0.09129876 
#> iter: 1  value: 1614.047 mgc: 0.1162833 ustep: 1 
#> iter: 2  value: 1614.047 mgc: 6.189653e-06 ustep: 1 
#> iter: 3  mgc: 5.775736e-13 
#> iter: 1  value: 1614.056 mgc: 0.02090679 ustep: 1 
#> iter: 2  value: 1614.056 mgc: 1.890208e-07 ustep: 1 
#> iter: 3  mgc: 1.413953e-13 
#> iter: 1  mgc: 1.413953e-13 
#> outer mgc:  0.02016466 
#> iter: 1  value: 1614.06 mgc: 0.02812634 ustep: 1 
#> iter: 2  value: 1614.06 mgc: 2.093897e-06 ustep: 1 
#> iter: 3  mgc: 1.162721e-13 
#> iter: 1  mgc: 1.162721e-13 
#> outer mgc:  0.08076828 
#> iter: 1  value: 1614.062 mgc: 0.009623661 ustep: 1 
#> iter: 2  value: 1614.062 mgc: 2.070019e-07 ustep: 1 
#> iter: 3  mgc: 1.207659e-13 
#> iter: 1  mgc: 1.207659e-13 
#> outer mgc:  0.08569161 
#> iter: 1  value: 1614.061 mgc: 0.02829466 ustep: 1 
#> iter: 2  value: 1614.061 mgc: 6.163568e-07 ustep: 1 
#> iter: 3  mgc: 9.83871e-14 
#> iter: 1  mgc: 9.83871e-14 
#> outer mgc:  0.02097428 
#> iter: 1  value: 1614.062 mgc: 0.01354775 ustep: 1 
#> iter: 2  value: 1614.062 mgc: 2.203652e-07 ustep: 1 
#> iter: 3  mgc: 1.857134e-13 
#> iter: 1  mgc: 1.857134e-13 
#> outer mgc:  0.1016028 
#> iter: 1  value: 1614.06 mgc: 0.02591805 ustep: 1 
#> iter: 2  value: 1614.06 mgc: 2.241895e-06 ustep: 1 
#> iter: 3  mgc: 1.183105e-13 
#> iter: 1  value: 1614.061 mgc: 0.01280334 ustep: 1 
#> iter: 2  value: 1614.061 mgc: 5.470731e-07 ustep: 1 
#> iter: 3  mgc: 1.029255e-13 
#> iter: 1  mgc: 1.029255e-13 
#> outer mgc:  0.09859426 
#> iter: 1  value: 1614.056 mgc: 0.02134829 ustep: 1 
#> iter: 2  value: 1614.056 mgc: 4.40676e-07 ustep: 1 
#> iter: 3  mgc: 1.251709e-13 
#> iter: 1  value: 1614.059 mgc: 0.01060214 ustep: 1 
#> iter: 2  value: 1614.059 mgc: 8.267486e-08 ustep: 1 
#> iter: 3  mgc: 1.994397e-13 
#> iter: 1  mgc: 1.994397e-13 
#> outer mgc:  0.01003301 
#> iter: 1  value: 1614.06 mgc: 0.01011325 ustep: 1 
#> iter: 2  value: 1614.06 mgc: 2.306859e-07 ustep: 1 
#> iter: 3  mgc: 1.347303e-13 
#> iter: 1  mgc: 1.347303e-13 
#> outer mgc:  0.03530474 
#> iter: 1  value: 1614.06 mgc: 0.007067646 ustep: 1 
#> iter: 2  value: 1614.06 mgc: 7.959323e-08 ustep: 1 
#> mgc: 1.454361e-13 
#> iter: 1  mgc: 1.454361e-13 
#> outer mgc:  0.03785864 
#> iter: 1  value: 1614.06 mgc: 0.001382708 ustep: 1 
#> iter: 2  mgc: 2.785879e-09 
#> iter: 1  value: 1614.06 mgc: 0.002933668 ustep: 1 
#> iter: 2  mgc: 6.160355e-09 
#> iter: 1  mgc: 6.160355e-09 
#> outer mgc:  0.02495388 
#> iter: 1  value: 1614.058 mgc: 0.01411416 ustep: 1 
#> iter: 2  value: 1614.058 mgc: 2.834599e-07 ustep: 1 
#> iter: 3  mgc: 1.866637e-13 
#> iter: 1  value: 1614.059 mgc: 0.005700018 ustep: 1 
#> iter: 2  value: 1614.059 mgc: 4.641769e-08 ustep: 1 
#> mgc: 1.366045e-13 
#> iter: 1  mgc: 1.366045e-13 
#> outer mgc:  0.007098823 
#> iter: 1  value: 1614.059 mgc: 0.002016309 ustep: 1 
#> iter: 2  mgc: 8.135752e-09 
#> iter: 1  mgc: 8.135752e-09 
#> outer mgc:  0.02088733 
#> iter: 1  value: 1614.06 mgc: 0.002118894 ustep: 1 
#> iter: 2  mgc: 6.45447e-09 
#> iter: 1  mgc: 6.45447e-09 
#> outer mgc:  0.00933426 
#> iter: 1  value: 1614.06 mgc: 0.003536252 ustep: 1 
#> iter: 2  value: 1614.06 mgc: 1.410383e-08 ustep: 1 
#> mgc: 1.218946e-13 
#> iter: 1  mgc: 1.218946e-13 
#> outer mgc:  0.01930622 
#> iter: 1  value: 1614.059 mgc: 0.003001299 ustep: 1 
#> iter: 2  value: 1614.059 mgc: 1.286718e-08 ustep: 1 
#> mgc: 1.719327e-13 
#> iter: 1  mgc: 1.719327e-13 
#> outer mgc:  0.02258169 
#> iter: 1  value: 1614.06 mgc: 0.003532735 ustep: 1 
#> iter: 2  value: 1614.06 mgc: 3.163432e-08 ustep: 1 
#> mgc: 1.361681e-13 
#> iter: 1  mgc: 1.361681e-13 
#> outer mgc:  0.01023574 
#> iter: 1  value: 1614.06 mgc: 0.006019856 ustep: 1 
#> iter: 2  value: 1614.06 mgc: 2.925762e-08 ustep: 1 
#> mgc: 1.452813e-13 
#> iter: 1  value: 1614.06 mgc: 0.001999503 ustep: 1 
#> iter: 2  mgc: 3.22743e-09 
#> iter: 1  mgc: 3.22743e-09 
#> outer mgc:  0.01059122 
#> iter: 1  value: 1614.06 mgc: 0.001755182 ustep: 1 
#> iter: 2  mgc: 6.715651e-09 
#> iter: 1  mgc: 6.715651e-09 
#> outer mgc:  0.008470149 
#> iter: 1  value: 1614.061 mgc: 0.005591367 ustep: 1 
#> iter: 2  value: 1614.061 mgc: 2.615277e-08 ustep: 1 
#> mgc: 1.001327e-13 
#> iter: 1  value: 1614.06 mgc: 0.001279451 ustep: 1 
#> iter: 2  mgc: 1.36939e-09 
#> iter: 1  mgc: 1.36939e-09 
#> outer mgc:  0.004782592 
#> iter: 1  value: 1614.06 mgc: 0.0003964001 ustep: 1 
#> iter: 2  mgc: 2.15897e-10 
#> iter: 1  value: 1614.06 mgc: 0.000571036 ustep: 1 
#> iter: 2  mgc: 6.52881e-10 
#> iter: 1  mgc: 6.52881e-10 
#> outer mgc:  0.005151132 
#> iter: 1  value: 1614.06 mgc: 0.002207468 ustep: 1 
#> iter: 2  value: 1614.06 mgc: 2.657116e-08 ustep: 1 
#> mgc: 1.247832e-13 
#> iter: 1  value: 1614.06 mgc: 0.0008529652 ustep: 1 
#> iter: 2  mgc: 3.967773e-09 
#> iter: 1  mgc: 3.967773e-09 
#> outer mgc:  0.003146232 
#> iter: 1  value: 1614.06 mgc: 0.0006089841 ustep: 1 
#> iter: 2  mgc: 1.25912e-09 
#> iter: 1  mgc: 1.25912e-09 
#> outer mgc:  0.000870969 
#> iter: 1  value: 1614.06 mgc: 0.0002836605 ustep: 1 
#> iter: 2  mgc: 1.216779e-10 
#> iter: 1  value: 1614.06 mgc: 0.0008283047 ustep: 1 
#> iter: 2  mgc: 1.006922e-09 
#> iter: 1  mgc: 1.006922e-09 
#> outer mgc:  0.00128074 
#> iter: 1  value: 1614.06 mgc: 0.001657506 ustep: 1 
#> iter: 2  mgc: 3.757201e-09 
#> iter: 1  mgc: 1.006922e-09 
#> iter: 1  mgc: 1.006922e-09 
#> outer mgc:  0.00128074 
#> iter: 1  mgc: 1.006922e-09 
#> outer mgc:  0.00128074 
#> Retro Peel: 5
#> Using stats::nlminb for optimization with these control parameters:
#> iter.max: 1000
#> eval.max: 1000
#> iter: 1  value: 1587.276 mgc: 1.645174 ustep: 1 
#> iter: 2  value: 1587.275 mgc: 0.07524493 ustep: 1 
#> iter: 3  value: 1587.275 mgc: 0.0001813551 ustep: 1 
#> iter: 4  mgc: 1.060953e-09 
#> iter: 1  mgc: 1.060953e-09 
#> Matching hessian patterns... Done
#> outer mgc:  15.5325 
#> iter: 1  value: 1864.233 mgc: 22.01379 ustep: 1 
#> iter: 2  value: 1864.013 mgc: 1.150604 ustep: 1 
#> iter: 3  value: 1864.013 mgc: 0.006958672 ustep: 1 
#> iter: 4  value: 1864.013 mgc: 7.005142e-07 ustep: 1 
#> iter: 5  mgc: 1.535246e-13 
#> iter: 1  value: 1587.399 mgc: 1.865596 ustep: 1 
#> iter: 2  value: 1587.399 mgc: 0.01020216 ustep: 1 
#> iter: 3  value: 1587.399 mgc: 8.634594e-07 ustep: 1 
#> iter: 4  mgc: 1.390773e-13 
#> iter: 1  mgc: 1.390773e-13 
#> outer mgc:  49.74289 
#> iter: 1  value: 1601.458 mgc: 11.03159 ustep: 1 
#> iter: 2  value: 1601.433 mgc: 0.2396149 ustep: 1 
#> iter: 3  value: 1601.433 mgc: 0.003378865 ustep: 1 
#> iter: 4  value: 1601.433 mgc: 6.728733e-07 ustep: 1 
#> iter: 5  mgc: 1.520853e-13 
#> iter: 1  value: 1586.737 mgc: 1.83494 ustep: 1 
#> iter: 2  value: 1586.737 mgc: 0.005700874 ustep: 1 
#> iter: 3  value: 1586.737 mgc: 1.53498e-06 ustep: 1 
#> iter: 4  mgc: 2.381665e-13 
#> iter: 1  mgc: 2.381665e-13 
#> outer mgc:  32.41291 
#> iter: 1  value: 1586.48 mgc: 0.7437807 ustep: 1 
#> iter: 2  value: 1586.48 mgc: 0.000658708 ustep: 1 
#> iter: 3  value: 1586.48 mgc: 2.340351e-08 ustep: 1 
#> iter: 4  mgc: 1.28466e-13 
#> iter: 1  mgc: 1.28466e-13 
#> outer mgc:  33.8438 
#> iter: 1  value: 1586.198 mgc: 1.696352 ustep: 1 
#> iter: 2  value: 1586.198 mgc: 0.005018766 ustep: 1 
#> iter: 3  value: 1586.198 mgc: 1.122937e-06 ustep: 1 
#> iter: 4  mgc: 1.391497e-13 
#> iter: 1  mgc: 1.391497e-13 
#> outer mgc:  26.60868 
#> iter: 1  value: 1585.719 mgc: 1.974849 ustep: 1 
#> iter: 2  value: 1585.719 mgc: 0.005435438 ustep: 1 
#> iter: 3  value: 1585.719 mgc: 1.35662e-06 ustep: 1 
#> iter: 4  mgc: 1.203095e-13 
#> iter: 1  mgc: 1.203095e-13 
#> outer mgc:  17.9691 
#> iter: 1  value: 1585.442 mgc: 0.8288661 ustep: 1 
#> iter: 2  value: 1585.442 mgc: 0.001357994 ustep: 1 
#> iter: 3  value: 1585.442 mgc: 8.584037e-08 ustep: 1 
#> iter: 4  mgc: 2.157255e-13 
#> iter: 1  mgc: 2.157255e-13 
#> outer mgc:  14.06058 
#> iter: 1  value: 1585.165 mgc: 1.006743 ustep: 1 
#> iter: 2  value: 1585.165 mgc: 0.001075362 ustep: 1 
#> iter: 3  value: 1585.165 mgc: 2.536529e-08 ustep: 1 
#> iter: 4  mgc: 1.167556e-13 
#> iter: 1  mgc: 1.167556e-13 
#> outer mgc:  22.62755 
#> iter: 1  value: 1584.931 mgc: 2.435992 ustep: 1 
#> iter: 2  value: 1584.931 mgc: 0.006531593 ustep: 1 
#> iter: 3  value: 1584.931 mgc: 3.645074e-06 ustep: 1 
#> iter: 4  mgc: 2.395684e-12 
#> iter: 1  mgc: 2.395684e-12 
#> outer mgc:  17.37382 
#> iter: 1  value: 1584.141 mgc: 2.975621 ustep: 1 
#> iter: 2  value: 1584.141 mgc: 0.007742179 ustep: 1 
#> iter: 3  value: 1584.141 mgc: 4.293975e-06 ustep: 1 
#> iter: 4  mgc: 3.510923e-12 
#> iter: 1  mgc: 3.510923e-12 
#> outer mgc:  16.30483 
#> iter: 1  value: 1591.038 mgc: 7.367296 ustep: 1 
#> iter: 2  value: 1591.033 mgc: 0.100458 ustep: 1 
#> iter: 3  value: 1591.033 mgc: 0.0005145526 ustep: 1 
#> iter: 4  value: 1591.033 mgc: 1.561649e-08 ustep: 1 
#> mgc: 1.527791e-13 
#> iter: 1  value: 1584.192 mgc: 1.660018 ustep: 1 
#> iter: 2  value: 1584.192 mgc: 0.004662379 ustep: 1 
#> iter: 3  value: 1584.192 mgc: 9.989905e-07 ustep: 1 
#> iter: 4  mgc: 1.264615e-13 
#> iter: 1  value: 1584.018 mgc: 0.6417835 ustep: 1 
#> iter: 2  value: 1584.018 mgc: 0.0006870869 ustep: 1 
#> iter: 3  value: 1584.018 mgc: 2.130077e-08 ustep: 1 
#> mgc: 1.499585e-13 
#> iter: 1  mgc: 1.499585e-13 
#> outer mgc:  5.974428 
#> iter: 1  value: 1583.979 mgc: 0.3707909 ustep: 1 
#> iter: 2  value: 1583.979 mgc: 0.0001966456 ustep: 1 
#> iter: 3  mgc: 1.645332e-09 
#> iter: 1  mgc: 1.645332e-09 
#> outer mgc:  5.881787 
#> iter: 1  value: 1583.9 mgc: 0.4808516 ustep: 1 
#> iter: 2  value: 1583.9 mgc: 0.0004021284 ustep: 1 
#> iter: 3  mgc: 7.490667e-09 
#> iter: 1  mgc: 7.490667e-09 
#> outer mgc:  7.471724 
#> iter: 1  value: 1583.824 mgc: 0.417227 ustep: 1 
#> iter: 2  value: 1583.824 mgc: 0.0002702676 ustep: 1 
#> iter: 3  mgc: 3.034496e-09 
#> iter: 1  mgc: 3.034496e-09 
#> outer mgc:  5.584954 
#> iter: 1  value: 1583.671 mgc: 0.3527971 ustep: 1 
#> iter: 2  value: 1583.671 mgc: 0.0002260026 ustep: 1 
#> iter: 3  mgc: 2.572637e-09 
#> iter: 1  value: 1583.472 mgc: 0.1575499 ustep: 1 
#> iter: 2  value: 1583.472 mgc: 4.528204e-05 ustep: 1 
#> iter: 3  mgc: 1.221941e-10 
#> iter: 1  value: 1582.83 mgc: 0.6393328 ustep: 1 
#> iter: 2  value: 1582.83 mgc: 0.0007452958 ustep: 1 
#> iter: 3  value: 1582.83 mgc: 3.222238e-08 ustep: 1 
#> mgc: 8.785033e-14 
#> iter: 1  mgc: 8.785033e-14 
#> outer mgc:  13.25885 
#> iter: 1  value: 1582.505 mgc: 0.7831688 ustep: 1 
#> iter: 2  value: 1582.505 mgc: 0.005925899 ustep: 1 
#> iter: 3  value: 1582.505 mgc: 1.383135e-06 ustep: 1 
#> iter: 4  mgc: 1.122827e-13 
#> iter: 1  value: 1582.192 mgc: 0.7557006 ustep: 1 
#> iter: 2  value: 1582.192 mgc: 0.002418785 ustep: 1 
#> iter: 3  value: 1582.192 mgc: 8.18608e-08 ustep: 1 
#> iter: 4  mgc: 1.711857e-13 
#> iter: 1  mgc: 1.711857e-13 
#> outer mgc:  21.86292 
#> iter: 1  value: 1582.717 mgc: 3.386025 ustep: 1 
#> iter: 2  value: 1582.716 mgc: 0.01325541 ustep: 1 
#> iter: 3  value: 1582.716 mgc: 1.052983e-05 ustep: 1 
#> iter: 4  mgc: 6.413916e-12 
#> iter: 1  value: 1582.36 mgc: 2.434936 ustep: 1 
#> iter: 2  value: 1582.36 mgc: 0.008082376 ustep: 1 
#> iter: 3  value: 1582.36 mgc: 3.810049e-06 ustep: 1 
#> iter: 4  mgc: 8.6541e-13 
#> iter: 1  value: 1582.279 mgc: 2.16728 ustep: 1 
#> iter: 2  value: 1582.279 mgc: 0.007093079 ustep: 1 
#> iter: 3  value: 1582.279 mgc: 2.492335e-06 ustep: 1 
#> iter: 4  mgc: 3.844863e-13 
#> iter: 1  value: 1582.029 mgc: 1.365837 ustep: 1 
#> iter: 2  value: 1582.029 mgc: 0.003039677 ustep: 1 
#> iter: 3  value: 1582.029 mgc: 3.93765e-07 ustep: 1 
#> iter: 4  mgc: 1.212336e-13 
#> iter: 1  mgc: 1.212336e-13 
#> outer mgc:  15.38572 
#> iter: 1  value: 1581.967 mgc: 0.8204384 ustep: 1 
#> iter: 2  value: 1581.967 mgc: 0.001187957 ustep: 1 
#> iter: 3  value: 1581.967 mgc: 4.527142e-08 ustep: 1 
#> iter: 4  mgc: 1.254207e-13 
#> iter: 1  mgc: 1.254207e-13 
#> outer mgc:  9.118484 
#> iter: 1  value: 1581.92 mgc: 0.8407738 ustep: 1 
#> iter: 2  value: 1581.92 mgc: 0.0005500608 ustep: 1 
#> iter: 3  value: 1581.92 mgc: 1.929094e-08 ustep: 1 
#> iter: 4  mgc: 1.405814e-13 
#> iter: 1  mgc: 1.405814e-13 
#> outer mgc:  10.57798 
#> iter: 1  value: 1581.877 mgc: 0.5943421 ustep: 1 
#> iter: 2  value: 1581.877 mgc: 0.0009156328 ustep: 1 
#> iter: 3  value: 1581.877 mgc: 2.363184e-08 ustep: 1 
#> iter: 4  mgc: 1.365715e-13 
#> iter: 1  mgc: 1.365715e-13 
#> outer mgc:  7.458686 
#> iter: 1  value: 1581.709 mgc: 0.4580262 ustep: 1 
#> iter: 2  value: 1581.709 mgc: 0.0003708431 ustep: 1 
#> iter: 3  mgc: 1.433748e-09 
#> iter: 1  value: 1581.485 mgc: 0.4514467 ustep: 1 
#> iter: 2  value: 1581.485 mgc: 0.000538973 ustep: 1 
#> iter: 3  mgc: 8.960611e-10 
#> iter: 1  value: 1581.175 mgc: 0.888888 ustep: 1 
#> iter: 2  value: 1581.175 mgc: 0.002095013 ustep: 1 
#> iter: 3  value: 1581.175 mgc: 1.351054e-08 ustep: 1 
#> mgc: 1.593134e-13 
#> iter: 1  mgc: 1.593134e-13 
#> outer mgc:  16.55595 
#> iter: 1  value: 1582.633 mgc: 3.561223 ustep: 1 
#> iter: 2  value: 1582.633 mgc: 0.01284681 ustep: 1 
#> iter: 3  value: 1582.633 mgc: 1.010251e-05 ustep: 1 
#> iter: 4  mgc: 5.831462e-12 
#> iter: 1  value: 1581.212 mgc: 1.029636 ustep: 1 
#> iter: 2  value: 1581.212 mgc: 0.002750918 ustep: 1 
#> iter: 3  value: 1581.212 mgc: 1.327681e-07 ustep: 1 
#> iter: 4  mgc: 1.371066e-13 
#> iter: 1  mgc: 1.371066e-13 
#> outer mgc:  16.97868 
#> iter: 1  value: 1581.226 mgc: 1.240195 ustep: 1 
#> iter: 2  value: 1581.226 mgc: 0.001634393 ustep: 1 
#> iter: 3  value: 1581.226 mgc: 8.297382e-08 ustep: 1 
#> iter: 4  mgc: 1.473406e-13 
#> iter: 1  mgc: 1.473406e-13 
#> outer mgc:  10.46867 
#> iter: 1  value: 1581.507 mgc: 1.580385 ustep: 1 
#> iter: 2  value: 1581.507 mgc: 0.004289211 ustep: 1 
#> iter: 3  value: 1581.507 mgc: 1.013564e-06 ustep: 1 
#> iter: 4  mgc: 1.418409e-13 
#> iter: 1  value: 1581.303 mgc: 0.9204828 ustep: 1 
#> iter: 2  value: 1581.303 mgc: 0.001291847 ustep: 1 
#> iter: 3  value: 1581.303 mgc: 8.252679e-08 ustep: 1 
#> iter: 4  mgc: 1.874311e-13 
#> iter: 1  value: 1581.272 mgc: 0.8235783 ustep: 1 
#> iter: 2  value: 1581.272 mgc: 0.001037904 ustep: 1 
#> iter: 3  value: 1581.272 mgc: 4.389471e-08 ustep: 1 
#> iter: 4  mgc: 1.716973e-13 
#> iter: 1  value: 1581.223 mgc: 0.4810571 ustep: 1 
#> iter: 2  value: 1581.223 mgc: 0.0003561551 ustep: 1 
#> iter: 3  mgc: 5.083264e-09 
#> iter: 1  mgc: 5.083264e-09 
#> outer mgc:  8.276689 
#> iter: 1  value: 1581.186 mgc: 0.2745367 ustep: 1 
#> iter: 2  value: 1581.186 mgc: 0.0001327297 ustep: 1 
#> iter: 3  mgc: 6.900775e-10 
#> iter: 1  mgc: 6.900775e-10 
#> outer mgc:  7.763997 
#> iter: 1  value: 1581.168 mgc: 0.3232395 ustep: 1 
#> iter: 2  value: 1581.168 mgc: 0.0001471063 ustep: 1 
#> iter: 3  mgc: 9.53383e-10 
#> iter: 1  mgc: 9.53383e-10 
#> outer mgc:  6.94323 
#> iter: 1  value: 1581.127 mgc: 0.2626311 ustep: 1 
#> iter: 2  value: 1581.127 mgc: 0.0001155559 ustep: 1 
#> iter: 3  mgc: 5.592392e-10 
#> iter: 1  mgc: 5.592392e-10 
#> outer mgc:  6.407543 
#> iter: 1  value: 1581.121 mgc: 0.3039189 ustep: 1 
#> iter: 2  value: 1581.121 mgc: 0.0001484868 ustep: 1 
#> iter: 3  mgc: 7.975359e-10 
#> iter: 1  value: 1581.124 mgc: 0.1864443 ustep: 1 
#> iter: 2  value: 1581.124 mgc: 6.154738e-05 ustep: 1 
#> iter: 3  mgc: 9.978891e-11 
#> iter: 1  value: 1581.17 mgc: 0.4165274 ustep: 1 
#> iter: 2  value: 1581.17 mgc: 0.000306547 ustep: 1 
#> iter: 3  mgc: 2.484067e-09 
#> iter: 1  mgc: 2.484067e-09 
#> outer mgc:  4.652116 
#> iter: 1  value: 1581.623 mgc: 2.521121 ustep: 1 
#> iter: 2  value: 1581.623 mgc: 0.02242759 ustep: 1 
#> iter: 3  value: 1581.623 mgc: 5.949404e-06 ustep: 1 
#> iter: 4  mgc: 4.134977e-13 
#> iter: 1  value: 1581.119 mgc: 0.7093263 ustep: 1 
#> iter: 2  value: 1581.119 mgc: 0.002129211 ustep: 1 
#> iter: 3  value: 1581.119 mgc: 5.438671e-08 ustep: 1 
#> iter: 4  mgc: 1.46905e-13 
#> iter: 1  mgc: 1.46905e-13 
#> outer mgc:  6.124379 
#> iter: 1  value: 1581.236 mgc: 2.055108 ustep: 1 
#> iter: 2  value: 1581.236 mgc: 0.003954741 ustep: 1 
#> iter: 3  value: 1581.236 mgc: 6.695131e-07 ustep: 1 
#> iter: 4  mgc: 1.827863e-13 
#> iter: 1  value: 1581.267 mgc: 1.400296 ustep: 1 
#> iter: 2  value: 1581.267 mgc: 0.002923021 ustep: 1 
#> iter: 3  value: 1581.267 mgc: 3.716652e-07 ustep: 1 
#> iter: 4  mgc: 1.159491e-13 
#> iter: 1  value: 1581.176 mgc: 0.8341392 ustep: 1 
#> iter: 2  value: 1581.176 mgc: 0.001094239 ustep: 1 
#> iter: 3  value: 1581.176 mgc: 5.144763e-08 ustep: 1 
#> iter: 4  mgc: 1.359853e-13 
#> iter: 1  value: 1581.109 mgc: 0.2555445 ustep: 1 
#> iter: 2  value: 1581.109 mgc: 0.0001035448 ustep: 1 
#> iter: 3  mgc: 4.6524e-10 
#> iter: 1  mgc: 4.6524e-10 
#> outer mgc:  6.288324 
#> iter: 1  value: 1581.08 mgc: 0.2751572 ustep: 1 
#> iter: 2  value: 1581.08 mgc: 0.0001206637 ustep: 1 
#> iter: 3  mgc: 5.421937e-10 
#> iter: 1  mgc: 5.421937e-10 
#> outer mgc:  5.253827 
#> iter: 1  value: 1581.067 mgc: 0.2466948 ustep: 1 
#> iter: 2  value: 1581.067 mgc: 9.315478e-05 ustep: 1 
#> iter: 3  mgc: 3.799387e-10 
#> iter: 1  mgc: 3.799387e-10 
#> outer mgc:  4.826603 
#> iter: 1  value: 1581.019 mgc: 0.2358301 ustep: 1 
#> iter: 2  value: 1581.019 mgc: 0.0001004252 ustep: 1 
#> iter: 3  mgc: 2.496654e-10 
#> iter: 1  mgc: 2.496654e-10 
#> outer mgc:  3.847128 
#> iter: 1  value: 1580.946 mgc: 0.4543508 ustep: 1 
#> iter: 2  value: 1580.946 mgc: 0.000323353 ustep: 1 
#> iter: 3  mgc: 1.160182e-09 
#> iter: 1  mgc: 1.160182e-09 
#> outer mgc:  2.227555 
#> iter: 1  value: 1580.926 mgc: 0.4432953 ustep: 1 
#> iter: 2  value: 1580.926 mgc: 0.0001160645 ustep: 1 
#> iter: 3  mgc: 4.957793e-10 
#> iter: 1  mgc: 4.957793e-10 
#> outer mgc:  1.896471 
#> iter: 1  value: 1581.125 mgc: 0.7717878 ustep: 1 
#> iter: 2  value: 1581.125 mgc: 0.001102056 ustep: 1 
#> iter: 3  value: 1581.125 mgc: 6.075844e-08 ustep: 1 
#> iter: 4  mgc: 1.260006e-13 
#> iter: 1  value: 1580.983 mgc: 0.6158653 ustep: 1 
#> iter: 2  value: 1580.983 mgc: 0.0005262931 ustep: 1 
#> iter: 3  value: 1580.983 mgc: 1.451137e-08 ustep: 1 
#> mgc: 2.432976e-13 
#> iter: 1  value: 1580.922 mgc: 0.3134551 ustep: 1 
#> iter: 2  value: 1580.922 mgc: 0.000135721 ustep: 1 
#> iter: 3  mgc: 1.175081e-09 
#> iter: 1  value: 1580.92 mgc: 0.152157 ustep: 1 
#> iter: 2  value: 1580.92 mgc: 3.202852e-05 ustep: 1 
#> iter: 3  mgc: 6.566701e-11 
#> iter: 1  mgc: 6.566701e-11 
#> outer mgc:  3.266205 
#> iter: 1  value: 1580.911 mgc: 0.1620927 ustep: 1 
#> iter: 2  value: 1580.911 mgc: 3.836169e-05 ustep: 1 
#> iter: 3  mgc: 5.742436e-11 
#> iter: 1  mgc: 5.742436e-11 
#> outer mgc:  1.844185 
#> iter: 1  value: 1580.915 mgc: 0.1251156 ustep: 1 
#> iter: 2  value: 1580.915 mgc: 3.024366e-05 ustep: 1 
#> iter: 3  mgc: 4.235785e-11 
#> iter: 1  mgc: 4.235785e-11 
#> outer mgc:  2.565353 
#> iter: 1  value: 1580.9 mgc: 0.2022081 ustep: 1 
#> iter: 2  value: 1580.9 mgc: 3.971075e-05 ustep: 1 
#> iter: 3  mgc: 8.336425e-11 
#> iter: 1  mgc: 8.336425e-11 
#> outer mgc:  2.062131 
#> iter: 1  value: 1580.83 mgc: 0.2175372 ustep: 1 
#> iter: 2  value: 1580.83 mgc: 0.0001102253 ustep: 1 
#> iter: 3  mgc: 8.063546e-10 
#> iter: 1  mgc: 8.063546e-10 
#> outer mgc:  1.928412 
#> iter: 1  value: 1580.723 mgc: 0.4046496 ustep: 1 
#> iter: 2  value: 1580.723 mgc: 0.0001451409 ustep: 1 
#> iter: 3  mgc: 1.098716e-09 
#> iter: 1  mgc: 1.098716e-09 
#> outer mgc:  3.788756 
#> iter: 1  value: 1580.81 mgc: 0.3737773 ustep: 1 
#> iter: 2  value: 1580.81 mgc: 0.000409838 ustep: 1 
#> iter: 3  mgc: 8.291801e-09 
#> iter: 1  mgc: 8.291801e-09 
#> outer mgc:  5.752619 
#> iter: 1  value: 1580.917 mgc: 0.3597538 ustep: 1 
#> iter: 2  value: 1580.917 mgc: 0.0002722156 ustep: 1 
#> iter: 3  mgc: 4.199562e-09 
#> iter: 1  mgc: 4.199562e-09 
#> outer mgc:  2.767175 
#> iter: 1  value: 1580.449 mgc: 0.3478299 ustep: 1 
#> iter: 2  value: 1580.449 mgc: 0.0003620954 ustep: 1 
#> iter: 3  mgc: 3.894204e-09 
#> iter: 1  mgc: 3.894204e-09 
#> outer mgc:  1.89849 
#> iter: 1  value: 1580.629 mgc: 1.132866 ustep: 1 
#> iter: 2  value: 1580.629 mgc: 0.001990366 ustep: 1 
#> iter: 3  value: 1580.629 mgc: 2.673679e-07 ustep: 1 
#> iter: 4  mgc: 1.382329e-13 
#> iter: 1  value: 1580.457 mgc: 0.341507 ustep: 1 
#> iter: 2  value: 1580.457 mgc: 0.0001839544 ustep: 1 
#> iter: 3  mgc: 2.279793e-09 
#> iter: 1  mgc: 2.279793e-09 
#> outer mgc:  4.349848 
#> iter: 1  value: 1580.438 mgc: 0.4571376 ustep: 1 
#> iter: 2  value: 1580.438 mgc: 0.000260726 ustep: 1 
#> iter: 3  mgc: 1.777486e-09 
#> iter: 1  mgc: 1.777486e-09 
#> outer mgc:  4.666428 
#> iter: 1  value: 1580.442 mgc: 0.1395692 ustep: 1 
#> iter: 2  value: 1580.442 mgc: 5.879488e-05 ustep: 1 
#> iter: 3  mgc: 2.104287e-10 
#> iter: 1  mgc: 2.104287e-10 
#> outer mgc:  1.438685 
#> iter: 1  value: 1580.438 mgc: 0.1805669 ustep: 1 
#> iter: 2  value: 1580.438 mgc: 5.845115e-05 ustep: 1 
#> iter: 3  mgc: 1.025701e-10 
#> iter: 1  value: 1580.441 mgc: 0.1301288 ustep: 1 
#> iter: 2  value: 1580.441 mgc: 0.0001441547 ustep: 1 
#> iter: 3  mgc: 6.237481e-10 
#> iter: 1  value: 1580.461 mgc: 0.2290869 ustep: 1 
#> iter: 2  value: 1580.461 mgc: 0.0004442352 ustep: 1 
#> iter: 3  mgc: 5.921354e-09 
#> iter: 1  mgc: 5.921354e-09 
#> outer mgc:  3.942095 
#> iter: 1  value: 1580.669 mgc: 0.8580031 ustep: 1 
#> iter: 2  value: 1580.669 mgc: 0.0008610215 ustep: 1 
#> iter: 3  mgc: 3.06812e-09 
#> iter: 1  mgc: 3.06812e-09 
#> outer mgc:  3.670379 
#> iter: 1  value: 1581.013 mgc: 2.474151 ustep: 1 
#> iter: 2  value: 1581.013 mgc: 0.01192936 ustep: 1 
#> iter: 3  value: 1581.013 mgc: 5.806424e-06 ustep: 1 
#> iter: 4  mgc: 1.651471e-12 
#> iter: 1  value: 1580.968 mgc: 1.661589 ustep: 1 
#> iter: 2  value: 1580.968 mgc: 0.004085018 ustep: 1 
#> iter: 3  value: 1580.968 mgc: 6.202061e-07 ustep: 1 
#> iter: 4  mgc: 1.751348e-13 
#> iter: 1  value: 1580.66 mgc: 0.2334084 ustep: 1 
#> iter: 2  value: 1580.66 mgc: 8.210817e-05 ustep: 1 
#> iter: 3  mgc: 2.562819e-10 
#> iter: 1  mgc: 2.562819e-10 
#> outer mgc:  4.619564 
#> iter: 1  value: 1580.634 mgc: 0.1986704 ustep: 1 
#> iter: 2  value: 1580.634 mgc: 5.967807e-05 ustep: 1 
#> iter: 3  mgc: 1.407368e-10 
#> iter: 1  mgc: 1.407368e-10 
#> outer mgc:  3.366107 
#> iter: 1  value: 1580.614 mgc: 0.1827201 ustep: 1 
#> iter: 2  value: 1580.614 mgc: 5.869581e-05 ustep: 1 
#> iter: 3  mgc: 1.53765e-10 
#> iter: 1  mgc: 1.53765e-10 
#> outer mgc:  3.555479 
#> iter: 1  value: 1580.584 mgc: 0.1721998 ustep: 1 
#> iter: 2  value: 1580.584 mgc: 4.162134e-05 ustep: 1 
#> iter: 3  mgc: 5.105614e-11 
#> iter: 1  mgc: 5.105614e-11 
#> outer mgc:  2.945313 
#> iter: 1  value: 1580.633 mgc: 0.2316421 ustep: 1 
#> iter: 2  value: 1580.633 mgc: 0.0001170506 ustep: 1 
#> iter: 3  mgc: 2.590627e-10 
#> iter: 1  mgc: 2.590627e-10 
#> outer mgc:  3.032964 
#> iter: 1  value: 1580.248 mgc: 0.3052068 ustep: 1 
#> iter: 2  value: 1580.248 mgc: 0.0004082469 ustep: 1 
#> iter: 3  mgc: 5.533161e-09 
#> iter: 1  mgc: 5.533161e-09 
#> outer mgc:  2.363199 
#> iter: 1  value: 1581.476 mgc: 2.395658 ustep: 1 
#> iter: 2  value: 1581.475 mgc: 0.01078552 ustep: 1 
#> iter: 3  value: 1581.475 mgc: 4.919037e-06 ustep: 1 
#> iter: 4  mgc: 1.346228e-12 
#> iter: 1  value: 1580.242 mgc: 0.2467267 ustep: 1 
#> iter: 2  value: 1580.242 mgc: 0.0001125551 ustep: 1 
#> iter: 3  mgc: 5.18115e-10 
#> iter: 1  mgc: 5.18115e-10 
#> outer mgc:  5.089911 
#> iter: 1  value: 1580.243 mgc: 0.2570009 ustep: 1 
#> iter: 2  value: 1580.243 mgc: 0.0001020525 ustep: 1 
#> iter: 3  mgc: 4.123674e-10 
#> iter: 1  mgc: 4.123674e-10 
#> outer mgc:  1.493673 
#> iter: 1  value: 1580.245 mgc: 0.1466035 ustep: 1 
#> iter: 2  value: 1580.245 mgc: 3.440849e-05 ustep: 1 
#> iter: 3  mgc: 4.558735e-11 
#> iter: 1  mgc: 4.558735e-11 
#> outer mgc:  2.624751 
#> iter: 1  value: 1580.256 mgc: 0.1131077 ustep: 1 
#> iter: 2  value: 1580.256 mgc: 3.555184e-05 ustep: 1 
#> iter: 3  mgc: 5.735973e-11 
#> iter: 1  mgc: 5.735973e-11 
#> outer mgc:  1.314377 
#> iter: 1  value: 1580.204 mgc: 0.2252297 ustep: 1 
#> iter: 2  value: 1580.204 mgc: 6.269419e-05 ustep: 1 
#> iter: 3  mgc: 3.408046e-11 
#> iter: 1  value: 1580.134 mgc: 0.1803915 ustep: 1 
#> iter: 2  value: 1580.134 mgc: 3.631101e-05 ustep: 1 
#> iter: 3  mgc: 2.512049e-11 
#> iter: 1  mgc: 2.512049e-11 
#> outer mgc:  1.553413 
#> iter: 1  value: 1580.219 mgc: 0.2014511 ustep: 1 
#> iter: 2  value: 1580.219 mgc: 0.000571323 ustep: 1 
#> iter: 3  mgc: 7.228908e-09 
#> iter: 1  value: 1580.415 mgc: 0.26006 ustep: 1 
#> iter: 2  value: 1580.415 mgc: 0.0007756203 ustep: 1 
#> iter: 3  value: 1580.415 mgc: 1.441899e-08 ustep: 1 
#> mgc: 2.150686e-13 
#> iter: 1  mgc: 2.150686e-13 
#> outer mgc:  7.890315 
#> iter: 1  value: 1580.311 mgc: 1.475154 ustep: 1 
#> iter: 2  value: 1580.311 mgc: 0.006237603 ustep: 1 
#> iter: 3  value: 1580.311 mgc: 1.615729e-07 ustep: 1 
#> iter: 4  mgc: 1.337526e-13 
#> iter: 1  mgc: 1.337526e-13 
#> outer mgc:  5.54294 
#> iter: 1  value: 1580.476 mgc: 1.43 ustep: 1 
#> iter: 2  value: 1580.476 mgc: 0.004667591 ustep: 1 
#> iter: 3  value: 1580.476 mgc: 1.270928e-07 ustep: 1 
#> iter: 4  mgc: 1.374281e-13 
#> iter: 1  mgc: 1.374281e-13 
#> outer mgc:  5.190807 
#> iter: 1  value: 1583.435 mgc: 5.563233 ustep: 1 
#> iter: 2  value: 1583.434 mgc: 0.03654386 ustep: 1 
#> iter: 3  value: 1583.434 mgc: 5.582332e-05 ustep: 1 
#> iter: 4  mgc: 1.813861e-10 
#> iter: 1  value: 1580.631 mgc: 1.464144 ustep: 1 
#> iter: 2  value: 1580.631 mgc: 0.003282003 ustep: 1 
#> iter: 3  value: 1580.631 mgc: 4.312624e-07 ustep: 1 
#> iter: 4  mgc: 9.956743e-14 
#> iter: 1  value: 1580.45 mgc: 0.2675943 ustep: 1 
#> iter: 2  value: 1580.45 mgc: 0.0001078471 ustep: 1 
#> iter: 3  mgc: 4.57131e-10 
#> iter: 1  mgc: 4.57131e-10 
#> outer mgc:  3.315373 
#> iter: 1  value: 1580.438 mgc: 0.06186912 ustep: 1 
#> iter: 2  value: 1580.438 mgc: 5.421391e-06 ustep: 1 
#> iter: 3  mgc: 1.467194e-12 
#> iter: 1  value: 1580.424 mgc: 0.1004239 ustep: 1 
#> iter: 2  value: 1580.424 mgc: 1.429361e-05 ustep: 1 
#> iter: 3  mgc: 9.894345e-12 
#> iter: 1  mgc: 9.894345e-12 
#> outer mgc:  3.60707 
#> iter: 1  value: 1580.416 mgc: 0.3091841 ustep: 1 
#> iter: 2  value: 1580.416 mgc: 6.343719e-05 ustep: 1 
#> iter: 3  mgc: 1.41035e-10 
#> iter: 1  mgc: 1.41035e-10 
#> outer mgc:  2.168454 
#> iter: 1  value: 1580.411 mgc: 0.06896028 ustep: 1 
#> iter: 2  value: 1580.411 mgc: 2.735706e-05 ustep: 1 
#> iter: 3  mgc: 3.264223e-11 
#> iter: 1  mgc: 3.264223e-11 
#> outer mgc:  1.53508 
#> iter: 1  value: 1580.347 mgc: 0.1402833 ustep: 1 
#> iter: 2  value: 1580.347 mgc: 9.07852e-05 ustep: 1 
#> iter: 3  mgc: 5.172787e-11 
#> iter: 1  mgc: 5.172787e-11 
#> outer mgc:  1.870135 
#> iter: 1  value: 1580.367 mgc: 0.1499285 ustep: 1 
#> iter: 2  value: 1580.367 mgc: 2.795871e-05 ustep: 1 
#> iter: 3  mgc: 5.34615e-12 
#> iter: 1  mgc: 5.34615e-12 
#> outer mgc:  1.030517 
#> iter: 1  value: 1580.303 mgc: 0.5330643 ustep: 1 
#> iter: 2  value: 1580.303 mgc: 0.0003052119 ustep: 1 
#> iter: 3  mgc: 7.971566e-10 
#> iter: 1  value: 1580.296 mgc: 0.5245796 ustep: 1 
#> iter: 2  value: 1580.296 mgc: 0.0002963237 ustep: 1 
#> iter: 3  mgc: 2.319005e-09 
#> iter: 1  value: 1580.335 mgc: 0.1987598 ustep: 1 
#> iter: 2  value: 1580.335 mgc: 4.354733e-05 ustep: 1 
#> iter: 3  mgc: 6.650067e-11 
#> iter: 1  mgc: 6.650067e-11 
#> outer mgc:  2.262351 
#> iter: 1  value: 1580.334 mgc: 0.1187148 ustep: 1 
#> iter: 2  value: 1580.334 mgc: 2.887213e-05 ustep: 1 
#> iter: 3  mgc: 3.909441e-11 
#> iter: 1  mgc: 3.909441e-11 
#> outer mgc:  0.956837 
#> iter: 1  value: 1580.326 mgc: 0.08190019 ustep: 1 
#> iter: 2  value: 1580.326 mgc: 6.657297e-06 ustep: 1 
#> iter: 3  mgc: 3.932933e-13 
#> iter: 1  mgc: 3.932933e-13 
#> outer mgc:  0.9015934 
#> iter: 1  value: 1580.319 mgc: 0.03575026 ustep: 1 
#> iter: 2  value: 1580.319 mgc: 9.264935e-06 ustep: 1 
#> iter: 3  mgc: 1.90936e-12 
#> iter: 1  value: 1580.315 mgc: 0.05493151 ustep: 1 
#> iter: 2  value: 1580.315 mgc: 6.815668e-06 ustep: 1 
#> iter: 3  mgc: 1.921299e-13 
#> iter: 1  value: 1580.305 mgc: 0.1317674 ustep: 1 
#> iter: 2  value: 1580.305 mgc: 3.921717e-05 ustep: 1 
#> iter: 3  mgc: 6.633895e-12 
#> iter: 1  mgc: 6.633895e-12 
#> outer mgc:  1.180261 
#> iter: 1  value: 1580.255 mgc: 0.2088403 ustep: 1 
#> iter: 2  value: 1580.255 mgc: 9.598146e-05 ustep: 1 
#> iter: 3  mgc: 5.233104e-10 
#> iter: 1  mgc: 5.233104e-10 
#> outer mgc:  0.8708525 
#> iter: 1  value: 1580.148 mgc: 0.3830329 ustep: 1 
#> iter: 2  value: 1580.148 mgc: 0.0001872869 ustep: 1 
#> iter: 3  mgc: 3.095968e-09 
#> iter: 1  value: 1580.222 mgc: 0.4616504 ustep: 1 
#> iter: 2  value: 1580.222 mgc: 0.0001936251 ustep: 1 
#> iter: 3  mgc: 2.453464e-09 
#> iter: 1  value: 1580.246 mgc: 0.08455942 ustep: 1 
#> iter: 2  value: 1580.246 mgc: 6.730657e-06 ustep: 1 
#> iter: 3  mgc: 2.421359e-12 
#> iter: 1  mgc: 2.421359e-12 
#> outer mgc:  1.048893 
#> iter: 1  value: 1580.244 mgc: 0.07082067 ustep: 1 
#> iter: 2  value: 1580.244 mgc: 8.245876e-06 ustep: 1 
#> iter: 3  mgc: 2.967073e-12 
#> iter: 1  mgc: 2.967073e-12 
#> outer mgc:  0.6298153 
#> iter: 1  value: 1580.235 mgc: 0.03511622 ustep: 1 
#> iter: 2  value: 1580.235 mgc: 1.32844e-06 ustep: 1 
#> iter: 3  mgc: 1.859147e-13 
#> iter: 1  mgc: 1.859147e-13 
#> outer mgc:  0.6371217 
#> iter: 1  value: 1580.224 mgc: 0.04827777 ustep: 1 
#> iter: 2  value: 1580.224 mgc: 6.216473e-06 ustep: 1 
#> iter: 3  mgc: 1.154058e-12 
#> iter: 1  value: 1580.212 mgc: 0.01711895 ustep: 1 
#> iter: 2  value: 1580.212 mgc: 2.490763e-06 ustep: 1 
#> iter: 3  mgc: 1.712299e-13 
#> iter: 1  value: 1580.177 mgc: 0.05242249 ustep: 1 
#> iter: 2  value: 1580.177 mgc: 2.33516e-05 ustep: 1 
#> iter: 3  mgc: 1.504824e-11 
#> iter: 1  value: 1580.122 mgc: 0.08563073 ustep: 1 
#> iter: 2  value: 1580.122 mgc: 6.231589e-05 ustep: 1 
#> iter: 3  mgc: 1.073192e-10 
#> iter: 1  mgc: 1.073192e-10 
#> outer mgc:  1.071297 
#> iter: 1  value: 1580.097 mgc: 0.1448125 ustep: 1 
#> iter: 2  value: 1580.097 mgc: 3.649301e-05 ustep: 1 
#> iter: 3  mgc: 1.698427e-11 
#> iter: 1  mgc: 1.698427e-11 
#> outer mgc:  0.560235 
#> iter: 1  value: 1580.061 mgc: 0.2749075 ustep: 1 
#> iter: 2  value: 1580.061 mgc: 5.884583e-05 ustep: 1 
#> iter: 3  mgc: 9.450403e-11 
#> iter: 1  mgc: 9.450403e-11 
#> outer mgc:  2.952874 
#> iter: 1  value: 1580.102 mgc: 0.2356333 ustep: 1 
#> iter: 2  value: 1580.102 mgc: 0.0001476499 ustep: 1 
#> iter: 3  mgc: 1.079253e-09 
#> iter: 1  mgc: 1.079253e-09 
#> outer mgc:  1.516727 
#> iter: 1  value: 1580.143 mgc: 0.1818257 ustep: 1 
#> iter: 2  value: 1580.143 mgc: 1.812636e-05 ustep: 1 
#> iter: 3  mgc: 5.697824e-12 
#> iter: 1  mgc: 5.697824e-12 
#> outer mgc:  1.331272 
#> iter: 1  value: 1580.124 mgc: 0.0954535 ustep: 1 
#> iter: 2  value: 1580.124 mgc: 1.564237e-05 ustep: 1 
#> iter: 3  mgc: 8.662789e-12 
#> iter: 1  mgc: 8.662789e-12 
#> outer mgc:  0.2846467 
#> iter: 1  value: 1580.113 mgc: 0.1376676 ustep: 1 
#> iter: 2  value: 1580.113 mgc: 2.906013e-05 ustep: 1 
#> iter: 3  mgc: 2.563122e-11 
#> iter: 1  mgc: 2.563122e-11 
#> outer mgc:  0.7032961 
#> iter: 1  value: 1580.102 mgc: 0.1112389 ustep: 1 
#> iter: 2  value: 1580.102 mgc: 8.634059e-05 ustep: 1 
#> iter: 3  mgc: 1.677362e-10 
#> iter: 1  value: 1580.11 mgc: 0.109168 ustep: 1 
#> iter: 2  value: 1580.11 mgc: 5.41335e-05 ustep: 1 
#> iter: 3  mgc: 6.676019e-11 
#> iter: 1  value: 1580.114 mgc: 0.130504 ustep: 1 
#> iter: 2  value: 1580.114 mgc: 3.748239e-05 ustep: 1 
#> iter: 3  mgc: 4.655802e-11 
#> iter: 1  value: 1580.116 mgc: 0.1485281 ustep: 1 
#> iter: 2  value: 1580.116 mgc: 3.443626e-05 ustep: 1 
#> iter: 3  mgc: 5.096689e-11 
#> iter: 1  value: 1580.114 mgc: 0.07693597 ustep: 1 
#> iter: 2  value: 1580.114 mgc: 8.883971e-06 ustep: 1 
#> iter: 3  mgc: 2.980016e-12 
#> iter: 1  mgc: 2.980016e-12 
#> outer mgc:  1.100376 
#> iter: 1  value: 1580.111 mgc: 0.0431937 ustep: 1 
#> iter: 2  value: 1580.111 mgc: 2.657208e-06 ustep: 1 
#> iter: 3  mgc: 2.59833e-13 
#> iter: 1  mgc: 2.59833e-13 
#> outer mgc:  0.4179495 
#> iter: 1  value: 1580.111 mgc: 0.02938571 ustep: 1 
#> iter: 2  value: 1580.111 mgc: 1.903769e-06 ustep: 1 
#> iter: 3  mgc: 1.363429e-13 
#> iter: 1  mgc: 1.363429e-13 
#> outer mgc:  0.5493543 
#> iter: 1  value: 1580.109 mgc: 0.03195848 ustep: 1 
#> iter: 2  value: 1580.109 mgc: 1.085969e-06 ustep: 1 
#> iter: 3  mgc: 1.552056e-13 
#> iter: 1  mgc: 1.552056e-13 
#> outer mgc:  0.3071183 
#> iter: 1  value: 1580.11 mgc: 0.01760535 ustep: 1 
#> iter: 2  value: 1580.11 mgc: 1.225784e-06 ustep: 1 
#> iter: 3  mgc: 1.642658e-13 
#> iter: 1  value: 1580.111 mgc: 0.03140519 ustep: 1 
#> iter: 2  value: 1580.111 mgc: 1.932099e-06 ustep: 1 
#> iter: 3  mgc: 1.969403e-13 
#> iter: 1  value: 1580.115 mgc: 0.06254898 ustep: 1 
#> iter: 2  value: 1580.115 mgc: 7.660905e-06 ustep: 1 
#> iter: 3  mgc: 1.587533e-12 
#> iter: 1  mgc: 1.587533e-12 
#> outer mgc:  0.4746458 
#> iter: 1  value: 1580.118 mgc: 0.07200637 ustep: 1 
#> iter: 2  value: 1580.118 mgc: 6.430464e-05 ustep: 1 
#> iter: 3  mgc: 8.89907e-11 
#> iter: 1  mgc: 8.89907e-11 
#> outer mgc:  0.300396 
#> iter: 1  value: 1580.143 mgc: 0.1446339 ustep: 1 
#> iter: 2  value: 1580.143 mgc: 7.617656e-05 ustep: 1 
#> iter: 3  mgc: 1.242536e-10 
#> iter: 1  value: 1580.128 mgc: 0.05941005 ustep: 1 
#> iter: 2  value: 1580.128 mgc: 1.596713e-05 ustep: 1 
#> iter: 3  mgc: 6.934664e-12 
#> iter: 1  mgc: 6.934664e-12 
#> outer mgc:  0.2049157 
#> iter: 1  value: 1580.129 mgc: 0.136408 ustep: 1 
#> iter: 2  value: 1580.129 mgc: 3.20367e-05 ustep: 1 
#> iter: 3  mgc: 1.295904e-11 
#> iter: 1  value: 1580.129 mgc: 0.08340344 ustep: 1 
#> iter: 2  value: 1580.129 mgc: 1.049677e-05 ustep: 1 
#> iter: 3  mgc: 2.029042e-12 
#> iter: 1  value: 1580.13 mgc: 0.06197702 ustep: 1 
#> iter: 2  value: 1580.13 mgc: 5.001659e-06 ustep: 1 
#> iter: 3  mgc: 1.143013e-12 
#> iter: 1  value: 1580.13 mgc: 0.05508082 ustep: 1 
#> iter: 2  value: 1580.13 mgc: 4.338325e-06 ustep: 1 
#> iter: 3  mgc: 8.354426e-13 
#> iter: 1  value: 1580.128 mgc: 0.01779829 ustep: 1 
#> iter: 2  value: 1580.128 mgc: 4.735365e-07 ustep: 1 
#> iter: 3  mgc: 1.598956e-13 
#> iter: 1  value: 1580.128 mgc: 0.008022748 ustep: 1 
#> iter: 2  value: 1580.128 mgc: 9.622825e-08 ustep: 1 
#> iter: 3  mgc: 1.332244e-13 
#> iter: 1  mgc: 1.332244e-13 
#> outer mgc:  0.1452781 
#> iter: 1  value: 1580.128 mgc: 0.004586401 ustep: 1 
#> iter: 2  value: 1580.128 mgc: 3.639536e-08 ustep: 1 
#> iter: 3  mgc: 1.038734e-13 
#> iter: 1  mgc: 1.038734e-13 
#> outer mgc:  0.1146081 
#> iter: 1  value: 1580.128 mgc: 0.006927566 ustep: 1 
#> iter: 2  value: 1580.128 mgc: 5.359524e-08 ustep: 1 
#> iter: 3  mgc: 1.463608e-13 
#> iter: 1  mgc: 1.463608e-13 
#> outer mgc:  0.1557494 
#> iter: 1  value: 1580.127 mgc: 0.005246627 ustep: 1 
#> iter: 2  value: 1580.127 mgc: 5.291177e-08 ustep: 1 
#> iter: 3  mgc: 1.432766e-13 
#> iter: 1  mgc: 1.432766e-13 
#> outer mgc:  0.1065086 
#> iter: 1  value: 1580.127 mgc: 0.007239521 ustep: 1 
#> iter: 2  value: 1580.127 mgc: 6.923104e-08 ustep: 1 
#> iter: 3  mgc: 1.098457e-13 
#> iter: 1  value: 1580.127 mgc: 0.002590312 ustep: 1 
#> iter: 2  value: 1580.127 mgc: 1.348415e-08 ustep: 1 
#> mgc: 1.09321e-13 
#> iter: 1  value: 1580.127 mgc: 0.01019692 ustep: 1 
#> iter: 2  value: 1580.127 mgc: 2.088894e-07 ustep: 1 
#> iter: 3  mgc: 1.297796e-13 
#> iter: 1  value: 1580.127 mgc: 0.02604183 ustep: 1 
#> iter: 2  value: 1580.127 mgc: 1.361302e-06 ustep: 1 
#> iter: 3  mgc: 1.491484e-13 
#> iter: 1  mgc: 1.491484e-13 
#> outer mgc:  0.1800129 
#> iter: 1  value: 1580.125 mgc: 0.03604156 ustep: 1 
#> iter: 2  value: 1580.125 mgc: 2.25534e-06 ustep: 1 
#> iter: 3  mgc: 1.295078e-13 
#> iter: 1  mgc: 1.295078e-13 
#> outer mgc:  0.1430783 
#> iter: 1  value: 1580.122 mgc: 0.03827742 ustep: 1 
#> iter: 2  value: 1580.122 mgc: 1.358617e-06 ustep: 1 
#> iter: 3  mgc: 1.393144e-13 
#> iter: 1  mgc: 1.393144e-13 
#> outer mgc:  0.3383863 
#> iter: 1  value: 1580.121 mgc: 0.03191204 ustep: 1 
#> iter: 2  value: 1580.121 mgc: 4.204103e-06 ustep: 1 
#> iter: 3  mgc: 3.828676e-13 
#> iter: 1  mgc: 3.828676e-13 
#> outer mgc:  0.1084597 
#> iter: 1  value: 1580.12 mgc: 0.02666544 ustep: 1 
#> iter: 2  value: 1580.12 mgc: 1.429448e-06 ustep: 1 
#> iter: 3  mgc: 1.593201e-13 
#> iter: 1  mgc: 1.593201e-13 
#> outer mgc:  0.1407966 
#> iter: 1  value: 1580.117 mgc: 0.02019797 ustep: 1 
#> iter: 2  value: 1580.117 mgc: 4.287925e-07 ustep: 1 
#> iter: 3  mgc: 1.690419e-13 
#> iter: 1  mgc: 1.690419e-13 
#> outer mgc:  0.3018873 
#> iter: 1  value: 1580.118 mgc: 0.02085543 ustep: 1 
#> iter: 2  value: 1580.118 mgc: 1.558215e-06 ustep: 1 
#> iter: 3  mgc: 1.67209e-13 
#> iter: 1  value: 1580.119 mgc: 0.0172798 ustep: 1 
#> iter: 2  value: 1580.119 mgc: 9.490809e-07 ustep: 1 
#> iter: 3  mgc: 1.347504e-13 
#> iter: 1  mgc: 1.347504e-13 
#> outer mgc:  0.1664069 
#> iter: 1  value: 1580.122 mgc: 0.03341399 ustep: 1 
#> iter: 2  value: 1580.122 mgc: 4.586278e-06 ustep: 1 
#> iter: 3  mgc: 8.073195e-13 
#> iter: 1  value: 1580.124 mgc: 0.02824103 ustep: 1 
#> iter: 2  value: 1580.124 mgc: 4.060412e-06 ustep: 1 
#> iter: 3  mgc: 4.20598e-13 
#> iter: 1  mgc: 4.20598e-13 
#> outer mgc:  0.5190698 
#> iter: 1  value: 1580.128 mgc: 0.08666156 ustep: 1 
#> iter: 2  value: 1580.128 mgc: 1.778413e-05 ustep: 1 
#> iter: 3  mgc: 3.27687e-12 
#> iter: 1  mgc: 3.27687e-12 
#> outer mgc:  0.2277815 
#> iter: 1  value: 1580.116 mgc: 0.1094803 ustep: 1 
#> iter: 2  value: 1580.116 mgc: 2.270782e-05 ustep: 1 
#> iter: 3  mgc: 1.377862e-11 
#> iter: 1  mgc: 1.377862e-11 
#> outer mgc:  0.2479595 
#> iter: 1  value: 1580.118 mgc: 0.08867878 ustep: 1 
#> iter: 2  value: 1580.118 mgc: 1.223614e-05 ustep: 1 
#> iter: 3  mgc: 4.692846e-13 
#> iter: 1  value: 1580.119 mgc: 0.1083451 ustep: 1 
#> iter: 2  value: 1580.119 mgc: 1.351154e-05 ustep: 1 
#> iter: 3  mgc: 4.868991e-13 
#> iter: 1  mgc: 4.692846e-13 
#> outer mgc:  0.309133 
#> iter: 1  value: 1580.117 mgc: 0.07687289 ustep: 1 
#> iter: 2  value: 1580.117 mgc: 4.871385e-06 ustep: 1 
#> iter: 3  mgc: 4.414217e-13 
#> iter: 1  mgc: 4.414217e-13 
#> outer mgc:  0.7629805 
#> iter: 1  value: 1580.12 mgc: 0.1459704 ustep: 1 
#> iter: 2  value: 1580.12 mgc: 2.868409e-05 ustep: 1 
#> iter: 3  mgc: 4.075153e-12 
#> iter: 1  mgc: 4.075153e-12 
#> outer mgc:  0.2104342 
#> iter: 1  value: 1580.118 mgc: 0.08848794 ustep: 1 
#> iter: 2  value: 1580.118 mgc: 1.675931e-05 ustep: 1 
#> iter: 3  mgc: 8.934328e-13 
#> iter: 1  mgc: 8.934328e-13 
#> outer mgc:  0.4325107 
#> iter: 1  value: 1580.113 mgc: 0.09110191 ustep: 1 
#> iter: 2  value: 1580.113 mgc: 1.142543e-05 ustep: 1 
#> iter: 3  mgc: 2.813987e-12 
#> iter: 1  mgc: 2.813987e-12 
#> outer mgc:  0.652926 
#> iter: 1  value: 1580.108 mgc: 0.07663279 ustep: 1 
#> iter: 2  value: 1580.108 mgc: 1.083945e-05 ustep: 1 
#> iter: 3  mgc: 4.727267e-13 
#> iter: 1  mgc: 4.727267e-13 
#> outer mgc:  0.2147305 
#> iter: 1  value: 1580.108 mgc: 0.08138902 ustep: 1 
#> iter: 2  value: 1580.108 mgc: 1.918945e-05 ustep: 1 
#> iter: 3  mgc: 1.130392e-12 
#> iter: 1  value: 1580.107 mgc: 0.09097249 ustep: 1 
#> iter: 2  value: 1580.107 mgc: 2.194487e-05 ustep: 1 
#> iter: 3  mgc: 1.613638e-12 
#> iter: 1  mgc: 1.613638e-12 
#> outer mgc:  0.5762062 
#> iter: 1  value: 1580.117 mgc: 0.09217513 ustep: 1 
#> iter: 2  value: 1580.117 mgc: 1.851655e-05 ustep: 1 
#> iter: 3  mgc: 2.847386e-12 
#> iter: 1  mgc: 2.847386e-12 
#> outer mgc:  0.2209992 
#> iter: 1  value: 1580.126 mgc: 0.05162422 ustep: 1 
#> iter: 2  value: 1580.126 mgc: 1.123302e-05 ustep: 1 
#> iter: 3  mgc: 2.526534e-12 
#> iter: 1  mgc: 2.526534e-12 
#> outer mgc:  0.4288448 
#> iter: 1  value: 1580.138 mgc: 0.06770335 ustep: 1 
#> iter: 2  value: 1580.138 mgc: 1.971487e-05 ustep: 1 
#> iter: 3  mgc: 1.690813e-12 
#> iter: 1  mgc: 1.690813e-12 
#> outer mgc:  0.4196069 
#> iter: 1  value: 1580.121 mgc: 0.1422318 ustep: 1 
#> iter: 2  value: 1580.121 mgc: 5.683362e-05 ustep: 1 
#> iter: 3  mgc: 9.928153e-12 
#> iter: 1  mgc: 9.928153e-12 
#> outer mgc:  0.208913 
#> iter: 1  value: 1580.118 mgc: 0.04532923 ustep: 1 
#> iter: 2  value: 1580.118 mgc: 5.718849e-06 ustep: 1 
#> iter: 3  mgc: 6.552926e-13 
#> iter: 1  value: 1580.113 mgc: 0.05146407 ustep: 1 
#> iter: 2  value: 1580.113 mgc: 9.86584e-06 ustep: 1 
#> iter: 3  mgc: 2.061602e-12 
#> iter: 1  mgc: 2.061602e-12 
#> outer mgc:  0.4197611 
#> iter: 1  value: 1580.093 mgc: 0.1615437 ustep: 1 
#> iter: 2  value: 1580.093 mgc: 7.235311e-05 ustep: 1 
#> iter: 3  mgc: 1.129094e-10 
#> iter: 1  mgc: 1.129094e-10 
#> outer mgc:  0.1871847 
#> iter: 1  value: 1580.107 mgc: 0.2501874 ustep: 1 
#> iter: 2  value: 1580.107 mgc: 0.0001659045 ustep: 1 
#> iter: 3  mgc: 8.450853e-11 
#> iter: 1  mgc: 8.450853e-11 
#> outer mgc:  0.2654423 
#> iter: 1  value: 1580.121 mgc: 0.1160487 ustep: 1 
#> iter: 2  value: 1580.121 mgc: 2.384773e-05 ustep: 1 
#> iter: 3  mgc: 1.116507e-11 
#> iter: 1  mgc: 1.116507e-11 
#> outer mgc:  0.495943 
#> iter: 1  value: 1580.141 mgc: 0.1839296 ustep: 1 
#> iter: 2  value: 1580.141 mgc: 1.898498e-05 ustep: 1 
#> iter: 3  mgc: 3.748328e-12 
#> iter: 1  mgc: 3.748328e-12 
#> outer mgc:  0.4685998 
#> iter: 1  value: 1580.154 mgc: 0.1306061 ustep: 1 
#> iter: 2  value: 1580.154 mgc: 6.423664e-05 ustep: 1 
#> iter: 3  mgc: 2.058834e-11 
#> iter: 1  mgc: 2.058834e-11 
#> outer mgc:  0.5912539 
#> iter: 1  value: 1580.134 mgc: 0.1800989 ustep: 1 
#> iter: 2  value: 1580.134 mgc: 4.885012e-05 ustep: 1 
#> iter: 3  mgc: 2.032687e-11 
#> iter: 1  value: 1580.142 mgc: 0.07129228 ustep: 1 
#> iter: 2  value: 1580.142 mgc: 1.40923e-05 ustep: 1 
#> iter: 3  mgc: 1.457186e-12 
#> iter: 1  mgc: 1.457186e-12 
#> outer mgc:  0.4917323 
#> iter: 1  value: 1580.121 mgc: 0.1281226 ustep: 1 
#> iter: 2  value: 1580.121 mgc: 9.682785e-06 ustep: 1 
#> iter: 3  mgc: 5.651639e-13 
#> iter: 1  mgc: 5.651639e-13 
#> outer mgc:  0.0992252 
#> iter: 1  value: 1580.111 mgc: 0.08843799 ustep: 1 
#> iter: 2  value: 1580.111 mgc: 1.718361e-05 ustep: 1 
#> iter: 3  mgc: 2.001721e-12 
#> iter: 1  mgc: 2.001721e-12 
#> outer mgc:  0.1268137 
#> iter: 1  value: 1580.11 mgc: 0.05709529 ustep: 1 
#> iter: 2  value: 1580.11 mgc: 1.195725e-05 ustep: 1 
#> iter: 3  mgc: 6.946241e-13 
#> iter: 1  mgc: 6.946241e-13 
#> outer mgc:  0.38068 
#> iter: 1  value: 1580.111 mgc: 0.03267004 ustep: 1 
#> iter: 2  value: 1580.111 mgc: 2.503729e-06 ustep: 1 
#> iter: 3  mgc: 1.281914e-13 
#> iter: 1  mgc: 1.281914e-13 
#> outer mgc:  0.07680508 
#> iter: 1  value: 1580.107 mgc: 0.04420505 ustep: 1 
#> iter: 2  value: 1580.107 mgc: 2.028416e-06 ustep: 1 
#> iter: 3  mgc: 1.477889e-13 
#> iter: 1  mgc: 1.477889e-13 
#> outer mgc:  0.15989 
#> iter: 1  value: 1580.106 mgc: 0.05462595 ustep: 1 
#> iter: 2  value: 1580.106 mgc: 1.742004e-06 ustep: 1 
#> iter: 3  mgc: 1.117482e-13 
#> iter: 1  mgc: 1.117482e-13 
#> outer mgc:  0.08417672 
#> iter: 1  value: 1580.104 mgc: 0.01798351 ustep: 1 
#> iter: 2  value: 1580.104 mgc: 6.794082e-07 ustep: 1 
#> iter: 3  mgc: 1.517278e-13 
#> iter: 1  value: 1580.099 mgc: 0.03945377 ustep: 1 
#> iter: 2  value: 1580.099 mgc: 3.731856e-06 ustep: 1 
#> iter: 3  mgc: 3.232068e-13 
#> iter: 1  mgc: 3.232068e-13 
#> outer mgc:  0.2372797 
#> iter: 1  value: 1580.108 mgc: 0.07383616 ustep: 1 
#> iter: 2  value: 1580.108 mgc: 6.639549e-06 ustep: 1 
#> iter: 3  mgc: 3.913897e-13 
#> iter: 1  mgc: 3.913897e-13 
#> outer mgc:  0.4061533 
#> iter: 1  value: 1580.107 mgc: 0.04677592 ustep: 1 
#> iter: 2  value: 1580.107 mgc: 4.129928e-06 ustep: 1 
#> iter: 3  mgc: 1.793541e-13 
#> iter: 1  mgc: 1.793541e-13 
#> outer mgc:  0.08787817 
#> iter: 1  value: 1580.11 mgc: 0.09382731 ustep: 1 
#> iter: 2  value: 1580.11 mgc: 5.999135e-06 ustep: 1 
#> iter: 3  mgc: 9.231791e-13 
#> iter: 1  value: 1580.108 mgc: 0.03709572 ustep: 1 
#> iter: 2  value: 1580.108 mgc: 9.999443e-07 ustep: 1 
#> iter: 3  mgc: 1.761733e-13 
#> iter: 1  mgc: 1.761733e-13 
#> outer mgc:  0.09341798 
#> iter: 1  value: 1580.107 mgc: 0.03300886 ustep: 1 
#> iter: 2  value: 1580.107 mgc: 9.306267e-07 ustep: 1 
#> iter: 3  mgc: 1.308041e-13 
#> iter: 1  mgc: 1.308041e-13 
#> outer mgc:  0.1128123 
#> iter: 1  value: 1580.11 mgc: 0.03672723 ustep: 1 
#> iter: 2  value: 1580.11 mgc: 1.268418e-06 ustep: 1 
#> iter: 3  mgc: 1.260354e-13 
#> iter: 1  mgc: 1.260354e-13 
#> outer mgc:  0.1165177 
#> iter: 1  value: 1580.108 mgc: 0.02294593 ustep: 1 
#> iter: 2  value: 1580.108 mgc: 7.179777e-07 ustep: 1 
#> iter: 3  mgc: 1.347408e-13 
#> iter: 1  mgc: 1.347408e-13 
#> outer mgc:  0.06702889 
#> iter: 1  value: 1580.108 mgc: 0.02399448 ustep: 1 
#> iter: 2  value: 1580.108 mgc: 1.277649e-06 ustep: 1 
#> iter: 3  mgc: 1.201529e-13 
#> iter: 1  mgc: 1.201529e-13 
#> outer mgc:  0.04540589 
#> iter: 1  value: 1580.106 mgc: 0.04127871 ustep: 1 
#> iter: 2  value: 1580.106 mgc: 1.085544e-06 ustep: 1 
#> iter: 3  mgc: 9.611452e-14 
#> iter: 1  mgc: 9.611452e-14 
#> outer mgc:  0.08499193 
#> iter: 1  value: 1580.108 mgc: 0.01536846 ustep: 1 
#> iter: 2  value: 1580.108 mgc: 5.375506e-07 ustep: 1 
#> iter: 3  mgc: 1.39833e-13 
#> iter: 1  mgc: 1.39833e-13 
#> outer mgc:  0.05867122 
#> iter: 1  value: 1580.108 mgc: 0.02138386 ustep: 1 
#> iter: 2  value: 1580.108 mgc: 9.389077e-07 ustep: 1 
#> iter: 3  mgc: 1.744977e-13 
#> iter: 1  mgc: 1.744977e-13 
#> outer mgc:  0.09052023 
#> iter: 1  value: 1580.107 mgc: 0.01280654 ustep: 1 
#> iter: 2  value: 1580.107 mgc: 4.515321e-07 ustep: 1 
#> iter: 3  mgc: 1.086331e-13 
#> iter: 1  mgc: 1.086331e-13 
#> outer mgc:  0.0523974 
#> iter: 1  value: 1580.109 mgc: 0.01519804 ustep: 1 
#> iter: 2  value: 1580.109 mgc: 4.684627e-07 ustep: 1 
#> iter: 3  mgc: 1.022028e-13 
#> iter: 1  mgc: 1.022028e-13 
#> outer mgc:  0.04848837 
#> iter: 1  value: 1580.109 mgc: 0.01144687 ustep: 1 
#> iter: 2  value: 1580.109 mgc: 2.814178e-07 ustep: 1 
#> iter: 3  mgc: 9.959124e-14 
#> iter: 1  mgc: 9.959124e-14 
#> outer mgc:  0.03083772 
#> iter: 1  value: 1580.11 mgc: 0.01831356 ustep: 1 
#> iter: 2  value: 1580.11 mgc: 2.271162e-07 ustep: 1 
#> iter: 3  mgc: 1.357939e-13 
#> iter: 1  mgc: 1.357939e-13 
#> outer mgc:  0.06546044 
#> iter: 1  value: 1580.109 mgc: 0.01061347 ustep: 1 
#> iter: 2  value: 1580.109 mgc: 2.28234e-07 ustep: 1 
#> iter: 3  mgc: 1.482857e-13 
#> iter: 1  mgc: 1.482857e-13 
#> outer mgc:  0.04387878 
#> iter: 1  value: 1580.108 mgc: 0.02389023 ustep: 1 
#> iter: 2  value: 1580.108 mgc: 6.703349e-07 ustep: 1 
#> iter: 3  mgc: 1.353592e-13 
#> iter: 1  mgc: 1.353592e-13 
#> outer mgc:  0.01813657 
#> iter: 1  value: 1580.108 mgc: 0.02349521 ustep: 1 
#> iter: 2  value: 1580.108 mgc: 1.486641e-06 ustep: 1 
#> iter: 3  mgc: 1.619054e-13 
#> iter: 1  mgc: 1.619054e-13 
#> outer mgc:  0.04020329 
#> iter: 1  value: 1580.109 mgc: 0.01670978 ustep: 1 
#> iter: 2  value: 1580.109 mgc: 7.556865e-07 ustep: 1 
#> iter: 3  mgc: 1.791363e-13 
#> iter: 1  value: 1580.108 mgc: 0.00877017 ustep: 1 
#> iter: 2  value: 1580.108 mgc: 1.98588e-07 ustep: 1 
#> iter: 3  mgc: 1.366059e-13 
#> iter: 1  mgc: 1.366059e-13 
#> outer mgc:  0.03414236 
#> iter: 1  value: 1580.107 mgc: 0.03437256 ustep: 1 
#> iter: 2  value: 1580.107 mgc: 3.106294e-06 ustep: 1 
#> iter: 3  mgc: 1.897745e-13 
#> iter: 1  value: 1580.107 mgc: 0.01691201 ustep: 1 
#> iter: 2  value: 1580.107 mgc: 7.521124e-07 ustep: 1 
#> iter: 3  mgc: 1.768219e-13 
#> iter: 1  mgc: 1.768219e-13 
#> outer mgc:  0.01963737 
#> iter: 1  value: 1580.109 mgc: 0.0171462 ustep: 1 
#> iter: 2  value: 1580.109 mgc: 8.237565e-07 ustep: 1 
#> iter: 3  mgc: 1.116406e-13 
#> iter: 1  value: 1580.108 mgc: 0.005821502 ustep: 1 
#> iter: 2  value: 1580.108 mgc: 9.590387e-08 ustep: 1 
#> mgc: 1.386856e-13 
#> iter: 1  mgc: 1.386856e-13 
#> outer mgc:  0.01818918 
#> iter: 1  value: 1580.108 mgc: 0.004537942 ustep: 1 
#> iter: 2  value: 1580.108 mgc: 6.115984e-08 ustep: 1 
#> mgc: 1.403777e-13 
#> iter: 1  mgc: 1.403777e-13 
#> outer mgc:  0.02367428 
#> iter: 1  value: 1580.108 mgc: 0.002179143 ustep: 1 
#> iter: 2  value: 1580.108 mgc: 1.090601e-08 ustep: 1 
#> mgc: 1.391501e-13 
#> iter: 1  mgc: 1.391501e-13 
#> outer mgc:  0.01266357 
#> iter: 1  value: 1580.108 mgc: 0.001700113 ustep: 1 
#> iter: 2  mgc: 4.358361e-09 
#> iter: 1  mgc: 4.358361e-09 
#> outer mgc:  0.01952649 
#> iter: 1  value: 1580.108 mgc: 0.003384432 ustep: 1 
#> iter: 2  value: 1580.108 mgc: 1.469984e-08 ustep: 1 
#> mgc: 1.534974e-13 
#> iter: 1  mgc: 1.534974e-13 
#> outer mgc:  0.02812634 
#> iter: 1  value: 1580.108 mgc: 0.002051229 ustep: 1 
#> iter: 2  value: 1580.108 mgc: 1.107803e-08 ustep: 1 
#> mgc: 1.322365e-13 
#> iter: 1  value: 1580.108 mgc: 0.001645037 ustep: 1 
#> iter: 2  mgc: 7.999017e-09 
#> iter: 1  mgc: 7.999017e-09 
#> outer mgc:  0.01040119 
#> iter: 1  value: 1580.107 mgc: 0.01091624 ustep: 1 
#> iter: 2  value: 1580.107 mgc: 1.196929e-07 ustep: 1 
#> iter: 3  mgc: 1.368561e-13 
#> iter: 1  value: 1580.108 mgc: 0.003166503 ustep: 1 
#> iter: 2  mgc: 9.829352e-09 
#> iter: 1  mgc: 9.829352e-09 
#> outer mgc:  0.01166554 
#> iter: 1  value: 1580.108 mgc: 0.004481115 ustep: 1 
#> iter: 2  value: 1580.108 mgc: 3.362306e-08 ustep: 1 
#> mgc: 1.995958e-13 
#> iter: 1  mgc: 1.995958e-13 
#> outer mgc:  0.00525045 
#> iter: 1  value: 1580.108 mgc: 0.00392934 ustep: 1 
#> iter: 2  value: 1580.108 mgc: 3.539257e-08 ustep: 1 
#> mgc: 2.193099e-13 
#> iter: 1  mgc: 2.193099e-13 
#> outer mgc:  0.003541607 
#> iter: 1  value: 1580.108 mgc: 0.002005092 ustep: 1 
#> iter: 2  value: 1580.108 mgc: 1.7652e-08 ustep: 1 
#> mgc: 1.276103e-13 
#> iter: 1  mgc: 1.276103e-13 
#> outer mgc:  0.003233384 
#> iter: 1  value: 1580.108 mgc: 0.001705586 ustep: 1 
#> iter: 2  mgc: 4.03167e-09 
#> iter: 1  mgc: 4.03167e-09 
#> outer mgc:  0.003158307 
#> iter: 1  value: 1580.108 mgc: 0.001899892 ustep: 1 
#> iter: 2  mgc: 2.636207e-09 
#> iter: 1  mgc: 2.636207e-09 
#> outer mgc:  0.005800196 
#> iter: 1  value: 1580.108 mgc: 0.0008511993 ustep: 1 
#> iter: 2  mgc: 1.203014e-09 
#> iter: 1  mgc: 1.203014e-09 
#> outer mgc:  0.001190837 
#> iter: 1  value: 1580.108 mgc: 0.0003675151 ustep: 1 
#> iter: 2  mgc: 1.798925e-10 
#> iter: 1  mgc: 1.798925e-10 
#> outer mgc:  0.0003962758 
#> iter: 1  mgc: 1.798925e-10 
#> iter: 1  mgc: 1.798925e-10 
#> outer mgc:  0.0003962758 
#> iter: 1  mgc: 1.798925e-10 
#> outer mgc:  0.0003962758 
#> Retro Peel: 6
#> Using stats::nlminb for optimization with these control parameters:
#> iter.max: 1000
#> eval.max: 1000
#> iter: 1  value: 1548.483 mgc: 1.563108 ustep: 1 
#> iter: 2  value: 1548.483 mgc: 0.012741 ustep: 1 
#> iter: 3  value: 1548.483 mgc: 3.836886e-06 ustep: 1 
#> iter: 4  mgc: 1.459343e-12 
#> iter: 1  mgc: 1.459343e-12 
#> Matching hessian patterns... Done
#> outer mgc:  17.33372 
#> iter: 1  value: 4561.474 mgc: 132.4805 ustep: 0.01407033 
#> iter: 2  value: 3816.558 mgc: 58.70817 ustep: 0.06122157 
#> iter: 3  value: 3486.655 mgc: 123.075 ustep: 0.2475052 
#> iter: 4  value: 3447.017 mgc: 32.46288 ustep: 0.4975492 
#> iter: 5  value: 3443.696 mgc: 6.571332 ustep: 0.7054011 
#> iter: 6  value: 3443.578 mgc: 0.6721669 ustep: 0.8398976 
#> iter: 7  value: 3443.576 mgc: 0.03707728 ustep: 0.9164676 
#> iter: 8  value: 3443.576 mgc: 0.003056224 ustep: 0.9573274 
#> iter: 9  value: 3443.576 mgc: 0.0001310532 ustep: 0.9784333 
#> iter: 10  value: 3443.576 mgc: 2.813787e-06 ustep: 0.9891589 
#> iter: 11  value: 3443.576 mgc: 3.014102e-08 ustep: 0.9945652 
#> iter: 12  mgc: 1.613185e-10 
#> iter: 1  value: 1564.876 mgc: 12.16936 ustep: 1 
#> iter: 2  value: 1564.847 mgc: 0.2970048 ustep: 1 
#> iter: 3  value: 1564.847 mgc: 0.005455828 ustep: 1 
#> iter: 4  value: 1564.847 mgc: 1.745359e-06 ustep: 1 
#> iter: 5  mgc: 1.968593e-13 
#> iter: 1  value: 1548.304 mgc: 1.19459 ustep: 1 
#> iter: 2  value: 1548.304 mgc: 0.002409817 ustep: 1 
#> iter: 3  value: 1548.304 mgc: 2.898754e-07 ustep: 1 
#> iter: 4  mgc: 1.619888e-13 
#> iter: 1  mgc: 1.619888e-13 
#> outer mgc:  17.3975 
#> iter: 1  value: 1548.053 mgc: 1.105617 ustep: 1 
#> iter: 2  value: 1548.053 mgc: 0.001794776 ustep: 1 
#> iter: 3  value: 1548.053 mgc: 1.641455e-07 ustep: 1 
#> iter: 4  mgc: 8.872583e-14 
#> iter: 1  mgc: 8.872583e-14 
#> outer mgc:  15.85452 
#> iter: 1  value: 1547.818 mgc: 0.8481612 ustep: 1 
#> iter: 2  value: 1547.818 mgc: 0.001268022 ustep: 1 
#> iter: 3  value: 1547.818 mgc: 8.094604e-08 ustep: 1 
#> iter: 4  mgc: 1.44908e-13 
#> iter: 1  mgc: 1.44908e-13 
#> outer mgc:  15.24269 
#> iter: 1  value: 1547.507 mgc: 0.8363322 ustep: 1 
#> iter: 2  value: 1547.507 mgc: 0.0009176808 ustep: 1 
#> iter: 3  value: 1547.507 mgc: 4.850952e-08 ustep: 1 
#> iter: 4  mgc: 1.894963e-13 
#> iter: 1  mgc: 1.894963e-13 
#> outer mgc:  13.62635 
#> iter: 1  value: 1547.082 mgc: 0.6588357 ustep: 1 
#> iter: 2  value: 1547.082 mgc: 0.001146349 ustep: 1 
#> iter: 3  value: 1547.082 mgc: 8.215218e-08 ustep: 1 
#> iter: 4  mgc: 1.299079e-13 
#> iter: 1  value: 1546.732 mgc: 0.2797587 ustep: 1 
#> iter: 2  value: 1546.732 mgc: 7.190506e-05 ustep: 1 
#> iter: 3  mgc: 9.772922e-11 
#> iter: 1  value: 1545.967 mgc: 0.8367649 ustep: 1 
#> iter: 2  value: 1545.967 mgc: 0.0006331299 ustep: 1 
#> iter: 3  mgc: 7.564182e-09 
#> iter: 1  mgc: 7.564182e-09 
#> outer mgc:  14.53304 
#> iter: 1  value: 1547.381 mgc: 4.761274 ustep: 1 
#> iter: 2  value: 1547.381 mgc: 0.02494173 ustep: 1 
#> iter: 3  value: 1547.381 mgc: 2.58295e-05 ustep: 1 
#> iter: 4  mgc: 7.003572e-11 
#> iter: 1  value: 1546.12 mgc: 2.853942 ustep: 1 
#> iter: 2  value: 1546.12 mgc: 0.005979132 ustep: 1 
#> iter: 3  value: 1546.12 mgc: 2.218863e-06 ustep: 1 
#> iter: 4  mgc: 8.607102e-13 
#> iter: 1  mgc: 8.607102e-13 
#> outer mgc:  38.66774 
#> iter: 1  value: 1545.638 mgc: 2.854518 ustep: 1 
#> iter: 2  value: 1545.638 mgc: 0.01405247 ustep: 1 
#> iter: 3  value: 1545.638 mgc: 9.777827e-06 ustep: 1 
#> iter: 4  mgc: 5.674594e-12 
#> iter: 1  mgc: 5.674594e-12 
#> outer mgc:  23.09286 
#> iter: 1  value: 1544.996 mgc: 2.094416 ustep: 1 
#> iter: 2  value: 1544.996 mgc: 0.007070864 ustep: 1 
#> iter: 3  value: 1544.996 mgc: 1.716261e-06 ustep: 1 
#> iter: 4  mgc: 1.825115e-13 
#> iter: 1  mgc: 1.825115e-13 
#> outer mgc:  6.295243 
#> iter: 1  value: 1544.705 mgc: 0.7861911 ustep: 1 
#> iter: 2  value: 1544.705 mgc: 0.001066463 ustep: 1 
#> iter: 3  value: 1544.705 mgc: 7.01046e-08 ustep: 1 
#> iter: 4  mgc: 2.289732e-13 
#> iter: 1  mgc: 2.289732e-13 
#> outer mgc:  9.400908 
#> iter: 1  value: 1544.384 mgc: 0.9113203 ustep: 1 
#> iter: 2  value: 1544.384 mgc: 0.000978636 ustep: 1 
#> iter: 3  value: 1544.384 mgc: 6.491669e-08 ustep: 1 
#> iter: 4  mgc: 1.16602e-13 
#> iter: 1  mgc: 1.16602e-13 
#> outer mgc:  11.35084 
#> iter: 1  value: 1543.912 mgc: 2.224688 ustep: 1 
#> iter: 2  value: 1543.912 mgc: 0.002926173 ustep: 1 
#> iter: 3  value: 1543.912 mgc: 6.214669e-07 ustep: 1 
#> iter: 4  mgc: 1.248998e-13 
#> iter: 1  mgc: 1.248998e-13 
#> outer mgc:  9.453118 
#> iter: 1  value: 1544.614 mgc: 2.988266 ustep: 1 
#> iter: 2  value: 1544.614 mgc: 0.01730126 ustep: 1 
#> iter: 3  value: 1544.614 mgc: 8.178286e-06 ustep: 1 
#> iter: 4  mgc: 5.01882e-12 
#> iter: 1  value: 1544.293 mgc: 2.193008 ustep: 1 
#> iter: 2  value: 1544.293 mgc: 0.008141927 ustep: 1 
#> iter: 3  value: 1544.293 mgc: 2.520773e-06 ustep: 1 
#> iter: 4  mgc: 5.292834e-13 
#> iter: 1  value: 1544.247 mgc: 1.866456 ustep: 1 
#> iter: 2  value: 1544.247 mgc: 0.005378018 ustep: 1 
#> iter: 3  value: 1544.247 mgc: 1.12021e-06 ustep: 1 
#> iter: 4  mgc: 1.60264e-13 
#> iter: 1  value: 1543.89 mgc: 0.5168758 ustep: 1 
#> iter: 2  value: 1543.89 mgc: 0.0004201168 ustep: 1 
#> iter: 3  mgc: 6.98391e-09 
#> iter: 1  mgc: 6.98391e-09 
#> outer mgc:  5.890071 
#> iter: 1  value: 1543.853 mgc: 0.2151853 ustep: 1 
#> iter: 2  value: 1543.853 mgc: 7.220754e-05 ustep: 1 
#> iter: 3  mgc: 2.480256e-10 
#> iter: 1  mgc: 2.480256e-10 
#> outer mgc:  5.882992 
#> iter: 1  value: 1543.845 mgc: 0.4910437 ustep: 1 
#> iter: 2  value: 1543.845 mgc: 0.0003886744 ustep: 1 
#> iter: 3  mgc: 6.511461e-09 
#> iter: 1  mgc: 6.511461e-09 
#> outer mgc:  8.70266 
#> iter: 1  value: 1543.769 mgc: 0.5908697 ustep: 1 
#> iter: 2  value: 1543.769 mgc: 0.0004735824 ustep: 1 
#> iter: 3  value: 1543.769 mgc: 1.393404e-08 ustep: 1 
#> mgc: 1.780884e-13 
#> iter: 1  mgc: 1.780884e-13 
#> outer mgc:  5.715404 
#> iter: 1  value: 1543.7 mgc: 0.3537391 ustep: 1 
#> iter: 2  value: 1543.7 mgc: 0.0001400385 ustep: 1 
#> iter: 3  mgc: 6.412113e-10 
#> iter: 1  value: 1543.642 mgc: 0.3918042 ustep: 1 
#> iter: 2  value: 1543.642 mgc: 0.0001878217 ustep: 1 
#> iter: 3  mgc: 1.936419e-09 
#> iter: 1  mgc: 1.936419e-09 
#> outer mgc:  10.18376 
#> iter: 1  value: 1543.04 mgc: 1.19637 ustep: 1 
#> iter: 2  value: 1543.04 mgc: 0.001554013 ustep: 1 
#> iter: 3  value: 1543.04 mgc: 8.234016e-08 ustep: 1 
#> iter: 4  mgc: 1.250593e-13 
#> iter: 1  mgc: 1.250593e-13 
#> outer mgc:  7.857989 
#> iter: 1  value: 1542.929 mgc: 1.693001 ustep: 1 
#> iter: 2  value: 1542.929 mgc: 0.005953938 ustep: 1 
#> iter: 3  value: 1542.929 mgc: 1.091904e-06 ustep: 1 
#> iter: 4  mgc: 1.147423e-13 
#> iter: 1  value: 1543.028 mgc: 1.393166 ustep: 1 
#> iter: 2  value: 1543.028 mgc: 0.003642538 ustep: 1 
#> iter: 3  value: 1543.028 mgc: 3.963106e-07 ustep: 1 
#> iter: 4  mgc: 1.821248e-13 
#> iter: 1  value: 1543.098 mgc: 1.407971 ustep: 1 
#> iter: 2  value: 1543.098 mgc: 0.003283819 ustep: 1 
#> iter: 3  value: 1543.098 mgc: 3.683004e-07 ustep: 1 
#> iter: 4  mgc: 1.407835e-13 
#> iter: 1  value: 1543.087 mgc: 1.162865 ustep: 1 
#> iter: 2  value: 1543.087 mgc: 0.002062875 ustep: 1 
#> iter: 3  value: 1543.087 mgc: 1.610017e-07 ustep: 1 
#> iter: 4  mgc: 1.09708e-13 
#> iter: 1  value: 1543.014 mgc: 0.559916 ustep: 1 
#> iter: 2  value: 1543.014 mgc: 0.0004821674 ustep: 1 
#> iter: 3  mgc: 8.8774e-09 
#> iter: 1  mgc: 8.8774e-09 
#> outer mgc:  5.049547 
#> iter: 1  value: 1542.98 mgc: 0.2897143 ustep: 1 
#> iter: 2  value: 1542.98 mgc: 6.845845e-05 ustep: 1 
#> iter: 3  mgc: 2.506641e-10 
#> iter: 1  mgc: 2.506641e-10 
#> outer mgc:  5.037527 
#> iter: 1  value: 1542.958 mgc: 0.4966308 ustep: 1 
#> iter: 2  value: 1542.958 mgc: 0.000403967 ustep: 1 
#> iter: 3  mgc: 5.617394e-09 
#> iter: 1  mgc: 5.617394e-09 
#> outer mgc:  6.854393 
#> iter: 1  value: 1542.831 mgc: 0.4265532 ustep: 1 
#> iter: 2  value: 1542.831 mgc: 0.0002150177 ustep: 1 
#> iter: 3  mgc: 2.601445e-09 
#> iter: 1  mgc: 2.601445e-09 
#> outer mgc:  5.053152 
#> iter: 1  value: 1542.682 mgc: 0.2557715 ustep: 1 
#> iter: 2  value: 1542.682 mgc: 0.0002581003 ustep: 1 
#> iter: 3  mgc: 8.419733e-10 
#> iter: 1  value: 1542.447 mgc: 0.3973181 ustep: 1 
#> iter: 2  value: 1542.447 mgc: 0.0004195107 ustep: 1 
#> iter: 3  mgc: 5.523835e-10 
#> iter: 1  value: 1542.019 mgc: 0.8896398 ustep: 1 
#> iter: 2  value: 1542.019 mgc: 0.002110112 ustep: 1 
#> iter: 3  value: 1542.019 mgc: 1.39583e-08 ustep: 1 
#> mgc: 1.228702e-13 
#> iter: 1  mgc: 1.228702e-13 
#> outer mgc:  13.28134 
#> iter: 1  value: 1541.817 mgc: 2.44871 ustep: 1 
#> iter: 2  value: 1541.817 mgc: 0.006422566 ustep: 1 
#> iter: 3  value: 1541.817 mgc: 3.963159e-07 ustep: 1 
#> iter: 4  mgc: 2.121552e-13 
#> iter: 1  mgc: 2.121552e-13 
#> outer mgc:  10.51303 
#> iter: 1  value: 1542.038 mgc: 2.061539 ustep: 1 
#> iter: 2  value: 1542.038 mgc: 0.0065995 ustep: 1 
#> iter: 3  value: 1542.038 mgc: 1.537606e-06 ustep: 1 
#> iter: 4  mgc: 3.22381e-13 
#> iter: 1  mgc: 3.22381e-13 
#> outer mgc:  8.070833 
#> iter: 1  value: 1542.454 mgc: 2.665668 ustep: 1 
#> iter: 2  value: 1542.454 mgc: 0.01072131 ustep: 1 
#> iter: 3  value: 1542.454 mgc: 2.811615e-06 ustep: 1 
#> iter: 4  mgc: 6.394698e-13 
#> iter: 1  value: 1542.132 mgc: 1.840618 ustep: 1 
#> iter: 2  value: 1542.132 mgc: 0.004056897 ustep: 1 
#> iter: 3  value: 1542.132 mgc: 1.342306e-06 ustep: 1 
#> iter: 4  mgc: 1.789356e-13 
#> iter: 1  value: 1542.278 mgc: 2.023595 ustep: 1 
#> iter: 2  value: 1542.278 mgc: 0.00726135 ustep: 1 
#> iter: 3  value: 1542.278 mgc: 2.052981e-06 ustep: 1 
#> iter: 4  mgc: 2.629467e-13 
#> iter: 1  value: 1541.993 mgc: 0.7607043 ustep: 1 
#> iter: 2  value: 1541.993 mgc: 0.0009336448 ustep: 1 
#> iter: 3  value: 1541.993 mgc: 3.360016e-08 ustep: 1 
#> mgc: 1.47092e-13 
#> iter: 1  mgc: 1.47092e-13 
#> outer mgc:  10.61675 
#> iter: 1  value: 1541.963 mgc: 0.5717887 ustep: 1 
#> iter: 2  value: 1541.963 mgc: 0.0004441852 ustep: 1 
#> iter: 3  mgc: 7.453063e-09 
#> iter: 1  mgc: 7.453063e-09 
#> outer mgc:  6.594834 
#> iter: 1  value: 1541.891 mgc: 0.442469 ustep: 1 
#> iter: 2  value: 1541.891 mgc: 0.0002402233 ustep: 1 
#> iter: 3  mgc: 3.919816e-09 
#> iter: 1  mgc: 3.919816e-09 
#> outer mgc:  6.621681 
#> iter: 1  value: 1541.8 mgc: 0.6208632 ustep: 1 
#> iter: 2  value: 1541.8 mgc: 0.0006112925 ustep: 1 
#> iter: 3  value: 1541.8 mgc: 1.096506e-08 ustep: 1 
#> mgc: 1.256714e-13 
#> iter: 1  mgc: 1.256714e-13 
#> outer mgc:  5.859981 
#> iter: 1  value: 1541.761 mgc: 0.626982 ustep: 1 
#> iter: 2  value: 1541.761 mgc: 0.0006252987 ustep: 1 
#> iter: 3  mgc: 9.363056e-09 
#> iter: 1  mgc: 9.363056e-09 
#> outer mgc:  5.519522 
#> iter: 1  value: 1541.763 mgc: 0.7044041 ustep: 1 
#> iter: 2  value: 1541.763 mgc: 0.0003933938 ustep: 1 
#> iter: 3  mgc: 5.007515e-09 
#> iter: 1  mgc: 5.007515e-09 
#> outer mgc:  5.330751 
#> iter: 1  value: 1541.648 mgc: 0.8418288 ustep: 1 
#> iter: 2  value: 1541.648 mgc: 0.0007894413 ustep: 1 
#> iter: 3  value: 1541.648 mgc: 4.204374e-08 ustep: 1 
#> iter: 4  mgc: 1.779642e-13 
#> iter: 1  mgc: 1.779642e-13 
#> outer mgc:  8.558844 
#> iter: 1  value: 1541.553 mgc: 1.193669 ustep: 1 
#> iter: 2  value: 1541.553 mgc: 0.002586078 ustep: 1 
#> iter: 3  value: 1541.553 mgc: 3.312925e-07 ustep: 1 
#> iter: 4  mgc: 1.51926e-13 
#> iter: 1  mgc: 1.51926e-13 
#> outer mgc:  15.56459 
#> iter: 1  value: 1541.222 mgc: 0.9396819 ustep: 1 
#> iter: 2  value: 1541.222 mgc: 0.001574233 ustep: 1 
#> iter: 3  value: 1541.222 mgc: 6.246101e-08 ustep: 1 
#> iter: 4  mgc: 1.326593e-13 
#> iter: 1  mgc: 1.326593e-13 
#> outer mgc:  4.427389 
#> iter: 1  value: 1541.636 mgc: 0.6438124 ustep: 1 
#> iter: 2  value: 1541.636 mgc: 0.0008841652 ustep: 1 
#> iter: 3  value: 1541.636 mgc: 3.334433e-08 ustep: 1 
#> iter: 4  mgc: 1.489489e-13 
#> iter: 1  mgc: 1.489489e-13 
#> outer mgc:  3.749082 
#> iter: 1  value: 1541.284 mgc: 0.5334189 ustep: 1 
#> iter: 2  value: 1541.284 mgc: 0.0005272438 ustep: 1 
#> iter: 3  mgc: 4.335469e-09 
#> iter: 1  mgc: 4.335469e-09 
#> outer mgc:  5.718274 
#> iter: 1  value: 1541.229 mgc: 0.9040974 ustep: 1 
#> iter: 2  value: 1541.229 mgc: 0.001725324 ustep: 1 
#> iter: 3  value: 1541.229 mgc: 8.373486e-08 ustep: 1 
#> iter: 4  mgc: 1.274574e-13 
#> iter: 1  mgc: 1.274574e-13 
#> outer mgc:  8.356719 
#> iter: 1  value: 1541.281 mgc: 0.9273492 ustep: 1 
#> iter: 2  value: 1541.281 mgc: 0.0006586901 ustep: 1 
#> iter: 3  value: 1541.281 mgc: 1.605265e-08 ustep: 1 
#> mgc: 1.159508e-13 
#> iter: 1  mgc: 1.159508e-13 
#> outer mgc:  7.837379 
#> iter: 1  value: 1541.219 mgc: 0.9203226 ustep: 1 
#> iter: 2  value: 1541.219 mgc: 0.0009292661 ustep: 1 
#> iter: 3  value: 1541.219 mgc: 4.186882e-08 ustep: 1 
#> iter: 4  mgc: 2.417634e-13 
#> iter: 1  mgc: 2.417634e-13 
#> outer mgc:  4.580925 
#> iter: 1  value: 1541.13 mgc: 0.8616585 ustep: 1 
#> iter: 2  value: 1541.13 mgc: 0.00078426 ustep: 1 
#> iter: 3  mgc: 2.372658e-09 
#> iter: 1  mgc: 2.372658e-09 
#> outer mgc:  4.792364 
#> iter: 1  value: 1541.218 mgc: 0.5450891 ustep: 1 
#> iter: 2  value: 1541.218 mgc: 0.0007223814 ustep: 1 
#> iter: 3  mgc: 5.350207e-09 
#> iter: 1  mgc: 5.350207e-09 
#> outer mgc:  11.00764 
#> iter: 1  value: 1540.857 mgc: 0.73258 ustep: 1 
#> iter: 2  value: 1540.857 mgc: 0.0007631543 ustep: 1 
#> iter: 3  value: 1540.857 mgc: 2.416496e-08 ustep: 1 
#> iter: 4  mgc: 1.401467e-13 
#> iter: 1  mgc: 1.401467e-13 
#> outer mgc:  3.936226 
#> iter: 1  value: 1540.861 mgc: 1.085748 ustep: 1 
#> iter: 2  value: 1540.861 mgc: 0.001196617 ustep: 1 
#> iter: 3  value: 1540.861 mgc: 6.129632e-08 ustep: 1 
#> iter: 4  mgc: 1.228904e-13 
#> iter: 1  mgc: 1.228904e-13 
#> outer mgc:  2.18728 
#> iter: 1  value: 1540.97 mgc: 0.5706925 ustep: 1 
#> iter: 2  value: 1540.97 mgc: 0.0003760422 ustep: 1 
#> iter: 3  mgc: 8.558283e-10 
#> iter: 1  mgc: 8.558283e-10 
#> outer mgc:  3.359641 
#> iter: 1  value: 1540.703 mgc: 0.2821602 ustep: 1 
#> iter: 2  value: 1540.703 mgc: 0.0004122308 ustep: 1 
#> iter: 3  mgc: 5.991392e-09 
#> iter: 1  mgc: 5.991392e-09 
#> outer mgc:  6.341819 
#> iter: 1  value: 1540.907 mgc: 0.8524876 ustep: 1 
#> iter: 2  value: 1540.907 mgc: 0.001470361 ustep: 1 
#> iter: 3  mgc: 6.000797e-09 
#> iter: 1  mgc: 6.000797e-09 
#> outer mgc:  4.668292 
#> iter: 1  value: 1541.005 mgc: 0.400672 ustep: 1 
#> iter: 2  value: 1541.005 mgc: 0.0004181917 ustep: 1 
#> iter: 3  mgc: 2.172079e-09 
#> iter: 1  mgc: 2.172079e-09 
#> outer mgc:  3.780603 
#> iter: 1  value: 1540.404 mgc: 0.2831392 ustep: 1 
#> iter: 2  value: 1540.404 mgc: 0.0004318634 ustep: 1 
#> iter: 3  value: 1540.404 mgc: 1.61862e-08 ustep: 1 
#> mgc: 1.168528e-13 
#> iter: 1  mgc: 1.168528e-13 
#> outer mgc:  2.46001 
#> iter: 1  value: 1540.379 mgc: 0.4514583 ustep: 1 
#> iter: 2  value: 1540.379 mgc: 0.001128114 ustep: 1 
#> iter: 3  value: 1540.379 mgc: 5.292723e-08 ustep: 1 
#> iter: 4  mgc: 1.503177e-13 
#> iter: 1  mgc: 1.503177e-13 
#> outer mgc:  7.802314 
#> iter: 1  value: 1540.456 mgc: 0.5716526 ustep: 1 
#> iter: 2  value: 1540.456 mgc: 0.0003447361 ustep: 1 
#> iter: 3  mgc: 2.959622e-09 
#> iter: 1  mgc: 2.959622e-09 
#> outer mgc:  4.229246 
#> iter: 1  value: 1540.599 mgc: 0.2966299 ustep: 1 
#> iter: 2  value: 1540.599 mgc: 0.000415746 ustep: 1 
#> iter: 3  value: 1540.599 mgc: 1.019756e-08 ustep: 1 
#> mgc: 1.677541e-13 
#> iter: 1  mgc: 1.677541e-13 
#> outer mgc:  1.354626 
#> iter: 1  value: 1540.75 mgc: 0.4453401 ustep: 1 
#> iter: 2  value: 1540.75 mgc: 0.0002643453 ustep: 1 
#> iter: 3  mgc: 1.561048e-09 
#> iter: 1  mgc: 1.561048e-09 
#> outer mgc:  1.377543 
#> iter: 1  value: 1540.795 mgc: 0.463023 ustep: 1 
#> iter: 2  value: 1540.795 mgc: 0.0008118101 ustep: 1 
#> iter: 3  value: 1540.795 mgc: 2.224727e-08 ustep: 1 
#> mgc: 1.587425e-13 
#> iter: 1  mgc: 1.587425e-13 
#> outer mgc:  4.388897 
#> iter: 1  value: 1540.648 mgc: 0.6256758 ustep: 1 
#> iter: 2  value: 1540.648 mgc: 0.001371873 ustep: 1 
#> iter: 3  value: 1540.648 mgc: 3.864977e-08 ustep: 1 
#> mgc: 1.343859e-13 
#> iter: 1  value: 1540.72 mgc: 0.5553506 ustep: 1 
#> iter: 2  value: 1540.72 mgc: 0.0005368971 ustep: 1 
#> iter: 3  mgc: 2.903734e-09 
#> iter: 1  value: 1540.766 mgc: 0.6145041 ustep: 1 
#> iter: 2  value: 1540.766 mgc: 0.0005039167 ustep: 1 
#> iter: 3  mgc: 4.229191e-09 
#> iter: 1  value: 1540.792 mgc: 0.6453891 ustep: 1 
#> iter: 2  value: 1540.792 mgc: 0.0006206741 ustep: 1 
#> iter: 3  value: 1540.792 mgc: 1.088608e-08 ustep: 1 
#> mgc: 8.338334e-14 
#> iter: 1  value: 1540.777 mgc: 0.2617168 ustep: 1 
#> iter: 2  value: 1540.777 mgc: 0.0001045744 ustep: 1 
#> iter: 3  mgc: 3.76544e-10 
#> iter: 1  mgc: 3.76544e-10 
#> outer mgc:  1.98867 
#> iter: 1  value: 1540.775 mgc: 0.1262957 ustep: 1 
#> iter: 2  value: 1540.775 mgc: 2.568165e-05 ustep: 1 
#> iter: 3  mgc: 3.272684e-11 
#> iter: 1  mgc: 3.272684e-11 
#> outer mgc:  1.815585 
#> iter: 1  value: 1540.757 mgc: 0.07498217 ustep: 1 
#> iter: 2  value: 1540.757 mgc: 1.186811e-05 ustep: 1 
#> iter: 3  mgc: 2.072276e-12 
#> iter: 1  mgc: 2.072276e-12 
#> outer mgc:  2.046422 
#> iter: 1  value: 1540.745 mgc: 0.1283313 ustep: 1 
#> iter: 2  value: 1540.745 mgc: 4.749455e-05 ustep: 1 
#> iter: 3  mgc: 1.39537e-10 
#> iter: 1  mgc: 1.39537e-10 
#> outer mgc:  1.620924 
#> iter: 1  value: 1540.712 mgc: 0.09577915 ustep: 1 
#> iter: 2  value: 1540.712 mgc: 3.923964e-05 ustep: 1 
#> iter: 3  mgc: 4.413178e-11 
#> iter: 1  value: 1540.666 mgc: 0.08199414 ustep: 1 
#> iter: 2  value: 1540.666 mgc: 0.0001755945 ustep: 1 
#> iter: 3  mgc: 8.907335e-10 
#> iter: 1  value: 1540.565 mgc: 0.1986889 ustep: 1 
#> iter: 2  value: 1540.565 mgc: 0.001021312 ustep: 1 
#> iter: 3  value: 1540.565 mgc: 3.011613e-08 ustep: 1 
#> mgc: 1.53671e-13 
#> iter: 1  mgc: 1.53671e-13 
#> outer mgc:  2.466329 
#> iter: 1  value: 1540.505 mgc: 0.3529984 ustep: 1 
#> iter: 2  value: 1540.505 mgc: 0.002116057 ustep: 1 
#> iter: 3  value: 1540.505 mgc: 1.791395e-07 ustep: 1 
#> iter: 4  mgc: 1.425617e-13 
#> iter: 1  mgc: 1.425617e-13 
#> outer mgc:  1.597328 
#> iter: 1  value: 1540.611 mgc: 1.154447 ustep: 1 
#> iter: 2  value: 1540.611 mgc: 0.003264725 ustep: 1 
#> iter: 3  value: 1540.611 mgc: 5.265506e-07 ustep: 1 
#> iter: 4  mgc: 1.91528e-13 
#> iter: 1  mgc: 1.91528e-13 
#> outer mgc:  7.357748 
#> iter: 1  value: 1540.568 mgc: 0.4522606 ustep: 1 
#> iter: 2  value: 1540.568 mgc: 0.00153263 ustep: 1 
#> iter: 3  value: 1540.568 mgc: 7.131462e-08 ustep: 1 
#> iter: 4  mgc: 1.248673e-13 
#> iter: 1  mgc: 1.248673e-13 
#> outer mgc:  6.067387 
#> iter: 1  value: 1540.659 mgc: 0.4266288 ustep: 1 
#> iter: 2  value: 1540.659 mgc: 0.001307643 ustep: 1 
#> iter: 3  value: 1540.659 mgc: 7.147052e-08 ustep: 1 
#> mgc: 1.513892e-13 
#> iter: 1  mgc: 1.513892e-13 
#> outer mgc:  3.85616 
#> iter: 1  value: 1540.467 mgc: 0.9601971 ustep: 1 
#> iter: 2  value: 1540.467 mgc: 0.005656157 ustep: 1 
#> iter: 3  value: 1540.467 mgc: 7.049839e-07 ustep: 1 
#> iter: 4  mgc: 1.262307e-13 
#> iter: 1  value: 1540.566 mgc: 0.3648549 ustep: 1 
#> iter: 2  value: 1540.566 mgc: 0.001139078 ustep: 1 
#> iter: 3  value: 1540.566 mgc: 6.525652e-08 ustep: 1 
#> iter: 4  mgc: 1.523242e-13 
#> iter: 1  mgc: 1.523242e-13 
#> outer mgc:  5.092477 
#> iter: 1  value: 1540.44 mgc: 0.6880564 ustep: 1 
#> iter: 2  value: 1540.44 mgc: 0.001245246 ustep: 1 
#> iter: 3  value: 1540.44 mgc: 2.991767e-08 ustep: 1 
#> mgc: 1.047615e-13 
#> iter: 1  value: 1540.483 mgc: 0.5125975 ustep: 1 
#> iter: 2  value: 1540.483 mgc: 0.0004100247 ustep: 1 
#> iter: 3  mgc: 1.304816e-09 
#> iter: 1  mgc: 1.304816e-09 
#> outer mgc:  1.107012 
#> iter: 1  value: 1540.571 mgc: 0.5152146 ustep: 1 
#> iter: 2  value: 1540.571 mgc: 0.0002305305 ustep: 1 
#> iter: 3  mgc: 8.360791e-10 
#> iter: 1  mgc: 8.360791e-10 
#> outer mgc:  2.708328 
#> iter: 1  value: 1540.633 mgc: 0.4606379 ustep: 1 
#> iter: 2  value: 1540.633 mgc: 0.0007864408 ustep: 1 
#> iter: 3  value: 1540.633 mgc: 3.698383e-08 ustep: 1 
#> iter: 4  mgc: 2.138906e-13 
#> iter: 1  mgc: 2.138906e-13 
#> outer mgc:  3.931374 
#> iter: 1  value: 1540.644 mgc: 0.511006 ustep: 1 
#> iter: 2  value: 1540.644 mgc: 0.0002493753 ustep: 1 
#> iter: 3  mgc: 6.476719e-10 
#> iter: 1  mgc: 6.476719e-10 
#> outer mgc:  1.459358 
#> iter: 1  value: 1540.612 mgc: 0.279262 ustep: 1 
#> iter: 2  value: 1540.612 mgc: 0.0001822413 ustep: 1 
#> iter: 3  mgc: 9.654967e-10 
#> iter: 1  mgc: 9.654967e-10 
#> outer mgc:  0.6972362 
#> iter: 1  value: 1540.596 mgc: 0.2531989 ustep: 1 
#> iter: 2  value: 1540.596 mgc: 0.0001798827 ustep: 1 
#> iter: 3  mgc: 1.3312e-09 
#> iter: 1  mgc: 1.3312e-09 
#> outer mgc:  2.310701 
#> iter: 1  value: 1540.559 mgc: 0.2761941 ustep: 1 
#> iter: 2  value: 1540.559 mgc: 0.0002022961 ustep: 1 
#> iter: 3  mgc: 9.948088e-10 
#> iter: 1  value: 1540.572 mgc: 0.2495961 ustep: 1 
#> iter: 2  value: 1540.572 mgc: 6.139878e-05 ustep: 1 
#> iter: 3  mgc: 6.703038e-11 
#> iter: 1  value: 1540.579 mgc: 0.252347 ustep: 1 
#> iter: 2  value: 1540.579 mgc: 7.290323e-05 ustep: 1 
#> iter: 3  mgc: 1.286055e-10 
#> iter: 1  value: 1540.584 mgc: 0.2697862 ustep: 1 
#> iter: 2  value: 1540.584 mgc: 0.0001008346 ustep: 1 
#> iter: 3  mgc: 3.117054e-10 
#> iter: 1  value: 1540.587 mgc: 0.2549214 ustep: 1 
#> iter: 2  value: 1540.587 mgc: 9.78134e-05 ustep: 1 
#> iter: 3  mgc: 3.30232e-10 
#> iter: 1  value: 1540.59 mgc: 0.1034795 ustep: 1 
#> iter: 2  value: 1540.59 mgc: 1.608424e-05 ustep: 1 
#> iter: 3  mgc: 8.975417e-12 
#> iter: 1  mgc: 8.975417e-12 
#> outer mgc:  0.3268373 
#> iter: 1  value: 1540.588 mgc: 0.02451057 ustep: 1 
#> iter: 2  value: 1540.588 mgc: 9.603925e-07 ustep: 1 
#> iter: 3  mgc: 1.827718e-13 
#> iter: 1  value: 1540.585 mgc: 0.0396434 ustep: 1 
#> iter: 2  value: 1540.585 mgc: 2.512438e-06 ustep: 1 
#> iter: 3  mgc: 2.010955e-13 
#> iter: 1  mgc: 2.010955e-13 
#> outer mgc:  0.6912881 
#> iter: 1  value: 1540.58 mgc: 0.05386596 ustep: 1 
#> iter: 2  value: 1540.58 mgc: 2.353844e-06 ustep: 1 
#> iter: 3  mgc: 1.697378e-13 
#> iter: 1  mgc: 1.697378e-13 
#> outer mgc:  0.6363633 
#> iter: 1  value: 1540.578 mgc: 0.04520223 ustep: 1 
#> iter: 2  value: 1540.578 mgc: 3.633286e-06 ustep: 1 
#> iter: 3  mgc: 6.850097e-13 
#> iter: 1  mgc: 6.850097e-13 
#> outer mgc:  0.4046531 
#> iter: 1  value: 1540.575 mgc: 0.02328942 ustep: 1 
#> iter: 2  value: 1540.575 mgc: 1.302969e-06 ustep: 1 
#> iter: 3  mgc: 1.485748e-13 
#> iter: 1  value: 1540.573 mgc: 0.01208981 ustep: 1 
#> iter: 2  value: 1540.573 mgc: 1.998762e-06 ustep: 1 
#> iter: 3  mgc: 1.268554e-13 
#> iter: 1  value: 1540.567 mgc: 0.03698315 ustep: 1 
#> iter: 2  value: 1540.567 mgc: 1.86951e-05 ustep: 1 
#> iter: 3  mgc: 1.020572e-11 
#> iter: 1  value: 1540.558 mgc: 0.06230768 ustep: 1 
#> iter: 2  value: 1540.558 mgc: 5.303347e-05 ustep: 1 
#> iter: 3  mgc: 8.243858e-11 
#> iter: 1  mgc: 8.243858e-11 
#> outer mgc:  0.7122878 
#> iter: 1  value: 1540.561 mgc: 0.1030704 ustep: 1 
#> iter: 2  value: 1540.561 mgc: 0.0001943533 ustep: 1 
#> iter: 3  mgc: 9.021382e-10 
#> iter: 1  mgc: 9.021382e-10 
#> outer mgc:  0.4807467 
#> iter: 1  value: 1540.566 mgc: 0.08349314 ustep: 1 
#> iter: 2  value: 1540.566 mgc: 0.0001155982 ustep: 1 
#> iter: 3  mgc: 6.361362e-10 
#> iter: 1  mgc: 6.361362e-10 
#> outer mgc:  0.9313742 
#> iter: 1  value: 1540.556 mgc: 0.1013537 ustep: 1 
#> iter: 2  value: 1540.556 mgc: 0.0001263272 ustep: 1 
#> iter: 3  mgc: 6.18364e-10 
#> iter: 1  mgc: 6.18364e-10 
#> outer mgc:  0.6354676 
#> iter: 1  value: 1540.547 mgc: 0.3352973 ustep: 1 
#> iter: 2  value: 1540.547 mgc: 0.0001068807 ustep: 1 
#> iter: 3  mgc: 1.410188e-10 
#> iter: 1  value: 1540.55 mgc: 0.2108892 ustep: 1 
#> iter: 2  value: 1540.55 mgc: 3.051245e-05 ustep: 1 
#> iter: 3  mgc: 1.402721e-11 
#> iter: 1  value: 1540.552 mgc: 0.1508576 ustep: 1 
#> iter: 2  value: 1540.552 mgc: 1.991385e-05 ustep: 1 
#> iter: 3  mgc: 9.742797e-12 
#> iter: 1  value: 1540.552 mgc: 0.1438832 ustep: 1 
#> iter: 2  value: 1540.552 mgc: 2.742816e-05 ustep: 1 
#> iter: 3  mgc: 2.38248e-11 
#> iter: 1  value: 1540.554 mgc: 0.05348891 ustep: 1 
#> iter: 2  value: 1540.554 mgc: 4.264079e-06 ustep: 1 
#> iter: 3  mgc: 6.653836e-13 
#> iter: 1  mgc: 6.653836e-13 
#> outer mgc:  0.5400231 
#> iter: 1  value: 1540.555 mgc: 0.02611608 ustep: 1 
#> iter: 2  value: 1540.555 mgc: 1.040888e-06 ustep: 1 
#> iter: 3  mgc: 1.334763e-13 
#> iter: 1  mgc: 1.334763e-13 
#> outer mgc:  0.2072348 
#> iter: 1  value: 1540.555 mgc: 0.01241673 ustep: 1 
#> iter: 2  value: 1540.555 mgc: 2.022558e-07 ustep: 1 
#> iter: 3  mgc: 1.301304e-13 
#> iter: 1  value: 1540.556 mgc: 0.01900037 ustep: 1 
#> iter: 2  value: 1540.556 mgc: 4.735929e-07 ustep: 1 
#> iter: 3  mgc: 1.688157e-13 
#> iter: 1  mgc: 1.688157e-13 
#> outer mgc:  0.5015708 
#> iter: 1  value: 1540.557 mgc: 0.02421969 ustep: 1 
#> iter: 2  value: 1540.557 mgc: 1.698257e-06 ustep: 1 
#> iter: 3  mgc: 1.537926e-13 
#> iter: 1  mgc: 1.537926e-13 
#> outer mgc:  0.1885744 
#> iter: 1  value: 1540.556 mgc: 0.01679768 ustep: 1 
#> iter: 2  value: 1540.556 mgc: 2.313116e-07 ustep: 1 
#> iter: 3  mgc: 1.18324e-13 
#> iter: 1  value: 1540.555 mgc: 0.02923689 ustep: 1 
#> iter: 2  value: 1540.555 mgc: 1.375333e-06 ustep: 1 
#> iter: 3  mgc: 1.247133e-13 
#> iter: 1  value: 1540.551 mgc: 0.1169551 ustep: 1 
#> iter: 2  value: 1540.551 mgc: 2.201457e-05 ustep: 1 
#> iter: 3  mgc: 1.461634e-11 
#> iter: 1  mgc: 1.461634e-11 
#> outer mgc:  0.2564583 
#> iter: 1  value: 1540.556 mgc: 0.0432245 ustep: 1 
#> iter: 2  value: 1540.556 mgc: 3.177978e-05 ustep: 1 
#> iter: 3  mgc: 2.368063e-11 
#> iter: 1  value: 1540.561 mgc: 0.04544111 ustep: 1 
#> iter: 2  value: 1540.561 mgc: 4.393958e-05 ustep: 1 
#> iter: 3  mgc: 4.510817e-11 
#> iter: 1  mgc: 4.510817e-11 
#> outer mgc:  0.4169681 
#> iter: 1  value: 1540.574 mgc: 0.05780013 ustep: 1 
#> iter: 2  value: 1540.574 mgc: 6.798022e-05 ustep: 1 
#> iter: 3  mgc: 1.737005e-10 
#> iter: 1  mgc: 1.737005e-10 
#> outer mgc:  0.3329653 
#> iter: 1  value: 1540.572 mgc: 0.09338963 ustep: 1 
#> iter: 2  value: 1540.572 mgc: 0.0001608391 ustep: 1 
#> iter: 3  mgc: 7.430838e-10 
#> iter: 1  mgc: 7.430838e-10 
#> outer mgc:  0.129021 
#> iter: 1  value: 1540.57 mgc: 0.09152539 ustep: 1 
#> iter: 2  value: 1540.57 mgc: 3.093408e-05 ustep: 1 
#> iter: 3  mgc: 3.962491e-11 
#> iter: 1  mgc: 3.962491e-11 
#> outer mgc:  0.2434734 
#> iter: 1  value: 1540.565 mgc: 0.08587928 ustep: 1 
#> iter: 2  value: 1540.565 mgc: 2.301097e-05 ustep: 1 
#> iter: 3  mgc: 9.3912e-12 
#> iter: 1  mgc: 9.3912e-12 
#> outer mgc:  0.2100995 
#> iter: 1  value: 1540.565 mgc: 0.1063948 ustep: 1 
#> iter: 2  value: 1540.565 mgc: 1.090961e-05 ustep: 1 
#> iter: 3  mgc: 8.449688e-13 
#> iter: 1  mgc: 8.449688e-13 
#> outer mgc:  0.4014802 
#> iter: 1  value: 1540.561 mgc: 0.1915212 ustep: 1 
#> iter: 2  value: 1540.561 mgc: 4.225658e-05 ustep: 1 
#> iter: 3  mgc: 6.452566e-12 
#> iter: 1  mgc: 6.452566e-12 
#> outer mgc:  0.3568506 
#> iter: 1  value: 1540.567 mgc: 0.1203727 ustep: 1 
#> iter: 2  value: 1540.567 mgc: 2.258447e-05 ustep: 1 
#> iter: 3  mgc: 1.931869e-11 
#> iter: 1  mgc: 1.931869e-11 
#> outer mgc:  1.323469 
#> iter: 1  value: 1540.574 mgc: 0.1675625 ustep: 1 
#> iter: 2  value: 1540.574 mgc: 1.500753e-05 ustep: 1 
#> iter: 3  mgc: 6.101419e-12 
#> iter: 1  mgc: 6.101419e-12 
#> outer mgc:  0.4446626 
#> iter: 1  value: 1540.553 mgc: 0.1734507 ustep: 1 
#> iter: 2  value: 1540.553 mgc: 3.50251e-05 ustep: 1 
#> iter: 3  mgc: 8.071238e-11 
#> iter: 1  mgc: 8.071238e-11 
#> outer mgc:  1.604477 
#> iter: 1  value: 1540.529 mgc: 0.1213536 ustep: 1 
#> iter: 2  value: 1540.529 mgc: 1.731875e-05 ustep: 1 
#> iter: 3  mgc: 5.601765e-12 
#> iter: 1  mgc: 5.601765e-12 
#> outer mgc:  0.4442598 
#> iter: 1  value: 1540.555 mgc: 0.2223221 ustep: 1 
#> iter: 2  value: 1540.555 mgc: 4.978965e-05 ustep: 1 
#> iter: 3  mgc: 1.655272e-10 
#> iter: 1  mgc: 1.655272e-10 
#> outer mgc:  1.189551 
#> iter: 1  value: 1540.568 mgc: 0.1991152 ustep: 1 
#> iter: 2  value: 1540.568 mgc: 8.582036e-05 ustep: 1 
#> iter: 3  mgc: 8.529139e-11 
#> iter: 1  mgc: 8.529139e-11 
#> outer mgc:  0.4559004 
#> iter: 1  value: 1540.596 mgc: 0.1785391 ustep: 1 
#> iter: 2  value: 1540.596 mgc: 6.504422e-05 ustep: 1 
#> iter: 3  mgc: 2.639195e-10 
#> iter: 1  mgc: 2.639195e-10 
#> outer mgc:  1.169216 
#> iter: 1  value: 1540.598 mgc: 0.1950972 ustep: 1 
#> iter: 2  value: 1540.598 mgc: 2.349224e-05 ustep: 1 
#> iter: 3  mgc: 1.286109e-11 
#> iter: 1  mgc: 1.286109e-11 
#> outer mgc:  0.5186548 
#> iter: 1  value: 1540.569 mgc: 0.1489258 ustep: 1 
#> iter: 2  value: 1540.569 mgc: 4.606875e-05 ustep: 1 
#> iter: 3  mgc: 1.602649e-11 
#> iter: 1  mgc: 1.602649e-11 
#> outer mgc:  0.8122973 
#> iter: 1  value: 1540.531 mgc: 0.1045296 ustep: 1 
#> iter: 2  value: 1540.531 mgc: 3.457492e-05 ustep: 1 
#> iter: 3  mgc: 2.106446e-11 
#> iter: 1  mgc: 2.106446e-11 
#> outer mgc:  0.3855716 
#> iter: 1  value: 1540.483 mgc: 0.1451672 ustep: 1 
#> iter: 2  value: 1540.483 mgc: 8.418202e-05 ustep: 1 
#> iter: 3  mgc: 1.964253e-10 
#> iter: 1  mgc: 1.964253e-10 
#> outer mgc:  0.3586984 
#> iter: 1  value: 1540.498 mgc: 0.06169304 ustep: 1 
#> iter: 2  value: 1540.498 mgc: 7.181954e-06 ustep: 1 
#> iter: 3  mgc: 5.804798e-13 
#> iter: 1  mgc: 5.804798e-13 
#> outer mgc:  0.2009154 
#> iter: 1  value: 1540.502 mgc: 0.1438276 ustep: 1 
#> iter: 2  value: 1540.502 mgc: 2.56345e-05 ustep: 1 
#> iter: 3  mgc: 4.01531e-12 
#> iter: 1  mgc: 4.01531e-12 
#> outer mgc:  0.3897005 
#> iter: 1  value: 1540.537 mgc: 0.08802762 ustep: 1 
#> iter: 2  value: 1540.537 mgc: 1.667102e-05 ustep: 1 
#> iter: 3  mgc: 8.871931e-13 
#> iter: 1  mgc: 8.871931e-13 
#> outer mgc:  0.1867671 
#> iter: 1  value: 1540.538 mgc: 0.05186165 ustep: 1 
#> iter: 2  value: 1540.538 mgc: 1.364281e-05 ustep: 1 
#> iter: 3  mgc: 4.476647e-12 
#> iter: 1  mgc: 4.476647e-12 
#> outer mgc:  0.3196255 
#> iter: 1  value: 1540.529 mgc: 0.09180527 ustep: 1 
#> iter: 2  value: 1540.529 mgc: 1.611115e-05 ustep: 1 
#> iter: 3  mgc: 4.822373e-12 
#> iter: 1  mgc: 4.822373e-12 
#> outer mgc:  0.2565409 
#> iter: 1  value: 1540.522 mgc: 0.06676667 ustep: 1 
#> iter: 2  value: 1540.522 mgc: 5.591754e-06 ustep: 1 
#> iter: 3  mgc: 7.838138e-13 
#> iter: 1  mgc: 7.838138e-13 
#> outer mgc:  0.2625251 
#> iter: 1  value: 1540.523 mgc: 0.08938688 ustep: 1 
#> iter: 2  value: 1540.523 mgc: 1.846414e-05 ustep: 1 
#> iter: 3  mgc: 3.68363e-12 
#> iter: 1  mgc: 3.68363e-12 
#> outer mgc:  0.5208377 
#> iter: 1  value: 1540.529 mgc: 0.07804466 ustep: 1 
#> iter: 2  value: 1540.529 mgc: 1.711118e-05 ustep: 1 
#> iter: 3  mgc: 4.825219e-12 
#> iter: 1  mgc: 4.825219e-12 
#> outer mgc:  0.2566215 
#> iter: 1  value: 1540.515 mgc: 0.03692695 ustep: 1 
#> iter: 2  value: 1540.515 mgc: 3.36921e-06 ustep: 1 
#> iter: 3  mgc: 4.078925e-13 
#> iter: 1  mgc: 4.078925e-13 
#> outer mgc:  0.1872298 
#> iter: 1  value: 1540.519 mgc: 0.07741648 ustep: 1 
#> iter: 2  value: 1540.519 mgc: 7.139412e-06 ustep: 1 
#> iter: 3  mgc: 1.831286e-12 
#> iter: 1  mgc: 1.831286e-12 
#> outer mgc:  0.2153993 
#> iter: 1  value: 1540.524 mgc: 0.08712658 ustep: 1 
#> iter: 2  value: 1540.524 mgc: 1.637894e-05 ustep: 1 
#> iter: 3  mgc: 1.36801e-11 
#> iter: 1  mgc: 1.36801e-11 
#> outer mgc:  0.2876342 
#> iter: 1  value: 1540.571 mgc: 0.2220763 ustep: 1 
#> iter: 2  value: 1540.571 mgc: 7.684566e-05 ustep: 1 
#> iter: 3  mgc: 1.282501e-10 
#> iter: 1  value: 1540.538 mgc: 0.06558474 ustep: 1 
#> iter: 2  value: 1540.538 mgc: 8.452786e-06 ustep: 1 
#> iter: 3  mgc: 1.309568e-12 
#> iter: 1  mgc: 1.309568e-12 
#> outer mgc:  0.1519961 
#> iter: 1  value: 1540.533 mgc: 0.06049321 ustep: 1 
#> iter: 2  value: 1540.533 mgc: 5.882284e-06 ustep: 1 
#> iter: 3  mgc: 1.908097e-13 
#> iter: 1  mgc: 1.908097e-13 
#> outer mgc:  0.2542028 
#> iter: 1  value: 1540.526 mgc: 0.05090035 ustep: 1 
#> iter: 2  value: 1540.526 mgc: 2.25922e-06 ustep: 1 
#> iter: 3  mgc: 2.256369e-13 
#> iter: 1  mgc: 2.256369e-13 
#> outer mgc:  0.1509034 
#> iter: 1  value: 1540.519 mgc: 0.07327251 ustep: 1 
#> iter: 2  value: 1540.519 mgc: 1.090456e-05 ustep: 1 
#> iter: 3  mgc: 3.051589e-12 
#> iter: 1  mgc: 3.051589e-12 
#> outer mgc:  0.163867 
#> iter: 1  value: 1540.52 mgc: 0.03784985 ustep: 1 
#> iter: 2  value: 1540.52 mgc: 2.23564e-06 ustep: 1 
#> iter: 3  mgc: 1.208747e-13 
#> iter: 1  mgc: 1.208747e-13 
#> outer mgc:  0.1207511 
#> iter: 1  value: 1540.516 mgc: 0.03283415 ustep: 1 
#> iter: 2  value: 1540.516 mgc: 4.129248e-06 ustep: 1 
#> iter: 3  mgc: 4.421653e-13 
#> iter: 1  mgc: 4.421653e-13 
#> outer mgc:  0.105723 
#> iter: 1  value: 1540.513 mgc: 0.04227213 ustep: 1 
#> iter: 2  value: 1540.513 mgc: 2.403234e-06 ustep: 1 
#> iter: 3  mgc: 1.211893e-13 
#> iter: 1  mgc: 1.211893e-13 
#> outer mgc:  0.1030712 
#> iter: 1  value: 1540.518 mgc: 0.03537415 ustep: 1 
#> iter: 2  value: 1540.518 mgc: 2.910924e-06 ustep: 1 
#> iter: 3  mgc: 1.990701e-13 
#> iter: 1  mgc: 1.990701e-13 
#> outer mgc:  0.1644156 
#> iter: 1  value: 1540.522 mgc: 0.03809032 ustep: 1 
#> iter: 2  value: 1540.522 mgc: 2.240939e-06 ustep: 1 
#> iter: 3  mgc: 1.403306e-13 
#> iter: 1  mgc: 1.403306e-13 
#> outer mgc:  0.07513752 
#> iter: 1  value: 1540.521 mgc: 0.03372834 ustep: 1 
#> iter: 2  value: 1540.521 mgc: 1.027368e-06 ustep: 1 
#> iter: 3  mgc: 1.636887e-13 
#> iter: 1  mgc: 1.636887e-13 
#> outer mgc:  0.06483881 
#> iter: 1  value: 1540.517 mgc: 0.05804625 ustep: 1 
#> iter: 2  value: 1540.517 mgc: 7.209859e-06 ustep: 1 
#> iter: 3  mgc: 7.33523e-13 
#> iter: 1  value: 1540.519 mgc: 0.02107726 ustep: 1 
#> iter: 2  value: 1540.519 mgc: 9.506838e-07 ustep: 1 
#> iter: 3  mgc: 1.002126e-13 
#> iter: 1  mgc: 1.002126e-13 
#> outer mgc:  0.05325636 
#> iter: 1  value: 1540.522 mgc: 0.04653099 ustep: 1 
#> iter: 2  value: 1540.522 mgc: 6.405485e-06 ustep: 1 
#> iter: 3  mgc: 2.242537e-13 
#> iter: 1  mgc: 2.242537e-13 
#> outer mgc:  0.05874642 
#> iter: 1  value: 1540.522 mgc: 0.03143549 ustep: 1 
#> iter: 2  value: 1540.522 mgc: 1.374013e-06 ustep: 1 
#> iter: 3  mgc: 1.558406e-13 
#> iter: 1  mgc: 1.558406e-13 
#> outer mgc:  0.07193887 
#> iter: 1  value: 1540.52 mgc: 0.02937296 ustep: 1 
#> iter: 2  value: 1540.52 mgc: 2.599484e-06 ustep: 1 
#> iter: 3  mgc: 1.845687e-13 
#> iter: 1  mgc: 1.845687e-13 
#> outer mgc:  0.04610762 
#> iter: 1  value: 1540.52 mgc: 0.03907315 ustep: 1 
#> iter: 2  value: 1540.52 mgc: 2.502669e-06 ustep: 1 
#> iter: 3  mgc: 1.325826e-13 
#> iter: 1  mgc: 1.325826e-13 
#> outer mgc:  0.07082905 
#> iter: 1  value: 1540.519 mgc: 0.01957235 ustep: 1 
#> iter: 2  value: 1540.519 mgc: 6.46639e-07 ustep: 1 
#> iter: 3  mgc: 1.523424e-13 
#> iter: 1  mgc: 1.523424e-13 
#> outer mgc:  0.0426529 
#> iter: 1  value: 1540.52 mgc: 0.03618948 ustep: 1 
#> iter: 2  value: 1540.52 mgc: 7.838915e-07 ustep: 1 
#> iter: 3  mgc: 1.635164e-13 
#> iter: 1  mgc: 1.635164e-13 
#> outer mgc:  0.0465363 
#> iter: 1  value: 1540.519 mgc: 0.01298783 ustep: 1 
#> iter: 2  value: 1540.519 mgc: 4.098941e-07 ustep: 1 
#> iter: 3  mgc: 1.93375e-13 
#> iter: 1  mgc: 1.93375e-13 
#> outer mgc:  0.02935276 
#> iter: 1  value: 1540.521 mgc: 0.01347953 ustep: 1 
#> iter: 2  value: 1540.521 mgc: 9.030946e-07 ustep: 1 
#> iter: 3  mgc: 1.376983e-13 
#> iter: 1  mgc: 1.376983e-13 
#> outer mgc:  0.03124451 
#> iter: 1  value: 1540.519 mgc: 0.03194314 ustep: 1 
#> iter: 2  value: 1540.519 mgc: 2.863476e-06 ustep: 1 
#> iter: 3  mgc: 1.775368e-13 
#> iter: 1  mgc: 1.775368e-13 
#> outer mgc:  0.03800708 
#> iter: 1  value: 1540.519 mgc: 0.02128468 ustep: 1 
#> iter: 2  value: 1540.519 mgc: 1.243772e-06 ustep: 1 
#> iter: 3  mgc: 1.488241e-13 
#> iter: 1  mgc: 1.488241e-13 
#> outer mgc:  0.02821498 
#> iter: 1  value: 1540.52 mgc: 0.01410883 ustep: 1 
#> iter: 2  value: 1540.52 mgc: 5.645127e-07 ustep: 1 
#> iter: 3  mgc: 1.455913e-13 
#> iter: 1  value: 1540.521 mgc: 0.01449513 ustep: 1 
#> iter: 2  value: 1540.521 mgc: 5.933576e-07 ustep: 1 
#> iter: 3  mgc: 9.75481e-14 
#> iter: 1  mgc: 9.75481e-14 
#> outer mgc:  0.01688005 
#> iter: 1  value: 1540.52 mgc: 0.01958254 ustep: 1 
#> iter: 2  value: 1540.52 mgc: 5.188067e-07 ustep: 1 
#> iter: 3  mgc: 1.313864e-13 
#> iter: 1  mgc: 1.313864e-13 
#> outer mgc:  0.02094265 
#> iter: 1  value: 1540.519 mgc: 0.009122875 ustep: 1 
#> iter: 2  value: 1540.519 mgc: 3.232392e-07 ustep: 1 
#> iter: 3  mgc: 1.246134e-13 
#> iter: 1  mgc: 1.246134e-13 
#> outer mgc:  0.01385052 
#> iter: 1  value: 1540.519 mgc: 0.01035409 ustep: 1 
#> iter: 2  value: 1540.519 mgc: 1.499923e-07 ustep: 1 
#> iter: 3  mgc: 2.172278e-13 
#> iter: 1  mgc: 2.172278e-13 
#> outer mgc:  0.01072753 
#> iter: 1  value: 1540.519 mgc: 0.006603622 ustep: 1 
#> iter: 2  value: 1540.519 mgc: 1.544729e-07 ustep: 1 
#> iter: 3  mgc: 1.514398e-13 
#> iter: 1  mgc: 1.514398e-13 
#> outer mgc:  0.01464834 
#> iter: 1  value: 1540.52 mgc: 0.006896831 ustep: 1 
#> iter: 2  value: 1540.52 mgc: 3.758076e-08 ustep: 1 
#> mgc: 1.421647e-13 
#> iter: 1  mgc: 1.421647e-13 
#> outer mgc:  0.01239297 
#> iter: 1  value: 1540.521 mgc: 0.003603659 ustep: 1 
#> iter: 2  value: 1540.521 mgc: 3.766559e-08 ustep: 1 
#> mgc: 1.25145e-13 
#> iter: 1  mgc: 1.25145e-13 
#> outer mgc:  0.006920792 
#> iter: 1  value: 1540.521 mgc: 0.002651879 ustep: 1 
#> iter: 2  value: 1540.521 mgc: 1.079356e-08 ustep: 1 
#> mgc: 1.143935e-13 
#> iter: 1  mgc: 1.143935e-13 
#> outer mgc:  0.005412093 
#> iter: 1  value: 1540.519 mgc: 0.001756429 ustep: 1 
#> iter: 2  mgc: 5.214405e-09 
#> iter: 1  mgc: 5.214405e-09 
#> outer mgc:  0.006465728 
#> iter: 1  value: 1540.52 mgc: 0.002788971 ustep: 1 
#> iter: 2  mgc: 6.426703e-09 
#> iter: 1  mgc: 6.426703e-09 
#> outer mgc:  0.0034086 
#> iter: 1  value: 1540.52 mgc: 0.002706672 ustep: 1 
#> iter: 2  value: 1540.52 mgc: 1.237832e-08 ustep: 1 
#> mgc: 1.467079e-13 
#> iter: 1  mgc: 1.467079e-13 
#> outer mgc:  0.002396305 
#> iter: 1  value: 1540.52 mgc: 0.002409134 ustep: 1 
#> iter: 2  value: 1540.52 mgc: 1.440783e-08 ustep: 1 
#> mgc: 2.164157e-13 
#> iter: 1  mgc: 2.164157e-13 
#> outer mgc:  0.004878518 
#> iter: 1  value: 1540.52 mgc: 0.002848311 ustep: 1 
#> iter: 2  value: 1540.52 mgc: 1.006228e-08 ustep: 1 
#> mgc: 1.754864e-13 
#> iter: 1  mgc: 1.754864e-13 
#> outer mgc:  0.006555884 
#> iter: 1  value: 1540.52 mgc: 0.001370179 ustep: 1 
#> iter: 2  mgc: 2.589549e-09 
#> iter: 1  value: 1540.52 mgc: 0.001418562 ustep: 1 
#> iter: 2  mgc: 2.77459e-09 
#> iter: 1  mgc: 2.77459e-09 
#> outer mgc:  0.006797284 
#> iter: 1  value: 1540.519 mgc: 0.003406031 ustep: 1 
#> iter: 2  mgc: 7.321168e-09 
#> iter: 1  mgc: 7.321168e-09 
#> outer mgc:  0.006788784 
#> iter: 1  value: 1540.519 mgc: 0.004764849 ustep: 1 
#> iter: 2  mgc: 9.03223e-09 
#> iter: 1  mgc: 9.03223e-09 
#> outer mgc:  0.002817804 
#> iter: 1  value: 1540.52 mgc: 0.002885321 ustep: 1 
#> iter: 2  value: 1540.52 mgc: 1.324811e-08 ustep: 1 
#> mgc: 1.924915e-13 
#> iter: 1  value: 1540.52 mgc: 0.001422571 ustep: 1 
#> iter: 2  mgc: 3.208421e-09 
#> iter: 1  mgc: 3.208421e-09 
#> outer mgc:  0.002825734 
#> iter: 1  value: 1540.52 mgc: 0.002476694 ustep: 1 
#> iter: 2  value: 1540.52 mgc: 1.235988e-08 ustep: 1 
#> mgc: 1.316504e-13 
#> iter: 1  mgc: 1.316504e-13 
#> outer mgc:  0.002045087 
#> iter: 1  value: 1540.52 mgc: 0.0005143001 ustep: 1 
#> iter: 2  mgc: 2.746535e-10 
#> iter: 1  mgc: 2.746535e-10 
#> outer mgc:  0.0007376433 
#> iter: 1  mgc: 2.746535e-10 
#> iter: 1  mgc: 2.746535e-10 
#> outer mgc:  0.0007376433 
#> iter: 1  mgc: 2.746535e-10 
#> outer mgc:  0.0007376433 
#> Retro Peel: 7
#> Using stats::nlminb for optimization with these control parameters:
#> iter.max: 1000
#> eval.max: 1000
#> iter: 1  value: 1491.368 mgc: 1.822871 ustep: 1 
#> iter: 2  value: 1491.368 mgc: 0.05968741 ustep: 1 
#> iter: 3  value: 1491.368 mgc: 9.891347e-05 ustep: 1 
#> iter: 4  mgc: 2.741569e-10 
#> iter: 1  mgc: 2.741569e-10 
#> Matching hessian patterns... Done
#> outer mgc:  18.00373 
#> iter: 1  value: 15516.5 mgc: 112.5273 ustep: 1 
#> iter: 2  value: 14661.66 mgc: 36.55139 ustep: 1 
#> iter: 3  value: 14455.13 mgc: 20.50171 ustep: 1 
#> iter: 4  value: 14449.83 mgc: 3.908055 ustep: 1 
#> iter: 5  value: 14449.81 mgc: 0.1686311 ustep: 1 
#> iter: 6  value: 14449.81 mgc: 0.008204331 ustep: 1 
#> iter: 7  value: 14449.81 mgc: 0.002829102 ustep: 0.4367355 
#> iter: 8  value: 14449.8 mgc: 0.001172117 ustep: 0.1263996 
#> iter: 9  value: 14449.8 mgc: 0.002068253 ustep: 0.02707895 
#> iter: 10  value: 14449.8 mgc: 0.002546653 ustep: 0.002805174 
#> iter: 11  value: 14449.8 mgc: 0.002618015 ustep: 0.0008509692 
#> iter: 12  value: 14449.8 mgc: 0.002641813 ustep: 0.0001159367 
#> iter: 13  value: 14449.8 mgc: 0.002643859 ustep: 8.529428e-05 
#> iter: 14  value: 14449.8 mgc: 0.002645365 ustep: 9.1494e-06 
#> iter: 15  value: 14449.8 mgc: 0.002645526 ustep: 1.907353e-07 
#> mgc: 0.00264553 
#> iter: 1  value: 1532.46 mgc: 17.86072 ustep: 1 
#> iter: 2  value: 1532.409 mgc: 0.3805504 ustep: 1 
#> iter: 3  value: 1532.409 mgc: 0.00510551 ustep: 1 
#> iter: 4  value: 1532.409 mgc: 1.592377e-06 ustep: 1 
#> iter: 5  mgc: 1.491708e-13 
#> iter: 1  value: 1491.43 mgc: 1.861211 ustep: 1 
#> iter: 2  value: 1491.43 mgc: 0.004852128 ustep: 1 
#> iter: 3  value: 1491.43 mgc: 1.099333e-06 ustep: 1 
#> iter: 4  mgc: 1.587691e-13 
#> iter: 1  mgc: 1.587691e-13 
#> outer mgc:  25.53972 
#> iter: 1  value: 1491.07 mgc: 1.10989 ustep: 1 
#> iter: 2  value: 1491.07 mgc: 0.002127004 ustep: 1 
#> iter: 3  value: 1491.07 mgc: 2.176667e-07 ustep: 1 
#> iter: 4  mgc: 1.164363e-13 
#> iter: 1  mgc: 1.164363e-13 
#> outer mgc:  17.58528 
#> iter: 1  value: 1490.775 mgc: 0.4539018 ustep: 1 
#> iter: 2  value: 1490.775 mgc: 0.0004765307 ustep: 1 
#> iter: 3  value: 1490.775 mgc: 1.354601e-08 ustep: 1 
#> mgc: 1.319805e-13 
#> iter: 1  value: 1490.405 mgc: 1.222101 ustep: 1 
#> iter: 2  value: 1490.405 mgc: 0.003450934 ustep: 1 
#> iter: 3  value: 1490.405 mgc: 7.12637e-07 ustep: 1 
#> iter: 4  mgc: 1.160735e-13 
#> iter: 1  mgc: 1.160735e-13 
#> outer mgc:  33.15559 
#> iter: 1  value: 1488.503 mgc: 3.148101 ustep: 1 
#> iter: 2  value: 1488.503 mgc: 0.01807655 ustep: 1 
#> iter: 3  value: 1488.503 mgc: 4.120256e-06 ustep: 1 
#> iter: 4  mgc: 3.405997e-13 
#> iter: 1  mgc: 3.405997e-13 
#> outer mgc:  17.11289 
#> iter: 1  value: 1488.022 mgc: 5.011704 ustep: 1 
#> iter: 2  value: 1488.022 mgc: 0.04551705 ustep: 1 
#> iter: 3  value: 1488.022 mgc: 0.0001112993 ustep: 1 
#> iter: 4  mgc: 7.162275e-10 
#> iter: 1  mgc: 7.162275e-10 
#> outer mgc:  19.73755 
#> iter: 1  value: 1486.711 mgc: 3.06701 ustep: 1 
#> iter: 2  value: 1486.711 mgc: 0.01538796 ustep: 1 
#> iter: 3  value: 1486.711 mgc: 3.985267e-06 ustep: 1 
#> iter: 4  mgc: 1.326967e-12 
#> iter: 1  mgc: 1.326967e-12 
#> outer mgc:  23.41472 
#> iter: 1  value: 1488.343 mgc: 6.431528 ustep: 1 
#> iter: 2  value: 1488.341 mgc: 0.06149614 ustep: 1 
#> iter: 3  value: 1488.341 mgc: 0.0002421372 ustep: 1 
#> iter: 4  mgc: 3.517481e-09 
#> iter: 1  value: 1488.843 mgc: 5.06536 ustep: 1 
#> iter: 2  value: 1488.843 mgc: 0.04470396 ustep: 1 
#> iter: 3  value: 1488.843 mgc: 9.195169e-05 ustep: 1 
#> iter: 4  mgc: 4.985584e-10 
#> iter: 1  value: 1486.446 mgc: 1.207372 ustep: 1 
#> iter: 2  value: 1486.446 mgc: 0.002409744 ustep: 1 
#> iter: 3  value: 1486.446 mgc: 2.464156e-07 ustep: 1 
#> iter: 4  mgc: 1.248859e-13 
#> iter: 1  mgc: 1.248859e-13 
#> outer mgc:  20.89533 
#> iter: 1  value: 1486.348 mgc: 0.7136258 ustep: 1 
#> iter: 2  value: 1486.348 mgc: 0.0007590249 ustep: 1 
#> iter: 3  value: 1486.348 mgc: 2.355191e-08 ustep: 1 
#> iter: 4  mgc: 1.232935e-13 
#> iter: 1  mgc: 1.232935e-13 
#> outer mgc:  17.34178 
#> iter: 1  value: 1486.024 mgc: 1.074244 ustep: 1 
#> iter: 2  value: 1486.024 mgc: 0.001820085 ustep: 1 
#> iter: 3  value: 1486.024 mgc: 2.165899e-07 ustep: 1 
#> iter: 4  mgc: 1.058607e-13 
#> iter: 1  mgc: 1.058607e-13 
#> outer mgc:  17.09293 
#> iter: 1  value: 1485.61 mgc: 1.285076 ustep: 1 
#> iter: 2  value: 1485.61 mgc: 0.003289757 ustep: 1 
#> iter: 3  value: 1485.61 mgc: 1.77962e-07 ustep: 1 
#> iter: 4  mgc: 1.027617e-13 
#> iter: 1  mgc: 1.027617e-13 
#> outer mgc:  11.28409 
#> iter: 1  value: 1485.202 mgc: 0.9778893 ustep: 1 
#> iter: 2  value: 1485.202 mgc: 0.001687637 ustep: 1 
#> iter: 3  value: 1485.202 mgc: 1.842371e-07 ustep: 1 
#> iter: 4  mgc: 1.610835e-13 
#> iter: 1  mgc: 1.610835e-13 
#> outer mgc:  12.03234 
#> iter: 1  value: 1484.93 mgc: 0.4586302 ustep: 1 
#> iter: 2  value: 1484.93 mgc: 0.0005143759 ustep: 1 
#> iter: 3  mgc: 4.991647e-09 
#> iter: 1  value: 1484.637 mgc: 0.4211854 ustep: 1 
#> iter: 2  value: 1484.637 mgc: 0.000337214 ustep: 1 
#> iter: 3  mgc: 3.50325e-09 
#> iter: 1  value: 1484.027 mgc: 1.164583 ustep: 1 
#> iter: 2  value: 1484.027 mgc: 0.002528056 ustep: 1 
#> iter: 3  value: 1484.027 mgc: 1.969677e-07 ustep: 1 
#> iter: 4  mgc: 1.129136e-13 
#> iter: 1  mgc: 1.129136e-13 
#> outer mgc:  13.75158 
#> iter: 1  value: 1484.146 mgc: 5.420393 ustep: 1 
#> iter: 2  value: 1484.146 mgc: 0.0231164 ustep: 1 
#> iter: 3  value: 1484.146 mgc: 3.373091e-05 ustep: 1 
#> iter: 4  mgc: 8.656349e-11 
#> iter: 1  value: 1483.676 mgc: 3.280788 ustep: 1 
#> iter: 2  value: 1483.676 mgc: 0.01099559 ustep: 1 
#> iter: 3  value: 1483.676 mgc: 7.756518e-06 ustep: 1 
#> iter: 4  mgc: 3.677608e-12 
#> iter: 1  mgc: 3.677608e-12 
#> outer mgc:  26.22733 
#> iter: 1  value: 1483.022 mgc: 1.498191 ustep: 1 
#> iter: 2  value: 1483.022 mgc: 0.005686693 ustep: 1 
#> iter: 3  value: 1483.022 mgc: 3.931097e-07 ustep: 1 
#> iter: 4  mgc: 1.352318e-13 
#> iter: 1  mgc: 1.352318e-13 
#> outer mgc:  15.74363 
#> iter: 1  value: 1483.381 mgc: 2.781538 ustep: 1 
#> iter: 2  value: 1483.381 mgc: 0.007979879 ustep: 1 
#> iter: 3  value: 1483.381 mgc: 7.408821e-06 ustep: 1 
#> iter: 4  mgc: 5.50153e-12 
#> iter: 1  value: 1482.875 mgc: 1.694261 ustep: 1 
#> iter: 2  value: 1482.875 mgc: 0.003688816 ustep: 1 
#> iter: 3  value: 1482.875 mgc: 8.37771e-07 ustep: 1 
#> iter: 4  mgc: 1.284404e-13 
#> iter: 1  mgc: 1.284404e-13 
#> outer mgc:  19.99543 
#> iter: 1  value: 1482.667 mgc: 1.22866 ustep: 1 
#> iter: 2  value: 1482.667 mgc: 0.00236316 ustep: 1 
#> iter: 3  value: 1482.667 mgc: 1.480278e-07 ustep: 1 
#> iter: 4  mgc: 1.298904e-13 
#> iter: 1  mgc: 1.298904e-13 
#> outer mgc:  3.725176 
#> iter: 1  value: 1482.576 mgc: 0.422281 ustep: 1 
#> iter: 2  value: 1482.576 mgc: 0.000195157 ustep: 1 
#> iter: 3  mgc: 2.51076e-09 
#> iter: 1  mgc: 2.51076e-09 
#> outer mgc:  4.649527 
#> iter: 1  value: 1482.54 mgc: 0.449115 ustep: 1 
#> iter: 2  value: 1482.54 mgc: 0.0003309445 ustep: 1 
#> iter: 3  mgc: 5.090264e-09 
#> iter: 1  mgc: 5.090264e-09 
#> outer mgc:  3.97857 
#> iter: 1  value: 1482.463 mgc: 0.4806927 ustep: 1 
#> iter: 2  value: 1482.463 mgc: 0.0003100809 ustep: 1 
#> iter: 3  mgc: 6.020296e-09 
#> iter: 1  value: 1482.427 mgc: 0.5187965 ustep: 1 
#> iter: 2  value: 1482.427 mgc: 0.0009785313 ustep: 1 
#> iter: 3  value: 1482.427 mgc: 1.176203e-08 ustep: 1 
#> mgc: 1.48897e-13 
#> iter: 1  mgc: 1.48897e-13 
#> outer mgc:  3.620662 
#> iter: 1  value: 1489.06 mgc: 8.350134 ustep: 1 
#> iter: 2  value: 1489.057 mgc: 0.1084577 ustep: 1 
#> iter: 3  value: 1489.057 mgc: 0.0005139248 ustep: 1 
#> iter: 4  value: 1489.057 mgc: 1.562475e-08 ustep: 1 
#> mgc: 1.349131e-13 
#> iter: 1  value: 1482.452 mgc: 1.15392 ustep: 1 
#> iter: 2  value: 1482.452 mgc: 0.001974495 ustep: 1 
#> iter: 3  value: 1482.452 mgc: 1.405868e-07 ustep: 1 
#> iter: 4  mgc: 1.409349e-13 
#> iter: 1  value: 1482.402 mgc: 0.3272833 ustep: 1 
#> iter: 2  value: 1482.402 mgc: 0.0001570601 ustep: 1 
#> iter: 3  mgc: 8.819405e-10 
#> iter: 1  mgc: 8.819405e-10 
#> outer mgc:  4.673897 
#> iter: 1  value: 1482.398 mgc: 0.2348311 ustep: 1 
#> iter: 2  value: 1482.398 mgc: 8.995916e-05 ustep: 1 
#> iter: 3  mgc: 3.176147e-10 
#> iter: 1  mgc: 3.176147e-10 
#> outer mgc:  3.457985 
#> iter: 1  value: 1482.359 mgc: 0.3072644 ustep: 1 
#> iter: 2  value: 1482.359 mgc: 0.000123359 ustep: 1 
#> iter: 3  mgc: 6.918811e-10 
#> iter: 1  mgc: 6.918811e-10 
#> outer mgc:  4.284231 
#> iter: 1  value: 1482.319 mgc: 0.1970246 ustep: 1 
#> iter: 2  value: 1482.319 mgc: 0.0001013927 ustep: 1 
#> iter: 3  mgc: 2.109849e-10 
#> iter: 1  mgc: 2.109849e-10 
#> outer mgc:  3.417262 
#> iter: 1  value: 1482.231 mgc: 0.3599432 ustep: 1 
#> iter: 2  value: 1482.231 mgc: 0.0001463267 ustep: 1 
#> iter: 3  mgc: 1.030634e-09 
#> iter: 1  value: 1482.139 mgc: 0.2604523 ustep: 1 
#> iter: 2  value: 1482.139 mgc: 0.0001057136 ustep: 1 
#> iter: 3  mgc: 8.250267e-11 
#> iter: 1  value: 1481.95 mgc: 0.6761641 ustep: 1 
#> iter: 2  value: 1481.95 mgc: 0.0007108039 ustep: 1 
#> iter: 3  mgc: 3.737436e-09 
#> iter: 1  mgc: 3.737436e-09 
#> outer mgc:  11.96194 
#> iter: 1  value: 1481.95 mgc: 0.9938168 ustep: 1 
#> iter: 2  value: 1481.95 mgc: 0.001656523 ustep: 1 
#> iter: 3  value: 1481.95 mgc: 7.88273e-08 ustep: 1 
#> iter: 4  mgc: 1.778807e-13 
#> iter: 1  mgc: 1.778807e-13 
#> outer mgc:  9.51638 
#> iter: 1  value: 1481.621 mgc: 2.150961 ustep: 1 
#> iter: 2  value: 1481.621 mgc: 0.003525262 ustep: 1 
#> iter: 3  value: 1481.621 mgc: 1.380718e-06 ustep: 1 
#> iter: 4  mgc: 1.782529e-13 
#> iter: 1  mgc: 1.782529e-13 
#> outer mgc:  10.57905 
#> iter: 1  value: 1481.362 mgc: 1.904292 ustep: 1 
#> iter: 2  value: 1481.362 mgc: 0.002636543 ustep: 1 
#> iter: 3  value: 1481.362 mgc: 2.454071e-07 ustep: 1 
#> iter: 4  mgc: 1.437481e-13 
#> iter: 1  mgc: 1.437481e-13 
#> outer mgc:  14.73725 
#> iter: 1  value: 1481.952 mgc: 1.933435 ustep: 1 
#> iter: 2  value: 1481.952 mgc: 0.007563283 ustep: 1 
#> iter: 3  value: 1481.952 mgc: 2.906435e-06 ustep: 1 
#> iter: 4  mgc: 4.729706e-13 
#> iter: 1  mgc: 4.729706e-13 
#> outer mgc:  11.92185 
#> iter: 1  value: 1481.812 mgc: 2.657457 ustep: 1 
#> iter: 2  value: 1481.812 mgc: 0.01311171 ustep: 1 
#> iter: 3  value: 1481.812 mgc: 2.133853e-06 ustep: 1 
#> iter: 4  mgc: 1.102591e-13 
#> iter: 1  mgc: 1.102591e-13 
#> outer mgc:  11.23098 
#> iter: 1  value: 1480.313 mgc: 2.449264 ustep: 1 
#> iter: 2  value: 1480.313 mgc: 0.007835328 ustep: 1 
#> iter: 3  value: 1480.313 mgc: 3.198329e-06 ustep: 1 
#> iter: 4  mgc: 1.203834e-12 
#> iter: 1  mgc: 1.203834e-12 
#> outer mgc:  11.46321 
#> iter: 1  value: 1481.611 mgc: 0.9543484 ustep: 1 
#> iter: 2  value: 1481.61 mgc: 0.005366925 ustep: 1 
#> iter: 3  value: 1481.61 mgc: 1.761968e-06 ustep: 1 
#> iter: 4  mgc: 4.411874e-13 
#> iter: 1  mgc: 4.411874e-13 
#> outer mgc:  9.527116 
#> iter: 1  value: 1481.333 mgc: 1.600286 ustep: 1 
#> iter: 2  value: 1481.333 mgc: 0.006601822 ustep: 1 
#> iter: 3  value: 1481.333 mgc: 7.607407e-07 ustep: 1 
#> iter: 4  mgc: 1.217217e-13 
#> iter: 1  mgc: 1.217217e-13 
#> outer mgc:  15.13245 
#> iter: 1  value: 1481.218 mgc: 1.605235 ustep: 1 
#> iter: 2  value: 1481.218 mgc: 0.004609048 ustep: 1 
#> iter: 3  value: 1481.218 mgc: 3.277141e-07 ustep: 1 
#> iter: 4  mgc: 1.353497e-13 
#> iter: 1  value: 1481.23 mgc: 1.516237 ustep: 1 
#> iter: 2  value: 1481.23 mgc: 0.00384419 ustep: 1 
#> iter: 3  value: 1481.23 mgc: 3.082671e-07 ustep: 1 
#> iter: 4  mgc: 1.331155e-13 
#> iter: 1  mgc: 1.331155e-13 
#> outer mgc:  18.75247 
#> iter: 1  value: 1481.185 mgc: 1.632182 ustep: 1 
#> iter: 2  value: 1481.185 mgc: 0.003020201 ustep: 1 
#> iter: 3  value: 1481.185 mgc: 2.999233e-07 ustep: 1 
#> iter: 4  mgc: 1.330147e-13 
#> iter: 1  mgc: 1.330147e-13 
#> outer mgc:  4.827255 
#> iter: 1  value: 1480.576 mgc: 0.6956258 ustep: 1 
#> iter: 2  value: 1480.576 mgc: 0.001355706 ustep: 1 
#> iter: 3  mgc: 9.431534e-09 
#> iter: 1  mgc: 9.431534e-09 
#> outer mgc:  2.874137 
#> iter: 1  value: 1481.976 mgc: 3.722496 ustep: 1 
#> iter: 2  value: 1481.975 mgc: 0.01970957 ustep: 1 
#> iter: 3  value: 1481.975 mgc: 6.72158e-06 ustep: 1 
#> iter: 4  mgc: 3.123952e-12 
#> iter: 1  value: 1480.593 mgc: 0.5740801 ustep: 1 
#> iter: 2  value: 1480.593 mgc: 0.0004930071 ustep: 1 
#> iter: 3  mgc: 5.607087e-09 
#> iter: 1  value: 1480.576 mgc: 0.2587967 ustep: 1 
#> iter: 2  value: 1480.576 mgc: 0.0001005901 ustep: 1 
#> iter: 3  mgc: 2.343888e-10 
#> iter: 1  mgc: 2.343888e-10 
#> outer mgc:  3.886759 
#> iter: 1  value: 1480.565 mgc: 0.1999935 ustep: 1 
#> iter: 2  value: 1480.565 mgc: 5.754166e-05 ustep: 1 
#> iter: 3  mgc: 1.45752e-10 
#> iter: 1  mgc: 1.45752e-10 
#> outer mgc:  2.023225 
#> iter: 1  value: 1480.575 mgc: 0.1860483 ustep: 1 
#> iter: 2  value: 1480.575 mgc: 4.213844e-05 ustep: 1 
#> iter: 3  mgc: 5.152712e-11 
#> iter: 1  mgc: 5.152712e-11 
#> outer mgc:  2.835513 
#> iter: 1  value: 1480.576 mgc: 0.1380132 ustep: 1 
#> iter: 2  value: 1480.576 mgc: 4.744316e-05 ustep: 1 
#> iter: 3  mgc: 1.368847e-10 
#> iter: 1  mgc: 1.368847e-10 
#> outer mgc:  2.008524 
#> iter: 1  value: 1480.544 mgc: 0.36791 ustep: 1 
#> iter: 2  value: 1480.544 mgc: 0.0001678498 ustep: 1 
#> iter: 3  mgc: 2.985776e-10 
#> iter: 1  value: 1480.51 mgc: 0.2516838 ustep: 1 
#> iter: 2  value: 1480.51 mgc: 7.99369e-05 ustep: 1 
#> iter: 3  mgc: 9.980118e-11 
#> iter: 1  mgc: 9.980118e-11 
#> outer mgc:  2.151082 
#> iter: 1  value: 1480.62 mgc: 0.5651091 ustep: 1 
#> iter: 2  value: 1480.62 mgc: 0.0005465567 ustep: 1 
#> iter: 3  value: 1480.62 mgc: 1.776454e-08 ustep: 1 
#> mgc: 1.864867e-13 
#> iter: 1  value: 1480.786 mgc: 0.4102145 ustep: 1 
#> iter: 2  value: 1480.786 mgc: 0.0003042526 ustep: 1 
#> iter: 3  mgc: 5.349851e-09 
#> iter: 1  mgc: 5.349851e-09 
#> outer mgc:  4.093134 
#> iter: 1  value: 1480.755 mgc: 2.483024 ustep: 1 
#> iter: 2  value: 1480.755 mgc: 0.007739861 ustep: 1 
#> iter: 3  value: 1480.755 mgc: 2.676646e-06 ustep: 1 
#> iter: 4  mgc: 3.332549e-13 
#> iter: 1  value: 1481.063 mgc: 1.653116 ustep: 1 
#> iter: 2  value: 1481.063 mgc: 0.00306698 ustep: 1 
#> iter: 3  value: 1481.063 mgc: 3.953225e-07 ustep: 1 
#> iter: 4  mgc: 1.465422e-13 
#> iter: 1  value: 1480.766 mgc: 0.1950136 ustep: 1 
#> iter: 2  value: 1480.766 mgc: 4.467749e-05 ustep: 1 
#> iter: 3  mgc: 7.964628e-11 
#> iter: 1  mgc: 7.964628e-11 
#> outer mgc:  2.454487 
#> iter: 1  value: 1480.759 mgc: 0.2026971 ustep: 1 
#> iter: 2  value: 1480.759 mgc: 5.631357e-05 ustep: 1 
#> iter: 3  mgc: 7.954128e-11 
#> iter: 1  mgc: 7.954128e-11 
#> outer mgc:  3.012803 
#> iter: 1  value: 1480.733 mgc: 0.1581135 ustep: 1 
#> iter: 2  value: 1480.733 mgc: 2.993649e-05 ustep: 1 
#> iter: 3  mgc: 4.057499e-11 
#> iter: 1  mgc: 4.057499e-11 
#> outer mgc:  1.915349 
#> iter: 1  value: 1480.707 mgc: 0.1524483 ustep: 1 
#> iter: 2  value: 1480.707 mgc: 4.795877e-05 ustep: 1 
#> iter: 3  mgc: 2.826216e-11 
#> iter: 1  mgc: 2.826216e-11 
#> outer mgc:  2.782223 
#> iter: 1  value: 1480.662 mgc: 0.1685434 ustep: 1 
#> iter: 2  value: 1480.662 mgc: 6.731408e-05 ustep: 1 
#> iter: 3  mgc: 1.972486e-10 
#> iter: 1  value: 1480.619 mgc: 0.09979414 ustep: 1 
#> iter: 2  value: 1480.619 mgc: 4.909044e-05 ustep: 1 
#> iter: 3  mgc: 5.421083e-11 
#> iter: 1  value: 1480.504 mgc: 0.2928363 ustep: 1 
#> iter: 2  value: 1480.504 mgc: 0.0004239851 ustep: 1 
#> iter: 3  mgc: 4.04006e-09 
#> iter: 1  mgc: 4.04006e-09 
#> outer mgc:  2.279047 
#> iter: 1  value: 1480.408 mgc: 0.4391415 ustep: 1 
#> iter: 2  value: 1480.408 mgc: 0.0005124681 ustep: 1 
#> iter: 3  mgc: 8.064346e-09 
#> iter: 1  mgc: 8.064346e-09 
#> outer mgc:  3.424182 
#> iter: 1  value: 1480.325 mgc: 1.019313 ustep: 1 
#> iter: 2  value: 1480.325 mgc: 0.00255479 ustep: 1 
#> iter: 3  value: 1480.325 mgc: 2.028547e-07 ustep: 1 
#> iter: 4  mgc: 1.371914e-13 
#> iter: 1  mgc: 1.371914e-13 
#> outer mgc:  6.417647 
#> iter: 1  value: 1480.54 mgc: 1.050054 ustep: 1 
#> iter: 2  value: 1480.54 mgc: 0.002752835 ustep: 1 
#> iter: 3  value: 1480.54 mgc: 9.329752e-08 ustep: 1 
#> iter: 4  mgc: 1.736078e-13 
#> iter: 1  mgc: 1.736078e-13 
#> outer mgc:  5.551258 
#> iter: 1  value: 1480.807 mgc: 0.6917315 ustep: 1 
#> iter: 2  value: 1480.807 mgc: 0.003896867 ustep: 1 
#> iter: 3  value: 1480.807 mgc: 5.812306e-07 ustep: 1 
#> iter: 4  mgc: 1.03765e-13 
#> iter: 1  value: 1480.633 mgc: 0.5127942 ustep: 1 
#> iter: 2  value: 1480.633 mgc: 0.001479467 ustep: 1 
#> iter: 3  value: 1480.633 mgc: 8.568943e-08 ustep: 1 
#> iter: 4  mgc: 2.245871e-13 
#> iter: 1  mgc: 2.245871e-13 
#> outer mgc:  3.988956 
#> iter: 1  value: 1480.692 mgc: 0.4443249 ustep: 1 
#> iter: 2  value: 1480.692 mgc: 0.0003523136 ustep: 1 
#> iter: 3  mgc: 2.401408e-09 
#> iter: 1  mgc: 2.401408e-09 
#> outer mgc:  1.745591 
#> iter: 1  value: 1480.327 mgc: 1.344813 ustep: 1 
#> iter: 2  value: 1480.327 mgc: 0.002130243 ustep: 1 
#> iter: 3  value: 1480.327 mgc: 2.952014e-07 ustep: 1 
#> iter: 4  mgc: 1.597086e-13 
#> iter: 1  value: 1480.591 mgc: 0.2350093 ustep: 1 
#> iter: 2  value: 1480.591 mgc: 6.388907e-05 ustep: 1 
#> iter: 3  mgc: 2.679771e-10 
#> iter: 1  mgc: 2.679771e-10 
#> outer mgc:  1.258222 
#> iter: 1  value: 1480.543 mgc: 0.3892988 ustep: 1 
#> iter: 2  value: 1480.543 mgc: 0.0002472189 ustep: 1 
#> iter: 3  mgc: 9.803654e-10 
#> iter: 1  mgc: 9.803654e-10 
#> outer mgc:  3.260739 
#> iter: 1  value: 1480.518 mgc: 0.2873954 ustep: 1 
#> iter: 2  value: 1480.518 mgc: 8.757387e-05 ustep: 1 
#> iter: 3  mgc: 5.105365e-10 
#> iter: 1  mgc: 5.105365e-10 
#> outer mgc:  2.146938 
#> iter: 1  value: 1480.513 mgc: 0.07480507 ustep: 1 
#> iter: 2  value: 1480.513 mgc: 1.35821e-05 ustep: 1 
#> iter: 3  mgc: 3.451259e-12 
#> iter: 1  mgc: 3.451259e-12 
#> outer mgc:  1.044692 
#> iter: 1  value: 1480.494 mgc: 0.1704351 ustep: 1 
#> iter: 2  value: 1480.494 mgc: 8.31167e-05 ustep: 1 
#> iter: 3  mgc: 3.283791e-10 
#> iter: 1  value: 1480.499 mgc: 0.08107152 ustep: 1 
#> iter: 2  value: 1480.499 mgc: 3.317422e-05 ustep: 1 
#> iter: 3  mgc: 5.373237e-11 
#> iter: 1  value: 1480.519 mgc: 0.2325496 ustep: 1 
#> iter: 2  value: 1480.519 mgc: 0.0002740008 ustep: 1 
#> iter: 3  mgc: 3.675374e-09 
#> iter: 1  mgc: 3.675374e-09 
#> outer mgc:  1.598186 
#> iter: 1  value: 1480.498 mgc: 0.6971712 ustep: 1 
#> iter: 2  value: 1480.498 mgc: 0.0005026323 ustep: 1 
#> iter: 3  mgc: 6.956094e-09 
#> iter: 1  mgc: 6.956094e-09 
#> outer mgc:  1.056585 
#> iter: 1  value: 1480.338 mgc: 0.2742891 ustep: 1 
#> iter: 2  value: 1480.338 mgc: 0.001014891 ustep: 1 
#> iter: 3  value: 1480.338 mgc: 3.158804e-08 ustep: 1 
#> mgc: 8.628544e-14 
#> iter: 1  mgc: 8.628544e-14 
#> outer mgc:  1.933678 
#> iter: 1  value: 1480.389 mgc: 0.9881291 ustep: 1 
#> iter: 2  value: 1480.389 mgc: 0.00189766 ustep: 1 
#> iter: 3  value: 1480.389 mgc: 6.730747e-08 ustep: 1 
#> iter: 4  mgc: 1.128162e-13 
#> iter: 1  value: 1480.349 mgc: 0.5230625 ustep: 1 
#> iter: 2  value: 1480.349 mgc: 0.0004888137 ustep: 1 
#> iter: 3  mgc: 4.756848e-09 
#> iter: 1  mgc: 4.756848e-09 
#> outer mgc:  4.132387 
#> iter: 1  value: 1480.347 mgc: 0.7534395 ustep: 1 
#> iter: 2  value: 1480.347 mgc: 0.0007450596 ustep: 1 
#> iter: 3  value: 1480.347 mgc: 2.265578e-08 ustep: 1 
#> mgc: 1.947617e-13 
#> iter: 1  mgc: 1.947617e-13 
#> outer mgc:  3.680418 
#> iter: 1  value: 1480.302 mgc: 0.4155585 ustep: 1 
#> iter: 2  value: 1480.302 mgc: 0.000622128 ustep: 1 
#> iter: 3  value: 1480.302 mgc: 1.631258e-08 ustep: 1 
#> mgc: 1.199613e-13 
#> iter: 1  mgc: 1.199613e-13 
#> outer mgc:  1.204428 
#> iter: 1  value: 1480.226 mgc: 0.2740807 ustep: 1 
#> iter: 2  value: 1480.226 mgc: 0.0001274734 ustep: 1 
#> iter: 3  mgc: 1.009286e-09 
#> iter: 1  mgc: 1.009286e-09 
#> outer mgc:  1.029404 
#> iter: 1  value: 1480.23 mgc: 0.2152312 ustep: 1 
#> iter: 2  value: 1480.23 mgc: 0.0003420824 ustep: 1 
#> iter: 3  mgc: 2.184203e-09 
#> iter: 1  mgc: 2.184203e-09 
#> outer mgc:  1.804884 
#> iter: 1  value: 1480.269 mgc: 0.4933803 ustep: 1 
#> iter: 2  value: 1480.269 mgc: 0.0002760244 ustep: 1 
#> iter: 3  mgc: 2.592163e-09 
#> iter: 1  mgc: 2.592163e-09 
#> outer mgc:  1.159007 
#> iter: 1  value: 1480.197 mgc: 0.2514098 ustep: 1 
#> iter: 2  value: 1480.197 mgc: 0.0003831805 ustep: 1 
#> iter: 3  mgc: 2.766063e-09 
#> iter: 1  mgc: 2.766063e-09 
#> outer mgc:  2.43547 
#> iter: 1  value: 1480.122 mgc: 0.4195116 ustep: 1 
#> iter: 2  value: 1480.122 mgc: 0.0005821322 ustep: 1 
#> iter: 3  value: 1480.122 mgc: 1.676527e-08 ustep: 1 
#> mgc: 1.606811e-13 
#> iter: 1  value: 1480.154 mgc: 0.3517261 ustep: 1 
#> iter: 2  value: 1480.154 mgc: 0.0002856099 ustep: 1 
#> iter: 3  mgc: 4.118976e-09 
#> iter: 1  mgc: 4.118976e-09 
#> outer mgc:  3.610521 
#> iter: 1  value: 1480.178 mgc: 0.226887 ustep: 1 
#> iter: 2  value: 1480.178 mgc: 0.000144511 ustep: 1 
#> iter: 3  mgc: 2.602891e-10 
#> iter: 1  mgc: 2.602891e-10 
#> outer mgc:  0.4416795 
#> iter: 1  value: 1480.226 mgc: 0.203553 ustep: 1 
#> iter: 2  value: 1480.226 mgc: 4.872469e-05 ustep: 1 
#> iter: 3  mgc: 8.550415e-11 
#> iter: 1  mgc: 8.550415e-11 
#> outer mgc:  0.3958739 
#> iter: 1  value: 1480.207 mgc: 0.1964812 ustep: 1 
#> iter: 2  value: 1480.207 mgc: 7.627051e-05 ustep: 1 
#> iter: 3  mgc: 1.098879e-10 
#> iter: 1  mgc: 1.098879e-10 
#> outer mgc:  0.6032935 
#> iter: 1  value: 1480.201 mgc: 0.2028428 ustep: 1 
#> iter: 2  value: 1480.201 mgc: 0.0002156831 ustep: 1 
#> iter: 3  mgc: 2.391098e-09 
#> iter: 1  mgc: 2.391098e-09 
#> outer mgc:  0.4713319 
#> iter: 1  value: 1480.213 mgc: 0.1721886 ustep: 1 
#> iter: 2  value: 1480.213 mgc: 0.0001761286 ustep: 1 
#> iter: 3  mgc: 6.782917e-10 
#> iter: 1  mgc: 6.782917e-10 
#> outer mgc:  0.4841169 
#> iter: 1  value: 1480.208 mgc: 0.2462303 ustep: 1 
#> iter: 2  value: 1480.208 mgc: 0.000103184 ustep: 1 
#> iter: 3  mgc: 2.549395e-11 
#> iter: 1  mgc: 2.549395e-11 
#> outer mgc:  1.71064 
#> iter: 1  value: 1480.191 mgc: 0.2402307 ustep: 1 
#> iter: 2  value: 1480.191 mgc: 8.001818e-05 ustep: 1 
#> iter: 3  mgc: 3.881903e-10 
#> iter: 1  mgc: 3.881903e-10 
#> outer mgc:  0.4078351 
#> iter: 1  value: 1480.188 mgc: 0.2767952 ustep: 1 
#> iter: 2  value: 1480.188 mgc: 7.847871e-05 ustep: 1 
#> iter: 3  mgc: 1.708479e-10 
#> iter: 1  mgc: 1.708479e-10 
#> outer mgc:  0.3009992 
#> iter: 1  value: 1480.196 mgc: 0.1484601 ustep: 1 
#> iter: 2  value: 1480.196 mgc: 6.238203e-05 ustep: 1 
#> iter: 3  mgc: 3.96794e-11 
#> iter: 1  mgc: 3.96794e-11 
#> outer mgc:  0.4575057 
#> iter: 1  value: 1480.204 mgc: 0.2061473 ustep: 1 
#> iter: 2  value: 1480.204 mgc: 0.0001149949 ustep: 1 
#> iter: 3  mgc: 4.020413e-11 
#> iter: 1  mgc: 4.020413e-11 
#> outer mgc:  1.091744 
#> iter: 1  value: 1480.195 mgc: 0.2067316 ustep: 1 
#> iter: 2  value: 1480.195 mgc: 3.846943e-05 ustep: 1 
#> iter: 3  mgc: 1.836407e-11 
#> iter: 1  mgc: 1.836407e-11 
#> outer mgc:  0.9720897 
#> iter: 1  value: 1480.198 mgc: 0.1652615 ustep: 1 
#> iter: 2  value: 1480.198 mgc: 6.038532e-05 ustep: 1 
#> iter: 3  mgc: 1.106912e-11 
#> iter: 1  mgc: 1.106912e-11 
#> outer mgc:  0.1879934 
#> iter: 1  value: 1480.178 mgc: 0.1170643 ustep: 1 
#> iter: 2  value: 1480.178 mgc: 2.022839e-05 ustep: 1 
#> iter: 3  mgc: 1.342885e-12 
#> iter: 1  mgc: 1.342885e-12 
#> outer mgc:  0.7445399 
#> iter: 1  value: 1480.203 mgc: 0.07742144 ustep: 1 
#> iter: 2  value: 1480.203 mgc: 4.295285e-05 ustep: 1 
#> iter: 3  mgc: 3.999078e-11 
#> iter: 1  mgc: 3.999078e-11 
#> outer mgc:  0.5120758 
#> iter: 1  value: 1480.204 mgc: 0.1328167 ustep: 1 
#> iter: 2  value: 1480.204 mgc: 2.492835e-05 ustep: 1 
#> iter: 3  mgc: 1.572703e-11 
#> iter: 1  mgc: 1.572703e-11 
#> outer mgc:  0.4812711 
#> iter: 1  value: 1480.205 mgc: 0.09325414 ustep: 1 
#> iter: 2  value: 1480.205 mgc: 2.610858e-05 ustep: 1 
#> iter: 3  mgc: 2.755143e-11 
#> iter: 1  mgc: 2.755143e-11 
#> outer mgc:  0.8162939 
#> iter: 1  value: 1480.193 mgc: 0.09120275 ustep: 1 
#> iter: 2  value: 1480.193 mgc: 1.221221e-05 ustep: 1 
#> iter: 3  mgc: 2.938942e-12 
#> iter: 1  mgc: 2.938942e-12 
#> outer mgc:  0.4043535 
#> iter: 1  value: 1480.187 mgc: 0.06669411 ustep: 1 
#> iter: 2  value: 1480.187 mgc: 7.594693e-06 ustep: 1 
#> iter: 3  mgc: 1.495871e-12 
#> iter: 1  value: 1480.172 mgc: 0.2059195 ustep: 1 
#> iter: 2  value: 1480.172 mgc: 3.563519e-05 ustep: 1 
#> iter: 3  mgc: 2.433375e-11 
#> iter: 1  mgc: 2.433375e-11 
#> outer mgc:  0.6133419 
#> iter: 1  value: 1480.141 mgc: 0.1512849 ustep: 1 
#> iter: 2  value: 1480.141 mgc: 0.0001595206 ustep: 1 
#> iter: 3  mgc: 6.206321e-10 
#> iter: 1  mgc: 6.206321e-10 
#> outer mgc:  0.6377674 
#> iter: 1  value: 1480.13 mgc: 0.1740391 ustep: 1 
#> iter: 2  value: 1480.13 mgc: 5.815903e-05 ustep: 1 
#> iter: 3  mgc: 1.235872e-10 
#> iter: 1  mgc: 1.235872e-10 
#> outer mgc:  0.3357017 
#> iter: 1  value: 1480.165 mgc: 0.2430705 ustep: 1 
#> iter: 2  value: 1480.165 mgc: 0.0002332536 ustep: 1 
#> iter: 3  mgc: 3.572448e-09 
#> iter: 1  mgc: 3.572448e-09 
#> outer mgc:  1.581374 
#> iter: 1  value: 1480.174 mgc: 0.1043481 ustep: 1 
#> iter: 2  value: 1480.174 mgc: 0.0001309403 ustep: 1 
#> iter: 3  mgc: 3.971396e-10 
#> iter: 1  mgc: 3.971396e-10 
#> outer mgc:  0.9085527 
#> iter: 1  value: 1480.163 mgc: 0.162355 ustep: 1 
#> iter: 2  value: 1480.163 mgc: 0.0002706277 ustep: 1 
#> iter: 3  mgc: 2.11849e-09 
#> iter: 1  mgc: 2.11849e-09 
#> outer mgc:  0.6737985 
#> iter: 1  value: 1480.147 mgc: 0.1632254 ustep: 1 
#> iter: 2  value: 1480.147 mgc: 0.0005510305 ustep: 1 
#> iter: 3  mgc: 6.810211e-09 
#> iter: 1  mgc: 6.810211e-09 
#> outer mgc:  0.6379296 
#> iter: 1  value: 1480.153 mgc: 0.2961021 ustep: 1 
#> iter: 2  value: 1480.153 mgc: 0.0002593853 ustep: 1 
#> iter: 3  mgc: 1.859187e-09 
#> iter: 1  mgc: 1.859187e-09 
#> outer mgc:  1.648297 
#> iter: 1  value: 1480.148 mgc: 0.1959037 ustep: 1 
#> iter: 2  value: 1480.148 mgc: 0.0001387879 ustep: 1 
#> iter: 3  mgc: 4.712229e-10 
#> iter: 1  mgc: 4.712229e-10 
#> outer mgc:  1.392841 
#> iter: 1  value: 1480.188 mgc: 0.1841252 ustep: 1 
#> iter: 2  value: 1480.188 mgc: 3.304829e-05 ustep: 1 
#> iter: 3  mgc: 5.683794e-11 
#> iter: 1  mgc: 5.683794e-11 
#> outer mgc:  0.6816964 
#> iter: 1  value: 1480.159 mgc: 0.2443019 ustep: 1 
#> iter: 2  value: 1480.159 mgc: 0.000344305 ustep: 1 
#> iter: 3  mgc: 2.700598e-09 
#> iter: 1  value: 1480.174 mgc: 0.1288384 ustep: 1 
#> iter: 2  value: 1480.174 mgc: 6.004442e-05 ustep: 1 
#> iter: 3  mgc: 8.181479e-11 
#> iter: 1  mgc: 8.181479e-11 
#> outer mgc:  0.6221611 
#> iter: 1  value: 1480.178 mgc: 0.1072898 ustep: 1 
#> iter: 2  value: 1480.178 mgc: 2.914791e-05 ustep: 1 
#> iter: 3  mgc: 1.816793e-11 
#> iter: 1  mgc: 1.816793e-11 
#> outer mgc:  0.8607478 
#> iter: 1  value: 1480.153 mgc: 0.1265159 ustep: 1 
#> iter: 2  value: 1480.153 mgc: 5.387563e-05 ustep: 1 
#> iter: 3  mgc: 1.339732e-10 
#> iter: 1  mgc: 1.339732e-10 
#> outer mgc:  0.247591 
#> iter: 1  value: 1480.162 mgc: 0.09867355 ustep: 1 
#> iter: 2  value: 1480.162 mgc: 7.909146e-05 ustep: 1 
#> iter: 3  mgc: 1.664445e-10 
#> iter: 1  mgc: 1.664445e-10 
#> outer mgc:  0.5566738 
#> iter: 1  value: 1480.17 mgc: 0.07711617 ustep: 1 
#> iter: 2  value: 1480.17 mgc: 2.938788e-05 ustep: 1 
#> iter: 3  mgc: 1.672164e-11 
#> iter: 1  mgc: 1.672164e-11 
#> outer mgc:  0.1527908 
#> iter: 1  value: 1480.152 mgc: 0.03309766 ustep: 1 
#> iter: 2  value: 1480.152 mgc: 2.512237e-05 ustep: 1 
#> iter: 3  mgc: 1.442793e-11 
#> iter: 1  value: 1480.128 mgc: 0.06115799 ustep: 1 
#> iter: 2  value: 1480.128 mgc: 5.311882e-05 ustep: 1 
#> iter: 3  mgc: 6.462279e-11 
#> iter: 1  mgc: 6.462279e-11 
#> outer mgc:  0.4390389 
#> iter: 1  value: 1480.116 mgc: 0.113141 ustep: 1 
#> iter: 2  value: 1480.116 mgc: 0.0002151486 ustep: 1 
#> iter: 3  mgc: 1.026093e-09 
#> iter: 1  mgc: 1.026093e-09 
#> outer mgc:  0.2281602 
#> iter: 1  value: 1480.125 mgc: 0.192936 ustep: 1 
#> iter: 2  value: 1480.125 mgc: 0.0002876924 ustep: 1 
#> iter: 3  mgc: 2.69121e-09 
#> iter: 1  mgc: 2.69121e-09 
#> outer mgc:  0.2224913 
#> iter: 1  value: 1480.158 mgc: 0.1113079 ustep: 1 
#> iter: 2  value: 1480.158 mgc: 0.0001345339 ustep: 1 
#> iter: 3  mgc: 4.671088e-10 
#> iter: 1  mgc: 4.671088e-10 
#> outer mgc:  0.3220598 
#> iter: 1  value: 1480.189 mgc: 0.1340004 ustep: 1 
#> iter: 2  value: 1480.189 mgc: 7.070735e-05 ustep: 1 
#> iter: 3  mgc: 2.256272e-10 
#> iter: 1  value: 1480.173 mgc: 0.08587348 ustep: 1 
#> iter: 2  value: 1480.173 mgc: 1.875345e-05 ustep: 1 
#> iter: 3  mgc: 9.428669e-12 
#> iter: 1  mgc: 9.428669e-12 
#> outer mgc:  0.3473143 
#> iter: 1  value: 1480.16 mgc: 0.07011624 ustep: 1 
#> iter: 2  value: 1480.16 mgc: 9.092769e-06 ustep: 1 
#> iter: 3  mgc: 1.756703e-12 
#> iter: 1  mgc: 1.756703e-12 
#> outer mgc:  0.1039736 
#> iter: 1  value: 1480.129 mgc: 0.07302291 ustep: 1 
#> iter: 2  value: 1480.129 mgc: 1.017137e-05 ustep: 1 
#> iter: 3  mgc: 1.295186e-12 
#> iter: 1  mgc: 1.295186e-12 
#> outer mgc:  0.169023 
#> iter: 1  value: 1480.129 mgc: 0.03819546 ustep: 1 
#> iter: 2  value: 1480.129 mgc: 4.351893e-06 ustep: 1 
#> iter: 3  mgc: 4.287087e-13 
#> iter: 1  mgc: 4.287087e-13 
#> outer mgc:  0.2372897 
#> iter: 1  value: 1480.131 mgc: 0.05206463 ustep: 1 
#> iter: 2  value: 1480.131 mgc: 3.470354e-06 ustep: 1 
#> iter: 3  mgc: 1.936793e-13 
#> iter: 1  mgc: 1.936793e-13 
#> outer mgc:  0.07273308 
#> iter: 1  value: 1480.136 mgc: 0.0342343 ustep: 1 
#> iter: 2  value: 1480.136 mgc: 2.379581e-06 ustep: 1 
#> iter: 3  mgc: 1.672756e-13 
#> iter: 1  mgc: 1.672756e-13 
#> outer mgc:  0.2695736 
#> iter: 1  value: 1480.146 mgc: 0.04792437 ustep: 1 
#> iter: 2  value: 1480.146 mgc: 2.077332e-06 ustep: 1 
#> iter: 3  mgc: 1.829786e-13 
#> iter: 1  mgc: 1.829786e-13 
#> outer mgc:  0.157986 
#> iter: 1  value: 1480.154 mgc: 0.06479902 ustep: 1 
#> iter: 2  value: 1480.154 mgc: 1.869643e-05 ustep: 1 
#> iter: 3  mgc: 2.003761e-11 
#> iter: 1  mgc: 2.003761e-11 
#> outer mgc:  0.3175579 
#> iter: 1  value: 1480.148 mgc: 0.04825717 ustep: 1 
#> iter: 2  value: 1480.148 mgc: 1.156696e-05 ustep: 1 
#> iter: 3  mgc: 2.970094e-12 
#> iter: 1  mgc: 2.970094e-12 
#> outer mgc:  0.3742325 
#> iter: 1  value: 1480.142 mgc: 0.1491029 ustep: 1 
#> iter: 2  value: 1480.142 mgc: 9.890277e-05 ustep: 1 
#> iter: 3  mgc: 2.062159e-10 
#> iter: 1  value: 1480.145 mgc: 0.07480849 ustep: 1 
#> iter: 2  value: 1480.145 mgc: 1.94785e-05 ustep: 1 
#> iter: 3  mgc: 7.969499e-12 
#> iter: 1  mgc: 7.969499e-12 
#> outer mgc:  0.06194798 
#> iter: 1  value: 1480.139 mgc: 0.06338834 ustep: 1 
#> iter: 2  value: 1480.139 mgc: 8.74572e-06 ustep: 1 
#> iter: 3  mgc: 1.968258e-12 
#> iter: 1  mgc: 1.968258e-12 
#> outer mgc:  0.1670961 
#> iter: 1  value: 1480.141 mgc: 0.05279193 ustep: 1 
#> iter: 2  value: 1480.141 mgc: 5.827525e-06 ustep: 1 
#> iter: 3  mgc: 2.444181e-13 
#> iter: 1  mgc: 2.444181e-13 
#> outer mgc:  0.306703 
#> iter: 1  value: 1480.146 mgc: 0.05581747 ustep: 1 
#> iter: 2  value: 1480.146 mgc: 4.367623e-06 ustep: 1 
#> iter: 3  mgc: 1.686242e-13 
#> iter: 1  mgc: 1.686242e-13 
#> outer mgc:  0.09105139 
#> iter: 1  value: 1480.152 mgc: 0.01927061 ustep: 1 
#> iter: 2  value: 1480.152 mgc: 2.305752e-06 ustep: 1 
#> iter: 3  mgc: 1.430645e-13 
#> iter: 1  mgc: 1.430645e-13 
#> outer mgc:  0.1468607 
#> iter: 1  value: 1480.15 mgc: 0.04529849 ustep: 1 
#> iter: 2  value: 1480.15 mgc: 2.803687e-06 ustep: 1 
#> iter: 3  mgc: 1.623284e-13 
#> iter: 1  mgc: 1.623284e-13 
#> outer mgc:  0.1489981 
#> iter: 1  value: 1480.143 mgc: 0.01384337 ustep: 1 
#> iter: 2  value: 1480.143 mgc: 5.497168e-07 ustep: 1 
#> iter: 3  mgc: 1.40868e-13 
#> iter: 1  mgc: 1.40868e-13 
#> outer mgc:  0.1131006 
#> iter: 1  value: 1480.14 mgc: 0.03087595 ustep: 1 
#> iter: 2  value: 1480.14 mgc: 1.617853e-06 ustep: 1 
#> iter: 3  mgc: 1.667411e-13 
#> iter: 1  mgc: 1.667411e-13 
#> outer mgc:  0.03298068 
#> iter: 1  value: 1480.145 mgc: 0.03936253 ustep: 1 
#> iter: 2  value: 1480.145 mgc: 3.41323e-06 ustep: 1 
#> iter: 3  mgc: 1.562297e-13 
#> iter: 1  mgc: 1.562297e-13 
#> outer mgc:  0.2148113 
#> iter: 1  value: 1480.147 mgc: 0.0565331 ustep: 1 
#> iter: 2  value: 1480.147 mgc: 9.55984e-06 ustep: 1 
#> iter: 3  mgc: 3.076991e-13 
#> iter: 1  mgc: 3.076991e-13 
#> outer mgc:  0.1254535 
#> iter: 1  value: 1480.127 mgc: 0.06120194 ustep: 1 
#> iter: 2  value: 1480.127 mgc: 1.03877e-05 ustep: 1 
#> iter: 3  mgc: 3.718446e-13 
#> iter: 1  value: 1480.141 mgc: 0.01652442 ustep: 1 
#> iter: 2  value: 1480.141 mgc: 7.566323e-07 ustep: 1 
#> iter: 3  mgc: 1.050824e-13 
#> iter: 1  mgc: 1.050824e-13 
#> outer mgc:  0.02677116 
#> iter: 1  value: 1480.143 mgc: 0.05989632 ustep: 1 
#> iter: 2  value: 1480.143 mgc: 1.032971e-05 ustep: 1 
#> iter: 3  mgc: 3.268442e-13 
#> iter: 1  mgc: 3.268442e-13 
#> outer mgc:  0.04998623 
#> iter: 1  value: 1480.143 mgc: 0.0201617 ustep: 1 
#> iter: 2  value: 1480.143 mgc: 8.057299e-07 ustep: 1 
#> iter: 3  mgc: 1.090281e-13 
#> iter: 1  mgc: 1.090281e-13 
#> outer mgc:  0.07930762 
#> iter: 1  value: 1480.143 mgc: 0.0130915 ustep: 1 
#> iter: 2  value: 1480.143 mgc: 3.700529e-07 ustep: 1 
#> iter: 3  mgc: 1.356371e-13 
#> iter: 1  mgc: 1.356371e-13 
#> outer mgc:  0.05908456 
#> iter: 1  value: 1480.143 mgc: 0.02172328 ustep: 1 
#> iter: 2  value: 1480.143 mgc: 1.176564e-06 ustep: 1 
#> iter: 3  mgc: 1.410717e-13 
#> iter: 1  mgc: 1.410717e-13 
#> outer mgc:  0.03246201 
#> iter: 1  value: 1480.141 mgc: 0.02406119 ustep: 1 
#> iter: 2  value: 1480.141 mgc: 1.14421e-06 ustep: 1 
#> iter: 3  mgc: 1.663585e-13 
#> iter: 1  value: 1480.142 mgc: 0.0115567 ustep: 1 
#> iter: 2  value: 1480.142 mgc: 2.639823e-07 ustep: 1 
#> iter: 3  mgc: 1.486307e-13 
#> iter: 1  mgc: 1.486307e-13 
#> outer mgc:  0.02704945 
#> iter: 1  value: 1480.144 mgc: 0.017273 ustep: 1 
#> iter: 2  value: 1480.144 mgc: 2.576245e-07 ustep: 1 
#> iter: 3  mgc: 9.854425e-14 
#> iter: 1  mgc: 9.854425e-14 
#> outer mgc:  0.03445316 
#> iter: 1  value: 1480.143 mgc: 0.01000968 ustep: 1 
#> iter: 2  value: 1480.143 mgc: 4.761093e-08 ustep: 1 
#> iter: 3  mgc: 1.510148e-13 
#> iter: 1  mgc: 1.510148e-13 
#> outer mgc:  0.0185883 
#> iter: 1  value: 1480.145 mgc: 0.02400609 ustep: 1 
#> iter: 2  value: 1480.145 mgc: 9.361352e-07 ustep: 1 
#> iter: 3  mgc: 1.740964e-13 
#> iter: 1  value: 1480.144 mgc: 0.006562191 ustep: 1 
#> iter: 2  value: 1480.144 mgc: 7.744706e-08 ustep: 1 
#> mgc: 1.355103e-13 
#> iter: 1  mgc: 1.355103e-13 
#> outer mgc:  0.01544897 
#> iter: 1  value: 1480.144 mgc: 0.002845358 ustep: 1 
#> iter: 2  value: 1480.144 mgc: 1.362173e-08 ustep: 1 
#> mgc: 1.294977e-13 
#> iter: 1  value: 1480.144 mgc: 0.0044475 ustep: 1 
#> iter: 2  value: 1480.144 mgc: 3.328705e-08 ustep: 1 
#> mgc: 1.697708e-13 
#> iter: 1  mgc: 1.697708e-13 
#> outer mgc:  0.02549827 
#> iter: 1  value: 1480.142 mgc: 0.007787693 ustep: 1 
#> iter: 2  value: 1480.142 mgc: 9.677998e-08 ustep: 1 
#> mgc: 1.978161e-13 
#> iter: 1  mgc: 1.978161e-13 
#> outer mgc:  0.04557033 
#> iter: 1  value: 1480.143 mgc: 0.007262248 ustep: 1 
#> iter: 2  value: 1480.143 mgc: 1.207655e-07 ustep: 1 
#> mgc: 1.5066e-13 
#> iter: 1  mgc: 1.5066e-13 
#> outer mgc:  0.02391267 
#> iter: 1  value: 1480.145 mgc: 0.01740201 ustep: 1 
#> iter: 2  value: 1480.145 mgc: 4.869962e-07 ustep: 1 
#> iter: 3  mgc: 1.239968e-13 
#> iter: 1  value: 1480.143 mgc: 0.00437046 ustep: 1 
#> iter: 2  value: 1480.143 mgc: 3.127877e-08 ustep: 1 
#> mgc: 1.38578e-13 
#> iter: 1  mgc: 1.38578e-13 
#> outer mgc:  0.009977765 
#> iter: 1  value: 1480.143 mgc: 0.002976993 ustep: 1 
#> iter: 2  value: 1480.143 mgc: 1.108922e-08 ustep: 1 
#> mgc: 1.312787e-13 
#> iter: 1  mgc: 1.312787e-13 
#> outer mgc:  0.01202429 
#> iter: 1  value: 1480.143 mgc: 0.0008449236 ustep: 1 
#> iter: 2  mgc: 1.611447e-09 
#> iter: 1  value: 1480.144 mgc: 0.002376587 ustep: 1 
#> iter: 2  value: 1480.144 mgc: 1.275166e-08 ustep: 1 
#> mgc: 1.545394e-13 
#> iter: 1  mgc: 1.545394e-13 
#> outer mgc:  0.009668075 
#> iter: 1  value: 1480.145 mgc: 0.007253548 ustep: 1 
#> iter: 2  value: 1480.145 mgc: 7.698629e-08 ustep: 1 
#> mgc: 1.278341e-13 
#> iter: 1  value: 1480.144 mgc: 0.002710157 ustep: 1 
#> iter: 2  mgc: 8.857797e-09 
#> iter: 1  mgc: 8.857797e-09 
#> outer mgc:  0.004439967 
#> iter: 1  value: 1480.144 mgc: 0.0008817317 ustep: 1 
#> iter: 2  mgc: 2.226323e-09 
#> iter: 1  value: 1480.144 mgc: 0.001436944 ustep: 1 
#> iter: 2  mgc: 5.912493e-09 
#> iter: 1  mgc: 5.912493e-09 
#> outer mgc:  0.01051654 
#> iter: 1  value: 1480.143 mgc: 0.001265884 ustep: 1 
#> iter: 2  mgc: 2.882773e-09 
#> iter: 1  mgc: 2.882773e-09 
#> outer mgc:  0.002875987 
#> iter: 1  value: 1480.144 mgc: 0.001781123 ustep: 1 
#> iter: 2  mgc: 3.604954e-09 
#> iter: 1  mgc: 3.604954e-09 
#> outer mgc:  0.002873761 
#> iter: 1  value: 1480.144 mgc: 0.0009472973 ustep: 1 
#> iter: 2  mgc: 1.154696e-09 
#> iter: 1  mgc: 1.154696e-09 
#> outer mgc:  0.0006834342 
#> iter: 1  mgc: 1.154696e-09 
#> iter: 1  mgc: 1.154696e-09 
#> outer mgc:  0.0006834342 
#> iter: 1  mgc: 1.154696e-09 
#> outer mgc:  0.0006834342 
mohns_rho(mod) # calculate Mohn's rho
#> $Fbar
#> [1] -0.02297575
#> 
#> $SSB
#> [1] 0.02896454
#> 
#> $naa
#> , , R
#> 
#>           [,1]
#> [1,] 0.1271331
#> 
#> , , N2
#> 
#>            [,1]
#> [1,] 0.02847806
#> 
#> , , N3
#> 
#>            [,1]
#> [1,] 0.03853523
#> 
#> , , N4
#> 
#>            [,1]
#> [1,] 0.03078133
#> 
#> , , N5
#> 
#>             [,1]
#> [1,] -0.01831377
#> 
#> , , N6+
#> 
#>            [,1]
#> [1,] 0.05937164
#> 
#> 
# }
```
