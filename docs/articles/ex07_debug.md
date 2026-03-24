# Ex 7: Debugging WHAM models

In this vignette we show how to use the `do.fit = FALSE` flag to debug a
WHAM model.

``` r

# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
```

Load the ASAP3 data file.

``` r

wham.dir <- find.package("wham")
asap3 <- read_asap3_dat(file.path(wham.dir,"extdata","ex7_SNEMAYT.dat"))
```

Prepare the WHAM model (`m3` from [example
1](https://timjmiller.github.io/wham/articles/ex1_SNEMA_yellowtail_flounder.html)).

``` r

input <- prepare_wham_input(asap3, recruit_model=2, model_name="Ex 1: SNEMA Yellowtail Flounder",
                              NAA_re = list(sigma="rec+1", cor="iid"))
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
#> --NAA--------------------------------------------------------------------------------------------------------------------------------
#> 
#>  Same NAA_re$sigma being used for all stocks (rec+1).
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
```

Try to fit the model… uh oh.

``` r

mod <- fit_wham(input, do.osa = F, do.retro = F)
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
#> iter: 1  Error in iterate(par) : 
#>   Newton dropout because inner gradient had non-finite components.
#> iter: 1  Error in iterate(par) : 
#>   Newton dropout because inner gradient had non-finite components.
#> Using stats::nlminb for optimization with these control parameters:
#> iter.max: 1000
#> eval.max: 1000
#> iter: 1  Error in iterate(par) : 
#>   Newton dropout because inner gradient had non-finite components.
#> Warning in stats::nlminb(model$par, model$fn, model$gr, control = opt.control):
#> NA/NaN function evaluation
#> iter: 1  Error in iterate(par) : 
#>   Newton dropout because inner gradient had non-finite components.
#> Error in ff(x, order = 1) : 
#>   inner newton optimization failed during gradient calculation
#> outer mgc:  NaN
#> Error in `mod$opt$par`:
#> ! $ operator is invalid for atomic vectors
```

What’s wrong? It looks like the likelihood function is returning `NaN`.
It is often easier to diagnose problems like this using the
*unoptimized* model, i.e. look at everything using the initial parameter
values. WHAM includes a `do.fit = F` flag in `fit_wham` to return the
unoptimized model object returned by
[`TMB::MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html). Let’s
see how it works.

``` r

mod <- fit_wham(input, do.fit = F)
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
#> iter: 1  Error in iterate(par) : 
#>   Newton dropout because inner gradient had non-finite components.
#> iter: 1  Error in iterate(par) : 
#>   Newton dropout because inner gradient had non-finite components.
#> iter: 1  Error in iterate(par) : 
#>   Newton dropout because inner gradient had non-finite components.
#> iter: 1  Error in iterate(par) : 
#>   Newton dropout because inner gradient had non-finite components.
#> iter: 1  Error in iterate(par) : 
#>   Newton dropout because inner gradient had non-finite components.
#> iter: 1  Error in iterate(par) : 
#>   Newton dropout because inner gradient had non-finite components.
#> iter: 1  Error in iterate(par) : 
#>   Newton dropout because inner gradient had non-finite components.
```

This runs without a fatal error. The optimization failed because the
likelihood was `NaN`, and now we can see which of the likelihood
components was responsible. To do so, we need to look at the `REPORT`ed
objects with `"nll"` in their name. We can use the `$report()` function
from TMB.

``` r

therep = mod$report()
```

See all of the `REPORT`ed objects.

``` r

names(therep)
#>  [1] "L"                         "selAA"                    
#>  [3] "FAA_static"                "waa_catch_static"         
#>  [5] "NAA_index"                 "index_Neff_out"           
#>  [7] "waa_catch"                 "pred_stock_CAA"           
#>  [9] "trans_mu_base"             "N1"                       
#> [11] "NAAPR0_static"             "log_SPR_FXSPR_static"     
#> [13] "log_index_resid"           "log_FXSPR_iter"           
#> [15] "log_F_tot"                 "log_FAA_tot"              
#> [17] "nll_M"                     "Ecov_obs_sigma"           
#> [19] "log_YPR_FXSPR"             "mu_static"                
#> [21] "mature_static"             "log_FAA_XSPR_static"      
#> [23] "YPR_srf_FXSPR_static"      "log_M_static"             
#> [25] "mature_all"                "Ecov_x"                   
#> [27] "fracyr_SSB_all"            "pred_log_catch"           
#> [29] "pred_IAA"                  "NAA"                      
#> [31] "log_FAA_by_region"         "log_SPR0"                 
#> [33] "FAA"                       "pred_NAA"                 
#> [35] "seasonal_Ps_terminal_year" "annual_Ps"                
#> [37] "nll_sel"                   "nll_agg_catch"            
#> [39] "QAA"                       "log_M"                    
#> [41] "log_Fbar_XSPR"             "NAAPR_FXSPR_static"       
#> [43] "log_Y_FXSPR_static"        "log_FXSPR_iter_static"    
#> [45] "FAA_by_region"             "selpars"                  
#> [47] "pred_log_indices"          "nll_Ecov_obs"             
#> [49] "log_YPR_FXSPR_static"      "waa_ssb_static"           
#> [51] "annual_SPR0AA"             "q"                        
#> [53] "all_NAA"                   "mu"                       
#> [55] "pred_N1"                   "log_FXSPR_static"         
#> [57] "log_SSB_FXSPR"             "log_FAA_XSPR"             
#> [59] "waa_ssb"                   "R_XSPR"                   
#> [61] "sel_static"                "pred_CAA"                 
#> [63] "log_FXSPR"                 "log_Y_FXSPR"              
#> [65] "NAA_spawn"                 "NAA_devs"                 
#> [67] "pred_catch_paa"            "pred_index_paa"           
#> [69] "Fbar"                      "nll_agg_indices"          
#> [71] "log_catch_resid"           "nll"                      
#> [73] "annual_SAA_spawn"          "nll_Ecov_obs_sig"         
#> [75] "nll_NAA"                   "catch_Neff_out"           
#> [77] "logit_q_mat"               "pred_catch"               
#> [79] "log_SSB_FXSPR_static"      "log_Fbar_XSPR_static"     
#> [81] "pred_stock_catch"          "marg_NAA_sigma"           
#> [83] "nll_catch_acomp"           "MAA"                      
#> [85] "SSB"                       "selpars_re_mats"          
#> [87] "nll_index_acomp"           "log_SPR_FXSPR"            
#> [89] "pred_indices"              "log_SPR0_static"
```

Now just get the objects with `"nll"` in their name, and sum over all
individual values.

``` r

sapply(grep("nll",names(therep),value=T), function(x) sum(therep[[x]]))
#>            nll_M          nll_sel    nll_agg_catch     nll_Ecov_obs 
#>           0.0000           0.0000              NaN           0.0000 
#>  nll_agg_indices              nll nll_Ecov_obs_sig          nll_NAA 
#>        1565.7938              NaN           0.0000         325.8628 
#>  nll_catch_acomp  nll_index_acomp 
#>        2844.2380        3921.8316
```

The likelihood components that are equal to 0 are not used in the model
(no random effects on `M`, `Ecov`, `selectivity`, etc.). `nll` is the
total likelihood and is `NaN`. The one troublesome component is
`nll_agg_catch`. This is the total (aggregate) catch from the fishery in
each year. It could be an issue with the catch data, the model-predicted
catch, or the input CVs for the annual catches, since they all are used
in the likelihood calculation. Let’s take a closer look at all of the
annual values in `nll_agg_catch`

``` r

therep$nll_agg_catch
#>              [,1]
#>  [1,] 103.3840617
#>  [2,] 122.5550305
#>  [3,]   8.4553933
#>  [4,]   0.3642440
#>  [5,]   1.7656695
#>  [6,]  14.0697924
#>  [7,]  32.3760427
#>  [8,]   9.0788070
#>  [9,]   8.8434262
#> [10,]  69.4381045
#> [11,] 119.1118361
#> [12,]  57.8351840
#> [13,]   0.7705734
#> [14,]  -0.9505520
#> [15,]   7.4019297
#> [16,]  13.9873964
#> [17,]  13.7551021
#> [18,]  96.2002607
#> [19,]   7.2293432
#> [20,]  15.4822178
#> [21,] 179.5247068
#> [22,] 128.8416770
#> [23,] 307.9921432
#> [24,] 156.5912411
#> [25,]  83.0194806
#> [26,]  97.6718012
#> [27,] 117.3983487
#> [28,]  79.6474016
#> [29,]  78.8584557
#> [30,] 143.3232225
#> [31,] 215.3673272
#> [32,] 208.3245273
#> [33,] 343.2713739
#> [34,] 294.8287577
#> [35,] 246.4181361
#> [36,] 232.2362516
#> [37,] 263.3032268
#> [38,] 367.0891293
#> [39,] 296.4944154
#> [40,] 203.8510507
#> [41,] 174.3395764
#> [42,]         NaN
#> [43,] 318.8528146
#> [44,] 524.1093361
```

We see that the problem is in year 42. The user will not necessarily
know what names are used in WHAM to denote the different inputs to
`nll_agg_catch`, so we can search the [WHAM .cpp
file](https://github.com/timjmiller/wham/blob/d3370a0f82d0f2012a1a19afdfb4b2c29da73720/src/multi_wham.cpp#L915)
for “nll_agg_catch”. We find that it depends on `agg_catch` (the catch
data), `pred_log_catch` (model-predicted log catch), `agg_catch_sigma`,
and `log_catch_sig_scale`. The aggregate catch CVs from ASAP are
transformed to standard deviations for the log-normal assumption on
aggregate catch and supplied to WHAM as `input$data$agg_catch_sigma`. It
is possible to attempt to estimate a scalar multiple of the annual SDs
via `log_catch_sig_scale`, but it is not used by default.

The catch data looks ok.

``` r

input$data$agg_catch
#>             [,1]
#>  [1,] 14549.0000
#>  [2,] 17088.0000
#>  [3,]  5732.0000
#>  [4,]  3436.0000
#>  [5,]  5223.0000
#>  [6,]  8085.0000
#>  [7,]  9883.0000
#>  [8,]  8021.0000
#>  [9,]  6607.0000
#> [10,] 15764.0000
#> [11,] 22211.0000
#> [12,] 11225.0000
#> [13,]  4817.0000
#> [14,]  4620.0000
#> [15,]  2652.0000
#> [16,]  2782.0000
#> [17,]  8349.0000
#> [18,] 17916.0000
#> [19,]  6430.0000
#> [20,]  2695.0000
#> [21,]   771.0008
#> [22,]   735.0004
#> [23,]   343.0004
#> [24,]   759.0001
#> [25,]  1222.0000
#> [26,]  1087.0000
#> [27,]  1403.0000
#> [28,]  1397.0000
#> [29,]  1449.0000
#> [30,]   945.0005
#> [31,]   666.0005
#> [32,]   619.0002
#> [33,]   346.0002
#> [34,]   396.0000
#> [35,]   502.0004
#> [36,]   583.0006
#> [37,]   453.0006
#> [38,]   290.9995
#> [39,]   390.0002
#> [40,]   563.0002
#> [41,]   645.9999
#> [42,]   625.0000
#> [43,]   337.0005
#> [44,]   151.9999
```

The catch SDs looks ok.

``` r

input$data$agg_catch_sigma
#>             [,1]
#>  [1,] 0.09975135
#>  [2,] 0.09975135
#>  [3,] 0.09975135
#>  [4,] 0.09975135
#>  [5,] 0.09975135
#>  [6,] 0.09975135
#>  [7,] 0.09975135
#>  [8,] 0.09975135
#>  [9,] 0.09975135
#> [10,] 0.09975135
#> [11,] 0.09975135
#> [12,] 0.09975135
#> [13,] 0.09975135
#> [14,] 0.09975135
#> [15,] 0.09975135
#> [16,] 0.09975135
#> [17,] 0.09975135
#> [18,] 0.09975135
#> [19,] 0.09975135
#> [20,] 0.09975135
#> [21,] 0.09975135
#> [22,] 0.09975135
#> [23,] 0.09975135
#> [24,] 0.09975135
#> [25,] 0.09975135
#> [26,] 0.09975135
#> [27,] 0.09975135
#> [28,] 0.09975135
#> [29,] 0.09975135
#> [30,] 0.09975135
#> [31,] 0.09975135
#> [32,] 0.09975135
#> [33,] 0.09975135
#> [34,] 0.09975135
#> [35,] 0.09975135
#> [36,] 0.09975135
#> [37,] 0.09975135
#> [38,] 0.09975135
#> [39,] 0.09975135
#> [40,] 0.09975135
#> [41,] 0.09975135
#> [42,] 0.09975135
#> [43,] 0.09975135
#> [44,] 0.09975135
```

The predicted log catch in year 42 is the issue.

``` r

therep$pred_log_catch
#>           [,1]
#>  [1,] 8.141326
#>  [2,] 8.175619
#>  [3,] 8.211267
#>  [4,] 8.328701
#>  [5,] 8.310382
#>  [6,] 8.443164
#>  [7,] 8.378883
#>  [8,] 8.533464
#>  [9,] 8.344692
#> [10,] 8.478282
#> [11,] 8.459799
#> [12,] 8.240291
#> [13,] 8.272735
#> [14,] 8.345046
#> [15,] 8.301266
#> [16,] 8.484047
#> [17,] 8.480970
#> [18,] 8.399881
#> [19,] 8.354660
#> [20,] 8.478542
#> [21,] 8.545121
#> [22,] 8.209722
#> [23,] 8.319028
#> [24,] 8.405094
#> [25,] 8.404287
#> [26,] 8.395213
#> [27,] 8.783863
#> [28,] 8.511973
#> [29,] 8.542322
#> [30,] 8.548186
#> [31,] 8.578196
#> [32,] 8.470991
#> [33,] 8.465392
#> [34,] 8.409350
#> [35,] 8.439292
#> [36,] 8.524398
#> [37,] 8.410996
#> [38,] 8.381255
#> [39,] 8.400900
#> [40,] 8.354259
#> [41,] 8.340842
#> [42,]      NaN
#> [43,] 8.344558
#> [44,] 8.257718
```

The predicted log catch is just a
[log-transformation](https://github.com/timjmiller/wham/blob/d3370a0f82d0f2012a1a19afdfb4b2c29da73720/src/multi_wham.cpp#L912)
of `pred_catch` which [aggregates over stock-specific
catches](https://github.com/timjmiller/wham/blob/d3370a0f82d0f2012a1a19afdfb4b2c29da73720/src/multi_wham.cpp#L910)
(if more than 1). The stock-specific catches are a
[function](https://github.com/timjmiller/wham/blob/d3370a0f82d0f2012a1a19afdfb4b2c29da73720/src/multi_wham.cpp#L908)
of stock-specific catch-at-age and weight-at-age for each fleet. These
objects are reported out by wham, so first let’s look at
`pred_stock_CAA` which is an array (n_fleetx x n_stocks x n_years x
n_ages). There is only 1 fleet and 1 stock so:

``` r

therep$pred_stock_CAA[1,1,,]
#>            [,1]     [,2]     [,3]     [,4]     [,5]     [,6]
#>  [1,] 5435.8017  725.069 1948.968 1417.205  958.388 1175.115
#>  [2,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#>  [3,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#>  [4,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#>  [5,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#>  [6,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#>  [7,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#>  [8,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#>  [9,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [10,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [11,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [12,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [13,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [14,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [15,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [16,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [17,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [18,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [19,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [20,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [21,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [22,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [23,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [24,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [25,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [26,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [27,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [28,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [29,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [30,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [31,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [32,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [33,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [34,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [35,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [36,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [37,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [38,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [39,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [40,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [41,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [42,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [43,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
#> [44,]  888.1566 1363.739 1654.605 1782.931 1839.482 1875.122
```

The catch in numbers at age looks ok. So, lets look at `waa_catch` which
is an array (n_fleet x n_years x n_ages):

``` r

therep$waa_catch[1,,] 
#>          [,1]  [,2]  [,3]  [,4]  [,5]  [,6]
#>  [1,]   0.210 0.296 0.348 0.374 0.382 0.428
#>  [2,]   0.203 0.308 0.352 0.396 0.439 0.457
#>  [3,]   0.218 0.289 0.376 0.432 0.435 0.481
#>  [4,]   0.228 0.303 0.408 0.498 0.499 0.557
#>  [5,]   0.215 0.283 0.381 0.504 0.513 0.542
#>  [6,]   0.234 0.292 0.383 0.536 0.662 0.656
#>  [7,]   0.189 0.301 0.364 0.475 0.590 0.662
#>  [8,]   0.206 0.281 0.384 0.500 0.682 0.925
#>  [9,]   0.140 0.262 0.342 0.474 0.596 0.650
#> [10,]   0.226 0.263 0.353 0.499 0.660 0.833
#> [11,]   0.175 0.261 0.339 0.496 0.668 0.819
#> [12,]   0.182 0.237 0.295 0.388 0.487 0.656
#> [13,]   0.183 0.260 0.365 0.408 0.504 0.608
#> [14,]   0.186 0.284 0.331 0.463 0.587 0.642
#> [15,]   0.247 0.268 0.353 0.404 0.520 0.631
#> [16,]   0.270 0.293 0.396 0.493 0.611 0.821
#> [17,]   0.061 0.216 0.275 0.489 0.735 0.957
#> [18,]   0.204 0.255 0.290 0.366 0.613 0.884
#> [19,]   0.090 0.214 0.294 0.378 0.664 0.798
#> [20,]   0.110 0.303 0.375 0.447 0.631 0.918
#> [21,]   0.122 0.314 0.420 0.439 0.640 1.040
#> [22,]   0.078 0.247 0.321 0.387 0.480 0.622
#> [23,]   0.076 0.216 0.325 0.401 0.579 0.758
#> [24,]   0.102 0.335 0.368 0.457 0.604 0.740
#> [25,]   0.139 0.251 0.396 0.466 0.584 0.768
#> [26,]   0.160 0.287 0.367 0.494 0.567 0.726
#> [27,]   0.131 0.309 0.400 0.557 0.629 1.695
#> [28,]   0.185 0.321 0.444 0.561 0.667 0.752
#> [29,]   0.145 0.360 0.419 0.567 0.684 0.824
#> [30,]   0.164 0.330 0.438 0.574 0.764 0.751
#> [31,]   0.095 0.313 0.413 0.572 0.722 0.945
#> [32,]   0.136 0.295 0.436 0.540 0.581 0.799
#> [33,]   0.102 0.295 0.415 0.511 0.634 0.795
#> [34,]   0.110 0.251 0.373 0.475 0.607 0.783
#> [35,]   0.111 0.268 0.363 0.472 0.628 0.834
#> [36,]   0.151 0.266 0.388 0.461 0.574 1.077
#> [37,]   0.105 0.281 0.367 0.502 0.601 0.753
#> [38,]   0.099 0.281 0.414 0.470 0.573 0.702
#> [39,]   0.130 0.280 0.412 0.491 0.572 0.717
#> [40,]   0.120 0.251 0.365 0.475 0.554 0.709
#> [41,]   0.154 0.254 0.351 0.453 0.565 0.683
#> [42,] -99.000 0.345 0.396 0.459 0.565 0.689
#> [43,]   0.046 0.281 0.388 0.476 0.548 0.685
#> [44,]   0.000 0.270 0.349 0.434 0.521 0.629
```

Ah, here is the problem. The weight at age data has an entry of `-99` in
year 42. This means that `pred_catch` is negative, and we take the log
of a negative number for `pred_log_catch`.

If we fix this issue with the data file, the model runs fine!
