wham 1.0.9.9000
=========================

### Major improvements

These features have been developed previously on the lab branch and not attributed to a specific commit, but are now included as the main version of WHAM
* 1 or more stocks and regions can now be modeled with movement among regions
* User can define seasonal intervals within years (possibly varying in length)
* Large changes in structure of input data, parameters, and reported output
* Seasonal operation of fleets
* movement parameters can incorporate time and age varying random effects and effects of environmental covariates like natural mortality
* priors for mean movement parameters are allowed because it is currently not possible to include tagging data observations.
* functions to configure various aspects of an input (e.g., set_NAA) are now exported and can be used to modify existing inputs
* greatly increased options to functions to allow user to specify initial parameter values and mapping of parameter estimation.
* random effects options for initial abundance at age
* functions to perform jittering (jitter_wham) and self tests (self_test) **(considered beta currently)**
* functions to add reference point estimation, TMB::sdreport objects and retrospective peels to a previously fitted model
* new Rmarkdown based generation of html by plot_wham_output that is self-contained (can be opened from any location).

### Minor improvements

* Can now specify fleet-specific catch or F in projection years. [f2a298e](https://github.com/timjmiller/wham/commit/f2a298e891d1608713115e9e32ff60ce10264433)
* add some more input options for M, catch, indices. [d50acc7](https://github.com/timjmiller/wham/commit/d50acc7305be99c6c6d5280313b7ad66dd85d48d)

### Bug fixes

* fix bugs in reporting of NAA re cor parameters in table of estimates and bug in setting covariate effect maps. [0ab78d8](https://github.com/timjmiller/wham/commit/0ab78d863597b3c6380c0c7cf9e199141e07d2a1)
*  fix mapping of selectivity RE when selectivity block does not span all years, correct marginal variance used for `basic_info$bias_correct_BRPs=TRUE` in the SSB/R and Y/R calculations when there is AR1 correlation of NAA RE, and correct dimension of recruitment random effects for projection of SCAA models. [b031806](https://github.com/timjmiller/wham/commit/b0318061e24b1f5820b8dfc0b369624742a2f69d)
* fix bug in use of user-defined waa, maturity in projection years. [24dd1ab](https://github.com/timjmiller/wham/commit/24dd1ab92d90aad2d9bb04dcc8e58f1a155def19)
* fix some indexing errors for retros with Ecovs or selectivity RE, fix some setup problems for movement. [2371fc3](https://github.com/timjmiller/wham/commit/2371fc3d5eb8dfe20b37ba9c9e81a0587aa863cb)
* fix specification of initial/fixed M when input M is both age and time varying. [1d637e3](https://github.com/timjmiller/wham/commit/1d637e39e900e5ac09fae3898ea1f61a0e89c1b0)
* fix specification of shared fleet selectivity blocks provided by ASAP input. [fcbdfca](https://github.com/timjmiller/wham/commit/fcbdfcacccce103ec1a6c55b000d2ad26a6d073a)

wham 1.0.9 (2024-06-13)
=========================

### IMPORTANT

* **This will be the last release for models allowing only 1 stock and 1 region.**
* The single-stock version will live on as the **single_wham** branch of the WHAM repository, but no further development will occur. Only bug fixes will be addressed.

### Bug fixes

* Fix bug in reporting of NAA autocorrelation parameters in parameter table. [ddad93c](https://github.com/timjmiller/wham/commit/ddad93c968c5986ebb2cd7aea089fbc4864fc7b8)
* Fix bug in projections for models with q random effets. [503fefc](https://github.com/timjmiller/wham/commit/503fefc86d2cb7f50ede53ce89e43c6d02c6b23a)
* Fix recruitment decoupling configuration bug. [18063fa](https://github.com/timjmiller/wham/commit/18063fad1fdd0d8787970cde9af303e0e24b9b95)


wham 1.0.8 (2024-04-24)
=========================

### Minor improvements

* Added option to bias-correct SSB/R and Y/R calculations for BRPs. Long-term projections can then be consistent with static BRPs. [4dade43](https://github.com/timjmiller/wham/commit/4dade43485ad35055c7bfc4ee982cf73dadac700)
* Added option to make recruitment in projection years consistent with that used for prevailing SPR-based BRPs [4edf29a](https://github.com/timjmiller/wham/commit/4edf29a97964b7cbacf830c25effd9100697f4ad)
* Added option to specify WAA and maturity in projection years. [d6f82f1](https://github.com/timjmiller/wham/commit/d6f82f12c525a06ef5f319c477e78a4cfed10dfb)
* Added option to decouple recruitment random effects from those for older ages when NAA_re$sigma = "rec+1". [a7880bf](https://github.com/timjmiller/wham/commit/a7880bf56dc89dd5c60c7f4832027d508d194f6c)
* Linear parameterization of Dirichlet-multinomial dispersion parameter option for age composition likelihood added [Thorson et al. 2017](https://doi.org/10.1016/j.fishres.2016.06.005) [e3ab77b](https://github.com/timjmiller/wham/commit/e3ab77b8e42c0e68e364f58a0a7855e04a262604). Effective sample size is reported for both Dirichlet-multinomial likelihood options.

### Bug fixes

* Fixed bug in make_osa_residuals when checking for a data frame. [9f65ad6](https://github.com/timjmiller/wham/commit/9f65ad62b366fc08e35df48178e87cd39d330cd3)
* Fix bug in NAA reporting simulated models with projection years when holding process errors constant in data years. [22633c](https://github.com/timjmiller/wham/commit/22633cd126a4d55e34f9d3cc39e3a4ed8229ea17) 
* Fixed bug in fit_peel function when modeling random effects on catchability. [a7880bf](https://github.com/timjmiller/wham/commit/a7880bf56dc89dd5c60c7f4832027d508d194f6c)
* Fixed issue in multinomial, D-M, MVTweedie likelihoods when selectivity = 0 for some ages with age-specific selectivity configuration [issue 80](https://github.com/timjmiller/wham/issues/80) [ad22b8a](https://github.com/timjmiller/wham/commit/ad22b8a2dec86d22525f58f5f7db9bc6e3d18fc6)
* Fixed NaN issue in multinomial for extremely small (e.g. 1e-17) predicted proportions [issue 76](https://github.com/timjmiller/wham/issues/76) [99ec44f](https://github.com/timjmiller/wham/commit/99ec44f7f3c4647d1e882d3102d3da5d31afbb8b)
* Fix incorrect flag for simulating catch age composition data and incorrect indexing for years to average over for SPR inputs [4f06dc7](https://github.com/timjmiller/wham/commit/4f06dc7baff3b6bef3e7c8eb262dcdedefbb33dc)
* Fixed bug in reporting correlation parameter estimates in pdf/html tables [31ea94b](https://github.com/timjmiller/wham/commit/31ea94b611844a557947b11ac758531abb83005b)
* Fixed bug in setting up catch age comp observations when ages are omitted due to selectivity=0 [6667efa](https://github.com/timjmiller/wham/commit/6667efac36df17215b5733e13866a3a78ed4f3b8)
* Fixed reference point estimation and use in projections when there are multiple fleets [65c0130](https://github.com/timjmiller/wham/commit/65c0130f56e9346a98c1a1a957efe6dd1441c4c1) [b4c1ca3](https://github.com/timjmiller/wham/commit/b4c1ca379476042ef3ae376179df369c13e1fc0e). 

wham 1.0.7 (2022-11-3)
=========================

### Major improvements

* Improved stability and reliability of one-step-ahead residuals for age compositon observations for most log-likelihoods [037b714](https://github.com/timjmiller/wham/commit/037b7145927824359cdba00ecdffc1e4fbaceee6).
* Multivariate Tweedie age composition likelihood option added (Thorson et al. in press) [6f77c16](https://github.com/timjmiller/wham/commit/6f77c164289b029725ecc0882ec50865832b696f).

### Minor improvements

* Added exported make_osa_residuals function that can make OSA residuals from an object returned by fit_wham where it is also used internally with do.osa=TRUE [1616ade](https://github.com/timjmiller/wham/commit/1616aded076a63c5ac375abf706ea710de6e9d0e). 
* Added ability to specify different F options for each year of projections to project_wham and prepare_projection [f763059](https://github.com/timjmiller/wham/commit/f763059a61fde0814f2de83ef2d08045a3ed59e0)
* Revised vignette on simulation studies [f763059](https://github.com/timjmiller/wham/commit/f763059a61fde0814f2de83ef2d08045a3ed59e0)

### Bug fixes

* Fix some plotting errors when there are multiple fleets and when there is no age comp for some fleets or indices [83f23ff](https://github.com/timjmiller/wham/commit/83f23ff2c2c676577be709135505cb5348a1c632)
* Fix age comp observations when selectivity is assumed 0 for one or more age classes (needed for osa obs) [7bba974](https://github.com/timjmiller/wham/commit/7bba974c4d9ea8d772902be4d62512f20274a3ad)
* Fix bug in simulation of selectivity random effects with "ar1_y" or "ar1" option [190000c](https://github.com/timjmiller/wham/commit/190000ccfa3ed2fb7ba30ee8d4af41329fdca655)
* Fix bug in simulation of M random effects with "ar1_y" option [77bbd94](https://github.com/timjmiller/wham/commit/77bbd946e4881216a439933473d1c58b21c270c3)
* Fix bug in check_projF which tests whether the F in projections is being specified correctly when FXSPR or F at catch is specified [c643d4b](https://github.com/timjmiller/wham/commit/c643d4ba8339c13dd4b9e3662aaa29f26d309624)
* Fix index proportions at age data specification in basic_info [issue 64](https://github.com/timjmiller/wham/issues/64) [a3e3afc](https://github.com/timjmiller/wham/commit/a3e3afc9b23e2ca4e3d369581dfcf2b33732686c)
* Fix age of full selection for N1_model = 1 [issue 56](https://github.com/timjmiller/wham/issues/56) [865ca3d](https://github.com/timjmiller/wham/commit/865ca3d304c2d5eddcf044fe3d3776ada42b07aa).

wham 1.0.6 (2022-04-08)
=========================

### Major improvements

* One-step-ahead prediction quantile residuals for age compositon observations for most log-likelihoods (Trijoulet et al. 2023) [80e3bba](https://github.com/timjmiller/wham/commit/80e3bbae0244bd1199f018dca53f09b08f5fd203).
* Can specify a prior distribution on fully-selected catchability (logit-scale) which is then estimated as a random effect [fcbc604](https://github.com/timjmiller/wham/commit/fcbc604068005cd7cce6c3f329376c2b4ef7b540).
* Auto-regressive random effects for fully-selected catchability [fcbc604](https://github.com/timjmiller/wham/commit/fcbc604068005cd7cce6c3f329376c2b4ef7b540).
* Environmental covariate effects on fully-selected catchability [fcbc604](https://github.com/timjmiller/wham/commit/fcbc604068005cd7cce6c3f329376c2b4ef7b540).
* One or more environmental covariates can have multiple effects on recruitment, M, and catchability [fcbc604](https://github.com/timjmiller/wham/commit/fcbc604068005cd7cce6c3f329376c2b4ef7b540).

### Minor improvements

* Expanded set of age composition models including allowing AR1 correlation for logistic-normal [80e3bba](https://github.com/timjmiller/wham/commit/80e3bbae0244bd1199f018dca53f09b08f5fd203).
* Example 11 script and vignette (and tests) added to demonstrated extensions for modeling catchability [173555d](https://github.com/timjmiller/wham/commit/173555d1485d52ca58aa77ae77fad856fa8b3bbb).
* New plots added to plot_wham_output.R and wham_plots_tables.R to include new catchability features [173555d](https://github.com/timjmiller/wham/commit/173555d1485d52ca58aa77ae77fad856fa8b3bbb).
* Added reporting of SPR-based reference points that use SPR inputs averged over the same years as they would for FXSPR in projection years [5d1cd9d](https://github.com/timjmiller/wham/commit/5d1cd9dd6b4ec8214e988be1fb4c5e6f4b51a661).
* Fixed and generalized configuration for estimation of observation error variances for environmental covariates for `prepare_wham_input` (`ecov$logsigma`) [a3328a9](https://github.com/timjmiller/wham/commit/a3328a9a9652a7fc3f0eb45ee4688fcdaa4334f7).
* Added tables (html or pdf) of estimated parameters, abundance at age, and fishing mortality at age to results generated by `plot_wham_output` [13094c6](https://github.com/timjmiller/wham/commit/13094c674a53d495021711531c083184e3c6b892).

### Bug fixes

* Fix age labeling for selectivity and numbers at age retro plots.
* Fix typo in set_NAA.r for N1_model option (issue [issue 55](https://github.com/timjmiller/wham/issues/55) caught by @Cole-Monnahan-NOAA).
* Fix errors in `prepare_wham_input` introduced in 1.0.5 when not all index data in `.dat` are used (aggregate or proportion-at-age, all years or a subset) [12abdef](https://github.com/timjmiller/wham/commit/12abdefb22e7b8b5f0640b0389b5e07fc7a00877)
* Updated version of `TMBhelper::check_estimability` is now internal [issue 47](https://github.com/timjmiller/wham/issues/47)
* Peel `which_F_age` (fix error when projecting a peel doing retrospective predictions, ex9)
* Ecov AR1 mean parameterization
* simulation of Dirichlet-multinomial age composition ([issue 49](https://github.com/timjmiller/wham/issues/49), caught by @seananderson).
* error reported in [issue 48](https://github.com/timjmiller/wham/issues/48). There is no `opt` element when do.fit = FALSE.
* `plot_wham_output` now removes any existing plot files in the directories where plots are saved before writing files so that any orphaned files from previous fits are not included in the results and html rendering.
* `plot_wham_output` now creates `.png` output by default because of 1) issues opening the `.html` in Firefox on Windows, and 2) it was annoying to close dozens of browser tabs during testing. issues [42](https://github.com/timjmiller/wham/issues/42) and [50](https://github.com/timjmiller/wham/issues/50), thanks to @ejosymart and @JDeroba.
* `set_catch` bug when catch age comp data exist for some years and not others

wham 1.0.5 (2021-08-25)
=========================

### Major improvements

* `prepare_wham_input` is now modularized and can take a "basic_info" argument or generate a dummy input for fit_wham without an asap3 object. This provides better organization and helps fix the operating model/MSE features that have been broken since version 1.0.0. Options have been added to NAA_re argument to configure how initial numbers at age and recruitment are treated. Including initial parameter values. A tenth example and vignette illustrate the usage to generate an assessment model without an asap3 file and a simple management strategy evaluation [9bc161d](https://github.com/timjmiller/wham/commit/87d87c9ded5c3fa58db779cbd1ac8eb911950e4e). 
  
### Minor improvements

* Options have been added to NAA_re argument to configure how initial numbers at age and recruitment are treated. Including initial parameter values.
* NOAA pkgdown website template, thanks to @Bai-Li-NOAA [7618bb](https://github.com/timjmiller/wham/commit/7618bb0427743fffbb6db49465fd8c2f73993719)
* Add `ecov$ages`, allows user to specify that an ecov only affects a subset of ages [e6518f](https://github.com/timjmiller/wham/commit/e6518f8f09f1c252517e6dc89b98776d687d417b) 

### Bug fixes

* Don't plot age comp residuals for fleets/indices/years that are not used [6876ee](https://github.com/timjmiller/wham/commit/6876eed1dd4bbdc7b26295356756bb7a30053887)
* Fix `Ecov_re` padding in projections with multiple ecovs that have different lags [a57f3b](https://github.com/timjmiller/wham/commit/a57f3b339881a39dbee3ab0ab190c732f855dec1)
* Default to use static reference points in projections (so ref pts are consistent across multiple projections) [74a6c3](https://github.com/timjmiller/wham/commit/74a6c3cda74e5f2f14a6471e410b227c407e0790)
* Diagnostic plots for log catch and indices did not use bias correction [issue 46](https://github.com/timjmiller/wham/issues/46). Now they do.
* Fixed bug in index diagnostic plot when index is in biomass and age comp is in numbers.

wham 1.0.4 (2021-05-03)
=========================

WHAM description + simulation test paper published: [https://doi.org/10.1016/j.fishres.2021.105967](https://doi.org/10.1016/j.fishres.2021.105967)

### Major improvements

* New ability to plot results from multiple ASAP3 and WHAM models together for comparison. Thanks to @liz-brooks for contributing `read_asap3_fit`! See `?compare_wham_models` and new [vignette 8](https://timjmiller.github.io/wham/articles/ex8_compare.html).

### Minor improvements

* New `save.sdrep = F` option to only save `summary(sdreport)` instead of `sdreport`. Can make saved models MUCH smaller (e.g. 2 MB vs. 150 MB). ([2f8875](https://github.com/timjmiller/wham/commit/2f8875323c0d6845a92444a9e7d4aaa92fe29d8d)).
* Added `proj.opts$percentFXSPR` option, percent of F_XSPR to use for calculating catch in projections. For example, GOM cod uses F = 75% F_XSPR, so `proj.opts$percentFXSPR = 75`.
* Added `proj.opts$useFMSY` and `proj.opts$percentFMSY` options, to project population and catches at (a percentage of) F_MSY. ([467532](https://github.com/timjmiller/wham/commit/46753202499cbf256fb79796f8570e045dddb9db))
* Allow peels to be projected so can do retrospective predictions, as in Fig. 4 of [Xu et al. 2018](https://onlinelibrary.wiley.com/doi/abs/10.1111/fog.12236). See new [vignette 9](https://timjmiller.github.io/wham/articles/ex9_retro_pred.html). ([e1e00c](https://github.com/timjmiller/wham/commit/e1e00c712fcba0d0631dd9321e3bb83b87a1485f#diff-7f3a4433e0573f6fd6787d47614e0303f714e2d71d9e618ecfb832195987b732))

### Bug fixes

* Broken links in vignettes (thanks to @tcarruth)
* Plotting issue with Ecov OSA residuals and peels (thanks to @h-du-pontavice, [1163df](https://github.com/timjmiller/wham/commit/1163df4290387174a5336ada98a13dc0f5a9644c))
* Double logistic selectivity setup in `prepare_wham_input` (thanks to @tcarruth, [f270dd](https://github.com/timjmiller/wham/commit/f270ddb66d253ac2aaf9a0109631b89036ddcd5e))
* `tryCatch` error assignment issues ([9d5c87](https://github.com/timjmiller/wham/commit/9d5c8792769dae9e640da3ad44b7f3e6b74e4a87))
* Remove error in default projection options (`percentSPR`)
* Setting up projections with multiple Ecovs (again thanks to @h-du-pontavice, [42a6a4](https://github.com/timjmiller/wham/commit/42a6a4950e85613219525e89c5590c37f3a6369f))
* Selectivity parameter initial values set to middle of range if unspecified in `prepare_wham_input`
* Fixed `fit_tmb` to make `$final_gradient` reported by wham equivalent to `sdreport()$gradient.fixed`. Issue was that the `model$env$last.par.best` is not updated by the newton steps after optimization so `model$opt$par` was slightly different (10^(-7) or smaller). [e20bd8](https://github.com/timjmiller/wham/commit/e20bd8d01a32b1cbdb826905257555f9a8b55c75)
* Factor order problem caused indices to be out of order when > 9 indices present. [9e0141](https://github.com/timjmiller/wham/commit/9e0141b59bb80d17b9fc79662e3f38b24b8991f1)
* Fixed default determination of age corresponding to fully selected (total) fishing mortality. Now `which_F_age` is a vector of annual values (including projection years) and the values are reassigned after fitting before reporting. ([189b6c](https://github.com/timjmiller/wham/commit/189b6cede494a53763afcda5b0bf6854fed095d5))
* Added check for F_XSPR to verify percent of SPR0 is correct and, in projections, catch at F from catch is verified. Initial values for these algorithms can be specified by the user. ([189b6c](https://github.com/timjmiller/wham/commit/189b6cede494a53763afcda5b0bf6854fed095d5))
* If ecov years differ from model years, peel to same year in `fit_peel` ([e1e00c](https://github.com/timjmiller/wham/commit/e1e00c712fcba0d0631dd9321e3bb83b87a1485f))

wham 1.0.3 (2021-02-05)
=========================

### Minor improvements

* allow for multiple ecov effects on recruitment and M
* add WHAM version to model object (version and GitHub commit). Requires `sessioninfo` package (within `devtools`).
* for SCAA (recruitment as fixed effects), treat recruitment in projections as random effects with mean and SD calculated from model years (previously was fixed at `exp(10)`). This approximates ASAP + AgePro, allows for simulating projected recruitment according to the ECDF.

### Bug fixes

* clean up error messages (remove global `err` object) [53f5ec6](https://github.com/timjmiller/wham/commit/53f5ec67a6f60ea85b928debc4a57d0ee6673e78)
* adjust age composition variance (`tau`) by effective sample size (`Neff`) for option 7 [4c13331](https://github.com/timjmiller/wham/commit/4c133312b2ecaab3972be13a96b4b456d9c1f6b0)

wham 1.0.2 (2020-12-14)
=========================

### Minor improvements

* add `1e-15` to predicted proportions to make age composition likelihoods robust to 0 predictions when selAA is fixed at 0. This affects the multinomial, Dirichlet, and Dirichlet-multinomial (options 1-3), since the logistic normal (options 4-7) already did this. [88f15d4](https://github.com/timjmiller/wham/commit/88f15d4a51f69a3d649d76bcac0a8cf299c3135e)
* specify age composition model using `age_comp` argument to `prepare_wham_input`. See [`?prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.html) for details. [fd94b3d](https://github.com/timjmiller/wham/commit/fd94b3dcaf189482e10a6750c2f1b8350837fd48)
* estimate age-specific M for any age(s) ([60f1358](https://github.com/timjmiller/wham/commit/60f13584e515950bf43f358c5a0cfdb0e8a30241))

### Bug fixes

* check for sel par inits outside lower/upper bounds ([4ed394b](https://github.com/timjmiller/wham/commit/4ed394bac38d6054727bc4d3e17c8f4452ae8289))
* fleet weight-at-age pointers fixed ([059e66a](https://github.com/timjmiller/wham/commit/059e66a1cc6e90554862b3e30e815cf0d5cd18ab))
* couple small fixes to handle multiple Ecovs ([9e46b48](https://github.com/timjmiller/wham/commit/9e46b4818e01b2bfd12e500bd4376869001f265e))

wham 1.0.1 (2020-11-12)
=========================

### Minor improvements

* `run_all_examples.R` script to test all examples work
* new vignette on debugging WHAM models

### Bug fixes

* update `tidyr::gather` to `tidyr::pivot_longer` ([80bf89a](https://github.com/timjmiller/wham/commit/80bf89a69975a4ac5c037127d094afc3b00b7f59))
* if all selpars fixed then use ASAP file inits ([2887e87](https://github.com/timjmiller/wham/commit/2887e87aa4af08530765fbf04c002a00afd8efa1))
* update example scripts 4, 5, 6 ([31cb9a9](https://github.com/timjmiller/wham/commit/31cb9a931a83e9c76e5d63fd52243de917fe7c99))
* avoid font family warnings on Windows ([af42be4](https://github.com/timjmiller/wham/commit/af42be4d940785539f00e3ec29a85c605ef2a0de))

wham 1.0.0 (2020-10-07)
=========================

* This is the first release of WHAM.

