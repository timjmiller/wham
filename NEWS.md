wham 1.0.3.9000
=========================

WHAM description + simulation test paper published: [https://doi.org/10.1016/j.fishres.2021.105967](https://doi.org/10.1016/j.fishres.2021.105967)

### Major improvements

* New ability to plot results from multiple ASAP3 and WHAM models together for comparison. Thanks to @liz-brooks for contributing `read_asap3_fit`! See `?compare_wham_models` and updated vignettes.

### Minor improvements

* New `save.sdrep = F` option to only save `summary(sdreport)` instead of `sdreport`. Can make saved models MUCH smaller (e.g. 2 MB vs. 150 MB). ([2f8875](https://github.com/timjmiller/wham/commit/2f8875323c0d6845a92444a9e7d4aaa92fe29d8d)).
* Added `proj.opts$percentFXSPR` option, percent of F_XSPR to use for calculating catch in projections. For example, GOM cod uses F = 75% F_XSPR, so `proj.opts$percentFXSPR = 75`.
* Added `proj.opts$useFMSY` and `proj.opts$percentFMSY` options, to project population and catches at (a percentage of) F_MSY.

### Bug fixes

* Broken links in vignettes (thanks to @tcarruth)
* Plotting issue with Ecov OSA residuals and peels (thanks to @h-du-pontavice, [1163df](https://github.com/timjmiller/wham/commit/1163df4290387174a5336ada98a13dc0f5a9644c))
* Double logistic selectivity setup in `prepare_wham_input` (thanks to @tcarruth, [f270dd](https://github.com/timjmiller/wham/commit/f270ddb66d253ac2aaf9a0109631b89036ddcd5e))
* `tryCatch` error assignment issues ([9d5c87](https://github.com/timjmiller/wham/commit/9d5c8792769dae9e640da3ad44b7f3e6b74e4a87))
* Remove error in default projection options (`percentSPR`)
* Setting up projections with multiple Ecovs (again thanks to @h-du-pontavice, [42a6a4](https://github.com/timjmiller/wham/commit/42a6a4950e85613219525e89c5590c37f3a6369f))
* Selectivity parameter initial values set to middle of range if unspecified in `prepare_wham_input`
* Fixed `fit_tmb` to make `$final_gradient` reported by wham equivalent to `sdreport()$gradient.fixed`. Issue was that the `model$env$last.par.best` is not updated by the newton steps after optimization so `model$opt$par` was slightly different (10^(-7) or smaller). [e20bd8](https://github.com/timjmiller/wham/commit/e20bd8d01a32b1cbdb826905257555f9a8b55c75)
* Fixed default determination of age corresponding to fully selected (total) fishing mortality. Now which_F_age is a vector of annual values (including projection years) and the values are reassigned after fitting before reporting.
* Added check for F_XSPR to verify percent of SPR0 is correct and, in projections, catch at F from catch is verified. Initial values for these algorithms can be specified by the user. 

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

