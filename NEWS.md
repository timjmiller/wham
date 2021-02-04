wham 1.0.2.9000
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

