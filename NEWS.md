wham 1.0.1.9000
=========================

### Minor improvements

* add `1e-15` to predicted proportions to make age composition likelihoods robust to 0 predictions when selAA is fixed at 0. This affects the multinomial, Dirichlet, and Dirichlet-multinomial (options 1-3), since the logistic normal (options 4-7) already did this. [88f15d4](https://github.com/timjmiller/wham/commit/88f15d4a51f69a3d649d76bcac0a8cf299c3135e)
* can now specify age composition model using `age_comp` argument to `prepare_wham_input`. See [`?prepare_wham_input`](https://timjmiller.github.io/wham/reference/prepare_wham_input.html) for details. [fd94b3d](https://github.com/timjmiller/wham/commit/fd94b3dcaf189482e10a6750c2f1b8350837fd48)

### Bug fixes

* check for sel par inits outside lower/upper bounds ([4ed394b](https://github.com/timjmiller/wham/commit/4ed394bac38d6054727bc4d3e17c8f4452ae8289))

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

