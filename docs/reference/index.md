# Package index

## Exported functions that user can call directly

- [`check_convergence()`](https://timjmiller.github.io/wham/reference/check_convergence.md)
  : Check convergence of WHAM model

- [`check_estimability()`](https://timjmiller.github.io/wham/reference/check_estimability.md)
  :

  Check for identifiability of fixed effects Originally provided by the
  [TMBhelper](https://github.com/kaskr/TMB_contrib_R/tree/master/TMBhelper)
  package. Internal function called by
  [`fit_tmb`](https://timjmiller.github.io/wham/reference/fit_tmb.html).

- [`compare_wham_models()`](https://timjmiller.github.io/wham/reference/compare_wham_models.md)
  : Compare multiple WHAM (or ASAP) models

- [`do_reference_points()`](https://timjmiller.github.io/wham/reference/do_reference_points.md)
  : Add reporting of biological reference points to WHAM model

- [`do_retro_peels()`](https://timjmiller.github.io/wham/reference/do_retro_peels.md)
  : Fit retrospective peels and add them to the fitted model object

- [`do_sdreport()`](https://timjmiller.github.io/wham/reference/do_sdreport.md)
  : Add TMB sdreport object to WHAM model

- [`fit_peel()`](https://timjmiller.github.io/wham/reference/fit_peel.md)
  :

  Fit model peeling off *i* years of data

- [`fit_tmb()`](https://timjmiller.github.io/wham/reference/fit_tmb.md)
  : Fit TMB model using nlminb

- [`fit_wham()`](https://timjmiller.github.io/wham/reference/fit_wham.md)
  : Fit WHAM model

- [`jitter_wham()`](https://timjmiller.github.io/wham/reference/jitter_wham.md)
  : Jitter starting values of a fitted WHAM model

- [`make_osa_residuals()`](https://timjmiller.github.io/wham/reference/make_osa_residuals.md)
  : Calculate one-step-ahead residuals

- [`mohns_rho()`](https://timjmiller.github.io/wham/reference/mohns_rho.md)
  : Calculate Mohn's rho for a WHAM model with peels

- [`plot_wham_output()`](https://timjmiller.github.io/wham/reference/plot_wham_output.md)
  : Plot WHAM output

- [`prepare_projection()`](https://timjmiller.github.io/wham/reference/prepare_projection.md)
  : Prepare input data and parameters to project an already fit WHAM
  model

- [`prepare_wham_input()`](https://timjmiller.github.io/wham/reference/prepare_wham_input.md)
  : Prepare input data and parameters for WHAM model

- [`project_wham()`](https://timjmiller.github.io/wham/reference/project_wham.md)
  : Project a fit WHAM model

- [`read_asap3_dat()`](https://timjmiller.github.io/wham/reference/read_asap3_dat.md)
  : Read an ASAP3 .dat file into R

- [`read_asap3_fit()`](https://timjmiller.github.io/wham/reference/read_asap3_fit.md)
  : Read ASAP3 fit

- [`read_wham_fit()`](https://timjmiller.github.io/wham/reference/read_wham_fit.md)
  : Read WHAM fit

- [`retro()`](https://timjmiller.github.io/wham/reference/retro.md) :
  Run retrospective analysis

- [`self_test()`](https://timjmiller.github.io/wham/reference/self_test.md)
  : Perform self-test simulation and estimation of a fitted WHAM model

- [`set_age_comp()`](https://timjmiller.github.io/wham/reference/set_age_comp.md)
  : Specify the age composition models for fleet(s) and indices.

- [`set_catch()`](https://timjmiller.github.io/wham/reference/set_catch.md)
  : Specify catch selectivity blocks and aggregate and age composition
  observations for catch

- [`set_ecov()`](https://timjmiller.github.io/wham/reference/set_ecov.md)
  : Specify configuration for environmental covariates, effects on the
  population, and parameter values

- [`set_F()`](https://timjmiller.github.io/wham/reference/set_F.md) :
  Specify configuration for fully-selected fishing mortality

- [`set_indices()`](https://timjmiller.github.io/wham/reference/set_indices.md)
  : Specify index selectivity blocks and aggregate and age composition
  observations for indices

- [`set_L()`](https://timjmiller.github.io/wham/reference/set_L.md) :
  Specify model and parameter configuration for "extra" mortality not
  directly attributed to natural mortality

- [`set_M()`](https://timjmiller.github.io/wham/reference/set_M.md) :
  Specify model and parameter configuration for natural mortality

- [`set_move()`](https://timjmiller.github.io/wham/reference/set_move.md)
  : Specify model and parameter configuration for movement when
  input\$data\$n_regions \> 1

- [`set_NAA()`](https://timjmiller.github.io/wham/reference/set_NAA.md)
  : Specify model and parameter configuration for numbers at age

- [`set_osa_obs()`](https://timjmiller.github.io/wham/reference/set_osa_obs.md)
  : Set up observation vector that is used by the model for likelihood
  calculations and one-step-ahead residuals.

- [`set_q()`](https://timjmiller.github.io/wham/reference/set_q.md) :
  Specify model and parameter configuration for catchability

- [`set_selectivity()`](https://timjmiller.github.io/wham/reference/set_selectivity.md)
  : Specify model and parameter configuration for selectivity

## Unexported functions

- [`extract_fixed()`](https://timjmiller.github.io/wham/reference/extract_fixed.md)
  :

  Extract fixed effects Originally provided by the
  [TMBhelper](https://github.com/kaskr/TMB_contrib_R/tree/master/TMBhelper)
  package. Internal function called by
  [`check_estimability`](https://timjmiller.github.io/wham/reference/check_estimability.html).

- [`reduce_input()`](https://timjmiller.github.io/wham/reference/reduce_input.md)
  : Reduce the years of the model

- [`wham_html()`](https://timjmiller.github.io/wham/reference/wham_html.md)
  : Create HTML file to view output plots in browser
