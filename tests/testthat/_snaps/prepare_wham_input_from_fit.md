# refitting saved code 3 uses current AR1 without changing the saved input

    Code
      capture.output(result <- suppressMessages(refit(fit)))
    Condition
      Warning:
      Historical environmental AR1 is no longer supported; refitting with the current AR1 model. Likelihoods and estimates may change.
    Output
      [1] "  nll_original nll_refit nll_diff max_abs_par_diff conv_original conv_refit" "1            0         0        0                0             0          0"

# historical environmental AR1 migrates to the current process

    Code
      input <- suppressMessages(migrate_legacy_refit_input(fit))
    Condition
      Warning:
      Historical environmental AR1 is no longer supported; refitting with the current AR1 model. Likelihoods and estimates may change.

# legacy steepness is not silently treated as alpha and beta

    Code
      migrate_legacy_refit_input(fit)
    Condition
      Error:
      ! Cannot safely migrate this legacy WHAM input: steepness parameterization.

---

    Code
      migrate_legacy_refit_input(fit)
    Condition
      Error:
      ! Cannot safely migrate this legacy WHAM input: steepness parameterization.

