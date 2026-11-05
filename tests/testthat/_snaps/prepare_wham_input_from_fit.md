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

# legacy decoupled Cod migration uses corrected variances

    Code
      input <- suppressMessages(migrate_legacy_refit_input(fit))
    Condition
      Warning:
      WHAM refit correction: historical decoupled recruitment used the wrong correlations in NAA marginal variances. The original fit used an affected historical calculation. Refitting uses the current implementation; the corrected NLL and parameter estimates may differ from the original fit.

---

    Code
      output <- capture.output(result <- suppressMessages(refit(fit)))
    Condition
      Warning:
      WHAM refit correction: historical decoupled recruitment used the wrong correlations in NAA marginal variances. The original fit used an affected historical calculation. Refitting uses the current implementation; the corrected NLL and parameter estimates may differ from the original fit.

---

    Code
      inferred <- suppressMessages(migrate_legacy_refit_input(fit))
    Condition
      Warning:
      WHAM refit correction: historical decoupled recruitment used the wrong correlations in NAA marginal variances. The original fit used an affected historical calculation. Refitting uses the current implementation; the corrected NLL and parameter estimates may differ from the original fit.

