# equilibrium refits with selectivity random effects use the correction

    Code
      expect_identical(refit_initial_equilibrium(input, result$rep$N1 * 2), input)
    Condition
      Warning:
      Initial NAA equilibrium for stock(s) 1: stored initial abundance differs from the corrected calculation and may reflect the historical equilibrium bug. The model will be refit using the corrected equilibrium calculation; parameter estimates and likelihoods may change.
    Code
      expect_identical(refit_initial_equilibrium(input, array(0, c(1L, 1L, 1L))), input)
    Condition
      Warning:
      Initial NAA equilibrium for stock(s) 1: stored initial abundance has incompatible dimensions, so exposure to the historical equilibrium bug cannot be checked. The model will be refit using the corrected equilibrium calculation; parameter estimates and likelihoods may change.

# equilibrium warnings handle saved compatibility modes and missing reports

    Code
      expect_identical(refit_initial_equilibrium(input, NULL), corrected)
    Condition
      Warning:
      Initial NAA equilibrium for stock(s) 2: the saved input used the historical equilibrium bug. The model will be refit using the corrected equilibrium calculation; parameter estimates and likelihoods may change.
    Code
      expect_identical(refit_initial_equilibrium(corrected, NULL), corrected)
    Condition
      Warning:
      Initial NAA equilibrium for stock(s) 2: stored initial abundance is unavailable, so exposure to the historical equilibrium bug cannot be checked. The model will be refit using the corrected equilibrium calculation; parameter estimates and likelihoods may change.

