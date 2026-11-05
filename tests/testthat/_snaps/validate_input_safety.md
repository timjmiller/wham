# R guards cover vector lengths, random-walk predecessors and projections

    Code
      validate_input_safety(broken)
    Condition
      Error:
      ! Unsafe WHAM input: NAA_re_model must have one entry per stock.

---

    Code
      validate_selpars_re(broken)
    Condition
      Error:
      ! Unsafe WHAM input: selblock_models must have one entry per selectivity block.

---

    Code
      validate_input_safety(broken)
    Condition
      Error:
      ! Unsafe WHAM input: Ecov random walk requires a preceding year.

---

    Code
      validate_input_safety(broken)
    Condition
      Error:
      ! Unsafe WHAM input: Ecov_model contains invalid option codes.

---

    Code
      validate_input_safety(broken)
    Condition
      Error:
      ! Unsafe WHAM input: n_years_proj must be a nonnegative integer within the supported year range.

---

    Code
      validate_input_safety(broken)
    Condition
      Error:
      ! Unsafe WHAM input: proj_Ecov_opt must have one valid option per Ecov.

---

    Code
      validate_input_safety(broken)
    Condition
      Error:
      ! Unsafe WHAM input: Ecov process does not cover projection years.

---

    Code
      validate_input_safety(broken)
    Condition
      Error:
      ! Unsafe WHAM input: Ecov projection averaging years are missing or out of bounds.

---

    Code
      validate_input_safety(broken)
    Condition
      Error:
      ! Unsafe WHAM input: Ecov_use_proj must cover projection years and covariates.

