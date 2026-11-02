# Historical Refit Checks

The verification counts later in this document describe earlier runs, not a
guarantee for the current workspace. The November 2026 rerun exposed steepness
rejections and cached v1.0.0 fixtures whose steepness option was inactive.
The fixes preserve the original steepness parameterization in the native model
and regenerate missing or inactive fixtures using the release's array layout.
Consult the current rerun reports before claiming release-by-feature support.

Rebuild WHAM before testing the native steepness fix. Local validation uses
`oldfitobjects/refit-fix-library` without replacing the user's installed package:

```text
R CMD INSTALL --library=oldfitobjects/refit-fix-library .
```

Set `WHAM_CURRENT_LIBRARY` to that directory for the commands below. The native
compatibility modes use internal recruitment codes 5/6 and retain the saved
steepness/R0 values and maps. Their coefficients are recalculated each objective
evaluation using annual mortality, maturity, spawning weights and spawn timing.
They remain restricted to one stock, region and season without projections.
The native regressions independently check the historical BH/Ricker formulas,
including perturbations to mortality and year-specific spawning inputs.

The shared optimizer and its Newton refinement are unchanged. Weakly identified
historical models can still report singular-Hessian warnings; likelihood
agreement does not establish identifiability or valid standard errors.

The local `oldfitobjects/test_refit_releases.R` runner checks all semantic-version
Git tags, reusing historical fixtures when available and generating missing
ones. It runs the nine baseline cases and thirteen feature cases independently
for each release, continues after failed stages, and writes `rerun_status.csv`
and `rerun_coverage.csv` under `oldfitobjects/release_refits/`:

```text
Rscript oldfitobjects/test_refit_releases.R
Rscript oldfitobjects/test_refit_releases.R v1.0.6
Rscript oldfitobjects/test_refit_releases.R --regenerate
```

The runner exits unsuccessfully if any stage or required case fails. Missing
fixtures, inactive steepness and incomplete refits are not counted as passes. This tests a defined
matrix across tagged releases, not every historical development commit or
every feature combination. Migration rejections, numerical failures, and
historical generation failures must be distinguished in the per-stage logs.
Use a current native build consistent with the workspace R code; an older
installed DLL can invalidate comparisons.

The 22-case matrix is not a full WHAM feature inventory. In particular it does
not establish cross-version support for catchability effects, movement,
multi-stock/multi-region models, length models, every age-composition
likelihood, stored projections, OSA, posterior sampling, retrospectives, or
uncertainty/reference-point workflows. Separate saved-fit regressions may cover
some of these, but must not be counted as release-by-feature coverage.

Run from the repository root in separate Rscript processes:

```text
Rscript tests/manual/refit-releases.R generate v1.0.0
Rscript tests/manual/refit-releases.R refit v1.0.0
Rscript tests/manual/refit-releases.R refit-defaults v1.0.0
Rscript tests/manual/refit-releases.R run v1.0.1 v1.0.2 v1.0.3 v1.0.4 v1.0.5 v1.0.6
Rscript tests/manual/refit-releases.R summarize
```

Set `WHAM_CURRENT_LIBRARY` to the library directory containing the current WHAM
installation if it is not on R's library path. The refit stages load the installed
native model and the workspace migration, refit and fitting functions; the
installed R functions are not modified.

`generate` requires Git, the local release tag, TMB, dplyr, tidyr, and a working R
C++ toolchain. It extracts the release, compiles an unchanged build copy of its
model source, and sources its historical R functions. It does not install an old
WHAM package or replace the current one. The toolchain is current, not a
recreation of the release's original R/TMB environment.

The harness evaluates the input-construction expressions from each tag's own
example-1 test, retaining that release's selectivity and likelihood settings.
The four original example-1 configurations cover recruitment-only and full
numbers-at-age random effects, each with multinomial and logistic-normal
age compositions. Generated likelihoods are checked against the release's
archived test results. Five additional cases cover SCAA, recruitment AR1 over
years, abundance AR1 over ages, abundance AR1 over years, and two-dimensional
abundance AR1. These additional cases have no archived reference optimum.
`refit` checks the migrated likelihood at the saved fixed-parameter estimates
before optimization, as well as after refitting. It disables retrospectives, standard errors, OSA,
reference points, and posterior sampling to isolate the migration and likelihood
check. `refit-defaults` leaves the saved call unchanged to exercise default
behavior. Warnings are retained separately; NLL agreement does not imply a
well-conditioned Hessian or successful uncertainty calculations.

Generation also evaluates each historical DLL at the saved fixed and latent
parameter values and saves a `_reference.rds` report for biological comparisons.
Some early fits cache reports that are inconsistent with their parameter lists,
so these checks use the historical DLL as the reference rather than assuming all
cached report fields describe the same parameter state. `reference <tag>` can
create these reports for existing fixtures without refitting. The original saved
NLL remains the benchmark for the marginal-likelihood checks.

Outputs are written under `oldfitobjects/release_refits/<tag>/`, configurable
with `WHAM_RELEASE_RESULTS`: provenance, original fits, fresh refits, comparison
text, and CSV summaries. No refit results are reused. The NLL tolerance is `1e-6`
absolute. The regression tests in `test_prepare_wham_input_from_fit.R` also check
the v1.0.0 through v1.0.6 fixtures when present; they skip absent fixtures.
The `run` stage launches separate generation and refit processes, saves their
logs, and stops on the first failure. Generation is repeatable but not cached;
use the `refit` stage alone to rerun migration against already-generated fits.

Other tags may be supplied, but a successful v1.0.0 run alone does not establish
coverage of every release or model configuration. Stored projections and
estimated early mortality baselines with age-specific offsets remain unsupported
by the legacy migration. Environmental cross-effect coefficient sharing,
ambiguous process conventions, and early observation-error random effects
are explicitly rejected.
For the early packed age-composition parameter layout, supported likelihood codes
are 1 (multinomial), 2 (Dirichlet-multinomial), and 7 (logistic-normal with zeros
missing). Historical code 7 becomes current code 5. Stored likelihood components
determine whether to retain the original effective sample sizes (weighted
variance) or use unit effective sample sizes (unweighted variance). The originals
are retained in `input$legacy_age_comp_Neff`. Missing or ambiguous evidence and
other early likelihood codes are explicitly rejected.
The later matrix parameter layout retains its likelihood codes and sample-size
weights independently of whether observation indices need rebuilding.

`summarize` checks all six requested releases by default and writes
`oldfitobjects/release_refits/release_summary.csv`. It requires all nine cases per
release, four archived NLL matches, matching migrated and refitted NLLs, optimizer
convergence, and no migration errors. Warnings are counted separately.

Verified with WHAM 2.1.0.9012 and R 4.4.3: all 54 cases across v1.0.1 through
v1.0.6 pass the `1e-6` absolute NLL tolerance before and after refitting. The
largest refit difference was `1.65e-8`. Five v1.0.1 cases retain singular-Hessian
warnings; v1.0.2 through v1.0.6 produced no refit warnings in this matrix.

## Legacy Feature Checks

```text
Rscript tests/manual/refit-releases.R generate-features v1.0.6
Rscript tests/manual/refit-releases.R refit-features v1.0.6
Rscript tests/manual/refit-releases.R summarize-features
```

Feature stages write separate results under `<tag>/features/`. Optional trailing
case names restrict generation without refitting the other fixtures. Run refits
against a rebuilt current package: native compatibility code is needed for the
early environmental AR1 equation and for BH/Ricker steepness. Loading only the
updated R functions is insufficient.

The feature matrix covers IID, correlated and shared-annual selectivity effects;
age-specific and constant mortality with annual or two-dimensional deviations;
environmental RW/AR1 effects on recruitment and mortality; allometric mortality
with and without shared annual intercept deviations; and BH/Ricker steepness.
Allometric fixtures use a separate, controlled positive mortality-weight matrix,
without changing catch or survey weights. Their annual process variance and
correlation are fixed to provide stable historical test fits. The original nine
feature cases are tested on v1.0.0, v1.0.5 and v1.0.6, spanning the main coefficient
layouts. These examples do not establish support for every possible combination.

Steepness retains the original parameterization rather than converting only its
starting values to alpha/beta. Its yearly unfished SPR changes with mortality.
The legacy mode is restricted to single-stock, single-region, single-season
models without projections. The biological regression checks use independent
reports evaluated by the historical DLL at saved fixed and latent parameters.

Verified on v1.0.0, v1.0.5 and v1.0.6: all 39 feature cases pass both the saved-
parameter and refitted NLL checks (`1e-6` absolute tolerance). Maximum differences
are `5.18e-7` at saved parameters and `2.82e-8` after refitting. All refits report
optimizer convergence. Eight v1.0.0 cases retain numerical warnings, and one
original v1.0.0 steepness fit did not report convergence; NLL agreement does not
remove these conditioning concerns. The later two releases produced no warnings
in this feature matrix. Results are saved in
`oldfitobjects/release_refits/feature_summary.csv`.

The focused migration regression run passes 1,292 assertions, including the
original release fixtures, the 39 feature fixtures, historical AR1 simulation,
and steepness formula checks that perturb mortality. Both v1.0.6 steepness models
also complete reference-point and standard-error calls successfully. OSA,
posterior sampling, stored projections, and every possible combination of these
features are not covered by this matrix.