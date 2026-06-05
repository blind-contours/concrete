# concrete 1.1.1.9000

## New RCT estimands and reporting

* Added `getRMST()` for restricted mean survival time (event-free) and
  cause-specific life-years lost. These are collapsible, clinically
  interpretable estimands obtained by integrating the targeted cumulative
  incidence over the fitted target-time grid; their influence functions are the
  time-integrals of the absolute-risk influence functions, so standard errors,
  differences, and ratios reuse the existing machinery.
* Added `targetRMST()`, which targets the RMST / life-years-lost estimand
  *directly* with the integrated clever covariate rather than integrating
  pointwise-targeted risks. It fluctuates the fitted hazards on a rescaled time
  axis until the RMST estimating equation is solved, which is better conditioned
  and tends to converge and cover better than the pointwise approach for sparse
  grids, rare events, and long horizons.
* `getOutput()` now reports a two-sided Wald `pValue` for the comparative
  (risk-difference and risk-ratio) estimands, and gains optional `NIMargin` /
  `NIDirection` arguments for one-sided non-inferiority assessment.
* Added `getRelativeEfficiency()` to quantify the precision gain from covariate
  adjustment versus an unadjusted analysis (relative efficiency, percentage
  variance reduction, and effective sample-size multiplier), as in the FDA 2023
  covariate-adjustment guidance.

## Trialist beta updates

* Updated convergence guidance for rare-event analyses. The default stopping
  rule remains the original relative empirical EIC rule, but the recommended
  rare-event sensitivity is now `EICStopRule = "absolute"` with
  `EICStopAbsTol = 0.02 / sqrt(nrow(data))`.
* Updated README, vignettes, help files, and the installed trialist smoke test
  to use the sample-size-scaled absolute stopping sensitivity.
* Expanded the convergence issue template so trialist beta testers report the
  learner library and, when possible, a relative-vs-absolute stopping comparison.

## Simulation support

* Added project-relative referee simulation scripts for event-1, failed-seed,
  rare-event stopping, and alternative convergence-method comparisons.
* Added notes summarizing the full survival learner rare-event validation. In
  the primary hard-seed screen with Cox/Coxnet, random survival forests,
  additive hazards, and HAL, absolute stopping converged on all hard seeds and
  avoided a relative-rule failure driven by a near-zero-variance rare-event
  component.

# concrete 1.1.0

## New features

* Added user-facing TMLE convergence controls through `formatArguments()`:
  `UpdateMethod`, `EICStopRule`, and `EICStopAbsTol`.
* Added hybrid and absolute empirical EIC stopping rules. The default remains
  the original relative rule.
* Added `getTmleDiagnostics()` for inspecting final component-wise empirical EIC
  diagnostics, update traces, and norm trajectories from fitted `ConcreteEst`
  objects.
* Added documented hazard learner aliases for Coxnet, random survival forests,
  additive hazards, and HAL.
* Added trialist-facing documentation, a pkgdown site, GitHub issue templates,
  and an installed smoke-test script for checking Cox-only and optional learner
  workflows.

## Stability

* Added trace fields for relative and absolute EIC ratios, maximum absolute
  empirical EIC, failing components, and stopping-rule metadata.
* Improved rare-event convergence behavior with the hybrid stopping rule while
  preserving the original default behavior.
* Restricted convergence stopping checks to the requested target event/time
  components, excluding internal complement rows used in EIC summaries.
* Aligned adaptive update acceptance with the active stopping rule: relative
  stopping uses the target EIC norm, while absolute and hybrid stopping use the
  component-wise stopping ratio.
* Added a zero-norm hazard-update guard and default smoke-test execution in the
  test suite.
