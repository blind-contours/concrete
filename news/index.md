# Changelog

## concrete 1.1.1.9000

### Win ratio

- [`getWinRatio()`](https://blind-contours.github.io/concrete/reference/getWinRatio.md)
  estimates the covariate-adjusted **restricted win ratio**, **win
  odds**, and **net benefit** for a single terminal time-to-event
  outcome, as functionals of the targeted counterfactual survival
  curves. Because the win/loss probabilities are smooth functionals of
  those curves, inference uses the influence functions already produced
  by the estimator (delta method), giving doubly-robust,
  covariate-adjusted, censoring-corrected win statistics – unlike the
  standard unadjusted, censoring-sensitive win ratio. A hierarchical /
  competing-risk win ratio is planned.

### Hazard learner: ensemble Super Learner and a baseline-hazard fix

- **Bug fix (Cox hazards):** the Cox hazard learner took its baseline
  hazard from a separate treatment-only model while using the full
  model’s linear predictor, which mis-scaled every covariate-adjusted
  Cox hazard (the misspecified treatment-only candidate could even win
  cross-validation). The baseline is now taken from the fitted model
  itself, so the conditional hazard is reconstructed consistently. This
  removes a finite-sample bias in covariate-adjusted estimates and lets
  the correct hazard model be selected.
- `formatArguments(..., HazEnsemble = TRUE)` combines the candidate
  hazard learners into a cross-validated convex-combination ensemble
  (Super Learner), minimizing the counting-process negative
  log-likelihood of the weighted hazard over the simplex, instead of
  discrete winner-take-all selection – matching how the treatment
  propensity is already estimated.

### Censoring sensitivity analysis

- [`senseCensoring()`](https://blind-contours.github.io/concrete/reference/senseCensoring.md)
  adds a tipping-point sensitivity analysis for the
  independent-censoring (conditional MAR) assumption: a fraction `delta`
  of the subjects censored before the target time are imputed as having
  experienced the event of interest and the analysis is re-fit, tracing
  the estimate and CI from the optimistic (`delta = 0`) to the
  pessimistic (`delta = 1`) bound and reporting the tipping point.
  (Scaling the censoring weight leaves a doubly-robust target unchanged,
  so imputation – which changes the estimand – is used instead.)

### ICH E9(R1) estimand framework and intercurrent events

- [`makeEstimand()`](https://blind-contours.github.io/concrete/reference/makeEstimand.md)
  records the analysis target by the five ICH E9(R1) attributes
  (treatment, population, endpoint, intercurrent-event strategy,
  summary), and travels with the results for a statistical analysis
  plan.
- [`applyIntercurrentEvent()`](https://blind-contours.github.io/concrete/reference/makeEstimand.md)
  implements the data handling for the intercurrent-event strategy:
  `"treatment policy"` (intent-to-treat, the default), `"hypothetical"`
  (recode the intercurrent event as censoring so the existing IPCW
  targets the no-intercurrent-event risk), and `"composite"` (merge the
  intercurrent event into the event of interest).

### Cross-fitting (CV-TMLE)

- `formatArguments(..., CrossFit = TRUE)` estimates the propensity and
  the cause-specific / censoring hazards by cross-fitting: each
  subject’s nuisances are predicted from models fit on the other folds
  (with an inner cross-validation for hazard-learner selection). This
  supports valid influence-function inference when flexible
  machine-learning learners are used – the regime where the in-sample
  fit can otherwise undercover – and is the principled answer to
  regulatory caution around ML-based covariate adjustment. The targeting
  and inference are unchanged; only how the nuisances are produced
  differs. The out-of-fold construction is verified in the test suite.

### New RCT estimands and reporting

- Added
  [`getRMST()`](https://blind-contours.github.io/concrete/reference/getRMST.md)
  for restricted mean survival time (event-free) and cause-specific
  life-years lost. These are collapsible, clinically interpretable
  estimands obtained by integrating the targeted cumulative incidence
  over the fitted target-time grid; their influence functions are the
  time-integrals of the absolute-risk influence functions, so standard
  errors, differences, and ratios reuse the existing machinery.
- Added
  [`targetRMST()`](https://blind-contours.github.io/concrete/reference/targetRMST.md),
  which targets the RMST / life-years-lost estimand *directly* with the
  integrated clever covariate rather than integrating pointwise-targeted
  risks. It fluctuates the fitted hazards on a rescaled time axis until
  the RMST estimating equation is solved, which is better conditioned
  and tends to converge and cover better than the pointwise approach for
  sparse grids, rare events, and long horizons.
- [`getOutput()`](https://blind-contours.github.io/concrete/reference/getOutput.md)
  now reports a two-sided Wald `pValue` for the comparative
  (risk-difference and risk-ratio) estimands, and gains optional
  `NIMargin` / `NIDirection` arguments for one-sided non-inferiority
  assessment.
- Added
  [`getRelativeEfficiency()`](https://blind-contours.github.io/concrete/reference/getRelativeEfficiency.md)
  to quantify the precision gain from covariate adjustment versus an
  unadjusted analysis (relative efficiency, percentage variance
  reduction, and effective sample-size multiplier), as in the FDA 2023
  covariate-adjustment guidance.

### Trialist beta updates

- Updated convergence guidance for rare-event analyses. The default
  stopping rule remains the original relative empirical EIC rule, but
  the recommended rare-event sensitivity is now
  `EICStopRule = "absolute"` with
  `EICStopAbsTol = 0.02 / sqrt(nrow(data))`.
- Updated README, vignettes, help files, and the installed trialist
  smoke test to use the sample-size-scaled absolute stopping
  sensitivity.
- Expanded the convergence issue template so trialist beta testers
  report the learner library and, when possible, a relative-vs-absolute
  stopping comparison.

### Simulation support

- Added project-relative referee simulation scripts for event-1,
  failed-seed, rare-event stopping, and alternative convergence-method
  comparisons.
- Added notes summarizing the full survival learner rare-event
  validation. In the primary hard-seed screen with Cox/Coxnet, random
  survival forests, additive hazards, and HAL, absolute stopping
  converged on all hard seeds and avoided a relative-rule failure driven
  by a near-zero-variance rare-event component.

## concrete 1.1.0

### New features

- Added user-facing TMLE convergence controls through
  [`formatArguments()`](https://blind-contours.github.io/concrete/reference/formatArguments.md):
  `UpdateMethod`, `EICStopRule`, and `EICStopAbsTol`.
- Added hybrid and absolute empirical EIC stopping rules. The default
  remains the original relative rule.
- Added
  [`getTmleDiagnostics()`](https://blind-contours.github.io/concrete/reference/getTmleDiagnostics.md)
  for inspecting final component-wise empirical EIC diagnostics, update
  traces, and norm trajectories from fitted `ConcreteEst` objects.
- Added documented hazard learner aliases for Coxnet, random survival
  forests, additive hazards, and HAL.
- Added trialist-facing documentation, a pkgdown site, GitHub issue
  templates, and an installed smoke-test script for checking Cox-only
  and optional learner workflows.

### Stability

- Added trace fields for relative and absolute EIC ratios, maximum
  absolute empirical EIC, failing components, and stopping-rule
  metadata.
- Improved rare-event convergence behavior with the hybrid stopping rule
  while preserving the original default behavior.
- Restricted convergence stopping checks to the requested target
  event/time components, excluding internal complement rows used in EIC
  summaries.
- Aligned adaptive update acceptance with the active stopping rule:
  relative stopping uses the target EIC norm, while absolute and hybrid
  stopping use the component-wise stopping ratio.
- Added a zero-norm hazard-update guard and default smoke-test execution
  in the test suite.
