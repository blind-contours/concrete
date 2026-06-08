# concrete 1.1.1.9000

## Clinical (death-priority) win ratio [experimental]

* New experimental `clinicalWinRatio()`: the death-priority win ratio that counts
  death **even when it follows a non-fatal event** (e.g. CV death > HF
  hospitalization) -- the estimand the competing-risks/first-event `getWinRatio()`
  cannot produce. It is built on a Markov illness-death model (alive -> non-fatal
  -> death, plus alive -> death), each transition estimated by a Super Learner
  (the post-non-fatal death hazard on a left-truncated risk set), with
  doubly-robust, covariate-adjusted, censoring-corrected (IPCW) influence-function
  inference. Returns the win ratio, win odds, net benefit, and P(win/loss/tie).
  Takes a one-row-per-subject multistate data frame (its own interface for now,
  not the formatArguments pipeline); assumes a single non-fatal event type and
  conditionally-independent censoring.
* Validated end-to-end against ground truth (closed-form path probabilities and a
  brute-force pairwise win ratio on full simulated histories): point estimate
  recovers truth (~1%), and 95% CI coverage is nominal with and without ~36%
  censoring. See the new "Win ratios for trialists" article and
  `scripts/make-clinical-wr-*.R`.
* Internal `R/getTransitionHazard.R`: a Super Learner for illness-death transition
  intensities on left-truncated risk sets (the building block), plus
  `multistateCurves()` (midpoint-quadrature path probabilities).

## Documentation: win ratios for trialists

* New "Win ratios for trialists" vignette explaining the single-event,
  hierarchical (prioritized), and clinical (death-priority) win ratios -- when to
  use each, how they work, and the simulation coverage evidence (validated vs
  ground truth) so the inference can be trusted.

## Documentation: RMST methods comparison

* New "Restricted mean survival time: how concrete compares to other packages"
  vignette. On a competing-risks toy trial with a closed-form true RMST, it puts
  `concrete` (`getRMST`/`targetRMST`) head-to-head with `survRM2` (unadjusted and
  Tian-adjusted) and `eventglm`, and describes `riskRegression::ate`. It reports,
  truthfully, parity on the event-free RMST difference; the genuine differentiators
  (direct RMST-equation targeting, competing-risks years-lost that `survRM2`
  cannot do, ML + cross-fitted double robustness, and the integrated estimand
  toolkit); where the established tools are better (maturity, speed, simplicity);
  and the grid-discretization caveat for `getRMST` (use a dense `TargetTime` grid
  and/or `targetRMST`). Adds `survRM2` and `eventglm` to Suggests.

## Bug fixes from package audit

* **`formatArguments()` no longer mutates the caller's `data.table`.**
  `formatDataTable()` assigned the ID column and reordered columns in place, so a
  user's own `data.table` was silently modified by reference on every call
  (column order changed; ID column added). It now operates on a copy.
* **Additive-hazards (`aareg`) learner: corrected baseline accumulation.**
  `survival::aareg()` returns the per-event-time coefficient *increments*
  `dB(t)`, but `predictAaregHazLearner()` treated them as the cumulative
  coefficients `B(t)`, producing near-zero, oscillating hazards. The
  coefficients are now cumulated (`cumsum`) before forming the cumulative hazard,
  so `aareg`-based cause-specific and censoring hazards (and the IPCW weights and
  estimates that depend on them) are correct. Other learners were unaffected.

## Hierarchical (prioritized) win ratio

* `getWinRatio()` now accepts `TargetEvent` as an ordered vector of event codes
  (highest priority first) and computes the **prioritized win ratio, win odds,
  and net benefit** over a hierarchy of competing events (e.g.
  `TargetEvent = c(1, 2)` for death > hospitalization). The win/loss
  probabilities remain smooth functionals of the per-arm cause-specific
  cumulative incidence curves, so inference is still covariate-adjusted,
  doubly-robust, and censoring-corrected via the influence-function delta method.
  The single-event call (one `TargetEvent`) is unchanged and is the exact K=1
  special case. The prioritized rule compares each patient's *first* event
  (treating the listed events as competing risks); events after a patient's first
  event are not used (the fully semi-competing version is future work). Validated
  against a brute-force pairwise simulation of the same rule.

## Documentation: trial-design and regulatory toolkit

* New "Trial-design and regulatory toolkit" vignette and README section that
  surface the regulatory-facing features in one place: cross-fitting
  (`CrossFit = TRUE`), the ensemble hazard Super Learner (`HazEnsemble = TRUE`),
  the win ratio / win odds / net benefit (`getWinRatio()`), the ICH E9(R1)
  estimand framework and intercurrent-event strategies (`makeEstimand()`,
  `applyIntercurrentEvent()`), censoring sensitivity (`senseCensoring()`), and
  adjustment efficiency (`getRelativeEfficiency()`). These were documented at the
  function level but were not discoverable from the README or the article index.

## Coxnet baseline-hazard fix

* The Coxnet (penalized Cox) hazard learner had the same baseline-hazard
  mismatch as the Cox learner: it took the baseline from a separate
  treatment-only model while using the glmnet linear predictor, mis-scaling the
  covariate-adjusted Coxnet hazard. The baseline is now the Breslow estimator
  computed from the glmnet linear predictor itself (centered), so the conditional
  hazard is reconstructed consistently. With the fix, a lightly-penalized Coxnet
  reproduces the Cox estimate. (The random-survival-forest, additive-hazards, and
  HAL learners estimate the conditional hazard directly and were unaffected.)

## Win ratio

* `getWinRatio()` estimates the covariate-adjusted **restricted win ratio**,
  **win odds**, and **net benefit** for a single terminal time-to-event outcome,
  as functionals of the targeted counterfactual survival curves. Because the
  win/loss probabilities are smooth functionals of those curves, inference uses
  the influence functions already produced by the estimator (delta method),
  giving doubly-robust, covariate-adjusted, censoring-corrected win statistics --
  unlike the standard unadjusted, censoring-sensitive win ratio. A hierarchical /
  competing-risk win ratio is planned.

## Hazard learner: ensemble Super Learner and a baseline-hazard fix

* **Bug fix (Cox hazards):** the Cox hazard learner took its baseline hazard
  from a separate treatment-only model while using the full model's linear
  predictor, which mis-scaled every covariate-adjusted Cox hazard (the
  misspecified treatment-only candidate could even win cross-validation). The
  baseline is now taken from the fitted model itself, so the conditional hazard
  is reconstructed consistently. This removes a finite-sample bias in
  covariate-adjusted estimates and lets the correct hazard model be selected.
* `formatArguments(..., HazEnsemble = TRUE)` combines the candidate hazard
  learners into a cross-validated convex-combination ensemble (Super Learner),
  minimizing the counting-process negative log-likelihood of the weighted hazard
  over the simplex, instead of discrete winner-take-all selection -- matching how
  the treatment propensity is already estimated.

## Censoring sensitivity analysis

* `senseCensoring()` adds a tipping-point sensitivity analysis for the
  independent-censoring (conditional MAR) assumption: a fraction `delta` of the
  subjects censored before the target time are imputed as having experienced the
  event of interest and the analysis is re-fit, tracing the estimate and CI from
  the optimistic (`delta = 0`) to the pessimistic (`delta = 1`) bound and
  reporting the tipping point. (Scaling the censoring weight leaves a
  doubly-robust target unchanged, so imputation -- which changes the estimand --
  is used instead.)

## ICH E9(R1) estimand framework and intercurrent events

* `makeEstimand()` records the analysis target by the five ICH E9(R1) attributes
  (treatment, population, endpoint, intercurrent-event strategy, summary), and
  travels with the results for a statistical analysis plan.
* `applyIntercurrentEvent()` implements the data handling for the
  intercurrent-event strategy: `"treatment policy"` (intent-to-treat, the
  default), `"hypothetical"` (recode the intercurrent event as censoring so the
  existing IPCW targets the no-intercurrent-event risk), and `"composite"`
  (merge the intercurrent event into the event of interest).

## Cross-fitting (CV-TMLE)

* `formatArguments(..., CrossFit = TRUE)` estimates the propensity and the
  cause-specific / censoring hazards by cross-fitting: each subject's nuisances
  are predicted from models fit on the other folds (with an inner
  cross-validation for hazard-learner selection). This supports valid
  influence-function inference when flexible machine-learning learners are used
  -- the regime where the in-sample fit can otherwise undercover -- and is the
  principled answer to regulatory caution around ML-based covariate adjustment.
  The targeting and inference are unchanged; only how the nuisances are produced
  differs. The out-of-fold construction is verified in the test suite.

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
