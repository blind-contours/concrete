# concrete 1.1.1.9000

## Estimator audit fixes

* **PRO tiers now inverse-probability-weight pre-horizon censoring.** The PRO
  block previously corrected only for landmark-visit non-attendance, not for
  censoring before the horizon, so its win/loss components were attenuated
  (~`G(tau)^2`) while the hard-tier reach was censoring-corrected --- biasing
  `clinicalWinRatio` / `clinicalPSNB` PRO tiers toward the null under meaningful
  censoring. Reachers are now weighted by `1/G(tau | W)` (the engine's `Ginv` at
  the horizon). Validated: the conditional PRO net benefit is now invariant to the
  censoring rate (0.25 across 0/18/45% pre-horizon censoring vs a drift to 0.08
  before), and the TRISCEND II 7-tier win-ratio bias falls from −0.10 to −0.05
  with coverage 0.967.
* **Ratio confidence intervals (`Rel Risk`) are now on the log scale** in
  `getOutput()` --- symmetric in log, always positive, better small-sample coverage
  (the win ratio already did this). Point estimates, SEs and Wald p-values are
  unchanged. Validated: null RR (=1) log-scale CI coverage 0.956. The core TMLE
  was re-confirmed well-calibrated in the same run (RD coverage 0.950, type-I
  0.050; per-arm risk 0.96–0.98).
* `getOutput()` simultaneous bands use 1e4 (was 1e3) multiplier draws for a stabler
  max-statistic quantile, matching `getSimultaneousFamily()`.

## Continuous / ordinal PRO tiers for the win statistics (TRISCEND II-style)

* `clinicalWinRatio()` and `clinicalPSNB()` gain a **`pro`** argument: continuous
  or ordinal patient-reported-outcome tiers (e.g. KCCQ, NYHA, 6-minute walk)
  measured at the final visit, appended at the **bottom** of the hierarchy (the
  clinical norm for soft markers; multiple PRO tiers may be stacked in priority
  order). This reproduces the **TRISCEND II** primary endpoint (Hahn et al., NEJM
  2025): `death > RV-assist/transplant > tricuspid reintervention > HF
  hospitalization > KCCQ >= 10 > NYHA >= 1 class > 6-min walk >= 30 m`.
* A pair reaches the PRO block iff tied on every hard-event tier (both event-free
  and alive at the horizon); within reach the PRO tiers are compared
  **sequentially**, each with a win margin `delta` and a direction
  (`higher.better` / `lower.better`). The block is estimated by a **reach-weighted,
  IPCW-corrected two-sample generalized pairwise comparison** on the joint marker
  vectors (pairs restricted to reachers; missing final visits inverse-probability
  weighted), so the markers may be arbitrarily correlated and the sequential
  tie-passing among PRO tiers is exact. Inference is the two-sample U-statistic
  (Hajek) influence function; hard-event tiers above keep the engine's
  covariate-adjusted, doubly-robust influence-function inference. Validated
  end-to-end against a brute-force pairwise win-ratio truth over the full 7-tier
  hierarchy. Experimental; the PRO landmark must equal the horizon (final-visit
  design), and a PRO ranked *above* a hard event is not supported.
* The multistate engine now supports a **death-only** hard-event hierarchy
  (`Kev = 1`), e.g. `death > KCCQ` with no non-fatal event tier.

## Restricted mean time in favor of treatment (RMT-IF)

* New exported **`getRMTIF()`** (first-event / competing risks) and
  **`clinicalRMTIF()`** (death-priority multistate, credits death after a
  non-fatal event): the average time over `[0, tau]` a treated patient spends in
  a more favorable state than a control --- the time-units companion to the win
  ratio (Mao 2023). `getRMTIF()` reduces **exactly** to the RMST difference for a
  single event (so RMST and the win statistics are two ends of one favor-time
  functional) and has closed-form influence-function inference; `clinicalRMTIF()`
  carries an analytic adjoint-value efficient influence function. Both validated
  against brute-force and closed-form truth (K=2/K=3 coverage 0.96--0.97,
  calibrated SEs, unbiased under cross-fitting). Feed `getSimultaneousFamily()`.

## Charter-weighted win ratio (priority-standardized net benefit)

* New exported **`clinicalPSNB()`**: the priority-standardized net benefit (PSNB)
  and win ratio (PSWR). Standard hierarchical win statistics aggregate the
  per-layer effects using implicit *reach* weights (set by the outcome/censoring
  distribution, not clinical priority), so a frequently-reached low-priority
  layer can dominate. `clinicalPSNB()` lets the user supply a prespecified
  clinical **charter** and recombines the layer-specific effects with it:
  `PSNB = sum_k alpha_k Delta_k`, `PSWR = sum_k alpha_k w_k / sum_k alpha_k l_k`.
  It reuses `clinicalWinRatio()`'s validated multistate engine (which already
  produces the per-tier components `W^(k)=r_k w_k`), divides the reach back out,
  and recombines with the charter, propagating influence functions by the delta
  method --- giving the covariate-adjusted, doubly-robust, IPCW-corrected
  estimator of these estimands. Reports each layer's reach and stage-conditional
  net benefit; `charter = "reach"` reproduces the standard net benefit / win
  ratio. Validated against a brute-force charter-weighted truth (coverage 0.967,
  SE ratio ~0.91 for both PSNB and PSWR).

## Joint simultaneous inference across the estimand family

* New exported **`getSimultaneousFamily()`**: family-wise simultaneous confidence
  bands across an arbitrary set of `concrete` estimands computed from the same
  fit --- e.g. a risk difference at several horizons together with the RMST
  difference and the win ratio. It stacks the per-subject efficient influence
  functions the estimators already produce into one influence matrix, estimates
  their joint correlation, and takes the simultaneous critical value from a
  Gaussian-multiplier bootstrap (`max_j |Z_j|`), so the bands control the
  family-wise error without the conservativeness of a Bonferroni correction
  (the estimands are strongly correlated, being functions of the same targeted
  curves). Ratio estimands are handled on the log scale; the bands inherit the
  covariate adjustment, censoring correction, cross-fitting, and
  stratified-randomization variance correction of the inputs. Every estimator
  (`getOutput`, `getRMST`, `targetRMST`, `getWinRatio`, `targetWinRatio`) now
  attaches the per-subject influence functions needed to feed it.

## Directly targeted win ratio

* New exported **`targetWinRatio()`**: a TMLE that fluctuates \strong{both arms'}
  cause-specific hazards to solve the efficient-influence-function estimating
  equations of the win and loss probabilities directly, instead of plugging
  pointwise-targeted risk curves into the win functional (`getWinRatio()`). The
  clever covariate is the win-gradient-weighted combination of the pointwise
  clever covariates over the \strong{full event-time grid} -- the chain-rule
  analogue of `targetRMST()`'s time-integrated clever covariate -- so the win
  integral no longer depends on the `TargetTime` grid. Because the gradient
  coefficients depend on both arms' current curves they are recomputed at every
  fluctuation step and the arms are updated jointly. Reports the win ratio, win
  odds, net benefit, and win/loss/tie probabilities with influence-function
  (and, with `Strata`, design-corrected) inference. Validated against
  brute-force pairwise ground truth (120 reps/cell, n = 500/arm): on a sparse
  2-time grid the plug-in win-ratio bias/coverage was +0.078 / 0.933 vs the
  direct method's +0.015 / 0.958; on a dense grid +0.042 / 0.967 vs
  +0.009 / 0.975, converging in 100% of replicates
  (`scripts/dev-targetwr-validation.R`). Derivation: `notes/target-win-ratio.md`.

## Variance correction for stratified / covariate-adaptive randomization

* `formatArguments()` gains a **`Strata`** argument naming the column(s) the
  trial actually randomized within (site, disease severity, ...). When supplied,
  **every** reported standard error --- absolute risk, risk difference and ratio,
  RMST / life-years-lost (both `getRMST()` and `targetRMST()`), and the win
  ratio / win odds / net benefit --- is corrected for covariate-adaptive
  randomization following Bugni--Canay--Shaikh / Ye--Shao: under permuted blocks
  or a stratified biased coin the usual iid influence-function variance is
  generically conservative, and ICH E9 / the FDA covariate-adjustment guidance
  ask the analysis to reflect the randomization scheme. The strata columns stay
  in the data as adjustment covariates; when the working models adjust for them
  the correction is approximately zero. Degenerate strata (fewer than 2 subjects
  in an arm) warn and fall back to the conservative iid SEs. Point estimates are
  unchanged. Validated by simulation (150 reps/cell, permuted blocks within 4
  prognostic strata): the correction shrinks SEs when the stratum is unadjusted,
  agrees with the iid SE to 4 decimals when it is adjusted, and never went
  anti-conservative (`scripts/dev-strata-validation.R`,
  `scripts/dev-strata-control.R`).

## Missing baseline covariates are handled

* `formatArguments()` previously stopped on any `NA` in the data. Baseline
  covariates with missing values are now **imputed** (numeric: median,
  categorical: mode) with a `<column>_missing` indicator added per affected
  column --- the handling endorsed by the FDA covariate-adjustment guidance for
  pre-randomization covariates. The event-time, event-type, treatment, and ID
  columns must still be complete.

## Treatment switching (crossover) and the hypothetical no-switching estimand

* `formatArguments()` gains a **`Crossover`** argument naming a per-subject
  switch-time column (`NA`/`Inf` = never switched). Each switcher's outcome is
  re-censored at the switch time and a **separate crossover hazard** is fit and
  multiplied into the IPCW alongside dropout (`1 / (S_dropout * S_crossover)`),
  targeting the hypothetical "no-switching" (per-protocol-type) estimand instead
  of naive per-protocol censoring. The crossover hazard inherits the censoring
  model's covariates, including `CensoringTV`.
* `senseCensoring()` gains **`mechanism = c("all", "dropout", "crossover")`** so
  the ICH E9(R1) tipping-point analysis can probe the dropout (MAR) and
  switching (no-switching counterfactual) assumptions individually or jointly
  --- each intercurrent-event assumption gets its own tipping point.

## Positivity / effective-sample-size diagnostics

* New exported **`getPositivityDx()`**: per-arm effective sample size (overall
  and worst-time), maximum inverse weight, minimum observation probability, and
  the share of weights pinned at the truncation bound, with a `CAUTION` flag ---
  surfacing when an IPCW estimand (especially the hypothetical no-switching one)
  rests on heavy extrapolation.

## New validation evidence

* **Null (type-I) simulation for the core absolute-risk TMLE** (160 reps,
  n = 500/arm, two competing events + censoring, no treatment effect): type-I
  error 0.050 at alpha = 0.05, risk-difference coverage at zero 0.950, SE
  calibration ratio 0.98 (`scripts/dev-null-typeI.R`).
* Crossover estimand construction checked against oracle IPCW weights
  (estimated ~ oracle) and a g-comp-vs-TMLE cut (TMLE removes ~half the g-comp
  bias); residual heavy-crossover bias is at-risk depletion that
  `getPositivityDx()` is designed to flag.

## Time-varying covariates in the censoring model [experimental]

* **Threaded through the core.** `formatArguments()` gains an optional `CensoringTV`
  argument (a long `data.frame` of time-varying covariates for the censoring model,
  with the `ID` column, a `time` column, and value columns). When supplied, the
  inverse-probability-of-censoring weight used by **every** estimand --- the
  targeted survival / cumulative-incidence curves, `getRMST`/`targetRMST`, and
  `getWinRatio` --- is conditioned on the LOCF value and change-from-baseline of
  each. Implemented as an override of the lagged censoring survival in
  `get(CV)InitialEstimate()` (the single object the censoring hazard feeds), so the
  outcome hazards and the rest of the TMLE are untouched; no `CensoringTV`
  reproduces the previous behavior exactly. This matters most for the **marginal**
  (non-ratio) estimands: in simulation a targeted absolute risk biased by
  informative censoring under baseline-only IPCW is recovered with the time-varying
  censoring model (`scripts/dev-tvcens-core.R`).
* `clinicalWinRatio()` gains optional `id` + `censoring.tv` arguments: a long
  `data.frame` of post-randomization, time-varying covariates (e.g. echo / KCCQ /
  6-minute-walk measured at follow-up visits) for the **censoring** model. When
  dropout is driven by such measurements, conditioning on baseline covariates alone
  makes censoring informative and biases the inverse-probability-of-censoring
  weights; conditioning the censoring hazard on the last-observation-carried-forward
  value and change-from-baseline of each restores conditional independence (CAR).
  They enter \strong{only} the censoring model (never the outcome hazards), so the
  marginal/ITT estimand is preserved (they are post-treatment mediators). Omitting
  them leaves the result unchanged.
* New shared internal machinery (`R/tvCensoring.R`): `.tvLOCF()` (LOCF value +
  change-from-baseline on the hazard grid) and `.tvCensoringInc()` (cross-fitted
  discrete-time censoring-hazard SL with time-varying covariates). Validated to
  remove informative-censoring bias: a survival probability biased +5% under
  baseline-only IPCW is recovered to within ~0.4% with the time-varying censoring
  model (`scripts/dev-tvcens*.R`). (Note: the win ratio, being a between-arm ratio,
  is far less sensitive to arm-symmetric informative censoring -- the larger payoff
  is for the marginal survival / cumulative-incidence curves and RMST, where the
  same machinery will be threaded next.)

## Hierarchical (death-priority) win ratio is now the flagship [experimental]

* `clinicalWinRatio()` is generalized from the two-tier illness-death case to an
  **arbitrary ordered hierarchy** of a terminal event (death) plus one or more
  non-fatal events: pass `illness.time` as an ordered vector of non-fatal event
  columns (highest priority first). A single column reproduces the previous
  two-tier behavior (backward compatible). It counts a **higher-priority event
  even when it follows a lower-priority one** (death after a hospitalization, a
  stroke after a hospitalization) -- the clinically intended hierarchy, and the
  one the first-event `getWinRatio()` cannot produce. This is now the recommended
  win ratio for most trials (e.g. device or oncology trials, or wherever a
  hierarchical composite endpoint is used).
* Internally it is a Markov multistate model whose states are the subsets of
  non-fatal events experienced; every transition (each non-fatal event out of each
  reachable state, and death out of every state) is a Super Learner, with
  doubly-robust, covariate-adjusted, censoring-corrected (IPCW) influence-function
  inference via adjoint-value clever covariates, and cross-fitting. The estimand
  and its inference are validated against a brute-force pairwise win ratio for
  hierarchies up to four time-to-event tiers (`scripts/genwr-*.R`).
* `getWinRatio()` documentation now flags its first-event limitation and points to
  `clinicalWinRatio()` for hierarchies where a higher-priority event can follow a
  lower-priority one.

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
* `clinicalWinRatio()` now cross-fits the transition and censoring hazards by
  default (`n.folds = 5`): nuisances are fit out-of-fold so each subject's
  influence-function contribution uses learners trained without them. This gives
  honest inference when the `SL.library` contains flexible learners that could
  over-fit; set `n.folds = 1` for faster in-sample fits with simple learners.
* Characterized and documented small-sample behavior. Like the win ratio in
  general (including the unadjusted Pocock win ratio), the estimate is a ratio and
  is mildly biased/anti-conservative at small n: in a null simulation it is biased
  downward ~1% at ~400/arm with coverage ~0.93-0.94 / type-I ~0.06-0.07, becoming
  nominal (0.95-0.97) by ~800/arm. This is a finite-sample property of the
  win-ratio functional, **not** an over-fitting artifact -- cross-fitting does not
  change it. Documented in the function help (`@section Small-sample behavior`),
  the vignette, and `scripts/make-clinical-wr-smalln.R`.
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
