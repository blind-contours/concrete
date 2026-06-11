# Trial-design and regulatory toolkit

This article collects the features of `concrete` that go beyond a single
covariate-adjusted absolute risk and that map onto how randomized trials
are designed, analyzed, and reviewed. The framing follows two regulatory
documents:

- **FDA (2023), *Adjusting for Covariates in Randomized Clinical Trials
  for Drugs and Biological Products*.** Covariate adjustment is
  encouraged because it can improve precision without changing the
  estimand in a randomized trial. The guidance favors collapsible
  marginal summaries (risk difference, risk ratio, restricted mean) over
  the hazard ratio, which is non-collapsible and hard to interpret under
  non-proportional hazards.
- **ICH E9(R1), the estimand and sensitivity-analysis addendum.** A
  trial analysis should name its estimand through five attributes and
  pre-specify how intercurrent events are handled, accompanied by
  sensitivity analyses for the assumptions that cannot be checked from
  the data.

Each section below shows the corresponding `concrete` function. The code
is illustrative; substitute your own trial `data.table`. Throughout,
`trial` is a one-row-per-subject table with an observed `time`, an
event-type `event` (`0` = censored, positive integers for events), a
binary `arm`, an `id`, and baseline covariates.

``` r

library(concrete)
library(data.table)

trial <- as.data.table(your_trial_data)
```

## 1. Name the estimand (ICH E9(R1))

[`makeEstimand()`](https://blind-contours.github.io/concrete/reference/makeEstimand.md)
records the analysis target using the five E9(R1) attributes. It does
not change the estimation — it documents the target so it can travel
with the results into a statistical analysis plan and the final report.

``` r

est_def <- makeEstimand(
  treatment  = "as randomized (intent-to-treat)",
  population  = "all randomized participants",
  endpoint    = "cause-specific absolute risk of relapse at 730 days",
  summary     = "risk difference",
  strategy    = "treatment policy"
)
est_def
```

## 2. Handle intercurrent events

An intercurrent event (ICE) — treatment discontinuation, rescue
medication, a competing terminal event — is coded as one of the
event-type values in the data.
[`applyIntercurrentEvent()`](https://blind-contours.github.io/concrete/reference/makeEstimand.md)
recodes the outcome to match the chosen E9(R1) strategy and returns a
copy of the data carrying the estimand as an attribute:

- **Treatment policy** (default ITT): the ICE is ignored, follow-up
  continues, data returned unchanged.
- **Hypothetical**: the ICE is recoded as censoring (`0`). `concrete`’s
  inverse-probability-of-censoring weighting then targets the risk that
  would be seen had the ICE not occurred, **under the strong, untestable
  assumption that the ICE acts as conditionally independent censoring**
  — pair it with the censoring sensitivity analysis in Section 6.
- **Composite**: the ICE is merged into the event of interest, defining
  a composite endpoint.

``` r

# Suppose event type 3 codes treatment discontinuation (the intercurrent event).
# Hypothetical strategy: estimate the risk had discontinuation not occurred.
trial_hyp <- applyIntercurrentEvent(
  trial, EventTime = "time", EventType = "event",
  Intercurrent = 3, strategy = "hypothetical", intercurrent = "discontinuation"
)

# Composite strategy: count discontinuation as an occurrence of event 1.
trial_comp <- applyIntercurrentEvent(
  trial, EventTime = "time", EventType = "event",
  Intercurrent = 3, strategy = "composite", TargetEvent = 1
)
```

The `"while on treatment"` and `"principal stratum"` strategies are
intentionally not supported: the former is a different (descriptive)
estimand that an IPCW estimator does not target, and the latter needs
latent stratification beyond the current scope.

### Treatment switching (crossover): ITT vs per-protocol

Treatment switching is the intercurrent event of device and oncology
trials: control participants cross over to the active arm (or vice
versa) at a time that is *not* fixed in advance, so it cannot be coded
as a single event type. Two estimands answer two different questions:

- **Intent-to-treat (treatment policy).** Analyze as randomized;
  switching is part of the strategy being estimated. This is the default
  — do nothing special.
- **Hypothetical (per-protocol-type).** Estimate the cumulative
  incidence that *would* have been seen had no one switched. `concrete`
  targets this when you pass a column of per-subject switch times to
  `Crossover`. Each switcher’s outcome is re-censored at their switch
  time, and a **separate crossover hazard** is fit and multiplied into
  the inverse-probability-of-censoring weight alongside the ordinary
  dropout hazard, so the weight denominator becomes
  `S_dropout(t) * S_crossover(t)`. This removes the selection bias that
  naive per-protocol censoring would introduce, because the crossover
  hazard reweights the participants who *would* have switched back into
  the at-risk set.

``` r

# `switch_time` holds each participant's crossover time (NA / Inf if they never
# switched). The crossover hazard inherits the censoring model's covariates.
args_pp <- formatArguments(
  DataTable = trial, EventTime = "time", EventType = "event", Treatment = "arm",
  ID = "id", Intervention = makeITT(), TargetTime = c(365, 730), TargetEvent = 1,
  Crossover = "switch_time"
)
est_pp <- doConcrete(args_pp)

# Report ITT and the hypothetical no-switching estimand side by side:
getOutput(doConcrete(args_itt), Estimand = "RD", TargetTime = 730)  # treatment policy
getOutput(est_pp,               Estimand = "RD", TargetTime = 730)  # no-switching
```

The no-switching estimand rests on the assumption that switching acts as
conditionally independent censoring given the covariates — the same
MAR-type assumption as ordinary censoring. Because it multiplies two
weights, the effective sample size can shrink in heavy-crossover trials;
check it with
[`getPositivityDx()`](https://blind-contours.github.io/concrete/reference/getPositivityDx.md)
(Section 8) and probe the assumption with the tipping-point analysis
(Section 6).

### Time-varying covariates in the censoring (and crossover) model

When dropout — or the decision to switch — is driven by
*post-randomization* measurements (echocardiography, KCCQ, six-minute
walk distance recorded at follow-up visits), conditioning the censoring
hazard on baseline covariates alone leaves it informative, biasing the
IPCW. Pass those measurements as `CensoringTV`, a long `data.frame` with
the `ID` column, a `time` column, and one value column per measurement.
`concrete` conditions the censoring hazard — and, when `Crossover` is
supplied, the crossover hazard, which inherits these covariates by
default — on the last-observation-carried-forward value and the change
from baseline of each, restoring conditional independence (CAR).

``` r

# Long-form post-randomization measurements: one row per (id, visit), with the
# ID column, a `time` column, and one column per measurement.
tv <- data.frame(id = visits$id, time = visits$visit_day,
                 kccq = visits$kccq, six_min_walk = visits$six_min_walk)

args_tv <- formatArguments(
  DataTable = trial, EventTime = "time", EventType = "event", Treatment = "arm",
  ID = "id", Intervention = makeITT(), TargetTime = c(365, 730), TargetEvent = 1,
  CensoringTV = tv,            # used ONLY in the censoring / crossover hazards
  Crossover   = "switch_time"  # optional; inherits the same covariates
)
est_tv <- doConcrete(args_tv)
```

These covariates enter **only** the censoring and crossover hazards,
never the outcome hazards: they are post-treatment mediators, so keeping
them out of the outcome model preserves the marginal / intent-to-treat
target while still correcting informative dropout. The correction flows
through every estimand that uses the IPCW — survival and
cumulative-incidence curves, RMST, and the win ratio — not just one of
them.

## 3. Cross-fitting (CV-TMLE) for machine-learning nuisances

When the propensity and hazard nuisances are flexible machine-learning
fits, the plug-in carries an empirical-process (overfitting) term that
can invalidate influence-function inference. Cross-fitting removes it by
predicting each subject’s nuisances from folds that did not include
them. Turn it on with a single argument:

``` r

args_cf <- formatArguments(
  DataTable = trial, EventTime = "time", EventType = "event", Treatment = "arm",
  ID = "id", Intervention = makeITT(), TargetTime = c(365, 730), TargetEvent = 1,
  CVArg = list(V = 5),
  CrossFit = TRUE
)
est_cf <- doConcrete(args_cf)
```

Cross-fitting is the continuous-time survival analogue of CV-TMLE /
double machine learning, and is what licenses the use of the
random-forest, HAL, and ensemble hazard learners with honest standard
errors.

## 4. Ensemble hazard Super Learner

By default `concrete` performs discrete Super Learner: cross-validation
*selects* the single best hazard learner. `HazEnsemble = TRUE` instead
fits a convex *combination* of the library, with weights chosen to
minimize the cross-validated counting-process loss. This is more robust
when no single learner dominates across the time axis.

``` r

args_ens <- formatArguments(
  DataTable = trial, EventTime = "time", EventType = "event", Treatment = "arm",
  ID = "id", Intervention = makeITT(), TargetTime = c(365, 730), TargetEvent = 1,
  CVArg = list(V = 5),
  CrossFit = TRUE, HazEnsemble = TRUE,
  Model = list(
    arm = c("SL.glm", "SL.glmnet"),
    "0" = list(Cox = survival::Surv(time, event == 0) ~ ., RSF = "rsf"),
    "1" = list(Cox = survival::Surv(time, event == 1) ~ ., HAL = "hal")
  )
)
est_ens <- doConcrete(args_ens)
```

## 5. Win ratio, win odds, and net benefit

For hierarchical or composite comparisons,
[`getWinRatio()`](https://blind-contours.github.io/concrete/reference/getWinRatio.md)
turns the targeted survival curves into a covariate-adjusted,
doubly-robust win ratio, win odds, and net benefit, each with an
influence-function confidence interval. Because it is built from the
targeted curves rather than from raw pairwise comparisons, it inherits
`concrete`’s covariate adjustment and censoring handling. As a plug-in
it evaluates the win integral over the `TargetTime` grid, so use a
reasonably dense grid up to the horizon — or use
[`targetWinRatio()`](https://blind-contours.github.io/concrete/reference/targetWinRatio.md)
below, which removes that grid sensitivity entirely.

``` r

args_wr <- formatArguments(
  DataTable = trial, EventTime = "time", EventType = "event", Treatment = "arm",
  ID = "id", Intervention = makeITT(),
  TargetTime = c(180, 365, 545, 730), TargetEvent = 1, CVArg = list(V = 5)
)
est_wr <- doConcrete(args_wr)

getWinRatio(est_wr, Horizon = 730, Intervention = c(1, 2))
#> reports P(win), P(loss), P(tie), Win Ratio, Win Odds, and Net Benefit
```

A win ratio above 1 (and a net benefit above 0) favors the active arm;
the confidence interval for the win ratio crossing 1 is the
corresponding test of no difference.

For a **prioritized hierarchy** of competing events, pass `TargetEvent`
as an ordered vector (highest priority first). For example, with cause 1
= death and cause 2 = hospitalization, `TargetEvent = c(1, 2)` decides
each pairwise comparison on death first and breaks ties on
hospitalization (the Pocock / Finkelstein–Schoenfeld rule). The events
must have been targeted in
[`doConcrete()`](https://blind-contours.github.io/concrete/reference/doConcrete.md),
and the comparison uses each patient’s first event (the competing- risks
structure `concrete` models).

``` r

args_h <- formatArguments(
  DataTable = trial, EventTime = "time", EventType = "event", Treatment = "arm",
  ID = "id", Intervention = makeITT(),
  TargetTime = c(180, 365, 545, 730), TargetEvent = c(1, 2), CVArg = list(V = 5)
)
getWinRatio(doConcrete(args_h), Horizon = 730, Intervention = c(1, 2),
            TargetEvent = c(1, 2))   # death > hospitalization
```

For the primary analysis, prefer the **directly targeted** version:
[`targetWinRatio()`](https://blind-contours.github.io/concrete/reference/targetWinRatio.md)
fluctuates both arms’ hazards over the full event-time grid until the
win and loss probabilities’ own estimating equations are solved, rather
than plugging the pointwise-targeted curves into the win functional. In
validation it cut the residual win-ratio bias about five-fold and
restored nominal coverage on sparse target grids (see the [win
ratio](https://blind-contours.github.io/concrete/articles/win-ratio.md)
article for the numbers):

``` r

targetWinRatio(doConcrete(args_h), Horizon = 730, Intervention = c(1, 2),
               TargetEvent = c(1, 2))
#> same six rows as getWinRatio(); attr(., "WRConverged") reports convergence
```

## 6. Censoring sensitivity (tipping point)

The hypothetical strategy and the core estimator both assume censoring
is independent given the measured covariates (MAR).
[`senseCensoring()`](https://blind-contours.github.io/concrete/reference/senseCensoring.md)
stress-tests that assumption: it imputes an increasing fraction `delta`
of censored subjects as having had the event of interest, re-fits the
estimator at each `delta`, and reports the **tipping point** — the
smallest fraction at which a significant conclusion would be overturned.

``` r

sc <- senseCensoring(
  args_wr,
  deltas = c(0, 0.05, 0.10, 0.15, 0.20),
  Estimand = "RD", Intervention = c(1, 2)
)
sc                          # estimate / CI / p-value at each delta
attr(sc, "tippingPoint")    # smallest delta that overturns the conclusion
```

A conclusion that survives a large `delta` is robust to departures from
independent censoring; one that tips at a small `delta` should be
reported with caution.

When the fit carries a crossover model (Section 2), the censored
subjects pool two intercurrent events with two *different* untestable
assumptions — ordinary dropout (MAR censoring) and switching (the
no-switching counterfactual). The `mechanism` argument lets you tip each
pool individually or jointly, so each assumption gets its own tipping
point:

``` r

sc_drop  <- senseCensoring(args_pp, deltas = c(0, 0.1, 0.2),
                           Estimand = "RD", mechanism = "dropout")    # MAR only
sc_xover <- senseCensoring(args_pp, deltas = c(0, 0.1, 0.2),
                           Estimand = "RD", mechanism = "crossover")  # no-switching only
sc_both  <- senseCensoring(args_pp, deltas = c(0, 0.1, 0.2),
                           Estimand = "RD", mechanism = "all")        # jointly (default)
rbind(sc_drop, sc_xover, sc_both)   # leading `mechanism` column distinguishes them
```

For a per-protocol (no-switching) analysis, `mechanism = "crossover"`
answers the reviewer’s question directly: *how sensitive is the
conclusion to the assumption that switchers would not have had the event
had they stayed on their assigned arm?*

## 7. How much did covariate adjustment buy you?

Covariate adjustment does not change the estimand in a randomized trial,
but it can tighten the confidence intervals. Fit an unadjusted
(treatment-only) version and compare with
[`getRelativeEfficiency()`](https://blind-contours.github.io/concrete/reference/getRelativeEfficiency.md),
which reports the variance ratio.

``` r

unadj_args <- formatArguments(
  DataTable = trial, EventTime = "time", EventType = "event", Treatment = "arm",
  ID = "id", Intervention = makeITT(), TargetTime = c(365, 730), TargetEvent = 1,
  CVArg = list(V = 5),
  Model = list(arm = "SL.mean",
               "0" = list(Cox = survival::Surv(time, event == 0) ~ arm),
               "1" = list(Cox = survival::Surv(time, event == 1) ~ arm))
)
unadj_est <- doConcrete(unadj_args)

getRelativeEfficiency(
  Adjusted   = getOutput(est_cf,    Estimand = "RD", Intervention = c(1, 2)),
  Unadjusted = getOutput(unadj_est, Estimand = "RD", Intervention = c(1, 2))
)
```

A `RelEfficiency` above 1 means adjustment was worth it: a value of 1.25
says the adjusted analysis has the precision of an unadjusted analysis
on 25% more participants.

## 8. Positivity and effective sample size

Every IPCW-based estimand depends on participants having a non-trivial
probability of remaining uncensored (and, under the hypothetical
no-switching estimand, un-switched) through the target time. When that
probability is small for some participants, their inverse weights blow
up, the effective sample size (ESS) drops, and the variance inflates —
the price of the hypothetical estimand.
[`getPositivityDx()`](https://blind-contours.github.io/concrete/reference/getPositivityDx.md)
reports these diagnostics per arm alongside the estimates so the
trialist can see when an extrapolation is fragile.

``` r

getPositivityDx(est_pp)
#> per-arm ESS (overall and worst-time), max weight, minimum observation
#> probability, and the % of weights pinned at the truncation bound, with a
#> CAUTION flag when ESS is low or truncation is heavy.
```

A low ESS or a large fraction of truncated weights is the signal that a
heavy-crossover (or heavy-dropout) trial cannot support the no-switching
estimand without strong extrapolation; report it, and lean on the
tipping-point analysis (Section 6) to bound the conclusion.

## 9. Variance under stratified randomization

Nearly every phase-3 trial randomizes within strata (site, disease
severity, biomarker) using permuted blocks or a stratified biased coin,
and ICH E9 and the FDA covariate-adjustment guidance ask the analysis to
reflect that design. The usual influence-function standard errors assume
*simple* randomization; under covariate-adaptive schemes they are
generically **conservative** — they ignore the
between-arm-within-stratum variance the design removes — which gives
away exactly the precision covariate adjustment is meant to buy.

Pass the randomization strata to
[`formatArguments()`](https://blind-contours.github.io/concrete/reference/formatArguments.md)
and every reported standard error — absolute risk, risk difference and
ratio, RMST and life-years lost, and the win ratio — is corrected
following Bugni–Canay–Shaikh / Ye–Shao:

``` r

args_strat <- formatArguments(
  DataTable = trial, EventTime = "time", EventType = "event", Treatment = "arm",
  ID = "id", Intervention = makeITT(), TargetTime = c(365, 730), TargetEvent = 1,
  Strata = c("site", "severity")   # the variables randomization was stratified on
)
est_strat <- doConcrete(args_strat)
getOutput(est_strat, Estimand = "RD")   # SEs reflect the stratified design
```

Three things to know:

- The strata columns stay in the data as adjustment covariates
  (recommended). When the working models adjust for them well, the
  correction is approximately zero — the iid variance is then already
  correct. The correction matters exactly when adjustment for the
  stratification variables is absent or imperfect.
- Only supply `Strata` when randomization truly was stratified: applying
  the correction under simple randomization understates the variance.
- Each stratum needs both arms represented (at least 2 subjects per
  arm); if not, `concrete` warns and reports the conservative iid
  standard errors instead. Pool very small strata before analysis.

Missing baseline covariates, incidentally, no longer stop the pipeline:
NA values in baseline covariates are imputed (median / mode) with a
missingness indicator added per affected column — the handling endorsed
for pre-randomization covariates by the FDA guidance. The outcome,
treatment, and ID columns must still be complete.

## Putting it together: an analysis-plan checklist

A reproducible, E9(R1)-aware analysis with `concrete` typically records:

1.  The estimand, via
    [`makeEstimand()`](https://blind-contours.github.io/concrete/reference/makeEstimand.md),
    including the intercurrent-event strategy — and, for treatment
    switching, whether it is intent-to-treat or the hypothetical
    no-switching estimand (`Crossover`).
2.  The data handling for that strategy, via
    [`applyIntercurrentEvent()`](https://blind-contours.github.io/concrete/reference/makeEstimand.md)
    for event-coded ICEs, or `Crossover` / `CensoringTV` for switching
    and informative dropout driven by post-randomization measurements.
3.  A cross-fitted (`CrossFit = TRUE`), optionally ensemble
    (`HazEnsemble = TRUE`) fit, so machine-learning nuisances keep
    honest inference.
4.  The primary summary — risk difference / ratio
    ([`getOutput()`](https://blind-contours.github.io/concrete/reference/getOutput.md)),
    restricted mean survival time
    ([`getRMST()`](https://blind-contours.github.io/concrete/reference/getRMST.md)
    /
    [`targetRMST()`](https://blind-contours.github.io/concrete/reference/targetRMST.md)),
    or win ratio
    ([`getWinRatio()`](https://blind-contours.github.io/concrete/reference/getWinRatio.md)
    /
    [`targetWinRatio()`](https://blind-contours.github.io/concrete/reference/targetWinRatio.md),
    preferring the directly targeted version).
5.  Pre-specified sensitivity analyses —
    [`senseCensoring()`](https://blind-contours.github.io/concrete/reference/senseCensoring.md)
    for the independent-censoring assumption — and positivity
    diagnostics via
    [`getPositivityDx()`](https://blind-contours.github.io/concrete/reference/getPositivityDx.md).
6.  The randomization design reflected in the variance: pass `Strata`
    when randomization was stratified (permuted blocks, biased coin,
    minimization).
7.  The precision gain from adjustment, via
    [`getRelativeEfficiency()`](https://blind-contours.github.io/concrete/reference/getRelativeEfficiency.md).

See the [Trialist
quickstart](https://blind-contours.github.io/concrete/articles/trialist-quickstart.md)
for the core workflow and [How concrete
works](https://blind-contours.github.io/concrete/articles/how-concrete-works.md)
for the estimator itself.
