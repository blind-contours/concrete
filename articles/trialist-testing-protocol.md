# Testing protocol and current limitations

This article gives trial analysts a short protocol for testing
`concrete` on their own data and reporting useful feedback.

## Step 1: Install and run the package smoke test

Install the development version from GitHub:

``` r

install.packages("remotes")
remotes::install_github("blind-contours/concrete")
```

Then run the built-in smoke test:

``` r

library(concrete)

source(system.file("examples", "trialist-smoke-test.R", package = "concrete"))
```

The smoke test runs a small competing-risk analysis using a Cox-only
hazard library, prints event counts, prints absolute risks, risk
differences, and risk ratios, and returns a `smoke_summary` table.

The first table should look like this:

| arm | event |   N |
|----:|------:|----:|
|   0 |     0 |  38 |
|   0 |     1 |  32 |
|   0 |     2 |   3 |
|   1 |     0 |  49 |
|   1 |     1 |  31 |
|   1 |     2 |   7 |

The smoke-test summary should show an `ok` status, convergence, and no
failing components:

| analysis | status | elapsed_sec | converged | step | max_ratio | failing_components |
|----------|--------|------------:|-----------|-----:|----------:|-------------------:|
| cox_only | ok     |         1.4 | TRUE      |    4 |     0.743 |                  0 |

To also try optional hazard learners that are installed on your machine:

``` r

Sys.setenv(CONCRETE_RUN_OPTIONAL_LEARNERS = "true")
source(system.file("examples", "trialist-smoke-test.R", package = "concrete"))
```

The optional pass attempts additive hazards, Coxnet, random survival
forests, and HAL when the required packages are installed.

## Step 2: Run the conservative analysis on your trial

Use a simple learner library first. This makes it easier to identify
whether any problem is due to the data structure, the estimand, or a
flexible learner.

``` r

Model <- list(
  arm = c("SL.mean", "SL.glm"),
  "0" = list(Censor = survival::Surv(time, event == 0) ~ arm + age + sex),
  "1" = list(Event = survival::Surv(time, event == 1) ~ arm + age + sex)
)

ConcreteArgs <- formatArguments(
  DataTable = trial,
  EventTime = "time",
  EventType = "event",
  Treatment = "arm",
  ID = "id",
  Intervention = makeITT(),
  TargetTime = c(365, 730),
  TargetEvent = 1,
  CVArg = list(V = 5),
  Model = Model,
  UpdateMethod = "adaptive",
  EICStopRule = "absolute",
  EICStopAbsTol = 0.02 / sqrt(nrow(trial)),
  Verbose = FALSE
)

ConcreteEst <- doConcrete(ConcreteArgs)
ConcreteOut <- getOutput(
  ConcreteEst,
  Estimand = c("Risk", "RD", "RR"),
  Intervention = c(1, 2)
)

ConcreteOut
getTmleDiagnostics(ConcreteEst, type = "components")
```

For a successful first trial run, expect three linked outputs:

1.  Event counts by arm showing enough events near the target time.
2.  A
    [`getOutput()`](https://blind-contours.github.io/concrete/reference/getOutput.md)
    table with absolute risks, risk differences, and risk ratios.
3.  A diagnostics table with `check = TRUE` for all targeted components.

Example trial-output rows:

| Time | Event | Estimand  | Intervention      | Estimator | Pt Est |   se |
|-----:|------:|-----------|-------------------|-----------|-------:|-----:|
|  730 |     1 | Abs Risk  | A=0               | tmle      |   0.12 | 0.03 |
|  730 |     1 | Abs Risk  | A=1               | tmle      |   0.10 | 0.03 |
|  730 |     1 | Risk Diff | \[A=1\] - \[A=0\] | tmle      |  -0.02 | 0.04 |
|  730 |     1 | Rel Risk  | \[A=1\] / \[A=0\] | tmle      |   0.83 | 0.24 |

For competing risks, add one model entry for each positive event code.

## Step 3: Compare against your usual analysis

Record the standard analysis that your trial team would usually report.

Useful comparisons:

- event and censoring counts by randomized arm
- unadjusted cumulative incidence or Kaplan-Meier estimates at the same
  target times
- cause-specific Cox model output, when relevant
- adjusted `concrete` absolute risks
- adjusted `concrete` risk differences and risk ratios
- g-computation plug-in estimates from `getOutput(..., GComp = TRUE)`

The Cox hazard ratio is not the same estimand as the marginal risk ratio
from `concrete`, so use it as context rather than as a direct equality
check.

## Step 4: Escalate learners

Use the same data and estimand while changing only the learner library.

``` r

model_cox <- list(
  arm = c("SL.mean", "SL.glm"),
  "0" = list(Cox = survival::Surv(time, event == 0) ~ .),
  "1" = list(Cox = survival::Surv(time, event == 1) ~ .)
)

model_coxnet <- list(
  arm = c("SL.mean", "SL.glm", "SL.glmnet"),
  "0" = list(Cox = survival::Surv(time, event == 0) ~ ., Coxnet = "coxnet"),
  "1" = list(Cox = survival::Surv(time, event == 1) ~ ., Coxnet = "coxnet")
)

model_flexible <- list(
  arm = c("SL.mean", "SL.glm", "SL.glmnet"),
  "0" = list(Cox = survival::Surv(time, event == 0) ~ ., Aalen = "aareg"),
  "1" = list(
    Cox = survival::Surv(time, event == 1) ~ .,
    Coxnet = "coxnet",
    RSF = "rsf",
    Aalen = "aareg",
    HAL = "hal"
  )
)
```

Compare point estimates, runtime, selected learners, and convergence
diagnostics.

## Step 5: Exercise the trial-design features

These are the newest parts of the package, so feedback here is the most
valuable. Try whichever apply to your trial:

- **Missing baseline covariates.** Run your data as-is — NA baseline
  covariates are imputed (median / mode) with a `<column>_missing`
  indicator added. Check the message lists what you expect, and tell us
  if the imputation behavior surprises you.
- **Stratified randomization.** If your trial randomized within strata
  (permuted blocks, biased coin), pass `Strata = c(...)` and compare the
  standard errors with and without it. The corrected SEs should be the
  same or tighter; the point estimates must not change.
- **Treatment switching.** If participants crossed over, pass a
  switch-time column as `Crossover` and compare the hypothetical
  no-switching estimand against your ITT run. Always check
  [`getPositivityDx()`](https://blind-contours.github.io/concrete/reference/getPositivityDx.md)
  afterwards — heavy switching shrinks the effective sample size.
- **Informative dropout.** If post-randomization measurements (labs, QoL
  scores, functional tests) drive dropout in your trial, supply them as
  `CensoringTV` and see how much the estimates move.
- **Sensitivity.** Run
  [`senseCensoring()`](https://blind-contours.github.io/concrete/reference/senseCensoring.md)
  — with `mechanism = "dropout"` / `"crossover"` separately if you used
  `Crossover` — and report whether the tipping point lands where your
  clinical intuition says it should.

## Current limitations

The current public testing target is intentionally narrow.

Supported:

- one row per participant
- right-censored event time outcome
- optional competing risks
- baseline binary treatment coded `0` and `1`
- static interventions such as everyone assigned `A = 1` versus everyone
  assigned `A = 0`
- baseline covariate adjustment, with missing baseline values imputed
  automatically (missingness indicators added)
- stratified / covariate-adaptive randomization via `Strata` (corrected
  standard errors)
- treatment switching via `Crossover` (hypothetical no-switching
  estimand)
- post-baseline time-varying covariates in the **censoring** model via
  `CensoringTV` (informative dropout)
- target absolute risks, risk differences, and risk ratios at
  prespecified times, plus RMST / life-years lost and the win ratio
  family

Not currently supported in the main trialist workflow:

- longitudinal treatment regimes
- recurrent events
- delayed entry or left truncation
- multi-arm or continuous treatment without custom intervention work
- post-baseline time-varying covariates in the **outcome** model (they
  are post-treatment mediators; only the censoring model uses them)
- clustered trial designs requiring special variance handling

Use caution when:

- one arm has very few events by the target time
- censoring is highly imbalanced by arm or covariates
- optional machine-learning learners dominate a small trial
- the flexible learner results differ sharply from Cox-only results
- the TMLE update has large absolute `PnEIC` values after the adaptive
  update

## What to send back

For useful feedback, include:

- package version from `packageVersion("concrete")`
- [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html)
- event and censoring counts by treatment arm
- target event and target time
- exact `Model` list
- `UpdateMethod`, `EICStopRule`, and `EICStopAbsTol`
- [`getOutput()`](https://blind-contours.github.io/concrete/reference/getOutput.md)
  table
- `getTmleDiagnostics(ConcreteEst, type = "components")`
- whether the issue also occurs with a Cox-only learner library

GitHub issue templates are available for convergence issues, learner
failures, and estimand questions.
