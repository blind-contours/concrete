# Trialist quickstart

This article is written for trial analysts who want to test `concrete`
on a randomized trial or trial-like observational data set. The core
output is a covariate-adjusted marginal absolute risk at prespecified
follow-up times. From those risks, `concrete` reports risk differences
and risk ratios.

## What question is being answered?

For a binary baseline treatment, `concrete` estimates counterfactual
absolute risks such as:

- risk of event 1 by day 365 if everyone followed the active arm
- risk of event 1 by day 365 if everyone followed the control arm
- active minus control risk difference
- active divided by control risk ratio

For competing risks, the estimate is the cumulative incidence of the
target event in the presence of the competing events. For a standard
right-censored survival endpoint with no competing risks, use one
positive event code and code censoring as `0`.

## Data checklist

Your analysis data set should have one row per participant for the
current package interface.

| Item       | Requirement                                         |
|------------|-----------------------------------------------------|
| Subject id | Unique participant id column                        |
| Time       | Observed event or censoring time, numeric, positive |
| Event type | `0` for censoring, positive integers for events     |
| Treatment  | Binary baseline treatment coded `0` and `1`         |
| Covariates | Baseline covariates measured before treatment       |

Check these before fitting:

``` r

trial[, .N, by = event]
trial[, .N, by = .(arm, event)]
summary(trial$time)
stopifnot(all(trial$event >= 0))
stopifnot(all(trial$arm %in% c(0, 1)))
```

Current scope:

- baseline binary treatment
- one row per participant
- right-censored event time data
- optional competing risks
- baseline covariate adjustment

Currently out of scope:

- longitudinal treatment regimes
- recurrent events
- left truncation or delayed entry
- non-binary treatment without custom intervention work

## Example data

The examples below use
[`survival::pbc`](https://rdrr.io/pkg/survival/man/pbc.html) only to
provide a complete data shape. Replace this block with your trial data.

``` r

library(concrete)
library(data.table)

trial <- as.data.table(survival::pbc)
trial <- trial[!is.na(trt), .(id, time, status, trt, age, sex, albumin, bili)]

# Example treatment coding: 1 = active arm, 0 = control arm.
trial[, arm := as.integer(trt == 2)]

# Example competing-risk coding:
# 0 = censored, 1 = death, 2 = transplant.
trial[, event := data.table::fifelse(
  status == 2, 1L,
  data.table::fifelse(status == 1, 2L, 0L)
)]

trial <- trial[, .(id, time, event, arm, age, sex, albumin, bili)]
```

Before fitting your own data, you can run the installed smoke test:

``` r

source(system.file("examples", "trialist-smoke-test.R", package = "concrete"))
```

This verifies that the package can run a small competing-risk analysis
and prints the same kinds of outputs and diagnostics you should inspect
in your trial.

## First intent-to-treat analysis

Start with a conservative analysis before adding flexible hazard
learners. This uses the default treatment Super Learner and default Cox
candidate hazards.

``` r

ConcreteArgs <- formatArguments(
  DataTable = trial,
  EventTime = "time",
  EventType = "event",
  Treatment = "arm",
  ID = "id",
  Intervention = makeITT(),
  TargetTime = c(365, 730, 1095),
  TargetEvent = 1,
  CVArg = list(V = 5),
  Verbose = FALSE
)

ConcreteEst <- doConcrete(ConcreteArgs)
print(ConcreteEst)
```

[`makeITT()`](https://blind-contours.github.io/concrete/reference/formatArguments.md)
creates two static interventions:

- `A=1`: everyone assigned to the active arm
- `A=0`: everyone assigned to the control arm

This is the usual starting point for an intent-to-treat comparison when
the treatment column is the randomized treatment assignment.

## Output for a trial report

``` r

ConcreteOut <- getOutput(
  ConcreteEst,
  Estimand = c("Risk", "RD", "RR"),
  Intervention = c(1, 2),
  GComp = TRUE,
  Simultaneous = TRUE
)

ConcreteOut[Event == 1]
```

Example output will look like this. The exact values will change with
your data, target times, learner library, and confidence-interval
settings.

| Time | Event | Estimand  | Intervention      | Estimator | Pt Est |    se | CI Low | CI Hi |
|-----:|------:|-----------|-------------------|-----------|-------:|------:|-------:|------:|
| 1000 |     1 | Abs Risk  | A=0               | tmle      |   0.19 | 0.036 |   0.11 |  0.26 |
| 1000 |     1 | Abs Risk  | A=1               | tmle      |   0.21 | 0.042 |   0.12 |  0.29 |
| 1000 |     1 | Risk Diff | \[A=1\] - \[A=0\] | tmle      |  0.021 | 0.052 | -0.082 |  0.12 |
| 1000 |     1 | Rel Risk  | \[A=1\] / \[A=0\] | tmle      |   1.10 |  0.30 |   0.53 |  1.70 |
| 2000 |     1 | Abs Risk  | A=0               | tmle      |   0.34 | 0.053 |   0.24 |  0.44 |
| 2000 |     1 | Abs Risk  | A=1               | tmle      |   0.33 | 0.047 |   0.24 |  0.43 |
| 2000 |     1 | Risk Diff | \[A=1\] - \[A=0\] | tmle      | -0.007 | 0.065 |  -0.14 |  0.12 |
| 2000 |     1 | Rel Risk  | \[A=1\] / \[A=0\] | tmle      |   0.98 |  0.19 |   0.61 |  1.40 |

Important columns:

- `Time`: target follow-up time
- `Event`: target event code
- `Estimand`: absolute risk, risk difference, or relative risk
- `Intervention`: intervention being estimated or compared
- `Estimator`: `tmle` or optional `gcomp` plug-in estimate
- `Pt Est`: point estimate
- `CI Low`, `CI Hi`: confidence interval

For risk differences, positive values mean the first intervention listed
in `Intervention = c(1, 2)` has higher risk than the second. With
[`makeITT()`](https://blind-contours.github.io/concrete/reference/formatArguments.md),
that is `A=1` minus `A=0`.

## Compare with familiar analyses

Use standard analyses as context, not as exact substitutes. A Cox hazard
ratio is not the same estimand as a marginal risk difference or
cumulative incidence risk ratio.

``` r

# Event counts by arm are always the first check.
trial[, .N, by = .(arm, event)]

# Cause-specific Cox model for the event of interest.
cox_event1 <- survival::coxph(
  survival::Surv(time, event == 1) ~ arm + age + sex + albumin + bili,
  data = trial
)
summary(cox_event1)

# Kaplan-Meier style context for a non-competing-risk endpoint.
# For competing risks, prefer a cumulative incidence estimator in your usual
# trial reporting workflow and compare target times to ConcreteOut.
km_event1 <- survival::survfit(
  survival::Surv(time, event == 1) ~ arm,
  data = trial
)
summary(km_event1, times = c(365, 730, 1095))
```

Useful comparison questions:

- Are event counts and censoring patterns plausible by treatment arm?
- Are `concrete` absolute risks close to unadjusted estimates in a
  balanced trial?
- Do adjusted estimates move in the expected direction when important
  prognostic covariates are included?
- Are the TMLE estimates and g-computation plug-in estimates close?
- Are confidence intervals wider in rare-event settings?

## Check convergence

``` r

components <- getTmleDiagnostics(ConcreteEst, type = "components")
components[order(ratio, decreasing = TRUE)][1:10]

trace <- getTmleDiagnostics(ConcreteEst, type = "trace")
trace
```

For a clean fit, the component diagnostics should have `check = TRUE`
for all targeted components under the selected stopping rule. A compact
success summary from the smoke test looks like:

| analysis | status | converged | step | max_ratio | failing_components |
|----------|--------|-----------|-----:|----------:|-------------------:|
| cox_only | ok     | TRUE      |    4 |     0.743 |                  0 |

If the fit does not converge, try the settings in the convergence
diagnostics article before changing the estimand.

``` r

ConcreteArgs$UpdateMethod <- "adaptive"
ConcreteArgs$EICStopRule <- "hybrid"
ConcreteArgs$EICStopAbsTol <- 1e-3
ConcreteArgs <- formatArguments(ConcreteArgs)

ConcreteEst <- doConcrete(ConcreteArgs)
getTmleDiagnostics(ConcreteEst, type = "components")
```

## What to share when testing

When sharing results with collaborators or opening a GitHub issue,
include:

- package version from `packageVersion("concrete")`
- event counts by treatment arm
- target event and target time
- learner library used in `Model`
- TMLE controls: `UpdateMethod`, `EICStopRule`, `EICStopAbsTol`
- `getTmleDiagnostics(ConcreteEst, type = "components")`
- whether the issue reproduces with the conservative Cox-only analysis
