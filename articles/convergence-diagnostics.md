# Convergence diagnostics

This article explains how to evaluate the TMLE update when `concrete`
reports slow convergence or non-convergence. The goal is to distinguish
a meaningful targeting problem from a numerically tiny empirical
efficient influence curve (EIC) component in a rare-event setting.

## What convergence means here

The TMLE update tries to make the empirical mean of each requested
target EIC component small. Stopping is evaluated for the requested
intervention, target event, and target time components; internal
complement rows used to form survival contrasts are not used as
additional stopping equations. The original stopping rule checks:

``` text
abs(PnEIC) <= seEIC / (sqrt(n) * log(n))
```

This is a relative, component-specific threshold. It can become
extremely strict when an event is rare or a component has near-zero
variability.

`concrete` now exposes three stopping rules:

| Rule | Meaning | Typical use |
|----|----|----|
| `relative` | Original component-specific rule | Default, best first check |
| `absolute` | `abs(PnEIC) <= EICStopAbsTol` | Sensitivity analyses on absolute-risk scale |
| `hybrid` | `abs(PnEIC) <= max(relative threshold, EICStopAbsTol)` | Rare events or near-zero EIC variance |

A setting such as `EICStopRule = "hybrid"` and `EICStopAbsTol = 1e-3`
means: keep the original relative criterion unless it is stricter than
an absolute empirical EIC tolerance of 0.001 on the risk scale.

## Inspect the diagnostics

``` r

components <- getTmleDiagnostics(ConcreteEst, type = "components")
components[order(ratio, decreasing = TRUE)]

trace <- getTmleDiagnostics(ConcreteEst, type = "trace")
trace

norm <- getTmleDiagnostics(ConcreteEst, type = "norm")
norm
```

Example component output for a converged fit:

| Intervention | Time | Event | PnEIC | RelativeCriteria | AbsoluteCriteria | StopCriteria | ratio | check |
|----|---:|---:|---:|---:|---:|---:|---:|----|
| A=1 | 1000 | 1 | 0.00057 | 0.00077 | 0.001 | 0.001 | 0.57 | TRUE |
| A=0 | 1000 | 1 | -0.00074 | 0.00089 | 0.001 | 0.001 | 0.74 | TRUE |
| A=1 | 2000 | 1 | 0.00022 | 0.00110 | 0.001 | 0.00110 | 0.20 | TRUE |
| A=0 | 2000 | 1 | -0.00031 | 0.00104 | 0.001 | 0.00104 | 0.30 | TRUE |

Example trace output:

| Step | NormPnEIC | MaxRatio | FailingComponents | MaxAbsPnEIC |
|-----:|----------:|---------:|------------------:|------------:|
|    0 |    0.0204 |     8.12 |                 4 |      0.0073 |
|    1 |    0.0068 |     2.45 |                 2 |      0.0022 |
|    2 |    0.0021 |     1.08 |                 1 |      0.0011 |
|    3 |    0.0013 |     0.74 |                 0 |      0.0007 |

The component table is usually the most useful first view.

Key columns:

- `PnEIC`: empirical mean EIC for the component
- `RelativeCriteria`: original relative stopping threshold
- `AbsoluteCriteria`: absolute threshold supplied by `EICStopAbsTol`
- `StopCriteria`: threshold used by the selected rule
- `ratio`: `abs(PnEIC) / StopCriteria`
- `check`: whether that component passed the selected rule
- `Converged`: whether the overall update converged
- `ConvergenceStep`: update step where convergence was reached

Focus first on rows with `check == FALSE`, sorted by `ratio`.

``` r

components[check == FALSE][order(ratio, decreasing = TRUE)]
```

If there are failures, the output will identify which
event/time/intervention combination is driving the problem:

| Intervention | Time | Event | AbsPnEIC | StopCriteria | ratio | check |
|--------------|-----:|------:|---------:|-------------:|------:|-------|
| A=1          |  730 |     1 |   0.0048 |       0.0010 |   4.8 | FALSE |
| A=0          |  730 |     1 |   0.0017 |       0.0010 |   1.7 | FALSE |

This pattern says the empirical EIC is still materially larger than the
chosen threshold. Start with adaptive updating and a simpler learner
library before loosening the stopping rule.

## Recommended escalation ladder

Start simple and add flexibility only after the conservative analysis
behaves as expected.

### 1. Run a conservative baseline

Use default Cox hazards and a small treatment Super Learner library.

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
  Verbose = FALSE
)

ConcreteEst <- doConcrete(ConcreteArgs)
```

### 2. Use adaptive updating

The adaptive method uses a line search with rollback and is the
recommended first convergence fix. With `EICStopRule = "relative"` it
accepts updates that reduce the target empirical EIC norm. With
`EICStopRule = "absolute"` or `"hybrid"` it accepts updates that reduce
the active component-wise stopping ratio.

``` r

ConcreteArgs$UpdateMethod <- "adaptive"
ConcreteArgs <- formatArguments(ConcreteArgs)
ConcreteEst <- doConcrete(ConcreteArgs)
```

### 3. Use the hybrid stopping rule for rare events

``` r

ConcreteArgs$UpdateMethod <- "adaptive"
ConcreteArgs$EICStopRule <- "hybrid"
ConcreteArgs$EICStopAbsTol <- 1e-3
ConcreteArgs <- formatArguments(ConcreteArgs)
ConcreteEst <- doConcrete(ConcreteArgs)
```

Use this when the largest failing components have very small absolute
`PnEIC` values but large ratios because the relative threshold is tiny.

### 4. Increase iterations only if progress is continuing

``` r

ConcreteArgs$MaxUpdateIter <- 1000
ConcreteArgs <- formatArguments(ConcreteArgs)
ConcreteEst <- doConcrete(ConcreteArgs)
```

If the trace has flattened and the same tiny components remain,
increasing iterations may not change the practical estimate.

### 5. Simplify or stabilize nuisance estimation

If `abs(PnEIC)` remains large, the issue may be nuisance instability
rather than only the stopping threshold.

Try:

- simpler hazard learner libraries
- fewer high-variance learners for small trials
- stronger propensity-score truncation through `MinNuisance`
- fewer target times for initial debugging
- checking for arms with no events near the target time

## Interpreting common patterns

| Pattern | Likely meaning | Next step |
|----|----|----|
| Large `ratio`, tiny `AbsPnEIC` | Relative rule is too strict on a near-zero component | Try `hybrid` with `1e-3` |
| Large `ratio`, large `AbsPnEIC` | Targeting problem remains meaningful | Use adaptive update and inspect learners |
| Many failing components for censoring | Censoring or positivity instability | Check censoring by arm and covariates |
| Failure only for one rare event/time | Sparse target component | Report event counts and try hybrid rule |
| Norm decreases then rebounds | Update overshooting | Use `UpdateMethod = "adaptive"` |

## A reporting template

For trial reports or issue reports, record:

``` r

list(
  package_version = as.character(packageVersion("concrete")),
  update_method = ConcreteArgs$UpdateMethod,
  eic_stop_rule = ConcreteArgs$EICStopRule,
  eic_stop_abs_tol = ConcreteArgs$EICStopAbsTol,
  max_update_iter = ConcreteArgs$MaxUpdateIter,
  target_time = ConcreteArgs$TargetTime,
  target_event = ConcreteArgs$TargetEvent,
  event_counts = trial[, .N, by = .(arm, event)],
  components = getTmleDiagnostics(ConcreteEst, type = "components")
)
```

This is usually enough to understand whether a convergence issue is
numerical, data-sparsity related, or learner related.
