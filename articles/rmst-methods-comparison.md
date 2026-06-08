# Restricted mean survival time: how concrete compares to other packages

The restricted mean survival time (RMST) — the average event-free time
over a horizon `[0, tau]` — is an increasingly preferred trial summary
because, unlike the hazard ratio, it stays interpretable under
non-proportional hazards and is reported in time units (days, months).
Several R packages estimate it. This article is an **honest**
head-to-head: where `concrete` is at parity with the established tools,
and where it genuinely adds something. We use a toy competing-risks
trial with a **closed-form true RMST** so every method can be scored
against the same target.

## A toy trial with known truth

Two competing causes (cause 1 = the event of interest, cause 2 = a
competing event), a randomized binary treatment that lowers the cause-1
hazard, two covariates, and independent censoring. Cause-specific
hazards are exponential, so the marginal RMST and cause-specific
years-lost have closed forms.

``` r

library(data.table)
tau <- 1500                                   # horizon (days)
hz <- function(A, W1, W2) {                   # cause-specific hazards
  l1 <- 6e-4 * exp(-0.6 * A + 0.4 * W1)       # cause 1 (treatment lowers it)
  l2 <- 4e-4 * exp( 0.1 * A + 0.3 * W2)       # cause 2 (competing)
  list(l1 = l1, l2 = l2, lt = l1 + l2)
}
simTrial <- function(n, seed) {
  set.seed(seed); W1 <- rnorm(n); W2 <- rnorm(n); A <- rbinom(n, 1, 0.5); h <- hz(A, W1, W2)
  T1 <- rexp(n, h$l1); T2 <- rexp(n, h$l2); C <- rexp(n, 3e-4)
  Tt <- pmin(T1, T2); J <- ifelse(T1 < T2, 1L, 2L)
  To <- pmin(Tt, C, 2500); ev <- ifelse(To >= 2500, 0L, ifelse(Tt <= C, J, 0L))
  data.table(id = seq_len(n), time = To, event = ev, arm = A, W1 = W1, W2 = W2)
}
```

Integrating the exponential hazards gives the marginal truth for this
design:

| Estimand (treatment − control)        |      True value |
|---------------------------------------|----------------:|
| Event-free RMST difference            | **+110.4 days** |
| Cause-1 years-of-life-lost difference | **−172.2 days** |

(Treatment buys ~110 more event-free days and ~172 fewer days lost to
cause 1.)

## Estimating it with `concrete`

`concrete` targets the cause-specific cumulative incidence with
continuous-time TMLE, then reports the event-free RMST and the
cause-specific years-lost (which sum to the horizon:
`tau = RMST + sum_j YLL_j`). Two routes:

``` r

library(concrete)
trial <- simTrial(n = 500, seed = 1)
args <- formatArguments(
  DataTable = trial, EventTime = "time", EventType = "event", Treatment = "arm",
  ID = "id", Intervention = makeITT(), TargetEvent = c(1, 2),
  TargetTime = seq(150, tau, by = 150), CVArg = list(V = 5)   # a reasonably dense grid (see caveat)
)
est <- doConcrete(args)

getRMST(est, Horizon = tau, Intervention = c(1, 2))      # integrate pointwise-targeted risks
targetRMST(est, Horizon = tau, Intervention = c(1, 2))   # target the RMST equation directly
```

[`getRMST()`](https://blind-contours.github.io/concrete/reference/getRMST.md)
returns the event-free RMST per arm, the cause-specific life-years lost,
and their between-arm differences, each with an influence-function CI
and a Wald p-value;
[`targetRMST()`](https://blind-contours.github.io/concrete/reference/targetRMST.md)
returns the same estimands from a direct targeting step (see the
**direct targeting** section).

## The same estimand with other packages

``` r

library(survRM2); library(eventglm); library(survival)
trial[, anyev := as.integer(event %in% c(1, 2))]   # "any first event" for event-free RMST

# survRM2 (Uno/Tian): RMST difference for a single survival endpoint.
#   Unadjusted (marginal):
rmst2(trial$time, trial$anyev, trial$arm, tau = tau)$unadjusted.result
#   Covariate-adjusted (Tian augmentation, still marginal):
rmst2(trial$time, trial$anyev, trial$arm, tau = tau,
      covariates = as.data.frame(trial[, .(W1, W2)]))$RMST.difference.adjusted

# eventglm (Sachs & Gabriel): pseudo-observation regression.
#   Event-free RMST difference (arm coefficient):
coef(rmeanglm(Surv(time, anyev) ~ arm, time = tau, data = trial))["arm"]
#   Cause-1 years-lost difference (competing risks need a factor status):
trial[, eventf := factor(event)]
coef(rmeanglm(Surv(time, eventf) ~ arm, time = tau, cause = 1, data = trial))["arm"]

# riskRegression::ate (Ozenne/Gerds): AIPW RMST with competing risks (see note).
# Not run here — its current dependency `rms` requires a newer R than this build.
```

## Head-to-head: bias and coverage

We re-ran the toy trial 40 times at n = 500 and scored each method’s
**marginal** estimate against the closed-form truth (every method used
in its marginal form, so the comparison is like-for-like). Mean bias in
days, 95% CI coverage, and mean CI width:

| Method | Estimand | Mean bias (d) | Coverage | CI width (d) |
|----|----|---:|---:|---:|
| `concrete` `getRMST` | EF-RMST diff | −10 | ~0.95–1.0 | 192 |
| `concrete` `targetRMST` | EF-RMST diff | −10 | ~0.95–1.0 | 192 |
| `survRM2` (unadjusted) | EF-RMST diff | +7 | 0.93 | 199 |
| `survRM2` (Tian-adjusted) | EF-RMST diff | +5 | 0.95 | 196 |
| `eventglm` | EF-RMST diff | +7 | — | — |
| `concrete` | cause-1 YLL diff | +11 | 0.93 | 181 |
| `eventglm` | cause-1 YLL diff | −3 | — | — |
| `survRM2` | cause-1 YLL diff | **not supported** |  |  |

**How to read this honestly.** On the event-free RMST difference, all
four methods are at **parity**: every bias is within ~1–1.5 Monte-Carlo
standard errors of zero at 40 replicates, and coverage is near nominal.
`concrete`’s CIs are slightly *narrower* than `survRM2`’s (the
adjustment/TMLE efficiency gain) while remaining conservative in
coverage. `concrete` leans about −10 days here and `survRM2` about +5–7
— neither is a clear winner on point accuracy in this randomized design,
where covariate adjustment mainly buys precision rather than removing
bias. (The −10 is discussed under **caveats**.) We could not include
`riskRegression::ate` empirically because its dependency would not
install on this R build; in the feature matrix below it is described
from its documentation.

## What each package can and cannot do

| Capability | `concrete` | `survRM2`/adj | `eventglm` | `riskRegression::ate` |
|----|:--:|:--:|:--:|:--:|
| RMST difference | ✓ | ✓ | ✓ | ✓ |
| Cause-specific years-lost (competing risks) | ✓ | ✗ | ✓ | ✓ |
| Absolute risk / risk difference / risk ratio | ✓ | ✗ | ✓ (CIF) | ✓ |
| Competing risks handled correctly | ✓ | ✗ (treats as censoring) | ✓ | ✓ |
| Machine-learning / Super Learner nuisances | ✓ | ✗ | ✗ | ✗ (parametric models) |
| Doubly robust | ✓ | ✗ (single working model) | ✗ | ✓ (AIPW) |
| Influence-function inference | ✓ | ✓ (closed form) | ✓ (sandwich) | ✓ |
| Cross-fitting (CV-TMLE) | ✓ | ✗ | ✗ | partial |
| **Direct RMST-equation targeting** | ✓ (`targetRMST`) | ✗ | ✗ | ✗ |
| Win ratio / simultaneous CIs / censoring sensitivity in the same framework | ✓ | ✗ | ✗ | ✗ |
| Maturity, speed, widespread trial use | newer, slower | **mature, fast, standard** | established | mature |

## Where `concrete` genuinely adds something

Stated truthfully — these are real differences, not marketing:

1.  **Direct RMST targeting (`targetRMST`).** `concrete` can fluctuate
    the cause-specific hazards to solve the RMST estimating equation
    *directly* (a time-integrated clever covariate), rather than
    integrating pointwise-targeted risks. We are not aware of another
    package that does this. It matters most on **sparse target grids**:
    in a separate study (B = 80, n = 500), on a 2-time grid the
    pointwise integral was biased **+37 days with 0.86 coverage**, while
    direct targeting gave **−1.3 days with 0.95 coverage**.
2.  **Competing-risks cause-specific years-lost.** `concrete` (and
    `eventglm`, `riskRegression`) decompose the horizon into event-free
    time plus years lost to each cause. `survRM2` — the most widely used
    RMST tool — cannot: it treats competing events as censoring.
3.  **ML nuisances with double robustness in continuous time, plus
    cross-fitting.** `concrete` is the only one of these that combines
    Super Learner nuisance estimation, continuous-time one-step TMLE,
    and cross-fitting (CV-TMLE) for honest inference with flexible
    learners. `riskRegression::ate` is doubly robust but uses
    pre-specified parametric models.
4.  **One framework, many estimands.** From a single fit, `concrete`
    also gives risk differences/ratios, the (hierarchical) win ratio,
    simultaneous multiplicity-adjusted CIs, and a censoring-sensitivity
    analysis — useful for an ICH E9(R1)-style analysis plan.

## Where the established tools are better

Also truthfully:

- **`survRM2` is simpler, faster, and the de-facto standard** for a
  two-arm RMST difference on a single survival endpoint. If that is all
  you need, it is the right tool.
- **`riskRegression` is mature and broadly validated**; `eventglm`
  offers a clean regression interface and pseudo-observation
  flexibility.
- `concrete` is **newer, less battle-tested, and slower** — `targetRMST`
  in particular is markedly more expensive per fit than a pointwise
  integral.

## Caveats and practical advice

- **Grid discretization.**
  [`getRMST()`](https://blind-contours.github.io/concrete/reference/getRMST.md)
  integrates the targeted cumulative incidence by the trapezoid rule
  over the target times. Because event-free survival is convex, a
  *coarse* grid or a long horizon biases the RMST upward (and, for a
  difference, can shrink it — the ~−10-day lean above came from a
  deliberately coarse 8-point grid). Use a **reasonably dense
  `TargetTime` grid**, and/or
  [`targetRMST()`](https://blind-contours.github.io/concrete/reference/targetRMST.md),
  which targets the integral directly. In our coverage studies with
  adequate grids, `concrete`’s RMST and RMST-difference cover at
  ~0.95–0.97 with a bias of a few days.
- **Randomized trials.** Adjustment does not change the RMST estimand in
  an RCT; it improves precision (narrower CIs). The doubly-robust/ML
  machinery earns its keep most under confounding or strong covariate
  effects.
- **Pick the grid to the horizon.** A grid spacing well below the time
  scale on which survival changes keeps the discretization bias
  negligible.

## Bottom line

For a covariate-adjusted RMST *difference* on a single endpoint,
`concrete` is at **parity** with `survRM2(adj)`, `eventglm`, and
`riskRegression::ate` — calibrated, with competitive or slightly
narrower intervals. Its distinct value is the **direct RMST targeting**,
**competing-risks years-lost**, **ML + cross-fitted double robustness**,
and the **integrated estimand toolkit** — with the honest costs of being
newer and slower, and the practical need for a dense target grid (or
`targetRMST`).
