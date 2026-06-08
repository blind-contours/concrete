# Win ratios for trialists: single, hierarchical, and clinical

The **win ratio** compares two treatment arms by, in effect, playing
every treated patient against every control patient and asking who does
better. It is popular for composite and hierarchical endpoints
(e.g. cardiovascular death and heart-failure hospitalization) and under
non-proportional hazards, where the hazard ratio is hard to interpret.
`concrete` provides three versions, all built as smooth functionals of
the **targeted counterfactual survival / cumulative- incidence curves**
it already estimates — so unlike the standard pairwise win ratio, they
are **covariate-adjusted, doubly-robust, and censoring-corrected**, with
influence-function confidence intervals.

This article explains the three, when to use each, and shows the
simulation evidence that they are calibrated.

## 1. Single endpoint

For one time-to-event outcome,
[`getWinRatio()`](https://blind-contours.github.io/concrete/reference/getWinRatio.md)
reports the restricted win ratio, win odds, and net benefit over
`[0, Horizon]`: a treated patient *wins* if the control has the event
first (within the horizon), *loses* if the treated has it first, and
*ties* otherwise.

``` r

library(concrete)
ConcreteEst <- doConcrete(formatArguments(
  DataTable = trial, EventTime = "time", EventType = "event", Treatment = "arm",
  ID = "id", Intervention = makeITT(), TargetTime = c(180, 365, 545, 730),
  TargetEvent = 1, CVArg = list(V = 5)))

getWinRatio(ConcreteEst, Horizon = 730, Intervention = c(1, 2))
#> Win Ratio, Win Odds, Net Benefit, P(win), P(loss), P(tie) -- each with a CI
```

## 2. Hierarchical (prioritized) endpoints

Give `TargetEvent` an **ordered** vector (highest priority first) and
the comparison becomes hierarchical (Pocock / Finkelstein–Schoenfeld):
decide each pair on the top-priority event; only if that ties, move to
the next.

``` r

# event 1 (e.g. death) outranks event 2 (e.g. hospitalization)
getWinRatio(ConcreteEst, Horizon = 730, Intervention = c(1, 2), TargetEvent = c(1, 2))
```

This version compares each patient’s **first** event, treating the
listed events as competing risks (the structure `concrete` models). It
is exact and fast, but note the consequence: a patient hospitalized
**and then** dying is recorded by their first event (the
hospitalization), so their later death is not used. When death after a
non-fatal event matters, use the clinical win ratio below.

## 3. Clinical (death-priority) win ratio — experimental

[`clinicalWinRatio()`](https://blind-contours.github.io/concrete/reference/clinicalWinRatio.md)
targets the win ratio cardiologists usually mean: **compare on death
first — counting death even when it follows the non-fatal event — then
break ties on the non-fatal event.** This requires modelling the full
illness–death structure (alive → non-fatal → death, plus alive → death),
which `concrete` does with three transition hazards, each a Super
Learner; the post-non-fatal death hazard is fit on a left-truncated risk
set (subjects enter at their non-fatal-event time). Inference is
doubly-robust, covariate-adjusted, and censoring-corrected (IPCW).

It currently takes its own one-row-per-subject multistate data frame:

``` r

# columns: arm (1=active), t_hfh (non-fatal-event time; NA if none),
#          t_term (death or censoring time), died (1/0), plus covariates
clinicalWinRatio(
  trial,
  arm = "arm", illness.time = "t_hfh",
  terminal.time = "t_term", terminal.status = "died",
  covariates = c("age", "sex"), horizon = 1460)
#> Win Ratio / Win Odds / Net Benefit / P(win,loss,tie), each with an IF CI
```

## Can you trust it? Simulation evidence

Every version is validated against **ground truth** — closed-form path
probabilities and a brute-force pairwise win ratio computed on full
simulated histories (which counts death-after-non-fatal correctly).
Across the win-ratio family, the 95% influence-function confidence
intervals cover at the nominal rate:

![win ratio CI coverage near 0.95](figures/winratio-coverage.png)

And the clinical (death-priority) win ratio recovers the brute-force
pairwise truth essentially exactly, with and without censoring:

![clinical win ratio recovers ground
truth](figures/clinical-winratio-truth.png)

(Reproduce with `scripts/make-winratio-coverage.R`,
`scripts/make-clinical-wr-censoring-validation.R`, and the other
`scripts/make-clinical-wr-*.R`.)

## Which one should I use?

| Situation | Use |
|----|----|
| One time-to-event endpoint | [`getWinRatio()`](https://blind-contours.github.io/concrete/reference/getWinRatio.md) (single `TargetEvent`) |
| Prioritized endpoints, mutually-exclusive first events | [`getWinRatio()`](https://blind-contours.github.io/concrete/reference/getWinRatio.md) (ordered `TargetEvent`) |
| Death must count *even after* a non-fatal event (e.g. CV death \> HF hospitalization) | [`clinicalWinRatio()`](https://blind-contours.github.io/concrete/reference/clinicalWinRatio.md) |

In all three, a win ratio above 1 (and net benefit above 0) favors the
active arm, and the win ratio’s CI crossing 1 is the test of no
difference.

## Caveats

- **[`clinicalWinRatio()`](https://blind-contours.github.io/concrete/reference/clinicalWinRatio.md)
  is experimental.** It uses its own multistate data interface (not yet
  the \[formatArguments()\] pipeline), assumes a **single** non-fatal
  event type, **conditionally-independent censoring** (CAR), and a
  **Markov** illness–death model (the post-non-fatal death hazard
  depends on calendar time, not time-since-event). A semi-Markov option
  and tighter integration are planned.
- Use a reasonably dense time grid (`n.grid`); like RMST, the
  path-probability integrals carry a small grid-discretization bias that
  shrinks with the grid.
- For the single/hierarchical versions, use a reasonably dense
  `TargetTime` grid up to the horizon. \`\`\`
