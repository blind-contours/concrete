# Covariate-adjusted restricted win ratio, win odds, and net benefit

`getWinRatio()` estimates the **restricted win ratio** and its relatives
(win odds, net benefit) from the covariate-adjusted, censoring-corrected
counterfactual cumulative-incidence curves that `concrete` targets. It
supports a single time-to-event outcome or a **prioritized hierarchy**
of competing events (e.g.\\ death \\\>\\ hospitalization \\\>\\ stroke).

**Single event.** Comparing a random treated patient to a random control
patient over \\\[0,\tau\]\\: the treated patient *wins* if the control
has the event first (and within the horizon), *loses* if the treated has
it first, and *ties* if both are event-free at \\\tau\\. By independence
of the two patients these are functionals of the marginal survival
\\\bar S_a\\ and cumulative incidence \\\bar F_a = 1 - \bar S_a\\:
\$\$P(\text{win}) = \int_0^\tau \bar S_1(t)\\ d\bar F_0(t), \qquad
P(\text{loss}) = \int_0^\tau \bar S_0(t)\\ d\bar F_1(t).\$\$

**Hierarchy.** When `TargetEvent` lists several event codes in priority
order (highest priority first), the comparison is the prioritized
(Pocock / Finkelstein–Schoenfeld) rule applied to the patients' *first*
events: a patient who is event-free beats one who had any event; between
two patients with events of *different* priority, the one whose event is
*lower* priority (less severe) wins; between two with the *same* event,
the one whose event is *later* wins. Writing the per-arm cause-specific
cumulative incidences \\F_a^{(k)}\\ for priority \\k\\ (\\k=1\\
highest), the win probability is again a smooth functional of those
marginal curves, \$\$P(\text{win}) = S_1(\tau)\bigl(1 -
S_0(\tau)\bigr) + \sum\_{a\>b} F_1^{(a)}(\tau)F_0^{(b)}(\tau) + \sum_k
\int_0^\tau \bigl\[F_1^{(k)}(\tau) - F_1^{(k)}(t)\bigr\]\\
dF_0^{(k)}(t),\$\$ with \\S_a(\tau) = 1 - \sum_k F_a^{(k)}(\tau)\\, and
\\P(\text{loss})\\ the mirror image. This reduces exactly to the
single-event formula when one event is given.

In all cases the win/loss probabilities are smooth functionals of the
targeted curves, so their influence functions are weighted combinations
of the per-subject curve influence functions and the win ratio, win
odds, and net benefit follow by the delta method — giving doubly-robust,
covariate-adjusted, censoring-corrected inference, unlike the standard
unadjusted, censoring- sensitive win ratio. The integral is taken over
the fitted target times, so use a reasonably dense `TargetTime` grid.

**First-event limitation — use
[`clinicalWinRatio()`](https://blind-contours.github.io/concrete/reference/clinicalWinRatio.md)
for most trials.** The prioritized comparison here uses each patient's
*first* observed event, treating the listed events as competing risks. A
**higher-priority event that follows a lower-priority one is therefore
not counted** — death after a non-fatal hospitalization, or a stroke
after a hospitalization. For the clinically intended hierarchy (compare
on the most serious event first, *whenever* it occurs), which is what
most trials mean, use
[`clinicalWinRatio()`](https://blind-contours.github.io/concrete/reference/clinicalWinRatio.md):
it estimates the death-priority win ratio over an ordered hierarchy of
fatal and non-fatal time-to-event tiers via a multistate model.
`getWinRatio()` is the right choice only when events are genuinely
mutually exclusive (a higher-priority event can never follow a
lower-priority one).

## Usage

``` r
getWinRatio(
  ConcreteEst,
  Horizon = NULL,
  Intervention = c(1, 2),
  TargetEvent = NULL,
  Signif = 0.05
)
```

## Arguments

- ConcreteEst:

  a `"ConcreteEst"` object from
  [`doConcrete()`](https://blind-contours.github.io/concrete/reference/doConcrete.md).

- Horizon:

  numeric: the restriction horizon \\\tau\\ (default: the largest target
  time).

- Intervention:

  length-2 numeric: treatment and control indices.

- TargetEvent:

  numeric: the event code, or an ordered vector of event codes giving
  the priority hierarchy from highest to lowest (default: the first
  targeted event).

- Signif:

  numeric (default 0.05): alpha for confidence intervals and p-values.
  Win ratio and win odds are inferred on the log scale.

## Value

a `data.table` of class `"ConcreteOut"` with the win ratio, win odds,
net benefit, and the win/loss/tie probabilities, each with a CI and (for
the comparative statistics) a p-value against the null of no difference.
If the fit was built with `Strata` (see
[`formatArguments()`](https://blind-contours.github.io/concrete/reference/formatArguments.md)),
the standard errors are corrected for the stratified /
covariate-adaptive randomization design.

## See also

[`getRMST()`](https://blind-contours.github.io/concrete/reference/getRMST.md),
[`getOutput()`](https://blind-contours.github.io/concrete/reference/getOutput.md)
