# Covariate-adjusted restricted win ratio, win odds, and net benefit

For a single terminal time-to-event outcome, `getWinRatio()` estimates
the **restricted win ratio** and its relatives (win odds, net benefit)
from the covariate-adjusted, censoring-corrected counterfactual survival
curves that `concrete` targets. Comparing a random treated patient to a
random control patient over `[0, Horizon]`:

- a treated patient **wins** if the control patient has the event first
  (before the treated patient and before the horizon),

- **loses** if the treated patient has the event first,

- **ties** otherwise (both event-free at the horizon).

These probabilities are functionals of the marginal counterfactual
survival \\\bar S_a\\ and cumulative incidence \\\bar F_a = 1 - \bar
S_a\\: \$\$P(\text{win}) = \sum\_{t_k \le \tau} \bar S_1(t_k)\\ d\bar
F_0(t_k), \qquad P(\text{loss}) = \sum\_{t_k \le \tau} \bar S_0(t_k)\\
d\bar F_1(t_k),\$\$ and the win statistics are \$\$\text{win ratio} =
\frac{P(\text{win})}{P(\text{loss})}, \quad \text{win odds} =
\frac{P(\text{win}) + P(\text{tie})/2}{P(\text{loss}) +
P(\text{tie})/2}, \quad \text{net benefit} = P(\text{win}) -
P(\text{loss}).\$\$ Because the win/loss probabilities are smooth
functionals of the targeted curves, their influence functions are
weighted combinations of the per-subject curve influence functions; the
win ratio, win odds, and net benefit then follow by the delta method,
giving doubly-robust, covariate-adjusted inference. Unlike the standard
(unadjusted, censoring-sensitive) win ratio, this is restricted to the
horizon and corrects for censoring through the same inverse-probability
machinery as the rest of the package.

This is the single-terminal-event version; a hierarchical /
competing-risk win ratio is planned. The integral is taken over the
fitted target times, so use a reasonably dense `TargetTime` grid.

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

  numeric: the single terminal event code (default: the first targeted
  event).

- Signif:

  numeric (default 0.05): alpha for confidence intervals and p-values.
  Win ratio and win odds are inferred on the log scale.

## Value

a `data.table` of class `"ConcreteOut"` with the win ratio, win odds,
net benefit, and the win/loss/tie probabilities, each with a CI and (for
the comparative statistics) a p-value against the null of no difference.

## See also

[`getRMST()`](https://blind-contours.github.io/concrete/reference/getRMST.md),
[`getOutput()`](https://blind-contours.github.io/concrete/reference/getOutput.md)
