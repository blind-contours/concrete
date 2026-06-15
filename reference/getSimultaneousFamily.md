# Joint simultaneous inference across a family of estimands

Trial reports rarely contain a single number: a typical analysis
presents a risk difference at several horizons, the RMST difference, and
a win ratio together. Pointwise 95% intervals do not control the
family-wise error across such a set, and a Bonferroni correction is
needlessly conservative because the estimands are computed from the
*same* targeted curves and are therefore strongly correlated.
`getSimultaneousFamily()` builds **simultaneous** (family-wise)
confidence bands across an arbitrary collection of `concrete` estimands
by exploiting that correlation.

Every `concrete` estimate carries its per-subject efficient influence
values. Stacking them into one \\n \times q\\ matrix gives the joint
influence function of the whole family; its empirical correlation \\R\\
is the correlation of the (asymptotically normal) estimators. The
simultaneous critical value is the \\1-\alpha\\ quantile of \\\max_j
\|Z_j\|\\ for \\Z \sim N(0, R)\\ (a multiplier/Gaussian-multiplier
bootstrap), and each band is \\\hat\psi \pm q\\\widehat{\mathrm{se}}\\
(or, for ratio estimands, the same on the log scale). When the family is
a single estimand at one time this reduces to the usual Wald interval.

This composes with everything else in the package: the influence values
used are exactly those reported by
[`getOutput()`](https://blind-contours.github.io/concrete/reference/getOutput.md),
[`getRMST()`](https://blind-contours.github.io/concrete/reference/getRMST.md),
[`targetRMST()`](https://blind-contours.github.io/concrete/reference/targetRMST.md),
[`getWinRatio()`](https://blind-contours.github.io/concrete/reference/getWinRatio.md),
and
[`targetWinRatio()`](https://blind-contours.github.io/concrete/reference/targetWinRatio.md),
so the simultaneous bands inherit the covariate adjustment, censoring
correction, cross-fitting, and (when present) the
stratified-randomization variance correction.

## Usage

``` r
getSimultaneousFamily(..., Signif = 0.05, nSim = 10000L)
```

## Arguments

- ...:

  two or more `"ConcreteOut"` objects produced from the **same** fitted
  [`doConcrete()`](https://blind-contours.github.io/concrete/reference/doConcrete.md)
  object (so that subjects align). Each may be named; names are carried
  into the output `family` column.

- Signif:

  numeric (default 0.05): family-wise alpha.

- nSim:

  integer (default 1e4): Monte Carlo draws for the multiplier bootstrap
  critical value.

## Value

a `data.table` with one row per scalar estimand in the family: `family`,
`Estimand`, `Event`, `Time`, `Intervention`, `Pt Est`, `se`, the
pointwise `CI Low` / `CI Hi`, and the simultaneous `SimCI Low` /
`SimCI Hi`, plus the shared simultaneous critical value in
`attr(., "critValue")`.

## See also

[`getOutput()`](https://blind-contours.github.io/concrete/reference/getOutput.md),
[`getRMST()`](https://blind-contours.github.io/concrete/reference/getRMST.md),
[`targetRMST()`](https://blind-contours.github.io/concrete/reference/targetRMST.md),
[`getWinRatio()`](https://blind-contours.github.io/concrete/reference/getWinRatio.md),
[`targetWinRatio()`](https://blind-contours.github.io/concrete/reference/targetWinRatio.md).
