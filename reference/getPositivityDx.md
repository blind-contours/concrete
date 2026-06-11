# Positivity / inverse-weight diagnostics for a fitted estimate

Reports the practical-positivity health of the inverse-probability
weights that every `concrete` estimand relies on. The nuisance weight is
\\1/(g(A\mid W)\\ S_C(t)\\ S_X(t))\\ — the inverse of the probability of
being assigned the regime's treatment *and* remaining uncensored
(\\S_C\\) *and*, when a crossover model is used, not yet switched
(\\S_X\\). Because these are multiplied, the denominator can become very
small at later times, which (i) inflates the influence-function variance
and (ii) triggers truncation that can bias the estimate. This is exactly
the regime to watch with informative censoring and crossover.

For each intervention it returns the **effective sample size**
\\\mathrm{ESS}(t) = (\sum_i w\_{it})^2 / \sum_i w\_{it}^2\\ (as a
fraction of \\n\\), the largest weight, the smallest observation
probability (the positivity floor), and the share of weights sitting at
the truncation bound — overall and at the worst time point. Read it
alongside any estimate to judge whether the inference is trustworthy or
weight-limited.

## Usage

``` r
getPositivityDx(ConcreteEst, Verbose = TRUE)
```

## Arguments

- ConcreteEst:

  a `"ConcreteEst"` object from
  [`doConcrete()`](https://blind-contours.github.io/concrete/reference/doConcrete.md).

- Verbose:

  logical (default TRUE): print a short interpreted summary.

## Value

invisibly, a list with `summary` (one row per intervention) and `byTime`
(the per-evaluation-time ESS fraction, max weight, and minimum
observation probability for each intervention).

## Examples

``` r
if (FALSE) { # \dontrun{
est <- doConcrete(formatArguments(...))
getPositivityDx(est)
} # }
```
