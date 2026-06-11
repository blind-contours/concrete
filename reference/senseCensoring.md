# Tipping-point sensitivity analysis for informative censoring

`concrete`'s primary analysis assumes censoring is independent of the
event given the measured covariates. `senseCensoring()` probes
robustness to departures from that assumption with a transparent
tipping-point (bounds) analysis: a fraction `delta` of the subjects who
are censored before the target time are assumed to have actually
experienced the event of interest, and the analysis is re-fit for each
`delta`. `delta = 0` is the optimistic bound (censored subjects never
have the event), `delta = 1` the pessimistic bound (all do), and the
primary inverse-probability-of-censoring analysis sits between them. The
**tipping point** is the smallest `delta` at which the conclusion
changes – a sensitivity artifact recommended by ICH E9(R1).

Unlike scaling the censoring weight (which leaves a doubly-robust
estimator's target unchanged), imputing the event status of the censored
changes the estimand and so produces a genuine, interpretable
sensitivity curve. Because it re-fits the estimator for every `delta`,
it is computationally heavier than the primary analysis.

When the fit carries a treatment-switching (crossover) model (i.e.
[`formatArguments()`](https://blind-contours.github.io/concrete/reference/formatArguments.md)
was called with `Crossover`), the censored subjects are a mix of two
intercurrent events with two *different* untestable assumptions:
ordinary **dropout** (conditionally-independent censoring) and
**crossover** (the no-switching counterfactual). `mechanism` selects
which pool is imputed, so the two assumptions can be probed individually
or jointly:

- `"dropout"` – tip only the genuinely-censored (drop-out) subjects,
  holding the switching handled by the crossover hazard;

- `"crossover"` – tip only the subjects re-censored at their switch
  time, i.e. probe “what if switchers would have had the event had they
  not switched”, holding dropout handled by the censoring weight;

- `"all"` (default) – tip both pools jointly (the original behaviour).

With no crossover model, all censored subjects are dropout and the three
modes coincide.

## Usage

``` r
senseCensoring(
  ConcreteArgs,
  deltas = seq(0, 1, by = 0.25),
  Estimand = c("RD", "RR", "Risk"),
  Intervention = c(1, 2),
  mechanism = c("all", "dropout", "crossover"),
  Signif = 0.05,
  Verbose = FALSE
)
```

## Arguments

- ConcreteArgs:

  a `"ConcreteArgs"` object from
  [`formatArguments()`](https://blind-contours.github.io/concrete/reference/formatArguments.md).

- deltas:

  numeric in `[0, 1]`: fractions of pre-target-time censored subjects
  imputed as having the event of interest. Should include 0.

- Estimand:

  one of `"RD"` (default), `"RR"`, `"Risk"`.

- Intervention:

  length-2 numeric: treatment and control indices.

- mechanism:

  one of `"all"` (default), `"dropout"`, `"crossover"`: which pool of
  censored subjects to impute (see Details). `"crossover"` requires a
  fit built with `Crossover`.

- Signif:

  numeric (default 0.05): two-sided alpha.

- Verbose:

  logical.

## Value

a `data.table` of estimate / CI / p-value by `delta` x event x time,
with a leading `mechanism` column, the tipping point in
`attr(., "tippingPoint")`, and the mechanism in `attr(., "mechanism")`.
