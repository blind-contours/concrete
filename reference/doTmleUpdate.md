# TMLE update with multiple implementation options

TMLE update with multiple implementation options

## Usage

``` r
doTmleUpdate(
  Estimates,
  SummEIC,
  Data,
  TargetEvent,
  TargetTime,
  MaxUpdateIter,
  OneStepEps,
  NormPnEIC,
  Verbose,
  Method = "standard",
  EICStopRule = "relative",
  EICStopAbsTol = 0
)
```

## Arguments

- Estimates:

  list of estimates for each arm/treatment, containing Hazards, Surv,
  SummEIC, etc.

- SummEIC:

  data.table summarizing the current empirical mean EIC (PnEIC), etc.

- Data:

  data.table containing the original data

- TargetEvent:

  numeric vector of event types

- TargetTime:

  numeric vector of target times

- MaxUpdateIter:

  maximum number of steps/iterations

- OneStepEps:

  initial step size for the incremental updates

- NormPnEIC:

  numeric, initial \|\|PnEIC\|\| for reference

- Verbose:

  boolean, whether to print debugging output

- Method:

  character - one of

  - "standard"

  - "adaptive"

  - "coordinated"

- EICStopRule:

  character stopping rule for empirical mean EIC checks. Supported
  values are `"relative"`, `"absolute"`, and `"hybrid"`.

- EICStopAbsTol:

  numeric absolute `|PnEIC|` tolerance used by the `"absolute"` and
  `"hybrid"` stopping rules.

## Value

Updated `Estimates` object, with updated hazards, SummEIC, etc.
