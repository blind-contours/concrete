# updateHazardsWithEps

Illustrative helper function that modifies local hazards based on
epsVec. You must define how each eps entry maps to a hazard shift or
scale.

## Usage

``` r
updateHazardsWithEps(
  localEst,
  epsVec,
  Data,
  TargetEvent,
  TargetTime,
  Verbose = FALSE
)
```

## Arguments

- localEst:

  ephemeral copy of the Estimates list

- epsVec:

  numeric vector (length = k)

- Data:

  data.table

- TargetEvent:

  numeric

- TargetTime:

  numeric

- Verbose:

  boolean

## Value

The updated localEst (ephemeral)
