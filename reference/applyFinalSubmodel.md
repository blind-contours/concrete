# applyFinalSubmodel

Commits the final epsVec solution from nleqslv onto the real
'Estimates'.

## Usage

``` r
applyFinalSubmodel(
  Estimates,
  final_eps,
  Data,
  TargetEvent,
  TargetTime,
  Verbose = FALSE
)
```

## Arguments

- Estimates:

  the real estimates list to be updated in-place

- final_eps:

  numeric vector from root solver

- Data:

  data.table

- TargetEvent:

  numeric

- TargetTime:

  numeric

- Verbose:

  boolean

## Value

Updated `Estimates` object with final hazards + EIC
