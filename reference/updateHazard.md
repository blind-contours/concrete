# Update hazards based on clever covariate and PnEIC

Update hazards based on clever covariate and PnEIC

## Usage

``` r
updateHazard(
  GStar,
  Hazards,
  TotalSurv,
  NuisanceWeight,
  EvalTimes,
  T.tilde,
  Delta,
  PnEIC,
  NormPnEIC,
  OneStepEps,
  TargetEvent,
  TargetTime
)
```

## Arguments

- GStar:

  numeric vector of star probabilities

- Hazards:

  list of hazard matrices

- TotalSurv:

  matrix of survival probabilities

- NuisanceWeight:

  matrix of nuisance weights

- EvalTimes:

  numeric vector of evaluation times

- T.tilde:

  numeric vector of observed event times

- Delta:

  numeric vector of observed event types

- PnEIC:

  data.table with PnEIC values

- NormPnEIC:

  numeric norm of PnEIC

- OneStepEps:

  numeric step size for updates

- TargetEvent:

  numeric vector of target events

- TargetTime:

  numeric vector of target times

## Value

list of updated hazard matrices
