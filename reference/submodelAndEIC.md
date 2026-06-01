# submodelAndEIC

Evaluate the empirical mean EIC (length k) for a guessed epsVec, without
modifying your real 'Estimates'.

Applies a multi-dimensional universal least-favorable submodel
parameterized by epsVec. Then returns the k-dimensional vector of
empirical means of the efficient influence curve for each coordinate.

## Usage

``` r
submodelAndEIC(
  Estimates,
  epsVec,
  Data,
  TargetEvent,
  TargetTime,
  Verbose = FALSE
)

submodelAndEIC(
  Estimates,
  epsVec,
  Data,
  TargetEvent,
  TargetTime,
  Verbose = FALSE
)
```

## Arguments

- Estimates:

  A list of arms, each containing:

  - Hazards: a named list of hazard matrices (rows = times, columns =
    subjects)

  - EvntFreeSurv: a matrix of survival probabilities

  - NuisanceWeight: a matrix of 1/(g\*(censorSurv)) or similar

  - PropScore: stored g-star?

- epsVec:

  numeric vector of length k = (#arms \* \#events \* \#times).

- Data:

  data.table with T.tilde, Delta, etc. (as in doTmleUpdate)

- TargetEvent:

  numeric vector of event types (e.g. c(1,2))

- TargetTime:

  numeric vector of target times (e.g. c(365, 730))

- Verbose:

  logical, if TRUE prints some progress

## Value

numeric vector of length k, the Pn(EIC) values

A numeric vector of length k. This is the stack of Pn(EIC_j) for j=1..k.
The root solver tries to make all of them ~0.
