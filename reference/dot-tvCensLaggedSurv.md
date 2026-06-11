# Build the lagged censoring-survival matrix with time-varying covariates, in the core's `LaggedCensSurv` convention (matches `getHazSurvPred()`), for overriding the IPCW throughout the pipeline. `CensoringTV` is a long data.frame with the data's id column (`attr(Data, "ID")`), a `time` column, and value columns.

Build the lagged censoring-survival matrix with time-varying covariates,
in the core's `LaggedCensSurv` convention (matches `getHazSurvPred()`),
for overriding the IPCW throughout the pipeline. `CensoringTV` is a long
data.frame with the data's id column (`attr(Data, "ID")`), a `time`
column, and value columns.

## Usage

``` r
.tvCensLaggedSurv(
  Data,
  CensoringTV = NULL,
  times,
  Crossover = NULL,
  SL.library = c("SL.mean", "SL.glm"),
  n.folds = 5L,
  nGridCens = 40L
)
```

## Value

a `length(times) x n` matrix of lagged censoring survival.
