# Specify an ICH E9(R1) estimand and apply an intercurrent-event strategy

`makeEstimand()` records the analysis target using the five attributes
of the ICH E9(R1) estimand framework (treatment condition, population,
variable / endpoint, intercurrent-event handling, and population-level
summary). It does not change the estimation; it documents the target for
a statistical analysis plan and travels with the results for
reproducibility and regulatory clarity.

`applyIntercurrentEvent()` implements the data handling for the
intercurrent-event (ICE) strategy. An ICE is coded as one of the
event-type values in the data (the observed time is the earliest of the
outcome, competing, censoring, and ICE times). Supported strategies:

- `"treatment policy"`:

  The ICE is ignored and follow-up continues to the outcome – the
  default intent-to-treat target. The data are returned unchanged (the
  ICE code, if any, is treated as a competing/censoring event exactly as
  supplied).

- `"hypothetical"`:

  The ICE is recoded as censoring (event type `0`). `concrete`'s
  inverse-probability-of-censoring weighting then targets the risk that
  would be observed if the ICE had not occurred, **under the assumption
  that the ICE acts as conditionally independent censoring given the
  measured covariates** – a strong, untestable assumption that should be
  accompanied by a censoring sensitivity analysis.

- `"composite"`:

  The ICE is recoded as an occurrence of the event of interest
  `TargetEvent`, defining a composite endpoint.

The `"while on treatment"` and `"principal stratum"` strategies are not
supported: the former is a different (descriptive, non-extrapolated)
estimand that an IPCW estimator does not target, and the latter requires
latent stratification beyond the current scope.

## Usage

``` r
makeEstimand(
  treatment = "as randomized (intent-to-treat)",
  population = "as enrolled",
  endpoint = "cause-specific absolute risk",
  summary = "risk difference",
  strategy = c("treatment policy", "hypothetical", "composite"),
  intercurrent = NULL
)

# S3 method for class 'ConcreteEstimand'
print(x, ...)

applyIntercurrentEvent(
  Data,
  EventTime,
  EventType,
  Intercurrent,
  strategy = c("treatment policy", "hypothetical", "composite"),
  TargetEvent = 1,
  intercurrent = "intercurrent event",
  Verbose = TRUE
)
```

## Arguments

- treatment:

  character: the treatment condition / intervention contrast.

- population:

  character: the target population.

- endpoint:

  character: the variable / endpoint (e.g. cause-specific absolute risk,
  RMST) and its horizon.

- summary:

  character: the population-level summary (e.g. risk difference, risk
  ratio, RMST difference).

- strategy:

  one of `"treatment policy"`, `"hypothetical"`, `"composite"`.

- intercurrent:

  character or NULL: a label for the intercurrent event.

- x:

  a `"ConcreteEstimand"` object.

- ...:

  ignored.

- Data:

  a `data.table`/`data.frame` of one row per subject.

- EventTime:

  character: name of the observed-time column.

- EventType:

  character: name of the event-type column (`0` = censoring, positive
  integers = events).

- Intercurrent:

  numeric: the event-type value in `EventType` that codes the
  intercurrent event.

- TargetEvent:

  numeric: for the `"composite"` strategy, the event-of- interest code
  that the intercurrent event is merged into.

- Verbose:

  logical: report what was recoded.

## Value

`makeEstimand()` returns a `"ConcreteEstimand"` object.

`applyIntercurrentEvent()` returns a copy of `Data` with `EventType`
recoded for the chosen `strategy`, carrying the estimand in
`attr(., "Estimand")`.
