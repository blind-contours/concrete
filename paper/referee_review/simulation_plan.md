# Simulation Plan for Referee Response

This plan targets the simulation evidence requested in the referee review copied to
`paper/referee_review/concrete_referee_report.docx`.

## Review-Driven Questions

1. How often does the one-step TMLE update converge in realistic finite samples?
2. When convergence is slow, what problem features predict it: rare early events,
   many target times/events, practical positivity, censoring, or small step sizes?
3. What guidance should users follow when convergence fails or requires many updates?
4. Does a continuous-time estimator with Cox/Coxnet hazard learners perform well
   relative to discrete-time TMLE benchmarks across discretization choices?
5. How much does the current Cox-only hazard learner limitation matter under
   proportional hazards and non-proportional hazards data-generating mechanisms?
6. What are the runtime and memory costs users should expect?

## Estimators to Compare

Primary concrete estimators:

- `concrete` TMLE with `UpdateMethod = "standard"`.
- `concrete` TMLE with `UpdateMethod = "adaptive"` as a convergence sensitivity.
- `concrete` g-formula plug-in estimates from the same initial hazards.

Concrete nuisance libraries:

- Minimal hazard library: treatment-only and main-effects Cox models.
- Rich Cox hazard library: main effects plus targeted interactions motivated by
  the data-generating mechanism.
- Coxnet hazard library where supported by the current package.
- Propensity score SuperLearner library: start with `SL.glm` and `SL.glmnet`;
  add `SL.xgboost`/`SL.ranger` only in the main run if pilot runtime is acceptable.

Benchmarks:

- Aalen-Johansen, unadjusted.
- `survtmle` with 1-month, 3-month, and 6-month discretization grids.
- Optional extension benchmark: if we implement a broader hazard learner
  interface, add survivalSL/survSuperLearner/survML/additive-hazard candidates.

## Data-Generating Scenarios

Use the existing `simConCR()`/`getTrueRisks()` scaffold, but run from a new driver
with only project-relative paths.

1. Proportional hazards, correctly specified:
   Cox models are well specified. This isolates baseline performance and
   convergence of the current implementation.
2. Proportional hazards with nonlinear covariate effects/interactions:
   Tests whether richer Cox/Coxnet libraries reduce bias and improve coverage.
3. Non-proportional hazards:
   Include time-varying treatment or covariate effects to directly address the
   reviewers' Cox-only concern.
4. Rare-event/early-target stress:
   Include early target times where few events have occurred to quantify the
   slow-convergence issue noted in the manuscript.
5. Positivity/censoring stress:
   Increase treatment imbalance and informative censoring to quantify truncation,
   bias, coverage, and convergence behavior.

## Targets and Truth

Target static interventions `A=0` and `A=1`, censoring prevented.

Target times:

- Pilot: `c(365, 730, 1095, 1460)`.
- Main: `c(180, 365, 730, 1095, 1460)`.

Target events:

- Primary: event `1`.
- Secondary: events `1` and `2` jointly, to stress multidimensional targeting.

Truth:

- Estimate true risks by large Monte Carlo under each intervention using
  censoring disabled, with at least `5e6` simulated subjects per scenario if
  runtime permits.
- Cache truth as scenario-specific RDS/CSV files.

## Simulation Size

Pilot:

- `B = 25` replicates per scenario.
- `n = 500` and `n = 1000`.
- `CVArg = list(V = 2)` for quick debugging.

Main:

- `B = 500` replicates per scenario, increasing to `B = 1000` for final paper
  tables if runtime is acceptable.
- `n = 500`, `n = 1000`, and `n = 2500`.
- `CVArg = list(V = 5)`.

## Metrics

Estimation:

- Bias and percent bias.
- Empirical standard deviation.
- Mean estimated standard error.
- RMSE.
- Pointwise 95% confidence interval coverage.
- Simultaneous band coverage for multi-time targets where available.

Convergence and diagnostics:

- Convergence rate.
- Number of accepted update steps.
- Final norm of empirical mean EIC.
- Final maximum `|PnEIC| / stopping threshold`.
- Number of target components still failing the stopping criterion.
- Runtime per replicate.
- Error rate.
- Fraction of observations affected by nuisance truncation.
- SuperLearner selected hazards and propensity learners.

## Implementation Steps

1. Create a new project-relative simulation directory under
   `scripts/sim-data/referee-sims/`.
2. Replace the old absolute-path simulation driver with a clean driver split into:
   data generation, truth generation, estimator runners, replicate execution,
   result aggregation, and figures/tables.
3. Add a small pilot script that can complete locally and writes structured RDS
   output per replicate.
4. Run the pilot first to validate convergence metadata, target-event handling,
   output schema, and runtime.
5. Review pilot diagnostics and prune any estimator/library combination that is
   too unstable or too slow for the main run.
6. Run the main simulation grid.
7. Produce referee-response tables:
   convergence summary, runtime summary, bias/RMSE/coverage summary, and
   discretization sensitivity for `survtmle`.
8. Use the results to write practical guidance for users:
   when to increase `MaxUpdateIter`, when to reduce target dimensionality, when
   slow convergence signals sparse/rare-event targets, and when nuisance
   truncation or bootstrap sensitivity should be considered.

## Immediate Recommendation

Run the first pilot on scenarios 1, 3, and 4 only. That gives direct evidence for
the strongest reviewer concerns before spending compute on the full grid.
