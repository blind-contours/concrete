# Survival SL pilot notes

Date: 2026-05-29

## Implementation

The package now supports a general cause-specific survival hazard library with a discrete Super Learner selector. Candidate learners can include:

- Existing Cox formula learners.
- `"coxnet"` via `glmnet::cv.glmnet`.
- `"rsf"` / `"randomForestSRC"` via `randomForestSRC::rfsrc`.
- `"aareg"` / `"additive_hazards"` via `survival::aareg`.
- `"hal"` / `"hal9001"` via a pooled discrete-time binomial HAL hazard on the continuous-time evaluation grid.

All candidates are compared with the same held-out counting-process negative log-likelihood based on predicted hazard increments:

`sum_i Lambda_j(T_i) - sum_i 1(Delta_i = j) log dLambda_j(T_i)`.

This avoids using Cox partial likelihood to compare non-Cox learners.

## Smoke checks

- Package install into `/private/tmp/concrete-lib`: passed.
- `test-formatArguments.R`: passed, 71 pass / 0 fail.
- `test-doConcrete.R`, with `data.table` attached for the existing test helper assumptions: passed, 5 pass / 0 fail.
- Small HAL-inclusive end-to-end smoke: passed.
- Small `run_one_replicate()` survival-SL no-HAL harness smoke: passed.

The full `testthat::test_dir()` invocation still fails in this ad hoc runner because the existing test files assume `data.table` is attached and the local test environment lacks `waldo`; this is separate from the survival-SL code path.

## Tiny saved pilot

Output directory:

`scripts/sim-data/referee-sims/output/survsl_pilot_nonph_b2_n150`

Design:

- Scenario: `nonph`
- Reps: `B = 2`
- Sample size: `n = 150`
- Target times: 365, 730
- Target events: 1, 2
- Survival SL: Cox + Coxnet + random survival forest + additive hazards; HAL excluded for this bulk pilot.

Top-line diagnostics:

- All concrete estimators converged in both reps.
- `concrete_standard_survsl_nohal`: convergence rate 1.00, median step 5.0, mean runtime about 2.1 sec.
- `concrete_adaptive_survsl_nohal`: convergence rate 1.00, median step 3.5, mean runtime about 2.1 sec.
- Adaptive survival-SL reduced update steps similarly to the earlier adaptive runs.

Interpretation:

This pilot only verifies wiring and runtime. It is too small for bias, coverage, or learner-selection conclusions. The next meaningful run should use the no-HAL survival-SL library over all referee scenarios, then run a smaller HAL sensitivity because pooled HAL scales with `n * number_of_grid_times`.

## B20 no-HAL referee grid

Output directory:

`scripts/sim-data/referee-sims/output/survsl_nohal_b20_n500_all_events`

Design:

- Scenarios: `ph_correct`, `nonph`, `rare_early`, `positivity`
- Reps: `B = 20`
- Sample size: `n = 500`
- Target times: 180, 365, 730, 1095, 1460
- Target events: 1, 2
- Survival SL: Cox + Coxnet + random survival forest + additive hazards; HAL excluded
- Update cap: 200 iterations
- CV folds for nuisance fitting: 2

Files written:

- `diagnostic_metrics.csv`
- `diagnostics_long.csv`
- `metrics.csv`
- `estimates_long.csv`
- `summary.rds`
- `row_status.csv`

Convergence and runtime:

| Scenario | Estimator | Convergence | Median step | Mean runtime sec |
|---|---:|---:|---:|---:|
| nonph | standard minimal | 0.95 | 55.5 | 103.5 |
| nonph | standard rich | 0.95 | 56.0 | 99.3 |
| nonph | adaptive rich | 0.95 | 51.5 | 117.1 |
| nonph | standard survSL no-HAL | 0.95 | 55.5 | 107.2 |
| nonph | adaptive survSL no-HAL | 0.95 | 48.0 | 125.4 |
| ph_correct | standard minimal | 1.00 | 32.0 | 67.6 |
| ph_correct | standard rich | 1.00 | 32.0 | 65.1 |
| ph_correct | adaptive rich | 1.00 | 30.5 | 76.3 |
| ph_correct | standard survSL no-HAL | 1.00 | 32.0 | 68.2 |
| ph_correct | adaptive survSL no-HAL | 1.00 | 30.5 | 81.4 |
| positivity | standard minimal | 0.90 | 82.0 | 197.7 |
| positivity | standard rich | 0.90 | 79.0 | 195.5 |
| positivity | adaptive rich | 0.90 | 71.5 | 216.3 |
| positivity | standard survSL no-HAL | 0.90 | 77.5 | 184.1 |
| positivity | adaptive survSL no-HAL | 0.90 | 73.0 | 211.7 |
| rare_early | standard minimal | 0.75 | 60.0 | 143.0 |
| rare_early | standard rich | 0.75 | 60.5 | 139.4 |
| rare_early | adaptive rich | 0.75 | 57.0 | 153.5 |
| rare_early | standard survSL no-HAL | 0.75 | 60.0 | 134.6 |
| rare_early | adaptive survSL no-HAL | 0.75 | 59.0 | 142.5 |

Nonconvergence was seed-level, not learner-specific. Every failed replicate failed for all five concrete estimators:

- `nonph`: rep 4
- `positivity`: reps 15 and 18
- `rare_early`: reps 8, 10, 11, 13, and 17

Mean absolute-risk performance across interventions, events, and target times:

| Scenario | Best concrete TMLE by mean RMSE | Mean RMSE | Mean absolute bias | Aalen-Johansen mean RMSE |
|---|---|---:|---:|---:|
| nonph | standard rich | 0.0218 | 0.0052 | 0.0252 |
| ph_correct | standard survSL no-HAL | 0.0242 | 0.0057 | 0.0270 |
| positivity | standard survSL no-HAL | 0.0250 | 0.0140 | 0.0283 |
| rare_early | standard rich | 0.0190 | 0.0063 | 0.0227 |

Mean risk-difference performance across interventions, events, and target times:

- `nonph`: all concrete TMLE variants were effectively tied; mean RMSE about 0.0260 to 0.0262.
- `ph_correct`: all concrete TMLE variants were effectively tied; mean RMSE about 0.0315.
- `positivity`: standard minimal/adaptive no-HAL were slightly better than standard rich; standard survSL no-HAL reduced mean absolute bias but not RMSE.
- `rare_early`: adaptive survSL no-HAL had the lowest mean RMSE, 0.0227, and lowest mean absolute bias, 0.0108.

Interpretation:

The no-HAL survival SL is stable and does not introduce errors. It also modestly improves update behavior in the hardest scenarios, especially `positivity` runtime for the standard update and `rare_early` runtime for both survival-SL variants. However, convergence failures are not fixed by changing the hazard learner library: failures are shared across all concrete estimators on the same replicate seeds. That points to the targeting/update path and stress-scenario positivity, not the initial hazard model alone, as the main convergence bottleneck.

For the referee response, this supports presenting survival-SL as a flexible nuisance-estimation extension, while treating slow/nonconvergence as a separate update-stability issue.

## HAL sensitivity

Output directory:

`scripts/sim-data/referee-sims/output/hal_sensitivity_nonph_b3_n300`

Design:

- Scenario: `nonph`
- Reps: `B = 3`
- Sample size: `n = 300`
- Target times: 365, 730
- Target events: 1, 2
- Survival SL: no-HAL variants plus a HAL-inclusive standard survival-SL variant
- Update cap: 100 iterations

Diagnostics:

| Estimator | Convergence | Median step | Max step | Mean runtime sec |
|---|---:|---:|---:|---:|
| standard minimal | 1.00 | 6 | 8 | 4.51 |
| standard rich | 1.00 | 5 | 6 | 2.46 |
| adaptive rich | 1.00 | 5 | 5 | 3.03 |
| standard survSL no-HAL | 1.00 | 6 | 6 | 3.96 |
| adaptive survSL no-HAL | 1.00 | 5 | 5 | 4.13 |
| standard survSL with HAL | 1.00 | 6 | 8 | 4.75 |

Mean absolute-risk performance across interventions, events, and target times:

| Estimator | Mean RMSE | Mean absolute bias | Mean coverage |
|---|---:|---:|---:|
| Aalen-Johansen | 0.0254 | 0.0154 | 0.958 |
| adaptive rich | 0.0238 | 0.0154 | 1.000 |
| adaptive survSL no-HAL | 0.0231 | 0.0146 | 1.000 |
| standard minimal | 0.0239 | 0.0151 | 1.000 |
| standard rich | 0.0234 | 0.0152 | 1.000 |
| standard survSL with HAL | 0.0239 | 0.0151 | 1.000 |
| standard survSL no-HAL | 0.0237 | 0.0151 | 1.000 |

Interpretation:

The HAL-inclusive survival-SL path runs and converges in this small sensitivity, but it does not show a clear accuracy gain at this size. It should remain in the implementation as an available library member, but the main simulation grid should default to no-HAL unless we allocate more compute or restrict the grid. A focused HAL sensitivity is more defensible than putting HAL into every full stress-scenario replicate.
