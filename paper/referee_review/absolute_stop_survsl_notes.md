# Absolute stopping rule survival-SL comparison

Date: 2026-06-01

## Objective

Evaluate whether a sample-size-scaled absolute empirical EIC stopping rule is a
better practical fallback than the original relative rule for rare-event or
rich-learner convergence problems.

Primary comparison:

- `relative_adaptive_e0.01_min0.05`: `UpdateMethod = "adaptive"`,
  `EICStopRule = "relative"`
- `absolute_nsqrt_0.02`: `UpdateMethod = "adaptive"`,
  `EICStopRule = "absolute"`, `EICStopAbsTol = 0.02 / sqrt(n)`

Hazard library: full survival learner stack from the referee simulation helper:
Cox formulas, Coxnet, random survival forest, additive hazards, and HAL.

## Main run

Output directory:

`scripts/sim-data/referee-sims/output/alt_convergence_primary_survsl_all_n150`

Design:

- 8 hard seeds: non-PH, positivity, and rare-early cases
- `n = 150`
- target event 1
- target times 365 and 1095
- `MaxUpdateIter = 120`
- treatment library `SL.glm`

Summary:

| Config | Jobs | Convergence | Median step | Max step | Median max abs PnEIC | Max abs PnEIC |
|---|---:|---:|---:|---:|---:|---:|
| absolute_nsqrt_0.02 | 8 | 8/8 | 12.0 | 65 | 0.00154 | 0.00163 |
| relative_adaptive_e0.01_min0.05 | 8 | 7/8 | 7.5 | 121 | 0.00318 | 0.00809 |

The relative rule failed on `rare_early_rep10` at the iteration cap:

| Case | Rule | Step | Max ratio | Max abs PnEIC | Failing components |
|---|---|---:|---:|---:|---:|
| rare_early_rep10 | relative | 121 | 11.05 | 0.00106 | 1 |

For the same seed, the absolute rule converged because the largest absolute
empirical EIC was still on the risk scale:

| Rule | Component | Abs PnEIC | Relative ratio | Absolute ratio | Passed |
|---|---|---:|---:|---:|---|
| absolute | A=1, 365d, event 1 | 0.00163 | 16.05 | 0.996 | TRUE |
| relative | A=1, 365d, event 1 | 0.00106 | 11.05 | Inf | FALSE |

Estimate differences versus the relative rule were small on absolute-risk and
risk-difference scales:

| Estimand class | Max absolute difference | Median absolute difference | 90th percentile |
|---|---:|---:|---:|
| Risk/RD | 0.00298 | 0.00056 | 0.00149 |
| RR | 0.07660 | 0.00649 | 0.04712 |

RR differences were larger because the event risks were small; this supports
reporting absolute-risk and risk-difference sensitivity results alongside risk
ratios in rare-event analyses.

## Interpretation

The absolute n-scaled rule is the best documented fallback from these runs. It
directly controls empirical EIC imbalance on the absolute-risk scale and avoids
spending many update iterations on components whose relative threshold is tiny
because the EIC variance is near zero.

Recommended user guidance:

1. Start with the default relative rule and `UpdateMethod = "adaptive"`.
2. If convergence fails with tiny `AbsPnEIC` but large relative ratios, run a
   sensitivity with `EICStopRule = "absolute"` and
   `EICStopAbsTol = 0.02 / sqrt(nrow(data))`.
3. Report the stopping rule and diagnostics from `getTmleDiagnostics()`.
4. Compare absolute risks and risk differences first; interpret rare-event risk
   ratios with caution.
