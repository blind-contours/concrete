# Event-1 Pilot Notes

Run:

```sh
R_LIBS_USER=/private/tmp/concrete-lib Rscript scripts/sim-data/referee-sims/run_pilot.R \
  --B=2 \
  --n=300 \
  --truth_n=30000 \
  --cv_v=2 \
  --target_times=180,365,730,1460 \
  --target_events=1 \
  --include_survtmle=true \
  --survtmle_months=6 \
  --out=scripts/sim-data/referee-sims/output/pilot_event1
```

Outputs:

- `scripts/sim-data/referee-sims/output/pilot_event1/summary.rds`
- `scripts/sim-data/referee-sims/output/pilot_event1/metrics.csv`
- `scripts/sim-data/referee-sims/output/pilot_event1/diagnostic_metrics.csv`
- `scripts/sim-data/referee-sims/output/pilot_event1/estimates_long.csv`

This was a smoke pilot only; `B=2` is not enough for inferential conclusions.

## Convergence

All concrete estimators completed without errors and converged in all six
scenario-replicates.

| Scenario | Estimator | Median step | Median max ratio | Mean runtime sec |
|---|---:|---:|---:|---:|
| nonph | concrete_standard_minimal | 9.0 | 0.841 | 3.65 |
| nonph | concrete_standard_rich | 14.5 | 0.957 | 4.55 |
| nonph | concrete_adaptive_rich | 11.5 | 0.978 | 5.21 |
| ph_correct | concrete_standard_minimal | 17.5 | 0.941 | 4.64 |
| ph_correct | concrete_standard_rich | 18.0 | 0.976 | 4.89 |
| ph_correct | concrete_adaptive_rich | 17.0 | 0.936 | 6.34 |
| rare_early | concrete_standard_minimal | 14.0 | 0.770 | 3.61 |
| rare_early | concrete_standard_rich | 15.5 | 0.772 | 4.05 |
| rare_early | concrete_adaptive_rich | 14.0 | 0.784 | 5.30 |

## Immediate Read

- The simulation harness works for `concrete`, Aalen-Johansen, and 6-month
  `survtmle`.
- `adaptive` now returns complete output after rollback; the previous
  `GCompEst`-stripping bug is fixed.
- For this small event-1 pilot, adaptive sometimes reduces accepted steps, but
  line-search overhead means runtime is not lower.
- The heavier joint event-1/event-2 pilot is feasible but much slower. One
  non-PH replicate with five target times and two events took about two minutes.

## Next Run

Use event 1 only for the first meaningful run:

- `B = 100`
- `n = 500`
- `truth_n = 500000`
- scenarios: `ph_correct`, `nonph`, `rare_early`
- target times: `180,365,730,1460`
- estimators: concrete standard-minimal, standard-rich, adaptive-rich,
  Aalen-Johansen, `survtmle_6mo`

After that run, decide whether to add event 2 jointly or to reserve the
multi-event target for a smaller convergence-focused supplement.
