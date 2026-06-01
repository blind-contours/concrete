# Event-1 B100 N500 Simulation Notes

Run:

```sh
env LC_ALL=C LANG=C R_LIBS_USER=/private/tmp/concrete-lib \
  Rscript scripts/sim-data/referee-sims/run_pilot.R \
  --B=100 \
  --n=500 \
  --truth_n=500000 \
  --cv_v=2 \
  --target_times=180,365,730,1460 \
  --target_events=1 \
  --scenarios=ph_correct,nonph,rare_early \
  --include_survtmle=true \
  --survtmle_months=6 \
  --max_update_iter=200 \
  --cores=4 \
  --out=scripts/sim-data/referee-sims/output/event1_b100_n500
```

All 300 scenario-replicates completed. The row-status file has 300 `done`
rows and no worker errors.

Outputs:

- `scripts/sim-data/referee-sims/output/event1_b100_n500/summary.rds`
- `scripts/sim-data/referee-sims/output/event1_b100_n500/metrics.csv`
- `scripts/sim-data/referee-sims/output/event1_b100_n500/diagnostic_metrics.csv`
- `scripts/sim-data/referee-sims/output/event1_b100_n500/diagnostics_long.csv`
- `scripts/sim-data/referee-sims/output/event1_b100_n500/estimates_long.csv`

## Convergence

| Scenario | Estimator | Reps | Errors | Convergence | Median step | Median max ratio | Mean runtime sec |
|---|---:|---:|---:|---:|---:|---:|---:|
| nonph | concrete_standard_minimal | 100 | 0 | 1.00 | 14.0 | 0.958 | 18.7 |
| nonph | concrete_standard_rich | 100 | 0 | 1.00 | 14.5 | 0.956 | 16.3 |
| nonph | concrete_adaptive_rich | 100 | 0 | 1.00 | 12.0 | 0.953 | 20.5 |
| ph_correct | concrete_standard_minimal | 100 | 0 | 1.00 | 14.0 | 0.956 | 15.9 |
| ph_correct | concrete_standard_rich | 100 | 0 | 1.00 | 13.0 | 0.948 | 13.5 |
| ph_correct | concrete_adaptive_rich | 100 | 0 | 1.00 | 11.0 | 0.940 | 17.3 |
| rare_early | concrete_standard_minimal | 100 | 0 | 0.96 | 25.0 | 0.980 | 27.2 |
| rare_early | concrete_standard_rich | 100 | 0 | 0.96 | 25.5 | 0.974 | 24.9 |
| rare_early | concrete_adaptive_rich | 100 | 0 | 0.96 | 23.0 | 0.979 | 29.7 |

Non-convergence occurred only in `rare_early`, in 4 of 100 simulated datasets.
For each failed dataset all three concrete variants failed one target component
at `MaxUpdateIter = 200`.

Failed rare-early seeds:

| Rep | Seed | Max ratio range across concrete variants |
|---:|---:|---:|
| 3 | 1394360859 | 9.81 to 21.49 |
| 31 | 1449804648 | 18.94 to 20.18 |
| 78 | 1756082870 | 17.86 to 23.92 |
| 85 | 1775030832 | 21.23 to 36.46 |

## High-Iteration Sensitivity

For the representative failed dataset `rare_early`, rep 3, standard-rich was
rerun with `MaxUpdateIter = 500`.

| Estimator | Converged | Step | Final norm | Max ratio | Failing components | Runtime sec |
|---|---:|---:|---:|---:|---:|---:|
| standard_rich_max500 | FALSE | 501 | 0.000290 | 4.25 | 1 | 189.8 |

Increasing the iteration cap materially reduced the residual violation but did
not solve convergence for this rare early target. This supports guidance that
more iterations can help, but persistent failure at rare/early target times is
a practical warning rather than just a software failure.

## Performance Averaged Across Target Times

The table below averages across target times `180`, `365`, `730`, and `1460`.

| Scenario | Estimand | Estimator | Mean abs bias | Mean RMSE | Mean coverage |
|---|---|---:|---:|---:|---:|
| nonph | Risk A=1 | Aalen-Johansen | 0.0093 | 0.0213 | 0.88 |
| nonph | Risk A=1 | concrete_adaptive_rich | 0.0277 | 0.0356 | 0.70 |
| nonph | Risk A=1 | concrete_standard_minimal | 0.0277 | 0.0356 | 0.70 |
| nonph | Risk A=1 | concrete_standard_rich | 0.0273 | 0.0355 | 0.69 |
| nonph | Risk A=1 | survtmle_6mo | 0.0170 | 0.0262 | 0.81 |
| nonph | Risk Diff | Aalen-Johansen | 0.0228 | 0.0380 | 0.88 |
| nonph | Risk Diff | concrete_adaptive_rich | 0.0153 | 0.0291 | 0.96 |
| nonph | Risk Diff | concrete_standard_minimal | 0.0153 | 0.0290 | 0.96 |
| nonph | Risk Diff | concrete_standard_rich | 0.0154 | 0.0291 | 0.96 |
| nonph | Risk Diff | survtmle_6mo | 0.0227 | 0.0375 | 0.88 |
| ph_correct | Risk A=1 | Aalen-Johansen | 0.0070 | 0.0198 | 0.93 |
| ph_correct | Risk A=1 | concrete_adaptive_rich | 0.0133 | 0.0240 | 0.92 |
| ph_correct | Risk A=1 | concrete_standard_minimal | 0.0138 | 0.0242 | 0.92 |
| ph_correct | Risk A=1 | concrete_standard_rich | 0.0135 | 0.0248 | 0.92 |
| ph_correct | Risk A=1 | survtmle_6mo | 0.0030 | 0.0155 | 1.00 |
| ph_correct | Risk Diff | Aalen-Johansen | 0.0160 | 0.0328 | 0.94 |
| ph_correct | Risk Diff | concrete_adaptive_rich | 0.0036 | 0.0241 | 1.00 |
| ph_correct | Risk Diff | concrete_standard_minimal | 0.0038 | 0.0242 | 1.00 |
| ph_correct | Risk Diff | concrete_standard_rich | 0.0038 | 0.0242 | 1.00 |
| ph_correct | Risk Diff | survtmle_6mo | 0.0051 | 0.0208 | 1.00 |
| rare_early | Risk A=1 | Aalen-Johansen | 0.0049 | 0.0151 | 0.88 |
| rare_early | Risk A=1 | concrete_adaptive_rich | 0.0059 | 0.0168 | 0.91 |
| rare_early | Risk A=1 | concrete_standard_minimal | 0.0048 | 0.0164 | 0.91 |
| rare_early | Risk A=1 | concrete_standard_rich | 0.0053 | 0.0167 | 0.91 |
| rare_early | Risk A=1 | survtmle_6mo | 0.0077 | 0.0140 | 0.95 |
| rare_early | Risk Diff | Aalen-Johansen | 0.0118 | 0.0247 | 0.93 |
| rare_early | Risk Diff | concrete_adaptive_rich | 0.0070 | 0.0183 | 0.99 |
| rare_early | Risk Diff | concrete_standard_minimal | 0.0071 | 0.0185 | 0.99 |
| rare_early | Risk Diff | concrete_standard_rich | 0.0068 | 0.0182 | 0.99 |
| rare_early | Risk Diff | survtmle_6mo | 0.0027 | 0.0151 | 0.99 |

## Read

- Convergence is reliable in the proportional-hazards and non-proportional
  hazards scenarios at this target dimension.
- Rare early targets cause the practical convergence issue reviewers asked
  about: 4% non-convergence at `MaxUpdateIter = 200`, always one component.
- Adaptive reduces median accepted update steps, but is slower in wall time.
  It is a useful sensitivity method, not a default runtime improvement.
- Cox-only hazards show a real limitation under non-proportional hazards for
  treatment-specific absolute risks. Risk differences remain well behaved in
  this simulation, which is useful nuance for the response.
- For the manuscript, report convergence/runtime separately from estimator
  performance and use the rare-early scenario as practical guidance on when
  iteration increases are not enough.
