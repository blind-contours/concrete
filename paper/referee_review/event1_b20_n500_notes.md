# Event-1 B20 N500 Simulation Notes

Run:

```sh
env LC_ALL=C LANG=C R_LIBS_USER=/private/tmp/concrete-lib \
  Rscript scripts/sim-data/referee-sims/run_pilot.R \
  --B=20 \
  --n=500 \
  --truth_n=300000 \
  --cv_v=2 \
  --target_times=180,365,730,1460 \
  --target_events=1 \
  --scenarios=ph_correct,nonph,rare_early \
  --include_survtmle=true \
  --survtmle_months=6 \
  --cores=4 \
  --out=scripts/sim-data/referee-sims/output/event1_b20_n500
```

All 60 scenario-replicates completed.

Outputs:

- `scripts/sim-data/referee-sims/output/event1_b20_n500/summary.rds`
- `scripts/sim-data/referee-sims/output/event1_b20_n500/metrics.csv`
- `scripts/sim-data/referee-sims/output/event1_b20_n500/diagnostic_metrics.csv`
- `scripts/sim-data/referee-sims/output/event1_b20_n500/diagnostics_long.csv`
- `scripts/sim-data/referee-sims/output/event1_b20_n500/estimates_long.csv`

## Convergence

| Scenario | Estimator | Reps | Errors | Convergence | Median step | Median max ratio | Mean runtime sec |
|---|---:|---:|---:|---:|---:|---:|---:|
| nonph | concrete_standard_minimal | 20 | 0 | 1.00 | 19.0 | 0.956 | 22.8 |
| nonph | concrete_standard_rich | 20 | 0 | 1.00 | 19.0 | 0.955 | 19.9 |
| nonph | concrete_adaptive_rich | 20 | 0 | 1.00 | 14.5 | 0.967 | 24.8 |
| ph_correct | concrete_standard_minimal | 20 | 0 | 1.00 | 11.0 | 0.950 | 16.5 |
| ph_correct | concrete_standard_rich | 20 | 0 | 1.00 | 12.0 | 0.909 | 13.9 |
| ph_correct | concrete_adaptive_rich | 20 | 0 | 1.00 | 11.0 | 0.899 | 17.7 |
| rare_early | concrete_standard_minimal | 20 | 0 | 0.95 | 22.5 | 0.984 | 31.6 |
| rare_early | concrete_standard_rich | 20 | 0 | 0.95 | 21.5 | 0.969 | 28.5 |
| rare_early | concrete_adaptive_rich | 20 | 0 | 0.95 | 19.0 | 0.979 | 35.4 |

The only non-converged case was `rare_early`, replicate 5. It failed one target
component for all three concrete variants at `MaxUpdateIter = 200`.

## Performance Averaged Across Target Times

The table below averages across the four target times. These are exploratory
because `B = 20` is still small.

| Scenario | Estimand | Estimator | Mean abs bias | Mean RMSE | Mean coverage |
|---|---|---:|---:|---:|---:|
| nonph | Risk A=1 | Aalen-Johansen | 0.0066 | 0.0231 | 0.86 |
| nonph | Risk A=1 | concrete_adaptive_rich | 0.0318 | 0.0414 | 0.56 |
| nonph | Risk A=1 | concrete_standard_minimal | 0.0319 | 0.0415 | 0.56 |
| nonph | Risk A=1 | concrete_standard_rich | 0.0312 | 0.0404 | 0.58 |
| nonph | Risk A=1 | survtmle_6mo | 0.0195 | 0.0294 | 0.75 |
| nonph | Risk Diff | Aalen-Johansen | 0.0217 | 0.0382 | 0.89 |
| nonph | Risk Diff | concrete_adaptive_rich | 0.0141 | 0.0279 | 0.95 |
| nonph | Risk Diff | concrete_standard_minimal | 0.0143 | 0.0280 | 0.95 |
| nonph | Risk Diff | concrete_standard_rich | 0.0142 | 0.0280 | 0.95 |
| nonph | Risk Diff | survtmle_6mo | 0.0242 | 0.0366 | 0.90 |
| ph_correct | Risk A=1 | Aalen-Johansen | 0.0049 | 0.0230 | 0.91 |
| ph_correct | Risk A=1 | concrete_adaptive_rich | 0.0155 | 0.0294 | 0.84 |
| ph_correct | Risk A=1 | concrete_standard_minimal | 0.0151 | 0.0285 | 0.86 |
| ph_correct | Risk A=1 | concrete_standard_rich | 0.0149 | 0.0288 | 0.86 |
| ph_correct | Risk A=1 | survtmle_6mo | 0.0036 | 0.0182 | 0.99 |
| ph_correct | Risk Diff | Aalen-Johansen | 0.0117 | 0.0335 | 0.92 |
| ph_correct | Risk Diff | concrete_adaptive_rich | 0.0059 | 0.0282 | 0.95 |
| ph_correct | Risk Diff | concrete_standard_minimal | 0.0064 | 0.0287 | 0.95 |
| ph_correct | Risk Diff | concrete_standard_rich | 0.0065 | 0.0287 | 0.95 |
| ph_correct | Risk Diff | survtmle_6mo | 0.0019 | 0.0222 | 0.98 |
| rare_early | Risk A=1 | Aalen-Johansen | 0.0039 | 0.0131 | 0.88 |
| rare_early | Risk A=1 | concrete_adaptive_rich | 0.0085 | 0.0169 | 0.94 |
| rare_early | Risk A=1 | concrete_standard_minimal | 0.0071 | 0.0167 | 0.94 |
| rare_early | Risk A=1 | concrete_standard_rich | 0.0079 | 0.0170 | 0.92 |
| rare_early | Risk A=1 | survtmle_6mo | 0.0085 | 0.0129 | 0.99 |
| rare_early | Risk Diff | Aalen-Johansen | 0.0068 | 0.0243 | 0.95 |
| rare_early | Risk Diff | concrete_adaptive_rich | 0.0028 | 0.0208 | 0.99 |
| rare_early | Risk Diff | concrete_standard_minimal | 0.0034 | 0.0197 | 0.99 |
| rare_early | Risk Diff | concrete_standard_rich | 0.0030 | 0.0202 | 0.99 |
| rare_early | Risk Diff | survtmle_6mo | 0.0026 | 0.0170 | 1.00 |

## Read

- The convergence issue appears targeted: only the rare-early stress scenario
  had non-convergence, and only one of 20 replicates.
- Adaptive reduced median accepted steps in `nonph` and `rare_early`, but it
  remained slower because line search adds work.
- Under non-proportional hazards, concrete absolute risk estimates for `A=1`
  show noticeable positive bias and low pointwise coverage, while risk
  differences remain better behaved. This directly supports the reviewer-facing
  caveat that Cox-only hazard learners can matter.
- For the next evidence-building run, increase `B` for event 1 before widening
  to joint events. Joint event-1/event-2 targeting should be a separate
  convergence-stress supplement.
