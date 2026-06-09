# Time-varying covariates in the censoring model (shared internal machinery).

Post-randomization measurements (e.g.\\ echo / KCCQ / 6-minute walk at
follow-up visits) can drive dropout, making censoring informative
conditional on baseline covariates alone. Conditioning the censoring
hazard on the time-varying covariate history restores conditional
independence (CAR) and removes the
inverse-probability-of-censoring-weight bias. These covariates are used
**only** in the censoring model — never the outcome hazards — so the
intent-to-treat / marginal estimand is preserved (they are
post-treatment mediators).

## Details

`.tvLOCF()` builds last-observation-carried-forward value and
change-from-baseline matrices on the hazard grid; `.tvCensoringInc()`
fits the (cross-fitted) discrete-time censoring hazard with those
covariates and returns per-interval cumulative-hazard increments. Both
the core IPCW (via an override of the lagged censoring survival) and
[`clinicalWinRatio()`](https://blind-contours.github.io/concrete/reference/clinicalWinRatio.md)
consume these.
