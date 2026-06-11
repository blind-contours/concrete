# Variance under stratified / covariate-adaptive randomization

Influence-function standard errors assume simple (iid) randomization.
Under covariate-adaptive schemes that achieve strong balance –
stratified permuted blocks, stratified biased coin, minimization (under
conditions) – the iid variance is generically conservative: it includes
a between-arm-within-stratum component that the randomization scheme
removes (Bugni, Canay & Shaikh 2018, JASA; Ye, Yi & Shao 2022,
Biometrika; Ye, Shao, Yi & Zhao 2023, JASA).

For an asymptotically linear estimator with influence function psi,
target allocation pi, strata S with frequencies p_s, and \\\Delta(s) =
E(\psi \mid A=1, S=s) - E(\psi \mid A=0, S=s)\\: \$\$\sigma^2\_{CAR} =
\sigma^2\_{iid} - \pi(1-\pi)\sum_s p_s \Delta(s)^2.\$\$ The estimator
below computes the equivalent residual + between-stratum form
\$\$\hat\sigma^2 = P_n\[(\psi - \bar\psi\_{A,S})^2\] + \sum_s \hat p_s
(\hat m(s) - \bar m)^2,\$\$ with \\\hat m(s) = \hat\pi \bar\psi\_{1,s} +
(1-\hat\pi)\bar\psi\_{0,s}\\, which is nonnegative by construction.
Derivation: notes/strata-variance.md.

Because every reported estimand (absolute risk, RD, RR, RMST/LYL, win
ratio) stores per-subject influence values, supplying `Strata` to
[`formatArguments()`](https://blind-contours.github.io/concrete/reference/formatArguments.md)
corrects the standard errors of all of them through this one helper.
When the working models adjust for the stratification variables
(recommended; the strata columns stay in the data as covariates),
Delta(s) is approximately 0 and the correction is approximately 0 – the
iid variance is then already asymptotically correct.

## Usage

``` r
.strataAdjSigma2(IC, A, S)
```

## Arguments

- IC:

  numeric: per-subject influence values (mean approximately 0).

- A:

  numeric/integer: per-subject randomized arm (binary, \> 0 = treated).

- S:

  per-subject randomization-stratum labels.

## Value

the adjusted variance of the influence function (so `se = sqrt(. / n)`),
or `NULL` when a stratum-arm cell has fewer than 2 subjects (caller
falls back to the iid variance).
