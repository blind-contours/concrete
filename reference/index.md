# Package index

## Main workflow

- [`formatArguments()`](https://blind-contours.github.io/concrete/reference/formatArguments.md)
  [`makeITT()`](https://blind-contours.github.io/concrete/reference/formatArguments.md)
  [`print(`*`<ConcreteArgs>`*`)`](https://blind-contours.github.io/concrete/reference/formatArguments.md)
  : formatArguments
- [`doConcrete()`](https://blind-contours.github.io/concrete/reference/doConcrete.md)
  [`print(`*`<ConcreteEst>`*`)`](https://blind-contours.github.io/concrete/reference/doConcrete.md)
  [`plot(`*`<ConcreteEst>`*`)`](https://blind-contours.github.io/concrete/reference/doConcrete.md)
  [`print(`*`<ConcreteOut>`*`)`](https://blind-contours.github.io/concrete/reference/doConcrete.md)
  : doConcrete
- [`getOutput()`](https://blind-contours.github.io/concrete/reference/getOutput.md)
  [`plot(`*`<ConcreteOut>`*`)`](https://blind-contours.github.io/concrete/reference/getOutput.md)
  : getOutput
- [`getTmleDiagnostics()`](https://blind-contours.github.io/concrete/reference/getTmleDiagnostics.md)
  : Extract TMLE convergence diagnostics

## RCT estimands

- [`getRMST()`](https://blind-contours.github.io/concrete/reference/getRMST.md)
  : Restricted mean survival time and cause-specific life-years lost
- [`targetRMST()`](https://blind-contours.github.io/concrete/reference/targetRMST.md)
  : Directly target the restricted mean survival time
- [`getRelativeEfficiency()`](https://blind-contours.github.io/concrete/reference/getRelativeEfficiency.md)
  : Relative efficiency of a covariate-adjusted vs unadjusted analysis
- [`getWinRatio()`](https://blind-contours.github.io/concrete/reference/getWinRatio.md)
  : Covariate-adjusted restricted win ratio, win odds, and net benefit
- [`clinicalWinRatio()`](https://blind-contours.github.io/concrete/reference/clinicalWinRatio.md)
  : Clinical (death-priority) win ratio for an illness-death outcome
  (experimental)

## Estimand framework

- [`makeEstimand()`](https://blind-contours.github.io/concrete/reference/makeEstimand.md)
  [`print(`*`<ConcreteEstimand>`*`)`](https://blind-contours.github.io/concrete/reference/makeEstimand.md)
  [`applyIntercurrentEvent()`](https://blind-contours.github.io/concrete/reference/makeEstimand.md)
  : Specify an ICH E9(R1) estimand and apply an intercurrent-event
  strategy
- [`senseCensoring()`](https://blind-contours.github.io/concrete/reference/senseCensoring.md)
  : Tipping-point sensitivity analysis for informative censoring

## Package overview

- [`concrete-package`](https://blind-contours.github.io/concrete/reference/concrete-package.md)
  [`concrete`](https://blind-contours.github.io/concrete/reference/concrete-package.md)
  : One-step continuous-time Targeted Minimum Loss-Based Estimator
  (TMLE) for outcome-specific absolute risk estimands in right-censored
  survival settings with or without competing risks
