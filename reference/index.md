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

## Lower-level helpers

- [`getInitialEstimate()`](https://blind-contours.github.io/concrete/reference/getInitialEstimate.md)
  : getInitialEstimate
- [`getPropScore()`](https://blind-contours.github.io/concrete/reference/getPropScore.md)
  : getPropScore
- [`getHazFit()`](https://blind-contours.github.io/concrete/reference/getHazFit.md)
  : Title
- [`getEIC()`](https://blind-contours.github.io/concrete/reference/getEIC.md)
  : get EICs
- [`doTmleUpdate()`](https://blind-contours.github.io/concrete/reference/doTmleUpdate.md)
  : TMLE update with multiple implementation options

## Developer and internal helpers

- [`concrete-package`](https://blind-contours.github.io/concrete/reference/concrete-package.md)
  [`concrete`](https://blind-contours.github.io/concrete/reference/concrete-package.md)
  : One-step continuous-time Targeted Minimum Loss-Based Estimator
  (TMLE) for outcome-specific absolute risk estimands in right-censored
  survival settings with or without competing risks
- [`printOneStepDiagnostics()`](https://blind-contours.github.io/concrete/reference/printOneStepDiagnostics.md)
  : Print diagnostic information about TMLE convergence
- [`deepcopyEstimates()`](https://blind-contours.github.io/concrete/reference/deepcopyEstimates.md)
  : Deep-copy an Estimates list so we can do parallel updates safely
  (Jacobi).
- [`submodelAndEIC()`](https://blind-contours.github.io/concrete/reference/submodelAndEIC.md)
  : submodelAndEIC
- [`applyFinalSubmodel()`](https://blind-contours.github.io/concrete/reference/applyFinalSubmodel.md)
  : applyFinalSubmodel
- [`updateHazardsWithEps()`](https://blind-contours.github.io/concrete/reference/updateHazardsWithEps.md)
  : updateHazardsWithEps
- [`updateHazard()`](https://blind-contours.github.io/concrete/reference/updateHazard.md)
  : Update hazards based on clever covariate and PnEIC
