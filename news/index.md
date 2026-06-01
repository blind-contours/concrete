# Changelog

## concrete 1.1.0

### New features

- Added user-facing TMLE convergence controls through
  [`formatArguments()`](https://blind-contours.github.io/concrete/reference/formatArguments.md):
  `UpdateMethod`, `EICStopRule`, and `EICStopAbsTol`.
- Added hybrid and absolute empirical EIC stopping rules. The default
  remains the original relative rule.
- Added
  [`getTmleDiagnostics()`](https://blind-contours.github.io/concrete/reference/getTmleDiagnostics.md)
  for inspecting final component-wise empirical EIC diagnostics, update
  traces, and norm trajectories from fitted `ConcreteEst` objects.
- Added documented hazard learner aliases for Coxnet, random survival
  forests, additive hazards, and HAL.
- Added trialist-facing documentation, a pkgdown site, GitHub issue
  templates, and an installed smoke-test script for checking Cox-only
  and optional learner workflows.

### Stability

- Added trace fields for relative and absolute EIC ratios, maximum
  absolute empirical EIC, failing components, and stopping-rule
  metadata.
- Improved rare-event convergence behavior with the hybrid stopping rule
  while preserving the original default behavior.
- Restricted convergence stopping checks to the requested target
  event/time components, excluding internal complement rows used in EIC
  summaries.
- Aligned adaptive update acceptance with the active stopping rule:
  relative stopping uses the target EIC norm, while absolute and hybrid
  stopping use the component-wise stopping ratio.
- Added a zero-norm hazard-update guard and default smoke-test execution
  in the test suite.
