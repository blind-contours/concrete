# formatArguments

formatArguments() checks and reformats inputs into a form that can be
interpreted by doConcrete(). makeITT() returns an Intervention list for
a single, binary, point-treatment variable

## Usage

``` r
formatArguments(
  DataTable,
  EventTime,
  EventType,
  Treatment,
  ID = NULL,
  TargetTime = NULL,
  TargetEvent = NULL,
  Intervention,
  CVArg = NULL,
  Model = NULL,
  MaxUpdateIter = 500,
  OneStepEps = 0.1,
  MinNuisance = 5/sqrt(nrow(DataTable))/log(nrow(DataTable)),
  Verbose = TRUE,
  GComp = TRUE,
  ReturnModels = TRUE,
  ConcreteArgs = NULL,
  RenameCovs = TRUE,
  UpdateMethod = c("standard", "adaptive"),
  EICStopRule = c("relative", "absolute", "hybrid"),
  EICStopAbsTol = 0,
  CrossFit = FALSE,
  HazEnsemble = FALSE,
  ...
)

makeITT(...)

# S3 method for class 'ConcreteArgs'
print(x, ...)
```

## Arguments

- DataTable:

  data.table (n x (d + (3:5)); data.table of the observed data, with
  rows n = the number of observations and d = the number of baseline
  covariates. DataTable must include the following columns:

  - "EventTime": numeric; real numbers \> 0, the observed event or
    censoring time

  - "EventType": numeric; the observed event type, censoring events
    indicated by integers \<= 0

  - "Treatment": numeric; the observed treatment value. Binary
    treatments must be coded as 0, 1

  - "Treatment": numeric; the observed treatment

  May include

  - "ID": factor, character, or numeric; unique subject id. If ID column
    is missing, row numbers will be used as ID. For longitudinal data,
    ID must be provided

  - "Baseline Covariates": factor, character, or numeric;

- EventTime:

  character: the column name of the observed event or censoring time

- EventType:

  character: the column name of the observed event type. (0 indicating
  censoring)

- Treatment:

  character: the column name of the observed treatment assignment

- ID:

  character (default: NULL): the column name of the observed subject id
  longitudinal data structures

- TargetTime:

  numeric: vector of target times. If NULL, the last observed
  non-censoring event time will be targeted.

- TargetEvent:

  numeric: vector of target events - some subset of unique EventTypes.
  If NULL, all non-censoring observed event types will be targeted.

- Intervention:

  list: a list of desired interventions on the treatment variable. Each
  intervention must be a list containing two named functions:
  'intervention' = function(treatment vector, covariate data) and
  'gstar' = function(treatment vector, covariate data)
  concrete::makeITT() can be used to specify an intent-to-treat analysis
  for a binary intervention variable

- CVArg:

  list: arguments to be passed into do.call(origami::make_folds). If
  NULL, the default is list(n = nrow(DataTable), fold_fun = folds_vfold,
  cluster_ids = NULL, strata_ids = NULL)

- Model:

  list (default: NULL): named list of models, one for each failure or
  censoring event and one for the `Treatment` variable. Treatment models
  are passed to `SuperLearner`. Hazard models can be Cox formulas or
  one-word learner aliases. Supported hazard aliases are `"coxnet"`,
  `"rsf"`/`"randomForestSRC"`, `"aareg"`/`"additive_hazards"`, and
  `"hal"`/`"hal9001"`. These optional learners require `glmnet`,
  `randomForestSRC`, or `hal9001` when selected. If `Model = NULL`, a
  Cox model template is generated for the user to amend.

- MaxUpdateIter:

  numeric (default: 500): the number of one-step update steps

- OneStepEps:

  numeric (default: 0.1): the one-step TMLE step size

- MinNuisance:

  numeric (default: 5/log(n)/sqrt(n)): value between (0, 1) for
  truncating the g-related denominator of the clever covariate

- Verbose:

  boolean

- GComp:

  boolean (default: TRUE): return g-computation formula plug-in
  estimates

- ReturnModels:

  boolean (default: TRUE): return fitted models from the initial
  estimation stage

- ConcreteArgs:

  list (default: NULL, not yet ready) : Use to recheck amended output
  from previous formatArguments() calls. A non-NULL input will cause all
  other arguments to be ignored.

- RenameCovs:

  boolean (default: TRUE): whether or not to rename covariates

- UpdateMethod:

  character (default: "standard"): TMLE update method. Supported values
  are `"standard"` and `"adaptive"`. `"adaptive"` uses a line search
  with rollback and is usually the most stable choice for difficult
  convergence cases. `"accelerated"` is accepted as a legacy alias for
  `"adaptive"`.

- EICStopRule:

  character (default: "relative"): TMLE stopping rule. Supported values
  are `"relative"`, `"absolute"`, and `"hybrid"`. `"relative"` preserves
  the original criterion, `|PnEIC| <= seEIC / (sqrt(n) log(n))`.
  `"absolute"` checks `|PnEIC| <= EICStopAbsTol`. `"hybrid"` checks
  `|PnEIC| <= max(seEIC / (sqrt(n) log(n)), EICStopAbsTol)`. The
  absolute rule is useful when rare-event or near-zero-variance
  components make the relative criterion numerically too strict.

- EICStopAbsTol:

  numeric (default: 0): absolute \|PnEIC\| tolerance used when
  `EICStopRule` is `"absolute"` or `"hybrid"`. A value such as
  `0.02 / sqrt(nrow(DataTable))` gives a small sample-size-scaled
  absolute-risk tolerance while the default `0` leaves the original
  relative rule unchanged.

- CrossFit:

  logical (default: FALSE): if TRUE, estimate the propensity and hazard
  nuisances by cross-fitting (CV-TMLE) – each subject's nuisances are
  predicted from models fit on the other folds, which supports valid
  influence-function inference when flexible machine-learning learners
  are used. Adds compute (the nuisance library is refit once per fold).

- HazEnsemble:

  logical (default: FALSE): if TRUE, combine the candidate hazard
  learners into a cross-validated convex-combination ensemble (Super
  Learner) by minimizing the counting-process negative log-likelihood of
  the weighted hazard, instead of the default discrete (winner-take-all)
  selection. The treatment propensity already uses an ensemble Super
  Learner.

- ...:

  additional arguments to be passed into print methods

- x:

  a ConcreteArgs object

## Value

a list of class "ConcreteArgs"

- Data: data.table containing EventTime, EventType, Treatment, and
  potentially ID and baseline covariates. Has the following attributes

  - EventTime: the column name of the observed event or censoring time

  - EventType: the column name of the observed event type. (0 indicating
    censoring)

  - Treatment: the column name of the observed treatment assignment

  - ID: the column name of the observed subject id

  - RenameCovs: boolean whether or not covariates are renamed

- TargetTime: numeric vector of target times to evaluate risk/survival

- TargetEvent: numeric vector of target events

- Regime: named list of desired regimes, each tagged with a 'g.star'
  attribute function

  - Regime\[\[i\]\]: a vector of desired treatment assignments

  - attr(Regime\[\[i\]\], "g.star"): function of Treatment and
    Covariates, outputting a vector of desired treatment assignment
    probabilities

- CVFolds: list of cross-validation fold assignments in the structure as
  output by origami::make_folds()

- Model: named list of model specifications, one for each unique
  'EventType' and one for the 'Treatment' variable.

- MaxUpdateIter: the number of one-step update steps

- OneStepEps: initial one-step TMLE step size

- MinNuisance: numeric lower bound for the propensity score denominator
  in the efficient influence function

- Verbose: boolean to print additional information

- GComp: boolean to return g-computation formula plug-in estimates

- ReturnModels: boolean to return fitted models from the initial
  estimation stage

- EICStopRule: TMLE stopping rule for empirical mean EIC checks

- EICStopAbsTol: absolute empirical mean EIC tolerance for absolute or
  hybrid stopping

## Functions

- `makeITT()`: makeITT ...

- `print(ConcreteArgs)`: print.ConcreteArgs print method for
  "ConcreteArgs" class

## Examples

``` r
library(data.table)
library(concrete)

data <- as.data.table(survival::pbc)
data <- data[1:200, .SD, .SDcols = c("id", "time", "status", "trt", "age", "sex")]
data[, trt := sample(0:1, nrow(data), TRUE)]
#>         id  time status   trt      age    sex
#>      <int> <int>  <int> <int>    <num> <fctr>
#>   1:     1   400      2     0 58.76523      f
#>   2:     2  4500      0     0 56.44627      f
#>   3:     3  1012      2     0 70.07255      m
#>   4:     4  1925      2     1 54.74059      f
#>   5:     5  1504      1     1 38.10541      f
#>  ---                                         
#> 196:   196  2363      0     1 57.04038      f
#> 197:   197  2365      0     0 44.62697      f
#> 198:   198  2357      0     1 35.79740      f
#> 199:   199  1592      0     1 40.71732      f
#> 200:   200  2318      0     0 32.23272      f

# makeITT() creates a list of functions to specify intent-to-treat
#   regimes for a binary, single, point treatment variable
intervention <- makeITT()

# formatArguments() returns correctly formatted arguments for doConcrete()
#   If no input is provided for the Model argument, a default will be generated
concrete.args <- formatArguments(DataTable = data,
                                 EventTime = "time",
                                 EventType = "status",
                                 Treatment = "trt",
                                 ID = "id",
                                 TargetTime = 2500,
                                 TargetEvent = c(1, 2),
                                 Intervention = intervention,
                                 CVArg = list(V = 2))

# Alternatively, estimation algorithms can be provided as a named list
model <- list("trt" = c("SL.glm", "SL.glmnet"),
              "0" = list(Surv(time, status == 0) ~ .),
              "1" = list(Surv(time, status == 1) ~ .),
              "2" = list(Surv(time, status == 2) ~ .))
concrete.args <- formatArguments(DataTable = data,
                                 EventTime = "time",
                                 EventType = "status",
                                 Treatment = "trt",
                                 ID = "id",
                                 TargetTime = 2500,
                                 TargetEvent = c(1, 2),
                                 Intervention = intervention,
                                 CVArg = list(V = 2),
                                 Model = model)

# 'ConcreteArgs' output can be modified and passed back through formatArguments()
# examples of modifying the censoring and failure event candidate regressions
concrete.args[["Model"]][["0"]] <-
    list(Surv(time, status == 0) ~ trt:sex + age)
concrete.args[["Model"]][["1"]] <-
    list("mod1" = Surv(time, status == 1) ~ trt,
         "mod2" = Surv(time, status == 1) ~ .)
formatArguments(concrete.args)
#> 
#> Observed Data (200 rows x 6 cols)
#> Unique IDs: "id" (n=200),  Time-to-Event: "time",  Event Type: "status",  Treatment: "trt"
#> 
#> - - - - - - - - - - - - - - - - - - - - 
#> Estimand Specification:
#> Target Events: 1, 2
#> 
#> Target Time (n at risk): 2500 (100/200)
#> 
#> Interventions
#>   A=1: ("trt" = [1,1,1,1,1,1,1,1,1,1,...])  -  Observed Prevalence = 0.52
#>   A=0: ("trt" = [0,0,0,0,0,0,0,0,0,0,...])  -  Observed Prevalence = 0.48
#> 
#> - - - - - - - - - - - - - - - - - - - - 
#> Estimation Specification:
#> Stratified 2-Fold Cross Validation 
#> "trt" Propensity Score Estimation (SuperLearner): Default SL Selector, Default Loss Fn, 2 candidates - SL.glm, SL.glmnet
#> Cens. 0 Estimation (coxph): Discrete SL Selector, Log Partial-LL Loss, 1 candidate - model1
#> Event 1 Estimation (coxph): Discrete SL Selector, Log Partial-LL Loss, 2 candidates - mod1, mod2
#> Event 2 Estimation (coxph): Discrete SL Selector, Log Partial-LL Loss, 1 candidate - model1
#> 
#> One-step TMLE (finite sum approx.) simultaneously targeting all cause-specific Absolute Risks
#> g nuisance bounds = [0.06673, 1],  max update steps = 500,  starting one-step epsilon = 0.1,  EIC stop rule = relative
#> 
#> ****
#> Cox model specifications have been renamed where necessary to reflect changed covariate names. Model specifications in .[['Model']] can be checked against the covariate names in attr(.[['DataTable']], 'CovNames')
#> ****

# For difficult rare-event settings, use an adaptive update with an
# absolute risk-scale stopping rule:
# concrete.args$UpdateMethod <- "adaptive"
# concrete.args$EICStopRule <- "absolute"
# concrete.args$EICStopAbsTol <- 0.02 / sqrt(nrow(data))
# concrete.args <- formatArguments(concrete.args)
```
