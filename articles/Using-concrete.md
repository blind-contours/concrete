# Using-concrete

## PBC competing risks data example

``` r

library(concrete)
library(data.table)
set.seed(12345)
data <- as.data.table(survival::pbc)
data <- data[!is.na(trt), ][, trt := trt - 1]
data <- data[, c("time", "status", "trt", "age", "sex", "albumin")]
```

## Baseline Intervention

``` r

ConcreteArgs <- formatArguments(DataTable = data,
                                EventTime = "time",
                                EventType = "status",
                                Treatment = "trt",
                                Intervention = 0:1,
                                TargetTime = 1500,
                                TargetEvent = 1:2,
                                MaxUpdateIter = 250,
                                CVArg = list(V = 10),
                                Verbose = FALSE)
ConcreteEst <- doConcrete(ConcreteArgs)
ConcreteOut <- getOutput(ConcreteEst)
```

## Convergence controls

``` r

ConcreteArgs$UpdateMethod <- "adaptive"
ConcreteArgs$EICStopRule <- "absolute"
ConcreteArgs$EICStopAbsTol <- 0.02 / sqrt(nrow(ConcreteArgs$Data))
ConcreteArgs <- formatArguments(ConcreteArgs)
ConcreteEst <- doConcrete(ConcreteArgs)
getTmleDiagnostics(ConcreteEst, type = "components")
```

## Survival learner library

For hazard learners, `concrete` currently uses cross-validated discrete
selection within each event-specific library.

``` r

Model <- list(
  trt = c("SL.glm", "SL.glmnet"),
  "0" = list(Cox = survival::Surv(time, status == 0) ~ .,
             RSF = "rsf",
             Aareg = "aareg"),
  "1" = list(Cox = survival::Surv(time, status == 1) ~ .,
             HAL = "hal")
)

ConcreteArgs <- formatArguments(DataTable = data,
                                EventTime = "time",
                                EventType = "status",
                                Treatment = "trt",
                                Intervention = 0:1,
                                TargetTime = 1500,
                                TargetEvent = 1:2,
                                MaxUpdateIter = 250,
                                CVArg = list(V = 10),
                                Model = Model,
                                Verbose = FALSE)
```

## Joint Intervention

``` r

data <- data[, trt2 := sample(0:1, .N, replace = TRUE, prob = c(0.3, .7))]
Intervention <- makeITT("A1" = data.frame(trt = rep_len(1, nrow(data)),
                                          trt2 = rep_len(1, nrow(data))),
                        "A0" = data.frame(trt = rep_len(0, nrow(data)),
                                          trt2 = rep_len(0, nrow(data))))

ConcreteArgs <- formatArguments(DataTable = data,
                                EventTime = "time",
                                EventType = "status",
                                Treatment = c("trt", "trt2"),
                                Intervention = Intervention,
                                TargetTime = 2000,
                                TargetEvent = 1:2,
                                MaxUpdateIter = 250,
                                CVArg = list(V = 10),
                                Verbose = FALSE)
ConcreteEst <- doConcrete(ConcreteArgs)
ConcreteOut <- getOutput(ConcreteEst)
```
