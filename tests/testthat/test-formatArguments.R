set.seed(0)
data <- data.table::as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
data[, trt := sample(0:1, length(trt), replace = TRUE)]

test_that("formatArguments works ", {
    concrete.args <- formatArguments(Data = data,
                                     EventTime = "time",
                                     EventType = "status",
                                     Treatment = "trt",
                                     ID = 'id',
                                     Intervention = 0:1,
                                     TargetTime = quantile(data[["time"]], probs = seq(.1, .9, .05)),
                                     TargetEvent = setdiff(unique(data[["status"]]), 0)
    )
    expect_s3_class(concrete.args, class = "ConcreteArgs")
    expect_s3_class(formatArguments(concrete.args), class = "ConcreteArgs")
    expect_s3_class(formatArguments(ConcreteArgs = concrete.args), class = "ConcreteArgs")
    expect_error(formatArguments(ConcreteArgs = data))
})

test_that("Data with missingness or incorrect type throw errors", {
    require(data.table)
    DataWithMissing <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
    expect_error(formatArguments(Data = DataWithMissing,
                                 EventTime = "time",
                                 EventType = "status",
                                 Treatment = "trt",
                                 ID = 'id',
                                 Intervention = 0:1,
                                 TargetTime = mean(data[["time"]]),
                                 TargetEvent = unique(data[["status"]])))
    expect_error(formatArguments(Data = as.data.frame(DataWithMissing),
                                 EventTime = "time",
                                 EventType = "status",
                                 Treatment = "trt",
                                 ID = 'id',
                                 Intervention = 0:1,
                                 TargetTime = mean(data[["time"]]),
                                 TargetEvent = unique(data[["status"]])))
    expect_error(formatArguments(Data = as.numeric(data$time),
                                 EventTime = "time",
                                 EventType = "status",
                                 Treatment = "trt",
                                 ID = 'id',
                                 Intervention = 0:1,
                                 TargetTime = mean(data[["time"]]),
                                 TargetEvent = unique(data[["status"]])))
    expect_error(formatArguments(Data = "foo",
                                 EventTime = "time",
                                 EventType = "status",
                                 Treatment = "trt",
                                 ID = 'id',
                                 Intervention = 0:1,
                                 TargetTime = mean(data[["time"]]),
                                 TargetEvent = unique(data[["status"]])))
})

test_that("EventTime is a positive, finite numeric vector", {
    test_vals <- list(NaN, NA, Inf, TRUE, "a", 0, -1)
    for (value in test_vals) {
        expect_error(concrete:::checkEventTime(value, data.frame("x" = value)))
        expect_error(concrete:::checkEventTime("x", data.frame("x" = value)))
    }
})

test_that("EventType is a non-negative numeric vector", {
    test_vals <- list(NaN, NA, TRUE, "a", -1)
    for (value in test_vals) {
        expect_error(concrete:::checkEventType(value, data.frame("x" = value)))
        expect_error(concrete:::checkEventType("x", data.frame("x" = value)))
    }
})

test_that("Treatment is a numeric vector", {
    test_vals <- list(NaN, NA, Inf, TRUE, "a")
    for (value in test_vals) {
        expect_error(concrete:::checkTreatment(value, data.frame("x" = value)))
        expect_error(concrete:::checkTreatment("x", data.frame("x" = value)))
    }
})

test_that("Intervention specifications", {
    test_vals <- list(NaN, NA, Inf, "a", matrix(1, 3, 3),
                      function(...) return(list(...)),
                      list(function(x) x, function(y) 1))
    for (value in test_vals) {
        expect_error(formatArguments(DataTable = data, EventTime = "time",
                                     EventType = "status", Treatment = "trt",
                                     ID = "id", Intervention = value))
    }
})

test_that("ID is a vector with non-\'null\'-type values", {
    require(data.table)
    data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]

    set.seed(0)
    data[, trt := sample(0:1, length(trt), replace = TRUE)]

    expect_error(concrete:::getID(NULL, data), regexp = NA)

    test_vals <- list(NaN, NA)
    for (value in test_vals) {
        data <- data.frame("x" = value)
        expect_error(concrete:::getID(value, data))
        expect_error(concrete:::getID("x", data))
    }
})

test_that("Boolean cheecks for non-boolean values and resets values to FALSE", {
    require(data.table)
    data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]

    set.seed(0)
    data[, trt := sample(0:1, length(trt), replace = TRUE)]
    concrete.args <- formatArguments(Data = data,
                                     EventTime = "time",
                                     EventType = "status",
                                     Treatment = "trt",
                                     ID = 'id',
                                     Intervention = 0:1,
                                     TargetTime = mean(data[["time"]]),
                                     TargetEvent = setdiff(unique(data[["status"]]), 0),
                                     Verbose = 2,
                                     GComp = NA,
                                     ReturnModels = "c",
                                     RenameCovs = Inf)
    for (bool in c("Verbose", "GComp", "ReturnModels", "RenameCovs")) {
        expect_equal(concrete.args[[bool]], FALSE)
    }
})

test_that("RenameCovs = FALSE gets processed correctly", {
    require(data.table)
    data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]

    set.seed(0)
    data[, trt := sample(0:1, length(trt), replace = TRUE)]
    concrete.args <- formatArguments(Data = data,
                                     EventTime = "time",
                                     EventType = "status",
                                     Treatment = "trt",
                                     ID = 'id',
                                     Intervention = 0:1,
                                     TargetTime = mean(data[["time"]]),
                                     TargetEvent = setdiff(unique(data[["status"]]), 0),
                                     RenameCovs = FALSE)
    expect_equal(colnames(concrete.args$DataTable), colnames(data))
})

test_that("TargetEvent preserves requested non-censoring subset", {
    concrete.args <- formatArguments(Data = data,
                                     EventTime = "time",
                                     EventType = "status",
                                     Treatment = "trt",
                                     ID = "id",
                                     Intervention = 0:1,
                                     TargetTime = 2500,
                                     TargetEvent = 1,
                                     CVArg = list(V = 2),
                                     Verbose = FALSE)
    expect_true(identical(concrete.args$TargetEvent, 1))

    expect_error(formatArguments(Data = data,
                                 EventTime = "time",
                                 EventType = "status",
                                 Treatment = "trt",
                                 ID = "id",
                                 Intervention = 0:1,
                                 TargetTime = 2500,
                                 TargetEvent = 0,
                                 CVArg = list(V = 2),
                                 Verbose = FALSE),
                 regexp = "non-censoring")
})

test_that("UpdateMethod is scalar and simulation-safe", {
    concrete.args <- formatArguments(Data = data,
                                     EventTime = "time",
                                     EventType = "status",
                                     Treatment = "trt",
                                     ID = "id",
                                     Intervention = 0:1,
                                     TargetTime = 2500,
                                     TargetEvent = 1,
                                     CVArg = list(V = 2),
                                     Verbose = FALSE)
    expect_true(identical(concrete.args$UpdateMethod, "standard"))

    concrete.args <- formatArguments(Data = data,
                                     EventTime = "time",
                                     EventType = "status",
                                     Treatment = "trt",
                                     ID = "id",
                                     Intervention = 0:1,
                                     TargetTime = 2500,
                                     TargetEvent = 1,
                                     CVArg = list(V = 2),
                                     Verbose = FALSE,
                                     UpdateMethod = "accelerated")
    expect_true(identical(concrete.args$UpdateMethod, "adaptive"))

    concrete.args <- formatArguments(Data = data,
                                     EventTime = "time",
                                     EventType = "status",
                                     Treatment = "trt",
                                     ID = "id",
                                     Intervention = 0:1,
                                     TargetTime = 2500,
                                     TargetEvent = 1,
                                     CVArg = list(V = 2),
                                     Verbose = FALSE,
                                     UpdateMethod = " Adaptive ")
    expect_true(identical(concrete.args$UpdateMethod, "adaptive"))

    expect_error(formatArguments(Data = data,
                                 EventTime = "time",
                                 EventType = "status",
                                 Treatment = "trt",
                                 ID = "id",
                                 Intervention = 0:1,
                                 TargetTime = 2500,
                                 TargetEvent = 1,
                                 CVArg = list(V = 2),
                                 Verbose = FALSE,
                                 UpdateMethod = "coordinated"),
                 regexp = "removed")

    expect_error(formatArguments(Data = data,
                                 EventTime = "time",
                                 EventType = "status",
                                 Treatment = "trt",
                                 ID = "id",
                                 Intervention = 0:1,
                                 TargetTime = 2500,
                                 TargetEvent = 1,
                                 CVArg = list(V = 2),
                                 Verbose = FALSE,
                                 UpdateMethod = "rootSolve"),
                 regexp = "removed")
})

test_that("EIC stopping rules are parsed and evaluated", {
    concrete.args <- formatArguments(Data = data,
                                     EventTime = "time",
                                     EventType = "status",
                                     Treatment = "trt",
                                     ID = "id",
                                     Intervention = 0:1,
                                     TargetTime = 2500,
                                     TargetEvent = 1,
                                     CVArg = list(V = 2),
                                     Verbose = FALSE,
                                     EICStopRule = "hybrid",
                                     EICStopAbsTol = 1e-3)
    expect_true(identical(concrete.args$EICStopRule, "hybrid"))
    expect_equal(concrete.args$EICStopAbsTol, 1e-3)

    expect_error(formatArguments(Data = data,
                                 EventTime = "time",
                                 EventType = "status",
                                 Treatment = "trt",
                                 ID = "id",
                                 Intervention = 0:1,
                                 TargetTime = 2500,
                                 TargetEvent = 1,
                                 CVArg = list(V = 2),
                                 Verbose = FALSE,
                                 EICStopRule = "loose"),
                 regexp = "EICStopRule")

    stop_dt <- data.table::data.table(
        Trt = "A=1",
        Time = 180,
        Event = 1,
        PnEIC = 7e-4,
        `seEIC/(sqrt(n)log(n))` = 2e-5
    )
    expect_false(concrete:::makeOneStepStop(stop_dt, "relative", 0)$check)
    expect_true(concrete:::makeOneStepStop(stop_dt, "HYBRID", 1e-3)$check)
})

test_that("Target convergence helpers ignore internal complement rows", {
    stop_dt <- data.table::data.table(
        Trt = c("A=1", "A=1"),
        Time = c(180, 180),
        Event = c(1, -1),
        PnEIC = c(7e-4, 1),
        `seEIC/(sqrt(n)log(n))` = c(1e-3, 1e-3)
    )

    expect_false(all(concrete:::makeOneStepStop(stop_dt, "relative", 0)$check))
    target_stop <- concrete:::targetOneStepStop(
        stop_dt,
        TargetTime = 180,
        TargetEvent = 1,
        EICStopRule = "relative",
        EICStopAbsTol = 0
    )
    expect_equal(nrow(target_stop), 1)
    expect_true(all(target_stop$check))
})

test_that("Hazard update is a no-op when target norm is numerically zero", {
    haz <- matrix(0.1, nrow = 2, ncol = 2)
    attr(haz, "j") <- 1
    hazards <- list("1" = haz)
    out <- concrete:::updateHazard(
        GStar = c(1, 1),
        Hazards = hazards,
        TotalSurv = matrix(0.9, nrow = 2, ncol = 2),
        NuisanceWeight = matrix(1, nrow = 2, ncol = 2),
        EvalTimes = c(0, 1),
        T.tilde = c(1, 1),
        Delta = c(1, 1),
        PnEIC = data.table::data.table(Time = 1, Event = 1, PnEIC = 0),
        NormPnEIC = 0,
        OneStepEps = 0.1,
        TargetEvent = 1,
        TargetTime = 1
    )

    expect_equal(out[["1"]], haz)
    expect_equal(attr(out[["1"]], "j"), 1)
})

test_that("Survival hazard learner aliases are parsed", {
    model <- list(
        trt = c("SL.glm"),
        "0" = list(RSF = "rsf", Aareg = "aareg", HAL = "hal", Coxnet = "coxnet"),
        "1" = list(RSF = "randomForestSRC", Aareg = "additive_hazards", HAL = "hal9001")
    )
    concrete.args <- formatArguments(Data = data,
                                     EventTime = "time",
                                     EventType = "status",
                                     Treatment = "trt",
                                     ID = "id",
                                     Intervention = 0:1,
                                     TargetTime = 2500,
                                     TargetEvent = 1,
                                     CVArg = list(V = 2),
                                     Model = model,
                                     Verbose = FALSE)
    expect_s3_class(concrete.args$Model[["0"]][["RSF"]], "Lrnr.RSF")
    expect_s3_class(concrete.args$Model[["0"]][["Aareg"]], "Lrnr.Aareg")
    expect_s3_class(concrete.args$Model[["0"]][["HAL"]], "Lrnr.HAL")
    expect_s3_class(concrete.args$Model[["0"]][["Coxnet"]], "Lrnr.Coxnet")
    expect_s3_class(concrete.args$Model[["1"]][["RSF"]], "Lrnr.RSF")
    expect_s3_class(concrete.args$Model[["1"]][["Aareg"]], "Lrnr.Aareg")
    expect_s3_class(concrete.args$Model[["1"]][["HAL"]], "Lrnr.HAL")
})
