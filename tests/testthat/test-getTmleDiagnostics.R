test_that("getTmleDiagnostics extracts component, trace, and norm diagnostics", {
    est <- list(
        "A=1" = list(
            SummEIC = data.table::data.table(
                Time = c(180, 365),
                Event = c(1, 1),
                PnEIC = c(7e-4, 2e-3),
                `seEIC/(sqrt(n)log(n))` = c(2e-5, 1e-3)
            )
        )
    )
    attr(est, "TargetTime") <- c(180, 365)
    attr(est, "TargetEvent") <- 1
    attr(est, "EICStopRule") <- "hybrid"
    attr(est, "EICStopAbsTol") <- 1e-3
    attr(est, "TmleConverged") <- list(converged = FALSE, step = 5)
    attr(est, "NormPnEICs") <- c(0.5, 0.25)
    attr(est, "TmleUpdateTrace") <- data.table::data.table(Step = 0, Status = "initial")
    class(est) <- "ConcreteEst"

    components <- getTmleDiagnostics(est, type = "components")
    expect_s3_class(components, "data.table")
    expect_true(components[Time == 180, check])
    expect_false(components[Time == 365, check])
    expect_equal(components[Time == 180, StopCriteria], 1e-3)
    expect_equal(components[Time == 365, ConvergenceStep], 5)
    expect_false(unique(components$Converged))

    trace <- getTmleDiagnostics(est, type = "trace")
    expect_equal(trace$Status, "initial")

    norm <- getTmleDiagnostics(est, type = "norm")
    expect_equal(norm$Step, c(0, 1))
    expect_equal(norm$NormPnEIC, c(0.5, 0.25))

    expect_error(getTmleDiagnostics(list()), regexp = "ConcreteEst")
})
