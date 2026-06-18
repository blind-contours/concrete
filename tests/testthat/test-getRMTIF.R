test_that("getRMTIF reduces exactly to the RMST difference when K = 1", {
  skip_on_cran()
  set.seed(7)
  n <- 400; W <- rnorm(2*n); A <- rep(0:1, each = n)
  T1 <- rexp(2*n, 0.18 * exp(0.3*W - 0.5*A)); C <- rexp(2*n, 0.08)
  dat <- data.table::data.table(id = 1:(2*n), time = pmin(T1, C, 5),
    status = as.integer(T1 <= pmin(C, 5)), trt = A, W = W, W2 = rnorm(2*n))
  a <- suppressMessages(formatArguments(DataTable = dat, EventTime = "time",
    EventType = "status", Treatment = "trt", ID = "id", Intervention = 0:1,
    TargetTime = c(1, 2, 3, 4), TargetEvent = 1, CVArg = list(V = 2),
    MaxUpdateIter = 8, Verbose = FALSE,
    Model = list(trt = "SL.mean", "0" = list(Cox = survival::Surv(time, status == 0) ~ .),
                 "1" = list(Cox = survival::Surv(time, status == 1) ~ .))))
  e <- suppressMessages(suppressWarnings(doConcrete(a)))
  rmst  <- suppressMessages(getRMST(e, Horizon = 4, Intervention = c(1, 2)))
  rmtif <- suppressMessages(getRMTIF(e, Horizon = 4, Intervention = c(1, 2), TargetEvent = 1))
  rd <- as.data.frame(rmst)[as.data.frame(rmst)$Estimand == "RMST Diff", ]
  ri <- as.data.frame(rmtif)[as.data.frame(rmtif)$Estimand == "RMT-IF", ]
  expect_equal(ri[["Pt Est"]], rd[["Pt Est"]], tolerance = 1e-6)
  expect_equal(ri[["se"]],     rd[["se"]],     tolerance = 1e-6)
})

test_that("getRMTIF structure, decomposition, and family attachment", {
  skip_on_cran()
  set.seed(8)
  n <- 350; W <- rnorm(2*n); A <- rep(0:1, each = n)
  l1 <- 0.14*exp(0.3*W - 0.5*A); l2 <- 0.16*exp(0.2*W - 0.3*A)
  Te <- rexp(2*n, l1+l2); cause <- 1L + rbinom(2*n, 1, l2/(l1+l2)); C <- rexp(2*n, 0.08)
  dat <- data.table::data.table(id = 1:(2*n), time = pmin(Te, C, 5),
    status = ifelse(Te <= pmin(C, 5), cause, 0L), trt = A, W = W, W2 = rnorm(2*n))
  a <- suppressMessages(formatArguments(DataTable = dat, EventTime = "time",
    EventType = "status", Treatment = "trt", ID = "id", Intervention = 0:1,
    TargetTime = c(1, 2, 3, 4), TargetEvent = c(1, 2), CVArg = list(V = 2),
    MaxUpdateIter = 8, Verbose = FALSE,
    Model = list(trt = "SL.mean", "0" = list(Cox = survival::Surv(time, status == 0) ~ .),
                 "1" = list(Cox = survival::Surv(time, status == 1) ~ .),
                 "2" = list(Cox = survival::Surv(time, status == 2) ~ .))))
  e <- suppressMessages(suppressWarnings(doConcrete(a)))
  o <- suppressMessages(getRMTIF(e, Horizon = 4, Intervention = c(1, 2), TargetEvent = c(1, 2)))
  expect_setequal(o$Estimand, c("RMT-IF", "Time in favor", "Time against"))
  ## net = favor - against, and all SEs positive
  net <- o[Estimand == "RMT-IF", `Pt Est`]
  fav <- o[Estimand == "Time in favor", `Pt Est`]
  agn <- o[Estimand == "Time against", `Pt Est`]
  expect_equal(net, fav - agn, tolerance = 1e-8)
  expect_true(all(o$se > 0))
  ## feeds getSimultaneousFamily
  expect_s3_class(attr(o, "famEst"), "data.table")
  fam <- getSimultaneousFamily(RMTIF = o)
  expect_true("RMT-IF" %in% fam$Estimand)
})

test_that("getRMTIF validates inputs", {
  expect_error(getRMTIF(list()), regexp = "ConcreteEst")
})
