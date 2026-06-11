test_that(".strataAdjSigma2 matches the closed form and never exceeds iid", {
  set.seed(1)
  n <- 400
  S <- sample(c("a", "b", "c"), n, TRUE)
  A <- unlist(lapply(split(seq_len(n), S), function(ix) sample(rep(0:1, length.out = length(ix)))))
  A <- A[order(unlist(split(seq_len(n), S)))]
  IC <- rnorm(n) + 0.5 * (S == "a") * (2 * A - 1)   # stratum-arm structure
  IC <- IC - mean(IC)

  s2 <- concrete:::.strataAdjSigma2(IC, A, S)
  expect_true(is.numeric(s2) && s2 > 0)
  expect_lte(s2, mean(IC^2) + 1e-12)                # correction can only reduce

  ## equals iid - pi(1-pi) sum_s p_s Delta(s)^2 (up to the centering of m(s))
  pihat <- mean(A)
  ps <- prop.table(table(S))
  m1 <- tapply(IC[A == 1], S[A == 1], mean)[names(ps)]
  m0 <- tapply(IC[A == 0], S[A == 0], mean)[names(ps)]
  ms <- pihat * m1 + (1 - pihat) * m0
  direct <- mean((IC - ifelse(A == 1, m1[S], m0[S]))^2) +
    sum(ps * (ms - sum(ps * ms))^2)
  expect_equal(s2, as.numeric(direct), tolerance = 1e-10)
})

test_that(".strataAdjSigma2 returns NULL on degenerate strata", {
  IC <- rnorm(20); A <- rep(0:1, 10); S <- c(rep("a", 19), "b")  # stratum b: 1 subject
  expect_null(concrete:::.strataAdjSigma2(IC, A, S))
  expect_null(concrete:::.strataSE(IC, seq_len(20),
                                   data.table::data.table(ID = 1:19, A = A[1:19], S = S[1:19])))
})

test_that("Strata threads through to corrected SEs across estimands", {
  skip_on_cran()
  set.seed(2)
  n <- 400
  Z <- sample(1:3, n, TRUE)
  A <- integer(n)
  for (s in 1:3) {
    ix <- which(Z == s)
    for (b in split(ix, ceiling(seq_along(ix) / 4))) A[b] <- sample(rep(0:1, length.out = length(b)))
  }
  T1 <- rexp(n, 0.25 * exp(c(-0.7, 0, 0.7)[Z] - 0.3 * A))
  C <- rexp(n, 0.1)
  dat <- data.table::data.table(id = seq_len(n), time = pmin(T1, C, 4),
                                status = as.integer(T1 <= pmin(C, 4)),
                                trt = A, Z = factor(Z), W2 = rnorm(n))
  Mdl <- list(trt = "SL.mean",
              "0" = list(Cox = survival::Surv(time, status == 0) ~ trt),
              "1" = list(Cox = survival::Surv(time, status == 1) ~ trt))
  fit <- function(st) {
    a <- suppressMessages(formatArguments(DataTable = data.table::copy(dat),
          EventTime = "time", EventType = "status", Treatment = "trt", ID = "id",
          Intervention = 0:1, TargetTime = 3, TargetEvent = 1, CVArg = list(V = 2),
          MaxUpdateIter = 8, Model = Mdl, Verbose = FALSE, Strata = st))
    suppressMessages(suppressWarnings(doConcrete(a)))
  }
  e0 <- fit(NULL); e1 <- fit("Z")
  expect_null(attr(e0, "StrataDT"))
  expect_s3_class(attr(e1, "StrataDT"), "data.table")

  o0 <- suppressMessages(getOutput(e0, Estimand = c("Risk", "RD"), Simultaneous = FALSE))
  o1 <- suppressMessages(getOutput(e1, Estimand = c("Risk", "RD"), Simultaneous = FALSE))
  s0 <- as.data.frame(o0)[o0$Estimator == "tmle", c("Estimand", "se")]
  s1 <- as.data.frame(o1)[o1$Estimator == "tmle", c("Estimand", "se")]
  ## prognostic strata + unadjusted models: corrected SEs strictly smaller
  expect_true(all(s1$se <= s0$se + 1e-12))
  expect_lt(s1$se[s1$Estimand == "Risk Diff"], s0$se[s0$Estimand == "Risk Diff"])
  ## point estimates untouched
  expect_equal(as.data.frame(o0)[["Pt Est"]], as.data.frame(o1)[["Pt Est"]], tolerance = 1e-12)
})

test_that("missing baseline covariates are imputed with indicators", {
  set.seed(3)
  n <- 120
  dat <- data.table::data.table(
    id = seq_len(n), time = rexp(n, 0.3) + 0.01,
    status = rbinom(n, 1, 0.7), trt = rep(0:1, n / 2),
    W = rnorm(n), G = sample(c("x", "y"), n, TRUE))
  dat$W[1:10] <- NA
  dat$G[5:12] <- NA
  expect_message(
    a <- formatArguments(DataTable = data.table::copy(dat), EventTime = "time",
                         EventType = "status", Treatment = "trt", ID = "id",
                         Intervention = 0:1, TargetTime = 2, TargetEvent = 1,
                         CVArg = list(V = 2)),
    regexp = "Imputed|missingness indicator", ignore.case = TRUE)
  DT <- a$DataTable
  expect_false(anyNA(DT))
  cn <- attr(DT, "CovNames")
  expect_true(any(grepl("W_missing", cn$CovName)) || "W_missing" %in% colnames(DT))
  expect_true(any(grepl("G_missing", cn$CovName)) || "G_missing" %in% colnames(DT))
})

test_that("missing required columns still error", {
  n <- 50
  dat <- data.table::data.table(id = seq_len(n), time = rexp(n, 0.3) + 0.01,
                                status = rbinom(n, 1, 0.7), trt = rep(0:1, n / 2),
                                W = rnorm(n))
  dat$trt[3] <- NA
  expect_error(
    formatArguments(DataTable = dat, EventTime = "time", EventType = "status",
                    Treatment = "trt", ID = "id", Intervention = 0:1,
                    TargetTime = 2, TargetEvent = 1, CVArg = list(V = 2)),
    regexp = "complete|missing")
})
