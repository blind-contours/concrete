test_that("clinicalRMTIF point estimate runs and decomposes", {
  skip_on_cran()
  set.seed(11)
  n <- 600; A <- rep(0:1, each = n); W <- rnorm(2*n)
  lD <- 0.10 * exp(0.3*W - 0.5*A); lH <- 0.20 * exp(0.2*W - 0.3*A)
  tH <- rexp(2*n, lH); tD <- rexp(2*n, lD); C <- rexp(2*n, 0.05)
  dat <- data.frame(arm = A, t_hosp = ifelse(tH <= pmin(tD, C, 5), tH, NA),
                    t_term = pmin(tD, C, 5), died = as.integer(tD <= pmin(C, 5)),
                    W = W, W2 = rnorm(2*n))
  o <- suppressMessages(suppressWarnings(clinicalRMTIF(dat, arm = "arm",
        illness.time = "t_hosp", terminal.time = "t_term", terminal.status = "died",
        covariates = c("W", "W2"), horizon = 4, n.grid = 30, n.folds = 1)))
  expect_s3_class(o, "ConcreteOut")
  expect_setequal(o$Estimand, c("RMT-IF", "Time in favor", "Time against"))
  net <- o[Estimand == "RMT-IF", `Pt Est`]
  fav <- o[Estimand == "Time in favor", `Pt Est`]
  agn <- o[Estimand == "Time against", `Pt Est`]
  expect_equal(net, fav - agn, tolerance = 1e-8)
  expect_true(fav > 0 && agn > 0)                 # both arms spend time in each role
  expect_true(is.na(o[Estimand == "RMT-IF", se]))  # no SE without nBoot
  expect_identical(attr(o, "Tiers"), 2L)
})

test_that("clinicalRMTIF validates inputs", {
  dat <- data.frame(arm = rep(0:1, 5), t_hosp = NA, t_term = runif(10),
                    died = rep(0:1, 5), W = rnorm(10))
  expect_error(clinicalRMTIF(dat, arm = "arm", illness.time = "t_hosp",
    terminal.time = "t_term", terminal.status = "died", covariates = "nope"),
    regexp = "not found")
})
