## Covariate-adjusted crossover (treatment-switching) IPCW in the win-ratio paths.

xoFit <- function(seed = 4, n = 600, tau = 4) {
  set.seed(seed); A <- rep(0:1, each = n); W <- rnorm(2 * n)
  tD <- rexp(2 * n, 0.15 * exp(0.3 * W - 0.5 * A))
  tH <- rexp(2 * n, 0.30 * exp(0.2 * W - 0.3 * A)); C <- rexp(2 * n, 0.05)
  ## OMT subjects switch to the device (informative: sicker switch sooner)
  sw <- ifelse(A == 0 & runif(2 * n) < plogis(-0.5 + 0.8 * W), runif(2 * n, 0.5, tau), NA)
  obst <- pmin(tD, C, tau)
  data.frame(arm = A, t_hosp = ifelse(tH < obst, tH, NA), t_term = obst,
             died = as.integer(tD <= pmin(C, tau)), W = W, W2 = rnorm(2 * n), xover = sw)
}

test_that("clinicalWinRatio accepts a crossover column and runs with finite inference", {
  skip_on_cran()
  dat <- xoFit()
  o <- suppressMessages(suppressWarnings(clinicalWinRatio(
    dat, arm = "arm", illness.time = "t_hosp", terminal.time = "t_term",
    terminal.status = "died", covariates = c("W", "W2"), horizon = 4, n.grid = 20,
    n.folds = 1, crossover = "xover")))
  wr <- as.data.frame(o)[o$Estimand == "Win Ratio", ]
  expect_true(is.finite(wr$`Pt Est`) && wr$`Pt Est` > 0 && is.finite(wr$se) && wr$se > 0)
})

test_that("crossover IPCW changes the estimand vs ITT (no-switching != treatment-policy)", {
  skip_on_cran()
  dat <- xoFit(seed = 9)
  ar <- list(data = dat, arm = "arm", illness.time = "t_hosp", terminal.time = "t_term",
             terminal.status = "died", covariates = c("W", "W2"), horizon = 4, n.grid = 20, n.folds = 1)
  itt <- suppressMessages(suppressWarnings(do.call(clinicalWinRatio, ar)))
  nos <- suppressMessages(suppressWarnings(do.call(clinicalWinRatio, c(ar, list(crossover = "xover")))))
  g <- function(o) as.data.frame(o)[o$Estimand == "Win Ratio", "Pt Est"]
  expect_false(isTRUE(all.equal(g(itt), g(nos))))          # the adjustment moves the estimate
})

test_that("min.cens.surv truncation floor is honored (looser floor -> larger weights)", {
  skip_on_cran()
  dat <- xoFit(seed = 2)
  ar <- list(data = dat, arm = "arm", illness.time = "t_hosp", terminal.time = "t_term",
             terminal.status = "died", covariates = c("W", "W2"), horizon = 4, n.grid = 20,
             n.folds = 1, crossover = "xover")
  o1 <- suppressMessages(suppressWarnings(do.call(clinicalWinRatio, c(ar, list(min.cens.surv = 0.05)))))
  o2 <- suppressMessages(suppressWarnings(do.call(clinicalWinRatio, c(ar, list(min.cens.surv = 0.01)))))
  g <- function(o) as.data.frame(o)[o$Estimand == "Win Ratio", "Pt Est"]
  expect_true(is.finite(g(o1)) && is.finite(g(o2)))
  expect_false(isTRUE(all.equal(g(o1), g(o2))))            # the floor is actually used
})

test_that("an invalid crossover column errors clearly", {
  expect_error(suppressWarnings(clinicalWinRatio(
    xoFit(n = 80), arm = "arm", illness.time = "t_hosp", terminal.time = "t_term",
    terminal.status = "died", covariates = c("W", "W2"), horizon = 4, n.grid = 10,
    n.folds = 1, crossover = "not_a_column")), regexp = "crossover column")
})
