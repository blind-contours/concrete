test_that("clinicalWinRatio runs on multistate data and returns coherent win statistics", {
  skip_on_cran()
  skip_if_not_installed("SuperLearner")
  set.seed(1); n <- 400; tau <- 1500
  W <- stats::rnorm(n); arm <- stats::rbinom(n, 1, 0.5)
  a01 <- 6e-4 * exp(-0.4 * arm + 0.3 * W); a02 <- 4e-4 * exp(0.2 * W); a12 <- 1e-3 * exp(0.3 * W)
  T01 <- stats::rexp(n, a01); T02 <- stats::rexp(n, a02); u0 <- pmin(T01, T02); hfh <- T01 < T02
  s <- ifelse(hfh, u0, Inf); dpost <- ifelse(hfh, u0 + stats::rexp(n, a12), Inf)
  D <- ifelse(hfh, dpost, u0)
  C <- stats::rexp(n, 3e-4); term <- pmin(D, C, tau); died <- as.integer(D < C & D < tau)
  t_hfh <- ifelse(hfh & s < term, s, NA)
  d <- data.frame(arm = arm, t_hfh = t_hfh, t_term = term, died = died, W = W)

  r <- clinicalWinRatio(d, arm = "arm", illness.time = "t_hfh", terminal.time = "t_term",
                        terminal.status = "died", covariates = "W", horizon = tau, n.grid = 25)

  expect_s3_class(r, "ConcreteOut")
  expect_true(all(c("Win Ratio", "Win Odds", "Net Benefit", "P(win)", "P(loss)", "P(tie)") %in% r$Estimand))
  # win/loss/tie form a valid distribution
  probs <- r[r$Estimand %in% c("P(win)", "P(loss)", "P(tie)"), `Pt Est`]
  expect_equal(sum(probs), 1, tolerance = 1e-6)
  expect_true(all(probs >= -1e-8 & probs <= 1 + 1e-8))
  # comparative statistics have finite positive SEs
  comp <- r[r$Estimand %in% c("Win Ratio", "Win Odds", "Net Benefit"), ]
  expect_true(all(is.finite(comp$se) & comp$se > 0))
  # internal consistency: Win Ratio = P(win)/P(loss)
  pw <- r$`Pt Est`[r$Estimand == "P(win)"]; pl <- r$`Pt Est`[r$Estimand == "P(loss)"]
  expect_equal(r$`Pt Est`[r$Estimand == "Win Ratio"], pw / pl, tolerance = 1e-6)
})

test_that("clinicalWinRatio validates its inputs", {
  d <- data.frame(arm = 0:1, t_hfh = c(NA, 1), t_term = c(2, 3), died = c(1, 0), W = c(0.1, -0.2))
  expect_error(clinicalWinRatio(d, arm = "nope", illness.time = "t_hfh",
                                terminal.time = "t_term", terminal.status = "died", covariates = "W"))
})

test_that("clinicalWinRatio handles a 3-tier hierarchy (death > E2 > E3)", {
  skip_on_cran()
  skip_if_not_installed("SuperLearner")
  set.seed(2); n <- 500; tau <- 3
  W <- stats::rnorm(n); arm <- stats::rbinom(n, 1, 0.5)
  fW <- exp(0.3 * W); fT <- exp(-0.4 * arm); fE <- exp(-0.3 * arm)
  # latent first-event times for two non-fatal events + death (illustrative)
  tE2 <- stats::rexp(n, 0.12 * fW * fE)   # higher-priority non-fatal
  tE3 <- stats::rexp(n, 0.20 * fW * fE)   # lower-priority non-fatal
  tD  <- stats::rexp(n, 0.06 * fW * fT)
  C   <- stats::rexp(n, 0.12)
  termT <- pmin(tD, C, tau); died <- as.integer(tD <= C & tD <= tau)
  ill <- function(tt) ifelse(tt <= termT, tt, NA)
  d <- data.frame(arm = arm, t_e2 = ill(tE2), t_e3 = ill(tE3),
                  t_term = termT, died = ifelse(termT >= tau, 0L, died), W = W)

  r <- clinicalWinRatio(d, arm = "arm", illness.time = c("t_e2", "t_e3"),
                        terminal.time = "t_term", terminal.status = "died",
                        covariates = "W", horizon = tau, n.grid = 20, n.folds = 1)

  expect_s3_class(r, "ConcreteOut")
  expect_identical(attr(r, "Tiers"), 3L)
  probs <- r[r$Estimand %in% c("P(win)", "P(loss)", "P(tie)"), `Pt Est`]
  expect_equal(sum(probs), 1, tolerance = 1e-6)
  expect_true(all(probs >= -1e-8 & probs <= 1 + 1e-8))
  comp <- r[r$Estimand %in% c("Win Ratio", "Win Odds", "Net Benefit"), ]
  expect_true(all(is.finite(comp$se) & comp$se > 0))
  pw <- r$`Pt Est`[r$Estimand == "P(win)"]; pl <- r$`Pt Est`[r$Estimand == "P(loss)"]
  expect_equal(r$`Pt Est`[r$Estimand == "Win Ratio"], pw / pl, tolerance = 1e-6)
})
