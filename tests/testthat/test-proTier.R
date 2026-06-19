## Continuous / ordinal PRO bottom tiers for the hierarchical win statistics.

## PRO tiers use a final-visit (landmark = horizon) design; the marker is observed
## among reachers (event-free, alive, in follow-up at the horizon).
proFit <- function(seed = 5, n = 500, tau = 4) {
  set.seed(seed); A <- rep(0:1, each = n); W <- rnorm(2 * n)
  tD <- rexp(2 * n, 0.12 * exp(0.3 * W - 0.6 * A))
  tH <- rexp(2 * n, 0.25 * exp(0.2 * W - 0.2 * A)); C <- rexp(2 * n, 0.04)
  obst <- pmin(tD, C, tau); reach <- tD > tau & C >= tau & tH > tau   # at the horizon
  kccq <- ifelse(reach & runif(2 * n) < plogis(0.9 + 0.4 * W), 58 + 6 * W + 7 * A + rnorm(2 * n, 0, 14), NA)
  nyha <- ifelse(reach & runif(2 * n) < plogis(0.9 + 0.4 * W),
                 pmin(4, pmax(1, round(2.5 - 0.6 * A + 0.4 * W + rnorm(2 * n)))), NA)
  data.frame(arm = A, t_hosp = ifelse(tH < obst, tH, NA), t_term = obst,
             died = as.integer(tD <= pmin(C, tau)), W = W, W2 = rnorm(2 * n),
             kccq = kccq, nyha = nyha)
}
kccqSpec <- list(marker = "kccq", landmark = 4, margin = 5, direction = "higher.better",
                 type = "continuous", label = "KCCQ")

test_that("a continuous PRO bottom tier runs in clinicalWinRatio with finite inference", {
  skip_on_cran()
  o <- suppressMessages(suppressWarnings(clinicalWinRatio(
    proFit(), arm = "arm", illness.time = "t_hosp", terminal.time = "t_term",
    terminal.status = "died", covariates = c("W", "W2"), horizon = 4, n.grid = 20,
    n.folds = 1, pro = kccqSpec)))
  expect_s3_class(o, "ConcreteOut")
  expect_equal(attr(o, "Tiers"), 3L)                         # death > hosp > KCCQ
  wr <- as.data.frame(o)[o$Estimand == "Win Ratio", ]
  expect_true(is.finite(wr$`Pt Est`) && wr$`Pt Est` > 0 && is.finite(wr$se) && wr$se > 0)
})

test_that("PRO tier appears in clinicalPSNB and charter length must count it", {
  skip_on_cran()
  ar <- list(data = proFit(), arm = "arm", illness.time = "t_hosp", terminal.time = "t_term",
             terminal.status = "died", covariates = c("W", "W2"), horizon = 4, n.grid = 20,
             n.folds = 1, pro = kccqSpec)
  o <- suppressMessages(suppressWarnings(do.call(clinicalPSNB, c(ar, list(charter = c(0.5, 0.3, 0.2))))))
  expect_true(any(grepl("KCCQ", o$Estimand)))                # PRO tier reported
  expect_equal(o[Estimand == "Reach[D]", `Pt Est`], 1)
  rk <- o[grepl("^Reach", Estimand), `Pt Est`]
  expect_true(all(rk[is.finite(rk)] <= 1 + 1e-6))            # reach is a probability
  ## a 2-vector charter (ignoring the PRO tier) must error: K = 3
  expect_error(suppressWarnings(do.call(clinicalPSNB, c(ar, list(charter = c(0.5, 0.5))))),
               regexp = "length K")
})

test_that("charter='reach' with a PRO tier reproduces sum reach_k * NB_k", {
  skip_on_cran()
  o <- suppressMessages(suppressWarnings(clinicalPSNB(
    proFit(seed = 8), arm = "arm", illness.time = "t_hosp", terminal.time = "t_term",
    terminal.status = "died", covariates = c("W", "W2"), charter = "reach",
    horizon = 4, n.grid = 20, n.folds = 1, pro = kccqSpec)))
  rk <- o[grepl("^Reach", Estimand), `Pt Est`]
  dk <- o[grepl("^NetBenefit", Estimand), `Pt Est`]
  expect_equal(o[Estimand == "PSNB", `Pt Est`], sum(rk * dk), tolerance = 1e-8)
})

test_that("an ordinal (lower-is-better) PRO tier runs", {
  skip_on_cran()
  nyhaSpec <- list(marker = "nyha", landmark = 4, margin = 0, direction = "lower.better",
                   type = "ordinal", label = "NYHA")
  o <- suppressMessages(suppressWarnings(clinicalPSNB(
    proFit(seed = 3), arm = "arm", illness.time = "t_hosp", terminal.time = "t_term",
    terminal.status = "died", covariates = c("W", "W2"), charter = c(0.5, 0.3, 0.2),
    horizon = 4, n.grid = 20, n.folds = 1, pro = nyhaSpec)))
  expect_true(any(grepl("NYHA", o$Estimand)))
  expect_true(is.finite(o[Estimand == "PSNB", `Pt Est`]))
})

test_that("stacked PRO tiers are compared sequentially (win+loss+tie sums to 1)", {
  skip_on_cran()
  nyhaSpec <- list(marker = "nyha", landmark = 4, margin = 0, direction = "higher.better",
                   type = "ordinal", label = "NYHA")
  o <- suppressMessages(suppressWarnings(clinicalWinRatio(
    proFit(seed = 6), arm = "arm", illness.time = "t_hosp", terminal.time = "t_term",
    terminal.status = "died", covariates = c("W", "W2"), horizon = 4, n.grid = 20,
    n.folds = 1, pro = list(kccqSpec, nyhaSpec))))      # death > hosp > KCCQ > NYHA
  expect_equal(attr(o, "Tiers"), 4L)
  d <- as.data.frame(o)
  pw <- d[d$Estimand == "P(win)", "Pt Est"]; pl <- d[d$Estimand == "P(loss)", "Pt Est"]
  pt <- d[d$Estimand == "P(tie)", "Pt Est"]
  expect_true(pw >= 0 && pw <= 1 && pl >= 0 && pl <= 1)   # proper probabilities (no double-count)
  expect_equal(pw + pl + pt, 1, tolerance = 1e-6)
})

test_that("the PRO block requires landmark = horizon", {
  skip_on_cran()
  bad <- list(marker = "kccq", landmark = 2, margin = 5, direction = "higher.better", type = "continuous")
  expect_error(suppressWarnings(clinicalWinRatio(
    proFit(), arm = "arm", illness.time = "t_hosp", terminal.time = "t_term",
    terminal.status = "died", covariates = c("W", "W2"), horizon = 4, n.grid = 20,
    n.folds = 1, pro = bad)), regexp = "landmark = horizon")
})

test_that("death > PRO with no non-fatal tier is supported (Kev = 1)", {
  skip_on_cran()
  o <- suppressMessages(suppressWarnings(clinicalPSNB(
    proFit(seed = 2), arm = "arm", illness.time = character(0), terminal.time = "t_term",
    terminal.status = "died", covariates = c("W", "W2"), charter = c(0.5, 0.5),
    horizon = 4, n.grid = 20, n.folds = 1, pro = kccqSpec)))
  expect_equal(attr(o, "Tiers"), 2L)                         # death > KCCQ
  expect_true(is.finite(o[Estimand == "PSWR", `Pt Est`]))
})
