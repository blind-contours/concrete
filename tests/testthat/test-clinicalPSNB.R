psnbFit <- function(seed = 7, n = 700) {
  set.seed(seed); A <- rep(0:1, each = n); W <- rnorm(2*n)
  tD <- rexp(2*n, 0.08*exp(0.3*W - 0.6*A)); tH <- rexp(2*n, 0.25*exp(0.2*W - 0.2*A)); C <- rexp(2*n, 0.05)
  obs <- pmin(tD, C, 5)
  data.frame(arm = A, t_hosp = ifelse(tH < obs, tH, NA), t_term = obs,
             died = as.integer(tD <= pmin(C, 5)), W = W, W2 = rnorm(2*n))
}
psnbArgs <- function(dat) list(data = dat, arm = "arm", illness.time = "t_hosp",
  terminal.time = "t_term", terminal.status = "died", covariates = c("W", "W2"),
  horizon = 4, n.grid = 30, n.folds = 1)

test_that("PSNB = sum of charter-weighted stage-conditional net benefits (exact, within call)", {
  skip_on_cran()
  ar <- psnbArgs(psnbFit())
  for (ch in list(c(0.5, 0.5), c(0.8, 0.2), "reach")) {
    o <- suppressMessages(suppressWarnings(do.call(clinicalPSNB, c(ar, list(charter = ch)))))
    rk <- o[grepl("^Reach", Estimand), `Pt Est`]
    dk <- o[grepl("^NetBenefit", Estimand), `Pt Est`]
    a  <- if (identical(ch, "reach")) rk else ch / sum(ch)
    expect_equal(o[Estimand == "PSNB", `Pt Est`], sum(a * dk), tolerance = 1e-8)
    ## PSWR > 1 iff PSNB > 0 (directional concordance)
    expect_equal(o[Estimand == "PSWR", `Pt Est`] > 1, o[Estimand == "PSNB", `Pt Est`] > 0)
    expect_true(all(o$se[!is.na(o$se)] > 0))
  }
})

test_that("charter='reach' reproduces the standard net benefit / win ratio", {
  skip_on_cran()
  dat <- psnbFit(seed = 11); ar <- psnbArgs(dat)
  wr  <- suppressMessages(suppressWarnings(do.call(clinicalWinRatio, ar)))
  ps  <- suppressMessages(suppressWarnings(do.call(clinicalPSNB, c(ar, list(charter = "reach")))))
  ## within one PSNB call, reach-charter PSNB = sum reach_k * NB_k = standard NB exactly
  rk <- ps[grepl("^Reach", Estimand), `Pt Est`]; dk <- ps[grepl("^NetBenefit", Estimand), `Pt Est`]
  expect_equal(ps[Estimand == "PSNB", `Pt Est`], sum(rk * dk), tolerance = 1e-8)
  ## across independent fits, matches clinicalWinRatio up to nuisance-fit noise
  expect_equal(ps[Estimand == "PSWR", `Pt Est`],
               as.data.frame(wr)[wr$Estimand == "Win Ratio", "Pt Est"], tolerance = 0.05)
  expect_equal(ps[Estimand == "PSNB", `Pt Est`],
               as.data.frame(wr)[wr$Estimand == "Net Benefit", "Pt Est"], tolerance = 0.02)
  expect_identical(ps[Estimand == "Reach[D]", `Pt Est`], 1)   # death tier always reached
})

test_that("an extreme charter isolates a single layer's stage-conditional effect", {
  skip_on_cran()
  ar <- psnbArgs(psnbFit(seed = 3))
  o <- suppressMessages(suppressWarnings(do.call(clinicalPSNB, c(ar, list(charter = c(1, 0))))))
  expect_equal(o[Estimand == "PSNB", `Pt Est`], o[Estimand == "NetBenefit[D]", `Pt Est`], tolerance = 1e-8)
})

test_that("clinicalPSNB validates the charter", {
  ar <- psnbArgs(psnbFit(n = 100))
  expect_error(suppressWarnings(do.call(clinicalPSNB, c(ar, list(charter = c(1, 1, 1))))),
               regexp = "length K")
  expect_error(suppressWarnings(do.call(clinicalPSNB, c(ar, list(charter = c(-1, 2))))),
               regexp = "non-negative")
})
