#!/usr/bin/env Rscript
# Cross-package RMST benchmark on a competing-risks RCT DGP with closed-form truth.
# Estimand = the MARGINAL (population-averaged) treatment difference; every method
# is used in its marginal form so the comparison is apples-to-apples:
#   - event-free RMST difference  (status = any first event)
#   - cause-1 years-lost difference (competing-risks specific)
# Methods: concrete getRMST / targetRMST; survRM2 (unadj + Tian-adjusted);
#          eventglm (pseudo-observation rmean); riskRegression::ate (if available).
# Because the trial is randomized, all marginal-form methods should be ~unbiased;
# the comparison then shows (a) calibration parity, (b) efficiency (CI width from
# adjustment), and (c) where the competing-risks / sparse-grid differences appear.
#   Rscript scripts/make-rmst-benchmark.R [B] [n]
suppressWarnings(suppressMessages({
  if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION"))
    devtools::load_all(".", quiet = TRUE) else library(concrete)
  library(data.table); library(survival)
}))
has <- function(p) requireNamespace(p, quietly = TRUE)
cli <- commandArgs(trailingOnly = TRUE)
B <- if (length(cli) >= 1) as.integer(cli[1]) else 60L
n <- if (length(cli) >= 2) as.integer(cli[2]) else 500L
tau <- 1500; tau_max <- 2500
grid <- seq(200, tau, by = 200)          # target grid for concrete

# --- competing-risks RCT DGP (two causes); treatment lowers cause 1 -----------
b10 <- 6e-4; b20 <- 4e-4; lc <- 3e-4
hz <- function(A, W1, W2) {
  l1 <- b10 * exp(-0.6 * A + 0.4 * W1); l2 <- b20 * exp(0.1 * A + 0.3 * W2)
  list(l1 = l1, l2 = l2, lt = l1 + l2)
}
truth <- function(Wn = 4e5) {
  W1 <- rnorm(Wn); W2 <- rnorm(Wn)
  efrmst <- function(a) { h <- hz(a, W1, W2); mean((1 - exp(-h$lt * tau)) / h$lt) }
  yll1   <- function(a) { h <- hz(a, W1, W2);
    mean((h$l1 / h$lt) * (tau - (1 - exp(-h$lt * tau)) / h$lt)) }
  list(efrmst = efrmst(1) - efrmst(0), yll1 = yll1(1) - yll1(0))
}
set.seed(20260607); TR <- truth()
cat(sprintf("TRUE marginal diffs: event-free RMST = %.2f | cause-1 YLL = %.2f\n",
            TR$efrmst, TR$yll1))

simData <- function(s) {
  set.seed(s)
  W1 <- rnorm(n); W2 <- rnorm(n); A <- rbinom(n, 1, 0.5); h <- hz(A, W1, W2)
  T1 <- rexp(n, h$l1); T2 <- rexp(n, h$l2); C <- rexp(n, lc)
  Tt <- pmin(T1, T2); J <- ifelse(T1 < T2, 1L, 2L)
  To <- pmin(Tt, C, tau_max); ev <- ifelse(To >= tau_max, 0L, ifelse(Tt <= C, J, 0L))
  data.table(id = seq_len(n), time = To, event = ev, arm = A, W1 = W1, W2 = W2,
             anyev = as.integer(ev %in% c(1, 2)))
}
trip <- function(est, lo, hi) c(est = est, lo = lo, hi = hi)
NA3 <- trip(NA, NA, NA)

# ---- method wrappers: each returns c(est, lo, hi) for the named estimand ------
pick <- function(out, estimand) {
  r <- if (estimand == "efrmst") out[grepl("RMST Diff", Estimand)]
       else out[grepl("LYL Diff", Estimand) & Event == 1]
  if (!nrow(r)) return(NA3)
  trip(r$`Pt Est`[1], r$`CI Low`[1], r$`CI Hi`[1])
}
# fit concrete ONCE per rep, return getRMST and targetRMST tables for both estimands
m_concrete_fit <- function(d) {
  a <- formatArguments(DataTable = d, EventTime = "time", EventType = "event", Treatment = "arm",
    ID = "id", Intervention = makeITT(), TargetTime = grid, TargetEvent = c(1, 2),
    CVArg = list(V = 5), UpdateMethod = "adaptive", EICStopRule = "absolute",
    MaxUpdateIter = 25, Verbose = FALSE)
  est <- suppressMessages(doConcrete(a))
  # targetRMST is excluded here (its iterative direct-targeting loop is ~50x slower
  # per fit); the getRMST-vs-targetRMST sparse-grid comparison lives in
  # make-rmst-comparison.R. This benchmark uses getRMST for the concrete entry.
  list(get = as.data.table(getRMST(est, Horizon = tau, Intervention = c(1, 2))))
}
m_survRM2 <- function(d, tau, adj = FALSE) {
  if (!has("survRM2")) return(NA3)
  r <- tryCatch(
    if (adj) survRM2::rmst2(d$time, d$anyev, d$arm, tau = tau, covariates = as.data.frame(d[, .(W1, W2)]))
    else     survRM2::rmst2(d$time, d$anyev, d$arm, tau = tau),
    error = function(e) NULL)
  if (is.null(r)) return(NA3)
  if (adj) { a <- r$RMST.difference.adjusted; trip(a["arm","coef"], a["arm","lower .95"], a["arm","upper .95"]) }
  else     { u <- r$unadjusted.result;        trip(u[1,1], u[1,2], u[1,3]) }
}
m_eventglm <- function(d, tau, estimand = c("efrmst","yll1")) {
  if (!has("eventglm")) return(NA3)
  estimand <- match.arg(estimand)
  g <- tryCatch({
    if (estimand == "efrmst") eventglm::rmeanglm(Surv(time, anyev) ~ arm, time = tau, data = d)
    else { d2 <- data.table::copy(d); d2[, eventf := factor(event)]   # competing risks needs factor status
           eventglm::rmeanglm(Surv(time, eventf) ~ arm, time = tau, cause = 1, data = d2) }
  }, error = function(e) NULL)
  if (is.null(g)) return(NA3)
  est <- unname(coef(g)["arm"])
  se <- tryCatch(sqrt(diag(vcov(g))["arm"]), error = function(e) NA_real_)  # sandwich SE
  z <- stats::qnorm(0.975)
  trip(est, est - z * se, est + z * se)
}
# riskRegression wrapper is filled in after inspecting its RMST API (see TODO).
m_riskreg <- function(d, tau, estimand) NA3

cover <- function(v, tr) is.finite(v["lo"]) & v["lo"] <= tr & tr <= v["hi"]
width <- function(v) unname(v["hi"] - v["lo"])

one <- function(s) {
  d <- simData(s)
  rows <- list()
  add <- function(method, estimand, tr, v) rows[[length(rows)+1]] <<- data.table(
    method = method, estimand = estimand, bias = unname(v["est"]) - tr,
    cover = cover(v, tr), width = width(v))
  tryCatch({
    cc <- m_concrete_fit(d)
    add("concrete",      "efrmst", TR$efrmst, pick(cc$get, "efrmst"))
    add("survRM2:unadj", "efrmst", TR$efrmst, m_survRM2(d, tau, adj = FALSE))
    add("survRM2:adj",   "efrmst", TR$efrmst, m_survRM2(d, tau, adj = TRUE))
    add("eventglm",      "efrmst", TR$efrmst, m_eventglm(d, tau, "efrmst"))
    add("concrete",      "yll1", TR$yll1, pick(cc$get, "yll1"))
    add("eventglm",      "yll1", TR$yll1, m_eventglm(d, tau, "yll1"))
    data.table::rbindlist(rows)
  }, error = function(e) NULL)
}
seeds <- 20260607L + seq_len(B)
envCores <- suppressWarnings(as.integer(Sys.getenv("VAL_CORES", unset = NA)))
uc <- if (!is.na(envCores) && envCores >= 1L) envCores else
  if (.Platform$OS.type == "unix") min(4L, max(1L, parallel::detectCores() - 1L)) else 1L
res <- data.table::rbindlist(Filter(Negate(is.null), parallel::mclapply(seeds, one, mc.cores = uc)))
summ <- res[, .(Reps = sum(!is.na(bias)), Bias = round(mean(bias, na.rm = TRUE), 2),
                Coverage = round(mean(cover, na.rm = TRUE), 3),
                CIwidth = round(mean(width, na.rm = TRUE), 1)),
            by = .(estimand, method)][order(estimand, method)]
print(summ)
data.table::fwrite(summ, "scripts/rmst-benchmark-summary.csv")
cat("\nNote: RCT DGP -> marginal-form methods should be ~unbiased; compare CIwidth",
    "(efficiency) and note survRM2 cannot do cause-1 YLL.\n")
