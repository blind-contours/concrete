## K-general multistate clinical win-ratio engine + configurable DGP + brute force
## + oracle point estimate. Hierarchy death E1 > E2 > ... > EK (K-1 non-fatal,
## non-recurrent, Markov). States = bitmask over (K-1) non-fatal events.
## configK(K) sets the state space / transitions / tiers; configDGP(K) the rates.
suppressWarnings(suppressMessages(library(data.table)))
tau <- 3.0
popcount <- function(s) sum(bitwAnd(bitwShiftR(as.integer(s), 0:30), 1L))

## ---- configuration ----
configK <- function(K) {
  KK <<- K
  NF <<- as.character(seq_len(K - 1L))                      # non-fatal labels, priority order
  EVB <<- setNames(2L^(seq_len(K - 1L) - 1L), NF)           # label -> bit
  ALIVE <<- 0:(2L^(K - 1L) - 1L)
  st <- character(0)
  for (s in ALIVE) { st <- c(st, paste0(s, "|D"))
    for (e in NF) if (bitwAnd(s, EVB[[e]]) == 0L) st <- c(st, paste0(s, "|", e)) }
  structTrans <<- st
  TIERS <<- c(list(list(self = "D", higher = character(0), lower = NF)),
              lapply(2:K, function(k) list(self = NF[k - 1L],
                higher = NF[seq_len(k - 2L)],
                lower  = if (k < K) NF[k:(K - 1L)] else character(0))))
}
selfBit <- function(ev) if (ev == "D") -1L else EVB[[ev]]
inSet   <- function(s, set) s %in% set
hasEvB  <- function(s, ev) bitwAnd(s, EVB[[ev]]) > 0L
outTrans <- function(s) {
  tr <- list(list(ev = "D", to = -1L))
  for (e in NF) if (!hasEvB(s, e)) tr <- c(tr, list(list(ev = e, to = bitwOr(s, EVB[[e]]))))
  tr
}

## ---- configurable DGP: death rises with accumulated morbidity; lower-priority
## events more frequent so they tend to precede higher-priority ones ----
DGP <- new.env()
configDGP <- function(K) {
  DGP$d0 <- 0.04; DGP$morb <- 1.5; DGP$bWd <- 0.30; DGP$bAd <- -0.40
  DGP$bWe <- 0.30; DGP$bAe <- -0.30
  DGP$nfbase <- setNames(seq(0.10, 0.22, length.out = K - 1L), NF)  # NF[1]=highest prio=rarest
}
rateState <- function(W, a, s, ev) {
  if (ev == "D") DGP$d0 * DGP$morb^popcount(s) * exp(DGP$bWd * W + DGP$bAd * a)
  else DGP$nfbase[[ev]] * exp(DGP$bWe * W + DGP$bAe * a)
}

## ---- Gillespie full-history sim; records first time of each event ----
simOne <- function(W, a, cens = 0) {
  s <- 0L; t <- 0; tev <- setNames(rep(Inf, length(NF)), NF); tD <- Inf
  C <- if (cens > 0) rexp(1, cens) else Inf
  repeat {
    trs <- outTrans(s); evs <- vapply(trs, function(x) x$ev, character(1))
    rs <- vapply(evs, function(e) rateState(W, a, s, e), numeric(1))
    tot <- sum(rs); t <- t + rexp(1, tot); if (t > tau) break
    pick <- evs[which(runif(1) * tot <= cumsum(rs))[1]]
    if (pick == "D") { tD <- t; break }
    tev[[pick]] <- min(tev[[pick]], t); s <- bitwOr(s, EVB[[pick]])
  }
  c(tD = tD, setNames(tev, paste0("t", NF)), C = C)
}
simArm <- function(n, a, cens = 0) {
  W <- rnorm(n)
  h <- t(vapply(seq_len(n), function(i) simOne(W[i], a, cens), numeric(length(NF) + 2L)))
  data.table(W = W, a = a, h)
}

## ---- brute-force pairwise K-tier hierarchical win ratio (uncensored truth) ----
cap <- function(x) ifelse(x <= tau, x, Inf)
bruteWR <- function(trt, ctl, npairs = 3e6) {
  it <- sample(nrow(trt), npairs, TRUE); ic <- sample(nrow(ctl), npairs, TRUE)
  cols <- c("tD", paste0("t", NF))                          # priority order: death, then NF
  xt <- lapply(cols, function(cc) cap(trt[[cc]][it])); yt <- lapply(cols, function(cc) cap(ctl[[cc]][ic]))
  win <- loss <- dec <- logical(npairs)
  for (tier in seq_along(cols)) { u <- !dec; x <- xt[[tier]]; y <- yt[[tier]]
    win <- win | (u & x > y); loss <- loss | (u & x < y); dec <- dec | (u & x != y) }
  c(Pwin = mean(win), Ploss = mean(loss), Ptie = mean(!dec), WR = mean(win) / mean(loss))
}

## ---- engine: occupancy from per-subject Lambda-increment matrices ----
M <- 60L; dt <- tau / M
setGrid <- function(MM) { M <<- as.integer(MM); dt <<- tau / M }
buildRateMats <- function(D) {                              # oracle: true rates -> Lambda increments
  n <- nrow(D); rmat <- list()
  for (key in structTrans) { sp <- strsplit(key, "\\|")[[1]]; s <- as.integer(sp[1]); ev <- sp[2]
    rmat[[key]] <- matrix(vapply(seq_len(n), function(i) rateState(D$W[i], D$a[i], s, ev), numeric(1)) * dt,
                          M, n, byrow = TRUE) }
  rmat
}
getR   <- function(rmat, s, ev) { k <- paste0(s, "|", ev); if (is.null(rmat[[k]])) NULL else rmat[[k]] }
getRow <- function(rmat, s, ev, j) { m <- getR(rmat, s, ev); if (is.null(m)) NULL else m[j, ] }
totRow <- function(rmat, s, j) { tot <- 0; for (tr in outTrans(s)) tot <- tot + getRow(rmat, s, tr$ev, j); tot }

tierNonFatal <- function(rmat, selfEv, higher, lower) {
  n <- ncol(rmat[[1]])
  isPre <- function(s) (!hasEvB(s, selfEv)) && all(vapply(higher, function(e) !hasEvB(s, e), logical(1)))
  pre  <- ALIVE[vapply(ALIVE, isPre, logical(1))]
  post <- bitwOr(pre, EVB[[selfEv]])
  A <- setNames(lapply(post, function(s) matrix(0, M + 1L, n)), as.character(post))
  for (s in post) A[[as.character(s)]][M + 1L, ] <- 1
  for (j in M:1) for (s in post) {
    Lam <- totRow(rmat, s, j); stay <- exp(-Lam); val <- stay * A[[as.character(s)]][j + 1L, ]
    for (e in lower) if (!hasEvB(s, e)) val <- val +
      (getRow(rmat, s, e, j) / Lam) * (1 - stay) * A[[as.character(bitwOr(s, EVB[[e]]))]][j + 1L, ]
    A[[as.character(s)]][j, ] <- val
  }
  p <- setNames(lapply(pre, function(s) matrix(0, M + 1L, n)), as.character(pre))
  p[["0"]][1L, ] <- 1; q <- matrix(0, M, n)
  for (j in 1:M) {
    pnext <- setNames(lapply(pre, function(s) numeric(n)), as.character(pre))
    for (s in pre) {
      ps <- p[[as.character(s)]][j, ]; if (all(ps == 0)) next
      Lam <- totRow(rmat, s, j); stay <- exp(-Lam)
      pnext[[as.character(s)]] <- pnext[[as.character(s)]] + ps * stay
      for (e in lower) if (!hasEvB(s, e)) { to <- bitwOr(s, EVB[[e]])
        pnext[[as.character(to)]] <- pnext[[as.character(to)]] + ps * (getRow(rmat, s, e, j) / Lam) * (1 - stay) }
      rs <- getRow(rmat, s, selfEv, j)
      if (!is.null(rs)) q[j, ] <- q[j, ] + ps * (rs / Lam) * (1 - stay) *
          A[[as.character(bitwOr(s, EVB[[selfEv]]))]][j + 1L, ]
    }
    for (s in pre) p[[as.character(s)]][j + 1L, ] <- pnext[[as.character(s)]]
  }
  list(q = q, m = Reduce(`+`, lapply(pre, function(s) p[[as.character(s)]][M + 1L, ])) + colSums(q))
}

tierDeath <- function(rmat) {
  n <- ncol(rmat[[1]])
  p <- setNames(lapply(ALIVE, function(s) matrix(0, M + 1L, n)), as.character(ALIVE))
  p[["0"]][1L, ] <- 1; q <- matrix(0, M, n)
  for (j in 1:M) {
    pnext <- setNames(lapply(ALIVE, function(s) numeric(n)), as.character(ALIVE))
    for (s in ALIVE) {
      ps <- p[[as.character(s)]][j, ]; if (all(ps == 0)) next
      Lam <- totRow(rmat, s, j); stay <- exp(-Lam)
      pnext[[as.character(s)]] <- pnext[[as.character(s)]] + ps * stay
      for (e in NF) if (!hasEvB(s, e)) { to <- bitwOr(s, EVB[[e]])
        pnext[[as.character(to)]] <- pnext[[as.character(to)]] + ps * (getRow(rmat, s, e, j) / Lam) * (1 - stay) }
      q[j, ] <- q[j, ] + ps * (getRow(rmat, s, "D", j) / Lam) * (1 - stay)
    }
    for (s in ALIVE) p[[as.character(s)]][j + 1L, ] <- pnext[[as.character(s)]]
  }
  list(q = q, m = rep(1, n))
}

armTiers <- function(D, rmat = buildRateMats(D)) lapply(TIERS, function(ti)
  if (ti$self == "D") tierDeath(rmat) else tierNonFatal(rmat, ti$self, ti$higher, ti$lower))
Wtier <- function(winner, loser, k) {
  qW <- rowMeans(winner[[k]]$q); qL <- rowMeans(loser[[k]]$q); mW <- mean(winner[[k]]$m)
  sum((mW - c(0, cumsum(qW)[-length(qW)])) * qL)
}
assembleWR <- function(TT, CC) {
  Pwin  <- sum(vapply(seq_along(TIERS), function(k) Wtier(TT, CC, k), numeric(1)))
  Ploss <- sum(vapply(seq_along(TIERS), function(k) Wtier(CC, TT, k), numeric(1)))
  c(Pwin = Pwin, Ploss = Ploss, Ptie = 1 - Pwin - Ploss, WR = Pwin / Ploss)
}
engineWR <- function(trtD, ctlD) assembleWR(armTiers(trtD), armTiers(ctlD))
