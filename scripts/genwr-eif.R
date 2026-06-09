## K-general analytic EIF assembly for the multistate win ratio. Sources the
## K-general engine (genwr-engine.R) and generalizes the tier loop / state-set
## logic from the validated K=3 assembly. armSetupFrom takes prebuilt (rmat, Ginv)
## so oracle and SL/cross-fit callers share it.
suppressWarnings(suppressMessages(library(data.table)))
source("scripts/genwr-engine.R")

## ---- per-subject at-risk Y_s and event N_{e,s} (n x M) from observed segments ----
buildYN <- function(D, grid) {
  n <- nrow(D); starts <- grid[-(M + 1L)]
  Y   <- setNames(lapply(ALIVE, function(s) matrix(0, n, M)), as.character(ALIVE))
  Nev <- setNames(lapply(structTrans, function(k) matrix(0, n, M)), structTrans)
  for (i in seq_len(n)) {
    obsEnd <- min(D$C[i], tau); et <- list()
    for (e in NF) { ti <- D[[paste0("t", e)]][i]
      if (is.finite(ti) && ti <= obsEnd) et[[length(et) + 1L]] <- c(ti, as.numeric(e)) }
    deathObs <- is.finite(D$tD[i]) && D$tD[i] <= obsEnd
    em <- if (length(et)) do.call(rbind, et) else matrix(numeric(0), 0, 2)
    if (nrow(em)) em <- em[order(em[, 1]), , drop = FALSE]
    s <- 0L; entry <- 0; segs <- list()
    if (nrow(em)) for (r in seq_len(nrow(em))) {
      te <- em[r, 1]; evc <- as.character(em[r, 2])
      segs[[length(segs) + 1L]] <- list(s = s, entry = entry, exit = te, ev = evc)
      s <- bitwOr(s, EVB[[evc]]); entry <- te
    }
    segs[[length(segs) + 1L]] <- if (deathObs) list(s = s, entry = entry, exit = D$tD[i], ev = "D")
                                 else           list(s = s, entry = entry, exit = obsEnd, ev = NA_character_)
    for (sg in segs) {
      jr <- which(starts >= sg$entry - 1e-12 & starts < sg$exit - 1e-12)
      if (length(jr)) Y[[as.character(sg$s)]][i, jr] <- 1
      if (!is.na(sg$ev)) { jx <- min(max(findInterval(sg$exit, grid), 1L), M)
        Nev[[paste0(sg$s, "|", sg$ev)]][i, jx] <- 1 }
    }
  }
  list(Y = Y, Nev = Nev)
}

## survival-to-tau adjoint over a state SET (terminal b; leaving set or death -> 0)
survKilled <- function(rmat, states, b = 1) {
  n <- ncol(rmat[[1]]); V <- setNames(lapply(states, function(s) matrix(0, M + 1L, n)), as.character(states))
  for (s in states) V[[as.character(s)]][M + 1L, ] <- b
  for (j in M:1) for (s in states) {
    Lam <- totRow(rmat, s, j); stay <- exp(-Lam); val <- stay * V[[as.character(s)]][j + 1L, ]
    for (tr in outTrans(s)) if (tr$to >= 0 && inSet(tr$to, states))
      val <- val + (getRow(rmat, s, tr$ev, j) / Lam) * (1 - stay) * V[[as.character(tr$to)]][j + 1L, ]
    V[[as.character(s)]][j, ] <- val
  }
  V
}

## pre adjoint U with self-transition reward w(j)*A_post(post(s))(j+1)
preAdjoint <- function(rmat, preStates, selfEv, Apost, w) {
  n <- ncol(rmat[[1]]); sb <- selfBit(selfEv)
  U <- setNames(lapply(preStates, function(s) matrix(0, M + 1L, n)), as.character(preStates))
  for (j in M:1) for (s in preStates) {
    Lam <- totRow(rmat, s, j); stay <- exp(-Lam); val <- stay * U[[as.character(s)]][j + 1L, ]
    for (tr in outTrans(s)) {
      mv <- (getRow(rmat, s, tr$ev, j) / Lam) * (1 - stay)
      if (tr$ev == selfEv) {
        Apj <- if (sb < 0) 1 else Apost[[as.character(bitwOr(s, sb))]][j + 1L, ]
        val <- val + mv * w[j] * Apj
      } else if (tr$to >= 0 && inSet(tr$to, preStates)) {
        val <- val + mv * U[[as.character(tr$to)]][j + 1L, ]
      }
    }
    U[[as.character(s)]][j, ] <- val
  }
  U
}

## generic per-subject survival-functional IF over a state set
survIFset <- function(rmat, YN, Ginv, states, V) {
  n <- ncol(rmat[[1]]); IF <- numeric(n)
  for (s in states) { Vs <- V[[as.character(s)]][-1L, , drop = FALSE]
    for (tr in outTrans(s)) {
      Vto <- if (tr$to >= 0 && inSet(tr$to, states)) V[[as.character(tr$to)]][-1L, , drop = FALSE] else 0
      clev <- t(Vto - Vs)
      dM <- YN$Nev[[paste0(s, "|", tr$ev)]] - YN$Y[[as.character(s)]] * t(rmat[[paste0(s, "|", tr$ev)]])
      IF <- IF + rowSums(clev * dM * Ginv)
    }
  }
  IF
}

## per-subject IF of Phi_w = int q^{(k)}(v) w(v) dv (weighted incidence)
phiIF <- function(rmat, YN, Ginv, preStates, postStates, selfEv, Apost, w) {
  n <- ncol(rmat[[1]]); sb <- selfBit(selfEv)
  U <- preAdjoint(rmat, preStates, selfEv, Apost, w); gi <- U[["0"]][1L, ]
  IF <- numeric(n)
  for (s in preStates) { Us <- U[[as.character(s)]][-1L, , drop = FALSE]
    for (tr in outTrans(s)) {
      if (tr$ev == selfEv) {
        Apj <- if (sb < 0) matrix(1, M, n) else Apost[[as.character(bitwOr(s, sb))]][-1L, , drop = FALSE]
        clev <- t(sweep(Apj, 1, w, "*") - Us)
      } else if (tr$to >= 0 && inSet(tr$to, preStates)) {
        clev <- t(U[[as.character(tr$to)]][-1L, , drop = FALSE] - Us)
      } else clev <- t(-Us)
      dM <- YN$Nev[[paste0(s, "|", tr$ev)]] - YN$Y[[as.character(s)]] * t(rmat[[paste0(s, "|", tr$ev)]])
      IF <- IF + rowSums(clev * dM * Ginv)
    }
  }
  if (length(postStates) && sb >= 0) {                       # post-part, weighted by self-firing time
    bracket <- survIFset(rmat, YN, Ginv, postStates, Apost)
    wfire <- numeric(n)
    for (s in preStates) wfire <- wfire + rowSums(YN$Nev[[paste0(s, "|", selfEv)]] * matrix(w, n, M, byrow = TRUE))
    IF <- IF + wfire * bracket
  }
  list(IF = (gi - mean(gi)) + IF, g = gi)
}

## ---- tier state sets, per-arm setup, tier IF blocks ----
tiedStates <- function(higher) ALIVE[vapply(ALIVE, function(s)
  all(vapply(higher, function(e) bitwAnd(s, EVB[[e]]) == 0L, logical(1))), logical(1))]
preStatesOf <- function(self, higher) { sb <- selfBit(self)
  tiedStates(higher)[vapply(tiedStates(higher), function(s) bitwAnd(s, sb) == 0L, logical(1))] }

armSetupFrom <- function(D, rmat, Ginv) {
  grid <- seq(0, tau, length.out = M + 1L)
  YN <- buildYN(D, grid); eng <- armTiers(D, rmat)
  list(D = D, rmat = rmat, YN = YN, Ginv = Ginv, n = nrow(D),
       q = lapply(eng, function(t) rowMeans(t$q)),
       m = lapply(eng, function(t) mean(t$m)))
}
armSetup <- function(D, censRate) {                          # oracle wrapper
  grid <- seq(0, tau, length.out = M + 1L)
  Ginv <- matrix(1 / exp(-censRate * grid[-1L]), nrow(D), M, byrow = TRUE)
  armSetupFrom(D, buildRateMats(D), Ginv)
}

DmIF <- function(arm, k) {
  ti <- TIERS[[k]]; if (ti$self == "D") return(numeric(arm$n))     # m = 1 for death tier
  states <- tiedStates(ti$higher); V <- survKilled(arm$rmat, states, 1); gi <- V[["0"]][1L, ]
  (gi - mean(gi)) + survIFset(arm$rmat, arm$YN, arm$Ginv, states, V)
}
PhiwIF <- function(arm, k, w) {
  ti <- TIERS[[k]]; pre <- preStatesOf(ti$self, ti$higher)
  post <- if (ti$self == "D") integer(0) else bitwOr(pre, selfBit(ti$self))
  Apost <- if (ti$self == "D") NULL else survKilled(arm$rmat, post, 1)
  phiIF(arm$rmat, arm$YN, arm$Ginv, pre, post, ti$self, Apost, w)$IF
}

## ---- cross-arm assembly over all K tiers ----
Qbar <- function(q) sum(q) - cumsum(q)
Qlag <- function(q) c(0, cumsum(q)[-length(q)])
assembleDP <- function(win, los) {
  DPt_win <- numeric(win$n); DPt_los <- numeric(los$n)
  for (k in seq_along(TIERS)) {
    DPt_win <- DPt_win + sum(los$q[[k]]) * DmIF(win, k) - PhiwIF(win, k, Qbar(los$q[[k]]))
    DPt_los <- DPt_los + PhiwIF(los, k, win$m[[k]] - Qlag(win$q[[k]]))
  }
  list(winnerIF = DPt_win, loserIF = DPt_los)
}

winRatioEIF <- function(trt, ctl, Signif = 0.05) {
  eng <- assembleWR(armTiers(trt$D, trt$rmat), armTiers(ctl$D, ctl$rmat))
  Pwin <- eng["Pwin"]; Ploss <- eng["Ploss"]
  win <- assembleDP(trt, ctl); los <- assembleDP(ctl, trt)
  DPwin_T <- win$winnerIF; DPwin_C <- win$loserIF; DPloss_C <- los$winnerIF; DPloss_T <- los$loserIF
  Ntot <- trt$n + ctl$n; piT <- trt$n / Ntot; piC <- ctl$n / Ntot; z <- qnorm(1 - Signif / 2)
  Pwin  <- Pwin  + mean(DPwin_T)  + mean(DPwin_C)            # one-step correction
  Ploss <- Ploss + mean(DPloss_T) + mean(DPloss_C)
  seGrad <- function(gw, gl) {
    Dt <- (1 / piT) * (gw * DPwin_T + gl * DPloss_T); Dc <- (1 / piC) * (gw * DPwin_C + gl * DPloss_C)
    sqrt((sum(Dt^2) + sum(Dc^2)) / Ntot^2)
  }
  WR <- Pwin / Ploss; seWR <- seGrad(1 / Ploss, -Pwin / Ploss^2); slWR <- seWR / WR
  c(WR = unname(WR), seWR = unname(seWR), lo = unname(WR * exp(-z * slWR)), hi = unname(WR * exp(z * slWR)))
}
