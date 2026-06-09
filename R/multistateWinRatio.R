#' Multistate engine + analytic EIF for the general hierarchical win ratio.
#'
#' Internal machinery for the death-priority (clinical) win ratio over an ordered
#' hierarchy of fatal + non-fatal time-to-event tiers, where a higher-priority
#' event may occur after a lower-priority one. States are subsets of the
#' \eqn{K-1} non-fatal events (bitmask); transitions are auto-generated. The
#' occupancy engine and the adjoint-value influence functions are validated
#' against a brute-force pairwise win ratio (see `scripts/genwr-*.R`).
#'
#' `.msEngine(K, grid)` returns an environment holding the configuration and all
#' engine/EIF functions (defined in one frame so they close over the config -- no
#' package globals).
#' @keywords internal
#' @noRd
.msEngine <- function(K, grid) {
  tau <- max(grid); M <- length(grid) - 1L
  NF <- as.character(seq_len(K - 1L))                       # non-fatal labels, priority order
  EVB <- stats::setNames(2L^(seq_len(K - 1L) - 1L), NF)     # label -> bit
  ALIVE <- 0:(2L^(K - 1L) - 1L)
  structTrans <- character(0)
  for (s in ALIVE) { structTrans <- c(structTrans, paste0(s, "|D"))
    for (e in NF) if (bitwAnd(s, EVB[[e]]) == 0L) structTrans <- c(structTrans, paste0(s, "|", e)) }
  TIERS <- c(list(list(self = "D", higher = character(0), lower = NF)),
             lapply(2:K, function(k) list(self = NF[k - 1L],
               higher = NF[seq_len(k - 2L)], lower = if (k < K) NF[k:(K - 1L)] else character(0))))
  selfBit <- function(ev) if (ev == "D") -1L else EVB[[ev]]
  inSet   <- function(s, set) s %in% set
  hasEvB  <- function(s, ev) bitwAnd(s, EVB[[ev]]) > 0L
  outTrans <- function(s) {
    tr <- list(list(ev = "D", to = -1L))
    for (e in NF) if (!hasEvB(s, e)) tr <- c(tr, list(list(ev = e, to = bitwOr(s, EVB[[e]]))))
    tr
  }
  getRow <- function(rmat, s, ev, j) { m <- rmat[[paste0(s, "|", ev)]]; if (is.null(m)) NULL else m[j, ] }
  totRow <- function(rmat, s, j) { tot <- 0; for (tr in outTrans(s)) tot <- tot + getRow(rmat, s, tr$ev, j); tot }

  ## ---- point-estimate occupancy curves ----
  tierNonFatal <- function(rmat, selfEv, higher, lower) {
    n <- ncol(rmat[[1]])
    isPre <- function(s) (!hasEvB(s, selfEv)) && all(vapply(higher, function(e) !hasEvB(s, e), logical(1)))
    pre  <- ALIVE[vapply(ALIVE, isPre, logical(1))]; post <- bitwOr(pre, EVB[[selfEv]])
    A <- stats::setNames(lapply(post, function(s) matrix(0, M + 1L, n)), as.character(post))
    for (s in post) A[[as.character(s)]][M + 1L, ] <- 1
    for (j in M:1) for (s in post) {
      Lam <- totRow(rmat, s, j); stay <- exp(-Lam); val <- stay * A[[as.character(s)]][j + 1L, ]
      for (e in lower) if (!hasEvB(s, e)) val <- val +
        (getRow(rmat, s, e, j) / Lam) * (1 - stay) * A[[as.character(bitwOr(s, EVB[[e]]))]][j + 1L, ]
      A[[as.character(s)]][j, ] <- val
    }
    p <- stats::setNames(lapply(pre, function(s) matrix(0, M + 1L, n)), as.character(pre))
    p[["0"]][1L, ] <- 1; q <- matrix(0, M, n)
    for (j in 1:M) {
      pnext <- stats::setNames(lapply(pre, function(s) numeric(n)), as.character(pre))
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
    p <- stats::setNames(lapply(ALIVE, function(s) matrix(0, M + 1L, n)), as.character(ALIVE))
    p[["0"]][1L, ] <- 1; q <- matrix(0, M, n)
    for (j in 1:M) {
      pnext <- stats::setNames(lapply(ALIVE, function(s) numeric(n)), as.character(ALIVE))
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
  armTiers <- function(rmat) lapply(TIERS, function(ti)
    if (ti$self == "D") tierDeath(rmat) else tierNonFatal(rmat, ti$self, ti$higher, ti$lower))
  Wtier <- function(winner, loser, k) {
    qW <- rowMeans(winner[[k]]$q); qL <- rowMeans(loser[[k]]$q); mW <- mean(winner[[k]]$m)
    sum((mW - c(0, cumsum(qW)[-length(qW)])) * qL)
  }
  assembleWR <- function(TT, CC) {
    Pwin  <- sum(vapply(seq_along(TIERS), function(k) Wtier(TT, CC, k), numeric(1)))
    Ploss <- sum(vapply(seq_along(TIERS), function(k) Wtier(CC, TT, k), numeric(1)))
    c(Pwin = Pwin, Ploss = Ploss)
  }

  ## ---- per-subject at-risk / event arrays from observed multistate segments ----
  ## Grid-indexed: loop over the M intervals with vector ops over subjects (O(M)
  ## not O(n)). State at grid node j = the non-fatal events observed by then; a
  ## subject is at risk in interval j while alive and uncensored at its start.
  buildYN <- function(D) {
    n <- nrow(D); starts <- grid[-(M + 1L)]; tcols <- paste0("t", NF); nbit <- as.integer(EVB[NF])
    Y   <- stats::setNames(lapply(ALIVE, function(s) matrix(0, n, M)), as.character(ALIVE))
    Nev <- stats::setNames(lapply(structTrans, function(k) matrix(0, n, M)), structTrans)
    obsEnd <- pmin(D$C, tau)
    TT <- matrix(Inf, n, length(NF))                          # observed non-fatal event times
    for (ei in seq_along(NF)) { ti <- D[[tcols[ei]]]; TT[, ei] <- ifelse(ti <= obsEnd, ti, Inf) }
    tdeath  <- ifelse(D$tD <= obsEnd, D$tD, Inf)
    exitAll <- pmin(tdeath, obsEnd)                           # end of the alive period
    for (j in 1:M) {
      gj <- starts[j]; state <- integer(n)
      for (ei in seq_along(NF)) state <- bitwOr(state, ifelse(TT[, ei] <= gj, nbit[ei], 0L))
      atrisk <- gj < exitAll
      if (any(atrisk)) for (s in unique(state[atrisk])) {
        idx <- which(atrisk & state == s); Y[[as.character(s)]][cbind(idx, j)] <- 1
      }
    }
    evspec <- c(lapply(seq_along(NF), function(ei) list(te = TT[, ei], lab = NF[ei], excl = ei)),
                list(list(te = tdeath, lab = "D", excl = 0L)))
    for (sp in evspec) {                                      # event fires from state = events strictly earlier
      obs <- is.finite(sp$te); if (!any(obs)) next
      jx <- pmin(pmax(findInterval(sp$te, grid), 1L), M); fromS <- integer(n)
      for (ej in seq_along(NF)) if (ej != sp$excl) fromS <- bitwOr(fromS, ifelse(TT[, ej] < sp$te, nbit[ej], 0L))
      for (f in unique(fromS[obs])) { ii <- which(obs & fromS == f)
        Nev[[paste0(f, "|", sp$lab)]][cbind(ii, jx[ii])] <- 1 }
    }
    list(Y = Y, Nev = Nev)
  }

  ## ---- adjoint-value influence functions ----
  survKilled <- function(rmat, states, b = 1) {
    n <- ncol(rmat[[1]]); V <- stats::setNames(lapply(states, function(s) matrix(0, M + 1L, n)), as.character(states))
    for (s in states) V[[as.character(s)]][M + 1L, ] <- b
    for (j in M:1) for (s in states) {
      Lam <- totRow(rmat, s, j); stay <- exp(-Lam); val <- stay * V[[as.character(s)]][j + 1L, ]
      for (tr in outTrans(s)) if (tr$to >= 0 && inSet(tr$to, states))
        val <- val + (getRow(rmat, s, tr$ev, j) / Lam) * (1 - stay) * V[[as.character(tr$to)]][j + 1L, ]
      V[[as.character(s)]][j, ] <- val
    }
    V
  }
  preAdjoint <- function(rmat, preStates, selfEv, Apost, w) {
    n <- ncol(rmat[[1]]); sb <- selfBit(selfEv)
    U <- stats::setNames(lapply(preStates, function(s) matrix(0, M + 1L, n)), as.character(preStates))
    for (j in M:1) for (s in preStates) {
      Lam <- totRow(rmat, s, j); stay <- exp(-Lam); val <- stay * U[[as.character(s)]][j + 1L, ]
      for (tr in outTrans(s)) {
        mv <- (getRow(rmat, s, tr$ev, j) / Lam) * (1 - stay)
        if (tr$ev == selfEv) {
          Apj <- if (sb < 0) 1 else Apost[[as.character(bitwOr(s, sb))]][j + 1L, ]
          val <- val + mv * w[j] * Apj
        } else if (tr$to >= 0 && inSet(tr$to, preStates)) val <- val + mv * U[[as.character(tr$to)]][j + 1L, ]
      }
      U[[as.character(s)]][j, ] <- val
    }
    U
  }
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
  phiIF <- function(rmat, YN, Ginv, preStates, postStates, selfEv, Apost, w) {
    n <- ncol(rmat[[1]]); sb <- selfBit(selfEv)
    U <- preAdjoint(rmat, preStates, selfEv, Apost, w); gi <- U[["0"]][1L, ]; IF <- numeric(n)
    for (s in preStates) { Us <- U[[as.character(s)]][-1L, , drop = FALSE]
      for (tr in outTrans(s)) {
        if (tr$ev == selfEv) {
          Apj <- if (sb < 0) matrix(1, M, n) else Apost[[as.character(bitwOr(s, sb))]][-1L, , drop = FALSE]
          clev <- t(sweep(Apj, 1, w, "*") - Us)
        } else if (tr$to >= 0 && inSet(tr$to, preStates)) clev <- t(U[[as.character(tr$to)]][-1L, , drop = FALSE] - Us)
        else clev <- t(-Us)
        dM <- YN$Nev[[paste0(s, "|", tr$ev)]] - YN$Y[[as.character(s)]] * t(rmat[[paste0(s, "|", tr$ev)]])
        IF <- IF + rowSums(clev * dM * Ginv)
      }
    }
    if (length(postStates) && sb >= 0) {
      bracket <- survIFset(rmat, YN, Ginv, postStates, Apost); wfire <- numeric(n)
      for (s in preStates) wfire <- wfire + rowSums(YN$Nev[[paste0(s, "|", selfEv)]] * matrix(w, n, M, byrow = TRUE))
      IF <- IF + wfire * bracket
    }
    (gi - mean(gi)) + IF
  }
  tiedStates <- function(higher) ALIVE[vapply(ALIVE, function(s)
    all(vapply(higher, function(e) bitwAnd(s, EVB[[e]]) == 0L, logical(1))), logical(1))]
  preStatesOf <- function(self, higher) { sb <- selfBit(self)
    tiedStates(higher)[vapply(tiedStates(higher), function(s) bitwAnd(s, sb) == 0L, logical(1))] }

  armSetup <- function(D, rmat, Ginv) {
    YN <- buildYN(D); eng <- armTiers(rmat)
    list(rmat = rmat, YN = YN, Ginv = Ginv, n = nrow(D),
         q = lapply(eng, function(t) rowMeans(t$q)), m = lapply(eng, function(t) mean(t$m)))
  }
  DmIF <- function(arm, k) {
    ti <- TIERS[[k]]; if (ti$self == "D") return(numeric(arm$n))
    states <- tiedStates(ti$higher); V <- survKilled(arm$rmat, states, 1); gi <- V[["0"]][1L, ]
    (gi - mean(gi)) + survIFset(arm$rmat, arm$YN, arm$Ginv, states, V)
  }
  PhiwIF <- function(arm, k, w) {
    ti <- TIERS[[k]]; pre <- preStatesOf(ti$self, ti$higher)
    post <- if (ti$self == "D") integer(0) else bitwOr(pre, selfBit(ti$self))
    Apost <- if (ti$self == "D") NULL else survKilled(arm$rmat, post, 1)
    phiIF(arm$rmat, arm$YN, arm$Ginv, pre, post, ti$self, Apost, w)
  }
  Qbar <- function(q) sum(q) - cumsum(q); Qlag <- function(q) c(0, cumsum(q)[-length(q)])
  assembleDP <- function(win, los) {
    DPt_win <- numeric(win$n); DPt_los <- numeric(los$n)
    for (k in seq_along(TIERS)) {
      DPt_win <- DPt_win + sum(los$q[[k]]) * DmIF(win, k) - PhiwIF(win, k, Qbar(los$q[[k]]))
      DPt_los <- DPt_los + PhiwIF(los, k, win$m[[k]] - Qlag(win$q[[k]]))
    }
    list(winnerIF = DPt_win, loserIF = DPt_los)
  }
  environment()
}

#' Observed multistate segments (long) for transition-hazard fitting.
#' @keywords internal
#' @noRd
.msSegments <- function(eng, D, covariates) {
  NF <- eng$NF; EVB <- eng$EVB; tau <- eng$tau; tcols <- paste0("t", NF)
  rows <- vector("list", nrow(D))
  for (i in seq_len(nrow(D))) {
    obsEnd <- min(D$C[i], tau); et <- list()
    for (ei in seq_along(NF)) { ti <- D[[tcols[ei]]][i]
      if (is.finite(ti) && ti <= obsEnd) et[[length(et) + 1L]] <- c(ti, as.numeric(NF[ei])) }
    deathObs <- is.finite(D$tD[i]) && D$tD[i] <= obsEnd
    em <- if (length(et)) do.call(rbind, et) else matrix(numeric(0), 0, 2)
    if (nrow(em)) em <- em[order(em[, 1]), , drop = FALSE]
    s <- 0L; entry <- 0; segs <- list()
    if (nrow(em)) for (r in seq_len(nrow(em))) { te <- em[r, 1]; evc <- as.character(em[r, 2])
      segs[[length(segs) + 1L]] <- data.table::data.table(s = s, entry = entry, exit = te, ev = evc)
      s <- bitwOr(s, EVB[[evc]]); entry <- te }
    segs[[length(segs) + 1L]] <- if (deathObs) data.table::data.table(s = s, entry = entry, exit = D$tD[i], ev = "D")
                                 else           data.table::data.table(s = s, entry = entry, exit = obsEnd, ev = NA_character_)
    seg <- as.data.frame(data.table::rbindlist(segs))
    rows[[i]] <- cbind(seg, D[rep(i, nrow(seg)), covariates, drop = FALSE])
  }
  as.data.frame(data.table::rbindlist(rows))   # data.frame so [, cov, drop=FALSE] selects columns
}

#' Cross-fitted Super Learner transition hazards + censoring G for one arm.
#' Returns `rmat` (list of M x n cumulative-hazard increments) and `Ginv` (n x M,
#' inverse lagged censoring survival). `n.folds <= 1` gives in-sample fits.
#' @keywords internal
#' @noRd
.msNuisances <- function(eng, D, covariates, SL.library, n.folds) {
  M <- eng$M; grid <- eng$grid; structTrans <- eng$structTrans; tau <- eng$tau; n <- nrow(D)
  V <- max(1L, min(as.integer(n.folds), floor(n / 30)))
  rmat <- stats::setNames(lapply(structTrans, function(k) matrix(1e-10, M, n)), structTrans)
  incC <- matrix(0, M, n)
  fold <- if (V <= 1L) rep(1L, n) else sample(rep(seq_len(V), length.out = n))
  obsT <- pmin(D$C, D$tD, tau); censInd <- as.integer(D$C < pmin(D$tD, tau))
  for (v in seq_len(max(fold))) {
    tr <- if (V <= 1L) seq_len(n) else which(fold != v)
    te <- if (V <= 1L) seq_len(n) else which(fold == v)
    seg <- .msSegments(eng, D[tr, , drop = FALSE], covariates); CovTe <- D[te, covariates, drop = FALSE]
    for (key in structTrans) {
      sp <- strsplit(key, "\\|")[[1]]; s <- as.integer(sp[1]); ev <- sp[2]
      sub <- seg[seg$s == s, ]
      rmat[[key]][, te] <- tryCatch(suppressWarnings({   # SL CV-weights any non-converged glm; preds bounded
        fit <- fitTransitionSL(sub$entry, sub$exit, as.integer(sub$ev == ev & !is.na(sub$ev)),
                               sub[, covariates, drop = FALSE], grid, SL.library = SL.library)
        predictTransitionSL(fit, CovTe)
      }), error = function(e) matrix(1e-10, M, length(te)))
    }
    cFit <- tryCatch(suppressWarnings(fitTransitionSL(rep(0, length(tr)), obsT[tr], censInd[tr],
              D[tr, covariates, drop = FALSE], grid, SL.library = SL.library)), error = function(e) NULL)
    incC[, te] <- if (is.null(cFit)) matrix(1e-8, M, length(te)) else predictTransitionSL(cFit, CovTe)
  }
  Glag <- pmax(rbind(1, exp(-apply(incC, 2, cumsum)))[1:M, , drop = FALSE], 0.05)
  list(rmat = rmat, Ginv = 1 / t(Glag))
}

#' Assemble the win ratio / win odds / net benefit table with IF inference.
#' @keywords internal
#' @noRd
.msWinRatioOut <- function(eng, trt, ctl, Signif) {
  base <- eng$assembleWR(eng$armTiers(trt$rmat), eng$armTiers(ctl$rmat))
  win <- eng$assembleDP(trt, ctl); los <- eng$assembleDP(ctl, trt)
  DPwin_T <- win$winnerIF; DPwin_C <- win$loserIF; DPloss_C <- los$winnerIF; DPloss_T <- los$loserIF
  Ntot <- trt$n + ctl$n; piT <- trt$n / Ntot; piC <- ctl$n / Ntot; z <- stats::qnorm(1 - Signif / 2)
  Pwin  <- unname(base["Pwin"]  + mean(DPwin_T)  + mean(DPwin_C))     # one-step
  Ploss <- unname(base["Ploss"] + mean(DPloss_T) + mean(DPloss_C))
  Ptie  <- max(0, 1 - Pwin - Ploss)
  seGrad <- function(gw, gl) {
    Dt <- (1 / piT) * (gw * DPwin_T + gl * DPloss_T); Dc <- (1 / piC) * (gw * DPwin_C + gl * DPloss_C)
    sqrt((sum(Dt^2) + sum(Dc^2)) / Ntot^2)
  }
  ratioRow <- function(label, val, gw, gl) { se <- seGrad(gw, gl); sl <- se / val
    data.table::data.table(Estimand = label, `Pt Est` = val, se = se,
      `CI Low` = val * exp(-z * sl), `CI Hi` = val * exp(z * sl),
      pValue = 2 * stats::pnorm(-abs(log(val) / sl))) }
  diffRow <- function(label, val, gw, gl) { se <- seGrad(gw, gl)
    data.table::data.table(Estimand = label, `Pt Est` = val, se = se,
      `CI Low` = val - z * se, `CI Hi` = val + z * se, pValue = 2 * stats::pnorm(-abs(val / se))) }
  probRow <- function(label, val, gw, gl) { se <- seGrad(gw, gl)
    data.table::data.table(Estimand = label, `Pt Est` = val, se = se,
      `CI Low` = val - z * se, `CI Hi` = val + z * se, pValue = NA_real_) }
  NB <- Pwin - Ploss; dWO <- 2 / (1 - NB)^2
  data.table::rbindlist(list(
    ratioRow("Win Ratio", Pwin / Ploss, 1 / Ploss, -Pwin / Ploss^2),
    ratioRow("Win Odds", (1 + NB) / (1 - NB), dWO, -dWO),
    diffRow("Net Benefit", NB, 1, -1),
    probRow("P(win)", Pwin, 1, 0), probRow("P(loss)", Ploss, 0, 1), probRow("P(tie)", Ptie, -1, -1)))
}
