#' submodelAndEIC
#'
#' Applies a multi-dimensional universal least-favorable submodel parameterized
#' by epsVec. Then returns the k-dimensional vector of empirical means of the
#' efficient influence curve for each coordinate.
#'
#' @param Estimates A list of arms, each containing:
#'   - Hazards: a named list of hazard matrices (rows = times, columns = subjects)
#'   - EvntFreeSurv: a matrix of survival probabilities
#'   - NuisanceWeight: a matrix of 1/(g*(censorSurv)) or similar
#'   - PropScore: stored g-star?
#' @param epsVec numeric vector of length k = (#arms * #events * #times).
#' @param Data data.table with T.tilde, Delta, etc. (as in doTmleUpdate)
#' @param TargetEvent numeric vector of event types (e.g. c(1,2))
#' @param TargetTime numeric vector of target times (e.g. c(365, 730))
#' @param Verbose logical, if TRUE prints some progress
#'
#' @return A numeric vector of length k. This is the stack of Pn(EIC_j) for j=1..k.
#'   The root solver tries to make all of them ~0.
#'
#' @import data.table
submodelAndEIC <- function(Estimates, epsVec, Data,
                           TargetEvent, TargetTime,
                           Verbose = FALSE) {
  Time <- Event <- NULL

  # ------------------------------------------------------------------------
  # 1. Setup basic dimension info
  # ------------------------------------------------------------------------
  arms         <- names(Estimates)                     # e.g. c("A=0","A=1") or dynamic regimes
  nA           <- length(arms)
  nE           <- length(TargetEvent)
  nT           <- length(TargetTime)
  k_dim        <- nA * nE * nT
  EvalTimes    <- attr(Estimates, "Times")             # row index of hazards
  if (length(epsVec) != k_dim) {
    stop("epsVec length does not match #arms * #events * #times = ", k_dim)
  }

  # We'll map each coordinate i in 1..k_dim to (armIndex, eventIndex, timeIndex).
  coordMap <- expand.grid(
    arm   = seq_len(nA),        # 1..nA
    event = seq_len(nE),        # 1..nE
    time  = seq_len(nT)         # 1..nT
  )
  coordMap$eps <- epsVec

  # Pull out time/event from Data
  T.tilde <- Data[[attr(Data, "EventTime")]]
  Delta   <- Data[[attr(Data, "EventType")]]

  # ------------------------------------------------------------------------
  # 2. Copy the existing hazards so we can exponentiate them by epsVec
  # ------------------------------------------------------------------------
  # We'll store new hazards in a structure similar to "updateHazard",
  # but in a single pass for all coords. We'll later exponentiate them
  # based on sum-of-eps*CleverCov at each (row, col).
  UpdatedEstimates <- lapply(Estimates, function(estA) {
    # Deep-copy hazard list
    newHazList <- lapply(estA[["Hazards"]], function(H) {
      Hcopy <- H
      attr(Hcopy, "j") <- attr(H, "j")
      Hcopy
    })
    # We'll store them in a new list
    out <- list(
      Hazards         = newHazList,
      EvntFreeSurv    = estA[["EvntFreeSurv"]],    # will recalc below
      SummEIC         = data.table::copy(estA[["SummEIC"]]),
      IC              = estA[["IC"]],
      PropScore       = estA[["PropScore"]],
      NuisanceWeight  = estA[["NuisanceWeight"]]
    )
    out
  })
  names(UpdatedEstimates) <- arms
  attr(UpdatedEstimates, "Times") <- EvalTimes

  # We'll build an "exponent increment matrix" for each arm x event,
  # so that hazard_{event}(t) <- hazard_{event}(t)*exp( sum_of_... ).
  # We'll initialize them to all 0, then fill them in from each coordinate i.
  # hazardExponent[[armName]][[eventName]] is a matrix same size as hazards.
  hazardExponent <- list()
  for (aName in arms) {
    eventList <- UpdatedEstimates[[aName]][["Hazards"]]
    eventExponentList <- lapply(eventList, function(hazMatrix) {
      matrix(0, nrow = nrow(hazMatrix), ncol = ncol(hazMatrix))
    })
    hazardExponent[[aName]] <- eventExponentList
  }

  # ------------------------------------------------------------------------
  # 3. For each coordinate i, build the increment in the universal submodel
  # ------------------------------------------------------------------------
  for (i in seq_len(k_dim)) {
    rowi        <- coordMap[i, ]
    aIdx        <- rowi$arm
    eIdx        <- rowi$event
    tIdx        <- rowi$time
    eps_i       <- rowi$eps

    aName       <- arms[aIdx]
    thisEvent   <- TargetEvent[eIdx]
    thisTime    <- TargetTime[tIdx]

    # We replicate the logic from your single-step "updateHazard" function, but
    # only for the single coordinate i. Then we add it to hazardExponent[arm][event].
    # We'll do something like:
    #    F.j.t = cumsum(hazard_j(t)*survival)
    #    For times <= thisTime, define the "hFS" = ...
    #    Then multiply by getCleverCovariate(...).
    #    Then multiply by eps_i, and add it to hazardExponent.

    evStr <- as.character(thisEvent)
    # If there's no hazard in this arm for that event name, skip
    if (!evStr %in% names(hazardExponent[[aName]])) {
      if (Verbose) cat("Warning: event ", thisEvent,
                       " not found in arm ", aName, "\n")
      next
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~
    # (A) Build F.j.t = \int hazard_{j} * survival_{j} up to each time row
    Hazards_j  <- UpdatedEstimates[[aName]][["Hazards"]][[evStr]]
    Surv_a     <- UpdatedEstimates[[aName]][["EvntFreeSurv"]]
    F.j.t      <- apply(Hazards_j * Surv_a, 2, cumsum)

    # (B) Subset to times <= thisTime
    idx_le_time <- (EvalTimes <= thisTime)
    if (!any(idx_le_time)) {
      # if no times <= thisTime, skip
      next
    }

    # hFS: for rows <= thisTime
    # hFS(s,:) = [F.j.t(tau) - F.j.t(s)] / Surv(s), for s <= tau
    # with tau = thisTime.
    # Similar to your single-step code:
    FtauRow <- F.j.t[EvalTimes == thisTime, , drop=FALSE]  # 1 x ncol
    # Replicate that row over all s that are <= thisTime
    nrows_le  <- sum(idx_le_time)
    Ftau_mat  <- matrix(FtauRow, nrow = nrows_le, ncol = ncol(FtauRow), byrow = TRUE)
    # for the same indices
    hFS <- Ftau_mat - F.j.t[idx_le_time, , drop=FALSE]
    hFS <- hFS / Surv_a[idx_le_time, , drop=FALSE]

    # (C) getCleverCovariate(...) with LeqJ=1 if this hazard is indeed j.
    # We'll build a zero matrix for the entire row count, fill the rows <= time
    ClevCov <- matrix(0, nrow = nrow(Hazards_j), ncol = ncol(Hazards_j))
    subClev <- getCleverCovariate(
      GStar          = as.numeric(unlist(attr(UpdatedEstimates[[aName]][["PropScore"]], "g.star.obs"))),
      NuisanceWeight = UpdatedEstimates[[aName]][["NuisanceWeight"]][idx_le_time, , drop=FALSE],
      hFS            = hFS,
      LeqJ           = 1L   # because we're definitely event == thisEvent
    )
    # place subClev in the first sum(idx_le_time) rows
    ClevCov[idx_le_time, ] <- subClev

    # (D) Multiply by eps_i, and add to hazardExponent
    hazardExponent[[aName]][[evStr]] <- hazardExponent[[aName]][[evStr]] +
      (ClevCov * eps_i)
  }

  # ------------------------------------------------------------------------
  # 4. Exponentiate each hazard by exp( hazardExponent[..][..] ), then
  #    recompute survival, EIC, SummEIC for each arm.
  # ------------------------------------------------------------------------
  for (aName in arms) {
    for (evNm in names(hazardExponent[[aName]])) {
      oldHaz <- UpdatedEstimates[[aName]][["Hazards"]][[evNm]]
      incMat <- hazardExponent[[aName]][[evNm]]  # sum of all eps * cleverCov
      newHaz <- oldHaz * exp(incMat)
      attr(newHaz, "j") <- attr(oldHaz, "j")
      UpdatedEstimates[[aName]][["Hazards"]][[evNm]] <- newHaz
    }

    # Recompute survival for this arm
    # sum up all event hazards (not censoring) to get total event hazard
    # or if you have separate code for censoring hazard, you might exclude it from the sum
    # For example, if Delta=0 indicates censoring, you might skip that in the sum.
    # We'll do something like:
    eventHlist <- UpdatedEstimates[[aName]][["Hazards"]]
    # If your naming scheme is c("1","2","0") for events/censoring, skip "0" if you want event-only sum:
    # eventHlist <- eventHlist[names(eventHlist) %in% as.character(TargetEvent)]
    # Or if you do want to include all events in Surv calc:
    Summed <- Reduce(`+`, eventHlist)
    newSurv <- apply(Summed, 2, function(hz) exp(-cumsum(hz)))
    newSurv[newSurv < 1e-12 | is.na(newSurv) | is.nan(newSurv)] <- 1e-12
    UpdatedEstimates[[aName]][["EvntFreeSurv"]] <- newSurv

    # Now compute EIC via your existing getIC()
    newIC <- getIC(
      GStar          = attr(UpdatedEstimates[[aName]][["PropScore"]], "g.star.obs"),
      Hazards        = UpdatedEstimates[[aName]][["Hazards"]],
      TotalSurv      = UpdatedEstimates[[aName]][["EvntFreeSurv"]],
      NuisanceWeight = UpdatedEstimates[[aName]][["NuisanceWeight"]],
      TargetEvent    = TargetEvent,
      TargetTime     = TargetTime,
      T.tilde        = T.tilde,
      Delta          = Delta,
      EvalTimes      = EvalTimes,
      GComp          = FALSE
    )

    UpdatedEstimates[[aName]][["IC"]]      <- newIC
    UpdatedEstimates[[aName]][["SummEIC"]] <- summarizeIC(newIC)
  }

  # ------------------------------------------------------------------------
  # 5. Extract the k = (#arms * #events * #times) means of EIC
  #    We loop over i in 1..k_dim, matching how we enumerated epsVec,
  #    and read off SummEIC's "PnEIC".
  # ------------------------------------------------------------------------
  outVec <- numeric(k_dim)
  for (i in seq_len(k_dim)) {
    rowi       <- coordMap[i, ]
    aIdx       <- rowi$arm
    eIdx       <- rowi$event
    tIdx       <- rowi$time
    aName      <- arms[aIdx]
    thisEvent  <- TargetEvent[eIdx]
    thisTime   <- TargetTime[tIdx]

    # SummEIC is a data.table with columns: Time, Event, PnEIC, ...
    # We'll pick out the row that matches Time==thisTime, Event==thisEvent
    sE <- UpdatedEstimates[[aName]][["SummEIC"]]
    rowMatch <- sE[Time == thisTime & Event == thisEvent]
    if (nrow(rowMatch) < 1) {
      # no row found => set zero or NA
      outVec[i] <- 0.0
    } else {
      outVec[i] <- rowMatch$PnEIC
    }
  }

  if (Verbose) {
    cat("submodelAndEIC returning vector of length ", length(outVec),
        "; norm=", round(sqrt(sum(outVec^2)), 4), "\n")
  }
  return(outVec)
}
