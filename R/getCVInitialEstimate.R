#' Cross-fitted (CV-TMLE) initial nuisance estimation
#'
#' @description
#' Produces out-of-fold nuisance predictions for the propensity score and the
#' cause-specific / censoring hazards, so that no subject's nuisance estimate
#' uses that subject's own data. This is the cross-fitting (CV-TMLE / double
#' machine-learning) construction that weakens the conditions for valid
#' influence-function inference when flexible machine-learning learners are used,
#' which is the relevant regime for the FDA's caution around ML-based covariate
#' adjustment. The downstream targeting (`getEIC`, `doTmleUpdate`) is unchanged;
#' only how the per-subject nuisances are produced differs.
#'
#' For each outer fold the propensity Super Learner and the hazard library are
#' fit on the training rows (with an inner cross-validation for hazard-learner
#' selection) and predicted on the held-out validation rows; the per-subject
#' columns are then assembled into the same structure returned by
#' [getInitialEstimate()].
#'
#' @inheritParams getInitialEstimate
#' @keywords internal
#' @importFrom origami make_folds
getCVInitialEstimate <- function(Data, Model, CVFolds, MinNuisance, TargetEvent, TargetTime,
                                 Regime, ReturnModels, HazEnsemble = FALSE) {
    Time <- NULL
    n <- nrow(Data)
    Censored <- any(Data[[attr(Data, "EventType")]] <= 0)
    TimeVal <- Data[[attr(Data, "EventTime")]]
    TrtCol <- attr(Data, "Treatment")
    CovCol <- attr(Data, "CovNames")[["ColName"]]

    HazTimes <- sort(unique(c(TargetTime, TimeVal)))
    HazTimes <- HazTimes[HazTimes <= max(TargetTime)]
    Hazards <- data.table("Time" = c(0, HazTimes))
    EvalTimes <- Hazards[["Time"]]
    nT <- length(EvalTimes)
    arms <- names(Regime)

    copyAttrs <- function(sub, full) {
        for (an in c("EventTime", "EventType", "Treatment", "ID", "CovNames"))
            attr(sub, an) <- attr(full, an)
        sub
    }
    subRegime <- function(idx) lapply(Regime, function(r) r[idx, , drop = FALSE])

    PropFull <- stats::setNames(lapply(arms, function(a) numeric(n)), arms)
    GStarFull <- stats::setNames(lapply(arms, function(a) numeric(n)), arms)
    SurvFull <- stats::setNames(lapply(arms, function(a) matrix(NA_real_, nT, n)), arms)
    LagCensFull <- stats::setNames(lapply(arms, function(a) matrix(NA_real_, nT, n)), arms)
    HazFull <- stats::setNames(vector("list", length(arms)), arms)

    message("\nCross-fitting nuisance parameters over ", length(CVFolds), " folds:\n")
    for (v in seq_along(CVFolds)) {
        tr <- CVFolds[[v]][["training_set"]]
        va <- CVFolds[[v]][["validation_set"]]
        trainData <- copyAttrs(data.table::copy(Data[tr]), Data)
        validData <- copyAttrs(data.table::copy(Data[va]), Data)
        innerV <- max(2L, min(length(CVFolds), floor(length(tr) / 30)))
        innerFolds <- origami::make_folds(n = length(tr), V = innerV)

        ## propensity: fit on train, predict held-out validation rows
        PS <- getPropScore(
            TrtVal = trainData[, .SD, .SDcols = TrtCol],
            CovDT = subset(trainData, select = CovCol),
            TrtModel = Model[which(names(Model) %in% TrtCol)],
            MinNuisance = MinNuisance, Regime = subRegime(tr),
            CVFolds = innerFolds, TrtLoss = NULL, ReturnModels = FALSE,
            PredTrtVal = validData[, .SD, .SDcols = TrtCol],
            PredCovDT = subset(validData, select = CovCol),
            PredRegime = subRegime(va))

        ## hazards: fit (with inner CV selection) on train, predict validation rows
        HazFits <- getHazFit(Data = trainData, Model = Model, CVFolds = innerFolds,
                             Hazards = Hazards, ReturnModels = FALSE, Ensemble = HazEnsemble)
        HSP <- getHazSurvPred(Data = validData, HazFits = HazFits, MinNuisance = MinNuisance,
                              TargetEvent = TargetEvent, TargetTime = TargetTime,
                              Regime = subRegime(va), Censored = Censored)

        for (a in arms) {
            PropFull[[a]][va] <- as.numeric(PS[[a]])
            GStarFull[[a]][va] <- as.numeric(unlist(attr(PS[[a]], "g.star.obs")))
            SurvFull[[a]][, va] <- HSP[[a]][["Survival"]][["TotalSurv"]]
            if (Censored) LagCensFull[[a]][, va] <- HSP[[a]][["Survival"]][["LaggedCensSurv"]]
            if (is.null(HazFull[[a]])) HazFull[[a]] <- list()
            for (j in names(HSP[[a]][["Hazards"]])) {
                if (is.null(HazFull[[a]][[j]])) {
                    HazFull[[a]][[j]] <- matrix(NA_real_, nT, n)
                    attr(HazFull[[a]][[j]], "j") <- attr(HSP[[a]][["Hazards"]][[j]], "j")
                }
                HazFull[[a]][[j]][, va] <- HSP[[a]][["Hazards"]][[j]]
            }
        }
    }

    ## Guard against an event type absent from a training fold (rare with the
    ## EventType-stratified outer folds), which would leave NA hazard columns:
    ## treat the unseen event's hazard as 0 for those subjects rather than NA.
    naFilled <- FALSE
    for (a in arms) for (j in names(HazFull[[a]])) {
        if (anyNA(HazFull[[a]][[j]])) {
            HazFull[[a]][[j]][is.na(HazFull[[a]][[j]])] <- 0
            naFilled <- TRUE
        }
    }
    if (naFilled)
        warning("Some event types were absent from a cross-fitting training fold; ",
                "their hazards were set to 0 for the affected subjects. Consider ",
                "fewer folds or pooling rare competing events.")

    InitialEstimates <- lapply(arms, function(a) {
        PropScore <- PropFull[[a]]
        attr(PropScore, "g.star.obs") <- GStarFull[[a]]
        if (Censored) {
            NuisanceDenom <- sapply(seq_len(n), function(i) PropScore[i] * LagCensFull[[a]][, i])
        } else {
            NuisanceDenom <- matrix(PropScore, nrow = nT, ncol = n, byrow = TRUE)
        }
        NuisanceWeight <- 1 / truncNuisanceWeight(NuisanceDenom = NuisanceDenom,
                                                  MinNuisance = MinNuisance, RegimeName = a)
        list("PropScore" = PropScore,
             "Hazards" = HazFull[[a]],
             "EvntFreeSurv" = SurvFull[[a]],
             "NuisanceWeight" = NuisanceWeight)
    })
    names(InitialEstimates) <- arms
    attr(InitialEstimates, "Times") <- EvalTimes
    attr(InitialEstimates, "InitFits") <- "Cross-fitted nuisances (per-fold fits not stored)"
    return(InitialEstimates)
}
