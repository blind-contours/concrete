#' Fit candidate hazard learners
#'
#' @param Data data.table
#' @param Model list
#' @param CVFolds list
#' @param Hazards list
#' @param ReturnModels boolean
# #' @param HazFits list
# #' @param MinNuisance numeric (in the future a function)
# #' @param TargetEvent numeric vector
# #' @param TargetTime numeric vector
# #' @param Regime list
# #' @param Censored boolean
#' @import survival
#' @importFrom stats predict
#' @keywords internal
#'

getHazFit <- function(Data, Model, CVFolds, Hazards, ReturnModels) {
    IDCol <- attr(Data, "ID")
    TimeCol <- attr(Data, "EventTime")
    TypeCol <- attr(Data, "EventType")
    HazModel <- Model[which(names(Model) %in% unique(Data[[TypeCol]]))]

    SupLrnModel <- lapply(HazModel, function(ModelJ) {
        CandidateRisk <- stats::setNames(rep(0, length(ModelJ)), names(ModelJ))
        ModelFits <- vector("list", length(ModelJ))
        names(ModelFits) <- names(ModelJ)
        j <- attr(ModelJ, "j")

        for (Fold_v in CVFolds) {
            TrainIndices <- Fold_v[["training_set"]]
            ValidIndices <- Fold_v[["validation_set"]]
            TrainData <- Data[TrainIndices, .SD, .SDcols = setdiff(colnames(Data), IDCol)]
            ValidData <- Data[ValidIndices, .SD, .SDcols = setdiff(colnames(Data), IDCol)]

            for (i in seq_along(ModelJ)) {
                Fit <- try(
                    fitHazLearner(ModelSpec = ModelJ[[i]],
                                  Data = TrainData,
                                  j = j,
                                  TimeCol = TimeCol,
                                  TypeCol = TypeCol,
                                  TrtCol = attr(Data, "Treatment"),
                                  IDCol = IDCol,
                                  Hazards = Hazards),
                    silent = TRUE
                )
                if (inherits(Fit, "try-error")) {
                    CandidateRisk[i] <- Inf
                    next
                }
                if (ReturnModels) {
                    ModelFits[[i]][[length(ModelFits[[i]]) + 1L]] <- Fit
                }

                HazMat <- try(predictHazLearner(Fit, ValidData), silent = TRUE)
                if (inherits(HazMat, "try-error")) {
                    CandidateRisk[i] <- Inf
                    next
                }
                FoldRisk <- hazardValidationLoss(HazMat = HazMat,
                                                 ValidData = ValidData,
                                                 j = j,
                                                 TimeCol = TimeCol,
                                                 TypeCol = TypeCol,
                                                 EvalTimes = Fit[["Hazards"]][["Time"]])
                if (!is.finite(FoldRisk)) {
                    FoldRisk <- Inf
                }
                CandidateRisk[i] <- CandidateRisk[i] + FoldRisk
            }
        }

        if (all(!is.finite(CandidateRisk))) {
            stop("All hazard learner candidates failed for event type ", j, ".")
        }
        SLModel <- ModelJ[[which.min(CandidateRisk)]]
        SLCoef <- as.numeric(seq_along(CandidateRisk) == which.min(CandidateRisk))
        names(SLCoef) <- names(CandidateRisk)

        if (ReturnModels) {
            return(list("SupLrnCVRisks" = CandidateRisk,
                        "SupLrnModel" = SLModel,
                        "j" = j,
                        "SLCoef" = SLCoef,
                        "ModelFits" = ModelFits))
        } else {
            return(list("SupLrnCVRisks" = CandidateRisk,
                        "SupLrnModel" = SLModel,
                        "j" = j,
                        "SLCoef" = SLCoef))
        }
    })
    names(SupLrnModel) <- sapply(SupLrnModel, function(sl) sl[["j"]])

    HazFits <- lapply(SupLrnModel, function(SLMod) {
        HazFitOut <- fitHazLearner(ModelSpec = SLMod$SupLrnModel,
                                   Data = Data[, .SD, .SDcols = setdiff(colnames(Data), IDCol)],
                                   j = SLMod[["j"]],
                                   TimeCol = TimeCol,
                                   TypeCol = TypeCol,
                                   TrtCol = attr(Data, "Treatment"),
                                   IDCol = IDCol,
                                   Hazards = Hazards)
        attr(HazFitOut, "j") <- SLMod[["j"]]
        SLMod[["SupLrnModel"]] <- NULL
        attr(HazFitOut, "HazSL") <- SLMod
        return(HazFitOut)
    })
    names(HazFits) <- names(SupLrnModel)
    return(HazFits)
}

getHazSurvPred <- function(Data, HazFits, MinNuisance, TargetEvent, TargetTime, Regime,
                           Censored = NULL) {
    # Censored can be forced (e.g. during cross-fitting a validation fold may
    # contain no censored subjects even though the censoring hazard was fit).
    if (is.null(Censored)) Censored <- any(Data[[attr(Data, "EventType")]] <= 0)
    IDCol <- attr(Data, "ID")
    TrtCol <- attr(Data, "Treatment")
    TimeCol <- attr(Data, "EventTime")
    TypeCol <- attr(Data, "EventType")
    CovCols <- getHazCovCols(Data = Data,
                             TimeCol = TimeCol,
                             TypeCol = TypeCol,
                             IDCol = IDCol,
                             TrtCol = TrtCol)

    PredHazSurv <- lapply(Regime, function(Reg) {
        PredData <- as.data.table(Data)[, .SD, .SDcols = CovCols]
        setcolorder(PredData, neworder = CovCols)
        TrtNames <- colnames(Reg)
        PredData[, (TrtNames) := Reg[, .SD, .SDcols = TrtNames]]

        PredHaz <- lapply(HazFits, function(HazFit) {
            haz <- predictHazLearner(HazFit, PredData)
            attr(haz, "j") <- attr(HazFit, "j")
            return(haz)
        })
        names(PredHaz) <- names(HazFits)

        CensInd <- which(sapply(PredHaz, function(haz) attr(haz, "j") <= 0))
        HazInd <- setdiff(seq_along(PredHaz), CensInd)

        TotalSurv <- apply(Reduce(`+`, PredHaz[HazInd]), 2, function(haz) exp(-cumsum(haz)))
        TotalSurv[TotalSurv < 1e-12] <- 1e-12
        if (Censored && length(CensInd) >= 1L) {
            LaggedCensSurv <- apply(PredHaz[[CensInd]], 2, function(haz) c(1, utils::head(exp(-cumsum(haz)), -1)))
        } else {
            # no censoring hazard available (e.g. a cross-fitting training fold with
            # no censored subjects); fall back to no censoring weighting for these rows
            if (isTRUE(Censored))
                warning("No censoring hazard was fit for a prediction set; ",
                        "censoring weights set to 1 for those subjects.")
            LaggedCensSurv <- 1
        }
        PredHaz <- PredHaz[HazInd]

        Survival <- list("TotalSurv" = TotalSurv, "LaggedCensSurv" = LaggedCensSurv)
        return(list("Hazards" = PredHaz, "Survival" = Survival))
    })
    return(PredHazSurv)
}

fitHazLearner <- function(ModelSpec, Data, j, TimeCol, TypeCol, TrtCol, IDCol, Hazards) {
    if (inherits(ModelSpec, "Lrnr.Coxnet")) {
        Fit <- fitCoxnetHazLearner(Data = Data, j = j, TimeCol = TimeCol,
                                   TypeCol = TypeCol, TrtCol = TrtCol,
                                   IDCol = IDCol, Hazards = Hazards)
    } else if (inherits(ModelSpec, "Lrnr.RSF")) {
        Fit <- fitRsfHazLearner(Data = Data, j = j, TimeCol = TimeCol,
                                TypeCol = TypeCol, TrtCol = TrtCol,
                                IDCol = IDCol, Hazards = Hazards)
    } else if (inherits(ModelSpec, "Lrnr.Aareg")) {
        Fit <- fitAaregHazLearner(Data = Data, j = j, TimeCol = TimeCol,
                                  TypeCol = TypeCol, TrtCol = TrtCol,
                                  IDCol = IDCol, Hazards = Hazards)
    } else if (inherits(ModelSpec, "Lrnr.HAL")) {
        Fit <- fitHalHazLearner(Data = Data, j = j, TimeCol = TimeCol,
                                TypeCol = TypeCol, TrtCol = TrtCol,
                                IDCol = IDCol, Hazards = Hazards)
    } else {
        Fit <- fitCoxHazLearner(ModelSpec = ModelSpec, Data = Data, j = j,
                                TimeCol = TimeCol, TypeCol = TypeCol,
                                TrtCol = TrtCol, IDCol = IDCol,
                                Hazards = Hazards)
    }
    class(Fit) <- union("ConcreteHazFit", class(Fit))
    return(Fit)
}

predictHazLearner <- function(Fit, PredData) {
    Learner <- Fit[["Learner"]]
    if (identical(Learner, "cox")) {
        haz <- predictCoxHazLearner(Fit, PredData)
    } else if (identical(Learner, "coxnet")) {
        haz <- predictCoxnetHazLearner(Fit, PredData)
    } else if (identical(Learner, "rsf")) {
        haz <- predictRsfHazLearner(Fit, PredData)
    } else if (identical(Learner, "aareg")) {
        haz <- predictAaregHazLearner(Fit, PredData)
    } else if (identical(Learner, "hal")) {
        haz <- predictHalHazLearner(Fit, PredData)
    } else {
        stop("Unrecognized hazard learner fit.")
    }
    haz <- sanitizeHazardMatrix(haz, floor = 0)
    attr(haz, "j") <- attr(Fit, "j")
    return(haz)
}

getHazCovCols <- function(Data, TimeCol, TypeCol, IDCol, TrtCol) {
    c(TrtCol, setdiff(colnames(Data), c(TimeCol, TypeCol, TrtCol, IDCol)))
}

fitCoxHazLearner <- function(ModelSpec, Data, j, TimeCol, TypeCol, TrtCol, IDCol, Hazards) {
    FitData <- Data[, .SD, .SDcols = setdiff(colnames(Data), IDCol)]
    ModelFit <- do.call(survival::coxph, list("formula" = ModelSpec, "data" = FitData))
    BaseHazJ <- fitTreatmentBaseHazard(Data = FitData, j = j, TimeCol = TimeCol,
                                       TypeCol = TypeCol, TrtCol = TrtCol,
                                       Hazards = Hazards)
    Fit <- list("Learner" = "cox",
                "HazFit" = ModelFit,
                "BaseHaz" = BaseHazJ,
                "Hazards" = Hazards)
    attr(Fit, "j") <- j
    return(Fit)
}

fitCoxnetHazLearner <- function(Data, j, TimeCol, TypeCol, TrtCol, IDCol, Hazards) {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Coxnet hazard learning requires the 'glmnet' package.")
    }
    CovCols <- getHazCovCols(Data = Data, TimeCol = TimeCol, TypeCol = TypeCol,
                             IDCol = IDCol, TrtCol = TrtCol)
    x <- as.matrix(Data[, .SD, .SDcols = CovCols])
    y <- survival::Surv(time = Data[[TimeCol]], event = (Data[[TypeCol]] == j),
                        type = "right")
    nfolds <- max(3L, min(10L, floor(nrow(Data) / 10L)))
    ModelFit <- glmnet::cv.glmnet(x = x,
                                  y = y,
                                  family = "cox",
                                  nfolds = nfolds,
                                  penalty.factor = c(rep(0, length(TrtCol)),
                                                     rep(1, length(CovCols) - length(TrtCol))))
    BaseHazJ <- fitTreatmentBaseHazard(Data = Data, j = j, TimeCol = TimeCol,
                                       TypeCol = TypeCol, TrtCol = TrtCol,
                                       Hazards = Hazards)
    Fit <- list("Learner" = "coxnet",
                "HazFit" = ModelFit,
                "BaseHaz" = BaseHazJ,
                "CovCols" = CovCols,
                "Hazards" = Hazards)
    attr(Fit, "j") <- j
    return(Fit)
}

fitTreatmentBaseHazard <- function(Data, j, TimeCol, TypeCol, TrtCol, Hazards) {
    BaseHazCox <- paste0("Surv(time=", TimeCol, ", event=", TypeCol, "==", j, ")~",
                         paste0(TrtCol, collapse = "+"))
    BaseHazCox <- survival::coxph(stats::as.formula(BaseHazCox), data = Data)
    BaseHazJ <- rbind(data.table(time = 0, hazard = 0),
                      suppressWarnings(data.table::setDT(survival::basehaz(BaseHazCox,
                                                                            centered = TRUE))))
    colnames(BaseHazJ) <- c("Time", "BaseHaz")
    BaseHazJ <- merge(Hazards, BaseHazJ, by = "Time", all.x = TRUE)
    BaseHazJ <- BaseHazJ[order(Time)]
    CumBaseHaz <- zoo::na.locf(BaseHazJ[["BaseHaz"]], na.rm = FALSE)
    CumBaseHaz[is.na(CumBaseHaz)] <- 0
    BaseHazJ[, BaseHaz := c(0, diff(CumBaseHaz))]
    BaseHazJ[BaseHaz < 0 | is.na(BaseHaz) | !is.finite(BaseHaz), BaseHaz := 0]
    return(BaseHazJ)
}

predictCoxHazLearner <- function(Fit, PredData) {
    exp.coef <- stats::predict(Fit[["HazFit"]], newdata = PredData, type = "risk")
    sapply(as.numeric(exp.coef), function(expLP) Fit[["BaseHaz"]][["BaseHaz"]] * expLP)
}

predictCoxnetHazLearner <- function(Fit, PredData) {
    z <- as.matrix(PredData[, .SD, .SDcols = Fit[["CovCols"]]])
    exp.coef <- stats::predict(Fit[["HazFit"]], newx = z,
                               s = Fit[["HazFit"]]$lambda.min,
                               type = "response")
    sapply(as.numeric(exp.coef), function(expLP) Fit[["BaseHaz"]][["BaseHaz"]] * expLP)
}

fitRsfHazLearner <- function(Data, j, TimeCol, TypeCol, TrtCol, IDCol, Hazards) {
    if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
        stop("Random survival forest hazard learning requires the 'randomForestSRC' package.")
    }
    CovCols <- getHazCovCols(Data = Data, TimeCol = TimeCol, TypeCol = TypeCol,
                             IDCol = IDCol, TrtCol = TrtCol)
    FitData <- as.data.frame(Data[, .SD, .SDcols = c(TimeCol, TypeCol, CovCols)])
    FitData[[".event_j"]] <- as.integer(FitData[[TypeCol]] == j)
    FitData[[TypeCol]] <- NULL
    fit_formula <- stats::as.formula(paste0("survival::Surv(", TimeCol, ", .event_j) ~ ."))
    ntime <- max(10L, min(150L, length(unique(Data[[TimeCol]]))))
    ModelFit <- randomForestSRC::rfsrc(formula = fit_formula,
                                       data = FitData,
                                       ntree = 150L,
                                       ntime = ntime,
                                       forest = TRUE)
    Fit <- list("Learner" = "rsf",
                "HazFit" = ModelFit,
                "CovCols" = CovCols,
                "Hazards" = Hazards)
    attr(Fit, "j") <- j
    return(Fit)
}

predictRsfHazLearner <- function(Fit, PredData) {
    NewData <- as.data.frame(PredData[, .SD, .SDcols = Fit[["CovCols"]]])
    Pred <- stats::predict(Fit[["HazFit"]], newdata = NewData)
    chf <- Pred[["chf"]]
    if (is.null(dim(chf))) {
        chf <- matrix(chf, nrow = 1L)
    }
    EvalTimes <- Fit[["Hazards"]][["Time"]]
    CumHazBySubject <- t(apply(chf, 1, function(chf_i) {
        approxStepCumulative(x = Pred[["time.interest"]],
                             y = chf_i,
                             xout = EvalTimes)
    }))
    apply(CumHazBySubject, 1, function(ch) c(0, diff(ch)))
}

fitAaregHazLearner <- function(Data, j, TimeCol, TypeCol, TrtCol, IDCol, Hazards) {
    CovCols <- getHazCovCols(Data = Data, TimeCol = TimeCol, TypeCol = TypeCol,
                             IDCol = IDCol, TrtCol = TrtCol)
    FitData <- as.data.frame(Data[, .SD, .SDcols = c(TimeCol, TypeCol, CovCols)])
    FitData[[".event_j"]] <- as.integer(FitData[[TypeCol]] == j)
    FitData[[TypeCol]] <- NULL
    fit_formula <- stats::as.formula(paste0("survival::Surv(", TimeCol, ", .event_j) ~ ",
                                            paste0(CovCols, collapse = " + ")))
    ModelFit <- survival::aareg(formula = fit_formula, data = FitData)
    Fit <- list("Learner" = "aareg",
                "HazFit" = ModelFit,
                "Formula" = fit_formula,
                "CovCols" = CovCols,
                "Hazards" = Hazards)
    attr(Fit, "j") <- j
    return(Fit)
}

predictAaregHazLearner <- function(Fit, PredData) {
    Terms <- stats::delete.response(stats::terms(Fit[["Formula"]]))
    x <- stats::model.matrix(Terms, data = as.data.frame(PredData))
    colnames(x)[colnames(x) == "(Intercept)"] <- "Intercept"
    Coefs <- Fit[["HazFit"]][["coefficient"]]
    x <- x[, colnames(Coefs), drop = FALSE]
    CumHazFit <- Coefs %*% t(x)
    EvalTimes <- Fit[["Hazards"]][["Time"]]
    CumHaz <- apply(CumHazFit, 2, function(chf_i) {
        approxStepCumulative(x = Fit[["HazFit"]][["times"]],
                             y = chf_i,
                             xout = EvalTimes)
    })
    apply(CumHaz, 2, function(ch) c(0, diff(ch)))
}

fitHalHazLearner <- function(Data, j, TimeCol, TypeCol, TrtCol, IDCol, Hazards) {
    if (!requireNamespace("hal9001", quietly = TRUE)) {
        stop("HAL hazard learning requires the 'hal9001' package.")
    }
    CovCols <- getHazCovCols(Data = Data, TimeCol = TimeCol, TypeCol = TypeCol,
                             IDCol = IDCol, TrtCol = TrtCol)
    LongData <- makeHalHazardLongData(Data = Data,
                                      j = j,
                                      TimeCol = TimeCol,
                                      TypeCol = TypeCol,
                                      CovCols = CovCols,
                                      Hazards = Hazards)
    if (length(unique(LongData[[".Y"]])) < 2L) {
        stop("HAL hazard learner requires at least one event and one non-event interval.")
    }
    X <- as.data.frame(LongData[, .SD, .SDcols = c(".HazTime", ".IntervalWidth", CovCols)])
    ModelFit <- hal9001::fit_hal(X = X,
                                 Y = LongData[[".Y"]],
                                 family = "binomial",
                                 max_degree = 2,
                                 yolo = TRUE)
    Fit <- list("Learner" = "hal",
                "HazFit" = ModelFit,
                "CovCols" = CovCols,
                "Hazards" = Hazards,
                "TimeScale" = max(Hazards[["Time"]]))
    attr(Fit, "j") <- j
    return(Fit)
}

makeHalHazardLongData <- function(Data, j, TimeCol, TypeCol, CovCols, Hazards) {
    EvalTimes <- Hazards[["Time"]]
    Ends <- EvalTimes[-1L]
    Starts <- utils::head(EvalTimes, -1L)
    Widths <- pmax(Ends - Starts, .Machine$double.eps)
    TimeScale <- max(EvalTimes)
    MaxEval <- max(EvalTimes)

    LongData <- data.table::rbindlist(lapply(seq_len(nrow(Data)), function(i) {
        idx <- findInterval(min(Data[[TimeCol]][i], MaxEval), EvalTimes)
        if (idx < 2L) {
            return(NULL)
        }
        keep <- seq_len(idx - 1L)
        out <- Data[rep(i, length(keep)), .SD, .SDcols = CovCols]
        out[, `:=`(".HazTime" = Ends[keep] / TimeScale,
                   ".IntervalWidth" = Widths[keep] / TimeScale,
                   ".Y" = as.integer(Data[[TypeCol]][i] == j &
                                     Data[[TimeCol]][i] == Ends[keep]))]
        return(out)
    }), fill = TRUE)
    return(LongData)
}

predictHalHazLearner <- function(Fit, PredData) {
    EvalTimes <- Fit[["Hazards"]][["Time"]]
    Ends <- EvalTimes[-1L]
    Starts <- utils::head(EvalTimes, -1L)
    Widths <- pmax(Ends - Starts, .Machine$double.eps)
    TimeScale <- Fit[["TimeScale"]]
    CovCols <- Fit[["CovCols"]]
    n <- nrow(PredData)

    LongPred <- data.table::rbindlist(lapply(seq_len(n), function(i) {
        out <- PredData[rep(i, length(Ends)), .SD, .SDcols = CovCols]
        out[, `:=`(".HazTime" = Ends / TimeScale,
                   ".IntervalWidth" = Widths / TimeScale)]
        return(out)
    }), fill = TRUE)
    X <- as.data.frame(LongPred[, .SD, .SDcols = c(".HazTime", ".IntervalWidth", CovCols)])
    Prob <- stats::predict(Fit[["HazFit"]], new_data = X, type = "response")
    Prob <- pmin(pmax(as.numeric(Prob), 1e-12), 1 - 1e-8)
    Haz <- -log1p(-Prob)
    rbind(0, matrix(Haz, nrow = length(Ends), ncol = n))
}

hazardValidationLoss <- function(HazMat, ValidData, j, TimeCol, TypeCol, EvalTimes) {
    HazMat <- sanitizeHazardMatrix(HazMat, floor = 0)
    CumHaz <- apply(HazMat, 2, cumsum)
    if (is.null(dim(CumHaz))) {
        CumHaz <- matrix(CumHaz, ncol = 1L)
    }
    obs_idx <- findInterval(pmin(ValidData[[TimeCol]], max(EvalTimes)), EvalTimes)
    obs_idx[obs_idx < 1L] <- 1L
    obs_idx[obs_idx > length(EvalTimes)] <- length(EvalTimes)
    loss <- sum(CumHaz[cbind(obs_idx, seq_len(nrow(ValidData)))], na.rm = TRUE)

    event_rows <- which(ValidData[[TypeCol]] == j & ValidData[[TimeCol]] <= max(EvalTimes))
    if (length(event_rows)) {
        event_idx <- match(ValidData[[TimeCol]][event_rows], EvalTimes)
        keep <- !is.na(event_idx)
        if (any(keep)) {
            event_haz <- HazMat[cbind(event_idx[keep], event_rows[keep])]
            loss <- loss - sum(log(pmax(event_haz, 1e-12)))
        }
    }
    return(loss)
}

approxStepCumulative <- function(x, y, xout) {
    if (!length(x) || !length(y)) {
        return(rep(0, length(xout)))
    }
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]
    y <- y[keep]
    if (!length(x) || !length(y)) {
        return(rep(0, length(xout)))
    }
    ord <- order(x)
    x <- x[ord]
    y <- y[ord]
    x <- c(0, x)
    y <- c(0, y)
    keep <- !duplicated(x, fromLast = TRUE)
    x <- x[keep]
    y <- y[keep]
    stats::approx(x = x, y = y, xout = xout, method = "constant",
                  f = 0, yleft = 0, rule = 2)$y
}

sanitizeHazardMatrix <- function(haz, floor = 0, ceiling = Inf) {
    haz <- as.matrix(haz)
    if (nrow(haz) == 0L || ncol(haz) == 0L) {
        return(haz)
    }
    haz[!is.finite(haz)] <- floor
    haz <- pmax(haz, floor)
    if (is.finite(ceiling)) {
        haz <- pmin(haz, ceiling)
    }
    haz[1L, ] <- 0
    return(haz)
}
