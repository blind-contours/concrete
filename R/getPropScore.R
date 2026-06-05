#' getPropScore
#'
#' @param TrtVal numeric vector
#' @param CovDT data.table
#' @param TrtModel list or fitted object
#' @param MinNuisance numeric
#' @param Regime list
#' @param CVFolds list
#' @param TrtLoss character or function(A, g.A)
#' @param ReturnModels boolean
#' 
#' @import SuperLearner
#' @importFrom stats binomial gaussian
#' @keywords internal

getPropScore <- function(TrtVal, CovDT, TrtModel, MinNuisance, Regime,
                         CVFolds, TrtLoss = NULL, ReturnModels,
                         PredTrtVal = NULL, PredCovDT = NULL, PredRegime = NULL) {
    old <- options()
    on.exit(options(old))
    options(warn = 0)
    # Cross-fitting: when Pred* are supplied, the propensity model is FIT on
    # (TrtVal, CovDT) but the scores are PREDICTED on the held-out (PredTrtVal,
    # PredCovDT, PredRegime). Default (NULL) predicts on the fitting data.
    if (is.null(PredTrtVal)) PredTrtVal <- TrtVal
    if (is.null(PredCovDT)) PredCovDT <- CovDT
    if (is.null(PredRegime)) PredRegime <- Regime
    
    if (all(sapply(TrtModel, function(a) inherits(a, "SuperLearner")))) {
        TrtFit <- TrtModel
    } else {
        TrtFit <- vector("list", ncol(TrtVal))
        names(TrtFit) <- colnames(TrtVal)
        for (a_i in 1:ncol(TrtVal)) {
            if (attr(TrtModel[[a_i]], "Backend") == "SuperLearner") {
                SLArgs <- list()
                SLArgs[["Y"]] <- unlist(subset(TrtVal, select = a_i))
                if (a_i > 1) {
                    SLArgs[["X"]] <- cbind(subset(TrtVal, select = 1:(a_i - 1)), CovDT)
                } else {
                    SLArgs[["X"]] <- CovDT
                }
                SLArgs[["family"]] <- ifelse(length(unique(unlist(TrtVal))) == 2, "binomial", "gaussian")
                SLArgs[["SL.library"]] <- TrtModel[[a_i]]
                SLArgs[["cvControl"]] <- list("V" = as.integer(length(CVFolds)), "stratifyCV" = FALSE, 
                                              "shuffle" = FALSE,
                                              "validRows" = lapply(CVFolds, function(v) v[["validation_set"]]))
                TrtFit[[a_i]] <- do.call(SuperLearner, SLArgs)
            # } else if (attr(TrtModel[[a_i]], "Backend") == "sl3") {
            #     data <- as.data.frame(cbind(subset(TrtVal, select = 1:a_i), CovDT))
            #     TrtTask <- sl3::make_sl3_Task(
            #         data = data,
            #         covariates = setdiff(colnames(data), colnames(TrtVal)[a_i]),
            #         outcome = colnames(TrtVal)[a_i]
            #     )
            #     if (is.null(TrtModel[[a_i]]$params$learners)) {
            #         TrtSL <- sl3::Lrnr_cv$new(learner = TrtModel[[a_i]], folds = CVFolds)
            #     } else {
            #         TrtSL <- sl3::Lrnr_sl$new(learners = TrtModel[[a_i]], folds = CVFolds)
            #     }
            #     TrtFit[[a_i]] <- TrtSL$train(TrtTask)
            } else {
                stop("functionality for propensity score estimation not using 'sl3' or ", 
                     "'SuperLearner' has not yet been implemented")
            }
        } 
    }
    
    PropScores <- lapply(PredRegime, function(a) {
        if (!all(dim(a) == dim(PredTrtVal)))
            stop("Regime dimensions don't match with observed treatment. Bugfix needed")

        PropScore <- rep_len(1, nrow(PredTrtVal))
        for (a_i in 1:ncol(PredTrtVal)) {
            a_vec <- unlist(subset(a, select = a_i))
            if (!all(a_vec %in% c(0, 1)))
                stop("support for non-binary intervention variables is not yet implemented")

            if (attr(TrtModel[[a_i]], "Backend") == "SuperLearner") {
                if (a_i > 1) {
                    newdata <- cbind(subset(PredTrtVal, select = 1:(a_i - 1)), PredCovDT)
                } else {
                    newdata <- PredCovDT
                }
                g.a <- as.numeric(predict.SuperLearner(object = TrtFit[[a_i]], newdata = newdata)$pred)
                g.a[a_vec == 0] <- 1 - g.a[a_vec == 0]
                PropScore <- PropScore * g.a
            # } else if (attr(TrtModel[[a_i]], "Backend") == "sl3") {
            #     g.a <- unlist(TrtFit[[a_i]]$predict())
            #     g.a[a_vec == 0] <- 1 - g.a[a_vec == 0]
            #     PropScore <- PropScore * g.a
            }
        }
        attr(PropScore, "g.star.intervention") <- attr(a, "g.star")(a, PredCovDT, PropScore, a)
        attr(PropScore, "g.star.obs") <- attr(a, "g.star")(PredTrtVal, PredCovDT, PropScore, a)
        return(PropScore)
    })
    
    if (ReturnModels) {
        attr(PropScores, "TrtFit") <- TrtFit
    } else {
        attr(PropScores, "TrtFit") <- "TrtFits not saved because `ReturnModels' was set = FALSE"
    }
    attr(PropScores, "warnings") <- summary(warnings())
    return(PropScores)
}
