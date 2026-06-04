#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# make-doc-figures.R
#
# Generates the *computed* figures embedded in the README and vignettes from a
# small, fast `concrete` fit. Run once after changing the plotting code or the
# example; the resulting PNGs are committed so the vignettes stay `eval = FALSE`
# and the CRAN/vignette build does no model fitting.
#
#   Rscript scripts/make-doc-figures.R
#
# Output: vignettes/figures/*.png and man/figures/*.png (README hero figures).
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  library(concrete)
  library(data.table)
  library(ggplot2)
}))

set.seed(20260601)

fig_dir <- "vignettes/figures"
man_dir <- "man/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(man_dir, showWarnings = FALSE, recursive = TRUE)

save_png <- function(plot, file, width = 7, height = 4, dpi = 150, also_man = FALSE) {
  ggplot2::ggsave(file.path(fig_dir, file), plot = plot, width = width,
                  height = height, dpi = dpi, bg = "white")
  if (also_man) {
    ggplot2::ggsave(file.path(man_dir, file), plot = plot, width = width,
                    height = height, dpi = dpi, bg = "white")
  }
  message("wrote ", file)
}

# --- Small competing-risks trial example (mirrors the installed smoke test) ---
trial <- as.data.table(survival::pbc)
trial <- trial[!is.na(trt), .(id, time, status, trt, age, sex, albumin, bili)]
trial <- trial[stats::complete.cases(trial)]
trial[, arm := as.integer(trt == 2)]
# 0 = censored, 1 = death (event of interest), 2 = transplant (competing event)
trial[, event := data.table::fifelse(status == 2L, 1L,
        data.table::fifelse(status == 1L, 2L, 0L))]
trial <- trial[, .(id, time, event, arm, age, sex, albumin, bili)]

TargetTime <- c(365, 730, 1460, 2190)

# A converged fit uses the adaptive update with the absolute stopping rule,
# which is the recommended setting for rare-event / competing-risk targets.
abs_tol <- 0.02 / sqrt(nrow(trial))

fit_args <- function(events, rule = "absolute") {
  formatArguments(
    DataTable     = trial,
    EventTime     = "time",
    EventType     = "event",
    Treatment     = "arm",
    ID            = "id",
    Intervention  = makeITT(),
    TargetTime    = TargetTime,
    TargetEvent   = events,
    CVArg         = list(V = 2),
    UpdateMethod  = "adaptive",
    EICStopRule   = rule,
    EICStopAbsTol = if (rule == "relative") 0 else abs_tol,
    MaxUpdateIter = 50,
    Verbose       = FALSE
  )
}

# === (A) Quickstart: single event of interest, converged ====================
message("Fitting (A) quickstart single-event model ...")
EstA <- doConcrete(fit_args(events = 1))
OutA <- getOutput(EstA, Estimand = c("Risk", "RD", "RR"),
                  Intervention = c(1, 2), GComp = TRUE, Simultaneous = TRUE)
figsA <- plot(OutA, ask = FALSE, NullLine = TRUE)

save_png(figsA[["risk"]] +
           labs(subtitle = "Adjusted cumulative incidence of the event of interest, PBC example (n = 276)"),
         "quickstart-risk.png", width = 7.5, height = 3.8, also_man = TRUE)
save_png(figsA[["risk"]], "readme-hero.png", width = 7.5, height = 3.4, also_man = TRUE)
save_png(figsA[["rd"]] +
           labs(subtitle = "Active minus control, PBC example"),
         "quickstart-rd.png", width = 7.5, height = 3.8)

out_tbl <- OutA[, .(Time, Event, Estimand, Intervention, Estimator,
    `Pt Est` = round(`Pt Est`, 3), se = round(se, 3),
    `CI Low` = round(`CI Low`, 3), `CI Hi` = round(`CI Hi`, 3))]
data.table::fwrite(out_tbl, file.path(fig_dir, "quickstart-output.csv"))
saveRDS(out_tbl, file.path(fig_dir, "quickstart-output.rds"))
saveRDS(getTmleDiagnostics(EstA, type = "components"),
        file.path(fig_dir, "quickstart-components.rds"))

# === (B) Competing risks: both events, for the "how it works" illustration ==
message("Fitting (B) competing-risks model ...")
EstB <- doConcrete(fit_args(events = c(1, 2)))
OutB <- getOutput(EstB, Estimand = "Risk", Intervention = c(1, 2),
                  GComp = FALSE, Simultaneous = TRUE)
figsB <- plot(OutB, ask = FALSE)
save_png(figsB[["risk"]] +
           labs(subtitle = "Event 1 = death, Event 2 = transplant (competing), PBC example"),
         "competing-risks.png", width = 8, height = 3.8)

png(file.path(fig_dir, "diagnostics-convergence.png"), width = 7, height = 4,
    units = "in", res = 150, bg = "white")
plot(EstB, convergence = TRUE, gweights = FALSE, ask = FALSE)
dev.off()
message("wrote diagnostics-convergence.png")

png(file.path(fig_dir, "diagnostics-nuisance-weights.png"), width = 7, height = 4,
    units = "in", res = 150, bg = "white")
plot(EstB, convergence = FALSE, gweights = TRUE, ask = FALSE)
dev.off()
message("wrote diagnostics-nuisance-weights.png")

# === (C) Convergence teaching: relative rule struggles on the rare event ====
message("Fitting (C) relative-rule model (rare-event stress) ...")
EstC <- doConcrete(fit_args(events = c(1, 2), rule = "relative"))
saveRDS(getTmleDiagnostics(EstC, type = "components"),
        file.path(fig_dir, "convergence-relative-components.rds"))
saveRDS(getTmleDiagnostics(EstB, type = "components"),
        file.path(fig_dir, "convergence-absolute-components.rds"))

message("Done. Documentation figures written to ", fig_dir, " and ", man_dir, ".")
