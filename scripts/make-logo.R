#!/usr/bin/env Rscript
# Generates the package hex logo at man/figures/logo.png (pkgdown auto-detects it).
#   Rscript scripts/make-logo.R
#
# Motif: two covariate-adjusted cumulative-incidence curves (active vs control)
# with a confidence ribbon and a target-time marker -- the package's core output.
# Data are confined to the central region of the panel so no ink touches the
# hex border.
suppressWarnings(suppressMessages({
  library(ggplot2)
  library(hexSticker)
}))

t  <- seq(0, 1, length.out = 220)
Fa <- 0.55 * (1 - exp(-1.7 * t^1.25))     # active arm cumulative incidence
Fc <- 0.55 * (1 - exp(-2.6 * t^1.25))     # control arm (higher risk)
band <- 0.05 + 0.06 * t                    # widening CI half-width

curves <- rbind(
  data.frame(t = t, y = Fa, arm = "active"),
  data.frame(t = t, y = Fc, arm = "control")
)
ribbon <- rbind(
  data.frame(t = t, lo = Fa - band, hi = Fa + band, arm = "active"),
  data.frame(t = t, lo = Fc - band, hi = Fc + band, arm = "control")
)
tt <- 0.78
marks <- data.frame(t = tt,
                    y = c(0.55 * (1 - exp(-1.7 * tt^1.25)),
                          0.55 * (1 - exp(-2.6 * tt^1.25))),
                    arm = c("active", "control"))

pal <- c(active = "#34e7c8", control = "#ff8364")

motif <- ggplot() +
  geom_ribbon(data = ribbon, aes(t, ymin = lo, ymax = hi, fill = arm), alpha = 0.38) +
  geom_vline(xintercept = tt, linetype = "22", colour = "#bfe0da", linewidth = 0.45) +
  geom_line(data = curves, aes(t, y, colour = arm), linewidth = 1.7, lineend = "round") +
  geom_point(data = marks, aes(t, y), colour = "#15323f", fill = "#f3f7f6",
             shape = 21, size = 2.2, stroke = 0.7) +
  scale_colour_manual(values = pal) + scale_fill_manual(values = pal) +
  # pad the limits so the curves occupy the central ~70% of the panel
  coord_cartesian(xlim = c(-0.18, 1.18), ylim = c(-0.12, 0.92), expand = FALSE) +
  theme_void() + theme(legend.position = "none")

sticker(
  motif,
  package = "concrete",
  p_size = 16, p_y = 1.5, p_color = "#f3f7f6", p_family = "sans", p_fontface = "bold",
  s_x = 1, s_y = 0.86, s_width = 1.45, s_height = 1.08,
  h_fill = "#15323f", h_color = "#34e7c8", h_size = 1.5,
  url = "blind-contours.github.io/concrete", u_size = 3.0, u_color = "#9fc7c0",
  white_around_sticker = FALSE,
  filename = "man/figures/logo.png", dpi = 300
)
message("wrote man/figures/logo.png")
