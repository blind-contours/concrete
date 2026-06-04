#!/usr/bin/env Rscript
# Generates the package hex logo at man/figures/logo.png (pkgdown auto-detects it).
#   Rscript scripts/make-logo.R
suppressWarnings(suppressMessages({
  library(ggplot2)
  library(hexSticker)
}))

# Motif: two covariate-adjusted cumulative-incidence curves (the package output)
t <- seq(0, 1, length.out = 100)
cif <- function(rate) 1 - exp(-rate * t^1.3)
d <- rbind(
  data.frame(t = t, y = cif(0.95), arm = "A=1"),
  data.frame(t = t, y = cif(0.65), arm = "A=0")
)

motif <- ggplot(d, aes(t, y, colour = arm)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = c("A=1" = "#2ec4c6", "A=0" = "#f25c54")) +
  theme_void() + theme(legend.position = "none")

sticker(
  motif,
  package = "concrete",
  p_size = 17, p_y = 1.45, p_color = "#ffffff",
  s_x = 1, s_y = 0.78, s_width = 1.5, s_height = 0.95,
  h_fill = "#2b3a55", h_color = "#2ec4c6",
  url = "blind-contours.github.io/concrete", u_size = 3.2, u_color = "#cfd8e3",
  filename = "man/figures/logo.png", dpi = 300
)
message("wrote man/figures/logo.png")
