# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : plot_functions.R                                           ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.com                                 ║
# ║ Date          : 02-11-2026                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────

# Plotting helpers shared by the power and simulation-study figures, so that
# they look alike and are all written to Figures/ in the same way. They
# require {ggplot2}.

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Plotting helpers
# ─────────────────────────────────────────────────────────────────────────────

# Shared theme, so that the figures are visually consistent.
theme_usdt <- function(base_size = 9) {

  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(hjust = 0, face = "bold",
                                               size = base_size - 0.5),
      panel.border     = ggplot2::element_rect(fill = NA, colour = "grey80"),
      axis.title       = ggplot2::element_text(size = base_size - 0.5)
    )
}

# Single place where output format, units and resolution are set, so that all
# figures are written to Figures/ in the same way. file includes the subfolder
# of the document, e.g. "Manuscript/Figure_5.png" or "Supplemental/Figure_S1.png".
save_figure <- function(plot, file, width, height, dpi = 300) {

  path <- file.path("Figures", file)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = path, plot = plot,
                  width = width, height = height, units = "in", dpi = dpi)

  message("Saved: ", path)
  invisible(path)
}

# ─────────────────────────────────────────────────────────────────────────────
