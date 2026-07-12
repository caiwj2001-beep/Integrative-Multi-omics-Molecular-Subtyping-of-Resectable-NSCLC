# =============================================================
# config.R — shared configuration (portable)
# Edit PROJECT_ROOT to point to the directory that contains
# data/, results/, figures/, tables/. All other paths derive from it.
# =============================================================

# --- EDIT THIS ONE LINE to your local project root ---
PROJECT_ROOT <- Sys.getenv("NSCLC_PROJECT_ROOT",
                           unset = "C:/Users/Administrator/科研文件/NSCLC多组学分子分型")

# Optional: user R library (leave as-is or point to your library)
.user_lib <- Sys.getenv("NSCLC_R_LIB", unset = "")
if (nzchar(.user_lib)) .libPaths(c(.user_lib, .libPaths()))
options(stringsAsFactors = FALSE)

PROJ    <- PROJECT_ROOT
DATA    <- file.path(PROJ, "data")
MATDIR  <- file.path(DATA, "matrices")
TCGADIR <- file.path(DATA, "tcga")
RESULTS <- file.path(PROJ, "results")
FIGDIR  <- file.path(PROJ, "figures")
TABDIR  <- file.path(PROJ, "tables")
for (d in c(RESULTS, FIGDIR, TABDIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Integrated molecular subtype colors
IMS_COL <- c(IMS1 = "#1B5E9B", IMS2 = "#E8A73C", IMS3 = "#B23A48", IMS4 = "#3E8E7E")
PAL_MAIN <- "#003366"

suppressMessages(library(ggplot2))
theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, color = "grey90"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", size = base_size + 1),
      legend.key = element_blank()
    )
}

# Save a ggplot as PNG + TIFF at 300 dpi (submission requirement)
save_fig <- function(plot_obj, name, width = 7, height = 6, dpi = 300, res_dir = FIGDIR) {
  ggplot2::ggsave(file.path(res_dir, paste0(name, ".png")),  plot_obj,
                  width = width, height = height, dpi = dpi, bg = "white")
  ggplot2::ggsave(file.path(res_dir, paste0(name, ".tiff")), plot_obj,
                  width = width, height = height, dpi = dpi,
                  device = "tiff", compression = "lzw", bg = "white")
  cat(sprintf("Saved: %s (.png/.tiff)\n", name))
}

cat("R config loaded. PROJECT_ROOT =", PROJ, "\n")
