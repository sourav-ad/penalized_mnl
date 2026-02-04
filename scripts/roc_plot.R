## TRUE ROC PIPELINE (penalized alpha=0.5,1 + unpenalized vanilla)
## - Reads your *_interaction_detection_long.csv files
## - Drops the redundant Threshold dimension (keeps one Estimate per iter×interaction)
## - Uses score = |Estimate|, label = 1{TrueBeta != 0}
## - Computes ROC by sweeping all unique score cutpoints (pROC does this)
## - Outputs: combined ROC plot + AUC table
##
## Requirements: dplyr, ggplot2, pROC, readr (optional)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(pROC)
})

# -----------------------------
# 0) Paths
# -----------------------------
n_val <- 973
iterations <- 100

pen_dir <- file.path("data", "v3", paste0("n_", n_val, "_iter_100"))
van_dir <- file.path("data", "v3_vanilla", paste0("n_", n_val, "_iter_100"))

pen_files <- c(
  "0.5" = file.path(pen_dir, paste0(n_val, "n_", iterations, "_alpha_", 0.5, "_interaction_detection_long.csv")),
  "1"   = file.path(pen_dir, paste0(n_val, "n_", iterations, "_alpha_", 1,   "_interaction_detection_long.csv"))
)

# pick ONE vanilla run (change the filename to whichever exists)
# Example placeholders; set this to your actual vanilla filename.
van_file <- file.path(van_dir, paste0(n_val, "n_100_alpha_1_interaction_detection_long.csv"))

# -----------------------------
# 1) Read + normalize one long file to "base" rows
#    (one row per iteration × interaction)
# -----------------------------
read_long_to_base <- function(path, method_label) {
  df <- read.csv(path)

  # Expected cols: Alpha, Threshold, Iteration, Interaction, TrueBeta, Estimate, Detected
  stopifnot(all(c("Iteration", "Interaction", "TrueBeta", "Estimate") %in% names(df)))

  base <- df %>%
    mutate(
      Method = method_label,
      y = as.integer(TrueBeta != 0),
      score = abs(Estimate)
    ) %>%
    group_by(Method, Iteration, Interaction) %>%
    summarise(
      y = first(y),
      score = first(score),
      .groups = "drop"
    )

  base
}

# -----------------------------
# 2) Build base data for penalized methods
# -----------------------------
base_pen <- bind_rows(
  read_long_to_base(pen_files[["0.5"]], "Penalized (alpha=0.5)"),
  read_long_to_base(pen_files[["1"]],   "Penalized (alpha=1)")
)

# -----------------------------
# 3) Build base data for vanilla method
# -----------------------------
# If your vanilla long file has no Alpha column, that's fine; we ignore it anyway.
base_van <- read_long_to_base(van_file, "Unpenalized (vanilla)")

# Combine
base_all <- bind_rows(base_pen, base_van)

# Sanity checks
# You want both classes present (active and inactive) for ROC.
class_check <- base_all %>%
  group_by(Method) %>%
  summarise(n_active = sum(y == 1), n_inactive = sum(y == 0), .groups = "drop")
print(class_check)

# -----------------------------
# 4) Compute TRUE ROC + AUC per Method
#    (pROC sweeps all unique score cutpoints)
# -----------------------------
roc_df <- base_all %>%
  group_by(Method) %>%
  group_modify(~{
    roc_obj <- pROC::roc(
      response  = .x$y,
      predictor = .x$score,
      direction = "<",
      quiet = TRUE
    )
    tibble(
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities,
      AUC = as.numeric(pROC::auc(roc_obj))
    )
  }) %>%
  ungroup()

auc_tab <- roc_df %>%
  group_by(Method) %>%
  summarise(AUC = first(AUC), .groups = "drop") %>%
  arrange(desc(AUC))

print(auc_tab)

# -----------------------------
# 5) Plot combined TRUE ROC (step curves) + AUC in legend labels
# -----------------------------
auc_labels <- auc_tab %>%
  mutate(MethodLabel = paste0(Method, " | AUC=", sprintf("%.3f", AUC)))

roc_plot_df <- roc_df %>%
  left_join(auc_labels %>% select(Method, MethodLabel), by = "Method")

# Preserve ordering in legend: penalized then vanilla (or by AUC)
roc_plot_df$MethodLabel <- factor(
  roc_plot_df$MethodLabel,
  levels = auc_labels$MethodLabel
)

roc_plot_df <- roc_plot_df %>%
  mutate(
    MethodType = case_when(
      grepl("Unpenalized", Method) ~ "Unpenalized",
      grepl("alpha=0.5", Method)   ~ "Penalized",
      grepl("alpha=1", Method)     ~ "Penalized"
    ),
    AlphaLine = case_when(
      grepl("alpha=0.5", Method) ~ "alpha = 0.5",
      grepl("alpha=1", Method)   ~ "alpha = 1",
      grepl("Unpenalized", Method) ~ "unpenalized"
    )
  )


p <- ggplot(
  roc_plot_df,
  aes(x = FPR, y = TPR,
      color = MethodType,
      linetype = AlphaLine)
) +
  geom_step(linewidth = 1.1) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_color_manual(
    values = c(
      "Penalized"   = "red",
      "Unpenalized" = "#279FF5"
    ),
    name = "Method"
  ) +
  scale_linetype_manual(
    values = c(
      "alpha = 0.5"   = "dotted",
      "alpha = 1"     = "solid",
      "unpenalized"   = "dotdash"
    ),
    name = "Penalty"
  ) +
  labs(
    x = "FPR",
    y = "TPR",
    title = paste0("ROC (n = ", n_val, ")")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x  = element_text(size = 15),
    axis.text.y  = element_text(size = 15),
    plot.title   = element_text(size = 26)
  )

print(p)


ggsave(filename = file.path("plots", paste0("ROC_n", n_val, "_iter", iterations, ".png")),
       plot = p, width = 9, height = 6, dpi = 300)

#write.csv(auc_tab, file.path("data", "v3",
#  paste0("trueROC_AUC_n", n_val, "_iter", iterations, ".csv")), row.names = FALSE)
