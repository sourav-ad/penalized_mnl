# Threshold-based screening plot
# Penalized (alpha = 0.5, alpha = 1) + Unpenalized (vanilla)


suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

#set accordingly
n_val      <- 50
iterations <- 100

pen_dir <- file.path("data", "v3", paste0("n_", n_val, "_iter_", iterations))
van_dir <- file.path("data", "v3_vanilla", paste0("n_", n_val, "_iter_", iterations))

alphas <- c(0.5, 1)


#Read in penalized confusion rates

read_penalized_rates <- function(alpha_val) {
  fn <- file.path(
    pen_dir,
    paste0(n_val, "n_", iterations, "_alpha_", alpha_val, "_confusion_rates.csv")
  )

  read.csv(fn) %>%
    mutate(
      Model = paste0("pen_alpha_", alpha_val),
      Threshold = as.numeric(as.character(Threshold))
    )
}

pen_df <- bind_rows(lapply(alphas, read_penalized_rates))

#Read unpenalized (vanilla)
#(stored with alpha=1/0.5 in filename but actually unpenalized, just 2 iterations of the same thing)

van_fn <- file.path(
  van_dir,
  paste0(n_val, "n_", iterations, "_alpha_1_confusion_rates.csv")
)

van_df <- read.csv(van_fn) %>%
  mutate(
    Model = "unpenalized",
    Threshold = as.numeric(as.character(Threshold))
  )

#Combine
rates_df <- bind_rows(pen_df, van_df) %>%
  mutate(
    ThresholdLab = factor(Threshold, levels = sort(unique(Threshold))),
    ShapeGroup = factor(
      Model,
      levels = c("pen_alpha_0.5", "pen_alpha_1", "unpenalized")
    )
  )

#Plot

p <- ggplot(
  rates_df,
  aes(x = FPR, y = TPR, color = ThresholdLab, shape = ShapeGroup)
) +
  geom_point(size = 3.8,
             #position = position_jitter(width = 0.01, height = 0.01),
             alpha = 0.9) +

  scale_shape_manual(
    values = c(
      "pen_alpha_0.5" = 16,  # dot
      "pen_alpha_1"   = 4,   # cross
      "unpenalized"   = 17   # triangle
    ),
    labels = c(
      "Alpha = 0.5",
      "Alpha = 1",
      "Unpenalized"
    ),
    name = "Model"
  ) +

  scale_color_brewer(palette = "Dark2", name = "Threshold") +

  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +

  labs(
    x = "FPR",
    y = "TPR",
    title = paste0("n = ", n_val)
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

#save
ggsave(filename = file.path("plots", paste0("fprtpr_all", n_val, "_iter", iterations, ".png")),
        plot = p, width = 10, height = 7, dpi = 350)

#Sanity checks
#table(rates_df$Model)
# head(rates_df)
