library(ggplot2)
library(dplyr)
library(stringr)

#Histogram distribution of optimal lambdas

input_dir <- "data"
output_dir <- "plots"
if(!dir.exists(output_dir)) dir.create(output_dir)

files <- list.files(input_dir, pattern = "_lambda_values.csv", full.names = TRUE)

for (f in files) {
  df <- read.csv(f)
  meta <- str_extract(basename(f), "\\d+n_\\d+iter")
  
  # extract numbers for nicer title
  n_val    <- str_extract(meta, "^[0-9]+")
  iter_val <- str_extract(meta, "(?<=_)\\d+")
  clean_title <- paste0("(n = ", n_val, ", iter = ", iter_val, ")")
  
  p <- ggplot(df, aes(x = BestLambda)) +
    geom_histogram(
      fill = "#269dc7", color = "black",
      bins = 20,        # adjust bins as needed
      boundary = 0,
      #size = 1.2,
      binwidth = 0.001
    ) +
    ggtitle(clean_title) +
    xlab("Best Lambda") +
    ylab("Frequency") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme(
      text            = element_text(size = 12, face = "plain", color = "black"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid       = element_blank(),
      axis.line        = element_line(color = "black", linewidth = 0.4),
      axis.ticks       = element_line(color = "black", linewidth = 0.4),
      axis.text        = element_text(size = 15, color = "black"),
      axis.title.x     = element_text(size = 20),
      axis.title.y     = element_text(size = 20),
      plot.title       = element_text(size = 26, hjust = 0.5, face = "plain")
    )
  
  ggsave(filename = file.path(output_dir, paste0(meta, "_hist.png")),
         plot = p, width = 7, height = 5, dpi = 300)
}


#Boxplot for detections 

files <- list.files(input_dir, pattern = "^[0-9]+n_[0-9]+iter_results\\.csv$", full.names = TRUE)

true_values <- data.frame(
  Feature = c("ASC2_male", "ASC3_age", "cost_income",
              "prot25_child", "prot50_child", "invasive_child"),
  Forced = c(-0.4, -0.2, 0.3, 0.5, 0.3, -0.4)
)

for (f in files) {
  results_df_interactions <- read.csv(f)
  base <- tools::file_path_sans_ext(basename(f))
  sample_size <- sub("^([0-9]+)n_.*$", "\\1", base)
  
  plot_data <- tidyr::pivot_longer(results_df_interactions,
                                   cols = everything(),
                                   names_to = "Feature",
                                   values_to = "Coefficient")
  
  p <- ggplot(plot_data, aes(x = Feature, y = Coefficient, fill = Feature)) +
    geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
    geom_point(data = true_values,
               aes(x = Feature, y = Forced),
               color = "red", shape = 95, size = 6, inherit.aes = FALSE) +
    labs(
      title = paste0("(n = ", sample_size, ")"),
      #subtitle = "Red dash = induced interaction value",
      x = "Interaction Feature",
      y = "Coefficient Estimate"
    ) +
    scale_fill_brewer(palette = "Pastel1") +
    theme(
      text = element_text(size = 15),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text        = element_text(size = 15, color = "black"),
      axis.title.x     = element_text(size = 20),
      axis.title.y     = element_text(size = 20),
      plot.title       = element_text(size = 26, hjust = 0.5, face = "plain"),
      legend.position = "none"
    )
  
  # Save with same base name in output_dir
  base <- tools::file_path_sans_ext(basename(f))
  outfile <- file.path(output_dir, paste0(base, ".png"))
  
  ggsave(outfile, plot = p, width = 8, height = 6, dpi = 400)
}


#Violin for detections: can detect multimodality

files <- list.files(input_dir, pattern = "^[0-9]+n_[0-9]+iter_results\\.csv$", full.names = TRUE)

true_values <- data.frame(
  Feature = c("ASC2_male", "ASC3_age", "cost_income",
              "prot25_child", "prot50_child", "invasive_child"),
  Forced = c(-0.4, -0.2, 0.3, 0.5, 0.3, -0.4)
)

for (f in files) {
  results_df_interactions <- read.csv(f)
  base <- tools::file_path_sans_ext(basename(f))
  sample_size <- sub("^([0-9]+)n_.*$", "\\1", base)
  
  plot_data <- tidyr::pivot_longer(results_df_interactions,
                                   cols = everything(),
                                   names_to = "Feature",
                                   values_to = "Coefficient")
  
  p <- ggplot(plot_data, aes(x = Feature, y = Coefficient, fill = Feature)) +
    geom_violin(trim = FALSE, alpha = 0.8, fill = "grey", colour = "black", width = 1.7) +
    stat_summary(fun = median, geom = "point", size = 2.5, colour = "black") +
    stat_summary(fun.data = ggplot2::median_hilow, geom = "errorbar",
                 width = 0.2, colour = "black") +
    geom_point(data = true_values,
               aes(x = Feature, y = Forced),
               color = "red", shape = 95, size = 6, inherit.aes = FALSE) +
    labs(
      title = paste0("(n = ", sample_size, ")"),
      #subtitle = "Red dash = induced interaction value",
      x = "Interaction Feature",
      y = "Coefficient Estimate"
    ) +
    scale_fill_brewer(palette = "Pastel1") +
    theme(
      text = element_text(size = 15),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text        = element_text(size = 15, color = "black"),
      axis.title.x     = element_text(size = 20),
      axis.title.y     = element_text(size = 20),
      plot.title       = element_text(size = 26, hjust = 0.5, face = "plain"),
      legend.position = "none"
    )
  
  # Save with same base name in output_dir
  base <- tools::file_path_sans_ext(basename(f))
  outfile <- file.path(output_dir, paste0(base, "violin.png"))
  
  ggsave(outfile, plot = p, width = 8, height = 6, dpi = 400)
}
