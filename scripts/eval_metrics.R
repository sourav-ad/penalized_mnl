library(patchwork)

#Semi synthetic simulation evaluation metrics
res <- read.csv("data/500n_250iter_results.csv")

true_beta <- c(-0.4, -0.2, 0.3, 0.5, 0.3, -0.4)

#sign recovery
sign_rec <- colMeans(sign(res) == sign(true_beta))

#bias
bias <- colMeans(res) - true_beta

#rmse
rmse <- sqrt(colMeans((res - matrix(true_beta, nrow(res), 6, byrow = TRUE))^2))


metrics_table <- data.frame(
  Feature      = colnames(res),
  True         = true_beta,
  SignRecovery = sign_rec,
  Bias         = bias,
  RMSE         = rmse
)

metrics_table



library(tidyverse)

# Sample sizes
ns <- c(500, 600, 700, 800, 900, 973)

# True coefficients in correct order
true_beta <- c(
  ASC2_male      = -0.4,
  ASC3_age       = -0.2,
  cost_income    =  0.3,
  prot25_child   =  0.5,
  prot50_child   =  0.3,
  invasive_child = -0.4
)

# Function to compute metrics for one file
compute_metrics <- function(n) {
  file <- paste0("data/", n, "n_250iter_results.csv")
  res <- read.csv(file)
  
  # Ensure column order matches true_beta names
  res <- res[, names(true_beta)]
  
  # Sign recovery
  sign_rec <- colMeans(sign(res) == sign(true_beta))
  
  # Bias
  bias <- colMeans(res) - true_beta
  
  # RMSE
  rmse <- sqrt(colMeans((res - 
                           matrix(true_beta, nrow(res), length(true_beta), byrow = TRUE))^2))
  
  tibble(
    n = n,
    Feature = names(true_beta),
    SignRecovery = as.numeric(sign_rec),
    Bias = as.numeric(bias),
    RMSE = as.numeric(rmse)
  )
}

# Apply across all n values
all_metrics <- map_df(ns, compute_metrics)

# # Plot: Sign Recovery
# p_sign <- ggplot(all_metrics, aes(x = n, y = SignRecovery, color = Feature)) +
#   geom_line(size = 1.1) +
#   geom_point(color = "black") +
#   scale_x_continuous(breaks = ns) +
#   labs(y = "Sign Recovery", x = "Sample Size (n)") +
#   theme_bw()

# Plot: Bias
p_bias <- ggplot(all_metrics, aes(x = n, y = Bias, color = Feature)) +
  geom_line(size = 1.1) +
  geom_point(color = "black") +
  scale_x_continuous(breaks = ns) +
  labs(y = "Bias", x = "Sample Size (n)") +
  theme_bw()

# Plot: RMSE
p_rmse <- ggplot(all_metrics, aes(x = n, y = RMSE, color = Feature)) +
  geom_line(size = 1.1) +
  geom_point(color = "black") +
  scale_x_continuous(breaks = ns) +
  labs(y = "RMSE", x = "Sample Size (n)") +
  theme_bw()

# Display
#p_sign
p_bias + theme(legend.position = "none")
p_rmse + theme(legend.position = "none")

combined <- p_bias + p_rmse + 
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

combined
