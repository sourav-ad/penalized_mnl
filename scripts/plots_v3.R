## =========================================
## POWER CURVES: y = 1-FNR vs |TrueBeta|
##  - Reads:
##    data/n_<n>_iter_<iter>/alpha_<alpha>_interaction_detection_long.csv
##  - Power is computed ONLY on induced cases (TrueBeta != 0)
##  - Lines: alpha 0.5 dotted, alpha 1 solid
##  - Facets: Threshold (since power varies by threshold)
## =========================================

library(dplyr)
library(ggplot2)

n_val <- 250
iter  <- 100
alphas <- c(0.5, 1)
bin_width <- 0.1  # set 0 to skip binning

read_det <- function(a) {
  fn <- file.path("data", paste0("n_", n_val, "_iter_", iter),
                  paste0("alpha_", a, "_interaction_detection_long.csv"))
  read.csv(fn) %>% mutate(Alpha = a)
}

dl <- bind_rows(lapply(alphas, read_det))

curve_df <- dl %>%
  mutate(
    BetaAbs = abs(TrueBeta),
    BetaBin = if (bin_width > 0) round(BetaAbs / bin_width) * bin_width else BetaAbs
  ) %>%
  filter(BetaAbs > 0) %>%                                  # POWER = P(Detected | induced)
  group_by(Alpha, Threshold, BetaBin) %>%
  summarise(Power = mean(Detected), .groups = "drop")

curve_df$Alpha <- factor(curve_df$Alpha, levels = c(0.5, 1))

p <- ggplot(curve_df, aes(x = BetaBin, y = Power, linetype = Alpha)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Threshold, scales = "free_x") +
  scale_linetype_manual(values = c("0.5" = "dotted", "1" = "solid")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "|TrueBeta|", y = "Power (1 - FNR)", linetype = "Alpha") +
  theme_bw()

print(p)

ggsave(
  filename = file.path("plots", paste0("n_", n_val, "_iter_", iter, "_power_curve_by_threshold.png")),
  plot = p, width = 9, height = 5.5, dpi = 300
)


#lambdas
library(ggplot2)
library(scales)

n_val  <- 50
iter   <- 100
alpha  <- 1

infile <- file.path("data",
                    paste0("n_", n_val, "_iter_", iter),
                    paste0(n_val, "n_", iter, "_alpha_", alpha, "_optimal_lambdas.csv"))

df <- read.csv(infile)

## adjust this: smaller = more bins, larger = fewer bins
binwidth_log10 <- 0.25

df$log10_lambda <- log10(df$BestLambda)

p <- ggplot(df, aes(x = log10_lambda)) +
  geom_histogram(binwidth = binwidth_log10, boundary = 0, closed = "left") +
  scale_x_continuous(
    name = "Optimal lambda (log scale, labeled in lambda units)",
    labels = x
  ) +
  ylab("Count") +
  ggtitle(paste0("Optimal lambda histogram | n=", n_val, " iter=", iter, " alpha=", alpha)) +
  theme_bw()

print(p)

outfile <- file.path("plots",
                     paste0("n_", n_val, "_iter_", iter, "_alpha_", alpha,
                            "_lambda_hist_binw_", binwidth_log10, ".png"))

ggsave(outfile, plot = p, width = 8, height = 5.5, dpi = 300)



###Post hoc confusion matrices
library(dplyr)

## --------------------------
## inputs
## --------------------------
n_persons  <- 973
iterations <- 100
alpha      <- 0.5

th_grid <- c(0.1)

out_dir <- file.path("data/v3_1", paste0("n_", n_persons, "_iter_", iterations))

coef_file <- file.path(
  out_dir,
  paste0(n_persons, "n_", iterations, "_alpha_", alpha, "_interaction_coefs.csv")
)

true_file <- file.path(
  out_dir,
  paste0(n_persons, "n_", iterations, "_alpha_", alpha, "_true_interactions.csv")
)

## --------------------------
## load data
## --------------------------
coef_df <- read.csv(coef_file)
true_df <- read.csv(true_file)

## drop meta columns
coef_mat <- coef_df %>% select(-Alpha, -Iteration)
true_mat <- true_df %>% select(-Alpha, -Iteration)

interaction_names <- colnames(coef_mat)
K <- ncol(coef_mat)

## --------------------------
## build detection-long for all thresholds
## --------------------------
det_long_all <- vector("list", length(th_grid))
cm_all       <- vector("list", length(th_grid))
rates_all    <- vector("list", length(th_grid))

for (i in seq_along(th_grid)) {

  th <- th_grid[i]

  detected_mat <- abs(as.matrix(coef_mat)) >= th
  truth_mat    <- as.matrix(true_mat) != 0

  ## detection-long
  det_long_all[[i]] <- data.frame(
    Alpha       = alpha,
    Threshold   = th,
    Iteration   = rep(seq_len(iterations), each = K),
    Interaction = rep(interaction_names, times = iterations),
    TrueBeta    = as.vector(t(as.matrix(true_mat))),
    Estimate    = as.vector(t(as.matrix(coef_mat))),
    Detected    = as.integer(as.vector(t(detected_mat)))
  )

  ## confusion counts
  TP <- sum(detected_mat &  truth_mat)
  FN <- sum(!detected_mat &  truth_mat)
  FP <- sum(detected_mat & !truth_mat)
  TN <- sum(!detected_mat & !truth_mat)

  cm_all[[i]] <- data.frame(
    Alpha = alpha,
    Threshold = th,
    TrueStatus = c("Induced","Induced","Not Induced","Not Induced"),
    Predicted  = c("Detected","Not Detected","Detected","Not Detected"),
    Count      = c(TP, FN, FP, TN)
  )

  rates_all[[i]] <- data.frame(
    Alpha = alpha,
    Threshold = th,
    TP = TP, FN = FN, FP = FP, TN = TN,
    TPR = TP / (TP + FN),
    FPR = FP / (FP + TN),
    TNR = TN / (TN + FP),
    FNR = FN / (FN + TP),
    OneMinusFPR = 1 - (FP / (FP + TN))
  )
}

det_long_df <- bind_rows(det_long_all)
cm_df       <- bind_rows(cm_all)
rates_df    <- bind_rows(rates_all)

## --------------------------
## write outputs
## --------------------------
write.csv(
  det_long_df,
  file.path(
    out_dir,
    paste0(n_persons, "n_", iterations, "_alpha_", alpha,
           "_interaction_detection_long05.csv")
  ),
  row.names = FALSE
)

write.csv(
  cm_df,
  file.path(
    out_dir,
    paste0(n_persons, "n_", iterations, "_alpha_", alpha,
           "_confusion_matrix05.csv")
  ),
  row.names = FALSE
)

write.csv(
  rates_df,
  file.path(
    out_dir,
    paste0(n_persons, "n_", iterations, "_alpha_", alpha,
           "_confusion_rates05.csv")
  ),
  row.names = FALSE
)


####Mega plotting

library(dplyr)
library(ggplot2)


n_val      <- 973
iterations <- 100

alphas    <- c(0.5, 1)
thresholds <- c(1e-4, 1e-3, 0.01, 0.05, 0.1, 0.5)


read_detection <- function(n, alpha, iter){
  fn <- file.path(
    "data/v3",
    paste0("n_", n, "_iter_100"),
    paste0(n, "n_100_alpha_", alpha, "_interaction_detection_long.csv")
    )
  read.csv(fn)
}

dl <- bind_rows(
  lapply(alphas, function(a) read_detection(n_val, a, iterations))
)

dl_induced <- dl %>% filter(TrueBeta != 0)

#Abs val only
dl_induced <- dl_induced %>%
  mutate(TrueBetaAbs = abs(TrueBeta))

power_df <- dl_induced %>% 
  group_by(n = n_val, Alpha, Threshold, TrueBetaAbs) %>%
  summarise(
    Power = mean(Detected),
    N = n(),
    .groups = "drop"
  )

power_df$Alpha <- factor(power_df$Alpha, levels = c(0.5, 1))
power_df$Threshold <- factor(power_df$Threshold)

#plot

p<- ggplot(
  power_df,
  aes(
    x = TrueBetaAbs,
    y = Power, 
    color = Threshold, 
    linetype = Alpha
  )
) + 
  geom_line(linewidth = 0.6) + 
  geom_point(size = 1) +
  scale_linetype_manual(
    values = c("0.5" = "dotted", "1" = "solid")
  )+
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  coord_cartesian(ylim = c(0, 1))+

  labs(
    x = "Induced interactions",
    y = "1-FNR (power)",
    color = "Threshold",
    linetype = "Alpha",
    title = paste("n = ", n_val)
  ) +
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x  = element_text(size = 15),
    axis.text.y  = element_text(size = 15),
    plot.title   = element_text(size = 26)
    #legend.position = "none"
  )

p

outfile <- file.path("plots",
                      paste0("n_", n_val, "_iter_100_power", ".png"))
 
ggsave(outfile, plot = p, width = 8, height = 5.5, dpi = 350)




###Optimal lambda histograms 
library(scales)

n_val      <- 973
iterations <- 100

alphas    <- c(0.5, 1)

read_lambda <- function(n, alpha, iter){
  fn <- file.path(
    "data/v3",
    paste0("n_", n, "_iter_100"),
    paste0(n, "n_100_alpha_", alpha, "_optimal_lambdas.csv")
  )
  
  df <- read.csv(fn)
  df$Alpha <- factor(alpha)
  df$n     <- n
  df
}

lambda_df <- bind_rows(
  read_lambda(n_val, 0.5),
  read_lambda(n_val, 1)
)

lambda_df <- lambda_df %>%
  mutate(
    logLambda = log(BestLambda),
    LambdaLabel = round(BestLambda, 4)
  )

lambda_df <- lambda_df %>%
  filter(BestLambda > 0)


lambda_grid_global <- exp(seq(log(1e-4), log(5e-1), length.out = 10))
log_breaks_global  <- log(lambda_grid_global)
binw               <- diff(log_breaks_global)[1]


# lambda_breaks <- sort(unique(lambda_df$logLambda))
# lambda_labels <- round(10^lambda_breaks, 4)
lambda_breaks <- log_breaks_global
lambda_labels <- round(lambda_grid_global, 4)


p <- ggplot(lambda_df, 
            aes(x = logLambda, fill = Alpha)) + 
  geom_histogram(binwidth = binw, 
                 boundary = min(log_breaks_global) - binw/2, 
                 position = position_dodge2(width = 2, padding = -0.5), 
                 alpha = 0.5, 
                 color = "white" ) + 
  scale_fill_manual(values = c("0.5" = "#4F7CD1", "1" = "#8671D1"), 
                    name = "Alpha" ) + 
  scale_x_continuous(name = "Optimal lambda", 
                     breaks = lambda_breaks, 
                     labels = lambda_labels, 
                     limits = c(min(log_breaks_global) - binw/2, 
                                max(log_breaks_global) + binw/2) ) + 
  labs( title = paste("n =", n_val), 
        y = "Count" ) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        legend.position = "inside", 
        legend.position.inside = c(0.02, 0.98), 
        legend.justification = c("left", "top"), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), 
        axis.text.y = element_text(size = 15), 
        plot.title = element_text(size = 26) 
        ) 

print(p + coord_cartesian(ylim = c(0, 80), expand = FALSE))

outfile <- file.path("plots",
                     paste0("n_", n_val, "_lambda_hist.png"))

ggsave(outfile, plot = p, width = 8, height = 5.5, dpi = 300)



##TPR vs FPR
#ROC

n_val      <- 973
iterations <- 100
out_dir <- file.path("data/v3", paste0("n_", n_val, "_iter_", iterations))

alphas <- c(0.5, 1)

read_rates <- function(a) {
  fn <- file.path(out_dir, paste0(n_val, "n_", iterations, "_alpha_", a, "_confusion_rates.csv"))
  df <- read.csv(fn)
  df$Alpha <- factor(a, levels = c(0.5, 1))
  df
}

rates_df <- bind_rows(lapply(alphas, read_rates)) %>%
  mutate(
    Threshold = as.numeric(as.character(Threshold)),
    ThresholdLab = factor(Threshold, levels = sort(unique(Threshold)))
  )

# extra_pts <- data.frame(
#   Alpha = factor(c(0.5, 1), levels = c(0.5, 1)),
#   ThresholdLab = factor(0.1, levels = levels(rates_df$ThresholdLab)),
#   FPR = c(0.205, 0.192),
#   TPR = c(0.953, 0.95)
# )

p <- ggplot(rates_df, aes(x = FPR, y = TPR, color = ThresholdLab, shape = Alpha)) +
  geom_point(
    size = 3,
    #position = position_jitter(width = 0.01, height = 0.01),
    alpha = 0.9
  ) +

#   geom_point(
#     data = extra_pts,
#     aes(x = FPR, y = TPR, color = ThresholdLab, shape = Alpha),
#     #position = position_jitter(width = 0.01, height = 0.01),
#     color = "black",
#     size = 3
#   ) +
  scale_shape_manual(
    values = c("0.5" = 16, "1" = 4),
    name = "Alpha"
  ) +
  scale_color_brewer(palette = "Dark2") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + 
  labs(
    x = "FPR",
    y = "TPR",
    color = "Threshold",
    title = paste0("n=", n_val)
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 15),
        axis.text.y  = element_text(size = 15),
        plot.title   = element_text(size = 26)
        #legend.position = "none"
        )

print(p)

outfile <- file.path("plots",
                     paste0("n_", n_val, "_fprtpr.png"))

ggsave(outfile, plot = p, width = 8, height = 5.5, dpi = 300)



####TRUE ROC
library(dplyr)
library(tidyr)
library(ggplot2)
library(pROC)

n_val      <- 973
iterations <- 100
out_dir <- file.path("data/v3", paste0("n_", n_val, "_iter_", iterations))

alphas <- c(0.5, 1)

# You need two inputs per alpha:
# 1) Bhat: iterations x p_interactions matrix (estimated interaction coefficients)
# 2) y_true: length p_interactions vector in {0,1} indicating which interactions are truly active

read_bhat <- function(a) {
  # Example filenames — replace with your real stored outputs
  # Must load the *interaction coefficient* estimates per iteration.
  fn <- file.path(out_dir, paste0(n_val, "n_", iterations, "_alpha_", a, "_bhat_interactions.rds"))
  readRDS(fn)  # should be a numeric matrix: iterations x p_interactions
}

read_truth <- function() {
  # Truth vector used in the semi-synthetic DGP.
  fn <- file.path(out_dir, paste0(n_val, "n_", iterations, "_truth_active_interactions.rds"))
  readRDS(fn)  # numeric/integer vector length p_interactions with 0/1
}

y_true <- read_truth()

roc_df_for_alpha <- function(a) {
  Bhat <- read_bhat(a)
  
  stopifnot(is.matrix(Bhat))
  stopifnot(ncol(Bhat) == length(y_true))
  
  # pooled scores and labels across iterations
  scores <- as.vector(abs(Bhat))                         # length iterations*p
  labels <- rep(y_true, times = nrow(Bhat))              # align with as.vector column-major? see note below
  
  # IMPORTANT alignment note:
  # as.vector() in R is column-major, so it stacks columns.
  # rep(y_true, times=nrow(Bhat)) matches that (each column label repeated for all rows).
  # If you used c(t(abs(Bhat))) instead, you'd need rep(y_true, each=??) differently.
  
  roc_obj <- pROC::roc(response = labels, predictor = scores, quiet = TRUE, direction = ">")
  
  tibble(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    Alpha = factor(a, levels = c(0.5, 1)),
    AUC = as.numeric(pROC::auc(roc_obj))
  )
}

roc_df <- bind_rows(lapply(alphas, roc_df_for_alpha))

# AUC labels (one per alpha)
auc_labels <- roc_df %>%
  group_by(Alpha) %>%
  summarise(AUC = first(AUC), .groups = "drop") %>%
  mutate(label = paste0("AUC=", sprintf("%.3f", AUC)))

ggplot(roc_df, aes(x = FPR, y = TPR, linetype = Alpha)) +
  geom_step(linewidth = 1.0) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    x = "FPR",
    y = "TPR",
    linetype = "Alpha",
    title = paste0("n=", n_val, " (true ROC from |beta_hat|)")
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 15),
        axis.text.y  = element_text(size = 15),
        plot.title   = element_text(size = 26))





###Plots vanilla MNL

##FPR vs TPR
n_val      <- 50
iterations <- 100

out_dir <- file.path("data/v3_vanilla", paste0("n_", n_val, "_iter_", iterations))

rates_df <- read.csv(
  file.path(out_dir, paste0(n_val, "n_100_alpha_1_confusion_rates.csv"))
) %>%
  mutate(
    Threshold = as.numeric(as.character(Threshold)),
    ThresholdLab = factor(Threshold, levels = sort(unique(Threshold)))
  )

p <- ggplot(
  rates_df,
  aes(x = FPR, y = TPR, color = ThresholdLab)
) +
  geom_point(
    size = 4,
    alpha = 0.9
  ) +
  scale_color_brewer(palette = "Dark2") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + 
  labs(
    x = "FPR",
    y = "TPR",
    color = "Threshold",
    title = paste0("Unpenalized (n=", n_val, ")")
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

outfile <- file.path("plots",
                     paste0("n_", n_val, "_fprtpr_vanilla.png"))

ggsave(outfile, plot = p, width = 8, height = 5.5, dpi = 300)


###Power vanilla

n_val      <- 973
iterations <- 100

thresholds <- c(1e-4, 1e-3, 0.01, 0.05, 0.1, 0.5)

read_detection <- function(n, iter){
  fn <- file.path(
    "data/v3_vanilla",
    paste0("n_", n, "_iter_100"),
    paste0(n, "n_100_alpha_1_interaction_detection_long.csv")
  )
  read.csv(fn)
}

dl <- read_detection(n_val, iterations)

dl_induced <- dl %>%
  filter(TrueBeta != 0) %>%
  mutate(TrueBetaAbs = abs(TrueBeta))

power_df <- dl_induced %>% 
  group_by(n = n_val, Threshold, TrueBetaAbs) %>%
  summarise(
    Power = mean(Detected),
    N = n(),
    .groups = "drop"
  )

power_df$Threshold <- factor(power_df$Threshold)

p <- ggplot(
  power_df,
  aes(
    x = TrueBetaAbs,
    y = Power,
    color = Threshold
  )
) + 
  geom_line(linewidth = 0.6) + 
  geom_point(size = 1) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Induced interactions",
    y = "1-FNR (power)",
    color = "Threshold",
    title = paste("Unpenalized (n =", n_val, ")")
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x  = element_text(size = 15),
    axis.text.y  = element_text(size = 15),
    plot.title   = element_text(size = 26)
  )

p

outfile <- file.path(
  "plots",
  paste0("n_", n_val, "_iter_100_power_vanilla.png")
)

ggsave(outfile, plot = p, width = 8, height = 5.5, dpi = 350)


