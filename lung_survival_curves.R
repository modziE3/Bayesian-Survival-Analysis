# file-name = "lung_survival_curves.R"

library(survival)

samples <- fit$sims.matrix
k_samples <- samples[, "k"]
lambda_samples <- samples[, "lambda"]

time_grid <- seq(
  0.01,
  max(time_months),
  length.out = 300
)

survival_samples <- sapply(seq_along(k_samples), function(j) {
  exp(-lambda_samples[j] * time_grid^k_samples[j])
})

survival_mean <- rowMeans(survival_samples)

survival_lower <- apply(
  survival_samples, 1, quantile, probs = 0.025
)

survival_upper <- apply(
  survival_samples, 1, quantile, probs = 0.975
)

density_samples <- sapply(seq_along(k_samples), function(j) {
  k_samples[j] * lambda_samples[j] *
    time_grid^(k_samples[j] - 1) *
    exp(-lambda_samples[j] * time_grid^k_samples[j])
})

density_mean <- rowMeans(density_samples)

density_lower <- apply(
  density_samples, 1, quantile, probs = 0.025
)

density_upper <- apply(
  density_samples, 1, quantile, probs = 0.975
)

# kaplan meier estimate
km <- survfit(Surv(time_months, event) ~ 1)

mean_col <- "blue"
km_col <- "red"
ci_col <- adjustcolor("skyblue", alpha.f = 0.4)
grid_col <- "grey50"

pdf("weibull_survival.pdf", width = 8, height = 6)

par(mar = c(5, 4, 4, 2) + 0.1)

plot(
  time_grid, survival_mean,
  type = "n",
  ylim = c(0, 1),
  xlab = "Time (months)",
  ylab = "Survival probability",
  main = "Bayesian Weibull Survival Model"
)

# Grid behind data
grid(col = grid_col, lty = "dotted")

polygon(
  c(time_grid, rev(time_grid)),
  c(survival_lower, rev(survival_upper)),
  col = ci_col,
  border = NA
)

lines(
  time_grid, survival_mean,
  lwd = 2,
  col = mean_col
)

lines(
  km,
  conf.int = FALSE,
  col = km_col,
  lwd = 2,
  lty = 1
)

legend(
  "topright",
  legend = c(
    "Bayesian Weibull mean",
    "95% credible interval",
    "Kaplan-Meier estimate"
  ),
  col = c(mean_col, ci_col, km_col),
  lwd = c(2, NA, 2),
  lty = c(1, NA, 1),
  pch = c(NA, 15, NA),
  pt.cex = 1.5,
  pt.bg = ci_col,
  bty = "n",
  cex = 0.85
)

density_ylim <- c(
  0,
  max(density_upper, na.rm = TRUE) * 1.10
)

plot(
  time_grid, density_mean,
  type = "n",
  ylim = density_ylim,
  xlab = "Time (months)",
  ylab = "Probability density",
  main = "Weibull Survival-Time Density"
)

grid(col = grid_col, lty = "dotted")

polygon(
  c(time_grid, rev(time_grid)),
  c(density_lower, rev(density_upper)),
  col = ci_col,
  border = NA
)

lines(
  time_grid, density_mean,
  lwd = 2,
  col = mean_col
)

legend(
  "topright",
  legend = c(
    "Posterior mean density",
    "95% credible interval"
  ),
  col = c(mean_col, ci_col),
  lwd = c(2, NA),
  lty = c(1, NA),
  pch = c(NA, 15),
  pt.cex = 1.5,
  pt.bg = ci_col,
  bty = "n",
  cex = 0.85
)

dev.off()

cat(
  "PDF saved to:",
  file.path(getwd(), "weibull_survival.pdf"),
  "\n"
)