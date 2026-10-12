# file-name = "additional_plots.R"

hazard_samples <- sapply(seq_along(k_samples), function(j) {
  lambda_samples[j] * k_samples[j] *
    time_grid^(k_samples[j] - 1)
})

hazard_mean <- rowMeans(hazard_samples)

hazard_lower <- apply(
  hazard_samples, 1, quantile, probs = 0.025
)

hazard_upper <- apply(
  hazard_samples, 1, quantile, probs = 0.975
)

mean_col <- "blue"
ci_col <- adjustcolor("skyblue", alpha.f = 0.4)
grid_col <- "grey50"

pdf("weibull_additional_results.pdf", width = 8, height = 6)

par(mar = c(5, 4, 4, 2) + 0.1)

plot(
  time_grid, hazard_mean,
  type = "n",
  ylim = c(0, max(hazard_upper, na.rm = TRUE) * 1.10),
  xlab = "Time (months)",
  ylab = "Hazard rate (per month)",
  main = "Posterior Weibull Hazard Function"
)

# Grid behind data
grid(col = grid_col, lty = "dotted")

# 95% credible interval
polygon(
  c(time_grid, rev(time_grid)),
  c(hazard_lower, rev(hazard_upper)),
  col = ci_col,
  border = NA
)

lines(
  time_grid, hazard_mean,
  lwd = 2,
  col = mean_col
)

legend(
  "topright",
  legend = c(
    "Posterior mean hazard",
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

hist(
  median_samples,
  breaks = 35,
  probability = TRUE,
  xlab = "Median survival time (months)",
  ylab = "Posterior density",
  main = "Posterior Distribution of Median Survival",
  col = "lightblue",
  border = "white"
)

grid(col = grid_col, lty = "dotted")

# Add a smooth density estimate
density_fit <- density(median_samples)

lines(
  density_fit,
  col = "blue",
  lwd = 2
)

posterior_median <- median(median_samples)

abline(
  v = posterior_median,
  lwd = 2,
  col = "red",
  lty = 2
)

legend(
  "topright",
  legend = c(
    "Posterior density",
    paste0(
      "Posterior median = ",
      round(posterior_median, 2),
      " months"
    )
  ),
  col = c("blue", "red"),
  lwd = 2,
  lty = c(1, 2),
  bty = "n",
  cex = 0.85
)

dev.off()

cat(
  "Additional plots saved to:",
  file.path(getwd(), "weibull_additional_results.pdf"),
  "\n"
)