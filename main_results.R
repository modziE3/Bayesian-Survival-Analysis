# file-name = "main_results.R"

samples <- fit$sims.matrix

k_samples <- samples[, "k"]
lambda_samples <- samples[, "lambda"]
median_samples <- samples[, "median"]

# (a) Posterior distribution of survival probability at 5 months
S5_samples <- exp(-lambda_samples * 5^k_samples)

cat("\n--- Survival probability at 5 months ---\n")
print(quantile(S5_samples, c(0.025, 0.25, 0.5, 0.75, 0.975)))

cat("Posterior mean S(5):", mean(S5_samples), "\n")
cat("Posterior SD S(5):", sd(S5_samples), "\n")

# (b) Probability of surviving beyond 5 months
# Posterior predictive probability for a randomly selected
# individual under the fitted population survival model.
cat("Estimated probability of surviving beyond 5 months:",
    mean(S5_samples), "\n")

# (c) Posterior distribution of median survival time
cat("\n--- Median survival time (months) ---\n")
print(quantile(median_samples,
               c(0.025, 0.25, 0.5, 0.75, 0.975)))

cat("Posterior mean median survival:",
    mean(median_samples), "months\n")

# (d) Is the hazard increasing or decreasing?
cat("\n--- Weibull hazard behaviour ---\n")
cat("Posterior probability hazard increases (k > 1):",
    mean(k_samples > 1), "\n")
cat("Posterior probability hazard decreases (k < 1):",
    mean(k_samples < 1), "\n")
cat("Posterior probability hazard is constant (k = 1):",
    mean(k_samples == 1), "\n")

# (e) Posterior uncertainty in Weibull parameters
cat("\n--- Parameter summaries ---\n")
print(fit$summary[
  c("k", "lambda"),
  c("mean", "sd", "2.5%", "50%", "97.5%", "Rhat", "n.eff")
])

# Save numerical results for the report
results <- data.frame(
  Quantity = c("S(5)", "Median survival time (months)",
               "Weibull shape k", "Weibull rate lambda"),
  Mean = c(mean(S5_samples), mean(median_samples),
           mean(k_samples), mean(lambda_samples)),
  SD = c(sd(S5_samples), sd(median_samples),
         sd(k_samples), sd(lambda_samples)),
  Lower95 = c(
    quantile(S5_samples, 0.025),
    quantile(median_samples, 0.025),
    quantile(k_samples, 0.025),
    quantile(lambda_samples, 0.025)
  ),
  Upper95 = c(
    quantile(S5_samples, 0.975),
    quantile(median_samples, 0.975),
    quantile(k_samples, 0.975),
    quantile(lambda_samples, 0.975)
  )
)

print(results, digits = 4)
write.csv(results, "bayesian_results.csv", row.names = FALSE)
