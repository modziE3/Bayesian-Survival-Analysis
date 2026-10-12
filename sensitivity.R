# file-name = "sensitivity.R"

library(survival)
library(R2WinBUGS)

bugs_dir <- "C:/WinBUGS14"
stopifnot(file.exists(file.path(bugs_dir, "WinBUGS14.exe")))

data(lung)

N <- nrow(lung)
time_months <- lung$time / 30
event <- ifelse(lung$status == 2, 1, 0)

stopifnot(all(time_months > 0))
stopifnot(sum(event == 1) == 165)
stopifnot(sum(event == 0) == 63)

bugs_data <- list(
  N = N,
  t = time_months,
  event = event,
  zeros = rep(0, N)
)

prior_settings <- list(
  
  "Research informed" = list(
    k_prior = "k ~ dgamma(2.6, 2)",
    lambda_prior = "lambda ~ dgamma(2, 80)",
    log_k = FALSE,
    init_k = 1.3,
    init_lambda = 0.03
  ),
  
  "Sensible initial" = list(
    k_prior = "k ~ dgamma(1, 1)",
    lambda_prior = "lambda ~ dgamma(2, 20)",
    log_k = FALSE,
    init_k = 1,
    init_lambda = 0.1
  ),
  
  "Weakly informed" = list(
    k_prior = "log_k ~ dunif(-5, 5)\n    k <- exp(log_k)",
    lambda_prior = "lambda ~ dgamma(1, 10)",
    log_k = TRUE,
    init_k = 1,
    init_lambda = 0.1
  )
)

fit_prior_model <- function(prior, model_name) {
  
  model_code <- paste0("
  model {
    for (i in 1:N) {
      loglik[i] <- event[i] * (
        log(lambda) + log(k) +
        (k - 1) * log(t[i])
      ) - lambda * pow(t[i], k)

      # Zero-trick likelihood
      phi[i] <- 10000 - loglik[i]
      zeros[i] ~ dpois(phi[i])
    }

    ", prior$k_prior, "

    ", prior$lambda_prior, "

    # Posterior quantities
    S12 <- exp(-lambda * pow(12, k))
    median <- pow(log(2) / lambda, 1 / k)
  }
  ")
  
  model_file <- paste0(
    "weibull_",
    gsub(" ", "_", tolower(model_name)),
    ".txt"
  )
  
  writeLines(model_code, model_file)
  
  inits <- function() {
    if (prior$log_k) {
      list(
        log_k = 0,
        lambda = prior$init_lambda
      )
    } else {
      list(
        k = prior$init_k,
        lambda = prior$init_lambda
      )
    }
  }
  
  bugs(
    data = bugs_data,
    inits = inits,
    parameters.to.save = c("k", "lambda", "S12", "median"),
    model.file = model_file,
    n.chains = 2,
    n.iter = 10000,
    n.burnin = 2000,
    n.thin = 2,
    bugs.directory = bugs_dir,
    debug = FALSE
  )
}

fits <- lapply(names(prior_settings), function(name) {
  cat("\nFitting:", name, "\n")
  fit_prior_model(prior_settings[[name]], name)
})

names(fits) <- names(prior_settings)
parameters <- c("k", "lambda", "S12", "median")
summary_columns <- c("mean", "sd", "2.5%", "50%", "97.5%", "Rhat", "n.eff")

results <- do.call(rbind, lapply(names(fits), function(name) {
  s <- fits[[name]]$summary[
    parameters, summary_columns, drop = FALSE
  ]
  
  data.frame(
    Prior = name,
    Parameter = rownames(s),
    s,
    row.names = NULL,
    check.names = FALSE
  )
}))

cat("\n--- Posterior summaries and diagnostics ---\n")
print(results, digits = 4, row.names = FALSE)
comparison <- do.call(rbind, lapply(parameters, function(p) {
  
  do.call(rbind, lapply(names(fits), function(name) {
    
    s <- fits[[name]]$summary[p, ]
    
    data.frame(
      Parameter = p,
      Prior = name,
      Mean = unname(s["mean"]),
      Lower_95 = unname(s["2.5%"]),
      Upper_95 = unname(s["97.5%"])
    )
  }))
}))

cat("\n--- Posterior comparison ---\n")
print(comparison, digits = 4, row.names = FALSE)
baseline <- fits[["Sensible initial"]]$summary

sensitivity <- do.call(rbind, lapply(names(fits), function(name) {
  s <- fits[[name]]$summary[parameters, ]
  data.frame(
    Prior = name,
    Parameter = parameters,
    Baseline_mean = baseline[parameters, "mean"],
    Posterior_mean = s[, "mean"],
    Absolute_change = abs(
      s[, "mean"] - baseline[parameters, "mean"]
    ),
    row.names = NULL
  )
}))

sensitivity$Relative_change_percent <- with(
  sensitivity,
  ifelse(
    abs(Baseline_mean) > 1e-10,
    100 * (Posterior_mean - Baseline_mean) /
      abs(Baseline_mean),
    NA_real_
  )
)

cat("\n--- Prior sensitivity analysis ---\n")
print(sensitivity, digits = 4, row.names = FALSE)

pdf_file <- file.path(
  getwd(),
  "weibull_prior_sensitivity.pdf"
)

# colour blind friencly
prior_colours <- c(
  "Research informed" = "#0072B2",
  "Sensible initial"  = "#D55E00",
  "Weakly informed"   = "#009E73"
)

prior_lty <- c(
  "Research informed" = 1,
  "Sensible initial"  = 2,
  "Weakly informed"   = 3
)

pdf(
  file = pdf_file,
  width = 11.7,
  height = 8.3,
  pointsize = 11,
  family = "Helvetica"
)

par(
  mfrow = c(2, 2),
  mar = c(4.2, 4.2, 3.2, 1),
  oma = c(0, 0, 2, 0),
  mgp = c(2.4, 0.7, 0),
  las = 1
)

for (p in parameters) {
  densities <- lapply(fits, function(fit) {
    density(fit$sims.matrix[, p], n = 1024)
  })
  
  names(densities) <- names(fits)
  xlim <- range(unlist(lapply(densities, function(d) d$x)))
  ylim <- c(
    0,
    max(unlist(lapply(densities, function(d) d$y))) * 1.12
  )
  plot(
    NA,
    xlim = xlim,
    ylim = ylim,
    type = "n",
    main = paste("Posterior distribution:", p),
    xlab = p,
    ylab = "Density",
    cex.main = 1,
    cex.lab = 1,
    cex.axis = 0.85
  )
  grid(
    col = "grey85",
    lty = "dotted",
    lwd = 0.8
  )
  for (name in names(fits)) {
    lines(
      densities[[name]]$x,
      densities[[name]]$y,
      col = prior_colours[name],
      lty = prior_lty[name],
      lwd = 2.2
    )
  }
  legend(
    "topright",
    legend = names(fits),
    col = prior_colours[names(fits)],
    lty = prior_lty[names(fits)],
    lwd = 2.2,
    cex = 0.78,
    bg = "white",
    box.col = "grey75",
    inset = 0.02,
    x.intersp = 0.7,
    y.intersp = 1.1
  )
}

mtext(
  "Prior sensitivity analysis: Weibull survival model",
  outer = TRUE,
  cex = 1.2,
  font = 2
)

dev.off()

cat("\nPDF saved successfully to:\n", pdf_file, "\n")

