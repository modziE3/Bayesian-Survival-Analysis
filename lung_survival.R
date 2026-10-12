# file-name = "lung_survival.R"


library(survival)
library(R2WinBUGS)

bugs_dir <- "C:/WinBUGS14"
file.exists(file.path(bugs_dir, "WinBUGS14.exe"))
bugs.directory = "C:/WinBUGS14"



data(lung)
N <- nrow(lung)
time_months <- lung$time / 30

# 1 = death, 0 = censored
event <- ifelse(lung$status == 2, 1, 0)

stopifnot(all(time_months > 0))
stopifnot(sum(event == 1) == 165)
stopifnot(sum(event == 0) == 63)

bugs_data <- list(
  N = N,
  t = time_months,
  event = event
)

model_code <- "
model {
    for (i in 1:N) {
        loglik[i] <- event[i] * (
            log(lambda) + log(k) +
            (k - 1) * log(t[i])
        ) - lambda * pow(t[i], k)

        phi[i] <- 1000 - loglik[i]
        zeros[i] ~ dpois(phi[i])
    }

    # Priors (initial choices; assess sensitivity later)
    k ~ dgamma(1, 1)
    lambda ~ dgamma(2, 20)

    # Posterior quantities
    S12 <- exp(-lambda * pow(12, k))
    median <- pow(log(2) / lambda, 1 / k)
}
"

writeLines(model_code, "weibull_model.txt")

bugs_data$zeros <- rep(0, N)

inits <- function() {
  list(
    k = 1,
    lambda = 0.1
  )
}


fit <- bugs(
  data = bugs_data,
  inits = inits,
  parameters.to.save = c("k", "lambda", "S12", "median"),
  model.file = "weibull_model.txt",
  n.chains = 2,
  n.iter = 10000,
  n.burnin = 2000,
  n.thin = 2,
  bugs.directory = "C:/WinBUGS14",
  debug = TRUE
)

print(fit, digits = 3)

head(fit$sims.matrix)

fit$summary[c("k", "lambda", "S12", "median"), ]

fit$summary[c("k", "lambda", "S12", "median"),
            c("mean", "sd", "2.5%", "50%", "97.5%", "Rhat", "n.eff")]

table(lung$sex)