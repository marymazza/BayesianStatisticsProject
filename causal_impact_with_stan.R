#stan Model recreating Causal Impact

#Libraries----
library(rstan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(bayesplot)
library(patchwork)

#Upload Data----
causalImpact_data_interpolated <- read.csv("causalImpact_data_interpolated.csv")
causalImpact_data_interpolated <- causalImpact_data_interpolated[,-1]

causalImpact_data_interpolated$date <- as.Date(causalImpact_data_interpolated$date)

causalImpact_data_interpolated <- causalImpact_data_interpolated[order(causalImpact_data_interpolated$date),]

# Outcome
y <- causalImpact_data_interpolated$avg_val_inside

# Covariate
X <- as.matrix(causalImpact_data_interpolated[ ,c("T_max","T_min","rain","wind",
                                                 "avg_val_outside") ])

T_ <- length(y)
K   <- ncol(X)
S   <- 52   
t0 <- which(causalImpact_data_interpolated$date == as.Date("2012-01-16"))

#Create Stan data----
stan_data <- list(
  T = T_,
  K = K,
  X = X,
  y = y,
  t0 = t0,
  S = S
)

#Stan Model----
fit <- stan(
  file = "causal_impact_stan_model.stan",
  data = stan_data,
  chains = 4,
  iter = 4000,
  warmup = 1000,
  cores = 4,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)

#Posterior----
post <- as.matrix(fit)
mu_post     <- post[, grep("^mu\\[", colnames(post)), drop = FALSE]
season_post <- post[, grep("^season\\[", colnames(post)), drop = FALSE]
beta_post   <- post[, grep("^beta\\[", colnames(post)), drop = FALSE]
impact_post      <- post[, "impact"]
sigma_level_post <- post[, "sigma_level"]
sigma_y_post     <- post[, "sigma_y"]
T <- length(y)
K <- ncol(X)
S <- ncol(season_post)

X_beta <- matrix(0, nrow = nrow(post), ncol = T)
for(k in 1:K){
  X_beta <- X_beta + beta_post[,k] %*% t(X[,k])
}

season_full <- matrix(0, nrow = nrow(post), ncol = T)
for(t in 1:T){
  idx <- 1 + ((t-1) %% S)
  season_full[,t] <- season_post[, idx]
}

mu_cf <- mu_post   

for (i in 1:nrow(post)) {
  for (t in t0:T) {
    if (t == t0) {
      mu_cf[i, t] <- mu_cf[i, t-1] + rnorm(1, 0, sigma_level_post[i])
    } else {
      mu_cf[i, t] <- mu_cf[i, t-1] + rnorm(1, 0, sigma_level_post[i])
    }
  }
}

y_cf_post <- mu_cf + season_full + X_beta
y_fitted_post <- y_cf_post
y_fitted_post[, t0:T] <- y_fitted_post[, t0:T] + impact_post
summ_cf <- apply(y_cf_post, 2, function(z) {
  c(mean = mean(z),
    lo   = quantile(z, 0.05),
    hi   = quantile(z, 0.95))
})

summ_cf <- as.matrix(summ_cf)  
summ_fit <- apply(y_fitted_post, 2, function(z) {
  c(mean = mean(z),
    lo   = quantile(z, 0.05),
    hi   = quantile(z, 0.95))
})

summ_fit <- as.matrix(summ_fit)

df <- data.frame(
  time = 1:T,
  y = y,
  cf_mean = summ_cf["mean", ],
  cf_lo   = summ_cf["lo.5%", ],
  cf_hi   = summ_cf["hi.95%", ],
  diff_mean = y - summ_cf["mean", ],
  diff_lo   = y - summ_cf["hi.95%", ],
  diff_hi   = y - summ_cf["lo.5%", ],
  cum_mean = cumsum(y - summ_cf["mean", ]),
  cum_lo   = cumsum(y - summ_cf["hi.95%", ]),
  cum_hi   = cumsum(y - summ_cf["lo.5%", ])
)

#Plotting the results----
p1 <- ggplot(df, aes(time)) +
  geom_ribbon(aes(ymin = cf_lo, ymax = cf_hi), fill="skyblue", alpha=.4) +
  geom_line(aes(y = cf_mean), color="blue", size=0.6) +
  geom_line(aes(y = y), color="black") +
  geom_vline(xintercept = t0, linetype=2) +
  labs(title="Original (Observed vs Counterfactual)") +
  theme_bw()

p2 <- ggplot(df, aes(time)) +
  geom_ribbon(aes(ymin = diff_lo, ymax = diff_hi), fill="skyblue", alpha=.4) +
  geom_line(aes(y = diff_mean), color="blue", size=0.6) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = t0, linetype=2) +
  labs(title="Pointwise Effect") +
  theme_bw()

p3 <- ggplot(df, aes(time)) +
  geom_ribbon(aes(ymin = cum_lo, ymax = cum_hi), fill="skyblue", alpha=.4) +
  geom_line(aes(y = cum_mean), color="blue", size=0.6) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = t0, linetype=2) +
  labs(title="Cumulative Effect") +
  theme_bw()

p1 / p2 / p3

#Summary----
make_impact_summary <- function(y, y_cf_post, t0) {
  S <- nrow(y_cf_post)
  T <- length(y)
  post_idx <- t0:T
  n_post <- length(post_idx)
  cf_mean_post <- apply(y_cf_post[, post_idx, drop = FALSE], 2, mean)
  cf_lo_post   <- apply(y_cf_post[, post_idx, drop = FALSE], 2, quantile, 0.05)
  cf_hi_post   <- apply(y_cf_post[, post_idx, drop = FALSE], 2, quantile, 0.95)
  eff_post <- sweep(matrix(y[post_idx], S, n_post, byrow = TRUE),
                    2, y_cf_post[, post_idx], "-")
  pt_mean <- apply(eff_post, 2, mean)
  pt_lo   <- apply(eff_post, 2, quantile, 0.05)
  pt_hi   <- apply(eff_post, 2, quantile, 0.95)
  cum_post <- apply(eff_post, 1, cumsum)
  cum_mean <- apply(cum_post, 1, function(x) x[n_post]) |> mean()
  cum_lo   <- apply(cum_post, 1, function(x) x[n_post]) |> quantile(0.05)
  cum_hi   <- apply(cum_post, 1, function(x) x[n_post]) |> quantile(0.95)
  avg_mean <- mean(pt_mean)
  avg_lo   <- quantile(rowMeans(eff_post), 0.05)
  avg_hi   <- quantile(rowMeans(eff_post), 0.95)
  tail_prob <- mean(rowMeans(eff_post) < 0)
  if(avg_mean > 0) tail_prob <- 1 - tail_prob
  cf_sum <- apply(y_cf_post[, post_idx], 1, sum)
  rel_eff <- cum_post[n_post, ] / cf_sum
  rel_mean <- mean(rel_eff)
  rel_lo   <- quantile(rel_eff, 0.05)
  rel_hi   <- quantile(rel_eff, 0.95)
  
  summary_df <- data.frame(
    Metric = c("Average pointwise effect", 
               "Cumulative effect", 
               "Relative effect"),
    Mean = c(avg_mean, cum_mean, rel_mean),
    `Lower 95%` = c(avg_lo, cum_lo, rel_lo),
    `Upper 95%` = c(avg_hi, cum_hi, rel_hi)
  )
  list(
    summary = summary_df,
    tail_probability = tail_prob
  )
}

impact_summary <- make_impact_summary(y, y_cf_post, t0)
impact_summary$summary
impact_summary$tail_probability




