data {
  int<lower=1> T;
  int<lower=1> K;
  matrix[T, K] X;
  vector[T] y;
  int<lower=1, upper=T> t0;
  int<lower=1> S;
}

parameters {
  real mu_0;                    
  vector[T-1] mu_raw;           
  vector[S-1] season_raw;
  vector[K] beta;
  real<lower=0> sigma_y;
  real<lower=0> sigma_level;
  real<lower=0> sigma_season;
  real impact;
}

transformed parameters {
  vector[T] mu;
  vector[S] season;
  vector[T] y_hat;
  
  mu[1] = mu_0;
  for (t in 2:T) {
    mu[t] = mu[t-1] + sigma_level * mu_raw[t-1];
  }
  
  season[1:(S-1)] = season_raw;
  season[S] = -sum(season_raw);
  
  for (t in 1:T) {
    int s_idx = 1 + ((t - 1) % S);
    y_hat[t] = mu[t] + season[s_idx] + X[t] * beta;
    if (t >= t0) y_hat[t] += impact;
  }
}

model {
  sigma_y ~ inv_gamma(0.01, 0.01 * variance(y));
  sigma_level ~ inv_gamma(50, 50 * 0.01 * sd(y));
  sigma_season ~ inv_gamma(0.01, 0.01 * sd(y));
  
  mu_0 ~ normal(y[1], 10);
  mu_raw ~ std_normal();  
  season_raw ~ normal(0, sigma_season);
  beta ~ normal(0, 1);
  impact ~ normal(0, 5);
  
  y ~ normal(y_hat, sigma_y);
}