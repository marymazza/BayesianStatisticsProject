#Spatial Stan Model

#Libraries----
library(rstan)
library(MASS)
library(dplyr)
library(tidyr)
library(prettymapr)
library(curl)
library(httr2)
library(ggplot2)
library(giscoR)
library(sf)
library(ggspatial)

#Upload the Dataset----

prepare_multivariate_data <- function(data_imputed, intervention_date) {
  
    y_wide <- data_imputed %>%
        select(idSensore, date, avg_val) %>%
        pivot_wider(
        names_from = idSensore,
        values_from = avg_val,
        names_prefix = "station_"
      ) %>%
      arrange(date)
  
  station_cols <- grep("^station_", names(y_wide), value = TRUE)
  M <- length(station_cols)
  
  Y_raw <- as.matrix(y_wide[, station_cols])
  
  for (j in 1:ncol(Y_raw)) {
    Y_raw[is.na(Y_raw[, j]), j] <- mean(Y_raw[, j], na.rm = TRUE)
  }
  
  Y <- scale(t(Y_raw))
  T <- ncol(Y)
  dates <- y_wide$date
  
  station_info <- data_imputed %>%
    distinct(idSensore, .keep_all = TRUE) %>%
    select(idSensore, UTM_est, UTM_nord, lat, long) %>%
    mutate(station_col = paste0("station_", idSensore)) %>%
    slice(match(station_col, station_cols))
  
  coords <- as.matrix(station_info[, c("UTM_est", "UTM_nord")])
  D <- as.matrix(dist(coords, method = "euclidean"))
  
  covariate_names <- c("T_max", "T_min", "rain", "wind")
  K <- length(covariate_names)
  
  X <- array(NA, dim = c(T, M, K))
  
  for (i in 1:M) {
    
    st_id <- gsub("station_", "", station_cols[i])
    
    st_data <- data_imputed %>%
      filter(idSensore == st_id)
    
    for (k in 1:K) {
      
      x <- st_data %>%
        select(date, !!covariate_names[k]) %>%
        right_join(
          tibble(date = dates),
          by = "date"
        ) %>%
        arrange(date) %>%
        pull(!!covariate_names[k])
      
      x[is.na(x)] <- mean(x, na.rm = TRUE)
      X[, i, k] <- as.numeric(scale(x))
    }
  }
  
  t0 <- which(dates == intervention_date)
 
  list(
    Y = Y,
    X = X,
    D = D,
    T = T,
    M = M,
    P = K,
    t0 = t0,
    dates = dates,
    station_ids = gsub("^station_", "", station_cols),
    station_info = station_info
  )
}
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

data <- read.csv("data_region_interpolated.csv")
intervention_date <- as.Date("2012-01-16")
data_stan <- prepare_multivariate_data(data, intervention_date)

#Split data----

which(data$date == intervention_date)
intervention_time <- 107 
T_train <- intervention_time - 1

y_train <- data_stan$Y[, 1:T_train]
X_train <- data_stan$X[1:T_train, , ]

#Threshold on distances----
threshold <- 100000  # 100 km

adj <- matrix(0, nrow = 51, ncol = 51)
adj[data_stan$D < threshold] <- 1
diag(adj) <- 1
pairs <- which(adj == 1 & upper.tri(adj), arr.ind = TRUE)
N_pairs <- nrow(pairs)
pair_i <- pairs[, 1]
pair_j <- pairs[, 2]

#Stan data----
ell <- max(data_stan$D)/3

stan_data <- list(
  M = nrow(y_train),
  P = dim(X_train)[3],
  T = T_train,
  Y = y_train,
  X = lapply(1:T_train, function(t) X_train[t, , ]),
  D = data_stan$D,
  ell = ell,
  N_pairs = N_pairs,
  pair_i = pair_i,
  pair_j = pair_j
)

#Fit model----

init_fun <- function() {
  list(
    alpha = 0.3,
    beta = rep(0, data_stan$P),
    sigma = 0.5,
    sigma_w = 0.5,
    z_w = rep(0, data_stan$M),
    rho_raw = rnorm(N_pairs, 0, 0.05)
  )
}

fit <- stan(
  file = "spatial_stan_model.stan",
  data = stan_data,
  chains = 4,
  cores = 4,
  iter = 2000,
  warmup = 1500,
  init = init_fun,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  refresh = 50
)

print(fit, pars = c("alpha", "beta", "sigma", "sigma_w"))


#Posterior prediction----

post <- rstan::extract(fit)
alpha_draws <- post$alpha
beta_draws  <- post$beta
sigma_draws <- post$sigma
w_draws     <- post$w
L_draws     <- post$L_Omega

S <- length(alpha_draws)
M <- data_stan$M
T <- data_stan$
Y <- data_stan$Y
X <- data_stan$X
T0 <- 107

Y_cf <- array(NA, dim = c(S, M, T - T0))
for (s in 1:S) {
  alpha_s <- alpha_draws[s]
  beta_s  <- beta_draws[s, ]
  sigma_s <- sigma_draws[s]
  w_s     <- w_draws[s, ]
  L_s     <- L_draws[s,,]
  Y_prev <- Y[, T0]
  for (t in (T0 + 1):T) {
    mu <- alpha_s * Y_prev + X[t,,] %*% beta_s + w_s
    z  <- rnorm(M)
    eps <- sigma_s * (L_s %*% z)
    Y_new <- mu + eps
    Y_cf[s, , t - T0] <- Y_new
    Y_prev <- Y_new
  }
}
Y_obs_post <- Y[, (T0 + 1):T]

obs_mean <- colMeans(Y_obs_post)
cf_mean_draws <- apply(Y_cf, c(1,3), mean)
cf_mean <- apply(cf_mean_draws, 2, mean)
cf_lower <- apply(cf_mean_draws, 2, quantile, 0.025)
cf_upper <- apply(cf_mean_draws, 2, quantile, 0.975)

effect_draws <- sweep(cf_mean_draws, 2, obs_mean, "-")
avg_effect_draws <- rowMeans(effect_draws)
mean_effect  <- mean(avg_effect_draws)
ci_effect    <- quantile(avg_effect_draws, c(0.025, 0.975))
ci_effect

#Plots----
time_post <- (T0 + 1):T
plot(time_post, obs_mean, type="l", lwd=2,
     ylim=range(c(obs_mean, cf_lower, cf_upper)),
     ylab="PM10", xlab="Time")

lines(time_post, cf_mean, col="blue", lwd=2)
lines(time_post, cf_lower, col="blue", lty=2)
lines(time_post, cf_upper, col="blue", lty=2)

legend("topright",
       legend=c("Observed", "Counterfactual mean", "95% CI"),
       col=c("black","blue","blue"),
       lty=c(1,1,2),
       lwd=2)

Y_obs_post <- Y[, (T0 + 1):T]
effect_array <- array(NA, dim = c(S, M, T - T0))
for (s in 1:S) {
  effect_array[s,,] <- Y_obs_post - Y_cf[s,,] 
}
station_effect_draws <- apply(effect_array, c(2,3), mean)

lower_effect <- apply(effect_array, c(2,3), quantile, 0.025)
upper_effect <- apply(effect_array, c(2,3), quantile, 0.975)
cum_effect_station <- rowSums(station_effect_draws)
cum_effect_station

lombardia_sf <- gisco_get_nuts(
  year = "2021",
  epsg = "4326",      
  resolution = "1",
  nuts_level = 2,
  country = "IT"
) %>%
  filter(NAME_LATN == "Lombardia") %>%
  st_transform(32632)   

#Random Spatial Effects
w_draws<-post$w
w_mean<-colMeans(w_draws)
summary(w_mean)
df_plot <- data.frame(
  x = data_stan$station_info$UTM_est,
  y = data_stan$station_info$UTM_nord,
  w_values = w_mean
)

points_sf <- st_as_sf(df_plot,
                      coords = c("x", "y"),
                      crs = 32632)

ggplot() + 
  annotation_map_tile(type = "cartolight", zoomin = 0, progress = "none") +
  geom_sf(
    data = lombardia_sf,
    color = "black",
    fill = NA,
    linewidth = 0.8
  ) +
  geom_sf(
    data = points_sf,
    aes(color = w_values),
    size = 4
  ) +
  scale_color_gradient2(
    low = "forestgreen",
    mid = "yellow",
    high = "red",
    midpoint = 0
  ) +
  annotation_scale(location = "br", width_hint = 0.4) +
  theme_minimal() +
  labs(
    title = "Lombardia - Random Spatial Effects",
    color = "Value"
  ) +
  coord_sf()
W_diff <- as.matrix(dist(w_mean))
cor(data_stan$D[upper.tri(data_stan$D)], W_diff[upper.tri(W_diff)])

#Cumulative Effects
df_plot <- data.frame(
  x = data_stan$station_info$UTM_est,
  y = data_stan$station_info$UTM_nord,
  cum_effect = cum_effect_station
)

points_sf <- st_as_sf(df_plot,
                      coords = c("x", "y"),
                      crs = 32632)

ggplot() + 
  annotation_map_tile(type = "cartolight", zoomin = 0, progress = "none") +
  geom_sf(
    data = lombardia_sf,
    color = "black",
    fill = NA,
    linewidth = 0.8
  ) +
  geom_sf(
    data = points_sf,
    aes(color = cum_effect),
    size = 4
  ) +
  scale_color_gradient2(
    low = "forestgreen",
    mid = "yellow",
    high = "red",
    midpoint = 0
  ) +
  annotation_scale(location = "br", width_hint = 0.4) +
  theme_minimal() +
  labs(
    title = "Lombardia - Cumulative Effects",
    color = "Value"
  ) +
  coord_sf()

#Posterior Probability of an Effect
cum_effect_draws <- apply(effect_array, c(1,2), sum)
posterior_prob_station <- apply(cum_effect_draws, 2, function(x) {
  max(mean(x > 0), mean(x < 0))
})

df_plot <- data.frame(
  x = data_stan$station_info$UTM_est,
  y = data_stan$station_info$UTM_nord,
  post_percentage = posterior_prob_station
)

df_plot$evidence <- ifelse(df_plot$post_percentage >= 0.95,
                           "Strong evidence (≥ 0.95)",
                           "Weak/moderate evidence (< 0.95)")

df_plot$evidence <- factor(df_plot$evidence,
                           levels = c("Weak/moderate evidence (< 0.95)",
                                      "Strong evidence (≥ 0.95)"))

points_sf <- st_as_sf(df_plot,
                      coords = c("x", "y"),
                      crs = 32632)

ggplot() + 
  annotation_map_tile(type = "cartolight", zoomin = 0, progress = "none") +
  geom_sf(
    data = lombardia_sf,
    color = "black",
    fill = NA,
    linewidth = 0.8
  ) +
  geom_sf(
    data = points_sf,
    aes(color = evidence),
    size=4
  ) +
    scale_color_manual(
      values = c("gray50", "red")
    ) +
  annotation_scale(location = "br", width_hint = 0.4) +
  theme_minimal() +
  labs(
    title = "Lombardia - Posterior Probability of an Effect",
    color = "Value"
  ) +
  coord_sf()

#Stations
ggplot() + 
  annotation_map_tile(type = "cartolight", zoomin = 0, progress = "none") +
  geom_sf(
    data = lombardia_sf,
    color = "black",
    fill = NA,
    linewidth = 0.8
  ) +
  geom_sf(
    data = points_sf,
    size=3,
    color= 'red'
  ) +
  annotation_scale(location = "br", width_hint = 0.4) +
  theme_minimal() +
  labs(
    title = "Lombardia - Stations"
  ) +
  coord_sf()