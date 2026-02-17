data {
  int<lower=1> M;                 
  int<lower=1> P;                 
  int<lower=2> T;                 
  matrix[M, T] Y;                 
  array[T] matrix[M, P] X;        
  matrix[M, M] D;                 
  real<lower=0> ell;              
  int<lower=1> N_pairs;
  int pair_i[N_pairs];
  int pair_j[N_pairs];
}

parameters {
  real<lower=-1, upper=1> alpha;   
  vector[P] beta;

  real<lower=0> sigma;                   

  real<lower=0> sigma_w;                 
  vector[M] z_w;                          
  
  vector<lower=-1, upper=1>[N_pairs] rho_raw;
}

transformed parameters {
  vector[M] w;
  matrix[M, M] Sw;
  matrix[M, M] Omega;
  matrix[M, M] L_Omega;

  for (i in 1:M)
    for (j in 1:M)
      Sw[i, j] = exp(-D[i, j] / ell);

  w = sigma_w * cholesky_decompose(Sw) * z_w;

  Omega = diag_matrix(rep_vector(1.0 + 1e-6, M));

  for (n in 1:N_pairs) {
    int i = pair_i[n];
    int j = pair_j[n];
    Omega[i, j] = rho_raw[n];
    Omega[j, i] = rho_raw[n];
  }

  L_Omega = cholesky_decompose(Omega);
}

model {
  
  beta ~ normal(0, 1);
  alpha ~ uniform(-1, 1);

  sigma   ~ normal(0, 0.5);
  sigma_w ~ inv_gamma(2, 1);

  z_w ~ normal(0, 1);
  rho_raw ~ normal(0, 0.3); 

  
  for (t in 2:T) {
    vector[M] mu;
    mu = alpha * Y[, t - 1] + X[t] * beta + w;
    Y[, t] ~ multi_normal_cholesky(mu, sigma * L_Omega);
  }
}

