
data {
  
  int<lower=1> N;
  
  int<lower=1> K;
  
  int<lower=1> Z;
  
  matrix[N, K] Y;
  
  matrix<lower=0, upper=1>[K, K] W;
  
}

transformed data{
  
  vector[K] zeros;
  
  matrix<lower = 0>[K, K] D;
  
  vector[K] W_rowsums;
  
  for (i in 1:K) {
    
    W_rowsums[i] = sum(W[i, ]);
    
  }
  
  D = diag_matrix(W_rowsums);
  
  zeros = rep_vector(0, K);
  
}

parameters {
  
  array[N] vector[K] mu_err;
  
  array[N] vector[K] v_err;
  
  vector<lower=0>[K] s_mu;
  
  vector<lower=0>[K] s_v;
  
  vector<lower=-1, upper=1>[K] rho;
  
  vector<lower=0>[K] sigma_eta;
  
  vector[K] phi;
  
  real<lower=0> tau;
  
  real<lower=0, upper=1> alpha;
  
  vector[K] m;
  
  real M;
  
  real<lower=0> sigma_m;

  
}

transformed parameters {
  
  array[N] vector[K] eta_obs;
  
  array[N] vector[K] mu;
  
  array[N] vector[K] v;
  
  array[N] vector[K] nu;
  
  for (k in 1:K) {
    
    mu[1, k] = Y[1, k] + s_mu[k] * mu_err[1, k];
    
    v[1, k] = s_v[k] * v_err[1, k];
    
    nu[1, k] = mu[1, k];
    
    for(t in 2:N) {
      
      mu[t, k] = mu[t-1, k] + v[t-1, k] + s_mu[k] * mu_err[t, k];
      
      v[t, k] = (m[k] + phi[k]) + rho[k] * (v[t-1, k] - (m[k] + phi[k])) + s_v[k] * v_err[t, k];
      
      nu[t, k] = mu[t, k] + mu[t-1, k];
      
    }
    
  }
  
  eta_obs = nu;
  
}

model {
  
  phi ~ multi_normal_prec(zeros, tau * (D - alpha * W));
  
  sum(phi) ~ normal(0, 0.001);
  
  for (k in 1:K) {
    
    m[k] ~ normal(M, sigma_m);

    s_mu[k] ~ normal(0, 0.1);
    
    s_v[k] ~ normal(0, 0.1);
    
    sigma_eta[k] ~ normal(0, 5); 
    
    for (t in 1:N) {
      
      mu_err[t, k] ~ normal(0, 1);
      
      v_err[t, k] ~ normal(0, 1);
      
      if (t <= Z) {
        
        Y[t, k] ~ normal(eta_obs[t, k], sigma_eta);
        
      }
      
    }
    
  }
  
}
