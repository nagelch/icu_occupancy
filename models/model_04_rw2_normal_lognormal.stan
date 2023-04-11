
data {
  
  int<lower=1> N;
  
  int<lower=1> K;
  
  int<lower=1> Z;
  
  matrix[N, K] Y;
  
  matrix[N, K] X;
  
  matrix<lower=0, upper=1>[K, K] W;
  
  matrix[N, K] X1;
  
  matrix[N, K] X2;
  
  matrix[N, K] X3;
  
  matrix[N, K] X4;
  
  matrix[N, K] X5;
  
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
  
  array[3] matrix[N, K] O;
  
  O[1, , ] = X;
  
  O[2, , ] = Y;
  
  O[3, , ] = Y;

}

parameters {
  
  array[3, N] vector[K] mu_err;
  
  array[3, N] vector[K] v_err;
  
  array[N] vector[K] b_err;
  
  array[3] vector<lower=0>[K] s_mu;
  
  array[3] vector<lower=0>[K] s_v;
  
  vector<lower=0>[K] s_b;
  
  array[3] vector<lower=-1, upper=1>[K] rho;
  
  vector<lower=0>[K] sigma_eta_m;
  
  vector<lower=0>[K] sigma_eta_log;
  
  vector<lower=0>[K] sigma_eta_obs;
  
  vector<lower=0>[6] sigma_hyp_beta;
  
  array[3] vector<lower=0>[K] sigma_s;
  
  // array[3] vector[K] M;
  
  array[3] vector[K] phi;
  
  vector<lower=0>[3] tau;
  
  vector<lower=0, upper=1>[3] alpha;
  
  array[3, N] vector[K] s;
  
  array[5] vector[K] beta;
  
  array[5] real Beta;
  
  vector[K] zeta;
  
}

transformed parameters {
  
  array[N] vector[K] eta_m;
  
  array[N] vector[K] eta_log;
  
  array[N] vector[K] eta_obs;
  
  array[3, N] vector[K] mu;
  
  array[3, N] vector[K] v;
  
  array[3, N] vector[K] nu;
  
  array[N] vector[K] b;
  
  for (i in 1:3) {
    
    for (k in 1:K) {
      
      if (i < 3) {
        
        mu[i, 1, k] = log(O[i, 1, k]) + s_mu[i, k] * mu_err[i, 1, k];
        
      } else {
        
        mu[i, 1, k] = O[i, 1, k] + s_mu[i, k] * mu_err[i, 1, k];
        
      }
      
      v[i, 1, k] = s_v[i, k] * v_err[i, 1, k];
      
      nu[i, 1, k] = mu[i, 1, k];
      
      for(t in 2:N) {
        
        mu[i, t, k] = mu[i, t-1, k] + v[i, t-1, k] + s_mu[i, k] * mu_err[i, t, k];
        
        v[i, t, k] = phi[i, k] + rho[i, k] * (v[i, t-1, k] - phi[i, k]) + s_v[i, k] * v_err[i, t, k];
        
        nu[i, t, k] = mu[i, t, k] + mu[i, t-1, k];
        
      }
      
    }
    
  }
  
  for (k in 1:K) {
    
    b[1, k] = b_err[1, k];
    
    for (t in 1:N) {
      
      if (t > 1) {
        
        b[t, k] = b[t-1, k] + s_b[k] * b_err[t, k];
        
      }
      
      eta_m[t, k] = nu[1, t, k] + s[1, t, k];
      
      eta_log[t, k] = nu[2, t, k] + s[2, t, k] + 
      
      b[t, k] * nu[1, t, k] +
      
      beta[1, k] * X1[t, k] +
      
      beta[2, k] * X2[t, k] +
      
      beta[3, k] * X3[t, k] + 
      
      beta[4, k] * X4[t, k] +
      
      beta[5, k] * X5[t, k];
      
      eta_obs[t, k] = nu[3, t, k] + s[3, t, k] + zeta[k] * eta_log[t, k];
      
    }
    
  }
  
}

model {
  
  sigma_eta_m ~ gamma(1, 1); 
  
  sigma_eta_log ~ gamma(1, 1); 
  
  sigma_eta_obs ~ gamma(1, 1); 
  
  s_b ~ normal(0, 0.1);
  
  zeta ~ normal(0, 5);
  
  for (i in 1:5) {
    
    Beta[i] ~ normal(0, 5);
    
    sigma_hyp_beta[i] ~ gamma(1, 1);
    
    beta[i, ] ~ normal(Beta[i], sigma_hyp_beta[i]);
    
  } 
  
  for (t in 1:N) {
    
    for (k in 1:K) {
      
      b_err[t, k] ~ normal(0, 1);
      
    }
    
  }
  
  for (i in 1:3) {
    
    phi[i, ] ~ multi_normal_prec(zeros, tau[i] * (D - alpha[i] * W));
    
    for (k in 1:K) {
      
      s_mu[i, k] ~ normal(0, 0.1);
      
      s_v[i, k] ~ normal(0, 0.1);
      
      sigma_s[i, k] ~ gamma(1, 1);
      
      for (t in 1:N) {
        
        mu_err[i, t, k] ~ normal(0, 1);
        
        v_err[i, t, k] ~ normal(0, 1);
        
        if (t >= 7) {
          
          s[i, t, k] ~ normal(-sum(s[i, t-6:t-1, k]), sigma_s[i, k]);
          
        }
        
      }
      
    }
    
  }
  
  for (k in 1:K) {
    
    for (t in 1:N) {
      
      X[t, k] ~ lognormal(eta_m[t, k], sigma_eta_m);
      
      if (t <= Z) {
        
        Y[t, k] ~ lognormal(eta_log[t, k], sigma_eta_log);
        
        Y[t, k] ~ normal(eta_obs[t, k], sigma_eta_obs);
        
      }
      
    }
    
  }
  
}
