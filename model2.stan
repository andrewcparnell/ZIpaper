functions {
  // Define log probability mass function
  real ZIPo_lpmf(int y, real omega, real theta) {
    real kappa = exp(omega) - 1;
    real lpmf = -log(1 + kappa * exp(-theta)) - theta + 
                  y * log(theta) - lgamma(y + 1);
    if(y == 0) lpmf += omega;
    return(lpmf);
  }
}
data {
  int<lower = 0> N; // Number of observations
  int<lower = 0> K; // Number of covariates
  int<lower = 0> y[N]; // Vector response
  matrix[N, K] x; // Covariate matrix
  vector[N] tau1; // ZI parameter 1
  vector[N] tau2; // ZI parameter 2
}
parameters {
  real alpha; // ZI intercept
  vector[K] beta; // Vector of mean coefficients
}
transformed parameters {
  real<lower = 0> theta[N]; // Mean paraemters
  real omega[N]; // ZI parameters
  for (n in 1:N) {
    theta[n] = exp(x[n] * beta); // Log link on mean
    omega[n] = alpha - // ZI behaviour with fixed tau
                tau1[n] * theta[n] + 
                tau2[n] * log(1 + x[n] * beta);
  }
}
model {
  // Priors - mostly vague
  alpha ~ normal(0, 10); 
  for (k in 1:K) {
    beta[k] ~ normal(0, 10);
  }

  // Likelihood
  for (n in 1:N){
    target += ZIPo_lpmf(y[n] | omega[n], theta[n]);
  }
}
generated quantities {
  // Generate log-likelihood for later loo comparison
  real log_lik[N];
  for (n in 1:N){
    log_lik[n] = ZIPo_lpmf(y[n] | omega[n], theta[n]);
  }
}
