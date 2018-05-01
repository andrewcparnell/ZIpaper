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
}
parameters {
  real alpha; // ZI intercept
  vector[K] beta; // Vector of mean coefficients
  real tau1p1; // ZI parameter 1
  real tau2p1; // ZI parameter 2
}
transformed parameters {
  real<lower = 0> theta[N]; // Mean paraemters
  real omega[N]; // ZI parameters
  real tau1; // ZI parameter 1
  real tau2; // ZI parameter 2
  for (n in 1:N) {
    theta[n] = exp(x[n] * beta); // Log link on mean
    omega[n] = alpha - // ZI behaviour with fixed tau
                tau1 * theta[n] + 
                tau2 * log(1 + x[n] * beta);
  }
  tau1 = -tau1p1; 
  tau2 = -tau2p1;
}
model {
  // Priors - mostly vague
  alpha ~ normal(0, 10); 
  for (k in 1:K) {
    beta[k] ~ normal(0, 10);
  }
  tau1p1 ~ beta(1, 1);
  tau2p1 ~ beta(1, 1);

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
