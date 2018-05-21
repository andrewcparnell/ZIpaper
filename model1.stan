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
  real tau1; // ZI parameter 1
  real tau2; // ZI parameter 2
  int<lower = 0> N_cov;
  matrix[N_cov, K] x_cov;
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
                tau1 * theta[n] + 
                tau2 * log(1 + x[n] * beta);
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
  real<lower = 0> theta_cov[N_cov];
  real omega_cov[N_cov];
  real lpit0[N_cov];
  real lpi0[N_cov];
  for (n in 1:N){
    log_lik[n] = ZIPo_lpmf(y[n] | omega[n], theta[n]);
  }
  for (n in 1:N_cov) {
    theta_cov[n] = exp(x_cov[n] * beta);
    omega_cov[n] = alpha - 
                tau1 * theta_cov[n] + 
                tau2 * log(1 + x_cov[n] * beta);
    lpit0[n] = ZIPo_lpmf(0 | omega_cov[n], theta_cov[n]);
    lpi0[n] = -x_cov[n] * beta;
  }
}
