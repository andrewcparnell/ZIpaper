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
  int<lower = 0> N_tau; // Number of levels of tau covariate
  int<lower = 0> y[N]; // Vector response
  int<lower = 0> x_k[N]; // Vector response
  matrix[N, K] x; // Covariate matrix
}
parameters {
  real alpha; // ZI intercept
  vector[K] beta; // Vector of mean coefficients
  real tau1p1[N_tau]; // ZI parameter 1
  real tau2p1[N_tau]; // ZI parameter 2

}
transformed parameters {
  real<lower = 0> theta[N]; // Mean paraemters
  real omega[N]; // ZI parameters
  real tau1[N_tau]; // ZI parameter 1
  real tau2[N_tau]; // ZI parameter 2
  for (n in 1:N) {
    theta[n] = exp(x[n] * beta); // Log link on mean
    omega[n] = alpha - // ZI behaviour with fixed tau
                tau1[x_k[n]] * theta[n] + 
                tau2[x_k[n]] * log(1 + x[n] * beta);
  }
  for (j in 1:N_tau) {
    tau1[j] = tau1p1[j] - 1; 
    tau2[j] = tau2p1[j] - 1;
  }
}
model {
  // Priors - mostly vague
  alpha ~ normal(0, 10); 
  for (k in 1:K) {
    beta[k] ~ normal(0, 10);
  }
  for (j in 1:N_tau) {
    tau1p1[j] ~ beta(0.01, 0.01);
    tau2p1[j] ~ beta(0.01, 0.01);
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
