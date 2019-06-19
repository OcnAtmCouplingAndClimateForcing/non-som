data {
  int<lower=0> n_levels;
  int<lower=0> n;
  int variable[n];
  int era[n];
  vector[n] y;
  vector[n] x;
}
parameters {
  vector[n_levels] beta;
  real<lower=0> sigma_resid;
  real<lower=0> sigma_beta;
  real mu_beta;
}
transformed parameters {
  vector[n] pred;
  for(i in 1:n) {
    if(era[i]==1) {
      pred[i] = x[i]*beta[variable[i]];
    } else {
      pred[i] = x[i]*beta[variable[i]];
    }
  }
}
model {
  mu_beta ~ normal(0,1);
  sigma_beta ~ student_t(3,0,2);
  sigma_resid ~ student_t(3,0,2);
  beta ~ normal(mu_beta,sigma_beta);
  y ~ normal(pred, sigma_resid);
}
