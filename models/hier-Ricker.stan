//Dynamic Linearized Ricker Model
//  This model is fit to multiple populations and estimates separate time-varying coefficients each region. 

//Notes:
  //  matrix[3, 3] m[6, 7] - m to be a two-dimensional array of size 6 x 7,
//    containing values that are 3 x 3 matrices. 

data {
  
  int<lower=0> S; //number of stocks
  int N[S];  //number of SR observations for each stock
  int<lower=0> maxN;  //maximum number of observations across stocks
  int years[S, maxN]; //Pointer vector for brood year references.
  matrix[S,maxN] ln_rps;
  matrix[S,maxN] spawn;
  int era[S,maxN]; //Indicator for before/after 89
  
  // Covariates
  matrix[S,maxN] covar;
  
  int<lower=0> R; //number of regions
  int region[S]; //region pointer
  
}

parameters {
  // Coefficients: Group Level
  vector[R] mu_beta;
  real<lower=0> sigma_beta[R];
  // Coefficients: Stock Level
  vector[S] beta;
  
  // Coef. Ratios: Group Level
  vector[R] mu_ratio;
  real<lower=0> sigma_ratio[R];
  // Coefficients: Stock Level
  vector[S] ratio; 
  
  //Ricker Params
  real ricker_mu_alpha;
  real<lower=0> ricker_sigma_alpha;
  
  real ricker_alpha[S]; //Note, could treat as RE by drawing from common Normal() where mean/var are estimates...
  real<lower=0> ricker_beta[S]; //No sense in structuring as RE, given scale issues. 
  //Variances
  real<lower=0> sigma_resid[S];
  //Autoregressive term for errors
  real phi;
  vector[S] init_residual;
}

transformed parameters {
  // Predited Log-R/S
  vector[maxN] pred[S]; // [S,maxN] - vector of length S, each element of which is a vector of length maxN
  
  //Iniital Residual for autocorrelated error structure
  vector[maxN] residual[S];
  
  // Convert Ratios
  vector[S] exp_ratio;
  
  //Generate Predicted values
  for(s in 1:S) {
      exp_ratio[s] = exp(ratio[s]);
    
      //Alternatively, we may be able to do this with an if statement.... But it makes the JAGSer inside of me cry :(
        for(n in 1:maxN) {
          if(n<=N[s]) {
            // Calculate Residual for AR(1) Error
            if(n==1) {
              residual[s,n] = init_residual[s];
            }else {
              residual[s,n] = ln_rps[s,n-1] - pred[s,n-1];
            }
            
            // Calculate predicted ln(R/S)
            if(era[s,n]==1) {
              pred[s,n] = ricker_alpha[s] - ricker_beta[s]*spawn[s,n] + covar[s,n]*beta[s]  + 
                            phi * residual[s,n];
            }else {
              pred[s,n] = ricker_alpha[s] - ricker_beta[s]*spawn[s,n] + covar[s,n]*(beta[s]*exp_ratio[s])  + 
                            phi * residual[s,n];
            }

                                                         
          }else {
            pred[s,n] = 0;
            residual[s,n] = 0;
          }
        }// next n
    }//next s
  }
  
model {
  //Priors
  phi ~ normal(0,sqrt(10)); //uniform(-0.99,0.99); //Autoregressive term for errors
  
  // Ricker alpha hyperparameters
  ricker_mu_alpha ~ normal(0,5);
  ricker_sigma_alpha ~ normal(0,5);
  
  for(s in 1:S) { // Stocks
    ricker_alpha[s] ~ normal(ricker_mu_alpha,ricker_sigma_alpha);
    ricker_beta[s] ~ normal(0,0.001);
    sigma_resid[s] ~ normal(0,5);//cauchy(0,5);
    init_residual[s] ~ normal(0,sqrt(10));//normal(0,(sigma_oe[s]/(1-phi^2))); // normal(0,10); //Initial residual
    
    //Coefficients
    beta[s] ~ normal(mu_beta[region[s]], sigma_beta[region[s]]);
    ratio[s] ~ normal(mu_ratio[region[s]], sigma_ratio[region[s]]);
  }//next s
    

    
  for(r in 1:R) { //Regions
    mu_beta[r] ~ normal(0,1);
    mu_ratio[r] ~ normal(0,1);
    sigma_beta[r] ~ normal(0,5);//cauchy(0,5);
    sigma_ratio[r] ~ normal(0,5);//cauchy(0,5);
  }
    
    
  //Likelihood
  for(s in 1:S) {
    for(n in 1:N[s]) {
      ln_rps[s,n] ~ normal(pred[s,n], sigma_resid[s]);
      // ln_rps[s,n] ~ normal(pred[s,n] + phi * residual[s,n], sigma_resid[s]);
    }//next n
  }//next s
}
    
generated quantities {
  vector[maxN] log_lik[S];
  for(s in 1:S) {
    // for(n in 1:N[s]){
    for(n in 1:maxN) {
      if(n<=N[s]) {
        log_lik[s,n] = normal_lpdf(ln_rps[s,n] | pred[s,n], sigma_resid[s]);
        // log_lik[s,n] = normal_lpdf(ln_rps[s,n] | pred[s,n] + phi * residual[s,n], sigma_oe[s]);
      }else {
        log_lik[s,n] = 0;
      }
    }//next n
  }//next s
}