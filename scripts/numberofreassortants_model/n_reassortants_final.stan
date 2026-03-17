
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// Input data 
data {
  int<lower=0> N; // Number of observations
  array[N] int<lower=0> y; // Observed counts
  int<lower=0> K; // Upper bound of reassortants
  int<lower=1> C; // Number of continents
  int<lower=1> Y; // Number of years
  array[N] int<lower=1, upper=C> continent_index; // Continent index for each observation
  array[N] int<lower=1, upper=Y> year_index; // Year index for each observation
  array[N] real cases; // Additional data for cases
  array[N] real sequences; // Additional data for sequences
}

// Declared parameters
parameters {
  array[C] real<lower=0, upper=1> continent_specific_theta; // Zero-inflation probability stratified by continent
  array[C] real<lower=0> continent_specific_abundance; // Poisson rate stratified by continent
  array[C] real<lower=0, upper=1> continent_specific_detection; // Detection probability stratified by continent
  
  // 'fixed' effects
  real beta_cases; // Coefficient for additional cases data
  real beta_sequences; // Coefficient for additional sequences data
  real beta_continentcases; // Coefficient for interaction between continent and cases

  // 'random' effects (standardised)
  //real<lower=0> sigma_year_abundance;    
 // array[Y] real year_abundance;      // Random intercepts for each year
  //real<lower=0> sigma_abundance;         // Residual standard deviation
  //real<lower=0> sigma_year_detection;    // Standard deviation of year intercepts
  //array[Y] real year_detection;      // Random intercepts for each year
  //real<lower=0> sigma_detection;         // Residual standard deviation
  
  real<lower=0> sigma_year_abundance; // Standard deviation of year intercepts
  array[Y] real z_year_abundance;
  real<lower=0> sigma_year_detection; // Standard deviation of year intercepts
  array[Y] real z_year_detection;


  //cholesky_factor_corr[2] L_Omega_theta;
 // vector<lower=0>[2] sigma_theta;
}

transformed parameters {
  array[Y] real year_abundance;
  array[Y] real year_detection;

  for (j in 1:Y) {
    year_abundance[j] = sigma_year_abundance * z_year_abundance[j];
    year_detection[j] = sigma_year_detection * z_year_detection[j];
  }
}

// Model
model {
  // Priors
  // Abundance Model
  continent_specific_abundance ~ normal(3, 1.5);
  beta_cases ~ normal(0, 1);
  z_year_abundance ~ normal(0, 1);
  //year_abundance ~ normal(0, sigma_year_abundance);  // Prior for year random intercepts
  sigma_year_abundance ~ exponential(0.5);     // Prior for year intercept std deviation

  // Detection model
  continent_specific_detection ~ beta(1.5,2);//beta(1.5, 1.5);
  beta_sequences ~ normal(0, 1);
  z_year_detection ~ normal(0, 1);
  //year_detection ~ normal(0, sigma_year_detection);  // Prior for year random intercepts
  sigma_year_detection ~ exponential(0.5);     // Prior for year intercept std deviation
  

  // Zero Inflation Model
  continent_specific_theta ~ beta(2, 5); 
  //L_Omega_theta ~ lkj_corr_cholesky(1);
  //sigma_theta ~ normal(0, 1);
  //to_vector(z_theta) ~ normal(0, 1);
  

  // Loop over number of observations
  for (i in 1:N) {
    vector[K] lp;
    int c = continent_index[i]; // Current continent
    int yr = year_index[i]; //Current year

    
    // Linear predictors 
    real lambda = exp(continent_specific_abundance[c] + beta_cases * cases[i] + year_abundance[yr] ); 
    real p = inv_logit(continent_specific_detection[c] + beta_sequences * sequences[i] + year_detection[yr]); 

    // Loop over plausible values of K to marginalise out discrete latent variables
    for (j in 1:K) {
    
      int current_population = y[i] + j - 1;
      
      // Likelihood for y[i] = 0
      if (y[i] == 0){
      
      vector[3] components;
      components[1] = bernoulli_lpmf(1 | continent_specific_theta[c]);
      components[2] = bernoulli_lpmf(0 | continent_specific_theta[c]) + poisson_lpmf(current_population | lambda);
      components[3] = bernoulli_lpmf(0 | continent_specific_theta[c]) + poisson_lpmf(current_population | lambda) + binomial_lpmf(0 | current_population, p);
                      
      lp[j] = log_sum_exp(components);
      
      // Likelihood for y[i] > 0
      } else {
      
      lp[j] = bernoulli_lpmf(0 | continent_specific_theta[c]) + poisson_lpmf(current_population | lambda) + binomial_lpmf(y[i] | current_population, p);
      }
    }
    // Aggregate the probabilities 
    target += log_sum_exp(lp);
  }
}


// Replications for the posterior predictive distribution
generated quantities {
  array[N] int y_rep; 
  array[N] int N_rep;
  
  for (i in 1:N) {
    int c = continent_index[i]; // Current continent
    int yr = year_index[i]; //Current year
    
    // Draw latent abundance from zero-inflated Poisson
    int current_population;
    if (bernoulli_rng(continent_specific_theta[c]) == 1) {
      // Zero inflation: y_rep[i] = 0
      y_rep[i] = 0;
      N_rep[i] = 0;
      
    } else {
      // Simulate from Poisson
      current_population = poisson_rng(exp(continent_specific_abundance[c] + beta_cases * cases[i] + year_abundance[yr] ));
      N_rep[i] = current_population;
      
      // Simulate observed count from binomial
      y_rep[i] = binomial_rng(current_population, inv_logit(continent_specific_detection[c] + beta_sequences * sequences[i]  + year_detection[yr]));
    }
  }
}


