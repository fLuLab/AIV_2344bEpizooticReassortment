// This model describes the detection of novel reassortments during the 2020
// - 2024 H5 Avian influenza panzootic. It is loosely based on previously 
// described abundance models used in ecological studies. 
  
// Input data 
data {
  int<lower=0> N; // Number of observations
  array[N] int<lower=0, upper=12> y; // Observed counts
  int<lower=1> C; // Number of continents
  int<lower=1> Y; // Number of years
  array[N] int<lower=1, upper=C> continent_index; // Continent index for each observation
  array[N] int<lower=1, upper=Y> year_index; // Year index for each observation
  array[N] real cases; // Additional data for cases
  array[N] real sequences; // Additional data for sequences
}


transformed data {
  int K = 12;
}


// Declared parameters
parameters {
  array[C] real<lower=0, upper=1> continent_specific_theta; // Zero-inflation probability stratified by continent
  array[C] real continent_specific_abundance; // Poisson rate stratified by continent
  array[C] real continent_specific_detection; // Detection probability stratified by continent
  
  // 'fixed' effects
  real beta_cases; // Coefficient for additional cases data
  real beta_sequences; // Coefficient for additional sequences data

  // 'random' effects (standardised)
  real<lower=0> sigma_year_abundance; // Standard deviation of year intercepts
  array[Y] real z_year_abundance;
  real<lower=0> sigma_year_detection; // Standard deviation of year intercepts
  array[Y] real z_year_detection;
  
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
  continent_specific_abundance ~ normal(1.5, 1); // on log scale
  beta_cases ~ normal(0, 1);
  z_year_abundance ~ normal(0, 1);
  sigma_year_abundance ~ exponential(0.5);     // Prior for year intercept std deviation
  
  // Detection model
  continent_specific_detection ~ normal(-1, 1.5); // on logit scale
  beta_sequences ~ normal(0, 1);
  z_year_detection ~ normal(0, 1);
  sigma_year_detection ~ exponential(0.5);     // Prior for year intercept std deviation
  
  
  // Zero Inflation Model
  continent_specific_theta ~ beta(2, 5); 
  
  // Loop over number of observations
  for (i in 1:N) {
    int c = continent_index[i]; // Current continent
    int yr = year_index[i]; //Current year
    
    
    // Linear predictors 
    real lambda = exp(continent_specific_abundance[c] + beta_cases * cases[i] + year_abundance[yr] ); 
    real p = inv_logit(continent_specific_detection[c] + beta_sequences * sequences[i] + year_detection[yr]); 
    
    vector[K] lp_pos;
    
    // Likelihood for y[i] = 0
    if (y[i] == 0){
      
      for (j in 1:K) {
        int N_latent = j - 1; 
        lp_pos[j] = poisson_lpmf(N_latent | lambda) + binomial_lpmf(0 | N_latent, p);
      }
      
      target += log_sum_exp({
        bernoulli_lpmf(1 | continent_specific_theta[c]),
        bernoulli_lpmf(0 | continent_specific_theta[c]) + log_sum_exp(lp_pos)
      });
      
      // Likelihood for y[i] > 0
    } else {
      
      for (j in 1:K) {
        int N_latent = y[i] + j - 1;
        lp_pos[j] = poisson_lpmf(N_latent | lambda) + binomial_lpmf(y[i] | N_latent, p);
      }
      
      target += bernoulli_lpmf(0 | continent_specific_theta[c]) + log_sum_exp(lp_pos);
    }

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

