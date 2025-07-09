
// File: BYM2_model_covariates.stan
// The BYM2 model with integrated ICAR functions, covariates, and updated syntax
// Restored bodies for helper functions

functions {
  /**
   * Log probability of the intrinsic conditional autoregressive (ICAR) prior,
   * excluding additive constants.
   * ... (rest of comments) ...
   **/
  real icar_normal_lpdf(vector phi, real spatial_scale,
                        array[] int node1, array[] int node2,
                        int k, array[] int group_size, array[] int group_idx,
                        int has_theta) {
    real lp = -0.5 * dot_self(phi[node1] - phi[node2]);
    int pos = 1;
    if (has_theta == 1) { // BYM2 case
      for (j in 1:k) {
        if (group_size[j] > 1) { // Apply constraint only to non-singletons
          lp += normal_lpdf(sum(phi[segment(group_idx, pos, group_size[j])]) | 0, 0.001 * sqrt(group_size[j]));
        }
        pos += group_size[j];
      }
    } else { // Pure ICAR case (no theta)
      for (j in 1:k) {
        if (group_size[j] > 1) {
           lp += normal_lpdf(sum(phi[segment(group_idx, pos, group_size[j])]) | 0, 0.001 * sqrt(group_size[j]));
        } else { // Singleton in pure ICAR model needs a prior
          lp += normal_lpdf(phi[segment(group_idx, pos, group_size[j])] | 0, spatial_scale);
        }
        pos += group_size[j];
      }
    }
    return lp;
  }

  /**
   * Combine local and global partial-pooling components into the convolved BYM2 term.
   * ... (rest of comments) ...
   */
  vector convolve_bym2(vector phi_tilde, vector theta_tilde,
                       real spatial_scale,
                       int n, int k,
                       array[] int group_size, array[] int group_idx,
                       real rho, vector inv_sqrt_scale_factor) {
    vector[n] convolution;
    int pos = 1;
     for (j in 1:k) {
        int current_size = group_size[j];
        int end_pos = pos + current_size - 1;
        if (current_size == 1) { // Singleton: only scaled IID component
            convolution[pos : end_pos] = spatial_scale * segment(theta_tilde, pos, current_size);
        } else { // Non-singleton: combine scaled spatial and IID
            convolution[pos : end_pos] = spatial_scale * (
                sqrt(rho) * inv_sqrt_scale_factor[j] * segment(phi_tilde, pos, current_size) +
                sqrt(1 - rho) * segment(theta_tilde, pos, current_size)
            );
       }
       pos = end_pos + 1;
    }
    return convolution;
  }

  // Other helper functions (with corrected syntax and bodies restored)
  vector convolve_bym(vector phi, vector theta,
            int n, int k,
            array[] int group_size, array[] int group_idx
            ) {
    vector[n] convolution;
    int pos=1;
    for (j in 1:k) {
       int current_size = group_size[j];
       int end_pos = pos + current_size - 1;
       if (current_size == 1) {
          convolution[pos : end_pos] = segment(theta, pos, current_size);
      } else {
          convolution[pos : end_pos] = segment(phi, pos, current_size) + segment(theta, pos, current_size);
       }
        pos = end_pos + 1;
    }
    return convolution; // Added return
  }

  vector make_phi(vector phi_tilde, real phi_scale,
            vector inv_sqrt_scale_factor,
            int n, int k,
            array[] int group_size, array[] int group_idx
            ) {
    vector[n] phi;
    int pos=1;
    for (j in 1:k) {
      int current_size = group_size[j];
      int end_pos = pos + current_size - 1;
      phi[pos : end_pos] = phi_scale * inv_sqrt_scale_factor[j] * segment(phi_tilde, pos, current_size);
      pos = end_pos + 1;
    }
    return phi; // Added return
  }

  vector make_phi2(vector phi_tilde, vector phi_scale, // phi_scale is vector[k] here
            vector inv_sqrt_scale_factor,
            int n, int k,
            array[] int group_size, array[] int group_idx
            ) {
    vector[n] phi;
    int pos=1;
    for (j in 1:k) {
      int current_size = group_size[j];
      int end_pos = pos + current_size - 1;
      phi[pos : end_pos] = phi_scale[j] * inv_sqrt_scale_factor[j] * segment(phi_tilde, pos, current_size);
      pos = end_pos + 1;
    }
    return phi; // Added return
  }
} // end functions block

data {
  // Structure / Spatial Info
  int<lower=1> n;
  int<lower=1> k;
  array[k] int group_size;
  array[n] int group_idx;
  int<lower=1> n_edges;
  array[n_edges] int<lower=1, upper=n> node1;
  array[n_edges] int<lower=1, upper=n> node2;
  vector[k] inv_sqrt_scale_factor;

  // Covariates
  int<lower=0> p;          // Number of covariates (p=0 allowed)
  matrix[n, p] X;          // Covariate matrix

  // Model data
  array[n] int<lower=0> y;
  vector[n] log_expected; // RENAMED from offset

  // Control flag
  int<lower=0, upper=1> prior_only;
}

transformed data {
  int<lower=0, upper=1> has_theta = 1;
}

parameters {
  // Intercept & Random Effects
  real alpha;
  vector[n] phi_tilde;
  vector[n] theta_tilde;
  real<lower=0> spatial_scale;
  real<lower=0, upper=1> rho;

  // Covariate Effects
  vector[p] beta;          // Coefficients for covariates
}

transformed parameters {
  vector[n] convolution;
  vector[n] eta;

  convolution = convolve_bym2(phi_tilde, theta_tilde, spatial_scale, n, k, group_size, group_idx, rho, inv_sqrt_scale_factor);

  // Calculate linear predictor, including covariates
  eta = log_expected + alpha + convolution; // Base predictor
  if (p > 0) {
      eta += X * beta; // Add covariate effects if p > 0
  }
}

model {
  // Priors
  target += icar_normal_lpdf(phi_tilde | spatial_scale, node1, node2, k, group_size, group_idx, has_theta); // Use target +=
  theta_tilde ~ std_normal();
  spatial_scale ~ normal(0, 1);
  rho ~ beta(0.5, 0.5);
  alpha ~ normal(0, 5);

  // Prior for covariate coefficients
  if (p > 0) {
      beta ~ normal(0, 2.5); // Weakly informative prior, adjust scale as needed
  }

  // Likelihood
  if (!prior_only) {
    y ~ poisson_log(eta);
  }
}

generated quantities {
  vector[n] mu = exp(eta);
  vector[n] resid = to_vector(y) - mu; // Needs y converted to vector
  vector[n] log_lik;
  array[n] int y_rep;

  for (i in 1:n) {
    log_lik[i] = poisson_log_lpmf(y[i] | eta[i]);
    y_rep[i] = poisson_log_rng(eta[i]);
  }
}

