
// Stan code for BYM2 model - Manual ICAR implementation (NO sparse_car)
// Assumes Poisson likelihood with log link and exposure offset
// USES UPDATED STAN ARRAY SYNTAX

data {
  int<lower=1> N;                     // Number of spatial areas
  array[N] int<lower=0> y;            // Observed counts per area
  vector<lower=0>[N] E;               // Expected counts (exposure offset) per area

  // Need adjacency matrix and degrees for manual ICAR
  matrix<lower=0, upper=1>[N, N] adj_matrix; // Adjacency matrix (0/1)
  vector<lower=0>[N] node_degree;     // Number of neighbours for each node
}

transformed data {
  vector[N] log_E = log(E);
  // ICAR Precision matrix Q = D - W
  // D is diagonal matrix of node degrees
  // W is adjacency matrix
  // Note: This Q is singular. We handle constraints in the model block.
  matrix[N, N] Q = diag_matrix(node_degree) - adj_matrix;
}

parameters {
  real alpha;                       // Overall intercept (average log RR)
  real<lower=0> sigma_s;            // Total standard deviation of the spatial effect (theta)
  real<lower=0, upper=1> phi;       // Mixing parameter (proportion of variance for structured effect)
  vector[N] v_raw;                  // Unstructured component (IID Normal)

  // Raw structured component (will apply sum-to-zero implicitly via prior)
  vector[N] u_raw;
}

transformed parameters {
  vector[N] theta; // Combined spatial random effect (structured + unstructured)
  vector[N] eta;   // Linear predictor on the log scale

  // Apply sum-to-zero constraint *explicitly* for transformation:
  // This ensures u_centered has mean zero for the BYM2 combination
  vector[N] u_centered = u_raw - mean(u_raw);

  // Combine the two components according to BYM2 parameterization
  theta = sigma_s * (sqrt(phi) * u_centered + sqrt(1.0 - phi) * v_raw);

  // Linear predictor: log(E) + intercept + spatial random effects
  eta = log_E + alpha + theta;
}

model {
  // --- Priors --- //
  alpha ~ normal(0, 5);
  sigma_s ~ normal(0, 1); // Constrained lower=0 -> Half-Normal(0, 1)
  phi ~ beta(0.5, 0.5);   // Corresponds roughly to rho or proportion spatial

  // Prior for unstructured component
  v_raw ~ std_normal();

  // --- Manual ICAR prior for u_raw ---
  // The ICAR precision matrix Q is singular.
  // Using multi_normal_prec requires Q to be positive semi-definite (which it is).
  // The model implies a precision tau=1, absorbed into sigma_s and phi later.
  // We must handle the rank deficiency (singularity). A common approach
  // is a 'soft' sum-to-zero constraint on u_raw.

  target += multi_normal_prec_lpdf(u_raw | rep_vector(0, N), Q); // ICAR structure prior

  // Soft sum-to-zero constraint: enforces identifiability by anchoring the mean.
  // Adjust the standard deviation (0.001 * sqrt(N)) if needed for sampling stability.
  // Smaller SD enforces the constraint more strongly.
  sum(u_raw) ~ normal(0, 0.001 * sqrt(N));

  // --- Likelihood --- //
  y ~ poisson_log(eta); // Use array directly
}

generated quantities {
  vector[N] mu;               // Expected counts (rate * exposure)
  vector[N] log_rr;           // Log relative risk (intercept + spatial)
  vector[N] rr;               // Relative risk (exp(log_rr))
  vector[N] rr_spatial = exp(theta); // Spatial component RR
  vector[N] u_scaled;         // Contribution of structured effect to theta
  vector[N] v_scaled;         // Contribution of unstructured effect to theta
  array[N] int<lower=0> y_rep;      // Posterior predictive replications
  vector[N] log_lik;          // Pointwise log-likelihood for LOO/WAIC

  // Re-calculate centered u consistently from the sampled u_raw
  vector[N] u_centered_gen = u_raw - mean(u_raw);
  u_scaled = sigma_s * sqrt(phi) * u_centered_gen;
  v_scaled = sigma_s * sqrt(1.0 - phi) * v_raw;

  log_rr = alpha + theta;
  rr = exp(log_rr);
  mu = exp(eta);

  // Posterior predictive checks & Log-likelihood
  for (i in 1:N) {
    y_rep[i] = poisson_log_rng(eta[i]);
    log_lik[i] = poisson_log_lpmf(y[i] | eta[i]);
  }
}


