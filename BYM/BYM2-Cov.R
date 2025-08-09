# 1. Load Required Packages
pkgs <- c("sf", "spdep", "INLA", "rstan", "SpatialEpi", "geostan")
lapply(pkgs, require, character.only = TRUE)

# Configure rstan
rstan_options(javascript = FALSE, auto_write = TRUE)

# 2. Source Helper Functions
source("./BYM/icar-functions.R") # Contains prep_icar_data and scale_c

# 3. Load Data
data(scotland)
df <- scotland$data
sp <- scotland$spatial.polygon
sp_sf <- st_as_sf(sp)
sp_sf$cases <- df$cases
sp_sf$expected <- df$expected
sp_sf$AFF <- df$AFF # Agriculture, Forestry, Fishing covariate

# 4. Prepare Spatial Structure
nb <- spdep::poly2nb(sp, queen = TRUE)
C <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
n <- nrow(df)
stopifnot(nrow(C) == n)

# 5. Prepare Data for Stan ICAR Functions
icar_data <- prep_icar_data(C)

# 6. Calculate Component-Specific Scaling Factor
k <- icar_data$k
scale_factor <- vector(mode = "numeric", length = k)
if (!exists("scale_c")) stop("Error: 'scale_c' function not found.")
for (j in 1:k) {
  g_idx <- which(icar_data$comp_id == j)
  if (length(g_idx) == 1) {
    scale_factor[j] <- 1.0
  } else {
    Cg <- C[g_idx, g_idx, drop = FALSE]
    scale_factor[j] <- scale_c(Cg)
  }
}
if(any(scale_factor <= 0)) {
  warning("Zero/negative scale factors found. Setting inv_sqrt_scale_factor to 1.")
  icar_data$inv_sqrt_scale_factor <- ifelse(scale_factor > 0, 1 / sqrt(scale_factor), 1.0)
} else {
  icar_data$inv_sqrt_scale_factor <- 1 / sqrt(scale_factor)
}
print("Updated inv_sqrt_scale_factor (for Stan):")
print(icar_data$inv_sqrt_scale_factor)

# 7. Prepare Covariate Matrix <--------- NEW SECTION --------->
# Select covariates and create model matrix *without* intercept
# (because 'alpha' is already the intercept in Stan)
covar_formula <- ~ AFF - 1 # Use '- 1' or '+ 0' to exclude intercept
X_unscaled <- model.matrix(covar_formula, data = sp_sf)

# Standardize continuous covariates (highly recommended)
# scale() returns a matrix with attributes, convert back to plain matrix
X_scaled <- scale(X_unscaled)
# Check for columns with zero variance after scaling (e.g., if a covariate was constant)
col_sds <- attr(X_scaled, "scaled:scale")
if(any(col_sds == 0)){
  warning("Some columns in X had zero variance. Removing them.")
  X_scaled <- X_scaled[, col_sds > 0, drop = FALSE]
  if(ncol(X_scaled) == 0) {
    stop("No valid covariates remaining after removing zero-variance columns.")
  }
}
X <- as.matrix(X_scaled) # Ensure it's a plain matrix for Stan

# Get number of covariates
p <- ncol(X)

print(paste("Number of covariates (p):", p))
print("Head of scaled covariate matrix (X):")
print(head(X))
# <----------------------- END NEW SECTION ------------------->

# 8. Combine all data into a list for Stan
stan_data_list <- list(
  # ICAR structure data
  n = n,
  k = icar_data$k,
  group_size = icar_data$group_size,
  group_idx = icar_data$group_idx,
  n_edges = icar_data$n_edges,
  node1 = icar_data$node1,
  node2 = icar_data$node2,
  inv_sqrt_scale_factor = icar_data$inv_sqrt_scale_factor,
  # Covariates <--------- NEW --------->
  p = p,
  X = X,
  # Response and offset (renamed)
  y = df$cases,
  log_expected = log(df$expected), # Use the renamed variable
  # Control flag
  prior_only = 0
)

# Verify data types and structure (optional but recommended)
str(stan_data_list)

# 9. Compile the Stan Model
stan_model_file <- "./BYM/BYM2_model_covariates.stan" # Use the new file name
if (!file.exists(stan_model_file)) {
  stop("Stan model file not found: ", stan_model_file)
}
bym2_compiled <- stan_model(stan_model_file)

# 10. Run the Stan Sampler
fit <- sampling(bym2_compiled,
                data = stan_data_list,
                chains = 4,
                cores = 4,
                iter = 2000,
                warmup = 1000,
                seed = 456) # Use a seed

# 11. Check Diagnostics and View Results
# Include 'beta' in parameters to check
print(fit, pars = c("alpha", "beta", "spatial_scale", "rho"), digits = 3)
summary(fit)$summary[c("alpha", "beta[1]", "spatial_scale", "rho"), ]

# Check Rhat and Neff
stan_rhat(fit)
stan_ess(fit)

# Plot key parameters, including beta
plot(fit, pars = c("alpha", "beta[1]", "spatial_scale", "rho"))
plot(fit, plotfun = "trace", pars = c("alpha", "beta", "spatial_scale", "rho"), inc_warmup = FALSE)

# Posterior predictive checks etc. (remain the same conceptually)
# ...

# Posterior predictive checks
y_rep_samples <- as.matrix(fit, pars = "y_rep")
# Example: Compare mean/sd of observed vs. replicated data
hist(apply(y_rep_samples, 1, mean), main="Mean(y_rep)", xlab="Mean(y_rep)")
abline(v=mean(df$cases), col="red")

hist(apply(y_rep_samples, 1, sd), main="SD(y_rep)", xlab="SD(y_rep)")
abline(v=sd(df$cases), col="red")


print("Script finished.")