rm(list = ls(all = TRUE))

# 1. Load Required Packages
pkgs <- c("sf", "spdep", "INLA", "rstan", "SpatialEpi")
# Install missing packages
# installed_pkgs <- pkgs %in% rownames(installed.packages())
# if (any(installed_pkgs == FALSE)) {
#   install.packages(pkgs[!installed_pkgs])
# }
lapply(pkgs, require, character.only = TRUE)

# Configure rstan
rstan_options(javascript = FALSE, auto_write = TRUE) # auto_write=TRUE is helpful

# 2. Source Helper Functions for ICAR Data Prep
# --- IMPORTANT ---
# Make sure 'icar-functions.R' is in your working directory
# or provide the full path. This file should contain at least
# 'prep_icar_data' and 'scale_c' functions.
source("icar-functions.R")
# --- ---

# 3. Load Data
# Using Scotland lip cancer data from SpatialEpi package
data(scotland)
df <- scotland$data
sp <- scotland$spatial.polygon

# Ensure sf object for plotting later if needed
sp_sf <- st_as_sf(sp)
# Add data to sf object
sp_sf$cases <- df$cases
sp_sf$expected <- df$expected
sp_sf$AFF <- df$AFF # Example covariate (replace/remove if not used)

# Check data
print(head(df))
# plot(sp_sf["cases"])

# 4. Prepare Spatial Structure (Adjacency)
nb <- spdep::poly2nb(sp, queen = TRUE) # Queen contiguity neighbours
C <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE) # Binary adjacency matrix
n <- nrow(df)
stopifnot(nrow(C) == n) # Check matrix dimension

# 5. Prepare Data for Stan ICAR Functions using prep_icar_data
# This function should return: n, k, group_size, group_idx, n_edges, node1, node2, comp_id
# And initializes inv_sqrt_scale_factor to 1.
icar_data <- prep_icar_data(C)

# Check initial scaling factor (should be all 1s)
print("Initial inv_sqrt_scale_factor:")
print(icar_data$inv_sqrt_scale_factor)

# 6. Calculate Component-Specific Scaling Factor (using scale_c)
# This corrects the variance of the ICAR prior for disconnected graphs
# and makes the 'spatial_scale' parameter more interpretable.
k <- icar_data$k
scale_factor <- vector(mode = "numeric", length = k)

# Check if scale_c function exists
if (!exists("scale_c")) {
  stop("Error: The 'scale_c' function (expected from icar-functions.R or INLA context) was not found.")
}

for (j in 1:k) {
  g_idx <- which(icar_data$comp_id == j)
  if (length(g_idx) == 1) {
    # Singletons have scale factor 1 (no scaling needed for prior)
    scale_factor[j] <- 1.0
  } else {
    # Calculate scaling factor for components with size > 1
    Cg <- C[g_idx, g_idx, drop = FALSE] # Subgraph adjacency matrix
    scale_factor[j] <- scale_c(Cg) # This function needs to be defined/sourced
  }
}

# Update the inverse square root scaling factor in the data list
# Handle potential division by zero or issues if scale_factor is non-positive
if(any(scale_factor <= 0)) {
  warning("Some calculated scale factors are zero or negative. Setting inv_sqrt_scale_factor to 1 for these components.")
  icar_data$inv_sqrt_scale_factor <- ifelse(scale_factor > 0, 1 / sqrt(scale_factor), 1.0)
} else {
  icar_data$inv_sqrt_scale_factor <- 1 / sqrt(scale_factor)
}


print("Calculated component scale factors (sqrt(lambda_geom)):")
print(sqrt(scale_factor))
print("Updated inv_sqrt_scale_factor (for Stan):")
print(icar_data$inv_sqrt_scale_factor)

# 7. Combine all data into a list for Stan
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
  # Response and offset
  y = df$cases,
  #offset = log(df$expected),
  log_expected = log(df$expected), # <<< CORRECTED NAME
  # Control flag
  prior_only = 0 # Set to 1 to sample from prior only
)

# Verify data types and structure (optional but recommended)
str(stan_data_list)

# 8. Compile the Stan Model
# Make sure the Stan file 'BYM2_model.stan' is in the working directory
stan_model_file <- "BYM2_model.stan"
if (!file.exists(stan_model_file)) {
  stop("Stan model file not found: ", stan_model_file)
}
bym2_compiled <- stan_model(stan_model_file)

# 9. Run the Stan Sampler
# Adjust chains, cores, iter, warmup as needed
fit <- sampling(bym2_compiled,
                data = stan_data_list,
                chains = 4,
                cores = 4,
                iter = 2000,
                warmup = 1000,
                seed = 123) # for reproducibility

# 10. Check Diagnostics and View Results
print(fit, pars = c("alpha", "spatial_scale", "rho"))
summary(fit)$summary[c("alpha", "spatial_scale", "rho"), ]

# Check Rhat and Neff
stan_rhat(fit)
stan_ess(fit)

# Plot key parameters
plot(fit, pars = c("alpha", "spatial_scale", "rho"))
plot(fit, plotfun = "trace", pars = c("alpha", "spatial_scale", "rho"), inc_warmup = FALSE)
plot(fit, plotfun = "hist", pars = "rho")
plot(fit, plotfun = "hist", pars = "spatial_scale")

# Extract convolution term
convolution_samples <- as.matrix(fit, pars = "convolution")
convolution_mean <- apply(convolution_samples, 2, mean)
convolution_median <- apply(convolution_samples, 2, median)

# Add results to sf object for mapping
sp_sf$bym2_conv_mean <- convolution_mean
sp_sf$bym2_conv_median <- convolution_median

# Simple map of the posterior mean convolution term
plot(sp_sf["bym2_conv_mean"], main = "Posterior Mean BYM2 Convolution")

# Calculate degree of spatial autocorrelation (SA) in the convolution term
# Need Moran's I calculation function (e.g., from spdep or a custom one)
# Example using a simplified Moran's I calculation (requires 'ape' package for Moran.I)
# if (require(ape)) {
#    W_listw <- mat2listw(C) # Convert adjacency matrix to listw object
#    sa <- apply(convolution_samples, 1, function(conv_row) {
#      Moran.I(conv_row, W_listw$weights, W_listw$neighbours, length(W_listw$neighbours), S0=Szero(W_listw))$observed
#    })
#    hist(sa, main = "Moran's I of Convolution Term Samples")
#    abline(v = mean(sa), col = "red", lwd = 2)
# } else {
#    print("Install 'ape' package to calculate Moran's I.")
# }

# Posterior predictive checks
y_rep_samples <- as.matrix(fit, pars = "y_rep")
# Example: Compare mean/sd of observed vs. replicated data
hist(apply(y_rep_samples, 1, mean), main="Mean(y_rep)", xlab="Mean(y_rep)")
abline(v=mean(df$cases), col="red")

hist(apply(y_rep_samples, 1, sd), main="SD(y_rep)", xlab="SD(y_rep)")
abline(v=sd(df$cases), col="red")

print("Script finished.")

