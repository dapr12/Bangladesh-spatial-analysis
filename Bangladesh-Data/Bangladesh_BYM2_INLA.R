################################################################################
#
# Bayesian Spatial Analysis of Unemployment in Bangladesh
#
# Author:       Diego Perez Ruiz
# Date:         1st of July, 2025
# Institution:  The University of Manchester, Department of Statistics
#
#-------------------------------------------------------------------------------
#
# Description:
# This script performs a Bayesian spatial analysis of unemployment rates across
# the districts of Bangladesh using the 2013 Labour Force Survey (LFS) data.
# It implements a Besag-York-Mollié 2 (BYM2) model using R-INLA to account for
# spatial autocorrelation and estimate the effect of urbanization on unemployment.
#
# The script covers the full workflow:
# 1. Data Loading and Pre-processing
# 2. Aggregation of survey microdata to the district level
# 3. Creation of a spatial adjacency matrix
# 4. Definition and execution of the BYM2 model
# 5. Visualization of raw data, smoothed model estimates, and spatial effects.
#
################################################################################


### STEP 0: SETUP - INSTALL AND LOAD PACKAGES ###
required_packages <- c("tidyverse", "sf", "spdep", "haven", "geodata", "viridis")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}
if (!require("INLA", character.only = TRUE)) {
  install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
}
library(INLA)
library(spdep)

### STEP 1: LOAD SURVEY DATA ###
file_path <- "/Users/user/Desktop/Data Science Models/Bangladesh -Data/Bangladesh_LF_and_CL_Survey_2013/Dataset/LFS-2013-By Quarter.dta"
lfs_data_factored <- read_dta(file_path) %>% as_factor()

### STEP 2: AGGREGATE DATA (WITH CORRECT CASE AND NAME CLEANING) ###
district_model_data <- lfs_data_factored %>%
  
  # --- FIX #1: USE LOWERCASE FOR FILTERING ---
  # Filter for the relevant population using the correct lowercase label
  filter(lf15 == "1" & !is.na(wgt_final)) %>%
  
  group_by(zl_name) %>%
  summarise(
    # --- FIX #2: USE LOWERCASE IN CALCULATIONS ---
    unemployed_count = sum((unemp15 == "unemployed") * wgt_final, na.rm = TRUE),
    lf_count = sum(wgt_final, na.rm = TRUE),
    prop_urban = weighted.mean(urb == "urban", w = wgt_final, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    unemployed_count = round(unemployed_count),
    lf_count = round(lf_count),
    
    # --- FIX #3: KEEP THE ROBUST NAME CLEANING ---
    NAME_2 = str_remove(zl_name, "^[0-9]+\\.\\s"),
    NAME_2 = str_to_title(NAME_2),
    NAME_2 = case_when(
      NAME_2 == "Cox's Bazar"   ~ "Coxs Bazar",
      NAME_2 == "Brahmanbaria"  ~ "Brahamanbaria",
      NAME_2 == "Nawabganj"     ~ "Chapai Nawabganj",
      TRUE                      ~ NAME_2
    )
  ) %>%
  filter(lf_count > 0) %>%
  dplyr::select(-zl_name)

### STEP 3: PREPARE SPATIAL DATA AND ADJACENCY ###
bd_map_district <- gadm(country = "BGD", level = 2, path = tempdir()) %>% st_as_sf()

# The join will now succeed because the names and cases match.
model_sf <- bd_map_district %>%
  inner_join(district_model_data, by = "NAME_2") %>%
  mutate(district_id = 1:n())

print(paste("Successfully joined", nrow(model_sf), "districts. This should be close to 64."))

adj.mat <- poly2nb(model_sf)
nb2INLA("map.adj", adj.mat)
g <- inla.read.graph(filename = "map.adj")

### STEP 4: DEFINE AND RUN THE BYM2 MODEL ###
# This step remains the same
formula <- unemployed_count ~ prop_urban +
  f(district_id, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(
      phi = list(prior = "pc", param = c(0.5, 2/3)),
      prec = list(prior = "pc.prec", param = c(1, 0.01))
    ))

bym2_model <- inla(
  formula,
  family = "binomial",
  Ntrials = lf_count,
  data = as.data.frame(model_sf),
  control.predictor = list(compute = TRUE, link = 1),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
)

print("--- Model Summary ---")
summary(bym2_model)

### STEP 5: VISUALIZE MODEL RESULTS ###
# This step remains the same
model_sf$raw_rate <- (model_sf$unemployed_count / model_sf$lf_count) * 100
model_sf$fitted_rate <- bym2_model$summary.fitted.values$mean * 100
model_sf$spatial_structured <- bym2_model$summary.random$district_id$mean[1:g$n]

map1 <- ggplot(model_sf) +
  geom_sf(aes(fill = raw_rate), color = "white", size = 0.1) +
  scale_fill_viridis_c(option = "magma") +
  labs(title = "Raw Weighted Unemployment Rate by District", fill = "Rate (%)") +
  theme_void()

map2 <- ggplot(model_sf) +
  geom_sf(aes(fill = fitted_rate), color = "white", size = 0.1) +
  scale_fill_viridis_c(option = "magma") +
  labs(title = "Model-Smoothed Unemployment Rate (BYM2)", fill = "Fitted Rate (%)") +
  theme_void()

map3 <- ggplot(model_sf) +
  geom_sf(aes(fill = spatial_structured), color = "white", size = 0.1) +
  scale_fill_viridis_c(option = "cividis") +
  labs(title = "Structured Spatial Effect (Latent Geographic Risk)", fill = "Log-Odds Effect") +
  theme_void()

print(map1)
print(map2)
print(map3)
