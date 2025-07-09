#Exploratory data analysis script. 

#Author: Anthony Kwong

pacman::p_load("tidyverse", "sf", "spdep", "haven", "geodata", "viridis", "INLA")

### STEP 1: load survey data ----
 
lfs_data = haven::read_dta("./Bangladesh-Data/Bangladesh_LF_and_CL_Survey_2013/Dataset/LFS-2013-By Quarter.dta") %>%
  as_factor()

# Data wrangling
lfs_data$age <- as.numeric(lfs_data$age)
#set female=1, male=0
lfs_data$sex_binary <- plyr::revalue(lfs_data$sex,
                                              c("1. male" = 0, "2. female" = 1))
#add household head variable. 1 if household head, 0 otherwise
unique(lfs_data$rel)
lfs_data = lfs_data%>%
  dplyr::mutate(house_head = ifelse(rel == "1. household head", 1, 0))

#add unemployed variable. based on nul1.
lfs_data = lfs_data %>%
  mutate(unemployed = case_when(
    nul1 == 1 ~ 1,
    .default = 0
  ))
#add unemployed to q48 as a factor
levels(lfs_data$q48) <- c(levels(lfs_data$q48), "unemployed")
lfs_data$q48[lfs_data$unemployed == 1] <- "unemployed"
#check 
all(which(lfs_data$q48 == "unemployed") == which(lfs_data$unemployed == 1))
#rename variables to be more understandable
lfs_data = lfs_data %>%
  dplyr::rename(
    #Rename class as edu_high
    edu_high = class,
    #Rename q48 to labour_status
    labour_status_type = q48,
    in_labforce = lf,
    avail_to_work = nul1
  )

#create table of labourforce and unemployment 
lfs_data %>%
  count(in_labforce, unemployed)

#change in_labforce to 0 if they are unemployed, otherwise it is 1 (Wendy to check)
lfs_data = lfs_data %>%
  dplyr::mutate(in_labforce = case_when(
    unemployed == 1 ~ 0, 
    .default = 1
  ))

###Exploratory Data Analysis----
### Aims 1) what % of women are in the Labour Force (LF), and men? 

#Recall that female = 1
lfs_data$sex_binary
lfs_data %>%
  group_by(sex) %>%
  summarise(
    n_in_lf = sum(in_labforce == 1, na.rm = TRUE),
    n_total = n(),
    pct_in_lf = 100 * n_in_lf / n_total
  )

# 2) What is the FemHhouseholdHead prevalence by UpaZilla (districts) upz?----

#add fem_hh_head variable
lfs_data = lfs_data %>%
  dplyr::mutate(fem_hh_head = case_when(
    sex_binary == 1 & house_head == 1 ~ 1,
    .default = 0
  ))

#check discrepancy between upz_name(500 unique values) vs upz (95 unique values). 
lfs_data %>%
  group_by(upz_name) %>%
  summarise(
    total = n(),
    num_fhead = sum(fem_hh_head == 1),
    pct_fhead = num_fhead/total
  )


# 3) Breaking up agegroups, and adding a sex-specific grouping, what is the
# smallest Zilla count you arrive at in this sample? 

lfs_data = lfs_data %>%
  dplyr::mutate(age_group = case_when(
    age < 15 ~ "<15",
    age >= 15 & age <= 44 ~ "15–44",
    age > 44 ~ ">44"
  ))

lfs_data <- lfs_data %>%
  mutate(age_group = factor(age_group, levels = c("<15", "15–44", ">44"), ordered = TRUE))

zilla_counts = lfs_data %>%
  group_by(upz_name, age_group) %>%
  summarise(
    zilla_count = n()
  )

zilla_counts %>%
  arrange(zilla_count)

### 4) What's the fewest Unemployed sample members in a sex in an 
### age-group in a Zilla? 

uemp_tab = zilla_counts2 = lfs_data %>%
  dplyr::group_by(age_group, sex) %>%
  summarise(
    zilla_count = n(),
    unemp = sum(unemployed == 1),
    uemp_rate = unemp/zilla_count,
    .groups = "drop"
  )

uemp_tab %>% 
  arrange(zilla_count)

# 5) % Unemployed (Weighted) in each Upazilla, by sex.

uemp_counts = lfs_data %>%
  dplyr::group_by(sex) %>%
  summarise(
    zilla_count = n(),
    unemp = sum(unemployed == 1),
    uemp_rate = unemp/zilla_count,
    .groups = "drop"
  )


# 6) % Unemployed (Weighted) in each Upazilla, by age-group by sex.

#found here
uemp_tab

#old comments from Wendy, to be removed later.

### Rename class as educhigh
### Rename q48 as labourstatustype. 
### Rename lf as isinlforce
### Rename nul1 as availtowork


###Select the key groups for analysis
###      KEEP GROUP 1
###Employed and otherwise Labour-Force active people
### %>% lf==1
### ADD GROUP 2 the unemployed: 
##  add group those cases where nul1==1
### Note:  nul1 is an unemployment variable, type #1. 
### Then add a factor level to q48, for the combined list of cases, mutate
###  so that a factor level 10 unemployed is added to it. 
### This gives levels 1 to 10 plus 99, ie 11 levels. 
### Now you can create a vector for all, unemployed.
### unemployed = 0 casewhen q48 ==10 unemployed==1. 


### Tabulate unemployed by isinlforce and check raw N matches the current number of raw N. 
### If not, then alter isinlforce by mutate , casewhen unemployed==1, isinlforce==1. 

### If this is not clear then discuss it with the team. 


### STEP 2: AGGREGATE DATA (WITH CORRECT CASE AND NAME CLEANING) ----
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

### STEP 3: PREPARE SPATIAL DATA AND ADJACENCY ----
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
