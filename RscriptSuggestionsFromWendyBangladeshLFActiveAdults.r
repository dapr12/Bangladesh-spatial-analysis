#Here annotate fully, this part contributed by Wendy Olsen,
# Working with Diego Perez Ruiz and Anthony Shing Yan Kwong;
# July 2025, Univ. of Manchester
# We gratefully acknowledge that the  Univ. of Manchester has provided the 
# facilities for this work via the Dept of Social Statistics. 

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

### STEP 1: LOAD and examine SURVEY DATA ###

file_path <- "C:/data/Bangladesh_LF_and_CL_Survey_2013/Dataset/BGD_2013_LFS-By Quarter_STATA/LFS-2013-By Quarter.dta"
lfs_data_factored <- read_dta(file_path) %>% as_factor()

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

### Rename class as educhigh
### Rename q48 as labourstatustype. 
### Rename lf as isinlforce
### Rename nul1 as availtowork

### Tabulate unemployed by isinlforce and check raw N matches the current number of raw N. 
### If not, then alter isinlforce by mutate , casewhen unemployed==1, isinlforce==1. 

### If this is not clear then discuss it with the team. 

###Exploratory Data Analysis
### Aims 1) what % of women are in the Labour Force (LF), and men? 
### 2) What is the FemHhouseholdHead prevalence by UpaZilla upz? 
### 3) Breaking up agegroups, and adding a sex-specific grouping, what is the
### smallest Zilla count you arrive at in this sample? 
### 4) What's the fewest Unemployed sample members in a sex in an 
### age-group in a Zilla? 

# 5) % Unemployed (Weighted) in each Upazilla, by sex.
# 6) % Unemployed (Weighted) in each Upazilla, by age-group by sex.



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
etc. as seen in earlier code versions.
