### STEP 0: INSTALL AND LOAD PACKAGES ###

# This script will create a choropleth map of the unemployment rate in Bangladesh
# by division, using the 2013 Labour Force Survey data.

# Install packages if you don't have them already (uncomment the lines below)
# install.packages("tidyverse")
# install.packages("sf")
# install.packages("geodata")
# install.packages("haven")

# Load all the necessary libraries
library(tidyverse) # For data manipulation (dplyr, ggplot2)
library(sf)        # For working with spatial data (shapefiles)
library(geodata)   # To download administrative boundaries
library(haven)     # To read Stata .dta files


### STEP 1: LOAD AND PREPARE SURVEY DATA ###

# --- Load your data ---
# !!! IMPORTANT: Replace this with the correct path to your file !!!
file_path <- "./Bangladesh-Data/Bangladesh_LF_and_CL_Survey_2013/Dataset/LFS-2013-By Quarter.dta"
lfs_data <- read_dta(file_path)

#read in the cleaned dataset, produced by the Lf_tut.qmd script
lfs_data_factored <- readr::read_csv("./Bangladesh-Data/Bangladesh_LF_and_CL_Survey_2013/Dataset/cleaned_lfdata.csv")


# # Convert coded numeric variables (like division codes) into their text labels
# # This turns the numeric code for 'div' (e.g., 10) into its name (e.g., "Barisal")
# lfs_data_factored <- as_factor(lfs_data)


### STEP 2: AGGREGATE DATA TO THE DIVISION LEVEL ----

#run preprocessing script, output data is called division_summary
source("./Code/preprocessing.R")


### STEP 3: GET THE BANGLADESH MAP DATA ----

# Download Bangladesh administrative boundaries (level 1 = Divisions)
# The data will be cached, so it's fast on subsequent runs.
# 'sf' objects are the standard for spatial data in the tidyverse.
bd_map_sf <- gadm(country = "BGD", level = 2, path = tempdir()) %>%
  st_as_sf()

# Check the column names of the map data to confirm 'NAME_1' is the division name
check = sort(division_summary$NAME_2) == sort(bd_map_sf$NAME_2)
change = which(check == F)

#spelling mistake in bd_map_sf. Should be "Brahmanbaria".
bd_map_sf$NAME_2[which(bd_map_sf$NAME_2 == "Brahamanbaria")] <- "Brahmanbaria"
#spelling mistake in division_summary. Should be Kishoreganj.
division_summary$NAME_2[which(division_summary$NAME_2 == "Kishorgonj")] <- "Kishoreganj"
#variant spelling of Cox's Bazar. It's capitalised in bd_map
bd_map_sf$NAME_2[which(bd_map_sf$NAME_2 == "Cox'S Bazar")] <- "Cox's Bazar"

# sort(division_summary$NAME_2)[12]
# sort(bd_map_sf$NAME_2)[12]

#check all names match up now
all(sort(division_summary$NAME_2) == sort(bd_map_sf$NAME_2))


### STEP 4: JOIN YOUR DATA WITH THE MAP ---

# Merge the calculated unemployment rates with the map polygons
# using the division name ('NAME_1') as the common key.
map_data_final <- bd_map_sf %>%
  left_join(division_summary, by = "NAME_2")

#rename NAME_2 to upazilla

map_data_final <- map_data_final %>%
  dplyr::rename(upazilla = NAME_2)

### STEP 5: CREATE AND DISPLAY THE VISUALIZATION MAP ----

# Use ggplot2 and geom_sf to plot the final data
final_map <- ggplot(data = map_data_final) +
  # Color the map polygons based on the 'unemployment_rate' column
  geom_sf(aes(fill = unemployment_rate)) +
  
  # Add labels to the divisions for better readability
  geom_sf_text(aes(label = upazilla), size = 2.5, color = "black", fontface = "bold") +
  
  # Use a color-blind friendly and visually appealing color scale
  scale_fill_viridis_c(option = "magma", direction = -1, na.value = "grey80") +
  
  # Add a descriptive title, subtitle, and caption
  labs(
    title = "Unemployment Rate in Bangladesh (Ages 15+)",
    subtitle = "By Division, based on the 2013 Labour Force Survey",
    caption = "Data Source: Bangladesh LFS 2013 | Map: GADM",
    fill = "Unemployment\nRate (%)" # Legend title
  ) +
  
  # Use a minimal theme that is ideal for maps (removes axes, gridlines etc.)
  theme_void() +
  
  # Customize the theme for better aesthetics
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 5)),
    plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 15)),
    plot.caption = element_text(hjust = 0.5, size = 9, color = "grey40"),
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA) # Set a background color
  )

# Print the final map to the plot viewer
print(final_map)

# To save the map to a file (e.g., as a PNG)
# ggsave("bangladesh_unemployment_map_2013.png", plot = final_map, width = 8, height = 9, dpi = 300)
