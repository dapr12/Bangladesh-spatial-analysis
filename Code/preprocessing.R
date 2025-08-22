#Author: Anthony Kwong

#Data cleaning code for Bangladesh data. 

#Here, we process the lfs_data into a final data set called lf_data, which we will use for the rest of this study.

#1. Load in survey data ----

#every row is an individual
lfs_data = haven::read_dta("./Bangladesh-Data/Bangladesh_LF_and_CL_Survey_2013/Dataset/LFS-2013-By Quarter.dta") %>% as_factor()

#2. Data cleaning ----

#Change age to numeric
lfs_data$age <- as.numeric(lfs_data$age)
#set female=1, male=0
lfs_data$sex_binary <- plyr::revalue(lfs_data$sex,c("1. male" = 0, "2. female" = 1))
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

#filter lf = 1, rename to lf_data since we filtered out students, retired etc

#check the levels of q48
levels(lfs_data$q48)
#add unemployed to q48 as a factor
levels(lfs_data$q48) <- c(levels(lfs_data$q48), "unemployed")
#set unemployed rows for q48
lfs_data$q48[lfs_data$unemployed == 1] <- "unemployed"
#check 
all(which(lfs_data$q48 == "unemployed") == which(lfs_data$unemployed == 1))

#add female househould head variable
lfs_data = lfs_data %>%
  dplyr::mutate(fem_hh_head = case_when(
    sex_binary == 1 & house_head == 1 ~ 1,
    .default = 0
  ))

#rename variables to be more understandable
lfs_data = lfs_data %>%
  dplyr::rename(
    #Rename class as edu_high
    edu_high = class,
    #Rename q48 to labour_status
    labour_status_type = q48,
    avail_to_work = nul1
  )

#age and lf filter
lf_data = lfs_data %>% 
  dplyr::filter(age > 14) %>%
  dplyr::filter(age < 45) %>%
  dplyr::filter(lf == 1)

#3. Computing variables of interest on a zilla (district level) ----

division_summary <- lf_data %>%
  # Group the data by zilla (districts) and compute statistics
  group_by(zl) %>%
  summarise(
    # Calculate the weighted unemployment rate and convert it to a percentage
    unemployment_rate = weighted.mean(as.numeric(unemployed), w = wgt_final, na.rm = TRUE),
    # Percentage of female household head
    fhead = weighted.mean(fem_hh_head, w = wgt_final, na.rm = TRUE),
    #Percentage in rural vs urban
    rural = weighted.mean(urb == "1. Rural", w = wgt_final, na.rm = TRUE),
    #percentage with at least class-x (grade 10)
    junior_highschool = weighted.mean(ifelse(edu_high == "99. others", NA, as.numeric(edu_high) > 10), w = wgt_final, na.rm = TRUE),
    #Percentage with at least class-xii (grade 12 high school), we omit people in 99.others
    highschool = weighted.mean(ifelse(edu_high == "99. others", NA, as.numeric(edu_high) > 12), w = wgt_final, na.rm = TRUE),
    #at least undergrad
    undergrad = weighted.mean(ifelse(edu_high == "99. others", NA, as.numeric(edu_high) > 14), w = wgt_final, na.rm = TRUE)
  ) %>%
  # Rename the 'div' column to 'NAME_2' to match the map data's column name
  rename(NAME_2 = zl)

#Rename zl by removing number and capitalising to match output of gadm
division_summary$NAME_2 <- division_summary$NAME_2 %>%
  str_replace("^[0-9]+\\.\\s*", "") %>%  # remove leading number, dot, and spaces
  str_to_title()                       # capitalise first letter of each word

# Print the resulting summary table to check it
print("Aggregated Unemployment Rate by Division:")
print(division_summary)




