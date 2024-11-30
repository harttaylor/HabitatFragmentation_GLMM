# Data preparation 
library(tidyr)
library(dplyr)
library(lubridate)


# Load datasets (update file paths as needed)
visit_matrix <- read.csv("1_Data/2yearvisitmatrix.csv")
visits_for_det_covs <- read.csv("1_Data/visitsfordetcovs.csv")
# 500m scale 
habitat_metricsA500 <- read.csv("1_Data/500m/habitat_metrics_hyp1.csv")
habitat_metricsB500 <- read.csv("1_Data/500m/habitat_metrics_hyp2.csv")
habitat_metricsC500 <- read.csv("1_Data/500m/habitat_metrics_hyp3.csv")


# Step 1: Reshape detection history to long format
visit_matrix_long <- visit_matrix %>%
  pivot_longer(cols = starts_with("X"),  # Select columns starting with "X"
               names_to = "visit", 
               values_to = "occupancy") %>%
  mutate(
    year_id = as.numeric(substr(visit, 2, 2)),  # Extract year (1 or 2) from column name
    visit_num = as.numeric(substr(visit, 4, 4))  # Extract visit number (1, 2, or 3)
  ) %>%
  select(gisid, occupancy, year_id, visit_num)  # Keep necessary columns

# Step 2: Merge with detection covariates
visit_with_covs <- visit_matrix_long %>%
  inner_join(
    visits_for_det_covs %>%
      select(gisid, offset, year, surveyid, year_id, visit_num, duration, Easting, Northing, date_time), 
    by = c("gisid", "year_id", "visit_num")
  )
head(visit_with_covs)
str(visit_with_covs)

# Step 3: Duplicate habitat metrics for each visit within the corresponding year 
# so we can create visit specific habitat metrics for the long dataframe that icnludes entire det history
# First, create a mapping between surveyid and gisid using the visit dataframe
site_mapping <- visit_with_covs %>%
  select(surveyid, gisid, year) %>%
  distinct()

# Join this mapping to the habitat metrics to add gisid
habitat_with_sites <- habitat_metricsA500 %>%
  left_join(site_mapping, by = c("surveyid", "survey_yea" = "year"))

# Now create a version of habitat metrics without the surveyid
# This prevents duplicate surveyid values in the final join
habitat_for_join <- habitat_with_sites %>%
  select(-surveyid) %>%
  distinct()

# Finally, join the visit data with habitat metrics using gisid and year
methodA_500m <- visit_with_covs %>%
  left_join(habitat_for_join, by = c("gisid", "year")) %>% 
  select(c(-mean_patch_index)) # remove these covariates because they are riddled with NAs

# Step 4: Remove rows with missing values
methodA_500m <- na.omit(methodA_500m)

# Step 5: Check final dataset structure
print(head(methodA_500m))
print(dim(methodA_500m))
str(methodA_500m)
hist(methodA_500m$edge_density)
# Remove the outlier
methodA_500m <- methodA_500m %>% 
  filter(edge_density < 0.1)

# Function to check if standardization is needed
check_scale <- function(x) {
  mean_val <- mean(x, na.rm = TRUE)
  sd_val <- sd(x, na.rm = TRUE)
  if(abs(mean_val) > 1 || sd_val > 1) return(TRUE)
  return(FALSE)
}

# Prepare data
methodA_500m <- methodA_500m %>%
  mutate(
    # Convert datetime
    date_time = as.POSIXct(date_time, format = "%Y-%m-%d %H:%M"),
    
    # Factorize categorical variables
    survey_year = factor(survey_yea),
    duration = factor(duration),
    
    # Create temporal variables
    julian_date = as.numeric(yday(date_time)),
    time_of_day = as.numeric(hour(date_time) + minute(date_time)/60),
    
    # Standardize continuous landscape metrics
    # Using scale() with center and scale parameters for more control
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    #nlsi_std = scale(nlsi, center = TRUE, scale = TRUE)[,1],
    #core_mn_std = scale(core_mn, center = TRUE, scale = TRUE)[,1],
    
    # Standardize temporal variables
    julian_date_std = scale(julian_date, center = TRUE, scale = TRUE)[,1],
    time_of_day_std = scale(time_of_day, center = TRUE, scale = TRUE)[,1],
    offset = offset
  )

# Check distributions and scaling
summary_stats <- methodA_500m %>%
  summarise(across(ends_with('_std'), 
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE))))

print(summary_stats)
range(methodA_500m$edge_density_std)
range(methodA_500m$prop_habitat_std)
# Check distributions 
hist(methodA_500m$duration)
hist(methodA_500m$julian_date)
hist(methodA_500m$time_of_day)
range(methodA_500m$julian_date) # ranges from 126 to 199 day of year
range(methodA_500m$time_of_day) # ranges from 3am to 10:30am, mainly between 4:30am and 8:30 am 

# Summary statistics for original variables
summary(methodA_500m$edge_density)



# Method B 
# Join this mapping to the habitat metrics to add gisid
habitat_with_sites <- habitat_metricsB500 %>%
  left_join(site_mapping, by = c("surveyid", "survey_yea" = "year"))

# Now create a version of habitat metrics without the surveyid
# This prevents duplicate surveyid values in the final join
habitat_for_join <- habitat_with_sites %>%
  select(-surveyid) %>%
  distinct()

# Finally, join the visit data with habitat metrics using gisid and year
methodB_500m <- visit_with_covs %>%
  left_join(habitat_for_join, by = c("gisid", "year")) %>% 
  select(c(-mean_patch_index)) # remove these covariates because they are riddled with NAs

# Step 4: Remove rows with missing values
methodB_500m <- na.omit(methodB_500m)

# Remove the outlier
methodB_500m <- methodB_500m %>% 
  filter(edge_density < 0.1)

# Step 5: Check final dataset structure
print(head(methodB_500m))
print(dim(methodB_500m))
str(methodB_500m)


# Prepare data
methodB_500m <- methodB_500m %>%
  mutate(
    # Convert datetime
    date_time = as.POSIXct(date_time, format = "%Y-%m-%d %H:%M"),
    
    # Factorize categorical variables
    survey_year = factor(survey_yea),
    duration = factor(duration),
    
    # Create temporal variables
    julian_date = as.numeric(yday(date_time)),
    time_of_day = as.numeric(hour(date_time) + minute(date_time)/60),
    
    # Standardize continuous landscape metrics
    # Using scale() with center and scale parameters for more control
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    #nlsi_std = scale(nlsi, center = TRUE, scale = TRUE)[,1],
    #core_mn_std = scale(core_mn, center = TRUE, scale = TRUE)[,1],
    
    # Standardize temporal variables
    julian_date_std = scale(julian_date, center = TRUE, scale = TRUE)[,1],
    time_of_day_std = scale(time_of_day, center = TRUE, scale = TRUE)[,1],
    offset = offset
    # autocov = scale(autocov, center = TRUE, scale = TRUE)[,1]
  )

# Check distributions and scaling
summary_stats <- methodB_500m %>%
  summarise(across(ends_with('_std'), 
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE))))

print(summary_stats)



# Method C 
# Join this mapping to the habitat metrics to add gisid
habitat_with_sites <- habitat_metricsC500 %>%
  left_join(site_mapping, by = c("surveyid", "survey_yea" = "year"))

# Now create a version of habitat metrics without the surveyid
# This prevents duplicate surveyid values in the final join
habitat_for_join <- habitat_with_sites %>%
  select(-surveyid) %>%
  distinct()

# Finally, join the visit data with habitat metrics using gisid and year
methodC_500m <- visit_with_covs %>%
  left_join(habitat_for_join, by = c("gisid", "year")) %>% 
  select(c(-mean_patch_index)) # remove these covariates because they are riddled with NAs

# Step 4: Remove rows with missing values
methodC_500m <- na.omit(methodC_500m)

# Remove the outlier
methodC_500m <- methodC_500m %>% 
  filter(edge_density < 0.1)

# Step 5: Check final dataset structure
print(head(methodC_500m))
print(dim(methodC_500m))
str(methodC_500m)


# Prepare data
methodC_500m <- methodC_500m %>%
  mutate(
    # Convert datetime
    date_time = as.POSIXct(date_time, format = "%Y-%m-%d %H:%M"),
    
    # Factorize categorical variables
    survey_year = factor(survey_yea),
    duration = factor(duration),
    
    # Create temporal variables
    julian_date = as.numeric(yday(date_time)),
    time_of_day = as.numeric(hour(date_time) + minute(date_time)/60),
    
    # Standardize continuous landscape metrics
    # Using scale() with center and scale parameters for more control
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    #nlsi_std = scale(nlsi, center = TRUE, scale = TRUE)[,1],
    #core_mn_std = scale(core_mn, center = TRUE, scale = TRUE)[,1],
    
    # Standardize temporal variables
    julian_date_std = scale(julian_date, center = TRUE, scale = TRUE)[,1],
    time_of_day_std = scale(time_of_day, center = TRUE, scale = TRUE)[,1],
    offset = offset
    # autocov = scale(autocov, center = TRUE, scale = TRUE)[,1]
  )

# Check distributions and scaling
summary_stats <- methodC_500m %>%
  summarise(across(ends_with('_std'), 
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE))))

print(summary_stats)

# save all those datasets
write.csv(methodA_500m, "1_Data/500m/methodA_glmmdata.csv")
write.csv(methodB_500m, "1_Data/500m/methodB_glmmdata.csv")
write.csv(methodC_500m, "1_Data/500m/methodC_glmmdata.csv")


# 150m scale -------------------------------------------------------------------
habitat_metricsA150 <- read.csv("1_Data/150m/habitat_metrics_hyp1.csv")
habitat_metricsB150 <- read.csv("1_Data/150m/habitat_metrics_hyp2.csv")
habitat_metricsC150 <- read.csv("1_Data/150m/habitat_metrics_hyp3.csv")

# Join this mapping to the habitat metrics to add gisid
habitat_with_sites <- habitat_metricsA150 %>%
  left_join(site_mapping, by = c("surveyid", "survey_yea" = "year"))

# Now create a version of habitat metrics without the surveyid
# This prevents duplicate surveyid values in the final join
habitat_for_join <- habitat_with_sites %>%
  select(-surveyid) %>%
  distinct()

# Finally, join the visit data with habitat metrics using gisid and year
methodA_150m <- visit_with_covs %>%
  left_join(habitat_for_join, by = c("gisid", "year")) %>% 
  select(c(-mean_patch_index)) # remove these covariates because they are riddled with NAs

# Step 4: Remove rows with missing values
methodA_150m <- na.omit(methodA_150m)

# Step 5: Check final dataset structure
print(head(methodA_150m))
print(dim(methodA_150m))

# Check distributions 
hist(methodA_150m$edge_density)
range(methodA_150m$edge_density)
summary(methodA_150m$edge_density)
# Remove the outliers if needed
#methodC_1000m <- methodC_1000m %>% 
# filter(edge_density < 0.1)

# Prepare data
methodA_150m <- methodA_150m %>%
  mutate(
    # Convert datetime
    date_time = as.POSIXct(date_time, format = "%Y-%m-%d %H:%M"),
    
    # Factorize categorical variables
    survey_year = factor(survey_yea),
    duration = factor(duration),
    
    # Create temporal variables
    julian_date = as.numeric(yday(date_time)),
    time_of_day = as.numeric(hour(date_time) + minute(date_time)/60),
    
    # Standardize continuous landscape metrics
    # Using scale() with center and scale parameters for more control
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    #nlsi_std = scale(nlsi, center = TRUE, scale = TRUE)[,1],
    #core_mn_std = scale(core_mn, center = TRUE, scale = TRUE)[,1],
    
    # Standardize temporal variables
    julian_date_std = scale(julian_date, center = TRUE, scale = TRUE)[,1],
    time_of_day_std = scale(time_of_day, center = TRUE, scale = TRUE)[,1],
    offset = offset
  )

# Check distributions and scaling
summary_stats <- methodA_150m %>%
  summarise(across(ends_with('_std'), 
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE))))

print(summary_stats)



# Method B 
# Join this mapping to the habitat metrics to add gisid
habitat_with_sites <- habitat_metricsB150 %>%
  left_join(site_mapping, by = c("surveyid", "survey_yea" = "year"))

# Now create a version of habitat metrics without the surveyid
# This prevents duplicate surveyid values in the final join
habitat_for_join <- habitat_with_sites %>%
  select(-surveyid) %>%
  distinct()

# Finally, join the visit data with habitat metrics using gisid and year
methodB_150m <- visit_with_covs %>%
  left_join(habitat_for_join, by = c("gisid", "year")) %>% 
  select(c(-mean_patch_index)) # remove these covariates because they are riddled with NAs

# Step 4: Remove rows with missing values
methodB_150m <- na.omit(methodB_150m)

# Step 5: Check final dataset structure
print(head(methodB_150m))
print(dim(methodB_150m))
str(methodB_150m)

# Check distributions 
hist(methodB_150m$edge_density)
range(methodB_150m$edge_density)
summary(methodB_150m$edge_density)
# Remove the outliers if needed
#methodC_1000m <- methodC_1000m %>% 
# filter(edge_density < 0.1)

# Prepare data
methodB_150m <- methodB_150m %>%
  mutate(
    # Convert datetime
    date_time = as.POSIXct(date_time, format = "%Y-%m-%d %H:%M"),
    
    # Factorize categorical variables
    survey_year = factor(survey_yea),
    duration = factor(duration),
    
    # Create temporal variables
    julian_date = as.numeric(yday(date_time)),
    time_of_day = as.numeric(hour(date_time) + minute(date_time)/60),
    
    # Standardize continuous landscape metrics
    # Using scale() with center and scale parameters for more control
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    #nlsi_std = scale(nlsi, center = TRUE, scale = TRUE)[,1],
    #core_mn_std = scale(core_mn, center = TRUE, scale = TRUE)[,1],
    
    # Standardize temporal variables
    julian_date_std = scale(julian_date, center = TRUE, scale = TRUE)[,1],
    time_of_day_std = scale(time_of_day, center = TRUE, scale = TRUE)[,1],
    offset = offset
    # autocov = scale(autocov, center = TRUE, scale = TRUE)[,1]
  )

# Check distributions and scaling
summary_stats <- methodB_150m %>%
  summarise(across(ends_with('_std'), 
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE))))

print(summary_stats)



# Method C 
# Join this mapping to the habitat metrics to add gisid
habitat_with_sites <- habitat_metricsC150 %>%
  left_join(site_mapping, by = c("surveyid", "survey_yea" = "year"))

# Now create a version of habitat metrics without the surveyid
# This prevents duplicate surveyid values in the final join
habitat_for_join <- habitat_with_sites %>%
  select(-surveyid) %>%
  distinct()

# Finally, join the visit data with habitat metrics using gisid and year
methodC_150m <- visit_with_covs %>%
  left_join(habitat_for_join, by = c("gisid", "year")) %>% 
  select(c(-mean_patch_index)) # remove these covariates because they are riddled with NAs

# Step 4: Remove rows with missing values
methodC_150m <- na.omit(methodC_150m)

# Step 5: Check final dataset structure
print(head(methodC_150m))
print(dim(methodC_150m))
str(methodC_150m)

# Check distributions 
hist(methodC_150m$edge_density)
range(methodC_150m$edge_density)
summary(methodC_150m$edge_density)
# Remove the outliers if needed
#methodC_1000m <- methodC_1000m %>% 
# filter(edge_density < 0.1)


# Prepare data
methodC_150m <- methodC_150m %>%
  mutate(
    # Convert datetime
    date_time = as.POSIXct(date_time, format = "%Y-%m-%d %H:%M"),
    
    # Factorize categorical variables
    survey_year = factor(survey_yea),
    duration = factor(duration),
    
    # Create temporal variables
    julian_date = as.numeric(yday(date_time)),
    time_of_day = as.numeric(hour(date_time) + minute(date_time)/60),
    
    # Standardize continuous landscape metrics
    # Using scale() with center and scale parameters for more control
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    #nlsi_std = scale(nlsi, center = TRUE, scale = TRUE)[,1],
    #core_mn_std = scale(core_mn, center = TRUE, scale = TRUE)[,1],
    
    # Standardize temporal variables
    julian_date_std = scale(julian_date, center = TRUE, scale = TRUE)[,1],
    time_of_day_std = scale(time_of_day, center = TRUE, scale = TRUE)[,1],
    offset = offset
    # autocov = scale(autocov, center = TRUE, scale = TRUE)[,1]
  )

# Check distributions and scaling
summary_stats <- methodC_150m %>%
  summarise(across(ends_with('_std'), 
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE))))

print(summary_stats)

# save all those datasets
write.csv(methodA_150m, "1_Data/150m/methodA_glmmdata.csv")
write.csv(methodB_150m, "1_Data/150m/methodB_glmmdata.csv")
write.csv(methodC_150m, "1_Data/150m/methodC_glmmdata.csv")



# 1000m scale ------------------------------------------------------------------
habitat_metricsA1000 <- read.csv("1_Data/1000m/habitat_metrics_hyp1.csv")
habitat_metricsB1000 <- read.csv("1_Data/1000m/habitat_metrics_hyp2.csv")
habitat_metricsC1000 <- read.csv("1_Data/1000m/habitat_metrics_hyp3.csv")

# Join this mapping to the habitat metrics to add gisid
habitat_with_sites <- habitat_metricsA1000 %>%
  left_join(site_mapping, by = c("surveyid", "survey_yea" = "year"))

# Now create a version of habitat metrics without the surveyid
# This prevents duplicate surveyid values in the final join
habitat_for_join <- habitat_with_sites %>%
  select(-surveyid) %>%
  distinct()

# Finally, join the visit data with habitat metrics using gisid and year
methodA_1000m <- visit_with_covs %>%
  left_join(habitat_for_join, by = c("gisid", "year")) %>% 
  select(c(-mean_patch_index)) # remove these covariates because they are riddled with NAs

# Step 4: Remove rows with missing values
methodA_1000m <- na.omit(methodA_1000m)

# Step 5: Check final dataset structure
print(head(methodA_1000m))
print(dim(methodA_1000m))
str(methodA_1000m)

# Check distributions 
hist(methodA_1000m$edge_density)
range(methodA_1000m$edge_density)
# Remove the outliers if needed
#methodC_1000m <- methodC_1000m %>% 
# filter(edge_density < 0.1)


# Prepare data
methodA_1000m <- methodA_1000m %>%
  mutate(
    # Convert datetime
    date_time = as.POSIXct(date_time, format = "%Y-%m-%d %H:%M"),
    
    # Factorize categorical variables
    survey_year = factor(survey_yea),
    duration = factor(duration),
    
    # Create temporal variables
    julian_date = as.numeric(yday(date_time)),
    time_of_day = as.numeric(hour(date_time) + minute(date_time)/60),
    
    # Standardize continuous landscape metrics
    # Using scale() with center and scale parameters for more control
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    #nlsi_std = scale(nlsi, center = TRUE, scale = TRUE)[,1],
    #core_mn_std = scale(core_mn, center = TRUE, scale = TRUE)[,1],
    
    # Standardize temporal variables
    julian_date_std = scale(julian_date, center = TRUE, scale = TRUE)[,1],
    time_of_day_std = scale(time_of_day, center = TRUE, scale = TRUE)[,1],
    offset = offset
  )

# Check distributions and scaling
summary_stats <- methodA_1000m %>%
  summarise(across(ends_with('_std'), 
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE))))

print(summary_stats)
range(methodA_1000m$edge_density_std)
range(methodA_1000m$prop_habitat_std)

# Summary statistics for original variables
summary(methodA_1000m$edge_density)


# Method B 
# Join this mapping to the habitat metrics to add gisid
habitat_with_sites <- habitat_metricsB1000 %>%
  left_join(site_mapping, by = c("surveyid", "survey_yea" = "year"))

# Now create a version of habitat metrics without the surveyid
# This prevents duplicate surveyid values in the final join
habitat_for_join <- habitat_with_sites %>%
  select(-surveyid) %>%
  distinct()

# Finally, join the visit data with habitat metrics using gisid and year
methodB_1000m <- visit_with_covs %>%
  left_join(habitat_for_join, by = c("gisid", "year")) %>% 
  select(c(-mean_patch_index)) # remove these covariates because they are riddled with NAs


# Step 4: Remove rows with missing values
methodB_1000m <- na.omit(methodB_1000m)

# Check distributions 
hist(methodB_1000m$edge_density)
range(methodB_1000m$edge_density)
# Remove the outliers if needed
#methodC_1000m <- methodC_1000m %>% 
# filter(edge_density < 0.1)


# Step 5: Check final dataset structure
print(head(methodB_1000m))
print(dim(methodB_1000m))
str(methodB_1000m)


# Prepare data
methodB_1000m <- methodB_1000m %>%
  mutate(
    # Convert datetime
    date_time = as.POSIXct(date_time, format = "%Y-%m-%d %H:%M"),
    
    # Factorize categorical variables
    survey_year = factor(survey_yea),
    duration = factor(duration),
    
    # Create temporal variables
    julian_date = as.numeric(yday(date_time)),
    time_of_day = as.numeric(hour(date_time) + minute(date_time)/60),
    
    # Standardize continuous landscape metrics
    # Using scale() with center and scale parameters for more control
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    #nlsi_std = scale(nlsi, center = TRUE, scale = TRUE)[,1],
    #core_mn_std = scale(core_mn, center = TRUE, scale = TRUE)[,1],
    
    # Standardize temporal variables
    julian_date_std = scale(julian_date, center = TRUE, scale = TRUE)[,1],
    time_of_day_std = scale(time_of_day, center = TRUE, scale = TRUE)[,1],
    offset = offset
    # autocov = scale(autocov, center = TRUE, scale = TRUE)[,1]
  )

# Check distributions and scaling
summary_stats <- methodB_1000m %>%
  summarise(across(ends_with('_std'), 
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE))))

print(summary_stats)


# Method C 
# Join this mapping to the habitat metrics to add gisid
habitat_with_sites <- habitat_metricsC1000 %>%
  left_join(site_mapping, by = c("surveyid", "survey_yea" = "year"))

# Now create a version of habitat metrics without the surveyid
# This prevents duplicate surveyid values in the final join
habitat_for_join <- habitat_with_sites %>%
  select(-surveyid) %>%
  distinct()

# Finally, join the visit data with habitat metrics using gisid and year
methodC_1000m <- visit_with_covs %>%
  left_join(habitat_for_join, by = c("gisid", "year")) %>% 
  select(c(-mean_patch_index)) # remove these covariates because they are riddled with NAs

# Step 4: Remove rows with missing values
methodC_1000m <- na.omit(methodC_1000m)

# Check distributions 
hist(methodC_1000m$edge_density)
range(methodC_1000m$edge_density)
# Remove the outliers if needed
#methodC_1000m <- methodC_1000m %>% 
 # filter(edge_density < 0.1)

# Step 5: Check final dataset structure
print(head(methodC_1000m))
print(dim(methodC_1000m))
str(methodC_1000m)


# Prepare data
methodC_1000m <- methodC_1000m %>%
  mutate(
    # Convert datetime
    date_time = as.POSIXct(date_time, format = "%Y-%m-%d %H:%M"),
    
    # Factorize categorical variables
    survey_year = factor(survey_yea),
    duration = factor(duration),
    
    # Create temporal variables
    julian_date = as.numeric(yday(date_time)),
    time_of_day = as.numeric(hour(date_time) + minute(date_time)/60),
    
    # Standardize continuous landscape metrics
    # Using scale() with center and scale parameters for more control
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    #nlsi_std = scale(nlsi, center = TRUE, scale = TRUE)[,1],
    #core_mn_std = scale(core_mn, center = TRUE, scale = TRUE)[,1],
    
    # Standardize temporal variables
    julian_date_std = scale(julian_date, center = TRUE, scale = TRUE)[,1],
    time_of_day_std = scale(time_of_day, center = TRUE, scale = TRUE)[,1],
    offset = offset
    # autocov = scale(autocov, center = TRUE, scale = TRUE)[,1]
  )

# Check distributions and scaling
summary_stats <- methodC_1000m %>%
  summarise(across(ends_with('_std'), 
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE))))

print(summary_stats)

# save all those datasets
write.csv(methodA_1000m, "1_Data/1000m/methodA_glmmdata.csv")
write.csv(methodB_1000m, "1_Data/1000m/methodB_glmmdata.csv")
write.csv(methodC_1000m, "1_Data/1000m/methodC_glmmdata.csv")
