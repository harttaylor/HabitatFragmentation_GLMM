# GLM Data preparation --------------------------------------------------------------
library(tidyverse)
library(brms)
library(bayesplot)
library(loo)
library(viridis)
library(sf)
library(patchwork)
library(mgcv)  # for GAM fitting
library(sjPlot)

# Load datasets (update file paths as needed)
visit_matrix <- read.csv("1_Data/2yearvisitmatrix.csv")
visits_for_det_covs <- read.csv("1_Data/visitsfordetcovs.csv")
habitat_metricsA500 <- read.csv("~/Chapter 2/Chapter 2 Analysis/Extract_patch_metrics/0_data/2_combined/500habitat_metrics_hyp1.csv")
habitat_metricsB500 <- read.csv("~/Chapter 2/Chapter 2 Analysis/Extract_patch_metrics/0_data/2_combined/500habitat_metrics_hyp2.csv")
habitat_metricsC500 <- read.csv("~/Chapter 2/Chapter 2 Analysis/Extract_patch_metrics/0_data/2_combined/500habitat_metrics_hyp3.csv")
#head(visit_matrix)
#head(visits_for_det_covs)
#head(habitat_metrics)
#str(visit_matrix)
#str(visits_for_det_covs)
#str(habitat_metrics)

# First, let's merge the QPAD offsets into our main dataset
# Step 1: Clean up visit matrix and get maximum detection per site-year
visit_matrix_long <- visit_matrix %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "visit_id", 
    values_to = "occupancy"
  ) %>%
  # Extract year and visit number using substr
  mutate(
    year = as.numeric(substr(visit_id, 2, 2)),
    visit_num = as.numeric(substr(visit_id, 4, 4))
  ) %>%
  group_by(gisid, year) %>%
  summarize(occupancy = max(occupancy)) %>%  # Take maximum detection per site-year
  ungroup()

# Step 2: Get detection covariates and offset for first visit of each site-year
det_covs <- visits_for_det_covs %>%
  mutate(
    date_time = as.POSIXct(date_time, format = "%Y-%m-%d %H:%M"),
    julian_date = yday(date_time),
    time_of_day = hour(date_time) + minute(date_time)/60
  ) %>%
  # Keep first visit info and offset for each site-year
  group_by(gisid, year_id) %>%
  slice(1) %>%
  ungroup() %>%
  select(gisid, offset, surveyid, year_id, duration, 
         Easting, Northing, date_time, julian_date, time_of_day, offset)

# Repeat for each measurement method: 
# Step 3: Clean habitat metrics
habitat_dataA500 <- habitat_metricsA500 %>%
  select(surveyid, habitat_amount, num_patches, edge_density, survey_yea) %>%
  distinct()

habitat_dataB500 <- habitat_metricsB500 %>%
  select(surveyid, habitat_amount, num_patches, edge_density, survey_yea) %>%
  distinct()

habitat_dataC500 <- habitat_metricsC500 %>%
  select(surveyid, habitat_amount, num_patches, edge_density, survey_yea) %>%
  distinct()

# Step 4: Merge datasets
full_dataA500 <- visit_matrix_long %>%
  left_join(
    det_covs,
    by = c("gisid", "year" = "year_id")
  ) %>%
  left_join(
    habitat_dataA500,
    by = c("surveyid")
  )

full_dataB500 <- visit_matrix_long %>%
  left_join(
    det_covs,
    by = c("gisid", "year" = "year_id")
  ) %>%
  left_join(
    habitat_dataB500,
    by = c("surveyid")
  )

full_dataC500 <- visit_matrix_long %>%
  left_join(
    det_covs,
    by = c("gisid", "year" = "year_id")
  ) %>%
  left_join(
    habitat_dataC500,
    by = c("surveyid")
  )


# Step 5: Clean and standardize covariates
method_A_data500 <- full_dataA500 %>%
  # Standardize continuous covariates and create factors
  mutate(
    # Standardize continuous variables
    habitat_amount_std = as.vector(scale(habitat_amount)),
    num_patches_std = as.vector(scale(num_patches)),
    edge_density_std = as.vector(scale(edge_density)),
    julian_date_std = as.vector(scale(julian_date)),
    time_of_day_std = as.vector(scale(time_of_day)),
    
    # Create factors
    survey_duration = factor(duration),
    survey_year = factor(survey_yea),
    gisid = factor(gisid)
  )

method_B_data500 <- full_dataB500 %>%
  # Standardize continuous covariates and create factors
  mutate(
    # Standardize continuous variables
    habitat_amount_std = as.vector(scale(habitat_amount)),
    num_patches_std = as.vector(scale(num_patches)),
    edge_density_std = as.vector(scale(edge_density)),
    julian_date_std = as.vector(scale(julian_date)),
    time_of_day_std = as.vector(scale(time_of_day)),
    
    # Create factors
    survey_duration = factor(duration),
    survey_year = factor(survey_yea),
    gisid = factor(gisid)
  )

method_C_data500 <- full_dataC500 %>%
  # Standardize continuous covariates and create factors
  mutate(
    # Standardize continuous variables
    habitat_amount_std = as.vector(scale(habitat_amount)),
    num_patches_std = as.vector(scale(num_patches)),
    edge_density_std = as.vector(scale(edge_density)),
    julian_date_std = as.vector(scale(julian_date)),
    time_of_day_std = as.vector(scale(time_of_day)),
    
    # Create factors
    survey_duration = factor(duration),
    survey_year = factor(survey_yea),
    gisid = factor(gisid)
  )
write.csv(method_C_data500, "1_Data/multinomial_model/methodCdata500.csv")
head(method_C_data500)
str(method_C_data500)

# Summary statistics for original variables
summary(method_C_data500$habitat_amount)
summary(method_C_data500$edge_density)

# Find and examine the outlier
outlier_row <- method_C_data500[method_C_data500$edge_density_std > 5, ]
print(outlier_row)

# Look at original edge density distribution
summary(method_C_data500$edge_density)
hist(method_C_data500$edge_density)

# Log transform then standardize (using log1p since you have zeros)
method_C_data500$edge_density_std <- scale(log1p(method_C_data500$edge_density))

# Check the new range
range(method_C_data500$edge_density_std)

# Remove outliers from dataset for now # Create new dataset without the outliers
method_C_data500 <- method_C_data500[method_C_data500$edge_density_std <= 5, ]

# Verify the removal
range(method_C_data500_clean$edge_density_std)

# Look at the distribution again
summary(method_C_data500_clean$edge_density)
hist(method_C_data500_clean$edge_density)
summary(method_C_data500_clean$habitat_amount_std)




# Step 1: Map the 'year' column in detection data to actual years using 'survey_year'
detection_data <- detection_data %>%
  mutate(actual_year = survey_year)

# Step 2: Rename the 'year' column in landscape metrics to 'actual_year' for alignment
landscape_metrics <- habitat_metrics %>%
  rename(actual_year = year)

# Step 3: Perform the join using 'surveyid' and 'actual_year'
matched_data <- detection_data %>%
  left_join(landscape_metrics, by = c("surveyid", "actual_year"))

# Step 4: Ensure landscape metrics are duplicated for each visit without altering occupancy values
final_data <- matched_data %>%
  select(-visit_num) %>% # Remove existing `visit_num` temporarily
  distinct() %>%         # Ensure no duplicated rows
  crossing(visit_num = 1:3) %>% # Add back `visit_num` with 3 duplicates
  left_join(select(detection_data, surveyid, visit_num, occupancy), 
            by = c("surveyid", "visit_num")) # Reattach the correct occupancy


# Check the number of rows and save the final dataset
print(paste("Number of rows in the final dataset:", nrow(final_data)))
write.csv(final_data, "final_dataset.csv", row.names = FALSE)
