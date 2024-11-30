# Load necessary libraries
library(dplyr)
library(DHARMa)
library(spdep)
library(lme4)
library(glmmTMB)
library(tidyverse)
library(brms)
library(bayesplot)
library(loo)
library(viridis)
library(sf)
library(patchwork)
library(mgcv)  # for GAM fitting
library(sjPlot)
library(glmm.hp)
library(MuMIn)

# Load datasets (update file paths as needed)
visit_matrix <- read.csv("1_Data/2yearvisitmatrix.csv")
visits_for_det_covs <- read.csv("1_Data/visitsfordetcovs.csv")
habitat_metrics <- read.csv("~/Chapter 2/Chapter 2 Analysis/Extract_patch_metrics/0_data/2_combined/landscape_metrics/500landscape_metrics_hyp3.csv")
head(visit_matrix)
head(visits_for_det_covs)
head(habitat_metrics)
str(visit_matrix)
str(visits_for_det_covs)
str(habitat_metrics)


# Prepare data

# Year 1 only models --------------------------------------------------------------
# Step 1: Filter only first-year visits (1_1, 1_2, 1_3) and create occupancy variable
visit_matrix <- visit_matrix %>%
  mutate(occupancy = pmax(`X1_1`, `X1_2`, `X1_3`, na.rm = TRUE)) %>%
  select(gisid, occupancy)  # Keep only necessary columns

# Step 2: Merge with visits_for_det_covs to get surveyid and first-year records only
visit_with_covs <- visit_matrix %>%
  inner_join(
    visits_for_det_covs %>%
      filter(year_id == 1) %>%  # Filter for the first year
      select(gisid, surveyid),
    by = "gisid"
  )

# Step 3: Merge with habitat_metrics to get habitat covariates and first-year survey data
full_data <- visit_with_covs %>%
  inner_join(
    habitat_metrics %>% 
      select(surveyid, habitat_amount, num_patches, edge_density, survey_yea), 
    by = "surveyid"
  )

# Remove missing values
full_data <- na.omit(full_data)

# Step 4: Split gisid into Easting and Northing for spatial analysis
full_data <- full_data %>%
  mutate(Easting = as.numeric(sub("_.*", "", gisid)),
         Northing = as.numeric(sub(".*_", "", gisid)))

# Check the filtered dataset
print(head(full_data))

# Fit the logistic regression model

# Step 1: Fit the GLM model with just habitat amount (logistic regression)
glm_model <- glm(occupancy ~ habitat_amount, data = full_data, family = binomial)

# Simulate residuals using DHARMa

# Step 1: Simulate residuals
residuals_sim <- simulateResiduals(fittedModel = glm_model)

# Step 2: Plot residual diagnostics
plot(residuals_sim)

# Step 3: Test for spatial autocorrelation 
# Create a spatial weights matrix
coords <- full_data %>% select(Easting, Northing)
neighbors <- dnearneigh(as.matrix(coords), d1 = 0, d2 = quantile(dist(coords), 0.1)) # 10th percentile distance
weights <- nb2listw(neighbors, style = "W")

# Perform Moran's I test for residuals
moran_test <- moran.test(residuals_sim$scaledResiduals, weights)
print(moran_test)


# Standardize habitat covariates
full_data$year <- as.factor(full_data$survey_yea)

full_data <- full_data %>%
  mutate(
    habitat_amount_std = scale(habitat_amount),
    num_patches_std = scale(num_patches),
    edge_density_std = scale(edge_density),
    year = droplevels(full_data$year)
  )
str(full_data$year)


table(full_data$occupancy, full_data$year)


# Fit Bayesian logistic regression with a random effect for year
bayesian_model <- brm(
  occupancy ~ habitat_amount_std + edge_density_std + Easting + (1 | year),
  data = full_data,
  family = bernoulli(),
  prior = c(
    prior(normal(0, 10), class = "b"),      # Priors for fixed effects
    prior(student_t(3, 0, 10), class = "sd") # Priors for random effect variance
  ),
  chains = 4, cores = 4, iter = 2000, control = list(adapt_delta = 0.95)
)

# Summarize the Bayesian model
summary(bayesian_model)


# try modelling the temporal trend with splines 
library(mgcv)
gam_model <- gam(
  occupancy ~ s(year, bs = "re") + habitat_amount_std + edge_density_std + Easting,
  family = binomial,
  data = full_data
)


# Fit a mixed-effects logistic regression model
mixed_model <- glm(
  occupancy ~ habitat_amount_std + edge_density_std + Easting,
  data = full_data,
  family = binomial
)

# Summary of the model
summary(mixed_model)

# Simulate residuals using DHARMa
residuals_sim <- simulateResiduals(fittedModel = mixed_model)

# Plot DHARMa residuals for diagnostics
plot(residuals_sim)

# Create a spatial weights matrix
coords <- full_data %>% select(Easting, Northing)
neighbors <- dnearneigh(as.matrix(coords), d1 = 0, d2 = quantile(dist(coords), 0.1))
weights <- nb2listw(neighbors, style = "W")


# Run GAM with residuals to add as autocovariate in my model 
residauto <- residuals(mixed_model)
modelgam <- gam(resid ~ s(Easting, Northing), family = gaussian, data = full_data)
plot(modelgam)
summary(modelgam)
outofcov <- fitted(modelgam)

# Moran's I Test 
moran_test <- moran.test(residauto, weights)
print(moran_test)
plot(full_data$Easting, residauto)
plot(full_data$Northing, residauto)


