# Load necessary libraries
library(dplyr)
library(DHARMa)
library(spdep)
library(lme4)
library(glmer)
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
library(performance)
library(ggplot2)
library(lattice)


# 500m Scale -----------------------------------------------------------------------------------------------------------------
# Load datasets 
methodA_500m <- read.csv("1_Data/500m/methodA_glmmdata.csv", header = T)
methodB_500m <- read.csv("1_Data/500m/methodB_glmmdata.csv", header = T)
methodC_500m <- read.csv("1_Data/500m/methodC_glmmdata.csv", header = T)

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



# A Method--------------------------------------------------------------------------------------------
# Required packages
library(lme4)
library(mgcv)
library(dplyr)
library(ggplot2)
library(MuMIn)
library(sjPlot)
library(partR2)
library(ncf)
library(spdep)
library(lmtest)

## 1. DATA PREPARATION -------------------------------------------------------
methodA_500m <- methodA_500m %>%
  mutate(
    season = paste(gisid, year_id, sep = "_"),
    survey_year = factor(survey_yea),
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    qpad = offset
  )

## 2. MODEL FITTING WITH MODEL-SPECIFIC SPATIAL TERMS ------------------------

# A. Full model with interaction
modA500int_base <- glmer(occupancy ~ offset(qpad) + prop_habitat_std * edge_density_std +
                           (1 | gisid / season) + (1 | survey_year),
                         data = methodA_500m, family = binomial)
residauto_int <- residuals(modA500int_base)
modelgam_int <- gam(residauto_int ~ s(Easting, Northing), family = gaussian, data = methodA_500m)
methodA_500m$spatial_autocov_int <- scale(fitted(modelgam_int))[,1]
modA500int <- glmer(occupancy ~ offset(qpad) + prop_habitat_std * edge_density_std + 
                      spatial_autocov_int + (1 | gisid / season) + (1 | survey_year), 
                    data = methodA_500m, family = binomial)
summary(modA500int)
summary(modA500_A)
summary(modA500_F)
summary(modA500_AF)
model_selectionA500 <- model.sel(modA500int, modA500_A, modA500_AF, modA500_F)
print(model_selectionA500)
print("R2 values:")
r.squaredGLMM(modA500int)
r.squaredGLMM(modA500_A)
r.squaredGLMM(modA500_F)
r.squaredGLMM(modA500_AF)

# B. Habitat amount only
modA500_A <- glmer(occupancy ~ offset(qpad) + prop_habitat_std + spatial_autocov_int +
                     (1 | gisid / season) + (1 | survey_year), 
                   data = methodA_500m, family = binomial)

# C. Fragmentation only
modA500_F <- glmer(occupancy ~ offset(qpad) + edge_density_std + spatial_autocov_int +
                     (1 | gisid / season) + (1 | survey_year), 
                   data = methodA_500m, family = binomial)

# D. Amount + fragmentation (no interaction)
modA500_AF <- glmer(occupancy ~ offset(qpad) + prop_habitat_std + edge_density_std + 
                      spatial_autocov_int + (1 | gisid / season) + (1 | survey_year), 
                    data = methodA_500m, family = binomial)

## 3. MODEL SELECTION AND COMPARISON ----------------------------------------

# AIC model selection
model_selectionA500 <- model.sel(modA500int, modA500_A, modA500_AF, modA500_F)
print(model_selectionA500)

# Likelihood ratio test
lrtest(modA500int, modA500_AF)

## 4. MODEL DIAGNOSTICS ----------------------------------------------------

# Spatial autocorrelation checks for each model
coord <- cbind(methodA_500m$Easting, methodA_500m$Northing)

# Full model
resid_int <- residuals(modA500int)
spline_corr_int <- ncf::spline.correlog(methodA_500m$Easting, methodA_500m$Northing,
                                        resid_int, resamp = 10, type = "boot") # change to 10 from 100 to make it faster 

# Amount only
resid_A <- residuals(modA500_A)
spline_corr_A <- ncf::spline.correlog(methodA_500m$Easting, methodA_500m$Northing,
                                      resid_A, resamp = 10, type = "boot")

# Fragmentation only
resid_F <- residuals(modA500_F)
spline_corr_F <- ncf::spline.correlog(methodA_500m$Easting, methodA_500m$Northing,
                                      resid_F, resamp = 10, type = "boot")

# Amount + Fragmentation
resid_AF <- residuals(modA500_AF)
spline_corr_AF <- ncf::spline.correlog(methodA_500m$Easting, methodA_500m$Northing,
                                       resid_AF, resamp = 10, type = "boot")

# Plot spatial correlograms
par(mfrow = c(2,2))
plot(spline_corr_int, ylim = c(-0.2, 0.2), main = "Interaction Model")
plot(spline_corr_A, ylim = c(-0.2, 0.2), main = "Amount Only")
plot(spline_corr_F, ylim = c(-0.2, 0.2), main = "Fragmentation Only")
plot(spline_corr_AF, ylim = c(-0.2, 0.2), main = "Amount + Fragmentation")
par(mfrow = c(1,1))

# Temporal autocorrelation checks
par(mfrow = c(2,2))
pacf(residuals(modA500int), main = "Interaction Model")
pacf(residuals(modA500_A), main = "Amount Only")
pacf(residuals(modA500_F), main = "Fragmentation Only")
pacf(residuals(modA500_AF), main = "Amount + Fragmentation")
par(mfrow = c(1,1))


## 6. VARIANCE PARTITIONING ------------------------------------------------

# Set up parallel processing
plan(multisession, workers = 18)

# Partition variance for full model
r2part <- partR2(
  modA500int,
  partvars = c("prop_habitat_std", "edge_density_std", "prop_habitat_std:edge_density_std"), # have to manually plot these results to only include the four model combinations 
  data = methodA_500m,
  R2_type = "marginal",
  nboot = 99,
  CI = 0.95,
  parallel = TRUE
)

# Plot variance partitioning results
forestplot(r2part, type = c("R2"))

## 7. VISUALIZATION -------------------------------------------------------

# Effect plots for best model (assuming interaction model)
plot_model(modA500int, type = "int", mdrt.values = "quart")
plot_model(modA500int, type = "pred")

# Coefficient plots for model comparison
plot_models(modA500int, modA500_A, modA500_AF, modA500_F, transform = NULL)

# Calculate R2 for all models
print("R2 values:")
r.squaredGLMM(modA500int)
r.squaredGLMM(modA500_A)
r.squaredGLMM(modA500_F)
r.squaredGLMM(modA500_AF)


# B Method--------------------------------------------------------------------------------------------

## 1. DATA PREPARATION -------------------------------------------------------
methodB_500m <- methodB_500m %>%
  mutate(
    season = paste(gisid, year_id, sep = "_"),
    survey_year = factor(survey_yea),
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    qpad = offset
  )

## 2. MODEL FITTING WITH MODEL-SPECIFIC SPATIAL TERMS ------------------------

# A. Full model with interaction
modB500int_base <- glmer(occupancy ~ offset(qpad) + prop_habitat_std * edge_density_std +
                           (1 | gisid / season) + (1 | survey_year),
                         data = methodB_500m, family = binomial)
residauto_int <- residuals(modB500int_base)
modelgam_int <- gam(residauto_int ~ s(Easting, Northing), family = gaussian, data = methodB_500m)
methodB_500m$spatial_autocov_int <- scale(fitted(modelgam_int))[,1]
modB500int <- glmer(occupancy ~ offset(qpad) + prop_habitat_std * edge_density_std + 
                      spatial_autocov_int + (1 | gisid / season) + (1 | survey_year), 
                    data = methodB_500m, family = binomial)

# B. Habitat amount only
modB500_A <- glmer(occupancy ~ offset(qpad) + prop_habitat_std + spatial_autocov_int +
                     (1 | gisid / season) + (1 | survey_year), 
                   data = methodB_500m, family = binomial)

# C. Fragmentation only
modB500_F <- glmer(occupancy ~ offset(qpad) + edge_density_std + spatial_autocov_int +
                     (1 | gisid / season) + (1 | survey_year), 
                   data = methodB_500m, family = binomial)

# D. Amount + fragmentation (no interaction)
modB500_AF <- glmer(occupancy ~ offset(qpad) + prop_habitat_std + edge_density_std + 
                      spatial_autocov_int + (1 | gisid / season) + (1 | survey_year), 
                    data = methodB_500m, family = binomial)

## 3. MODEL SELECTION AND COMPARISON ----------------------------------------

# AIC model selection
model_selectionB500 <- model.sel(modB500int, modB500_A, modB500_AF, modB500_F)
print(model_selectionB500)

# Likelihood ratio test
lrtest(modB500int, modB500_AF)

## 4. MODEL DIAGNOSTICS ----------------------------------------------------

# Spatial autocorrelation checks for each model
coord <- cbind(methodB_500m$Easting, methodB_500m$Northing)

# Full model
resid_int <- residuals(modB500int)
spline_corr_int <- ncf::spline.correlog(methodB_500m$Easting, methodB_500m$Northing,
                                        resid_int, resamp = 10, type = "boot") # change to 10 from 100 to make it faster 

# Amount only
resid_A <- residuals(modB500_A)
spline_corr_A <- ncf::spline.correlog(methodB_500m$Easting, methodB_500m$Northing,
                                      resid_A, resamp = 10, type = "boot")

# Fragmentation only
resid_F <- residuals(modB500_F)
spline_corr_F <- ncf::spline.correlog(methodB_500m$Easting, methodB_500m$Northing,
                                      resid_F, resamp = 10, type = "boot")

# Amount + Fragmentation
resid_AF <- residuals(modB500_AF)
spline_corr_AF <- ncf::spline.correlog(methodB_500m$Easting, methodB_500m$Northing,
                                       resid_AF, resamp = 10, type = "boot")

# Plot spatial correlograms
par(mfrow = c(2,2))
plot(spline_corr_int, ylim = c(-0.2, 0.2), main = "Interaction Model")
plot(spline_corr_A, ylim = c(-0.2, 0.2), main = "Amount Only")
plot(spline_corr_F, ylim = c(-0.2, 0.2), main = "Fragmentation Only")
plot(spline_corr_AF, ylim = c(-0.2, 0.2), main = "Amount + Fragmentation")
par(mfrow = c(1,1))

# Temporal autocorrelation checks
par(mfrow = c(2,2))
pacf(residuals(modB500int), main = "Interaction Model")
pacf(residuals(modB500_A), main = "Amount Only")
pacf(residuals(modB500_F), main = "Fragmentation Only")
pacf(residuals(modB500_AF), main = "Amount + Fragmentation")
par(mfrow = c(1,1))


## 6. VARIANCE PARTITIONING ------------------------------------------------

# Set up parallel processing
plan(multisession, workers = 18)

# Partition variance for full model
r2part <- partR2(
  modB500int,
  partvars = c("prop_habitat_std", "edge_density_std", "prop_habitat_std:edge_density_std"), # have to manually plot these results to only include the four model combinations 
  data = methodB_500m,
  R2_type = "marginal",
  nboot = 99,
  CI = 0.95,
  parallel = TRUE
)

# Plot variance partitioning results
forestplot(r2part, type = c("R2"))

## 7. VISUALIZATION -------------------------------------------------------

# Effect plots for best model (assuming interaction model)
plot_model(modB500int, type = "int", mdrt.values = "quart")
plot_model(modB500int, type = "pred")

# Coefficient plots for model comparison
plot_models(modB500int, modB500_A, modB500_AF, modB500_F, transform = NULL)

# Calculate R2 for all models
print("R2 values:")
r.squaredGLMM(modB500int)
r.squaredGLMM(modB500_A)
r.squaredGLMM(modB500_F)
r.squaredGLMM(modB500_AF)



# C Method------------------------------------------------------------------------------------------------- 

## 1. DATA PREPARATION -------------------------------------------------------
methodC_500m <- methodC_500m %>%
  mutate(
    season = paste(gisid, year_id, sep = "_"),
    survey_year = factor(survey_yea),
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    qpad = offset
  )

## 2. MODEL FITTING WITH MODEL-SPECIFIC SPATIAL TERMS ------------------------

# A. Full model with interaction
modC500int_base <- glmer(occupancy ~ offset(qpad) + prop_habitat_std * edge_density_std +
                           (1 | gisid / season) + (1 | survey_year),
                         data = methodC_500m, family = binomial)
residauto_int <- residuals(modC500int_base)
modelgam_int <- gam(residauto_int ~ s(Easting, Northing), family = gaussian, data = methodC_500m)
methodC_500m$spatial_autocov_int <- scale(fitted(modelgam_int))[,1]
modC500int <- glmer(occupancy ~ offset(qpad) + prop_habitat_std * edge_density_std + 
                      spatial_autocov_int + (1 | gisid / season) + (1 | survey_year), 
                    data = methodC_500m, family = binomial)

# B. Habitat amount only
modC500_A <- glmer(occupancy ~ offset(qpad) + prop_habitat_std + spatial_autocov_int +
                     (1 | gisid / season) + (1 | survey_year), 
                   data = methodC_500m, family = binomial)

# C. Fragmentation only
modC500_F <- glmer(occupancy ~ offset(qpad) + edge_density_std + spatial_autocov_int +
                     (1 | gisid / season) + (1 | survey_year), 
                   data = methodC_500m, family = binomial)

# D. Amount + fragmentation (no interaction)
modC500_AF <- glmer(occupancy ~ offset(qpad) + prop_habitat_std + edge_density_std + 
                      spatial_autocov_int + (1 | gisid / season) + (1 | survey_year), 
                    data = methodC_500m, family = binomial)

## 3. MODEL SELECTION AND COMPARISON ----------------------------------------

# AIC model selection
model_selectionC500 <- model.sel(modC500int, modC500_A, modC500_AF, modC500_F)
print(model_selectionC500)

# Likelihood ratio test
lrtest(modC500int, modC500_AF)

## 4. MODEL DIAGNOSTICS ----------------------------------------------------

# Spatial autocorrelation checks for each model
coord <- cbind(methodC_500m$Easting, methodC_500m$Northing)

# Full model
resid_int <- residuals(modC500int)
spline_corr_int <- ncf::spline.correlog(methodC_500m$Easting, methodC_500m$Northing,
                                        resid_int, resamp = 10, type = "boot") # change to 10 from 100 to make it faster 

# Amount only
resid_A <- residuals(modC500_A)
spline_corr_A <- ncf::spline.correlog(methodC_500m$Easting, methodC_500m$Northing,
                                      resid_A, resamp = 10, type = "boot")

# Fragmentation only
resid_F <- residuals(modC500_F)
spline_corr_F <- ncf::spline.correlog(methodC_500m$Easting, methodC_500m$Northing,
                                      resid_F, resamp = 10, type = "boot")

# Amount + Fragmentation
resid_AF <- residuals(modC500_AF)
spline_corr_AF <- ncf::spline.correlog(methodC_500m$Easting, methodC_500m$Northing,
                                       resid_AF, resamp = 10, type = "boot")

# Plot spatial correlograms
par(mfrow = c(2,2))
plot(spline_corr_int, ylim = c(-0.2, 0.2), main = "Interaction Model")
plot(spline_corr_A, ylim = c(-0.2, 0.2), main = "Amount Only")
plot(spline_corr_F, ylim = c(-0.2, 0.2), main = "Fragmentation Only")
plot(spline_corr_AF, ylim = c(-0.2, 0.2), main = "Amount + Fragmentation")
par(mfrow = c(1,1))

# Temporal autocorrelation checks
par(mfrow = c(2,2))
pacf(residuals(modC500int), main = "Interaction Model")
pacf(residuals(modC500_A), main = "Amount Only")
pacf(residuals(modC500_F), main = "Fragmentation Only")
pacf(residuals(modC500_AF), main = "Amount + Fragmentation")
par(mfrow = c(1,1))


## 6. VARIANCE PARTITIONING ------------------------------------------------

# Set up parallel processing
plan(multisession, workers = 18)

# Partition variance for full model
r2part <- partR2(
  modC500int,
  partvars = c("prop_habitat_std", "edge_density_std", "prop_habitat_std:edge_density_std"), # have to manually plot these results to only include the four model combinations 
  data = methodC_500m,
  R2_type = "marginal",
  nboot = 99,
  CI = 0.95,
  parallel = TRUE
)

# Plot variance partitioning results
forestplot(r2part, type = c("R2"))

## 7. VISUALIZATION -------------------------------------------------------

# Effect plots for best model (assuming interaction model)
plot_model(modC500int, type = "int", mdrt.values = "quart")
plot_model(modC500int, type = "pred")

# Coefficient plots for model comparison
plot_models(modC500int, modC500_A, modC500_AF, modC500_F, transform = NULL)

# Calculate R2 for all models
print("R2 values:")
r.squaredGLMM(modC500int)
r.squaredGLMM(modC500_A)
r.squaredGLMM(modC500_F)
r.squaredGLMM(modC500_AF)



# 1000m Scale -----------------------------------------------------------------------------------------------------------------
# Load datasets 
methodA_1000m <- read.csv("1_Data/1000m/methodA_glmmdata.csv", header = T)
methodB_1000m <- read.csv("1_Data/1000m/methodB_glmmdata.csv", header = T)
methodC_1000m <- read.csv("1_Data/1000m/methodC_glmmdata.csv", header = T)

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



# A Method--------------------------------------------------------------------------------------------
# Required packages
library(lme4)
library(mgcv)
library(dplyr)
library(ggplot2)
library(MuMIn)
library(sjPlot)
library(partR2)
library(ncf)
library(spdep)
library(lmtest)

## 1. DATA PREPARATION -------------------------------------------------------
methodA_1000m <- methodA_1000m %>%
  mutate(
    season = paste(gisid, year_id, sep = "_"),
    survey_year = factor(survey_yea),
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    qpad = offset
  )

## 2. MODEL FITTING WITH MODEL-SPECIFIC SPATIAL TERMS ------------------------

# A. Full model with interaction
modA1000int_base <- glmer(occupancy ~ offset(qpad) + prop_habitat_std * edge_density_std +
                           (1 | gisid / season) + (1 | survey_year),
                         data = methodA_1000m, family = binomial)
residauto_int <- residuals(modA1000int_base)
modelgam_int <- gam(residauto_int ~ s(Easting, Northing), family = gaussian, data = methodA_1000m)
methodA_1000m$spatial_autocov_int <- scale(fitted(modelgam_int))[,1]
modA1000int <- glmer(occupancy ~ offset(qpad) + prop_habitat_std * edge_density_std + 
                      spatial_autocov_int + (1 | gisid / season) + (1 | survey_year), 
                    data = methodA_1000m, family = binomial)
summary(modA1000int)
# B. Habitat amount only
modA1000_A <- glmer(occupancy ~ offset(qpad) + prop_habitat_std + spatial_autocov_int +
                     (1 | gisid / season) + (1 | survey_year), 
                   data = methodA_1000m, family = binomial)
summary(modA1000_A)
# C. Fragmentation only
modA1000_F <- glmer(occupancy ~ offset(qpad) + edge_density_std + spatial_autocov_int +
                     (1 | gisid / season) + (1 | survey_year), 
                   data = methodA_1000m, family = binomial)
summary(modA1000_F)
# D. Amount + fragmentation (no interaction)
modA1000_AF <- glmer(occupancy ~ offset(qpad) + prop_habitat_std + edge_density_std + 
                      spatial_autocov_int + (1 | gisid / season) + (1 | survey_year), 
                    data = methodA_1000m, family = binomial)
summary(modA1000_AF)

## 3. MODEL SELECTION AND COMPARISON ----------------------------------------

# AIC model selection
model_selectionA1000 <- model.sel(modA1000int, modA1000_A, modA1000_AF, modA1000_F)
print(model_selectionA1000)

# Likelihood ratio test
lrtest(modA1000int, modA1000_AF)

## 4. MODEL DIAGNOSTICS ----------------------------------------------------

# Spatial autocorrelation checks for each model
coord <- cbind(methodA_1000m$Easting, methodA_1000m$Northing)

# Full model
resid_int <- residuals(modA1000int)
spline_corr_int <- ncf::spline.correlog(methodA_1000m$Easting, methodA_1000m$Northing,
                                        resid_int, resamp = 10, type = "boot") # change to 10 from 100 to make it faster 

# Amount only
resid_A <- residuals(modA1000_A)
spline_corr_A <- ncf::spline.correlog(methodA_1000m$Easting, methodA_1000m$Northing,
                                      resid_A, resamp = 10, type = "boot")

# Fragmentation only
resid_F <- residuals(modA1000_F)
spline_corr_F <- ncf::spline.correlog(methodA_1000m$Easting, methodA_1000m$Northing,
                                      resid_F, resamp = 10, type = "boot")

# Amount + Fragmentation
resid_AF <- residuals(modA1000_AF)
spline_corr_AF <- ncf::spline.correlog(methodA_1000m$Easting, methodA_1000m$Northing,
                                       resid_AF, resamp = 10, type = "boot")

# Plot spatial correlograms
par(mfrow = c(2,2))
plot(spline_corr_int, ylim = c(-0.2, 0.2), main = "Interaction Model")
plot(spline_corr_A, ylim = c(-0.2, 0.2), main = "Amount Only")
plot(spline_corr_F, ylim = c(-0.2, 0.2), main = "Fragmentation Only")
plot(spline_corr_AF, ylim = c(-0.2, 0.2), main = "Amount + Fragmentation")
par(mfrow = c(1,1))

# Temporal autocorrelation checks
par(mfrow = c(2,2))
pacf(residuals(modA1000int), main = "Interaction Model")
pacf(residuals(modA1000_A), main = "Amount Only")
pacf(residuals(modA1000_F), main = "Fragmentation Only")
pacf(residuals(modA1000_AF), main = "Amount + Fragmentation")
par(mfrow = c(1,1))


## 6. VARIANCE PARTITIONING ------------------------------------------------

# Set up parallel processing
plan(multisession, workers = 18)

# Partition variance for full model
r2part <- partR2(
  modA1000int,
  partvars = c("prop_habitat_std", "edge_density_std", "prop_habitat_std:edge_density_std"), # have to manually plot these results to only include the four model combinations 
  data = methodA_1000m,
  R2_type = "marginal",
  nboot = 99,
  CI = 0.95,
  parallel = TRUE
)

# Plot variance partitioning results
forestplot(r2part, type = c("R2"))

## 7. VISUALIZATION -------------------------------------------------------

# Effect plots for best model (assuming interaction model)
plot_model(modA1000int, type = "int", mdrt.values = "quart")
plot_model(modA1000int, type = "pred")

# Coefficient plots for model comparison
plot_models(modA1000int, modA1000_A, modA1000_AF, modA1000_F, transform = NULL)

# Calculate R2 for all models
print("R2 values:")
r.squaredGLMM(modA1000int)
r.squaredGLMM(modA1000_A)
r.squaredGLMM(modA1000_F)
r.squaredGLMM(modA1000_AF)

# B Method--------------------------------------------------------------------------------------------

## 1. DATA PREPARATION -------------------------------------------------------
methodB_1000m <- methodB_1000m %>%
  mutate(
    season = paste(gisid, year_id, sep = "_"),
    survey_year = factor(survey_yea),
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    qpad = offset
  )

## 2. MODEL FITTING WITH MODEL-SPECIFIC SPATIAL TERMS ------------------------

# A. Full model with interaction
modB1000int_base <- glmer(occupancy ~ offset(qpad) + prop_habitat_std * edge_density_std +
                           (1 | gisid / season) + (1 | survey_year),
                         data = methodB_1000m, family = binomial)
residauto_int <- residuals(modB1000int_base)
modelgam_int <- gam(residauto_int ~ s(Easting, Northing), family = gaussian, data = methodB_1000m)
methodB_1000m$spatial_autocov_int <- scale(fitted(modelgam_int))[,1]
modB1000int <- glmer(occupancy ~ offset(qpad) + prop_habitat_std * edge_density_std + 
                      spatial_autocov_int + (1 | gisid / season) + (1 | survey_year), 
                    data = methodB_1000m, family = binomial)
summary(modB1000int)
# B. Habitat amount only
modB1000_A <- glmer(occupancy ~ offset(qpad) + prop_habitat_std + spatial_autocov_int +
                     (1 | gisid / season) + (1 | survey_year), 
                   data = methodB_1000m, family = binomial)

# C. Fragmentation only
modB1000_F <- glmer(occupancy ~ offset(qpad) + edge_density_std + spatial_autocov_int +
                     (1 | gisid / season) + (1 | survey_year), 
                   data = methodB_1000m, family = binomial)

# D. Amount + fragmentation (no interaction)
modB1000_AF <- glmer(occupancy ~ offset(qpad) + prop_habitat_std + edge_density_std + 
                      spatial_autocov_int + (1 | gisid / season) + (1 | survey_year), 
                    data = methodB_1000m, family = binomial)

## 3. MODEL SELECTION AND COMPARISON ----------------------------------------

# AIC model selection
model_selectionB1000 <- model.sel(modB1000int, modB1000_A, modB1000_AF, modB1000_F)
print(model_selectionB1000)

# Likelihood ratio test
lrtest(modB1000int, modB1000_AF)

## 4. MODEL DIAGNOSTICS ----------------------------------------------------

# Spatial autocorrelation checks for each model
coord <- cbind(methodB_1000m$Easting, methodB_1000m$Northing)

# Full model
resid_int <- residuals(modB1000int)
spline_corr_int <- ncf::spline.correlog(methodB_1000m$Easting, methodB_1000m$Northing,
                                        resid_int, resamp = 10, type = "boot") # change to 10 from 100 to make it faster 

# Amount only
resid_A <- residuals(modB1000_A)
spline_corr_A <- ncf::spline.correlog(methodB_1000m$Easting, methodB_1000m$Northing,
                                      resid_A, resamp = 10, type = "boot")

# Fragmentation only
resid_F <- residuals(modB1000_F)
spline_corr_F <- ncf::spline.correlog(methodB_1000m$Easting, methodB_1000m$Northing,
                                      resid_F, resamp = 10, type = "boot")

# Amount + Fragmentation
resid_AF <- residuals(modB1000_AF)
spline_corr_AF <- ncf::spline.correlog(methodB_1000m$Easting, methodB_1000m$Northing,
                                       resid_AF, resamp = 10, type = "boot")

# Plot spatial correlograms
par(mfrow = c(2,2))
plot(spline_corr_int, ylim = c(-0.2, 0.2), main = "Interaction Model")
plot(spline_corr_A, ylim = c(-0.2, 0.2), main = "Amount Only")
plot(spline_corr_F, ylim = c(-0.2, 0.2), main = "Fragmentation Only")
plot(spline_corr_AF, ylim = c(-0.2, 0.2), main = "Amount + Fragmentation")
par(mfrow = c(1,1))

# Temporal autocorrelation checks
par(mfrow = c(2,2))
pacf(residuals(modB1000int), main = "Interaction Model")
pacf(residuals(modB1000_A), main = "Amount Only")
pacf(residuals(modB1000_F), main = "Fragmentation Only")
pacf(residuals(modB1000_AF), main = "Amount + Fragmentation")
par(mfrow = c(1,1))


## 6. VARIANCE PARTITIONING ------------------------------------------------

# Set up parallel processing
plan(multisession, workers = 18)

# Partition variance for full model
r2part <- partR2(
  modB1000int,
  partvars = c("prop_habitat_std", "edge_density_std", "prop_habitat_std:edge_density_std"), # have to manually plot these results to only include the four model combinations 
  data = methodB_1000m,
  R2_type = "marginal",
  nboot = 99,
  CI = 0.95,
  parallel = TRUE
)

# Plot variance partitioning results
forestplot(r2part, type = c("R2"))

## 7. VISUALIZATION -------------------------------------------------------

# Effect plots for best model (assuming interaction model)
plot_model(modB1000int, type = "int", mdrt.values = "quart")
plot_model(modB1000int, type = "pred")

# Coefficient plots for model comparison
plot_models(modB1000int, modB1000_A, modB1000_AF, modB1000_F, transform = NULL)

# Calculate R2 for all models
print("R2 values:")
r.squaredGLMM(modB1000int)
r.squaredGLMM(modB1000_A)
r.squaredGLMM(modB1000_F)
r.squaredGLMM(modB1000_AF)



# C Method------------------------------------------------------------------------------------------------- 

## 1. DATA PREPARATION -------------------------------------------------------
methodC_1000m <- methodC_1000m %>%
  mutate(
    season = paste(gisid, year_id, sep = "_"),
    survey_year = factor(survey_yea),
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    qpad = offset
  )

## 2. MODEL FITTING WITH MODEL-SPECIFIC SPATIAL TERMS ------------------------

# A. Full model with interaction
modC1000int_base <- glmer(occupancy ~ offset(qpad) + prop_habitat_std * edge_density_std +
                           (1 | gisid / season) + (1 | survey_year),
                         data = methodC_1000m, family = binomial)
residauto_int <- residuals(modC1000int_base)
modelgam_int <- gam(residauto_int ~ s(Easting, Northing), family = gaussian, data = methodC_1000m)
methodC_1000m$spatial_autocov_int <- scale(fitted(modelgam_int))[,1]
modC1000int <- glmer(occupancy ~ offset(qpad) + prop_habitat_std * edge_density_std + 
                      spatial_autocov_int + (1 | gisid / season) + (1 | survey_year), 
                    data = methodC_1000m, family = binomial)
summary(modC1000int)
# B. Habitat amount only
modC1000_A <- glmer(occupancy ~ offset(qpad) + prop_habitat_std + spatial_autocov_int +
                     (1 | gisid / season) + (1 | survey_year), 
                   data = methodC_1000m, family = binomial)

# C. Fragmentation only
modC1000_F <- glmer(occupancy ~ offset(qpad) + edge_density_std + spatial_autocov_int +
                     (1 | gisid / season) + (1 | survey_year), 
                   data = methodC_1000m, family = binomial)

# D. Amount + fragmentation (no interaction)
modC1000_AF <- glmer(occupancy ~ offset(qpad) + prop_habitat_std + edge_density_std + 
                      spatial_autocov_int + (1 | gisid / season) + (1 | survey_year), 
                    data = methodC_1000m, family = binomial)
summary(modC1000_AF)
## 3. MODEL SELECTION AND COMPARISON ----------------------------------------

# AIC model selection
model_selectionC1000 <- model.sel(modC1000int, modC1000_A, modC1000_AF, modC1000_F)
print(model_selectionC1000)

# Likelihood ratio test
lrtest(modC1000int, modC1000_AF)

## 4. MODEL DIAGNOSTICS ----------------------------------------------------

# Spatial autocorrelation checks for each model
coord <- cbind(methodC_1000m$Easting, methodC_1000m$Northing)

# Full model
resid_int <- residuals(modC1000int)
spline_corr_int <- ncf::spline.correlog(methodC_1000m$Easting, methodC_1000m$Northing,
                                        resid_int, resamp = 10, type = "boot") # change to 10 from 100 to make it faster 

# Amount only
resid_A <- residuals(modC1000_A)
spline_corr_A <- ncf::spline.correlog(methodC_1000m$Easting, methodC_1000m$Northing,
                                      resid_A, resamp = 10, type = "boot")

# Fragmentation only
resid_F <- residuals(modC1000_F)
spline_corr_F <- ncf::spline.correlog(methodC_1000m$Easting, methodC_1000m$Northing,
                                      resid_F, resamp = 10, type = "boot")

# Amount + Fragmentation
resid_AF <- residuals(modC1000_AF)
spline_corr_AF <- ncf::spline.correlog(methodC_1000m$Easting, methodC_1000m$Northing,
                                       resid_AF, resamp = 10, type = "boot")

# Plot spatial correlograms
par(mfrow = c(2,2))
plot(spline_corr_int, ylim = c(-0.2, 0.2), main = "Interaction Model")
plot(spline_corr_A, ylim = c(-0.2, 0.2), main = "Amount Only")
plot(spline_corr_F, ylim = c(-0.2, 0.2), main = "Fragmentation Only")
plot(spline_corr_AF, ylim = c(-0.2, 0.2), main = "Amount + Fragmentation")
par(mfrow = c(1,1))

# Temporal autocorrelation checks
par(mfrow = c(2,2))
pacf(residuals(modC1000int), main = "Interaction Model")
pacf(residuals(modC1000_A), main = "Amount Only")
pacf(residuals(modC1000_F), main = "Fragmentation Only")
pacf(residuals(modC1000_AF), main = "Amount + Fragmentation")
par(mfrow = c(1,1))


## 6. VARIANCE PARTITIONING ------------------------------------------------

# Set up parallel processing
plan(multisession, workers = 18)

# Partition variance for full model
r2part <- partR2(
  modC1000int,
  partvars = c("prop_habitat_std", "edge_density_std", "prop_habitat_std:edge_density_std"), # have to manually plot these results to only include the four model combinations 
  data = methodC_1000m,
  R2_type = "marginal",
  nboot = 99,
  CI = 0.95,
  parallel = TRUE
)

# Plot variance partitioning results
forestplot(r2part, type = c("R2"))

## 7. VISUALIZATION -------------------------------------------------------

# Effect plots for best model (assuming interaction model)
plot_model(modC1000int, type = "int", mdrt.values = "quart")
plot_model(modC1000int, type = "pred")

# Coefficient plots for model comparison
plot_models(modC1000int, modC1000_A, modC1000_AF, modC1000_F, transform = NULL)

# Calculate R2 for all models
print("R2 values:")
r.squaredGLMM(modC1000int)
r.squaredGLMM(modC1000_A)
r.squaredGLMM(modC1000_F)
r.squaredGLMM(modC1000_AF)




# 150m Scale --------------------------------------------------------------------------------------------
# A Method--------------------------------------------------------------------------------------------

## 1. DATA PREPARATION -------------------------------------------------------
methodA_150m <- methodA_150m %>%
  mutate(
    season = paste(gisid, year_id, sep = "_"),
    survey_year = factor(survey_yea),
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    qpad = offset
  )

## 2. MODEL FITTING WITH MODEL-SPECIFIC SPATIAL TERMS ------------------------

# A. Full model with interaction
modA150int_base <- glmer(occupancy ~ offset(qpad) + prop_habitat_std * edge_density_std +
                            (1 | gisid / season) + (1 | survey_year),
                          data = methodA_150m, family = binomial)
residauto_int <- residuals(modA150int_base)
modelgam_int <- gam(residauto_int ~ s(Easting, Northing), family = gaussian, data = methodA_150m)
methodA_150m$spatial_autocov_int <- scale(fitted(modelgam_int))[,1]
modA150int <- glmer(occupancy ~ offset(qpad) + prop_habitat_std * edge_density_std + 
                       spatial_autocov_int + (1 | gisid / season) + (1 | survey_year), 
                     data = methodA_150m, family = binomial)
summary(modA150int)
# B. Habitat amount only
modA150_A <- glmer(occupancy ~ offset(qpad) + prop_habitat_std + spatial_autocov_int +
                      (1 | gisid / season) + (1 | survey_year), 
                    data = methodA_150m, family = binomial)
summary(modA150_A)
# C. Fragmentation only
modA150_F <- glmer(occupancy ~ offset(qpad) + edge_density_std + spatial_autocov_int +
                      (1 | gisid / season) + (1 | survey_year), 
                    data = methodA_150m, family = binomial)
summary(modA150_F)
# D. Amount + fragmentation (no interaction)
modA150_AF <- glmer(occupancy ~ offset(qpad) + prop_habitat_std + edge_density_std + 
                       spatial_autocov_int + (1 | gisid / season) + (1 | survey_year), 
                     data = methodA_150m, family = binomial)
summary(modA150_AF)

## 3. MODEL SELECTION AND COMPARISON ----------------------------------------

# AIC model selection
model_selectionA150 <- model.sel(modA150int, modA150_A, modA150_AF, modA150_F)
print(model_selectionA150)

# Likelihood ratio test
lrtest(modA150int, modA150_AF)

## 4. MODEL DIAGNOSTICS ----------------------------------------------------

# Spatial autocorrelation checks for each model
coord <- cbind(methodA_150m$Easting, methodA_150m$Northing)

# Full model
resid_int <- residuals(modA150int)
spline_corr_int <- ncf::spline.correlog(methodA_150m$Easting, methodA_150m$Northing,
                                        resid_int, resamp = 10, type = "boot") # change to 10 from 100 to make it faster 

# Amount only
resid_A <- residuals(modA150_A)
spline_corr_A <- ncf::spline.correlog(methodA_150m$Easting, methodA_150m$Northing,
                                      resid_A, resamp = 10, type = "boot")

# Fragmentation only
resid_F <- residuals(modA150_F)
spline_corr_F <- ncf::spline.correlog(methodA_150m$Easting, methodA_150m$Northing,
                                      resid_F, resamp = 10, type = "boot")

# Amount + Fragmentation
resid_AF <- residuals(modA150_AF)
spline_corr_AF <- ncf::spline.correlog(methodA_150m$Easting, methodA_150m$Northing,
                                       resid_AF, resamp = 10, type = "boot")

# Plot spatial correlograms
par(mfrow = c(2,2))
plot(spline_corr_int, ylim = c(-0.2, 0.2), main = "Interaction Model")
plot(spline_corr_A, ylim = c(-0.2, 0.2), main = "Amount Only")
plot(spline_corr_F, ylim = c(-0.2, 0.2), main = "Fragmentation Only")
plot(spline_corr_AF, ylim = c(-0.2, 0.2), main = "Amount + Fragmentation")
par(mfrow = c(1,1))

# Temporal autocorrelation checks
par(mfrow = c(2,2))
pacf(residuals(modA150int), main = "Interaction Model")
pacf(residuals(modA150_A), main = "Amount Only")
pacf(residuals(modA150_F), main = "Fragmentation Only")
pacf(residuals(modA150_AF), main = "Amount + Fragmentation")
par(mfrow = c(1,1))


## 6. VARIANCE PARTITIONING ------------------------------------------------

# Set up parallel processing
plan(multisession, workers = 18)

# Partition variance for full model
r2part <- partR2(
  modA150int,
  partvars = c("prop_habitat_std", "edge_density_std", "prop_habitat_std:edge_density_std"), # have to manually plot these results to only include the four model combinations 
  data = methodA_150m,
  R2_type = "marginal",
  nboot = 99,
  CI = 0.95,
  parallel = TRUE
)

# Plot variance partitioning results
forestplot(r2part, type = c("R2"))

## 7. VISUALIZATION -------------------------------------------------------

# Effect plots for best model (assuming interaction model)
plot_model(modA150int, type = "int", mdrt.values = "quart")
plot_model(modA150int, type = "pred")

# Coefficient plots for model comparison
plot_models(modA150int, modA150_A, modA150_AF, modA150_F, transform = NULL)

# Calculate R2 for all models
print("R2 values:")
r.squaredGLMM(modA150int)
r.squaredGLMM(modA150_A)
r.squaredGLMM(modA150_F)
r.squaredGLMM(modA150_AF)


# B Method--------------------------------------------------------------------------------------------

## 1. DATA PREPARATION -------------------------------------------------------
methodB_150m <- methodB_150m %>%
  mutate(
    season = paste(gisid, year_id, sep = "_"),
    survey_year = factor(survey_yea),
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    qpad = offset
  )

## 2. MODEL FITTING WITH MODEL-SPECIFIC SPATIAL TERMS ------------------------

# A. Full model with interaction
modB150int_base <- glmer(occupancy ~ offset(qpad) + prop_habitat_std * edge_density_std +
                            (1 | gisid / season) + (1 | survey_year),
                          data = methodB_150m, family = binomial)
residauto_int <- residuals(modB150int_base)
modelgam_int <- gam(residauto_int ~ s(Easting, Northing), family = gaussian, data = methodB_150m)
methodB_150m$spatial_autocov_int <- scale(fitted(modelgam_int))[,1]
modB150int <- glmer(occupancy ~ offset(qpad) + prop_habitat_std * edge_density_std + 
                       spatial_autocov_int + (1 | gisid / season) + (1 | survey_year), 
                     data = methodB_150m, family = binomial)
summary(modB150int)
# B. Habitat amount only
modB150_A <- glmer(occupancy ~ offset(qpad) + prop_habitat_std + spatial_autocov_int +
                      (1 | gisid / season) + (1 | survey_year), 
                    data = methodB_150m, family = binomial)

# C. Fragmentation only
modB150_F <- glmer(occupancy ~ offset(qpad) + edge_density_std + spatial_autocov_int +
                      (1 | gisid / season) + (1 | survey_year), 
                    data = methodB_150m, family = binomial)

# D. Amount + fragmentation (no interaction)
modB150_AF <- glmer(occupancy ~ offset(qpad) + prop_habitat_std + edge_density_std + 
                       spatial_autocov_int + (1 | gisid / season) + (1 | survey_year), 
                     data = methodB_150m, family = binomial)
summary(modB150_AF)

## 3. MODEL SELECTION AND COMPARISON ----------------------------------------

# AIC model selection
model_selectionB150 <- model.sel(modB150int, modB150_A, modB150_AF, modB150_F)
print(model_selectionB150)

# Likelihood ratio test
lrtest(modB150int, modB150_AF)

## 4. MODEL DIAGNOSTICS ----------------------------------------------------

# Spatial autocorrelation checks for each model
coord <- cbind(methodB_150m$Easting, methodB_150m$Northing)

# Full model
resid_int <- residuals(modB150int)
spline_corr_int <- ncf::spline.correlog(methodB_150m$Easting, methodB_150m$Northing,
                                        resid_int, resamp = 10, type = "boot") # change to 10 from 100 to make it faster 

# Amount only
resid_A <- residuals(modB150_A)
spline_corr_A <- ncf::spline.correlog(methodB_150m$Easting, methodB_150m$Northing,
                                      resid_A, resamp = 10, type = "boot")

# Fragmentation only
resid_F <- residuals(modB150_F)
spline_corr_F <- ncf::spline.correlog(methodB_150m$Easting, methodB_150m$Northing,
                                      resid_F, resamp = 10, type = "boot")

# Amount + Fragmentation
resid_AF <- residuals(modB150_AF)
spline_corr_AF <- ncf::spline.correlog(methodB_150m$Easting, methodB_150m$Northing,
                                       resid_AF, resamp = 10, type = "boot")

# Plot spatial correlograms
par(mfrow = c(2,2))
plot(spline_corr_int, ylim = c(-0.2, 0.2), main = "Interaction Model")
plot(spline_corr_A, ylim = c(-0.2, 0.2), main = "Amount Only")
plot(spline_corr_F, ylim = c(-0.2, 0.2), main = "Fragmentation Only")
plot(spline_corr_AF, ylim = c(-0.2, 0.2), main = "Amount + Fragmentation")
par(mfrow = c(1,1))

# Temporal autocorrelation checks
par(mfrow = c(2,2))
pacf(residuals(modB150int), main = "Interaction Model")
pacf(residuals(modB150_A), main = "Amount Only")
pacf(residuals(modB150_F), main = "Fragmentation Only")
pacf(residuals(modB150_AF), main = "Amount + Fragmentation")
par(mfrow = c(1,1))


## 6. VARIANCE PARTITIONING ------------------------------------------------

# Set up parallel processing
plan(multisession, workers = 18)

# Partition variance for full model
r2part <- partR2(
  modB150int,
  partvars = c("prop_habitat_std", "edge_density_std", "prop_habitat_std:edge_density_std"), # have to manually plot these results to only include the four model combinations 
  data = methodB_150m,
  R2_type = "marginal",
  nboot = 99,
  CI = 0.95,
  parallel = TRUE
)

# Plot variance partitioning results
forestplot(r2part, type = c("R2"))

## 7. VISUALIZATION -------------------------------------------------------

# Effect plots for best model (assuming interaction model)
plot_model(modB150int, type = "int", mdrt.values = "quart")
plot_model(modB150int, type = "pred")

# Coefficient plots for model comparison
plot_models(modB150int, modB150_A, modB150_AF, modB150_F, transform = NULL)

# Calculate R2 for all models
print("R2 values:")
r.squaredGLMM(modB150int)
r.squaredGLMM(modB150_A)
r.squaredGLMM(modB150_F)
r.squaredGLMM(modB150_AF)



# C Method------------------------------------------------------------------------------------------------- 

## 1. DATA PREPARATION -------------------------------------------------------
methodC_150m <- methodC_150m %>%
  mutate(
    season = paste(gisid, year_id, sep = "_"),
    survey_year = factor(survey_yea),
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    qpad = offset
  )

## 2. MODEL FITTING WITH MODEL-SPECIFIC SPATIAL TERMS ------------------------

# A. Full model with interaction
modC150int_base <- glmer(occupancy ~ offset(qpad) + prop_habitat_std * edge_density_std +
                            (1 | gisid / season) + (1 | survey_year),
                          data = methodC_150m, family = binomial)
residauto_int <- residuals(modC150int_base)
modelgam_int <- gam(residauto_int ~ s(Easting, Northing), family = gaussian, data = methodC_150m)
methodC_150m$spatial_autocov_int <- scale(fitted(modelgam_int))[,1]
modC150int <- glmer(occupancy ~ offset(qpad) + prop_habitat_std * edge_density_std + 
                       spatial_autocov_int + (1 | gisid / season) + (1 | survey_year), 
                     data = methodC_150m, family = binomial)
summary(modC150int)

# B. Habitat amount only
modC150_A <- glmer(occupancy ~ offset(qpad) + prop_habitat_std + spatial_autocov_int +
                      (1 | gisid / season) + (1 | survey_year), 
                    data = methodC_150m, family = binomial)
summary(modC150_A)

# C. Fragmentation only
modC150_F <- glmer(occupancy ~ offset(qpad) + edge_density_std + spatial_autocov_int +
                      (1 | gisid / season) + (1 | survey_year), 
                    data = methodC_150m, family = binomial)
summary(modC150_F)

# D. Amount + fragmentation (no interaction)
modC150_AF <- glmer(occupancy ~ offset(qpad) + prop_habitat_std + edge_density_std + 
                       spatial_autocov_int + (1 | gisid / season) + (1 | survey_year), 
                     data = methodC_150m, family = binomial)
summary(modC150_AF)
## 3. MODEL SELECTION AND COMPARISON ----------------------------------------

# AIC model selection
model_selectionC150 <- model.sel(modC150int, modC150_A, modC150_AF, modC150_F)
print(model_selectionC150)

# Likelihood ratio test
lrtest(modC150int, modC150_AF)

## 4. MODEL DIAGNOSTICS ----------------------------------------------------

# Spatial autocorrelation checks for each model
coord <- cbind(methodC_150m$Easting, methodC_150m$Northing)

# Full model
resid_int <- residuals(modC150int)
spline_corr_int <- ncf::spline.correlog(methodC_150m$Easting, methodC_150m$Northing,
                                        resid_int, resamp = 10, type = "boot") # change to 10 from 100 to make it faster 

# Amount only
resid_A <- residuals(modC150_A)
spline_corr_A <- ncf::spline.correlog(methodC_150m$Easting, methodC_150m$Northing,
                                      resid_A, resamp = 10, type = "boot")

# Fragmentation only
resid_F <- residuals(modC150_F)
spline_corr_F <- ncf::spline.correlog(methodC_150m$Easting, methodC_150m$Northing,
                                      resid_F, resamp = 10, type = "boot")

# Amount + Fragmentation
resid_AF <- residuals(modC150_AF)
spline_corr_AF <- ncf::spline.correlog(methodC_150m$Easting, methodC_150m$Northing,
                                       resid_AF, resamp = 10, type = "boot")

# Plot spatial correlograms
par(mfrow = c(2,2))
plot(spline_corr_int, ylim = c(-0.2, 0.2), main = "Interaction Model")
plot(spline_corr_A, ylim = c(-0.2, 0.2), main = "Amount Only")
plot(spline_corr_F, ylim = c(-0.2, 0.2), main = "Fragmentation Only")
plot(spline_corr_AF, ylim = c(-0.2, 0.2), main = "Amount + Fragmentation")
par(mfrow = c(1,1))

# Temporal autocorrelation checks
par(mfrow = c(2,2))
pacf(residuals(modC150int), main = "Interaction Model")
pacf(residuals(modC150_A), main = "Amount Only")
pacf(residuals(modC150_F), main = "Fragmentation Only")
pacf(residuals(modC150_AF), main = "Amount + Fragmentation")
par(mfrow = c(1,1))


## 6. VARIANCE PARTITIONING ------------------------------------------------

# Set up parallel processing
plan(multisession, workers = 18)

# Partition variance for full model
r2part <- partR2(
  modC150int,
  partvars = c("prop_habitat_std", "edge_density_std", "prop_habitat_std:edge_density_std"), # have to manually plot these results to only include the four model combinations 
  data = methodC_150m,
  R2_type = "marginal",
  nboot = 99,
  CI = 0.95,
  parallel = TRUE
)

# Plot variance partitioning results
forestplot(r2part, type = c("R2"))

## 7. VISUALIZATION -------------------------------------------------------

# Effect plots for best model (assuming interaction model)
plot_model(modC150int, type = "int", mdrt.values = "quart")
plot_model(modC150int, type = "pred")

# Coefficient plots for model comparison
plot_models(modC150int, modC150_A, modC150_AF, modC150_F, transform = NULL)

# Calculate R2 for all models
print("R2 values:")
r.squaredGLMM(modC150int)
r.squaredGLMM(modC150_A)
r.squaredGLMM(modC150_F)
r.squaredGLMM(modC150_AF)



# Compare all models-----------------------------------------------------------------------
# 500 m
# Method A
summary(modA500int)
summary(modA500_A)
summary(modA500_F)
summary(modA500_AF)
model_selectionA500 <- model.sel(modA500int, modA500_A, modA500_AF, modA500_F)
print(model_selectionA500)
print("R2 values:")
r.squaredGLMM(modA500int)
r.squaredGLMM(modA500_A)
r.squaredGLMM(modA500_F)
r.squaredGLMM(modA500_AF)
# Method B
summary(modB500int)
summary(modB500_A)
summary(modB500_F)
summary(modB500_AF)
model_selectionB500 <- model.sel(modB500int, modB500_A, modB500_AF, modB500_F)
print(model_selectionB500)
print("R2 values:")
r.squaredGLMM(modB500int)
r.squaredGLMM(modB500_A)
r.squaredGLMM(modB500_F)
r.squaredGLMM(modB500_AF)
# Method C
summary(modC500int)
summary(modC500_A)
summary(modC500_F)
summary(modC500_AF)
model_selectionC500 <- model.sel(modC500int, modC500_A, modC500_AF, modC500_F)
print(model_selectionC500)
print("R2 values:")
r.squaredGLMM(modC500int)
r.squaredGLMM(modC500_A)
r.squaredGLMM(modC500_F)
r.squaredGLMM(modC500_AF)

# 1000 m
# Method A
summary(modA1000int)
summary(modA1000_A)
summary(modA1000_F)
summary(modA1000_AF)
model_selectionA1000 <- model.sel(modA1000int, modA1000_A, modA1000_AF, modA1000_F)
print(model_selectionA1000)
print("R2 values:")
r.squaredGLMM(modA1000int)
r.squaredGLMM(modA1000_A)
r.squaredGLMM(modA1000_F)
r.squaredGLMM(modA1000_AF)
# Method B
summary(modB1000int)
summary(modB1000_A)
summary(modB1000_F)
summary(modB1000_AF)
model_selectionB1000 <- model.sel(modB1000int, modB1000_A, modB1000_AF, modB1000_F)
print(model_selectionB1000)
print("R2 values:")
r.squaredGLMM(modB1000int)
r.squaredGLMM(modB1000_A)
r.squaredGLMM(modB1000_F)
r.squaredGLMM(modB1000_AF)
# Method C
summary(modC1000int)
summary(modC1000_A)
summary(modC1000_F)
summary(modC1000_AF)
model_selectionC1000 <- model.sel(modC1000int, modC1000_A, modC1000_AF, modC1000_F)
print(model_selectionC1000)
print("R2 values:")
r.squaredGLMM(modC1000int)
r.squaredGLMM(modC1000_A)
r.squaredGLMM(modC1000_F)
r.squaredGLMM(modC1000_AF)

# 150 m
# Method A
summary(modA150int)
summary(modA150_A)
summary(modA150_F)
summary(modA150_AF)
model_selectionA150 <- model.sel(modA150int, modA150_A, modA150_AF, modA150_F)
print(model_selectionA150)
print("R2 values:")
r.squaredGLMM(modA150int)
r.squaredGLMM(modA150_A)
r.squaredGLMM(modA150_F)
r.squaredGLMM(modA150_AF)
# Method B
summary(modB150int)
summary(modB150_A)
summary(modB150_F)
summary(modB150_AF)
model_selectionB150 <- model.sel(modB150int, modB150_A, modB150_AF, modB150_F)
print(model_selectionB150)
print("R2 values:")
r.squaredGLMM(modB150int)
r.squaredGLMM(modB150_A)
r.squaredGLMM(modB150_F)
r.squaredGLMM(modB150_AF)
# Method C
summary(modC150int)
summary(modC150_A)
summary(modC150_F)
summary(modC150_AF)
model_selectionC150 <- model.sel(modC150int, modC150_A, modC150_AF, modC150_F)
print(model_selectionC150)
print("R2 values:")
r.squaredGLMM(modC150int)
r.squaredGLMM(modC150_A)
r.squaredGLMM(modC150_F)
r.squaredGLMM(modC150_AF)











# 1. Create prediction plots
# For habitat amount (holding edge density at mean)
new_data_amount <- data.frame(
  habitat_amount_std = seq(min(full_data_clean$habitat_amount_std),
                           max(full_data_clean$habitat_amount_std),
                           length.out = 100),
  edge_density_std = mean(full_data_clean$edge_density_std),
  annual_detection_rate = mean(full_data_clean$annual_detection_rate),
  autocov = mean(full_data_clean$autocov),
  offset = mean(log(full_data_clean$offset))
)

# For edge density (holding habitat amount at mean)
new_data_edge <- data.frame(
  habitat_amount_std = mean(full_data_clean$habitat_amount_std),
  edge_density_std = seq(min(full_data_clean$edge_density_std),
                         max(full_data_clean$edge_density_std),
                         length.out = 100),
  annual_detection_rate = mean(full_data_clean$annual_detection_rate),
  autocov = mean(full_data_clean$autocov),
  offset = mean(log(full_data_clean$offset))
)

# Get predictions and confidence intervals
pred_amount <- predict(model_full, newdata = new_data_amount, type = "link", se.fit = TRUE)
pred_edge <- predict(model_full, newdata = new_data_edge, type = "link", se.fit = TRUE)

# Convert to probabilities with CI
new_data_amount$pred_prob <- plogis(pred_amount$fit)
new_data_amount$lower_ci <- plogis(pred_amount$fit - 1.96 * pred_amount$se.fit)
new_data_amount$upper_ci <- plogis(pred_amount$fit + 1.96 * pred_amount$se.fit)

new_data_edge$pred_prob <- plogis(pred_edge$fit)
new_data_edge$lower_ci <- plogis(pred_edge$fit - 1.96 * pred_edge$se.fit)
new_data_edge$upper_ci <- plogis(pred_edge$fit + 1.96 * pred_edge$se.fit)

# 2. Calculate relative importance
# Using standardized coefficients and their standard errors
importance <- data.frame(
  Variable = c("Habitat Amount", "Edge Density"),
  Std_Coef = abs(coef(model_full)[c("habitat_amount_std", "edge_density_std")]),
  SE = summary(model_full)$coefficients[c("habitat_amount_std", "edge_density_std"), "Std. Error"]
)
importance$Relative_Importance <- importance$Std_Coef / sum(importance$Std_Coef)

# 3. Find thresholds using quantiles of raw data
thresholds <- data.frame(
  habitat_amount = quantile(full_data_clean$habitat_amount, probs = c(0.25, 0.5, 0.75)),
  edge_density = quantile(full_data_clean$edge_density, probs = c(0.25, 0.5, 0.75))
)

# Calculate occupancy probabilities at these thresholds
threshold_pred <- expand.grid(
  habitat_amount_std = scale(thresholds$habitat_amount)[,1],
  edge_density_std = scale(thresholds$edge_density)[,1]
)
threshold_pred$annual_detection_rate <- mean(full_data_clean$annual_detection_rate)
threshold_pred$autocov <- mean(full_data_clean$autocov)
threshold_pred$offset <- mean(log(full_data_clean$offset))

threshold_pred$pred_prob <- predict(model_full, newdata = threshold_pred, type = "response")

# Print results
print("Relative Importance of Variables:")
print(importance)

print("\nPredicted Occupancy at Different Thresholds:")
print(threshold_pred)

# Create plots using base R
par(mfrow = c(1,2))

# Habitat Amount Plot
plot(new_data_amount$habitat_amount_std, new_data_amount$pred_prob, type = "l",
     xlab = "Standardized Habitat Amount", ylab = "Predicted Occupancy Probability",
     main = "Effect of Habitat Amount")
lines(new_data_amount$habitat_amount_std, new_data_amount$lower_ci, lty = 2)
lines(new_data_amount$habitat_amount_std, new_data_amount$upper_ci, lty = 2)

# Edge Density Plot
plot(new_data_edge$edge_density_std, new_data_edge$pred_prob, type = "l",
     xlab = "Standardized Edge Density", ylab = "Predicted Occupancy Probability",
     main = "Effect of Edge Density")
lines(new_data_edge$edge_density_std, new_data_edge$lower_ci, lty = 2)
lines(new_data_edge$edge_density_std, new_data_edge$upper_ci, lty = 2)