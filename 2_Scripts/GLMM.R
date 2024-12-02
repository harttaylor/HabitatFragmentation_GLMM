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



# Start with basic GLM + spatial autocorrelation
# Step 1: Initial GLM
# Fit the GLMM
## A Method---------------------------------------------------------------------------------------------------
glm_long <- glm(occupancy ~ offset(log(qpad)) + prop_habitat_std + edge_density_std, data = methodA_500m, family = binomial)
summary(glm_long)

# Step 2: Handle spatial autocorrelation
# Run GAM with residuals to add as autocovariate in my model 
residauto <- residuals(glm_long)
modelgam <- gam(residauto ~ s(Easting, Northing), family = gaussian, data = methodA_500m)
methodA_500m$autocov <- fitted(modelgam)

# Step 3: Add autocovariate term back into dataset and standardize 
methodA_500m <- methodA_500m %>%
  mutate(
    survey_year = factor(survey_yea),
    duration = factor(duration),
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    #nlsi_std = scale(nlsi, center = TRUE, scale = TRUE)[,1],
    #core_mn_std = scale(core_mn, center = TRUE, scale = TRUE)[,1],
    qpad = offset,
    autocov = scale(autocov, center = TRUE, scale = TRUE)[,1]
  )


modA500 <- glmer(occupancy ~ offset(log(qpad)) + prop_habitat_std + edge_density_std + autocov +
                     (1 | gisid) + (1 | survey_year), 
                   data = methodA_500m, family = binomial)
summary(modA500)

modA500_A <- glmer(occupancy ~ offset(log(qpad)) + prop_habitat_std + autocov +
                     (1 | gisid) + (1 | survey_year), 
                   data = methodA_500m, family = binomial)
summary(modA500_A)

modA500_F <- glmer(occupancy ~ offset(log(qpad)) + edge_density_std + autocov +
                       (1 | gisid) + (1 | survey_year), 
                     data = methodA_500m, family = binomial)
summary(modA500_F)

modA500_int <- glmer(occupancy ~ offset(log(qpad)) + prop_habitat_std*edge_density_std + autocov +
                       (1 | gisid) + (1 | survey_year), 
                     data = methodA_500m, family = binomial)
summary(modA500_int)

aic_table <- AIC(modA500, modA500_A, modA500_F, modA500_int)
aic_sorted <- aic_table[order(aic_table$AIC), ]
print(aic_sorted)

## B Method--------------------------------------------------------------------------------------------
# Step 1: Initial GLM
# Fit the GLMM
glm_long <- glm(occupancy ~ prop_habitat_std + edge_density_std + duration +
                  survey_year, data = methodB_500m, family = binomial)
summary(glm_long)

# Step 2: Handle spatial autocorrelation
# Run GAM with residuals to add as autocovariate in my model 
residauto <- residuals(glm_long)
modelgam <- gam(residauto ~ s(Easting, Northing), family = gaussian, data = methodB_500m)
methodB_500m$autocov <- fitted(modelgam)


# Step 3: Add autocovariate term back into dataset and standardize 
methodB_500m <- methodB_500m %>%
  mutate(
    # Factorize categorical variables
    survey_year = factor(survey_yea),
    duration = factor(duration),
    # Standardize continuous landscape metrics
    # Using scale() with center and scale parameters for more control
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    #nlsi_std = scale(nlsi, center = TRUE, scale = TRUE)[,1],
    #core_mn_std = scale(core_mn, center = TRUE, scale = TRUE)[,1],
    qpad = offset,
    autocov = scale(autocov, center = TRUE, scale = TRUE)[,1]
  )

modB500 <- glmer(occupancy ~ offset(log(qpad)) + prop_habitat_std + edge_density_std + autocov +
                     (1 | gisid) + (1 | survey_year), 
                   data = methodB_500m, family = binomial)
summary(modB500)

modB500_A <- glmer(occupancy ~ offset(log(qpad)) + prop_habitat_std + autocov +
                       (1 | gisid) + (1 | survey_year), 
                     data = methodB_500m, family = binomial)
summary(modB500_A)

modB500_F <- glmer(occupancy ~ offset(log(qpad)) + edge_density_std + autocov +
                       (1 | gisid) + (1 | survey_year), 
                     data = methodB_500m, family = binomial)
summary(modB500_F)

modB500_int <- glmer(occupancy ~ offset(log(qpad)) + prop_habitat_std*edge_density_std + autocov +
                         (1 | gisid) + (1 | survey_year), 
                       data = methodB_500m, family = binomial)
summary(modB500_int)

aic_table <- AIC(modB500, modB500_A, modB500_F, modB500_int)
aic_sortedB <- aic_table[order(aic_table$AIC), ]
print(aic_sortedB)


# C Method------------------------------------------------------------------------------------------------- 
# Step 1: Initial GLM
# Fit the GLMM
modC500int_NOxy <- glmer(occupancy ~ offset(qpad) + prop_habitat_std * edge_density_std +
                      (1 | gisid / season) + (1 | survey_year),
                    data = methodC_500m, family = binomial)
summary(modC500int_NOxy)

# Step 2: Handle spatial autocorrelation
# Run GAM with residuals to add as autocovariate in my model 
residauto <- residuals(modC500int_NOxy)
modelgam <- gam(residauto ~ s(Easting, Northing), family = gaussian, data = methodC_500m)
methodC_500m$autocov <- fitted(modelgam)


# Step 3: Add autocovariate term back into dataset and standardize 
methodC_500m <- methodC_500m %>%
  mutate(
    survey_year = factor(survey_yea),
    season = paste(gisid, year_id, sep = "_"),
    duration = factor(duration),
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    #nlsi_std = scale(nlsi, center = TRUE, scale = TRUE)[,1],
    #core_mn_std = scale(core_mn, center = TRUE, scale = TRUE)[,1],
    qpad = offset,
    autocov = scale(autocov, center = TRUE, scale = TRUE)[,1]
  )

modC500int <- glmer(occupancy ~ offset(qpad) + prop_habitat_std * edge_density_std + autocov + scale(as.numeric(survey_year)) +
                     (1 | gisid / season) + (1 | survey_year), 
                   data = methodC_500m, family = binomial)
summary(modC500int)


# Habitat amt only 
modC500_A <- glmer(occupancy ~ offset(qpad) + prop_habitat_std + autocov +
                     (1 | gisid / season) + (1 | survey_year), 
                   data = methodC_500m, family = binomial)
summary(modC500_A)


# Frag only 
modC500_F <- glmer(occupancy ~ offset(qpad) + edge_density_std + autocov +
                       (1 | gisid / season) + (1 | survey_year), 
                     data = methodC_500m, family = binomial)
summary(modC500_F)

# amt + frag 
modC500_AF <- glmer(occupancy ~ offset(qpad) + prop_habitat_std + edge_density_std + autocov +
                         (1 | gisid / season) + (1 | survey_year), 
                       data = methodC_500m, family = binomial)
summary(modC500_AF)

model.sel(modC500int, modC500_A, modC500_AF, modC500_F)
  
plot_models(modC500int, modC500_A, modC500_AF, modC500_F, transform = NULL) # can make this pretty and use it as supp mat 
# Forest coefficient plots for the four model varainats showing standardized effect sizes in the link scale 

# Marginal model predictions (response curves)
plot_model(modC500int, type = "int", mdrt.values = "quart")

plot_model(modC500int, type = "pred")

library(lmtest)

lrtest(modC500int, modC500_AF)

# Variance partitioning 
library(partR2)
library(parallel)
library(future)
library(furrr)
parallel::detectCores()
plan(multisession, workers=18)
glmm.hp(modC500_AF)
?partR2
summary(modC500int)
r2part <- partR2(
  modC500int,
  partvars = c("prop_habitat_std", "edge_density_std", "prop_habitat_std:edge_density_std"),
  data = methodC_500m,
  R2_type = "marginal",
  max_level = NULL,
  nboot = 99,
  CI = 0.95,
  parallel = T,
  expct = "meanobs",
  olre = F,
  partbatch = NULL,
  allow_neg_r2 = FALSE
)
(y <- forestplot(r2part, type = c("R2")))

# Calculate R2 for mixed models (both marginal and conditional)
r.squaredGLMM(modA500)
r.squaredGLMM(modB500)
r.squaredGLMM(modC500)



# spatial autocorrelation test 
library(ncf)
install.packages("pgirmess")
library(pgirmess)
library(spdep) 
resid <- residuals(modC500int)
residnospace <- residuals(modC500_NOxy)
coord <- cbind(methodC_500m$Easting, methodC_500m$Northing)

correlog_noXY <- pgirmess::correlog(coords = coord, z = residnospace, method = "Moran", nbclass = 10)
spline_corr_noXY <- ncf::spline.correlog(methodC_500m$Easting, methodC_500m$Northing, 
                                    residnospace, resamp = 100, type = "boot")
plot(spline_corr_noXY, ylim = c(-0.2, 0.2))
correlog <- pgirmess::correlog(coords = coord, z = resid, method = "Moran", nbclass = 10)
plot(correlog)
spline_corr <- ncf::spline.correlog(methodC_500m$Easting, methodC_500m$Northing, 
                                         resid, resamp = 100, type = "boot")

plot(spline_corr, ylim = c(-0.2, 0.2))


# temporal autocorrelation checks 
pacf(residuals(modC500int))

# Diagnostics ---------------------------------------------------------------------------------------- 
# Function to perform comprehensive GLMM diagnostics
check_glmm <- function(model, modelname = "Model") {
  
  # 1. Create DHARMa residuals
  sim_resid <- simulateResiduals(model)
  
  # 2. Basic DHARMa diagnostic plots
  plot(sim_resid)
  
  # 3. Check for zero-inflation
  testZeroInflation(sim_resid)
  
  # 4. Check for dispersion
  testDispersion(sim_resid)
  
  # 5. Spatial autocorrelation in residuals
  acf(residuals(model), main = paste(modelname, "- ACF of Residuals"))
  
  # 7. Check variance inflation factors for fixed effects
  print("Variance Inflation Factors:")
  print(check_collinearity(model))
  
  # 8. Check model convergence
  print("Model convergence check:")
  print(check_convergence(model))
  
  # 9. Check singularity
  print("Singularity check:")
  print(check_singularity(model))
  
  # 10. Performance metrics
  print("Performance metrics:")
  print(model_performance(model))
  
}

# Usage example for your models:
# Individual model checks
check_glmm(modA500, "Model A")
check_glmm(modB500, "Model B")
check_glmm(modC500, "Model C")











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