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
remove.packages("glmmTMB")
install.packages('TMB', type = 'source')

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
glm_long <- glm(occupancy ~ prop_habitat_std + edge_density_std + duration +
                  survey_year, data = methodA_500m, family = binomial)
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


modA500 <- glmmTMB(occupancy ~ offset(log(qpad)) + prop_habitat_std + edge_density_std + autocov +
                     (1 | gisid) + (1 | survey_year), 
                   data = methodA_500m, family = binomial)
summary(modA500)




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

modB500 <- glmmTMB(occupancy ~ offset(log(qpad)) + prop_habitat_std + edge_density_std + autocov +
                     (1 | gisid) + (1 | survey_year), 
                   data = methodB_500m, family = binomial)
summary(modB500)

# C Method------------------------------------------------------------------------------------------------- 
# Step 1: Initial GLM
# Fit the GLMM
glm_long <- glm(occupancy ~ prop_habitat_std + edge_density_std + duration +
                  survey_year, data = methodC_500m, family = binomial)
summary(glm_long)

# Step 2: Handle spatial autocorrelation
# Run GAM with residuals to add as autocovariate in my model 
residauto <- residuals(glm_long)
modelgam <- gam(residauto ~ s(Easting, Northing), family = gaussian, data = methodC_500m)
methodC_500m$autocov <- fitted(modelgam)


# Step 3: Add autocovariate term back into dataset and standardize 
methodC_500m <- methodC_500m %>%
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

modC500 <- glmmTMB(occupancy ~ offset(log(qpad)) + prop_habitat_std + edge_density_std + autocov +
                     (1 | gisid) + (1 | survey_year), 
                   data = methodC_500m, family = binomial)
summary(modC500)
  
# Check model summaries
summary(modA500)
summary(modB500)
summary(modC500)


# Model Comparison
# Compare AIC
print("AIC comparison:")
print(AIC(modA500, modB500, modC500))

# Calculate R2 for mixed models (both marginal and conditional)
r.squaredGLMM(modA500)
r.squaredGLMM(modB500)
r.squaredGLMM(modC500)










# Retain detection histories 
# Load data
visit_matrix <- read.csv("1_Data/2yearvisitmatrix.csv")
detection_covariates <- read.csv("1_Data/visitsfordetcovs.csv")
habitat_metricsA <- read.csv("~/Chapter 2/Chapter 2 Analysis/Extract_patch_metrics/0_data/2_combined/564habitat_metrics_hyp1.csv")
habitat_metricsB <- read.csv("~/Chapter 2/Chapter 2 Analysis/Extract_patch_metrics/0_data/2_combined/564habitat_metrics_hyp2.csv")
habitat_metrics <- read.csv("~/Chapter 2/Chapter 2 Analysis/Extract_patch_metrics/0_data/2_combined/landscape_metrics/500landscape_metrics_hyp3.csv")

# Step 1: Pivot Visit Matrix to Long Format
visit_long <- visit_matrix %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "visit_id", 
    values_to = "occupancy"
  ) %>%
  # Extract year and visit number using substr
  mutate(
    year = as.numeric(substr(visit_id, 2, 2)),
    visit_num = as.numeric(substr(visit_id, 4, 4)))
    
    
# Step 2: Join Detection Covariates
detection_data <- detection_covariates %>%
  rename(survey_year = year) %>%
  rename(year = year_id) %>%
  select(gisid, surveyid, survey_year, year, visit_num, offset, BTNW_binary) %>%
  inner_join(visit_long, by = c("gisid", "year", "visit_num"))
head(detection_data)
head(habitat_metrics_expanded)

expanded_habitat <- habitat_metrics %>%
  group_by(surveyid, survey_yea) %>%
  expand(visit_num = 1:3) %>%  # Add rows for three visits
  left_join(habitat_metrics, by = c("surveyid", "survey_yea"))
head(expanded_habitat)
# Step 3: Merge Habitat Metrics (Site × Year Level)
final_data <- detection_data %>%
  inner_join(habitat_metrics_expanded, by = c("surveyid", "visit_num")) %>% 
  arrange(surveyid)

write.csv(habitat_metrics_expanded, "3_Outputs/habitat_metrics_expanded.csv")
write.csv(detection_data, "3_Outputs/detection_data.csv")
# Step 4: Ensure Habitat Metrics Are Site × Year-Specific
# Add standardized habitat covariates
final_data <- final_data %>%
  # Standardize continuous covariates and create factors
  mutate(
    # Standardize continuous variables
    habitat_amount_std = as.vector(scale(habitat_amount)),
    num_patches_std = as.vector(scale(num_patches)),
    edge_density_std = as.vector(scale(edge_density)),
    
    # Create factors
    survey_year = factor(survey_year),
    gisid = factor(gisid)
  )
write.csv(final_data, "1_Data/glmm/methodC.csv")
glm_detC <- brm(
  occupancy ~ offset(log(offset)) + habitat_amount_std + edge_density_std + autocov +
    (1 | gisid) + (1 | survey_year),
  family = bernoulli(),
  data = final_data,
  prior = c(
    prior(normal(0, 0.5), class = "b"),
    prior(normal(0, 0.5), class = "sd", lb = 0)
  ))#,
  control = list(adapt_delta = 0.95)  # Increased from default 0.8
)

summary(glm_detC)

glm_detCint <- brm(
  occupancy ~ offset(log(offset)) + habitat_amount_std:edge_density_std + autocov +
    (1 | gisid) + (1 | survey_year),
  family = bernoulli(),
  data = final_data,
  prior = c(
    prior(normal(0, 0.5), class = "b"),
    prior(normal(0, 0.5), class = "sd", lb = 0)
  )
)
summary(glm_detCint)

glm_detCamt <- brm(
  occupancy ~ offset(log(offset)) + habitat_amount_std + autocov +
    (1 | gisid) + (1 | survey_year),
  family = bernoulli(),
  data = final_data,
  prior = c(
    prior(normal(0, 0.5), class = "b"),
    prior(normal(0, 0.5), class = "sd", lb = 0)
  )
)
summary(glm_detCamt)

glm_detCfrag <- brm(
  occupancy ~ offset(log(offset)) + edge_density_std + autocov +
    (1 | gisid) + (1 | survey_year),
  family = bernoulli(),
  data = final_data,
  prior = c(
    prior(normal(0, 0.5), class = "b"),
    prior(normal(0, 0.5), class = "sd", lb = 0)
  )
)
summary(glm_detCfrag)


# Add LOO criterion to models
glm_detC <- add_criterion(glm_detC, "loo")
glm_detCint <- add_criterion(glm_detCint, "loo")
glm_detCamt <- add_criterion(glm_detCamt, "loo")
glm_detCfrag <- add_criterion(glm_detCfrag, "loo")

# Now compare
model_comparison <- loo_compare(
  glm_detC,
  glm_detCint,
  glm_detCamt, 
  glm_detCfrag
)

# Variance partitioning for full model
R2_full <- bayes_R2(glm_detC)
var_partition <- posterior_summary(VarCorr(glm_detC))

# Get random effects variance
ranef_summary <- VarCorr(glm_detC)
print(ranef_summary)

# Get fixed effects variance
fixef_posterior <- fixef(glm_detC)
var_habitat <- fixef_posterior["habitat_amount_std", "Estimate"]^2
var_frag <- fixef_posterior["edge_density_std", "Estimate"]^2

# Compare magnitudes
effect_comparison <- data.frame(
  habitat = abs(posterior$b_habitat_amount_std),
  fragmentation = abs(posterior$b_edge_density_std)
)
prob_hab_larger <- mean(effect_comparison$habitat > effect_comparison$fragmentation)


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