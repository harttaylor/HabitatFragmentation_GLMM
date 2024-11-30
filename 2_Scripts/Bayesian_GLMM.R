# Load required packages
library(brms)
library(bayesplot)
library(tidybayes)
library(tidyverse)
library(ggplot2)

# Model with QPAD offset
model_qpad <- brm(
  bf(occupancy ~ offset(log(offset)) + prop_habitat_std + edge_density_std + 
       autocov + (1 | gisid) + (1 | survey_year),
     family = bernoulli()),
  data = long_dataA500,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.95),
  prior = c(
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 2), class = "sd")
  ),
  save_pars = save_pars(all = TRUE)  # Important for LOO calculation
)

# Model with raw duration
model_duration <- brm(
  bf(occupancy ~ duration + prop_habitat_std + edge_density_std + 
       autocov + (1 | gisid) + (1 | survey_year),
     family = bernoulli()),
  data = long_dataA500,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.95),
  prior = c(
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 2), class = "sd")
  ),
  save_pars = save_pars(all = TRUE)  # Important for LOO calculation
)

# Model diagnostics function with corrected function names
check_model_fit <- function(model) {
  # Convergence diagnostics
  print("Checking Rhat values...")
  rhat_vals <- brms::rhat(model)
  print(summary(as.numeric(rhat_vals)))
  
  # Effective sample sizes
  print("Checking effective sample sizes...")
  neff_vals <- brms::neff_ratio(model)
  print(summary(as.numeric(neff_vals)))
  
  # LOO calculation
  print("Computing LOO...")
  loo_result <- brms::loo(model)
  print(loo_result)
  
  # Posterior predictive checks
  print("Creating posterior predictive check plot...")
  pp_check(model, ndraws = 50)
  
  # Extract random effects
  print("Extracting random effects...")
  re_checks <- VarCorr(model)
  print(re_checks)
  
  # Return full diagnostics
  list(
    rhat = rhat_vals,
    neff = neff_vals,
    loo = loo_result,
    ranef = re_checks
  )
}

# Get detailed parameter summaries
parameter_summary <- get_parameter_summary(model_qpad)
parameter_summary <- get_parameter_summary(model_duration)

# Run basic diagnostics
qpad_diagnostics <- check_model_fit(model_qpad) # loo = 1607.6 +- 62.2
duration_diagnostics <- check_model_fit(model_duration) # loo = 1580.3 +- 60.8


# Now proceed with duration model 

# 500 m scale 
# Method A
# Define model settings
model_settings <- list(
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.95),
  prior = c(
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 2), class = "sd")
  ),
  save_pars = save_pars(all = TRUE)
)
# Fit models for Method A at 500m
modelA500_amt <- brm(
  bf(occupancy ~  prop_habitat_std + duration +
       autocov + (1 | gisid) + (1 | survey_year),
     family = bernoulli()),
  data = methodA_500m,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.95),
  prior = c(
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 2), class = "sd")
  ),
  save_pars = save_pars(all = TRUE)  # Important for LOO calculation
)


modelA500_frag <- brm(
  bf(occupancy ~ edge_density_std + duration +
       autocov + (1 | gisid) + (1 | survey_year),
     family = bernoulli()),
  data = methodA_500m,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.95),
  prior = c(
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 2), class = "sd")
  ),
  save_pars = save_pars(all = TRUE)  # Important for LOO calculation
)

modelA500_both <- brm(
  bf(occupancy ~ prop_habitat_std + edge_density_std + duration + 
       autocov + (1 | gisid) + (1 | survey_year),
     family = bernoulli()),
  data = methodA_500m,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.95),
  prior = c(
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 2), class = "sd")
  ),
  save_pars = save_pars(all = TRUE)  # Important for LOO calculation
)

modelA500_int <- brm(
  bf(occupancy ~ prop_habitat_std * edge_density_std + duration + 
       autocov + (1 | gisid) + (1 | survey_year),
     family = bernoulli()),
  data = methodA_500m,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.95),
  prior = c(
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 2), class = "sd")
  ),
  save_pars = save_pars(all = TRUE)  # Important for LOO calculation
)

# Save models
#saveRDS(modelA500_amt, " /methodA_500m_amount.rds")
#saveRDS(modelA500_frag, "models/methodA_500m_frag.rds")
#saveRDS(modelA500_both, "models/methodA_500m_both.rds")
#saveRDS(modelA500_int, "models/methodA_500m_int.rds")

# Compare model fits
models_A500 <- list(
  amount = modelA500_amt,
  frag = modelA500_frag,
  both = modelA500_both,
  interaction = modelA500_int
)

# Add LOO to all models and compare
models_loo <- lapply(models_A500, add_criterion, criterion = "loo")
model_comparison <- loo_compare(models_loo)

# Calculate R² for all models
models_R2 <- lapply(models_A500, bayes_R2)

# Extract and compare effect sizes
extract_effects <- function(model, model_name) {
  fixef_summary <- fixef(model)
  data.frame(
    model = model_name,
    parameter = rownames(fixef_summary),
    estimate = fixef_summary[,"Estimate"],
    sd = fixef_summary[,"Est.Error"],
    lower = fixef_summary[,"Q2.5"],
    upper = fixef_summary[,"Q97.5"]
  )
}

effects_comparison <- do.call(rbind, mapply(
  extract_effects,
  models_A500,
  names(models_A500),
  SIMPLIFY = FALSE
))

# Calculate variance explained by fixed effects
calc_var_explained <- function(model) {
  # Extract posterior samples
  post <- posterior_samples(model)
  # Get columns for fixed effects (excluding Intercept)
  fx_cols <- grep("^b_(?!Intercept)", names(post), value = TRUE, perl = TRUE)
  # Calculate variance explained by each fixed effect
  sapply(fx_cols, function(col) var(post[[col]]))
}

var_explained <- lapply(models_A500, calc_var_explained)

# Create summary results
results_A500 <- list(
  model_comparison = model_comparison,
  R2 = models_R2,
  effects = effects_comparison,
  variance_explained = var_explained
)

# Save results
#saveRDS(results_A500, "results/methodA_500m_results.rds")

# Print summary table
summary_table <- data.frame(
  Model = names(models_A500),
  LOOIC = sapply(models_loo, function(x) x$criteria$loo$estimates["looic","Estimate"]),
  R2_mean = sapply(models_R2, mean),
  R2_sd = sapply(models_R2, sd)
)

# Create effect size comparison plot
ggplot(effects_comparison[effects_comparison$parameter %in% 
                            c("prop_habitat_std", "edge_density_std"),], 
       aes(x = model, y = estimate, ymin = lower, ymax = upper)) +
  geom_pointrange() +
  facet_wrap(~parameter) +
  coord_flip() +
  theme_bw() +
  labs(title = "Effect Sizes Across Models - Method A 500m",
       y = "Standardized Effect Size",
       x = "Model")

# Print results
print("Model Comparison Results:")
print(summary_table)
print("\nEffect Sizes and Variance Explained:")
print(effects_comparison)
print("\nVariance Explained by Fixed Effects:")
print(var_explained)


# Variance partitioning for full model
R2_full <- bayes_R2(modelA500_both)
R2_amt <- bayes_R2(modelA500_amt)
R2_frag <- bayes_R2(modelA500_frag)
R2_int <- bayes_R2(modelA500_int)
var_partition <- posterior_summary(VarCorr(modelA500_both))
print(var_partition)
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




# Method B 

# Define model settings
model_settings <- list(
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.95),
  prior = c(
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 2), class = "sd")
  ),
  save_pars = save_pars(all = TRUE)
)
# Fit models for Method B at 500m
modelB500_amt <- brm(
  bf(occupancy ~  prop_habitat_std + duration +
       autocov + (1 | gisid) + (1 | survey_year),
     family = bernoulli()),
  data = methodB_500m,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.95),
  prior = c(
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 2), class = "sd")
  ),
  save_pars = save_pars(all = TRUE)  # Important for LOO calculation
)


modelB500_frag <- brm(
  bf(occupancy ~ edge_density_std + duration +
       autocov + (1 | gisid) + (1 | survey_year),
     family = bernoulli()),
  data = methodB_500m,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.95),
  prior = c(
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 2), class = "sd")
  ),
  save_pars = save_pars(all = TRUE)  # Important for LOO calculation
)

modelB500_both <- brm(
  bf(occupancy ~ prop_habitat_std + edge_density_std + duration + 
       autocov + (1 | gisid) + (1 | survey_year),
     family = bernoulli()),
  data = methodB_500m,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.95),
  prior = c(
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 2), class = "sd")
  ),
  save_pars = save_pars(all = TRUE)  # Important for LOO calculation
)

modelB500_int <- brm(
  bf(occupancy ~ prop_habitat_std * edge_density_std + duration + 
       autocov + (1 | gisid) + (1 | survey_year),
     family = bernoulli()),
  data = methodB_500m,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.95),
  prior = c(
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 2), class = "sd")
  ),
  save_pars = save_pars(all = TRUE)  # Important for LOO calculation
)

# Save models
#saveRDS(modelB500_amt, " /methodB_500m_amount.rds")
#saveRDS(modelB500_frag, "models/methodB_500m_frag.rds")
#saveRDS(modelB500_both, "models/methodB_500m_both.rds")
#saveRDS(modelB500_int, "models/methodB_500m_int.rds")

# Compare model fits
models_B500 <- list(
  amount = modelB500_amt,
  frag = modelB500_frag,
  both = modelB500_both,
  interaction = modelB500_int
)

# Add LOO to all models and compare
models_loo <- lapply(models_B500, add_criterion, criterion = "loo")
model_comparison <- loo_compare(models_loo)

# Calculate R² for all models
models_R2 <- lapply(models_B500, bayes_R2)

# Extract and compare effect sizes
extract_effects <- function(model, model_name) {
  fixef_summary <- fixef(model)
  data.frame(
    model = model_name,
    parameter = rownames(fixef_summary),
    estimate = fixef_summary[,"Estimate"],
    sd = fixef_summary[,"Est.Error"],
    lower = fixef_summary[,"Q2.5"],
    upper = fixef_summary[,"Q97.5"]
  )
}

effects_comparison <- do.call(rbind, mapply(
  extract_effects,
  models_B500,
  names(models_B500),
  SIMPLIFY = FALSE
))

# Calculate variance explained by fixed effects
calc_var_explained <- function(model) {
  # Extract posterior samples
  post <- posterior_samples(model)
  # Get columns for fixed effects (excluding Intercept)
  fx_cols <- grep("^b_(?!Intercept)", names(post), value = TRUE, perl = TRUE)
  # Calculate variance explained by each fixed effect
  sapply(fx_cols, function(col) var(post[[col]]))
}

var_explained <- lapply(models_B500, calc_var_explained)

# Create summary results
results_B500 <- list(
  model_comparison = model_comparison,
  R2 = models_R2,
  #effects = effects_comparison,
  variance_explained = var_explained
)

# Save results
#saveRDS(results_B500, "results/methodB_500m_results.rds")

# Print summary table
summary_table <- data.frame(
  Model = names(models_B500),
  LOOIC = sapply(models_loo, function(x) x$criteria$loo$estimates["looic","Estimate"]),
  R2_mean = sapply(models_R2, mean),
  R2_sd = sapply(models_R2, sd)
)

# Create effect size comparison plot
ggplot(effects_comparison[effects_comparison$parameter %in% 
                            c("prop_habitat_std", "edge_density_std"),], 
       aes(x = model, y = estimate, ymin = lower, ymax = upper)) +
  geom_pointrange() +
  facet_wrap(~parameter) +
  coord_flip() +
  theme_bw() +
  labs(title = "Effect Sizes Across Models - Method B 500m",
       y = "Standardized Effect Size",
       x = "Model")

# Print results
print("Model Comparison Results:")
print(summary_table)
print("\nEffect Sizes and Variance Explained:")
print(effects_comparison)
print("\nVariance Explained by Fixed Effects:")
print(var_explained)


# Variance partitioning for full model
R2_full <- bayes_R2(modelB500_both)
R2_amt <- bayes_R2(modelB500_amt)
R2_frag <- bayes_R2(modelB500_frag)
R2_int <- bayes_R2(modelB500_int)




# Method C

# Define model settings
model_settings <- list(
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.95),
  prior = c(
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 2), class = "sd")
  ),
  save_pars = save_pars(all = TRUE)
)
# Fit models for Method C at 500m
modelC500_amt <- brm(
  bf(occupancy ~  prop_habitat_std + duration +
       autocov + (1 | gisid) + (1 | survey_year),
     family = bernoulli()),
  data = methodC_500m,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.95),
  prior = c(
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 2), class = "sd")
  ),
  save_pars = save_pars(all = TRUE)  # Important for LOO calculation
)


modelC500_frag <- brm(
  bf(occupancy ~ edge_density_std + duration +
       autocov + (1 | gisid) + (1 | survey_year),
     family = bernoulli()),
  data = methodC_500m,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.95),
  prior = c(
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 2), class = "sd")
  ),
  save_pars = save_pars(all = TRUE)  # Important for LOO calculation
)

modelC500_both <- brm(
  bf(occupancy ~ prop_habitat_std + edge_density_std + duration + 
       autocov + (1 | gisid) + (1 | survey_year),
     family = bernoulli()),
  data = methodC_500m,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.95),
  prior = c(
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 2), class = "sd")
  ),
  save_pars = save_pars(all = TRUE)  # Important for LOO calculation
)

modelC500_int <- brm(
  bf(occupancy ~ prop_habitat_std * edge_density_std + duration + 
       autocov + (1 | gisid) + (1 | survey_year),
     family = bernoulli()),
  data = methodC_500m,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.95),
  prior = c(
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 2), class = "sd")
  ),
  save_pars = save_pars(all = TRUE)  # Important for LOO calculation
)

# Save models
#saveRDS(modelC500_amt, " /methodC_500m_amount.rds")
#saveRDS(modelC500_frag, "models/methodC_500m_frag.rds")
#saveRDS(modelC500_both, "models/methodC_500m_both.rds")
#saveRDS(modelC500_int, "models/methodC_500m_int.rds")

# Compare model fits
models_C500 <- list(
  #amount = modelC500_amt,
  #frag = modelC500_frag,
  both = modelC500_both
  #interaction = modelC500_int
)

# Add LOO to all models and compare
models_loo <- lapply(models_C500, add_criterion, criterion = "loo")
model_comparison <- loo_compare(models_loo)

# Calculate R² for all models
models_R2 <- lapply(models_C500, bayes_R2)

# Extract and compare effect sizes
extract_effects <- function(model, model_name) {
  fixef_summary <- fixef(model)
  data.frame(
    model = model_name,
    parameter = rownames(fixef_summary),
    estimate = fixef_summary[,"Estimate"],
    sd = fixef_summary[,"Est.Error"],
    lower = fixef_summary[,"Q2.5"],
    upper = fixef_summary[,"Q97.5"]
  )
}

effects_comparison <- do.call(rbind, mapply(
  extract_effects,
  models_C500,
  names(models_C500),
  SIMPLIFY = FALSE
))

# Calculate variance explained by fixed effects
calc_var_explained <- function(model) {
  # Extract posterior samples
  post <- posterior_samples(model)
  # Get columns for fixed effects (excluding Intercept)
  fx_cols <- grep("^b_(?!Intercept)", names(post), value = TRUE, perl = TRUE)
  # Calculate variance explained by each fixed effect
  sapply(fx_cols, function(col) var(post[[col]]))
}
posterior_samples(modelC500_both)
var_explained <- lapply(models_C500, calc_var_explained)

# Create summary results
results_C500 <- list(
  model_comparison = model_comparison,
  R2 = models_R2,
  #effects = effects_comparison,
  variance_explained = var_explained
)

# Save results
#saveRDS(results_C500, "results/methodC_500m_results.rds")

# Print summary table
summary_table <- data.frame(
  Model = names(models_C500),
  LOOIC = sapply(models_loo, function(x) x$criteria$loo$estimates["looic","Estimate"]),
  R2_mean = sapply(models_R2, mean),
  R2_sd = sapply(models_R2, sd)
)

# Create effect size comparison plot
ggplot(effects_comparison[effects_comparison$parameter %in% 
                            c("prop_habitat_std", "edge_density_std"),], 
       aes(x = model, y = estimate, ymin = lower, ymax = upper)) +
  geom_pointrange() +
  facet_wrap(~parameter) +
  coord_flip() +
  theme_bw() +
  labs(title = "Effect Sizes Across Models - Method C 500m",
       y = "Standardized Effect Size",
       x = "Model")

# Print results
print("Model Comparison Results:")
print(summary_table)
print("\nEffect Sizes and Variance Explained:")
print(effects_comparison)
print("\nVariance Explained by Fixed Effects:")
print(var_explained)


# Variance partitioning for full model
R2_full <- bayes_R2(modelC500_both)
R2_amt <- bayes_R2(modelC500_amt)
R2_frag <- bayes_R2(modelC500_frag)
R2_int <- bayes_R2(modelC500_int)

