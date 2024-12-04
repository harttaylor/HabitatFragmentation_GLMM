# Supplementary Analysis - Correlations and Robust Effects 

## Method C ---------------------------------------------------------------------
### 1. ASSESS CORRELATION BETWEEN PREDICTORS ------------------------------------

# Calculate correlation and VIF
cor_val <- cor(methodC_500m$prop_habitat_std, methodC_500m$edge_density_std)

# Calculate VIF (using basic model without spatial term)
mod_vif <- lm(occupancy ~ prop_habitat_std * edge_density_std, data = methodC_500m)
vif_vals <- car::vif(mod_vif)

# Visualize correlation
correlation_plot <- ggplot(methodC_500m, 
                           aes(x = edge_density_std, y = prop_habitat_std)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "red") +
  labs(x = "Edge Density (standardized)",
       y = "Habitat Amount (standardized)",
       title = "Correlation between Habitat Amount and Fragmentation") +
  theme_bw() +
  annotate("text", x = max(methodC_500m$edge_density_std), 
           y = min(methodC_500m$prop_habitat_std),
           label = paste("r =", round(cor_val, 2)))

### 2. SEPARATE EFFECTS USING GAM ---------------------------------------------

# Fit GAMs
amount_gam <- gam(prop_habitat_std ~ s(edge_density_std), 
                  family = gaussian, data = methodC_500m)
fragment_gam <- gam(edge_density_std ~ s(prop_habitat_std), 
                    family = gaussian, data = methodC_500m)

# Store residuals (uncorrelated components)
methodC_500m$fragres <- residuals(fragment_gam)
methodC_500m$amountres <- residuals(amount_gam)

### 3. REFIT MODELS WITH UNCORRELATED PREDICTORS ------------------------------

# Base models with residualized predictors
modC500_res_base <- glmer(occupancy ~ offset(qpad) + amountres * fragres +
                            (1 | gisid / season),
                          data = methodC_500m, family = binomial)

# Create spatial autocovariate
residauto_res <- residuals(modC500_res_base)
modelgam_res <- gam(residauto_res ~ s(Easting, Northing), 
                    family = gaussian, data = methodC_500m)
methodC_500m$spatial_autocov_res <- scale(fitted(modelgam_res))[,1]

# Final model with residualized predictors
modC500_res <- glmer(occupancy ~ offset(qpad) + amountres * fragres + 
                       spatial_autocov_res + (1 | gisid / season),
                     data = methodC_500m, family = binomial)

# Compare coefficients
coef_comparison <- bind_rows(
  "Original" = as.data.frame(t(fixef(modC500int))),
  "Residualized" = as.data.frame(t(fixef(modC500_res))),
  .id = "Model"
)

### 4. VISUALIZATION OF RESULTS ----------------------------------------------

# Plot coefficient comparison
coef_plot <- plot_models(modC500int, modC500_res,
                         transform = NULL,
                         m.labels = c("Original Model", "Residualized Model")) +
  labs(title = "Comparison of Effect Sizes") +
  theme_bw()

# Effect plots for residualized model
effect_plot_res <- plot_model(modC500_res, type = "pred") +
  theme_bw()