# Correlation between metrics -----------------------
hyp1_150 <- read.csv("1_Data/150m/habitat_metrics_hyp1.csv")
hyp2_150 <- read.csv("1_Data/150m/habitat_metrics_hyp2.csv")
hyp3_150 <- read.csv("1_Data/150m/habitat_metrics_hyp3.csv")

hyp1_500 <- read.csv("1_Data/500m/habitat_metrics_hyp1.csv")
hyp2_500 <- read.csv("1_Data/500m/habitat_metrics_hyp2.csv")
hyp3_500 <- read.csv("1_Data/500m/habitat_metrics_hyp3.csv")

hyp1_1000 <- read.csv("1_Data/1000m/habitat_metrics_hyp1.csv")
hyp2_1000 <- read.csv("1_Data/1000m/habitat_metrics_hyp2.csv")
hyp3_1000 <- read.csv("1_Data/1000m/habitat_metrics_hyp3.csv")


# Load required libraries
library(corrplot)
library(car)
library(dplyr)

# Function to analyze correlations
analyze_correlations <- function(data) {
  # Select only numeric columns for correlation analysis
  numeric_cols <- sapply(data, is.numeric)
  metrics_data <- data[, numeric_cols]
  
  # Remove any columns that are all NA or constant
  metrics_data <- metrics_data[, apply(metrics_data, 2, function(x) {
    !all(is.na(x)) && length(unique(x)) > 1
  })]
  
  # Calculate correlation matrix
  cor_matrix <- cor(metrics_data, use = "pairwise.complete.obs")
  
  # Create correlation plot
  corrplot(cor_matrix, 
           method = "color",
           type = "upper",
           tl.col = "black",
           tl.srt = 45,
           addCoef.col = "black",
           number.cex = 0.7,
           diag = FALSE)
  
  # Find highly correlated pairs (|r| > 0.7)
  high_cors <- which(abs(cor_matrix) > 0.7 & abs(cor_matrix) < 1, arr.ind = TRUE)
  high_cor_pairs <- data.frame(
    var1 = rownames(cor_matrix)[high_cors[,1]],
    var2 = colnames(cor_matrix)[high_cors[,2]],
    correlation = cor_matrix[high_cors]
  )
  
  return(list(
    correlation_matrix = cor_matrix,
    high_correlations = high_cor_pairs
  ))
}

# Analyze both datasets
hyp1_cors <- analyze_correlations(hyp1_1000)
hyp2_cors <- analyze_correlations(hyp2_1000)
hyp3_cors <- analyze_correlations(hyp3_1000)

# Print highly correlated pairs for both hypotheses
cat("\nHighly correlated metrics in Hypothesis 1 dataset (|r| > 0.7):\n")
print(hyp1_cors$high_correlations)

cat("\nHighly correlated metrics in Hypothesis 2 dataset (|r| > 0.7):\n")
print(hyp2_cors$high_correlations)

cat("\nHighly correlated metrics in Hypothesis 3 dataset (|r| > 0.7):\n")
print(hyp3_cors$high_correlations)

