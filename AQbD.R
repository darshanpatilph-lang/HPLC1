
# Install packages (run once)
packages <- c(
  "readxl",
  "dplyr",
  "ggplot2",
  "rsm",
  "car",
  "patchwork",
  "writexl",
  "tidyr",
  "corrplot"
)

installed <- rownames(installed.packages())

for (p in packages) {
  if (!p %in% installed) {
    install.packages(p)
  }
}

# Load packages
library(readxl)
library(dplyr)
library(ggplot2)
library(rsm)
library(car)
library(patchwork)
library(writexl)
library(tidyr)
library(corrplot)

install.packages("writexl")
install.packages("corrplot")



# Import Excel master table
getwd()
setwd("D:/Academics/research article/HPLC MID TRE/")


data <- read_excel(
  "AQbD_BBD_Numerical_Master_Table.xlsx",
  sheet = "BBD_Master_Data"
)

# Inspect dataset
View(data)

# Check dimensions
dim(data)

# Check column names
colnames(data)

# Check missing values
colSums(is.na(data))

# Display dataset
print(data)


# Calculate experimental resolution
data <- data %>%
  mutate(
    Rs_calculated =
      2 * (TRE_RT_min - MID_RT_min) /
      (MID_baseline_width_min + TRE_baseline_width_min)
  )

# Compare reported and calculated resolution
data <- data %>%
  mutate(
    Rs_difference =
      Rs_calculated - Rs_reported
  )

# Display resolution verification
data %>%
  select(
    Run,
    Rs_reported,
    Rs_calculated,
    Rs_difference
  ) %>%
  print()

# Check discrepancies greater than 0.02
data %>%
  filter(abs(Rs_difference) > 0.02)

# Code experimental factors
data <- data %>%
  mutate(
    x1 = (Flow_mL_min - 1.00) / 0.10,
    x2 = (ACN_percent - 70.0) / 2.0,
    x3 = (Temperature_C - 40.0) / 2.0
  )



# Descriptive statistics
descriptive_stats <- data %>%
  summarise(
    Mean_Rs = mean(Rs_calculated),
    SD_Rs = sd(Rs_calculated),
    Min_Rs = min(Rs_calculated),
    Max_Rs = max(Rs_calculated),
    
    Mean_MID_RT = mean(MID_RT_min),
    SD_MID_RT = sd(MID_RT_min),
    
    Mean_TRE_RT = mean(TRE_RT_min),
    SD_TRE_RT = sd(TRE_RT_min),
    
    Mean_MID_N = mean(MID_theoretical_plates),
    Mean_TRE_N = mean(TRE_theoretical_plates),
    
    Mean_MID_Tf = mean(MID_tailing_factor),
    Mean_TRE_Tf = mean(TRE_tailing_factor)
  )

print(descriptive_stats)

# Save descriptive statistics
write.csv(
  descriptive_stats,
  "AQbD_Descriptive_Statistics.csv",
  row.names = FALSE
)



# Fit quadratic response surface model
model <- rsm(
  Rs_calculated ~ SO(x1, x2, x3),
  data = data
)

# Model summary
summary(model)

# ANOVA
anova_results <- anova(model)

print(anova_results)

# Regression coefficients
coefficients <- coef(model)

print(coefficients)

# Model performance
model_summary <- summary(model)

R_squared <- model_summary$r.squared
Adjusted_R_squared <- model_summary$adj.r.squared

Predicted_Rs <- predict(model, newdata = data)

Residuals <- data$Rs_calculated - Predicted_Rs

RMSE <- sqrt(mean(Residuals^2))

cat("R-squared:", R_squared, "\n")
cat("Adjusted R-squared:", Adjusted_R_squared, "\n")
cat("RMSE:", RMSE, "\n")



# Generate predicted values
data$Predicted_Rs <- predict(model)

# Experimental vs predicted plot
p1 <- ggplot(
  data,
  aes(
    x = Rs_calculated,
    y = Predicted_Rs
  )
) +
  geom_point(size = 3) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  labs(
    title = "Experimental versus predicted resolution",
    x = "Experimental resolution (Rs)",
    y = "Predicted resolution (Rs)"
  ) +
  theme_classic(base_size = 13)

print(p1)

ggsave(
  "Figure_9c_Experimental_vs_Predicted.png",
  p1,
  width = 7,
  height = 5,
  dpi = 600
)



# Install if required
install.packages("plotly")

library(plotly)
# Create prediction grid
grid_p2 <- expand.grid(
  x1 = seq(-1, 1, length.out = 50),
  x2 = seq(-1, 1, length.out = 50),
  x3 = 0
)

# Predict resolution
grid_p2$Predicted_Rs <- predict(
  model,
  newdata = grid_p2
)

# Convert coded variables into actual values
grid_p2$Flow <- 1.00 + 0.10 * grid_p2$x1
grid_p2$ACN <- 70.0 + 2.0 * grid_p2$x2

# Convert predicted values into matrix
z_matrix <- matrix(
  grid_p2$Predicted_Rs,
  nrow = 50,
  ncol = 50
)

# Generate 3D response surface
p2 <- plot_ly(
  x = seq(0.90, 1.10, length.out = 50),
  y = seq(68, 72, length.out = 50),
  z = z_matrix,
  type = "surface"
) %>%
  layout(
    title = "Effect of Flow Rate and ACN on Resolution",
    scene = list(
      xaxis = list(title = "Flow rate (mL/min)"),
      yaxis = list(title = "ACN concentration (%)"),
      zaxis = list(title = "Predicted resolution (Rs)")
    )
  )

p2
