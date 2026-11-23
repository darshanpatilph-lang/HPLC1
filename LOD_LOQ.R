# ==============================================================================
# STEP 1: Load Libraries and Create Data Frame
# ==============================================================================

# Load required analytical & statistical packages
library(lmtest)
library(car)
library(ggplot2)
library(gridExtra)

# Define concentration points (6 levels in triplicate)
conc <- c(0.06, 0.06, 0.06, 0.10, 0.10, 0.10, 0.50, 0.50, 0.50, 
          1.00, 1.00, 1.00, 1.50, 1.50, 1.50, 2.00, 2.00, 2.00)

# Define Peak Areas for Midostaurin
midostaurin <- c(1021.29, 936.42, 927.57, 1688.05, 1726.82, 1686.04, 
                 8450.38, 8426.75, 8453.75, 16763.79, 16983.58, 16873.69, 
                 25124.92, 25361.16, 25260.00, 33621.06, 33747.36, 33570.54)

# Define Peak Areas for Tretinoin
tretinoin <- c(949.14, 906.42, 915.35, 1593.15, 1580.98, 1591.92, 
               7953.00, 7983.27, 7938.67, 15829.71, 16210.76, 15750.12, 
               23546.43, 23798.43, 23499.18, 31230.41, 31441.90, 31277.41)

# Assemble into a clean data frame
df <- data.frame(Conc = conc, Midostaurin = midostaurin, Tretinoin = tretinoin)

# Add 1/x^2 weighting vector for later Weighted Least Squares analysis
df$Weight <- 1 / (df$Conc^2)

# Display the first few rows of the dataset
head(df, 6)

# ==============================================================================
# STEP 2: Fit OLS Model for Midostaurin and Extract Parameters
# ==============================================================================

# Fit the linear model: Peak Area = Slope * Conc + Intercept
mido_ols <- lm(Midostaurin ~ Conc, data = df)

# Display the summary of the linear model
summary(mido_ols)

# ==============================================================================
# STEP 3: Extract Intercept, Slope, SE, R^2, and Calculate OLS LOD/LOQ
# ==============================================================================

# Extract summary object
mido_summary <- summary(mido_ols)

# Extract parameters
slope_mido <- mido_summary$coefficients["Conc", "Estimate"]
intercept_mido <- mido_summary$coefficients["(Intercept)", "Estimate"]
intercept_se_mido <- mido_summary$coefficients["(Intercept)", "Std. Error"]
r2_mido <- mido_summary$r.squared

# Calculate LOD and LOQ via ICH Q2 formula
lod_mido_ols <- (3.3 * intercept_se_mido) / slope_mido
loq_mido_ols <- (10.0 * intercept_se_mido) / slope_mido

# Print extracted metrics
cat(sprintf("Slope (S): %.3f\n", slope_mido))
cat(sprintf("Intercept (a): %.3f (SE = %.3f)\n", intercept_mido, intercept_se_mido))
cat(sprintf("R-squared (R2): %.6f\n", r2_mido))
cat(sprintf("Calculated OLS LOD: %.4f ug/mL\n", lod_mido_ols))
cat(sprintf("Calculated OLS LOQ: %.4f ug/mL\n", loq_mido_ols))

# ==============================================================================
# MASTER METHOD VALIDATION TABLE GENERATOR FOR R
# Analytes: Midostaurin & Tretinoin (0.06 - 2.00 ug/mL)
# ==============================================================================

# 1. Input Experimental Dataset
conc <- c(0.06, 0.06, 0.06, 0.10, 0.10, 0.10, 0.50, 0.50, 0.50, 
          1.00, 1.00, 1.00, 1.50, 1.50, 1.50, 2.00, 2.00, 2.00)

midostaurin <- c(1021.29, 936.42, 927.57, 1688.05, 1726.82, 1686.04, 
                 8450.38, 8426.75, 8453.75, 16763.79, 16983.58, 16873.69, 
                 25124.92, 25361.16, 25260.00, 33621.06, 33747.36, 33570.54)

tretinoin <- c(949.14, 906.42, 915.35, 1593.15, 1580.98, 1591.92, 
               7953.00, 7983.27, 7938.67, 15829.71, 16210.76, 15750.12, 
               23546.43, 23798.43, 23499.18, 31230.41, 31441.90, 31277.41)

df <- data.frame(Conc = conc, Midostaurin = midostaurin, Tretinoin = tretinoin)
df$Weight <- 1 / (df$Conc^2)

# 2. Fit Ordinary Least Squares (OLS) and Weighted Least Squares (WLS 1/x^2)
mido_ols <- summary(lm(Midostaurin ~ Conc, data = df))
tret_ols <- summary(lm(Tretinoin ~ Conc, data = df))

mido_wls <- summary(lm(Midostaurin ~ Conc, weights = Weight, data = df))
tret_wls <- summary(lm(Tretinoin ~ Conc, weights = Weight, data = df))

# 3. Construct Master Validation Table
master_table <- data.frame(
  Parameter_Category = c(
    "--- I. CHROMATOGRAPHIC BEHAVIOR ---",
    "Target Retention Time (Rt, min)",
    "Evaluated Concentration Range (ug/mL)",
    "Total Calibration Injections (N)",
    
    "--- II. EXPERIMENTAL PEAK AREA RESPONSES (Mean ± SD) ---",
    "Conc 0.06 ug/mL (LOQ floor)",
    "Conc 0.10 ug/mL (Trace point)",
    "Conc 0.50 ug/mL (Low-mid standard)",
    "Conc 1.00 ug/mL (Geometric midpoint)",
    "Conc 1.50 ug/mL (Mid-high standard)",
    "Conc 2.00 ug/mL (Upper linear limit)",
    
    "--- III. REPEATABILITY & PRECISION (CV%) ---",
    "Repeatability at LOQ (0.06 ug/mL)",
    "Repeatability Range Across Assay (0.10 - 2.00 ug/mL)",
    
    "--- IV. SIGNAL-TO-NOISE (S/N) METRICS ---",
    "Mean S/N at 0.06 ug/mL (LOQ Threshold)",
    "Mean S/N at 2.00 ug/mL (Upper Limit)",
    
    "--- V. REGRESSION METRICS: ORDINARY LEAST SQUARES (OLS) ---",
    "OLS Linear Equation",
    "Slope (S, Area per ug/mL)",
    "Intercept (a)",
    "Standard Error of Intercept (SE_a)",
    "Coefficient of Determination (R²)",
    
    "--- VI. REGRESSION METRICS: WEIGHTED LEAST SQUARES (WLS 1/x²) ---",
    "WLS Linear Equation",
    "Weighted Slope (S_w)",
    "Weighted Intercept Standard Error (SE_aw)",
    "Weighted Coefficient of Determination (R²_w)",
    
    "--- VII. DERIVED SENSITIVITY LIMITS (LOD & LOQ) ---",
    "OLS Calculated LOD (3.3 * SE_a / S)",
    "OLS Calculated LOQ (10 * SE_a / S)",
    "WLS Calculated LOD (3.3 * SE_aw / S_w)",
    "WLS Calculated LOQ (10 * SE_aw / S_w)",
    "Empirical Operational LOQ Floor (S/N ≥ 10:1)"
  ),
  
  Midostaurin = c(
    "",
    "3.31 ± 0.03",
    "0.06 - 2.00",
    "18 (6 levels × 3 reps)",
    
    "",
    "961.76 ± 51.74",
    "1700.30 ± 22.99",
    "8443.63 ± 14.71",
    "16873.69 ± 109.90",
    "25248.69 ± 118.53",
    "33646.32 ± 91.08",
    
    "",
    "5.38%",
    "0.17% - 1.35%",
    
    "",
    "10.47",
    "119.43",
    
    "",
    sprintf("y = %.2fx + %.2f", mido_ols$coefficients["Conc","Estimate"], mido_ols$coefficients["(Intercept)","Estimate"]),
    sprintf("%.2f", mido_ols$coefficients["Conc","Estimate"]),
    sprintf("%.2f", mido_ols$coefficients["(Intercept)","Estimate"]),
    sprintf("%.2f", mido_ols$coefficients["(Intercept)","Std. Error"]),
    sprintf("%.5f", mido_ols$r.squared),
    
    "",
    sprintf("y = %.2fx + %.2f", mido_wls$coefficients["Conc","Estimate"], mido_wls$coefficients["(Intercept)","Estimate"]),
    sprintf("%.2f", mido_wls$coefficients["Conc","Estimate"]),
    sprintf("%.2f", mido_wls$coefficients["(Intercept)","Std. Error"]),
    sprintf("%.5f", mido_wls$r.squared),
    
    "",
    sprintf("%.4f ug/mL", (3.3 * mido_ols$coefficients["(Intercept)","Std. Error"]) / mido_ols$coefficients["Conc","Estimate"]),
    sprintf("%.4f ug/mL", (10.0 * mido_ols$coefficients["(Intercept)","Std. Error"]) / mido_ols$coefficients["Conc","Estimate"]),
    sprintf("%.4f ug/mL", (3.3 * mido_wls$coefficients["(Intercept)","Std. Error"]) / mido_wls$coefficients["Conc","Estimate"]),
    sprintf("%.4f ug/mL", (10.0 * mido_wls$coefficients["(Intercept)","Std. Error"]) / mido_wls$coefficients["Conc","Estimate"]),
    "0.0600 ug/mL"
  ),
  
  Tretinoin = c(
    "",
    "9.21 ± 0.05",
    "0.06 - 2.00",
    "18 (6 levels × 3 reps)",
    
    "",
    "923.64 ± 22.53",
    "1588.68 ± 6.70",
    "7958.31 ± 22.77",
    "15930.20 ± 246.21",
    "23614.68 ± 160.88",
    "31316.57 ± 111.05",
    
    "",
    "2.44%",
    "0.29% - 1.55%",
    
    "",
    "10.66",
    "111.29",
    
    "",
    sprintf("y = %.2fx + %.2f", tret_ols$coefficients["Conc","Estimate"], tret_ols$coefficients["(Intercept)","Estimate"]),
    sprintf("%.2f", tret_ols$coefficients["Conc","Estimate"]),
    sprintf("%.2f", tret_ols$coefficients["(Intercept)","Estimate"]),
    sprintf("%.2f", tret_ols$coefficients["(Intercept)","Std. Error"]),
    sprintf("%.5f", tret_ols$r.squared),
    
    "",
    sprintf("y = %.2fx + %.2f", tret_wls$coefficients["Conc","Estimate"], tret_wls$coefficients["(Intercept)","Estimate"]),
    sprintf("%.2f", tret_ols$coefficients["Conc","Estimate"]),
    sprintf("%.2f", tret_wls$coefficients["(Intercept)","Std. Error"]),
    sprintf("%.5f", tret_wls$r.squared),
    
    "",
    sprintf("%.4f ug/mL", (3.3 * tret_ols$coefficients["(Intercept)","Std. Error"]) / tret_ols$coefficients["Conc","Estimate"]),
    sprintf("%.4f ug/mL", (10.0 * tret_ols$coefficients["(Intercept)","Std. Error"]) / tret_ols$coefficients["Conc","Estimate"]),
    sprintf("%.4f ug/mL", (3.3 * tret_wls$coefficients["(Intercept)","Std. Error"]) / tret_wls$coefficients["Conc","Estimate"]),
    sprintf("%.4f ug/mL", (10.0 * tret_wls$coefficients["(Intercept)","Std. Error"]) / tret_wls$coefficients["Conc","Estimate"]),
    "0.0600 ug/mL"
  ),
  
  Regulatory_Acceptance_Criteria = c(
    "",
    "Baseline peak resolution",
    "33-fold operational window",
    "ICH Q2(R2) minimum requirement",
    
    "",
    "Empirical LOQ floor",
    "Trace quantitative standard",
    "Low-range linear point",
    "Nominal midpoint",
    "High-range linear point",
    "Upper boundary (non-saturating)",
    
    "",
    "Target CV < 10.0% at LOQ",
    "Target CV < 2.0% across range",
    
    "",
    "Target S/N ≥ 10:1 for LOQ",
    "Target S/N < 300 (No optical saturation)",
    
    "",
    "Unweighted first-order fit",
    "High detector sensitivity",
    "p > 0.05 (Not different from 0)",
    "Residual baseline noise variance",
    "Target R² ≥ 0.9990",
    
    "",
    "1/x² weighting factor",
    "Optimized for trace range",
    "Minimized intercept uncertainty",
    "Target R²_w ≥ 0.9990",
    
    "",
    "Unweighted mathematical limit",
    "Unweighted mathematical limit",
    "Weighted trace-refined LOD",
    "Weighted trace-refined LOQ",
    "Empirical operational LOQ floor"
  )
)

# 4. Print Master Table in RStudio Console
print(master_table, row.names = FALSE)

# 5. Optional: Save Master Table to CSV File
write.csv(master_table, "Master_HPLC_Validation_Table.csv", row.names = FALSE)


# ==============================================================================
# STEP 4: Complete Remaining Calculations & Generate Updated Master Table
# ==============================================================================

library(lmtest)

# 1. Dataset Setup
conc <- c(0.06, 0.06, 0.06, 0.10, 0.10, 0.10, 0.50, 0.50, 0.50, 
          1.00, 1.00, 1.00, 1.50, 1.50, 1.50, 2.00, 2.00, 2.00)

midostaurin <- c(1021.29, 936.42, 927.57, 1688.05, 1726.82, 1686.04, 
                 8450.38, 8426.75, 8453.75, 16763.79, 16983.58, 16873.69, 
                 25124.92, 25361.16, 25260.00, 33621.06, 33747.36, 33570.54)

tretinoin <- c(949.14, 906.42, 915.35, 1593.15, 1580.98, 1591.92, 
               7953.00, 7983.27, 7938.67, 15829.71, 16210.76, 15750.12, 
               23546.43, 23798.43, 23499.18, 31230.41, 31441.90, 31277.41)

df <- data.frame(Conc = conc, Midostaurin = midostaurin, Tretinoin = tretinoin)
df$Weight <- 1 / (df$Conc^2)

# 2. Fit Models
mido_ols_fit <- lm(Midostaurin ~ Conc, data = df)
tret_ols_fit <- lm(Tretinoin ~ Conc, data = df)

mido_wls_fit <- lm(Midostaurin ~ Conc, weights = Weight, data = df)
tret_wls_fit <- lm(Tretinoin ~ Conc, weights = Weight, data = df)

mido_ols <- summary(mido_ols_fit)
tret_ols <- summary(tret_ols_fit)

mido_wls <- summary(mido_wls_fit)
tret_wls <- summary(tret_wls_fit)
print(mido_wls)
print(tret_wls)
print(tret_ols)
# 3. Extract 95% Confidence Intervals
mido_ci <- confint(mido_ols_fit)
tret_ci <- confint(tret_ols_fit)

# 4. Generate Dynamically Updated Master Table
master_table_updated <- data.frame(
  Parameter_Category = c(
    "--- I. CHROMATOGRAPHIC BEHAVIOR ---",
    "Target Retention Time (Rt, min)",
    "Evaluated Concentration Range (ug/mL)",
    "Total Calibration Injections (N)",
    
    "--- II. EXPERIMENTAL PEAK AREA RESPONSES (Mean ± SD) ---",
    "Conc 0.06 ug/mL (LOQ Floor)",
    "Conc 0.10 ug/mL (Trace Point)",
    "Conc 0.50 ug/mL (Low-Mid Standard)",
    "Conc 1.00 ug/mL (Geometric Midpoint)",
    "Conc 1.50 ug/mL (Mid-High Standard)",
    "Conc 2.00 ug/mL (Upper Linear Limit)",
    
    "--- III. REPEATABILITY & PRECISION (CV%) ---",
    "Repeatability at LOQ (0.06 ug/mL)",
    "Repeatability Range Across Assay (0.10 - 2.00 ug/mL)",
    
    "--- IV. SIGNAL-TO-NOISE (S/N) METRICS ---",
    "Mean S/N at 0.06 ug/mL (LOQ Floor)",
    "Mean S/N at 2.00 ug/mL (Upper Limit)",
    
    "--- V. REGRESSION PARAMETERS: ORDINARY LEAST SQUARES (OLS) ---",
    "OLS Linear Equation",
    "Slope (S, Area per ug/mL) ± SE",
    "Intercept (a) ± SE",
    "95% CI of Intercept",
    "95% CI of Slope",
    "Residual Standard Error (Sy/x)",
    "Coefficient of Determination (R²)",
    
    "--- VI. REGRESSION PARAMETERS: WEIGHTED LEAST SQUARES (WLS 1/x²) ---",
    "WLS Linear Equation",
    "Weighted Slope (S_w) ± SE",
    "Weighted Intercept (a_w) ± SE",
    "Weighted Coefficient of Determination (R²_w)",
    
    "--- VII. DERIVED SENSITIVITY LIMITS (LOD & LOQ) ---",
    "OLS Calculated LOD (3.3 * SE_a / S)",
    "OLS Calculated LOQ (10 * SE_a / S)",
    "WLS Calculated LOD (3.3 * SE_aw / S_w)",
    "WLS Calculated LOQ (10 * SE_aw / S_w)",
    "Empirical Operational LOQ Floor (S/N ≥ 10:1)"
  ),
  
  Midostaurin = c(
    "", "3.31 ± 0.03", "0.06 - 2.00", "18 (6 levels × 3 reps)",
    "", "961.76 ± 51.74", "1700.30 ± 22.99", "8443.63 ± 14.71", 
    "16873.69 ± 109.90", "25248.69 ± 118.53", "33646.32 ± 91.08",
    "", "5.38%", "0.17% - 1.35%",
    "", "10.47", "119.43",
    "",
    sprintf("y = %.2fx + %.2f", mido_ols$coefficients["Conc","Estimate"], mido_ols$coefficients["(Intercept)","Estimate"]),
    sprintf("%.2f ± %.2f", mido_ols$coefficients["Conc","Estimate"], mido_ols$coefficients["Conc","Std. Error"]),
    sprintf("%.2f ± %.2f", mido_ols$coefficients["(Intercept)","Estimate"], mido_ols$coefficients["(Intercept)","Std. Error"]),
    sprintf("[%.2f to %.2f]", mido_ci["(Intercept)",1], mido_ci["(Intercept)",2]),
    sprintf("[%.2f to %.2f]", mido_ci["Conc",1], mido_ci["Conc",2]),
    sprintf("%.2f", mido_ols$sigma),
    sprintf("%.5f", mido_ols$r.squared),
    
    "",
    sprintf("y = %.2fx + %.2f", mido_wls$coefficients["Conc","Estimate"], mido_wls$coefficients["(Intercept)","Estimate"]),
    sprintf("%.2f ± %.2f", mido_wls$coefficients["Conc","Estimate"], mido_wls$coefficients["Conc","Std. Error"]),
    sprintf("%.2f ± %.2f", mido_wls$coefficients["(Intercept)","Estimate"], mido_wls$coefficients["(Intercept)","Std. Error"]),
    sprintf("%.5f", mido_wls$r.squared),
    
    "",
    sprintf("%.4f ug/mL", (3.3 * mido_ols$coefficients["(Intercept)","Std. Error"]) / mido_ols$coefficients["Conc","Estimate"]),
    sprintf("%.4f ug/mL", (10.0 * mido_ols$coefficients["(Intercept)","Std. Error"]) / mido_ols$coefficients["Conc","Estimate"]),
    sprintf("%.4f ug/mL", (3.3 * mido_wls$coefficients["(Intercept)","Std. Error"]) / mido_wls$coefficients["Conc","Estimate"]),
    sprintf("%.4f ug/mL", (10.0 * mido_wls$coefficients["(Intercept)","Std. Error"]) / mido_wls$coefficients["Conc","Estimate"]),
    "0.0600 ug/mL"
  ),
  
  Tretinoin = c(
    "", "9.21 ± 0.05", "0.06 - 2.00", "18 (6 levels × 3 reps)",
    "", "923.64 ± 22.53", "1588.68 ± 6.70", "7958.31 ± 22.77", 
    "15930.20 ± 246.21", "23614.68 ± 160.88", "31316.57 ± 111.05",
    "", "2.44%", "0.29% - 1.55%",
    "", "10.66", "111.29",
    "",
    sprintf("y = %.2fx + %.2f", tret_ols$coefficients["Conc","Estimate"], tret_ols$coefficients["(Intercept)","Estimate"]),
    sprintf("%.2f ± %.2f", tret_ols$coefficients["Conc","Estimate"], tret_ols$coefficients["Conc","Std. Error"]),
    sprintf("%.2f ± %.2f", tret_ols$coefficients["(Intercept)","Estimate"], tret_ols$coefficients["(Intercept)","Std. Error"]),
    sprintf("[%.2f to %.2f]", tret_ci["(Intercept)",1], tret_ci["(Intercept)",2]),
    sprintf("[%.2f to %.2f]", tret_ci["Conc",1], tret_ci["Conc",2]),
    sprintf("%.2f", tret_ols$sigma),
    sprintf("%.5f", tret_ols$r.squared),
    
    "",
    sprintf("y = %.2fx + %.2f", tret_wls$coefficients["Conc","Estimate"], tret_wls$coefficients["(Intercept)","Estimate"]),
    sprintf("%.2f ± %.2f", tret_wls$coefficients["Conc","Estimate"], tret_wls$coefficients["Conc","Std. Error"]),
    sprintf("%.2f ± %.2f", tret_wls$coefficients["(Intercept)","Estimate"], tret_wls$coefficients["(Intercept)","Std. Error"]),
    sprintf("%.5f", tret_wls$r.squared),
    
    "",
    sprintf("%.4f ug/mL", (3.3 * tret_ols$coefficients["(Intercept)","Std. Error"]) / tret_ols$coefficients["Conc","Estimate"]),
    sprintf("%.4f ug/mL", (10.0 * tret_ols$coefficients["(Intercept)","Std. Error"]) / tret_ols$coefficients["Conc","Estimate"]),
    sprintf("%.4f ug/mL", (3.3 * tret_wls$coefficients["(Intercept)","Std. Error"]) / tret_wls$coefficients["Conc","Estimate"]),
    sprintf("%.4f ug/mL", (10.0 * tret_wls$coefficients["(Intercept)","Std. Error"]) / tret_wls$coefficients["Conc","Estimate"]),
    "0.0600 ug/mL"
  )
)

# Print updated master table
print(master_table_updated, row.names = FALSE)
View(master_table_updated)

# ==============================================================================
# STEP 5: Advanced Linearity & Residual Assumption Diagnostics
# ==============================================================================

library(lmtest)

# 1. ANOVA Lack-of-Fit (LOF) Test
mido_factor <- lm(Midostaurin ~ as.factor(Conc), data = df)
tret_factor <- lm(Tretinoin ~ as.factor(Conc), data = df)

mido_lof <- anova(mido_ols_fit, mido_factor)
tret_lof <- anova(tret_ols_fit, tret_factor)

# 2. Mandel's Fitting Test (Linear vs. Quadratic)
mido_quad <- lm(Midostaurin ~ Conc + I(Conc^2), data = df)
tret_quad <- lm(Tretinoin ~ Conc + I(Conc^2), data = df)

mido_mandel <- anova(mido_ols_fit, mido_quad)
tret_mandel <- anova(tret_ols_fit, tret_quad)

# 3. Shapiro-Wilk Test for Normality of Residuals
mido_sw <- shapiro.test(residuals(mido_ols_fit))
tret_sw <- shapiro.test(residuals(tret_ols_fit))

# 4. Breusch-Pagan Test for Homoscedasticity
mido_bp <- bptest(mido_ols_fit)
tret_bp <- bptest(tret_ols_fit)

# Print Summary of Diagnostic Tests
cat("==================================================================\n")
cat("            DIAGNOSTIC HYPOTHESIS TEST RESULTS                   \n")
cat("==================================================================\n\n")

cat("--- MIDOSTAURIN ---\n")
cat(sprintf("ANOVA Lack-of-Fit    : F = %.3f, p-value = %.4f (Pass: p > 0.05)\n", mido_lof$F[2], mido_lof$`Pr(>F)`[2]))
cat(sprintf("Mandel's Fitting Test : F = %.3f, p-value = %.4f (Linear optimal: p > 0.05)\n", mido_mandel$F[2], mido_mandel$`Pr(>F)`[2]))
cat(sprintf("Shapiro-Wilk Normality: W = %.3f, p-value = %.4f (Normal: p > 0.05)\n", mido_sw$statistic, mido_sw$p.value))
cat(sprintf("Breusch-Pagan Homosc. : BP = %.3f, p-value = %.4f (Homoscedastic: p > 0.05)\n\n", mido_bp$statistic, mido_bp$p.value))

cat("--- TRETINOIN ---\n")
cat(sprintf("ANOVA Lack-of-Fit    : F = %.3f, p-value = %.4f (Pass: p > 0.05)\n", tret_lof$F[2], tret_lof$`Pr(>F)`[2]))
cat(sprintf("Mandel's Fitting Test : F = %.3f, p-value = %.4f\n", tret_mandel$F[2], tret_mandel$`Pr(>F)`[2]))
cat(sprintf("Shapiro-Wilk Normality: W = %.3f, p-value = %.4f\n", tret_sw$statistic, tret_sw$p.value))
cat(sprintf("Breusch-Pagan Homosc. : BP = %.3f, p-value = %.4f (Homoscedastic: p > 0.05)\n", tret_bp$statistic, tret_bp$p.value))

# ==============================================================================
# NATURE-STYLE PUBLICATION FIGURE: LINEARITY & RESIDUAL DIAGNOSTICS
# Analytes: Midostaurin & Tretinoin (0.06 - 2.00 ug/mL)
# ==============================================================================

# 1. Install & Load Required Packages
required_packages <- c("ggplot2", "gridExtra", "grid", "scales")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(ggplot2)
library(gridExtra)
library(grid)
library(scales)

# 2. Input Dataset (Triplicate Injections across 6 Concentrations)
conc <- c(0.06, 0.06, 0.06, 0.10, 0.10, 0.10, 0.50, 0.50, 0.50, 
          1.00, 1.00, 1.00, 1.50, 1.50, 1.50, 2.00, 2.00, 2.00)

midostaurin <- c(1021.29, 936.42, 927.57, 1688.05, 1726.82, 1686.04, 
                 8450.38, 8426.75, 8453.75, 16763.79, 16983.58, 16873.69, 
                 25124.92, 25361.16, 25260.00, 33621.06, 33747.36, 33570.54)

tretinoin <- c(949.14, 906.42, 915.35, 1593.15, 1580.98, 1591.92, 
               7953.00, 7983.27, 7938.67, 15829.71, 16210.76, 15750.12, 
               23546.43, 23798.43, 23499.18, 31230.41, 31441.90, 31277.41)

df <- data.frame(Conc = conc, Midostaurin = midostaurin, Tretinoin = tretinoin)

# 3. Fit Models and Calculate Residuals
mido_lm <- lm(Midostaurin ~ Conc, data = df)
tret_lm <- lm(Tretinoin ~ Conc, data = df)

df$Mido_Resid <- residuals(mido_lm)
df$Tret_Resid <- residuals(tret_lm)

# 4. Define Nature Publication Theme
theme_nature <- function() {
  theme_classic(base_size = 11, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0, color = "#111111", margin = margin(b = 4)),
      plot.subtitle = element_text(size = 9.5, color = "#555555", margin = margin(b = 8)),
      axis.title = element_text(face = "bold", size = 10, color = "#222222"),
      axis.text = element_text(size = 9, color = "#333333"),
      axis.line = element_line(color = "#222222", linewidth = 0.6),
      axis.ticks = element_line(color = "#222222", linewidth = 0.6),
      panel.grid.major = element_line(color = "#EAEAEA", linewidth = 0.4, linetype = "dashed"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 12, 10, 12)
    )
}

# Color Palettes
col_mido_point <- "#008080" # Nature Teal
col_mido_line  <- "#004D40"
col_mido_fill  <- "#80CBC4"

col_tret_point <- "#8E24AA" # Nature Deep Purple
col_tret_line  <- "#4A148C"
col_tret_fill  <- "#E1BEE7"

# 5. Panel A: Midostaurin Calibration Curve
p1 <- ggplot(df, aes(x = Conc, y = Midostaurin)) +
  geom_smooth(method = "lm", formula = y ~ x, fill = col_mido_fill, color = col_mido_line, 
              linewidth = 1.1, alpha = 0.35) +
  geom_point(shape = 21, fill = col_mido_point, color = "#00332C", size = 3, stroke = 0.8) +
  scale_y_continuous(labels = comma, limits = c(-500, 36000)) +
  scale_x_continuous(breaks = c(0.06, 0.10, 0.50, 1.00, 1.50, 2.00)) +
  annotate("text", x = 0.15, y = 32000, label = "y = 16832.56x + 3.07\nR² = 0.99996", 
           hjust = 0, size = 3.3, fontface = "bold", color = col_mido_line) +
  labs(title = "a  Midostaurin Calibration Profile",
       subtitle = "Linear response across 0.06–2.00 µg/mL (n = 18)",
       x = "Concentration (µg/mL)", y = "Peak Area (µV·sec)") +
  theme_nature()

# 6. Panel B: Tretinoin Calibration Curve
p2 <- ggplot(df, aes(x = Conc, y = Tretinoin)) +
  geom_smooth(method = "lm", formula = y ~ x, fill = col_tret_fill, color = col_tret_line, 
              linewidth = 1.1, alpha = 0.35) +
  geom_point(shape = 21, fill = col_tret_point, color = "#2A0052", size = 3, stroke = 0.8) +
  scale_y_continuous(labels = comma, limits = c(-500, 34000)) +
  scale_x_continuous(breaks = c(0.06, 0.10, 0.50, 1.00, 1.50, 2.00)) +
  annotate("text", x = 0.15, y = 30000, label = "y = 15680.50x + 70.12\nR² = 0.99983", 
           hjust = 0, size = 3.3, fontface = "bold", color = col_tret_line) +
  labs(title = "b  Tretinoin Calibration Profile",
       subtitle = "Linear response across 0.06–2.00 µg/mL (n = 18)",
       x = "Concentration (µg/mL)", y = "Peak Area (µV·sec)") +
  theme_nature()

# 7. Panel C: Midostaurin Residuals Plot
p3 <- ggplot(df, aes(x = Conc, y = Mido_Resid)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#D32F2F", linewidth = 0.8) +
  geom_point(shape = 21, fill = col_mido_point, color = "#00332C", size = 3, stroke = 0.8) +
  scale_x_continuous(breaks = c(0.06, 0.10, 0.50, 1.00, 1.50, 2.00)) +
  scale_y_continuous(limits = c(-180, 180)) +
  labs(title = "c  Midostaurin Residual Distribution",
       subtitle = "Homoscedastic random scatter around zero-residual line",
       x = "Concentration (µg/mL)", y = "Residual Area (µV·sec)") +
  theme_nature()

# 8. Panel D: Tretinoin Residuals Plot
p4 <- ggplot(df, aes(x = Conc, y = Tret_Resid)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#D32F2F", linewidth = 0.8) +
  geom_point(shape = 21, fill = col_tret_point, color = "#2A0052", size = 3, stroke = 0.8) +
  scale_x_continuous(breaks = c(0.06, 0.10, 0.50, 1.00, 1.50, 2.00)) +
  scale_y_continuous(limits = c(-250, 500)) +
  labs(title = "d  Tretinoin Residual Distribution",
       subtitle = "Variance evaluation across 33-fold concentration range",
       x = "Concentration (µg/mL)", y = "Residual Area (µV·sec)") +
  theme_nature()

# 9. Arrange All Panels and Export High-Resolution PNG & Vector PDF
nature_figure <- grid.arrange(p1, p2, p3, p4, ncol = 2)

# Save to your current working directory
ggsave("Figure_1_Nature_Style.pdf", plot = nature_figure, width = 9.5, height = 8.5, dpi = 600)
ggsave("Figure_1_Nature_Style.png", plot = nature_figure, width = 9.5, height = 8.5, dpi = 300)
