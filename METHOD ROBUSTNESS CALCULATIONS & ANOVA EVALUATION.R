# ==============================================================================
# SECTION 3.6: METHOD ROBUSTNESS CALCULATIONS & ANOVA EVALUATION IN R
# ==============================================================================

# 1. Create Robustness Evaluation Dataset
df_rob <- data.frame(
  Parameter = factor(rep(c("Flow Rate (±0.1 mL/min)", 
                           "Organic Ratio (±2%)", 
                           "Column Temp (±2°C)", 
                           "Wavelength (±2 nm)"), each = 3)),
  Variation = factor(c("0.9 mL/min", "1.0 mL/min (Nominal)", "1.1 mL/min",
                       "78% Organic", "80% Organic (Nominal)", "82% Organic",
                       "38°C", "40°C (Nominal)", "42°C",
                       "308 nm", "310 nm (Nominal)", "312 nm")),
  Mido_Rec = c(100.15, 100.23, 99.88, 99.72, 100.23, 100.41, 100.08, 100.23, 99.92, 100.31, 100.23, 99.85),
  Mido_RSD = c(0.58, 0.65, 0.71, 0.82, 0.65, 0.74, 0.61, 0.65, 0.69, 0.52, 0.65, 0.60),
  Tret_Rec = c(100.82, 101.15, 100.45, 100.25, 101.15, 101.32, 100.95, 101.15, 100.62, 101.02, 101.15, 100.78),
  Tret_RSD = c(1.25, 1.55, 1.38, 1.48, 1.55, 1.62, 1.31, 1.55, 1.40, 1.28, 1.55, 1.35)
)

# 2. Perform ANOVA Across Parameter Variations
mido_anova <- aov(Mido_Rec ~ Parameter, data = df_rob)
tret_anova <- aov(Tret_Rec ~ Parameter, data = df_rob)

cat("=== MIDOSTAURIN ROBUSTNESS ANOVA ===\n")
print(summary(mido_anova))

cat("\n=== TRETINOIN ROBUSTNESS ANOVA ===\n")
print(summary(tret_anova))


# ==============================================================================
# SECTION 3.6: REMAINING ROBUSTNESS CALCULATIONS & YOUDEN EFFECT MAGNITUDES
# ==============================================================================

# 1. Input Robustness Recovery Data
mido_rec <- c(100.15, 100.23, 99.88, 99.72, 100.23, 100.41, 100.08, 100.23, 99.92, 100.31, 100.23, 99.85)
tret_rec <- c(100.82, 101.15, 100.45, 100.25, 101.15, 101.32, 100.95, 101.15, 100.62, 101.02, 101.15, 100.78)

# 2. Overall Summary Statistics Across All 12 Perturbation Conditions
mido_overall_mean <- mean(mido_rec)
mido_overall_sd   <- sd(mido_rec)
mido_overall_rsd  <- (mido_overall_sd / mido_overall_mean) * 100

tret_overall_mean <- mean(tret_rec)
tret_overall_sd   <- sd(tret_rec)
tret_overall_rsd  <- (tret_overall_sd / tret_overall_mean) * 100

cat("=== OVERALL ROBUSTNESS SUMMARY METRICS ===\n")
cat("Midostaurin Overall Recovery:", round(mido_overall_mean, 2), "% ±", round(mido_overall_sd, 2), "%\n")
cat("Midostaurin Overall % RSD:  ", round(mido_overall_rsd, 2), "%\n\n")

cat("Tretinoin Overall Recovery:  ", round(tret_overall_mean, 2), "% ±", round(tret_overall_sd, 2), "%\n")
cat("Tretinoin Overall % RSD:    ", round(tret_overall_rsd, 2), "%\n\n")

# 3. Youden's Parameter Effect Magnitudes (Delta E = High Level - Low Level)
# Flow Rate Effect (1.1 mL/min - 0.9 mL/min)
mido_eff_flow <- 99.88 - 100.15
tret_eff_flow <- 100.45 - 100.82

# Organic Ratio Effect (82% Organic - 78% Organic)
mido_eff_org <- 100.41 - 99.72
tret_eff_org <- 101.32 - 100.25

# Column Temperature Effect (42°C - 38°C)
mido_eff_temp <- 99.92 - 100.08
tret_eff_temp <- 100.62 - 100.95

# Detection Wavelength Effect (312 nm - 308 nm)
mido_eff_wave <- 99.85 - 100.31
tret_eff_wave <- 100.78 - 101.02

cat("=== YOUDEN PARAMETER EFFECT MAGNITUDES (Delta E) ===\n")
cat("Flow Rate Effect (Delta E):     Midostaurin =", round(mido_eff_flow, 2), "%, Tretinoin =", round(tret_eff_flow, 2), "%\n")
cat("Organic Ratio Effect (Delta E): Midostaurin =", round(mido_eff_org, 2), "%, Tretinoin =", round(tret_eff_org, 2), "%\n")
cat("Column Temp Effect (Delta E):    Midostaurin =", round(mido_eff_temp, 2), "%, Tretinoin =", round(tret_eff_temp, 2), "%\n")
cat("Wavelength Effect (Delta E):    Midostaurin =", round(mido_eff_wave, 2), "%, Tretinoin =", round(tret_eff_wave, 2), "%\n")
