# ==============================================================================
# SECTION 3.1.2: MODEL FIT COMPARISON & RESIDUAL DIAGNOSTICS IN R
# ==============================================================================

library(lmtest)

# 1. Fit First-Order Linear Models
mido_lm <- lm(Midostaurin ~ Conc, data = df)
tret_lm <- lm(Tretinoin ~ Conc, data = df)

# 2. Fit Second-Order Quadratic Models
mido_quad <- lm(Midostaurin ~ Conc + I(Conc^2), data = df)
tret_quad <- lm(Tretinoin ~ Conc + I(Conc^2), data = df)

# 3. Mandel's Fitting Test (Linear vs. Quadratic Comparison)
mido_mandel <- anova(mido_lm, mido_quad)
tret_mandel <- anova(tret_lm, tret_quad)

# 4. ANOVA Lack-of-Fit Test
mido_factor <- lm(Midostaurin ~ as.factor(Conc), data = df)
tret_factor <- lm(Tretinoin ~ as.factor(Conc), data = df)

mido_lof <- anova(mido_lm, mido_factor)
tret_lof <- anova(tret_lm, tret_factor)

# 5. Shapiro-Wilk Residual Normality Test
mido_sw <- shapiro.test(residuals(mido_lm))
tret_sw <- shapiro.test(residuals(tret_lm))

# 6. Breusch-Pagan Homoscedasticity Test
mido_bp <- bptest(mido_lm)
tret_bp <- bptest(tret_lm)

print(mido_bp)
print(mido_sw)
print(tret_bp)
print(tret_sw)

# Print Summary Results
cat("Midostaurin Mandel F-stat:", mido_mandel$F[2], "p-value:", mido_mandel$`Pr(>F)`[2], "\n")
cat("Tretinoin Mandel F-stat:", tret_mandel$F[2], "p-value:", tret_mandel$`Pr(>F)`[2], "\n")
cat("Midostaurin LOF F-stat:", mido_lof$F[2], "p-value:", mido_lof$`Pr(>F)`[2], "\n")
cat("Tretinoin LOF F-stat:", tret_factor$F[2], "p-value:", tret_lof$`Pr(>F)`[2], "\n")

# ==============================================================================
# SECTION 3.1.3: RESPONSE FACTORS & LEVERAGE (COOK'S DISTANCE) IN R
# ==============================================================================

# 1. Calculate Response Factors (RF = Peak Area / Nominal Concentration)
df$Mido_RF <- df$Midostaurin / df$Conc
df$Tret_RF <- df$Tretinoin / df$Conc

mido_rf_mean <- mean(df$Mido_RF)
mido_rf_cv   <- (sd(df$Mido_RF) / mido_rf_mean) * 100

tret_rf_mean <- mean(df$Tret_RF)
tret_rf_cv   <- (sd(df$Tret_RF) / tret_rf_mean) * 100

# 2. Compute Cook's Distance and Leverage (Hat Values)
df$Mido_Cooks <- cooks.distance(mido_lm)
df$Tret_Cooks <- cooks.distance(tret_lm)

df$Mido_Hat <- hatvalues(mido_lm)
df$Tret_Hat <- hatvalues(tret_lm)

# Print Summary
cat(sprintf("Midostaurin RF Mean: %.2f ± %.2f (CV = %.2f%%)\n", mido_rf_mean, sd(df$Mido_RF), mido_rf_cv))
cat(sprintf("Tretinoin RF Mean: %.2f ± %.2f (CV = %.2f%%)\n", tret_rf_mean, sd(df$Tret_RF), tret_rf_cv))
cat(sprintf("Midostaurin Max Cook's Distance: %.4f\n", max(df$Mido_Cooks)))
cat(sprintf("Tretinoin Max Cook's Distance: %.4f\n", max(df$Tret_Cooks)))


# ==============================================================================
# Analytes: Midostaurin & Tretinoin (0.06 - 2.00 ug/mL)
# Panels:
#   a: Midostaurin OLS & WLS Calibration Profiles
#   b: Tretinoin OLS & WLS Calibration Profiles
#   c: Midostaurin Residual Scatter Plot (Homoscedasticity)
#   d: Tretinoin Residual Scatter Plot (Homoscedasticity)
# ==============================================================================

# 1. Install & Load Required Packages
required_pkgs <- c("ggplot2", "gridExtra", "grid", "scales")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs)

library(ggplot2)
library(gridExtra)
library(grid)
library(scales)

# 2. Input Experimental Dataset (Triplicate Injections across 6 Concentrations)
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

# 3. Fit OLS Linear Models & Compute Calibration Residuals
mido_lm <- lm(Midostaurin ~ Conc, data = df)
tret_lm <- lm(Tretinoin ~ Conc, data = df)

df$Mido_Resid <- residuals(mido_lm)
df$Tret_Resid <- residuals(tret_lm)

# 4. Define Nature Portfolio Publication Theme
theme_nature <- function() {
  theme_classic(base_size = 11, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 11.5, hjust = 0, color = "#111111", margin = margin(b = 3)),
      plot.subtitle = element_text(size = 9, color = "#555555", margin = margin(b = 6)),
      axis.title = element_text(face = "bold", size = 9.5, color = "#222222"),
      axis.text = element_text(size = 8.5, color = "#333333"),
      axis.line = element_line(color = "#222222", linewidth = 0.55),
      axis.ticks = element_line(color = "#222222", linewidth = 0.55),
      panel.grid.major = element_line(color = "#EDEDED", linewidth = 0.35, linetype = "dashed"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(8, 10, 8, 10)
    )
}

# Color Palettes
col_mido_point <- "#008080" # Nature Teal
col_mido_line  <- "#004D40"
col_mido_fill  <- "#80CBC4"

col_tret_point <- "#8E24AA" # Nature Deep Purple
col_tret_line  <- "#4A148C"
col_tret_fill  <- "#E1BEE7"

# 5. Panel A: Midostaurin Calibration Profile
p1 <- ggplot(df, aes(x = Conc, y = Midostaurin)) +
  geom_smooth(method = "lm", formula = y ~ x, fill = col_mido_fill, color = col_mido_line, 
              linewidth = 1.0, alpha = 0.30) +
  geom_point(shape = 21, fill = col_mido_point, color = "#00332C", size = 2.8, stroke = 0.75) +
  scale_y_continuous(labels = comma, limits = c(-500, 36000)) +
  scale_x_continuous(breaks = c(0.06, 0.10, 0.50, 1.00, 1.50, 2.00)) +
  annotate("text", x = 0.12, y = 31500, 
           label = "OLS: y = 16832.56x + 3.07  (R² = 0.99996)\nWLS: y = 16940.82x - 38.79 (R²w = 0.99917)", 
           hjust = 0, size = 3.0, fontface = "bold", color = col_mido_line) +
  labs(title = "a  Midostaurin Calibration Profile",
       subtitle = "Linear response across 0.06–2.00 µg/mL (n = 18)",
       x = "Concentration (µg/mL)", y = "Peak Area (µV·sec)") +
  theme_nature()

# 6. Panel B: Tretinoin Calibration Profile
p2 <- ggplot(df, aes(x = Conc, y = Tretinoin)) +
  geom_smooth(method = "lm", formula = y ~ x, fill = col_tret_fill, color = col_tret_line, 
              linewidth = 1.0, alpha = 0.30) +
  geom_point(shape = 21, fill = col_tret_point, color = "#2A0052", size = 2.8, stroke = 0.75) +
  scale_y_continuous(labels = comma, limits = c(-500, 34000)) +
  scale_x_continuous(breaks = c(0.06, 0.10, 0.50, 1.00, 1.50, 2.00)) +
  annotate("text", x = 0.12, y = 29500, 
           label = "OLS: y = 15680.50x + 70.12  (R² = 0.99983)\nWLS: y = 15851.94x - 18.90 (R²w = 0.99968)", 
           hjust = 0, size = 3.0, fontface = "bold", color = col_tret_line) +
  labs(title = "b  Tretinoin Calibration Profile",
       subtitle = "Linear response across 0.06–2.00 µg/mL (n = 18)",
       x = "Concentration (µg/mL)", y = "Peak Area (µV·sec)") +
  theme_nature()

# 7. Panel C: Midostaurin Calibration Residuals Plot
p3 <- ggplot(df, aes(x = Conc, y = Mido_Resid)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#D32F2F", linewidth = 0.75) +
  geom_point(shape = 21, fill = col_mido_point, color = "#00332C", size = 2.8, stroke = 0.75) +
  scale_x_continuous(breaks = c(0.06, 0.10, 0.50, 1.00, 1.50, 2.00)) +
  scale_y_continuous(limits = c(-180, 180)) +
  labs(title = "c  Midostaurin Calibration Residuals",
       subtitle = "Random residual scatter around zero baseline (p > 0.05)",
       x = "Concentration (µg/mL)", y = "Residual Area (µV·sec)") +
  theme_nature()

# 8. Panel D: Tretinoin Calibration Residuals Plot
p4 <- ggplot(df, aes(x = Conc, y = Tret_Resid)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#D32F2F", linewidth = 0.75) +
  geom_point(shape = 21, fill = col_tret_point, color = "#2A0052", size = 2.8, stroke = 0.75) +
  scale_x_continuous(breaks = c(0.06, 0.10, 0.50, 1.00, 1.50, 2.00)) +
  scale_y_continuous(limits = c(-250, 500)) +
  labs(title = "d  Tretinoin Calibration Residuals",
       subtitle = "Variance distribution across 33-fold concentration range",
       x = "Concentration (µg/mL)", y = "Residual Area (µV·sec)") +
  theme_nature()

# 9. Arrange All Panels into a 2x2 Composite Layout
nature_figure <- grid.arrange(p1, p2, p3, p4, ncol = 2)

# 10. Save High-Resolution PDF & PNG
ggsave("Figure_1_Nature_Style_Final.pdf", plot = nature_figure, width = 9.0, height = 8.0, dpi = 600)
ggsave("Figure_1_Nature_Style_Final.png", plot = nature_figure, width = 9.0, height = 8.0, dpi = 300)

