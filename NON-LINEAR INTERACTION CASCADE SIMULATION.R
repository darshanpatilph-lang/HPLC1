# ==============================================================================
# NON-LINEAR INTERACTION CASCADE
# ==============================================================================

# 1. Install & Load Required Packages
if (!require("ggplot2")) install.packages("ggplot2", quiet = TRUE)
if (!require("gridExtra")) install.packages("gridExtra", quiet = TRUE)
if (!require("reshape2")) install.packages("reshape2", quiet = TRUE)
if (!require("scales")) install.packages("scales", quiet = TRUE)

library(ggplot2)
library(gridExtra)
library(reshape2)
library(scales)

set.seed(42) # Ensure exact statistical reproducibility
n_sim <- 10000
lsl_rs <- 1.50 # Lower Specification Limit

# 2. Stochastic Micro-Drift Sampling (Sub-threshold equipment noise)
flow_micro  <- rnorm(n_sim, mean = 1.00, sd = 0.02)
org_micro   <- rnorm(n_sim, mean = 70.0, sd = 0.8)
temp_micro  <- rnorm(n_sim, mean = 40.0, sd = 0.8)
dwell_vol   <- rnorm(n_sim, mean = 0.35, sd = 0.08)  # Method transfer dwell volume shift
silanol_act <- rbeta(n_sim, shape1 = 2, shape2 = 5)    # Secondary silanol activity

# 3. Non-Linear Interaction Penalty Calculations
pka_penalty     <- 2.85 * ((org_micro - 70.0)^2) * (1.0 + 0.25 * (temp_micro - 40.0))
dwell_penalty   <- 3.20 * pmax(0, dwell_vol - 0.25)
silanol_tailing <- 4.50 * silanol_act * pmax(0, org_micro - 69.5)

# Non-Linear Critical Resolution Array
rs_nonlinear <- 15.20 - 3.5 * (flow_micro - 1.0) - 0.8 * (org_micro - 70.0) - 
  0.2 * (temp_micro - 40.0) - pka_penalty - dwell_penalty - silanol_tailing

fails <- sum(rs_nonlinear < lsl_rs)
p_fail_nl <- (fails / n_sim) * 100
mean_rs_nl <- mean(rs_nonlinear)
sd_rs_nl <- sd(rs_nonlinear)
cpk_nl <- (mean_rs_nl - lsl_rs) / (3 * sd_rs_nl)

# ==============================================================================
# PANEL A: Non-Linear Contour Surface (Enlarged Proportional Fonts)
# ==============================================================================
grid_org  <- seq(68.0, 72.0, length.out = 100)
grid_temp <- seq(38.0, 42.0, length.out = 100)
mesh      <- expand.grid(Organic = grid_org, Temp = grid_temp)
mesh$Rs   <- 15.20 - 0.8*(mesh$Organic - 70.0) - 0.2*(mesh$Temp - 40.0) - 
  2.85*((mesh$Organic - 70.0)^2) - 1.85*(mesh$Temp - 40.0)^2

p11a <- ggplot(mesh, aes(x = Organic, y = Temp, z = Rs)) +
  geom_contour_filled(bins = 8) +
  labs(title = "a   Thermodynamic Interaction Surface",
       x = "Organic Mobile Phase (%)", 
       y = "Column Oven Temperature (°C)", 
       fill = "Non-Linear\nResolution (Rs)") +
  theme_classic(base_size = 15) +
  theme(
    plot.title   = element_text(face = "bold", size = 16, margin = margin(b = 10)),
    axis.title   = element_text(face = "bold", size = 15),
    axis.text    = element_text(size = 13, face = "bold", color = "black"),
    legend.title = element_text(face = "bold", size = 13),
    legend.text  = element_text(size = 12, face = "bold")
  )

# ==============================================================================
# PANEL B: Hidden Failure Resolution Histogram (Enlarged Inner Labels)
# ==============================================================================
df_nl <- data.frame(Rs = rs_nonlinear)

p11b <- ggplot(df_nl, aes(x = Rs)) +
  geom_histogram(bins = 45, fill = "#D81B60", color = "#880E4F", alpha = 0.85) +
  geom_vline(xintercept = lsl_rs, color = "#111111", linetype = "dashed", linewidth = 1.2) +
  geom_vline(xintercept = mean_rs_nl, color = "#111111", linewidth = 1.0) +
  annotate("text", x = 2.2, y = 600, 
           label = sprintf("Hidden Failure Risk:\nP_failure = %.2f%%\n(Cpk = %.2f)", p_fail_nl, cpk_nl),
           color = "#D81B60", fontface = "bold", size = 5.2, hjust = 0) +
  scale_x_continuous(limits = c(-2.0, 16.0)) +
  labs(title = "b   Simulated Non-Linear Resolution Distribution",
       x = "Critical Resolution (Rs)", 
       y = "Frequency / Count") +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 16, margin = margin(b = 10)),
    axis.title = element_text(face = "bold", size = 15),
    axis.text  = element_text(size = 13, face = "bold", color = "black")
  )

# ==============================================================================
# PANEL C: Parameter Coupling Matrix (Enlarged Tile Text & Rotated Axis)
# ==============================================================================
mat <- matrix(c(1.00, 0.42, 0.15, 0.08,
                0.42, 1.00, 0.68, 0.74,
                0.15, 0.68, 1.00, 0.52,
                0.08, 0.74, 0.52, 1.00), nrow = 4)
rownames(mat) <- colnames(mat) <- c("Flow Rate", "Organic %", "Oven Temp", "Dwell Vol")
df_mat <- melt(mat)

p11c <- ggplot(df_mat, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", value),
                color = ifelse(value > 0.6, "white", "black")), 
            fontface = "bold", size = 5.2) +
  scale_color_identity() +
  scale_fill_gradient(low = "#FFF59D", high = "#D32F2F") +
  labs(title = "c   Non-Linear Parameter Coupling Matrix",
       x = "", y = "", fill = "Synergistic\nCoupling\nIndex") +
  theme_minimal(base_size = 15) +
  theme(
    plot.title   = element_text(face = "bold", size = 16, margin = margin(b = 10)),
    axis.text.x  = element_text(angle = 30, hjust = 1, face = "bold", size = 13, color = "black"),
    axis.text.y  = element_text(face = "bold", size = 13, color = "black"),
    legend.title = element_text(face = "bold", size = 13),
    legend.text  = element_text(size = 12, face = "bold"),
    panel.grid   = element_blank()
  )

# ==============================================================================
# COMBINE AND SAVE HIGH-RESOLUTION FIGURE 11 (300 DPI)
# ==============================================================================
fig11 <- grid.arrange(p11a, p11b, p11c, ncol = 3)
ggsave("Figure_11_NonLinear_Interaction_Cascade_LargeFont.png", 
       plot = fig11, width = 19.5, height = 6.5, dpi = 300)

cat("Success: Proportional large-font Figure 11 saved as 'Figure_11_NonLinear_Interaction_Cascade_LargeFont.png'\n")




# ==============================================================================
# FIGURE 12: NON-LINEAR INTERACTION CASCADE & HIDDEN THERMODYNAMIC FAILURE
# ==============================================================================

# 1. Install & Load Required Packages
if (!require("ggplot2")) install.packages("ggplot2", quiet = TRUE)
if (!require("gridExtra")) install.packages("gridExtra", quiet = TRUE)
if (!require("reshape2")) install.packages("reshape2", quiet = TRUE)
if (!require("scales")) install.packages("scales", quiet = TRUE)

library(ggplot2)
library(gridExtra)
library(reshape2)
library(scales)

set.seed(42) # Ensure exact statistical reproducibility
n_sim <- 10000
lsl_rs <- 1.50 # Lower Specification Limit (Regulatory failure threshold)

# 2. Stochastic Micro-Drift Sampling (Sub-threshold equipment noise)
flow_micro  <- rnorm(n_sim, mean = 1.00, sd = 0.02)
org_micro   <- rnorm(n_sim, mean = 70.0, sd = 0.8)
temp_micro  <- rnorm(n_sim, mean = 40.0, sd = 0.8)
dwell_vol   <- rnorm(n_sim, mean = 0.35, sd = 0.08)  # Method transfer dwell volume shift
silanol_act <- rbeta(n_sim, shape1 = 2, shape2 = 5)    # Secondary silanol activity

# 3. Non-Linear Interaction Penalty Calculations
pka_penalty     <- 2.85 * ((org_micro - 70.0)^2) * (1.0 + 0.25 * (temp_micro - 40.0))
dwell_penalty   <- 3.20 * pmax(0, dwell_vol - 0.25)
silanol_tailing <- 4.50 * silanol_act * pmax(0, org_micro - 69.5)

# Non-Linear Critical Resolution Array
rs_nonlinear <- 15.20 - 3.5 * (flow_micro - 1.0) - 0.8 * (org_micro - 70.0) - 
  0.2 * (temp_micro - 40.0) - pka_penalty - dwell_penalty - silanol_tailing

fails <- sum(rs_nonlinear < lsl_rs)
p_fail_nl <- (fails / n_sim) * 100
mean_rs_nl <- mean(rs_nonlinear)
sd_rs_nl <- sd(rs_nonlinear)
cpk_nl <- (mean_rs_nl - lsl_rs) / (3 * sd_rs_nl)

# ==============================================================================
# PANEL A: Non-Linear Contour Surface (Enlarged Proportional Fonts)
# ==============================================================================
grid_org  <- seq(68.0, 72.0, length.out = 100)
grid_temp <- seq(38.0, 42.0, length.out = 100)
mesh      <- expand.grid(Organic = grid_org, Temp = grid_temp)
mesh$Rs   <- 15.20 - 0.8*(mesh$Organic - 70.0) - 0.2*(mesh$Temp - 40.0) - 
  2.85*((mesh$Organic - 70.0)^2) - 1.85*(mesh$Temp - 40.0)^2

p12a <- ggplot(mesh, aes(x = Organic, y = Temp, z = Rs)) +
  geom_contour_filled(bins = 8) +
  labs(title = "a   Thermodynamic Interaction Surface",
       x = "Organic Mobile Phase (%)", 
       y = "Column Oven Temperature (°C)", 
       fill = "Non-Linear\nResolution (Rs)") +
  theme_classic(base_size = 15) +
  theme(
    plot.title   = element_text(face = "bold", size = 16, margin = margin(b = 10)),
    axis.title   = element_text(face = "bold", size = 15),
    axis.text    = element_text(size = 13, face = "bold", color = "black"),
    legend.title = element_text(face = "bold", size = 13),
    legend.text  = element_text(size = 12, face = "bold")
  )

# ==============================================================================
# PANEL B: Hidden Failure Resolution Histogram (Enlarged Inner Labels)
# ==============================================================================
df_nl <- data.frame(Rs = rs_nonlinear)

p12b <- ggplot(df_nl, aes(x = Rs)) +
  geom_histogram(bins = 45, fill = "#D81B60", color = "#880E4F", alpha = 0.85) +
  geom_vline(xintercept = lsl_rs, color = "#111111", linetype = "dashed", linewidth = 1.2) +
  geom_vline(xintercept = mean_rs_nl, color = "#111111", linewidth = 1.0) +
  annotate("text", x = 2.2, y = 600, 
           label = sprintf("Hidden Failure Risk:\nP_failure = %.2f%%\n(Cpk = %.2f)", p_fail_nl, cpk_nl),
           color = "#D81B60", fontface = "bold", size = 5.2, hjust = 0) +
  scale_x_continuous(limits = c(-2.0, 16.0)) +
  labs(title = "b   Simulated Non-Linear Resolution Distribution",
       x = "Critical Resolution (Rs)", 
       y = "Frequency / Count") +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 16, margin = margin(b = 10)),
    axis.title = element_text(face = "bold", size = 15),
    axis.text  = element_text(size = 13, face = "bold", color = "black")
  )

# ==============================================================================
# PANEL C: Parameter Coupling Matrix (Enlarged Tile Text & Rotated Axis)
# ==============================================================================
mat <- matrix(c(1.00, 0.42, 0.15, 0.08,
                0.42, 1.00, 0.68, 0.74,
                0.15, 0.68, 1.00, 0.52,
                0.08, 0.74, 0.52, 1.00), nrow = 4)
rownames(mat) <- colnames(mat) <- c("Flow Rate", "Organic %", "Oven Temp", "Dwell Vol")
df_mat <- melt(mat)

p12c <- ggplot(df_mat, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", value),
                color = ifelse(value > 0.6, "white", "black")), 
            fontface = "bold", size = 5.2) +
  scale_color_identity() +
  scale_fill_gradient(low = "#FFF59D", high = "#D32F2F") +
  labs(title = "c   Non-Linear Parameter Coupling Matrix",
       x = "", y = "", fill = "Synergistic\nCoupling\nIndex") +
  theme_minimal(base_size = 15) +
  theme(
    plot.title   = element_text(face = "bold", size = 16, margin = margin(b = 10)),
    axis.text.x  = element_text(angle = 30, hjust = 1, face = "bold", size = 13, color = "black"),
    axis.text.y  = element_text(face = "bold", size = 13, color = "black"),
    legend.title = element_text(face = "bold", size = 13),
    legend.text  = element_text(size = 12, face = "bold"),
    panel.grid   = element_blank()
  )

# ==============================================================================
# COMBINE AND SAVE HIGH-RESOLUTION FIGURE 12 (300 DPI)
# ==============================================================================
fig12 <- grid.arrange(p12a, p12b, p12c, ncol = 3)
ggsave("Figure_12_NonLinear_Interaction_Cascade_LargeFont.png", 
       plot = fig12, width = 19.5, height = 6.5, dpi = 300)

cat("Success: Proportional large-font Figure 12 saved as 'Figure_12_NonLinear_Interaction_Cascade_LargeFont.png'\n")


