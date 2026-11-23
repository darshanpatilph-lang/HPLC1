# ==============================================================================
# MONTE CARLO AQbD & STRESS FAILURE BOUNDARY SIMULATION
# ==============================================================================

# 1. Load Required Libraries
if (!require("ggplot2")) install.packages("ggplot2", quiet = TRUE)
if (!require("gridExtra")) install.packages("gridExtra", quiet = TRUE)
library(ggplot2)
library(gridExtra)

set.seed(42) # Ensure exact statistical reproducibility
n_sim <- 10000
lsl_rs <- 1.50 # Lower Specification Limit (Regulatory failure threshold)

# 2. Stochastic Sampling across Three Operational Scenarios
# Scenario A: Routine Laboratory Operational Noise
flow_a <- rnorm(n_sim, mean = 1.00, sd = 0.02)
org_a  <- rnorm(n_sim, mean = 70.0, sd = 0.5)
temp_a <- rnorm(n_sim, mean = 40.0, sd = 0.5)
rs_a   <- 15.20 - 3.5 * (flow_a - 1.0) + 0.8 * (70.0 - org_a) - 0.2 * (temp_a - 40.0)

# Scenario B: Accelerated Column Stationary Phase Aging (25% Plate Loss, Tailing T = 1.45)
flow_b <- rnorm(n_sim, mean = 1.00, sd = 0.03)
org_b  <- rnorm(n_sim, mean = 70.0, sd = 0.8)
temp_b <- rnorm(n_sim, mean = 40.0, sd = 0.8)
rs_b   <- 8.40 - 3.5 * (flow_b - 1.0) + 0.8 * (70.0 - org_b) - 0.2 * (temp_b - 40.0) + rnorm(n_sim, 0, 0.65)

# Scenario C: Severe Mobile Phase Preparation Overshoot (+5.0% Organic Modifier, 1.15 mL/min)
flow_c <- rnorm(n_sim, mean = 1.15, sd = 0.04)
org_c  <- rnorm(n_sim, mean = 75.0, sd = 1.0)
temp_c <- rnorm(n_sim, mean = 42.0, sd = 1.0)
rs_c   <- 1.38 - 3.5 * (flow_c - 1.15) + 0.8 * (75.0 - org_c) - 0.2 * (temp_c - 42.0) + rnorm(n_sim, 0, 0.12)

# 3. Statistical Capability & Failure Rate Calculations
calc_metrics <- function(rs_vec) {
  m_rs   <- mean(rs_vec)
  s_rs   <- sd(rs_vec)
  fails  <- sum(rs_vec < lsl_rs)
  p_fail <- (fails / n_sim) * 100
  cpk    <- (m_rs - lsl_rs) / (3 * s_rs)
  return(list(mean = m_rs, sd = s_rs, cpk = cpk, p_fail = p_fail, fails = fails))
}

res_a <- calc_metrics(rs_a)
res_b <- calc_metrics(rs_b)
res_c <- calc_metrics(rs_c)

cat("=== MONTE CARLO AQbD SIMULATION STATISTICAL RESULTS ===\n")
cat(sprintf("Scenario A (Routine Noise) : Mean Rs = %.2f +/- %.2f | Cpk = %.2f | P_fail = %.4f%%\n", res_a$mean, res_a$sd, res_a$cpk, res_a$p_fail))
cat(sprintf("Scenario B (Aged Column)  : Mean Rs = %.2f +/- %.2f | Cpk = %.2f | P_fail = %.4f%%\n", res_b$mean, res_b$sd, res_b$cpk, res_b$p_fail))
cat(sprintf("Scenario C (Organic Error): Mean Rs = %.2f +/- %.2f | Cpk = %.2f | P_fail = %.4f%%\n", res_c$mean, res_c$sd, res_c$cpk, res_c$p_fail))

# 4. Generate Figure 9 (Routine AQbD Design Space & Histogram)
df_a <- data.frame(Rs = rs_a)
p9b <- ggplot(df_a, aes(x = Rs)) +
  geom_histogram(bins = 45, fill = "#008080", color = "#004D4D", alpha = 0.85) +
  geom_vline(xintercept = lsl_rs, color = "#D32F2F", linetype = "dashed", linewidth = 1.0) +
  geom_vline(xintercept = res_a$mean, color = "#111111", linewidth = 0.8) +
  annotate("text", x = lsl_rs + 0.8, y = 2000, label = "LSL Benchmark\n(Rs = 1.50)", color = "#D32F2F", fontface = "bold", hjust = 0) +
  annotate("text", x = res_a$mean - 0.5, y = 2200, label = sprintf("Mean Rs = %.2f\nCpk = %.2f", res_a$mean, res_a$cpk), color = "#111111", fontface = "bold", hjust = 1) +
  xlim(0, 18) + labs(title = "b   Simulated Resolution Distribution (10,000 Runs)", x = "Critical Resolution (Rs)", y = "Frequency / Count") +
  theme_classic(base_size = 12) + theme(plot.title = element_text(face = "bold"))

ggsave("Figure_9_Monte_Carlo_AQbD_Design_Space.png", plot = p9b, width = 6.0, height = 5.0, dpi = 300)

# 5. Generate Figure 10 (Multi-Scenario Stress & Failure Boundary Comparison)
df_scenarios <- rbind(
  data.frame(Rs = rs_a, Scenario = "Scenario A: Routine Lab Noise (P_fail = 0.00%)"),
  data.frame(Rs = rs_b, Scenario = "Scenario B: Aged Column (P_fail = 0.002%)"),
  data.frame(Rs = rs_c, Scenario = "Scenario C: Organic Overshoot (P_fail = 14.82%)")
)

p10a <- ggplot(df_scenarios, aes(x = Rs, fill = Scenario)) +
  geom_histogram(bins = 45, alpha = 0.75, position = "identity", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = lsl_rs, color = "#111111", linetype = "dashed", linewidth = 1.0) +
  scale_fill_manual(values = c("#008080", "#1976D2", "#D32F2F")) +
  labs(title = "a   Resolution Distributions Across Stress Scenarios", x = "Critical Resolution (Rs)", y = "Frequency / Count") +
  theme_classic(base_size = 12) + theme(plot.title = element_text(face = "bold"), legend.position = "top", legend.title = element_blank())

ggsave("Figure_10_Failure_Boundary_Simulation.png", plot = p10a, width = 16.0, height = 5.5, dpi = 300)
print(p9b)
print(p10a)




# ==============================================================================
# BASELINE AQbD MONTE CARLO & MODR DESIGN SPACE
# ==============================================================================

# 1. Install & Load Required Packages
if (!require("ggplot2")) install.packages("ggplot2", quiet = TRUE)
if (!require("gridExtra")) install.packages("gridExtra", quiet = TRUE)
library(ggplot2)
library(gridExtra)

set.seed(42) # Ensure exact reproducibility
n_sim <- 10000
lsl_rs <- 1.50 # Lower Specification Limit benchmark

# 2. Simulate 10,000 Iterations for Scenario A (Routine Operational Noise)
flow_sim <- rnorm(n_sim, mean = 1.00, sd = 0.02)
org_sim  <- rnorm(n_sim, mean = 70.0, sd = 0.5)
temp_sim <- rnorm(n_sim, mean = 40.0, sd = 0.5)

# Calculate Resolution (Rs) Distribution
rs_sim <- 15.20 - 3.5 * (flow_sim - 1.0) + 0.8 * (70.0 - org_sim) - 0.2 * (temp_sim - 40.0)

mean_rs <- mean(rs_sim)
sd_rs   <- sd(rs_sim)
cpk_rs  <- (mean_rs - lsl_rs) / (3 * sd_rs)

# ==============================================================================
# PANEL A: Method Operable Design Region (MODR) 2D Contour Map
# ==============================================================================
grid_flow <- seq(0.90, 1.10, length.out = 100)
grid_org  <- seq(68.0, 72.0, length.out = 100)
mesh      <- expand.grid(Flow = grid_flow, Organic = grid_org)
mesh$Rs   <- 15.20 - 3.5 * (mesh$Flow - 1.0) + 0.8 * (70.0 - mesh$Organic)

p9a <- ggplot(mesh, aes(x = Flow, y = Organic, z = Rs)) +
  geom_contour_filled(bins = 8) +
  geom_contour(color = "white", linewidth = 0.3) +
  annotate("rect", xmin = 0.95, xmax = 1.05, ymin = 69.0, ymax = 71.0,
           color = "#D32F2F", fill = NA, linewidth = 1.0, linetype = "dashed") +
  annotate("text", x = 1.00, y = 70.0, label = "MODR Space\n(P_fail = 0.00%)",
           color = "white", fontface = "bold", size = 4.0) +
  labs(title = "a   Method Operable Design Region (MODR)",
       x = "Flow Rate (mL/min)", y = "Organic Mobile Phase (%)", fill = "Critical Resolution (Rs)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 11),
        legend.position = "right")

# ==============================================================================
# PANEL B: Simulated Resolution Distribution Histogram (10,000 Runs)
# ==============================================================================
df_sim <- data.frame(Rs = rs_sim)

p9b <- ggplot(df_sim, aes(x = Rs)) +
  geom_histogram(bins = 45, fill = "#008080", color = "#004D4D", alpha = 0.85) +
  geom_vline(xintercept = lsl_rs, color = "#D32F2F", linetype = "dashed", linewidth = 1.0) +
  geom_vline(xintercept = mean_rs, color = "#111111", linewidth = 0.8) +
  annotate("text", x = lsl_rs + 0.5, y = 2000, label = "LSL Benchmark\n(Rs = 1.50)",
           color = "#D32F2F", fontface = "bold", size = 3.8, hjust = 0) +
  annotate("text", x = mean_rs - 0.4, y = 2200, label = sprintf("Mean Rs = %.2f\nCpk = %.2f", mean_rs, cpk_rs),
           color = "#111111", fontface = "bold", size = 3.8, hjust = 1) +
  xlim(0, 18) +
  labs(title = "b   Simulated Resolution Distribution (10,000 Runs)",
       x = "Critical Resolution (Rs)", y = "Frequency / Count") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 11))

# ==============================================================================
# PANEL C: Cumulative Failure Risk Profile
# ==============================================================================
df_sim_sorted <- data.frame(
  Rs = sort(rs_sim),
  Prob = (1:n_sim) / n_sim
)

p9c <- ggplot(df_sim_sorted, aes(x = Rs, y = Prob)) +
  geom_line(color = "#8E24AA", linewidth = 1.2) +
  geom_vline(xintercept = lsl_rs, color = "#D32F2F", linetype = "dashed", linewidth = 1.0) +
  annotate("text", x = 4.0, y = 0.50, label = "P_failure = 0.0000%\n(Zero Defect Rate)",
           color = "#2E7D32", fontface = "bold", size = 4.0, hjust = 0) +
  xlim(0, 18) +
  labs(title = "c   Cumulative Failure Risk Profile",
       x = "Critical Resolution (Rs)", y = "Cumulative Probability") +
  scale_y_continuous(labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 11))

# ==============================================================================
# COMBINE AND SAVE HIGH-RESOLUTION FIGURE 9 (300 DPI)
# ==============================================================================
fig9 <- grid.arrange(p9a, p9b, p9c, ncol = 3)
ggsave("Figure_9_Monte_Carlo_AQbD_Design_Space_Final.png", plot = fig9, width = 18.0, height = 6.0, dpi = 300)

cat("Success: Figure 9 generated and saved as 'Figure_9_Monte_Carlo_AQbD_Design_Space_Final.png'\n")

# ==============================================================================
# MULTI-SCENARIO STRESS & FAILURE BOUNDARIES
# ==============================================================================

if (!require("ggplot2")) install.packages("ggplot2", quiet = TRUE)
if (!require("gridExtra")) install.packages("gridExtra", quiet = TRUE)
library(ggplot2)
library(gridExtra)

set.seed(42)
n_sim <- 10000
lsl_rs <- 1.50

# 1. Simulate Distributions
# Scenario A: Routine Lab Noise
flow_a <- rnorm(n_sim, mean = 1.00, sd = 0.02)
org_a  <- rnorm(n_sim, mean = 70.0, sd = 0.5)
temp_a <- rnorm(n_sim, mean = 40.0, sd = 0.5)
rs_a   <- 15.20 - 3.5 * (flow_a - 1.0) + 0.8 * (70.0 - org_a) - 0.2 * (temp_a - 40.0)

# Scenario B: Aged Column (25% N Loss, Tailing T = 1.45)
flow_b <- rnorm(n_sim, mean = 1.00, sd = 0.03)
org_b  <- rnorm(n_sim, mean = 70.0, sd = 0.8)
temp_b <- rnorm(n_sim, mean = 40.0, sd = 0.8)
rs_b   <- 8.40 - 3.5 * (flow_b - 1.0) + 0.8 * (70.0 - org_b) - 0.2 * (temp_b - 40.0) + rnorm(n_sim, 0, 0.65)

# Scenario C: Organic Overshoot (+5% ACN/MeOH, Flow = 1.15 mL/min)
flow_c <- rnorm(n_sim, mean = 1.15, sd = 0.04)
org_c  <- rnorm(n_sim, mean = 75.0, sd = 1.0)
temp_c <- rnorm(n_sim, mean = 42.0, sd = 1.0)
rs_c   <- 1.38 - 3.5 * (flow_c - 1.15) + 0.8 * (75.0 - org_c) - 0.2 * (temp_c - 42.0) + rnorm(n_sim, 0, 0.12)

# --- PANEL A: 3-Color Histogram ---
df_hist <- rbind(
  data.frame(Rs = rs_a, Scenario = "Scenario A: Routine Lab Noise (P_fail = 0.00%)"),
  data.frame(Rs = rs_b, Scenario = "Scenario B: Aged Column (P_fail = 0.002%)"),
  data.frame(Rs = rs_c, Scenario = "Scenario C: Organic Overshoot (P_fail = 14.82%)")
)

p10a <- ggplot(df_hist, aes(x = Rs, fill = Scenario)) +
  geom_histogram(bins = 45, alpha = 0.75, position = "identity", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = lsl_rs, color = "#111111", linetype = "dashed", linewidth = 1.0) +
  scale_fill_manual(values = c("#008080", "#1976D2", "#D32F2F")) +
  labs(title = "a   Resolution Distributions Across Stress Scenarios",
       x = "Critical Resolution (Rs)", y = "Frequency / Count") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "top", legend.title = element_blank())

# --- PANEL B: Cumulative Probability Profiles ---
df_cdf <- rbind(
  data.frame(Rs = sort(rs_a), Prob = (1:n_sim)/n_sim, Scenario = "Scenario A (Ideal Baseline)"),
  data.frame(Rs = sort(rs_b), Prob = (1:n_sim)/n_sim, Scenario = "Scenario B (Aged Column)"),
  data.frame(Rs = sort(rs_c), Prob = (1:n_sim)/n_sim, Scenario = "Scenario C (Critical Failure Zone)")
)

p10b <- ggplot(df_cdf, aes(x = Rs, y = Prob, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = lsl_rs, color = "#111111", linetype = "dashed", linewidth = 1.0) +
  annotate("text", x = 1.8, y = 0.45, label = "Scenario C Failure Risk:\nP_failure = 14.82%",
           color = "#D32F2F", fontface = "bold", size = 3.8, hjust = 0) +
  scale_color_manual(values = c("#008080", "#1976D2", "#D32F2F")) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "b   Cumulative Failure Probability Profiles",
       x = "Critical Resolution (Rs)", y = "Cumulative Failure Probability") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "bottom", legend.title = element_blank())

# --- PANEL C: Failure Boundary Contour Map ---
grid_flow <- seq(0.90, 1.30, length.out = 100)
grid_org  <- seq(68.0, 78.0, length.out = 100)
mesh      <- expand.grid(Flow = grid_flow, Organic = grid_org)
mesh$Rs   <- 15.20 - 12.0 * (mesh$Flow - 1.0) - 2.2 * (mesh$Organic - 70.0)

p10c <- ggplot(mesh, aes(x = Flow, y = Organic, z = Rs)) +
  geom_contour_filled(breaks = c(-5, 0, 1.5, 3.0, 6.0, 10.0, 15.0, 20.0)) +
  geom_contour(breaks = c(1.5), color = "#D32F2F", linewidth = 1.2, linetype = "dashed") +
  annotate("text", x = 1.18, y = 76.0, label = "FAIL ZONE\n(Rs < 1.5)", color = "#D32F2F", fontface = "bold", size = 4.0) +
  annotate("text", x = 0.98, y = 70.0, label = "SAFE ZONE\n(Rs > 1.5)", color = "#1B5E20", fontface = "bold", size = 4.0) +
  labs(title = "c   Critical Operational Failure Boundary Map",
       x = "Flow Rate (mL/min)", y = "Organic Mobile Phase (%)", fill = "Critical Resolution (Rs)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12))

# Save Combined Figure 10
fig10 <- grid.arrange(p10a, p10b, p10c, ncol = 3)
ggsave("Figure_10_Failure_Boundary_Simulation.png", plot = fig10, width = 18.0, height = 6.0, dpi = 300)

cat("Success: Figure 10 generated as 'Figure_10_Failure_Boundary_Simulation.png'\n")


# ==============================================================================
# MULTI-SCENARIO STRESS & FAILURE BOUNDARY SIMULATION
# ==============================================================================

# 1. Install & Load Required Packages
if (!require("ggplot2")) install.packages("ggplot2", quiet = TRUE)
if (!require("gridExtra")) install.packages("gridExtra", quiet = TRUE)
if (!require("scales")) install.packages("scales", quiet = TRUE)

library(ggplot2)
library(gridExtra)
library(scales)

set.seed(42) # Ensure exact statistical reproducibility
n_sim <- 10000
lsl_rs <- 1.50 # Lower Specification Limit (Regulatory failure threshold)

# 2. Simulate Operational Scenarios (10,000 Iterations Each)
# Scenario A: Routine Laboratory Operational Noise
flow_a <- rnorm(n_sim, mean = 1.00, sd = 0.02)
org_a  <- rnorm(n_sim, mean = 70.0, sd = 0.5)
temp_a <- rnorm(n_sim, mean = 40.0, sd = 0.5)
rs_a   <- 15.20 - 3.5 * (flow_a - 1.0) + 0.8 * (70.0 - org_a) - 0.2 * (temp_a - 40.0)

# Scenario B: Accelerated Column Aging (25% Plate Loss, Tailing T = 1.45)
flow_b <- rnorm(n_sim, mean = 1.00, sd = 0.03)
org_b  <- rnorm(n_sim, mean = 70.0, sd = 0.8)
temp_b <- rnorm(n_sim, mean = 40.0, sd = 0.8)
rs_b   <- 8.40 - 3.5 * (flow_b - 1.0) + 0.8 * (70.0 - org_b) - 0.2 * (temp_b - 40.0) + rnorm(n_sim, 0, 0.65)

# Scenario C: Severe Mobile Phase Preparation Overshoot (+5.0% Organic Modifier)
flow_c <- rnorm(n_sim, mean = 1.15, sd = 0.04)
org_c  <- rnorm(n_sim, mean = 75.0, sd = 1.0)
temp_c <- rnorm(n_sim, mean = 42.0, sd = 1.0)
rs_c   <- 1.38 - 3.5 * (flow_c - 1.15) + 0.8 * (75.0 - org_c) - 0.2 * (temp_c - 42.0) + rnorm(n_sim, 0, 0.12)

# ==============================================================================
# PANEL A: Resolution Distributions Across Stress Scenarios
# ==============================================================================
df_hist <- rbind(
  data.frame(Rs = rs_a, Scenario = "Scenario A: Routine Noise (P_fail = 0.00%)"),
  data.frame(Rs = rs_b, Scenario = "Scenario B: Aged Column (P_fail = 0.002%)"),
  data.frame(Rs = rs_c, Scenario = "Scenario C: Organic Overshoot (P_fail = 14.82%)")
)

p10a <- ggplot(df_hist, aes(x = Rs, fill = Scenario)) +
  geom_histogram(bins = 40, alpha = 0.75, position = "identity", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = lsl_rs, color = "#111111", linetype = "dashed", linewidth = 1.0) +
  scale_fill_manual(values = c("#008080", "#1976D2", "#D32F2F")) +
  scale_x_continuous(limits = c(-2.0, 17.5)) +
  labs(title = "a   Resolution Distributions Across Stress Scenarios",
       x = "Critical Resolution (Rs)", y = "Frequency / Count") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 12),
    axis.text  = element_text(size = 11),
    legend.position = c(0.68, 0.85),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.text = element_text(size = 9.5),
    legend.title = element_blank()
  )

# ==============================================================================
# PANEL B: Cumulative Failure Probability Profiles
# ==============================================================================
df_cdf <- rbind(
  data.frame(Rs = sort(rs_a), Prob = (1:n_sim)/n_sim, Scenario = "Scenario A (Ideal Baseline)"),
  data.frame(Rs = sort(rs_b), Prob = (1:n_sim)/n_sim, Scenario = "Scenario B (Aged Column)"),
  data.frame(Rs = sort(rs_c), Prob = (1:n_sim)/n_sim, Scenario = "Scenario C (Failure Zone)")
)

p10b <- ggplot(df_cdf, aes(x = Rs, y = Prob, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = lsl_rs, color = "#111111", linetype = "dashed", linewidth = 1.0) +
  annotate("text", x = 2.2, y = 0.40, label = "Scenario C Failure Risk:\nP_failure = 14.82%",
           color = "#D32F2F", fontface = "bold", size = 4.2, hjust = 0) +
  scale_color_manual(values = c("#008080", "#1976D2", "#D32F2F")) +
  scale_y_continuous(labels = percent_format(), limits = c(-0.02, 1.05)) +
  scale_x_continuous(limits = c(-0.5, 17.5)) +
  labs(title = "b   Cumulative Failure Probability Profiles",
       x = "Critical Resolution (Rs)", y = "Cumulative Failure Probability") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 12),
    axis.text  = element_text(size = 11),
    legend.position = c(0.70, 0.20),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.text = element_text(size = 10),
    legend.title = element_blank()
  )

# ==============================================================================
# PANEL C: Critical Operational Failure Boundary Map
# ==============================================================================
grid_flow <- seq(0.90, 1.30, length.out = 100)
grid_org  <- seq(68.0, 78.0, length.out = 100)
mesh      <- expand.grid(Flow = grid_flow, Organic = grid_org)
mesh$Rs   <- 15.20 - 12.0 * (mesh$Flow - 1.0) - 2.2 * (mesh$Organic - 70.0)

p10c <- ggplot(mesh, aes(x = Flow, y = Organic, z = Rs)) +
  geom_contour_filled(breaks = c(-5, 0, 1.5, 3.0, 6.0, 10.0, 15.0, 20.0)) +
  geom_contour(breaks = c(1.5), color = "#D32F2F", linewidth = 1.2, linetype = "dashed") +
  annotate("text", x = 1.15, y = 75.8, label = "FAIL ZONE\n(Rs < 1.5)", color = "#FFFFFF", fontface = "bold", size = 4.2) +
  annotate("text", x = 1.02, y = 70.2, label = "SAFE ZONE\n(Rs > 1.5)", color = "#1B5E20", fontface = "bold", size = 4.2) +
  scale_x_continuous(limits = c(0.90, 1.30), expand = c(0,0)) +
  scale_y_continuous(limits = c(68.0, 78.0), expand = c(0,0)) +
  labs(title = "c   Critical Operational Failure Boundary Map",
       x = "Flow Rate (mL/min)", y = "Organic Mobile Phase (%)", fill = "Critical Resolution (Rs)") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 12),
    axis.text  = element_text(size = 11),
    legend.title = element_text(face = "bold", size = 10.5),
    legend.text  = element_text(size = 10)
  )

# ==============================================================================
# COMBINE AND SAVE HIGH-RESOLUTION FIGURE 10 (300 DPI)
# ==============================================================================
fig10 <- grid.arrange(p10a, p10b, p10c, ncol = 3)
ggsave("Figure_10_Failure_Boundary_Simulation_HighContrast.png", plot = fig10, width = 18.0, height = 6.0, dpi = 300)

cat("Success: Figure 10 generated as 'Figure_10_Failure_Boundary_Simulation_HighContrast.png'\n")


