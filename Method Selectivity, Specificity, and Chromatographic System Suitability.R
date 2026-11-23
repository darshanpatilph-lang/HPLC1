# ==============================================================================
# SECTION 3.5: CHROMATOGRAPHIC SYSTEM SUITABILITY CALCULATIONS IN R
# ==============================================================================

# 1. Input Experimental System Suitability Data (n = 12 Injections)
inj <- 1:12

mido_rt <- c(3.915, 3.911, 3.915, 3.916, 3.917, 3.918, 3.918, 3.922, 3.923, 3.921, 3.921, 3.922)
tret_rt <- c(9.175, 9.159, 9.178, 9.180, 9.192, 9.186, 9.186, 9.189, 9.194, 9.197, 9.200, 9.203)

mido_plates <- c(8564, 8501, 8560, 8520, 8510, 8540, 8550, 8530, 8525, 8540, 8515, 8505)
tret_plates <- c(15612, 15580, 15600, 15550, 15570, 15590, 15610, 15540, 15530, 15580, 15560, 15575)

mido_tf <- c(1.185, 1.189, 1.186, 1.188, 1.187, 1.185, 1.188, 1.189, 1.186, 1.187, 1.188, 1.186)
tret_tf <- c(1.073, 1.070, 1.075, 1.078, 1.072, 1.076, 1.074, 1.077, 1.075, 1.073, 1.076, 1.074)

# 2. Compute Summary Statistics
mido_rt_mean <- mean(mido_rt); mido_rt_rsd <- (sd(mido_rt)/mido_rt_mean)*100
tret_rt_mean <- mean(tret_rt); tret_rt_rsd <- (sd(tret_rt)/tret_rt_mean)*100

cat("Midostaurin tR:", mido_rt_mean, "min (RSD =", mido_rt_rsd, "%)\n")
cat("Midostaurin N:", mean(mido_plates), "plates\n")
cat("Midostaurin Tf:", mean(mido_tf), "\n\n")

cat("Tretinoin tR:", tret_rt_mean, "min (RSD =", tret_rt_rsd, "%)\n")
cat("Tretinoin N:", mean(tret_plates), "plates\n")
cat("Tretinoin Tf:", mean(tret_tf), "\n")
cat("Chromatographic Resolution (Rs): 6.78\n")


# ==============================================================================
# FIGURE 4: CHROMATOGRAPHIC SYSTEM SUITABILITY & SEQUENCE STABILITY IN R
# Analytes: Midostaurin & Tretinoin (n = 12 Injections)
# Panels:
#   a: Retention Time Consistency (t_R) Across 12 Sequence Injections
#   b: System Performance Metrics (Plates, Tailing, Resolution) vs. USP Limits
# ==============================================================================

# 1. Install & Load Packages
required_pkgs <- c("ggplot2", "gridExtra", "grid")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs)

library(ggplot2)
library(gridExtra)
library(grid)

# 2. Panel A Dataset: Retention Times across 12 injections
inj_seq <- 1:12
mido_rt <- c(3.915, 3.911, 3.915, 3.916, 3.917, 3.918, 3.918, 3.922, 3.923, 3.921, 3.921, 3.922)
tret_rt <- c(9.175, 9.159, 9.178, 9.180, 9.192, 9.186, 9.186, 9.189, 9.194, 9.197, 9.200, 9.203)

df_rt <- data.frame(
  Injection = rep(inj_seq, 2),
  Retention_Time = c(mido_rt, tret_rt),
  Analyte = factor(rep(c("Midostaurin (t_R = 3.918 min)", "Tretinoin (t_R = 9.187 min)"), each = 12))
)

# 3. Panel B Dataset: System Suitability Parameters
df_metrics <- data.frame(
  Parameter = factor(rep(c("Theoretical\nPlates (N × 10³)", "Tailing\nFactor (T_f)", "Chromatographic\nResolution (R_s)"), each = 2),
                     levels = c("Theoretical\nPlates (N × 10³)", "Tailing\nFactor (T_f)", "Chromatographic\nResolution (R_s)")),
  Analyte = factor(rep(c("Midostaurin", "Tretinoin"), times = 3), levels = c("Midostaurin", "Tretinoin")),
  Value = c(8.530, 15.574,   # Theoretical Plates (in thousands)
            1.187, 1.075,    # Tailing Factor
            0, 6.786)        # Resolution
)

# 4. Nature Portfolio Theme
theme_nature <- function() {
  theme_classic(base_size = 11, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 11.5, hjust = 0, color = "#111111", margin = margin(b = 4)),
      axis.title = element_text(face = "bold", size = 9.5, color = "#222222"),
      axis.text = element_text(size = 8.5, color = "#333333"),
      axis.line = element_line(color = "#222222", linewidth = 0.55),
      axis.ticks = element_line(color = "#222222", linewidth = 0.55),
      panel.grid.major = element_line(color = "#EDEDED", linewidth = 0.35, linetype = "dashed"),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 8.5, color = "#222222"),
      plot.margin = margin(8, 10, 8, 10)
    )
}

# 5. Panel A Plot
p1 <- ggplot(df_rt, aes(x = Injection, y = Retention_Time, color = Analyte, shape = Analyte)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("Midostaurin (t_R = 3.918 min)" = "#008080", "Tretinoin (t_R = 9.187 min)" = "#8E24AA")) +
  scale_x_continuous(breaks = 1:12) +
  scale_y_continuous(limits = c(3.0, 10.0), breaks = seq(3, 10, 1)) +
  labs(title = "a  Retention Time Consistency (n = 12)",
       x = "Injection Sequence Number", y = "Retention Time (t_R, min)") +
  theme_nature() +
  theme(legend.position = "right")

# 6. Panel B Plot
p2 <- ggplot(df_metrics, aes(x = Parameter, y = Value, fill = Analyte)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.65, color = "#222222", linewidth = 0.4) +
  geom_hline(yintercept = 2.0, linetype = "dashed", color = "#D32F2F", linewidth = 0.75) +
  scale_fill_manual(values = c("Midostaurin" = "#008080", "Tretinoin" = "#8E24AA")) +
  scale_y_continuous(limits = c(0, 18), breaks = seq(0, 18, 2), expand = c(0, 0)) +
  labs(title = "b  Chromatographic Performance Parameters",
       x = "", y = "Parameter Value") +
  theme_nature() +
  theme(legend.position = "right")

# 7. Arrange Panels and Export
fig4 <- grid.arrange(p1, p2, ncol = 2)

ggsave("Figure_4_System_Suitability_Nature.pdf", plot = fig4, width = 10.0, height = 4.2, dpi = 600)
ggsave("Figure_4_System_Suitability_Nature.png", plot = fig4, width = 10.0, height = 4.2, dpi = 300)


library(ggplot2)

df_acc <- data.frame(
  Level = factor(
    rep(c("50% Level\n(0.50 µg/mL)", "100% Level\n(1.00 µg/mL)", "150% Level\n(1.50 µg/mL)"), 2),
    levels = c("50% Level\n(0.50 µg/mL)", "100% Level\n(1.00 µg/mL)", "150% Level\n(1.50 µg/mL)")
  ),
  Analyte = rep(c("Midostaurin", "Tretinoin"), each = 3),
  Recovery = c(100.29, 100.23, 99.99, 100.61, 101.15, 100.10)
)

p3 <- ggplot(df_acc, aes(x = Level, y = Recovery, fill = Analyte)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65, color = "black", linewidth = 0.4) +
  geom_hline(yintercept = 100, linetype = "solid", color = "black", linewidth = 0.6) +
  geom_hline(yintercept = c(98, 102), linetype = "dashed", color = "#B22222", linewidth = 0.6) +
  scale_fill_manual(values = c("Midostaurin" = "#008080", "Tretinoin" = "#8B1A8B")) +
  coord_cartesian(ylim = c(95, 105)) +
  scale_y_continuous(breaks = 95:105) +
  labs(
    title = "Method Spike-Recovery Accuracy Profile",
    subtitle = "Percentage recoveries across 50%, 100%, and 150% nominal levels (n = 9)",
    x = "Spike Addition Level",
    y = "Mean Percentage Recovery (%)",
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "top",
    axis.title = element_text(face = "bold")
  )

print(p3)
ggsave("Figure_3_Accuracy_Profile.png", p3, width = 8.5, height = 5.5, dpi = 300)




















