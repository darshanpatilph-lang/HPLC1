# ==============================================================================
# Corrected R Script for Figure 7 (Reflecting 60.64% Intact Tretinoin)
# ==============================================================================

library(ggplot2)
library(gridExtra)

# 1. Dataset Construction with Updated Photolytic Values
deg_data <- data.frame(
  Stress_Condition = factor(
    rep(c("Acid Hydrolysis\n(0.1 M HCl, 60°C)", 
          "Base Hydrolysis\n(0.1 M NaOH, 60°C)", 
          "Oxidative Stress\n(3% H2O2, 25°C)", 
          "Thermal Stress\n(80°C Dry Heat)", 
          "Photolytic Stress\n(UV/Vis Light)", 
          "Untreated Control\n(Baseline Ref)"), each = 2),
    levels = c("Acid Hydrolysis\n(0.1 M HCl, 60°C)", 
               "Base Hydrolysis\n(0.1 M NaOH, 60°C)", 
               "Oxidative Stress\n(3% H2O2, 25°C)", 
               "Thermal Stress\n(80°C Dry Heat)", 
               "Photolytic Stress\n(UV/Vis Light)", 
               "Untreated Control\n(Baseline Ref)")
  ),
  Analyte = factor(rep(c("Midostaurin", "Tretinoin"), times = 6), 
                   levels = c("Midostaurin", "Tretinoin")),
  Intact_Remaining_Pct = c(71.27, 78.47, 
                           79.87, 59.74, 
                           90.64, 68.19, 
                           92.41, 70.60, 
                           93.89, 60.64,   # <-- Updated Photolytic: Mido = 93.89%, Tre = 60.64%
                           100.00, 100.00),
  Rate_Constant_k = c(0.1693, 0.1212, 
                      0.1124, 0.2576, 
                      0.0491, 0.1915, 
                      0.0395, 0.1741, 
                      0.0315, 0.2501,     # <-- Updated Photolytic k_obs: Mido = 0.0315, Tre = 0.2501
                      0.0000, 0.0000)
)

teal_color   <- "#008080"  # Midostaurin
purple_color <- "#7B1FA2"  # Tretinoin

# 2. Panel A: Grouped Bar Chart of Intact Remaining (%)
p7a <- ggplot(deg_data, aes(x = Stress_Condition, y = Intact_Remaining_Pct, fill = Analyte)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 100.0, color = "black", linetype = "solid", linewidth = 0.8) +
  geom_hline(yintercept = 80.0, color = "#B71C1C", linetype = "dashed", linewidth = 0.7) +
  scale_fill_manual(values = c("Midostaurin" = teal_color, "Tretinoin" = purple_color)) +
  scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, by = 20), expand = c(0, 0)) +
  geom_text(aes(label = sprintf("%.1f%%", Intact_Remaining_Pct)), 
            position = position_dodge(width = 0.8), vjust = -0.5, size = 3.2, fontface = "bold") +
  labs(
    title = "(a) Forced Degradation Profiles Relative to Untreated Reference Control",
    x = "ICH Q1A(R2) Stress Degradation Condition",
    y = "Intact Drug Remaining (%)",
    fill = "Tested Analyte"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 10, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(face = "bold", size = 11),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 10)
  )

# 3. Panel B: Apparent First-Order Rate Constants (k_obs, h⁻¹)
kinetic_data <- subset(deg_data, Rate_Constant_k > 0)

p7b <- ggplot(kinetic_data, aes(x = Stress_Condition, y = Rate_Constant_k, fill = Analyte)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("Midostaurin" = teal_color, "Tretinoin" = purple_color)) +
  scale_y_continuous(limits = c(0, 0.32), breaks = seq(0, 0.30, by = 0.05), expand = c(0, 0)) +
  geom_text(aes(label = sprintf("%.4f", Rate_Constant_k)), 
            position = position_dodge(width = 0.8), vjust = -0.5, size = 3.0, fontface = "bold") +
  labs(
    title = "(b) Apparent First-Order Degradation Rate Constants (k_obs, h⁻¹)",
    x = "Stress Degradation Pathway",
    y = "Rate Constant k_obs (h⁻¹)",
    fill = "Tested Analyte"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 9.5, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(face = "bold", size = 11),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 10)
  )

# 4. Render Figure 7
grid.arrange(p7a, p7b, nrow = 2, heights = c(1, 1))

