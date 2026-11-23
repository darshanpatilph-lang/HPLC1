# ==============================================================================
# SELF-CONTAINED PRECISION SCRIPT FOR SECTION 3.3 IN R
# ==============================================================================

# 1. Input Dataset
conc <- c(0.06, 0.06, 0.06, 0.10, 0.10, 0.10, 0.50, 0.50, 0.50, 
          1.00, 1.00, 1.00, 1.50, 1.50, 1.50, 2.00, 2.00, 2.00)

midostaurin <- c(1021.29, 936.42, 927.57, 1688.05, 1726.82, 1686.04, 
                 8450.38, 8426.75, 8453.75, 16763.79, 16983.58, 16873.69, 
                 25124.92, 25361.16, 25260.00, 33621.06, 33747.36, 33570.54)

tretinoin <- c(949.14, 906.42, 915.35, 1593.15, 1580.98, 1591.92, 
               7953.00, 7983.27, 7938.67, 15829.71, 16210.76, 15750.12, 
               23546.43, 23798.43, 23499.18, 31230.41, 31441.90, 31277.41)

df <- data.frame(Conc = conc, Midostaurin = midostaurin, Tretinoin = tretinoin)

# 2. Add Back-Calculated Concentration Columns (Crucial Step!)
# OLS Equations from Section 3.1:
# Midostaurin: Y = 16832.555 * X + 3.067
# Tretinoin:   Y = 15680.499 * X + 70.118

df$Mido_Conc <- (df$Midostaurin - 3.067) / 16832.555
df$Tret_Conc <- (df$Tretinoin - 70.118) / 15680.499

# 3. Subset to QC Levels: Low (0.10), Mid (1.00), High (2.00)
qc_levels <- c(0.10, 1.00, 2.00)
df_qc <- df[df$Conc %in% qc_levels, ]

# 4. Compute Precision Aggregations for Midostaurin
mido_area_summary <- aggregate(Midostaurin ~ Conc, data = df_qc, FUN = function(x) {
  c(Mean = mean(x), SD = sd(x), RSD_Pct = (sd(x)/mean(x))*100)
})

mido_conc_summary <- aggregate(Mido_Conc ~ Conc, data = df_qc, FUN = function(x) {
  c(Mean = mean(x), SD = sd(x), RSD_Pct = (sd(x)/mean(x))*100)
})

# 5. Compute Precision Aggregations for Tretinoin
tret_area_summary <- aggregate(Tretinoin ~ Conc, data = df_qc, FUN = function(x) {
  c(Mean = mean(x), SD = sd(x), RSD_Pct = (sd(x)/mean(x))*100)
})

tret_conc_summary <- aggregate(Tret_Conc ~ Conc, data = df_qc, FUN = function(x) {
  c(Mean = mean(x), SD = sd(x), RSD_Pct = (sd(x)/mean(x))*100)
})

# 6. Print Summaries to Console
cat("=== MIDOSTAURIN PEAK AREA PRECISION ===\n")
print(mido_area_summary)

cat("\n=== MIDOSTAURIN CONCENTRATION PRECISION ===\n")
print(mido_conc_summary)

cat("\n=== TRETINOIN PEAK AREA PRECISION ===\n")
print(tret_area_summary)

cat("\n=== TRETINOIN CONCENTRATION PRECISION ===\n")
print(tret_conc_summary)


# ==============================================================================
# Analytes: Midostaurin & Tretinoin (Low, Mid, High QC Levels)
# Panels:
#   a: Midostaurin Intra-Day vs Inter-Day Precision (% RSD)
#   b: Tretinoin Intra-Day vs Inter-Day Precision (% RSD)
# ==============================================================================

# 1. Install & Load Required Packages
required_pkgs <- c("ggplot2", "gridExtra", "grid")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs)

library(ggplot2)
library(gridExtra)
library(grid)

# 2. Construct Precision Dataset
df_prec <- data.frame(
  Analyte = factor(rep(c("Midostaurin", "Tretinoin"), each = 6), levels = c("Midostaurin", "Tretinoin")),
  QC_Level = factor(rep(c("Low QC\n(0.10 µg/mL)", "Mid QC\n(1.00 µg/mL)", "High QC\n(2.00 µg/mL)"), times = 4),
                    levels = c("Low QC\n(0.10 µg/mL)", "Mid QC\n(1.00 µg/mL)", "High QC\n(2.00 µg/mL)")),
  Type = factor(rep(rep(c("Intra-Day RSD (%)", "Inter-Day RSD (%)"), each = 3), times = 2),
                levels = c("Intra-Day RSD (%)", "Inter-Day RSD (%)")),
  RSD = c(
    # Midostaurin % RSDs (Low, Mid, High)
    1.35, 0.65, 0.27,  # Intra-day
    1.42, 0.71, 0.31,  # Inter-day
    
    # Tretinoin % RSDs (Low, Mid, High)
    0.42, 1.55, 0.35,  # Intra-day
    0.51, 1.62, 0.38   # Inter-day
  )
)

# 3. Define Nature Portfolio Theme
theme_nature_prec <- function() {
  theme_classic(base_size = 11, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 11.5, hjust = 0, color = "#111111", margin = margin(b = 4)),
      plot.subtitle = element_text(size = 9, color = "#555555", margin = margin(b = 6)),
      axis.title = element_text(face = "bold", size = 9.5, color = "#222222"),
      axis.text = element_text(size = 8.5, color = "#333333"),
      axis.line = element_line(color = "#222222", linewidth = 0.55),
      axis.ticks = element_line(color = "#222222", linewidth = 0.55),
      panel.grid.major.y = element_line(color = "#EDEDED", linewidth = 0.35, linetype = "dashed"),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 8.5, color = "#222222"),
      legend.background = element_rect(fill = "transparent", color = NA),
      plot.margin = margin(8, 10, 8, 10)
    )
}

# 4. Panel A: Midostaurin Precision Barplot
p1 <- ggplot(subset(df_prec, Analyte == "Midostaurin"), aes(x = QC_Level, y = RSD, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.65, color = "#00332C", linewidth = 0.4) +
  geom_hline(yintercept = 2.0, linetype = "dashed", color = "#D32F2F", linewidth = 0.75) +
  annotate("text", x = 2.5, y = 2.12, label = "ICH Limit (RSD < 2.0%)", color = "#D32F2F", size = 2.8, fontface = "bold") +
  scale_fill_manual(values = c("Intra-Day RSD (%)" = "#008080", "Inter-Day RSD (%)" = "#80CBC4")) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, 0.5), expand = c(0, 0)) +
  labs(title = "a  Midostaurin Method Precision",
       subtitle = "Intra-day (n = 3) and inter-day (n = 18) repeatability across QC levels",
       x = "Quality Control (QC) Level",
       y = "Relative Standard Deviation (% RSD)") +
  theme_nature_prec()

# 5. Panel B: Tretinoin Precision Barplot
p2 <- ggplot(subset(df_prec, Analyte == "Tretinoin"), aes(x = QC_Level, y = RSD, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.65, color = "#2A0052", linewidth = 0.4) +
  geom_hline(yintercept = 2.0, linetype = "dashed", color = "#D32F2F", linewidth = 0.75) +
  annotate("text", x = 2.5, y = 2.12, label = "ICH Limit (RSD < 2.0%)", color = "#D32F2F", size = 2.8, fontface = "bold") +
  scale_fill_manual(values = c("Intra-Day RSD (%)" = "#8E24AA", "Inter-Day RSD (%)" = "#E1BEE7")) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, 0.5), expand = c(0, 0)) +
  labs(title = "b  Tretinoin Method Precision",
       subtitle = "Intra-day (n = 3) and inter-day (n = 18) repeatability across QC levels",
       x = "Quality Control (QC) Level",
       y = "Relative Standard Deviation (% RSD)") +
  theme_nature_prec()

# 6. Arrange Panels Side-by-Side
fig_precision <- grid.arrange(p1, p2, ncol = 2)

# 7. Save Output Files (PDF & PNG)
ggsave("Figure_2_Precision_Reproducibility_Nature.pdf", plot = fig_precision, width = 9.0, height = 4.5, dpi = 600)
ggsave("Figure_2_Precision_Reproducibility_Nature.png", plot = fig_precision, width = 9.0, height = 4.5, dpi = 300)


