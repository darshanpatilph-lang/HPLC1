# ==============================================================================
# SECTION 3.4: ACCURACY & SPIKE-RECOVERY CALCULATIONS IN R
# ==============================================================================

# 1. Input Standard Addition Spike Levels (0.50, 1.00, 1.50 ug/mL)
conc <- c(0.50, 0.50, 0.50, 1.00, 1.00, 1.00, 1.50, 1.50, 1.50)

midostaurin <- c(8450.38, 8426.75, 8453.75, 
                 16763.79, 16983.58, 16873.69, 
                 25124.92, 25361.16, 25260.00)

tretinoin   <- c(7953.00, 7983.27, 7938.67, 
                 15829.71, 16210.76, 15750.12, 
                 23546.43, 23798.43, 23499.18)

df_acc <- data.frame(Conc = conc, Midostaurin = midostaurin, Tretinoin = tretinoin)

# 2. Back-calculate Concentrations via Inverse Regression
df_acc$Mido_Conc <- (df_acc$Midostaurin - 3.067) / 16832.555
df_acc$Tret_Conc <- (df_acc$Tretinoin - 70.118) / 15680.499

# 3. Calculate Spike Recoveries (% Recovery = (Found / Nominal) * 100)
df_acc$Mido_Recovery <- (df_acc$Mido_Conc / df_acc$Conc) * 100
df_acc$Tret_Recovery <- (df_acc$Tret_Conc / df_acc$Conc) * 100

# 4. Aggregate Recovery Statistics
mido_rec_summary <- aggregate(Mido_Recovery ~ Conc, data = df_acc, FUN = function(x) {
  c(Mean = mean(x), SD = sd(x), RSD_Pct = (sd(x)/mean(x))*100)
})

tret_rec_summary <- aggregate(Tret_Recovery ~ Conc, data = df_acc, FUN = function(x) {
  c(Mean = mean(x), SD = sd(x), RSD_Pct = (sd(x)/mean(x))*100)
})

# 5. One-Sample Student's t-Test Against 100.0% Nominal Recovery
mido_ttest <- t.test(df_acc$Mido_Recovery, mu = 100.0)
tret_ttest <- t.test(df_acc$Tret_Recovery, mu = 100.0)

# Print Summary Results
cat("=== MIDOSTAURIN SPIKE-RECOVERY ===\n")
print(mido_rec_summary)

cat("\n=== TRETINOIN SPIKE-RECOVERY ===\n")
print(tret_rec_summary)

cat("\n=== MIDOSTAURIN t-TEST FOR SYSTEMATIC BIAS ===\n")
print(mido_ttest)

cat("\n=== TRETINOIN t-TEST FOR SYSTEMATIC BIAS ===\n")
print(tret_ttest)


# ==============================================================================
# FIGURE 3: SPIKE-RECOVERY ACCURACY BARPLOT IN R (CORRECTED)
# ==============================================================================

library(ggplot2)

df_fig3 <- data.frame(
  Spike_Level = factor(rep(c("80% Level\n(0.50 µg/mL)", "100% Level\n(1.00 µg/mL)", "120% Level\n(1.50 µg/mL)"), times = 2),
                       levels = c("80% Level\n(0.50 µg/mL)", "100% Level\n(1.00 µg/mL)", "120% Level\n(1.50 µg/mL)")),
  Analyte = factor(rep(c("Midostaurin", "Tretinoin"), each = 3)),
  Recovery = c(100.29, 100.23, 99.99, 100.61, 101.15, 100.10)
)

p_fig3 <- ggplot(df_fig3, aes(x = Spike_Level, y = Recovery, fill = Analyte)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.65, color = "#222222", linewidth = 0.4) +
  geom_hline(yintercept = 100.0, linetype = "solid", color = "#111111", linewidth = 0.8) +
  geom_hline(yintercept = c(98.0, 102.0), linetype = "dashed", color = "#D32F2F", linewidth = 0.75) +
  scale_fill_manual(values = c("Midostaurin" = "#008080", "Tretinoin" = "#8E24AA")) +
  scale_y_continuous(breaks = seq(95, 105, 1), expand = c(0, 0)) +
  coord_cartesian(ylim = c(95, 105)) +  # FIX: Zooms view without dropping bars
  labs(title = "Method Spike-Recovery Accuracy Profile",
       subtitle = "Percentage recoveries across 80%, 100%, and 120% nominal levels (n = 9)",
       x = "Spike Addition Level", y = "Mean Percentage Recovery (%)") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11.5),
        legend.position = "top", legend.title = element_blank())

# Save High-Resolution PNG & PDF Outputs
ggsave("Figure_3_Accuracy_Recovery_Nature.png", plot = p_fig3, width = 7.5, height = 4.5, dpi = 300)
ggsave("Figure_3_Accuracy_Recovery_Nature.pdf", plot = p_fig3, width = 7.5, height = 4.5, dpi = 600)
print(p_fig3)



