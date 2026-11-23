#Midostaurin + Tretinoin HPLC Method Validation 
#Script: Data Cleaning 
#Purpose: Import and quality check raw HPLC data 
#Author: Darshan Patil 
#Project: Simultaneous RP-HPLC-PDA method development

install.packages("readxl")
install.packages("dplyr")
install.packages("tidyr")

library(readxl)
library(dplyr)
library(tidyr)

linearity.raw<-read_excel(
  "D:/Academics/research article/HPLC MID TRE/Linearity.xlsx")
print(linearity.raw)
str(linearity.raw)
names(linearity.raw)
dim(linearity.raw)
colSums(linearity.raw)

#Quality checks 
colSums(is.na(linearity.raw))
sum(duplicated(linearity.raw))
sum(duplicated(linearity.raw$Conc))
summary(linearity.raw)

#Clean Analysis Dataset
linearity<-linearity.raw %>%
  rename(
    concentration=Conc,
    Midostaurin_Area=Mid,
    Tretinoin_Area=Tre
  )
print(linearity)

str(linearity)
summary(linearity)

#Response factor 
linearity <- linearity %>%
  mutate(
    Midostaurin_RF = ifelse(
      concentration > 0,
      Midostaurin_Area / concentration,
      NA
    ),
    Tretinoin_RF = ifelse(
      concentration > 0,
      Tretinoin_Area / concentration,
      NA
    )
  )


print(linearity)

linearity %>%
  summarise(
    Midostaurin_RF_mean=mean(Midostaurin_RF, na.rm=TRUE),
    Midostaurin_RF_SD=sd(Midostaurin_RF, na.rm=TRUE),
    Midostaurin_RF_RSD=100*Midostaurin_RF_SD/Midostaurin_RF_mean,
    
    Tretinoin_RF_mean=mean(Tretinoin_RF, na.rm=TRUE),
    Tretinoin_RF_SD=sd(Tretinoin_RF, na.rm=TRUE),
    Tretinoin_RF_RSD=100*Tretinoin_RF_SD/Tretinoin_RF_mean
  )

#Midostaurin calibration response 
library(ggplot2)
ggplot(linearity%>% filter(concentration>0),
       aes(x=concentration, y=Midostaurin_Area))+
  geom_point(size=3)+
  labs(
    title = "Midostaurin calibration response",
    x = "Concentration (ug/mL)",
    y = "Peak area"
  )+
  theme_minimal()

#Tretinoin calibration response 
ggplot(linearity %>% filter(concentration>0),
       aes(x = concentration, y = Tretinoin_Area))+
  geom_point(size = 3)+
  labs(
    title = "Tretinoin calibration response",
    x = "Concentration (ug/mL)",
    y = "Peak area"
  )+
  theme_minimal()

#Linearity: Calibration Regression 
mid_model<- lm(Midostaurin_Area~concentration, data=linearity)
tre_model<- lm(Tretinoin_Area~concentration, data = linearity)
summary(mid_model)
summary(tre_model)

#Regression Statistics 
regression_stats<-function(model){
  s<-summary(model)
  data.frame(
    Slope=coef(model)[2],
    Intercept=coef(model)[1],
    R_squared= s$r.squared, 
    Adjusted_R_squared=s$adj.r.squared,
    Residual_SE=sigma(model)
  )
}
mid_stats<-regression_stats(mid_model)
tre_stats<-regression_stats(tre_model)
mid_stats
tre_stats
confint(mid_model, level=0.95)

confint(tre_model, level=0.95)
par(mfrow=c(2,2))
dev.off()
plot(mid_model)
par(mfrow=c(2,2))
plot(tre_model)

mid_residuals<-data.frame(
  concentration=linearity$concentration,
  observed=linearity$Midostaurin_Area,
  fitted=fitted(mid_model),
  residual=residuals(mid_model)
)

tre_residuals<-data.frame(
  concentration=linearity$concentration,
  observed=linearity$Tretinoin_Area,
  fitted=fitted(tre_model),
  residual=residuals(tre_model)
)
mid_residuals
tre_residuals

#Residual Analysis 
mid_residuals<-data.frame(
  concentration=linearity$concentration,
  observed=linearity$Midostaurin_Area,
  fitted=fitted(mid_model),
  residual=residuals(mid_model)
)
 tre_residuals<-data.frame(
   concentration=linearity$concentration,
   observed=linearity$Tretinoin_Area,
   fitted=fitted(tre_model),
   residual=residuals(tre_model)
 )
mid_residuals
tre_residuals

#Percentage residuals 
mid_residuals$percent_residual <- ifelse(
  mid_residuals$observed == 0,
  NA,
  100 * mid_residuals$residual / mid_residuals$observed
)

tre_residuals$percent_residual <- ifelse(
  tre_residuals$observed == 0,
  NA,
  100 * tre_residuals$residual / tre_residuals$observed
)

mid_residuals
tre_residuals

#Residual plots
library(ggplot2)
ggplot(mid_residuals,
       aes(x=concentration, y=residual))+
  geom_point(size=3)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  labs(
    title = "Midostaurin: Residual Analysis",
    x= "Concentration (ug/mL)",
    y= "Residual peak area"
  )+
  theme_minimal()

ggplot(tre_residuals,
       aes(x=concentration, y=residual))+
  geom_point(size = 3)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  labs(
    title = "Tretinoin: Residual Analysis",
    x= "Concentration (ug/mL)",
    y= "Residual peak area"
  )+
  theme_minimal()

#Formal regression diagnostics 
par(mfrow=c(2,2))
dev.off()
plot(mid_model)

par(mfrow=c(2,2))
plot(tre_model)

#ANOVA for the regression 
anova(mid_model)
anova(tre_model)
summary(mid_model)
summary(tre_model)

#Complete calibration statistics 
#Midostaurin 
mid_coef<-summary(mid_model)$coefficients
mid_r2<-summary(mid_model)$r.squared
mid_adj_r2<-summary(mid_model)$adj.r.squared

#Tretinoin 
tre_coef<-summary(tre_model)$coefficients
tre_r2<-summary(tre_model)$r.squared
tre_adj_r2<-summary(tre_model)$adj.r.squared

mid_coef
mid_r2
mid_adj_r2

tre_coef
tre_r2
tre_adj_r2

#Predicted concentrations and residuals
linearity_stats <- linearity %>%
  mutate(
    Mid_predicted = predict(mid_model),
    Tre_predicted = predict(tre_model),
    Mid_residual = Midostaurin_Area - Mid_predicted,
    Tre_residual = Tretinoin_Area - Tre_predicted
  )

linearity_stats

#% deviation from the regression 
linearity_stats<-linearity_stats %>%
  mutate(
    Mid_Deviation=ifelse(
      concentration==0,
      NA,
      (Mid_residual/Mid_predicted)*100
    ),
    Tre_Deviation=ifelse(
      concentration==0,
      NA,
      (Tre_residual/Tre_predicted)*100
    )
  )

linearity_stats
View(linearity_stats)

mid_model_noblank<-lm(
  Midostaurin_Area~concentration,
  data=linearity%>%filter(concentration>0)
)
tre_model_noblank<-lm(
  Tretinoin_Area~concentration,
  data=linearity%>%filter(concentration>0)
)
summary(mid_model_noblank)
summary(tre_model_noblank)
View(mid_model_noblank)

linearity_backcalc <- linearity %>%
  filter(concentration > 0) %>%
  mutate(
    Mid_predicted_area = predict(mid_model_noblank),
    Tre_predicted_area = predict(tre_model_noblank),
    
    Mid_back_calculated =
      (Midostaurin_Area - coef(mid_model_noblank)[1]) /
      coef(mid_model_noblank)[2],
    
    Tre_back_calculated =
      (Tretinoin_Area - coef(tre_model_noblank)[1]) /
      coef(tre_model_noblank)[2],
    
    Mid_deviation_percent =
      ((Mid_back_calculated - concentration) /
         concentration) * 100,
    
    Tre_deviation_percent =
      ((Tre_back_calculated - concentration) /
         concentration) * 100
  )

print(linearity_backcalc)

max(abs(linearity_backcalc$Mid_deviation_percent))
max(abs(linearity_backcalc$Tre_deviation_percent))

mean(abs(linearity_backcalc$Mid_deviation_percent))
mean(abs(linearity_backcalc$Tre_deviation_percent))

#Final linearity table 
linearity_table<-linearity_backcalc%>%
  select(
    concentration,
    Midostaurin_Area,
    Mid_back_calculated,
    Mid_deviation_percent,
    Tretinoin_Area,
    Tre_back_calculated,
    Tre_deviation_percent
  )
linearity_table

View(linearity_table)

#Regression Summary Table

regression_table <- data.frame(
  Parameter = c(
    "Concentration range",
    "Number of calibration levels",
    "Slope",
    "Intercept",
    "R-squared",
    "Adjusted R-squared",
    "Residual standard error",
    "Slope p-value",
    "Intercept p-value"
  ),
  
  Midostaurin = c(
    "50-150 µg/mL",
    5,
    unname(coef(mid_model_noblank)[2]),
    unname(coef(mid_model_noblank)[1]),
    summary(mid_model_noblank)$r.squared,
    summary(mid_model_noblank)$adj.r.squared,
    sigma(mid_model_noblank),
    summary(mid_model_noblank)$coefficients[2, 4],
    summary(mid_model_noblank)$coefficients[1, 4]
  ),
  
  Tretinoin = c(
    "50-150 µg/mL",
    5,
    unname(coef(tre_model_noblank)[2]),
    unname(coef(tre_model_noblank)[1]),
    summary(tre_model_noblank)$r.squared,
    summary(tre_model_noblank)$adj.r.squared,
    sigma(tre_model_noblank),
    summary(tre_model_noblank)$coefficients[2, 4],
    summary(tre_model_noblank)$coefficients[1, 4]
  )
)

print(regression_table)
View(regression_table)

#Midostaurin calibration graph
# ============================================================
# PUBLICATION-READY CALIBRATION CURVE: MIDOSTAURIN
# ============================================================

library(ggplot2)

mid_plot <- ggplot(
  linearity %>% filter(concentration > 0),
  aes(x = concentration, y = Midostaurin_Area)
) +
  
  # Experimental points
  geom_point(
    size = 3.2
  ) +
  
  # Linear regression
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.9
  ) +
  
  # Regression equation and statistics
  annotate(
    "text",
    x = 55,
    y = 3900000,
    label = paste0(
      "y = ",
      round(coef(mid_model_noblank)[2], 2),
      "x ",
      ifelse(coef(mid_model_noblank)[1] >= 0, "+ ", "- "),
      abs(round(coef(mid_model_noblank)[1], 2)),
      "\nR² = ",
      round(summary(mid_model_noblank)$r.squared, 4),
      "\nn = 5 (50–150 µg/mL)"
    ),
    hjust = 0,
    size = 5
  ) +
  
  labs(
    title = "Midostaurin calibration curve",
    x = "Concentration (µg/mL)",
    y = "Peak area"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(
      size = 14,
      face = "bold"
    ),
    axis.text = element_text(
      size = 12
    ),
    panel.border = element_rect(
      fill = NA,
      linewidth = 0.8
    )
  )

mid_plot

#Tretinoin calibration graph 
# ============================================================
# PUBLICATION-READY CALIBRATION CURVE: TRETINOIN
# ============================================================

tre_plot <- ggplot(
  linearity %>% filter(concentration > 0),
  aes(x = concentration, y = Tretinoin_Area)
) +
  
  # Experimental points
  geom_point(
    size = 3.2
  ) +
  
  # Linear regression
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.9
  ) +
  
  # Regression equation and statistics
  annotate(
    "text",
    x = 55,
    y = 7200000,
    label = paste0(
      "y = ",
      round(coef(tre_model_noblank)[2], 2),
      "x ",
      ifelse(coef(tre_model_noblank)[1] >= 0, "+ ", "- "),
      abs(round(coef(tre_model_noblank)[1], 2)),
      "\nR² = ",
      round(summary(tre_model_noblank)$r.squared, 4),
      "\nn = 5 (50–150 µg/mL)"
    ),
    hjust = 0,
    size = 5
  ) +
  
  labs(
    title = "Tretinoin calibration curve",
    x = "Concentration (µg/mL)",
    y = "Peak area"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(
      size = 14,
      face = "bold"
    ),
    axis.text = element_text(
      size = 12
    ),
    panel.border = element_rect(
      fill = NA,
      linewidth = 0.8
    )
  )

tre_plot

# for 600dpi
ggsave(
  "Midostaurin_calibration_Curve.tiff",
  plot = mid_plot,
  width = 7,
  height = 5.5,
  units = "in",
  dpi = 600,
  compression="lzw"
)

ggsave(
  "Tretinoin_Calibration_Curve.tiff",
  plot = tre_plot,
  width = 7,
  height = 5.5,
  units = "in",
  dpi = 600,
  compression="lzw"
)

getOption("repos")
download.file(
  "https://cloud.r-project.org/",
  tempfile(),
  quiet = FALSE
)


install.packages("patchwork")
library(patchwork)
combined_calibration<-mid_plot+tre_plot+
  plot_annotation(
    tag_levels = "A"
  )
combined_calibration
library(ggplot2)
ggsave(
  "Figure_Calibration_Curves.tiff",
  plot = combined_calibration,
  width = 14,
  height = 5.5,
  units = "in",
  dpi = 600,
  compression="lzw"
)

#Regression plot - Midostaurin 
library(ggplot2)
library(dplyr)


mid_regression_plot<-ggplot(
  linearity %>% filter(concentration>0),
  aes(x=concentration, y=Midostaurin_Area)
)+
  geom_point(size=3)+
  geom_smooth(
    method = "lm",
    formula = y~x,
    se=TRUE,
    linewidth = 1
  )+
  annotate(
    "text",
    x=52,
    y=3800000,
    hjust=0,
    size=4.5,
    label=paste0(
      "y=",
      format(round(coef(mid_model_noblank)[2],2),
             big.mark = ","),
      "x",
      ifelse(coef(mid_model_noblank)[1]<0,"-","+"),
      format(abs(round(coef(mid_model_noblank)[1],2)),
             big.mark = ","),
      "\nR2=",
      round(summary(mid_model_noblank)$r.squared,4)
      
    )
  )+
  scale_x_continuous(
    breaks = c(50,75,100,125,150)
  )+
  labs(
    title = "Regression analysis of Midostaurin",
    x= "Concentration (ug/mL)",
    y="Peak area"
  )+
  theme_classic(base_size = 14)+
  theme(
    plot.title = element_text(
      face="bold",
      hjust = 0.5
    ),
    axis.title = element_text(face = "bold"),
    panel.border = element_rect(
      fill = NA,
      linewidth = 0.8
    )
  )
mid_regression_plot

#Regression plot Tretinoin 
tre_regression_plot<-ggplot(
  linearity%>%filter(concentration>0),
  aes(x=concentration, y=Tretinoin_Area)
)+
  geom_point(size=3)+
  geom_smooth(
    method = "lm",
    formula = y~x,
    se=TRUE,
    linewidth = 1
  )+
  annotate(
    "text",
    x=52,
    y=7000000,
    hjust=0,
    size=4.5,
    label=paste0(
      "y=",
      format(round(coef(tre_model_noblank)[2],2),
             big.mark = ","),
      "x",
      ifelse(coef(tre_model_noblank)[1]<0,"-","+"),
      format(abs(round(coef(tre_model_noblank)[1],2)),
             big.mark = ","),
      "\nR2=",
      round(summary(tre_model_noblank)$r.squared,4)
    )
  )+
  scale_x_continuous(
    breaks = c(50,75,100,125,150)
  )+
  labs(
    title = "Regression analysis of Tretinoin",
    x= "Concentration (ug/mL)",
    y= "Peak area"
  )+
  theme_classic(base_size = 14)+
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(face = "bold"),
    panel.border = element_rect(
      fill = NA,
      linewidth = 0.8
    )
  )
tre_regression_plot

ggplot(
  mid_residuals,
  aes(x = concentration, y = residual)
) +
  geom_point(size = 3) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    linewidth = 0.8
  ) +
  labs(
    title = "Residual analysis of Midostaurin calibration",
    x = "Concentration (µg/mL)",
    y = "Residual peak area"
  ) +
  theme_classic(base_size = 14)

ggplot(
  tre_residuals,
  aes(x = concentration, y = residual)
) +
  geom_point(size = 3) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    linewidth = 0.8
  ) +
  labs(
    title = "Residual analysis of Tretinoin calibration",
    x = "Concentration (µg/mL)",
    y = "Residual peak area"
  ) +
  theme_classic(base_size = 14)

#Formal residual/variance assessment 
install.packages("lmtest")
library(lmtest)
#Midostaurin
bptest(mid_model_noblank)
#Tretinoin 
bptest(tre_model_noblank)

#The null hypothesis is 
#H0= Residual variance is constant (homoscedasticity).
# If p>0.05, there is no statistically significant evidence of heteroscedasticity
# If p<0.05, there is evidence that the variance changes with concentration. 

# Shapiro Wilk test of residuals 
#Midostaurin 
shapiro.test(residuals(mid_model_noblank))
#Tretinoin 
shapiro.test((residuals(tre_model_noblank)))
#Here H0= residuals are consistent with a normal distribution 

#Influence diagnostics 
mid_influence<-data.frame(
  concentration=linearity%>%
    filter(concentration>0)%>%
    pull(concentration),
  leverage=hatvalues(mid_model_noblank),
  cooks_distance= cooks.distance(mid_model_noblank),
  standardized_residual=rstandard(mid_model_noblank)
)
tre_influence<-data.frame(
  concentration=linearity%>%
    filter(concentration>0)%>%
    pull(concentration),
  leverage=hatvalues(tre_model_noblank),
  cooks_distance=cooks.distance(tre_model_noblank),
  standardized_residual=rstandard(tre_model_noblank)
)
mid_influence
tre_influence

#Residual vs. concentration figure 
#Midostaurin
mid_residual_plot<-ggplot(
  mid_residuals%>% filter(concentration>0),
  aes(x=concentration, y=residual)
)+
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  )+
  geom_point(size = 3)+
  geom_smooth(
    method="loess",
    formula = y~x,
    se=TRUE,
    linewidth = 0.8
  )+
  labs(
    title = "Residual analysis:Midostaurin",
    x= "Concentration (ug/mL)",
    y= "Residual peak area"
  )+
  theme_classic(base_size = 14)
mid_residual_plot

#Tretinoin 
tre_residual_plot<-ggplot(
  tre_residuals%>% filter(concentration>0),
  aes(x=concentration, y=residual)
)+
  geom_hline(
    yintercept =0,
    linetype = "dashed"
  )+
  geom_point(size = 3)+
  geom_smooth(
    method = "loess",
    formula = y~x,
    se=TRUE,
    linewidth = 0.8
  )+
  labs(
    title = "Residual analysis: Tretinoin",
    x= "Concentration (ug/mL)",
    y= "Residual peak area"
  )+
  theme_classic(base_size = 14)
tre_residual_plot


bptest(mid_model_noblank)
bptest(tre_model_noblank)
shapiro.test(residuals(mid_model_noblank))
shapiro.test((residuals(tre_model_noblank)))
mid_influence
tre_influence

#Sensitivity analysis influence of 150 ug/mL
#Midostaurin
mid_model_125<-lm(
  Midostaurin_Area~concentration,
  data=linearity%>%
    filter(concentration>0, concentration<=125)
)
summary(mid_model_125)
#Tretinoin 
tre_model_125<-lm(
  Tretinoin_Area~concentration, 
  data=linearity%>%
    filter(concentration>0, concentration<=125)
)
summary(tre_model_125)
#Compare slopes
coef(mid_model_noblank)
coef(mid_model_125)
coef(tre_model_noblank)
coef(tre_model_125)
#Compare R2
summary(mid_model_noblank)$r.squared
summary(mid_model_125)$r.squared
summary(tre_model_noblank)$r.squared
summary(tre_model_125)$r.squared

#150ug/mL has appreciable influence on the fitted regression, particularly on
#the residual error and intercept, but the fundamental concentration response 
#relationship remains strongly linear. 

#Conclusion Midostaurin 
#Linear response was observed over 50-150 ug/mL, with R2=0.9973. the 150 ug/mL
#level exerted measurable influence on the regression, but sensitivity analysis 
#excluding this level retained a strong linear relationship (R2=0.9992). The
#slope changed by approximately 4.8%, supporting retention of the 150 ug/mL level
#in the primary calibration model. 

#Conclusion Tretinoin 
#Linear response was observed over 50-150ug/mL, R2=0.9992, while the slope
#changed by approximately 2.8%, including that the overall concentration response 
#relationship remained stable. 
#Darshan 





