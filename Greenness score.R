
# ================================================================
# FIGURE 7: GREEN ANALYTICAL CHEMISTRY ASSESSMENT
# Three coordinated circular panels:
# A. AGREE (12 principles)
# B. AGREEprep (10 criteria)
# C. ComplexMoGAPI (6 lifecycle stages)
#
# Method: RP-HPLC-PDA
# Mobile phase: ACN:aqueous 0.1% OPA (85:15, v/v)
# Flow rate: 1.00 mL/min
# Runtime: 12 min
#
# Output:
# 1. Figure_7A_AGREE.png / PDF
# 2. Figure_7B_AGREEprep.png / PDF
# 3. Figure_7C_ComplexMoGAPI.png / PDF
# 4. Figure_7_Combined.png / PDF
# ================================================================


# 1. PACKAGES -----------------------------------------------------

packages <- c("ggplot2", "patchwork", "dplyr")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(ggplot2)
library(patchwork)
library(dplyr)


# 2. COLOUR PALETTE -----------------------------------------------

GREEN_DARK <- "#176B35"
GREEN_LIGHT <- "#65A83F"
YELLOW <- "#F5B51B"
RED <- "#D92D27"

LOW_GREEN <- "#198344"
MODERATE_YELLOW <- "#F4B51C"

# Score colour classification for AGREE and AGREEprep
score_colour <- function(x) {
  ifelse(
    x >= 0.80, GREEN_DARK,
    ifelse(x >= 0.60, GREEN_LIGHT,
           ifelse(x >= 0.40, YELLOW, RED))
  )
}


# 3. AGREE DATA ---------------------------------------------------

agree <- data.frame(
  Code = paste0("P", 1:12),
  
  Principle = c(
    "Direct analysis",
    "Minimal sample",
    "In-situ analysis",
    "Process integration",
    "Automation",
    "No derivatization",
    "Waste generation",
    "Multi-analyte",
    "Energy",
    "Renewable reagents",
    "Toxic reagents",
    "Operator safety"
  ),
  
  Score = c(
    0.85, 0.95, 0.40, 0.88,
    0.80, 1.00, 0.62, 0.83,
    0.75, 0.38, 0.72, 0.90
  )
)

agree_overall <- 0.78


# 4. AGREEprep DATA -----------------------------------------------

agreeprep <- data.frame(
  Code = paste0("C", 1:10),
  
  Criterion = c(
    "Preparation strategy",
    "Hazardous materials",
    "Renewability",
    "Waste generation",
    "Sample consumption",
    "Throughput",
    "Automation",
    "Energy consumption",
    "Post-preparation compatibility",
    "Operator safety"
  ),
  
  Score = c(
    0.90, 0.70, 0.45, 0.85, 0.95,
    0.90, 0.80, 0.85, 1.00, 0.95
  )
)

agreeprep_overall <- 0.84


# 5. COMPLEXMoGAPI DATA -------------------------------------------

complex_data <- data.frame(
  Code = c("I", "II", "III", "IV", "V", "VI"),
  
  Stage = c(
    "Sample collection",
    "Sample preparation",
    "Reagents and solvents",
    "Instrumentation",
    "Waste management",
    "Pre-synthetic lifecycle"
  ),
  
  Short = c(
    "SAMPLE\nCOLLECTION",
    "SAMPLE\nPREPARATION",
    "REAGENTS &\nSOLVENTS",
    "INSTRUMENTATION",
    "WASTE\nMANAGEMENT",
    "PRE-SYNTHETIC\nLIFECYCLE"
  ),
  
  Impact = c(
    "Low impact",
    "Low impact",
    "Moderate impact",
    "Moderate impact",
    "Moderate impact",
    "Low impact"
  )
)

complex_colours <- c(
  "Low impact" = LOW_GREEN,
  "Moderate impact" = MODERATE_YELLOW
)


# 6. FUNCTION TO GENERATE CIRCULAR SCORE FIGURES ------------------

make_score_circle <- function(
    data,
    overall,
    panel_title,
    centre_name,
    score_column = "Score",
    code_column = "Code",
    inner_radius = 0.48,
    max_radius = 1.00,
    gap = 0.025
) {
  
  n <- nrow(data)
  
  # Start at 12 o'clock, proceeding clockwise
  start_angles <- pi / 2 - (0:(n - 1)) * 2 * pi / n
  end_angles <- pi / 2 - (1:n) * 2 * pi / n
  
  polygons <- list()
  labels <- list()
  
  for (i in seq_len(n)) {
    
    score <- data[[score_column]][i]
    
    # Radial length reflects the individual score
    outer_radius <- inner_radius +
      score * (max_radius - inner_radius)
    
    # Segment gap
    a1 <- start_angles[i] - gap
    a2 <- end_angles[i] + gap
    
    theta <- seq(a1, a2, length.out = 70)
    
    # Outer arc
    outer <- data.frame(
      x = outer_radius * cos(theta),
      y = outer_radius * sin(theta)
    )
    
    # Inner arc
    inner <- data.frame(
      x = inner_radius * cos(rev(theta)),
      y = inner_radius * sin(rev(theta))
    )
    
    poly <- rbind(outer, inner)
    
    poly$segment <- i
    poly$fill <- score_colour(score)
    
    polygons[[i]] <- poly
    
    # Score label positioned inside the segment
    middle_angle <- (a1 + a2) / 2
    
    label_radius <- inner_radius +
      (outer_radius - inner_radius) * 0.57
    
    labels[[i]] <- data.frame(
      x = label_radius * cos(middle_angle),
      y = label_radius * sin(middle_angle),
      label = sprintf("%.2f", score)
    )
  }
  
  polygon_data <- bind_rows(polygons)
  score_labels <- bind_rows(labels)
  
  # Principle / criterion codes around the outside
  code_data <- data.frame(
    Code = data[[code_column]],
    angle = (start_angles + end_angles) / 2
  )
  
  code_data$x <- 1.12 * cos(code_data$angle)
  code_data$y <- 1.12 * sin(code_data$angle)
  
  # Circular score chart
  p <- ggplot() +
    
    geom_polygon(
      data = polygon_data,
      aes(
        x = x,
        y = y,
        group = segment,
        fill = fill
      ),
      colour = "white",
      linewidth = 1.2
    ) +
    
    scale_fill_identity() +
    
    # Numeric scores
    geom_text(
      data = score_labels,
      aes(x = x, y = y, label = label),
      size = 3.3,
      fontface = "bold",
      colour = "white"
    ) +
    
    # Principle/criterion labels outside the circle
    geom_text(
      data = code_data,
      aes(x = x, y = y, label = Code),
      size = 3.8,
      fontface = "bold",
      colour = "black"
    ) +
    
    # Centre circle
    annotate(
      "point",
      x = 0,
      y = 0,
      size = 31,
      shape = 21,
      fill = "white",
      colour = "grey75",
      stroke = 0.8
    ) +
    
    annotate(
      "text",
      x = 0,
      y = 0.15,
      label = centre_name,
      size = 5.2,
      fontface = "bold",
      lineheight = 1.0
    ) +
    
    annotate(
      "text",
      x = 0,
      y = -0.10,
      label = sprintf("Overall\n%.2f / 1.00", overall),
      size = 4.0,
      fontface = "bold",
      lineheight = 1.05
    ) +
    
    coord_fixed(
      xlim = c(-1.28, 1.28),
      ylim = c(-1.28, 1.28),
      clip = "off"
    ) +
    
    labs(
      title = panel_title
    ) +
    
    theme_void() +
    
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5,
        colour = "black",
        margin = margin(b = 8)
      ),
      
      plot.margin = margin(
        t = 15, r = 15, b = 15, l = 15
      )
    )
  
  return(p)
}


# 7. FIGURE 7A: AGREE ---------------------------------------------

fig7A <- make_score_circle(
  data = agree,
  overall = agree_overall,
  panel_title = "(A) AGREE Metric",
  centre_name = "AGREE\n12 PRINCIPLES"
)

print(fig7A)


# 8. FIGURE 7B: AGREEprep -----------------------------------------

fig7B <- make_score_circle(
  data = agreeprep,
  overall = agreeprep_overall,
  panel_title = "(B) AGREEprep Metric",
  centre_name = "AGREEprep\n10 CRITERIA"
)

print(fig7B)


# 9. FIGURE 7C: COMPLEXMoGAPI CIRCULAR LIFECYCLE ------------------

make_complex_circle <- function(data) {
  
  n <- nrow(data)
  
  start_angles <- pi / 2 - (0:(n - 1)) * 2 * pi / n
  end_angles <- pi / 2 - (1:n) * 2 * pi / n
  
  polygons <- list()
  labels <- list()
  
  for (i in seq_len(n)) {
    
    theta <- seq(
      start_angles[i] - 0.025,
      end_angles[i] + 0.025,
      length.out = 90
    )
    
    outer_radius <- 1.00
    inner_radius <- 0.42
    
    outer <- data.frame(
      x = outer_radius * cos(theta),
      y = outer_radius * sin(theta)
    )
    
    inner <- data.frame(
      x = inner_radius * cos(rev(theta)),
      y = inner_radius * sin(rev(theta))
    )
    
    poly <- rbind(outer, inner)
    
    poly$segment <- i
    poly$Impact <- data$Impact[i]
    
    polygons[[i]] <- poly
    
    mid_angle <- (start_angles[i] + end_angles[i]) / 2
    
    labels[[i]] <- data.frame(
      x = 0.70 * cos(mid_angle),
      y = 0.70 * sin(mid_angle),
      label = paste0(
        data$Code[i], "\n",
        data$Short[i], "\n",
        toupper(data$Impact[i])
      )
    )
  }
  
  polygon_data <- bind_rows(polygons)
  label_data <- bind_rows(labels)
  
  ggplot() +
    
    geom_polygon(
      data = polygon_data,
      aes(
        x = x,
        y = y,
        group = segment,
        fill = Impact
      ),
      colour = "white",
      linewidth = 1.5
    ) +
    
    scale_fill_manual(
      values = complex_colours,
      name = "Environmental classification"
    ) +
    
    geom_text(
      data = label_data,
      aes(
        x = x,
        y = y,
        label = label
      ),
      size = 3.3,
      fontface = "bold",
      lineheight = 1.0,
      colour = "black"
    ) +
    
    # White central circle
    annotate(
      "point",
      x = 0,
      y = 0,
      size = 27,
      shape = 21,
      fill = "white",
      colour = "grey75",
      stroke = 0.8
    ) +
    
    annotate(
      "text",
      x = 0,
      y = 0.12,
      label = "COMPLEX\nMoGAPI",
      size = 5.0,
      fontface = "bold",
      lineheight = 1.0
    ) +
    
    annotate(
      "text",
      x = 0,
      y = -0.14,
      label = "12-min\n85:15 ACN:OPA",
      size = 3.2,
      fontface = "bold",
      lineheight = 1.0
    ) +
    
    coord_fixed(
      xlim = c(-1.12, 1.12),
      ylim = c(-1.12, 1.12),
      clip = "off"
    ) +
    
    labs(
      title = "(C) ComplexMoGAPI Lifecycle Assessment",
      subtitle = "Six analytical lifecycle stages"
    ) +
    
    theme_void() +
    
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5,
        colour = "black"
      ),
      
      plot.subtitle = element_text(
        size = 11,
        face = "bold",
        hjust = 0.5,
        colour = "grey25"
      ),
      
      legend.position = "bottom",
      
      legend.title = element_text(
        size = 11,
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 10,
        face = "bold"
      ),
      
      legend.key.size = grid::unit(0.5, "cm"),
      
      plot.margin = margin(
        t = 15, r = 15, b = 15, l = 15
      )
    )
}

fig7C <- make_complex_circle(complex_data)

print(fig7C)


# 10. COMBINE THREE FIGURES ---------------------------------------

# Equal-sized panels, arranged horizontally.
# Use a wide output canvas to retain readable labels.

combined_figure <- (
  fig7A | fig7B | fig7C
) +
  
  plot_layout(
    nrow = 1,
    widths = c(1, 1, 1),
    guides = "collect"
  ) +
  
  plot_annotation(
    title = "Figure 7. Green Analytical Chemistry Assessment",
    subtitle = paste(
      "Simultaneous RP-HPLC-PDA determination of",
      "midostaurin and tretinoin"
    ),
    theme = theme(
      plot.title = element_text(
        size = 21,
        face = "bold",
        hjust = 0.5
      ),
      
      plot.subtitle = element_text(
        size = 13,
        face = "bold",
        hjust = 0.5
      )
    )
  )

print(combined_figure)


# 11. EXPORT INDIVIDUAL PANELS ------------------------------------

ggsave(
  "Figure_7A_AGREE.png",
  fig7A,
  width = 7,
  height = 7,
  dpi = 600,
  bg = "white"
)

ggsave(
  "Figure_7B_AGREEprep.png",
  fig7B,
  width = 7,
  height = 7,
  dpi = 600,
  bg = "white"
)

ggsave(
  "Figure_7C_ComplexMoGAPI.png",
  fig7C,
  width = 7,
  height = 7,
  dpi = 600,
  bg = "white"
)


# 12. EXPORT COMBINED FIGURE --------------------------------------

ggsave(
  "Figure_7_Combined.png",
  combined_figure,
  width = 21,
  height = 8,
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

# Vector PDF for journal submission
ggsave(
  "Figure_7_Combined.pdf",
  combined_figure,
  width = 21,
  height = 8,
  device = grDevices::cairo_pdf,
  bg = "white",
  limitsize = FALSE
)


# 13. EXPORT SOURCE DATA ------------------------------------------

write.csv(
  agree,
  "AGREE_Principle_Scores.csv",
  row.names = FALSE
)

write.csv(
  agreeprep,
  "AGREEprep_Criteria_Scores.csv",
  row.names = FALSE
)

write.csv(
  complex_data,
  "ComplexMoGAPI_Lifecycle_Classifications.csv",
  row.names = FALSE
)

cat(
  "\nSUCCESS: All three circular panels generated.\n",
  "Individual panels: 600 DPI PNG\n",
  "Combined figure: 600 DPI PNG and vector PDF\n"
)
