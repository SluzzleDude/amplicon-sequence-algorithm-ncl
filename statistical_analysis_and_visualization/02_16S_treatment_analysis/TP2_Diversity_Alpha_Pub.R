# Load necessary libraries
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(agricolae)
library(ggpubr)

# Load data
data <- read.csv("TP2_SeqTab_Meta.csv")

# Calculate Shannon diversity index
data$Shannon <- diversity(data[,8:3852], index = "shannon")

# Perform ANOVA for Shannon diversity
anova_shannon <- aov(Shannon ~ Pesticide * Pathogen + Block, data = data)
summary(anova_shannon)

# Perform Tukey HSD for Shannon diversity
tukey_shannon <- HSD.test(anova_shannon, "Pesticide", group = TRUE)
print(tukey_shannon$groups)

# Manually input the p-value from ANOVA
manual_p_value <- 0.046255   # Replace with your actual p-value

# Visualize the expression profile
ggboxplot(data, x = "Pesticide", y = "Shannon", color = "Pesticide", 
          add = "jitter", legend = "none") +
  rotate_x_text(angle = 45) +
  annotate("text", x = 1, y = max(data$Shannon) + 0.5, label = paste("ANOVA, p =", manual_p_value), hjust = 0) +  # Manually add p-value
  theme(legend.position = "none",
        plot.subtitle = element_text(hjust = 0.5)) + 
  labs(subtitle = "Bacterial Species Diversity at Disease Mid-Point",
       x = "Pesticide",
       y = "Shannon Diversity") +
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE)  # Remove NS or *** above the whiskers
