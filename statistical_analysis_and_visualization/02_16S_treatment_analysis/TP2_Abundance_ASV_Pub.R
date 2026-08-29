# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(multcomp)
library(agricolae)
library(ggpubr)

# Set the path to the CSV file
file_path <- "TP2_SeqTab_Meta.csv"

# Read the data
data <- read.csv(file_path)

# Calculate the total abundance for each sample
data <- data %>%
  mutate(Total_Abundance = rowSums(across(starts_with("ASV"))))

# ANOVA for Pesticide and Pathogen
anova_model <- aov(Total_Abundance ~ Pesticide * Pathogen + Block, data = data)
summary(anova_model)

# Tukey HSD for mean comparisons
tukey_results <- HSD.test(anova_model, "Pesticide", group = TRUE)
print(tukey_results)

# Manually input the p-value from ANOVA
manual_p_value <- 0.306  # Replace with your actual p-value



# Visualize the expression profile
ggboxplot(data, x = "Pesticide", y = "Total_Abundance", color = "Pesticide", 
          add = "jitter", legend = "none") +
  rotate_x_text(angle = 45) +
  annotate("text", x = 1, y = 35000, label = paste("ANOVA, p =", manual_p_value), hjust = 0) +  # Manually add p-value
  theme(legend.position = "none",
        plot.subtitle = element_text(hjust = 0.5)) + 
  labs(subtitle = "Bacterial Abundance at Disease Mid-Point",
       x = "Pesticide",
       y = "Total Abundance") +
  ylim(0, 35000) +  # Set y-axis range
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE)  # Remove NS or *** above the whiskers
