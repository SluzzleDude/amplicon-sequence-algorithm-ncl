# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(multcomp)
library(agricolae)
library(ggpubr)

# Set the path to the CSV file
file_path <- "ITS_Seqab_Metadata_TP3.csv"

# Read the data
data <- read.csv(file_path)

# Calculate the total abundance for each sample
# Assuming ASVs start from column 8 and there are 7949 ASV columns
data <- data %>%
  mutate(Total_Abundance = rowSums(across(8:(8 + 7949 - 1))))

# ANOVA for Pesticide and Pathogen
anova_model <- aov(Total_Abundance ~ Pesticide * Pathogen + Block, data = data)
summary(anova_model)

# Tukey HSD for mean comparisons
tukey_results <- HSD.test(anova_model, "Pesticide", group = TRUE)
print(tukey_results)

# Manually input the p-value from ANOVA
manual_p_value <- 0.263  # Replace with your actual p-value

# Visualize the expression profile
ggboxplot(data, x = "Pesticide", y = "Total_Abundance", color = "Pesticide", 
          add = "jitter", legend = "none") +
  rotate_x_text(angle = 45) +
  annotate("text", x = 1, y = 60000, label = paste("ANOVA, p =", manual_p_value), hjust = 0) +  # Manually add p-value
  theme(legend.position = "none",
        plot.subtitle = element_text(hjust = 0.5)) + 
  labs(subtitle = "Fungal Abundance at Disease End-Point",
       x = "Pesticide",
       y = "Total Abundance") +
  ylim(0, 20000) +  # Set y-axis range
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE)  # Remove NS or *** above the whiskers
