# NOTE: this script's original input file, Seqtab_Linked_TP2.csv, was not
# among the recovered project files and is not included in this repository.
# The expected columns are SampleID/Treatment/Block plus ASV_* count columns.
#
# Merged with the near-duplicate Absulote_Abundance.R (also fixed the
# filename typo), which was a strict subset of this script's ANOVA/Tukey/plot
# steps, missing the homogeneity-of-variance, normality, and Kruskal-Wallis
# checks below.

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(agricolae)
library(car)

# Set the file path
file_path <- "Seqtab_Linked_TP2.csv"

# Load the data
data <- read.csv(file_path)

# Calculate the total absolute abundance of all ASVs per sample
data <- data %>%
  mutate(Total_Abundance = rowSums(select(., starts_with("ASV"))))

# Perform ANOVA for total absolute abundance
anova_total_abundance <- aov(Total_Abundance ~ Treatment + Block, data = data)
summary(anova_total_abundance)

# Perform Tukey HSD for total absolute abundance
tukey_total_abundance <- HSD.test(anova_total_abundance, "Treatment", group = TRUE)

# Create table for total absolute abundance
total_abundance_table <- data %>%
  group_by(Treatment) %>%
  summarise(avg_total_abundance = mean(Total_Abundance)) %>%
  left_join(tibble(Treatment = rownames(tukey_total_abundance$groups), letter = tukey_total_abundance$groups$groups), by = "Treatment")

print(total_abundance_table)

# Plot total absolute abundance with Tukey groups
ggplot(data, aes(x = Treatment, y = Total_Abundance, fill = Treatment)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge()) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2, position = position_dodge(0.9)) +
  theme_minimal() +
  labs(title = "Total Absolute Abundance by Treatment", y = "Mean Total Absolute Abundance", x = "Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data = total_abundance_table, aes(label = letter, y = avg_total_abundance + 0.5), vjust = -0.5)

# Save the plot
ggsave("total_abundance_tukey_plot.png")

# Test for homogeneity of variances (Levene's Test)
levene_test <- leveneTest(Total_Abundance ~ Treatment, data = data)
print(levene_test)

# Test for normality of residuals (Shapiro-Wilk Test)
shapiro_test <- shapiro.test(residuals(anova_total_abundance))
print(shapiro_test)


##################
# Perform Kruskal-Wallis test for total absolute abundance
kruskal_test <- kruskal.test(Total_Abundance ~ Treatment, data = data)
print(kruskal_test)

# Perform pairwise comparisons using Dunn's test
library(FSA)
dunn_test <- dunnTest(Total_Abundance ~ Treatment, data = data, method = "bonferroni")
print(dunn_test)
