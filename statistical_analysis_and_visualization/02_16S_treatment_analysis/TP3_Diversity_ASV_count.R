# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(agricolae)

# Load data
data <- read.csv("TP3_SeqTab_Meta.csv")

# Calculate the count of different ASVs
data$ASV_count <- rowSums(data[, 8:3852] > 0)

# Perform ANOVA for ASV count
anova_asv <- aov(ASV_count ~ Pesticide * Pathogen + Block, data = data)
summary(anova_asv)

# Perform Tukey HSD for ASV count
tukey_asv <- HSD.test(anova_asv, "Pesticide", group = TRUE)
print(tukey_asv$groups)

# Create a summary table with mean ASV count and letters
summary_table_asv <- data %>%
  group_by(Pesticide) %>%
  summarise(Average_ASV_count = mean(ASV_count), .groups = 'drop')

# Ensure column names are correct in tukey_asv$groups
tukey_asv_groups <- data.frame(Pesticide = rownames(tukey_asv$groups), ASV_Letter = tukey_asv$groups[, "groups"])

summary_table_asv <- summary_table_asv %>%
  left_join(tukey_asv_groups, by = "Pesticide")

# Save summary table
write.csv(summary_table_asv, "ASV_Count_Summary_Table.csv", row.names = FALSE)

# Create plots for ASV count
ggplot(data, aes(x = Pesticide, y = ASV_count)) +
  geom_boxplot() +
  geom_text(data = summary_table_asv, aes(x = Pesticide, y = max(data$ASV_count) + 5, label = ASV_Letter), size = 5) +
  labs(title = "ASV Count by Pesticide", y = "ASV Count") +
  theme_minimal()

# Create separate plots for Pathogen
ggplot(data, aes(x = Pathogen, y = ASV_count)) +
  geom_boxplot() +
  labs(title = "ASV Count by Pathogen", y = "ASV Count") +
  theme_minimal()

