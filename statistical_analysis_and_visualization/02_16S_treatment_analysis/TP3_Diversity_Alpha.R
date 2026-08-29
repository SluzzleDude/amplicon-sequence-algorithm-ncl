# Alpha diversity (Shannon, Simpson) by pesticide and pathogen treatment,
# TP3 (16S). Merges in the presentation-styled Shannon boxplot from
# TP3_Diversity_Alpha_beau.R (its second plot duplicated the styling already
# covered by TP3_Diversity_Alpha_Pub.R, so only the distinct one is kept).

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(agricolae)

data <- read.csv("TP3_SeqTab_Meta.csv")

data$Shannon <- diversity(data[, 8:3852], index = "shannon")
data$Simpson <- diversity(data[, 8:3852], index = "simpson")

anova_shannon <- aov(Shannon ~ Pesticide * Pathogen + Block, data = data)
summary(anova_shannon)
tukey_shannon <- HSD.test(anova_shannon, "Pesticide", group = TRUE)
print(tukey_shannon$groups)

anova_simpson <- aov(Simpson ~ Pesticide * Pathogen + Block, data = data)
summary(anova_simpson)
tukey_simpson <- HSD.test(anova_simpson, "Pesticide", group = TRUE)
print(tukey_simpson$groups)

summary_table_shannon <- data %>%
  group_by(Pesticide) %>%
  summarise(Average_Shannon = mean(Shannon), .groups = "drop") %>%
  left_join(data.frame(Pesticide = rownames(tukey_shannon$groups), Shannon_Letter = tukey_shannon$groups$groups),
            by = "Pesticide")

summary_table_simpson <- data %>%
  group_by(Pesticide) %>%
  summarise(Average_Simpson = mean(Simpson), .groups = "drop") %>%
  left_join(data.frame(Pesticide = rownames(tukey_simpson$groups), Simpson_Letter = tukey_simpson$groups$groups),
            by = "Pesticide")

write.csv(summary_table_shannon, "Shannon_Summary_Table.csv", row.names = FALSE)
write.csv(summary_table_simpson, "Simpson_Summary_Table.csv", row.names = FALSE)

# Shannon and Simpson by pesticide, with Tukey group letters
ggplot(data, aes(x = Pesticide, y = Shannon)) +
  geom_boxplot() +
  geom_text(data = summary_table_shannon, aes(x = Pesticide, y = max(data$Shannon) + 0.2, label = Shannon_Letter), size = 5) +
  labs(title = "Shannon Diversity Index by Pesticide", y = "Shannon Diversity Index") +
  theme_minimal()

ggplot(data, aes(x = Pesticide, y = Simpson)) +
  geom_boxplot() +
  geom_text(data = summary_table_simpson, aes(x = Pesticide, y = max(data$Simpson) + 0.02, label = Simpson_Letter), size = 5) +
  labs(title = "Simpson Diversity Index by Pesticide", y = "Simpson Diversity Index") +
  theme_minimal()

# Shannon and Simpson by pathogen presence
ggplot(data, aes(x = Pathogen, y = Shannon)) +
  geom_boxplot() +
  labs(title = "Shannon Diversity Index by Pathogen", y = "Shannon Diversity Index") +
  theme_minimal()

ggplot(data, aes(x = Pathogen, y = Simpson)) +
  geom_boxplot() +
  labs(title = "Simpson Diversity Index by Pathogen", y = "Simpson Diversity Index") +
  theme_minimal()

# Presentation-styled Shannon boxplot, by pesticide
ggplot(data, aes(x = Pesticide, y = Shannon, fill = Pesticide)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 1.5, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, fill = "white") +
  labs(title = "Shannon Diversity Index by Pesticide Treatment", x = "Pesticide Treatment", y = "Shannon Diversity Index") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")
