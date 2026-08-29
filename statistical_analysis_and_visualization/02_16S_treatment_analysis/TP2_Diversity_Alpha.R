# Alpha diversity (Shannon, Simpson) and total ASV abundance by pesticide
# and pathogen treatment, TP2 (16S).
#
# Consolidated from TP2_Diversity_Alpha.R, TP2_Diversity_Alpha_beau.R,
# TP2_Abundance_ASV_New.R and TP2_Abundance_ASV_beau.R, which had all
# accumulated overlapping/duplicate content (the "_New"/"_beau" abundance
# scripts turned out to contain the same Shannon-diversity plot as this one,
# not abundance content, likely a copy-paste mixup during drafting).

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(agricolae)

data <- read.csv("TP2_SeqTab_Meta.csv")

# Diversity indices
data$Shannon <- diversity(data[, 8:3852], index = "shannon")
data$Simpson <- diversity(data[, 8:3852], index = "simpson")
data$Total_Abundance <- rowSums(data[, 8:3852])

# ANOVA + Tukey HSD for each metric, by pesticide treatment
anova_shannon <- aov(Shannon ~ Pesticide * Pathogen + Block, data = data)
summary(anova_shannon)
tukey_shannon <- HSD.test(anova_shannon, "Pesticide", group = TRUE)
print(tukey_shannon$groups)

anova_simpson <- aov(Simpson ~ Pesticide * Pathogen + Block, data = data)
summary(anova_simpson)
tukey_simpson <- HSD.test(anova_simpson, "Pesticide", group = TRUE)
print(tukey_simpson$groups)

anova_total_abundance <- aov(Total_Abundance ~ Pesticide * Pathogen + Block, data = data)
summary(anova_total_abundance)
tukey_total_abundance <- HSD.test(anova_total_abundance, "Pesticide", group = TRUE)
print(tukey_total_abundance$groups)

# Summary tables with Tukey group letters
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

# Total abundance by pesticide
ggplot(data, aes(x = Pesticide, y = Total_Abundance)) +
  geom_boxplot(aes(fill = Pesticide), outlier.shape = NA, width = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Effect of Pesticide on Total Abundance", x = "Pesticide Treatment", y = "Total Abundance") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(0, 20000), expand = c(0, 0)) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, color = "black")

# Shannon diversity by pesticide, presentation styling (jitter + mean CI errorbar)
ggplot(data, aes(x = Pesticide, y = Shannon)) +
  geom_boxplot(aes(fill = Pesticide), outlier.shape = NA, width = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Effect of Pesticide on Shannon Diversity", x = "Pesticide Treatment", y = "Shannon Diversity Index") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, color = "black")
