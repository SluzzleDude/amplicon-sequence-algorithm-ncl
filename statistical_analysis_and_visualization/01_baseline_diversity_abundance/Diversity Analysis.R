# NOTE: this script's original input file, Seqtab_Linked_TP2.csv, was not
# among the recovered project files and is not included in this repository.
# The expected columns are SampleID/Treatment/Block plus ASV_* count columns.
#
# Shannon and Simpson diversity by treatment: ANOVA, Tukey HSD, and boxplot +
# bar chart outputs. Merged with "Mean Comparison - Diversity.R", which
# repeated the same diversity calculation just to add Tukey HSD group
# letters to the plots.

library(vegan)
library(ggplot2)
library(dplyr)
library(agricolae)

file_path <- "Seqtab_Linked_TP2.csv"
data <- read.csv(file_path)

data$Shannon <- diversity(data[, grep("ASV", names(data))], index = "shannon")
data$Simpson <- diversity(data[, grep("ASV", names(data))], index = "simpson")

anova_shannon <- aov(Shannon ~ Treatment + Block, data = data)
summary(anova_shannon)
tukey_shannon <- HSD.test(anova_shannon, "Treatment", group = TRUE)

anova_simpson <- aov(Simpson ~ Treatment + Block, data = data)
summary(anova_simpson)
tukey_simpson <- HSD.test(anova_simpson, "Treatment", group = TRUE)

shannon_table <- data %>%
  group_by(Treatment) %>%
  summarise(avg_shannon = mean(Shannon)) %>%
  left_join(tibble(Treatment = rownames(tukey_shannon$groups), letter = tukey_shannon$groups$groups), by = "Treatment")
print(shannon_table)

simpson_table <- data %>%
  group_by(Treatment) %>%
  summarise(avg_simpson = mean(Simpson)) %>%
  left_join(tibble(Treatment = rownames(tukey_simpson$groups), letter = tukey_simpson$groups$groups), by = "Treatment")
print(simpson_table)

# Boxplots (no statistical grouping)
ggplot(data, aes(x = Treatment, y = Shannon, fill = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Shannon Diversity Index by Treatment", y = "Shannon Diversity Index", x = "Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("shannon_diversity_plot.png")

ggplot(data, aes(x = Treatment, y = Simpson, fill = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Simpson Diversity Index by Treatment", y = "Simpson Diversity Index", x = "Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("simpson_diversity_plot.png")

# Bar charts with Tukey HSD group letters
ggplot(data, aes(x = Treatment, y = Shannon)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge()) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2, position = position_dodge(0.9)) +
  geom_text(data = shannon_table, aes(y = avg_shannon + 0.5, label = letter), vjust = -0.5) +
  theme_minimal() +
  labs(title = "Shannon Diversity Index by Treatment", y = "Mean Shannon Diversity Index", x = "Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("shannon_diversity_bar_chart.png")

ggplot(data, aes(x = Treatment, y = Simpson)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge()) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2, position = position_dodge(0.9)) +
  geom_text(data = simpson_table, aes(y = avg_simpson + 0.5, label = letter), vjust = -0.5) +
  theme_minimal() +
  labs(title = "Simpson Diversity Index by Treatment", y = "Mean Simpson Diversity Index", x = "Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("simpson_diversity_bar_chart.png")
