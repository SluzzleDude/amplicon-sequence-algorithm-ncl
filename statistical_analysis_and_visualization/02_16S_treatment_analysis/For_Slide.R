# Shannon diversity by pesticide treatment, TP3 (16S).
# Presentation-ready plot iterations: boxplot, then a bar chart with error
# bars and Tukey HSD group letters.

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(agricolae)

# Load data
data <- read.csv("TP3_SeqTab_Meta.csv")

# Calculate Shannon and Simpson diversity indices
data$Shannon <- diversity(data[, 8:3852], index = "shannon")
data$Simpson <- diversity(data[, 8:3852], index = "simpson")

# ANOVA and Tukey HSD for Shannon diversity by pesticide treatment
anova_shannon <- aov(Shannon ~ Pesticide * Pathogen + Block, data = data)
summary(anova_shannon)

tukey_shannon <- HSD.test(anova_shannon, "Pesticide", group = TRUE)
print(tukey_shannon$groups)

tukey_letters <- data.frame(Pesticide = rownames(tukey_shannon$groups),
                             Letter = tukey_shannon$groups$groups)

fill_colors <- c("Control" = "lightblue", "Metalaxyl" = "blue",
                  "K61" = "dodgerblue", "Root_Shield" = "darkblue")

plot_theme <- theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

# Boxplot with per-point jitter and Tukey group letters
p_boxplot <- ggplot(data, aes(x = Pesticide, y = Shannon)) +
  geom_boxplot(aes(fill = Pesticide), outlier.shape = NA, width = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Effect of Pesticide on Shannon Diversity",
       x = "Pesticide Treatment", y = "Shannon Diversity Index") +
  plot_theme +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = fill_colors) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, fill = "white") +
  geom_text(data = tukey_letters, aes(x = Pesticide, y = max(data$Shannon) + 0.1, label = Letter), vjust = 0)

print(p_boxplot)

# Bar chart with mean +/- standard error and Tukey group letters
data_summary <- data %>%
  group_by(Pesticide) %>%
  summarise(mean_Shannon = mean(Shannon), se_Shannon = sd(Shannon) / sqrt(n())) %>%
  merge(tukey_letters, by = "Pesticide")

p_barchart <- ggplot(data_summary, aes(x = Pesticide, y = mean_Shannon, fill = Pesticide)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = mean_Shannon - se_Shannon, ymax = mean_Shannon + se_Shannon), width = 0.2, color = "black") +
  geom_text(aes(y = mean_Shannon + se_Shannon + 0.2, label = Letter), vjust = 0) +
  labs(title = "Effect of Pesticide on Shannon Diversity",
       x = "Pesticide Treatment", y = "Shannon Diversity Index") +
  plot_theme +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = fill_colors)

print(p_barchart)
