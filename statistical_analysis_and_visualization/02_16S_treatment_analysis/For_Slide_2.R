# Shannon diversity by pesticide treatment, TP3 (16S) - alternate presentation
# variant to For_Slide.R: pairwise significance brackets, then error bars.

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(agricolae)
library(ggsignif)

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

data_summary <- data %>%
  group_by(Pesticide) %>%
  summarise(mean_Shannon = mean(Shannon), se_Shannon = sd(Shannon) / sqrt(n()))

fill_colors <- c("Control" = "lightblue", "Metalaxyl" = "blue",
                  "K61" = "dodgerblue", "Root_Shield" = "darkblue")

plot_theme <- theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

# Boxplot with pairwise significance brackets
p_signif <- ggplot(data, aes(x = Pesticide, y = Shannon, fill = Pesticide)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Effect of Pesticide on Shannon Diversity",
       x = "Pesticide Treatment", y = "Shannon Diversity Index") +
  plot_theme +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = fill_colors) +
  geom_signif(comparisons = list(c("Metalaxyl", "Control"), c("Metalaxyl", "K61"), c("Metalaxyl", "Root_Shield"),
                                  c("Control", "K61"), c("Control", "Root_Shield"), c("K61", "Root_Shield")),
              map_signif_level = TRUE, y_position = c(3.0, 3.2, 3.4, 3.6, 3.8, 4.0))

print(p_signif)

# Boxplot with mean +/- standard error bars
p_errorbar <- ggplot(data, aes(x = Pesticide, y = Shannon, fill = Pesticide)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  geom_errorbar(data = data_summary, aes(x = Pesticide, y = mean_Shannon,
                                          ymin = mean_Shannon - se_Shannon, ymax = mean_Shannon + se_Shannon),
                width = 0.2, color = "black") +
  labs(title = "Effect of Pesticide on Shannon Diversity",
       x = "Pesticide Treatment", y = "Shannon Diversity Index") +
  plot_theme +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = fill_colors)

print(p_errorbar)
