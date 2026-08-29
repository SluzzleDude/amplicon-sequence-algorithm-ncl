# Total ASV abundance by pesticide and pathogen treatment, TP3 (16S).
#
# Trimmed from the original TP3_Abundance_ASV.R, which had accumulated two
# throwaway early boxplots and two trailing blocks that recreated the
# summary stats from hardcoded numbers instead of computing them (leftover
# drafting artifacts). Also merges in the presentation-styled boxplot from
# TP3_Abundance_ASV_beau.R.

library(dplyr)
library(ggplot2)
library(tidyr)
library(multcomp)
library(agricolae)

data <- read.csv("TP3_SeqTab_Meta.csv")
data$Total_Abundance <- rowSums(data[, grep("^ASV", names(data))])

# ANOVA + Tukey HSD by pesticide treatment
anova_model <- aov(Total_Abundance ~ Pesticide * Pathogen + Block, data = data)
summary(anova_model)
tukey_results <- HSD.test(anova_model, "Pesticide", group = TRUE)
print(tukey_results)

summary_table <- tukey_results$groups %>%
  rownames_to_column(var = "Pesticide") %>%
  rename(Average_Abundance = means, Tukey_Group = groups)
write.csv(summary_table, "TP3_Pesticide_Abundance_Tukey_Summary.csv", row.names = FALSE)

# Boxplot with Tukey group letters, by pesticide
ggplot(data, aes(x = Pesticide, y = Total_Abundance)) +
  geom_boxplot() +
  geom_text(data = summary_table, aes(x = Pesticide, y = Average_Abundance + 5000, label = Tukey_Group), vjust = -0.5, color = "red") +
  labs(title = "Total Abundance by Pesticide", x = "Pesticide", y = "Total Abundance") +
  theme_minimal()
ggsave("TP3_Pesticide_Abundance_Boxplot.png")

# Boxplot by pathogen presence
ggplot(data, aes(x = Pathogen, y = Total_Abundance)) +
  geom_boxplot() +
  labs(title = "Total Abundance by Pathogen Presence", x = "Pathogen", y = "Total Abundance") +
  theme_minimal()
ggsave("TP3_Pathogen_Abundance_Boxplot.png")

# Bar chart with error bars and Tukey group letters, by pesticide
pesticide_summary <- data %>%
  group_by(Pesticide) %>%
  summarise(Average_Abundance = mean(Total_Abundance), SD = sd(Total_Abundance)) %>%
  left_join(summary_table, by = "Pesticide")

ggplot(pesticide_summary, aes(x = Pesticide, y = Average_Abundance)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Average_Abundance - SD, ymax = Average_Abundance + SD), width = 0.2) +
  geom_text(aes(label = Tukey_Group, y = Average_Abundance + SD + 5000), vjust = -0.5, color = "red") +
  labs(title = "Average Total Abundance by Pesticide", x = "Pesticide", y = "Average Total Abundance") +
  theme_minimal()
ggsave("TP3_Pesticide_Abundance_BarChart.png")

# Bar chart with error bars, by pathogen presence
pathogen_summary <- data %>%
  group_by(Pathogen) %>%
  summarise(Average_Abundance = mean(Total_Abundance), SD = sd(Total_Abundance))

ggplot(pathogen_summary, aes(x = Pathogen, y = Average_Abundance)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Average_Abundance - SD, ymax = Average_Abundance + SD), width = 0.2) +
  labs(title = "Average Total Abundance by Pathogen Presence", x = "Pathogen", y = "Average Total Abundance") +
  theme_minimal()
ggsave("TP3_Pathogen_Abundance_BarChart.png")

# Presentation-styled boxplot, by pesticide (y-axis capped at 25000 to keep
# outliers from compressing the comparison)
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
  ylim(0, 25000) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, color = "black")
