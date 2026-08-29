# NOTE: this script's original input file, Seqtab_Pareto.csv, was not among
# the recovered project files and is not included in this repository. It is
# expected to be a sequence table with samples as rows and ASVs as columns
# (ASV counts starting from column 2).
#
# Companion to Pareto.R: visualizes ASV count distribution as a heatmap, a
# Pareto chart, a cumulative-distribution curve, and top-N bar charts.

library(dplyr)
library(ggplot2)
library(reshape2)

file_path <- "Seqtab_Pareto.csv"
df <- read.csv(file_path)

asv_matrix <- as.matrix(df[, 2:ncol(df)])

sorted_asv_sums_df <- data.frame(ASV = names(colSums(asv_matrix)), TotalSum = colSums(asv_matrix)) %>%
  arrange(desc(TotalSum)) %>%
  mutate(CumulativeSum = cumsum(TotalSum))

asv_means_df <- data.frame(ASV = names(colMeans(asv_matrix)), AverageCount = colMeans(asv_matrix)) %>%
  arrange(desc(AverageCount))

top_n <- 50

# Heatmap of ASV counts across samples
heatmap(asv_matrix, scale = "row", col = heat.colors(256), margins = c(5, 10))

# Pareto chart: ASV counts with cumulative-sum overlay
ggplot(sorted_asv_sums_df, aes(x = reorder(ASV, -TotalSum), y = TotalSum)) +
  geom_bar(stat = "identity") +
  geom_line(aes(y = CumulativeSum), group = 1, color = "blue") +
  coord_flip() +
  xlab("ASV") + ylab("Total Sum") +
  ggtitle("Pareto Chart of ASV Counts")

# Cumulative distribution of ASV counts by rank
ggplot(sorted_asv_sums_df, aes(x = seq_len(nrow(sorted_asv_sums_df)), y = CumulativeSum)) +
  geom_line() +
  xlab("ASV Rank") + ylab("Cumulative Sum of Counts") +
  ggtitle("Cumulative Distribution of ASV Counts")

# Top-N ASVs by total count
ggplot(head(sorted_asv_sums_df, n = top_n), aes(x = reorder(ASV, -TotalSum), y = TotalSum)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  xlab("ASV") + ylab("Total Sum") +
  ggtitle(paste("Top", top_n, "ASVs by Total Sum"))

# Top-N ASVs by average count per sample
ggplot(head(asv_means_df, n = top_n), aes(x = reorder(ASV, -AverageCount), y = AverageCount)) +
  geom_bar(stat = "identity", fill = "blue") +
  scale_fill_manual(values = colorRampPalette(c("lightblue", "blue"))(top_n)) +
  xlab("ASV") + ylab("Average Count per Sample") +
  ggtitle(paste("Top", top_n, "ASVs by Average Count per Sample")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Same, with standard-error bars
asv_se <- apply(asv_matrix, 2, function(x) sd(x) / sqrt(length(x)))
asv_stats_df <- data.frame(ASV = names(asv_se), AverageCount = colMeans(asv_matrix), SE = asv_se) %>%
  arrange(desc(AverageCount))

ggplot(head(asv_stats_df, n = top_n), aes(x = reorder(ASV, -AverageCount), y = AverageCount)) +
  geom_bar(stat = "identity", fill = "blue") +
  geom_errorbar(aes(ymin = AverageCount - SE, ymax = AverageCount + SE), width = 0.2) +
  scale_fill_manual(values = colorRampPalette(c("lightblue", "blue"))(top_n)) +
  xlab("ASV") + ylab("Average Count per Sample") +
  ggtitle(paste("Top", top_n, "ASVs by Average Count per Sample (+/- SE)")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
