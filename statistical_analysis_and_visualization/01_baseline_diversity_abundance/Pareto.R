# NOTE: this script's original input file, Seqtab_Pareto.csv, was not among
# the recovered project files and is not included in this repository. It is
# expected to be a sequence table with samples as rows and ASVs as columns
# (ASV counts starting from column 2).

library(dplyr)

file_path <- "Seqtab_Pareto.csv"
df <- read.csv(file_path)

# Total count per ASV across all samples, sorted descending
asv_sums <- colSums(df[, 2:ncol(df)])
sorted_asv_sums_df <- data.frame(ASV = names(asv_sums), TotalSum = asv_sums) %>%
  arrange(desc(TotalSum))
write.csv(sorted_asv_sums_df, "sorted_asv_sums.csv", row.names = FALSE)

# ASVs making up the top 80% of total counts (Pareto cutoff)
sorted_asv_sums_df <- sorted_asv_sums_df %>%
  mutate(CumulativeSum = cumsum(TotalSum),
         CumulativePercentage = CumulativeSum / sum(TotalSum) * 100)

important_asvs <- sorted_asv_sums_df %>% filter(CumulativePercentage <= 80)
print(important_asvs)
write.csv(important_asvs, "important_asvs_80_percent.csv", row.names = FALSE)
