# Load necessary libraries
library(dplyr)
library(ggplot2)

# Read the data from the CSV file
file_path <- "ITS_Seqab_Metadata_TP3.csv"
df <- read.csv(file_path)

# Define the range of columns that contain ASV data
asv_columns <- 8:(8 + 7949 - 1)

# Calculate the sum for each ASV across all samples
asv_sums <- colSums(df[ , asv_columns])

# Convert to a data frame for easier manipulation
asv_sums_df <- data.frame(ASV = names(asv_sums), TotalSum = asv_sums)

# Sort the ASVs by the total sum in decreasing order
sorted_asv_sums_df <- asv_sums_df %>% arrange(desc(TotalSum))

# Plot the top N ASVs, for example, top 50
top_n <- 50
top_asvs <- head(sorted_asv_sums_df, n = top_n)

ggplot(top_asvs, aes(x = reorder(ASV, -TotalSum), y = TotalSum)) +
  geom_bar(stat = "identity", fill = "blue") +
  scale_fill_manual(values = colorRampPalette(c("lightblue", "blue"))(top_n)) +
  xlab("ASV") +
  ylab("Total Sum") +
  ggtitle(paste("Top", top_n, "ASVs by Total Sum")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

############################
############################
############################

# Calculate the average count for each ASV per sample
asv_means <- colMeans(df[ , asv_columns])

# Convert to a data frame for easier manipulation
asv_means_df <- data.frame(ASV = names(asv_means), AverageCount = asv_means)

# Sort the ASVs by the average count in decreasing order
sorted_asv_means_df <- asv_means_df %>% arrange(desc(AverageCount))

# Plot the top N ASVs by average count per sample, for example, top 50
top_n <- 50
top_asvs <- head(sorted_asv_means_df, n = top_n)

ggplot(top_asvs, aes(x = reorder(ASV, -AverageCount), y = AverageCount)) +
  geom_bar(stat = "identity", fill = "blue") +
  scale_fill_manual(values = colorRampPalette(c("lightblue", "blue"))(top_n)) +
  xlab("ASV") +
  ylab("Average Count per Sample") +
  ggtitle(paste("Top", top_n, "ASVs by Average Count per Sample")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

