# NOTE: this script's original input file, Shannon_Diversity_Index.csv, was
# not among the recovered project files and is not included in this
# repository. It is expected to have a Shannon_Diversity column with one row
# per sample.

library(dplyr)

file_path <- "Shannon_Diversity_Index.csv"
data <- read.csv(file_path)

# Calculate the mean Shannon Diversity Index
mean_H <- mean(data$Shannon_Diversity)

# Determine the number of unique ASVs (S)
# Assuming you have this information, for example, if S = 100:
S <- 3052  # Replace this with the actual number of unique ASVs if known

# Calculate the Shannon Evenness Index
E <- mean_H / log(S)

# Print the result
print(paste("The overall bacterial species evenness is:", E))
