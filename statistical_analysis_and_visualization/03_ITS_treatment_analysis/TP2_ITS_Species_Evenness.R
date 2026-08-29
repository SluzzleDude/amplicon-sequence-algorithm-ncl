# Load necessary libraries
library(dplyr)
library(vegan)  # Needed for diversity() function

# Set file path and file name
file_path <- "ITS_Seqab_Metadata_TP2.csv"

# Read the data from the specified file path
data <- read.csv(file_path)

# Define the columns containing the ASVs (columns 8 to 7956)
asv_data <- data[, 8:7956]  # 8 is the first ASV column, 7956 is the 8th column + 7949 ASVs

# Calculate the Shannon Diversity Index for each sample
shannon_diversity <- diversity(asv_data, index = "shannon")

# Calculate the mean Shannon Diversity Index
mean_H <- mean(shannon_diversity)

# Determine the number of unique ASVs (S)
S <- 7949  # Total number of ASVs

# Calculate the Shannon Evenness Index
E <- mean_H / log(S)

# Print the result
print(paste("The overall fungal species evenness is:", E))
