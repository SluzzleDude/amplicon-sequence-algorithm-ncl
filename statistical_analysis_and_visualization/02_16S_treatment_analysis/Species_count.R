# Load necessary libraries
library(vegan)
library(dplyr)

# Load data
data <- read.csv("TP2_SeqTab_Meta.csv")

# Calculate the count of species per sample
data$SpeciesCount <- specnumber(data[, 8:3852])

# Create a new dataframe with SpeciesCount and Pooled_number
new_data <- data %>%
  select(SpeciesCount, Pooled_number)

# Save the new dataframe to a CSV file
write.csv(new_data, "SpeciesCount_PooledNumber.csv", row.names = FALSE)

# View the first few rows of the new data
head(new_data)
