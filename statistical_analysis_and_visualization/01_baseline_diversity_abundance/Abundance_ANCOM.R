# NOTE: this script's original input file, Seqtab_Linked_TP2.csv, was not
# among the recovered project files and is not included in this repository.
# The expected columns are SampleID/Treatment/Block plus ASV_* count columns.

# Set the file path
file_path <- "Seqtab_Linked_TP2.csv"

# Load the data
data <- read.csv(file_path)

# Check the column names
print(colnames(data))

# Verify that columns for SampleID, Treatment, and Block exist
# If your dataset does not have a SampleID column, you can create one
if (!"SampleID" %in% colnames(data)) {
  data$SampleID <- paste0("Sample_", 1:nrow(data))
}

# Convert Treatment and Block to factors
data$Treatment <- as.factor(data$Treatment)
data$Block <- as.factor(data$Block)

# Create OTU table
otu_table <- data %>%
  select(starts_with("ASV")) %>%
  as.matrix()
rownames(otu_table) <- data$SampleID  # Assuming your dataset has a SampleID column

# Create sample metadata
sample_metadata <- data %>%
  select(SampleID, Treatment, Block) %>%
  as.data.frame()
rownames(sample_metadata) <- sample_metadata$SampleID
sample_metadata <- sample_metadata %>% select(-SampleID)

# Create phyloseq object
OTU = otu_table(otu_table, taxa_are_rows = FALSE)
SAMPLES = sample_data(sample_metadata)
physeq = phyloseq(OTU, SAMPLES)

###############################

# Perform ANCOMBC2 analysis
ancombc2_out <- ancombc2(data = physeq, 
                         fix_formula = "Treatment + Block", 
                         p_adj_method = "holm", 
                         struc_zero = TRUE, 
                         neg_lb = TRUE, 
                         alpha = 0.05, 
                         global = TRUE)

################################
# Extract results
res <- ancombc2_out$res
print(res)

# Extracting relevant results for ASVs
beta <- res$beta
p_val <- res$p_val
q_val <- res$q_val

# Combine results into a single dataframe
results <- cbind(beta, p_val, q_val)

# Filter significant ASVs based on adjusted p-values
sig_asvs <- results[results$q_val < 0.05, ]
print(sig_asvs)

# Plot significant ASVs
sig_asvs_df <- as.data.frame(sig_asvs)
ggplot(sig_asvs_df, aes(x = rownames(sig_asvs_df), y = beta, fill = rownames(sig_asvs_df))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(title = "Significant ASVs by Treatment (ANCOMBC2)", x = "ASV", y = "log2 Fold Change") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
