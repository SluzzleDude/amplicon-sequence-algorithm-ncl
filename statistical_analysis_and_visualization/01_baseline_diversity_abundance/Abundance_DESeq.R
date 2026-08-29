# NOTE: this script's original input file, Seqtab_Linked_TP2.csv, was not
# among the recovered project files and is not included in this repository.
# The expected columns are SampleID/Treatment/Block plus ASV_* count columns.

# Load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
BiocManager::install("DESeq2")

library(phyloseq)
library(DESeq2)
library(dplyr)

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

# Handle zero counts by adding a small constant (e.g., 1)
otu_table <- otu_table + 1

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

# Convert phyloseq object to DESeq2 object
dds <- phyloseq_to_deseq2(physeq, ~ Block + Treatment)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Set "Control_pyth" as the reference level for the Treatment factor
dds$Treatment <- relevel(dds$Treatment, ref = "Control_pyth")

# Extract results for the treatment effect
results_names <- resultsNames(dds)
for (i in results_names) {
  if (i != "Intercept") {
    res <- results(dds, name = i)
    # Order results by adjusted p-value
    res <- res[order(res$padj),]
    
    # Display the results
    print(paste("Results for comparison with:", i))
    print(res)
    
    # Summarize results
    summary_res <- summary(res)
    print(summary_res)
    
    # Plot results (MA plot)
    plotMA(res, main = paste("DESeq2: ", i), ylim = c(-2, 2))
    
    # Save results to CSV
    write.csv(as.data.frame(res), file = paste("DESeq2_results_", i, ".csv", sep = ""))
  }
}
