# 16S quality inspection: load paired reads, verify forward/reverse pairing,
# check for empty files, and plot quality profiles.
#
# The full dissertation dataset (191 paired samples) is not included in this
# repository. Set `path` below to input_files/16S_example for the bundled
# 4-file example, or to your own local copy of the full raw dataset.

library(ShortRead)
library(dada2)
library(BiocParallel)

path <- "input_files/16S_example"

# Forward and reverse fastq filenames
fnFs <- sort(list.files(path, pattern = "_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq", full.names = TRUE))

cat("Forward reads files:\n"); print(fnFs)
cat("Reverse reads files:\n"); print(fnRs)

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
cat("Sample names:\n"); print(sample.names)

# Function to check if a fastq file is empty
is_empty_fastq <- function(file) {
  fastqFile <- readFastq(file)
  length(sread(fastqFile)) == 0
}

# Filter out empty fastq files
non_empty_fnFs <- fnFs[!sapply(fnFs, is_empty_fastq)]
non_empty_fnRs <- fnRs[!sapply(fnRs, is_empty_fastq)]

# Ensure the paired files are correctly matched by sample name
paired_fnFs <- non_empty_fnFs[sample.names %in% sapply(strsplit(basename(non_empty_fnRs), "_"), `[`, 1)]
paired_fnRs <- non_empty_fnRs[sample.names %in% sapply(strsplit(basename(non_empty_fnFs), "_"), `[`, 1)]

if (length(paired_fnFs) != length(paired_fnRs)) {
  stop("The number of forward and reverse files does not match. Please check your file pairs.")
}

cat("Non-empty, paired forward fastq files:\n"); print(paired_fnFs)
cat("Non-empty, paired reverse fastq files:\n"); print(paired_fnRs)

# Inspect read quality profiles
register(SerialParam())
plotQualityProfile(paired_fnFs)
plotQualityProfile(paired_fnRs)
