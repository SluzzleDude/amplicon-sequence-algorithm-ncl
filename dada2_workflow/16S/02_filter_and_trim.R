# 16S filter and trim step.
#
# The full dissertation dataset (191 paired samples) is not included in this
# repository. Set `path` below to input_files/16S_example for the bundled
# 4-file example, or to your own local copy of the full raw dataset.

library(dada2)

path <- "input_files/16S_example"

fnFs <- sort(list.files(path, pattern = "_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq", full.names = TRUE))

if (length(fnFs) != length(fnRs)) {
  stop("The number of forward and reverse files do not match.")
}

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Place filtered files in a filtered/ subdirectory
filt_path <- file.path(path, "filtered")
if (!dir.exists(filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))

# Truncation lengths chosen from the quality profiles inspected in
# 01_quality_inspection.R
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(240, 180),
                      maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
                      compress = TRUE, multithread = FALSE) # set TRUE on Linux/macOS
head(out)
