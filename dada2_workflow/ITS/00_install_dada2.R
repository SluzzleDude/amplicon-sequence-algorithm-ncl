# Install DADA2 and the supporting packages used throughout the ITS workflow.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dada2")
BiocManager::install("ShortRead")
BiocManager::install("Biostrings")
install.packages("RcppParallel")
