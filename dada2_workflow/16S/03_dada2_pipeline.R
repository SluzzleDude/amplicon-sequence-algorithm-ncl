# 16S DADA2 pipeline: error learning, denoising, paired-read merging, chimera
# removal, and read tracking, following the standard DADA2 tutorial:
# https://benjjneb.github.io/dada2/tutorial.html
#
# Continues on from 02_filter_and_trim.R in the same R session (uses `out`
# from that step). Taxonomy assignment continues in 04_taxonomy_assignment.R.

library(dada2)

path <- "input_files/16S_example"

fnFs <- sort(list.files(path, pattern = "_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))

# Learn the error rates
errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

# Dereplication
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Sample inference
dadaFs <- dada(derepFs, err = errF, multithread = FALSE)
dadaRs <- dada(derepRs, err = errR, multithread = FALSE)

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE, verbose = TRUE)
dim(seqtab.nochim)

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

saveRDS(seqtab.nochim, "seqtab_nochim_16S.rds")
