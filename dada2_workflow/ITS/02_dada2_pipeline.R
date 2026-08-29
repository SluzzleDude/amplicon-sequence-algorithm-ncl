# ITS DADA2 pipeline: quality inspection, filter/trim, error learning,
# denoising, paired-read merging, chimera removal, and sequence-table export,
# following the standard DADA2 ITS workflow:
# https://benjjneb.github.io/dada2/ITS_workflow.html
#
# Continues on from 01_primer_check_and_trim.R (uses the primer-trimmed reads
# in path.cut). Note there is no truncLen here, unlike the 16S pipeline: ITS
# amplicon length varies naturally, so reads are filtered on quality
# (maxEE/truncQ) and a minimum length only.

library(dada2)

path <- "input_files/ITS_example"
path.cut <- file.path(path, "cutadapt")

cutFs <- sort(list.files(path.cut, pattern = "_1.fastq$", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_2.fastq$", full.names = TRUE))

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))

plotQualityProfile(cutFs)
plotQualityProfile(cutRs)

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2,
                      minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE) # set TRUE on Linux/macOS
head(out)

errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)
plotErrors(errF, nominalQ = TRUE)

dadaFs <- dada(filtFs, err = errF, multithread = FALSE)
dadaRs <- dada(filtRs, err = errR, multithread = FALSE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE, verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim) / sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

seqtab.nochim.df <- as.data.frame(seqtab.nochim)
write.csv(seqtab.nochim.df, "seqtab_nochim_ITS.csv", row.names = TRUE)
saveRDS(seqtab.nochim, "seqtab_nochim_ITS.rds")
