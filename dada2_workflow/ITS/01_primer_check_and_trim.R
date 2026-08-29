# ITS primer removal, following the standard DADA2 ITS workflow:
# https://benjjneb.github.io/dada2/ITS_workflow.html
#
# ITS amplicons vary in length (unlike 16S), so primers are removed with
# cutadapt rather than a fixed truncLen in filterAndTrim. This requires
# cutadapt to be installed separately (https://cutadapt.readthedocs.io) and
# available on the system PATH, or pointed at explicitly below.
#
# The full dissertation dataset (191 paired samples) is not included in this
# repository. Set `path` below to input_files/ITS_example for the bundled
# 4-file example, or to your own local copy of the full raw dataset.

library(dada2)
library(ShortRead)
library(Biostrings)

path <- "input_files/ITS_example"
cutadapt <- "cutadapt" # or the full path to the cutadapt executable

fnFs <- sort(list.files(path, pattern = "_1.fastq$", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq$", full.names = TRUE))

# ITS1/ITS2 primer pair used in the original sequencing run
FWD <- "CTTGGTCATTTAGAGGAAGTAA"
REV <- "GCTGCGTTCATCGATGC"

# Count how many reads contain each primer, in every orientation, before
# removal (sanity check recommended by the tutorial)
allOrients <- function(primer) {
  dna <- DNAString(primer)
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna),
               Reverse = Biostrings::reverse(dna), RevComp = Biostrings::reverseComplement(dna))
  sapply(orients, toString)
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  sum(nhits > 0)
}

# Pre-filter reads containing ambiguous bases, which cutadapt/dada2 cannot
# handle, before searching for primers
fnFs.filtN <- file.path(path, "filtN", basename(fnFs))
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = FALSE)

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

# Remove primers with cutadapt
path.cut <- file.path(path, "cutadapt")
if (!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
R2.flags <- paste("-G", REV, "-A", FWD.RC)

for (i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,
                              "-o", fnFs.cut[i], "-p", fnRs.cut[i],
                              fnFs.filtN[i], fnRs.filtN[i]))
}

# Confirm primers were removed
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
