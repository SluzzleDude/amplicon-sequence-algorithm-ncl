# 16S taxonomy assignment, following the standard DADA2 tutorial:
# https://benjjneb.github.io/dada2/tutorial.html
#
# assignTaxonomy() implements the RDP naive Bayesian classifier (Wang et al.
# 2007): it assigns each ASV a taxonomy by bootstrapping k-mer frequency
# comparisons against a reference training set. addSpecies() then does exact
# string matching against a species-level reference to add species-level
# calls where possible.
#
# Reference training sets are not included in this repository (they are
# several hundred MB). Download the Silva reference files from
# https://benjjneb.github.io/dada2/training.html and set the paths below.
# Continues on from 03_dada2_pipeline.R (uses seqtab.nochim from that step).

library(dada2)

silva_train_set <- "path/to/silva_nr99_v138.2_toGenus_trainset.fa.gz"
silva_species_assignment <- "path/to/silva_v138.2_assignSpecies.fa.gz"

taxa <- assignTaxonomy(seqtab.nochim, silva_train_set, multithread = FALSE)
taxa <- addSpecies(taxa, silva_species_assignment)

saveRDS(taxa, "taxa_16S.rds")
