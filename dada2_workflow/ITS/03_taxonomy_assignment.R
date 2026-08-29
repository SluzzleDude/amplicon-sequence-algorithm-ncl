# ITS taxonomy assignment, following the standard DADA2 ITS workflow:
# https://benjjneb.github.io/dada2/ITS_workflow.html
#
# assignTaxonomy() implements the RDP naive Bayesian classifier (Wang et al.
# 2007) to assign each ASV a taxonomy against a reference training set.
# tryRC = TRUE also checks the reverse-complement orientation, since ITS
# reads are not always oriented consistently.
#
# The UNITE reference database is not included in this repository. Download
# it from https://unite.ut.ee/repository.php and set the path below.
# Continues on from 02_dada2_pipeline.R (uses seqtab.nochim from that step).

library(dada2)

unite_ref <- "path/to/sh_general_release_dynamic_s_all.fasta"

taxa <- assignTaxonomy(seqtab.nochim, unite_ref, multithread = FALSE, tryRC = TRUE)

saveRDS(taxa, "taxa_ITS.rds")
