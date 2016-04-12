library(ape)
library(seqinr)

seq1 <- read.GenBank("NC_001477")

write.dna(seq1, file = "seq1.fasta", format="fasta")

seq1_format <- read.fasta(file = "seq1.fasta", seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE)
head(seq1_format)
