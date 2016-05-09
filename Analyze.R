source("Functions.R")
library(msa)

setwd("~/BioInf/")

x<-read.csv("pscores1.csv", row.names = 1)

y<-x[1]
z<-y[order(y, decreasing=TRUE),,drop=FALSE]

seqnames<-rownames(z)[2:11]
msa.seq.filenames<-sapply(seqnames, FUN=function(thing){paste(thing, ".fasta", sep="")},USE.NAMES = FALSE)
msa.seq.filenames<-c(msa.seq.filenames, "KU602083.fasta")

setwd("~/BioInf/VH Prot/")
msa.seq<- readAAStringSet(msa.seq.filenames)
msAlign<- msa(msa.seq)

print(msAlign, show = "complete")

msaPrettyPrint(msAlign, output="pdf", showNames="none",
               showLogo="none", askForOverwrite=FALSE, verbose=FALSE)
