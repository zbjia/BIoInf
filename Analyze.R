source("Functions.R")
library(msa)

setwd("~/BioInf/")

x<-read.csv("scoresVHProtBLOSUM100.csv", row.names = 1)
y<-x[237] #column number
z<-y[order(y, decreasing=FALSE),,drop=FALSE]

#Generate top 10 sequences to align
seqnames<-rownames(z)[2:11]
msa.seq.filenames<-sapply(seqnames, FUN=function(thing){paste(thing, ".fasta", sep="")},USE.NAMES = FALSE)
msa.seq.filenames<-c(msa.seq.filenames, "KU602555.fasta")

setwd("~/BioInf/VH Prot/")
msa.seq<- readAAStringSet(msa.seq.filenames)
msAlign<- msa(msa.seq, gapOpening = 5, gapExtension = 2) #MSA using ClustalW default parameters

print(msAlign, show = "complete")
setwd("~/BioInf/")

#write.table(msa.seq.filenames, file="topmatrices.txt", append = TRUE, row.names = FALSE, col.names = "BLOSUM62")
#head(z,11)