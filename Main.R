library(ape)
library(seqinr)

#seq1 <- read.GenBank("KU602083")

#write.dna(seq1, file = "seq1.fasta", format="fasta")

#seq1_format <- read.fasta(file = "seq1.fasta", seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE)
#head(seq1_format)

translate.fasta<- function(dfile){
  seqchar<- gsub(" ","",dfile, fixed = TRUE)
  seqvect<- substring(seqchar,seq(1,nchar(seqchar),1), seq(1,nchar(seqchar),1))
  protvect<- translate(seqvect)
  protchar<- paste(protvect,collapse = "")
}

wprotfile<- function(accnum, pchar, fname){
  txt<-c(paste(">",accnum,sep = ""),pchar)
  writeLines(txt, fname)
}



for(seq1 in paste("KU", seq(602083, 602085,2),sep="")){
  
  filename=paste(seq1,".fasta",sep = "")
  seqread<- read.GenBank(seq1)
  write.dna(seqread,file=filename, format = "fasta")
  seqfile<- read.fasta(filename, as.string = TRUE, seqonly = TRUE)
  file.rename(from = paste("~/BioInf/",filename,sep = ""), to = paste("~/BioInf/VH DNA/",filename,sep = ""))
  
  protstr<- translate.fasta(seqfile)
  wprotfile(seq1, protstr, filename)
  file.rename(from = paste("~/BioInf/",filename,sep = ""), to = paste("~/BioInf/VH Prot/",filename,sep = ""))

}

##test##602083, 602723

##test2##