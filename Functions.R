#read fasta sequence as a single string of uppercase letters
read.string<- function(fastafile){
  a<-toupper(read.fasta(fastafile, seqonly = TRUE))
  seqchar<- gsub(" ","",a, fixed = TRUE)
}

#read fasta sequence as a matrix of one character strings
read.single<- function(fastafile){
  seqchar<- read.string(fastafile)
  seqvect<- substring(seqchar,seq(1,nchar(seqchar),1), seq(1,nchar(seqchar),1))
}

#translates a single string of UPPERCASE letters to one letter AA codes
translate.fasta<- function(dfile){
  seqchar<- gsub(" ","",dfile, fixed = TRUE)
  seqvect<- substring(seqchar,seq(1,nchar(seqchar),1), seq(1,nchar(seqchar),1))
  protvect<- translate(seqvect)
  protchar<- paste(protvect,collapse = "")
}

#Writes protein sequence to file, AC# in first line
wprotfile<- function(accnum, pchar, fname){
  txt<-c(paste(">",accnum,sep = ""),pchar)
  writeLines(txt, fname)
}

directory.seq<- function(single.acc.num){
  paste("~/BioInf/VH Prot/", single.acc.num, ".fasta", sep="")
}
