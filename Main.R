library(Biostrings)
library(ape)
library(seqinr)
library(xlsx)

#seq1 <- read.GenBank("KU602083")

#write.dna(seq1, file = "seq1.fasta", format="fasta")

#seq1_format <- read.fasta(file = "seq1.fasta", seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE)
#head(seq1_format)

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

#initializes text file
sink("VHmasteracc.txt")
sink()

accession.num<-paste("KU", seq(602083, 602087,2),sep="")

#pulls sequences from NCBI from a matrix of accesion number strings
for(seq1 in accession.num){
  
  filename=paste(seq1,".fasta",sep = "")
  seqread<- read.GenBank(seq1)
  write.dna(seqread,file=filename, format = "fasta")
  #seqfile<- read.fasta(filename, as.string = TRUE, seqonly = TRUE)
  seqfile<- read.string(filename)
  file.rename(from = paste("~/BioInf/",filename,sep = ""), to = paste("~/BioInf/VH DNA/",filename,sep = ""))
  
  protstr<- translate.fasta(seqfile)
  wprotfile(seq1, protstr, filename)
  file.rename(from = paste("~/BioInf/",filename,sep = ""), to = paste("~/BioInf/VH Prot/",filename,sep = ""))
  
  write(seq1, file = "VHmasteracc.txt", append = TRUE)
}

#x<-readLines("VHmasteracc.txt") #read file as vector of chars
#x<- readLines("~/BioInf/VH DNA/KU602083.fasta")
#y<- readLines("~/BioInf/VH DNA/KU602085.fasta")

directory.seq<- function(single.acc.num){
  paste("~/BioInf/VH DNA/", single.acc.num, ".fasta", sep="")
}

masteracc.num <- readLines("VHmasteracc.txt")
xlsize=length(masteracc.num)+1 #number of sequences

xlwb    <- createWorkbook(type="xlsx")           # create an empty workbook
sheet <- createSheet(xlwb, sheetName="Sheet1")   # create an empty sheet 
rows  <- createRow(sheet, rowIndex=1:xlsize)      #rows
cells <- createCell(rows, colIndex=1:xlsize)      #columns

data("BLOSUM62")

for (i in 1:length(masteracc.num)){
  setCellValue(cells[[1+i,1]], masteracc.num[i])
  setCellValue(cells[[1,1+i]], masteracc.num[i])
  for (j in 2:length(masteracc.num)){
    s1=read.string(directory.seq(masteracc.num[i]))
    s2=read.string(directory.seq(masteracc.num[j]))
    localAlign <- pairwiseAlignment(s1,s2, substitutionMatrix=BLOSUM62, gapOpening=5, gapExtension=2, scoreOnly=TRUE)
    setCellValue(cells[[1+i,1+j]], localAlign)
    setCellValue(cells[[1+j,1+i]], localAlign)
  }
}

saveWorkbook(xlwb, "pscores.xlsx")

##test##602083, 602723

##test2##