


ComputeStatsFromFasta <- function(FastaFilePath){
  library(seqinr)
  
  
  GoIFasta=read.fasta(file=FastaFilePath)
  
  PreMat=as.data.frame(matrix(0,nrow = length(GoIFasta),ncol=6),stringsAsFactors = F)
  names(PreMat)=c("OriginalName","ProtName","Annot","Specie","Length","Type")
  PreMat$OriginalName=unlist(getName(GoIFasta))
  PreMat$ProtName=unlist(lapply(PreMat[,1],function(x) paste(strsplit(strsplit(x,"\\$")[[1]][1],"@")[[1]][-1],collapse ="-")))
  PreMat$Assembly=unlist(lapply(PreMat[,1],function(x) strsplit(x,split = "\\$")[[1]][-1]))
  PreMat$Taxid=unlist(lapply(PreMat[,1],function(x) strsplit(x,split = "_")[[1]][1]))
  #Stupidly complicated but ok
  PreMat$Annot=unlist(lapply(PreMat[,1],function(x) paste(strsplit(strsplit(x,"\\[")[[1]][1]," ")[[1]][-1],collapse ="-")))
  PreMat$Length=getLength(GoIFasta)
  #Need to have upper case letter for some reasons. 
  PreMat$Sequence=toupper(unlist(getSequence(GoIFasta,as.string=T)))
  #Vector of new names : 
  NamesToadd=names(unlist(AAstat(s2c(PreMat$Sequence[1]),plot = F)))
  #Creating new columns
  PreMat[,NamesToadd]=NA
  #Running this one liner, the transposition is necessary, I'm not sure why.
  PreMat[,NamesToadd]=t(apply(PreMat,1,function(x) as.numeric(as.character(unlist(AAstat(s2c(x[which(names(x)=="Sequence")]),plot = F))))))
  
  PreMatExport=PreMat[,-which(names(PreMat)=="Sequence")]
  
  #End of function
}




ComputeStatsFromFasta_reducedName <- function(FastaFilePath){
  library(seqinr)
  
  
  GoIFasta=read.fasta(file=FastaFilePath)
  
  PreMat=as.data.frame(matrix(0,nrow = length(GoIFasta),ncol=4),stringsAsFactors = F)
  names(PreMat)=c("OriginalName","ProtName","Length","Type")
  PreMat$OriginalName=unlist(getName(GoIFasta))
  PreMat$ProtName=unlist(lapply(PreMat[,1], function(x) strsplit(strsplit(x,"_")[[1]][2],"\\$")[[1]][1]))
  PreMat$Assembly=unlist(lapply(PreMat[,1],function(x) strsplit(x,split = "\\$")[[1]][-1]))
  PreMat$Taxid=unlist(lapply(PreMat[,1],function(x) strsplit(x,split = "_")[[1]][1]))
  #Stupidly complicated but ok
  PreMat$Length=getLength(GoIFasta)
  #Need to have upper case letter for some reasons. 
  PreMat$Sequence=toupper(unlist(getSequence(GoIFasta,as.string=T)))
  #Vector of new names : 
 
  PreMatExport=PreMat[,-which(names(PreMat)=="Sequence")]
  
  #End of function
}


ComputeStatsFromFasta_ProtName <- function(FastaFilePath){
  library(seqinr)
  
  
  GoIFasta=read.fasta(file=FastaFilePath)
  
  PreMat=as.data.frame(matrix(0,nrow = length(GoIFasta),ncol=4),stringsAsFactors = F)
  names(PreMat)=c("OriginalName","ProtName","Length","Type")
  PreMat$OriginalName=unlist(getName(GoIFasta))
  #Stupidly complicated but ok
  PreMat$Length=getLength(GoIFasta)
  #Need to have upper case letter for some reasons. 
  PreMat$Sequence=toupper(unlist(getSequence(GoIFasta,as.string=T)))
  #Vector of new names : 
  #Vector of new names : 
  NamesToadd=names(unlist(AAstat(s2c(PreMat$Sequence[1]),plot = F)))
  #Creating new columns
  PreMat[,NamesToadd]=NA
  #Running this one liner, the transposition is necessary, I'm not sure why.
  PreMat[,NamesToadd]=t(apply(PreMat,1,function(x) as.numeric(as.character(unlist(AAstat(s2c(x[which(names(x)=="Sequence")]),plot = F))))))
  
  PreMatExport=PreMat[,-which(names(PreMat)=="Sequence")]
  
}
