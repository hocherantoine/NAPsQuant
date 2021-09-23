#####################################################
## Project: diaforarchaea chromatin
## Analysis of operon mapper results
## https://doi.org/10.1038/s41586-020-2402-x
## Aim : 
## compute nb of tss / gene
## see if naps are never in an operon yes
## Associate proteomes to our local database
## 
## Date: feb 2021
## Author: Antoine Hocher
####################################################
source("/Users/ahocher/Dropbox/Scripts/PLOT_LIBRARY/AH_plot_library.R")
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.

library(rtracklayer)
DNAperSpecie=read.table("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/MASS_SPEC_PFAM/PerSpeciePFAM_abundancy.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
DNAperSpecie$PercentTSS=NA
setwd("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/OPERO_MAPPER/ARCHAEA/")
OperonFiles=list.files()
GFFDir="/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/DB_Archaea95_update012020/genomes_assemblies_gff/"
#initialise for merging : 
Opf=read.delim(OperonFiles[1],header=T,sep="\t",na.strings = "",stringsAsFactors = F)[,1:4]
Opf$OperonSize=NA
Opf$IdGeneReduced=NA
Opf$ProtName=NA
Opf$Specie=NA
Opf=Opf[1,]

for(FILE in OperonFiles){
  print(FILE)
NameFile=basename(FILE)
Assembly=paste("GCA_",strsplit(NameFile,"_")[[1]][2],sep="")

#there is a specific case for lumyniensis as I renamed proteins using their chromosme name
if(Assembly=="GCA_000308215.1"){
  GFFfile=list.files(GFFDir)[grep(Assembly,list.files(GFFDir))]
  GFF=readGFF(paste(GFFDir,GFFfile,sep=""))

  GFF$gene=gsub("cds-","",GFF$ID)
  GFF$protein_id=unlist(lapply(GFF$ID,function(x) strsplit(x,"_")[[1]][2]))
  GFF$protein_id=paste(GFF$seqid,"-",GFF$protein_id,sep="")
  
  
  }else{
  GFFfile=list.files(GFFDir)[grep(Assembly,list.files(GFFDir))]
GFF=readGFF(paste(GFFDir,GFFfile,sep=""))}



OutFile=paste("Processed/",NameFile,"_reformat_ah.txt",sep="")
if(file.exists(OutFile)==F){
  Op=read.delim(NameFile,header=T,sep="\t",na.strings = "",stringsAsFactors = F)[,1:4]
  Op$ProtName=NA
head(Op,n=20)

Liste=Op$Operon
for(i in 1:length(Liste)){
  if(is.na(Liste[i])==F){
    ValueToStore=Liste[i]}else{Liste[i]=ValueToStore}}
Liste
Op$Operon=Liste
Op=Op[which(is.na(Op$IdGene)==F),]

Tzar=table(Op$Operon)
Op$OperonSize=NA
for(j in unique(Op$Operon)){
  Op[which(Op$Operon==j),]$OperonSize=as.numeric(Tzar[j])
}


head(Op)

#Remove the cds- mark because it prevents good parsing
Op$IdGeneReduced=gsub("cds-","",Op$IdGene)

for(j in 1:dim(Op)[1]){
  print(j)
  Locus=Op[j,]$IdGeneReduced
  if(length(which(GFF$gene==Locus & GFF$type=="CDS"))>0){
    
  ProtName=GFF[which(GFF$gene==Locus & GFF$type=="CDS")[1],]$protein_id
  Op[j,]$ProtName=ProtName
  }else(Op[j,]$ProtName=Op[j,]$IdGeneReduced)
  }
  
Op$Specie=DNAperSpecie[which(DNAperSpecie$Assembly==Assembly),]$Specie
write.table(Op, file = OutFile,row.names = F,sep="\t",quote=F)
}else(Op=read.delim(OutFile,header=T,sep="\t",quote="",stringsAsFactors = F))
DNAperSpecie[which(DNAperSpecie$Assembly==Assembly),]$PercentTSS=100*length(table(Op$Operon))/dim(Op)[1]

Opf=rbind(Opf,Op)
}

head(Opf)
tail(Opf)
write.table(Opf, file = "Processed/AllOperons_MassSpec_Combined.txt",row.names = F,sep="\t",quote=F)


# Export
write.table(DNAperSpecie,file = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/MASS_SPEC_PFAM/PerSpeciePFAM_abundancy_and_TSSOperons.txt",row.names=F,quote=F,sep="\t")




#RELOAAAAD
DNAperSpecie=read.table(file = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/MASS_SPEC_PFAM/PerSpeciePFAM_abundancy_and_TSSOperons.txt",header=T,sep="\t",stringsAsFactors = F,quote="")


#Percent TSS vs NAPs
Correlationtest=cor.test(DNAperSpecie$PercentTSS,DNAperSpecie$PercentNAP,method="spearman")
options(scipen = 2)
PercentTSSNap=ggplot(data=DNAperSpecie,aes(y=PercentTSS,x=PercentNAP,color=Specie))+geom_point()+ggrepel::geom_text_repel(aes(label=Specie),size=2)+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("Spearman Cor ", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector)

ggsave(PercentTSSNap,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentTSS_OperonMapper_vsPFAMNAP.pdf",width=12,height = 6)


#Load #Plots with lines for CD-hit clusters : 

CandidateCount=read.table(file = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/NAP_PREDICTION/PerSpecie_CandidateNapCount_and_abundancy.txt",header=T,sep="\t",quote="",stringsAsFactors = F)


CandidateCountOperon=merge(CandidateCount,DNAperSpecie,by=c("Specie","PercentNAP"))
CandidateCountOperon$NAPs_w_candidates=CandidateCountOperon$PercentNAP+CandidateCountOperon$TotalProportionCandidates
Correlationtest=cor.test(CandidateCountOperon$PercentTSS,CandidateCountOperon$NAPs_w_candidates,method="spearman")
options(scipen = 2)
PercentTSSNAPs_w_candidates=ggplot(data=CandidateCountOperon,aes(y=PercentTSS,x=NAPs_w_candidates,color=Specie))+geom_point()+ggrepel::geom_text_repel(aes(label=Specie),size=2)+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("Spearman Cor ", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector)


ggsave(PercentTSSNAPs_w_candidates,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentTSS_OperonMapper_vsPFAMNAP_w_Candidates.pdf",width=12,height = 6)



