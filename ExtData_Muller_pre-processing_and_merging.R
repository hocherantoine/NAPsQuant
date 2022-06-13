#####################################################
## Project: diaforarchaea chromatin
## Analysis of external data from Muller et al. 2020 DOI : 
## https://doi.org/10.1038/s41586-020-2402-x
## Aim : 
## 
## Pre-process Muller et al. data
## To account for match to multiple proteins
## and to keep only proteins present in genbank complete genome
## 
## Date: January 2021
## Author: Antoine Hocher
## 
## 
## NOTE : all raw and processed data for individual species processed in house (M. alvus, M. luminyensis, C. divulgatum,.. ) are available from PRIDE
## they are NOT hosted on github to avoid any confusion regarding versions.
####################################################

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.

Sys.setenv(PATH="/Users/ahocher/opt/miniconda3/bin:/Users/ahocher/opt/miniconda3/condabin:/opt/local/bin:/opt/local/sbin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin")

#Libraries : 
library(rtracklayer)
library(metablastr)
library(ape)
library(seqinr)



#Loading the raw supplementary table (For Github user please download it first): 

MassSpec=read.csv("41586_2020_2402_MOESM3_ESM.csv",stringsAsFactors = F)

#Typing the genomes that should be retrieved (only archaea in our case): 

Genomes=read.table("Annotated_genomes_info_to_retrieve_mass_spec.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

#Subselecting mass spectrometry results from archaea : 
MassSpecArchaea=MassSpec[which(MassSpec$Organism%in%Genomes$Specie_name_orig_dataset),]

table(MassSpecArchaea$Organism)
#Verification : 
#hist(table(MassSpecArchaea$Organism))
#Many organism have surprisingly low number of proteins identified, I do not know why
#It is however consistent with the original publication


# First step, for each entry of the mass spectrometry table, we will attempt to retrieve the corresponding Uniprot ID from our genbank gff database : 
GFFDir="/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/DB_Archaea95_update012020/genomes_assemblies_gff"
ProteomeDir="/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/DB_Archaea95_update012020/genome_assemblies_prot_fasta/ncbi-genomes-2020-05-14/"
setwd(GFFDir)
GFF_available=list.files()
ProteinToID=c()
for (Organism in names(table(MassSpecArchaea$Organism))){
  
Assembly=Genomes[which(Genomes$Specie_name_orig_dataset==Organism),]$Assembly
Gff_to_retrieve=GFF_available[grep(Assembly,GFF_available)]
Gffoi=readGFF(Gff_to_retrieve)

ProteinToID=c(unique(Gffoi$protein_id),ProteinToID)
}
ProteinToID=ProteinToID[is.na(ProteinToID)==F]
write.table(t(as.matrix(ProteinToID)),"/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/MASS_SPEC_Muller2020/AllGenbankID_from_Species_inMuller.txt",row.names=F,col.names=F,sep=" ",quote=F)


#After we retrieve the corresponding Uniprot / genbank link between IDs, we add those to the archaeal MassSpectable : 

IDmatch=read.table(file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/MASS_SPEC_Muller2020/Genbank_uniprot_matching/uniprot-yourlist_M20210104A94466D2655679D1FD8953E075198DA82B009C6.tab",header=T,sep="\t",quote="",stringsAsFactors = F)



#Have to do things one by one because some values in the initial table match to two proteins, 
#in this case I'll adopt the strategy of multiple even matching
#

MassSpecArchaea$ID=NA
#Produre when non unique id : 1st flag the multiple id'ed to delete later
#Separately produce a table of entries to add with the corrected intensity
MassSpecArchaea$ToDelete=NA
MassSpecArchaea$NonUnique=NA
MassSpecArchaeaToadd=as.data.frame(matrix(ncol = 7,nrow=0))
names(MassSpecArchaeaToadd)=names(MassSpecArchaea)

for(i in 1:dim(MassSpecArchaea)[1]){
  print(i)
  IDs=unlist(strsplit(MassSpecArchaea[i,]$Protein_IDs,";")[[1]])
  MatchingID=which(IDmatch$Entry%in%IDs)
  if(length(MatchingID)==1){
    MassSpecArchaea[i,]$ID=IDmatch[MatchingID,]$ID
  }
  if(length(MatchingID)>1){
    for(j in 1:length(MatchingID)){
    MultiEntry=MassSpecArchaea[i,]
    MultiEntry$Intensity=MultiEntry$Intensity/length(MatchingID)
    MultiEntry$ID=IDmatch[MatchingID[j],]$ID
    MultiEntry$NonUnique=length(MatchingID)
    MassSpecArchaeaToadd=rbind(MassSpecArchaeaToadd,MultiEntry)}
    #Now we can flag the entries that we de-duplicated :
    MassSpecArchaea[i,]$ToDelete=1
    }
}

#We add the lines we stored : 
MassSpecArchaea=rbind(MassSpecArchaea,MassSpecArchaeaToadd)

table(MassSpecArchaea$ToDelete)

MassSpecArchaea=MassSpecArchaea[-which(MassSpecArchaea$ToDelete==1),]


#There are 7.7% of the initial dataset that we don't include. 
table(is.na(MassSpecArchaea$ID))[2]/dim(MassSpecArchaea)[1]*100


#Export pre-processed table:
write.table(MassSpecArchaea,file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/MASS_SPEC_Muller2020/41586_2020_2402_MOESM3_ESM_archaeaOnly.txt",row.names=F,quote=F,sep="\t")



#Merging with data produced in-house:
MassSpecArchaea=read.table(file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/MASS_SPEC_Muller2020/41586_2020_2402_MOESM3_ESM_archaeaOnly.txt",header=T,sep="\t",quote="",stringsAsFactors = F)[,c(2:5)]

#We add luminyensis and Alvus : 
MassSpecAlvus=read.table("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/b023p050_M_alvus_reanalysis_proteinGroups_20210204.txt",header=T,sep='\t',quote="",stringsAsFactors = F)
names(MassSpecAlvus)[1]="UniprotID"
MassSpecAlvus$ID=MassSpecAlvus$UniprotID
MassSpecAlvus$Organism="Methanomethylophilus alvus"

#We merge all replicates (after having checked that they are all very well correlated, not shown here)
MassSpecAlvus$mean.LFQ.intensity_BIOREPs_WCE1.2=(MassSpecAlvus$LFQ.intensity.M.alvus_WCE_1_1+MassSpecAlvus$LFQ.intensity.M.alvus_WCE_1_2+MassSpecAlvus$LFQ.intensity.M.alvus_WCE_2_1+MassSpecAlvus$LFQ.intensity.M.alvus_WCE_2_2)/4

MassSpecAlvusm=MassSpecAlvus[,c("UniprotID","mean.LFQ.intensity_BIOREPs_WCE1.2","Organism","ID")]

#Making names consistent accross tables
names(MassSpecAlvusm)=names(MassSpecArchaea)

head(MassSpecAlvusm)


################
#Load Luminiensis data : 
MassSpecLumi=read.table("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/b023p050_M_lumi_reanalysis_proteinGroups_20210204.txt",header=T,sep='\t',quote="",stringsAsFactors = F)
MassSpecLumi=MassSpecLumi[,c(1,44:47)]

#We merge all replicates (after having checked that they are all very well correlated, not shown here)
MassSpecLumi$Average_WCE_LFQ=(MassSpecLumi$LFQ.intensity.M.lumi_WCE1_1+MassSpecLumi$LFQ.intensity.M.lumi_WCE1_2+MassSpecLumi$LFQ.intensity.M.lumi_WCE2_1+MassSpecLumi$LFQ.intensity.M.lumi_WCE2_2)/4
MassSpecLumi$Organism="Methanomassiliicoccus luminyensis"
MassSpecLumi$ID=MassSpecLumi$Protein.IDs
MassSpecLumim=MassSpecLumi[,c("Protein.IDs","Average_WCE_LFQ","Organism","ID")]
names(MassSpecLumim)=names(MassSpecArchaea)
MassSpecArchaeam=rbind(MassSpecArchaea,MassSpecAlvusm,MassSpecLumim)

write.table(MassSpecArchaeam,"/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/MASS_SPEC_Muller2020/41586_2020_2402_MOESM3_ESM_archaeaOnly_Alvus_Lumi.txt",row.names = F,sep='\t',quote=F)


