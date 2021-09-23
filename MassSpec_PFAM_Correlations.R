#####################################################
## Project: diaforarchaea chromatin
## Analysis of external data from Muller et al. 2020 DOI : 
## https://doi.org/10.1038/s41586-020-2402-x
## in relationship with pFAM domains.
## Aim : 
## This scripts input is from MassSpec_NAPs_detection
## 
## 
## Date: January 2021
## Author: Antoine Hocher
####################################################

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
source("/Users/ahocher/Dropbox/Scripts/PLOT_LIBRARY/AH_plot_library.R")
# 
# #Load annotated mass spec data : 
# MassSpecArchaea=read.table(file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/MASS_SPEC_Muller2020/41586_2020_2402_MOESM3_ESM_archaeaOnly_Alvus_Lumi_with_all_PfamAnnotation.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
# 
# MassSpecArchaea[which(MassSpecArchaea$Organism=="Methanomethylophilus alvus" & is.na(MassSpecArchaea$Intensity)==T),c(1:3)]
# #I have a weird bug -somtimes R changing my header but not for all tables-  with names being added X or not so we fix them using the make.names function:
# names(MassSpecArchaea)=make.names(names(MassSpecArchaea))
# 
# #Merging with Kodakarensis data : 
# MassSpecKodak=read.table(file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/MASS_SPEC_Muller2020/SanChenetal2020_Kodakarensis_annotated_PFAM_Domains.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
# names(MassSpecKodak)[1]="Protein_IDs"
# #I have a weird but with names being added X or not so we fix them :
# names(MassSpecKodak)=make.names(names(MassSpecKodak))
# Namelist.a=names(MassSpecArchaea);Namelist.b=names(MassSpecKodak)
# 
# Allnames=unique(c(Namelist.a,Namelist.b))
# 
# Toadd.a=Allnames[which(!(Allnames%in%Namelist.a))];Toadd.b=Allnames[which(!(Allnames%in%Namelist.b))]
# 
# MassSpecArchaea[,Toadd.a]=0;MassSpecKodak[,Toadd.b]=0
# 
# MassSpecKodak=MassSpecKodak[,names(MassSpecArchaea)]
# 
# MassSpecArchaea=rbind(MassSpecArchaea,MassSpecKodak)
# 
# 
# #Merging with Natrialba data : 
# MassSpecToAdd=read.table(file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/DIVERSEMASSSPEC/Natrialba_magalii_annotated_PFAM_Domains.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
# 
# 
# #I have a weird but with names being added X or not so we fix them :
# names(MassSpecToAdd)=make.names(names(MassSpecToAdd))
# Namelist.a=names(MassSpecArchaea);Namelist.b=names(MassSpecToAdd)
# 
# Allnames=unique(c(Namelist.a,Namelist.b))
# 
# Toadd.a=Allnames[which(!(Allnames%in%Namelist.a))];Toadd.b=Allnames[which(!(Allnames%in%Namelist.b))]
# 
# MassSpecArchaea[,Toadd.a]=0;MassSpecToAdd[,Toadd.b]=0
# 
# MassSpecToAdd=MassSpecToAdd[,names(MassSpecArchaea)]
# 
# MassSpecArchaea=rbind(MassSpecArchaea,MassSpecToAdd)
# 
# 
# #Merging with Haloferax volcanii data : 
# MassSpecToAdd=read.table(file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/DIVERSEMASSSPEC/Haloferax_volcanii_annotated_PFAM_Domains.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
# 
# 
# #I have a weird but with names being added X or not so we fix them :
# names(MassSpecToAdd)=make.names(names(MassSpecToAdd))
# Namelist.a=names(MassSpecArchaea);Namelist.b=names(MassSpecToAdd)
# 
# Allnames=unique(c(Namelist.a,Namelist.b))
# 
# Toadd.a=Allnames[which(!(Allnames%in%Namelist.a))];Toadd.b=Allnames[which(!(Allnames%in%Namelist.b))]
# 
# MassSpecArchaea[,Toadd.a]=0;MassSpecToAdd[,Toadd.b]=0
# 
# MassSpecToAdd=MassSpecToAdd[,names(MassSpecArchaea)]
# 
# MassSpecArchaea=rbind(MassSpecArchaea,MassSpecToAdd)
# 
# 
# 
# #Merging with Nmaritimus data : 
# MassSpecToAdd=read.table(file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/MASS_SPEC_Muller2020/Qinetal2018_Nmaritimus_annotated_PFAM_Domains.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
# 
# #I have a weird but with names being added X or not so we fix them :
# names(MassSpecToAdd)=make.names(names(MassSpecToAdd))
# Namelist.a=names(MassSpecArchaea);Namelist.b=names(MassSpecToAdd)
# 
# Allnames=unique(c(Namelist.a,Namelist.b))
# 
# Toadd.a=Allnames[which(!(Allnames%in%Namelist.a))];Toadd.b=Allnames[which(!(Allnames%in%Namelist.b))]
# 
# MassSpecArchaea[,Toadd.a]=0;MassSpecToAdd[,Toadd.b]=0
# 
# MassSpecToAdd=MassSpecToAdd[,names(MassSpecArchaea)]
# 
# MassSpecArchaea=rbind(MassSpecArchaea,MassSpecToAdd)
# 
# 
# 
# 
# #normalizing to be able to compare mass spec archaea to the initial mass spec data : : 
# MassSpecArchaea$NormInt=NA
# for(i in names(table(MassSpecArchaea$Organism))){
#   
#   MassSpecArchaea[which(MassSpecArchaea$Organism==i),]$NormInt=100*MassSpecArchaea[which(MassSpecArchaea$Organism==i),]$Intensity/sum(MassSpecArchaea[which(MassSpecArchaea$Organism==i),]$Intensity)
# }
# 
# 
# write.table(MassSpecArchaea,"/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/MASS_SPEC_Muller2020/41586_2020_2402_MOESM3_ESM_archaeaOnly_Alvus_Lumi_Kodak_natrialba_Haloferax_Nmaritimus_With_PFam_annotations.txt",row.names = F,quote=F,sep="\t")
# 
# 



#Load annotated mass spec data : 
MassSpecArchaea=read.table(file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/MASS_SPEC_Muller2020/41586_2020_2402_MOESM3_ESM_archaeaOnly_Alvus_Lumi_Kodak_natrialba_Haloferax_Nmaritimus_With_PFam_annotations.txt",header=T,sep="\t",quote="",stringsAsFactors = F)






#Loading archaeal genomes ids :  
Genomes=read.table("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/MASS_SPEC_Muller2020/Annotated_genomes_info_to_retrieve_mass_spec.txt",header=T,sep="\t",quote="",stringsAsFactors = F)



#Adding all the external mass spec species  :
Genomes[which(Genomes$Specie=="Thermococcus^kodakarensis^KOD1"),]$To_retrieve_MassSpec=1
Genomes[which(Genomes$Specie=="Thermococcus^kodakarensis^KOD1"),]$Specie_name_orig_dataset="Thermococcus^kodakarensis^KOD1"

Genomes[which(Genomes$Specie=="Haloferax^volcanii^DS2"),]$To_retrieve_MassSpec=1
Genomes[which(Genomes$Specie=="Haloferax^volcanii^DS2"),]$Specie_name_orig_dataset="Haloferax volcanii"

#
Genomes[which(Genomes$Specie=="Natrialba^magadii^ATCC^43099"),]$To_retrieve_MassSpec=1
Genomes[which(Genomes$Specie=="Natrialba^magadii^ATCC^43099"),]$Specie_name_orig_dataset="Natrialba magadii"

Genomes[which(Genomes$Specie=="Nitrosopumilus^maritimus^SCM1"),]$To_retrieve_MassSpec=1
Genomes[which(Genomes$Specie=="Nitrosopumilus^maritimus^SCM1"),]$Specie_name_orig_dataset="Nitrosopumilus^maritimus^SCM1"

Genomes=Genomes[,-which(names(Genomes)=="Specie")]
Genomesppties=read.table("/Users/ahocher/Dropbox/Laboratory/ArchealChromatinProteom/Nucleotide_periodicity/Analysis/Archaea_EBMC2_genomes_properties_results.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
Genomesppties=Genomesppties[,-which(names(Genomesppties)=="Specie")]

OperonPred=read.table("/Users/ahocher/Dropbox/Laboratory/ArchealChromatinProteom/Nucleotide_periodicity/Analysis/Archaea_EBMC2_Operons_properties_results.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
OperonPred=OperonPred[,which(names(OperonPred)%in%c("Assembly","Tss.Proportion"))]

Genomesppties=merge(Genomesppties,OperonPred,by="Assembly")

#Exporting a list of assemblies : 
Assemblies=Genomes[which(Genomes$To_retrieve_MassSpec==1),]$Assembly



#Remove organism with less than 500 proteins : 
table(MassSpecArchaea$Organism)

#Annotate DNA binding Pfam and NAPs

#Reloading the manually annotated table :
PfamDNAbinding=read.table(file="/Users/ahocher/Dropbox/Laboratory/ArchealChromatinProteom/HMMs/DBPome/PFAM_DNA_binding_list_present_MassSpecDataSet_annot_AH.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

TFHmm=PfamDNAbinding[which(PfamDNAbinding$TF==1),]$PfamID
PfamDNAbinding=PfamDNAbinding$PfamID

NAPsHMM=c("Cren7","X7kD_DNA_binding","CBFD_NFYB_HMF","Alba","Bac_DNA_binding","MC1","Histone")
AllPfamID=make.names(read.table("/Users/ahocher/Dropbox/Laboratory/ArchealChromatinProteom/HMMs/PfamDb/All_pfamID.txt")$V1)
AllPfamID=AllPfamID[-which(AllPfamID=="ID")]

#There is a PFAM to GO Annotation -GO:0003700-, we load it.
Pfam2Go=read.table(file="/Users/ahocher/Dropbox/Laboratory/ArchealChromatinProteom/HMMs/PfamDb/Pfam2go.txt",skip=6,sep=";",stringsAsFactors = F)
head(Pfam2Go)
Pfam2Go$ID=unlist(lapply(Pfam2Go$V1,function(x) strsplit(x,split = " ")[[1]][2]))
Pfam2Go$Description=unlist(lapply(Pfam2Go$V1,function(x) strsplit(x,split = ">")[[1]][2]))
names(Pfam2Go)[2]="GO"
Pfam2Go$GO=gsub(" ","",Pfam2Go$GO)
Pfam2Go$GO=gsub(":",".",Pfam2Go$GO)

TFPFAMGO=Pfam2Go[which(Pfam2Go$GO=="GO.0003700"),]$ID

#The final TF HMM models are  :

TF_HMMs=unique(TFHmm,TFPFAMGO)


#annotating : DNA binding, and NAP : 
#Annotating all proteins containing a DNA binding domain : 
MassSpecArchaea$PFAM_DNA_binding=rowSums(MassSpecArchaea[,which(names(MassSpecArchaea)%in%PfamDNAbinding)])
#Binarize
MassSpecArchaea[which(MassSpecArchaea$PFAM_DNA_binding>0),]$PFAM_DNA_binding=1

MassSpecArchaea$NAP=rowSums(MassSpecArchaea[,which(names(MassSpecArchaea)%in%NAPsHMM)])
#Binarize
MassSpecArchaea[which(MassSpecArchaea$NAP>0),]$NAP=1

#TF : 
#Annotating Transcription factors : 
MassSpecArchaea$TFPfam=0
MassSpecArchaea[which(rowSums(MassSpecArchaea[,which(names(MassSpecArchaea)%in%TF_HMMs)])>0),]$TFPfam=1




#Adding CC1 at last : 
CC1=read.table(file = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/HMMSEARCH/CC1/Archaea_Distribution_CC1_Proteins_with_Length.txt",header=T,sep="\t",stringsAsFactors = F)
CC1=CC1[which(CC1$Length<200),]
MassSpecArchaea$CC1=0
MassSpecArchaea[which(MassSpecArchaea$ID%in%CC1$ProtName),]$CC1=1




write.table(MassSpecArchaea,file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/MASS_SPEC_Muller2020/41586_2020_2402_MOESM3_ESM_archaeaOnly_Alvus_Lumi_Kodak_natrialba_Haloferax_Nmaritimus_With_PFam_annotations.txt",row.names = F,sep="\t",quote=F)





#Building a table to host results : 
DNAperSpecie=as.data.frame(table(MassSpecArchaea$Organism))
names(DNAperSpecie)=c("Specie","NbDetectedProt")
DNAperSpecie=merge(DNAperSpecie,Genomes,by.x="Specie",by.y="Specie_name_orig_dataset")
DNAperSpecie=merge(DNAperSpecie,Genomesppties,by=c("Assembly","Taxid"))


#Creating new variables, 1 steps is for whole categories
DNAperSpecie$NbDNA=NA;DNAperSpecie$NbNAP=NA
DNAperSpecie$TotalProt=NA;DNAperSpecie$TotalProtDNA=NA;
DNAperSpecie$TotalProtGenbank=NA
DNAperSpecie$TotalProtNAP=NA;DNAperSpecie$TotalProtDNAnotNAP=NA
DNAperSpecie$TotalProtNoPfam=NA
DNAperSpecie$TotalTF=NA
for(i in 1:dim(DNAperSpecie)[1]){
  DNAperSpecie[i,]$NbDNA=length(which(MassSpecArchaea$Organism==DNAperSpecie[i,]$Specie & MassSpecArchaea$PFAM_DNA_binding==1))
  DNAperSpecie[i,]$NbNAP=length(which(MassSpecArchaea$Organism==DNAperSpecie[i,]$Specie & MassSpecArchaea$NAP==1))
  DNAperSpecie[i,]$TotalProt=sum(MassSpecArchaea[which(MassSpecArchaea$Organism==DNAperSpecie[i,]$Specie),]$Intensity)
  DNAperSpecie[i,]$TotalProtGenbank=sum(MassSpecArchaea[which(MassSpecArchaea$Organism==DNAperSpecie[i,]$Specie & is.na(MassSpecArchaea$ID)==F),]$Intensity)
  
  DNAperSpecie[i,]$TotalProtDNA=sum(MassSpecArchaea[which(MassSpecArchaea$Organism==DNAperSpecie[i,]$Specie & MassSpecArchaea$PFAM_DNA_binding==1),]$Intensity)
  DNAperSpecie[i,]$TotalProtNAP=sum(MassSpecArchaea[which(MassSpecArchaea$Organism==DNAperSpecie[i,]$Specie & MassSpecArchaea$NAP==1),]$Intensity)
  DNAperSpecie[i,]$TotalProtDNAnotNAP=sum(MassSpecArchaea[which(MassSpecArchaea$Organism==DNAperSpecie[i,]$Specie & MassSpecArchaea$PFAM_DNA_binding==1 & MassSpecArchaea$NAP==0),]$Intensity)
  DNAperSpecie[i,]$TotalProtNoPfam=sum(MassSpecArchaea[which(MassSpecArchaea$Organism==DNAperSpecie[i,]$Specie & MassSpecArchaea$NoPfam==1),]$Intensity)
  DNAperSpecie[i,]$TotalTF=sum(MassSpecArchaea[which(MassSpecArchaea$Organism==DNAperSpecie[i,]$Specie & MassSpecArchaea$TFPfam==1),]$Intensity)
  
}

#Converting into proportions : 
DNAperSpecie$PercentDNA=100*DNAperSpecie$TotalProtDNA/DNAperSpecie$TotalProt
DNAperSpecie$PercentNAP=100*DNAperSpecie$TotalProtNAP/DNAperSpecie$TotalProt
DNAperSpecie$PercentNAPDNA=100*DNAperSpecie$TotalProtNAP/DNAperSpecie$TotalProtDNA
DNAperSpecie$PercentDNAnotNAP=100*DNAperSpecie$TotalProtDNAnotNAP/DNAperSpecie$TotalProt
DNAperSpecie$PercentnoPfam=100*DNAperSpecie$TotalProtNoPfam/DNAperSpecie$TotalProt
DNAperSpecie$PercentTF=100*DNAperSpecie$TotalTF/DNAperSpecie$TotalProt



#Compute percent of proteome for each PFAM annotations : 
PFAMDomains=names(MassSpecArchaea)[7:c(dim(MassSpecArchaea)[2]-1)]

for( PF in PFAMDomains){
  print(PF)
  DNAperSpecie[,PF]=0
  for(i in 1:dim(DNAperSpecie)[1]){
    DNAperSpecie[i,PF]=100*sum(MassSpecArchaea[which(MassSpecArchaea$Organism==DNAperSpecie[i,]$Specie & MassSpecArchaea[,PF]==1),]$Intensity)/DNAperSpecie[i,]$TotalProt}}

# Cleanup
VarToExclude=names(which(colSums(DNAperSpecie[,PFAMDomains])==0))

if(length(which(names(DNAperSpecie)%in%VarToExclude))>0){
DNAperSpecie=DNAperSpecie[,-which(names(DNAperSpecie)%in%VarToExclude)]}

# Export
write.table(DNAperSpecie,file = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/MASS_SPEC_PFAM/PerSpeciePFAM_abundancy.txt",row.names=F,quote=F,sep="\t")




#Re-load
DNAperSpecie=read.table("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/MASS_SPEC_PFAM/PerSpeciePFAM_abundancy.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

#Add the NAP candidates ( NAP detection scheme has to be run 1st): 
CandidatePerSpecie=read.table("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/NAP_PREDICTION/PerSpecie_CandidateNapCount_and_abundancy.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
CandidatePerSpecie$PercentNAPwCandidates=CandidatePerSpecie$PercentNAP+CandidatePerSpecie$TotalProportionCandidates
CandidatePerSpecie=CandidatePerSpecie[-which(names(CandidatePerSpecie)=="PercentNAP")]

DNAperSpecie=merge(DNAperSpecie,CandidatePerSpecie,by="Specie")

DNAperSpecie$Specie[4]="Halobacterium salinarium"
DNAperSpecie$Specie[16]="Sulfolobus acidocaldarius"

DNAperSpecie$Specie[17]=gsub("\\^"," ",DNAperSpecie$Specie[17])
DNAperSpecie$Specie[14]=gsub("\\^"," ",DNAperSpecie$Specie[14])

DNAperSpecie$SpecieShort=unlist(lapply(DNAperSpecie$Specie,function(x) paste(strsplit(strsplit(x,'')[[1]][1],' ')[[1]][1],strsplit(x," ")[[1]][2],sep=".")))










#PCA on species to see how well PFAM proportions : describes species: 
##First on all PFAM domains 
#VERIFITY MAYBE BEST TO GO THROUGH CORRELATIONS


PCAInput=DNAperSpecie[,which(names(DNAperSpecie)%in%PFAMDomains)]
row.names(PCAInput)=DNAperSpecie$Specie
PCA=prcomp(PCAInput)
library(ggfortify)
pdf("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/PCA_on_PFAM_and_species.pdf",width=12,height=8)
A=autoplot(PCA, data = DNAperSpecie, colour = 'Specie',label=T,label.size=3)+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")
B=autoplot(PCA, data = DNAperSpecie,x=2,y=3, colour = 'Specie',label=T,label.size=3)+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")
ggarrange(A,B)
dev.off()





#Quick fix for the compatibility
MassSpecArchaea$Specie=MassSpecArchaea$Organism

#We use a tree to order the plot : 
ArchTree=read.tree("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/TREES/GTDB/Arch_120_pruned_to_MassSpec_dataSet.tree")
ArchTree$tip.label=gsub("_"," ",ArchTree$tip.label)
ArchTree$tip.label=gsub("\\^"," ",ArchTree$tip.label)
ArchTree$tip.label[17]="Sulfolobus acidocaldarius"
#Simpler, correlation matrix

CorRes=cor(t(as.matrix(PCAInput)),method="spearman")
NewOrder=colnames(CorRes)[match(ArchTree$tip.label,colnames(CorRes))]
NewOrder=c("Halobacterium salinarium",NewOrder[c(1:10)],"Picrophilus torridus",NewOrder[c(11:19)])

if(length(which(is.na(NewOrder)==T))>0){
  NewOrder=NewOrder[-which(is.na(NewOrder)==T)]}
#Re-ordering the matrix : 
CorRes=CorRes[NewOrder,NewOrder]

pdf("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/Correlation_on_PFAMDomains_and_species.pdf",width=12,height=8)
corrplot(CorRes, method="color", cl.lim = c(0,1),diag = F)
dev.off()





#DO THE SAME WITH DNA BINDING PFAM : 
#PCA on species to see how DNA binding PFAM describes species: 

#Restrict to Hmm corresponding to PFAM DNA binding , or to NAPs :
PfamDNAbinding=read.table(file="/Users/ahocher/Dropbox/Laboratory/ArchealChromatinProteom/HMMs/DBPome/PFAM_DNA_binding_list.txt",stringsAsFactors = F)$V1

PfamDNAbinding=c(PfamDNAbinding,"TFIIE_alpha")


PCAInput=DNAperSpecie[,which(names(DNAperSpecie)%in%PfamDNAbinding)]
row.names(PCAInput)=DNAperSpecie$Specie
PCA=prcomp(dist(cor(t(as.matrix(PCAInput)))))
pdf("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/PCA_on_PFAM_DNAbinding_and_species.pdf",width=12,height=8)
A=autoplot(PCA, data = DNAperSpecie, colour = 'Specie',label=T,label.size=3)+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")
B=autoplot(PCA, data = DNAperSpecie,x=2,y=3, colour = 'Specie',label=T,label.size=3)+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")
ggarrange(A,B)
dev.off()


#V IMPTT DOUBEL CHECK IM DOING PCA AS I SHOULD : on correlations

SpecieCorPfam=cor(t(as.matrix(PCAInput)),method="spearman")
pdf("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/PCA_on_Correaltions_PFAM_DNAbinding_and_species.pdf",width=12,height=8)
autoplot(prcomp(dist(SpecieCorPfam)), label = TRUE, label.size = 3)
dev.off()



library(corrplot)
pdf("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/Correlation_on_PFAM_DNABINDING_Domains_and_species.pdf",width=12,height=8)
corrplot(CorRes, method="color", cl.lim = c(0,1),diag = F)
dev.off()



#Now computing the correlations with being a Diaforarchaea : 
DNAperSpecie$Diaf=0
DNAperSpecie[which(DNAperSpecie$Specie%in%c("Methanomassiliicoccus luminyensis","Methanomethylophilus alvus","Thermoplasma acidophilum","Thermoplasma volcanium","Picrophilus torridus")),]$Diaf=1


#CORRELATIONS between all pfam domains : 
#Aim is to see what goes down with histones of TFIIEa for example
AA=as.data.frame(cor(DNAperSpecie[,which(names(DNAperSpecie)%in%c(PFAMDomains,"DivergentIntergeneProportion","ConvergentIntergeneProportion","AlignedIntergeneProportion","CompactionBp","GeneNb","NbDetectedProt","Diaf","Tss.Proportion","PercentNAP","PercentNAPwCandidates","GC.content.avg"))],use="complete.obs",method="spearman"))

write.table(AA,file = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/MASS_SPEC_PFAM/MASSSpec_PFAM_abundancy_correlation.txt",row.names=F,quote=F,sep="\t")




setwd("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/")


Selecta=order(AA$TFIIE_alpha)[c(1:25,c((dim(AA)[1]-25):dim(AA)[1]))]
postscript("Correlation_PFAM_Abundancy_Top50_TFIIEa.eps",width=12,height=12)
heatmap3::heatmap3(as.matrix(AA[Selecta,Selecta]),scale = "none",Rowv=NA,Colv = NA)
dev.off()


Selecta=order(AA$TFIIE_alpha)[c(1:50,c((dim(AA)[1]-50):dim(AA)[1]))]
postscript("Correlation_PFAM_Abundancy_Top100_TFIIEa.eps",width=12,height=12)
heatmap3::heatmap3(as.matrix(AA[Selecta,Selecta]),scale = "none",Rowv=NA,Colv = NA,cexRow = 0.35,cexCol  = 0.35)
dev.off()

Selecta=order(AA$CBFD_NFYB_HMF)[c(1:25,c((dim(AA)[1]-25):dim(AA)[1]))]
postscript("Correlation_PFAM_Abundancy_Top50_CBFD_NFYB_HMF.eps",width=12,height=12)
heatmap3::heatmap3(as.matrix(AA[Selecta,Selecta]),scale = "none",Rowv=NA,Colv = NA,cexRow = 0.35,cexCol  = 0.35)
dev.off()

Selecta=order(AA$Compaction)[c(1:25,c((dim(AA)[1]-25):dim(AA)[1]))]
postscript("Correlation_PFAM_Abundancy_Top50_Compaction.eps",width=12,height=12)
heatmap3::heatmap3(as.matrix(AA[Selecta,Selecta]),scale = "none",Rowv=NA,Colv = NA,cexRow = 0.35,cexCol  = 0.35)
dev.off()

Selecta=order(AA$CompactionBp)[c(1:25,c((dim(AA)[1]-25):dim(AA)[1]))]
postscript("Correlation_PFAM_Abundancy_Top50_CompactionBp.eps",width=12,height=12)
heatmap3::heatmap3(as.matrix(AA[Selecta,Selecta]),scale = "none",Rowv=NA,Colv = NA,cexRow = 0.35,cexCol  = 0.35)
dev.off()


Selecta=order(AA$Tss.Proportion)[c(1:50,c((dim(AA)[1]-50):dim(AA)[1]))]
postscript("Correlation_PFAM_Abundancy_Top50_Tss.Proportion.eps",width=12,height=12)
heatmap3::heatmap3(as.matrix(AA[Selecta,Selecta]),scale = "none",Rowv=NA,Colv = NA,cexRow = 0.35,cexCol  = 0.35)
dev.off()




library(viridis)
Selecta=order(AA$PercentNAP)[c(1:25,c((dim(AA)[1]-25):dim(AA)[1]))]
postscript("Correlation_PFAM_Abundancy_Top50_PercentNAP.eps",width=12,height=12)
heatmap3::heatmap3(as.matrix(AA[Selecta,Selecta]),scale = "none",cexRow = 0.35,cexCol  = 0.35,col=viridis(100))
dev.off()


Selecta=order(AA$PercentNAPwCandidates)[c(1:25,c((dim(AA)[1]-25):dim(AA)[1]))]
pdf("Correlation_PFAM_Abundancy_Top50_PercentNAPwcandidates.pdf",width=12,height=12)
heatmap3::heatmap3(as.matrix(AA[Selecta,Selecta]),scale = "none",cexRow = 0.35,cexCol  = 0.35,col=viridis(100))
dev.off()




Selecta=order(AA$Tss.Proportion)[c(1:25,c((dim(AA)[1]-25):dim(AA)[1]))]
postscript("Correlation_PFAM_Abundancy_Top50_Tss.Proportion.eps",width=12,height=12)
heatmap3::heatmap3(as.matrix(AA[Selecta,Selecta]),scale = "none",cexRow = 0.35,cexCol  = 0.35)
dev.off()



###################
#Specific plot on histones and other picked category: 


library(ggplot2);library(ggpubr)

NewOrder=rev(NewOrder)
DNAperSpecie$Specie=factor(DNAperSpecie$Specie,levels = NewOrder)

HistonePlot=ggplot()+geom_bar(data=DNAperSpecie[,c("CBFD_NFYB_HMF","Specie")],aes(x=Specie,y=CBFD_NFYB_HMF),stat = "identity", width=0.7,fill="darkblue")+geom_point(data=MassSpecArchaea[which(MassSpecArchaea$CBFD_NFYB_HMF==1),c("Specie","NormInt")],aes(x=Specie,y=NormInt),col="grey80",size=0.8)+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of proteome")+ggtitle("Histone abundancy accross species")

ggsave(HistonePlot,filename ="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CBFD_NFYB_HMF_Hits_abundancy_perSpecie.pdf",width = 6,height=4 )

#OtherFactors : 
for(FactoPlot in c("Alba","TFIIE_alpha","RNA_pol_Rpb1_3","PadR","Bac_DNA_binding","DUF1931","Regulator_TrmB","MC1","X7kD_DNA_binding")){
TheFactoPlot=ggplot()+geom_bar(data=DNAperSpecie[,c(FactoPlot,"Specie")],aes_string(x="Specie",y=FactoPlot),stat = "identity", width=0.7,fill="darkblue")+geom_point(data=MassSpecArchaea[which(MassSpecArchaea[,FactoPlot]==1),c("Specie","NormInt")],aes(x=Specie,y=NormInt),col="grey80",size=0.8)+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of proteome")+ggtitle(paste(FactoPlot," abundancy accross species",sep=""))
ggsave(TheFactoPlot,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/",FactoPlot,"_Hits_abundancy_perSpecie.pdf",sep=""),width = 6,height=4 )
}


#Plotting the % of protein that are NAPs
PercentNAPPLOT=ggplot()+geom_bar(data=DNAperSpecie[,c("PercentNAP","Specie")],aes_string(x="Specie",y="PercentNAP"),stat = "identity", width=0.7,fill="darkblue")+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of proteome")+ggtitle(paste("NAPs, accross species",sep=""))

ggsave(PercentNAPPLOT,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/Proportion_of_proteome_NAP_perSpecie.pdf",sep=""),width = 6,height=4)

#Plotting the % of protein that are NAPs
CompactionPLOT=ggplot()+geom_bar(data=DNAperSpecie[,c("CompactionBp","Specie")],aes_string(x="Specie",y="100*CompactionBp"),stat = "identity", width=0.7,fill="darkblue")+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of Genome")+ggtitle(paste("NAPs, accross species",sep=""))

ggsave(CompactionPLOT,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/Proportion_of_CodingIntergene_perSpecie.pdf",sep=""),width = 6,height=4)


#Plotting the % of genes with TSS
TSSPLOT=ggplot()+geom_bar(data=DNAperSpecie[,c("Tss.Proportion","Specie")],aes_string(x="Specie",y="Tss.Proportion"),stat = "identity", width=0.7,fill="darkblue")+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of Genome")+ggtitle(paste("NAPs, accross species",sep=""))

ggsave(TSSPLOT,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/Proportion_TSS_perSpecie.pdf",sep=""),width = 6,height=4)



#Plotting the % of protein that do not have PFAM domain hits.
NoPFAMPLOT=ggplot()+geom_bar(data=DNAperSpecie[,c("PercentnoPfam","Specie")],aes_string(x="Specie",y="PercentnoPfam"),stat = "identity", width=0.7,fill="darkblue")+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of proteome")+ggtitle(paste("no PFAM domains, accross species",sep=""))

ggsave(NoPFAMPLOT,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/Proportion_of_proteome_noPFAM_perSpecie.pdf",sep=""),width = 6,height=4)







#Plotting the % of proteome covered 
DNAperSpecie$Coverage=100*DNAperSpecie$NbDetectedProt/DNAperSpecie$GeneNb
DNAperSpecie=DNAperSpecie[order(-DNAperSpecie$Coverage),]
DNAperSpecie$Specie=factor(DNAperSpecie$Specie,levels = DNAperSpecie$Specie)
CoveragePlot=ggplot()+geom_bar(data=DNAperSpecie,aes(x=Specie,y=Coverage),stat = "identity", width=0.7,fill="darkblue")+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6),aspect.ratio = 1/3)+ylab("% of proteome measured")+ggtitle(paste("Proportion of proteome coverage",sep=""))+xlab("")

ggsave(CoveragePlot,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/Proportion_of_proteome_with_Coverage_perSpecie.pdf",sep=""),width = 6,height=4)

write.table(DNAperSpecie[,which(names(DNAperSpecie)%in%c("SpecieShort","Coverage"))],file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/Coverage_all_data_INhouse_and_muller.txt",row.names = F,sep="\t")


#Restart
#Focus on correlations with NAP abundancy : 

CorResT=data.frame(ID=AllPfamID[which(AllPfamID%in%names(DNAperSpecie))],Cor=NA,pval=NA,NbZero=NA)

for(i in 1: dim(CorResT)[1]){
  Test=cor.test(DNAperSpecie[,which(names(DNAperSpecie)==CorResT[i,]$ID)],DNAperSpecie$PercentNAPwCandidates,method="spearman")
  CorResT[i,]$Cor=Test$estimate
  CorResT[i,]$pval=Test$p.value
  CorResT[i,]$NbZero=length(which(DNAperSpecie[,which(names(DNAperSpecie)==CorResT[i,]$ID)]==0))
}

#To prevent spurious correlations and reduce the number of hits to an intelligible subset : 
CorResT=CorResT[which(CorResT$NbZero<14),]
CorResT=CorResT[which(CorResT$pval<0.01),]

CorResT=CorResT[order(CorResT$Cor),]

MatPlot=cbind(as.matrix(CorResT$Cor),as.matrix(CorResT$Cor))
row.names(MatPlot)=CorResT$ID
pdf("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/PFAM_vs_GENOME_PROP/PercentNAPwCandidates/Correlation_SignificantHits_SumNAPwcandidates_less14zeros.pdf",width=6,height=12)
heatmap3::heatmap3(MatPlot,scale = "none",cexRow = 0.35,cexCol  = 0.35,col=viridis(100),Rowv = NA,Colv = NA)
dev.off()


CorResT$ID=as.character(CorResT$ID)
#plotting those : 
for(IDPLOT in CorResT$ID){
Correlationtest=CorResT[which(CorResT$ID==IDPLOT),]
plotExport=addSmallLegend(ggplot(data=DNAperSpecie,aes_string(y=IDPLOT,x="PercentNAPwCandidates"))+geom_point(shape=1,size=1.8,stroke=0.8)+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25)+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$Cor,3), " p-value :", round(Correlationtest$pval,digits = 6),sep=""))+xlab("NAP abundancy (%)")+ylab(paste(IDPLOT," abundancy (%)"))+scale_color_manual(values=col_vector[10:30]))
ggsave(plotExport,filename = paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/PFAM_vs_GENOME_PROP/PercentNAPwCandidates/Formated_PercentNAPwCandidates_vs_",IDPLOT,".pdf",sep=""),width=4,height=4)
}


#Focus on eIF.5a : 
DNAperSpecie$RiboL=rowSums(DNAperSpecie[,grep("Ribosomal_L",names(DNAperSpecie))])
DNAperSpecie$RiboS=rowSums(DNAperSpecie[,grep("Ribosomal_S",names(DNAperSpecie))])



Correlationtest=cor.test(DNAperSpecie$eIF.5a,DNAperSpecie$PercentNAP,method="spearman")

plotExport=addSmallLegend(ggplot(data=DNAperSpecie,aes(y=eIF.5a,x=PercentNAP))+geom_point(shape=1,size=1.8,stroke=0.8)+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25)+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+xlab("NAP abundancy (%)")+ylab(paste("eIF5a abundancy (%)"))+scale_color_manual(values=col_vector[10:30]))

ggsave(plotExport,filename = paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/PFAM_vs_GENOME_PROP/PercentNAPwCandidates/Formated_PercentNAP_vs_","eIF5a",".pdf",sep=""),width=4,height=4)


DNAperSpecie$eIF.5aNorm=DNAperSpecie$eIF.5a/(DNAperSpecie$RiboL+DNAperSpecie$RiboS)
Correlationtest=cor.test(DNAperSpecie$eIF.5aNorm,DNAperSpecie$PercentNAPwCandidates,method="spearman")

plotExport=addSmallLegend(ggplot(data=DNAperSpecie,aes(y=eIF.5a,x=PercentNAPwCandidates))+geom_point(shape=1,size=1.8,stroke=0.8)+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25)+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+xlab("NAP abundancy (%)")+ylab(paste("eIF5a abundancy (%)"))+scale_color_manual(values=col_vector[10:30]))
ggsave(plotExport,filename = paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/PFAM_vs_GENOME_PROP/PercentNAPwCandidates/Formated_PercentNAPwCandidates_vs_","eIF5aNormalized_Ribosomes",".pdf",sep=""),width=4,height=4)




addSmallLegend(ggplot(data=DNAperSpecie,aes(y=RiboL,x=PercentNAPwCandidates))+geom_point(shape=1,size=1.8,stroke=0.8)+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25)+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+xlab("NAP abundancy (%)")+ylab(paste("eIF5a abundancy (%)"))+scale_color_manual(values=col_vector[10:30]))


 #Example : 
#Plot PFAM abundancy correlation between M.alvus and M.luminiensis :
#Transformation of the table for the plot is not elegant but has been double checked.
library(reshape2)
AlvLumComp=as.data.frame(t(as.matrix(DNAperSpecie[which(DNAperSpecie$Specie%in%c("Methanomassiliicoccus luminyensis","Methanomethylophilus alvus")),which(names(DNAperSpecie)%in%c("Specie",AllPfamID))])))
names(AlvLumComp)=unlist(AlvLumComp[1,])
AlvLumComp=AlvLumComp[-1,]
AlvLumComp$`Methanomethylophilus alvus`=as.numeric(as.character(AlvLumComp$`Methanomethylophilus alvus`))
AlvLumComp$`Methanomassiliicoccus luminyensis`=as.numeric(as.character(AlvLumComp$`Methanomassiliicoccus luminyensis`))

Correlationtest=cor.test(AlvLumComp$`Methanomethylophilus alvus`,AlvLumComp$`Methanomassiliicoccus luminyensis`,method="spearman")

AlvuLumComp=ggplot(data = AlvLumComp)+geom_point(aes(x=`Methanomethylophilus alvus`,y=`Methanomassiliicoccus luminyensis`),color="grey70",size=0.8,stroke=0.8)+scale_x_log10()+scale_y_log10()+geom_abline(intercept = 0)+xlab("M.alvus, % of proteome\nper PFAM domain")+ylab("M.luminyensis, % of proteome\nper PFAM domain")+theme_pubr()+theme(aspect.ratio = 1)+ggtitle(paste("Rho = ", round(Correlationtest$estimate,2), " p-value :", round(Correlationtest$p.value,3),sep=""))
ggsave(AlvuLumComp,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/Comparison_Lumi_Alvus_PFAMdomains.pdf",sep=""),width = 4,height=4)





library(ggrepel)

setwd("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/")


AA=as.data.frame(cor(DNAperSpecie[,which(names(DNAperSpecie)%in%c(PFAMDomains,"DivergentIntergeneProportion","ConvergentIntergeneProportion","AlignedIntergeneProportion","CompactionBp","GeneNb","NbDetectedProt","Diaf","Tss.Proportion","PercentNAP","PercentNAPwCandidates","GC.content.avg","Diaf"))],use="complete.obs",method="spearman"))

Selecta=order(AA$CBFD_NFYB_HMF)[c(1:30,c((dim(AA)[1]-30):dim(AA)[1]))]

PlotList=row.names(AA[Selecta,c("Alba","CBFD_NFYB_HMF")])

for( i in PlotList){

  PLOTOTOSAVE=ggplot(DNAperSpecie,aes_string(x="CBFD_NFYB_HMF",y=i,color="Specie",label="Specie"))+geom_point()+geom_text_repel(size=2)+theme_pubr()+theme(aspect.ratio = 1,legend.position="")
  
  ggsave(PLOTOTOSAVE,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/PFAM_vs_GENOME_PROP/CBFD_NFYB_HMF/CBFD_NFYB_HMV_vs_",i,".pdf",sep=""),width = 6,height=4)
  }


#On Alba : 


Selecta=order(AA$Alba)[c(1:30,c((dim(AA)[1]-30):dim(AA)[1]))]

PlotList=row.names(AA[Selecta,c("Alba","CBFD_NFYB_HMF")])

for( i in PlotList){
  
  PLOTOTOSAVE=ggplot(DNAperSpecie,aes_string(x="Alba",y=i,color="Specie",label="Specie"))+geom_point()+geom_text_repel(size=2)+theme_pubr()+theme(aspect.ratio = 1,legend.position="")
  
  ggsave(PLOTOTOSAVE,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/PFAM_vs_GENOME_PROP/ALBA/Alba_vs_",i,".pdf",sep=""),width = 6,height=4)
}

#On TFIIE_alpha
Selecta=order(AA$TFIIE_alpha)[c(1:30,c((dim(AA)[1]-30):dim(AA)[1]))]

PlotList=row.names(AA[Selecta,c("Alba","TFIIE_alpha")])

for( i in PlotList){
  
  PLOTOTOSAVE=ggplot(DNAperSpecie,aes_string(x="TFIIE_alpha",y=i,color="Specie",label="Specie"))+geom_point()+geom_text_repel(size=2)+theme_pubr()+theme(aspect.ratio = 1,legend.position="")
  
  ggsave(PLOTOTOSAVE,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/PFAM_vs_GENOME_PROP/TFIIE_alpha/TFIIE_alpha_vs_",i,".pdf",sep=""),width = 6,height=4)
}

#On tRNAsynth
Selecta=order(AA$tRNA.synt_1)[c(1:20,c((dim(AA)[1]-20):dim(AA)[1]))]

PlotList=row.names(AA[Selecta,c("Alba","tRNA.synt_1")])

for( i in PlotList){
  
  PLOTOTOSAVE=ggplot(DNAperSpecie,aes_string(x="tRNA.synt_1",y=i,color="Specie",label="Specie"))+geom_point()+geom_text_repel(size=2)+theme_pubr()+theme(aspect.ratio = 1,legend.position="")
  
  ggsave(PLOTOTOSAVE,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/PFAM_vs_GENOME_PROP/tRNA.synt_1/tRNA.synt_1_vs_",i,".pdf",sep=""),width = 6,height=4)
}





 #On DNA gyrase
Selecta=order(AA$DNA_gyraseB)[c(1:30,c((dim(AA)[1]-30):dim(AA)[1]))]

PlotList=row.names(AA[Selecta,c("Alba","DNA_gyraseB")])

for( i in PlotList){
  
  PLOTOTOSAVE=ggplot(DNAperSpecie,aes_string(x="DNA_gyraseB",y=i,color="Specie",label="Specie"))+geom_point()+geom_text_repel(size=2)+theme_pubr()+theme(aspect.ratio = 1,legend.position="")
  
  ggsave(PLOTOTOSAVE,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/PFAM_vs_GENOME_PROP/DNA_gyraseB/DNA_gyraseB_vs_",i,".pdf",sep=""),width = 6,height=4)
}



#On genome compaction
Selecta=order(AA$Compaction)[c(1:10,c((dim(AA)[1]-10):dim(AA)[1]))]

PlotList=row.names(AA[Selecta,c("Alba","Compaction")])

for( i in PlotList){
  
  PLOTOTOSAVE=ggplot(DNAperSpecie,aes_string(x="Compaction",y=i,color="Specie",label="Specie"))+geom_point()+geom_text_repel(size=2)+theme_pubr()+theme(aspect.ratio = 1,legend.position="")
  
  ggsave(PLOTOTOSAVE,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/PFAM_vs_GENOME_PROP/Compaction/Compaction_vs_",i,".pdf",sep=""),width = 6,height=4)
}




#On genome compaction
Selecta=order(AA$CompactionBp)[c(1:25,c((dim(AA)[1]-25):dim(AA)[1]))]

PlotList=row.names(AA[Selecta,c("Alba","CompactionBp")])

for( i in PlotList){
  
  PLOTOTOSAVE=ggplot(DNAperSpecie,aes_string(x="CompactionBp",y=i,color="Specie",label="Specie"))+geom_point()+geom_text_repel(size=2)+theme_pubr()+theme(aspect.ratio = 1,legend.position="")
  
  ggsave(PLOTOTOSAVE,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/PFAM_vs_GENOME_PROP/CompactionBp/CompactionBp_vs_",i,".pdf",sep=""),width = 6,height=4)
}


#On convergent intergene %
Selecta=order(AA$ConvergentIntergeneProportion)[c(1:25,c((dim(AA)[1]-25):dim(AA)[1]))]

PlotList=row.names(AA[Selecta,c("Alba","ConvergentIntergeneProportion")])

for( i in PlotList){
  
  PLOTOTOSAVE=ggplot(DNAperSpecie,aes_string(x="ConvergentIntergeneProportion",y=i,color="Specie",label="Specie"))+geom_point()+geom_text_repel(size=2)+theme_pubr()+theme(aspect.ratio = 1,legend.position="")
  
  ggsave(PLOTOTOSAVE,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/PFAM_vs_GENOME_PROP/ConvergentIntergeneProportion/ConvergentIntergeneProportion_vs_",i,".pdf",sep=""),width = 6,height=4)
}



#On NAP %
Selecta=order(AA$PercentNAP)[c(1:25,c((dim(AA)[1]-25):dim(AA)[1]))]

PlotList=row.names(AA[Selecta,c("Alba","PercentNAP")])

for( i in PlotList){
  
  PLOTOTOSAVE=ggplot(DNAperSpecie,aes_string(x="PercentNAP",y=i,color="Specie",label="Specie"))+geom_point()+geom_text_repel(size=2)+theme_pubr()+theme(aspect.ratio = 1,legend.position="")
  
  ggsave(PLOTOTOSAVE,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/PFAM_vs_GENOME_PROP/PFAM_vs_GENOME_PROP/PercentNAP/PercentNAP_vs_",i,".pdf",sep=""),width = 6,height=4)
}


#On NAP w candidate %
Selecta=order(AA$PercentNAPwCandidates)[c(1:50,c((dim(AA)[1]-50):dim(AA)[1]))]

PlotList=row.names(AA[Selecta,c("Alba","PercentNAPwCandidates")])

for( i in PlotList){
  Correlationtest=cor.test(DNAperSpecie$PercentNAPwCandidates,DNAperSpecie[,i],method="spearman")
  
  PLOTOTOSAVE=ggplot(DNAperSpecie,aes_string(x="PercentNAPwCandidates",y=i,color="SpecieShort",label="SpecieShort"))+geom_point()+geom_text_repel(size=2)+theme_pubr()+theme(aspect.ratio = 1,legend.position="")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))
  
  ggsave(PLOTOTOSAVE,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/PFAM_vs_GENOME_PROP/PercentNAPwCandidates/PercentNAPwCandidates_vs_",i,".pdf",sep=""),width = 6,height=4)
}





cor.test(DNAperSpecie$eIF.5_eIF.2B,DNAperSpecie$CBFD_NFYB_HMF,method="spearman")

ggplot(DNAperSpecie)+geom_point(aes(x=CBFD_NFYB_HMF,y=eIF.5_eIF.2B))


#Here we plot the most abundancy PFAM DNA binding categories,
#to see how easy we can spot / or not, NAPs for example

for(S in DNAperSpecie$Specie){
SubSpe=DNAperSpecie[which(DNAperSpecie$Specie==S),which(names(DNAperSpecie)%in%c("NormInt",PfamDNAbinding))]

SubSpet=as.data.frame(t(as.matrix(SubSpe)))

head(SubSpet)
names(SubSpet)="NormInt"
SubSpet$Domain=row.names(SubSpet)
SubSpet=SubSpet[order(-SubSpet$NormInt),]
head(SubSpet)
Top50=SubSpet[1:50,]
Top50Names=Top50$Domain
Top50$Domain=factor(Top50$Domain,levels = Top50$Domain)

SpecieAbundancy=ggplot()+geom_bar(data=Top50,aes_string(x="Domain",y="NormInt"),stat = "identity", width=0.7,fill="darkblue")+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6),aspect.ratio = 1/3)+ylab("% of proteome")+ggtitle(paste("Top50 DNA binding PFAM,\n",S,sep=""))
ggsave(SpecieAbundancy,filename = gsub(" ","_",paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/DNAbindingPFAM/","Top50 indiv.prot \ncontaining DNA binding PFAM",S,".pdf",sep="")))
}



#Same but instead of doing for PFAM motif, I do it for individual proteins

for(S in DNAperSpecie$Specie){
  
  SubSpe=MassSpecArchaea[which(MassSpecArchaea$Specie==S & MassSpecArchaea$PFAM_DNA_binding==T),c(1:4,which(names(MassSpecArchaea)%in%c(PfamDNAbinding,"NormInt")))]
  SubSpe=SubSpe[order(-SubSpe$NormInt),]
  PFAMCount=colSums(SubSpe[,which(names(SubSpe)%in%PfamDNAbinding)])
  
  PFAMToRemove=names(PFAMCount[which(PFAMCount==0)])
  SubSpe=SubSpe[,-which(names(SubSpe)%in%PFAMToRemove)]
  SubSpe$PFAMIDs=NA
  for(i in 1:dim(SubSpe)[1]){
    List=SubSpe[i,which(names(SubSpe)%in%PfamDNAbinding)]
    List=List[which(List>0)]
    IDs=paste(names(List),collapse ="-")
    IDs=paste(SubSpe[i,]$ID,IDs)
    SubSpe[i,]$PFAMIDs=IDs
  }
  SubSpe=SubSpe[order(-SubSpe$NormInt),]
  head(SubSpe)
  Top50=SubSpe[1:50,]
  Top50Names=Top50$PFAMIDs
  Top50$PFAMIDs=factor(Top50$PFAMIDs,levels = Top50Names)
  
  SpecieAbundancy=ggplot()+geom_bar(data=Top50,aes_string(x="PFAMIDs",y="NormInt"),stat = "identity", width=0.7,fill="darkblue")+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6),aspect.ratio = 1/3)+ylab("% of proteome")+ggtitle(paste("Top50 DNA binding PFAM,\n",S,sep=""))
  ggsave(SpecieAbundancy,filename = gsub(" ","_",paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/DNAbindingPFAM/INDIVIDUAL_PROTS/","Top50 DNA binding PFAM",S,".pdf",sep="")),height=12)
}
















#######################
#PFAM VS GO PART : 
#######################
#####################
##################
###############
############
#########
######
#An archive of GOCount can be found in my diaforarchaea script subfolder

#Loading the correspondance table : Pfam to G.O  (this is done by interpro as far as I understood)
Pfam2Go=read.table(file="/Users/ahocher/Dropbox/Laboratory/ArchealChromatinProteom/HMMs/PfamDb/Pfam2go.txt",skip=6,sep=";",stringsAsFactors = F)
head(Pfam2Go)
Pfam2Go$ID=unlist(lapply(Pfam2Go$V1,function(x) strsplit(x,split = " ")[[1]][2]))
Pfam2Go$Description=unlist(lapply(Pfam2Go$V1,function(x) strsplit(x,split = ">")[[1]][2]))
names(Pfam2Go)[2]="GO"



#For each GO compute the sum of detected proteins :  
GOCount=Pfam2Go[,c(2,4)]
GOCount=GOCount[-which(duplicated(GOCount$GO)),]
GOCount$GO=as.character(GOCount$GO)

#Table to host results : 
GOperSpecie=as.data.frame(table(MassSpecArchaea$Organism))
names(GOperSpecie)=c("Specie","NbDetectedProt")

for(MyGO in GOCount$GO){
  if(length(which(names(DNAperSpecie)%in%Pfam2Go[which(Pfam2Go$GO==MyGO),]$ID))>0){
    GOperSpecie[,MyGO]=0
    for(i in 1: dim(GOperSpecie)[1]){
      print(i)
      #First we Subselect the Specie of interest
      SubSpecie=MassSpecArchaea[which(MassSpecArchaea$Organism==GOperSpecie[i,]$Specie),]
      print(MyGO)
      #only enter the loop if there is PFAM corresponding to GO
      if(length(which(names(SubSpecie)%in%Pfam2Go[which(Pfam2Go$GO==MyGO),]$ID))>0){
        
        #Flaging all proteins in the GO : 
        #A special case is when there's only 1 PFAM ID in the GO.
        if(length(which(names(SubSpecie)%in%Pfam2Go[which(Pfam2Go$GO==MyGO),]$ID))==1){
          SubSpecie$InGO=0
          PFAMID=names(SubSpecie)[which(names(SubSpecie)%in%Pfam2Go[which(Pfam2Go$GO==MyGO),]$ID)]
          
            if(sum(SubSpecie[,PFAMID])>0){
              SubSpecie[which(SubSpecie[,PFAMID]>0),]$InGO=1}
      }else{
        #This is simply the sum of all flags for pfam that are included in the GO of interest.
        #If a protein has any of the PFAM domain, it will be included in the GOcount, but only once
        #This prevents polymerases for example to be counted multiple times in a GO category.
        SubSpecie$InGO=rowSums(SubSpecie[,which(names(SubSpecie)%in%Pfam2Go[which(Pfam2Go$GO==MyGO),]$ID)])
        }
        #Then we sum the normalized intensity of all proteins that contain a pfam in the GO
        #This counts all proteins only once 
        #If there is a protein containing a protein with a pfam in the go : 
        if(length(which(SubSpecie$InGO>0))>0){
          #We sum the normalised intensity of all those
        GOperSpecie[i,MyGO]=sum(SubSpecie[which(SubSpecie$InGO>0),]$NormInt)}
      }
    }
  }
}



#Now computing the correlations with being a Diaforarchaea : 
GOperSpecie$Diaf=0
GOperSpecie[which(GOperSpecie$Specie%in%c("Methanomassiliicoccus luminyensis","Methanomethylophilus alvus","Thermoplasma acidophilum","Thermoplasma volcanium","Picrophilus torridus")),]$Diaf=1


VarToExclude=names(which(colSums(GOperSpecie[,3:dim(GOperSpecie)[2]])==0))
if(length(VarToExclude)>0){
GOperSpecie=GOperSpecie[,-which(names(GOperSpecie)%in%VarToExclude)]}

head(GOperSpecie[order(GOperSpecie$Diaf),])

write.table(GOperSpecie,file = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/MASS_SPEC_PFAM/PerSpecie_GO_PFAM_abundancy.txt",row.names=F,quote=F,sep="\t")


#Reloaaad: 
GOperSpecie=read.table("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/MASS_SPEC_PFAM/PerSpecie_GO_PFAM_abundancy.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

Pfam2Go$GO=gsub(" ","",Pfam2Go$GO)
Pfam2Go$GO=gsub(":",".",Pfam2Go$GO)

GOCount$GO=gsub(" ","",GOCount$GO)
GOCount$GO=gsub(":",".",GOCount$GO)


GOperSpecie=merge(GOperSpecie,CandidatePerSpecie[,which(names(CandidatePerSpecie)%in%c("Specie","PercentNAPwCandidates"))],by="Specie")

#Re-prep the Species names : 
GOperSpecie$Specie[4]="Halobacterium salinarium"
GOperSpecie$Specie[16]="Sulfolobus acidocaldarius"

GOperSpecie$Specie[17]=gsub("\\^"," ",GOperSpecie$Specie[17])
GOperSpecie$Specie[14]=gsub("\\^"," ",GOperSpecie$Specie[14])

GOperSpecie$SpecieShort=unlist(lapply(GOperSpecie$Specie,function(x) paste(strsplit(strsplit(x,'')[[1]][1],' ')[[1]][1],strsplit(x," ")[[1]][2],sep=".")))

  
BB=cor(GOperSpecie[,3:608],method = "spearman")


DiafCorVal=BB[order(BB[,"Diaf"]),"Diaf"]



GOCount$Diaf=NA
for(i in GOCount$GO){
  GOCount[which(GOCount$GO==i),]$Diaf=DiafCorVal[i]}


GOCount=GOCount[order(GOCount$Diaf),]
if(length(which(is.na(GOCount$Diaf)==T))>0){
  GOCount=GOCount[-which(is.na(GOCount$Diaf)==T),]}
head(GOCount,n=25);tail(GOCount,n=25)



#We try to see how PCA on GOCount describes species: 

Goterms=Pfam2Go$GO
#Do a PCA on GOperSpecie and check it out : 
PCAInput=GOperSpecie[,which(names(GOperSpecie)%in%Goterms)]
row.names(PCAInput)=GOperSpecie$Specie
PCA=prcomp(cor(t(as.matrix(PCAInput)),method="spearman"))

library(ggfortify)
pdf("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/PCA_on_GOTerms_Correlation_and_species.pdf",width=12,height=8)
autoplot(PCA, data = GOperSpecie, colour = 'Specie',label=T,label.size=3)+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")
dev.off()


CorRes=cor(t(as.matrix(PCAInput)),method="spearman")
NewOrder=colnames(CorRes)[match(ArchTree$tip.label,colnames(CorRes))]
NewOrder=c("Halobacterium sp",NewOrder[c(1:10)],"Picrophilus torridus",NewOrder[c(11:19)])
if(length(which(is.na(NewOrder)==T))>0){
  NewOrder=NewOrder[-which(is.na(NewOrder)==T)]}
#Re-ordering the matrix : 
CorRes=CorRes[NewOrder,NewOrder]

pdf("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/Correlation_on_GOTerms_and_species.pdf",width=12,height=8)
corrplot(CorRes, method="color", cl.lim = c(0,1),diag = F)
dev.off()


#Example : 
#Plot PFAM abundancy correlation between M.alvus and M.luminiensis :
#Transformation of the table for the plot is not elegant but has been double checked.
library(reshape2)
AlvLumComp=as.data.frame(t(as.matrix(GOperSpecie[which(GOperSpecie$Specie%in%c("Methanomassiliicoccus luminyensis","Methanomethylophilus alvus")),c(1,3:592)])))
names(AlvLumComp)=unlist(AlvLumComp[1,])
AlvLumComp=AlvLumComp[-1,]
AlvLumComp$`Methanomethylophilus alvus`=as.numeric(as.character(AlvLumComp$`Methanomethylophilus alvus`))
AlvLumComp$`Methanomassiliicoccus luminyensis`=as.numeric(as.character(AlvLumComp$`Methanomassiliicoccus luminyensis`))

Correlationtest=cor.test(AlvLumComp$`Methanomethylophilus alvus`,AlvLumComp$`Methanomassiliicoccus luminyensis`,method="spearman")

AlvuLumComp=ggplot(data = AlvLumComp)+geom_point(aes(x=`Methanomethylophilus alvus`,y=`Methanomassiliicoccus luminyensis`),color="grey50")+scale_x_log10()+scale_y_log10()+geom_abline(intercept = 0)+xlab("M.alvus, % of proteome\nper G.O")+ylab("M.luminyensis, % of proteome\nper G.O")+theme_pubr()+theme(aspect.ratio = 1)+ggtitle(paste("Spearman Cor ", round(Correlationtest$estimate,2), " p-value :", round(Correlationtest$p.value,3),sep=""))
ggsave(AlvuLumComp,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/Comparison_Lumi_Alvus_PFAMtoGO.pdf",sep=""),width = 4,height=4)





GOCount$pvalue=NA

for(i in GOCount$GO){
  if(i%in%names(GOperSpecie)==T){
  Test=t.test(GOperSpecie[,i]~GOperSpecie$Diaf)
  GOCount[which(GOCount$GO==i),]$pvalue=Test$p.value}}


#ANOVA TO SEE WHICH GO DONT CHANGE : 
GOCount$Stdev=NA
GOCount$Average=NA

for(i in GOCount$GO){
  if(i%in%names(GOperSpecie)==T){
    Test=t.test(GOperSpecie[,i]~GOperSpecie$Diaf)
    GOCount[which(GOCount$GO==i),]$pvalue=Test$p.value
  GOCount[which(GOCount$GO==i),]$Stdev=sd(GOperSpecie[,i])
  GOCount[which(GOCount$GO==i),]$Average=mean(GOperSpecie[,i])
  
  }}


# Export data
write.table(GOCount,file = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/MASS_SPEC_PFAM/CorrelationDiafor_GO_PFAM_abundancy.txt",row.names=F,quote=F,sep="\t")


GOCount=read.delim(file = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/MASS_SPEC_PFAM/CorrelationDiafor_GO_PFAM_abundancy.txt",header=T,sep="\t",quote="",stringsAsFactors = F)



# GOCount=GOCount[-which(GOCount$pvalue>0.05),]

GOCount=GOCount[order(GOCount$pvalue),]
GOCount=GOCount[-which(GOCount$GO=="Pfam:PF02654 CobS "),]
#Here we produce ( for exploration purposes, various GO vs GO )
library(ggrepel)
Xvar="GO.0006351"
Yvar="GO.0005840"
for(i in GOCount$GO[1:10]){
  for(j in GOCount$GO[1:10]){
    if(i!=j){
    Xvar=i;Yvar=j
    
    if(file.exists(paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/PFAM_GO/",Yvar,"_vs_Yaxis_",Xvar,"_percent_proteome.pdf",sep=""))==F){
Comparison=ggplot(data=GOperSpecie,aes_string(x=Xvar,y=Yvar,color="as.factor(Diaf)",size="NbDetectedProt",label="Specie"))+geom_point()+geom_text_repel(size=1.5)+xlab(label = paste(GOCount[which(GOCount$GO==Xvar),]$Description," (%)",sep=""))+ylab(label = paste(GOCount[which(GOCount$GO==Yvar),]$Description," (%)",sep=""))+theme_pubr()+theme(aspect.ratio=1)+scale_color_manual(values=c("forestgreen","royalblue"))
ggsave(plot = Comparison,filename = paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/PFAM_GO/",Xvar,"_vs_Yaxis_",Yvar,"_percent_proteome.pdf",sep=""),width=5,height=5)}}}
}




#Methanogenesis : 
Xvar="GO.0006351"
#Membrane :
Yvar="GO.0005840"

Comparison=ggplot(data=GOperSpecie,aes_string(x=Xvar,y=Yvar,color="as.factor(Diaf)",size="NbDetectedProt",label="Specie"))+geom_point()+geom_text_repel(size=1.5)+xlab(label = paste(GOCount[which(GOCount$GO==Xvar),]$Description," (%)",sep=""))+ylab(label = paste(GOCount[which(GOCount$GO==Yvar),]$Description," (%)",sep=""))+theme_pubr()+theme(aspect.ratio=1)+scale_color_manual(values=c("forestgreen","royalblue"))

ggsave(plot = Comparison,filename = paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/PFAM_GO/",Xvar,"_vs_Yaxis_",Yvar,"_percent_proteome.pdf",sep=""),width=5,height=5)



#Focus on supercoiling : 
NN=names(BB[c(head(order(BB[,"GO.0000786"]),n=10),tail(order(BB[,"GO.0000786"])),n=10),"GO.0000786"])

GOCount[which(GOCount$GO%in%NN),]

#Nucleosome : 
Xvar="GO.0000786"
#Cell septum :
Yvar="GO.0055085"

Comparison=ggplot(data=GOperSpecie,aes_string(x=Xvar,y=Yvar,color="as.factor(Diaf)",size="NbDetectedProt",label="Specie"))+geom_point()+geom_text_repel(size=1.5)+xlab(label = paste(GOCount[which(GOCount$GO==Xvar),]$Description," (%)",sep=""))+ylab(label = paste(GOCount[which(GOCount$GO==Yvar),]$Description," (%)",sep=""))+theme_pubr()+theme(aspect.ratio=1)+scale_color_manual(values=c("forestgreen","royalblue"))

ggsave(plot = Comparison,filename = paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/",Xvar,"_vs_Yaxis_",Yvar,"_percent_proteome.pdf",sep=""),width=5,height=5)




TFPlot=ggplot()+geom_bar(data=GOperSpecie[,c("GO.0003700","Specie")],aes(x=Specie,y=GO.0003700),stat = "identity", width=0.7,fill="darkblue")+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of proteome")+ggtitle(paste(GOCount[which(GOCount$GO=="GO.0003700"),]$Description," abundancy accross species",sep=""))

ggsave( plot= TFPlot,filename = paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/GO_Transcription_factor_accross_species.pdf",sep=""),width=6,height=4)




NN=names(BB[c(head(order(BB[,"GO.0005840"])),tail(order(BB[,"GO.0005840"]))),"GO.0005840"])

#Nucleosome : 
Xvar="GO.0000786"
#Coiling :
Yvar="GO.0090529"

Comparison=ggplot(data=GOperSpecie,aes_string(x=Xvar,y=Yvar,color="as.factor(Diaf)",size="NbDetectedProt",label="Specie"))+geom_point()+geom_text_repel(size=1.5)+xlab(label = paste(GOCount[which(GOCount$GO==Xvar),]$Description," (%)",sep=""))+ylab(label = paste(GOCount[which(GOCount$GO==Yvar),]$Description," (%)",sep=""))+theme_pubr()+theme(aspect.ratio=1)+scale_color_manual(values=c("forestgreen","royalblue"))

ggsave(plot = Comparison,filename = paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/",Xvar,"_vs_Yaxis_",Yvar,"_percent_proteome.pdf",sep=""),width=5,height=5)

#Correlations between genoe properties and GOterms : 
GOperSpecie=merge(GOperSpecie,DNAperSpecie[,which(names(DNAperSpecie)%in%c("Specie","DivergentIntergeneProportion","ConvergentIntergeneProportion","AlignedIntergeneProportion","CompactionBp","GeneNb","NbDetectedProt","Diaf","Tss.Proportion","PercentNAP","GC.content.avg"))],by="Specie")

VariableList=c("DivergentIntergeneProportion","ConvergentIntergeneProportion","AlignedIntergeneProportion","CompactionBp","GeneNb","NbDetectedProt","Diaf","Tss.Proportion","PercentNAP","GC.content.avg","PercentNAPwCandidates",GOCount$GO)
BB=as.data.frame(cor(GOperSpecie[,which(names(GOperSpecie)%in%c(VariableList))],method = "spearman"))


GONAPCorVal=as.data.frame(BB[order(BB[,"PercentNAP"]),c("PercentNAP","CompactionBp")])


head(GONAPCorVal)
GONAPCorVal$GOID=row.names(GONAPCorVal)

GONAPCorVal$Description=NA
for(i in 1:dim(GONAPCorVal)[1]){
  if(length(which(GOCount$GO==GONAPCorVal[i,]$GOID))>0){
  GONAPCorVal[i,]$Description=GOCount[which(GOCount$GO==GONAPCorVal[i,]$GOID),]$Description}}
  


#
library(ggrepel)

for(i in c(head(GONAPCorVal$GOID,n=20),tail(GONAPCorVal$GOID,n=20))){
  CorTest=cor.test(GOperSpecie$PercentNAP,GOperSpecie[,which(names(GOperSpecie)==i)],method="spearman")
        Comparison=ggplot(data=GOperSpecie,aes_string(x="PercentNAP",y=i,color="Specie",label="Specie"))+geom_point()+geom_text_repel(size=1.5)+xlab(label = "PercentNAP")+ylab(label = paste(GONAPCorVal[which(GONAPCorVal$GOID==i),]$Description," (%)",sep=""))+theme_pubr()+theme(aspect.ratio=1,legend.position = "")+scale_color_manual(values = col_vector)+ggtitle(paste("Rho",round(CorTest$estimate,digits=2),"p=",round(CorTest$p.value,digits = 2)))
        ggsave(plot = Comparison,filename = paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/GO_VS_GENOMEPROP/PFAM_GO_NAPs/",i,"_vs_percentNAP_percent_proteome.pdf",sep=""),width=5,height=5)}


#With NAP candidates : 

GONAPCorVal=as.data.frame(BB[order(BB[,"PercentNAPwCandidates"]),c("PercentNAPwCandidates","CompactionBp")])


head(GONAPCorVal)
GONAPCorVal$GOID=row.names(GONAPCorVal)

GONAPCorVal$Description=NA
for(i in 1:dim(GONAPCorVal)[1]){
  if(length(which(GOCount$GO==GONAPCorVal[i,]$GOID))>0){
    GONAPCorVal[i,]$Description=GOCount[which(GOCount$GO==GONAPCorVal[i,]$GOID),]$Description}}



GONAPCorVal=GONAPCorVal[order(GONAPCorVal$CompactionBp),]
for(i in c(head(GONAPCorVal$GOID,n=40),tail(GONAPCorVal$GOID,n=40))){
  CorTest=cor.test(GOperSpecie$PercentNAPwCandidates,GOperSpecie[,which(names(GOperSpecie)==i)],method="spearman")
  options(scipen = 2)
  Comparison=ggplot(data=GOperSpecie,aes_string(x="PercentNAPwCandidates",y=i,color="Specie",label="Specie"))+geom_point()+geom_text_repel(size=1.5)+xlab(label = "PercentNAPwCandidates")+ylab(label = paste(GONAPCorVal[which(GONAPCorVal$GOID==i),]$Description," (%)",sep=""))+theme_pubr()+theme(aspect.ratio=1,legend.position = "")+scale_color_manual(values = col_vector)+ggtitle(paste("Rho",round(CorTest$estimate,digits=2),"p=",CorTest$p.value))
  
  if(CorTest$p.value<0.01){
  ggsave(plot = Comparison,filename = paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/GO_VS_GENOMEPROP/PFAM_GO_NAPs_wcandidates/",i,"_vs_PercentNAPwCandidates_percent_proteome.pdf",sep=""),width=5,height=5)}
    }


CC=BB
for(i in 1:dim(CC)[1]){
  if(length(which(GOCount$GO==row.names(CC)[i]))>0){
  row.names(CC)[i]=GOCount[which(GOCount$GO==row.names(CC)[i]),]$Description
}}
#HEatmap of hits :
Selecta=order(CC$PercentNAPwCandidates)[c(1:25,c((dim(CC)[1]-25):dim(CC)[1]))]
pdf("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/GO_VS_GENOMEPROP/Correlation_GO_Top50_SumNAPwcandidates.pdf",width=12,height=12)
heatmap3::heatmap3(as.matrix(CC[Selecta,Selecta]),scale = "none",cexRow = 0.35,cexCol  = 0.35,col=viridis(100))
dev.off()






i="GO.0006265"
CorTest=cor.test(GOperSpecie$PercentNAP,GOperSpecie[,which(names(GOperSpecie)==i)],method="spearman")
Comparison=ggplot(data=GOperSpecie,aes_string(x="PercentNAP",y=i,color="Specie",label="Specie"))+geom_point()+geom_text_repel(size=1.5)+xlab(label = "PercentNAP")+ylab(label = paste(GONAPCorVal[which(GONAPCorVal$GOID==i),]$Description," (%)",sep=""))+theme_pubr()+theme(aspect.ratio=1,legend.position = "")+scale_color_manual(values = col_vector)+ggtitle(paste("Rho",round(CorTest$estimate,digits=2),"p=",round(CorTest$p.value,digits = 2)))

GOCount[which(GOCount$Description==" GO:helicase activity"),]

GOCount[grep("helicase",GOCount$Description),]
GONAPCorVal[grep("helicase",GONAPCorVal$Description),]

Pfam2Go[which(Pfam2Go$GO=="GO.0004386"),]





#ARCHIVE : 

#Checking the percentage of proteins that I don't include because they don't correspond to genbank ids from local database.
table(is.na(MassSpecArchaea[which(MassSpecArchaea$Organism=="Halobacterium sp"),]$ID))
126/(1644+126)

for(S in DNAperSpecie$Specie){
Check=MassSpecArchaea[which(is.na(MassSpecArchaea$ID)==T & MassSpecArchaea$Organism==S),c("Protein_IDs","NormInt")]
print(S)
print(sum(Check$NormInt))}


MassSpecArchaea[which(MassSpecArchaea$Organism=="Halobacterium sp" & MassSpecArchaea$HTH_11==1),c("ID","NormInt")]


#Archive : 
MassSpecArchaea[which(MassSpecArchaea$CBFD_NFYB_HMF==1),1:3]
