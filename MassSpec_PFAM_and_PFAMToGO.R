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

#Load annotated mass spec data : 
#has to be unzipped first from Hocher_2022_MassSpecMeasurementsAccrossSpecies_wPFAM_Annotations.txt.zip
#
#Important note : This table contains all species, including the ones that were added during revisions (C. divulgatum)
#Some graphs in the manuscript are computed without Cuniculiplasma, which was added later as a confirmatory experiment
#I added C. divulgatum data for the sake of completion and reproducibility,
#However, to reproduce graphs from the manuscript one should first remove C. divulgatum data
MassSpecArchaea=read.table(file="Hocher_2022_MassSpecMeasurementsAccrossSpecies_wPFAM_Annotations.txt",header=T,sep="\t",quote="",stringsAsFactors = F)



#Loading archaeal genomes ids :  
Genomes=read.table("Annotated_genomes_info_to_retrieve_mass_spec.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

#Revisions : adding Cuniculiplasma divulgatum, which was not in my initial database : 
CdivGenomes=read.table("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/Single_Genomes/CUNICULIPLASMA_DIVULGATUM/AnnotatedGenomes/Annotated_genomes_info.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
CdivGenomes$To_retrieve_MassSpec=1;CdivGenomes$Specie_name_orig_dataset="Cuniculiplasma divulgatum"
Genomes=rbind(Genomes,CdivGenomes)

#Adding all the external mass spec data which were obtained not from Muller et al. :
Genomes[which(Genomes$Specie=="Thermococcus^kodakarensis^KOD1"),]$To_retrieve_MassSpec=1
Genomes[which(Genomes$Specie=="Thermococcus^kodakarensis^KOD1"),]$Specie_name_orig_dataset="Thermococcus^kodakarensis^KOD1"

Genomes[which(Genomes$Specie=="Haloferax^volcanii^DS2"),]$To_retrieve_MassSpec=1
Genomes[which(Genomes$Specie=="Haloferax^volcanii^DS2"),]$Specie_name_orig_dataset="Haloferax volcanii"

#
Genomes[which(Genomes$Specie=="Natrialba^magadii^ATCC^43099"),]$To_retrieve_MassSpec=1
Genomes[which(Genomes$Specie=="Natrialba^magadii^ATCC^43099"),]$Specie_name_orig_dataset="Natrialba magadii"

Genomes[which(Genomes$Specie=="Nitrosopumilus^maritimus^SCM1"),]$To_retrieve_MassSpec=1
Genomes[which(Genomes$Specie=="Nitrosopumilus^maritimus^SCM1"),]$Specie_name_orig_dataset="Nitrosopumilus^maritimus^SCM1"

#This is simply for compatibility later on
Genomes=Genomes[,-which(names(Genomes)=="Specie")]




#Adding genome properties :
Genomesppties=read.table("Archaea_EBMC2_genomes_properties_results.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
Genomesppties=Genomesppties[,-which(names(Genomesppties)=="Specie")]

#Adding the Operonic properties for each gene / protein
OperonPred=read.table("Archaea_EBMC2_Operons_properties_results.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
OperonPred=OperonPred[,which(names(OperonPred)%in%c("Assembly","Tss.Proportion"))]

Genomesppties=merge(Genomesppties,OperonPred,by="Assembly")

#Exporting a list of assemblies : 
Assemblies=Genomes[which(Genomes$To_retrieve_MassSpec==1),]$Assembly



#Annotate DNA binding Pfam and NAPs

#Reloading the manually annotated table :
PfamDNAbinding=read.table(file="PFAM_DNA_binding_list_present_MassSpecDataSet_annot_AH.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

TFHmm=PfamDNAbinding[which(PfamDNAbinding$TF==1),]$PfamID
PfamDNAbinding=PfamDNAbinding$PfamID

NAPsHMM=c("Cren7","X7kD_DNA_binding","CBFD_NFYB_HMF","Alba","Bac_DNA_binding","MC1","Histone")
AllPfamID=make.names(read.table("All_pfamID.txt")$V1)
AllPfamID=AllPfamID[-which(AllPfamID=="ID")]

#There is a PFAM to GO Annotation -GO:0003700-, we load it.
Pfam2Go=read.table(file="Pfam2go.txt",skip=6,sep=";",stringsAsFactors = F)
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


#Adding CC1 at last : 
CC1=read.table(file = "Archaea_Distribution_CC1_Proteins_with_Length.txt",header=T,sep="\t",stringsAsFactors = F)
#Cut-off on length to limit false positive (CC1 is typically very short)
CC1=CC1[which(CC1$Length<200),]
MassSpecArchaea$CC1=0
MassSpecArchaea[which(MassSpecArchaea$ID%in%CC1$ProtName),]$CC1=1


MassSpecArchaea$NAP=rowSums(MassSpecArchaea[,which(names(MassSpecArchaea)%in%NAPsHMM)])
#Binarize
MassSpecArchaea[which(MassSpecArchaea$NAP>0),]$NAP=1

#TF : 
#Annotating Transcription factors : 
MassSpecArchaea$TFPfam=0
MassSpecArchaea[which(rowSums(MassSpecArchaea[,which(names(MassSpecArchaea)%in%TF_HMMs)])>0),]$TFPfam=1



#PART II Counting all proteins with a given PFAM domains to have a metric to compare PFAM domains : 


#Building a table to host results : 
DNAperSpecie=as.data.frame(table(MassSpecArchaea$Organism))
names(DNAperSpecie)=c("Specie","NbDetectedProt")
DNAperSpecie=merge(DNAperSpecie,Genomes,by.x="Specie",by.y="Specie_name_orig_dataset",all.x=T)
DNAperSpecie=merge(DNAperSpecie,Genomesppties,by=c("Assembly","Taxid"),all.x=T)


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


DNAperSpecie[,which(names(DNAperSpecie)%in%c("Specie","PercentNAP"))]

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
write.table(DNAperSpecie,file = "Hocher_2022_Revisions_PerSpeciePFAM_abundancy.txt",row.names=F,quote=F,sep="\t")



##############
#Re-load
##############
DNAperSpecie=read.table("Hocher_2022_Revisions_PerSpeciePFAM_abundancy.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

#Add the NAP candidates ( NAP detection scheme (Predict_NAP_from_MassSpecdata.R) has to be run 1st): 
CandidatePerSpecie=read.table("PerSpecie_CandidateNapCount_and_abundancy.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
CandidatePerSpecie$PercentNAPwCandidates=CandidatePerSpecie$PercentNAP+CandidatePerSpecie$TotalProportionCandidates
CandidatePerSpecie=CandidatePerSpecie[-which(names(CandidatePerSpecie)=="PercentNAP")]

DNAperSpecie=merge(DNAperSpecie,CandidatePerSpecie,by="Specie")

#Re formating names (first due to Spelling mistakes in Muller et al. (I verified that Halo sp is indeed Halobacterium salinarium as stated in the sup and from the correspondence with unipro))
DNAperSpecie$Specie[which(DNAperSpecie$Specie=="Halobacterium sp")]="Halobacterium salinarium"
DNAperSpecie$Specie[which(DNAperSpecie$Specie=="Sulfolobus acidocaldicarius")]="Sulfolobus acidocaldarius"

#This is simply some formating to enable strsplit function bellow
DNAperSpecie$Specie[17]=gsub("\\^"," ",DNAperSpecie$Specie[17])
DNAperSpecie$Specie[14]=gsub("\\^"," ",DNAperSpecie$Specie[14])

DNAperSpecie$SpecieShort=unlist(lapply(DNAperSpecie$Specie,function(x) paste(strsplit(strsplit(x,'')[[1]][1],' ')[[1]][1],strsplit(x," ")[[1]][2],sep=".")))




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






###################
#Specific plot on histones and other single PFAM models: 


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


#On Alba and histones, exploratory plots: 


Selecta=order(AA$Alba)[c(1:30,c((dim(AA)[1]-30):dim(AA)[1]))]

PlotList=row.names(AA[Selecta,c("Alba","CBFD_NFYB_HMF")])

for( i in PlotList){
  
  PLOTOTOSAVE=ggplot(DNAperSpecie,aes_string(x="Alba",y=i,color="Specie",label="Specie"))+geom_point()+geom_text_repel(size=2)+theme_pubr()+theme(aspect.ratio = 1,legend.position="")
  
  ggsave(PLOTOTOSAVE,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/PFAM_vs_GENOME_PROP/ALBA/Alba_vs_",i,".pdf",sep=""),width = 6,height=4)
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
Pfam2Go=read.table(file="Pfam2go.txt",skip=6,sep=";",stringsAsFactors = F)
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


#Removing GO that are empty (without measurements)
VarToExclude=names(which(colSums(GOperSpecie[,3:dim(GOperSpecie)[2]])==0))
if(length(VarToExclude)>0){
GOperSpecie=GOperSpecie[,-which(names(GOperSpecie)%in%VarToExclude)]}


write.table(GOperSpecie,file = "PerSpecie_GO_PFAM_abundancy.txt",row.names=F,quote=F,sep="\t")





#Revisions : double checking what RNA polymerase subunits are used for normalisation : 
#
MassSpecArchaea[which(MassSpecArchaea$RNA_pol_Rpb2_3==1),which(names(MassSpecArchaea)%in%c("ID","Organism","NormInt"))]

MassSpecArchaea[which(MassSpecArchaea$RNA_pol_Rpb1_3==1),which(names(MassSpecArchaea)%in%c("ID","Organism","NormInt"))]

MassSpecArchaea[which(MassSpecArchaea$RNA_pol_A_bac==1),which(names(MassSpecArchaea)%in%c("ID","Organism","NormInt"))]

