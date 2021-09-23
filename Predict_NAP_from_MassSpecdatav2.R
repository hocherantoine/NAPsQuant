#####################################################
## Project: diaforarchaea chromatin
## Analysis of external data from Muller et al. 2020 DOI : 
## https://doi.org/10.1038/s41586-020-2402-x
## in relationship with pFAM domains.
## Aim : 
## See if there's anything special about Diaforarchaea
## This scripts input is from MassSpec_NAPs_detection
## 
## 
## Date: January 2021
## Author: Antoine Hocher
####################################################

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.

library(rhmmer)
source("/Users/ahocher/Dropbox/Scripts/PLOT_LIBRARY/AH_plot_library.R")

Sys.setenv(PATH="/Users/ahocher/opt/miniconda3/bin:/Users/ahocher/opt/miniconda3/condabin:/opt/local/bin:/opt/local/sbin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin")

#Load annotated mass spec data : 
MassSpecArchaea=read.table(file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/MASS_SPEC_Muller2020/41586_2020_2402_MOESM3_ESM_archaeaOnly_Alvus_Lumi_Kodak_natrialba_Haloferax_Nmaritimus_With_PFam_annotations.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
MassSpecArchaea$Specie=MassSpecArchaea$Organism
#AddingCC1 :
CC1=read.table(file = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/HMMSEARCH/CC1/Archaea_Distribution_CC1_Proteins_with_Length.txt",header=T,sep="\t",stringsAsFactors = F)
CC1=CC1[which(CC1$Length<200),]
MassSpecArchaea$CC1=0
MassSpecArchaea[which(MassSpecArchaea$ID%in%CC1$ProtName),]$CC1=1


OperonData=read.delim(file="/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/OPERO_MAPPER/ARCHAEA/Processed/AllOperons_MassSpec_Combined.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
names(OperonData)[3]="ORF_type"
OperonData=OperonData[which(is.na(OperonData$IdGeneReduced)==F),]
OperonData=OperonData[which(duplicated(OperonData$ProtName)==F),]
#not in operon :
NotInOperon=OperonData[which(OperonData$OperonSize==1),]$ProtName

#Load HMMnames of DNA binding Pfam and NAPs
PfamDNAbinding=read.table(file="/Users/ahocher/Dropbox/Laboratory/ArchealChromatinProteom/HMMs/DBPome/PFAM_DNA_binding_list.txt",stringsAsFactors = F)$V1

PfamDNAbinding=c(PfamDNAbinding,"TFIIE_alpha")


Restrict=PfamDNAbinding[which(make.names(PfamDNAbinding)%in%make.names(names(MassSpecArchaea)))]

#Exporting this to manually annotate everything
write.table(Restrict,"/Users/ahocher/Dropbox/Laboratory/ArchealChromatinProteom/HMMs/DBPome/PFAM_DNA_binding_list_present_MassSpecDataSet.txt",quote=F,sep="\t",row.names = F,col.names = F)



#Reloading the manually annotated table :
PfamDNAbinding=read.table(file="/Users/ahocher/Dropbox/Laboratory/ArchealChromatinProteom/HMMs/DBPome/PFAM_DNA_binding_list_present_MassSpecDataSet_annot_AH.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

TFHmm=PfamDNAbinding[which(PfamDNAbinding$TF==1),]$PfamID

NAPsHMM=c("Cren7","X7kD_DNA_binding","CBFD_NFYB_HMF","Alba","Bac_DNA_binding","MC1","Histone","CC1")











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

#Exporting all protein names to 
#retrieve their properties : 

#Loading Pfam HMM search :
Arch_Res=as.data.frame(read_tblout("/Users/ahocher/Dropbox/Laboratory/Final_analysis_Archaeal_chromatin/HMM_RESULTS/concatArchaeaEBMC2_vs_Pfam.txt"))

#Loading archaeal genomes ids to retieve assemblies :  
Genomes=read.table("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/MASS_SPEC_Muller2020/Annotated_genomes_info_to_retrieve_mass_spec.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

#Adding all the external mass spec species  :
Genomes[which(Genomes$Specie=="Thermococcus^kodakarensis^KOD1"),]$To_retrieve_MassSpec=1
Genomes[which(Genomes$Specie=="Thermococcus^kodakarensis^KOD1"),]$Specie_name_orig_dataset="Thermococcus^kodakarensis^KOD1"

Genomes[which(Genomes$Specie=="Haloferax^volcanii^DS2"),]$To_retrieve_MassSpec=1
Genomes[which(Genomes$Specie=="Haloferax^volcanii^DS2"),]$Specie_name_orig_dataset="Haloferax volcanii"

Genomes[which(Genomes$Specie=="Natrialba^magadii^ATCC^43099"),]$To_retrieve_MassSpec=1
Genomes[which(Genomes$Specie=="Natrialba^magadii^ATCC^43099"),]$Specie_name_orig_dataset="Natrialba magadii"


Genomes[which(Genomes$Specie=="Nitrosopumilus^maritimus^SCM1"),]$To_retrieve_MassSpec=1
Genomes[which(Genomes$Specie=="Nitrosopumilus^maritimus^SCM1"),]$Specie_name_orig_dataset="Nitrosopumilus^maritimus^SCM1"

Assemblies=Genomes[which(Genomes$Specie_name_orig_dataset%in%unique(MassSpecArchaea$Organism)),]$Assembly




#Because not all proteins are PFAM hits, we go back to original genomes : 
if(file.exists("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/NAP_PREDICTION/Archaea_all_species_MassSpec_proteins.fa")==F){
setwd("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/DB_Archaea95_update012020/genome_assemblies_prot_fasta/ncbi-genomes-2020-05-14")
Proteomes=list.files(pattern = "*.faa$")

ProteomesToRetrieve=c()
for(A in Assemblies){
ProteomesToRetrieve=c(ProteomesToRetrieve,Proteomes[grep(A,Proteomes)])
}

Prot=c()
for(P in ProteomesToRetrieve){
  Prot=c(Prot,read.fasta(P))
}

write.fasta(Prot,getName(Prot),file.out = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/NAP_PREDICTION/Archaea_all_species_MassSpec_proteins.fa")
}



#We compute the A.A properties of the whole genome.
source("/Users/ahocher/Dropbox/Scripts/ComputeProtStats_from_Fasta.R")
if(file.exists("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/NAP_PREDICTION/Archaea_all_species_MassSpec_proteinsPpties.txt")==F){
  MSProtStats=ComputeStatsFromFasta_ProtName(FastaFilePath = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/NAP_PREDICTION/Archaea_all_species_MassSpec_proteins.fa")
  
  
  MSProtStats=MSProtStats[,-which(names(MSProtStats)=="Type")]
  write.table(MSProtStats,"/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/NAP_PREDICTION/Archaea_all_species_MassSpec_proteinsPpties.txt",sep="\t",quote=F,row.names = F,col.names =T)
  
  
}





#Annotating proteins unnanotated by Pfam :
MassSpecArchaea$NoPfam=0
MassSpecArchaea[which(rowSums(MassSpecArchaea[,which(names(MassSpecArchaea)%in%AllPfamID)])==0),]$NoPfam=1


#Annotating Transcription factors : 
MassSpecArchaea$TFPfam=0
MassSpecArchaea[which(rowSums(MassSpecArchaea[,which(names(MassSpecArchaea)%in%TF_HMMs)])>0),]$TFPfam=1


MassSpecArchaea$Type="other"
MassSpecArchaea[which(MassSpecArchaea$TFPfam==1),]$Type="TF"
MassSpecArchaea[which(MassSpecArchaea$NAP==1),]$Type="NAP"
MassSpecArchaea$Type=factor(MassSpecArchaea$Type,levels = c("other","TF","NAP"))



#Reloading
MSProtStats=read.table("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/NAP_PREDICTION/Archaea_all_species_MassSpec_proteinsPpties.txt",sep="\t",header=T,quote="",stringsAsFactors = F)

names(MSProtStats)[1]="ID" 

MassSpecArchaeaGB=MassSpecArchaea[which(is.na(MassSpecArchaea$ID)==F),]

MassSpecArchaeaAnno=merge(MassSpecArchaeaGB,MSProtStats,by="ID",all.x=T)




#Have a sensible size cut-off for NAPs definition :
#The maximum of known NAP in archaea :
max(MassSpecArchaeaAnno[which(MassSpecArchaeaAnno$NAP==1),]$Length)
#This is 154 base pairs; however, TrmBl2, which is also a NAP in Kodakarensis
#but is not included in my list, is 264 a.a, thus, we will use 264 a.a


################START OF THE PREDICTION SCHEME : 
###
#Are all NAPs not in operon :

NAPsId=MassSpecArchaea[which(MassSpecArchaea$NAP==1),]$ID



MassSpecArchaea.m=merge(OperonData,MassSpecArchaea,by.x="ProtName",by.y="ID")




#Examplifying the few NAPs that are not in Operons : 
MassSpecArchaea.m[which(MassSpecArchaea.m$OperonSize>1 & MassSpecArchaea.m$NAP==1),c(1:3,which(names(MassSpecArchaea.m)%in%c("Organism",NAPsHMM)))]

dim(MassSpecArchaea.m)
table(MassSpecArchaea.m$OperonSize)
table(MassSpecArchaea.m[which(MassSpecArchaea.m$NAP==1),]$OperonSize)
length(which(MassSpecArchaea.m$NAP==1))

table(MassSpecArchaea.m[which(MassSpecArchaea.m$TFPfam==1),]$OperonSize)
length(which(MassSpecArchaea.m$TFPfam==1))
#I computed p-value of hypergeometric law for : NAPs are more oftent NOT in operon : 
#it gave 1e-9 https://stattrek.com/online-calculator/hypergeometric.aspx for NAP vs All
#0 for TF vs ALL
#and 0.00012 for TF vs NAP ( nap more oftent not in operon)


write.table(MassSpecArchaea.m[which(MassSpecArchaea.m$OperonSize>1 & MassSpecArchaea.m$NAP==1),c(1:3,which(names(MassSpecArchaea.m)%in%c("Assembly","Organism",NAPsHMM)))],file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/NAP_PREDICTION/Naps_that_are_in_Operons.txt",row.names=F,sep="\t",quote=F)

#Safety check : orf not in Operons, are not detected by Operon Mapper, it is a minor fraction and should not affect the result ( I manually checked this)


100*dim(MassSpecArchaea[which(!(MassSpecArchaea$ID%in%OperonData$ProtName)),c(1:5)])[1]/dim(MassSpecArchaea)[1]
#it's 0.3% of all data.


OperonDistrib=ggplot(MassSpecArchaea.m)+geom_boxplot(aes(x=as.factor(Type),y=OperonSize,color=as.factor(Type),group=as.factor(Type)),varwidth=F,outlier.size = 0.7)+scale_color_manual(values=col_vector)+theme_pubr()+theme(aspect.ratio = 1,axis.ticks.x = element_blank(),axis.line.x = element_blank(),legend.position = "")+xlab("")+ylab("Operon Size (# genes)")+ggtitle("p-value : 1e-9")

ggsave(OperonDistrib,file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/Comparison_Operon_TF_NAPs.pdf",width=3,height=3)

#PART I : Informed Size selection 
#Only takes proteins smaller than TrmBl2 +10 % : 290 a.a
MassSpecArchaeaAnno$Specie=MassSpecArchaeaAnno$Organism

MSmall=MassSpecArchaeaAnno[which(MassSpecArchaeaAnno$Length<290),]


#Storing the proteins without pfam domain for latter use : 
CandidateNapsNoPfam=MSmall[which(MSmall$NoPfam==1 ),]


#Part II compute the median abundancy of transcription factors, check if it's different from specie to specie.
TranscriptionFactors=MSmall[which(MSmall$TFPfam==1),]

TFcount=table(TranscriptionFactors$Organism)

# export TF sequences and candidate no PFAM for DNA binder analysis : 
library(seqinr);library(ape)
AllProteins=read.fasta("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/NAP_PREDICTION/Archaea_all_species_MassSpec_proteins.fa")

TFProteins=AllProteins[which(getName(AllProteins)%in%TranscriptionFactors$ID)]

write.fasta(sequences = TFProteins,names = getName(TFProteins),file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/NAP_PREDICTION/Archaea_all_species_MassSpec_proteins_PFAMTF_only.fa")

#Not in pfam : 
NoPFAMProteins=AllProteins[which(getName(AllProteins)%in%CandidateNapsNoPfam$ID)]

write.fasta(sequences = NoPFAMProteins,names = getName(NoPFAMProteins),file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/NAP_PREDICTION/Archaea_all_species_MassSpec_proteins_NoPFAMProteins.fa")


TFDNAbinder=read.table("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/NAP_PREDICTION/DNABinder_Archaea_all_species_MassSpec_proteins_PFAMTF_only.txt",header = T,sep="\t",stringsAsFactors = F,quote="")
TFquartile=quantile(TFDNAbinder$SVM.score,probs=seq(0,1,by=0.05))

NoPfamDNAbinder=read.table("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/NAP_PREDICTION/DNABinder_Archaea_all_species_MassSpec_proteins_NoPFAMProteins.txt",header = T,sep="\t",stringsAsFactors = F,quote="")
NoPfamquartile=quantile(NoPfamDNAbinder$SVM.score,probs=seq(0,1,by=0.05))

#Compare score : 
TFDNAbinder$TypeBinding="TF"
NoPfamDNAbinder$TypeBinding="NoPfam"
DNAbinderRes=rbind(TFDNAbinder,NoPfamDNAbinder)

library(ggplot2)

ggplot(data=as.data.frame(DNAbinderRes))+geom_histogram(aes(x=as.numeric(SVM.score),fill=as.factor(TypeBinding)),position="identity")




#Removing all proteins that do not qualify as DNA binding :
SVM.Threshold=TFquartile[5]
#TFToKeep=TFDNAbinder[which(TFDNAbinder$SVM.score>SVM.Threshold),]$ID
NoPfamToKeep=NoPfamDNAbinder[which(NoPfamDNAbinder$SVM.score>0),]$ID

#TranscriptionFactors=TranscriptionFactors[which(TranscriptionFactors$ID%in%TFToKeep),]
CandidateNapsNoPfam=CandidateNapsNoPfam[which(CandidateNapsNoPfam$ID%in%NoPfamToKeep),]








#Outlier detection on everyone : 
library(EnvStats)


AllCandidates=MSmall[which(MSmall$ID%in%c(TranscriptionFactors$ID,CandidateNapsNoPfam$ID)),]
AllCandidates$Specie=AllCandidates$Organism

ThresholdperSpecie=as.data.frame(table(AllCandidates$Organism))
names(ThresholdperSpecie)=c("Specie","Threshold")
ThresholdperSpecie$NbOutliers=NA
AllCandidates$IsOutlier=0
ThresholdperSpecie$NbOutlierMedian=0
ThresholdperSpecie$Shapiro=0
for(i in 1:dim( ThresholdperSpecie)[1]){
  print(i)
  ThresholdperSpecie[i,]$Threshold=15*median(AllCandidates[which(AllCandidates$Organism==as.character(ThresholdperSpecie[i,]$Specie)),]$NormInt)
  TF.spe=AllCandidates[which(AllCandidates$Organism==as.character(ThresholdperSpecie[i,]$Specie)),]
  #testing normality
  RosTest=rosnerTest(TF.spe$NormInt,k=5,alpha = 0.0001)
  ThresholdperSpecie[i,]$NbOutliers=RosTest$n.outliers

  RosTable=RosTest$all.stats
  ValuesToFlag=RosTable[which(RosTable$Outlier==T),]$Obs.Num
  if(length(ValuesToFlag)>0){
    AllCandidates[which(AllCandidates$Organism==ThresholdperSpecie[i,]$Specie)[ValuesToFlag],]$IsOutlier=1
    
    #To see if distrib is normal without outliers : 
    ShapTest=shapiro.test(TF.spe[-which(TF.spe$ID%in%AllCandidates[which(AllCandidates$Organism==ThresholdperSpecie[i,]$Specie & AllCandidates$IsOutlier==1),]$ID),]$NormInt)
    ThresholdperSpecie[i,]$Shapiro=ShapTest$p.value
    }
  
  
  
  ThresholdperSpecie[i,]$NbOutlierMedian=dim(AllCandidates[which(AllCandidates$Organism==ThresholdperSpecie[i,]$Specie & AllCandidates$NormInt>ThresholdperSpecie[i,]$Threshold),])[1]
}



#Graphical output : 
#Check if we have a meaningfull number of TF per genomes 
#and meaninfull number of noPfam+TF: 
DNAperSpecie=read.table("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/MASS_SPEC_PFAM/PerSpeciePFAM_abundancy.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

DNAperSpecie$NumberOfTF=NA
DNAperSpecie$NumberOfCandidatesBeforeThresh=NA
CandCount=table(AllCandidates$Organism)
for(i in DNAperSpecie$Specie){
  DNAperSpecie[which(DNAperSpecie$Specie==i),]$NumberOfTF=as.numeric(TFcount[i])
  DNAperSpecie[which(DNAperSpecie$Specie==i),]$NumberOfCandidatesBeforeThresh=as.numeric(CandCount[i])
}

#We just plot the Nb of candidates per specie : 
#We use a tree to order the plot : 

ArchTree=read.tree("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/TREES/GTDB/Arch_120_pruned_to_MassSpec_dataSet.tree")
ArchTree$tip.label=gsub("_"," ",ArchTree$tip.label)
NewOrder=TranscriptionFactors$Specie[match(ArchTree$tip.label,TranscriptionFactors$Specie)]
NewOrder=c("Halobacterium sp",NewOrder[c(1:10)],"Picrophilus torridus",NewOrder[c(11:19)])
if(length(which(is.na(NewOrder)==T))>0){
  NewOrder=NewOrder[-which(is.na(NewOrder)==T)]}

TranscriptionFactors$Specie=factor(TranscriptionFactors$Specie,levels = rev(NewOrder))


source("/Users/ahocher/Dropbox/Scripts/PLOT_LIBRARY/AH_plot_library.R")
Correlationtest=cor.test(DNAperSpecie$NumberOfTF,DNAperSpecie$GeneNb,method="spearman")
TFNbPlot=ggplot(data=DNAperSpecie,aes_string(x="NumberOfTF",y="GeneNb",color="Specie",label="Specie"))+geom_point()+geom_text_repel(size=2)+theme_pubr()+theme(aspect.ratio = 1,legend.position="")+ggtitle(paste("Spearman Cor ", round(Correlationtest$estimate,2), " p-value :", round(Correlationtest$p.value,3),sep=""))+xlab("Transcription factor (#)")+ylab("Total Genes (#)")+scale_color_manual(values=col_vector)
ggsave(TFNbPlot,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/Nb_Tf_vs_GeneNb.pdf",sep=""),width = 4,height=4)



Correlationtest=cor.test(DNAperSpecie$NumberOfCandidatesBeforeThresh,DNAperSpecie$GeneNb,method="spearman")
TFNbPlot=ggplot(data=DNAperSpecie,aes_string(x="NumberOfCandidatesBeforeThresh",y="GeneNb",color="Specie",label="Specie"))+geom_point()+geom_text_repel(size=2)+theme_pubr()+theme(aspect.ratio = 1,legend.position="")+ggtitle(paste("Spearman Cor ", round(Correlationtest$estimate,2), " p-value :", round(Correlationtest$p.value,3),sep=""))+xlab("Transcription factor & NoPfam (#)")+ylab("Total Genes (#)")+scale_color_manual(values=col_vector)
ggsave(TFNbPlot,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/Nb_Tf_and_no_Pfam_vs_GeneNb.pdf",sep=""),width = 4,height=4)


TFDistributionPlot=ggplot(data=TranscriptionFactors,aes(x=Specie,y=NormInt))+geom_violin(color="grey30")+geom_abline(slope=0,intercept = log10(median(TranscriptionFactors$NormInt)),linetype="dashed")+geom_jitter(color="grey60",size=0.75,shape=20,width = 0.2)+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% proteome (log scale)")+ggtitle("Transcription factor abundancy distribution")+xlab("")+scale_y_log10()
ggsave(TFDistributionPlot,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/TF_NormInt_distribution_perSpecie.pdf",sep=""),width = 6,height=4)



#Merging TF and putative TF/NAPs

AllCandidates$Specie=factor(AllCandidates$Specie,levels = rev(NewOrder))

CandidatesDistributionPlot=ggplot(data=AllCandidates,aes(x=Specie,y=NormInt,col=as.factor(IsOutlier)))+geom_violin(color="grey30")+geom_point(aes(shape=Type),position = position_jitter(width = 0.2),size=0.4)+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% proteome (sqrt scale)")+ggtitle("Transcription factor and NoPfam  abundancy distribution")+xlab("")+scale_y_continuous(trans="sqrt",breaks = c(0.0001,0.01,0.1,1))+scale_color_manual(values=c("grey60","red"))+scale_shape_manual(values=c(19,1))
ggsave(CandidatesDistributionPlot,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/Candidates_distribution_perSpecie_sqrt.pdf",sep=""),width = 6,height=4)


CandidatesDistributionPlot=ggplot(data=AllCandidates,aes(x=Specie,y=NormInt,col=as.factor(IsOutlier)))+geom_violin(color="grey30")+geom_point(aes(shape=Type),position = position_jitter(width = 0.2),size=0.4)+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% proteome (log2 scale)")+ggtitle("Transcription factor and NoPfam abundancy distribution")+xlab("")+scale_y_continuous(trans="log2",breaks = c(0.0001,0.01,0.1,1))+scale_color_manual(values=c("grey60","red"))+scale_shape_manual(values=c(19,1))
ggsave(CandidatesDistributionPlot,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/Candidates_distribution_perSpecie_log2.pdf",sep=""),width = 6,height=4)


CandidatesDistributionPlot=ggplot(data=AllCandidates,aes(x=Specie,y=NormInt,col=as.factor(IsOutlier)))+geom_violin(color="grey30")+geom_point(aes(shape=Type),position = position_jitter(width = 0.2),size=0.4)+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% proteome (log2 scale)")+ggtitle("Transcription factor and NoPfam abundancy distribution")+xlab("")+scale_y_continuous(breaks = c(0.01,0.1,1))+scale_color_manual(values=c("grey60","red"))+scale_shape_manual(values=c(19,1))
ggsave(CandidatesDistributionPlot,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/Candidates_distribution_perSpecie_linear_axis.pdf",sep=""),width = 6,height=4)

#Same for no pfam : 
#
NoPfamCandidatesDistrib=ggplot(data=CandidateNapsNoPfam,aes(x=Organism,y=NormInt))+geom_violin(color="grey30")+geom_jitter(color="grey60",size=0.75,shape=20,width = 0.2)+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% proteome (log scale)")+ggtitle("Transcription factor abundancy distribution")+xlab("")+scale_y_log10()

ggsave(NoPfamCandidatesDistrib,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/NoPfam_Subselection_distribution_perSpecie_log2.pdf",sep=""),width = 6,height=4)



#Step III : 

#Removing TF which are in an Operon : 

CandidateNaps=AllCandidates[which(AllCandidates$IsOutlier==1 & AllCandidates$ID%in%NotInOperon),]



#Candidate Table :
#For existing TF : 
ExportCandidates=CandidateNaps[,c(1:5,which(names(CandidateNaps)%in%c("NormInt","Pi","Length",TF_HMMs,"Type"))),]

CleanNames=ExportCandidates[,which(names(ExportCandidates)%in%TF_HMMs)]
CleanNames=names(CleanNames[,which(colSums(CleanNames)==0)])

ExportCandidates=ExportCandidates[,-which(names(ExportCandidates)%in%CleanNames)]

write.table(ExportCandidates,file ="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/NAP_PREDICTION/NAPs_candidates_all.tab" ,row.names=F,quote=F,sep="\t")

ExportCandidates=read.table(file ="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/NAP_PREDICTION/NAPs_candidates_all.tab",header=T,sep="\t",quote="",stringsAsFactors = F)
#Clustering of candidate by alignement and tree
CandidateIDs=c(ExportCandidates$ID)
SeqNames="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/NAP_PREDICTION/NAP_candidates_protein_IDsGenbank.txt"
write.table(CandidateIDs,file=SeqNames,sep="\t",quote=F,row.names = F,col.names = F)







#GRAPHICAL OUTPUT : 
#Annotate potential NAPS : 

MassSpecArchaea$CandidateNap=0
MassSpecArchaea[which(MassSpecArchaea$ID%in%CandidateIDs),]$CandidateNap=1





CandidateCount=data.frame(Specie=names(table(MassSpecArchaea$Organism)),CountTF=0,CountNoPfam=0,TotalProportionCandidates=0,PercentNAP=NA)
for(i in as.character(CandidateCount$Specie)){
  CandidateCount[which(CandidateCount$Specie==i),]$CountTF=length(which(CandidateNaps$Organism==i & CandidateNaps$TFPfam==1))
  CandidateCount[which(CandidateCount$Specie==i),]$CountNoPfam=length(which(CandidateNaps$Organism==i & CandidateNaps$NoPfam==1))
  
  
  
  CandidateCount[which(CandidateCount$Specie==i),]$TotalProportionCandidates=sum(MassSpecArchaea[which(MassSpecArchaea$Organism==i & MassSpecArchaea$CandidateNap==1),]$NormInt)
  CandidateCount[which(CandidateCount$Specie==i),]$PercentNAP=DNAperSpecie[which(DNAperSpecie$Specie==i),]$PercentNAP
}

write.table(CandidateCount,file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/NAP_PREDICTION/PerSpecie_CandidateNapCount_and_abundancy.txt",row.names=F,quote=F,sep="\t")


#Re-load the table
CandidateCount=read.table("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/NAP_PREDICTION/PerSpecie_CandidateNapCount_and_abundancy.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

library(reshape2)
CandidateCount.m=melt(CandidateCount[,which(names(CandidateCount)%in%c("Specie","TotalProportionCandidates","PercentNAP"))],id.vars = "Specie")
CandidateCount.m$Specie=as.character(CandidateCount.m$Specie)
#We just plot the Nb of candidates per specie : 
#We use a tree to order the plot : 

CandidateCount.m$Specie=factor(CandidateCount.m$Specie,levels = rev(NewOrder))

CandidateProtPlot=ggplot(data=CandidateCount.m,aes(x=Specie,y=value,fill=variable))+geom_bar(position = "stack",stat = "identity", width=0.7)+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% Proteome")+ggtitle("Nb of Candidate NAP per specie")+scale_fill_manual(values=c("royalblue","darkorange","green"))+scale_y_continuous(breaks = c(c(2,5,10)))+xlab("")


ggsave(CandidateProtPlot,filename ="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/Proportion_proteome_CandidateNaps_perSpecie.pdf",width = 6,height=4 )

#Just plotting the nb of candidates : 
CandidateCount.m=melt(CandidateCount[,which(names(CandidateCount)%in%c("Specie","CountTF","CountNoPfam"))],id.vars = "Specie")
CandidateCount.m$Specie=as.character(CandidateCount.m$Specie)
#We just plot the Nb of candidates per specie : 
#We use a tree to order the plot : 

CandidateCount.m$Specie=factor(CandidateCount.m$Specie,levels = rev(NewOrder))

CandidateProtPlot=ggplot(data=CandidateCount.m,aes(x=Specie,y=value,fill=variable))+geom_bar(position = "stack",stat = "identity", width=0.7)+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("Count")+ggtitle("Nb of Candidate NAP per specie")+scale_fill_manual(values=c("royalblue","darkorange","green"))+scale_y_continuous(breaks = c(c(1,3)))+xlab("")


ggsave(CandidateProtPlot,filename ="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/Count_proteome_CandidateNaps_perSpecie.pdf",width = 6,height=4 )











#
CandidateCount$Specie=as.character(CandidateCount$Specie)

CandidateCount$Specie=factor(CandidateCount$Specie,levels = rev(NewOrder))


CandidatePlot=ggplot()+geom_bar(data=CandidateCount[,c("TotalProportionCandidates","Specie")],aes(x=Specie,y=TotalProportionCandidates),stat = "identity", width=0.7,fill="darkblue")+geom_point(data=MassSpecArchaea[which(MassSpecArchaea$CandidateNap==1),c("Specie","NormInt")],aes(x=Specie,y=NormInt),col="grey80",size=0.8)+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of proteome")+ggtitle("NAP candidates abundancy accross species")+xlab("")

ggsave(CandidatePlot,filename ="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/Abundancy_CandidateNaps_perSpecie.pdf",width = 6,height=4 )


#MEGA PLOT : 

CandidateCount.m=melt(CandidateCount[,which(names(CandidateCount)%in%c("Specie","PercentNAP","TotalProportionCandidates"))],id.vars = "Specie")

CandidateCount.m$variable=factor(CandidateCount.m$variable,levels =c("PercentNAP","TotalProportionCandidates"))

CandidateCount.m$Specie=factor(CandidateCount.m$Specie,level=rev(NewOrder))
MassSpecArchaea$Specie=factor(MassSpecArchaea$Specie,level=rev(NewOrder))

CandidatePlot=ggplot()+geom_bar(data=CandidateCount.m,aes(x=Specie,y=value,fill=variable),stat = "identity", width=0.7)+geom_point(data=MassSpecArchaea[which(MassSpecArchaea$CandidateNap==1),c("Specie","NormInt")],aes(x=Specie,y=NormInt),col="grey80",size=0.8)+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of proteome")+ggtitle("NAP candidates abundancy accross species")+xlab("")+scale_fill_manual(values=c("#FFD28F","#65A8DB"))

ggsave(CandidatePlot,filename ="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/Abundancy_CandidateNaps_and_NAPs_perSpecie.pdf",width = 6,height=4 )




#Export sequences
PathToGenomes="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/NAP_PREDICTION/Archaea_all_species_MassSpec_proteins.fa"

setwd(dirname(PathToGenomes))

FastaNames=paste("NAP_candidate_Proteins.fa",sep="")

Command=paste("seqtk subseq ",dirname(PathToGenomes),"/",basename(PathToGenomes)," ",SeqNames," > ",FastaNames,sep="")
system(Command)





#Plots with lines for CD-hit clusters : 

CandidateCountClusters=read.table(file = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/NAP_PREDICTION/NAPs_candidates_all_CD_HIT_Annot.tab.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

MassSpecArchaea$ClusterNb=NA
for(i in 1:dim(CandidateCountClusters)[1]){
  if(is.na(CandidateCountClusters[i,]$Cluster_CDHIT)==F){
    
    MassSpecArchaea[which(MassSpecArchaea$ID==CandidateCountClusters[i,]$ID),]$ClusterNb=CandidateCountClusters[i,]$Cluster_CDHIT}}


CandidatePlot=ggplot()+geom_bar(data=CandidateCount.m,aes(x=Specie,y=value,fill=variable),stat = "identity", width=0.7)+geom_line(data=MassSpecArchaea[which(MassSpecArchaea$CandidateNap==1 & is.na(MassSpecArchaea$ClusterNb)==F),c("Specie","NormInt","ClusterNb")],aes(x=Specie,y=NormInt,group=ClusterNb,col=as.factor(ClusterNb)),size=0.5)+geom_point(data=MassSpecArchaea[which(MassSpecArchaea$CandidateNap==1),c("Specie","NormInt")],aes(x=Specie,y=NormInt),col="grey80",size=0.8)+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of proteome")+ggtitle("NAP candidates abundancy accross species")+xlab("")+scale_fill_manual(values=c("#FFD28F","#65A8DB"))+theme(legend.position="")+scale_color_manual(values=c(col_vector[5:10]))

ggsave(CandidatePlot,filename ="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/Abundancy_CandidateNaps_and_NAPs_perSpecie_with_clusters.pdf",width = 6,height=4 )




MassSpecArchaea$NAPorCandidate=NA
MassSpecArchaea[which(MassSpecArchaea$CandidateNap==1),]$NAPorCandidate="Candidate"
MassSpecArchaea[which(MassSpecArchaea$NAP==1),]$NAPorCandidate="NAP"

NAPs_vs_Candidates=ggplot(data=MassSpecArchaea[which(is.na(MassSpecArchaea$NAPorCandidate)==F),],aes(x=Specie,y=NormInt,col=as.factor(NAPorCandidate)))+geom_point(aes(col=NAPorCandidate),position = position_jitter(width = 0.4),size=1.5)+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of proteome")+ggtitle("NAP candidates abundancy accross species")+xlab("")+scale_fill_manual(values=c("#FFD28F","#65A8DB"))+theme(legend.position="")+scale_color_manual(values=c(col_vector[5:10]))
ggsave(NAPs_vs_Candidates,filename ="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/NAPs_vs_CandidatesNAPs.pdf",width = 6,height=4 )



#Make a nice and too complicated heatmap for results : 
library(reshape2)
CandidateCountClusters=CandidateCountClusters[order(CandidateCountClusters$Type,CandidateCountClusters$Organism),]

CandidateCountClusters$Organism[which(CandidateCountClusters$Organism=="Halobacterium sp")]="Halobacterium salinarium"
CandidateCountClusters$Organism[which(CandidateCountClusters$Organism=="Sulfolobus acidocaldicarius")]="Sulfolobus acidocaldarius"

CandidateCountClusters$Organism=gsub("\\^"," ",CandidateCountClusters$Organism)

CandidateCountClusters$SpecieShort=unlist(lapply(CandidateCountClusters$Organism,function(x) paste(strsplit(strsplit(x,'')[[1]][1],' ')[[1]][1],strsplit(x," ")[[1]][2],sep=".")))


CandidateCountClusters$SpecieShortName=paste(CandidateCountClusters$SpecieShort," (",CandidateCountClusters$ID,")",sep="")

IDorder=CandidateCountClusters$ID

Heat1=melt(CandidateCountClusters[,c(6:7,9:21,33)],id.vars=c("SpecieShortName"))
Heat1$SpecieShortName=factor(Heat1$SpecieShortName,levels = CandidateCountClusters$SpecieShortName)
H1=ggplot(data=Heat1)+ geom_tile(aes(fill = as.factor(value),x=variable,y=SpecieShortName),)+scale_fill_manual(drop=F,values=c("white","black"))+theme(axis.text.x=element_text(angle=90,hjust=1),legend.position = "")+coord_fixed(ratio = dim(CandidateCountClusters)[1]/dim(CandidateCountClusters)[2])+xlab("PFAM domain")+ylab("")


#Iso electric point : 
library(reshape2)
Heat2=melt(CandidateCountClusters[,c(1,26,24,4)],id.vars=c("ID","Type","Organism"))

Heat2=melt(CandidateCountClusters[,c(26,33)],id.vars=c("SpecieShortName"))
Heat2$SpecieShortName=factor(Heat2$SpecieShortName,levels = CandidateCountClusters$SpecieShortName)

H2=ggplot(data=Heat2)+ geom_tile(aes(fill = value,x=variable,y=SpecieShortName),)+scale_fill_viridis_c()+theme(axis.text.x=element_text(angle=90,hjust=1),legend.position = "")+coord_fixed(ratio = dim(CandidateCountClusters)[1]/dim(CandidateCountClusters)[2])+xlab("PFAM domain")+ylab("")


#Protein length : 

Heat3=melt(CandidateCountClusters[,c(25,33)],id.vars=c("SpecieShortName"))
Heat3$SpecieShortName=factor(Heat3$SpecieShortName,levels = CandidateCountClusters$SpecieShortName)

H3=ggplot(data=Heat3)+ geom_bar(aes(y = value,x=SpecieShortName),stat = "identity",fill="black",width = 0.7)+theme(axis.text.x=element_text(angle=90,hjust=1),legend.position = "right")+ylab("Length (a.a)")+coord_flip()+theme_pubr()+theme(aspect.ratio = 4)+xlab("")


#Norm Int

Heat4=melt(CandidateCountClusters[,c(23,33)],id.vars=c("SpecieShortName"))
Heat4$SpecieShortName=factor(Heat3$SpecieShortName,levels = CandidateCountClusters$SpecieShortName)

H4=ggplot(data=Heat4)+ geom_bar(aes(y = value,x=SpecieShortName),stat = "identity",fill="black",width = 0.7)+theme(axis.text.x=element_text(angle=90,hjust=1),legend.position = "right")+ylab("% of proteome")+coord_flip()+theme_pubr()+theme(aspect.ratio = 4)+xlab("")



ggsave(plot = H1,file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/SumUpFig/PfamHeatmap.pdf",height=4)

ggsave(plot = H2,file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/SumUpFig/PiHeatmap.pdf",height=4)
H2b=H2+theme(legend.position="top")
ggsave(plot = H2b,file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/SumUpFig/PiHeatmapwlegend.pdf",height=4)

ggsave(plot = H3,file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/SumUpFig/LengthBar.pdf",height=4)

ggsave(plot = H4+scale_y_continuous(position = "right"),file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/SumUpFig/NorMint.pdf",height=4)




#only on methanomassilicoccales : 

Heat1=melt(CandidateCountClusters[which(CandidateCountClusters$Organism%in%c("Methanomassiliicoccus luminyensis","Methanomethylophilus alvus")),c(6:7,9:21,33)],id.vars=c("SpecieShortName"))
Heat1$SpecieShortName=factor(Heat1$SpecieShortName,levels = CandidateCountClusters$SpecieShortName)
H1=ggplot(data=Heat1)+ geom_tile(aes(fill = as.factor(value),x=variable,y=SpecieShortName),)+scale_fill_manual(drop=F,values=c("white","black"))+theme(axis.text.x=element_text(angle=90,hjust=1),legend.position = "")+coord_fixed(ratio = dim(CandidateCountClusters)[1]/dim(CandidateCountClusters)[2])+xlab("PFAM domain")+ylab("")



Heat2=melt(CandidateCountClusters[which(CandidateCountClusters$Organism%in%c("Methanomassiliicoccus luminyensis","Methanomethylophilus alvus")),c(1,26,24,4)],id.vars=c("ID","Type","Organism"))

Heat2=melt(CandidateCountClusters[which(CandidateCountClusters$Organism%in%c("Methanomassiliicoccus luminyensis","Methanomethylophilus alvus")),c(26,33)],id.vars=c("SpecieShortName"))
Heat2$SpecieShortName=factor(Heat2$SpecieShortName,levels = CandidateCountClusters$SpecieShortName)

H2=ggplot(data=Heat2)+ geom_tile(aes(fill = value,x=variable,y=SpecieShortName),)+scale_fill_viridis_c()+theme(axis.text.x=element_text(angle=90,hjust=1),legend.position = "")+coord_fixed(ratio = dim(CandidateCountClusters)[1]/dim(CandidateCountClusters)[2])+xlab("PFAM domain")+ylab("")


#Protein length : 

Heat3=melt(CandidateCountClusters[which(CandidateCountClusters$Organism%in%c("Methanomassiliicoccus luminyensis","Methanomethylophilus alvus")),c(25,33)],id.vars=c("SpecieShortName"))
Heat3$SpecieShortName=factor(Heat3$SpecieShortName,levels = CandidateCountClusters$SpecieShortName)

H3=ggplot(data=Heat3)+ geom_bar(aes(y = value,x=SpecieShortName),stat = "identity",fill="black",width = 0.7)+theme(axis.text.x=element_text(angle=90,hjust=1),legend.position = "right")+ylab("Length (a.a)")+coord_flip()+theme_pubr()+theme(aspect.ratio = 4)+xlab("")


#Norm Int

Heat4=melt(CandidateCountClusters[which(CandidateCountClusters$Organism%in%c("Methanomassiliicoccus luminyensis","Methanomethylophilus alvus")),c(23,33)],id.vars=c("SpecieShortName"))
Heat4$SpecieShortName=factor(Heat3$SpecieShortName,levels = CandidateCountClusters$SpecieShortName)

H4=ggplot(data=Heat4)+ geom_bar(aes(y = value,x=SpecieShortName),stat = "identity",fill="black",width = 0.7)+theme(axis.text.x=element_text(angle=90,hjust=1),legend.position = "right")+ylab("% of proteome")+coord_flip()+theme_pubr()+theme(aspect.ratio = 4)+xlab("")



ggsave(plot = H1,file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/SumUpFig/PfamHeatmap_massilicoccales.pdf",height=4)

ggsave(plot = H2,file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/SumUpFig/PiHeatmap_massilicoccales.pdf",height=4)
H2b=H2+theme(legend.position="top")
ggsave(plot = H2b,file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/SumUpFig/PiHeatmapwlegend_massilicoccales.pdf",height=4)

ggsave(plot = H3,file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/SumUpFig/LengthBar_massilicoccales.pdf",height=4)

ggsave(plot = H4+scale_y_continuous(position = "right"),file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/SumUpFig/NorMint_massilicoccales.pdf",height=4)


#Only not methanomassilioccales : 




Heat1=melt(CandidateCountClusters[which(!(CandidateCountClusters$Organism%in%c("Methanomassiliicoccus luminyensis","Methanomethylophilus alvus"))),c(6:7,9:21,33)],id.vars=c("SpecieShortName"))
Heat1$SpecieShortName=factor(Heat1$SpecieShortName,levels = CandidateCountClusters$SpecieShortName)
H1=ggplot(data=Heat1)+ geom_tile(aes(fill = as.factor(value),x=variable,y=SpecieShortName),)+scale_fill_manual(drop=F,values=c("white","black"))+theme(axis.text.x=element_text(angle=90,hjust=1),legend.position = "")+coord_fixed(ratio = dim(CandidateCountClusters)[1]/dim(CandidateCountClusters)[2])+xlab("PFAM domain")+ylab("")



Heat2=melt(CandidateCountClusters[which(!(CandidateCountClusters$Organism%in%c("Methanomassiliicoccus luminyensis","Methanomethylophilus alvus"))),c(1,26,24,4)],id.vars=c("ID","Type","Organism"))

Heat2=melt(CandidateCountClusters[which(!(CandidateCountClusters$Organism%in%c("Methanomassiliicoccus luminyensis","Methanomethylophilus alvus"))),c(26,33)],id.vars=c("SpecieShortName"))
Heat2$SpecieShortName=factor(Heat2$SpecieShortName,levels = CandidateCountClusters$SpecieShortName)

H2=ggplot(data=Heat2)+ geom_tile(aes(fill = value,x=variable,y=SpecieShortName),)+scale_fill_viridis_c()+theme(axis.text.x=element_text(angle=90,hjust=1),legend.position = "")+coord_fixed(ratio = dim(CandidateCountClusters)[1]/dim(CandidateCountClusters)[2])+xlab("PFAM domain")+ylab("")


#Protein length : 

Heat3=melt(CandidateCountClusters[which(!(CandidateCountClusters$Organism%in%c("Methanomassiliicoccus luminyensis","Methanomethylophilus alvus"))),c(25,33)],id.vars=c("SpecieShortName"))
Heat3$SpecieShortName=factor(Heat3$SpecieShortName,levels = CandidateCountClusters$SpecieShortName)

H3=ggplot(data=Heat3)+ geom_bar(aes(y = value,x=SpecieShortName),stat = "identity",fill="black",width = 0.7)+theme(axis.text.x=element_text(angle=90,hjust=1),legend.position = "right")+ylab("Length (a.a)")+coord_flip()+theme_pubr()+theme(aspect.ratio = 4)+xlab("")


#Norm Int

Heat4=melt(CandidateCountClusters[which(!(CandidateCountClusters$Organism%in%c("Methanomassiliicoccus luminyensis","Methanomethylophilus alvus"))),c(23,33)],id.vars=c("SpecieShortName"))
Heat4$SpecieShortName=factor(Heat3$SpecieShortName,levels = CandidateCountClusters$SpecieShortName)

H4=ggplot(data=Heat4)+ geom_bar(aes(y = value,x=SpecieShortName),stat = "identity",fill="black",width = 0.7)+theme(axis.text.x=element_text(angle=90,hjust=1),legend.position = "right")+ylab("% of proteome")+coord_flip()+theme_pubr()+theme(aspect.ratio = 4)+xlab("")



ggsave(plot = H1,file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/SumUpFig/PfamHeatmap_NOTmassilicoccales.pdf",height=4)

ggsave(plot = H2,file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/SumUpFig/PiHeatmap_NOTmassilicoccales.pdf",height=4)
H2b=H2+theme(legend.position="top")
ggsave(plot = H2b,file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/SumUpFig/PiHeatmapwlegend_NOTmassilicoccales.pdf",height=4)

ggsave(plot = H3,file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/SumUpFig/LengthBar_NOTmassilicoccales.pdf",height=4)

ggsave(plot = H4+scale_y_continuous(position = "right"),file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/SumUpFig/NorMint_NOTmassilicoccales.pdf",height=4)































NAPsHMM=c("Cren7","X7kD_DNA_binding","CBFD_NFYB_HMF","Alba","Bac_DNA_binding","MC1","CC1")

MassSpecArchaea$NAPorCandidate=NA
MassSpecArchaea[which(MassSpecArchaea$CandidateNap==1),]$NAPorCandidate="Candidate"
for(j in NAPsHMM){
MassSpecArchaea[which(MassSpecArchaea[,j]==1),]$NAPorCandidate=j}


NAPs_vs_CandidatesDetailed=addSmallLegend(ggplot(data=MassSpecArchaea[which(is.na(MassSpecArchaea$NAPorCandidate)==F),],aes(x=Specie,y=NormInt,col=as.factor(NAPorCandidate)))+geom_point(aes(col=NAPorCandidate),position = position_jitter(width = 0.2),size=0.8)+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6),legend.position = "right")+ylab("% of proteome")+ggtitle("NAP candidates abundancy accross species")+xlab("")+scale_fill_manual(values=c("#FFD28F","#65A8DB"))+scale_color_manual(values=c(col_vector[10:20])))
ggsave(NAPs_vs_CandidatesDetailed,filename ="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/NAPs_vs_CandidatesNAPsDetailed.pdf",width = 10,height=4 )

MassSpecArchaea$NAPorCandidate=factor(MassSpecArchaea$NAPorCandidate,levels=rev(c("Alba","CBFD_NFYB_HMF","Bac_DNA_binding","CC1","X7kD_DNA_binding","Cren7","MC1","Candidate")))

NAPs_vs_CandidatesDetailed=ggplot(data=MassSpecArchaea[which(is.na(MassSpecArchaea$NAPorCandidate)==F),],aes(x=Specie,y=NormInt,fill=as.factor(NAPorCandidate)))+geom_bar(stat = "identity",width = 0.5)+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6),legend.position = "right")+ylab("% of proteome")+ggtitle("NAP candidates abundancy accross species")+xlab("")+scale_fill_manual(values=c(col_vector[30:40]))
ggsave(NAPs_vs_CandidatesDetailed,filename ="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/NAPs_vs_CandidatesNAPsDetailed_StackedBar.pdf",width = 10,height=4 )



#Without candidates : 
MassSpecArchaea$NAPorCandidate=factor(MassSpecArchaea$NAPorCandidate,levels=rev(c("Alba","CBFD_NFYB_HMF","Bac_DNA_binding","CC1","X7kD_DNA_binding","Cren7","MC1","Candidate")))

NAPs_Detailed=ggplot(data=MassSpecArchaea[which( !(MassSpecArchaea$ID%in%CandidateNaps$ID)& is.na(MassSpecArchaea$NAPorCandidate)==F),],aes(x=Specie,y=NormInt,fill=as.factor(NAPorCandidate)))+geom_bar(stat = "identity",width = 0.5)+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6),legend.position = "right")+ylab("% of proteome")+ggtitle("NAPs abundancy accross species")+xlab("")+scale_fill_manual(values=c(col_vector[28:36]))
ggsave(NAPs_Detailed,filename ="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/NAPs_Detailed_StackedBar.pdf",width = 10,height=4 )




#Exporting each PFAM hit and candidates abundancy
for(D in names(ExportCandidates)[6:20]){
TFPLOT=ggplot(data=MassSpecArchaea[which(MassSpecArchaea[,D]==1),],aes(x=Specie,y=NormInt,col=as.factor(CandidateNap)))+geom_point(position=position_jitter(width  = 0.2))+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6),legend.position = "")+ylab("% of proteome")+ggtitle(paste(D,", NAP candidates vs others ",sep=""))+xlab("")+scale_fill_manual(values=c(col_vector[30:40]))+scale_color_manual(values=c("grey80","red"))
ggsave(plot = TFPLOT,filename = paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/Comparison_per_TF/",D,",NAPcandidates_vs_others.pdf",sep="" ))
}

#Just to check Arsr : 
for(D in "ArsR"){
  TFPLOT=ggplot(data=MassSpecArchaea[which(MassSpecArchaea[,D]==1),],aes(x=Specie,y=NormInt,col=as.factor(CandidateNap)))+geom_point(position=position_jitter(width  = 0.2))+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6),legend.position = "")+ylab("% of proteome")+ggtitle(paste(D,", NAP candidates vs others ",sep=""))+xlab("")+scale_fill_manual(values=c(col_vector[30:40]))+scale_color_manual(values=c("grey80","red"))
  ggsave(plot = TFPLOT,filename = paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/Comparison_per_TF/",D,",NAPcandidates_vs_others.pdf",sep="" ))
}
#Conclusion : not very highly expressed. 



#Building a table to host results relative to transcription factors  : 
#
#
TranscriptionFactors$IsOutlier=0
TranscriptionFactors[which(TranscriptionFactors$ID%in%CandidateNaps$ID),]$IsOutlier=1

TFperSpecie=as.data.frame(table(TranscriptionFactors$Organism))
names(TFperSpecie)=c("Specie","NbDetectedProt")
TFperSpecie$NbTFmeasured=NA
TFperSpecie$TotalTF=NA
TFperSpecie$TotalTF_no_outliers=NA
TFperSpecie$MedianTF=NA
TFperSpecie$VarTF=NA

for(i in 1:dim(TFperSpecie)[1]){
  
TFperSpecie[i,]$NbTFmeasured=length(which(TranscriptionFactors$Organism==TFperSpecie[i,]$Specie))
TFperSpecie[i,]$TotalTF=sum(TranscriptionFactors[which(TranscriptionFactors$Organism==TFperSpecie[i,]$Specie),]$NormInt)
TFperSpecie[i,]$TotalTF_no_outliers=sum(TranscriptionFactors[which(TranscriptionFactors$Organism==TFperSpecie[i,]$Specie & TranscriptionFactors$IsOutlier==0),]$NormInt)
TFperSpecie[i,]$MedianTF=median(TranscriptionFactors[which(TranscriptionFactors$Organism==TFperSpecie[i,]$Specie ),]$NormInt)
TFperSpecie[i,]$VarTF=var(TranscriptionFactors[which(TranscriptionFactors$Organism==TFperSpecie[i,]$Specie ),]$NormInt)
}

write.table(TFperSpecie,file = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/MASS_SPEC_PFAM/PerSpecieTF_abundancy.txt",row.names = F,sep="\t",quote=F)








MassSpecArchaea$Type=""
MassSpecArchaea[which(MassSpecArchaea$CBFD_NFYB_HMF==1),]$Type="Histone"

HistoneTF=ggplot(data=MassSpecArchaea[which((MassSpecArchaea$ID%in%TranscriptionFactors$ID | MassSpecArchaea$CBFD_NFYB_HMF==1) & MassSpecArchaea$Specie%in%c("Methanosarcina barkeri","Natrialba magadii","Haloferax volcanii","Halobacterium sp","Nitrosopumilus^maritimus^SCM1")),])+geom_jitter(aes(y=NormInt,color=Type,x=Specie))+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of proteome")+xlab("")+scale_color_manual(values=c("grey80","#FFD28F","#65A8DB"))+scale_y_sqrt()


ggsave(addSmallLegend(HistoneTF),file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/Comparison_HistoneTF_in_Subspecies_NAPsnocandidates.pdf",width=4,height=4)




MassSpecArchaea$Candidate=0
MassSpecArchaea[which(MassSpecArchaea$ID%in%CandidateCountClusters$ID),]$Candidate=1

MassSpecArchaea$Type=""
MassSpecArchaea[which(MassSpecArchaea$ID%in%CandidateCountClusters$ID),]$Type="Candidate"
MassSpecArchaea[which(MassSpecArchaea$CBFD_NFYB_HMF==1),]$Type="Histone"

HistoneTF=ggplot(data=MassSpecArchaea[which((MassSpecArchaea$ID%in%TranscriptionFactors$ID | MassSpecArchaea$CBFD_NFYB_HMF==1) & MassSpecArchaea$Specie%in%c("Methanosarcina barkeri","Natrialba magadii","Haloferax volcanii","Halobacterium sp","Nitrosopumilus^maritimus^SCM1")),])+geom_jitter(aes(y=NormInt,color=Type,x=Specie))+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of proteome")+xlab("")+scale_color_manual(values=c("grey80","#FFD28F","#65A8DB"))+scale_y_sqrt()


ggsave(addSmallLegend(HistoneTF),file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/CANDIDATENAPS/Comparison_HistoneTF_in_Subspecies_NAPs.pdf",width=4,height=4)







