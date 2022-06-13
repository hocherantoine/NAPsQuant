#####################################################
## Project: Archaeal chromatin                                          #
## Analyse T. volcanium nucleoid enrichment data
## Date: january 2022                                                               #
## Author: Antoine Hocher                                                                  #
###################################################################
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
source("/Users/ahocher/Dropbox/Scripts/PLOT_LIBRARY/AH_plot_library.R")
#Creating colors for later : 
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))


#Loading all data :
#Load annotated mass spec data : 
MassSpecArchaea=read.table(file="Hocher_2022_MassSpecMeasurementsAccrossSpecies_wPFAM_Annotations.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
MassSpecArchaea$Specie=MassSpecArchaea$Organism


TvolcWCE=MassSpecArchaea[which(MassSpecArchaea$Specie=="Thermoplasma volcanium"),]



#Load Mass Spectrometry data : 
##NOTE GITHUB : data available from PRIDE mass spec database
Tvolcanium=read.delim(file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/Revisions/b023p058_Tvolc_proteinGroups_table.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
#Change name of ID column to merge later on:
names(Tvolcanium)[1]="ID"

# get LFQ column numbers
LFQ_columns <- grep("LFQ.", colnames(Tvolcanium))


#Computing Replicates correlations and average : 
PearsCorMat=cor(Tvolcanium[,LFQ_columns])
SpearCorMat=cor(Tvolcanium[,LFQ_columns],method="spearman")

library(corrplot)
corrplot(PearsCorMat)
corrplot.mixed(SpearCorMat)


#merging Public data with nucleoid enrichment data
Tvolc=merge(Tvolcanium,TvolcWCE,by="ID",all.x=T)

#Differential abundancy calculation
data=Tvolc

#Averaging technical duplicates
data$LFQ.intensity.Tvolc_nuc_1=(data$LFQ.intensity.Tvolc_Nuc_1_TR01+data$LFQ.intensity.Tvolc_Nuc_1_TR02)/2
data$LFQ.intensity.Tvolc_nuc_2=(data$LFQ.intensity.Tvolc_Nuc_2_TR01+data$LFQ.intensity.Tvolc_Nuc_2_TR02)/2

data$LFQ.intensity.Tvolc_top_1=(data$LFQ.intensity.Tvolc_Top_1_TR01+data$LFQ.intensity.Tvolc_Top_1_TR02)/2
data$LFQ.intensity.Tvolc_top_2=(data$LFQ.intensity.Tvolc_Top_2_TR01+data$LFQ.intensity.Tvolc_Top_2_TR02)/2


# get LFQ column numbers
LFQ_columns <- grep("LFQ.", colnames(data))

#Removed the pooled measurements and technical rep :
LFQ_columns_toRemove=LFQ_columns[1:12]
data=data[,-LFQ_columns_toRemove]

# get LFQ column numbers again (they changed after deletion)
LFQ_columns <- grep("LFQ.", colnames(data))


#This is required for dep
data$Gene.names=data$ID
data$Protein.IDs=data$ID



ShortNames=gsub("LFQ.intensity.","",names(data)[LFQ_columns])

#Exp.design : 
experimental_design=as.data.frame(matrix(nrow=length(LFQ_columns),ncol=3))
names(experimental_design)=c("label","condition","replicate")
experimental_design$label=names(data)[LFQ_columns]
experimental_design$condition=unlist(lapply(ShortNames,function(x) strsplit(x,"_")[[1]][2]))
experimental_design$replicate=unlist(lapply(ShortNames,function(x) strsplit(x,"_")[[1]][3]))



data_results <- LFQ(data, experimental_design,fun ="man" ,type = "all", alpha = 0.05, lfc = 0)

AA=data_results$results

#Computing another p-value correction, as the initial algorythm deems no gene sig, which is not correct and due to multiple testing pval correction :
AA$nuc_vs_top_p.val.FDR=p.adjust(AA$nuc_vs_top_p.val,method = "fdr")

#Annotating dna binding proteins
AA$PFAM_DNA_binding=0
AA[which(AA$ID%in%data[which(data$PFAM_DNA_binding==1),]$Protein.IDs),]$PFAM_DNA_binding=1

write.table(AA,file="b023p058_Tvolcanium_proteinGroups_DEP_output.txt",sep="\t",quote=F,row.names=F)



#Part II  : comparing abundancy and enrichment : 
#Merging Fold change and enrichment : 
#
Tvolc.temp=merge(AA,data,by="ID",all=T)

table(is.na(Tvolc.temp$ID))

#Adding localisation data obtained from psort db (should be re-done by user, using genbank IDs for proteins to enable compatibility:
#not essential for reproducibility, was used as control internaly 
LocaT=read.table("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PSORTDB/PSORTbTvolcanium.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
names(LocaT)[1]="ID"

#Formating ID properly
LocaT$ID=unlist(lapply(LocaT$ID,function(x) strsplit(x," ")[[1]][1]))

#Adding location data
Tvolc.temp2=merge(Tvolc.temp,LocaT,by="ID",all.x=T)

#loading operon mapper functional information 
#NOTE GITHUB : has to be redone by user, not essential at all, in fact it was mostly used for the nice protein annotation provided by operon mapper
OperonMapper=read.delim("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/OPERO_MAPPER/ARCHAEA/GCA_000011185.1_list_of_operons_27106",sep="\t",quote="",stringsAsFactors = F)
OperonMapper.annot=read.delim("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/OPERO_MAPPER/ARCHAEA/Processed/GCA_000011185.1_list_of_operons_27106_reformat_ah.txt",sep="\t",quote="",stringsAsFactors = F)
names(OperonMapper.annot)[which(names(OperonMapper.annot)=="ProtName")]="ID"
OperonMapper.f=merge(OperonMapper,OperonMapper.annot,by="IdGene")

OperonMapper.f$IdGene

OperonMapper.f$Category=unlist(lapply(OperonMapper.f$Function,function(x) strsplit(x,"\\]")[[1]][1]))
OperonMapper.f$Category=gsub("\\[","",OperonMapper.f$Category)

#####
Tvolcanium.final=as.data.frame(merge(Tvolc.temp2,OperonMapper.f,by="ID",all.x=T))




#Computing the Average LFQ intensity of the nucleoid fraction:
Tvolcanium.final$AvgNuc=(Tvolcanium.final$LFQ.intensity.Tvolc_nuc_1+Tvolcanium.final$LFQ.intensity.Tvolc_nuc_2)/2
Tvolcanium.final$AvgTop=(Tvolcanium.final$LFQ.intensity.Tvolc_top_1+Tvolcanium.final$LFQ.intensity.Tvolc_top_2)/2


EnrichmentPlot=ggplot(data=Tvolcanium.final)+geom_point(aes(y=nuc_vs_top_ratio,x=sqrt((AvgNuc)),color=as.factor(PFAM_DNA_binding.x)),size=0.5)+geom_text_repel(data=Tvolcanium.final[which(Tvolcanium.final$name%in%c("BAB59743.1","BAB59768.1")),],aes(y=nuc_vs_top_ratio,x=sqrt(AvgNuc),label=ID),color="orange",size=2)+geom_text_repel(data=Tvolcanium.final[which(Tvolcanium.final$name%in%c("BAB59303.1","BAB60227.1")),],aes(y=nuc_vs_top_ratio,x=sqrt(AvgNuc),label=ID),color="red",size=2)+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")+scale_color_manual(values=c("grey80","blue"))+ylab("Nucleoid fraction LFQ intensity (sqrt)")

ggsave(plot = EnrichmentPlot,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/REVISIONS/Tvolcanium_nucleoidEnrichment.pdf")



EnrichmentPlot=ggplot(data=Tvolcanium.final)+geom_point(aes(y=nuc_vs_top_ratio,x=NormInt,color=as.factor(PFAM_DNA_binding.x)),size=0.5)+geom_text_repel(data=Tvolcanium.final[which(Tvolcanium.final$name%in%c("BAB59743.1","BAB59768.1")),],aes(y=nuc_vs_top_ratio,x=NormInt,label=ID),color="orange",size=2)+geom_text_repel(data=Tvolcanium.final[which(Tvolcanium.final$name%in%c("BAB59303.1","BAB60227.1")),],aes(y=nuc_vs_top_ratio,x=NormInt,label=ID),color="red",size=2)+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")+scale_color_manual(values=c("grey80","blue"))+xlab("WCE %")

ggsave(plot = EnrichmentPlot,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/REVISIONS/Tvolcanium_nucleoidEnrichment_WCE_NormInt.pdf")


EnrichmentPlot=ggplot(data=Tvolcanium.final)+geom_point(aes(y=nuc_vs_top_ratio,x=sqrt(NormInt),color=as.factor(PFAM_DNA_binding.x)),size=0.5)+geom_point(data=Tvolcanium.final[which(Tvolcanium.final$Category=="K"),],aes(y=nuc_vs_top_ratio,x=sqrt(NormInt)),color="orange",size=0.5)+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")+scale_color_manual(values=c("grey80","blue"))+xlab("WCE % (sqrt)")

ggsave(plot = EnrichmentPlot,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/REVISIONS/Tvolcanium_nucleoidEnrichment_WCE_NormInt_K_orange.pdf")



EnrichmentPlot=ggplot(data=Tvolcanium.final)+geom_point(aes(y=nuc_vs_top_ratio,x=sqrt(NormInt),color=as.factor(PFAM_DNA_binding.x)),size=0.5)+geom_point(data=Tvolcanium.final[which(Tvolcanium.final$Localization%in%c("CytoplasmicMembrane","Extracellular","CellWall")),],aes(y=nuc_vs_top_ratio,x=sqrt(NormInt)),color="red",size=0.5)+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")+scale_color_manual(values=c("grey80","blue"))+ylab("WCE % (sqrt)")

ggsave(plot = EnrichmentPlot,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/REVISIONS/Tvolcanium_nucleoidEnrichment_WCE_NormInt_Membrane_etc_red.pdf")



#Plotting enrichment in function of localisation : 

Tvolcanium.final$Localization=factor(Tvolcanium.final$Localization,levels = rev(c("Cellwall","CytoplasmicMembrane","Extracellular","Unknown","Cytoplasmic")))

LocalisationPlot=ggplot(data=Tvolcanium.final)+geom_violin(aes(y=nuc_vs_top_ratio,x=as.factor(Localization),group=as.factor(Localization),fill=as.factor(Localization)),size=0,width=1)+geom_point(data=Tvolcanium.final[which(Tvolcanium.final$name=="B0SMR4"),],aes(y=nuc_vs_top_ratio,x=as.factor(Localization),group=as.factor(Localization)),color="black",size=2)+geom_point(data=Tvolcanium.final[which(Tvolcanium.final$name=="B0STP0"),],aes(y=nuc_vs_top_ratio,x=as.factor(Localization),group=as.factor(Localization)),color="black",size=2) +theme_pubclean()+theme(aspect.ratio=1,legend.position = "")+scale_fill_manual(values=col_vector[4:10])+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10))+xlab("")+coord_flip()

ggsave(plot = LocalisationPlot,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/REVISIONS/Tvolcanium_nucleoidEnrichment_Localisation.pdf")




Ktrans=Tvolcanium.final[which(Tvolcanium.final$Category=="K"),which(names(Tvolcanium.final)%in%c("ID","Function","PFAM_DNA_binding.x","AvgNuc","nuc_vs_top_ratio","Mol..weight..kDa.","NormInt"))]
Ktrans=Ktrans[order(-Ktrans$AvgNuc),]


#Annotating candidates
Tvolcanium.final$Candidate=0
Tvolcanium.final[which(Tvolcanium.final$ID%in%c("BAB59743.1","BAB59768.1")),]$Candidate=1

#Export data : 
write.table(Tvolcanium.final[,which(names(Tvolcanium.final)%in%c("ID","Function","PFAM_DNA_binding.x","AvgNuc","AvgTop","nuc_vs_top_ratio","Mol..weight..kDa.","NormInt","Category","NAP","Candidate"))],file = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/Revisions/ANALYSIS/b023p058_Tvolcanium_Nucleoid_enrichment_annotated.txt",sep="\t",quote=F,row.names = F)





WorkingGraph=ggplot(data=Tvolcanium.final[which(Tvolcanium.final$Category=="K"),])+geom_point(aes(y=nuc_vs_top_ratio,x=NormInt,color=as.factor(PFAM_DNA_binding.x)),size=0.5)+geom_point(data=Tvolcanium.final[which(Tvolcanium.final$name%in%c("CAJE01000015.1-144","CAJE01000001.1-33","CAJE01000005.1-32")),],aes(y=nuc_vs_top_ratio,x=NormInt),color="orange",size=2)+ggrepel::geom_text_repel(aes(y=nuc_vs_top_ratio,x=NormInt,label=Function),max.overlaps = 100,size=2)+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")+scale_color_manual(values=c("grey80","blue"))+ylab("NormInt")


ggsave(plot = WorkingGraph,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/REVISIONS/Tvolcanium_nucleoidEnrichment_K_categoryOnly_wLabels.pdf",width=10,height=10)
