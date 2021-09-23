#####################################################
## Project: diaforarchaea chromatin
## Analysis of growth temp and genome compactness
## https://doi.org/10.1038/s41586-020-2402-x
## Aim : 
## 
## see if naps are best correlated with temp or with genome compactness
## Associate proteomes to our local database
## 
## Date: feb 2021
## Author: Antoine Hocher
####################################################

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.

source("/Users/ahocher/Dropbox/Scripts/PLOT_LIBRARY/AH_plot_library.R")

#This table is fro Operon_Mapper_Results, it's computed using Operon mapper.
DNAperSpecie=read.table("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/MASS_SPEC_PFAM/PerSpeciePFAM_abundancy_and_TSSOperons.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

#This table contains G.O measurements : 
GOperSpecie=read.table("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/MASS_SPEC_PFAM/PerSpecie_GO_PFAM_abundancy.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

#And loading corresponding GO names :
Pfam2Go=read.table(file="/Users/ahocher/Dropbox/Laboratory/ArchealChromatinProteom/HMMs/PfamDb/Pfam2go.txt",skip=6,sep=";",stringsAsFactors = F)
head(Pfam2Go)
Pfam2Go$ID=unlist(lapply(Pfam2Go$V1,function(x) strsplit(x,split = " ")[[1]][2]))
Pfam2Go$Description=unlist(lapply(Pfam2Go$V1,function(x) strsplit(x,split = ">")[[1]][2]))
names(Pfam2Go)[2]="GO"
Pfam2Go$GO=gsub(" ","",Pfam2Go$GO)
Pfam2Go$GO=gsub(":",".",Pfam2Go$GO)


#Operons & purine bias: 
GenomePpties=read.table("/Users/ahocher/Dropbox/Laboratory/ArchealChromatinProteom/Nucleotide_periodicity/Analysis/Archaea_EBMC2_Operons_properties_results.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
GenomePpties=GenomePpties[,which(names(GenomePpties)%in%c("Assembly","Tss.Proportion","PurineBias"))]


CandidatePerSpecie=read.table("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/NAP_PREDICTION/PerSpecie_CandidateNapCount_and_abundancy.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
CandidatePerSpecie$PercentNAPwCandidates=CandidatePerSpecie$PercentNAP+CandidatePerSpecie$TotalProportionCandidates

CandidatePerSpecie$PercentNAPwCandidateswOperon=CandidatePerSpecie$PercentNAP+CandidatePerSpecie$TotalProportionCandidatesOperon

CandidatePerSpecie=CandidatePerSpecie[-which(names(CandidatePerSpecie)=="PercentNAP")]
write.table(CandidatePerSpecie,"/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/NAP_PREDICTION/PerSpecie_CandidateNapCount_and_abundancy_FinalNAPwCandidatevalues.txt",sep="\t",quote=F,row.names=F)
TFPerSpecie=read.table("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/MASS_SPEC_PFAM/PerSpecieTF_abundancy.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
names(TFPerSpecie)[4]="TotalTF_percentage"
#Loading archaeal genomes ids :  
Genomes=read.table("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/MASS_SPEC_Muller2020/Annotated_genomes_info_OperonsMapped_and_Phenotypes.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

Genomes=merge(Genomes, GenomePpties,by="Assembly")
Genomes=Genomes[-which(names(Genomes)=="Specie")]

#Get NCBI Taxonomy for colors : 
#NCBI Taxonomy file path is : 
source("/Users/ahocher/Dropbox/Scripts/tax_id_to_ranked_lineage.R")
TaxFilePath="/Users/ahocher/Dropbox/Laboratory/ArchealChromatinProteom/Database/ncbi_taxonomy_files/"


TaxonomyInfo=data.frame(tax_id_to_ranked_lineage(as.character(Genomes$Taxid),TaxFilePath),stringsAsFactors=T)
names(TaxonomyInfo)[which(names(TaxonomyInfo)=="tax_id")]="Taxid"
TaxonomyInfo[which(TaxonomyInfo$Taxid=="338192"),]$class="Nitrosopumilales"
Genomes.m=merge(Genomes,TaxonomyInfo,by="Taxid")



DNAperSpecie.a=merge(DNAperSpecie,Genomes.m,by=c("Assembly"))
DNAperSpecie.m1=merge(DNAperSpecie.a,CandidatePerSpecie,by="Specie")
DNAperSpecie.m=merge(DNAperSpecie.m1,TFPerSpecie,by="Specie")


#Craftin a nice specie name for later : 
#Fix kodak
DNAperSpecie.m$Specie[4]="Halobacterium salinarium"
DNAperSpecie.m$Specie[16]="Sulfolobus acidocaldarius"

DNAperSpecie.m$Specie[17]=gsub("\\^"," ",DNAperSpecie.m$Specie[17])
DNAperSpecie.m$Specie[14]=gsub("\\^"," ",DNAperSpecie.m$Specie[14])

DNAperSpecie.m$SpecieShort=unlist(lapply(DNAperSpecie.m$Specie,function(x) paste(strsplit(strsplit(x,'')[[1]][1],' ')[[1]][1],strsplit(x," ")[[1]][2],sep=". ")))



#First computing the PFAM domains that correlate best with temperature :
PFAMDomains=read.table("/Users/ahocher/Dropbox/Laboratory/ArchealChromatinProteom/HMMs/PfamDb/All_pfamID.txt",stringsAsFactors = F,quote="")$V1
PFAMDomains=make.names(PFAMDomains)
PFAMDomains=PFAMDomains[which(PFAMDomains%in%names(DNAperSpecie.m))]





#Subselect the PFAM domains for which > 75* data are populated : 
ListPFAM=apply(DNAperSpecie.m[,which(names(DNAperSpecie.m)%in%c(PFAMDomains,"PercentNAP","PercentNAPwCandidates"))],2, function(x) length(which(x==0)) )
ListPFAM=names(ListPFAM[which(ListPFAM<=10)])


CorTable=data.frame(ID=ListPFAM,Cor=NA,pval=NA)


for(i in 1:dim(CorTable)[1]){
  testres=cor.test(DNAperSpecie.m[,which(names(DNAperSpecie.m) ==  CorTable[i,]$ID)],DNAperSpecie.m$OGT_Manual,method="spearman")
  CorTable[i,]$Cor=testres$estimate
  CorTable[i,]$pval=testres$p.value}

CorTable$Adj.p=p.adjust(CorTable$pval,method="BH")
CorTable=CorTable[order(CorTable$Cor),]
head(CorTable,n=30);tail(CorTable,n=30)
write.table(CorTable,file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/MASS_SPEC_PFAM/Significant_PFAM_correlations_with_Growth_Temperature.txt",row.names = F,sep="\t",quote=F)

CorTable=read.table(file = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/MASS_SPEC_PFAM/Significant_PFAM_correlations_with_Growth_Temperature.txt",header=T,sep="\t",quote="",stringsAsFactors = F)



MinVal=-log10(min(CorTable$pval))
CorTable$logpval=-log10(CorTable$pval)
tail(CorTable)
Q=ggplot(data=CorTable,aes(x=Cor,y=logpval,color=Cor))+geom_rect(mapping=aes(xmin= -1.25, xmax= -1, ymin= 0, ymax=MinVal),col="white",fill="white")+geom_rect(mapping=aes(xmin= 1, xmax= 1.25, ymin= 0, ymax=MinVal),col="white",fill="white")+geom_rect(mapping=aes(xmin= -1, xmax= -1, ymin= 0, ymax=MinVal),col="grey",fill="grey")+geom_rect(mapping=aes(xmin= -1.25, xmax= 1.25, ymin= MinVal, ymax=MinVal),col="grey",fill="grey")+geom_point(aes(color=(Cor)))+ggrepel::geom_text_repel(data=CorTable[c(1:8,c((dim(CorTable)[1]-8):dim(CorTable)[1])),],aes(x=Cor,y=-log10(pval),label=ID),size=2.25,point.padding = 0.2,arrow = arrow(length = unit(0.01, "npc")))+xlim(-1.25,1.25)+coord_polar(theta = "x",start = 0,clip = "on",direction = -1)+scale_color_gradient2(low = "#000099",mid = "white",high = "red",midpoint = 0)+theme_bw()+theme(axis.line = element_blank())


ggsave(Q,file= "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/MASS_SPEC_PFAM/Cor_tempgraphtesttt.pdf",height=10,width = 10)


Q=ggplot(data=CorTable,aes(x=Cor,y=logpval,color=Cor))+geom_point(aes(color=(Cor)))+ggrepel::geom_text_repel(data=CorTable[c(1:8,c((dim(CorTable)[1]-8):dim(CorTable)[1])),],aes(x=Cor,y=-log10(pval),label=ID),size=2.25,point.padding = 0.2,arrow = arrow(length = unit(0.01, "npc")))+scale_color_gradient2(low = "#000099",mid = "white",high = "red",midpoint = 0)+theme_pubclean()+theme(axis.line = element_blank())+xlab("Spearman Rho")+ylab("-log10(p-value)")


ggsave(Q,file= "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/MASS_SPEC_PFAM/Cor_tempgraphtesLiner.pdf",height=4,width = 4)

CorTableSig=CorTable[which(CorTable$Adj.p<0.05),];CorTableSig=CorTableSig[order(CorTableSig$Cor),]
StripCorChart=ggplot(data=CorTable)+geom_abline(slope=0,intercept = 0)+geom_histogram(aes(x=Cor),bins = 100,fill="grey90")+geom_point(data=CorTableSig,aes(x=Cor,y=-3,color=(Cor)),shape=124,size=2)+ggrepel::geom_text_repel(data=CorTableSig[c(1:4,c((dim(CorTableSig)[1]-4):dim(CorTableSig)[1])),],aes(x=Cor,y=-3,label=ID),size=1,point.padding = 0.2,arrow = arrow(length = unit(0.01, "npc")))+scale_color_gradient2(low = "#000099",mid = "white",high = "red",midpoint = 0)+theme_pubr()+xlab("Spearman Rho")+ylab("count")+scale_y_reverse()+theme(aspect.ratio = 1/6)

ggsave(StripCorChart,file= "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/MASS_SPEC_PFAM/StripCorChart_Temperature_vs_PFAM.pdf",height=4,width = 6)


#Same for Gene ontology : 
GOperSpecie$Specie=gsub("\\^"," ",GOperSpecie$Specie)

GOperSpecie$Specie[16]="Sulfolobus acidocaldarius"
GOperSpecie$Specie[4]="Halobacterium salinarium"

GOmerge=merge(GOperSpecie,DNAperSpecie.m[,which(names(DNAperSpecie.m)%in%c("PercentNAPwCandidates","OGT_Manual","Specie"))],by="Specie")


GONAMES=names(GOperSpecie)[3:606]

#Subselect the PFAM domains for which > 75* data are populated : 
ListPFAM=apply(GOmerge[,which(names(GOmerge)%in%c(GONAMES,"PercentNAPwCandidates"))],2, function(x) length(which(x==0)) )
ListPFAM=names(ListPFAM[which(ListPFAM<=10)])


CorTable=data.frame(ID=ListPFAM,Cor=NA,pval=NA)


for(i in 1:dim(CorTable)[1]){
  testres=cor.test(GOmerge[,which(names(GOmerge) ==  CorTable[i,]$ID)],GOmerge$OGT_Manual,method="spearman")
  CorTable[i,]$Cor=testres$estimate
  CorTable[i,]$pval=testres$p.value}

CorTable$Adj.p=p.adjust(CorTable$pval,method="BH")
CorTable=CorTable[order(CorTable$Cor),]
head(CorTable,n=30);tail(CorTable,n=30)

#Merging with descrition of gene ontology
Pfam2Go=Pfam2Go[,c(2,4)]
names(Pfam2Go)[1]="ID"
Pfam2Go=Pfam2Go[which(duplicated(Pfam2Go$ID)==F),]
CorTable.m=merge(CorTable,Pfam2Go,by="ID",all.x=T)




write.table(CorTable.m,file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/MASS_SPEC_PFAM/Significant_GO_correlations_with_Growth_Temperature.txt",row.names = F,sep="\t",quote=F)

CorTable=read.delim(file = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/MASS_SPEC_PFAM/Significant_GO_correlations_with_Growth_Temperature.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

CorTable$Description=gsub("\\GO: ","",CorTable$Description)
CorTable=CorTable[order(CorTable$Cor),]
MinVal=-log10(min(CorTable$pval))
CorTable$logpval=-log10(CorTable$pval)
tail(CorTable)
CorTable[which(CorTable$ID=="PercentNAPwCandidates"),]$Description="% NAP with Candidates"
Q=ggplot(data=CorTable,aes(x=Cor,y=logpval,color=Cor))+geom_point(aes(color=(Cor)))+ggrepel::geom_text_repel(data=CorTable[c(1:8,c((dim(CorTable)[1]-8):dim(CorTable)[1])),],aes(x=Cor,y=-log10(pval),label=Description),size=2,point.padding = 0.2,arrow = arrow(length = unit(0.01, "npc")),max.overlaps = Inf)+scale_color_gradient2(low = "#000099",mid = "white",high = "red",midpoint = 0)+theme_pubclean()+theme(axis.line = element_blank())+xlab("Spearman Rho")+ylab("-log10(p-value)")


ggsave(Q,file= "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/MASS_SPEC/MASS_SPEC_PFAM/Cor_tempgraphtesLinear_GeneOntology.pdf",height=6,width = 6)






#Now plots on Genome properties :


#Graphical output : 

Correlationtest=cor.test(DNAperSpecie.m$GeneNb,DNAperSpecie.m$PercentDNA,method="spearman")
options(scipen = 2)
GNbDNA=addSmallLegend(ggplot(data=DNAperSpecie.m,aes(y=GeneNb,x=PercentDNA))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+xlab("Proteins containing DNA binding pfam (%)")+ylab("Number of genes")+scale_color_manual(values=col_vector[10:30]))

ggsave(GNbDNA,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_Percent_GeneNb_vsPFAMDNA.pdf",width=12,height = 6)



Correlationtest=cor.test(DNAperSpecie.m$PercentNAP,DNAperSpecie.m$PercentDNA,method="spearman")
options(scipen = 2)
DNA_NAP=addSmallLegend(ggplot(data=DNAperSpecie.m,aes(y=PercentNAP,x=PercentDNA))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+xlab("Proteins containing DNA binding pfam (%)")+ylab("NAP abundancy (%)")+scale_color_manual(values=col_vector[10:30]))

ggsave(DNA_NAP,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_Percent_PFAMDNA_vsPFAMNAP.pdf",width=12,height = 6)




Correlationtest=cor.test(DNAperSpecie.m$GeneNb,DNAperSpecie.m$PercentDNAnotNAP,method="spearman")
options(scipen = 2)
DNANONAP=addSmallLegend(ggplot(data=DNAperSpecie.m,aes(y=GeneNb,x=PercentDNAnotNAP))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+xlab("Proteins containing DNA binding pfam, not NAP (%)")+ylab("Gene Number")+scale_color_manual(values=col_vector[10:30]))
ggsave(DNANONAP,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_Percent_GeneNb_vsPFAMDNANONAP.pdf",width=12,height = 6)




Correlationtest=cor.test(DNAperSpecie.m$GeneNb,DNAperSpecie.m$PercentNAP,method="spearman")
options(scipen = 2)
GnbNap=addSmallLegend(ggplot(data=DNAperSpecie.m,aes(y=GeneNb,x=PercentNAP))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+xlab("NAP abundancy (%)")+ylab("Gene Number")+scale_color_manual(values=col_vector[10:30]))


ggsave(GnbNap,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_Percent_GeneNb_vsPFAMNAP.pdf",width=4,height = 4)




#Just to make sure the % of NAP is not influence by the quantity of peptide detected.
NApvsTotalIntensity=ggplot(data=DNAperSpecie[-which(DNAperSpecie$Specie=="Thermococcus^kodakarensis^KOD1"),],aes(x=PercentNAP,y=TotalProt,color=Specie))+geom_point()+ggrepel::geom_text_repel(aes(label=Specie),size=2)+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")

ggsave(NApvsTotalIntensity,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/Nap_vs_TotalInt_no_Kodak_bcause_dif_technique.pdf",width=6,height = 6)



Correlationtest=cor.test(DNAperSpecie.m$Compaction,DNAperSpecie.m$PercentNAP,method="spearman")
options(scipen = 2)
CompactionNap=addSmallLegend(ggplot(data=DNAperSpecie.m,aes(y=Compaction,x=PercentNAP))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+ylab("Compactness (Gene Nb / Size)")+xlab("NAP abundancy (% of proteome)")+scale_color_manual(values=col_vector[10:30]))

ggsave(CompactionNap,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_Compaction_geneNbSize_vsPFAMNAP.pdf",width=12,height = 6)





#Compaction Bp

Correlationtest=cor.test(DNAperSpecie.m$CompactionBp,DNAperSpecie.m$PercentNAP,method="spearman")
options(scipen = 2)
CompactionBpNap=addSmallLegend(ggplot(data=DNAperSpecie.m,aes(y=100*CompactionBp,x=PercentNAP))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+xlab("NAP abundancy (% of proteome)")+ylab("% of coding intergenes")+scale_color_manual(values=col_vector[10:30]))

ggsave(CompactionBpNap,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_Percent_CompactionBp_vsPFAMNAP.pdf",width=12,height = 6)










#Divergent Intergenes vs NAP
Correlationtest=cor.test(DNAperSpecie.m$DivergentIntergeneProportion,DNAperSpecie.m$PercentNAP,method="spearman")
options(scipen = 2)
DivergentIntergenesBpNap=addSmallLegend(ggplot(data=DNAperSpecie.m,aes(y=100*DivergentIntergeneProportion,x=PercentNAP))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+xlab("NAP abundancy (% of proteome)")+ylab("Divergent intergenes proportion (%)")+scale_color_manual(values=col_vector[10:30]))
ggsave(DivergentIntergenesBpNap,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_Proportion_Divergent_Intergenes_vsPFAMNAP.pdf",width=12,height = 6)


#Aligned Intergenes vs NAP
Correlationtest=cor.test(DNAperSpecie.m$AlignedIntergeneProportion,DNAperSpecie.m$PercentNAP,method="spearman")
options(scipen = 2)
AlignedIntergenesBpNap=addSmallLegend(ggplot(data=DNAperSpecie.m,aes(y=100*AlignedIntergeneProportion,x=PercentNAP))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+xlab("NAP abundancy (% of proteome)")+ylab("Aligned intergenes proportion (%)")+scale_color_manual(values=col_vector[10:30]))

ggsave(AlignedIntergenesBpNap,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_Proportion_Aligned_Intergenes_vsPFAMNAP.pdf",width=12,height = 6)



#Very different : convergent intergenes don't have promoters.
#Convergent Intergenes vs NAP
Correlationtest=cor.test(DNAperSpecie.m$ConvergentIntergeneProportion,DNAperSpecie.m$PercentNAP,method="spearman")
options(scipen = 2)
ConvergentIntergenesBpNap=addSmallLegend(ggplot(data=DNAperSpecie.m,aes(y=100*ConvergentIntergeneProportion,x=PercentNAP))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+xlab("NAP abundancy (% of proteome)")+ylab("Convergent intergenes proportion (%)")+scale_color_manual(values=col_vector[10:30]))


ggsave(ConvergentIntergenesBpNap,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_Proportion_Convergent_Intergenes_vsPFAMNAP.pdf",width=12,height = 6)








#Verifying if the results on convergent vs Divergent intergenes and all intergenes stands  with candidatse: 

#Very different : convergent intergenes don't have promoters.
#Convergent Intergenes vs NAP

Correlationtest=cor.test(DNAperSpecie.m$CompactionBp,DNAperSpecie.m$PercentNAPwCandidates,method="spearman")
options(scipen = 2)
CompactionBpNap=addSmallLegend(ggplot(data=DNAperSpecie.m,aes(y=100*CompactionBp,x=PercentNAPwCandidates))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+xlab("NAP abundancy (% of proteome)")+ylab("non convergent intergenes proportion (%)")+scale_color_manual(values=col_vector[10:30]))+scale_y_log10()


ggsave(CompactionBpNap,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_CompactionBp_vsPFAMNAPwcandidates.pdf",width=5,height = 5)



Correlationtest=cor.test(DNAperSpecie.m$ConvergentIntergeneProportion,DNAperSpecie.m$PercentNAPwCandidates,method="spearman")
options(scipen = 2)
ConvergentIntergenesBpNap=addSmallLegend(ggplot(data=DNAperSpecie.m,aes(y=100*ConvergentIntergeneProportion,x=PercentNAPwCandidates))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+xlab("NAP abundancy (% of proteome)")+ylab("Convergent intergenes proportion (%)")+scale_color_manual(values=col_vector[10:30]))+scale_y_log10()


ggsave(ConvergentIntergenesBpNap,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_ConvergentIntergeneProportions_vsPFAMNAPwcandidates.pdf",width=5,height = 5)


#Very different : divergent intergens have 
#Divergent Intergenes vs NAP
Correlationtest=cor.test(DNAperSpecie.m$DivergentIntergeneProportion,DNAperSpecie.m$PercentNAPwCandidates,method="spearman")
options(scipen = 2)
DivergentIntergenesBpNap=addSmallLegend(ggplot(data=DNAperSpecie.m,aes(y=100*DivergentIntergeneProportion,x=PercentNAPwCandidates))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+xlab("NAP abundancy (% of proteome)")+ylab("Divergent intergenes proportion (%)")+scale_color_manual(values=col_vector[10:30]))
ggsave(DivergentIntergenesBpNap,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_Proportion_Divergent_Intergenes_vsPFAMNAPwcandidates.pdf",width=12,height = 6)







#Percent TSS vs temperature
Correlationtest=cor.test(DNAperSpecie.m$PercentTSS,DNAperSpecie.m$OGT_Manual,method="spearman")
options(scipen = 2)
PercentTSSOGT_Manual=addSmallLegend(ggplot(data=DNAperSpecie.m,aes(y=PercentTSS,x=OGT_Manual))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+xlab("Growth temperature (°C)")+ylab("Transcriptional Units / Total Genes (%)")+scale_color_manual(values=col_vector[10:30]))

ggsave(PercentTSSOGT_Manual,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentTSS_OperonMapper_vsOGT_Manual.pdf",width=4,height = 4)



#Histones vs temp , Alba vs temps :

Correlationtest=cor.test(DNAperSpecie.m$CBFD_NFYB_HMF,DNAperSpecie.m$OGT_Manual,method="spearman")
options(scipen = 2)
CBFD_NFYB_HMFOGT_Manual=addSmallLegend(ggplot(data=DNAperSpecie.m,aes(y=CBFD_NFYB_HMF,x=OGT_Manual))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+xlab("Growth temperature (°C)")+ylab("Histones (% of proteome)")+scale_color_manual(values=col_vector[10:30]))

ggsave(CBFD_NFYB_HMFOGT_Manual,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_CBFD_NFYB_HMF_vsOGT_Manual.pdf",width=4,height = 4)



Correlationtest=cor.test(DNAperSpecie.m$Alba,DNAperSpecie.m$OGT_Manual,method="spearman")
options(scipen = 2)
AlbaOGT_Manual=addSmallLegend(ggplot(data=DNAperSpecie.m,aes(y=Alba,x=OGT_Manual))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+xlab("Growth temperature (°C)")+ylab("Alba (% of proteome)")+scale_color_manual(values=col_vector[10:30]))

ggsave(AlbaOGT_Manual,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_Alba_vsOGT_Manual.pdf",width=4,height = 4)











#Growth temp vs GC content : 

library(ppcor)
Correlationtest=pcor.test(DNAperSpecie.m$GC.content.avg,DNAperSpecie.m$OGT_Manual,DNAperSpecie.m$PercentNAPwCandidates,method="spearman")
options(scipen = 2)
PercentGCOGT_Manual=addSmallLegend(ggplot(data=DNAperSpecie.m,aes(y=GC.content.avg,x=OGT_Manual))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("partial rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+xlab("Growth temperature (°C)")+ylab("GC (%)")+scale_color_manual(values=col_vector[10:30]))

ggsave(PercentGCOGT_Manual,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentGCOGT_Manual.pdf",width=4,height = 4)


Correlationtest=pcor.test(DNAperSpecie.m$PurineBias,DNAperSpecie.m$PercentNAPwCandidates,DNAperSpecie.m$OGT_Manual,method="spearman")
options(scipen = 2)
PercentpurineOGT_Manual=addSmallLegend(ggplot(data=DNAperSpecie.m,aes(y=100*PurineBias,x=PercentNAPwCandidates))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("partial rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+xlab("NAPs (% of proteome)")+ylab("Purine Bias (%)")+scale_color_manual(values=col_vector[10:30]))

ggsave(PercentpurineOGT_Manual,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_Purine_bias_NAPs_color_Growth_temperature.pdf",width=4,height = 4)




#NAPs vs temperature
Correlationtest=cor.test(DNAperSpecie.m$PercentNAP,DNAperSpecie.m$OGT_Manual,method="spearman")
options(scipen = 2)
PercentNAPOGT_Manual=ggplot(data=DNAperSpecie.m,aes(y=PercentNAP,x=OGT_Manual))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+xlab("Growth temperature (°C)")+ylab("NAPs (% of measured proteome)")

ggsave(PercentNAPOGT_Manual,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAP_vsOGT_Manual.pdf",width=4,height = 4)


#NAPs vs temperature, including all candidates that are in operons
Correlationtest=cor.test(DNAperSpecie.m$PercentNAPwCandidateswOperon,DNAperSpecie.m$OGT_Manual,method="spearman")
options(scipen = 2)
PercentNAPOperonOGT_Manual=ggplot(data=DNAperSpecie.m,aes(y=PercentNAPwCandidateswOperon,x=OGT_Manual))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+xlab("Growth temperature (°C)")+ylab("NAPs (% of measured proteome)")

ggsave(PercentNAPOperonOGT_Manual,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAPwCandidates_includingOperons_vsOGT_Manual.pdf",width=4,height = 4)

#NAPs vs Size
Correlationtest=cor.test(DNAperSpecie.m$PercentNAP,DNAperSpecie.m$Size,method="spearman")
options(scipen = 2)
PercentNAPOGT_Manual=ggplot(data=DNAperSpecie.m,aes(y=PercentNAP,x=Size*1e-6))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+xlab("Genome Size (Mbp)")+ylab("NAPs (% of measured proteome)")
ggsave(PercentNAPOGT_Manual,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAP_vsSize.pdf",width=4,height = 4)

#NAPs vs TSS
Correlationtest=cor.test(DNAperSpecie.m$PercentNAP,DNAperSpecie.m$PercentTSS,method="spearman")
options(scipen = 2)
PercentNAPTSS=ggplot(data=DNAperSpecie.m,aes(y=PercentNAP,x=PercentTSS))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+ylab("NAPs (% of measured proteome)")+xlab("Transcriptional Units / Genes (%)")

ggsave(PercentNAPTSS,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAP_vsPercentTSS.pdf",width=4,height = 4)


#NAPs vs GC.content
Correlationtest=cor.test(DNAperSpecie.m$PercentNAP,DNAperSpecie.m$GC.content.avg,method="spearman")
options(scipen = 2)
PercentNAPGC=ggplot(data=DNAperSpecie.m,aes(x=PercentNAP,y=100*GC.content.avg))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+xlab("NAPs (% of measured proteome)")+ylab("GC content (%)")

ggsave(PercentNAPGC,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAP_vsGC.pdf",width=4,height = 4)


#NAPs vs GC.content
Correlationtest=cor.test(DNAperSpecie.m$PercentNAP,DNAperSpecie.m$GC.content.avg,method="spearman")
options(scipen = 2)
PercentNAPGC=addSmallLegend(ggplot(data=DNAperSpecie.m,aes(x=PercentNAP,y=100*GC.content.avg,color=OGT_Manual))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1)+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_viridis_c()+xlab("NAPs (% of measured proteome)")+ylab("GC content (%)"))

ggsave(PercentNAPGC,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAP_vsGC_Size_Temp.pdf",width=4,height = 4)



#GC vs NA vs Temp
Correlationtest=pcor.test(DNAperSpecie.m$PercentNAPwCandidates,DNAperSpecie.m$GC.content.avg,DNAperSpecie.m$OGT_Manual,method="spearman")
options(scipen = 2)
PercentNAPGC=addSmallLegend(ggplot(data=DNAperSpecie.m,aes(x=PercentNAPwCandidates,y=100*GC.content.avg,color=OGT_Manual))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1)+ggtitle(paste("partial rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_viridis_c()+xlab("NAPs (% of measured proteome)")+ylab("GC content (%)"))

ggsave(PercentNAPGC,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAPwCandidates_vsGC_Size_Temp.pdf",width=4,height = 4)



#NAPs vs GC.content
Correlationtest=cor.test(DNAperSpecie.m$PercentNAP,DNAperSpecie.m$NbNAP,method="spearman")
options(scipen = 2)
PercentNAP_vs_NAPNb=ggplot(data=DNAperSpecie.m,aes(x=PercentNAP,y=NbNAP))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[1:20])+xlab("NAPs (% of measured proteome)")+ylab("Number of NAP genes")

ggsave(PercentNAP_vs_NAPNb,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAP_vs_NAPNb.pdf",width=4,height = 4)



#Idem on NAP + candidates
#NAPs vs temperature
Correlationtest=cor.test(DNAperSpecie.m$PercentNAPwCandidates,DNAperSpecie.m$OGT_Manual,method="spearman")

PercentNAPOGT_Manual=ggplot(data=DNAperSpecie.m,aes(y=PercentNAPwCandidates,x=OGT_Manual))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+ylab("NAPs (% of measured proteome)")+xlab("Growth temperature (°C)")

ggsave(PercentNAPOGT_Manual,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAPwCandidates_vsOGT_Manual.pdf",width=4,height = 4)



Correlationtest=cor.test(DNAperSpecie.m$PercentNAP,DNAperSpecie.m$OGT_Manual,method="spearman")

PercentNAPOGT_Manual=ggplot(data=DNAperSpecie.m,aes(y=PercentNAP,x=OGT_Manual))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+ylab("NAPs (% of measured proteome)")+xlab("Growth temperature (°C)")
ggsave(PercentNAPOGT_Manual,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAP_vs_DSMZ_GrowthTemp.pdf",width=4,height = 4)


Correlationtest=cor.test(DNAperSpecie.m$PercentNAPwCandidates,DNAperSpecie.m$OGT_Manual,method="spearman")
PercentNAPwCandidatesOGT_Manual=ggplot(data=DNAperSpecie.m,aes(y=PercentNAPwCandidates,x=OGT_Manual))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+ylab("NAPs (% of measured proteome)")+xlab("Growth temperature (°C)")
ggsave(PercentNAPwCandidatesOGT_Manual,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAPwCandidates_vs_DSMZ_GrowthTemp.pdf",width=4,height = 4)




library(ppcor)
cor.test(DNAperSpecie.m$PercentNAPwCandidates,DNAperSpecie.m$doubling_h,method="spearman")
cor.test(DNAperSpecie.m$PercentNAPwCandidates,DNAperSpecie.m$optimum_ph,method="spearman")
pcor.test(DNAperSpecie.m$PercentNAPwCandidates,DNAperSpecie.m$optimum_ph,DNAperSpecie.m$OGT_Manual,method="spearman")
pcor.test(DNAperSpecie.m$PercentNAPwCandidates,DNAperSpecie.m$OGT_Manual,DNAperSpecie.m$optimum_ph,method="spearman")




#Same with Kodakarensis at 65 and 85°C : 
DNAperSpecie.restricted=DNAperSpecie.m[which(names(DNAperSpecie.m)%in%c("PercentNAPwCandidates","OGT_Manual","phylum","SpecieShort"))]
#Load Kodakarensis Data at 65°C : 
KodakNAPSu=read.table("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/MASS_SPEC_Muller2020/SanChenetal2020_KodakarensisSumNapTemperature.txt",header=T,sep="\t",stringsAsFactors = F)
DNAperSpecie.restricted=rbind(DNAperSpecie.restricted,c("65","Euryarchaeota",KodakNAPSu[which(KodakNAPSu$Specie=="Methanococcus Kodakarensis 65"),]$SumNap,"T.kodakarensis"))

DNAperSpecie.restricted$PercentNAPwCandidates=as.numeric(DNAperSpecie.restricted$PercentNAPwCandidates)
DNAperSpecie.restricted$OGT_Manual=as.numeric(DNAperSpecie.restricted$OGT_Manual)
#Idem on NAP + candidates
#NAPs vs temperature
Correlationtest=cor.test(DNAperSpecie.restricted$PercentNAPwCandidates,DNAperSpecie.restricted$OGT_Manual,method="spearman")

PercentNAPOGT_Manual=ggplot(data=DNAperSpecie.restricted,aes(y=PercentNAPwCandidates,x=OGT_Manual))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+ylab("NAPs (% of measured proteome)")+xlab("Growth temperature (°C)")

ggsave(PercentNAPOGT_Manual,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAPwCandidates_vsOGT_Manual_with_Kodak.pdf",width=4,height = 4)







#NAPs vs Size
Correlationtest=cor.test(DNAperSpecie.m$PercentNAPwCandidates,DNAperSpecie.m$Size,method="spearman")
options(scipen = 2)
PercentNAPSize=ggplot(data=DNAperSpecie.m,aes(y=PercentNAPwCandidates,x=Size*1e-6))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+xlab("Genome Size (Mbp)")+ylab("NAPs (% of measured proteome)")
ggsave(PercentNAPSize,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAPwCandidates_vsSize.pdf",width=4,height = 4)


#NAPs vs TSS
Correlationtest=cor.test(DNAperSpecie.m$PercentNAPwCandidates,DNAperSpecie.m$PercentTSS,method="spearman")
options(scipen = 2)
PercentNAPTSS=ggplot(data=DNAperSpecie.m,aes(y=PercentNAPwCandidates,x=PercentTSS))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+ylab("NAPs (% of measured proteome)")+xlab("Transcriptional Units / Genes (%)")
ggsave(PercentNAPTSS,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAPwCandidates_vsPercentTSS.pdf",width=4,height = 4)


#NAPs vs GC.content
Correlationtest=cor.test(DNAperSpecie.m$PercentNAPwCandidates,DNAperSpecie.m$GC.content.avg,method="spearman")
options(scipen = 2)
PercentNAPGC=ggplot(data=DNAperSpecie.m,aes(x=PercentNAPwCandidates,y=100*GC.content.avg))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+xlab("NAPs (% of measured proteome)")+ylab("GC content (%)")

ggsave(PercentNAPGC,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAPwcandidate_vsGC.pdf",width=4,height = 4)

#TFss vs GC.content
Correlationtest=cor.test(DNAperSpecie.m$TotalTF_no_outliers,DNAperSpecie.m$GC.content.avg,method="spearman")
options(scipen = 2)
PercentTFGC=ggplot(data=DNAperSpecie.m,aes(x=TotalTF_no_outliers,y=100*GC.content.avg))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+xlab("TFs (% of measured proteome)")+ylab("GC content (%)")

ggsave(PercentTFGC,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentTFnooutier_vsGC.pdf",width=4,height = 4)




#NAPs vs TF
Correlationtest=cor.test(DNAperSpecie.m$PercentNAPwCandidates,DNAperSpecie.m$TotalTF_no_outliers,method="spearman")
options(scipen = 2)
PercentNAPTF=ggplot(data=DNAperSpecie.m,aes(x=PercentNAPwCandidates,y=TotalTF_no_outliers))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+xlab("NAPs (% of measured proteome)")+ylab("TF (% of measured proteome)")

ggsave(PercentNAPTF,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAPwcandidate_vsTFnooutlier.pdf",width=4,height = 4)







#Computing another metric :
DNAperSpecie.m$PercentNAP_per_bp=DNAperSpecie.m$PercentNAP/DNAperSpecie.m$Size



#NAPS per Bp vs Size
Correlationtest=cor.test(DNAperSpecie.m$PercentNAP_per_bp,DNAperSpecie.m$PercentTSS,method="spearman")
options(scipen = 2)
PercentNAPperBp_vs_Percent_TSS=ggplot(data=DNAperSpecie.m,aes(y=PercentNAP_per_bp,x=PercentTSS))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+xlab("% of TSS")

ggsave(PercentNAPperBp_vs_Percent_TSS,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAPperBp_vs_Percent_TSS.pdf",width=4,height = 4)



#NAPS per Bp vs growth temp
Correlationtest=cor.test(DNAperSpecie.m$PercentNAP_per_bp,DNAperSpecie.m$OGT_Manual,method="spearman")
options(scipen = 2)
PercentNAPperBp_vs_growthTemp=ggplot(data=DNAperSpecie.m,aes(y=PercentNAP_per_bp,x=OGT_Manual))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+xlab("Growth temperature (°C)")

ggsave(PercentNAPperBp_vs_growthTemp,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAPperBp_vs_growthTemp.pdf",width=4,height = 4)





PercentNAPperBp_vs_PercentNAP=ggplot(data=DNAperSpecie.m,aes(y=scale(PercentNAP_per_bp),x=scale(PercentNAP)))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+scale_color_manual(values=col_vector[10:30])+ggtitle("zscore vs zscore")

ggsave(PercentNAPperBp_vs_PercentNAP,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAPperBp_vs_PercentNAP_zscore.pdf",width=4,height = 4)





#NAPs vs doubling
Correlationtest=cor.test(DNAperSpecie.m$PercentNAP,DNAperSpecie.m$doubling_h,method="spearman")
options(scipen = 2)
PercentNAPdoubling_h=ggplot(data=DNAperSpecie.m,aes(y=PercentNAP,x=doubling_h))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])

ggsave(PercentNAPdoubling_h,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAP_vsdoubling_h.pdf",width=4,height = 4)




#NAPs vs doubling
Correlationtest=cor.test(DNAperSpecie.m$PercentNAPwCandidates,DNAperSpecie.m$doubling_h,method="spearman")
options(scipen = 2)
PercentNAPdoubling_h=ggplot(data=DNAperSpecie.m,aes(y=PercentNAPwCandidates,x=doubling_h))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])

ggsave(PercentNAPdoubling_h,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAPwCandidates_vsdoubling_h.pdf",width=4,height = 4)




#NAPs vs pH
#Without candidates : 
#NAPs vs pH
Correlationtest=cor.test(DNAperSpecie.m$PercentNAP,DNAperSpecie.m$optimum_ph,method="spearman")
options(scipen = 2)
PercentNAPoptimum_ph=ggplot(data=DNAperSpecie.m,aes(y=PercentNAP,x=optimum_ph))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])

ggsave(PercentNAPoptimum_ph,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAPwCandidates_vsoptimum_ph.pdf",width=4,height = 4)

#With candidates
Correlationtest=cor.test(DNAperSpecie.m$PercentNAPwCandidates,DNAperSpecie.m$optimum_ph,method="spearman")
options(scipen = 2)
PercentNAPdoubling_h=ggplot(data=DNAperSpecie.m,aes(y=PercentNAPwCandidates,x=optimum_ph))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])

ggsave(PercentNAPoptimum_ph,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAPwCandidates_vsoptimum_ph.pdf",width=4,height = 4)





#Trancription factors : 


Correlationtest=cor.test(DNAperSpecie.m$TotalTF_percentage,DNAperSpecie.m$OGT_Manual,method="spearman")
TFOGT_Manual=ggplot(data=DNAperSpecie.m,aes(y=TotalTF_percentage,x=OGT_Manual))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+xlab("Growth temperature (°C)")
ggsave(TFOGT_Manual,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_TFvsOGT_Manual.pdf",width=4,height = 4)

Correlationtest=cor.test(DNAperSpecie.m$TotalTF_percentage,DNAperSpecie.m$PercentTSS,method="spearman")
TFvsPercentTSS=ggplot(data=DNAperSpecie.m,aes(y=TotalTF_percentage,x=PercentTSS))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+xlab("TSS (%)")
ggsave(TFvsPercentTSS,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_TFvsPercentTSS.pdf",width=4,height = 4)

Correlationtest=cor.test(DNAperSpecie.m$TotalTF_percentage,DNAperSpecie.m$Size,method="spearman")
TFvsGenomeSize=ggplot(data=DNAperSpecie.m,aes(y=TotalTF_percentage,x=Size))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+xlab("Genome Size (Mbp)")
ggsave(TFvsGenomeSize,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_TFvsGenomeSize.pdf",width=4,height = 4)

Correlationtest=cor.test(DNAperSpecie.m$NbTFmeasured,DNAperSpecie.m$Size,method="spearman")
NbTFmeasuredvsGenomeSize=ggplot(data=DNAperSpecie.m,aes(y=NbTFmeasured,x=Size/1e6))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+xlab("Genome Size (Mbp)")
ggsave(NbTFmeasuredvsGenomeSize,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_NbTFmeasuredvsGenomeSize.pdf",width=4,height = 4)


#Same calculation, without including TF that have been classified as NAP candidates: 
Correlationtest=cor.test(DNAperSpecie.m$TotalTF_no_outliers,DNAperSpecie.m$OGT_Manual,method="spearman")
TotalTF_no_outliersvs_growthTemp=ggplot(data=DNAperSpecie.m,aes(y=TotalTF_no_outliers,x=OGT_Manual))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+xlab("Growth temperature (°C)")+ylab("TFs (% measured proteome)")
ggsave(TotalTF_no_outliersvs_growthTemp,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_TotalTF_no_outliersvs_growthTemp.pdf",width=4,height = 4)

Correlationtest=cor.test(DNAperSpecie.m$TotalTF_no_outliers,DNAperSpecie.m$doubling_h,method="spearman")
TotalTF_no_outliersvs_doubling_h=ggplot(data=DNAperSpecie.m,aes(y=TotalTF_no_outliers,x=doubling_h))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+xlab("Growth temperature (°C)")+ylab("TFs (% measured proteome)")
ggsave(TotalTF_no_outliersvs_doubling_h,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_TotalTF_no_outliersvs_doublingtime.pdf",width=4,height = 4)

Correlationtest=cor.test(DNAperSpecie.m$TotalTF_no_outliers,DNAperSpecie.m$PercentTSS,method="spearman")
TotalTF_no_outliersvs_PercentTSS=ggplot(data=DNAperSpecie.m,aes(y=TotalTF_no_outliers,x=PercentTSS))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+xlab("Fraction of genes with a TSS (%)")+ylab("TFs (% measured proteome)")
ggsave(TotalTF_no_outliersvs_PercentTSS,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_TotalTF_no_outliersvs_PercentTSS.pdf",width=4,height = 4)

Correlationtest=cor.test(DNAperSpecie.m$TotalTF_no_outliers,DNAperSpecie.m$Size,method="spearman")
TotalTF_no_outliersvs_Size=ggplot(data=DNAperSpecie.m,aes(y=TotalTF_no_outliers,x=Size/1e6))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+xlab("Genome Size (Mbp)")+ylab("TFs (% measured proteome)")
ggsave(TotalTF_no_outliersvs_Size,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_TotalTF_no_outliersvs_Size.pdf",width=4,height = 4)


Correlationtest=cor.test(DNAperSpecie.m$TotalTF_no_outliers,DNAperSpecie.m$PercentNAP,method="spearman")
TotalTF_no_outliersvs_PercentNAP=ggplot(data=DNAperSpecie.m,aes(y=TotalTF_no_outliers,x=PercentNAP))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+ylab("TFs (% measured proteome)")
ggsave(TotalTF_no_outliersvs_PercentNAP,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_TotalTF_no_outliersvs_PercentNAP.pdf",width=4,height = 4)


#NAPs vs tRNA Nb : 

DNAperSpecie.m$tRNAProp=DNAperSpecie.m$tRNA/DNAperSpecie$GeneNb
Correlationtest=cor.test(DNAperSpecie.m$tRNAProp,DNAperSpecie.m$PercentNAP,method="spearman")
tRNAvs_PercentNAP=ggplot(data=DNAperSpecie.m,aes(y=tRNA/GeneNb,x=PercentNAPwCandidates))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])
ggsave(tRNAvs_PercentNAP,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_tRNAvs_vs_PercentNAP.pdf",width=4,height = 4)


#Normalizing NAP by TF abundancy, to circumvent the total proteome argument: 

DNAperSpecie.m$RatioNAP_TF=DNAperSpecie.m$PercentNAP/DNAperSpecie.m$TotalTF_percentage

DNAperSpecie.m$RatioNAP_TF_w_candidates=DNAperSpecie.m$PercentNAPwCandidates/DNAperSpecie.m$TotalTF_no_outliers



Correlationtest=cor.test(DNAperSpecie.m$RatioNAP_TF,DNAperSpecie.m$OGT_Manual,method="spearman")

NAPTF=ggplot(data=DNAperSpecie.m,aes(y=OGT_Manual,x=RatioNAP_TF))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])
ggsave(NAPTF,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentTF_vs_PercentNAP_vs_growth_Temp.pdf",width=4,height = 4)


Correlationtest=cor.test(DNAperSpecie.m$RatioNAP_TF_w_candidates,DNAperSpecie.m$OGT_Manual,method="spearman")

NAPTFcandidates=ggplot(data=DNAperSpecie.m,aes(y=OGT_Manual,x=RatioNAP_TF_w_candidates))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])
ggsave(NAPTFcandidates,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentTF_vs_PercentNAP_allcandidates_included_vs_growth_Temp.pdf",width=4,height = 4)

DNAperSpecie.m$RatioNAP_w_candidates_RNApol_Rpb1=DNAperSpecie.m$PercentNAPwCandidates/DNAperSpecie.m$RNA_pol_Rpb1_3


Correlationtest=cor.test(DNAperSpecie.m$RatioNAP_w_candidates_RNApol_Rpb1,DNAperSpecie.m$OGT_Manual,method="spearman")

NAPcandidatesPol=ggplot(data=DNAperSpecie.m,aes(y=OGT_Manual,x=RatioNAP_w_candidates_RNApol_Rpb1))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])

ggsave(NAPcandidatesPol,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAP_allcandidates_included_vs_RNA_PolRpb1_vs_growth_Temp.pdf",width=4,height = 4)

#Same per Bp : 
DNAperSpecie.m$RatioNAP_w_candidates_RNApol_Rpb1Size=DNAperSpecie.m$RatioNAP_w_candidates_RNApol_Rpb1/DNAperSpecie.m$Size

Correlationtest=cor.test(DNAperSpecie.m$RatioNAP_w_candidates_RNApol_Rpb1Size,DNAperSpecie.m$OGT_Manual,method="spearman")

NAPcandidatesPolSize=ggplot(data=DNAperSpecie.m,aes(y=OGT_Manual,x=RatioNAP_w_candidates_RNApol_Rpb1Size))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])

ggsave(NAPcandidatesPolSize,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAP_allcandidates_included_vs_RNA_PolRpb1Size_vs_growth_Temp.pdf",width=4,height = 4)








Correlationtest=cor.test(DNAperSpecie.m$PercentNAPwCandidates,DNAperSpecie.m$GeneNb,method="spearman")

NAP_geneNb=ggplot(data=DNAperSpecie.m,aes(y=GeneNb,x=PercentNAPwCandidates))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])
ggsave(NAP_geneNb,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_PercentNAP_allcandidates_included_vs_GeneNb.pdf",width=4,height = 4)


ggplot(data=DNAperSpecie.m,aes(y=OGT_Manual,x=TotalProt))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])


Correlationtest=cor.test(DNAperSpecie.m$PercentNAPwCandidates,DNAperSpecie.m$PurineBias,method="spearman")

ggplot(data=DNAperSpecie.m,aes(y=OGT_Manual,x=PurineBias))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])


###########
#Using the tRNA metric : 
DNAperSpecie.m$tRNA.syntAll=DNAperSpecie.m$tRNA.synt_1+DNAperSpecie.m$tRNA.synt_2
DNAperSpecie.m$PercentNAPwCandidatestRNA=DNAperSpecie.m$PercentNAPwCandidates/(DNAperSpecie.m$tRNA.synt_1+DNAperSpecie.m$tRNA.synt_2)
DNAperSpecie.m$PercentNAPwCandidatestRNASize=DNAperSpecie.m$PercentNAPwCandidatestRNA/DNAperSpecie.m$Size





Correlationtest=cor.test(DNAperSpecie.m$OGT_Manual,DNAperSpecie.m$tRNA.syntAll,method="spearman")

tRNATemp=ggplot(data=DNAperSpecie.m,aes(x=tRNA.syntAll,y=OGT_Manual))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+xlab("tRNA synthetases (% proteome)")+ylab("Growth temperature (°C)")
ggsave(tRNATemp,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_tRNAsynthetase_vsTemp.pdf",width=4,height = 4)


Correlationtest=cor.test(DNAperSpecie.m$tRNA.syntAll,DNAperSpecie.m$Size,method="spearman")

tRNASize=ggplot(data=DNAperSpecie.m,aes(x=tRNA.syntAll,y=Size/1000000))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+xlab("tRNA synthetases (% proteome)")+ylab("Size (Mbp)")
ggsave(tRNASize,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_tRNAsynthetase_vs_Size.pdf",width=4,height = 4)


Correlationtest=cor.test(DNAperSpecie.m$tRNA.syntAll,DNAperSpecie.m$TotalProt,method="spearman")

tRNATotalProt=ggplot(data=DNAperSpecie.m,aes(x=tRNA.syntAll,y=TotalProt))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+xlab("tRNA synthetases (% proteome)")+ylab("TotalProt")
ggsave(tRNATotalProt,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_tRNAsynthetase_vs_TotalProt.pdf",width=4,height = 4)


# PerNAP
Correlationtest=cor.test(DNAperSpecie.m$PercentNAPwCandidatestRNA,DNAperSpecie.m$OGT_Manual,method="spearman")

tRNANAPTemp=ggplot(data=DNAperSpecie.m,aes(x=PercentNAPwCandidatestRNA,y=OGT_Manual))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+xlab("NAP abundancy, tRNA synthetase normalized")+ylab("Growth temperature (°C)")
ggsave(tRNANAPTemp,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_NAPstRNAsynthetaseNormalized_vsTemp.pdf",width=4,height = 4)



Correlationtest=cor.test(DNAperSpecie.m$PercentNAPwCandidatestRNA,DNAperSpecie.m$PercentTSS,method="spearman")

tRNANAPTSS=ggplot(data=DNAperSpecie.m,aes(x=PercentNAPwCandidatestRNA,y=PercentTSS))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+xlab("NAP abundancy, tRNA synthetase normalized")+ylab("Percent of TSS")
ggsave(tRNANAPTSS,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_NAPstRNAsynthetaseNormalized_vsPercentTSS.pdf",width=4,height = 4)



Correlationtest=cor.test(DNAperSpecie.m$PercentNAPwCandidatestRNA,DNAperSpecie.m$CompactionBp,method="spearman")

tRNANAPCompactionBp=ggplot(data=DNAperSpecie.m,aes(x=PercentNAPwCandidatestRNA,y=100*CompactionBp))+geom_point()+ggrepel::geom_text_repel(aes(label=SpecieShort),size=2.25,color="grey65",fontface = "italic")+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")+ggtitle(paste("rho =", round(Correlationtest$estimate,3), " p-value :", round(Correlationtest$p.value,digits = 6),sep=""))+scale_color_manual(values=col_vector[10:30])+xlab("NAP abundancy, tRNA synthetase normalized")+ylab("fraction of coding intergenes (%)")
ggsave(tRNATemp,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_NAPstRNAsynthetaseNormalized_vsCompactionBp.pdf",width=4,height = 4)

DNAperSpecie.m$RiboL=rowSums(DNAperSpecie.m[,grep("Ribosomal_L",names(DNAperSpecie.m))])
DNAperSpecie.m$RiboS=rowSums(DNAperSpecie.m[,grep("Ribosomal_S",names(DNAperSpecie.m))])

ggplot(data=DNAperSpecie.m)+geom_point(aes(x=Spt5.NGN,y=PercentNAPwCandidates))
ggplot(data=DNAperSpecie.m)+geom_point(aes(x=Ribosomal_L5e,y=PercentNAPwCandidates))

ggplot(data=DNAperSpecie.m)+geom_point(aes(x=RtcB,y=PercentNAPwCandidates))
ggplot(data=DNAperSpecie.m)+geom_point(aes(x=Proteasome,y=PercentNAPwCandidates))

ggplot(data=DNAperSpecie.m)+geom_point(aes(x=eIF.5a,y=IF.2))

ggplot(data=DNAperSpecie.m)+geom_point(aes(x=eIF.5a,y=IF.2))
ggplot(data=DNAperSpecie.m)+geom_point(aes(x=eIF.5a-IF.2,y=PercentNAPwCandidates))

cor.test(DNAperSpecie.m$TFIIB,DNAperSpecie.m$PercentNAPwCandidates,method="spearman")
cor.test(DNAperSpecie.m$TFIIS_C,DNAperSpecie.m$PercentNAPwCandidates,method="spearman")
cor.test(DNAperSpecie.m$TFIIE_alpha,DNAperSpecie.m$PercentNAPwCandidates,method="spearman")
cor.test(DNAperSpecie.m$TFIID.18kDa,DNAperSpecie.m$PercentNAPwCandidates,method="spearman")
cor.test(DNAperSpecie.m$Proteasome,DNAperSpecie.m$PercentNAPwCandidates,method="spearman")



library(ppcor)

ppcor::pcor.test(DNAperSpecie.m$PercentNAPwCandidates,DNAperSpecie.m$PercentTSS,DNAperSpecie.m$OGT_Manual,method="spearman")
ppcor::pcor.test(DNAperSpecie.m$PercentNAPwCandidates,DNAperSpecie.m$ConvergentIntergeneProportion,DNAperSpecie.m$OGT_Manual,method="spearman")




ppcor::pcor.test(DNAperSpecie.m$PercentNAP,DNAperSpecie.m$Tss.Proportion,DNAperSpecie.m$OGT_Manual,method="pearson")
ppcor::pcor.test(DNAperSpecie.m$Tss.Proportion,DNAperSpecie.m$OGT_Manual,DNAperSpecie.m$Size,method="spearman")

ppcor::pcor.test(DNAperSpecie.m$PercentNAP,DNAperSpecie.m$OGT_Manual,DNAperSpecie.m$PercentTSS,method="spearman")

ppcor::pcor.test(DNAperSpecie.m$PercentNAP,DNAperSpecie.m$PercentTSS,DNAperSpecie.m$Size,method="pearson")





ggplot(data=DNAperSpecie.m)+geom_point(aes(x=ConvergentIntergeneProportion,y=OGT_Manual))
ggplot(data=DNAperSpecie.m)+geom_point(aes(x=DivergentIntergeneProportion,y=OGT_Manual))

ggplot(data=DNAperSpecie.m)+geom_point(aes(x=ConvergentIntergeneProportion,y=OGT_Manual))












#Phylogenetic linear odel : 
library(phylolm);library(car)



#For the Gtdb tree :
PhyloCor.gtdb.gassemblies <-function(TraitTable,TreePath,PredVar,TraitVar,Binarise,ModelType,Error,Bootstraps){
  #Removing Na values if there are any
  if(length(which(is.na(TraitTable[,which(names(TraitTable)==TraitVar)])==T))>0){
    TraitTable=TraitTable[-which(is.na(TraitTable[,which(names(TraitTable)==TraitVar)])==T),]}
  #Removing duplicated assemblies from Taxid if there are any
  if(length(which(duplicated(TraitTable$Taxid)==T))>0){
    TraitTable=TraitTable[-which(duplicated(TraitTable$Taxid)==T),]}
  TraitTable$Taxid=as.character(TraitTable$Taxid)
  #Logit transform the TraitVar :
  #Details : 
  TraitTable$Tr.logit=TraitTable[,which(names(TraitTable)==TraitVar)]
  #First transform between 0 and 1 : 
  #Set min to 0
  TraitTable$Tr.logit=TraitTable$Tr.logit-min(TraitTable$Tr.logit)
  #Set max to one :
  TraitTable$Tr.logit=TraitTable$Tr.logit/max(TraitTable$Tr.logit)
  TraitTable$Tr.logit=car::logit(TraitTable$Tr.logit,adjust = .025)
  #Binarise PredVar :
  if(Binarise==1){
    TraitTable[,which(names(TraitTable)==PredVar)][which(TraitTable[,which(names(TraitTable)==PredVar)]>0)]=1}
  
  #Loading tree
  MyTree=read.tree(TreePath)
  #List of the Assemblies we will work with : 
  species <- TraitTable$Assembly
  #Prune tree : keep only leaf for which we have info : 
  #Make sure data rows align w/ node names
  pruned.Tree <- drop.tip(MyTree,MyTree$tip.label[!(MyTree$tip.label %in% species)])
  sum(pruned.Tree$edge.length==0)
  min(pruned.Tree$edge.length[pruned.Tree$edge.length>0])
  pruned.Tree$edge.length <- pruned.Tree$edge.length + 1e-30
  
  #Trimming our table to keep only genomes present in the tree :
  
  TraitTable=TraitTable[which(TraitTable$Assembly%in%pruned.Tree$tip.label),]
  
  
  #ordering TraitTable to have the same order as the tree :
  order(pruned.Tree$tip.label)
  TraitTable.o=TraitTable[order(match(TraitTable$Assembly,pruned.Tree$tip.label)),]
  #To verify that both order correspond :
  pruned.Tree$tip.label[1:10]
  TraitTable.o$Assembly[1:10]
  tail(pruned.Tree$tip.label)
  tail(TraitTable.o$Assembly)
  table(row.names(TraitTable))
  #Specifying the model :
  Model=paste("Tr.logit~",PredVar,sep="")
  
  #Testing the model :
  Uncorrected=lm(Model,data=TraitTable.o)
  summary(Uncorrected)
  boxplot(TraitTable.o$Tr.logit~TraitTable.o[,which(names(TraitTable.o)==PredVar)],outline=F)
  
  rownames(TraitTable.o)=TraitTable.o$Assembly
  #Testing the model corrected for phylogeny :
  Corrected=phylolm(Model,data=TraitTable.o,model=ModelType,phy=pruned.Tree,measurement_error = Error,boot = Bootstraps)
  #Corrected=phylolm(Model,data=TraitTable.o,model="OUfixedRoot",phy=pruned.Tree,measurement_error = T)
  
  summary(Corrected)
  return(Corrected)
}


TraitTable=DNAperSpecie.m
TreePath="/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/TREES/GTDB/ar122_reducedEBMCDB2_GenbankAssemblyNames.tree"
names(TraitTable)[2]="Assembly"
names(TraitTable)[3]="Taxid"

ModelSum=summary(PhyloCor.gtdb.gassemblies(TraitTable = TraitTable ,TreePath = TreePath,PredVar = "OGT_Manual",TraitVar = "PercentNAPwCandidates",Binarise = 0,ModelType = "BM",Error = T,Bootstraps = 10000))

ModelSum2=summary(PhyloCor.gtdb.gassemblies(TraitTable = TraitTable ,TreePath = TreePath,PredVar = "OGT_Manual",TraitVar = "PercentNAPwCandidates",Binarise = 0,ModelType = "OUrandomRoot",Error = T,Bootstraps = F))

#pH
ModelSumpH=summary(PhyloCor.gtdb.gassemblies(TraitTable = TraitTable ,TreePath = TreePath,PredVar = "optimum_ph",TraitVar = "PercentNAPwCandidates",Binarise = 0,ModelType = "BM",Error = T,Bootstraps = 10000))

ModelSum2pH=summary(PhyloCor.gtdb.gassemblies(TraitTable = TraitTable ,TreePath = TreePath,PredVar = "optimum_ph",TraitVar = "PercentNAPwCandidates",Binarise = 0,ModelType = "OUrandomRoot",Error=T,Bootstraps = F))


ModelSumdh=summary(PhyloCor.gtdb.gassemblies(TraitTable = TraitTable ,TreePath = TreePath,PredVar = "doubling_h",TraitVar = "PercentNAPwCandidates",Binarise = 0,ModelType = "BM",Error = T,Bootstraps = 10000))

ModelSum2dh=summary(PhyloCor.gtdb.gassemblies(TraitTable = TraitTable ,TreePath = TreePath,PredVar = "doubling_h",TraitVar = "PercentNAPwCandidates",Binarise = 0,ModelType = "OUrandomRoot",Error=T,Bootstraps = F))




#With DSMZ values : 
ModelSum=summary(PhyloCor.gtdb.gassemblies(TraitTable = TraitTable ,TreePath = TreePath,PredVar = "OGT_Manual",TraitVar = "PercentNAPwCandidates",Binarise = 0,ModelType = "BM",Error = T,Bootstraps = 10000))

ModelSum2=summary(PhyloCor.gtdb.gassemblies(TraitTable = TraitTable ,TreePath = TreePath,PredVar = "OGT_Manual",TraitVar = "PercentNAPwCandidates",Binarise = 0,ModelType = "OUrandomRoot",Error = T,Bootstraps = F))


cor.test(DNAperSpecie.m$PercentNAPwCandidates,DNAperSpecie.m$OGT_Manual)

#At last probing polyamines : 


DNAperSpecie.m$Specie
Polyamines$Specie

Polyamines=read.table(file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/Kneifel_Polyamines_concentrations.txt",header=T,sep="\t",quote="",stringsAsFactors = F,)
Polyamines[which(is.na(Polyamines)==T),]=0
Polyamines$TotalPolyamines=rowSums(x = Polyamines[,which(names(Polyamines)%in%c("PUT","NSPD","SPD","HSPD","NSPM","SPM"))],na.rm = T)

Polyamines.m=merge(DNAperSpecie.m,Polyamines,by="Specie")

Polyamines.m$Specie
Polyamine_plot=ggplot(data=Polyamines.m,aes(x=TotalPolyamines,y=PercentNAPwCandidates))+geom_point(shape=21,stroke=0.1,size=3,color="grey20")+ggrepel::geom_text_repel(aes(label=Specie),size=2.5,segment.color = "grey80")+theme(aspect.ratio = 1)+theme_pubr()+theme(aspect.ratio = 1,panel.grid.major = element_line(color = "#EFEEF0"))+geom_point(aes(fill=OGT_Manual),shape=21,stroke=0.8,size=2,color="grey20")+scale_fill_viridis_c(option = "B")
ggsave(plot = Polyamine_plot,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/ARCHAEA_NAPSwcandidates_vs_POLYAMINES.pdf",height=5,width=5)



#Added 26 May :
#
#
#Showing that PercentNAPwcandidates doesn't change much when normalized by tRNA : 

#Reloading DNAperspecie.m so that it fits Archtree names
DNAperSpecie.m=merge(DNAperSpecie.m1,TFPerSpecie,by="Specie")

library(ape)
ArchTree=read.tree("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/TREES/GTDB/Arch_120_pruned_to_MassSpec_dataSet.tree")
ArchTree$tip.label=gsub("_"," ",ArchTree$tip.label)
NewOrder=DNAperSpecie.m$Specie[match(ArchTree$tip.label,DNAperSpecie.m$Specie)]
NewOrder=c("Halobacterium sp",NewOrder[c(1:10)],"Picrophilus torridus",NewOrder[c(11:19)])
if(length(which(is.na(NewOrder)==T))>0){
  NewOrder=NewOrder[-which(is.na(NewOrder)==T)]}

DNAperSpecie.m$Specie=factor(DNAperSpecie.m$Specie,levels = rev(NewOrder))

DNAperSpecie.m$PercentNAPwCandidatestRNA=DNAperSpecie.m$PercentNAPwCandidates/(DNAperSpecie.m$tRNA.synt_1+DNAperSpecie.m$tRNA.synt_2)


NAPAbudancytRNASynth=ggplot(DNAperSpecie.m)+geom_bar(aes(x=Specie,y=PercentNAPwCandidatestRNA),stat = "identity", width=0.7,fill="darkblue")+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("Abundancy \n(tRNA synthetase normalized)")+ggtitle("NAP abudancy")+xlab("")

ggsave(NAPAbudancytRNASynth,filename =paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/NAPwcandidates_Normalized_bytRNAsynth.pdf",sep=""),width = 5,height=4)
