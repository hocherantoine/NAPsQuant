#####################################################
## Project: diaforarchaea chromatin
## Analysis of published mass spec data : 
## https://doi.org/10.1038/s41586-020-2402-x
## 
## Aim : 
## Compute 1 to 1 correlation based on reciprocal blast hits
## This scripts input is agregated data from MAssSpecPFAM correlation
## 
## NOTE GITHUB : 
## In order to reproduce the figures produced by thsi script,
## user has to download all proteomes from NCBI assembly, 
## using the accession number (GCA_xxx) available from SupTable S1
## 
## Date: January 2021
## Author: Antoine Hocher
####################################################
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.

Sys.setenv(PATH="/Users/ahocher/opt/miniconda3/bin:/Users/ahocher/opt/miniconda3/condabin:/opt/local/bin:/opt/local/sbin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin")

#Libraries : 
library(rtracklayer)
library(metablastr)
library(ape)
library(seqinr)


#Load annotated mass spec data : 
MassSpecArchaea=read.table(file="Hocher_2022_MassSpecMeasurementsAccrossSpecies_wPFAM_Annotations.txt",header=T,sep="\t",quote="",stringsAsFactors = F)[,c(1:5)]
#To reproduce figure from main manuscript one has to remove Cuniculiplasma divulgatum data, which was only used for final validation.



#Loading archaeal genomes ids :  
Genomes=read.table("Annotated_genomes_info_to_retrieve_mass_spec.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

#Adding all the external mass spec species  :
Genomes[which(Genomes$Specie=="Thermococcus^kodakarensis^KOD1"),]$To_retrieve_MassSpec=1
Genomes[which(Genomes$Specie=="Thermococcus^kodakarensis^KOD1"),]$Specie_name_orig_dataset="Thermococcus^kodakarensis^KOD1"

Genomes[which(Genomes$Specie=="Haloferax^volcanii^DS2"),]$To_retrieve_MassSpec=1
Genomes[which(Genomes$Specie=="Haloferax^volcanii^DS2"),]$Specie_name_orig_dataset="Haloferax volcanii"

Genomes[which(Genomes$Specie=="Natrialba^magadii^ATCC^43099"),]$To_retrieve_MassSpec=1
Genomes[which(Genomes$Specie=="Natrialba^magadii^ATCC^43099"),]$Specie_name_orig_dataset="Natrialba magadii"

Genomes[which(Genomes$Specie=="Nitrosopumilus^maritimus^SCM1"),]$To_retrieve_MassSpec=1
Genomes[which(Genomes$Specie=="Nitrosopumilus^maritimus^SCM1"),]$Specie_name_orig_dataset="Nitrosopumilus^maritimus^SCM1"



#Exporting a list of assemblies : 
Assemblies=Genomes[which(Genomes$To_retrieve_MassSpec==1),]$Assembly









ProteomeDir="/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/DB_Archaea95_update012020/genome_assemblies_prot_fasta/ncbi-genomes-2020-05-14/"

#For each proteome, compute the best reciprocal hits : 
setwd(ProteomeDir)
Prot_available=list.files()
ProteinToID=c()

for (Assembly in Assemblies){
  
Organism=Genomes[which(Genomes$Assembly==Assembly),]$Specie_name_orig_dataset
  Prot_to_retrieve=Prot_available[grep(Assembly,Prot_available)]
  Prot_to_retrieve=Prot_to_retrieve[grep(".faa$",Prot_to_retrieve)]
  
  
  for(Assembly2 in Assemblies){
    Organism2=Genomes[which(Genomes$Assembly==Assembly2),]$Specie_name_orig_dataset
    
    Prot_to_retrieve2=Prot_available[grep(Assembly2,Prot_available)]
    Prot_to_retrieve2=Prot_to_retrieve2[grep(".faa$",Prot_to_retrieve2)]
    
    if(Organism!=Organism2){
      if(file.exists(paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/BLAST/RECIPROCALBLASTS/",Assembly,"_vs_",Assembly2,"__",Organism,"_vs_",Organism2,"_rblast_eval_1e3.txt",sep=""))==F){
        
        
        Rhits=blast_best_reciprocal_hit(query=Prot_to_retrieve,subject = Prot_to_retrieve2,evalue = 0.001 ,search_type = "protein_to_protein",task = "blastp",is.subject.db = T,cores = 4)
        
        write.table(Rhits,paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/BLAST/RECIPROCALBLASTS/",Assembly,"_vs_",Assembly2,"__",Organism,"_vs_",Organism2,"_rblast_eval_1e3.txt",sep=""),row.names = F,sep="\t",quote=F)
      }
    }
  }
}





#Annotating best r hits result files to retrieve intensity

setwd("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/BLAST/RECIPROCALBLASTS/")
RhitsFiles=list.files(pattern="*.txt")
for( R in RhitsFiles){
  Rhits=read.table(R,header=T,sep="\t",quote="",stringsAsFactors = F)
  Rhitsm1=merge(MassSpecArchaea[,c("ID","Protein_IDs","Intensity")],Rhits,by.x="ID",by.y="query_id",all.y=T)
  names(Rhitsm1)[which(names(Rhitsm1)=="Intensity")]="Intensity_query"
  names(Rhitsm1)[which(names(Rhitsm1)=="ID")]="id_query"
  names(Rhitsm1)[which(names(Rhitsm1)=="Protein_IDs")]="Protein_IDs_query"
  
  Rhitsm2=merge(MassSpecArchaea[,c("ID","Protein_IDs","Intensity")],Rhitsm1,by.x="ID",by.y="subject_id",all.y=T)
  names(Rhitsm2)[which(names(Rhitsm2)=="Intensity")]="Intensity_subject"
  names(Rhitsm2)[which(names(Rhitsm2)=="ID")]="id_subject"
  names(Rhitsm2)[which(names(Rhitsm2)=="Protein_IDs")]="Protein_IDs_subject"
  
  write.table(Rhitsm2,paste("ANNOTATED/",  sub(".txt","_MassSpec.txt",x=R),sep=""),row.names = F,sep="\t",quote=F)
}  




#Last part : Figures and Correlations :
library(ggplot2)
library(ggpubr)


#Here we chose to only work on genomes for which more than 500 protein have a measured abundancy ( otherwise results are meaningless I'd say) : 
#So I add back genomes I chose not to include : 
OlderMassSpec=read.table("41586_2020_2402_MOESM3_ESM_archaeaOnly.txt",header=T,sep="\t",quote="")
table(OlderMassSpec$Organism)

CountTablep=as.data.frame(table(OlderMassSpec$Organism))
names(CountTablep)=c("Specie","Protein_Count")
CountTablep$Specie=as.character(CountTablep$Specie)
CountTablep=CountTablep[which(CountTablep$Protein_Count<500),]
#I initially forgot T.tenax, but it ends up not meeting my criterion.
CountTablep=rbind(CountTablep,c("Thermoproteus tenax",491))

#Plot Nb of proteins : 
CountTable=as.data.frame(table(MassSpecArchaea$Organism))
names(CountTable)=c("Specie","Protein_Count")
CountTable=rbind(CountTable,CountTablep)
CountTable$Protein_Count=as.numeric(CountTable$Protein_Count)
CountTable=CountTable[order(-CountTable$Protein_Count),]

CountTable$Specie=factor(CountTable$Specie,levels = CountTable$Specie)

BarChart=ggplot()+geom_hline(yintercept = 500,color="grey20", linetype="dashed")+geom_bar(data=CountTable,aes(x=Specie,y=Protein_Count),stat = "identity", width=0.7,fill="darkblue")+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 8),aspect.ratio=1/3)+ylab("# proteins identified")+ggtitle("Species with n > 500 only")+xlab("")


ggsave(BarChart,filename = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/All_Proteomes_MassSpec_ProteinIdentified_per_specie.pdf",width = 6,height=4)


Assemblies=Genomes[which(Genomes$Specie_name_orig_dataset%in%unique(MassSpecArchaea$Organism)),]$Assembly
#From this we create two table : one for correlation estimate and one for pvalues
CorRes=as.data.frame(matrix(NA,nrow=length(Assemblies),ncol=length(Assemblies)))
row.names(CorRes)=Assemblies
names(CorRes)=Assemblies
CorRespval=CorRes
setwd("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/BLAST/RECIPROCALBLASTS/ANNOTATED/")
RhitsFiles=list.files(pattern="*.txt")

#This bit is to cleanup data from species I did not include because data quality was insufficient.
 AssembliesNotIncluded=Genomes[which(!(Genomes$Specie_name_orig_dataset%in%unique(MassSpecArchaea$Organism)) & Genomes$To_retrieve_MassSpec==1),]$Assembly
 
 for(i in AssembliesNotIncluded){
   FilesToRemove=RhitsFiles[grep(i,RhitsFiles)]
     file.remove(FilesToRemove)
     }
RhitsFiles=list.files(pattern="*.txt")



for( R in RhitsFiles){
  BaseR=basename(R)
  AssemblyQuery=paste(strsplit(BaseR,"_")[[1]][1:2],collapse="_")
  AssemblySubject=paste(strsplit(BaseR,"_")[[1]][4:5],collapse="_")
  QuerySpecie=Genomes[which(Genomes$Assembly==AssemblyQuery),]$Specie
  SubjectSpecie=Genomes[which(Genomes$Assembly==AssemblySubject),]$Specie
  
  Rhits=read.table(R,header=T,sep="\t",quote="",stringsAsFactors = F)
  
  Correlationtest=cor.test(Rhits$Intensity_query,Rhits$Intensity_subject,use="complete.obs",method="spearman")
  
  #Adding those to the result matrix
  CorRes[AssemblyQuery,AssemblySubject]=Correlationtest$estimate
  CorRespval[AssemblyQuery,AssemblySubject]=Correlationtest$p.value
  
  #Exporting a nice plot
  PlotToExport=ggplot(data=Rhits)+geom_point(aes(x=log2(Intensity_query),y=log2(Intensity_subject)),col="grey80",size=0.8)+xlab(paste("Intensity (log2) \n",QuerySpecie))+ylab(paste("Intensity (log2) \n",SubjectSpecie))+theme_pubr()+theme(aspect.ratio = 1)+ggtitle(paste("Spearman Cor ", round(Correlationtest$estimate,2), " p-value :", round(Correlationtest$p.value,3),sep=""))
  
  if(file.exists(paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/RECIPROCALBLAST_ABUNDANCY/",sub(".txt",".pdf",BaseR),collapse = ""))==F){
    ggsave(plot=PlotToExport,filename = paste("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/RECIPROCALBLAST_ABUNDANCY/",sub(".txt",".pdf",BaseR),collapse = ""))}
}  





write.table(CorRes,"/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/MASS_SPEC_Muller2020/MassSpec_Muller_and_Inhouse_Correlation_Table_Coef.txt",row.names = T,sep='\t',quote=F)
write.table(CorRespval,"/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/MASS_SPEC_Muller2020/MassSpec_Muller_and_Inhouse_Correlation_Table_pval.txt",row.names = T,sep='\t',quote=F)

CorRes=read.table(file = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/MASS_SPEC_Muller2020/MassSpec_Muller_and_Inhouse_Correlation_Table_Coef.txt",header=T,sep="\t")


median(as.matrix(CorRes),na.rm = T)
min(as.matrix(CorRes),na.rm = T)
max(as.matrix(CorRes),na.rm = T)

library(corrplot)

GenomesSub=Genomes[which(Genomes$Assembly%in%Assemblies),]
#Replacing assembly names by Specie name :
NewOrder=match(GenomesSub$Assembly,names(CorRes))
NewOrder=NewOrder[is.na(NewOrder)==F]

names(CorRes)=GenomesSub[NewOrder,]$Specie_name_orig_dataset
row.names(CorRes)=GenomesSub[NewOrder,]$Specie_name_orig_dataset



#Re-order according to phylogeny : 

#We use a tree to order the plot : 
ArchTree=read.tree("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/TREES/GTDB/Arch_120_pruned_to_MassSpec_dataSet.tree")
ArchTree$tip.label=gsub("_"," ",ArchTree$tip.label)
NewOrder=names(CorRes)[match(ArchTree$tip.label,names(CorRes))]
NewOrder=c("Halobacterium sp",NewOrder[c(1:10)],"Picrophilus torridus",NewOrder[c(11:19)])
if(length(which(is.na(NewOrder)==T))>0){
  NewOrder=NewOrder[-which(is.na(NewOrder)==T)]}


#Re-ordering the matrix : 
CorRes=CorRes[NewOrder,NewOrder]


#Preping the names

names(CorRes)=unlist(lapply(names(CorRes),function(x) paste(substr(strsplit(x," ")[[1]][1],1,1),".",strsplit(x," ")[[1]][2],sep="")))
row.names(CorRes)=unlist(lapply(row.names(CorRes),function(x) paste(substr(strsplit(x," ")[[1]][1],1,1),".",strsplit(x," ")[[1]][2],sep="")))


CorRes[is.na(CorRes)]=1
pdf("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/Correlation_matrix_Abundancy_of_bestreciprocalBlastHits_allSpecies_wKodak_Halof_Natria.pdf",height=18,width=18)
corrplot.mixed(as.matrix(CorRes),tl.pos="lt", cl.lim = c(0,1),upper = "square")
dev.off()

pdf("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/Correlation_matrix_Abundancy_of_bestreciprocalBlastHits_allSpecies_wKodak_Halof_Natria_no_numbers.pdf",height=18,width=18)
corrplot(as.matrix(CorRes),tl.pos="lt", method="color", cl.lim = c(0,1),diag = F)
dev.off()



library(ggfortify)

CorResT=CorRes
CorResT$Specie=names(CorResT)
PCA=prcomp(CorRes)
Panel1=autoplot(PCA,data=CorResT,colour="Specie",label=T,label.size=3)+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")
Panel2=autoplot(PCA,data=CorResT,x=2,y=3,colour="Specie",label=T,label.size=3)+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")
Panel3=autoplot(PCA,data=CorResT,x=3,y=4,colour="Specie",label=T,label.size=3)+theme_pubr()+theme(aspect.ratio = 1,legend.position = "")
ggsave(ggarrange(Panel1,Panel2,Panel3),file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PLOTS/MASS_SPEC/PCA_on_Spearman_Correlation_of_best_RHits.pdf")
