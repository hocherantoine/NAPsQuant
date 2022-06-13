#####################################################
## Project: diaforarchaea chromatin
## Analysis of external data from Muller et al. 2020 DOI : 
## https://doi.org/10.1038/s41586-020-2402-x
## Aim : 
## Detect known NAPs with unbiased criterion
## in all archaeal genomes with sufficient peptide ID
## Associate proteomes to our local database
## Compute all reciprocal blast matrices and link to abundancy
## Date: January 2021
## Author: Antoine Hocher
####################################################
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.

Sys.setenv(PATH="/Users/ahocher/opt/miniconda3/bin:/Users/ahocher/opt/miniconda3/condabin:/opt/local/bin:/opt/local/sbin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin")

library(rhmmer)


#Loading archaeal genomes ids :  
Genomes=read.table("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/MASS_SPEC_Muller2020/Annotated_genomes_info_to_retrieve_mass_spec.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

#Re-loading the complete table (output from ExtData_Muller_pre-processing_and_merging.R): 
MassSpecArchaea=read.table(file="/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/MASS_SPEC_Muller2020/41586_2020_2402_MOESM3_ESM_archaeaOnly_Alvus_Lumi.txt",header=T,sep="\t",quote="",stringsAsFactors = F)


#Sub-selection of correct genomes (genome with > 500 measurements, quality filterng): 
SpecieFocus=c(names(table(MassSpecArchaea$Organism)[which(table(MassSpecArchaea$Organism)>500)]),"Candidatus^Methanomethylophilus^alvus^Mx1201","Methanomassiliicoccus^luminyensis^B10")
SpecieFocus2=Genomes[which(Genomes$Specie_name_orig_dataset%in%SpecieFocus | Genomes$Specie%in%SpecieFocus),]$Specie

#Subsetting initial dataset to those species: 
MassSpecArchaea=MassSpecArchaea[which(MassSpecArchaea$Organism%in%SpecieFocus),]

#Verify 
table(MassSpecArchaea$Organism)
#Exporting a list of assemblies : 
Assemblies=Genomes[which(Genomes$Specie_name_orig_dataset%in%SpecieFocus | Genomes$Specie%in%SpecieFocus),]$Assembly


#Loading Pfam HMM search for all PFAM-A domains for those species:
#NOTE FOR GITHUB : this table is unfortunately too big for me to upload it 
#but is relatively simple to produce : 
#it's the results - tableout format- of hmmsearch of all PFAM-A versus the genomes of interest
#the only important option is the cut_ga one.
Arch_Res=as.data.frame(read_tblout("/Users/ahocher/Dropbox/Laboratory/Final_analysis_Archaeal_chromatin/HMM_RESULTS/concatArchaeaEBMC2_vs_Pfam.txt"))

#Restrict to Hmm corresponding to PFAM DNA binding , or to NAPs :
#PFAM domains corresponding to DNA binding proteins has been obtained from another pulication (see methods)
PfamDNAbinding=read.table(file="/Users/ahocher/Dropbox/Laboratory/ArchealChromatinProteom/HMMs/DBPome/PFAM_DNA_binding_list.txt",stringsAsFactors = F)$V1

#Adding TFIIE as it was absent from previous list
PfamDNAbinding=c(PfamDNAbinding,"TFIIE_alpha")

#listing all HMM models corresponding to NAPs
NAPsHMM=c("Cren7","7kD_DNA_binding","CBFD_NFYB_HMF","Alba","Bac_DNA_binding","MC1","Histone")


#Exporting all PFAM results for the sub-selection of species : 
Arch_Res$Assembly=unlist(lapply(Arch_Res$domain_name,function(x) strsplit(x,"\\$")[[1]][2]))
Arch_ResSub=Arch_Res[which(Arch_Res$Assembly%in%Assemblies),]
Arch_ResSub$ProteinID=unlist(lapply(Arch_ResSub$domain_name,function(x) paste(strsplit(strsplit(x,"\\$")[[1]][1],"@")[[1]][-1],collapse ="-")))
write.table(Arch_ResSub,file = "/Users/ahocher/Dropbox/Laboratory/Final_analysis_Archaeal_chromatin/HMM_RESULTS/Archaea_EBMCdb2_restricted_toMassSpecMuller2020.txt",row.names=F,col.names=F,sep="\t",quote=F)

#Subselect DNA binding proteins only
ArchDNAbinding=Arch_Res[which(Arch_Res$query_name%in%PfamDNAbinding),]

ArchDNAbinding$Assembly=unlist(lapply(ArchDNAbinding$domain_name,function(x) strsplit(x,"\\$")[[1]][2]))

#Further restricting tot he species for which mass spectrometry data is available : 
ArchSubSelDNA=ArchDNAbinding[which(ArchDNAbinding$Assembly%in%Assemblies),]

ArchSubSelDNA$ProteinID=unlist(lapply(ArchSubSelDNA$domain_name,function(x) paste(strsplit(strsplit(x,"\\$")[[1]][1],"@")[[1]][-1],collapse ="-")))

#Annotating all proteins containing a DNA binding domain : 
MassSpecArchaea$PFAM_DNA_binding=0
MassSpecArchaea[which(MassSpecArchaea$ID%in%ArchSubSelDNA$ProteinID),]$PFAM_DNA_binding=1

#Annotating all proteins containing a NAP binding domain : 
MassSpecArchaea$NAP=0
MassSpecArchaea[which(MassSpecArchaea$ID%in%ArchSubSelDNA[which(ArchSubSelDNA$query_name%in%NAPsHMM),]$ProteinID),]$NAP=1


#Annotating all proteins with all known Pfam Domains : 
for( PFAM in unique(Arch_ResSub$query_name)){
  if(length(which(MassSpecArchaea$ID%in%Arch_ResSub[which(Arch_ResSub$query_name==PFAM),]$ProteinID))>0){
    print(PFAM)
    MassSpecArchaea[,PFAM]=0
    MassSpecArchaea[which(MassSpecArchaea$ID%in%Arch_ResSub[which(Arch_ResSub$query_name==PFAM),]$ProteinID),PFAM]=1}
}

#Exporting data:
write.table(MassSpecArchaea,file = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/EXTERNAL_DATA/MASS_SPEC_Muller2020/41586_2020_2402_MOESM3_ESM_archaeaOnly_Alvus_Lumi_with_all_PfamAnnotation.txt",row.names=F,col.names=T,sep="\t",quote=F)
