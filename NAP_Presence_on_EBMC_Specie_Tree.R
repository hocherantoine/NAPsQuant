#####################################################################################################
## Project: Diaforarchaea chromatin analysis :                                                     #
## Plots NAP presence absence next to specie tree, iTol environment         #
## Date: Septembre 2020                                                                                 #
## Author: Antoine Hocher                                                                          #
####################################################################################################


library(ape)
source("/Users/ahocher/Dropbox/Scripts/table2itol-master/table2itol.R")
#Creating colors for later : 

library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))


#Loading the tree :

IF_Tree=read.tree(file = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PROTEIN_TREES/IF_2a_Subalignment_EBMCDb2_Only.tree")
#Loading PFam Results : 
IF_TreeName=IF_Tree$tip.label
IF_Taxid=unlist(lapply(IF_Tree$tip.label, function(x) strsplit(x,"_")[[1]][1]))

Genomes.arch=read.table(file = "/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/GENOMES_FULL_ANNOT/Archaea_EBMCdb2_complete_data.txt",header=T,sep="\t")
#add CC1 and 7kMk : 
CC1=read.table(file = "/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/HMMSEARCH/CC1/Archaea_Distribution_CC1_Proteins_with_Length.txt",header=T,sep="\t")
CC1=CC1[which(CC1$Length<200),]
DomainTable=table(CC1$Assembly)
AssemblyToFill=names(DomainTable)
Genomes.arch$CC1=0
Genomes.arch[which(Genomes.arch$Assembly%in%AssemblyToFill),"CC1"]=as.numeric(DomainTable)


#Adding a column : the Protein ID to join the table to the tree : 
Genomes.arch$TID=NA

for (i in 1 :length(IF_Taxid)){
  Genomes.arch[which(Genomes.arch$Taxid==IF_Taxid[i]),]$TID=IF_TreeName[i]
}
#In addition to NAP we plot variables associated to methanomassilicoccales : 

VarTotest=read.delim2("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/PFAM_CORRELATION/Correlation_PFam_With_Methanomassilicoccales_within_Diafor.txt",header=T,sep="\t")
VariablesToPlotBinary=c("Alba","Bac_DNA_binding","CBFD_NFYB_HMF","MC1","CC1","Alba","Hc1","Histone","Cren7","Regulator_TrmB","Lsr2","X7kD_DNA_binding",VarTotest$ID[1:10],tail(VarTotest$ID,n=10),"MutL","MutS_I","Ferritin","Histone_HNS","PadR","HMG-box","HTH_8","DUF1931","RelB","AsnC_trans_reg")
VariablesToPlot=c("assembly_level","GeneNb","GC.content.avg","tRNA","MaxQ","Size")
ExportTable=Genomes.arch[,c(which(names(Genomes.arch)=="TID"),which(names(Genomes.arch)%in%c(VariablesToPlotBinary,VariablesToPlot,"Specie")))]


#Binaryzing all data corresponding to Pfam : 
ExportTable[,which(names(ExportTable)%in%VariablesToPlotBinary)][ExportTable[,which(names(ExportTable)%in%VariablesToPlotBinary)]>1]=1

write.table(ExportTable,"/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/ITOLSPECIETREE/NAPTestinfoTable.txt",row.names=F,quote=F,sep="\t")

setwd("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/ITOLSPECIETREE/")
create_itol_files("/Users/ahocher/Dropbox/Laboratory/Diaforarchaea_Chromatin_Study/ITOLSPECIETREE/NAPTestinfoTable.txt",identifier = "TID",separator = "\t",label = "Specie")

BInaryFiles=list.files(pattern = "iTOL_binary")
                       
#Here we change all shapes to circles and color to darkblue: 
for(j in 1:length(BInaryFiles)){
i=BInaryFiles[j]
A=read.delim2(i,header=F,sep="\t",stringsAsFactors = F)
A[c(4,6,9),2]=col_vector[j]
A[c(8,11),2]=2
#A[A[,2]==0,2] = -1
write.table(A,file=i,quote = F,row.names = F,col.names = F,sep="\t")}
