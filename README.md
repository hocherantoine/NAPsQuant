Hi there,
Those scripts are relative to our quantitative analysis of NAPs in archaea. This folder contains most script and the processed data required to run it. In some instances, users will have to download publicly available data to run scripts entirely. The data folder contains the key data used in our manuscript. Raw and processed Mass spectrometry data that were produced in house are available from the PRIDE database. List of all NAPs detected as well as NCBI unique genome ID are available from Table S1 of the manuscript. 


1) ExtData_pre-processing_and_merging.R
This script run data pre-processing (mostly repartition of intensities when one entry corresponds to multiple proteins and matching between uniprot and genbank ids)
2) MassSpec_PFAM_and_NAPs_detection
This is to annotate for each protein measured by mass spec the PFAM domain that were found by hmmsearch
Subsequently, we compute a metric score for each PFAM and PFAM2GO models
3) Predict_NAP_from_MassSpecdata.R
To predict potential new NAPs based on proteins abundancy and DNA binding properties
4) OGT_vs_PFAM_NAPs_vs_Phenotype_Genome_ppties
To compute most of the correlations between OGT and PFAM, or NAPs and other variables. 

I tried to comment scripts in an intelligible manner however if that is not the case, 
I'm happy to answer any question by email (see corresponding manuscript)
