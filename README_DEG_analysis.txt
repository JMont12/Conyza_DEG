#Process of identifying differentially expressed genes in Conyza sumatrensis following treatment with 2,4-D
#The results of this pipeline are presented in the paper "LOSS OF 2,4-D TRANSLOCATION AND METABOLISM IN 2,4-D RESISTANT Conyza sumatrensis DUE TO RAPID PHYSIOLOGICAL RESPONSE"
#This file is meant to document the tools and general syntax used. Path names to files will need to be changed and some files are not stored in this github repository
#please contact me at jake.montgomery@colostate.edu with any questions
#Jake Montgomery
#02-16-2023

#raw sequence data can be accessed on through the Gene Expression Omnibus repository under accession XXX
#raw reads were trimmed using trimmomatic and a shell script was used to loop through the read files
#this script will trim reads and output them to files in the curent working directory. This should be run from a directory called trimmed_reads to contain the results
./trim_conyza_reads.sh

#use salmon installed in a conda environment to quantify the number of reads coming from each transcript
conda create -n salmon
conda activate salmon
conda install -c bioconda salmon

#build index of canadensis transcriptome for use by salmon. This should be run in the directory containing the transcriptome .fa file
salmon index -t canadensis_56884-CDS.fasta -i canadensis_transcriptome_index

#again, a shell script is used to loop through all the trimmed reads and run salmon on each sample
./quant_conyza_to_canadensis.sh

#I made a tx2gene file by opening one salmon quant.sf file and copying the list of transcript names into another column called geneid
#this file is important for importing the data into R

#I pulled out the quant.sf file from each sample result directory and renamed it with the sample name and the _quant.sf file extension.

#these raw read counts were imported into R and the Conyza_DESeq_to_canadensis_revision.R script was used to conduct each contrast and make plots
#significant DEG's from each contrast were exported to csv files 
#these csv files were filtered in excel to remove any DEG's with log2 fold change >-2 and <2

#The three_way_DEG_filter.py script was used to identify DEG's specific to the treated vs control contrast for the resistant population, 
#DEG's found in both sensitive treated vs control contrasts but not in the resistant treated vs untreated contrast, 
#and DEG's in common amongst all treated vs untreated contrasts
python ./three_way_DEG_filter.py ./RT_vs_RC_significant_DEG.csv ./ST_vs_SC_significant_DEG.csv TT_vs_TC_significant_DEG.csv ./

#protein sequence of each gene was blasted to a local copy of the uniprotKB database
blastp -query ./canadensis_56884-CDS-prot.fasta -db /path/to/uniprot.fa -evalue 0.00001 -outfmt 6 -out canadensis_protein_blast.tsv

#the results were filtered to keep the single best alignment for each gene using a custom python script
python ./blastqc_JSM.py ./canadensis_protein_blast.tsv ./canadensis_protein_blast_filtered.tsv

#The protein IDs were extracted from each alignment
cut -f 2 ./canadensis_protein_blast_filtered.tsv | cut -d '|' -f 2 > protein_IDs.txt

#this list of protein IDs was uploaded to uniprots ID mapping tool accessible at https://www.uniprot.org/id-mapping
#after the mapping is done, the results were downloaded after ensuring that GO term information is included (customize columns tab)

#I combined the mapping results and the blast results in one excel file (canadensis_protein_blast_uniprot_filtered.xlsx) and used VLOOKUP functions to assign GO terms to each DEG
#a list of DEGs can be copied onto the Assign_GO sheet and GO terms can be copied from the GO column back to the original file.

#GO terms were counted from each contrast and this data was used to build table 1 in the paper
sed 's/\[/\n/g' genes_in_common.csv | sed 's/\]/\n/g' | grep 'GO:' | sort | uniq -c | sort -n 
sed 's/\[/\n/g' genes_specific_to_sensitive_contrasts.csv | sed 's/\]/\n/g' | grep 'GO:' | sort | uniq -c | sort -n
sed 's/\[/\n/g' genes_specific_to_RTRC.csv | sed 's/\]/\n/g' | grep 'GO:' | sort | uniq -c | sort -n