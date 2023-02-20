#!/bin/bash

for folder in ~/Documents/PhD/Conyza_sumatrensis/trimmed_reads/*; do

name="$(basename $folder)"
cd $folder
echo "working on $name"
salmon quant -i /Users/jake/Documents/PhD/Conyza_sumatrensis/Genome/Canadensis_CDS/canadensis_transcriptome_index -l A \
         -1 *_1_paired.fq.gz \
         -2 *_2_paired.fq.gz \
         -p 7 --validateMappings -o /Users/jake/Documents/PhD/Conyza_sumatrensis/salmon_quant_to_canadensis/${name}_quant

done
