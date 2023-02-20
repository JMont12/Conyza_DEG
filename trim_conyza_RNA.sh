#!/bin/bash

for folder in /Users/jake/Documents/PhD/Conyza_sumatrensis/raw_reads/*; do

name="$(basename $folder)"
echo "working on $name"
cd $folder
java -jar ~/apps/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 7 ${folder}/${name}_1.fq.gz ${folder}/${name}_2.fq.gz ${name}/${name}_1_paired.fq.gz ${name}/${name}_1_unpaired.fq.gz ${name}/${name}_2_paired.fq.gz ${name}/${name}_2_unpaired.fq.gz ILLUMINACLIP:/Users/jake/apps/Trimmomatic-0.36/adapters/TruSeq2-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

done

