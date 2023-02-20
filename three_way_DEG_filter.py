#!/usr/bin/python

#this script will take in 3 csv files of DEGs from DEseq2 and output 4 files, one for genes specific to the RTRC contrast, genes specific to STSC, genes specific to TTTC, and genes in common.
#the 3 input files are a resistant treated vs resistant control contrast followed by two independant sensitive treated vs sensitive control contrasts. 
#usage /path/to/conyza_DEG_filter.py /data/RTvRC.csv /data/STvSC.csv /path/to/TTvTC /path/to/out_dir


from sys import argv, version_info
from os.path import realpath, splitext
import os

RTRC, STSC, TTTC, out_dir = realpath(argv[1]), realpath(argv[2]), realpath(argv[3]), realpath(argv[4])

fh1=open(RTRC, 'r')
fh2=open(STSC, 'r')
fh3=open(TTTC, 'r')
line_count=0
parts=[]
genes1={}
genes2={}
genes3={}

for line in fh1:
	if line_count>0:
			line=line.strip('\n')
			line=line.strip('\r')
			parts=line.split(',')
			genes1[parts[0]]=parts
	line_count+=1
fh1.close

line_count=0

for line in fh2:
	if line_count>0:
			line=line.strip('\n')
			line=line.strip('\r')
			parts=line.split(',')
			genes2[parts[0]]=parts
	line_count+=1
fh2.close

for line in fh3:
	if line_count>0:
			line=line.strip('\n')
			line=line.strip('\r')
			parts=line.split(',')
			genes3[parts[0]]=parts
	line_count+=1
fh3.close

out1_fh=open(out_dir+"/genes_specific_to_RTRC.csv",'w+')
out1_fh.write(',baseMean,log2FoldChange,lfcSE,pvalue,padj')
out2_fh=open(out_dir+"/genes_specific_to_sensitive_contrasts.csv", 'w+')
out2_fh.write('STSC,baseMean,log2FoldChange,lfcSE,pvalue,padj,TTTC,baseMean,log2FoldChange,lfcSE,pvalue,padj')
out_common_fh=open(out_dir+"/genes_in_common.csv", 'w+')
out_common_fh.write('RTRC,baseMean,log2FoldChange,lfcSE,pvalue,padj,STSC,baseMean,log2FoldChange,lfcSE,pvalue,padj,TTTC,baseMean,log2FoldChange,lfcSE,pvalue,padj')

for entry in sorted(genes1):
	if entry in genes2 and entry in genes3:
		out_common_fh.write('\n'+str(genes1[entry][0])+','+str(genes1[entry][1])+','+str(genes1[entry][2])+','+str(genes1[entry][3])+','+str(genes1[entry][4])+','+str(genes1[entry][5])+','+str(genes2[entry][0])+','+str(genes2[entry][1])+','+str(genes2[entry][2])+','+str(genes2[entry][3])+','+str(genes2[entry][4])+','+str(genes2[entry][5])+','+str(genes3[entry][0])+','+str(genes3[entry][1])+','+str(genes3[entry][2])+','+str(genes3[entry][3])+','+str(genes3[entry][4])+','+str(genes3[entry][5]))	
	elif entry not in genes2 and entry not in genes3:
		out1_fh.write('\n'+str(genes1[entry][0])+','+str(genes1[entry][1])+','+str(genes1[entry][2])+','+str(genes1[entry][3])+','+str(genes1[entry][4])+','+str(genes1[entry][5]))

for entry in sorted(genes2):
	if entry not in genes1 and entry in genes3:
		out2_fh.write('\n'+str(genes2[entry][0])+','+str(genes2[entry][1])+','+str(genes2[entry][2])+','+str(genes2[entry][3])+','+str(genes2[entry][4])+','+str(genes2[entry][5])+','+str(genes3[entry][0])+','+str(genes3[entry][1])+','+str(genes3[entry][2])+','+str(genes3[entry][3])+','+str(genes3[entry][4])+','+str(genes3[entry][5]))

