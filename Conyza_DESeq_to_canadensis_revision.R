#install.packages('BiocManager')
library(BiocManager)
#BiocManager::install("tximport")
library(tximport)
#install.packages("readr")
library(readr)
#BiocManager::install("DESeq2")
library(DESeq2)
#BiocManager::install('apeglm')
library(apeglm)

#import a list of trascrupt and gene names from a file
#to make this, I took the first column from a quant file and copied it to the second column as well
setwd("/Users/jake/Documents/PhD/Conyza_sumatrensis/salmon_quant_to_canadensis")
dir <- "/Users/jake/Documents/PhD/Conyza_sumatrensis/salmon_quant_to_canadensis"
tx2gene <- read.table("tx2gene.txt", header = TRUE)

#create an object that points to the quant files and give it meaningful names
files <- c(file.path(dir, "RC1_quant.sf"),
           file.path(dir, "RC2_quant.sf"),
           file.path(dir, "RC3_quant.sf"),
           file.path(dir, "RC4_quant.sf"),
           file.path(dir, "RT1_quant.sf"),
           file.path(dir, "RT2_quant.sf"),
           file.path(dir, "RT3_quant.sf"),
           file.path(dir, "SC1_quant.sf"),
           file.path(dir, "SC2_quant.sf"),
           file.path(dir, "SC3_quant.sf"),
           file.path(dir, "ST1_quant.sf"),
           file.path(dir, "ST2_quant.sf"),
           file.path(dir, "ST3_quant.sf"),
           file.path(dir, "ST4_quant.sf"),
           file.path(dir, "TC1_quant.sf"),
           file.path(dir, "TC2_quant.sf"),
           file.path(dir, "TC3_quant.sf"),
           file.path(dir, "TC4_quant.sf"),
           file.path(dir, "TT1_quant.sf"),
           file.path(dir, "TT2_quant.sf"),
           file.path(dir, "TT3_quant.sf"),
           file.path(dir, "TT4_quant.sf"))
names(files)<- c("RC1","RC2","RC3","RC4","RT1","RT2","RT3","SC1","SC2","SC3","ST1","ST2","ST3","ST4","TC1","TC2","TC3","TC4","TT1","TT2","TT3","TT4")

#import count data along with some metadata
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)
head(txi$counts)
samples <- read.csv("sample_info_revision.csv", header = TRUE)
ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)

#run the differential expression analysis and store the results for the resistant vs sensitive comparison
dds <- DESeq(ddsTxi)
dds$condition <- relevel(dds$condition, ref = "SC")
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, contrast = c("condition", "RC", "SC"))
resLFC_RC_SC <- lfcShrink(dds, coef="condition_RC_vs_SC", type="apeglm")
resLFC_ST_SC <- lfcShrink(dds, coef="condition_ST_vs_SC", type="apeglm")
dds$condition <- relevel(dds$condition, ref = "RC")
dds <- DESeq(dds)
resLFC_RT_RC <- lfcShrink(dds, coef="condition_RT_vs_RC", type="apeglm")
dds$condition <- relevel(dds$condition, ref = "TC")
dds <- DESeq(dds)
resLFC_TT_TC <- lfcShrink(dds, coef="condition_TT_vs_TC", type="apeglm")
resLFC_RC_TC <- lfcShrink(dds, coef="condition_RC_vs_TC", type="apeglm")

#change into a new directory to save the results to
setwd("/Users/jake/Documents/PhD/Conyza_sumatrensis/salmon_quant_to_canadensis/01_Revision_data")

#inspect the results of our contrast and then some summary statistics
resLFC_RC_SC
summary(resLFC_RC_SC)

#order the genes by p-value
resOrdered_RC_SC <- resLFC_RC_SC[order(resLFC_RC_SC$padj),]
head(resOrdered_RC_SC)
resOrdered_RT_RC["g29487.t1",]

#how many adjusted p-values are <0.05?
sum(resLFC_RC_SC$padj < 0.05, na.rm=TRUE)

#make a volcano plot
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT <- as.data.frame(resLFC_RC_SC)
#make the volcano plot and color significant DEGs (adj P<0.05 + LFC>2) 
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="RC vs. SC", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value), xlim=c(-10,10)))
with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
#this line will add gene names as labels, but it is way too many to show
#with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), text(log2FoldChange, -log10(padj), labels=subset(rownames(topT), topT$padj<0.05 & abs(topT$log2FoldChange)>2), cex=0.5, pos=3))
#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(0.05), col="black", lty=4, lwd=2.0)




#subset significant DEG's and output results to a csv file
resSig <- subset(resOrdered_RC_SC, padj < 0.05)
write.csv(as.data.frame(resSig), file = 'RC_vs_SC_significant_DEG.csv')

##############################################################################################
#inspect the results of our contrast and then some summary statistics
resLFC_RC_TC
summary(resLFC_RC_TC)

#order the genes by p-value
resOrdered_RC_TC <- resLFC_RC_TC[order(resLFC_RC_TC$padj),]
head(resOrdered_RC_TC)

#how many adjusted p-values are <0.05?
sum(resLFC_RC_TC$padj < 0.05, na.rm=TRUE)

#make a volcano plot
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT <- as.data.frame(resLFC_RC_TC)
#make the volcano plot and color significant DEGs (adj P<0.05 + LFC>2) 
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="RC vs. TC", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value), xlim=c(-10,10)))
with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
#this line will add gene names as labels, but it is way too many to show
#with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), text(log2FoldChange, -log10(padj), labels=subset(rownames(topT), topT$padj<0.05 & abs(topT$log2FoldChange)>2), cex=0.5, pos=3))
#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(0.05), col="black", lty=4, lwd=2.0)

#plot count data for the gene with the lowest p-value then for a specific gene
plotCounts(dds, gene=which.min(resLFC_RC_TC$padj), intgroup="condition")

#subset significant DEG's and output results to a csv file
resSig <- subset(resOrdered_RC_TC, padj < 0.05)
write.csv(as.data.frame(resSig), file = 'RC_vs_TC_significant_DEG.csv')

##############################################################################################
#inspect the results of our contrast and then some summary statistics
resLFC_ST_SC
summary(resLFC_ST_SC)

#order the genes by p-value
resOrdered_ST_SC <- resLFC_ST_SC[order(resLFC_ST_SC$padj),]
head(resOrdered_ST_SC)
resOrdered_ST_SC["g14746.t1",]
resOrdered_ST_SC["g29487.t1",]
resOrdered_ST_SC["g15488.t1",]

#how many adjusted p-values are <0.05?
sum(resLFC_ST_SC$padj < 0.05, na.rm=TRUE)

#make a volcano plot
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT <- as.data.frame(resLFC_ST_SC)
#make the volcano plot and color significant DEGs (adj P<0.05 + LFC>2) 
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="ST vs. SC", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value), xlim=c(-10,10)))
with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
#this line will add gene names as labels, but it is way too many to show
#with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), text(log2FoldChange, -log10(padj), labels=subset(rownames(topT), topT$padj<0.05 & abs(topT$log2FoldChange)>2), cex=0.5, pos=3))
#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(0.05), col="black", lty=4, lwd=2.0)

#plot count data for the gene with the lowest p-value then for a specific gene
plotCounts(dds, gene=which.min(resLFC_ST_SC$padj), intgroup="condition")

#subset significant DEG's and output results to a csv file
resSig <- subset(resOrdered_ST_SC, padj < 0.05)
write.csv(as.data.frame(resSig), file = 'ST_vs_SC_significant_DEG.csv')

##############################################################################################
#inspect the results of our contrast and then some summary statistics
resLFC_RT_RC
summary(resLFC_RT_RC)

#order the genes by p-value
resOrdered_RT_RC <- resLFC_RT_RC[order(resLFC_RT_RC$padj),]
head(resOrdered_RT_RC)

#how many adjusted p-values are <0.05?
sum(resLFC_RT_RC$padj < 0.05, na.rm=TRUE)

#make a volcano plot
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT <- as.data.frame(resLFC_RT_RC)
#make the volcano plot and color significant DEGs (adj P<0.05 + LFC>2) 
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="RT vs. RC", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value), xlim=c(-10,10)))
with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
#this line will add gene names as labels, but it is way too many to show
#with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), text(log2FoldChange, -log10(padj), labels=subset(rownames(topT), topT$padj<0.05 & abs(topT$log2FoldChange)>2), cex=0.5, pos=3))
#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(0.05), col="black", lty=4, lwd=2.0)

#plot count data for the gene with the lowest p-value then for a specific gene
plotCounts(dds, gene=which.min(resLFC_RT_RC$padj), intgroup="condition")

#subset significant DEG's and output results to a csv file
resSig <- subset(resOrdered_RT_RC, padj < 0.05)
write.csv(as.data.frame(resSig), file = 'RT_vs_RC_significant_DEG.csv')

##############################################################################################
#inspect the results of our contrast and then some summary statistics
resLFC_TT_TC
summary(resLFC_TT_TC)

#order the genes by p-value
resOrdered_TT_TC <- resLFC_TT_TC[order(resLFC_TT_TC$padj),]
head(resOrdered_TT_TC)
resOrdered_TT_TC["g14746.t1",]
resOrdered_TT_TC["g29487.t1",]
resOrdered_TT_TC["g15488.t1",]

#how many adjusted p-values are <0.05?
sum(resLFC_TT_TC$padj < 0.05, na.rm=TRUE)

#make a volcano plot
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT <- as.data.frame(resLFC_TT_TC)
#make the volcano plot and color significant DEGs (adj P<0.05 + LFC>2) 
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="TT vs. TC", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value), xlim=c(-10,10)))
with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
#this line will add gene names as labels, but it is way too many to show
#with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), text(log2FoldChange, -log10(padj), labels=subset(rownames(topT), topT$padj<0.05 & abs(topT$log2FoldChange)>2), cex=0.5, pos=3))
#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(0.05), col="black", lty=4, lwd=2.0)

#plot count data for the gene with the lowest p-value then for a specific gene
plotCounts(dds, gene=which.min(resLFC_TT_TC$padj), intgroup="condition")

#subset significant DEG's and output results to a csv file
resSig <- subset(resOrdered_TT_TC, padj < 0.05)
write.csv(as.data.frame(resSig), file = 'TT_vs_TC_significant_DEG.csv')

#transform the count data using a variance stabilizing transformation (VST) 
vsd <- vst(dds, blind=FALSE)
rlogd <- rlog(dds, blind = FALSE)
head(assay(vsd), 3)
assay(vsd)["g14746.t1",]

#plot the standard deviation of the transformed data across samples against the mean
#BiocManager::install("vsn")
library("vsn")
meanSdPlot(assay(dds))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rlogd))

#build a heatmap of the 20 genes with highest num of reads
#install.packages('pheatmap')
library(pheatmap)
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition", "pop")])
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

#calculate sample-to-sample distances to build a phylogeny
sampleDists <- dist(t(assay(dds)))
sampleDists <- dist(t(assay(vsd)))
sampleDists <- dist(t(assay(rlogd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$run
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#plot a PCA of the samples
plotPCA(vsd, intgroup="condition")
PCAdata <- plotPCA(vsd, intgroup="condition", returnData = TRUE)
plotPCA(rlogd, intgroup=c("condition", 'pop'))

library(ggplot2)
ggplot(PCAdata, aes(x=PC1, y=PC2, shape=group, color=group)) + 
  geom_point(size = 4) +
  labs(title="Principal Components of Transcript Expression Profiles \nBefore and After 2,4-D Treatment",
       x="PC1: 39% variance", y = "PC2: 25% variance") +
  scale_shape_manual(values=c(0,1,2,16, 17, 15)) +
  scale_color_manual(values=c('grey', 'pink', 'cyan', 'red', 'blue', 'black')) +
  theme_classic()
        