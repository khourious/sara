# load R packages
library(DESeq2)
library(dplyr)
library(calibrate)
library(genefilter)
library(gplots)
library(tidyr)
library(pcaExplorer)
library(PCAtools)
library(ggplot2)
library(EnhancedVolcano)
library(RColorBrewer)
library(cowplot)
library(pcaExplorer)
library(factoextra)



#Load sample annotation
chikv_metadata <- read.csv("chikv_metadata.txt", sep = "\t")
sampleTable
rownames(sampleTable) <- chikv_metadata$sampleID
sampleTable$Shera <- as.factor(chikv_metadata$Shera)
sampleTable$Crh <- as.factor(chikv_metadata$Crh)
sampleTable$DayDisease <- as.factor(chikv_metadata$DayDisease)
sampleTable$Type <- as.factor(chikv_metadata$Type)

#count table
counts <- chikv.counts
rownames(counts) <- counts$Geneid
counts <- counts[, 7:ncol(counts)]

# remove prefix and suffix from filename
colnames(counts) <- gsub("X.home.joycesilva.sara.analysis.ArboBA.align.", "", colnames(counts))
colnames(counts) <- gsub("_TM_ST_GRCh38.Aligned.out.bam", "", colnames(counts))
counts$ZK0214 <- NULL
counts$ZK0215 <- NULL
counts$ZK0219 <- NULL

sampleTable$group <- factor(paste(sampleTable$Type,sampleTable$DayDisease))

#chikv is interaction
dds <- DESeqDataSetFromMatrix(countData=counts,
                              colData=sampleTable,
                              design=~ Type + DayDisease)

#relevel to get EXP as reference
dds$Type <- relevel(dds$Type, ref="CTRL")

#Not required
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- dds[rowSums(counts(dds) > 0) > 6, ] 

#LTR test
dds <- DESeq(dds, test="LRT", full=~ Type + DayDisease, reduced=~DayDisease)

# plot dispersions
pdf("qc-dispersions.pdf", 50, 50, pointsize=100)
plotDispEsts(dds, main="dispersion plot")
dev.off()

# transform countdata to log2 scale to decrease the differences between samples with small counts
rld <- rlogTransformation(dds)

#Variance stabilizing transformation
vsd <- vst(dds, blind=FALSE)

# sample distance heatmap
pdf("qc-distance-heatmap.pdf", w=50, h=50, pointsize=100)
distsRL <- as.matrix(dist(t(assay(rld))))
hmcol <- colorRampPalette(brewer.pal(11,"RdYlGn"))(100)
rownames(distsRL) <- colnames(distsRL)
heatmap.2(distsRL,trace="none",col=rev(hmcol),margin=c(7,7),dendrogram="both",main="sample distance matrix")
dev.off()

# dot plot of the effect of the transformation to rlog
par(mfrow = c( 1, 2))
dds <- estimateSizeFactors(dds)
pdf("qc-rlog.pdf")
par(mfrow=c(1,2))
plot(log2(1+counts(dds)[,1:2]),col=rgb(0,0,0,.2),pch=16,cex=1.0,main="log2")
plot(assay(rld)[,1:2],col=rgb(0,0,0,.2),pch=16,cex=1.0,main="rlog")
dev.off()

# principal components analysis (PCA) RLD
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(rld$Type))])
rld_pca <- function (rld, intgroup = "Type", ntop = 500, colors=NULL, legendpos="bottomright", main="PCA", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop=FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")}
    else {colors = c("black", "red")}}
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)}
pdf("qc-pca-rld.pdf", 50, 50, pointsize=70)
rld_pca(rld, colors=mycols, intgroup="Type", xlim=c(-35, 45))
dev.off()

# principal components analysis (PCA) VSD
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(vsd$Type))])
vsd_pca <- function (vsd, intgroup = "Type", ntop = 500, colors=NULL, legendpos="bottomright", main="PCA", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(vsd))
  select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(vsd)[, intgroup, drop=FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")}
    else {colors = c("black", "red")}}
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)}
pdf("qc-pca-vsd.pdf", 50, 50, pointsize=80)
vsd_pca(vsd, colors=mycols, intgroup="Type", xlim=c(-35, 45))
dev.off()

# decrease the fold change noise with shrinkage function
resShrink <- lfcShrink(dds, coef="Type_CHIKV_vs_CTRL", type="normal")

# MA plot of the effect of the shrinkage correction
res <- results(dds, name= "Type_CHIKV_vs_CTRL")
resultsNames(dds)
mcols(res,use.names=TRUE)
pdf("qc-shrinkage-correction.pdf", 50, 50, pointsize=80)
par(mfrow=c(1,2))
plotMA(res,ylim=c(-7,7),main="unshurunken log2 fold change")
plotMA(resShrink,ylim=c(-7,7),main="shurunken log2 fold change")
dev.off()

# get differential expression results
table(resShrink$padj<0.05)

# order by adjusted p-value
resShrink <- resShrink[order(resShrink$padj), ]
resShrink <-na.omit(resShrink)

# merge with normalized count data
resShrinkdata <- merge(as.data.frame(resShrink),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
names(resShrinkdata)[1] <- "GeneId"
head(resShrinkdata)

# DEGs (volcano plot)
resShrink <- left_join(resShrinkdata,chikv.gs.results, by = "GeneId")
rownames(resShrink) <- resShrink$GeneSymbol
data <- resShrink
head(resShrink)
chikv.down <- data[(data$padj <= 0.05 & data$log2FoldChange <= -2),]
chikv.up <-data[(data$log2FoldChange >= 2 & data$padj <= 0.05),]
chikv.degs <- data[(data$padj <= 0.05 & data$log2FoldChange <= -2) | 
             (data$log2FoldChange >= 2 & data$padj <= 0.05),]
keyvals <- ifelse ( data$log2FoldChange < -2 & data$padj <= 0.05 , 'royalblue',
                    ifelse(data$log2FoldChange > 2 & data$padj <= 0.05, 'red', 'gray'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red'] <- 'high'
names(keyvals)[keyvals == 'gray'] <- 'mid'
names(keyvals)[keyvals == 'royalblue'] <- 'low'
pdf("diffexprShrinkage-volcanoplot.pdf", 10, 10, pointsize=20)
png("volcanoplot.png", width=3000, height=3000, units="px", pointsize=20, bg="white",res=300)
EnhancedVolcano(resShrink, lab=rownames(resShrink),title = 'CHIKV vs CTRL', x='log2FoldChange', y='padj', pCutoff=0.05, FCcutoff=2, pointSize = 4.0, labSize = 5.0) #colCustom = keyvals, xlim = c(-6, 7))
dev.off()

# write results
write.csv(resShrinkdata, file="chikv-results.csv")
write.csv(chikv.up, file="chikv-up-results.csv")
write.csv(chikv.down, file="chikv-down-results.csv")
write.csv(chikv.degs, file="chikv-degs.csv")
write.csv(data, file="chikv-results.csv")