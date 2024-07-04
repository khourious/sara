if(!interactive()) pdf(NULL)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("apeglm")
BiocManager::install('PCAtools')
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
denv_metadata <- read.csv("denv_metadata.txt", sep = "\t")
sampleTable
rownames(sampleTable) <- denv_metadata$sampleID
sampleTable$type <- as.factor(denv_metadata$type)
sampleTable$sorotype <- as.factor(denv_metadata$sorotype)
sampleTable$gender <- as.factor(denv_metadata$gender)

#count table
counts <- denv
rownames(counts) <- counts$Geneid
counts <- counts[, 7:ncol(counts)]

# remove prefix and suffix from filename
colnames(counts) <- gsub("X.home.joycesilva.sara.analysis.ArboBA.align.", "", colnames(counts))
colnames(counts) <- gsub("_TM_ST_GRCh38.Aligned.out.bam", "", colnames(counts))


#denv is factor
dds <- DESeqDataSetFromMatrix(countData=counts,
                              colData=sampleTable,
                              design=~gender + type)

#relevel to get control as reference
dds$type <- relevel(dds$type, ref="CTRL")

#Not required
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- dds[rowSums(counts(dds) > 0) > 6, ] 

#Very few replicates, LRT test
dds <- DESeq(dds,test="LRT",full=~gender + type, reduced=~gender)

#Variance stabilizing transformation
vsd <- vst(dds, blind=FALSE)

# transform countdata to log2 scale to decrease the differences between samples with small counts
rld <- rlog(dds, blind=FALSE)

# decrease the fold change noise with shrinkage function
resultsNames(dds)
resShrink <- lfcShrink(dds, coef="type_EXP_vs_CTRL", type = "normal")

# plot dispersions
pdf("qc-dispersions.pdf", 50, 50, pointsize=100)
plotDispEsts(dds, main="dispersion plot")
dev.off()

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

# principal components analysis (PCA) VSD
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(vsd$type))])
vsd_pca <- function (vsd, intgroup = "type", ntop = 500, colors=NULL, legendpos="bottomright", main="PCA", textcx=1, ...) {
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
pdf("qc-pca-vsd.pdf", 50, 50, pointsize=70)
vsd_pca(vsd, colors=mycols, intgroup="type", xlim=c(-35, 45))
dev.off()

# Analysis PCs
png("scree-vsd-denv.png", width=3000, height=3000, units="px", pointsize=20, bg="white",res=300)
pcaobj <- prcomp(t(assay(vsd)))
fviz_eig(pcaobj, 
         addlabels = TRUE, 
         ylim = c(0, 50),
         main="Proportion of explained proportion of variance - DENV")
dev.off()


# PCA
plotPCA(vsd, intgroup="type")

pca <- prcomp(t(assay(vsd)))
af <- cbind(pca$x[,1:3]) %>% as.data.frame()
af$PC1 <- as.numeric(af$PC1) / (pca$sdev[1] * sqrt(nrow(sampleTable)))
af$PC2 <- as.numeric(af$PC2) / (pca$sdev[2] * sqrt(nrow(sampleTable)))
af$PC3 <- as.numeric(af$PC3) / (pca$sdev[3] * sqrt(nrow(sampleTable)))
af[,"Type"] <- sampleTable$type
af$Type <- as.factor(af$Type)
af[,"SampleId"] <- sampleTable$sampleID
af$SampleId <- as.factor(af$SampleId)

# plot graphic
pca <- ggplot(af, aes(PC1, PC2, colour = `Type`)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = SampleId ), size = 2.5, box.padding = 0.5, point.padding = 0.3) +
  theme (panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(title = "PCA - DENV",
       x = "PC1 (40.3%)",
       y = "PC2 (15.6%)") 
# See
pca

# Add density curves to y and x axis
xdens <- 
  axis_canvas(pca, axis = "x") + 
  geom_density(data = af, aes(x = PC1, fill = `Type`, colour = `Type`), alpha = 0.3)
ydens <-
  axis_canvas(pca, axis = "y", coord_flip = TRUE) + 
  geom_density(data = af, aes(x = PC2, fill = `Type`, colour = `Type`), alpha = 0.3) +
  coord_flip()
pdf("PCA-DENV.pdf", w=8, h=6, pointsize=20)
pca %>%
  insert_xaxis_grob(xdens, grid::unit(1, "in"), position = "top") %>%
  insert_yaxis_grob(ydens, grid::unit(1, "in"), position = "right") %>%
  ggdraw()
dev.off()

# plot graphic
pca <- ggplot(af, aes(PC1, PC3, colour = `Type`)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = SampleId ), size = 2.5, box.padding = 0.5, point.padding = 0.3) +
  theme (panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(title = "PCA - DENV",
       x = "PC1 (40.3%)",
       y = "PC3 (13.1%)") 
# See
pca

# Add density curves to y and x axis
xdens <- 
  axis_canvas(pca, axis = "x") + 
  geom_density(data = af, aes(x = PC1, fill = `Type`, colour = `Type`), alpha = 0.3)
ydens <-
  axis_canvas(pca, axis = "y", coord_flip = TRUE) + 
  geom_density(data = af, aes(x = PC3, fill = `Type`, colour = `Type`), alpha = 0.3) +
  coord_flip()
pdf("PCA-DENV.2.pdf", w=8, h=6, pointsize=20)
pca %>%
  insert_xaxis_grob(xdens, grid::unit(1, "in"), position = "top") %>%
  insert_yaxis_grob(ydens, grid::unit(1, "in"), position = "right") %>%
  ggdraw()
dev.off()

# MA plot of the effect of the shrinkage correction
res <- results(dds, name= "type_EXP_vs_CTRL")
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

# Volcano plot
pdf("diffexprShrinkage-volcanoplot.pdf", 10, 10, pointsize=20)
EnhancedVolcano(resShrink, lab=rownames(resShrink), 
                x='log2FoldChange', y='padj', pCutoff=0.05, 
                FCcutoff=2, pointSize = 5.0, labSize = 5.0)
dev.off()

# write results
write.csv(resShrinkdata, file="denv-results.csv")

# DEGs (volcano plot)
data <- left_join(resShrinkdata,denv_gs, by = "GeneId")
rownames(data) <- data$GeneSymbol
denv.down <- data[(data$padj <= 0.05 & data$log2FoldChange <= -2),]
denv.up <-data[(data$log2FoldChange >= 2 & data$padj <= 0.05),]
denv.degs <- data[(data$padj <= 0.05 & data$log2FoldChange <= -2) | 
             (data$log2FoldChange >= 2 & data$padj <= 0.05),]
keyvals <- ifelse ( data$log2FoldChange < -2 & data$padj <= 0.05 , 'royalblue',
                    ifelse(data$log2FoldChange > 2 & data$padj <= 0.05, 'red', 'gray'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red'] <- 'Up Regulated'
names(keyvals)[keyvals == 'gray'] <- 'No Significance'
names(keyvals)[keyvals == 'royalblue'] <- 'Down Regulated'
pdf("volcanoplot.pdf", 10, 10, pointsize=20)
png("volcanoplot.png", width=3000, height=3000, units="px", pointsize=20, bg="white",res=300)
EnhancedVolcano(data, 
                lab=rownames(data),title = 'DENV vs CONTROL', 
                subtitle = ".",
                x='log2FoldChange', y='padj', pCutoff=0.05, 
                FCcutoff=2, pointSize = 4.0, labSize = 5.0, 
                colCustom = keyvals, xlim = c(-6, 6))
dev.off()


# write results
write.csv(denv.down, file="denv-down.csv")
write.csv(denv.up, file="denv-up.csv")
write.csv(denv.degs, file="denv-degs.csv")
write.csv(data, file="denv-results.csv")