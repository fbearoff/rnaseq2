library(tximport)
library(readr)
setwd("/home/frank/R_projects/rnaseq2/")
WD <- getwd()
samples <- read.table(file.path(WD, "samples.txt"), header = TRUE) #read in sample list file
files <- file.path(WD, "quants", samples$location, "quant.sf") #identify file locations
names(files) <- samples$sample #name samples 
all(file.exists(files)) # are all files present?
tx2gene <- read.csv(file.path(WD, "tx2gene_(Mus_musculus.GRCm38.cdna.all.fa.gz).csv")) #bring in tx2gene file
txi <- tximport(files, type = "salmon", tx2gene = tx2gene) # gene-level estimat
txi.tx <- tximport(files, type = "salmon", txOut = TRUE, tx2gene = tx2gene) #transcript-level estimate
txi.sum <- summarizeToGene(txi.tx, tx2gene) #summarixe transcripts to gene
all.equal(txi$counts, txi.sum$counts) #are txi and txi.sum equivalent?
write.csv(txi,"txi.csv")

#DESeq2
library(DESeq2)
sampleTable <- data.frame(condition = samples$condition) #make condition table
rownames(sampleTable) <- colnames(txi$counts) #add sample names to condition table
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition) #prepare DESeq2 file for import
dds <- DESeq(dds) #do differential expression
res <- results(dds, contrast=c("condition","10.3","FL10")) #store results of above as res, PICK CONDITIONS HERE


gene_synonym <- unique(tx2gene[, -1])
#head(rownames(res))
#head(gene_synonym$gene_id)
#match(rownames(res), gene_synonym$gene_id)
head(gene_synonym$gene_id[match(rownames(res), gene_synonym$gene_id)])

z <- data.frame(res)
z$gene_symbol <- gene_synonym$gene_symbol[match(rownames(res), gene_synonym$gene_id)]
z$gene_biotype <- gene_synonym$gene_biotype[match(rownames(res), gene_synonym$gene_id)]
z$mgi_id <- gene_synonym$mgi_id[match(rownames(res), gene_synonym$gene_id)]
z$chr <- gene_synonym$chr[match(rownames(res), gene_synonym$gene_id)]
z$start <- gene_synonym$start[match(rownames(res), gene_synonym$gene_id)]
z$end <- gene_synonym$end[match(rownames(res), gene_synonym$gene_id)]
z$description <- gene_synonym$description[match(rownames(res), gene_synonym$gene_id)]
write.csv(z, "res.csv")

# Show change in expression by mean expression level
# Would be better normalized to feature size
library(ggplot2)

p <- ggplot(z, aes(x = baseMean + 0.01, y = log2FoldChange, 
                   col = paste(padj < 0.05, is.na(padj)),
                   shape = paste(padj < 0.05, is.na(padj))))
p + geom_point() + 
  scale_x_log10() + 
  scale_y_continuous(limits = c(-2, 2)) +
  scale_color_manual(values = c("lightblue", "salmon", "green")) + 
  scale_shape_manual(values = c(0, 1, 2))


subset(z, padj < 0.05 & abs(log2FoldChange) > 1 )[ , -(3:5)]

p <- ggplot(z, aes(x = log2FoldChange, y = 1 - padj, col = abs(log2FoldChange)))
p + geom_point(size = 1) + scale_y_log10()

z[order((z$padj)), ]

top.count <- 50
top.genes <- rownames(z)[order(z$padj)][1:top.count]

abund <- txi$abundance
abund <- abund[which(rownames(abund) %in% top.genes), ]

#install.packages("gplots")
library(gplots)

#heatmap.2(abund)
#heatmap.2(log10(abund + 0.01))
#hist(log10(abund))

#abund.means <- rowMeans(abund, na.rm = TRUE)
abund.scale <- t(scale(t(abund), center = TRUE, scale = TRUE))
#heatmap.2(abund.scale)

library(RColorBrewer)

#display.brewer.all()
range(abund.scale)

heatmap.labels <- z$gene_symbol[which(rownames(z) %in% top.genes)]

no.symbol.ind <- which(is.na(heatmap.labels))

heatmap.labels[no.symbol.ind] <- z$gene_biotype[which(rownames(z) %in% top.genes)][no.symbol.ind]

samplecol <- c("forestgreen", "forestgreen", "forestgreen", "red", "red")

#par(mar=c(0.1, 0.1, 0.1, 0.1))
pdf(file = "heatmap.pdf")
heatmap.2(abund.scale, trace = "none", Colv = TRUE, margins = c(7, 10),
          breaks = seq(-2, 2, by = 4/11), col = brewer.pal(11, "BrBG"),
          labRow = heatmap.labels, cexRow = 0.8, cexCol = 1, ColSideColors = samplecol
)
dev.off()

?heatmap.2

#visualize results
library("RColorBrewer")
library("pheatmap")

rld <- rlog(dds, blind=FALSE) #get transformed values
sampleDists <- dist(t(assay(rld))) #get sample-to-sample distances
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste0(samples$sample, rld$type)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#PCA analysis
library(ggrepel)
pca_plot_data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pca_plot_data, "percentVar"))
ggplot(pca_plot_data, aes(PC1, PC2, color=condition)) + 
  geom_point(size=2) + geom_text_repel(aes(label=name)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + 
  scale_color_manual(values=c("green", "red", "purple")) #changes colors


######### THINGS TO DO ############

#head(tx2gene)

#head(txi)

#sapply(txi, dim)

# apply functions
# apply(matrix, MARGIN=1, function)
# tapply
# lapply(list, function)
# sapply
# mapply

#head(txi$abundance)
#head(txi$counts)
#head(txi$length)

#apply(head(txi$abundance, 20), MARGIN=1, FUN = mean)

#apply(txi$abundance, 2, sum)
#apply(txi$counts, 2, sum)

#head(sampleTable)
#head(dds)

#length(which(res$padj < 0.05 & abs(res$log2FoldChange) >= 1))

#bestest.ind <- which(res$padj < 0.05 & abs(res$log2FoldChange) >= 1)

#bestest.geneIDS <- rownames(res)[bestest.ind]

#unique(tx2gene$gene_symbol[which(tx2gene$gene_id %in% bestest.geneIDS)])

#which(res$padj < 0.05 & abs(res$log2FoldChange) >= 1)

#files[-6]

#head(tx2gene)

#gene2symb <- unique(tx2gene[-1])

#all(rownames(dds) == gene2symb$gene_id)