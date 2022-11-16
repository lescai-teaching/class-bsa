library(DESeq2)
library(tximport)
library(tidyverse)
library(pheatmap)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)


load("deseq2_analysis.RData")


############################################
## CLUSTERING ##############################
############################################

ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition")])

pdf("plots_heatmap.pdf")
pheatmap(assay(ntd)[select,],
         cluster_cols=FALSE, annotation_col=df$condition)
dev.off()

pdf("plots_pca.pdf")
plotPCA(ntd, intgroup=c("condition"))
dev.off()

save.image("deseq2_analysis.RData")