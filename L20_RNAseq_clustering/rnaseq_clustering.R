library(DESeq2)
library(tximport)
library(tidyverse)
library(pheatmap) ### this might have to be installed
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)

setwd("/workspaces/class-rnaseq/analysis")
load("deseq2_analysis.RData")


############################################
## CLUSTERING ##############################
############################################

ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition")])

pheatmap(assay(ntd)[select,],
         cluster_cols=FALSE, annotation_col=df$condition)

plotPCA(ntd, intgroup=c("condition"))

save.image("deseq2_analysis.RData")