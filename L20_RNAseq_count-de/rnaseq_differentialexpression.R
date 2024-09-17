library(DESeq2)
library(tximport)
library(tidyverse)
library(pheatmap)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)

###################################
## PREPARE DATASET CONDITIONS #####
###################################

setwd("/workspaces/class-rnaseq/analysis")

dataset <- tibble(
  sample = c("sample_01",
             "sample_02",
             "sample_03",
             "sample_04",
             "sample_05",
             "sample_06"),
  condition = c(rep("control", 3),
                rep("case", 3))
)
tx2gene <- read_tsv("/workspaces/class-rnaseq/datasets_reference_only/trascriptome/gencode.v29.transcripts_no-vers_chr21_tx2gene.txt")


###################################
#### READ LOCAL FILES IN ##########
###################################

files <- file.path("/workspaces/class-rnaseq/analysis/reads/", paste0(dataset$sample,".quant"), "quant.sf")
names(files) <- dataset$sample

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

colnames(txi$counts)
rownames(dataset) <- colnames(txi$counts)

dds <- DESeqDataSetFromTximport(txi, dataset, ~condition)

###################################
## PREFILTER MIN COUNTS >10 #####
###################################

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

### make sure base level is control
dds$condition <- relevel(dds$condition, ref = "control")


###################################
##### DIFFERENTIAL EXPRESSION #####
###################################

dds <- DESeq(dds)


###################################
## EXTRACT ANALYSIS RESULTS #####
###################################

res <- results(dds)
resOrdered <- res[order(res$pvalue),]

## writeLines(summary(res), "differential_expression_summary.txt")

plotMA(res, ylim=c(-3,3))

plotDispEsts(dds)

plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

###################################
## WRITE RESULTS OF ANALYSIS #####
###################################

resdata <- as_tibble(resOrdered)
resdata$gene <- rownames(resOrdered)
write_tsv(resdata, "analysis_results.tsv")

save.image("deseq2_analysis.RData")


#### on bash
# cd /workspaces/class-rnaseq/analysis
# git add deseq2_analysis.RData 
# git commit -m "saving analysis"
# git push