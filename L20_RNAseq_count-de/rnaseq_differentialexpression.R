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
tx2gene <- read_tsv("/home/rstudio/data/datasets_class/reference/trascriptome/gencode.v29.transcripts_no-vers_chr21_tx2gene.txt")


###################################
#### READ LOCAL FILES IN ##########
###################################

files <- file.path("/home/rstudio/data/rnaseq_exercise/reads/", paste0(dataset$sample,".quant"), "quant.sf")
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

pdf("plots_maplot.pdf")
plotMA(res, ylim=c(-3,3))
dev.off()

pdf("plots_dispersion.pdf")
plotDispEsts(dds)
dev.off()

pdf("plots_counts.pdf")
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
dev.off()

###################################
## WRITE RESULTS OF ANALYSIS #####
###################################

resdata <- as_tibble(resOrdered)
resdata$gene <- rownames(resOrdered)
write_tsv(resdata, "analysis_results.tsv")

save.image("deseq2_analysis.RData")